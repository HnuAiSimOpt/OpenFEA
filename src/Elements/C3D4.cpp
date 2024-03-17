/**************************************************************************

Copyright:  WH team

Author: YinJichao <jichaoyinyjc@163.com>、ChenBinXiang <chen_mech99@163.com>

Completion date:  XXX

Description: XXX

**************************************************************************/

#pragma once

#include "include/C3D4.h"

namespace CAE
{
    REGISTER(ele_base, tetra_ele_elastic, "C3D4");
    // 建立本构矩阵
    void tetra_ele_elastic::build_cons_mat()
    {
        double em = matrial_struc_.young_modulus;
        double nu = matrial_struc_.poisson_ratio;
        double fac = (em * (1.0 - nu)) / ((1.0 + nu) * (1.0 - 2.0 * nu));
        double fac_a = fac * nu / (1.0 - nu);
        double fac_b = fac * (1.0 - 2.0 * nu) / (2.0 * (1.0 - nu));
        // 本构矩阵赋值
        C_matrix_ << fac, fac_a, fac_a, 0., 0., 0.,
            fac_a, fac, fac_a, 0., 0., 0.,
            fac_a, fac_a, fac, 0., 0., 0.,
            0., 0., 0., fac_b, 0., 0.,
            0., 0., 0., 0., fac_b, 0.,
            0., 0., 0., 0., 0., fac_b;
    }

    // 建立应变矩阵
    void tetra_ele_elastic::build_strain_mat(Eigen::Ref<Eigen::MatrixXd> node_coords, Eigen::Ref<Eigen::MatrixXd> strain_mat)
    {
        // TODO 后期加入关键词或者参数，控制 积分方式选择
        det_jacobi_ = build_strain_mat_hammer(node_coords, strain_mat);
    }

    // 建立单元刚度矩阵
    void tetra_ele_elastic::build_ele_stiff_mat(Eigen::Ref<Eigen::MatrixXd> node_coords, Eigen::Ref<Eigen::MatrixXd> stiffness_matrix)
    {
        build_ele_stiff_mat_hammer(node_coords, stiffness_matrix, C_matrix_);
    }

    // 建立切线单元刚度矩阵（几何非线性）
    void tetra_ele_elastic::build_ele_nl_stiff_mat(Eigen::Ref<Eigen::MatrixXd> node_coords, Eigen::Ref<Eigen::MatrixXd> node_dis,
                                                   Eigen::Ref<Eigen::MatrixXd> stiffness_matrix, Eigen::Ref<Eigen::MatrixXd> inter_force)
    {
        build_ele_nl_stiff_mat_hammer(node_coords, node_dis, stiffness_matrix, inter_force, C_matrix_);
    }

    // 建立单元质量矩阵
    void tetra_ele_elastic::build_ele_mass(const vector<int> &node_topos, const vector<vector<double>> &coords, vector<double> &Mass)
    {
        // define & initialize variables
        double volume, mass;
        Eigen::Matrix4d nodes_coor;
        nodes_coor.setZero();
        int item_node;
        // trans nodes coor to a matrix
        for (int i = 0; i < 4; i++)
        {
            item_node = node_topos[i] - 1;
            nodes_coor(i, 0) = coords[item_node][0];
            nodes_coor(i, 1) = coords[item_node][1];
            nodes_coor(i, 2) = coords[item_node][2];
            nodes_coor(i, 3) = 1.0;
        }
        // calculate the element volume & mass = volums * density
        volume = std::abs(nodes_coor.determinant() / 6.0);
        mass = volume * this->matrial_struc_.density;
        // assemble the Global mass vector
        for (int i = 0; i < 4; i++)
        {
            Mass[node_topos[i] - 1] = Mass[node_topos[i] - 1] + mass / 4.0;
        }
    }

    // 计算单元内力
    void tetra_ele_elastic::cal_in_force(const vector<int> &node_topos, const vector<vector<double>> &real_coords, const vector<double> &disp_d,
                                         vector<double> &stress, vector<double> &strain, vector<double> &InFroce)
    {
        Matrix4d3 nodes_coor;
        nodes_coor.setZero();
        int item_node;
        // trans nodes coor to a matrix
        for (int i = 0; i < 4; i++)
        {
            item_node = node_topos[i] - 1;
            nodes_coor(i, 0) = real_coords[item_node][0];
            nodes_coor(i, 1) = real_coords[item_node][1];
            nodes_coor(i, 2) = real_coords[item_node][2];
        }
        Eigen::MatrixXd stiffness_matrix, disp, dstrain, inforce;
        stiffness_matrix.resize(12, 12);
        disp.resize(12, 1);
        int idx_temp;
        for (int i = 0; i < 4; i++)
        {
            idx_temp = node_topos[i] - 1;
            disp(3 * i, 0) = disp_d[3 * idx_temp];
            disp(3 * i + 1, 0) = disp_d[3 * idx_temp + 1];
            disp(3 * i + 2, 0) = disp_d[3 * idx_temp + 2];
        }
        build_ele_stiff_mat(nodes_coor, stiffness_matrix);
        inforce = stiffness_matrix * disp; // 12*12 12*1 = 12*1
        for (int i = 0; i < 4; i++)
        {
            idx_temp = node_topos[i] - 1;
            // 内力为负，故减去
            InFroce[3 * idx_temp] -= inforce(3 * i, 0);
            InFroce[3 * idx_temp + 1] -= inforce(3 * i + 1, 0);
            InFroce[3 * idx_temp + 2] -= inforce(3 * i + 2, 0);
        }
    }

    // 计算单元时间步长
    void tetra_ele_elastic::update_timestep(Eigen::Ref<Eigen::MatrixXd> node_coords, double &time_step)
    {
        Eigen::MatrixXd stiffness_matrix;
        stiffness_matrix.resize(12, 12);
        build_ele_stiff_mat(node_coords, stiffness_matrix);
        double max_v = DBL_MIN;
        double temp;
        for (int i = 0; i < 12; i++)
        {
            temp = 0;
            for (int j = 0; j < 12; j++)
            {
                temp = temp + abs(stiffness_matrix(i, j));
            }
            if (max_v < temp)
            {
                max_v = temp;
            }
        }
        if (time_step > 2 / sqrt(max_v))
        {
            time_step = 2 / sqrt(max_v);
        }
    }
}