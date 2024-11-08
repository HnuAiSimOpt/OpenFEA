/**************************************************************************

Copyright:  WH team

Author: YinJichao <jichaoyinyjc@163.com>、ChenBinXiang <chen_mech99@163.com>

Completion date:  XXX

Description: XXX

**************************************************************************/

#pragma once

#include "include/C3D8.h"

namespace CAE
{
    REGISTER(ele_base, hex_ele_elastic, "C3D8");
    // 建立本构矩阵
    void hex_ele_elastic::build_cons_mat()
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

    // 建立应变矩阵(积分点)
    void hex_ele_elastic::build_strain_mat(Eigen::Ref<Eigen::MatrixXd> node_coords, Eigen::Ref<Eigen::MatrixXd> strain_mat, vector<double> &gp_points, double *det_jacobi_point)
    {
        // TODO 后期加入关键词或者参数，控制 积分方式选择
        (*det_jacobi_point) = build_strain_mat_gauss(node_coords, strain_mat, gp_points);
    }

    // 建立单元刚度矩阵
    void hex_ele_elastic::build_ele_stiff_mat(Eigen::Ref<Eigen::MatrixXd> node_coords, Eigen::Ref<Eigen::MatrixXd> stiffness_matrix)
    {
        build_ele_stiff_mat_gauss(node_coords, stiffness_matrix, C_matrix_);
    }

    // 建立切线单元刚度矩阵（几何非线性）
    void hex_ele_elastic::build_ele_nl_stiff_mat(Eigen::Ref<Eigen::MatrixXd> node_coords, Eigen::Ref<Eigen::MatrixXd> node_dis,
                                                 Eigen::Ref<Eigen::MatrixXd> stiffness_matrix, Eigen::Ref<Eigen::MatrixXd> inter_force)
    {
        build_ele_nl_stiff_mat_gauss(node_coords, node_dis, stiffness_matrix, inter_force, C_matrix_);
    }

    // 建立质量矩阵
    void hex_ele_elastic::build_ele_mass(const vector<int> &node_topos, const vector<vector<double>> &coords, vector<double> &Mass)
    {
        Matrix24d24 eM;
        eM.setZero();
        Eigen::MatrixXd shape_fun;
        shape_fun.resize(3, 24);
        Matrix8d3 node_coords;
        int item_node;
        for (int i = 0; i < 8; i++)
        {
            item_node = node_topos[i] - 1;
            node_coords(i, 0) = coords[item_node][0];
            node_coords(i, 1) = coords[item_node][1];
            node_coords(i, 2) = coords[item_node][2];
        }
        // 初始化等参坐标
        double gp_values = 1. / sqrt(3.);
        Matrix8d3 gps;
        gps << -gp_values, -gp_values, -gp_values,
            gp_values, -gp_values, -gp_values,
            gp_values, gp_values, -gp_values,
            -gp_values, gp_values, -gp_values,
            -gp_values, -gp_values, gp_values,
            gp_values, -gp_values, gp_values,
            gp_values, gp_values, gp_values,
            -gp_values, gp_values, gp_values;
        //
        double weight = 1.; // 两点高斯积分 权重
        double ro = this->matrial_struc_.density;
        for (int i = 0; i < 8; i++)
        {
            vector<double> gp_points = {gps(i, 0), gps(i, 1), gps(i, 2)};
            double det_jacobi_point;
            build_shape_fun(node_coords, shape_fun, gp_points, det_jacobi_point);
            eM = eM + shape_fun.transpose() * shape_fun * det_jacobi_point * ro;
        }
        // 组装进整体质量荷载列阵，8个节点
        double temp;
        for (int i = 0; i < 8; i++)
        {
            temp = 0;
            for (int j = 0; j < 3; j++)
            {
                for (int k = 0; k < 24; k++)
                {
                    temp = temp + eM(i * 3 + j, k);
                } // 集中质量列阵
            }
            Mass[(node_topos[i] - 1)] = Mass[(node_topos[i] - 1)] + temp / 3;
        }
    }

    void hex_ele_elastic::build_shape_fun(Eigen::Ref<Eigen::MatrixXd> node_coords, Eigen::MatrixXd &shape_fun, vector<double> &gp_points, double &det_jacobi_point)
    {
        // 初始化
        shape_fun.setZero();
        Matrix3d8 dN_drst(3, 8);
        Matrix3d3 jacobi(3, 3);
        double r = gp_points[0], s = gp_points[1], t = gp_points[2];
        // 计算形函数
        vector<vector<double>> Parentnode{{-1, 1, 1, -1, -1, 1, 1, -1}, {-1, -1, 1, 1, -1, -1, 1, 1}, {-1, -1, -1, -1, 1, 1, 1, 1}};
        double temp;
        for (int i = 0; i < 8; i++)
        {
            temp = (1 + Parentnode[0][i] * r) * (1 + Parentnode[1][i] * s) * (1 + Parentnode[2][i] * t) / 8;
            shape_fun(0, 3 * i) = shape_fun(1, 3 * i + 1) = shape_fun(2, 3 * i + 2) = temp;
        }
        // 计算形函数对等参坐标系的偏导
        dN_drst(0, 0) = -0.125 * (1.0 - s) * (1.0 - t);
        dN_drst(0, 1) = 0.125 * (1.0 - s) * (1.0 - t);
        dN_drst(0, 2) = 0.125 * (1.0 + s) * (1.0 - t);
        dN_drst(0, 3) = -0.125 * (1.0 + s) * (1.0 - t);
        dN_drst(0, 4) = -0.125 * (1.0 - s) * (1.0 + t);
        dN_drst(0, 5) = 0.125 * (1.0 - s) * (1.0 + t);
        dN_drst(0, 6) = 0.125 * (1.0 + s) * (1.0 + t);
        dN_drst(0, 7) = -0.125 * (1.0 + s) * (1.0 + t);
        //
        dN_drst(1, 0) = -0.125 * (1.0 - r) * (1.0 - t);
        dN_drst(1, 1) = -0.125 * (1.0 + r) * (1.0 - t);
        dN_drst(1, 2) = 0.125 * (1.0 + r) * (1.0 - t);
        dN_drst(1, 3) = 0.125 * (1.0 - r) * (1.0 - t);
        dN_drst(1, 4) = -0.125 * (1.0 - r) * (1.0 + t);
        dN_drst(1, 5) = -0.125 * (1.0 + r) * (1.0 + t);
        dN_drst(1, 6) = 0.125 * (1.0 + r) * (1.0 + t);
        dN_drst(1, 7) = 0.125 * (1.0 - r) * (1.0 + t);
        //
        dN_drst(2, 0) = -0.125 * (1.0 - r) * (1.0 - s);
        dN_drst(2, 1) = -0.125 * (1.0 + r) * (1.0 - s);
        dN_drst(2, 2) = -0.125 * (1.0 + r) * (1.0 + s);
        dN_drst(2, 3) = -0.125 * (1.0 - r) * (1.0 + s);
        dN_drst(2, 4) = 0.125 * (1.0 - r) * (1.0 - s);
        dN_drst(2, 5) = 0.125 * (1.0 + r) * (1.0 - s);
        dN_drst(2, 6) = 0.125 * (1.0 + r) * (1.0 + s);
        dN_drst(2, 7) = 0.125 * (1.0 - r) * (1.0 + s);
        jacobi = dN_drst * node_coords; // 3x3 ：3x8 by 8x3
        // 计算应变转换矩阵
        det_jacobi_point = jacobi(0, 0) * jacobi(1, 1) * jacobi(2, 2) + jacobi(0, 1) * jacobi(1, 2) * jacobi(2, 0) + jacobi(0, 2) * jacobi(1, 0) * jacobi(2, 1) - jacobi(0, 2) * jacobi(1, 1) * jacobi(2, 0) - jacobi(1, 2) * jacobi(2, 1) * jacobi(0, 0) - jacobi(2, 2) * jacobi(0, 1) * jacobi(1, 0);
    }

    // 计算单元内力
    void hex_ele_elastic::cal_in_force(const vector<int> &node_topos, const vector<vector<double>> &real_coords, const vector<double> &disp_d, vector<double> &stress, vector<double> &strain, vector<double> &InFroce)
    {
        Matrix8d3 nodes_coor;
        nodes_coor.setZero();
        int item_node;
        // trans nodes coor to a matrix
        for (int i = 0; i < 8; i++)
        {
            item_node = node_topos[i] - 1;
            nodes_coor(i, 0) = real_coords[item_node][0];
            nodes_coor(i, 1) = real_coords[item_node][1];
            nodes_coor(i, 2) = real_coords[item_node][2];
        }
        Eigen::MatrixXd stiffness_matrix, disp, dstrain, inforce;
        stiffness_matrix.resize(24, 24);
        disp.resize(24, 1);
        int idx_temp;
        for (int i = 0; i < 8; i++)
        {
            idx_temp = node_topos[i] - 1;
            disp(3 * i, 0) = disp_d[3 * idx_temp];
            disp(3 * i + 1, 0) = disp_d[3 * idx_temp + 1];
            disp(3 * i + 2, 0) = disp_d[3 * idx_temp + 2];
        }
        build_ele_stiff_mat(nodes_coor, stiffness_matrix);
        inforce = stiffness_matrix * disp; // 24*24 24*1 = 24*1
        for (int i = 0; i < 8; i++)
        {
            idx_temp = node_topos[i] - 1;
            // 内力为负，故减去
            InFroce[3 * idx_temp] -= inforce(3 * i, 0);
            InFroce[3 * idx_temp + 1] -= inforce(3 * i + 1, 0);
            InFroce[3 * idx_temp + 2] -= inforce(3 * i + 2, 0);
        }
    }

    // 计算单元时间步长
    void hex_ele_elastic::update_timestep(Eigen::Ref<Eigen::MatrixXd> node_coords, double &time_step)
    {
        Eigen::MatrixXd stiffness_matrix;
        stiffness_matrix.resize(24, 24);
        build_ele_stiff_mat(node_coords, stiffness_matrix);
        double max_v = DBL_MIN;
        double temp;
        for (int i = 0; i < 24; i++)
        {
            temp = 0;
            for (int j = 0; j < 24; j++)
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

    // 交界面积分点物理空间坐标（非协调）
    void hex_ele_elastic::gps_phy_coords(Eigen::Ref<Eigen::MatrixXd> nodes1,
        Eigen::Ref<Eigen::MatrixXd> phy_gps, vector<double>& W_1,
        vector<Eigen::Vector3d>& Normal)
    {
        // 交界面(四边形)积分 权重为1
        double gp_values = 1. / sqrt(3.);
        Eigen::MatrixXd gps(4, 2);//高斯积分点
        gps << gp_values, gp_values,
            gp_values, -gp_values,
            -gp_values, gp_values,
            -gp_values, -gp_values;

        Eigen::MatrixXd dNdxi_1(2, 4);
        Eigen::MatrixXd N_1(1, 4);
        for (int q = 0; q < 4; q++)
        {
            dNdxi_1(0, 0) = -0.25 * (1 - gps(q, 1));
            dNdxi_1(0, 1) = 0.25 * (1 - gps(q, 1));
            dNdxi_1(0, 2) = 0.25 * (1 + gps(q, 1));
            dNdxi_1(0, 3) = -0.25 * (1 + gps(q, 1));
            dNdxi_1(1, 0) = -0.25 * (1 - gps(q, 0));
            dNdxi_1(1, 1) = -0.25 * (1 + gps(q, 0));
            dNdxi_1(1, 2) = 0.25 * (1 + gps(q, 0));
            dNdxi_1(1, 3) = 0.25 * (1 - gps(q, 0));
				
            N_1(0, 0) = 0.25 * (1 - gps(q, 0)) * (1 - gps(q, 1));
            N_1(0, 1) = 0.25 * (1 + gps(q, 0)) * (1 - gps(q, 1));
            N_1(0, 2) = 0.25 * (1 + gps(q, 0)) * (1 + gps(q, 1));
            N_1(0, 3) = 0.25 * (1 - gps(q, 0)) * (1 + gps(q, 1));

            Eigen::MatrixXd Jac = dNdxi_1 * nodes1;
            Eigen::Vector3d a1 = Jac.row(0);
            Eigen::Vector3d a2 = Jac.row(1);
            Eigen::Vector3d a3 = a1.cross(a2);
            double norm_a3 = a3.norm();
            Eigen::Vector3d unit_a3 = a3.normalized();

            phy_gps.row(q) = N_1 * nodes1;
            W_1[q]=norm_a3;
            Normal[q]=unit_a3;
            
        }

    }
    


}