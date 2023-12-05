/**************************************************************************

Copyright:  WH team

Author: YinJichao <jichaoyinyjc@163.com>��ChenBinXiang <chen_mech99@163.com>

Completion date:  XXX

Description: XXX

**************************************************************************/

#pragma once

#include "C3D4.h"

namespace CAE
{
    REGISTER(ele_base, tetra_ele_elastic, "C3D4");
    // ������������
    void tetra_ele_elastic::build_cons_mat()
    {
        double em = matrial_struc_.young_modulus;
        double nu = matrial_struc_.poisson_ratio;
        double fac = (em * (1.0 - nu)) / ((1.0 + nu) * (1.0 - 2.0 * nu));
        double fac_a = fac * nu / (1.0 - nu);
        double fac_b = fac * (1.0 - 2.0 * nu) / (2.0 * (1.0 - nu));
        // ��������ֵ
        C_matrix_ << fac, fac_a, fac_a, 0., 0., 0.,
            fac_a, fac, fac_a, 0., 0., 0.,
            fac_a, fac_a, fac, 0., 0., 0.,
            0., 0., 0., fac_b, 0., 0.,
            0., 0., 0., 0., fac_b, 0.,
            0., 0., 0., 0., 0., fac_b;
    }

    // ����Ӧ�����
    //void tetra_ele_elastic::build_strain_mat(Matrix4d3 &node_coords, Matrix6d12 &strain_mat)
    void tetra_ele_elastic::build_strain_mat(Eigen::Ref<Eigen::MatrixXd> node_coords, Eigen::Ref<Eigen::MatrixXd> strain_mat) {
        strain_mat.setZero();
        Matrix3d4 dN_drst(3, 4), dN_dxyz(3, 4);
        Matrix3d3 jacobi(3, 3), inv_jacobi(3, 3);
        // �����κ����ԵȲ�����ϵ��ƫ��
        dN_drst << -1., 1., 0., 0.,
            -1., 0., 1., 0.,
            -1., 0., 0., 1.;
        jacobi = dN_drst * node_coords; // 3x3 ��3x4 by 4x3
        inv_jacobi = jacobi.inverse();  // 3x3
        // �����κ�������Ȼ����ϵ��ƫ��
        dN_dxyz = inv_jacobi * dN_drst; // 3x4
        // ����Ӧ��ת������
        det_jacobi_ = jacobi(0, 0) * jacobi(1, 1) * jacobi(2, 2) + jacobi(0, 1) * jacobi(1, 2) * jacobi(2, 0) + jacobi(0, 2) * jacobi(1, 0) * jacobi(2, 1) - jacobi(0, 2) * jacobi(1, 1) * jacobi(2, 0) - jacobi(1, 2) * jacobi(2, 1) * jacobi(0, 0) - jacobi(2, 2) * jacobi(0, 1) * jacobi(1, 0);
        for (int i = 0; i < 4; i++)
        {
            int id_1 = 3 * i;
            int id_2 = 3 * i + 1;
            int id_3 = 3 * i + 2;
            strain_mat(0, id_1) = dN_dxyz(0, i);
            strain_mat(1, id_2) = dN_dxyz(1, i);
            strain_mat(2, id_3) = dN_dxyz(2, i);
            strain_mat(3, id_1) = dN_dxyz(1, i);
            strain_mat(3, id_2) = dN_dxyz(0, i);
            strain_mat(4, id_2) = dN_dxyz(2, i);
            strain_mat(4, id_3) = dN_dxyz(1, i);
            strain_mat(5, id_1) = dN_dxyz(2, i);
            strain_mat(5, id_3) = dN_dxyz(0, i);
        };
    }

    // ������Ԫ�նȾ���
    void tetra_ele_elastic::build_ele_stiff_mat(Eigen::Ref<Eigen::MatrixXd> node_coords, Eigen::Ref<Eigen::MatrixXd> stiffness_matrix)
    {
        // ���� ����Hammer���� ���������嵥Ԫ
        stiffness_matrix.setZero();
        double weight = 1. / 6.; // Hammer���� Ȩ��
        Matrix6d12 strain_mat;
        build_strain_mat(node_coords, strain_mat);
        Matrix12d6 item_temp = strain_mat.transpose() * C_matrix_; // 12 x 6
        stiffness_matrix = item_temp * strain_mat;                        // 12 x 12
        stiffness_matrix = weight * det_jacobi_ * stiffness_matrix;
    }

    // ������Ԫ��������
    void tetra_ele_elastic::build_ele_mass(const vector<int>& node_topos, const vector<vector<double>>& coords, vector<double>& Mass)
    {
        // define & initialize variables
        double volume, mass;
        Eigen::Matrix4d nodes_coor;
        nodes_coor.setZero();
        int item_node;
        // trans nodes coor to a matrix
        for (int i = 0; i < 4; i++) {
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
        for (int i = 0; i < 4;i++) {
            Mass[node_topos[i] - 1] = Mass[node_topos[i] - 1] + mass / 4.0;
        }
    }

    // ���㵥Ԫ����
    void tetra_ele_elastic::cal_in_force(const vector<int>& node_topos, const vector<vector<double>>& real_coords, const vector<double>& disp_d,
                                         vector<double>& stress, vector<double>& strain, vector<double>& InFroce)
    {
        Matrix4d3 nodes_coor;
        nodes_coor.setZero();
        int item_node;
        // trans nodes coor to a matrix
        for (int i = 0; i < 4; i++) {
            item_node = node_topos[i] - 1;
            nodes_coor(i, 0) = real_coords[item_node][0];
            nodes_coor(i, 1) = real_coords[item_node][1];
            nodes_coor(i, 2) = real_coords[item_node][2];
        }
        Eigen::MatrixXd stiffness_matrix, disp,dstrain, inforce;
        stiffness_matrix.resize(12, 12);
        disp.resize(12, 1);
        int idx_temp;
        for (int i = 0; i < 4; i++) {
            idx_temp = node_topos[i] - 1;
            disp(3 * i, 0) = disp_d[3 * idx_temp];
            disp(3 * i + 1, 0) = disp_d[3 * idx_temp + 1];
            disp(3 * i + 2, 0) = disp_d[3 * idx_temp + 2];
        }
        build_ele_stiff_mat(nodes_coor, stiffness_matrix);
        inforce = stiffness_matrix * disp;//12*12 12*1 = 12*1
        for (int i = 0; i < 4; i++) {
            idx_temp = node_topos[i] - 1;
            // ����Ϊ�����ʼ�ȥ
            InFroce[3 * idx_temp] -= inforce(3 * i, 0);
            InFroce[3 * idx_temp + 1] -= inforce(3 * i + 1, 0);
            InFroce[3 * idx_temp + 2] -= inforce(3 * i + 2, 0);
        }
    }

    // ���㵥Ԫʱ�䲽��
    void tetra_ele_elastic::update_timestep(Eigen::Ref<Eigen::MatrixXd> node_coords, double& time_step)
    {
        Eigen::MatrixXd stiffness_matrix;
        stiffness_matrix.resize(12, 12);
        build_ele_stiff_mat(node_coords, stiffness_matrix);
        double max_v = DBL_MIN;
        double temp;
        for (int i = 0; i < 12; i++) {
            temp = 0;
            for (int j = 0; j < 12; j++) {
                temp = temp + abs(stiffness_matrix(i, j));
            }
            if (max_v < temp) {
                max_v = temp;
            }
        }
        if (time_step > 2 / sqrt(max_v)) {
            time_step = 2 / sqrt(max_v);
        }




        //Eigen::MatrixXd stiffness_matrix;
        //stiffness_matrix.resize(12,12);
        //build_ele_stiff_mat(node_coords, stiffness_matrix);
        //Eigen::EigenSolver<Eigen::MatrixXd> es(stiffness_matrix);
        //Matrix12d12 D = es.pseudoEigenvalueMatrix();
        //double max_v = DBL_MIN;
        //for (int i = 0; i < 12; i++) {
        //    if (max_v < D(i, i)) {
        //        max_v = D(i, i);
        //    }
        //}
        //if (time_step > 2 / sqrt(max_v)) {
        //    time_step = 2 / sqrt(max_v);
        //}
    }
}