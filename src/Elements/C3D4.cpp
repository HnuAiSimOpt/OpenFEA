/**************************************************************************

Copyright:  WH team

Author: YinJichao <jichaoyinyjc@163.com>、ChenBinXiang <chen_mech99@163.com>

Completion date:  XXX

Description: XXX

**************************************************************************/

#pragma once

#include "C3D4.h"

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
    //void tetra_ele_elastic::build_strain_mat(Matrix4d3 &node_coords, Matrix6d12 &strain_mat)
    void tetra_ele_elastic::build_strain_mat(Eigen::Ref<Eigen::MatrixXd> node_coords, Eigen::Ref<Eigen::MatrixXd> strain_mat) {
        strain_mat.setZero();
        Matrix3d4 dN_drst(3, 4), dN_dxyz(3, 4);
        Matrix3d3 jacobi(3, 3), inv_jacobi(3, 3);
        // 计算形函数对等参坐标系的偏导
        dN_drst << -1., 1., 0., 0.,
            -1., 0., 1., 0.,
            -1., 0., 0., 1.;
        jacobi = dN_drst * node_coords; // 3x3 ：3x4 by 4x3
        inv_jacobi = jacobi.inverse();  // 3x3
        // 计算形函数对自然坐标系得偏导
        dN_dxyz = inv_jacobi * dN_drst; // 3x4
        // 计算应变转换矩阵
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

    // 建立单元刚度矩阵
    void tetra_ele_elastic::build_ele_stiff_mat(Eigen::Ref<Eigen::MatrixXd> node_coords, Eigen::Ref<Eigen::MatrixXd> stiffness_matrix)
    {
        // 基于 单点Hammer积分 计算四面体单元
        stiffness_matrix.setZero();
        double weight = 1. / 6.; // Hammer积分 权重
        Matrix6d12 strain_mat;
        build_strain_mat(node_coords, strain_mat);
        Matrix12d6 item_temp = strain_mat.transpose() * C_matrix_; // 12 x 6
        stiffness_matrix = item_temp * strain_mat;                        // 12 x 12
        stiffness_matrix = weight * det_jacobi_ * stiffness_matrix;
    }
}