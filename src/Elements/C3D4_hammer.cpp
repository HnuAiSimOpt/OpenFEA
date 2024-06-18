/**************************************************************************

Copyright:  WH team

Author: YinJichao <jichaoyinyjc@163.com>、ChenBinXiang <chen_mech99@163.com>

Completion date:  XXX

Description: XXX

**************************************************************************/

#pragma once

#include "include/C3D4_hammer.h"

namespace CAE
{
    // 建立应变矩阵(积分点)
    double build_strain_mat_hammer(Eigen::Ref<Eigen::MatrixXd> node_coords, Eigen::Ref<Eigen::MatrixXd> strain_mat)
    {
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
        double det_jacobi = jacobi(0, 0) * jacobi(1, 1) * jacobi(2, 2) + jacobi(0, 1) * jacobi(1, 2) * jacobi(2, 0) + jacobi(0, 2) * jacobi(1, 0) * jacobi(2, 1) -
                            jacobi(0, 2) * jacobi(1, 1) * jacobi(2, 0) - jacobi(1, 2) * jacobi(2, 1) * jacobi(0, 0) - jacobi(2, 2) * jacobi(0, 1) * jacobi(1, 0);
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
        return det_jacobi;
    }

    // 建立单元刚度矩阵
    void build_ele_stiff_mat_hammer(Eigen::Ref<Eigen::MatrixXd> node_coords, Eigen::Ref<Eigen::MatrixXd> stiffness_matrix, const Matrix6d6 &C_matrix)
    {
        // 基于 单点Hammer积分 计算四面体单元
        stiffness_matrix.setZero();
        double weight = 1. / 6.; // Hammer积分 权重
        Matrix6d12 strain_mat;
        double det_jacobi = build_strain_mat_hammer(node_coords, strain_mat);
        Matrix12d6 item_temp = strain_mat.transpose() * C_matrix; // 12 x 6
        stiffness_matrix = item_temp * strain_mat;                // 12 x 12
        stiffness_matrix = weight * det_jacobi * stiffness_matrix;
    }

    // 考虑几何非线性建立应变矩阵 Green-Lagrangian应变张量
    double build_GL_strain_hamer(Eigen::Ref<Eigen::MatrixXd> node_coords, Eigen::Ref<Eigen::MatrixXd> ele_dis,
                                 Eigen::Ref<Eigen::MatrixXd> strain_L, Eigen::Ref<Eigen::MatrixXd> strain_NL, Matrix6d1 &PK2_vec, const Matrix6d6 &C_matrix)
    {
        strain_L.setZero();
        strain_NL.setZero();
        PK2_vec.setZero();
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
        // 计算雅可比行列式
        double det_jacobi = jacobi(0, 0) * jacobi(1, 1) * jacobi(2, 2) + jacobi(0, 1) * jacobi(1, 2) * jacobi(2, 0) + jacobi(0, 2) * jacobi(1, 0) * jacobi(2, 1) -
                            jacobi(0, 2) * jacobi(1, 1) * jacobi(2, 0) - jacobi(1, 2) * jacobi(2, 1) * jacobi(0, 0) - jacobi(2, 2) * jacobi(0, 1) * jacobi(1, 0);
        // 计算形变矩阵
        Matrix3d3 eye_mat(3, 3);
        eye_mat << 1., 0., 0.,
            0., 1., 0.,
            0., 0., 1.;
        Matrix3d3 F_mat = ele_dis * dN_dxyz.transpose() + eye_mat;
        // 计算 Green-Lagrangian 应变张量
        Matrix3d3 GL_mat_temp = F_mat.transpose() * F_mat - eye_mat;
        Matrix3d3 GL_mat = 0.5 * GL_mat_temp;
        Matrix6d1 Gl_vec(6, 1);
        Gl_vec << GL_mat(0, 0), GL_mat(1, 1), GL_mat(2, 2),
            2. * GL_mat(0, 1), 2. * GL_mat(1, 2), 2. * GL_mat(0, 2);
        // 计算 PK-II 应力张量
        PK2_vec = C_matrix * Gl_vec;
        // 计算位移应变转移矩阵
        for (int i = 0; i < 4; i++)
        {
            int id_1 = 3 * i;
            int id_2 = 3 * i + 1;
            int id_3 = 3 * i + 2;
            // BL
            strain_L(0, id_1) = F_mat(0, 0) * dN_dxyz(0, i);
            strain_L(0, id_2) = F_mat(1, 0) * dN_dxyz(0, i);
            strain_L(0, id_3) = F_mat(2, 0) * dN_dxyz(0, i);

            strain_L(1, id_1) = F_mat(0, 1) * dN_dxyz(1, i);
            strain_L(1, id_2) = F_mat(1, 1) * dN_dxyz(1, i);
            strain_L(1, id_3) = F_mat(2, 1) * dN_dxyz(1, i);

            strain_L(2, id_1) = F_mat(0, 2) * dN_dxyz(2, i);
            strain_L(2, id_2) = F_mat(1, 2) * dN_dxyz(2, i);
            strain_L(2, id_3) = F_mat(2, 2) * dN_dxyz(2, i);

            strain_L(3, id_1) = F_mat(0, 0) * dN_dxyz(1, i) + F_mat(0, 1) * dN_dxyz(0, i);
            strain_L(3, id_2) = F_mat(1, 0) * dN_dxyz(1, i) + F_mat(1, 1) * dN_dxyz(0, i);
            strain_L(3, id_3) = F_mat(2, 0) * dN_dxyz(1, i) + F_mat(2, 1) * dN_dxyz(0, i);

            strain_L(4, id_1) = F_mat(0, 1) * dN_dxyz(2, i) + F_mat(0, 2) * dN_dxyz(1, i);
            strain_L(4, id_2) = F_mat(1, 1) * dN_dxyz(2, i) + F_mat(1, 2) * dN_dxyz(1, i);
            strain_L(4, id_3) = F_mat(2, 1) * dN_dxyz(2, i) + F_mat(2, 2) * dN_dxyz(1, i);

            strain_L(5, id_1) = F_mat(0, 0) * dN_dxyz(2, i) + F_mat(0, 2) * dN_dxyz(0, i);
            strain_L(5, id_2) = F_mat(1, 0) * dN_dxyz(2, i) + F_mat(1, 2) * dN_dxyz(0, i);
            strain_L(5, id_3) = F_mat(2, 0) * dN_dxyz(2, i) + F_mat(2, 2) * dN_dxyz(0, i);
            // BN
            strain_NL(0, id_1) = dN_dxyz(0, i);
            strain_NL(1, id_1) = dN_dxyz(1, i);
            strain_NL(2, id_1) = dN_dxyz(2, i);

            strain_NL(3, id_2) = dN_dxyz(0, i);
            strain_NL(4, id_2) = dN_dxyz(1, i);
            strain_NL(5, id_2) = dN_dxyz(2, i);

            strain_NL(6, id_3) = dN_dxyz(0, i);
            strain_NL(7, id_3) = dN_dxyz(1, i);
            strain_NL(8, id_3) = dN_dxyz(2, i);
        };
        return det_jacobi;
    }

    // 建立切线单元刚度矩阵
    void build_ele_nl_stiff_mat_hammer(Eigen::Ref<Eigen::MatrixXd> node_coords, Eigen::Ref<Eigen::MatrixXd> ele_dis,
                                       Eigen::Ref<Eigen::MatrixXd> stiffness_matrix, Eigen::Ref<Eigen::MatrixXd> inter_force, const Matrix6d6 &C_matrix)
    {
        Matrix6d12 strain_L;
        Matrix9d12 strain_NL;
        Matrix6d1 PK2_vec;
        double weight = 1. / 6.; // Hammer积分 权重
        double det_jacobi = build_GL_strain_hamer(node_coords, ele_dis, strain_L, strain_NL, PK2_vec, C_matrix);
        Matrix9d9 SIG(9, 9);
        SIG << PK2_vec(0, 0), PK2_vec(3, 0), PK2_vec(5, 0), 0., 0., 0., 0., 0., 0.,
            PK2_vec(3, 0), PK2_vec(1, 0), PK2_vec(4, 0), 0., 0., 0., 0., 0., 0.,
            PK2_vec(5, 0), PK2_vec(4, 0), PK2_vec(2, 0), 0., 0., 0., 0., 0., 0.,
            0., 0., 0., PK2_vec(0, 0), PK2_vec(3, 0), PK2_vec(5, 0), 0., 0., 0.,
            0., 0., 0., PK2_vec(3, 0), PK2_vec(1, 0), PK2_vec(4, 0), 0., 0., 0.,
            0., 0., 0., PK2_vec(5, 0), PK2_vec(4, 0), PK2_vec(2, 0), 0., 0., 0.,
            0., 0., 0., 0., 0., 0., PK2_vec(0, 0), PK2_vec(3, 0), PK2_vec(5, 0),
            0., 0., 0., 0., 0., 0., PK2_vec(3, 0), PK2_vec(1, 0), PK2_vec(4, 0),
            0., 0., 0., 0., 0., 0., PK2_vec(5, 0), PK2_vec(4, 0), PK2_vec(2, 0);
        // 内力
        inter_force = weight * det_jacobi * (strain_L.transpose() * PK2_vec);
        // 切线刚度矩阵 线性部分+位移引起部分
        Matrix12d6 item_KL_mat_temp = strain_L.transpose() * C_matrix; // 12 x 6
        Matrix12d12 KL = item_KL_mat_temp * strain_L;                  // 12 x 12
        // 切线刚度矩阵 应力引起部分
        Matrix12d9 item_KN_mat_temp = strain_NL.transpose() * SIG; // 12 x 9
        Matrix12d12 KN = item_KN_mat_temp * strain_NL;             // 12 x 12
        stiffness_matrix = weight * det_jacobi * (KL + KN);
    }
}