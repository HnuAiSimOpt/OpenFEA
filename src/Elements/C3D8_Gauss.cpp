/**************************************************************************

Copyright:  WH team

Author: YinJichao <jichaoyinyjc@163.com>

Completion date:  XXX

Description: XXX

**************************************************************************/

#pragma once
#include "include/C3D8_Gauss.h"

namespace CAE
{
    // 建立应变矩阵(积分点)
    double build_strain_mat_gauss(Eigen::Ref<Eigen::MatrixXd> node_coords, Eigen::Ref<Eigen::MatrixXd> strain_mat, vector<double> &gp_points)
    {
        // 初始化
        strain_mat.setZero();
        Matrix3d8 dN_drst(3, 8), dN_dxyz(3, 8);
        Matrix3d3 jacobi(3, 3), inv_jacobi(3, 3);
        double r = gp_points[0], s = gp_points[1], t = gp_points[2];
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
        inv_jacobi = jacobi.inverse();  // 3x3
        // 计算形函数对自然坐标系得偏导
        dN_dxyz = inv_jacobi * dN_drst; // 3x8
        // 计算应变转换矩阵
        double det_jacobi_point = jacobi(0, 0) * jacobi(1, 1) * jacobi(2, 2) +
                                  jacobi(0, 1) * jacobi(1, 2) * jacobi(2, 0) +
                                  jacobi(0, 2) * jacobi(1, 0) * jacobi(2, 1) -
                                  jacobi(0, 2) * jacobi(1, 1) * jacobi(2, 0) -
                                  jacobi(1, 2) * jacobi(2, 1) * jacobi(0, 0) -
                                  jacobi(2, 2) * jacobi(0, 1) * jacobi(1, 0);
        for (int i = 0; i < 8; i++)
        {
            strain_mat(0, 3 * i) = dN_dxyz(0, i);
            strain_mat(1, 3 * i + 1) = dN_dxyz(1, i);
            strain_mat(2, 3 * i + 2) = dN_dxyz(2, i);
            strain_mat(3, 3 * i) = dN_dxyz(1, i);
            strain_mat(3, 3 * i + 1) = dN_dxyz(0, i);
            strain_mat(4, 3 * i + 1) = dN_dxyz(2, i);
            strain_mat(4, 3 * i + 2) = dN_dxyz(1, i);
            strain_mat(5, 3 * i) = dN_dxyz(2, i);
            strain_mat(5, 3 * i + 2) = dN_dxyz(0, i);
        }
        return det_jacobi_point;
    }

    // 建立单元刚度矩阵
    void build_ele_stiff_mat_gauss(Eigen::Ref<Eigen::MatrixXd> node_coords, Eigen::Ref<Eigen::MatrixXd> stiffness_matrix, const Matrix6d6 &C_matrix)
    {
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
        double weight = 1.; // 两点高斯积分 权重
        stiffness_matrix.setZero();
        Matrix6d24 strain_mat;
        Matrix24d6 item_temp_1;
        Matrix24d24 item_temp_2;
        for (int i = 0; i < 8; i++)
        {
            vector<double> gp_points = {gps(i, 0), gps(i, 1), gps(i, 2)};
            double det_jacobi_point;
            det_jacobi_point = build_strain_mat_gauss(node_coords, strain_mat, gp_points);
            item_temp_1 = strain_mat.transpose() * C_matrix; // 24 x 6
            item_temp_2 = item_temp_1 * strain_mat;          // 24 x 24
            stiffness_matrix = stiffness_matrix + weight * det_jacobi_point * item_temp_2;
        }
    }

    // 考虑几何非线性建立应变矩阵 Green-Lagrangian应变张量
    double build_GL_strain_gauss(Eigen::Ref<Eigen::MatrixXd> node_coords, Eigen::Ref<Eigen::MatrixXd> ele_dis,
                                 Eigen::Ref<Eigen::MatrixXd> strain_L, Eigen::Ref<Eigen::MatrixXd> strain_NL, Matrix6d1 &PK2_vec,
                                 const Matrix6d6 &C_matrix, vector<double> &gp_points)
    {
        strain_L.setZero();
        strain_NL.setZero();
        PK2_vec.setZero();
        Matrix3d8 dN_drst, dN_dxyz;
        Matrix3d3 jacobi, inv_jacobi;
        double r = gp_points[0], s = gp_points[1], t = gp_points[2];
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
        inv_jacobi = jacobi.inverse();  // 3x3
        // 计算形函数对自然坐标系得偏导
        dN_dxyz = inv_jacobi * dN_drst; // 3x8
                                        // 计算应变转换矩阵
        double det_jacobi_point = jacobi(0, 0) * jacobi(1, 1) * jacobi(2, 2) +
                                  jacobi(0, 1) * jacobi(1, 2) * jacobi(2, 0) +
                                  jacobi(0, 2) * jacobi(1, 0) * jacobi(2, 1) -
                                  jacobi(0, 2) * jacobi(1, 1) * jacobi(2, 0) -
                                  jacobi(1, 2) * jacobi(2, 1) * jacobi(0, 0) -
                                  jacobi(2, 2) * jacobi(0, 1) * jacobi(1, 0);

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
        for (int i = 0; i < 8; i++)
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
        return det_jacobi_point;
    }

    // 建立切线单元刚度矩阵
    void build_ele_nl_stiff_mat_gauss(Eigen::Ref<Eigen::MatrixXd> node_coords, Eigen::Ref<Eigen::MatrixXd> ele_dis,
                                      Eigen::Ref<Eigen::MatrixXd> stiffness_matrix, Eigen::Ref<Eigen::MatrixXd> inter_force, const Matrix6d6 &C_matrix)
    {
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
        double weight = 1.; // 两点高斯积分 权重
        Matrix6d24 strain_L;
        Matrix9d24 strain_NL;
        Matrix6d1 PK2_vec;
        stiffness_matrix.setZero();
        inter_force.setZero();
        for (int i = 0; i < 8; i++)
        {
            vector<double> gp_points = {gps(i, 0), gps(i, 1), gps(i, 2)};
            double det_jacobi;
            det_jacobi = build_GL_strain_gauss(node_coords, ele_dis, strain_L, strain_NL, PK2_vec, C_matrix, gp_points);
            //  SIG
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
            inter_force = inter_force + weight * det_jacobi * (strain_L.transpose() * PK2_vec);
            // 切线刚度矩阵 线性部分+位移引起部分
            Matrix24d6 item_KL_mat_temp = strain_L.transpose() * C_matrix; // 24 x 6
            Matrix24d24 KL = item_KL_mat_temp * strain_L;                  // 24 x 24
            // 切线刚度矩阵 应力引起部分
            Matrix24d9 item_KN_mat_temp = strain_NL.transpose() * SIG; // 24 x 9
            Matrix24d24 KN = item_KN_mat_temp * strain_NL;             // 24 x 24
            stiffness_matrix = stiffness_matrix + weight * det_jacobi * (KL + KN);
        }
    }
}