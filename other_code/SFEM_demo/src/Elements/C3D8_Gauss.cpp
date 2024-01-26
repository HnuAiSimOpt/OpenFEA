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
        double det_jacobi_point = jacobi(0, 0) * jacobi(1, 1) * jacobi(2, 2) + jacobi(0, 1) * jacobi(1, 2) * jacobi(2, 0) + jacobi(0, 2) * jacobi(1, 0) * jacobi(2, 1) - jacobi(0, 2) * jacobi(1, 1) * jacobi(2, 0) - jacobi(1, 2) * jacobi(2, 1) * jacobi(0, 0) - jacobi(2, 2) * jacobi(0, 1) * jacobi(1, 0);
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
        //
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
}