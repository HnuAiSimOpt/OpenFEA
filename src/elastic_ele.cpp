/**************************************************************************

Copyright:  WH team

Author: YinJichao <jichaoyinyjc@163.com>

Completion date:  XXX

Description: XXX

**************************************************************************/

#pragma once

#include "include/elastic_ele.h"

namespace CAE
{
    //****************************************************************************//
    //                                四面体单元                                   //
    //****************************************************************************//
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
    void tetra_ele_elastic::build_strain_mat(Matrix4d3 &node_coords, Matrix6d12 &strain_mat)
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
    void tetra_ele_elastic::build_ele_stiff_mat(Matrix4d3 &node_coords, Matrix12d12 &ele_stiff)
    {
        // 基于 单点Hammer积分 计算四面体单元
        ele_stiff.setZero();
        double weight = 1. / 6.; // Hammer积分 权重
        Matrix6d12 strain_mat;
        build_strain_mat(node_coords, strain_mat);
        Matrix12d6 item_temp = strain_mat.transpose() * C_matrix_; // 12 x 6
        ele_stiff = item_temp * strain_mat;                        // 12 x 12
        ele_stiff = weight * det_jacobi_ * ele_stiff;
    }

    //****************************************************************************//
    //                                六面体单元                                   //
    //****************************************************************************//
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
    void hex_ele_elastic::build_strain_mat(Matrix8d3 &node_coords, Matrix6d24 &strain_mat, vector<double> &gp_points, double *det_jacobi_point)
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
        (*det_jacobi_point) = jacobi(0, 0) * jacobi(1, 1) * jacobi(2, 2) + jacobi(0, 1) * jacobi(1, 2) * jacobi(2, 0) + jacobi(0, 2) * jacobi(1, 0) * jacobi(2, 1) - jacobi(0, 2) * jacobi(1, 1) * jacobi(2, 0) - jacobi(1, 2) * jacobi(2, 1) * jacobi(0, 0) - jacobi(2, 2) * jacobi(0, 1) * jacobi(1, 0);
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
    }

    // 建立单元刚度矩阵
    void hex_ele_elastic::build_ele_stiff_mat(Matrix8d3 &node_coords, Matrix24d24 &ele_stiff)
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
        ele_stiff.setZero();
        Matrix6d24 strain_mat;
        Matrix24d6 item_temp_1;
        Matrix24d24 item_temp_2;
        for (int i = 0; i < 8; i++)
        {
            vector<double> gp_points = {gps(i, 0), gps(i, 1), gps(i, 2)};
            double det_jacobi_point;
            build_strain_mat(node_coords, strain_mat, gp_points, &det_jacobi_point);
            item_temp_1 = strain_mat.transpose() * C_matrix_; // 24 x 6
            item_temp_2 = item_temp_1 * strain_mat;           // 24 x 24
            ele_stiff = ele_stiff + weight * det_jacobi_point * item_temp_2;
        }
    }
}