/**************************************************************************

Copyright:  WH team

Author: YinJichao <jichaoyinyjc@163.com>

Completion date:  XXX

Description: XXX

**************************************************************************/

#pragma once
#include "ele_base.h"

typedef Eigen::Matrix<double, 6, 1> Matrix6d1;
typedef Eigen::Matrix<double, 9, 9> Matrix9d9;
typedef Eigen::Matrix<double, 9, 24> Matrix9d24;
typedef Eigen::Matrix<double, 24, 9> Matrix24d9;
typedef Eigen::Matrix<double, 24, 24> Matrix24d24;

namespace CAE
{
    // 建立应变矩阵(积分点)
    double build_strain_mat_gauss(Eigen::Ref<Eigen::MatrixXd> node_coords, Eigen::Ref<Eigen::MatrixXd> strain_mat, vector<double> &gp_points);

    // 建立单元刚度矩阵
    void build_ele_stiff_mat_gauss(Eigen::Ref<Eigen::MatrixXd> node_coords, Eigen::Ref<Eigen::MatrixXd> stiffness_matrix, const Matrix6d6 &C_matrix);

    // 考虑几何非线性建立应变矩阵 Green-Lagrangian应变张量
    double build_GL_strain_gauss(Eigen::Ref<Eigen::MatrixXd> node_coords, Eigen::Ref<Eigen::MatrixXd> ele_dis,
                                 Eigen::Ref<Eigen::MatrixXd> strain_L, Eigen::Ref<Eigen::MatrixXd> strain_NL, 
                                 Matrix6d1 &PK2_vec, const Matrix6d6 &C_matrix, vector<double> &gp_points);

    // 建立切线单元刚度矩阵
    void build_ele_nl_stiff_mat_gauss(Eigen::Ref<Eigen::MatrixXd> node_coords, Eigen::Ref<Eigen::MatrixXd> ele_dis,
                                      Eigen::Ref<Eigen::MatrixXd> stiffness_matrix, Eigen::Ref<Eigen::MatrixXd> inter_force, const Matrix6d6 &C_matrix);
}