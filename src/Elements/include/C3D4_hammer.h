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
typedef Eigen::Matrix<double, 9, 12> Matrix9d12;
typedef Eigen::Matrix<double, 12, 9> Matrix12d9;

namespace CAE
{
    // 建立应变矩阵(积分点)
    double build_strain_mat_hammer(Eigen::Ref<Eigen::MatrixXd> node_coords, Eigen::Ref<Eigen::MatrixXd> strain_mat);

    // 建立单元刚度矩阵
    void build_ele_stiff_mat_hammer(Eigen::Ref<Eigen::MatrixXd> node_coords, Eigen::Ref<Eigen::MatrixXd> stiffness_matrix, const Matrix6d6 &C_matrix);

    // 考虑几何非线性建立应变矩阵 Green-Lagrangian应变张量
    double build_GL_strain_hamer(Eigen::Ref<Eigen::MatrixXd> node_coords, Eigen::Ref<Eigen::MatrixXd> ele_dis,
                                 Eigen::Ref<Eigen::MatrixXd> strain_L, Eigen::Ref<Eigen::MatrixXd> strain_NL, Matrix6d1 &PK2_vec, const Matrix6d6 &C_matrix);

    // 建立切线单元刚度矩阵
    void build_ele_nl_stiff_mat_hammer(Eigen::Ref<Eigen::MatrixXd> node_coords, Eigen::Ref<Eigen::MatrixXd> ele_dis,
                                       Eigen::Ref<Eigen::MatrixXd> stiffness_matrix, Eigen::Ref<Eigen::MatrixXd> inter_force, const Matrix6d6 &C_matrix);

}