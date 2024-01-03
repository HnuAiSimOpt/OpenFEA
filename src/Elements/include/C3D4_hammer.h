/**************************************************************************

Copyright:  WH team

Author: YinJichao <jichaoyinyjc@163.com>

Completion date:  XXX

Description: XXX

**************************************************************************/

#pragma once
#include "ele_base.h"

namespace CAE
{
    // 建立应变矩阵(积分点)
    double build_strain_mat_hammer(Eigen::Ref<Eigen::MatrixXd> node_coords, Eigen::Ref<Eigen::MatrixXd> strain_mat);

    // 建立单元刚度矩阵
    void build_ele_stiff_mat_hammer(Eigen::Ref<Eigen::MatrixXd> node_coords, Eigen::Ref<Eigen::MatrixXd> stiffness_matrix, const Matrix6d6 &C_matrix);

}