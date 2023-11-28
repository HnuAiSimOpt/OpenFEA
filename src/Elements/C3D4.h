/**************************************************************************

Copyright:  WH team

Author: YinJichao <jichaoyinyjc@163.com>、ChenBinXiang <chen_mech99@163.com>

Completion date:  XXX

Description: XXX

**************************************************************************/

#pragma once

#include "ele_base.h"

namespace CAE
{
    class tetra_ele_elastic : public ele_base
    {
    public:
        double det_jacobi_; // 雅可比矩阵行列式(中心积分点)
        elastic_mat matrial_struc_;
        Matrix6d6 C_matrix_;
        Matrix6d12 strain_mat;

    public:
        // 构造函数，析构函数
        tetra_ele_elastic() { type_ = "C3D4"; node_dof_ = 3; nnode_ = 4; };
        tetra_ele_elastic(elastic_mat matrial_struc) : matrial_struc_(matrial_struc) { type_ = "C3D4"; node_dof_ = 3; nnode_ = 4; };

        // 材料赋属性
        void set_matrial(elastic_mat matrial_struc)override { matrial_struc_ = matrial_struc; }  ;

        // 建立本构矩阵
        virtual void build_cons_mat();

        // 建立应变矩阵(积分点)
        virtual void build_strain_mat(Eigen::Ref<Eigen::MatrixXd> node_coords, Eigen::Ref<Eigen::MatrixXd> strain_mat);


        // 建立单元刚度矩阵
        void build_ele_stiff_mat(Eigen::MatrixXd& node_coords, Eigen::MatrixXd& stiffness_matrix) override;

        // 建立单元密度矩阵
        // void build_ele_den_mat();
    };


}