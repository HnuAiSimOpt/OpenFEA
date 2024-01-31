/**************************************************************************

Copyright:  WH team

Author: YinJichao <jichaoyinyjc@163.com>、ChenBinXiang <chen_mech99@163.com>

Completion date:  XXX

Description: XXX

**************************************************************************/

#pragma once

#include "ele_base.h"
#include "C3D8_Gauss.h"

namespace CAE
{
    class hex_ele_elastic : public ele_base
    {
    public:
        double det_jacobi_; // 雅可比矩阵行列式（中心积分点）
        elastic_mat matrial_struc_;
        Matrix6d6 C_matrix_;

    public:
        // 构造函数，析构函数
        hex_ele_elastic()
        {
            type_ = "C3D8";
            nnode_ = 8;
            node_dof_ = 3;
            ngps_ = 8;
        };
        hex_ele_elastic(elastic_mat matrial_struc) : matrial_struc_(matrial_struc)
        {
            type_ = "C3D8";
            nnode_ = 8;
            node_dof_ = 3;
            ngps_ = 8;
        };

        // 材料赋属性
        void set_matrial(elastic_mat matrial_struc) override { matrial_struc_ = matrial_struc; };

        // 建立本构矩阵
        virtual void build_cons_mat();

        // 建立应变矩阵(积分点)
        virtual void build_strain_mat(Eigen::Ref<Eigen::MatrixXd> node_coords, Eigen::Ref<Eigen::MatrixXd> strain_mat, vector<double> &gp_points, double *det_jacobi_point);

        // 建立单元刚度矩阵
        void build_ele_stiff_mat(Eigen::Ref<Eigen::MatrixXd> node_coords, Eigen::Ref<Eigen::MatrixXd> stiffness_matrix) override;

        // 建立单元密度矩阵
        void build_ele_mass(const vector<int> &node_topos, const vector<vector<double>> &coords, vector<double> &Mass) override;

        // 建立形函数
        void build_shape_fun(Eigen::Ref<Eigen::MatrixXd> node_coords, Eigen::MatrixXd &shape_fun, vector<double> &gp_points, double &det_jacobi_point);

        // 计算单元内力
        void cal_in_force(const vector<int> &node_topos, const vector<vector<double>> &real_coords, const vector<double> &disp_d,
                          vector<double> &stress, vector<double> &strain, vector<double> &InFroce) override;

        // 计算单元时间步长
        void update_timestep(Eigen::Ref<Eigen::MatrixXd> node_coords, double &time_step) override;
    };
}