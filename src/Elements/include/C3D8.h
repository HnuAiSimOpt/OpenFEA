/**************************************************************************

Copyright:  WH team

Author: YinJichao <jichaoyinyjc@163.com>、ChenBinXiang <chen_mech99@163.com>

Completion date:  XXX

Description: XXX

**************************************************************************/

#pragma once

#include "ele_base.h"
#include "include/C3D8_Gauss.h"

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

        hex_ele_elastic() { type_ = "C3D8R"; nnode_ = 8; node_dof_ = 3; ngps_ = 8; face_node = 4; face_gps = 4; };
        hex_ele_elastic(elastic_mat matrial_struc) : matrial_struc_(matrial_struc) { type_ = "C3D8R"; nnode_ = 8; node_dof_ = 3; ngps_ = 8; face_node = 4; face_gps = 4;};

        // 材料赋属性
        void set_matrial(elastic_mat matrial_struc) override { matrial_struc_ = matrial_struc; };

        // 建立本构矩阵
        virtual void build_cons_mat();

        // 建立应变矩阵(积分点)
        virtual void build_strain_mat(Eigen::Ref<Eigen::MatrixXd> node_coords, Eigen::Ref<Eigen::MatrixXd> strain_mat, vector<double> &gp_points, double *det_jacobi_point);

        // 计算节点处位移应变转换矩阵
        // virtual void build_strain_node_mat(Eigen::Ref<Eigen::MatrixXd> node_coords, Eigen::Ref<Eigen::MatrixXd> strain_mat);

        // 计算节点处应力
        // virtual void get_stress_node(Eigen::Ref<Eigen::MatrixXd> node_coords, Eigen::Ref<Eigen::MatrixXd> stress_mat, Eigen::Ref<Eigen::MatrixXd> dis_vec);

        // 建立单元刚度矩阵 type 1：线弹性；type 2：几何非线性
        void build_ele_stiff_mat(Eigen::Ref<Eigen::MatrixXd> node_coords, Eigen::Ref<Eigen::MatrixXd> stiffness_matrix) override;

        // 建立切线单元刚度矩阵（几何非线性）
        void build_ele_nl_stiff_mat(Eigen::Ref<Eigen::MatrixXd> node_coords, Eigen::Ref<Eigen::MatrixXd> node_dis,
                                    Eigen::Ref<Eigen::MatrixXd> stiffness_matrix, Eigen::Ref<Eigen::MatrixXd> inter_force) override;

        // 建立单元密度矩阵
        void build_ele_mass(const vector<int> &node_topos, const vector<vector<double>> &coords, vector<double> &Mass) override;

        // 建立形函数
        void build_shape_fun(Eigen::Ref<Eigen::MatrixXd> node_coords, Eigen::MatrixXd &shape_fun, vector<double> &gp_points, double &det_jacobi_point);

        // 计算单元内力
        void cal_in_force(const vector<int> &node_topos, const vector<vector<double>> &real_coords, const vector<double> &disp_d,
                          vector<double> &stress, vector<double> &strain, vector<double> &InFroce) override;

        // 计算单元时间步长

        void update_timestep(Eigen::Ref<Eigen::MatrixXd> node_coords, double& time_step) override;

        // 交界面积分点物理空间坐标（非协调）
        void gps_phy_coords( Eigen::Ref<Eigen::MatrixXd> nodes1,
            Eigen::Ref<Eigen::MatrixXd> phy_gps, vector<double>& W_1,
            vector<Eigen::Vector3d>& Normal) override;
        

    };
}