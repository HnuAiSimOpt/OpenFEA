/**************************************************************************

Copyright:  WH team

Author: YinJichao <jichaoyinyjc@163.com>、ChenBinXiang <chen_mech99@163.com>

Completion date:  XXX

Description: XXX

**************************************************************************/

#pragma once

#include "Eigen/Dense"
#include "Eigen/SVD"
#include <vector>
#include "../../include/elastic_mat.h"
#include "../../include/Factory.h"

typedef Eigen::Matrix<double, 3, 3> Matrix3d3;
typedef Eigen::Matrix<double, 3, 4> Matrix3d4;
typedef Eigen::Matrix<double, 3, 8> Matrix3d8;

typedef Eigen::Matrix<double, 4, 3> Matrix4d3;

typedef Eigen::Matrix<double, 6, 6> Matrix6d6;
typedef Eigen::Matrix<double, 6, 12> Matrix6d12;
typedef Eigen::Matrix<double, 6, 24> Matrix6d24;

typedef Eigen::Matrix<double, 8, 3> Matrix8d3;

typedef Eigen::Matrix<double, 12, 6> Matrix12d6;
typedef Eigen::Matrix<double, 12, 12> Matrix12d12;

typedef Eigen::Matrix<double, 24, 24> Matrix24d24;
typedef Eigen::Matrix<double, 24, 6> Matrix24d6;

using std::vector;

namespace CAE
{
    class ele_base
    {
    public:
        // 构造函数，析构函数
        ele_base() {};
    public:
        std::string type_;//单元类型名字
        int nnode_;//该单元拥有节点数量
        int node_dof_;//该单元每个节点自由度数
        int ngps_;//积分点数量
        // 赋值材料属性
        virtual void set_matrial(elastic_mat data_mat) {};

        // 建立本构矩阵
        virtual void build_cons_mat() {};

        // 建立应变矩阵
        virtual void build_strain_mat() {};

        // 建立单元刚度矩阵
        virtual void build_ele_stiff_mat(Eigen::Ref<Eigen::MatrixXd> node_coords, Eigen::Ref<Eigen::MatrixXd> stiffness_matrix) {};
    
        // 建立质量列阵
        virtual void build_ele_mass(const vector<int>& node_topos, const vector<vector<double>>& coords, vector<double>& Mass) {};

        // 计算单元内力
        virtual void cal_in_force(const vector<int>& node_topos, const vector<vector<double>>& real_coords, const vector<double>& disp_d,
                                  vector<double>& stress, vector<double>& strain, vector<double>& InFroce) {};
    
        // 计算单元时间步长
        virtual void update_timestep(Eigen::Ref<Eigen::MatrixXd> node_coords, double& time_step) {};
    };
    CREAT_FACTORY(ele_base);
}