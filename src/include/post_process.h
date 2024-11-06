/**************************************************************************

Copyright:  WH team

Author: YinJichao <jichaoyinyjc@163.com>

Completion date:  XXX

Description: XXX

**************************************************************************/

#pragma once
#include <cmath>
#include "include/data_management.h"

using std::cout;
using std::endl;

namespace CAE
{
    using Eigen::MatrixXd;
    typedef Eigen::Matrix<double, 12, 1> Matrix6d1;

    class simulation_post
    {
    public:
        // 构造函数，析构函数
        simulation_post() {};

        // 重置全分析位移
        bool reset_displacement(data_management &data_cae);

        // 重置组合近似分析位移
        bool reset_ca_displacement(data_management &data_cae, vector<double> &dis);

        // 计算节点应力
        bool get_cauchy_stress_3d(data_management &data_cae);

        // 基于单元编号，单元类型和节点拓扑关系，返回自由度和节点坐标
        void build_ele_dofs_dis_coors(vector<int> &item_ele_dofs, Eigen::Ref<Eigen::MatrixXd> item_ele_coors,
                                      Eigen::Ref<Eigen::MatrixXd> item_ele_disp, data_management &data_cae, int ele_id, int num_nodes);
    };
}