/**************************************************************************

Copyright:  WH team

Author: YinJichao <jichaoyinyjc@163.com>

Completion date:  XXX

Description: XXX

**************************************************************************/

#pragma once

#include <iostream>
#include <vector>
#include <algorithm>

using std::string;
using std::vector;
using std::sort;

namespace CAE
{
    class data_management
    {
    public:
        int ne_, nd_;                        // 单元，节点总数
        vector<vector<double>> coords_;      // 节点坐标
        vector<vector<int>> node_topos_;     // 节点拓扑关系
        vector<string> ele_type_;            // 各单元类型
        vector<int> load_set_;               // 承载节点集合
        int load_dof_;                       // 载荷自由度，1:X   2:Y   3:Z
        double load_value_;                  // 承载幅值
        vector<int> dis_bc_set_;             // 约束节点集合
        vector<int> resort_free_nodes_;      // 重排无约束自由度索引
        vector<double> single_load_vec_;     // 基于重排自由度的单载荷向量
        vector<double> single_dis_vec_;      // 仅考虑无约束自由度的位移向量
        vector<double> single_full_dis_vec_; // 考虑所有自由度的位移向量

    public:
        // 单元种类 去重
        void filter_ele_type(vector<string> &ele_typr_sets);
    };
}
