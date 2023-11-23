/**************************************************************************

Copyright:  WH team

Author: YinJichao <jichaoyinyjc@163.com>

Completion date:  XXX

Description: XXX

**************************************************************************/

#pragma once

#include <iostream>
#include <vector>
using std::vector;
using std::string;

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
        int load_dof_;                       // 载荷自由度，0:X   1:Y   2:Z
        double load_value_;                  // 承载幅值
        vector<int> dis_bc_set_;                // 约束节点集合
    };
}
