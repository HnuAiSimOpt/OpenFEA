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
        int ne_, nd_;                     // 单元，节点总数
        vector<vector<double>> coords_;   // 节点坐标
        vector<vector<int>> node_topos_;  // 节点拓扑关系
        vector<string> ele_type_;         // 各单元类型
    };
}
