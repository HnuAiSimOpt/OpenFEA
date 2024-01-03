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
#include <map>
#include <set>
#include "include/Factory.h"
#include "elements/include/ele_base.h"
#include "solver/include/solver_pardiso.h"

using std::map;
using std::set;
using std::sort;
using std::string;
using std::vector;

namespace CAE
{
    class data_management
    {
    public:
        int ne_, nd_;                                            // 单元，节点总数
        vector<vector<double>> coords_;                          // 节点坐标
        vector<vector<int>> node_topos_;                         // 节点拓扑关系
        vector<int> BndMesh_F;                                   // 非协调细网格单元编号
        vector<int> BndMesh_C;                                   // 非协调粗网格单元编号
        vector<vector<int>> bndFace_finemesh;                    // 非协调细网格 面单元节点编号
        vector<vector<int>> bndFace_coarsemesh;                  // 非协调粗网格 面单元节点编号
        vector<int> load_set_;                                   // 承载节点集合
        int load_dof_;                                           // 载荷自由度，1:X   2:Y   3:Z
        double load_value_;                                      // 承载幅值
        vector<int> dis_bc_set_;                                 // 约束节点集合
        vector<int> resort_free_nodes_;                          // 重排无约束自由度索引
        vector<double> single_load_vec_;                         // 基于重排自由度的单载荷向量
        vector<double> single_dis_vec_;                          // 仅考虑无约束自由度的位移向量
        vector<double> single_full_dis_vec_;                     // 考虑所有自由度的位移向量
        vector<ele_base *> ele_list_;                            // 单元类型列表
        vector<int> ele_list_idx_;                               // 各单元对应的ele_list的索引
        set<int> ele_class_;                                     // 单元总类
        double time_total_;                                      // 计算总时间
        double time_step_;                                       // 计算时间步长，为0则采用自动步长
        map<string, int> ELE_TYPES = {{"C3D4", 0}, {"C3D8", 1}}; // 建立单元类型到整型的映射
        PardisoSolution item_pardiso;                            // paidiso求解器对象

    public:
        //  单元初始化
        void ele_inite(elastic_mat &data_mat);
        // 工厂注册单元
        int add_ele(string &ele_type);
    };
}
