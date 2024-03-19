/**************************************************************************

Copyright:  WH team

Author: YinJichao <jichaoyinyjc@163.com>

Completion date:  XXX

Description: XXX

**************************************************************************/

#pragma once

#include <iostream>
#include <vector>
#include <set>
#include <map>
#include "Eigen/Dense"
#include "include/data_management.h"

using std::cout;
using std::endl;
using std::map;
using std::set;
using std::vector;

namespace CAE
{
    using Eigen::MatrixXd;

    class assamble_nl_stiffness
    {
    public:
        int num_row;           // 稀疏矩阵行数
        int num_col;           // 稀疏矩阵列数
        int num_nz_val;        // 稀疏矩阵非零元素个数
        vector<double> nz_val; // 储存每个非零元素值的值
        vector<int> row_idx;   // 储存每个非零元素值的行索引
        vector<int> col_idx;   // 储存每行第一个非零元素值的列索引

    public:
        // 建立压缩稀疏行（Compressed Sparse Row，CSR）索引
        void build_CSR(data_management &data_cae);

        // 基于单元编号，单元类型和节点拓扑关系，返回自由度
        void build_ele_dofs(vector<int> &item_ele_dofs, data_management &data_cae, int ele_id, int num_nodes);

        // 基于CSR索引格式填充稀疏矩阵
        void fill_CSR_sparse_mat(data_management &data_cae, elastic_mat &data_mat, vector<double> &current_dis_vec, vector<double> &inter_force_vec);

        // 基于单元编号，单元类型和节点拓扑关系，返回自由度，节点坐标，和结点位移
        void build_ele_dofs_coors(vector<int> &item_ele_dofs, Eigen::Ref<Eigen::MatrixXd> item_ele_coors,
                                  data_management &data_cae, int ele_id, int num_nodes);
    };
}