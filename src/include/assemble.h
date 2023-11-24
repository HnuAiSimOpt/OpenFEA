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
    typedef Eigen::Matrix<double, 4, 3> Matrix4d3;
    typedef Eigen::Matrix<double, 8, 3> Matrix8d3;
    class assamble_stiffness
    {
    public:
        int num_row;                                              // 稀疏矩阵行数
        int num_col;                                              // 稀疏矩阵列数
        int num_nz_val;                                           // 稀疏矩阵非零元素个数
        vector<double> nz_val;                                    // 储存每个非零元素值的值
        vector<int> row_idx;                                      // 储存每个非零元素值的行索引
        vector<int> col_idx;                                      // 储存每行第一个非零元素值的列索引
        map<string, int> ELE_TYPES = {{"C3D4", 1}, {"C3D8R", 2}}; // 建立单元类型到整型的映射

    public:
        // 建立压缩稀疏行（Compressed Sparse Row，CSR）索引
        void build_CSR(data_management &data_cae);

        // 基于单元编号，单元类型和节点拓扑关系，返回自由度
        void build_ele_dofs(vector<int> &item_ele_dofs, data_management &data_cae, int ele_id, string ele_type);

        // 基于CSR索引格式填充稀疏矩阵
        void fill_CSR_sparse_mat(data_management &data_cae);

        // 基于单元编号，单元类型和节点拓扑关系，返回自由度和节点坐标
        void build_tetra_dofs_coors(vector<int> &item_ele_dofs, Matrix4d3 &item_ele_coors,
                                    data_management &data_cae, int ele_id, string ele_type);
        void build_hex_dofs_coors(vector<int> &item_ele_dofs, Matrix8d3 &item_ele_coors,
                                  data_management &data_cae, int ele_id, string ele_type);
    };
}