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
#include "./data_management.h"
#include "./NCF_map.h"


using std::cout;
using std::endl;
using std::map;
using std::set;
using std::vector;

namespace CAE
{
    using Eigen::MatrixXd;

    class SFEM3D;
    class assamble_stiffness
    {
    public:
        int num_row_;                                              // 稀疏矩阵行数
        int num_col_;                                              // 稀疏矩阵列数
        int num_nz_val_;                                           // 稀疏矩阵非零元素个数
        vector<double> nz_val_;                                    // 储存每个非零元素值的值
        vector<int> row_idx_;                                      // 储存每个非零元素值的行索引
        vector<int> col_idx_;                                      // 储存每行第一个非零元素值的列索引

    public:
        // 建立压缩稀疏行（Compressed Sparse Row，CSR）索引
        void build_CSR(data_management &data_cae);

        // 基于单元编号，单元类型和节点拓扑关系，返回自由度
        void build_ele_dofs(vector<int>& item_ele_dofs, data_management& data_cae, int ele_id, int num_nodes);

        // 基于CSR索引格式填充稀疏矩阵
        void fill_CSR_sparse_mat(data_management &data_cae, elastic_mat &data_mat);

        // 基于单元编号，单元类型和节点拓扑关系，返回自由度和节点坐标
        void build_ele_dofs_coors(vector<int>& item_ele_dofs, Eigen::Ref<Eigen::MatrixXd> item_ele_coors,
            data_management& data_cae, int ele_id, int num_nodes);

        //[*******非协调部分********]
        // 组装非协调刚度矩阵
        void NCF_assembleStiffness(data_management& data_cae, elastic_mat& data_mat);
        void NCF_build_CSR(data_management& data_cae);
        void Storematrix_columns(vector<set<int>>& columns, vector<int>& A_eper_dof,
            vector<int>& B_eper_dof);

        //[*******光滑有限元部分********]
        void SFEM_build_CSR(SFEM3D* sfemData);

        void iNs_dofs(vector<int>& item_ele_dofs, SFEM3D* sfemData, int iNs);

        void SFEM_fill_CSR_sparse_mat(SFEM3D* sfemData, elastic_mat& data_mat);
    };
}