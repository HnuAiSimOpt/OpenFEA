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
#include "./Factory.h"
#include "../Elements/include/ele_base.h"
#include "../solver/include/solver_superlu.h"
#include "../solver/include/solver_pardiso.h"

using std::map;
using std::set;
using std::sort;
using std::string;
using std::vector;

namespace CAE
{
    struct ca_mesh_map_info
    {
        vector<vector<double>> union_coords_;   // 修改前后模型前后所有节点（节点集合）坐标
        vector<vector<int>> change_node_topos_; // 修改的单元拓扑, [type, N1, N2, ... Nn] ，第1个为修改类型标记，0：删除单元，1：移动单元，2：增加单元
        vector<int> change_ele_idx_;            // 修改单元对应的ele_list的索引
        vector<int> node_idx_union_map_;        // 修改节点索引映射（节点合集在修改模型中节点的映射，如果为删除标记为-1）
    };

    struct ca_mesh_new_info
    {
        vector<vector<double>> mdf_coords_;  // 修改的模型的节点坐标
        vector<vector<int>> mdf_node_topos_; // 修改的模型的拓扑关系
        vector<int> mdf_ele_list_idx_;       // 修改后模型(所有单元), 各单元对应的ele_list的索引
        vector<double> ca_dis_vec_show_;     // 考虑所有自由度的位移向量(用于重分析)
    };

    struct NodeData                 // 定义一个结构体来表示节点的属性，包括节点编号、自由度
    {
        int node_id;                // 初始节点编号
        int free_node_id;           // 重排节点编号
        vector<int> full_dof;       // 初始全自由度
        vector<int> re_free_dof;    // 重排自由度
    };

    // 避免data_mana和assamble交叉引用
    class assamble_stiffness_save
    {
    public:
        int num_row_;               // 稀疏矩阵行数
        int num_col_;               // 稀疏矩阵列数
        int num_nz_val_;            // 稀疏矩阵非零元素个数
        vector<int> col_idx_;       // 储存每行第一个非零元素值的列索引
        vector<int> row_idx_;       // 储存每个非零元素值的行索引
        vector<double> nz_val_;     // 储存每个非零元素值的值
        vector<set<int>> col_data_; //  用于建立CSR索引
    };

    struct save_ref_info
    {
        assamble_stiffness_save stiff_ori_; // 完整隐式分析的总刚
        vector<double> dis_vec_ori_;        // 仅考虑无约束自由度的位移向量-原始模型
        vector<double> load_vec_ori_;       // 基于重排自由度的单载荷向量-原始模型
        vector<int> resort_free_nodes_ori_; // 重排无约束自由度索引-原始模型
    };

    class data_management
    {
    public:
        int ne_, nd_;                                            // 单元，节点总数
        vector<NodeData> nodes_;                                 // 存储节点数据
        int re_free_dof_num;                                     // 无约束节点自由度总数
        vector<vector<double>> coords_;                          // 节点坐标
        vector<vector<int>> node_topos_;                         // 节点拓扑关系
        vector<int> BndMesh_F;                                   // 非协调细网格单元编号(读inp）
        vector<int> BndMesh_C;                                   // 非协调粗网格单元编号(读inp）
        vector<vector<int>> bndFace_finemesh;                    // 非协调细网格 面单元节点编号(读inp）
        vector<vector<int>> bndFace_coarsemesh;                  // 非协调粗网格 面单元节点编号(读inp）
        vector<int> C_mesh;                                      // 与积分点对应的非协调粗网格单元编号(重排）
        vector<int> F_mesh;                                      // 与积分点对应的非协调细网格单元编号(重排）
        vector<int> load_set_;                                   // 承载节点集合
        int load_dof_;                                           // 载荷自由度，1:X   2:Y   3:Z
        double load_value_;                                      // 承载幅值
        vector<int> dis_bc_set_;                                 // 约束节点集合
        vector<int> resort_free_nodes_;                          // 重排无约束自由度索引
        vector<double> single_load_vec_;                         // 基于重排自由度的单载荷向量
        vector<double> single_dis_vec_;                          // 仅考虑无约束自由度的位移向量
        vector<double> single_full_dis_vec_;                     // 考虑所有自由度的位移向量
        vector<vector<double>> stress_mat_;                      // 针对每个"单元"的 柯西应力(前6列) + 冯米塞斯应力(最后一列)
        vector<vector<double>> stress_node_mat_;                 // 针对每个"节点"的 柯西应力(前6列) + 冯米塞斯应力(最后一列)
        vector<ele_base *> ele_list_;                            // 单元类型列表
        map<int, int> ele_map_list_;                             // 单元类型影射列表，即映射单元在ele_list_的索引
        vector<int> ele_list_idx_;                               // 各单元对应的ele_list的索引
        set<int> ele_class_;                                     // 单元总类
        double time_total_;                                      // 计算总时间
        double time_step_;                                       // 计算时间步长，为0则采用自动步长
        map<string, int> ELE_TYPES = {{"C3D4", 0}, {"C3D8", 1}}; // 建立单元类型到整型的映射
        PardisoSolution item_pardiso;                            // paidiso求解器对象
        PardisoSolution item_pardiso_ca;                            // paidiso求解器对象
        SuperLUSolution item_superlu;                            // superlu求解器对象
        bool NLFEA = false;                                      // 是否开启非线性计算
        double res_lmit = 1.0e-6;                                // 默认非线性的残差收敛标准

        // for CA
        ca_mesh_new_info mesh_ca_new_; // 修改后的模型
        ca_mesh_map_info mesh_ca_map_; // 网格映射数据
        save_ref_info save_ref_info_;  // 参考模型需要保存的数据
    public:
        //  单元初始化
        void ele_inite(elastic_mat &data_mat);
        // 工厂注册单元
        int add_ele(string &ele_type);
    };
}
