/**************************************************************************

Copyright:  WH team

Author: YinJichao <jichaoyinyjc@163.com>

Completion date:  XXX

Description: XXX

**************************************************************************/

#include <Eigen/Dense>
#include "include/assemble.h"

namespace CAE
{
    void assamble_stiffness::build_CSR(data_management &data_cae)
    {
        int num_free_node = data_cae.nd_ - data_cae.dis_bc_set_.size();
        vector<set<int>> col_data(3 * num_free_node);
        vector<int> item_ele_dofs;
        // 遍历单元，储存所有自由度
        for (int id_ele = 0; id_ele < data_cae.ne_; id_ele++)
        {
            // 识别单元类型
            int num_nodes = data_cae.ele_list_[data_cae.ele_list_idx_[id_ele]]->nnode_;
            // 基于单元类型和节点拓扑关系，计算单元包含的自由度
            build_ele_dofs(item_ele_dofs, data_cae, id_ele, num_nodes);
            // 删除负自由度，即被约束自由度
            for (auto it = item_ele_dofs.begin(); it != item_ele_dofs.end();)
            {
                if ((*it) < 0)
                {
                    it = item_ele_dofs.erase(it);
                }
                else
                {
                    ++it;
                }
            }
            // 压缩稀疏矩阵
            for (int id_dofs_row : item_ele_dofs)
            {
                for (int id_dofs_col : item_ele_dofs)
                {
                    col_data[id_dofs_col].insert(id_dofs_row);
                }
            }
        }
        // 建立CSR索引格式
        num_nz_val = 0;
        for (int id_dof = 0; id_dof < 3 * num_free_node; id_dof++)
        {
            num_nz_val += col_data[id_dof].size(); // 计算每个单元非零元数目
        }
        // 分配行索引容量
        row_idx.resize(num_nz_val);
        // 分配列索引容量
        col_idx.resize(3 * num_free_node + 1);
        // 建立列索引
        int item_idx_csr = 0;
        col_idx[0] = 0;
        for (int i = 0; i < 3 * num_free_node; i++)
        {
            for (int row : col_data[i])
            {
                row_idx[item_idx_csr] = row;
                item_idx_csr++;
            }
            col_idx[i + 1] = item_idx_csr;
        }
        cout << "The CSR index has been built" << endl;
    }

    // 基于单元类型和节点拓扑关系，返回自由度
    void assamble_stiffness::build_ele_dofs(vector<int> &item_ele_dofs, data_management &data_cae, int ele_id, int num_nodes)
    {
        item_ele_dofs.resize(3 * num_nodes);
        std::fill(item_ele_dofs.begin(), item_ele_dofs.end(), 0);
        int item_dof;
        for (int i = 0; i < num_nodes; i++)
        {
            item_dof = data_cae.resort_free_nodes_[data_cae.node_topos_[ele_id][i] - 1];
            item_ele_dofs[3 * i] = 3 * item_dof;
            item_ele_dofs[3 * i + 1] = 3 * item_dof + 1;
            item_ele_dofs[3 * i + 2] = 3 * item_dof + 2;
        }
    }

    // 基于CSR索引格式填充稀疏矩阵
    void assamble_stiffness::fill_CSR_sparse_mat(data_management &data_cae, elastic_mat &data_mat)
    {
        // 声明内存
        nz_val.resize(num_nz_val);
        // 初始化单元
        data_cae.ele_inite(data_mat);
        // 初始化单元坐标，刚度矩阵
        int num_nodes;
        vector<int> item_ele_dofs;
        MatrixXd item_ele_coors;
        MatrixXd stiffness_matrix;
        for (int id_ele = 0; id_ele < data_cae.ne_; id_ele++)
        {
            // 获取该单元的节点数量
            int node_num_ele = data_cae.ele_list_[data_cae.ele_list_idx_[id_ele]]->nnode_;
            // 查找节点自由度及坐标
            item_ele_coors.resize(node_num_ele, 3);
            build_ele_dofs_coors(item_ele_dofs, item_ele_coors, data_cae, id_ele, node_num_ele);
            // 计算单元刚度矩阵
            stiffness_matrix.resize(3 * node_num_ele, 3 * node_num_ele);
            data_cae.ele_list_[data_cae.ele_list_idx_[id_ele]]->build_ele_stiff_mat(item_ele_coors, stiffness_matrix);
            // 组装
            int ii_dof, jj_dof, loop_size = item_ele_dofs.size();
            for (int mm = 0; mm < loop_size; mm++)
            {
                jj_dof = item_ele_dofs[mm];
                if (jj_dof >= 0)
                {
                    int start = col_idx[jj_dof];
                    for (int nn = 0; nn < loop_size; nn++)
                    {
                        int t = start;
                        ii_dof = item_ele_dofs[nn];
                        if (ii_dof >= 0)
                        {
                            for (; row_idx[t] < ii_dof; t++)
                            {
                            } // 使得 t 对应的行索引 对应 ii_dof
                            nz_val[t] = nz_val[t] + stiffness_matrix(mm, nn);
                        }
                    }
                }
            }
        }
        cout << "The CSR index has been filled" << endl;
    }

    void assamble_stiffness::build_ele_dofs_coors(vector<int> &item_ele_dofs, Eigen::Ref<Eigen::MatrixXd> item_ele_coors,
                                                  data_management &data_cae, int ele_id, int num_nodes)
    {
        item_ele_dofs.resize(3 * num_nodes);
        std::fill(item_ele_dofs.begin(), item_ele_dofs.end(), 0);
        item_ele_coors.resize(num_nodes, 3);
        item_ele_coors.setZero();
        int item_dof, item_node;
        for (int i = 0; i < num_nodes; i++)
        {
            // 自由度
            item_dof = data_cae.resort_free_nodes_[data_cae.node_topos_[ele_id][i] - 1];
            item_ele_dofs[3 * i] = 3 * item_dof;
            item_ele_dofs[3 * i + 1] = 3 * item_dof + 1;
            item_ele_dofs[3 * i + 2] = 3 * item_dof + 2;
            // 坐标
            item_node = data_cae.node_topos_[ele_id][i] - 1;
            item_ele_coors(i, 0) = data_cae.coords_[item_node][0]; // X 坐标
            item_ele_coors(i, 1) = data_cae.coords_[item_node][1]; // Y 坐标
            item_ele_coors(i, 2) = data_cae.coords_[item_node][2]; // Z 坐标
        }
    };
}