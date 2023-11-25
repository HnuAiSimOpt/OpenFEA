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
            string item_ele_type = data_cae.ele_type_[id_ele];
            // 基于单元类型和节点拓扑关系，计算单元包含的自由度
            build_ele_dofs(item_ele_dofs, data_cae, id_ele, item_ele_type);
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
    void assamble_stiffness::build_ele_dofs(vector<int> &item_ele_dofs, data_management &data_cae, int ele_id, string ele_type)
    {
        int num_nodes = 0;
        switch (ELE_TYPES[ele_type])
        {
        case 1:
        {
            num_nodes = 4;
            break;
        }
        case 2:
        {
            num_nodes = 8;
            break;
        }
        default:
        {
            cout << "This type does not exist in the element library" << endl;
            break;
        }
        }
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
        // 初始化
        int num_nodes;
        vector<int> item_ele_dofs;
        Matrix4d3 item_tetra_coors;
        Matrix8d3 item_hex_coors;
        tetra_ele_elastic item_tetra_stiff(data_mat);
        item_tetra_stiff.build_cons_mat();
        hex_ele_elastic item_hex_stiff(data_mat);
        item_hex_stiff.build_cons_mat();
        for (int id_ele = 0; id_ele < data_cae.ne_; id_ele++)
        {
            // 识别单元类型
            string item_ele_type = data_cae.ele_type_[id_ele];
            // 基于单元类型和节点拓扑关系，计算单元包含的自由度
            switch (ELE_TYPES[item_ele_type])
            {
            case 1:
            {
                build_tetra_dofs_coors(item_ele_dofs, item_tetra_coors, data_cae, id_ele);
                Matrix12d12 ele_stiff;
                item_tetra_stiff.build_ele_stiff_mat(item_tetra_coors, ele_stiff);
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
                                } // 使用上三角矩阵
                                nz_val[t] = nz_val[t] + ele_stiff(mm, nn);
                            }
                        }
                    }
                }
                break;
            }
            case 2:
            {
                build_hex_dofs_coors(item_ele_dofs, item_hex_coors, data_cae, id_ele);
                Matrix24d24 ele_stiff;
                item_hex_stiff.build_ele_stiff_mat(item_hex_coors, ele_stiff);
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
                                }  // 使用上三角矩阵
                                nz_val[t] = nz_val[t] + ele_stiff(mm, nn);
                            }
                        }
                    }
                }
                break;
            }
            default:
            {
                cout << "This type does not exist in the element library" << endl;
                break;
            }
            }
        }
    }

    // 基于单元编号，单元类型和节点拓扑关系，返回自由度和节点坐标
    // 四面体
    void assamble_stiffness::build_tetra_dofs_coors(vector<int> &item_ele_dofs, Matrix4d3 &item_ele_coors,
                                                    data_management &data_cae, int ele_id)
    {
        item_ele_dofs.resize(12);
        std::fill(item_ele_dofs.begin(), item_ele_dofs.end(), 0);
        item_ele_coors.setZero();
        cout << item_ele_coors;
        int item_dof, item_node;
        for (int i = 0; i < 4; i++)
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
    }
    // 六面体
    void assamble_stiffness::build_hex_dofs_coors(vector<int> &item_ele_dofs, Matrix8d3 &item_ele_coors,
                                                  data_management &data_cae, int ele_id)
    {
        item_ele_dofs.resize(24);
        std::fill(item_ele_dofs.begin(), item_ele_dofs.end(), 0);
        item_ele_coors.setZero();
        int item_dof, item_node;
        for (int i = 0; i < 8; i++)
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
    }
}