/**************************************************************************

Copyright:  WH team

Author: YinJichao <jichaoyinyjc@163.com>

Completion date:  XXX

Description: XXX

**************************************************************************/

#pragma once

#include "include/set_bcs.h"

namespace CAE
{
    // 建立无约束索引，被约束自由度置为-1，无约束自由度重新顺序编号
    void set_BCs::build_free_index(data_management &data_cae)
    {
        // 按升阶排列被约束节点编号
        sort(data_cae.dis_bc_set_.begin(), data_cae.dis_bc_set_.end());
        // 逐被约束节点重拍自由度编号
        int id_resort_idx = 0;
        int id_constraint = 0;
        for (int i = 0; i < data_cae.nd_; i++)
        {
            if (i == data_cae.dis_bc_set_[id_constraint] - 1) // 因为 .inp文件 从 1 开始编号， 故 -1
            {
                data_cae.resort_free_nodes.push_back(-1);
                if (id_constraint < data_cae.dis_bc_set_.size() - 1)
                {
                    id_constraint++;
                }
            }
            else
            {
                data_cae.resort_free_nodes.push_back(id_resort_idx);
                id_resort_idx++;
            }
        }
    }

    // 基于重排的自由度建立载荷向量
    void set_BCs::build_single_load(data_management &data_cae)
    {
        // 初始化载荷向量
        int num_free_node = data_cae.nd_ - data_cae.dis_bc_set_.size();
        data_cae.single_load_vec.resize(3 * num_free_node);
        std::fill(data_cae.single_load_vec.begin(), data_cae.single_load_vec.end(), 0.0);

        // 建立载荷向量
        int num_load = data_cae.load_set_.size();
        int load_node_idx, load_dof_idx;
        for (int i = 0; i < num_load; i++)
        {
            load_node_idx = data_cae.resort_free_nodes[data_cae.load_set_[i] - 1];
            if (data_cae.load_dof_ == 1)
            {
                load_dof_idx = 3 * load_node_idx;
            }
            else if (data_cae.load_dof_ == 2)
            {
                load_dof_idx = 3 * load_node_idx + 1;
            }
            else if (data_cae.load_dof_ == 3)
            {
                load_dof_idx = 3 * load_node_idx + 2;
            }
            else
            {
                cout << "the DOF of load node is error !!!";
            }
            data_cae.single_load_vec[load_dof_idx] = data_cae.load_value_;
        }
    }
}