/**************************************************************************

Copyright:  WH team

Author: YinJichao <jichaoyinyjc@163.com>

Completion date:  XXX

Description: XXX

**************************************************************************/

#include "include/post_process.h"

namespace CAE
{
    // 重置全分析位移
    bool simulation_post::reset_displacement(data_management &data_cae)
    {
        data_cae.single_full_dis_vec_.resize(3 * data_cae.nd_);
        for (int i = 0; i < data_cae.nd_; i++)
        {
            int resort_node = data_cae.resort_free_nodes_[i];
            if (resort_node >= 0)
            {
                data_cae.single_full_dis_vec_[3 * i] = data_cae.single_dis_vec_[3 * resort_node];
                data_cae.single_full_dis_vec_[3 * i + 1] = data_cae.single_dis_vec_[3 * resort_node + 1];
                data_cae.single_full_dis_vec_[3 * i + 2] = data_cae.single_dis_vec_[3 * resort_node + 2];
            }
            else
            {
                data_cae.single_full_dis_vec_[3 * i] = 0.0;
                data_cae.single_full_dis_vec_[3 * i + 1] = 0.0;
                data_cae.single_full_dis_vec_[3 * i + 2] = 0.0;
            }
        }
        cout << "the full displacement has been filled." << endl;
        return true;
    }

    // 重置组合近似分析位移
    bool simulation_post::reset_ca_displacement(data_management &data_cae, vector<double> &dis)
    {
        // dis 未包含位移约束点，故还原到完整全位移
        data_cae.single_full_ca_dis_vec_.resize(3 * data_cae.coords_m_.size());
        for (int i = 0; i < data_cae.coords_m_.size(); i++)
        {
            int resort_node = data_cae.resort_free_nodes_o_[i];
            if (resort_node >= 0)
            {
                data_cae.single_full_ca_dis_vec_[3 * i] = dis[3 * resort_node];
                data_cae.single_full_ca_dis_vec_[3 * i + 1] = dis[3 * resort_node + 1];
                data_cae.single_full_ca_dis_vec_[3 * i + 2] = dis[3 * resort_node + 2];
            }
            else
            {
                data_cae.single_full_ca_dis_vec_[3 * i] = 0.0;
                data_cae.single_full_ca_dis_vec_[3 * i + 1] = 0.0;
                data_cae.single_full_ca_dis_vec_[3 * i + 2] = 0.0;
            }
        }
        cout << "the full displacement has been filled." << endl;

        // 删除修改后结构不包含的节点
        std::vector<double> new_full_dis_vec;
        new_full_dis_vec.resize(3 * data_cae.coords_mfull_.size());
        for (int i = 0; i < data_cae.node_idx_m_.size(); i++)
        {
            if (data_cae.node_idx_m_[i] == -1)
            {
                continue;
            }
            new_full_dis_vec[data_cae.node_idx_m_[i] * 3] = data_cae.single_full_ca_dis_vec_[i * 3];
            new_full_dis_vec[data_cae.node_idx_m_[i] * 3 + 1] = data_cae.single_full_ca_dis_vec_[i * 3 + 1];
            new_full_dis_vec[data_cae.node_idx_m_[i] * 3 + 2] = data_cae.single_full_ca_dis_vec_[i * 3 + 2];
        }
        data_cae.single_full_ca_dis_vec_ = std::move(new_full_dis_vec);
        cout << "the displacement has been mapped to the new model. the size is: " << data_cae.single_full_ca_dis_vec_.size() << endl;
        return true;
    }

    // 计算节点应力
    bool simulation_post::get_cauchy_stress_3d(data_management &data_cae)
    {
        data_cae.stress_mat_.resize(data_cae.ne_);
        data_cae.stress_node_mat_.resize(data_cae.nd_);
        MatrixXd item_ele_coors;
        MatrixXd item_ele_disp;
        MatrixXd item_ele_stress;
        vector<int> item_ele_dofs;
        for (int id_ele = 0; id_ele < data_cae.ne_; id_ele++)
        {
            int ele_type = data_cae.ele_list_idx_[id_ele];          // 单元类型
            int map_idx = data_cae.ele_map_list_[ele_type];         // 单元映射索引
            int node_num_ele = data_cae.ele_list_[map_idx]->nnode_; // 单元节点数
            // 查找节点自由度及坐标
            item_ele_disp.resize(3 * node_num_ele, 1);
            item_ele_disp.setZero();
            item_ele_coors.resize(node_num_ele, 3);
            item_ele_coors.setZero();
            item_ele_coors.resize(node_num_ele, 3);
            build_ele_dofs_dis_coors(item_ele_dofs, item_ele_disp, item_ele_coors, data_cae, id_ele, node_num_ele);
            // 计算柯西应力
            item_ele_stress.resize(6, 1);
            item_ele_stress.setZero();
            data_cae.ele_list_[map_idx]->get_stress_node(item_ele_coors, item_ele_stress, item_ele_disp);
            // 米塞斯应力
            double vm_s = std::sqrt(0.5 *
                                    ((item_ele_stress(0, 0) - item_ele_stress(1, 0)) * (item_ele_stress(0, 0) - item_ele_stress(1, 0)) +
                                     (item_ele_stress(1, 0) - item_ele_stress(2, 0)) * (item_ele_stress(1, 0) - item_ele_stress(2, 0)) +
                                     (item_ele_stress(2, 0) - item_ele_stress(0, 0)) * (item_ele_stress(2, 0) - item_ele_stress(0, 0)) +
                                     6.0 * (item_ele_stress(3, 0) * item_ele_stress(3, 0) +
                                            item_ele_stress(4, 0) * item_ele_stress(4, 0) +
                                            item_ele_stress(5, 0) * item_ele_stress(5, 0))));
            vector<double> stress = {item_ele_stress(0, 0),
                                     item_ele_stress(1, 0),
                                     item_ele_stress(2, 0),
                                     item_ele_stress(3, 0),
                                     item_ele_stress(4, 0),
                                     item_ele_stress(5, 0),
                                     vm_s};
            // cout << "S_VM ele"<<id_ele<< ":" << stress[2] << endl;
            data_cae.stress_mat_[id_ele] = stress;
        }

        cout << "the stress has been calculated." << endl;
        // 初始化节点应力矩阵
        data_cae.stress_node_mat_.resize(7);
        vector<double> node_count(data_cae.nd_);
        std::fill(node_count.begin(), node_count.end(), 0); // 统计节点累计次数
        for (int i = 0; i < 7; i++)
        {
            data_cae.stress_node_mat_[i].resize(data_cae.nd_);
            std::fill(data_cae.stress_node_mat_[i].begin(), data_cae.stress_node_mat_[i].end(), 0);
        }
        // 填充节点应力
        for (int id_ele = 0; id_ele < data_cae.ne_; id_ele++)
        {
            int ele_type = data_cae.ele_list_idx_[id_ele];          // 单元类型
            int map_idx = data_cae.ele_map_list_[ele_type];         // 单元映射索引
            int node_num_ele = data_cae.ele_list_[map_idx]->nnode_; // 单元节点数
            int item_node_id;
            // 自由度
            for (int i = 0; i < node_num_ele; i++)
            {
                item_node_id = data_cae.node_topos_[id_ele][i] - 1;
                data_cae.stress_node_mat_[0][item_node_id] = data_cae.stress_node_mat_[0][item_node_id] + data_cae.stress_mat_[id_ele][0];
                data_cae.stress_node_mat_[1][item_node_id] = data_cae.stress_node_mat_[1][item_node_id] + data_cae.stress_mat_[id_ele][1];
                data_cae.stress_node_mat_[2][item_node_id] = data_cae.stress_node_mat_[2][item_node_id] + data_cae.stress_mat_[id_ele][2];
                data_cae.stress_node_mat_[3][item_node_id] = data_cae.stress_node_mat_[3][item_node_id] + data_cae.stress_mat_[id_ele][3];
                data_cae.stress_node_mat_[4][item_node_id] = data_cae.stress_node_mat_[4][item_node_id] + data_cae.stress_mat_[id_ele][4];
                data_cae.stress_node_mat_[5][item_node_id] = data_cae.stress_node_mat_[5][item_node_id] + data_cae.stress_mat_[id_ele][5];
                data_cae.stress_node_mat_[6][item_node_id] = data_cae.stress_node_mat_[6][item_node_id] + data_cae.stress_mat_[id_ele][6];
                node_count[item_node_id] = node_count[item_node_id] + 1.;
            }
        }
        // 节点平均
        for (int id_node = 0; id_node < data_cae.nd_; id_node++)
        {
            if (node_count[id_node] == 0)
                continue;
            data_cae.stress_node_mat_[0][id_node] = data_cae.stress_node_mat_[0][id_node] / node_count[id_node];
            data_cae.stress_node_mat_[1][id_node] = data_cae.stress_node_mat_[1][id_node] / node_count[id_node];
            data_cae.stress_node_mat_[2][id_node] = data_cae.stress_node_mat_[2][id_node] / node_count[id_node];
            data_cae.stress_node_mat_[3][id_node] = data_cae.stress_node_mat_[3][id_node] / node_count[id_node];
            data_cae.stress_node_mat_[4][id_node] = data_cae.stress_node_mat_[4][id_node] / node_count[id_node];
            data_cae.stress_node_mat_[5][id_node] = data_cae.stress_node_mat_[5][id_node] / node_count[id_node];
            data_cae.stress_node_mat_[6][id_node] = data_cae.stress_node_mat_[6][id_node] / node_count[id_node];
            // cout << "S_VM node" << id_node << ":" << data_cae.stress_node_mat_[2][id_node] << endl;
        }
        cout << "the nodal stress has been smoothed." << endl;
        return true;
    }

    // 基于单元编号，单元类型和节点拓扑关系，返回自由度
    void simulation_post::build_ele_dofs_dis_coors(vector<int> &item_ele_dofs, Eigen::Ref<Eigen::MatrixXd> item_ele_disp,
                                                   Eigen::Ref<Eigen::MatrixXd> item_ele_coors, data_management &data_cae, int ele_id, int num_nodes)
    {
        item_ele_dofs.resize(3 * num_nodes);
        std::fill(item_ele_dofs.begin(), item_ele_dofs.end(), 0);
        int item_dof, item_node;
        for (int i = 0; i < num_nodes; i++)
        {
            // 自由度
            item_dof = data_cae.resort_free_nodes_[data_cae.node_topos_[ele_id][i] - 1];
            item_ele_dofs[3 * i] = 3 * item_dof;
            item_ele_dofs[3 * i + 1] = 3 * item_dof + 1;
            item_ele_dofs[3 * i + 2] = 3 * item_dof + 2;
            // 位移
            item_dof = data_cae.node_topos_[ele_id][i] - 1;
            item_ele_disp(3 * i, 0) = data_cae.single_full_dis_vec_[3 * item_dof];
            item_ele_disp(3 * i + 1, 0) = data_cae.single_full_dis_vec_[3 * item_dof + 1];
            item_ele_disp(3 * i + 2, 0) = data_cae.single_full_dis_vec_[3 * item_dof + 2];
            // 坐标
            item_node = data_cae.node_topos_[ele_id][i] - 1;
            item_ele_coors(i, 0) = data_cae.coords_[item_node][0]; // X 坐标
            item_ele_coors(i, 1) = data_cae.coords_[item_node][1]; // Y 坐标
            item_ele_coors(i, 2) = data_cae.coords_[item_node][2]; // Z 坐标
        }
    }
}