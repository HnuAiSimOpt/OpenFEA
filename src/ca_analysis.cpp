/**************************************************************************

Copyright:  WH team

Author: YinJichao <jichaoyinyjc@163.com>

Completion date:  XXX

Description: XXX

**************************************************************************/

#include "include/ca_reanalysis.h"
namespace CAE
{

    // 重分析信息预处理
    void ca_pre_process(data_management &data_cae)
    {
        // 拓扑关系
        for (auto part_topo : data_cae.parts_node_topos_)
        {
            for (auto topo : part_topo)
            {
                int nn_node = find_max_element(data_cae.node_topos_);
                for (int i = 0; i < topo.size(); i++)
                {
                    for (int j = 0; j < topo[i].size(); j++)
                    {
                        topo[i][j] += nn_node;
                    }
                }
                data_cae.node_topos_.insert(data_cae.node_topos_.end(), topo.begin(), topo.end());
            }
        }
        data_cae.ne_ = data_cae.node_topos_.size();
        // 节点坐标
        for (auto part_node : data_cae.parts_coods_)
            data_cae.coords_.insert(data_cae.coords_.end(), part_node.begin(), part_node.end());
        data_cae.nd_ = data_cae.coords_.size();
        // 单元类型
        for (auto e_type : data_cae.parts_ele_type_)
            data_cae.ele_list_idx_.insert(data_cae.ele_list_idx_.end(), e_type.begin(), e_type.end());
        // print
        std::cout << "Element:" << data_cae.ne_ << std::endl
                  << "Node: " << data_cae.nd_ << std::endl;
    }

    // 建立
    void ca_build_rom(data_management &data_cae)
    {
    }

    // 找二维向量最大值
    int find_max_element(vector<vector<int>> &vec)
    {
        int max_value = 0;
        for (int i = 0; i < vec.size(); i++)
        {
            for (int j = 0; j < vec[i].size(); j++)
            {
                if (vec[i][j] > max_value)
                {
                    max_value = vec[i][j];
                }
            }
        }
        return max_value;
    }
}