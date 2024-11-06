/**************************************************************************

Copyright:  WH team

Author: YinJichao <jichaoyinyjc@163.com>

Completion date:  XXX

Description: XXX

**************************************************************************/

#pragma once

#include "include/data_management.h"

namespace CAE
{

    void data_management::ele_inite(elastic_mat& data_mat)
    {
        for (auto ele : ele_list_) {
            ele->set_matrial(data_mat);
            ele->build_cons_mat();
        }
    }

    int data_management::add_ele(string& ele_type)
    {
        // 声明单元实例
        ele_base* ele = GetFactory<ele_baseFactory>::instance()->create_class(ele_type);
        if (ele == NULL) {
            std::cout << ele_type << " does not exist in the element library" << std::endl;
            exit(EXIT_FAILURE);
        }
        // 添加映射索引
        int map_idx = ELE_TYPES[ele_type];
        int ne = ele_list_.size();
        ele_map_list_.insert(std::pair<int, int>(map_idx, ne));
        // 添加单元
        ele_list_.push_back(ele);
        return ele->nnode_;
    }

}