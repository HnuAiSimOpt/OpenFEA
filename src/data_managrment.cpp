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
    // 筛选出单元种类
    void data_management::filter_ele_type(vector<string> &ele_type_sets)
    {
        ele_type_sets.resize(ele_type_.size());
        ele_type_sets.assign(ele_type_.begin(), ele_type_.end());
        sort(ele_type_sets.begin(), ele_type_sets.end());
        vector<string>::iterator pos = unique(ele_type_sets.begin(), ele_type_sets.end());
        ele_type_sets.erase(pos, ele_type_sets.end());
        for (int i = 0; i < ele_type_sets.size(); i++)
        {
            std::cout << ele_type_sets[i] << std::endl;
        }
    }


    void data_management::ele_inite(elastic_mat& data_mat)
    {
        for (auto ele : ele_list_) {
            ele->set_matrial(data_mat);
            ele->build_cons_mat();
        }
    }

    int data_management::add_ele(string& ele_type)
    {
        ele_base* ele = GetFactory<ele_baseFactory>::instance()->create_class(ele_type);
        if (ele == NULL) {
            std::cout << ele_type << " does not exist in the element library" << std::endl;
            exit(EXIT_FAILURE);
        }
        ele_list_.push_back(ele);
        return ele->nnode_;
    }

}