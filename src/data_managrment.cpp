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
        for(int i=0;i<ele_type_sets.size();i++)
        {
            std::cout<<ele_type_sets[i]<<std::endl;
        }
    }

    void data_management::ele_inite(vector<string>& ele_type_sets, map<string, int>& ELE_TYPES, elastic_mat& data_mat)
    {
        map<string, int> type_idx;
        for (int i = 0; i < ele_type_sets.size();i++) {
            ele_base* ele;
            auto type = ele_type_sets[i];
            switch (ELE_TYPES[type])
            {
            case 1:
            {
                ele = new tetra_ele_elastic(data_mat);
                ele->build_cons_mat();
                ele_list_.push_back(ele);
                type_idx[type] = i;
                break; 
            }
            case 2:
            {
                ele = new hex_ele_elastic(data_mat);
                ele->build_cons_mat();
                ele_list_.push_back(ele);
                type_idx[type] = i;
                break;
            }
            default:
            {
                std::cout << "This type does not exist in the element library" << std::endl;
                break;
            }
            }
        }

        //单元名字对应其在ele_list中的索引
        elelist_idx_.resize(ne_);
        for (int i = 0; i < ne_; i++) {
            elelist_idx_[i] = type_idx[ele_type_[i]];
        }
    }

}