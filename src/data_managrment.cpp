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
    };
}