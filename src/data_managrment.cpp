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
        ele_base* ele = GetFactory<ele_baseFactory>::instance()->create_class(ele_type);
        if (ele == NULL) {
            std::cout << ele_type << " does not exist in the element library" << std::endl;
            exit(EXIT_FAILURE);
        }
        ele_list_.push_back(ele);
        return ele->nnode_;
    }

}