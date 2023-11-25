/**************************************************************************

Copyright:  WH team

Author: YinJichao <jichaoyinyjc@163.com>

Completion date:  XXX

Description: XXX

**************************************************************************/

#pragma once
#include <iostream>
#include "include/cae.h"
#include "include/elastic_mat.h"

void main()
{
    // 材料属性赋值
    CAE::elastic_mat mat_item{2.1e5, 0.3, 0.};

    // 材料路径
    std::string path = "E:/CADCAE_project/OpenFEA/model/mix_ele_model.inp";

    // 关键字
    string load_set_keyword = "Set-load";
    string load_value_keyword = "Cload";
    string dis_set_keyword = "Set-constrain";

    // 建立CAE分析对象
    CAE::CAE_process cae_item(path, mat_item);

    // 读取计算文件
    cae_item.pre_info(load_set_keyword, load_value_keyword, dis_set_keyword);

    // 执行结构响应分析
    cae_item.implict_analysis();
}