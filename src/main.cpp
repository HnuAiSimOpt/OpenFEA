#pragma once
#include <iostream>
#include "include/cae.h"
#include "include/mat.h"

void main()
{
    // 材料属性赋值
    CAE::elastic_mat mat_item{2.1e5, 0.3, 0.};

    // 材料路径
    std::string path = "E:/CADCAE_project/OpenFEA/model/mix_ele_model.inp";

    // 建立CAE分析对象
    CAE::CAE_process cae_item(path, mat_item);

    // 读取计算文件
    cae_item.pre_info();

    // 
}
