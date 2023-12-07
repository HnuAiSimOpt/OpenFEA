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
    CAE::elastic_mat mat_item{2.1e9, 0.3, 7800};

    // 材料路径
    std::string path = "E:\\CADCAE_project\\model\\C3D4_C3D8R\\Job-1.inp";

    // 关键字
    string load_set_keyword = "Set-load";
    string load_value_keyword = "Cload";
    string dis_set_keyword = "Set-fix";

    // 建立CAE分析对象
    CAE::CAE_process cae_item(path, mat_item);

    // 读取计算文件
    cae_item.pre_info(load_set_keyword, load_value_keyword, dis_set_keyword);

    // 执行结构响应分析
    string result_path = "E:\\CADCAE_project\\model\\output\\dis_C3D4_C3D8R.vtk";
    string path_abaqus = "E:\\CADCAE_project\\model\\Abaqus_U_C3D4_C3D8R.txt";
    cae_item.implict_analysis(result_path, path_abaqus);

//     string result_path = "F:/OpenFEM/model/";
//     string path_abaqus = "F:/OpenFEM/model/Abaqus_U.txt";
//     cae_item.data_cae_.time_total_ = 0.1;
//     cae_item.data_cae_.time_step_ = 0.0;
//     cae_item.explicit_analysis(result_path, path_abaqus);
}
