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
#include "include/sample_eigen_superlu_mkl.h"

void main()
{
    int case_num = 2;
    if (case_num == 1)
    {
        // 材料属性赋值
        CAE::elastic_mat mat_item{2.1e9, 0.3, 7800};
        // 材料路径
        std::string path = "E:\\CADCAE_BY_ME\\model\\C3D4_C3D8\\Job-1.inp";
        // 关键字
        string load_set_keyword = "Set-load";
        string load_value_keyword = "Cload";
        string dis_set_keyword = "Set-fix";
        // 建立CAE分析对象
        CAE::CAE_process cae_item(path, mat_item);
        // 读取计算文件
        cae_item.pre_info(load_set_keyword, load_value_keyword, dis_set_keyword);
        // 执行结构响应分析
        string result_path = "E:\\CADCAE_BY_ME\\output\\dis_C3D4_C3D8.vtk";
        string path_abaqus = "E:\\CADCAE_BY_ME\\model\\Abaqus_U_C3D4_C3D8.txt";
        cae_item.implict_analysis(result_path, path_abaqus);
    }
    else if (case_num == 2)//非协调计算
    {
        // 材料属性赋值
        CAE::elastic_mat mat_item{2.1e5, 0.3, 7800};//此材料参数用于非协调计算
        // 材料路径
        //std::string path = "E:/project/CADCAE/bracket_nonconforming/inp/model.inp";//c3d8非协调路径
        //std::string path = "E:/project/CADCAE/bracket_nonconforming/inp/Job-c3d4-random.inp";//非协调路径
        //std::string path = "E:/project/CADCAE/bracket_nonconforming/inp/Job-c3d4-ncf.inp";//非协调路径
        string path = "E:/project/CADCAE/bracket_nonconforming/inp/Job-c3d4-ele12k.inp";//非协调路径
        //std::string path = "E:/project/CADCAE/bracket_nonconforming/inp/Job-c3d4_110k.inp";//非协调路径
        //std::string path = "E:/project/CADCAE/bracket_nonconforming/inp/Job-c3d4-1000k.inp";//非协调路径
        // 关键字
        string load_set_keyword = "Set-load";
        string load_value_keyword = "Cload";
        string dis_set_keyword = "Set-fix";
        // 建立CAE分析对象
        CAE::CAE_process cae_item(path, mat_item);
        // 读取计算文件
        cae_item.pre_info(load_set_keyword, load_value_keyword, dis_set_keyword);
        // 执行结构响应分析
        string result_path = "E:/project/CADCAE/bracket_nonconforming/output/dis_with_error.vtk";
        string path_abaqus = "E:/project/CADCAE/bracket_nonconforming/output/Abaqus_U.txt";
        cae_item.implict_analysis(result_path, path_abaqus);

    }
    else
    {
        std::cout << "Please check your code\n";
    }
}
