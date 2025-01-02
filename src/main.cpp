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
#include "include/sample_eigen_superlu_mkl_svd.h"

void code_test();

int main(int argc, char *argv[])
{
    bool test_flag = true; // 设置true为以前的启动模式
    if (test_flag)
    {
        code_test();
        return 0;
    }
    else
    {
        // 建立CAE分析对象
        CAE::CAE_process cae_item;
        // 初始化
        cae_item.Init(argc, argv);
        // 求解
        cae_item.Solve();
        return 0;
    }
}

void code_test()
{
    int case_num = 1;
    if (case_num == 1)
    {
        // 材料属性赋值
        CAE::elastic_mat mat_item{2.1e5, 0.3, 7800};
        // 材料路径
        // std::string path = "C:\\Users\\jicha\\Desktop\\test\\reduced\\Job-CASL.inp";
        // string result_path = "C:\\Users\\jicha\\Desktop\\test\\output\\verify_ca_reduced.vtk";
        std::string path = "C:\\Users\\jicha\\Desktop\\test\\modified\\modify_for_full.inp";
        string result_path = "C:\\Users\\jicha\\Desktop\\test\\output\\verify_ca_modify.vtk";
        // 关键字
        string load_set_keyword = "Set-load";
        string load_value_keyword = "Cload";
        string dis_set_keyword = "Set-fix";
        // 建立CAE分析对象
        CAE::CAE_process cae_item(path, mat_item);
        // 读取计算文件
        cae_item.pre_info(load_set_keyword, load_value_keyword, dis_set_keyword);
        // 执行结构响应分析
        cae_item.implict_analysis(result_path);
    }
    else if (case_num == 2) // 隐式完整分析+重分析
    {
        // -----------------------------------------------------------------------------------------------------------
        // 隐式 全分析
        // -----------------------------------------------------------------------------------------------------------
        // 材料属性赋值
        CAE::elastic_mat mat_item{2.1e5, 0.3, 7800};
        // 材料路径
        // string path = "C:\\Users\\jicha\\Desktop\\test\\reduced\\Job-SL.inp";
        // string result_path = "C:\\Users\\jicha\\Desktop\\test\\output\\Job-SL.vtk";
        string path = "C:\\Users\\jicha\\Desktop\\test\\modified\\original.inp";  // modify_for_full.inp, original.inp
        string result_path = "C:\\Users\\jicha\\Desktop\\test\\output\\original.vtk";
        // 关键字
        string load_set_keyword = "Set-load";
        string load_value_keyword = "Cload";
        string dis_set_keyword = "Set-fix";
        // 建立CAE分析对象
        CAE::CAE_process cae_item(path, mat_item);
        // 读取计算文件
        cae_item.pre_info(load_set_keyword, load_value_keyword, dis_set_keyword);
        // 执行结构响应分析
        
        bool is_save_stiffness = true;
        cae_item.implict_analysis(result_path, is_save_stiffness);
        // -----------------------------------------------------------------------------------------------------------
        // 重分析
        // -----------------------------------------------------------------------------------------------------------
        // 修改网格路径(仅包含1个part的网格)
        // string map_info = "C:\\Users\\jicha\\Desktop\\test\\reduced\\map_info.inp";
        // string mesh_path = "C:\\Users\\jicha\\Desktop\\test\\reduced\\m_model.inp";
        // string CA_result_path = "C:\\Users\\jicha\\Desktop\\test\\output\\ca_reduced.vtk";
        string map_info = "C:\\Users\\jicha\\Desktop\\test\\modified\\map_info.inp";
        string mesh_path = "C:\\Users\\jicha\\Desktop\\test\\modified\\modify.inp";
        string CA_result_path = "C:\\Users\\jicha\\Desktop\\test\\output\\ca_modify.vtk";
        cae_item.CA_pre_process(mesh_path, map_info);
        // 开始执行重分析
        int n_basis = 4;
        cae_item.CA_ReAnalysis(CA_result_path, n_basis);
    }
    else
    {
        std::cout << "Please check your code\n";
    }
}