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
        string path = "C:\\Users\\jicha\\Desktop\\new_model\\new_model.inp";
        string result_path = "C:\\Users\\jicha\\Desktop\\new_model\\ref_fem.vtk";
        string exepath = "E:\\WH_CAE\\CA_Topo_OpenFEA\\bin\\Debug\\OpenFEA.exe";
        // 关键字
        string load_set_keyword = "Set-load";
        string load_value_keyword = "Cload";
        string dis_set_keyword = "Set-fix";
        // 建立CAE分析对象
        CAE::CAE_process cae_item(path, mat_item);
        // 读取计算文件
        const char *args[] = {exepath.c_str(), path.c_str()};
        cae_item.Init(2, const_cast<char **>(args));
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
        string path = "C:\\Users\\jicha\\Desktop\\CA_Topo_OpenFEA\\OUTPUT\\ref_model.inp";
        string result_path = "C:\\Users\\jicha\\Desktop\\CA_Topo_OpenFEA\\OUTPUT\\ref_fem.vtk";
        string exepath = "C:\\Users\\jicha\\Desktop\\CA_Topo_OpenFEA\\bin\\Debug\\OpenFEA.exe";
        // 关键字
        string load_set_keyword = "Set-load";
        string load_value_keyword = "Cload";
        string dis_set_keyword = "Set-fix";
        // 建立CAE分析对象
        CAE::CAE_process cae_item(path, mat_item);
        // 读取计算文件
        const char *args[] = {exepath.c_str(), path.c_str()};
        cae_item.Init(2, const_cast<char **>(args));
        // 执行结构响应分析
        bool is_save_stiffness = true;
        cae_item.implict_analysis(result_path, is_save_stiffness);
        // -----------------------------------------------------------------------------------------------------------
        // 重分析
        // -----------------------------------------------------------------------------------------------------------
        // 修改网格路径(仅包含1个part的网格)
        string map_info = "C:\\Users\\jicha\\Desktop\\CA_Topo_OpenFEA\\OUTPUT\\mdf_info_map.inp";
        string mesh_path = "C:\\Users\\jicha\\Desktop\\CA_Topo_OpenFEA\\OUTPUT\\mdf_model.inp";
        string CA_result_path = "C:\\Users\\jicha\\Desktop\\CA_Topo_OpenFEA\\OUTPUT\\mdf_ca_2.vtk";
        cae_item.CA_pre_process(mesh_path, map_info);
        // 开始执行重分析
        int n_basis = 2;
        cae_item.CA_ReAnalysis(CA_result_path, n_basis);
    }
    else if (case_num == 3) // 拓扑优化
    {
        // 材料属性赋值
        CAE::elastic_mat mat_item{2.1e5, 0.3, 7800};
        // 材料路径
        string path = "C:\\Users\\jicha\\Desktop\\FangXiangPan.inp";
        string result_path = "C:\\Users\\jicha\\Desktop\\FangXiangPan.vtk";
        string exepath = "C:\\Users\\jicha\\Desktop\\CA_Topo_OpenFEA\\bin\\Debug\\OpenFEA.exe";
        // 关键字
        string load_set_keyword = "Set-load";
        string load_value_keyword = "Cload";
        string dis_set_keyword = "Set-fix";
        // 建立CAE分析对象
        CAE::CAE_process cae_item(path, mat_item);
        // 读取计算文件
        const char *args[] = {exepath.c_str(), path.c_str()};
        cae_item.Init(2, const_cast<char **>(args));
        // 拓扑优化
        double vol = 0.15;
        double rmin = 30;
        double penal = 3;
        double comp = cae_item.topo_process(result_path, vol, rmin, penal);
        
    }
    else
    {
        std::cout << "Please check your code\n";
    }
}