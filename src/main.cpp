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

void main()
{
    int case_num = 101;
    if (case_num == 1)
    {
        // 材料属性赋值
        CAE::elastic_mat mat_item{2.1e9, 0.3, 7800};
        // 材料路径
        std::string path = "E:\\WH_CAE\\for_CA\\sample_model\\CA\\full.inp";
        // 关键字
        string load_set_keyword = "Set-load";
        string load_value_keyword = "Cload";
        string dis_set_keyword = "Set-fix";
        // 建立CAE分析对象
        CAE::CAE_process cae_item(path, mat_item);
        // 读取计算文件
        cae_item.pre_info(load_set_keyword, load_value_keyword, dis_set_keyword);
        // 执行结构响应分析
        string result_path = "E:\\WH_CAE\\for_CA\\sample_model\\CA\\ca_solution.vtk";
        string path_abaqus = "E:\\WH_CAE\\for_CA\\sample_model\\CA\\ca_solution.txt";
        cae_item.implict_analysis(result_path, path_abaqus);
    }
    else if (case_num == 2)
    {
        // 执行自己的操作
    }
    else if(case_num == 101){// 隐式完整分析+重分析
        // 开始隐式分析
        // 材料属性赋值
        CAE::elastic_mat mat_item{2.1e5, 0.3, 7800};
        // 材料路径
        std::string path = "E:\\WH_CAE\\for_CA\\sample_model\\CA\\Job-cafull.inp";
        // 关键字
        string load_set_keyword = "Set-load";
        string load_value_keyword = "Cload";
        string dis_set_keyword = "Set-fix";
        // 建立CAE分析对象
        CAE::CAE_process cae_item(path, mat_item);
        // 读取计算文件
        cae_item.pre_info(load_set_keyword, load_value_keyword, dis_set_keyword);
        // 执行结构响应分析
        string result_path = "E:\\WH_CAE\\for_CA\\sample_model\\CA\\Im_Job-cafull.vtk";
        string path_abaqus = "E:\\WH_CAE\\for_CA\\sample_model\\CA\\Im_Job-cafull.txt";
        cae_item.implict_analysis(result_path, path_abaqus);

        // 开始重分析
        // 修改网格路径(仅包含1个part的网格)
        string mesh_path = "E:\\WH_CAE\\for_CA\\sample_model\\CA\\Job-cafull.inp";
        // 关键字
        string CA_del_set_keyword = "Set-del";
        // 建立重分析CAE分析对象
        cae_item.path_ = mesh_path;
        // 处理重复节点网格
        vector<int> del_topo;
        cae_item.CA_pre_process(CA_del_set_keyword, del_topo);
        // 执行重分析
        string CA_result_path = "E:\\WH_CAE\\for_CA\\sample_model\\CA\\CA_Job-cafull.vtk";
        string CA_path_abaqus = "E:\\WH_CAE\\for_CA\\sample_model\\CA\\CA_Job-cafull.txt";
        bool Is_Update = false;// 是否将原模型更新为修改后的模型
        cae_item.CA_ReAnalysis(CA_result_path, CA_path_abaqus, del_topo, Is_Update);// 提交完整隐式分析的data_managrment
    }
    else if(case_num == 102){
        // TODO 
        // 1、data_managrment输出为文件
        // 2、读文件+重分析
    }
    else if(case_num == 201){
        //显式动力学分析
    }
    else
    {
        std::cout << "Please check your code\n";
    }
}
