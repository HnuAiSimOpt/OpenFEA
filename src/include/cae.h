/**************************************************************************

Copyright:  WH team

Author: YinJichao <jichaoyinyjc@163.com>

Completion date:  XXX

Description: XXX

**************************************************************************/

#pragma once

#include <iostream>
#include <ctime>
#include "./elastic_mat.h"
#include "./data_management.h"
#include "./data2cae.h"
#include "./set_bcs.h"
#include "./assemble.h"
#include "./data2vtk.h"
#include "./explicit_tools.h"
#include "./linear_solution.h"
#include "./ca_reanalysis.h"

using namespace std;
namespace CAE
{
    class CAE_process
    {
    public:
        string path_;
        elastic_mat mat_;
        data_management data_cae_;

    public:
        // 构造函数，析构函数
        CAE_process(){};
        CAE_process(string path, elastic_mat mat) : path_(path), mat_(mat){};
        // ~CAE_process();

        // 读取计算文件
        void pre_info(string load_set_keyword, string load_value_keyword, string dis_set_keyword);

        // 执行结构响应分析
        void implict_analysis(string result_path, string path_abaqus);

        // 读取重分析网格并做重复节点处理
        void CA_pre_process(data_management& data_cae_Im);

        // 执行重分析
        void CA_ReAnalysis(data_management& data_cae_Im, string result_path, string path_abaqus);

        // 执行结构动态响应分析
        void explicit_analysis(string result_path, string path_abaqus);
    };
}
