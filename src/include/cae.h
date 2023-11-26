/**************************************************************************

Copyright:  WH team

Author: YinJichao <jichaoyinyjc@163.com>

Completion date:  XXX

Description: XXX

**************************************************************************/

#pragma once

#include <iostream>
#include "include/elastic_mat.h"
#include "include/data_management.h"
#include "include/data2cae.h"
#include "include/set_bcs.h"
#include "include/assemble.h"
#include "include/solver_linear.h"
#include "include/data2vtk.h"

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
    };
}
