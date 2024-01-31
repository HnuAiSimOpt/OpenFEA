/**************************************************************************

Copyright:  WH team

Author: YinJichao <jichaoyinyjc@163.com>

Completion date:  XXX

Description: XXX

**************************************************************************/

#pragma once

#include <iostream>
#include "elastic_mat.h"
#include "data_management.h"
#include "data2cae.h"
#include "set_bcs.h"
#include "assemble.h"
#include "data2vtk.h"
#include "explicit_tools.h"
#include "linear_solution.h"


using namespace std;
namespace CAE
{
        struct controlPara
    {
        string work_path;            //exe文件所在路径
        string file_name;           //文件路径
        int file_type = 1;          //文件类型  1:inp文件     2:fem文件

        int analysis_type = 1;      //求解类型  1:隐式        2:显示

        bool sfem_flag = true;      //光滑有限元标志
        bool reanalysis_flag = false;//重分析标志
    };

    class CAE_process
    {
    public:
        controlPara  option;
        string path_;
        elastic_mat mat_;
        data_management data_cae_;

    public:
        // 构造函数，析构函数
        CAE_process(){};
        CAE_process(string path, elastic_mat mat) : path_(path), mat_(mat){};
        // ~CAE_process();


    public://初始化
        void Init(int nargs, char* argv[]);
        private://初始化步骤
            void processCmdLine(int nargs, char* argv[]);
            void readFile();

    public://求解
        void Solve();
        private:
            void implict_analysis(string result_path, string path_abaqus);      //隐式分析
            void implict_SFEManalysis(string result_path, string path_abaqus);  //光滑有限元
            void explicit_analysis(string result_path, string path_abaqus);     //显式分析

    };
}
