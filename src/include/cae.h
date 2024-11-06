/**************************************************************************

Copyright:  WH team

Author: YinJichao <jichaoyinyjc@163.com>

Completion date:  XXX

Description: XXX

**************************************************************************/

#pragma once

#include <math.h>
#include <iostream>
#include <ctime>
#include "./elastic_mat.h"
#include "./data_management.h"
#include "./data2cae.h"
#include "./set_bcs.h"
#include "./assemble.h"
#include "./nl_assemble.h"
#include "./data2vtk.h"
#include "./post_process.h"
#include "./explicit_tools.h"
#include "./linear_solution.h"
#include "./ca_reanalysis.h"
#include "./SFEM3D.h"

using namespace std;
namespace CAE
{
    struct controlPara
    {
        string work_path;            //exe文件所在路径
        string ifile_name;           //输入文件路径
        string ofile_name;           //输出文件路径
        int file_type = 1;          //文件类型  1:inp文件     2:fem文件     3:re文件

        int analysis_type = 1;      //求解类型  1:隐式        2:显式

        bool sfem_flag = false;      //光滑有限元标志
        bool reanalysis_flag = false;//重分析标志
        bool is_save_stiffness = true; // 重分析--是否保存隐式完整分析刚度
    };

    class CAE_process
    {
    public:
        controlPara  option_;
        string path_;
        elastic_mat mat_;
        data_management data_cae_;

    public:
        // 构造函数，析构函数
        CAE_process() {};
        CAE_process(string path) : path_(path) {};
        CAE_process(string path, elastic_mat mat) : path_(path), mat_(mat) {};
        // ~CAE_process();

    public://初始化
        void Init(int nargs, char* argv[]);
    public://初始化步骤
        void processCmdLine(int nargs, char* argv[]);
        void readFile();
        // 读取计算文件
        void pre_info(string load_set_keyword, string load_value_keyword, string dis_set_keyword);

    public://求解
        void Solve();
    public:
        // 执行结构响应分析
        void implict_analysis(string result_path, bool is_save_stiffness=false);
        // 读取重分析网格并做重复节点处理
        void CA_pre_process(string mesh_path, string node_now_path);
        // 执行重分析
        void CA_ReAnalysis(string result_path, int n_basis = 4, bool is_Update = false);
        // 执行光滑有限元
        void implict_SFEManalysis(string result_path, string path_abaqus);  
        // 执行结构动态响应分析
        void explicit_analysis(string result_path, string path_abaqus);
        // 转存刚度矩阵
        void Save_stiffness(assamble_stiffness &item_assam);    // 线性
    };
}
