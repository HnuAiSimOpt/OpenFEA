/**************************************************************************

Copyright:  WH team

Author: YinJichao <jichaoyinyjc@163.com>

Completion date:  XXX

Description: XXX

**************************************************************************/

#pragma once

#include <fstream>
#include <iostream>
#include <sstream>
#include <map>
#include "data_management.h"

using std::cout;
using std::endl;
using std::map;

namespace CAE
{
    class data_process
    {
    public:
        map<string, int> ELE_TYPES = {{"C3D4", 1}, {"C3D8", 2}, {"C3D8R", 3}}; // 建立单元类型到整型的映射
    public:
        // 重排序位移，填充为完整自由度位移
        void fill_full_dis(data_management &data_cae);

        // 读取Abaqus位移场
        void read_abaqus_dis(string path, vector<double> &abaqus_dis);

        // 导出 位移场 VTK
        void export_dis_2_vtk(data_management &data_cae, string path, double scale_dis, string path_abaqus = " ", bool read_flag = false);
    };
}