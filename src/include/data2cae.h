/**************************************************************************

Copyright:  WH team

Author: YinJichao <jichaoyinyjc@163.com>

Completion date:  2023/11/23

Description: read a calculation file

**************************************************************************/

#pragma once

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include "include/data_management.h"

using std::string;
using std::vector;
using std::cout;
using std::endl;

namespace CAE
{
    class ReadInfo
    {
    public:
        string path_;

    public:
        // 构造函数，析构函数
        ReadInfo(){};
        ReadInfo(string path) : path_(path){};
        // ~ReadInfo();

        // 读取单元、节点总数
        void read_ele_node_num(data_management &data_cae);

        // 基于指定符号划分字符串
        vector<string> split_str(const string str, const string part);

        // 读取节点坐标，节点拓扑关系，单元类型
        void read_geo_mesh(data_management &data_cae);

        // 读取载荷信息 
        void read_load_bcs(string load_set_keyword, string load_value_keyword, data_management &data_cae);

        // 读取约束信息
        void read_dis_bcs(string dis_set_keyword, data_management &data_cae);
    };

}
