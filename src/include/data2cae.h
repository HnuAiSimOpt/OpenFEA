#pragma once

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include "include/data_management.h"

using std::string;
using std::vector;

namespace CAE
{
    class ReadInfo
    {
    public:
        string path_;

    public:
        // 构造函数，析构函数
        ReadInfo();
        ReadInfo(string path) : path_(path){};
        // ~ReadInfo();

        // 读取单元、节点总数
        void read_ele_node_num(data_management &data_cae);

        // 基于指定符号划分字符串
        vector<string> split_str(const string str, const string part);

        // 读取节点坐标，节点拓扑关系，单元类型
        void read_geo_mesh(data_management &data_cae);
    };

    // // read boundary and load  (data: 2023 / 5 / 11)
    // void ReadBoundaryLoad(string path, IntArray& constrain, IntArray& load);

    // // read Abaqus VM_s and displacement (data: 2023 / 6 / 30)
    // void ReadAbaqus(DoubleMatrix& dis, DoubleMatrix& VM_s, string path);

}
