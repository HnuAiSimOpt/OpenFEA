/**************************************************************************

Copyright:  WH team

Author: YinJichao <jichaoyinyjc@163.com>

Completion date:  2023/11/23

Description: read a calculation file

**************************************************************************/

#pragma once

#include <sstream>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
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
        string path_now_ndoe_;

    public:
        // 构造函数，析构函数
        ReadInfo(){};
        ReadInfo(string path1) : path_(path1){};
        ReadInfo(string path1, string path2) : path_(path1), path_now_ndoe_(path2){};
        // ~ReadInfo();

        // 读取单元、节点总数
        void read_ele_node_num(data_management &data_cae);

        // 基于指定符号划分字符串
        vector<string> split_str(const string str, const string part);

        // 读取节点坐标，节点拓扑关系，单元类型
        void read_geo_mesh(data_management &data_cae);

        //读取非协调信息
        void readNconformingMessage(data_management& data_cae);

        // 读取载荷信息 
        void read_load_bcs(string load_set_keyword, string load_value_keyword, data_management &data_cae);

        // 读取约束信息
        void read_dis_bcs(string dis_set_keyword, data_management &data_cae);

        // 读取材料类型
        void read_mat(elastic_mat& mat_);

        // 读取时间步长信息
        void read_time(data_management& data_cae_);

        void del_blank(string& str);
        int del_blank(string &str, map<string, int> &ele_map);

        // 重分析处理数据
        void CA_read_change(data_management &data_cae);
        void CA_read_now_node(data_management &data_cae);
        void CA_data_convert(data_management& data_cae, vector<int>& del_topo, bool* node_del_idx, bool* topo_del_idx);

        // 统一处理跨平台不同行尾序的getline
        std::istream& safeGetline(std::istream& is, std::string& t);
    };

}
