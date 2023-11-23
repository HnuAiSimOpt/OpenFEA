
#pragma once

#include <iostream>
#include "include/mat.h"
#include "include/data2cae.h"
#include "include/data_management.h"

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
        CAE_process();
        CAE_process(string path, elastic_mat mat) : path_(path), mat_(mat){};
        // ~CAE_process();

        // 读取计算文件
        void pre_info();
    };
}
