/**************************************************************************

Copyright:  WH team

Author: YinJichao <jichaoyinyjc@163.com>

Completion date:  XXX

Description: XXX

**************************************************************************/

#pragma once

#include <iostream>
#include <algorithm>
#include "include/data_management.h"

using std::sort;
using std::cout;

namespace CAE
{
    class set_BCs
    {
    public:
        // 构造函数，析构函数
        set_BCs(){};
        // ~CAE_process();

        // 建立无约束索引，被约束自由度置为-1，无约束自由度重新顺序编号
        void build_free_index(data_management &data_cae);

        // 基于重排的自由度建立载荷向量
        void build_single_load(data_management &data_cae);
    };
}