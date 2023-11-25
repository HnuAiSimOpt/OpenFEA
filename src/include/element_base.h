/**************************************************************************

Copyright:  WH team

Author: YinJichao <jichaoyinyjc@163.com>

Completion date:  XXX

Description: XXX

**************************************************************************/

#pragma once

#include "include/elastic_mat.h"

namespace CAE
{
    class ele_base
    {
    public:
        // 构造函数，析构函数
        ele_base(){};

        // 赋值材料属性
        virtual void set_matrial() = 0;

        // 建立本构矩阵
        virtual void build_cons_mat() = 0;

        // 建立应变矩阵
        virtual void build_strain_mat() = 0;
    };
}