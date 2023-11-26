/**************************************************************************

Copyright:  WH team

Author: YinJichao <jichaoyinyjc@163.com>

Completion date:  XXX

Description: XXX

**************************************************************************/

#pragma once

#include "Eigen/Dense"
#include "include/elastic_mat.h"

typedef Eigen::Matrix<double, 4, 3> Matrix4d3;
typedef Eigen::Matrix<double, 8, 3> Matrix8d3;
typedef Eigen::Matrix<double, 12, 12> Matrix12d12;
typedef Eigen::Matrix<double, 24, 24> Matrix24d24;

namespace CAE
{
    class ele_base
    {
    public:
        // 构造函数，析构函数
        ele_base(){};

        // 赋值材料属性
        virtual void set_matrial(){};

        // 建立本构矩阵
        virtual void build_cons_mat(){};

        // 建立应变矩阵
        virtual void build_strain_mat(){};

        // 建立单元刚度矩阵
        void build_ele_stiff_mat(){};
    };
}