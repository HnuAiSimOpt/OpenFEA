/**************************************************************************

Copyright:  WH team

Author: YinJichao <jichaoyinyjc@163.com>

Completion date:  XXX

Description: XXX

**************************************************************************/

#pragma once

namespace CAE
{
    class stress
    {
    public:
        // 构造函数，析构函数
        stress(){};
        // stress(string path, elastic_mat mat) : path_(path), mat_(mat){};
        void get_vms();
    };
}