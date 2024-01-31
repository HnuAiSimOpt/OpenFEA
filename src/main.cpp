/**************************************************************************

Copyright:  WH team

Author: YinJichao <jichaoyinyjc@163.com>

Completion date:  XXX

Description: XXX

**************************************************************************/

#pragma once
#include <iostream>
#include "include/cae.h"

void main(int nargs, char* argv[])
{
    int case_num = 1;
    if (case_num == 1)
    {
        // 建立CAE分析对象
        CAE::CAE_process cae_item;

        //初始化
        cae_item.Init(nargs, argv);

        //求解
        cae_item.Solve();

    }
    else if (case_num == 2)
    {
        // 执行自己的操作
    }
    else
    {
        std::cout << "Please check your code\n";
    }
}
