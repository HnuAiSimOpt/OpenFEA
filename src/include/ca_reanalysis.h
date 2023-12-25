/**************************************************************************

Copyright:  WH team

Author: YinJichao <jichaoyinyjc@163.com>

Completion date:  XXX

Description: XXX

**************************************************************************/

#include <algorithm>
#include "include/data_management.h"

namespace CAE
{
    // 重分析信息预处理
    void ca_pre_process(data_management &data_cae);

    // 建立
    void ca_build_rom(data_management &data_cae);

    // 找二维向量最大值
    int find_max_element(vector<vector<int>> &vec);
}