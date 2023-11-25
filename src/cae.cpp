/**************************************************************************

Copyright:  WH team

Author: YinJichao <jichaoyinyjc@163.com>

Completion date:  XXX

Description: XXX

**************************************************************************/

#include "include/cae.h"

namespace CAE
{
    void CAE_process::pre_info(string load_set_keyword, string load_value_keyword, string dis_set_keyword)
    {
        ReadInfo item_info(path_);

        // 读取单元、节点总数
        item_info.read_ele_node_num(data_cae_);

        // 读取几何信息
        item_info.read_geo_mesh(data_cae_);

        // 读取载荷边界信息
        item_info.read_load_bcs(load_set_keyword, load_value_keyword, data_cae_);

        // 读取位移边界信息
        item_info.read_dis_bcs(dis_set_keyword, data_cae_);
    }

    // 执行结构响应分析
    void CAE_process::implict_analysis()
    {
        set_BCs item_bcs;

        // 设置边界条件
        item_bcs.build_free_index(data_cae_);

        // 建立单载荷向量
        item_bcs.build_single_load(data_cae_);

        // 组装刚度矩阵
        assamble_stiffness item_assam;
        item_assam.build_CSR(data_cae_);
        item_assam.fill_CSR_sparse_mat(data_cae_);
    }
}