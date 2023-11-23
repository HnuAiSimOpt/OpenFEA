
#include "include/cae.h"

namespace CAE
{
    void CAE_process::pre_info()
    {
        ReadInfo item_info(path_);

        // 读取单元、节点总数
        item_info.read_ele_node_num(data_cae_);

        // 读取几何信息
        item_info.read_geo_mesh(data_cae_);
    }
}