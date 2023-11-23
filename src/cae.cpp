
#include "include/cae.h"

namespace CAE
{
    void CAE_process::pre_info()
    {
        ReadInfo item_info(path_);

        // 读取单元、节点总数
        item_info.read_ele_node_num(data_cae_);
        std::cout << "the number of element is: " << data_cae_.ne_ << std::endl
                  << "the number of node is: " << data_cae_.nd_ << std::endl;
    }
}