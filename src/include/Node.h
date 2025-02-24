#include <iostream>
#include <string>
#include <vector>
#include "include/data_management.h"

using namespace std;

namespace CAE
{
    // 基础节点类
    class node_base
    {
    public:

        // 构造函数，析构函数
        node_base(){};

        // 预处理节点编号、自由度
        void pre_nodes( data_management &data_cae );
        
        ////计算无约束节点自由度总数
        void count_re_free_dof_num(data_management &data_cae );
    };

} // namespace CAE
