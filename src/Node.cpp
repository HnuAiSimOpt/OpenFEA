#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include "include/Node.h"

namespace CAE
{
    //预处理节点编号和自由度
    void node_base::pre_nodes(data_management &data_cae)
    {
        data_cae.nodes_.resize(data_cae.nd_);
        for ( int i = 0 ; i < data_cae.nd_  ; i++ )
        {
            // node_id 是 int 类型，直接赋值
            data_cae.nodes_[i].node_id = (i + 1);
            // 设置 full_dof 
            data_cae.nodes_[i].full_dof = { 3*i , 3*i  + 1, 3*i + 2 };
        }
    }
    //计算无约束节点自由度总数
    void node_base::count_re_free_dof_num(data_management &data_cae)
    {
        int totalcount = 0;
        for(const auto& node : data_cae.nodes_ )
        {
            for (int value : node.re_free_dof)
            {
                if (value > -1){
                    totalcount++;
                }
            }
        }

        data_cae.re_free_dof_num = totalcount;
        std::cout << "the num of free nodes dofs is :" << data_cae.re_free_dof_num<< std::endl;
    }
}
