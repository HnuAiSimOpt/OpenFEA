/**************************************************************************

Copyright:  WH team

Author: GuoDaozhen <397908710@qq.com>

Completion date:  XXX

Description: XXX

**************************************************************************/

#pragma once

#include "ele_base.h"

namespace CAE
{
    class SFEM2D;
    class SFEM_T3 : public ele_base
    {
    public:
        elastic_mat matrial_struc_;
        Matrix3d3 C_matrix_;

        SFEM2D* sfemData;

	    //节点积分域局部数据
	    int n = 5;//节点积分域编号

    public:
        // 构造函数，析构函数
        SFEM_T3() { type_ = "SFEM_T3"; nnode_ = 0; node_dof_ = 2; ngps_ = 0; };
        SFEM_T3(elastic_mat matrial_struc) : matrial_struc_(matrial_struc) {

        }
        // 建立单元刚度矩阵
         void build_ele_stiff_mat();
    };
}