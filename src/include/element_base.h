///**************************************************************************
//
//Copyright:  WH team
//
//Author: YinJichao <jichaoyinyjc@163.com>
//
//Completion date:  XXX
//
//Description: XXX
//
//**************************************************************************/
//
//#pragma once
//
//#include "Eigen/Dense"
//#include "include/elastic_mat.h"
//
//using Eigen::MatrixXd;
//
//namespace CAE
//{
//    class ele_base
//    {
//    public:
//        // 构造函数，析构函数
//        ele_base(){};
//    public:
//        std::string type_;
//        // 赋值材料属性
//        virtual void set_matrial(){};
//
//        // 建立本构矩阵
//        virtual void build_cons_mat(){};
//
//        // 建立应变矩阵
//        virtual void build_strain_mat(){};
//
//        // 建立单元刚度矩阵
//        virtual void build_ele_stiff_mat(MatrixXd &node_coords, MatrixXd &stiffness_matrix){};
//    };
//}