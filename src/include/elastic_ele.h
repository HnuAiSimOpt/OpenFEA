/**************************************************************************

Copyright:  WH team

Author: YinJichao <jichaoyinyjc@163.com>

Completion date:  XXX

Description: XXX

**************************************************************************/

#pragma once

#include <vector>
#include "Eigen/Dense"
#include "include/element_base.h"

typedef Eigen::Matrix<double, 3, 3> Matrix3d3;
typedef Eigen::Matrix<double, 3, 4> Matrix3d4;
typedef Eigen::Matrix<double, 3, 8> Matrix3d8;

typedef Eigen::Matrix<double, 4, 3> Matrix4d3;

typedef Eigen::Matrix<double, 6, 6> Matrix6d6;
typedef Eigen::Matrix<double, 6, 12> Matrix6d12;
typedef Eigen::Matrix<double, 6, 24> Matrix6d24;

typedef Eigen::Matrix<double, 8, 3> Matrix8d3;

typedef Eigen::Matrix<double, 12, 6> Matrix12d6;
typedef Eigen::Matrix<double, 12, 12> Matrix12d12;
typedef Eigen::Matrix<double, 24, 6> Matrix24d6;
typedef Eigen::Matrix<double, 24, 24> Matrix24d24;

using std::vector;

namespace CAE
{
    //****************************************************************************//
    //                                四面体单元                                   //
    //****************************************************************************//
    class tetra_ele_elastic : public ele_base
    {
    public:
        double det_jacobi_; // 雅可比矩阵行列式(中心积分点)
        elastic_mat matrial_struc_;
        Matrix6d6 C_matrix_;
        Matrix12d12 ele_stiff_;

    public:
        // 构造函数，析构函数
        tetra_ele_elastic(){};
        tetra_ele_elastic(elastic_mat matrial_struc) : matrial_struc_(matrial_struc){};

        // 材料赋属性
        virtual void set_matrial(elastic_mat matrial_struc) { matrial_struc_ = matrial_struc; };

        // 建立本构矩阵
        virtual void build_cons_mat();

        // 建立应变矩阵(积分点)
        virtual void build_strain_mat(Matrix4d3 &node_coords, Matrix6d12 &strain_mat);

        // 建立单元刚度矩阵
        void build_ele_stiff_mat(Matrix4d3 &node_coords);

        // 建立单元密度矩阵
        // void build_ele_den_mat();
    };

    //****************************************************************************//
    //                                六面体单元                                   //
    //****************************************************************************//
    class hex_ele_elastic : public ele_base
    {
    public:
        double det_jacobi_; // 雅可比矩阵行列式（中心积分点）
        elastic_mat matrial_struc_;
        Matrix6d6 C_matrix_;
        Matrix24d24 ele_stiff_;

    public:
        // 构造函数，析构函数
        hex_ele_elastic(){};
        hex_ele_elastic(elastic_mat matrial_struc) : matrial_struc_(matrial_struc){};

        // 材料赋属性
        virtual void set_matrial(elastic_mat matrial_struc) { matrial_struc_ = matrial_struc; };

        // 建立本构矩阵
        virtual void build_cons_mat();

        // 建立应变矩阵(积分点)
        virtual void build_strain_mat(Matrix8d3 &node_coords, Matrix6d24 &strain_mat, vector<double> &gp_points, double *det_jacobi_point);

        // 建立单元刚度矩阵
        void build_ele_stiff_mat(Matrix8d3 &node_coords);

        // 建立单元密度矩阵
        // void build_ele_den_mat();
    };
}