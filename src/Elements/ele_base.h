/**************************************************************************

Copyright:  WH team

Author: YinJichao <jichaoyinyjc@163.com>��ChenBinXiang <chen_mech99@163.com>

Completion date:  XXX

Description: XXX

**************************************************************************/

#pragma once

#include "Eigen/Dense"
#include "Eigen/SVD"
#include <vector>
#include "../include/elastic_mat.h"
#include "../include/Factory.h"

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

typedef Eigen::Matrix<double, 24, 24> Matrix24d24;
typedef Eigen::Matrix<double, 24, 6> Matrix24d6;

using std::vector;

namespace CAE
{
    class ele_base
    {
    public:
        // ���캯������������
        ele_base() {};
    public:
        std::string type_;//��Ԫ��������
        int nnode_;//�õ�Ԫӵ�нڵ�����
        int node_dof_;//�õ�Ԫÿ���ڵ����ɶ���
        // ��ֵ��������
        virtual void set_matrial(elastic_mat data_mat) {};

        // ������������
        virtual void build_cons_mat() {};

        // ����Ӧ�����
        virtual void build_strain_mat() {};

        // ������Ԫ�նȾ���
        virtual void build_ele_stiff_mat(Eigen::MatrixXd& node_coords, Eigen::MatrixXd& stiffness_matrix) {};
    };
    CREAT_FACTORY(ele_base);
}