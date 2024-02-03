/**************************************************************************

Copyright:  WH team

Author: YinJichao <jichaoyinyjc@163.com>、ChenBinXiang <chen_mech99@163.com>

Completion date:  XXX

Description: XXX

**************************************************************************/

#pragma once

#include "ele_base.h"
#include "../../include/SupportForSFEM.h"
namespace CAE
{
    class CTRIA3 : public ele_base
    {
    public:
        // 构造函数，析构函数
        CTRIA3()
        {
            type_ = "CTRIA3";
            nnode_ = 3;
            node_dof_ = 6;
            ngps_ = 1;
        };
        CTRIA3(elastic_mat matrial_struc) : matrial_struc_(matrial_struc)
        {
            type_ = "CTRIA3";
            nnode_ = 3;
            node_dof_ = 6;
            ngps_ = 1;
        };

        // 材料赋属性
        void set_matrial(elastic_mat matrial_struc) override { matrial_struc_ = matrial_struc; };

        // 建立本构矩阵
        void build_cons_mat(double D[8][8]);

        // 建立单元刚度矩阵
        void build_ele_stiff_mat(Eigen::Ref<Eigen::MatrixXd> node_coords, Eigen::Ref<Eigen::MatrixXd> stiffness_matrix) override;

        void computeLocalSystems(Eigen::Ref<Eigen::MatrixXd> node_coords);

        void Ke(Eigen::Ref<Eigen::MatrixXd> node_coords, double D[8][8], double ke[18][18]);

        void Stress(vector<double>& UEL);

    public:
        double t = 1;//壳单元厚度，后续建议放置在Section类中

        elastic_mat matrial_struc_;
        Matrix6d6 C_matrix_;
    
        double BE1[3][5] = { 0. };/*应变计算矩阵(不存储零元素, 实际对应自由)
                                      |  1  |     |  7  |     |     |
                                      |     |  2  |     |  8  |  14 |
                                      |  1  |  2  |  7  |  8  |  13 | */

        double BE2[3][6] = { 0. };/*应变计算矩阵(不存储零元素，实际对应自由度
                                      |  4  |  5  |  10 |  11 |  16 |  17 |
                                      |  4  |  5  |  10 |  11 |  16 |  17 |
                                      |  4  |  5  |  10 |  11 |  16 |  17 |*/



        double BE3[2][9] = { 0. };/*应变计算矩阵(不存储零元素，实际对应自由度
                                      |  3  |  4  |  5  |  9  |  10 |  11 |  15 |  16 |  17 |
                                      |  3  |  4  |  5  |  9  |  10 |  11 |  15 |  16 |  17 |*/


        double PHI_SQ = 0.0;        

        double area;                //三角形面积

        vector<Point> lcoords;

        Mat3 TEG;//单元局部坐标转换矩阵
    };
}