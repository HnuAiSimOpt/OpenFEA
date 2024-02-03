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
    class CQUAD4 : public ele_base
    {
    public:
        // 构造函数，析构函数
        CQUAD4()
        {
            type_ = "CQUAD4";
            nnode_ = 4;
            node_dof_ = 6;
            ngps_ = 4;
        };
        CQUAD4(elastic_mat matrial_struc) : matrial_struc_(matrial_struc)
        {
            type_ = "CQUAD4";
            nnode_ = 4;
            node_dof_ = 6;
            ngps_ = 4;
        };

        // 材料赋属性
        void set_matrial(elastic_mat matrial_struc) override { matrial_struc_ = matrial_struc; };

        // 建立本构矩阵

        void build_cons_mat(double D[8][8]);

        // 建立单元刚度矩阵
        void build_ele_stiff_mat(Eigen::Ref<Eigen::MatrixXd> node_coords, Eigen::Ref<Eigen::MatrixXd> stiffness_matrix) override;

        void computeLocalSystems(Eigen::Ref<Eigen::MatrixXd> node_coords);

        void KmPart(Eigen::Ref<Eigen::MatrixXd> node_coords, double D[8][8], double ke[24][24]);

        void KpPart(double D[8][8], double k0[24][24]);

        void caculateKB(double* A, double* B, double AREA, double A4, double A42, double D[8][8], double KB[9][9], double KS[9][9], double& PHI_SQ);

        void STATIC_CONDENSATION(double KM_QQ_5[15][15], double B2M_QQ_5[3][15], double B3M_QQ_5[3][15], double k0[24][24]);

        void Stress(vector<double>& UEL);

   

    public:
        double t = 1;//壳单元厚度，后续建议放置在Section类中

        elastic_mat matrial_struc_;
        Matrix6d6 C_matrix_;
    
        //子三角形及其他索引（所有CQUAD4单元公用）
        static int IDI[4][9]; 
        static int IDM[12];
        static int IDV[9];
        static int kmIndex[8];

        double XSD[4];
        double YSD[4];


        double HBAR;                //判断膜部分刚度矩阵计算策略

        double BE1[2][4];           //应变计算矩阵(不存储零元素,实际对应自由度1,2,7,8,     13,14,19,20)
        double BE2[3][12];          //应变计算矩阵(不存储零元素，实际对应自由度3,4,5,9,10,11,15,16,17,21,22,23)
        double BE3[3][12];          //应变计算矩阵(不存储零元素，实际对应自由度3,4,5,9,10,11,15,16,17,21,22,23)

        double PHI_SQ_TRIA[4];      //四个三角形单元的PHI_SQ
        double PHI_SQ = 0.0;        //四个三角形单元的PHI_SQ平均值

        double AREA_TRIA[4];        //四个三角形的面积
        double AREA_QUAD = 0.0;     //四边形面积

        vector<Point> lcoords;

        Mat3 TEG;//单元局部坐标转换矩阵



    };
}