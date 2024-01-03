
/**************************************************************************

Copyright:  WH team

Author: Zhumingjun <1765380405@qq.com>

Completion date:  XXX

Description: XXX

**************************************************************************/

#pragma once

#include "Eigen/Dense"
#include "Eigen/SVD"
#include <vector>
#include "include/data_management.h"
#include <iostream>
using std::vector;
using Eigen::MatrixXd;


using namespace std;
namespace CAE
{
    class NCF_map
    {
    public:
        MatrixXd GP1, GP2;



    public:
        
        struct LagrangeBR
        {
            MatrixXd N_out;
            MatrixXd dNdxi_out;
            
        };

        void PhySpaceGPs(data_management& data_cae, elastic_mat& data_mat);
        void GetIntF_ele_Inform(data_management& data_cae, vector<int>& sctr1, vector<int>& sctr2,
            MatrixXd& pts1, MatrixXd& pts2, int& e);
        MatrixXd GlobalMap3D(MatrixXd gpoint, MatrixXd nodes);
        LagrangeBR lagrange_basis(MatrixXd& coord);
        //填充界面刚度矩阵
        void InterfacialStifMatrix(data_management& data_cae, elastic_mat& data_mat,
            vector<double>& nz_val, vector<int>& row_idx, vector<int>& col_idx);
        //罚参数计算
        double Get_Alpha(data_management& data_cae, elastic_mat& data_mat, MatrixXd& pts1);
        //计算Kp，Kd
        void Calculate_Kp_Kd(MatrixXd& Kp11, MatrixXd& Kp12, MatrixXd& Kp22,
            MatrixXd& Kd11, MatrixXd& Kd12, MatrixXd& Kd21, MatrixXd& Kd22,
            MatrixXd& Nm1, MatrixXd& Nm2, MatrixXd& B1, MatrixXd& B2,
            MatrixXd& n, MatrixXd& Ce, double& alpha, double& wt1);
        //计算总界面刚度矩阵
        void Calculate_InterFMatrix(MatrixXd& K11, MatrixXd& K12,
            MatrixXd& K21, MatrixXd& K22, vector<int>& F_eper_dof, vector<int>& C_eper_dof,
            MatrixXd& Kd11, MatrixXd& Kd12, MatrixXd& Kd21, MatrixXd& Kd22,
            MatrixXd& Kp11, MatrixXd& Kp12, MatrixXd& Kp22,
            data_management& data_cae, int& e);
        //储存界面刚度矩阵
        void Fill_InterFMatrix(vector<int> j_eper_dof, vector<int> i_eper_dof,
            vector<double>& nz_val, vector<int> row_idx, vector<int>& col_idx, MatrixXd& K_interface);
    };
}
