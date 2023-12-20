
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
        //MatrixXd GP1, GP2;

        //vector<MatrixXd> P_GP1;//���ռ���ֵ�
        //vector<MatrixXd> P_GP2;//���ռ���ֵ�
        //vector<double> W_1;//Ȩ��
        //vector<Eigen::Vector3d> Normal;//������






    public:
        
        struct LagrangeBR
        {
            MatrixXd N_out;
            MatrixXd dNdxi_out;
            
        };

       // void PhySpaceGPs(data_management& data_cae, elastic_mat& data_mat);
        void GetIntF_face_Inform(data_management& data_cae, MatrixXd& nodes1,
            int& face_nodes, int& e);
        void GetIntF_ele_Inform(MatrixXd& pts1, MatrixXd& pts2, vector<int>& F_eper_dof,
            vector<int>& C_eper_dof, data_management& data_cae, int& e, 
            int& n_node_F, int& n_node_C);
        MatrixXd GlobalMap3D(MatrixXd gpoint, MatrixXd nodes, int& n_node_mesh);
        LagrangeBR lagrange_basis(MatrixXd& coord, int& n_node_mesh);
        //������նȾ���
        void InterfacialStifMatrix(data_management& data_cae, elastic_mat& data_mat,
            vector<double>& nz_val, vector<int>& row_idx, vector<int>& col_idx);
        //Ce����
        void Get_Ce(elastic_mat& data_mat, MatrixXd& Ce);
        //����B��Nm
        void Calculate_B_Nm(MatrixXd& p_gps, MatrixXd& nodes, MatrixXd& B, MatrixXd& Nm, int& n_node_face);
        //����������
        double Get_Alpha(data_management& data_cae, elastic_mat& data_mat, MatrixXd& pts1);
        //����Kp��Kd
        void Calculate_Kp_Kd(MatrixXd& Kp11, MatrixXd& Kp12, MatrixXd& Kp22,
            MatrixXd& Kd11, MatrixXd& Kd12, MatrixXd& Kd21, MatrixXd& Kd22,
            MatrixXd& Nm1, MatrixXd& Nm2, MatrixXd& B1, MatrixXd& B2,
            MatrixXd& n, MatrixXd& Ce, double& alpha, double& wt1);
        //���㵥Ԫ�ڵ����ɶ�
       // void Get_ele_dof(vector<int>& F_eper_dof, vector<int>& C_eper_dof,
            //data_management& data_cae, int& e, int& n_node_F, int& n_node_C);
        //�����ܽ���նȾ���
        void Calculate_InterFMatrix(MatrixXd& K11, MatrixXd& K12,MatrixXd& K21, MatrixXd& K22,
            MatrixXd& Kd11, MatrixXd& Kd12, MatrixXd& Kd21, MatrixXd& Kd22,
            MatrixXd& Kp11, MatrixXd& Kp12, MatrixXd& Kp22, int& n_node_F, int& n_node_C);
        //�������նȾ���
        void Fill_InterFMatrix(vector<int>& j_eper_dof, vector<int>& i_eper_dof,
            vector<double>& nz_val, vector<int>& row_idx, vector<int>& col_idx,
            MatrixXd& K_interface);
    };
}
