/**************************************************************************

Copyright:  WH team

Author: YinJichao <jichaoyinyjc@163.com>��ChenBinXiang <chen_mech99@163.com>

Completion date:  XXX

Description: XXX

**************************************************************************/

#pragma once

#include "ele_base.h"

namespace CAE
{
    class tetra_ele_elastic : public ele_base
    {
    public:
        double det_jacobi_; // �ſɱȾ�������ʽ(���Ļ��ֵ�)
        elastic_mat matrial_struc_;
        Matrix6d6 C_matrix_;
        Matrix6d12 strain_mat;

    public:
        // ���캯������������
        tetra_ele_elastic() { type_ = "C3D4"; node_dof_ = 3; nnode_ = 4; ngps_ = 1;};
        tetra_ele_elastic(elastic_mat matrial_struc) : matrial_struc_(matrial_struc) { type_ = "C3D4"; node_dof_ = 3; nnode_ = 4; ngps_ = 1;};

        // ���ϸ�����
        void set_matrial(elastic_mat matrial_struc)override { matrial_struc_ = matrial_struc; }  ;

        // ������������
        virtual void build_cons_mat();

        // ����Ӧ�����(���ֵ�)
        virtual void build_strain_mat(Eigen::Ref<Eigen::MatrixXd> node_coords, Eigen::Ref<Eigen::MatrixXd> strain_mat);


        // ������Ԫ�նȾ���
        void build_ele_stiff_mat(Eigen::Ref<Eigen::MatrixXd> node_coords, Eigen::Ref<Eigen::MatrixXd> stiffness_matrix) override;

        // ������Ԫ��������
        void build_ele_mass(const vector<int>& node_topos, const vector<vector<double>>& coords, vector<double>& Mass)override;

        // ���㵥Ԫ����
        void cal_in_force(const vector<int>& node_topos, const vector<vector<double>>& real_coords, const vector<double>& disp_d,
                          vector<double>& stress, vector<double>& strain, vector<double>& InFroce)override;
        
        // ���㵥Ԫʱ�䲽��
        void update_timestep(Eigen::Ref<Eigen::MatrixXd> node_coords, double& time_step) override;
    };


}