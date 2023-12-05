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
    class hex_ele_elastic : public ele_base
    {
    public:
        double det_jacobi_; // �ſɱȾ�������ʽ�����Ļ��ֵ㣩
        elastic_mat matrial_struc_;
        Matrix6d6 C_matrix_;

    public:
        // ���캯������������
        hex_ele_elastic() { type_ = "C3D8R"; nnode_ = 8; node_dof_ = 3; ngps_ = 8; };
        hex_ele_elastic(elastic_mat matrial_struc) : matrial_struc_(matrial_struc) { type_ = "C3D8R"; nnode_ = 8; node_dof_ = 3; ngps_ = 8;};

        // ���ϸ�����
        void set_matrial(elastic_mat matrial_struc) override { matrial_struc_ = matrial_struc; };

        // ������������
        virtual void build_cons_mat();

        // ����Ӧ�����(���ֵ�)
        virtual void build_strain_mat(Eigen::Ref<Eigen::MatrixXd> node_coords, Matrix6d24& strain_mat, vector<double>& gp_points, double* det_jacobi_point);

        // ������Ԫ�նȾ���
        void build_ele_stiff_mat(Eigen::Ref<Eigen::MatrixXd> node_coords, Eigen::Ref<Eigen::MatrixXd> stiffness_matrix) override;

        // ������Ԫ�ܶȾ���
        void build_ele_mass(const vector<int>& node_topos, const vector<vector<double>>& coords, vector<double>& Mass)override;

        // �����κ���
        void build_shape_fun(Eigen::Ref<Eigen::MatrixXd> node_coords, Eigen::MatrixXd& shape_fun, vector<double>& gp_points, double& det_jacobi_point);
        
        // ���㵥Ԫ����
        void cal_in_force(const vector<int>& node_topos, const vector<vector<double>>& real_coords, const vector<double>& disp_d,
            vector<double>& stress, vector<double>& strain, vector<double>& InFroce)override;

        // ���㵥Ԫʱ�䲽��
        void update_timestep(Eigen::Ref<Eigen::MatrixXd> node_coords, double& time_step) override;
    };
}