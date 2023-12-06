/**************************************************************************

Copyright:  WH team

Author: ChenBinXiang <chen_mech99@163.com>

Completion date:  XXX

Description: XXX

**************************************************************************/

#pragma once
#include "Eigen/Dense"
#include "Eigen/SVD"
#include <vector>
#include "include/data_management.h"
#include "Elements/ele_base.h"
using std::vector;
using Eigen::MatrixXd;

namespace CAE {
	// ��ȡ��Ԫ�ڵ��������
	void get_ele_coords(const vector<int>& node_topos, const vector<vector<double>>& coords, MatrixXd& item_ele_coors, int nnode);

	// �����Զ�ʱ�䲽��
	void UpdateTimeStep(const vector<vector<int>>& node_topos, const vector<vector<double>>& coords, const vector<ele_base*>& ele_list,
						const vector<int>& ele_list_idx, double& time_step);

	// Ϊ��һ������������ر���
	void SwapData(vector<double>& disp_tp1, vector<double>& disp_t0, vector<double>& vel_tphalf, vector<double>& vel_thalf,
				  vector<double>& acc_t0, vector<double>& InFroce_, double& time_step, double& time_step_old, bool auto_time);

	// ���ʱ�䲽
	bool CheckTime(double& timestep, double limit);

	// ����ʵ������
	void update_coords(const vector<vector<double>>&  coords_, const vector<double>& disp_t0, vector<vector<double>>& real_coords_);

}// namespace
