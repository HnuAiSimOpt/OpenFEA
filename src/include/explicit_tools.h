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
#include "elements/include/ele_base.h"
using std::vector;
using Eigen::MatrixXd;

namespace CAE {
	// 获取单元节点坐标矩阵
	void get_ele_coords(const vector<int>& node_topos, const vector<vector<double>>& coords, MatrixXd& item_ele_coors, int nnode);

	// 更新自动时间步长
	void UpdateTimeStep(const vector<vector<int>>& node_topos, const vector<vector<double>>& coords, const vector<ele_base*>& ele_list,
						const vector<int>& ele_list_idx, map<int, int>& ele_map_list, double& time_step);

	// 为下一迭代步更新相关变量
	void SwapData(vector<double>& disp_tp1, vector<double>& disp_t0, vector<double>& vel_tphalf, vector<double>& vel_thalf,
				  vector<double>& acc_t0, vector<double>& InFroce_, double& time_step, double& time_step_old, bool auto_time);

	// 检查时间步
	bool CheckTime(double& timestep, double limit);

	// 更新实际坐标
	void update_coords(const vector<vector<double>>&  coords_, const vector<double>& disp_t0, vector<vector<double>>& real_coords_);

}// namespace
