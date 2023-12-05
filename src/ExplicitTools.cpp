#include "include/ExplicitTools.h"

namespace CAE {
    // 获取单元节点坐标矩阵
    void CAE::get_ele_coords(const vector<int>& node_topos, const vector<vector<double>>& coords, MatrixXd& item_ele_coors,int nnode)
    {
        item_ele_coors.resize(nnode, 3);
        item_ele_coors.setZero();
        int  item_node;
        for (int i = 0; i < nnode; i++)
        {
            // 坐标
            item_node = node_topos[i] - 1;
            item_ele_coors(i, 0) = coords[item_node][0]; // X 坐标
            item_ele_coors(i, 1) = coords[item_node][1]; // Y 坐标
            item_ele_coors(i, 2) = coords[item_node][2]; // Z 坐标
        }
    }

    // 更新自动时间步长
    void CAE::UpdateTimeStep(const vector<vector<int>>& node_topos, const vector<vector<double>>& coords,
        const vector<ele_base*>& ele_list, const vector<int>& ele_list_idx, double& time_step)
    {
        for (int i = 0; i < ele_list_idx.size(); i++) {
            int idx = ele_list_idx[i];
            MatrixXd item_ele_coors;
            get_ele_coords(node_topos[i], coords, item_ele_coors, ele_list[idx]->nnode_);
            ele_list[idx]->update_timestep(item_ele_coors, time_step);
        }
    }

    // 为下一迭代步更新相关变量
    void CAE::SwapData(vector<double>& disp_tp1, vector<double>& disp_t0, vector<double>& vel_tphalf, vector<double>& vel_thalf,
        vector<double>& acc_t0, vector<double>& InFroce_, double& time_step, double& time_step_old, bool auto_time)
    {
        disp_tp1 = disp_t0;
        std::fill(disp_t0.begin(), disp_t0.end(), 0);
        vel_tphalf = vel_thalf;
        std::fill(vel_thalf.begin(), vel_thalf.end(), 0);
        std::fill(acc_t0.begin(), acc_t0.end(), 0);
        //std::fill(InFroce_.begin(), InFroce_.end(), 0);
        time_step_old = time_step;
        if (auto_time) {
            time_step = DBL_MAX;
        }
        else {
            time_step = 0.000001;
        }

    }

    // 检查时间步
    bool CheckTime(double& timestep, double limit)
    {
        if (timestep > limit) {
            timestep = limit;
            return 1;
        }
        return 0;
    }

    // 更新实际坐标
    void CAE::update_coords(const vector<vector<double>>& coords_, const vector<double>& disp_t0, vector<vector<double>>& real_coords_)
    {
        for (int i = 0; i < coords_.size(); i++) {
            real_coords_[i][0] = coords_[i][0] + disp_t0[3 * i];
            real_coords_[i][1] = coords_[i][1] + disp_t0[3 * i + 1];
            real_coords_[i][2] = coords_[i][2] + disp_t0[3 * i + 2];
        }
    }
}// namespace

