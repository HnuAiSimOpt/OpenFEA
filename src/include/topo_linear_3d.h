/**************************************************************************

Copyright:  WH team

Author: YinJichao <jichaoyinyjc@163.com>

Completion date:  XXX

Description: XXX

**************************************************************************/
#pragma once

#define DMAX(a, b) (a > b ? a : b)
#define DMIN(a, b) (a < b ? a : b)

#include <cmath>
#include <vector>
#include <string>
#include <Eigen/Dense>
#include "./utils.h"
#include "assemble.h"
#include "./set_bcs.h"
#include "./elastic_mat.h"
#include "./data_management.h"
#include "./data2vtk.h"
#include "./post_process.h"
#include "solver/include/solver_superlu.h"

using namespace std;
namespace CAE
{
    class TopoLinear3D
    {
    private:
        elastic_mat mat_;
        double vol_;
        double rmin_;
        double penal_;
        double move_ = 0.1;
        double e_min_ = 1.0E-6;

    private:
        int nele_;

    public:
        // 构造函数，析构函数
        TopoLinear3D() {};
        TopoLinear3D(double vol, double rmin, double penal, elastic_mat mat) : vol_(vol), rmin_(rmin), penal_(penal), mat_(mat) {};

        // 拓扑优化循环
        double topo_linear_cycle_3d(data_management &data_cae, string result_path, int loop_max = 80);

        // 计算过滤矩阵
        void sen_filter(SpMatrix &H, vector<double> &Hs, vector<vector<double>> &coords, vector<vector<int>> &node_topos);

        // SIMP 密度插值
        void simp_density(const vector<double> &den, vector<double> &simp_den, bool sen);

        // 敏度分析
        double sen_analysis(vector<double> &sen_den, vector<double> &filter_den, vector<double> &dc_topo, data_management &data_cae);

        // 变量更新
        void oc_update(vector<double> &den, vector<double> &filter_den, SpMatrix &H, vector<double> &Hs, 
            vector<double> &dc_topo, vector<double> &dv_topo);

        // 返回自由度和坐标
        void TopoLinear3D::build_ele_coors(Eigen::Ref<Eigen::MatrixXd> item_ele_coors, data_management &data_cae, int ele_id, int num_nodes);

    protected:
        // 计算两点之间的距离
        double distance2points(vector<double> &node1, vector<double> &node2);
    };
}