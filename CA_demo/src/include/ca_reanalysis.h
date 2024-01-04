/**************************************************************************

Copyright:  WH team

Author: YinJichao <jichaoyinyjc@163.com>

Completion date:  XXX

Description: XXX

**************************************************************************/

#include <algorithm>
#include "include/assemble.h"
#include "include/data_management.h"
#include <Eigen/SVD>
#include <Eigen/Core>
#include "solver/include/solver_superlu.h"
#include <fstream>
#include <string>

using std::to_string;

namespace CAE
{
    // 重分析信息预处理
    void ca_pre_process(data_management &data_cae);

    // 基于单元类型和节点拓扑关系，返回自由度，节点坐标
    void ca_build_ele_dofs_coors(vector<int> &item_ele_dofs, Eigen::Ref<Eigen::MatrixXd> item_ele_coors, data_management &data_cae, int ele_id, int num_nodes);

    // 计算delt K
    void ca_get_delt_stiffness(data_management &data_cae, assamble_stiffness &item_k, assamble_stiffness &item_delt_k, elastic_mat &data_mat);

    // 基于CSR索引格式填充稀疏矩阵
    void ca_fill_CSR_sparse_mat(data_management &data_cae, elastic_mat &data_mat, assamble_stiffness &item_delt_k);

    // 建立组合近似降阶模型
    void ca_build_rom(data_management &data_cae, assamble_stiffness &item_delt_k, int n);

    // 基于 Eigen 求解 SVD
    void ca_SVD(vector<vector<double>> &ca_rom_n, vector<vector<double>> &ca_rom_SVD);

    // 求解降阶后的模型
    void ca_solve(data_management &data_cae, assamble_stiffness &item_k);

    // 稀疏矩阵与列向量 乘法
    void Sparese_dot_vector(assamble_stiffness &item_delt_k, vector<double> &rom_vec_i_0, vector<double> &rom_vec_1);

    // 计算约简后的系数矩阵
    void ca_reduced_K(assamble_stiffness &item_k, vector<vector<double>> &ca_rom_n, vector<vector<double>> &rk);

    // 向量与向量的乘法，返回标量
    double ca_vec_dot_vec(vector<double> &vec_1, vector<double> &vec_2);

    // 计算约简后的 载荷向量
    void ca_reduced_F(vector<double> &f, vector<vector<double>> &ca_rom_n, vector<double> &rf);

    // 找二维向量最大值
    int find_max_element(vector<vector<int>> &vec);
}