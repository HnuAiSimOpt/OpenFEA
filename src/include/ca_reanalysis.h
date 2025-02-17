/**************************************************************************

Copyright:  WH team

Author: YinJichao <jichaoyinyjc@163.com>

Completion date:  XXX

Description: XXX

**************************************************************************/

#include <cmath>
#include <algorithm>
#include "include/assemble.h"
#include "include/data_management.h"
#include <Eigen/SVD>
#include <Eigen/Core>
#include "solver/include/solver_superlu.h"
#include <fstream>
#include <string>
#include "./utils.h"

using std::to_string;

namespace CAE
{
    class CA_solution
    {
    private:
        int n_basis_ = 0;
        int num_fix_nodes_ = 0;
        int num_all_nodes_ = 0;  // 对应nodes_all 的所有节点数
        int num_free_nodes_ = 0; // 所有点减去被固定点，包含被删除的点，用以计算修改后的刚度矩阵自由度数
        int num_change_eles_ = 0;
        vector<vector<double>> ROM_;

    public:
        vector<double> ca_dis_vec_;      // 对应 nodes_all 的全位移 无约束位移
        vector<double> ca_full_dis_vec_; // 对应 nodes_all 的全位移
        vector<double> extract_dis_vec_; // 对应 new_model 的全位移

    public:
        // 重分析过程
        void ca_analysis(data_management &data_cae, elastic_mat &data_mat, int n_basis, bool is_Update);
        // 重分析重置索引
        void CA_re_CSR(data_management &data_cae, assamble_stiffness &current_K);
        // 重分析填充参考刚度矩阵
        void CA_copy_ref_stiff(assamble_stiffness &current_K, assamble_stiffness_save &stiff_ori);
        // 计算刚度矩阵变化量
        void CA_get_deltK(data_management &data_cae, assamble_stiffness &delt_K, assamble_stiffness &ref_K, elastic_mat &data_mat);
        // 获取自由度和坐标
        void CA_get_dof_coor(data_management &data_cae, vector<int> &item_ele_dofs, MatrixXd &item_ele_coors, int eid, int node_num_ele);
        // 获取坐标【for 单元变形】
        void CA_get_coor_for_mdf(data_management &data_cae, MatrixXd &item_ele_coors_new, int eid, int node_num_ele);
        // 构造组合近似降阶模型
        void ca_build_rom(data_management &data_cae, assamble_stiffness &delt_K, int n_basis);
        // 稀疏矩阵与列向量 乘法
        void Sparese_dot_vector(assamble_stiffness &sp_mat, vector<double> &vec_in, vector<double> &vec_out);
        // Schmidt 正交
        void ca_schmidt(vector<vector<double>> &ca_rom_n, vector<vector<double>> &ca_rom_schmidt);
        // 返回 norm
        double ca_vec_norm(const vector<double> vec);
        // SVD分解
        void ca_SVD(vector<vector<double>> &ca_rom_n, vector<vector<double>> &ca_rom_SVD);
        // 求解降阶后的模型
        void ca_solve(data_management &data_cae, assamble_stiffness &item_k, vector<double> &solution);
        // 计算约简后的系数矩阵
        void ca_reduced_K(assamble_stiffness &item_k, vector<vector<double>> &rk);
        // 计算约简后的 载荷向量
        void ca_reduced_F(vector<double> &f, vector<double> &rf);
        // 向量与向量的乘法，返回标量
        double ca_vec_dot_vec(vector<double> &vec_1, vector<double> &vec_2);
        // 计算位移进行处理
        void dis_process(data_management &data_cae);

        // 检查变形计算后的刚度矩阵
        void check_stiffness(data_management &data_cae, assamble_stiffness &current_K, assamble_stiffness &old_K);
        void check_stiffness(data_management &data_cae, assamble_stiffness &current_K, assamble_stiffness_save &old_K);
    };
}