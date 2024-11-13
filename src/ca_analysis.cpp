/**************************************************************************

Copyright:  WH team

Author: YinJichao <jichaoyinyjc@163.com>

Completion date:  XXX

Description: XXX

**************************************************************************/

#include "include/ca_reanalysis.h"
namespace CAE
{

    // 计算delt K
    void ca_get_delt_stiffness(data_management &data_cae, assamble_stiffness &item_delt_k, elastic_mat &data_mat)
    {
        // 填充索引
        item_delt_k.row_idx_.assign(data_cae.item_assam_implicit_.row_idx_.begin(), data_cae.item_assam_implicit_.row_idx_.end());
        item_delt_k.col_idx_.assign(data_cae.item_assam_implicit_.col_idx_.begin(), data_cae.item_assam_implicit_.col_idx_.end());
        // 填充 非零元素
        item_delt_k.nz_val_.resize(data_cae.item_assam_implicit_.num_nz_val_);
        std::fill(item_delt_k.nz_val_.begin(), item_delt_k.nz_val_.end(), 0.);
        ca_fill_CSR_sparse_mat(data_cae, data_mat, item_delt_k);
    }

    // 基于CSR索引格式填充稀疏矩阵
    void ca_fill_CSR_sparse_mat(data_management &data_cae, elastic_mat &data_mat, assamble_stiffness &item_delt_k)
    {
        // 初始化单元
        data_cae.ele_inite(data_mat);
        // 初始化单元坐标，刚度矩阵
        vector<int> item_ele_dofs;
        MatrixXd item_ele_coors;
        MatrixXd stiffness_matrix;
        //  开始计算
        int n_modify = data_cae.node_topos_m_.size();
        for (int i = 0; i < n_modify; i++)
        {
            // 获取该单元的节点数量
            int ele_type = data_cae.ele_list_idx_m_[i];
            int map_idx = data_cae.ele_map_list_[ele_type];
            int node_num_ele = data_cae.ele_list_[map_idx]->nnode_;
            // 查找节点自由度及坐标
            item_ele_dofs.resize(3 * node_num_ele);
            item_ele_coors.resize(node_num_ele, 3);
            int item_dof, item_node;
            for (int j = 0; j < node_num_ele; j++)
            {
                item_node = data_cae.node_topos_m_[i][j + 1];
                // 自由度
                item_dof = data_cae.resort_free_nodes_o_[item_node];
                item_ele_dofs[3 * j] = 3 * item_dof;
                item_ele_dofs[3 * j + 1] = 3 * item_dof + 1;
                item_ele_dofs[3 * j + 2] = 3 * item_dof + 2;
                // if (i == 0)
                // {
                //     cout << item_dof << ":  " << item_ele_dofs[3 * j] << "  " << item_ele_dofs[3 * j + 1] << "  " << item_ele_dofs[3 * j + 2] << endl;
                // }
                // 坐标
                item_ele_coors(j, 0) = data_cae.coords_m_[item_node][0]; // X 坐标
                item_ele_coors(j, 1) = data_cae.coords_m_[item_node][1]; // Y 坐标
                item_ele_coors(j, 2) = data_cae.coords_m_[item_node][2]; // Z 坐标
                // if (i == 0)
                // {
                //     cout << item_ele_coors(j, 0) << "  " << item_ele_coors(j, 1) << "  " << item_ele_coors(j, 2) << endl;
                // }
            }
            // 计算 单元刚度矩阵变化量 K-K0
            stiffness_matrix.resize(3 * node_num_ele, 3 * node_num_ele);
            if (data_cae.node_topos_m_[i][0] == 0)
            {
                // 0：删除单元
                data_cae.ele_list_[map_idx]->build_ele_stiff_mat(item_ele_coors, stiffness_matrix);
                stiffness_matrix = -1. * stiffness_matrix;
            }
            else if (data_cae.node_topos_m_[i][0] == 1)
            {
                // 1：移动单元（即仅改变单元形状）
                MatrixXd item_ele_coors_o;
                item_ele_coors_o.resize(node_num_ele, 3); // 单元修改前的节点坐标
                MatrixXd stiffness_matrix_o;
                stiffness_matrix_o.resize(3 * node_num_ele, 3 * node_num_ele);
                // MatrixXd stiffness_matrix_m;
                // stiffness_matrix_m.resize(3 * node_num_ele, 3 * node_num_ele);
                for (int j = 0; j < node_num_ele; j++)
                {
                    item_node = data_cae.node_topos_m_[i][j + 1];
                    // if (i==0)
                    //     cout<<item_node<<":  ";
                    item_ele_coors_o(j, 0) = data_cae.coords_o_[item_node][0]; // X 坐标
                    item_ele_coors_o(j, 1) = data_cae.coords_o_[item_node][1]; // Y 坐标
                    item_ele_coors_o(j, 2) = data_cae.coords_o_[item_node][2]; // Z 坐标
                    // if (i == 0)
                    //     cout << item_ele_coors_o(j, 0) << "  " << item_ele_coors_o(j, 1) << "  " << item_ele_coors_o(j, 2) << endl;
                }
                data_cae.ele_list_[map_idx]->build_ele_stiff_mat(item_ele_coors, stiffness_matrix);
                // if (i == 0)
                //     cout << stiffness_matrix << endl;
                data_cae.ele_list_[map_idx]->build_ele_stiff_mat(item_ele_coors_o, stiffness_matrix_o);
                // if (i == 0)
                //     cout << endl << stiffness_matrix_o << endl;
                stiffness_matrix = stiffness_matrix - stiffness_matrix_o;
                // stiffness_matrix.setZero();
                // if (i == 0)
                //     cout << endl << stiffness_matrix << endl;
            }
            else if (data_cae.node_topos_m_[i][0] == 2)
            {
                // 2：增加单元
                // data_cae.ele_list_[data_cae.ele_list_idx_m_[i]]->build_ele_stiff_mat(item_ele_coors, stiffness_matrix);
            }
            else
            {
                cout << "Error in type of modified element.";
            }

            // 组装
            int ii_dof, jj_dof, loop_size = item_ele_dofs.size();
            for (int mm = 0; mm < loop_size; mm++)
            {
                jj_dof = item_ele_dofs[mm];
                if (jj_dof >= 0)
                {
                    int start = item_delt_k.col_idx_[jj_dof];
                    for (int nn = 0; nn < loop_size; nn++)
                    {
                        int t = start;
                        ii_dof = item_ele_dofs[nn];
                        if (ii_dof >= 0)
                        {
                            for (; item_delt_k.row_idx_[t] < ii_dof; t++)
                            {
                            } // 使得 t 对应的行索引 对应 ii_dof
                            item_delt_k.nz_val_[t] = item_delt_k.nz_val_[t] + stiffness_matrix(mm, nn);
                        }
                    }
                }
            }
        }
        cout << "The CSR index of delt_K has been filled" << endl;
    }

    // 建立组合近似降阶模型
    void ca_build_rom(data_management &data_cae, assamble_stiffness &item_delt_k, int n)
    {
        // 声明 ca 模型
        int row = int(data_cae.single_dis_vec_o_.size());
        vector<vector<double>> ca_rom(n, vector<double>(row, 0.));
        // 初始化 第一列
        ca_rom[0].assign(data_cae.single_dis_vec_o_.begin(), data_cae.single_dis_vec_o_.end());
        bool flag_clear = false;
        for (int i = 1; i < n; i++)
        {
            //
            vector<double> item_temp(row, 0.);
            Sparese_dot_vector(item_delt_k, ca_rom[i - 1], item_temp);
            //
            if (i == n - 1)
                bool flag_clear = true;
            bool aaa = data_cae.item_superlu.superlu_solution_next(item_temp, ca_rom[i], flag_clear);
            //
            for (int j = 0; j < row; j++)
            {
                ca_rom[i][j] = -1. * ca_rom[i][j];
            }
        }
        // SVD
        vector<vector<double>> ca_rom_SVD;
        ca_SVD(ca_rom, ca_rom_SVD);
        data_cae.ca_rom_n_ = ca_rom_SVD;
        cout << "SVD has been completed !!!\n";
    }

    // 基于 Eigen 求解 SVD
    void ca_SVD(vector<vector<double>> &ca_rom, vector<vector<double>> &ca_rom_SVD)
    {
        int row = int(ca_rom[0].size());
        int col = int(ca_rom.size());
        // 二维 Vector 赋值给 Eigen
        Eigen::MatrixXf rom_svd(row, col);
        for (int i = 0; i < row; i++)
        {
            for (int j = 0; j < col; j++)
            {
                rom_svd(i, j) = ca_rom[j][i];
            }
        }
        // SVD 分解
        // Eigen::BDCSVD<Eigen::MatrixXf> svd_holder;
        Eigen::JacobiSVD<Eigen::MatrixXf> svd_holder;
        svd_holder.compute(rom_svd, Eigen::ComputeThinU | Eigen::ComputeThinV);
        Eigen::MatrixXf svd_u = svd_holder.matrixU();
        ca_rom_SVD.resize(col, vector<double>(row, 0.));
        for (int i = 0; i < row; i++)
        {
            for (int j = 0; j < col; j++)
            {
                ca_rom_SVD[j][i] = svd_u(i, j);
            }
        }
    }

    // 求解降阶后的模型
    void ca_solve(data_management &data_cae, assamble_stiffness &item_k, vector<double> &solution)
    {
        // 计算缩减后的 系数矩阵
        vector<vector<double>> rk;
        ca_reduced_K(item_k, data_cae.ca_rom_n_, rk);
        // 计算缩减后的载荷向量
        vector<double> rf;
        ca_reduced_F(data_cae.single_load_vec_o_, data_cae.ca_rom_n_, rf);
        // 建立约简后的系数矩阵的 SCR
        int nn = int(data_cae.ca_rom_n_.size());
        vector<double> nz_val(nn * nn, 0.);
        vector<int> row_idx(nn * nn, 0.);
        vector<int> col_idx(nn + 1, 0.);
        int id_;
        for (int i = 0; i < nn; i++)
        {
            for (int j = 0; j < nn; j++)
            {
                id_ = i + j * nn; // 列压缩，“对称，行/列压缩无差别”
                nz_val[id_] = rk[i][j];
                row_idx[id_] = i;
            }
            col_idx[i + 1] = (i + 1) * nn;
        }
        // 求解 降阶方程
        vector<double> rx(nn, 0.);
        superlu_solver_func(nz_val, row_idx, col_idx, rf, rx);

        // 线性映射 降阶解 到 完全解
        int row = int(data_cae.ca_rom_n_[0].size());
        solution.resize(row);
        for (int i = 0; i < row; i++)
        {
            double sum = 0.;
            for (int j = 0; j < nn; j++)
            {
                sum += data_cae.ca_rom_n_[j][i] * rx[j];
            }
            solution[i] = sum;
            // cout<< i+1<<": "<<solution[i]<<endl;
        }
    }

    // 稀疏矩阵与列向量 乘法
    void Sparese_dot_vector(assamble_stiffness &item_delt_k, vector<double> &rom_vec_i_0, vector<double> &rom_vec_1)
    {
        int row = int(rom_vec_i_0.size());
        int vec_row;
        for (int i = 0; i < row; i++)
        {
            double sum = 0.;
            for (int t = item_delt_k.col_idx_[i]; t < item_delt_k.col_idx_[i + 1]; t++)
            {
                vec_row = item_delt_k.row_idx_[t];
                sum += item_delt_k.nz_val_[t] * rom_vec_i_0[vec_row];
            }
            rom_vec_1[i] = sum;
        }
    }

    // 计算约简后的系数矩阵
    void ca_reduced_K(assamble_stiffness &item_k, vector<vector<double>> &ca_rom_svd, vector<vector<double>> &rk)
    {
        int col = int(ca_rom_svd.size());
        int row = int(ca_rom_svd[0].size());
        rk.resize(col, vector<double>(col, 0.));
        vector<vector<double>> k_dot_rom(col, vector<double>(row, 0.));
        // k_dot_rom = item_k * ca_rom_n
        for (int i = 0; i < col; i++)
        {
            vector<double> temp(row, 0.);
            Sparese_dot_vector(item_k, ca_rom_svd[i], temp);
            k_dot_rom[i].assign(temp.begin(), temp.end());
        }
        // rk = ca_rom_n.T * k_dot_rom
        for (int i = 0; i < col; i++)
        {
            for (int j = 0; j < col; j++)
            {
                rk[i][j] = ca_vec_dot_vec(ca_rom_svd[i], k_dot_rom[j]);
                // std::cout << rk[i][j] << " ";
            }
            // std::cout << std::endl;
        }
    }

    // 向量与向量的乘法，返回标量
    double ca_vec_dot_vec(vector<double> &vec_1, vector<double> &vec_2)
    {
        double answer = 0.;
        int n1 = int(vec_1.size()), n2 = int(vec_2.size());
        if (n1 != n2)
        {
            cout << "the sizes betwenn vec_1 and vec_2 is uneuqal !!!\n";
        }
        else
        {
            for (int i = 0; i < n1; i++)
            {
                answer += vec_1[i] * vec_2[i];
            }
        }
        return answer;
    }

    // 计算约简后的 载荷向量
    void ca_reduced_F(vector<double> &f, vector<vector<double>> &ca_rom_n, vector<double> &rf)
    {
        int col = int(ca_rom_n.size());
        rf.resize(col, 0.);
        for (int i = 0; i < col; i++)
        {
            rf[i] = ca_vec_dot_vec(ca_rom_n[i], f);
        }
    }

    // 找二维向量最大值
    int find_max_element(vector<vector<int>> &vec)
    {
        int max_value = 0;
        for (int i = 0; i < vec.size(); i++)
        {
            for (int j = 0; j < vec[i].size(); j++)
            {
                if (vec[i][j] > max_value)
                {
                    max_value = vec[i][j];
                }
            }
        }
        return max_value;
    }
}