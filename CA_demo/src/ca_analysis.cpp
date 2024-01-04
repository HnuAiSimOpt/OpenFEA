/**************************************************************************

Copyright:  WH team

Author: YinJichao <jichaoyinyjc@163.com>

Completion date:  XXX

Description: XXX

**************************************************************************/

#include "include/ca_reanalysis.h"
namespace CAE
{

    // 重分析信息预处理
    void ca_pre_process(data_management &data_cae)
    {
        // 节点坐标
        vector<vector<int>> node_overlap;
        int kkk = 0;
        for (auto part_node : data_cae.parts_coods_)
        {
            if (kkk == 0)
            {
                data_cae.coords_.insert(data_cae.coords_.end(), part_node.begin(), part_node.end());
                kkk += 1;
            }
            else
            {
                int nn = data_cae.coords_.size();
                int kk = data_cae.coords_.size() + 1;
                for (int i = 0; i < part_node.size(); i++)
                {
                    part_node[i][0] += 140.;
                    part_node[i][1] += 8.;
                    if (part_node[i][1] == 0. || part_node[i][1] == 42.)
                    {
                        for (int j = 0; j < nn; j++)
                        {
                            if (data_cae.coords_[j][0] == part_node[i][0])
                            {
                                if (data_cae.coords_[j][1] == part_node[i][1])
                                {
                                    if (data_cae.coords_[j][2] == part_node[i][2])
                                    {
                                        vector<int> temp = {i, j + 1};
                                        node_overlap.push_back(temp);
                                    }
                                }
                            }
                        }
                    }
                    else
                    {
                        vector<double> temp_ = {part_node[i][0], part_node[i][1], part_node[i][2]};
                        data_cae.coords_.push_back(temp_);
                        vector<int> temp = {i, kk};
                        node_overlap.push_back(temp);
                        kk += 1;
                    }
                }
                // std::cout<<part_node.size()<<"  "<<node_overlap.size()<<std::endl;
            }
        }
        data_cae.nd_ = data_cae.coords_.size();
        // print
        std::cout << "Node: " << data_cae.nd_ << std::endl;
        // 拓扑关系
        auto part_topo = data_cae.parts_node_topos_[0][0];
        data_cae.node_topos_.insert(data_cae.node_topos_.end(), part_topo.begin(), part_topo.end());
        //
        auto part_topo_ = data_cae.parts_node_topos_[1][0];
        for (int i = 0; i < part_topo_.size(); i++)
        {
            vector<int> temp_ = {node_overlap[part_topo_[i][0] - 1][1], node_overlap[part_topo_[i][1] - 1][1], node_overlap[part_topo_[i][2] - 1][1], node_overlap[part_topo_[i][3] - 1][1]};
            data_cae.node_topos_.push_back(temp_);
        }
        data_cae.ne_ = data_cae.node_topos_.size();
        // print
        std::cout << "Element:" << data_cae.ne_ << std::endl;
        // 单元类型
        for (auto e_type : data_cae.parts_ele_type_)
            data_cae.ele_list_idx_.insert(data_cae.ele_list_idx_.end(), e_type.begin(), e_type.end());
    }

    // 计算delt K
    void ca_get_delt_stiffness(data_management &data_cae, assamble_stiffness &item_k, assamble_stiffness &item_delt_k, elastic_mat &data_mat)
    {
        // 赋值 模板刚度矩阵的 CSR 索引
        item_delt_k.nz_val.resize(item_k.num_nz_val);
        std::fill(item_delt_k.nz_val.begin(), item_delt_k.nz_val.end(), 0.);
        item_delt_k.row_idx.assign(item_k.row_idx.begin(), item_k.row_idx.end());
        item_delt_k.col_idx.assign(item_k.col_idx.begin(), item_k.col_idx.end());
        // 填充 非零元素
        ca_fill_CSR_sparse_mat(data_cae, data_mat, item_delt_k);
        // delt_K = K - K0, 所以取负号
        for (int i = 0; i < item_k.num_nz_val; i++)
        {
            item_delt_k.nz_val[i] = -1. * item_delt_k.nz_val[i];
        }
    }

    // 基于单元类型和节点拓扑关系，返回自由度，节点坐标
    void ca_build_ele_dofs_coors(vector<int> &item_ele_dofs, Eigen::Ref<Eigen::MatrixXd> item_ele_coors, data_management &data_cae, int ele_id, int num_nodes)
    {
        item_ele_dofs.resize(3 * num_nodes);
        std::fill(item_ele_dofs.begin(), item_ele_dofs.end(), 0);
        item_ele_coors.resize(num_nodes, 3);
        item_ele_coors.setZero();
        int item_dof, item_node;
        for (int i = 0; i < num_nodes; i++)
        {
            // 自由度
            item_dof = data_cae.resort_free_nodes_[data_cae.node_topos_[ele_id][i] - 1];
            item_ele_dofs[3 * i] = 3 * item_dof;
            item_ele_dofs[3 * i + 1] = 3 * item_dof + 1;
            item_ele_dofs[3 * i + 2] = 3 * item_dof + 2;
            // 坐标
            item_node = data_cae.node_topos_[ele_id][i] - 1;
            item_ele_coors(i, 0) = data_cae.coords_[item_node][0]; // X 坐标
            item_ele_coors(i, 1) = data_cae.coords_[item_node][1]; // Y 坐标
            item_ele_coors(i, 2) = data_cae.coords_[item_node][2]; // Z 坐标
        }
    }

    // 基于CSR索引格式填充稀疏矩阵
    void ca_fill_CSR_sparse_mat(data_management &data_cae, elastic_mat &data_mat, assamble_stiffness &item_delt_k)
    {
        // 初始化单元
        data_cae.ele_inite(data_mat);
        // 初始化单元坐标，刚度矩阵
        int num_nodes;
        vector<int> item_ele_dofs;
        MatrixXd item_ele_coors;
        MatrixXd stiffness_matrix;
        for (int id_ele = data_cae.part_ne_[0]; id_ele < data_cae.ne_; id_ele++)
        {
            // 获取该单元的节点数量
            int node_num_ele = data_cae.ele_list_[data_cae.ele_list_idx_[id_ele] - 1]->nnode_;
            // 查找节点自由度及坐标
            item_ele_coors.resize(node_num_ele, 3);
            ca_build_ele_dofs_coors(item_ele_dofs, item_ele_coors, data_cae, id_ele, node_num_ele);
            // 计算单元刚度矩阵
            stiffness_matrix.resize(3 * node_num_ele, 3 * node_num_ele);
            data_cae.ele_list_[data_cae.ele_list_idx_[id_ele] - 1]->build_ele_stiff_mat(item_ele_coors, stiffness_matrix);
            // 组装
            int ii_dof, jj_dof, loop_size = item_ele_dofs.size();
            for (int mm = 0; mm < loop_size; mm++)
            {
                jj_dof = item_ele_dofs[mm];
                if (jj_dof >= 0)
                {
                    int start = item_delt_k.col_idx[jj_dof];
                    for (int nn = 0; nn < loop_size; nn++)
                    {
                        int t = start;
                        ii_dof = item_ele_dofs[nn];
                        if (ii_dof >= 0)
                        {
                            for (; item_delt_k.row_idx[t] < ii_dof; t++)
                            {
                            } // 使得 t 对应的行索引 对应 ii_dof
                            item_delt_k.nz_val[t] = item_delt_k.nz_val[t] + stiffness_matrix(mm, nn);
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
        int row = int(data_cae.single_dis_vec_.size());
        vector<vector<double>> ca_rom(n, vector<double>(row, 0.));
        ca_rom[0].assign(data_cae.single_dis_vec_.begin(), data_cae.single_dis_vec_.end());
        for (int i = 1; i < n; i++)
        {
            //
            vector<double> item_temp(row, 0.);
            Sparese_dot_vector(item_delt_k, ca_rom[i - 1], item_temp);
            //
            for (int j = 0; j < row; j++)
            {
                item_temp[j] = -1. * item_temp[j];
            }
            //
            int phase_33 = data_cae.item_pardiso.pardiso_solution(item_temp, ca_rom[i]);
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
    void ca_solve(data_management &data_cae, assamble_stiffness &item_k)
    {
        // 计算缩减后的 系数矩阵
        vector<vector<double>> rk;
        ca_reduced_K(item_k, data_cae.ca_rom_n_, rk);
        // 计算缩减后的载荷向量
        vector<double> rf;
        ca_reduced_F(data_cae.single_load_vec_, data_cae.ca_rom_n_, rf);
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
        // 求解
        vector<double> rx(nn, 0.);
        superlu_solver(nz_val, row_idx, col_idx, rf, rx);
        //
        int row = int(data_cae.ca_rom_n_[0].size());
        vector<double> solution(row, 0.);
        for (int i = 0; i < row; i++)
        {
            double sum = 0.;
            for (int j = 0; j < nn; j++)
            {
                sum += data_cae.ca_rom_n_[j][i] * rx[j];
            }
            solution[i] = sum;
        }
        // 写入 TXT
        fstream f;
        f.open("data.txt", ios::out);
        for (int i = 0; i < data_cae.part_nd_[0]; i++)
        {
            string s1 = to_string(solution[i]);
            f << s1 <<"\n";
        }
        f.close();
        cout<<"ending !!!\n";
    }

    // 稀疏矩阵与列向量 乘法
    void Sparese_dot_vector(assamble_stiffness &item_delt_k, vector<double> &rom_vec_i_0, vector<double> &rom_vec_1)
    {
        int row = int(rom_vec_i_0.size());
        int vec_row;
        for (int i = 0; i < row; i++)
        {
            double sum = 0.;
            for (int t = item_delt_k.col_idx[i]; t < item_delt_k.col_idx[i + 1]; t++)
            {
                vec_row = item_delt_k.row_idx[t];
                sum += item_delt_k.nz_val[t] * rom_vec_i_0[vec_row];
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
            }
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