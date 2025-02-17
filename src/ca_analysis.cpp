/**************************************************************************

Copyright:  WH team

Author: YinJichao <jichaoyinyjc@163.com>

Completion date:  XXX

Description: XXX

**************************************************************************/

#include "include/ca_reanalysis.h"
namespace CAE
{
    void CA_solution::ca_analysis(data_management &data_cae, elastic_mat &data_mat, int n_basis, bool is_Update)
    {
        clock_t start, end; // 定义clock_t变量
        this->n_basis_ = n_basis;
        this->num_fix_nodes_ = data_cae.dis_bc_set_.size();
        this->num_all_nodes_ = data_cae.mesh_ca_map_.union_coords_.size();
        this->num_free_nodes_ = num_all_nodes_ - this->num_fix_nodes_;
        this->num_change_eles_ = data_cae.mesh_ca_map_.change_node_topos_.size();

        // 建立索引, 填充稀疏刚度矩阵
        assamble_stiffness delt_K, current_K, ref_K;
        current_K.num_row_ = 3 * this->num_free_nodes_;
        current_K.num_col_ = 3 * this->num_free_nodes_;
        this->CA_re_CSR(data_cae, current_K);
        this->CA_copy_ref_stiff(current_K, data_cae.save_ref_info_.stiff_ori_);
        // 为了不修改原始的参考刚度矩阵，复制拷贝一份
        ref_K.num_col_ = current_K.num_col_;
        ref_K.num_row_ = current_K.num_row_;
        ref_K.num_nz_val_ = current_K.num_nz_val_;
        ref_K.col_idx_.assign(current_K.col_idx_.begin(), current_K.col_idx_.end());
        ref_K.row_idx_.assign(current_K.row_idx_.begin(), current_K.row_idx_.end());
        ref_K.nz_val_.assign(current_K.nz_val_.begin(), current_K.nz_val_.end());
        // write_spmat(current_K.row_idx_,
        //             current_K.col_idx_,
        //             current_K.nz_val_,
        //             "C:\\Users\\jicha\\Desktop\\bracketcubediff\\C\\current_K_row_value.txt",
        //             "C:\\Users\\jicha\\Desktop\\bracketcubediff\\C\\current_K_col.txt");
        // write_spmat(ref_K.row_idx_,
        //             ref_K.col_idx_,
        //             ref_K.nz_val_,
        //             "C:\\Users\\jicha\\Desktop\\bracketcubediff\\C\\ref_K_row_value.txt",
        //             "C:\\Users\\jicha\\Desktop\\bracketcubediff\\C\\ref_K_col.txt");
        // write_spmat(data_cae.save_ref_info_.stiff_ori_.row_idx_,
        //             data_cae.save_ref_info_.stiff_ori_.col_idx_,
        //             data_cae.save_ref_info_.stiff_ori_.nz_val_,
        //             "C:\\Users\\jicha\\Desktop\\bracketcubediff\\C\\old_K_row_value.txt",
        //             "C:\\Users\\jicha\\Desktop\\bracketcubediff\\C\\old_K_col.txt");
        // 计算刚度矩阵变化
        start = clock();
        this->CA_get_deltK(data_cae, delt_K, ref_K, data_mat);
        // 如果新增单元，参考刚度矩阵被修改，则需要重新分解参考刚度矩阵
        if (ref_K.num_col_ != data_cae.save_ref_info_.stiff_ori_.num_col_)
        {
            data_cae.item_pardiso.clear_data();
            int n = ref_K.num_col_;
            data_cae.item_pardiso.phase_00 = data_cae.item_pardiso_ca.pardiso_init(ref_K.nz_val_, ref_K.row_idx_, ref_K.col_idx_, n);
            data_cae.item_pardiso.phase_1122 = data_cae.item_pardiso_ca.pardiso_decomposition();
        }
        // 根据 delt_K 计算 current_K
        for (int i = 0; i < current_K.num_nz_val_; i++)
        {
            current_K.nz_val_[i] = current_K.nz_val_[i] + delt_K.nz_val_[i];
        }
        end = clock();
        cout << "It took " << double(end - start) / CLOCKS_PER_SEC << " s to compute the amount of change in the stiffness matrix" << endl;
        // write_spmat(current_K.row_idx_,
        //             current_K.col_idx_,
        //             current_K.nz_val_,
        //             "C:\\Users\\jicha\\Desktop\\bracketcubediff\\C\\current_K_row_value.txt",
        //             "C:\\Users\\jicha\\Desktop\\bracketcubediff\\C\\current_K_col.txt");
        // check_stiffness(data_cae, current_K, data_cae.save_ref_info_.stiff_ori_); // 由于删除单元，故current_K是奇异的，可将减少单元-1.0设置为-0.99999999以检查current_K

        // 构造组合近似降阶模型计时
        start = clock();
        ca_build_rom(data_cae, delt_K, n_basis); // 计算组合近似降阶模型
        end = clock();
        cout << "It took " << double(end - start) / CLOCKS_PER_SEC << " s to compute the CA model" << endl;

        // 求解
        ca_dis_vec_.resize(current_K.num_col_);
        std::fill(ca_dis_vec_.begin(), ca_dis_vec_.end(), 0.);
        this->ca_solve(data_cae, current_K, ca_dis_vec_);

        //
        dis_process(data_cae);
        // write_vec(ca_full_dis_vec_, "ca_full_dis_vec_.txt");
        // write_vec(extract_dis_vec_, "extract_dis_vec_.txt");
    };

    // 重构CSR索引
    void CA_solution::CA_re_CSR(data_management &data_cae, assamble_stiffness &current_K)
    {
        // 重新申明内存大小,并存入原始数据
        current_K.col_data_.resize(current_K.num_col_);
        for (int i = 0; i < data_cae.save_ref_info_.stiff_ori_.num_col_; i++)
        {
            current_K.col_data_[i] = data_cae.save_ref_info_.stiff_ori_.col_data_[i];
        }
        // 新增单元更新 col_data_ 内容
        for (int id_ele = 0; id_ele < this->num_change_eles_; id_ele++)
        {
            vector<int> item_ele_dofs;
            if (data_cae.mesh_ca_map_.change_node_topos_[id_ele][0] == 2)
            {
                // 识别单元类型
                int ele_type = data_cae.mesh_ca_map_.change_ele_idx_[id_ele];
                int map_idx = data_cae.ele_map_list_[ele_type];
                int num_nodes = data_cae.ele_list_[map_idx]->nnode_;
                //
                item_ele_dofs.resize(3 * num_nodes);
                std::fill(item_ele_dofs.begin(), item_ele_dofs.end(), 0);
                int item_coor_idx_use, item_coor_idx;
                for (int i = 0; i < num_nodes; i++)
                {
                    item_coor_idx = data_cae.mesh_ca_map_.change_node_topos_[id_ele][i + 1]; // data_cae_.node_topos_diff_map_中节点索引从 0 开始
                    if (item_coor_idx + 1 > data_cae.nd_)
                    {
                        item_coor_idx_use = item_coor_idx - this->num_fix_nodes_; // 新增节点 直接减去 无约束节点数 就是 当前自由节点索引
                    }
                    else
                    {
                        item_coor_idx_use = data_cae.save_ref_info_.resort_free_nodes_ori_[item_coor_idx]; // 非新增节点 在 重排节点顺序 索引
                    }
                    item_ele_dofs[3 * i] = 3 * item_coor_idx_use;
                    item_ele_dofs[3 * i + 1] = 3 * item_coor_idx_use + 1;
                    item_ele_dofs[3 * i + 2] = 3 * item_coor_idx_use + 2;
                }
                // 删除负自由度，即被约束自由度
                delete_negative(item_ele_dofs);
                // 压缩稀疏矩阵
                for (int id_dofs_row : item_ele_dofs)
                {
                    for (int id_dofs_col : item_ele_dofs)
                    {
                        current_K.col_data_[id_dofs_col].insert(id_dofs_row);
                    }
                }
            }
        }
        // 建立新增单元后的CSR索引格式
        int num_nz_val_ = 0;
        for (int id_dof = 0; id_dof < current_K.num_col_; id_dof++)
        {
            num_nz_val_ += current_K.col_data_[id_dof].size(); // 计算每个单元非零元数目
        }
        current_K.num_nz_val_ = num_nz_val_;
        // 分配行索引容量
        current_K.row_idx_.resize(num_nz_val_);
        // 分配列索引容量
        current_K.col_idx_.resize(current_K.num_col_ + 1);
        // 建立列索引
        int item_idx_csr = 0;
        current_K.col_idx_[0] = 0;
        for (int i = 0; i < current_K.num_col_; i++)
        {
            for (int row : current_K.col_data_[i])
            {
                current_K.row_idx_[item_idx_csr] = row;
                item_idx_csr++;
            }
            current_K.col_idx_[i + 1] = item_idx_csr;
        }
        cout << "The CSR index of increase CA mode has been built" << endl;
    }

    // 重分析填充参考刚度矩阵
    void CA_solution::CA_copy_ref_stiff(assamble_stiffness &current_K, assamble_stiffness_save &stiff_ori)
    {
        if (current_K.num_col_ == stiff_ori.num_col_) // 即没有新增单元， 直接拷贝非0元素
        {
            current_K.nz_val_.resize(stiff_ori.num_nz_val_);
            current_K.nz_val_.assign(stiff_ori.nz_val_.begin(), stiff_ori.nz_val_.end()); // 填充参考模型的值
        }
        else
        {
            current_K.nz_val_.resize(current_K.num_nz_val_);
            std::fill(current_K.nz_val_.begin(), current_K.nz_val_.end(), 0.); // 初始化并置零
            // 填充参考模型的值
            int stiff_ori_sta, stiff_ori_end, current_ptr;
            for (int col = 0; col < stiff_ori.num_col_; col++) // 参考模型自由度不超过num_col_，即非0元素在 [0, ref_K.num_col_]^2 区域
            {
                stiff_ori_sta = stiff_ori.col_idx_[col];                            // 参考矩阵当前列 第一个非0元素 行索引
                stiff_ori_end = stiff_ori.col_idx_[col + 1] - 1;                    // 参考矩阵当前列 最后一个非0元素 行索引
                current_ptr = current_K.col_idx_[col];                              // 参考矩阵每行元素超过该行索引最大值
                for (int nz_idx = stiff_ori_sta; nz_idx <= stiff_ori_end; nz_idx++) // 填充 值
                {
                    for (; current_K.row_idx_[current_ptr] < stiff_ori.row_idx_[nz_idx]; current_ptr++) // 保证插值到相同行
                    {
                    }
                    current_K.nz_val_[current_ptr] = stiff_ori.nz_val_[nz_idx];
                    current_ptr += 1;
                }
            }
        }
        cout << "The CSR index of ref stiffness matrix has been copied" << endl;
    };

    // 计算刚度矩阵变化量
    void CA_solution::CA_get_deltK(data_management &data_cae, assamble_stiffness &delt_K, assamble_stiffness &ref_K, elastic_mat &data_mat)
    {
        // 初始化
        delt_K.num_row_ = ref_K.num_row_;
        delt_K.num_col_ = ref_K.num_col_;
        delt_K.num_nz_val_ = ref_K.num_nz_val_;
        delt_K.col_idx_.assign(ref_K.col_idx_.begin(), ref_K.col_idx_.end());
        delt_K.row_idx_.assign(ref_K.row_idx_.begin(), ref_K.row_idx_.end());
        delt_K.nz_val_.resize(delt_K.num_nz_val_);
        std::fill(delt_K.nz_val_.begin(), delt_K.nz_val_.end(), 0.);
        // 开始计算形变量
        data_cae.ele_inite(data_mat);
        vector<int> item_ele_dofs;
        MatrixXd item_ele_coors, item_ele_coors_new;
        MatrixXd stiff_matrix_new, stiff_matrix, stiff_matrix_delt;
        int change_type;
        for (int i = 0; i < num_change_eles_; i++)
        {
            // 获取该单元的节点数量
            int ele_type = data_cae.mesh_ca_map_.change_ele_idx_[i];
            int map_idx = data_cae.ele_map_list_[ele_type];
            int node_num_ele = data_cae.ele_list_[map_idx]->nnode_;

            // 查找节点自由度及坐标
            item_ele_dofs.resize(3 * node_num_ele);
            item_ele_coors.resize(node_num_ele, 3);
            this->CA_get_dof_coor(data_cae, item_ele_dofs, item_ele_coors, i, node_num_ele);
            // 计算 单元刚度矩阵变化量 K-K0
            stiff_matrix_delt.resize(3 * node_num_ele, 3 * node_num_ele);
            stiff_matrix.resize(3 * node_num_ele, 3 * node_num_ele);
            data_cae.ele_list_[map_idx]->build_ele_stiff_mat(item_ele_coors, stiff_matrix);
            change_type = data_cae.mesh_ca_map_.change_node_topos_[i][0];
            if (change_type == 0) // 0：删除单元
            {
                stiff_matrix_delt = -1.0 * stiff_matrix; // 由于删除单元，故current_K是奇异的，可将1.0设置为0.99999999以检查current_K
            }
            else if (change_type == 1) // 1：变形单元（即仅改变单元形状）
            {
                item_ele_coors_new.resize(node_num_ele, 3);
                this->CA_get_coor_for_mdf(data_cae, item_ele_coors_new, i, node_num_ele);
                stiff_matrix_new.resize(3 * node_num_ele, 3 * node_num_ele);
                data_cae.ele_list_[map_idx]->build_ele_stiff_mat(item_ele_coors_new, stiff_matrix_new);
                stiff_matrix_delt = stiff_matrix_new - stiff_matrix;
            }
            else if (change_type == 2) // 2：增加单元
            {
                stiff_matrix_delt = 1.0 * stiff_matrix;
                // cout<<stiff_matrix_delt;
                // exit(0);
            }
            else
            {
                cout << "Error in type of modified element.";
            }
            // 组装
            int tt, ii_dof, jj_dof, loop_size = item_ele_dofs.size();
            for (int mm = 0; mm < loop_size; mm++)
            {
                jj_dof = item_ele_dofs[mm];
                int start = delt_K.col_idx_[jj_dof]; // 列指针
                for (int nn = 0; nn < loop_size; nn++)
                {
                    tt = start;
                    ii_dof = item_ele_dofs[nn]; // 行索引
                    for (; delt_K.row_idx_[tt] < ii_dof; tt++)
                    {
                    } // 使得 t 对应的行索引 对应 ii_dof
                    delt_K.nz_val_[tt] += stiff_matrix_delt(mm, nn);
                    if (change_type == 2)
                    {
                        ref_K.nz_val_[tt] += 1.0E-3 * stiff_matrix(mm, nn);
                    }
                }
            }
        }
        cout << "The deltK has been built" << endl;
    };

    // 获取自由度和坐标
    void CA_solution::CA_get_dof_coor(data_management &data_cae, vector<int> &item_ele_dofs, MatrixXd &item_ele_coors, int eid, int node_num_ele)
    {
        int item_dof, item_node_use, item_node;
        for (int j = 0; j < node_num_ele; j++)
        {
            item_node = data_cae.mesh_ca_map_.change_node_topos_[eid][j + 1]; // data_cae_.node_topos_diff_map_中节点索引从 0 开始
            if (item_node + 1 > data_cae.nd_)
            {
                item_node_use = item_node - num_fix_nodes_;
            }
            else
            {
                item_node_use = data_cae.save_ref_info_.resort_free_nodes_ori_[item_node];
            }
            item_ele_dofs[3 * j] = 3 * item_node_use;
            item_ele_dofs[3 * j + 1] = 3 * item_node_use + 1;
            item_ele_dofs[3 * j + 2] = 3 * item_node_use + 2;
            item_ele_coors(j, 0) = data_cae.mesh_ca_map_.union_coords_[item_node][0]; // X 坐标
            item_ele_coors(j, 1) = data_cae.mesh_ca_map_.union_coords_[item_node][1]; // Y 坐标
            item_ele_coors(j, 2) = data_cae.mesh_ca_map_.union_coords_[item_node][2]; // Z 坐标
        }
    };
    // 获取坐标【for 单元变形】
    void CA_solution::CA_get_coor_for_mdf(data_management &data_cae, MatrixXd &item_ele_coors_new, int eid, int node_num_ele)
    {
        int item_node;
        for (int j = 0; j < node_num_ele; j++)
        {
            int idx_ = data_cae.mesh_ca_map_.change_node_topos_[eid][j + 1];
            item_node = data_cae.mesh_ca_map_.node_idx_union_map_[idx_];
            item_ele_coors_new(j, 0) = data_cae.mesh_ca_new_.mdf_coords_[item_node][0]; // X 坐标
            item_ele_coors_new(j, 1) = data_cae.mesh_ca_new_.mdf_coords_[item_node][1]; // Y 坐标
            item_ele_coors_new(j, 2) = data_cae.mesh_ca_new_.mdf_coords_[item_node][2]; // Z 坐标
        }
    };

    // 构造组合近似降阶模型
    void CA_solution::ca_build_rom(data_management &data_cae, assamble_stiffness &delt_K, int n_basis)
    {
        // 声明 ca 模型
        int row = delt_K.num_row_;
        vector<vector<double>> init_ROM(n_basis, vector<double>(row, 0.));
        // 初始化 第一列

        for (int i = 0; i < data_cae.save_ref_info_.dis_vec_ori_.size(); i++)
        {
            init_ROM[0][i] = data_cae.save_ref_info_.dis_vec_ori_[i];
        }
        bool flag_clear = false;
        vector<double> item_temp(row, 0.);
        vector<double> item_temp_x(row, 0.);
        for (int i = 1; i < n_basis; i++)
        {
            std::fill(item_temp.begin(), item_temp.end(), 0.);
            Sparese_dot_vector(delt_K, init_ROM[i - 1], item_temp);
            // if (i == 1)
            // {
            //     write_vec(item_temp, "vec.txt");
            //     exit(0);
            // }
            if (i == n_basis - 1)
            {
                bool flag_clear = true;
            }
            std::fill(item_temp_x.begin(), item_temp_x.end(), 0.);
            // bool splu = data_cae.item_superlu.superlu_solution_next(item_temp, ca_rom[i], flag_clear);
            bool pds = data_cae.item_pardiso_ca.pardiso_solution(item_temp, item_temp_x);
            for (int j = 0; j < row; j++)
            {
                init_ROM[i][j] = -1.0 * item_temp_x[j];
            }
        }
        // SVD
        // vector<vector<double>> ca_rom_SVD;
        // this->ca_SVD(init_ROM, ca_rom_SVD);
        // this->ROM_ = ca_rom_SVD
        // cout << "SVD has been completed !!!\n";
        vector<vector<double>> ca_rom_schmidt;
        this->ca_schmidt(init_ROM, ca_rom_schmidt);
        this->ROM_ = ca_rom_schmidt;
        cout << "The Gram-Schmidt Orthogonalization has been completed !!!\n";
    };

    // 稀疏矩阵与列向量 乘法
    void CA_solution::Sparese_dot_vector(assamble_stiffness &sp_mat, vector<double> &vec_in, vector<double> &vec_out)
    {
        int row = int(vec_in.size());
        int vec_row;
        for (int i = 0; i < row; i++)
        {
            double sum = 0.;
            for (int t = sp_mat.col_idx_[i]; t < sp_mat.col_idx_[i + 1]; t++)
            {
                vec_row = sp_mat.row_idx_[t];
                sum += sp_mat.nz_val_[t] * vec_in[vec_row];
            }
            vec_out[i] = sum;
        }
    };

    // Schmidt 正交
    void CA_solution::ca_schmidt(vector<vector<double>> &ca_rom_n, vector<vector<double>> &ca_rom_schmidt)
    {
        int col = ca_rom_n.size();
        int row = ca_rom_n[0].size();
        ca_rom_schmidt.resize(col, vector<double>(row, 0.));
        double dot_1, dot_2, fac, norm_vec;
        // 
        ca_rom_schmidt[0] = ca_rom_n[0];
        norm_vec = this->ca_vec_norm(ca_rom_schmidt[0]);
        for (int k = 0; k < row; k++)
        {
            ca_rom_schmidt[ 0][k] = ca_rom_schmidt[0][k] / norm_vec;
        }
        for (int i = 1; i < col; i++)
        {
            ca_rom_schmidt[i] = ca_rom_n[i];
            for (int j = 0; j < i - 1; j++)
            {
                dot_1 = this->ca_vec_dot_vec(ca_rom_schmidt[j], ca_rom_n[i]);
                dot_2 = this->ca_vec_dot_vec(ca_rom_schmidt[j], ca_rom_schmidt[j]);
                fac = dot_1 / dot_2;
                for (int k = 0; k < row; k++)
                {
                    ca_rom_schmidt[i][k] = ca_rom_schmidt[i][k] - fac * ca_rom_schmidt[j][k];
                }
            }
            norm_vec = this->ca_vec_norm(ca_rom_schmidt[i]);
            if (norm_vec < 1.0E-10)
            {
                cout << "The vectors are linearly dependent.";
            }
            else
            {
                for (int k = 0; k < row; k++)
                {
                    ca_rom_schmidt[i][k] = ca_rom_schmidt[i][k] / norm_vec;
                }
            }
        }
    };

    // 返回 norm
    double CA_solution::ca_vec_norm(const vector<double> vec)
    {
        int nn = vec.size();
        double vec_sum = 0.;
        for (int i = 0; i < nn; i++)
        {
            vec_sum += vec[i] * vec[i];
        }
        double vec_norm = sqrt(vec_sum);
        return vec_norm;
    };

    // 基于 Eigen 求解 SVD
    void CA_solution::ca_SVD(vector<vector<double>> &ca_rom, vector<vector<double>> &ca_rom_SVD)
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
    };

    // 检查变形计算后的刚度矩阵
    void CA_solution::check_stiffness(data_management &data_cae, assamble_stiffness &current_K, assamble_stiffness &old_K)
    {
        vector<double> current_single_dis_vec(current_K.num_col_, 0.);
        vector<double> single_load_vec_1(current_K.num_col_, 0.);
        vector<double> old_single_dis_vec(old_K.num_col_, 0.);
        vector<double> single_load_vec_2(old_K.num_col_, 0.);
        for (int i = 0; i < data_cae.save_ref_info_.load_vec_ori_.size(); i++)
        {
            single_load_vec_1[i] = data_cae.save_ref_info_.load_vec_ori_[i];
            single_load_vec_2[i] = data_cae.save_ref_info_.load_vec_ori_[i];
        }
        superlu_solver_func(current_K.nz_val_, current_K.row_idx_, current_K.col_idx_, single_load_vec_1, current_single_dis_vec);
        superlu_solver_func(old_K.nz_val_, old_K.row_idx_, old_K.col_idx_, single_load_vec_2, old_single_dis_vec);
        write_vec(current_single_dis_vec, "dis_current.txt");
        write_vec(old_single_dis_vec, "dis_old.txt");
    };

    //
    void CA_solution::check_stiffness(data_management &data_cae, assamble_stiffness &current_K, assamble_stiffness_save &old_K)
    {
        vector<double> current_single_dis_vec(current_K.num_col_, 0.);
        vector<double> single_load_vec_1(current_K.num_col_, 0.);
        vector<double> old_single_dis_vec(old_K.num_col_, 0.);
        vector<double> single_load_vec_2(old_K.num_col_, 0.);
        for (int i = 0; i < data_cae.save_ref_info_.load_vec_ori_.size(); i++)
        {
            single_load_vec_1[i] = data_cae.save_ref_info_.load_vec_ori_[i];
            single_load_vec_2[i] = data_cae.save_ref_info_.load_vec_ori_[i];
        }
        superlu_solver_func(current_K.nz_val_, current_K.row_idx_, current_K.col_idx_, single_load_vec_1, current_single_dis_vec);
        superlu_solver_func(old_K.nz_val_, old_K.row_idx_, old_K.col_idx_, single_load_vec_2, old_single_dis_vec);
        write_vec(current_single_dis_vec, "dis_current.txt");
        write_vec(old_single_dis_vec, "dis_old.txt");
    };

    // 求解降阶后的模型
    void CA_solution::ca_solve(data_management &data_cae, assamble_stiffness &item_k, vector<double> &solution)
    {
        // 计算缩减后的 系数矩阵
        vector<vector<double>> rk;
        this->ca_reduced_K(item_k, rk);

        // 计算缩减后的载荷向量
        vector<double> rf;
        vector<double> new_f(item_k.num_row_, 0.);
        for (int i = 0; i < data_cae.save_ref_info_.load_vec_ori_.size(); i++)
        {
            new_f[i] = data_cae.save_ref_info_.load_vec_ori_[i];
        }
        this->ca_reduced_F(new_f, rf);
        // 建立约简后的系数矩阵的 SCR
        vector<double> nz_val(this->n_basis_ * this->n_basis_, 0.);
        vector<int> row_idx(this->n_basis_ * this->n_basis_, 0.);
        vector<int> col_idx(this->n_basis_ + 1, 0.);
        int id_;
        for (int i = 0; i < this->n_basis_; i++)
        {
            for (int j = 0; j < this->n_basis_; j++)
            {
                id_ = i + j * this->n_basis_; // 列压缩，“对称，行/列压缩无差别”
                nz_val[id_] = rk[i][j];
                row_idx[id_] = i;
            }
            col_idx[i + 1] = (i + 1) * this->n_basis_;
        }
        // 求解 降阶方程
        vector<double> rx(this->n_basis_, 0.);
        superlu_solver_func(nz_val, row_idx, col_idx, rf, rx);

        // 线性映射 降阶解 到 完全解
        int row = int(ROM_[0].size());
        solution.resize(row);
        for (int i = 0; i < row; i++)
        {
            double sum = 0.;
            for (int j = 0; j < this->n_basis_; j++)
            {
                sum += ROM_[j][i] * rx[j];
            }
            solution[i] = sum;
            // cout<< i+1<<": "<<solution[i]<<endl;
        }
        cout << "CA solution has been finished !!!" << endl;
    };

    // 计算约简后的系数矩阵
    void CA_solution::ca_reduced_K(assamble_stiffness &item_k, vector<vector<double>> &rk)
    {
        int col = int(ROM_.size());
        int row = int(ROM_[0].size());
        rk.resize(col, vector<double>(col, 0.));
        vector<vector<double>> k_dot_rom(col, vector<double>(row, 0.));
        // k_dot_rom = item_k * ca_rom_n
        for (int i = 0; i < col; i++)
        {
            vector<double> temp(row, 0.);
            Sparese_dot_vector(item_k, ROM_[i], temp);
            k_dot_rom[i].assign(temp.begin(), temp.end());
        }
        for (int i = 0; i < col; i++)
        {
            for (int j = 0; j < col; j++)
            {
                rk[i][j] = this->ca_vec_dot_vec(ROM_[i], k_dot_rom[j]);
                // std::cout << rk[i][j] << " ";
            }
            // std::cout << std::endl;
        }
    };

    // 计算约简后的 载荷向量
    void CA_solution::ca_reduced_F(vector<double> &f, vector<double> &rf)
    {
        int col = int(ROM_.size());
        rf.resize(col, 0.);
        for (int i = 0; i < col; i++)
        {
            rf[i] = ca_vec_dot_vec(ROM_[i], f);
        }
    };

    double CA_solution::ca_vec_dot_vec(vector<double> &vec_1, vector<double> &vec_2)
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
    };

    // 计算位移进行处理
    void CA_solution::dis_process(data_management &data_cae)
    {
        // 将位移填充为全位移
        this->ca_full_dis_vec_.resize(3 * this->num_all_nodes_);
        std::fill(this->ca_full_dis_vec_.begin(), this->ca_full_dis_vec_.end(), 0.);
        int dis_idx = 0;
        int ori_length = data_cae.save_ref_info_.resort_free_nodes_ori_.size();
        for (int i = 0; i < this->num_all_nodes_; i++)
        {
            if (i < ori_length)
            {
                if (data_cae.save_ref_info_.resort_free_nodes_ori_[i] >= 0)
                {
                    this->ca_full_dis_vec_[3 * i] = this->ca_dis_vec_[3 * dis_idx];
                    this->ca_full_dis_vec_[3 * i + 1] = this->ca_dis_vec_[3 * dis_idx + 1];
                    this->ca_full_dis_vec_[3 * i + 2] = this->ca_dis_vec_[3 * dis_idx + 2];
                    dis_idx += 1;
                }
            }
            else
            {
                this->ca_full_dis_vec_[3 * i] = this->ca_dis_vec_[3 * dis_idx];
                this->ca_full_dis_vec_[3 * i + 1] = this->ca_dis_vec_[3 * dis_idx + 1];
                this->ca_full_dis_vec_[3 * i + 2] = this->ca_dis_vec_[3 * dis_idx + 2];
                dis_idx += 1;
            }
        }
        // 删除被删除节点位移
        int num_node_new = data_cae.mesh_ca_new_.mdf_coords_.size();
        this->extract_dis_vec_.resize(3 * num_node_new);
        std::fill(this->extract_dis_vec_.begin(), this->extract_dis_vec_.end(), 0.);
        int node_map_id;
        for (int i = 0; i < this->num_all_nodes_; i++)
        {
            node_map_id = data_cae.mesh_ca_map_.node_idx_union_map_[i];
            if (node_map_id >= 0)
            {
                this->extract_dis_vec_[3 * node_map_id] = this->ca_full_dis_vec_[3 * i];
                this->extract_dis_vec_[3 * node_map_id + 1] = this->ca_full_dis_vec_[3 * i + 1];
                this->extract_dis_vec_[3 * node_map_id + 2] = this->ca_full_dis_vec_[3 * i + 2];
            }
        }
    };
}