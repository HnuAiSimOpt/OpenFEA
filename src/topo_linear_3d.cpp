/**************************************************************************

Copyright:  WH team

Author: YinJichao <jichaoyinyjc@163.com>

Completion date:  XXX

Description: XXX

**************************************************************************/

#include "include/topo_linear_3d.h"

namespace CAE
{

    double TopoLinear3D::topo_linear_cycle_3d(data_management &data_cae, string result_path, int loop_max)
    {
        this->nele_ = data_cae.node_topos_.size();
        vector<double> XDens(this->nele_);
        vector<double> XPhys(this->nele_);
        for (int i = 0; i < this->nele_; i++)
        {
            XDens[i] = vol_;
            XPhys[i] = vol_;
        }
        // 预处理过滤
        cout << "start density filter !!!" << endl;
        SpMatrix H;
        vector<double> Hs;
        this->sen_filter(H, Hs, data_cae.coords_, data_cae.node_topos_);
        cout << "finish density filter !!!" << endl;

        // 设置边界条件
        set_BCs item_bcs;
        item_bcs.build_free_index(data_cae);
        item_bcs.build_single_load(data_cae); // 建立单载荷向量

        // 预计算刚度矩阵索引
        assamble_stiffness topo_assam;
        topo_assam.build_CSR(data_cae);

        // 开始优化循环
        double obj = 0.;
        int loop_idx = 0, change = 1.;
        int num_free_nodes = data_cae.nd_ - data_cae.dis_bc_set_.size();
        data_cae.single_dis_vec_.clear();
        data_cae.single_dis_vec_.resize(3 * num_free_nodes);
        // 初始化
        vector<double> XPhys_simp(this->nele_);
        vector<double> XPhys_sen(this->nele_);
        vector<double> dc_topo(this->nele_);
        vector<double> dc_topo_filter(this->nele_);
        vector<double> dv_topo_filter(this->nele_);
        vector<double> dv_topo(this->nele_);
        vector<double> F_copy(data_cae.single_load_vec_.size(), 0.);
        data_process data2vtk;
        SuperLUSolution topo_solver_splu;
        cout << "The pre-processing has been finished, and the topology optimization start!" << endl;
        while ((loop_idx < loop_max) && (change > 0.01))
        {
            cout << "The <" << loop_idx + 1 << "-th> cycles, ";

            // 根据密度组装 刚度矩阵
            this->simp_density(XPhys, XPhys_simp, false);
            topo_assam.topo_fill_CSR_sparse_mat(data_cae, this->mat_, XPhys_simp);

            // 求解
            data_cae.single_dis_vec_.clear();
            // 初始化
            bool a1 = topo_solver_splu.superlu_init(topo_assam.nz_val_, topo_assam.row_idx_, topo_assam.col_idx_);
            F_copy.assign(data_cae.single_load_vec_.begin(), data_cae.single_load_vec_.end());
            bool a2 = topo_solver_splu.superlu_solution_1st(F_copy, data_cae.single_dis_vec_, true);
            cout << "the solution has been finished, ";
            // 填充位移_位移检查
            simulation_post post_item;
            post_item.reset_displacement(data_cae);
            // write_vec(data_cae.single_dis_vec_, "C:\\Users\\jicha\\Desktop\\CA_Topo_OpenFEA\\OUTPUT\\load.txt");
            // write_vec(data_cae.single_full_dis_vec_, "C:\\Users\\jicha\\Desktop\\CA_Topo_OpenFEA\\OUTPUT\\dis.txt");

            // 敏度分析
            this->simp_density(XPhys, XPhys_sen, true);
            obj = this->sen_analysis(XPhys_sen, XPhys_simp, dc_topo, data_cae);
            cout << "the sensitivity analysis has been finished, ";
            // write_vec(dc_topo, "C:\\Users\\jicha\\Desktop\\CA_Topo_OpenFEA\\OUTPUT\\dc_topo.txt");
            // write_vec(Hs, "C:\\Users\\jicha\\Desktop\\CA_Topo_OpenFEA\\OUTPUT\\Hs.txt");

            // 过滤
            fill(dv_topo.begin(), dv_topo.end(), 1.0);            
            for (int i = 0; i < this->nele_; i++)            
            {
                dc_topo[i] = dc_topo[i] / Hs[i];
                dv_topo[i] = dv_topo[i] / Hs[i];
            }
            H.this_dot_vector(dc_topo, dc_topo_filter);
            H.this_dot_vector(dv_topo, dv_topo_filter);
            write_vec(dc_topo_filter, "C:\\Users\\jicha\\Desktop\\\\dc_topo.txt");
            write_vec(dv_topo_filter, "C:\\Users\\jicha\\Desktop\\\\dv_topo.txt");

            // 变量更新
            this->oc_update(XDens, XPhys, H, Hs, dc_topo_filter, dv_topo_filter);
            cout << "the variables have been updated." << endl;
            // write_vec(XDens, "C:\\Users\\jicha\\Desktop\\CA_Topo_OpenFEA\\OUTPUT\\XDens.txt");
            // write_vec(XPhys, "C:\\Users\\jicha\\Desktop\\CA_Topo_OpenFEA\\OUTPUT\\XPhys.txt");

            // 打印信息
            loop_idx += 1;
            cout << obj << endl;
            std::string name_ = std::to_string(loop_idx);
            data2vtk.export_topo2vtk(data_cae, XPhys, 0.0, "C:\\Users\\jicha\\Desktop\\" + name_ + ".vtk");
        }
        // 保存优化结果，返回优化目标
        data2vtk.export_topo2vtk(data_cae, XPhys, 0.0, "C:\\Users\\jicha\\Desktop\\end.vtk");
        return obj;
    };

    // 计算过滤矩阵
    void TopoLinear3D::sen_filter(SpMatrix &H, vector<double> &Hs, vector<vector<double>> &coords, vector<vector<int>> &node_topos)
    {
        // 计算单元中心坐标
        vector<vector<double>> coor_center(this->nele_, vector<double>(3, 0.));
        for (int i = 0; i < this->nele_; i++)
        {
            double x = 0., y = 0., z = 0.;
            vector<int> node_per = node_topos[i];
            for (int j = 0; j < 4; j++) // 此处应该根据单元的节点数进行循环判断
            {
                x += coords[node_per[j] - 1][0];
                y += coords[node_per[j] - 1][1];
                z += coords[node_per[j] - 1][2];
            }
            coor_center[i][0] = x / 4.;
            coor_center[i][1] = y / 4.;
            coor_center[i][2] = z / 4.;
        }
        // 计算过滤矩阵
        Hs.clear();
        Hs.resize(this->nele_);
        Hs.assign(Hs.size(), 0);                    // 初始化Hs
        vector<std::set<int>> columns(this->nele_); // 稀疏过滤矩阵 H 的行和列
        vector<double> values;                      // 稀疏过滤矩阵 H 的值
        double r_points, weight = 0.;
        ;
        for (int i = 0; i < this->nele_; i++)
        {
            vector<double> node1 = coor_center[i];
            for (int j = 0; j < this->nele_; j++)
            {
                vector<double> node2 = coor_center[j];
                if (abs(node1[0] - node2[0]) < this->rmin_)
                {
                    if (abs(node1[1] - node2[1]) < this->rmin_)
                    {
                        if (abs(node1[2] - node2[2]) < this->rmin_)
                        {
                            r_points = this->distance2points(node1, node2);
                            weight = this->rmin_ - r_points;
                            if (weight > 0.)
                            {
                                columns[i].insert(j); // 因为j自身即位升序，因此set不会打乱排序
                                values.push_back(weight);
                                Hs[i] += weight;
                            }
                        }
                    }
                }
            }
            if(Hs[i]==0)
            {
                cout<<"the rmin is too small.";  // 避免过滤半径过小导致 Hs 为 0， Hs作为分布不能为0
                exit(1);
            }
        }
        H.build_sp(this->nele_, columns, values);
    }

    // 敏度分析
    double TopoLinear3D::sen_analysis(vector<double> &sen_den, vector<double> &filter_den, vector<double> &dc_topo, data_management &data_cae)
    {
        dc_topo.clear();
        dc_topo.resize(this->nele_);
        vector<int> node_topo_per;
        //
        int ele_type, map_idx, node_num_ele, node_idx;
        // 初始化单元
        data_cae.ele_inite(mat_);
        MatrixXd ele_dis;
        MatrixXd ele_coors;
        MatrixXd stiffness_matrix;
        double obj = 0.;
        for (int id_ele = 0; id_ele < this->nele_; id_ele++)
        {
            node_topo_per = data_cae.node_topos_[id_ele];
            // 获取该单元的节点数量
            ele_type = data_cae.ele_list_idx_[id_ele];
            map_idx = data_cae.ele_map_list_[ele_type];
            node_num_ele = data_cae.ele_list_[map_idx]->nnode_;
            // 获取单元位移
            ele_dis.resize(3 * node_num_ele, 1);
            ele_dis.setZero();
            for (int k = 0; k < node_num_ele; k++)
            {
                node_idx = data_cae.node_topos_[id_ele][k] - 1;
                ele_dis(3 * k, 0) = data_cae.single_full_dis_vec_[3 * node_idx];
                ele_dis(3 * k + 1, 0) = data_cae.single_full_dis_vec_[3 * node_idx + 1];
                ele_dis(3 * k + 2, 0) = data_cae.single_full_dis_vec_[3 * node_idx + 2];
                // cout << ele_dis(3 * k, 0)<<"  "<<ele_dis(3 * k + 1, 0)<<"  "<<ele_dis(3 * k + 2, 0)<<endl;
            }
            // 查找节点自由度及坐标, 计算单元刚度矩阵
            ele_coors.resize(node_num_ele, 3);
            build_ele_coors(ele_coors, data_cae, id_ele, node_num_ele);
            stiffness_matrix.resize(3 * node_num_ele, 3 * node_num_ele);
            data_cae.ele_list_[map_idx]->build_ele_stiff_mat(ele_coors, stiffness_matrix);
            // 计算单元的应变能
            MatrixXd ele_sen = (ele_dis.transpose() * stiffness_matrix) * ele_dis;
            dc_topo[id_ele] = sen_den[id_ele] * ele_sen(0, 0);
            obj += filter_den[id_ele] * ele_sen(0, 0);
            // cout<<dc_topo[id_ele]<<endl;
        }
        return obj;
    }

    // 变量更新
    void TopoLinear3D::oc_update(vector<double> &den, vector<double> &filter_den, SpMatrix &H, vector<double> &Hs,
                                 vector<double> &dc_topo, vector<double> &dv_topo)
    {
        vector<double> NewX(this->nele_, 0.);
        //
        filter_den.clear();
        filter_den.resize(this->nele_);
        double Xsum, Lmid, L1 = 0.0, L2 = 1E9, vol_all = this->vol_ * this->nele_;
        while (((L2 - L1) / (L2 + L1)) > 1e-3)
        {
            Lmid = 0.5 * (L1 + L2);
            for (int i = 0; i < this->nele_; i++)
            {
                NewX[i] = DMAX(0.0,
                               DMAX(den[i] - move_,
                                    DMIN(1.0, DMIN(den[i] + move_,
                                                   den[i] * sqrt(dc_topo[i] / dv_topo[i] / Lmid)))));
            }

            Xsum = 0.0;
            H.this_dot_vector(NewX, filter_den);
            for (int i = 0; i < this->nele_; i++)
            {
                filter_den[i] = filter_den[i] / Hs[i];
                Xsum += filter_den[i];
            }
            if (Xsum - vol_all > 0.0)
            {
                L1 = Lmid;
            }
            else
            {
                L2 = Lmid;
            }
            // write_vec(NewX, "aaa.txt");
            // write_vec(filter_den, "bbb.txt");
            // exit(0);
        }
        den.assign(NewX.begin(), NewX.end());
    }

    void TopoLinear3D::build_ele_coors(Eigen::Ref<Eigen::MatrixXd> ele_coors,
                                       data_management &data_cae, int ele_id, int num_nodes)
    {
        ele_coors.resize(num_nodes, 3);
        ele_coors.setZero();
        int node_idx;
        for (int i = 0; i < num_nodes; i++)
        {
            // 坐标
            node_idx = data_cae.node_topos_[ele_id][i] - 1;
            ele_coors(i, 0) = data_cae.coords_[node_idx][0]; // X 坐标
            ele_coors(i, 1) = data_cae.coords_[node_idx][1]; // Y 坐标
            ele_coors(i, 2) = data_cae.coords_[node_idx][2]; // Z 坐标
        }
    }

    // SIMP 密度插值
    void TopoLinear3D::simp_density(const vector<double> &den, vector<double> &simp_den, bool sen)
    {
        simp_den.clear();
        simp_den.resize(this->nele_);
        double one_minus_emin = 1. - e_min_;
        if (sen)
        {
            for (int i = 0; i < this->nele_; i++)
            {
                simp_den[i] = penal_ * one_minus_emin * pow(den[i], penal_ - 1.);
            }
        }
        else
        {
            for (int i = 0; i < this->nele_; i++)
            {
                simp_den[i] = e_min_ + one_minus_emin * pow(den[i], penal_);
            }
        }
    }

    // 计算两点之间的距离
    double TopoLinear3D::distance2points(vector<double> &node1, vector<double> &node2)
    {
        double r;
        double r2 = pow((node1[0] - node2[0]), 2) +
                    pow((node1[1] - node2[1]), 2) +
                    pow((node1[2] - node2[2]), 2);
        r = sqrt(r2);
        return r;
    }
}