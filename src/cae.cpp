/**************************************************************************

Copyright:  WH team

Author: YinJichao <jichaoyinyjc@163.com>

Completion date:  XXX

Description: XXX

**************************************************************************/

#include "include/cae.h"

namespace CAE
{
    void CAE_process::pre_info(string load_set_keyword, string load_value_keyword, string dis_set_keyword)
    {
        ReadInfo item_info(path_);

        // 读取单元、节点总数
        item_info.read_ele_node_num(data_cae_);

        // 读取几何信息
        item_info.read_geo_mesh(data_cae_);

        //读取非协调信息
        item_info.readNconformingMessage(data_cae_);

        // 读取载荷边界信息
        item_info.read_load_bcs(load_set_keyword, load_value_keyword, data_cae_);

        // 读取位移边界信息
        item_info.read_dis_bcs(dis_set_keyword, data_cae_);
    }

    // 执行结构响应分析
    void CAE_process::implict_analysis(string result_path, string path_abaqus)
    {
        set_BCs item_bcs;

        // 设置边界条件
        item_bcs.build_free_index(data_cae_);

        // 建立单载荷向量
        item_bcs.build_single_load(data_cae_);

        // 组装刚度矩阵
        assamble_stiffness item_assam;
        

        //判断是否为非协调
        if (data_cae_.BndMesh_F.empty())
        {  //不是非协调
            item_assam.build_CSR(data_cae_);
            item_assam.fill_CSR_sparse_mat(data_cae_, mat_);
        }
        else
        {  //是非协调
            item_assam.NCF_assembleStiffness(data_cae_, mat_);
        }
        
        // 求解
        int num_free_nodes = data_cae_.nd_ - data_cae_.dis_bc_set_.size();
        data_cae_.single_dis_vec_.resize(3 * num_free_nodes);
        string type_solver = "Pardiso_class"; // "SuperLU"; "Pardiso_func"; "Pardiso_class"; "CA"
        solution_api(item_assam, data_cae_, type_solver);

        // 输出物理场
        data_process item_output;
        item_output.fill_full_dis(data_cae_);
        double scale_dis = 1.0;
        item_output.export_dis_2_vtk(data_cae_, result_path, scale_dis, path_abaqus, false);

        // 存储刚度矩阵，用于重分析
        // TODO:存成文件？pum是否能访问上一次分析的data_managrement?
        data_cae_.item_assam_implicit.num_row = item_assam.num_row;
        data_cae_.item_assam_implicit.num_col = item_assam.num_col;
        data_cae_.item_assam_implicit.num_nz_val = item_assam.num_nz_val;
        data_cae_.item_assam_implicit.nz_val = std::move(item_assam.nz_val);
        data_cae_.item_assam_implicit.row_idx = std::move(item_assam.row_idx);
        data_cae_.item_assam_implicit.col_idx = std::move(item_assam.col_idx);
    }

    void CAE_process::CA_pre_process(string CA_del_set_keyword, vector<int>& del_topo)
    {
        // 读取网格信息
        ReadInfo item_info(path_);

        // 读取删除单元拓扑信息
        item_info.CA_read_del(CA_del_set_keyword, del_topo);
        // TODO：
        // 1、读取修改节点信息
        // 2、读取增加节点及单元拓扑信息

        // 处理CA数据
        // item_info.CA_data_convert(data_cae_, del_topo, node_del_idx, topo_del_idx);
    }

    void CAE_process::CA_ReAnalysis(string result_path, string path_abaqus, vector<int>& del_topo, bool Is_Update)
    {
        // 重分析
        clock_t start, end;     //定义clock_t变量
        int n_basis = 4;
        vector<double> ca_solution;
        assamble_stiffness delt_K;
        assamble_stiffness ref_K;
        assamble_stiffness current_K;
        // 参考刚度矩阵
        ref_K.num_row = data_cae_.item_assam_implicit.num_row;
        ref_K.num_col = data_cae_.item_assam_implicit.num_col;
        ref_K.num_nz_val = data_cae_.item_assam_implicit.num_nz_val;
        ref_K.nz_val = std::move(data_cae_.item_assam_implicit.nz_val);
        ref_K.row_idx = std::move(data_cae_.item_assam_implicit.row_idx);
        ref_K.col_idx = std::move(data_cae_.item_assam_implicit.col_idx);
        // 当前刚度矩阵
        current_K.num_row = data_cae_.item_assam_implicit.num_row;
        current_K.num_col = data_cae_.item_assam_implicit.num_col;
        current_K.num_nz_val = data_cae_.item_assam_implicit.num_nz_val;
        current_K.nz_val.resize(current_K.num_nz_val);
        std::fill(current_K.nz_val.begin(), current_K.nz_val.end(), 0.);
        current_K.row_idx.assign(ref_K.row_idx.begin(), ref_K.row_idx.end());
        current_K.col_idx.assign(ref_K.col_idx.begin(), ref_K.col_idx.end());
        // 刚度矩阵变化量 计时
        start = clock();
        ca_get_delt_stiffness(data_cae_, ref_K, delt_K, mat_, del_topo); // 计算delt_K
        // 填充当前刚度矩阵的值
        for(int i = 0; i < current_K.num_nz_val; i++)
        {
            current_K.nz_val[i] = ref_K.nz_val[i] + delt_K.nz_val[i];
        }
        end = clock();
        cout<<"It took "<<double(end-start)/CLOCKS_PER_SEC<<" s to compute the amount of change in the stiffness matrix"<<endl; 
        // 构造组合近似降阶模型计时
        start = clock();
        ca_build_rom(data_cae_, delt_K, n_basis);                   // 计算组合近似降阶模型
        end = clock();
        cout<<"It took "<<double(end-start)/CLOCKS_PER_SEC<<" s to compute the CA model"<<endl; 
        data_cae_.single_dis_vec_.clear();
        data_cae_.single_full_dis_vec_.clear();
        // 求解计时
        start = clock();
        ca_solve(data_cae_, current_K, ca_solution);                          // 求解降阶后的模型
        end = clock();
        cout<<"It took "<<double(end-start)/CLOCKS_PER_SEC<<" s to solve the reduced model"<<endl; 
        // 提取节点位移
        vector<double> full_dis;
        get_real_dis(data_cae_, ca_solution, full_dis);
        cout<<ca_solution.size()<<endl;
        // 写入 TXT
        // fstream f;
        // f.open(path_abaqus, ios::out);
        // for (int i = 0; i < 3 * data_cae_.nd_; i++)
        // {
        //     string s1 = to_string(full_dis[i]);
        //     f << s1 <<"\n";
        // }
        // f.close();
        // cout<<"ending !!!\n";

        // 输出物理场
        data_process item_output;
        double scale_dis = 1.0;
        item_output.CA_export_dis_2_vtk(data_cae_, del_topo, result_path, scale_dis, full_dis);

        if (Is_Update) {
            // TODO:留个接口，将原模型更新为修改后的模型（节点、拓扑、刚度、位移等）
            // 可能以后CAD/CAE会用到
        }
    }

    // 执行结构动态响应分析
    void CAE_process::explicit_analysis(string result_path, string path_abaqus)
    {
        std::cout << "Explicit solving ......\n";
        // prepare explicit data
        // 1.set material
        data_cae_.ele_inite(mat_);
        // 2. inite time vars
        bool auto_time_ = false;
        double time_now_ = 0;
        double time_step_ = DBL_MAX;
        double time_step_old_ = 0;
        double time_scale_ = 0.9;// use to scale the timestep
        int output_gap_ = 10;
        // start explicit solve 
        // 1.allocate and inite disp_tp1, disp_t0, vel_tphalf, vel_thalf, acc_t0, InFroce_, OutFroce_, Mass_, stress_, strain_, strain_p_, real_coords_;
        int nnode_ = data_cae_.nd_, nele_ = data_cae_.ne_;
        vector<double> disp_tp1(3 * nnode_), disp_t0(3 * nnode_), disp_d(3 * nnode_), vel_tphalf(3 * nnode_),
            vel_thalf(3 * nnode_), acc_t0(3 * nnode_), InFroce_(3 * nnode_), OutFroce_(3 * nnode_), Mass_(nnode_);
        vector<vector<double>> stress_(nele_), strain_(nele_), strain_p_(nele_), real_coords_(data_cae_.coords_);
        for (int i = 0; i < nele_; i++) {
            int idx = data_cae_.ele_list_idx_[i];
            vector<double> temp(6 * data_cae_.ele_list_[idx]->ngps_);
            stress_.push_back(temp);
            strain_.push_back(temp);
            strain_p_.push_back(temp);
            // fill value into Mass_ at the same time.
            data_cae_.ele_list_[idx]->build_ele_mass(data_cae_.node_topos_[i], real_coords_, Mass_);
        }
        // 2.fill value into OutFroce_(the external force is not considered to vary with time, so it is a constant.)
        int force_dof = data_cae_.load_dof_ - 1;
        int force_value = data_cae_.load_value_;
        for (int node_idx_:data_cae_.load_set_) {
            int node_idx = node_idx_ - 1;
            OutFroce_[3 * node_idx + force_dof] = force_value;
        }
        // 3.update timestep & calcul half time vel_thalf for solve loop
        // 3.1 update timestep(at this step, we can use time_step as parameter directly)
        if (data_cae_.time_step_ == 0.0) {
            auto_time_ = true;
            UpdateTimeStep(data_cae_.node_topos_, real_coords_, data_cae_.ele_list_, data_cae_.ele_list_idx_, time_step_);
            time_step_ *= time_scale_;
        }
        else {
            auto_time_ = false;
            time_step_ = data_cae_.time_step_;
        }
        // 3.2 calcul half time velocity(vel_thalf) 
        for (int i = 0; i < nnode_; i++) {
            // update acc_t0
            acc_t0[3 * i] = (InFroce_[3 * i] + OutFroce_[3 * i]) / Mass_[i];
            acc_t0[3 * i + 1] = (InFroce_[3 * i + 1] + OutFroce_[3 * i + 1]) / Mass_[i];
            acc_t0[3 * i + 2] = (InFroce_[3 * i + 2] + OutFroce_[3 * i + 2]) / Mass_[i];
            // vel_thalf = vel_tphalf + acc_t0 * (time_step_old + time_step) / 2;
            // at this time, we think of vel_tphalf and time_step_old_ as 0;
            vel_thalf[3 * i] = vel_tphalf[3 * i] + (acc_t0[3 * i] * (time_step_old_ + time_step_) / 2);
            vel_thalf[3 * i + 1] = vel_tphalf[3 * i + 1] + (acc_t0[3 * i + 1] * (time_step_old_ + time_step_) / 2);
            vel_thalf[3 * i + 2] = vel_tphalf[3 * i + 2] + (acc_t0[3 * i + 2] * (time_step_old_ + time_step_) / 2);
        }
        // 4.begin solve loop until time_now = time_total
        // 4.1 defind contor parameters
        int output_count = 0;       // save times of output 
        int step_count = 0;         // save times of step
        bool save_VTK = 0;          // parameter of save result
        bool finish_solve = 0;        // parameter of finish solve
        std::string save_file = result_path + "ExplicitCae_result_" + std::to_string(output_count) + ".vtk";
        // 4.2 save the origin struct  
        // 输出物理场
        data_process item_output;
        double scale_dis = 1.0;
        data_cae_.single_full_dis_vec_ = disp_t0;
        item_output.export_dis_2_vtk(data_cae_, save_file, scale_dis, path_abaqus, false);
        // 4.3 change time to first step, swap the result
        time_now_ += time_step_;
        SwapData(disp_tp1, disp_t0, vel_tphalf, vel_thalf, acc_t0, InFroce_, time_step_, time_step_old_, auto_time_);
        // 4.4 start loop
        while (time_now_ <= data_cae_.time_total_) {
            step_count++;
            // 4.4.1 update displacement by last step result
            // disp_d = vel_tphalf * time_step_old
            // disp_t0 = disp_tp1 + disp_d
            for (int i = 0; i < disp_t0.size(); i++) {
                disp_d[i] = vel_tphalf[i] * time_step_old_;
                disp_t0[i] = disp_tp1[i] + disp_d[i];
            }
            // 4.4.2 boundary condiction
             for (int node_idx_ : data_cae_.dis_bc_set_) {
                int node_idx = node_idx_ - 1;
                disp_d[3 * node_idx] = 0;
                disp_d[3 * node_idx + 1] = 0;
                disp_d[3 * node_idx + 2] = 0;
                disp_t0[3 * node_idx] = 0;
                disp_t0[3 * node_idx + 1] = 0;
                disp_t0[3 * node_idx + 2] = 0;
            }
            // 4.4.3 iterate over all elements to update Infroce
            for (int i = 0; i < nele_; i++) {
                int idx = data_cae_.ele_list_idx_[i];
                data_cae_.ele_list_[idx]->cal_in_force(data_cae_.node_topos_[i], data_cae_.coords_, disp_d, stress_[i], strain_[i], InFroce_);
            }
            // 4.4.4 update timestep
            // update real_coords
            update_coords(data_cae_.coords_, disp_t0, real_coords_);
            if (auto_time_) {
                UpdateTimeStep(data_cae_.node_topos_, real_coords_, data_cae_.ele_list_, data_cae_.ele_list_idx_, time_step_);
                time_step_ *= time_scale_;
            }
            // check wether if time_now + timestep > time_total or next output time
            //next output time = (time_total / output_gap_) * (output_count + 1)
            if (CheckTime(time_step_, (data_cae_.time_total_ / output_gap_) * (output_count + 1) - time_now_)) {
                output_count++;
            }
            CheckTime(time_step_, data_cae_.time_total_ - time_now_);
            // 4.4.5 calcul half time velocity(vel_thalf)
            for (int i = 0; i < nnode_; i++) {
                // update acc_t0
                acc_t0[3 * i] = (InFroce_[3 * i] + OutFroce_[3 * i]) / Mass_[i];
                acc_t0[3 * i + 1] = (InFroce_[3 * i + 1] + OutFroce_[3 * i + 1]) / Mass_[i];
                acc_t0[3 * i + 2] = (InFroce_[3 * i + 2] + OutFroce_[3 * i + 2]) / Mass_[i];
                // vel_thalf = vel_tphalf + acc_t0 * (time_step_old + time_step) / 2;
                vel_thalf[3 * i] = vel_tphalf[3 * i] + (acc_t0[3 * i] * (time_step_old_ + time_step_) / 2);
                vel_thalf[3 * i + 1] = vel_tphalf[3 * i + 1] + (acc_t0[3 * i + 1] * (time_step_old_ + time_step_) / 2);
                vel_thalf[3 * i + 2] = vel_tphalf[3 * i + 2] + (acc_t0[3 * i + 2] * (time_step_old_ + time_step_) / 2);
            }
            // 4.4.6 save result
            if (time_now_ == data_cae_.time_total_ || time_now_ == (data_cae_.time_total_ / output_gap_) * (output_count)) {
                save_file = result_path + "ExplicitCae_result_" + std::to_string(output_count) + ".vtk";
                data_cae_.single_full_dis_vec_ = disp_t0;
                item_output.export_dis_2_vtk(data_cae_, save_file, scale_dis, path_abaqus, false);
                if (time_now_ == data_cae_.time_total_) {
                    break;
                }
            }
            // 4.4.7 change time to first step, swap the result
            time_now_ += time_step_;
            SwapData(disp_tp1, disp_t0, vel_tphalf, vel_thalf, acc_t0, InFroce_, time_step_, time_step_old_, auto_time_);
        }
        // finish solve
        std::cout << "Explicit finished\n";
    }
}