/**************************************************************************

Copyright:  WH team

Author: YinJichao <jichaoyinyjc@163.com>

Completion date:  XXX

Description: XXX

**************************************************************************/

#include "include/cae.h"

namespace CAE
{

    void CAE_process::Init(int nargs, char *argv[])
    {
        // this->mat_.young_modulus = 21000.0;
        // this->mat_.poisson_ratio = 0.3;
        // 处理命令行
        processCmdLine(nargs, argv);

        // 读入求解文件
        readFile();
    }

    void CAE_process::processCmdLine(int nargs, char *argv[])
    {
        // 第一个命令为工作目录
        char *work_path = argv[0];
        option_.work_path = string(work_path);

        if (nargs == 2)
        {
            char *file_name = argv[1];
            option_.ifile_name = string(file_name);
            path_ = string(file_name);
        }
        else
        {
            // 依次处理后续命令符
            char *sz;
            bool ifile_on = false;
            bool ofile_on = false;
            bool cafile_on = false;
            for (int i = 1; i < nargs; i++)
            {
                sz = argv[i];
                if (ifile_on)
                {
                    option_.ifile_name = string(sz);
                    path_ = string(sz);
                    ifile_on = false;
                    continue;
                }
                if (ofile_on)
                {
                    option_.ofile_name = string(sz);
                    ofile_on = false;
                    continue;
                }
                if (cafile_on)
                {
                    option_.ofile_name = string(sz);
                    if (option_.ca1file_name == "None")
                        option_.ca1file_name = string(sz);
                    else
                        option_.ca2file_name = string(sz);
                    cafile_on = false;
                    continue;
                }
                if (strcmp(sz, "-imp") == 0)
                {
                    option_.analysis_type = 1;
                }
                else if (strcmp(sz, "-exp") == 0)
                {
                    option_.analysis_type = 2;
                }
                else if (strcmp(sz, "-nl") == 0)
                {
                    data_cae_.NLFEA = true;
                }
                else if (strcmp(sz, "-sfem") == 0)
                {
                    option_.sfem_flag = true;
                }
                else if (strcmp(sz, "-re") == 0)
                {
                    option_.reanalysis_flag = true;
                    // TODO
                    // 因为前处理相关功能不全，可根据后续项目需求完善重分析
                    // 比如找个变量存重分析文件路径等，检查是否文件给全了，用于后面读取。
                    // 以及是否需要加一个命令控制隐式分析是否保存刚度信息。
                }
                else if (strcmp(sz, "-i") == 0)
                {
                    ifile_on = true;
                }
                else if (strcmp(sz, "-o") == 0)
                {
                    ofile_on = true;
                }
                else if (strcmp(sz, "-ca") == 0)
                {
                    cafile_on = true;
                    option_.reanalysis_flag = true;
                }
            }
        }

        // 检查文件后缀和输出设置
        if (option_.ifile_name.empty())
        {
            cout << "Please enter the correct path to the input file ending in .inp" << option_.ifile_name << endl;
        }
        else
        {
            // 找到最后一个'/'或'\'的位置
            size_t last_slash_pos = option_.ifile_name.find_last_of("/\\");
            // 找到最后一个'.'的位置
            size_t last_dot_pos = option_.ifile_name.find_last_of('.');
            // 如果没有找到'/'或'\'，则文件名就是整个路径;如果没有找到'.'，则没有后缀名;或者后缀不是inp
            if (last_slash_pos == std::string::npos || last_dot_pos == std::string::npos || option_.ifile_name.substr(last_dot_pos + 1) != "inp")
            {
                cout << "Please enter the correct path to the input file ending in .inp" << option_.ifile_name << endl;
            }
            if (!option_.ofile_name.empty())
            {
                size_t last_dot_pos_out = option_.ofile_name.find_last_of('.');
                // 如果没有找到'.'，则没有后缀名;或者后缀不是vtk
                if (last_dot_pos == std::string::npos || option_.ofile_name.substr(last_dot_pos_out + 1) != "vtk")
                {
                    cout << "Please enter the correct path to the output file ending in .vtk" << option_.ofile_name << endl;
                }
            }
            else
            {
                option_.ofile_name = option_.ifile_name.substr(0, last_dot_pos) + ".vtk";
            }
        }

        if (option_.ifile_name.find(".inp") == option_.ifile_name.size() - 4)
        {
            option_.file_type = 1;
        }
        else if (option_.ifile_name.find(".fem") == option_.ifile_name.size() - 4)
        {
            option_.file_type = 2;
        }
    }

    void CAE_process::readFile()
    {
        if (option_.file_type == 1)
        { // inp文件读取
            // 关键字
            string load_set_keyword = "Set-load";
            string load_value_keyword = "Cload";
            string dis_set_keyword = "Set-fix";
            // 读取计算文件
            pre_info(load_set_keyword, load_value_keyword, dis_set_keyword);
        }
        else if (option_.file_type == 2)
        { // fem文件读取
        }
        else if (option_.file_type == 3)
        { // re文件读取(重分析文件)
        }
    }

    void CAE_process::Solve()
    {
        if (option_.sfem_flag)
        {
            // SFEM
            this->implict_SFEManalysis(option_.ofile_name, "");
        }
        else if (option_.analysis_type == 1)
        {
            // 隐式
            this->implict_analysis(option_.ofile_name, option_.is_save_stiffness);
            // 重分析
            if (option_.reanalysis_flag)
            {
                // 重分析
                string mesh_path1 = option_.ca1file_name;     // "E:\\WH_CAE\\test_model\\CA_cadcae\\CA_info.txt"; 
                string mesh_path2 = option_.ca2file_name;     // "E:\\WH_CAE\\test_model\\CA_cadcae\\mesh_m.inp";  
                this->CA_pre_process(mesh_path1, mesh_path2); 
                // 开始执行重分析                               
                string CA_result_path = option_.ofile_name;  // "E:\\WH_CAE\\test_model\\output\\cadcae_ca.vtk";
                int n_basis = 4;
                this->CA_ReAnalysis(CA_result_path, n_basis);
            }
        }
        else if (option_.analysis_type == 2)
        {
            // 显式
            option_.ofile_name = option_.ofile_name.substr(0, option_.ofile_name.size() - 4);
            this->explicit_analysis(option_.ofile_name, "");
        }
    }

    void CAE_process::pre_info(string load_set_keyword, string load_value_keyword, string dis_set_keyword)
    {
        ReadInfo item_info(option_.ifile_name);

        // 读取单元、节点总数
        item_info.read_ele_node_num(data_cae_);

        // 读取几何信息
        item_info.read_geo_mesh(data_cae_);

        // 读取非协调信息
        item_info.readNconformingMessage(data_cae_);

        // 读取载荷边界信息
        item_info.read_load_bcs(load_set_keyword, load_value_keyword, data_cae_);

        // 读取位移边界信息
        item_info.read_dis_bcs(dis_set_keyword, data_cae_);

        // 读取材料信息
        item_info.read_mat(mat_);

        // 读取时间信息--显式
        if (option_.analysis_type == 2)
        {
            item_info.read_time(data_cae_);
        }
    }

    // 执行结构响应分析
    void CAE_process::implict_analysis(string result_path, bool is_save_stiffness)
    {
        set_BCs item_bcs;

        // 设置边界条件
        item_bcs.build_free_index(data_cae_);

        // 建立单载荷向量
        item_bcs.build_single_load(data_cae_);

        // 几何非线性判断
        if (data_cae_.NLFEA)
        {
            double res_norm = 1.;
            int num_free_nodes = data_cae_.nd_ - data_cae_.dis_bc_set_.size();
            // 声明刚度矩阵组装对象, 刚度矩阵的索引不变，因此仅组装一次
            assamble_nl_stiffness item_nl_assam;
            item_nl_assam.build_CSR(data_cae_);
            vector<double> current_dis_vec(3 * num_free_nodes, 0.);
            vector<double> inter_force_vec(3 * num_free_nodes, 0.);
            vector<double> res_vec(3 * num_free_nodes, 0.);
            int loop = 0;
            while (res_norm > data_cae_.res_lmit && loop < 100)
            {
                std::fill(inter_force_vec.begin(), inter_force_vec.end(), 0.);
                item_nl_assam.fill_CSR_sparse_mat(data_cae_, mat_, current_dis_vec, inter_force_vec);
                // 外力-内力
                for (int i = 0; i < 3 * num_free_nodes; i++)
                {
                    res_vec[i] = data_cae_.single_load_vec_[i] - inter_force_vec[i];
                }
                // 求解
                vector<double> incre_dis_vec;
                incre_dis_vec.resize(3 * num_free_nodes);
                string type_solver = "Pardiso_class"; // "SuperLU_func"; "SuperLU_class"; "Pardiso_func"; "Pardiso_class"; "CA"
                solution_nl_api(item_nl_assam, data_cae_, res_vec, incre_dis_vec, type_solver);

                // 位移修正
                res_norm = 0.;
                for (int i = 0; i < 3 * num_free_nodes; i++)
                {
                    current_dis_vec[i] = current_dis_vec[i] + incre_dis_vec[i];
                    res_norm = res_norm + res_vec[i] * res_vec[i];
                }
                res_norm = sqrt(res_norm) / (3 * num_free_nodes);
                loop++;
                cout << "this is loop " << loop << ", and res is " << res_norm << endl;
            }
            // 将当前位移赋到data_cae_成员变量中
            cout << "the loop is over " << endl;
            data_cae_.single_dis_vec_.resize(3 * num_free_nodes);
            for (int i = 0; i < 3 * num_free_nodes; i++)
            {
                data_cae_.single_dis_vec_[i] = current_dis_vec[i];
            }
        }
        else
        {
            clock_t start, end; // 定义clock_t变量
            start = clock();

            // 声明刚度矩阵组装 对象
            assamble_stiffness item_assam;
            // 非协调判断
            if (data_cae_.BndMesh_F.empty())
            {
                item_assam.build_CSR(data_cae_);
                item_assam.fill_CSR_sparse_mat(data_cae_, mat_);
            }
            else
            {
                item_assam.NCF_assembleStiffness(data_cae_, mat_); // 12.8非协调面自由度索引有问题还未改完
            }
            // 求解
            int num_free_nodes = data_cae_.nd_ - data_cae_.dis_bc_set_.size();
            data_cae_.single_dis_vec_.resize(3 * num_free_nodes);
            string type_solver = "Pardiso_class"; // "SuperLU_class"; "Pardiso_class"; "SuperLU_func"; "Pardiso_func";  "CA"
            solution_api(item_assam, data_cae_, type_solver);
            if (is_save_stiffness)
            {
                Save_stiffness(item_assam);
            }
            end = clock();
            cout << "It took " << double(end - start) / CLOCKS_PER_SEC << " s to compute full analysis" << endl;
        }

        // 填充位移
        simulation_post post_item;
        post_item.reset_displacement(data_cae_);
        // 计算应力
        post_item.get_cauchy_stress_3d(data_cae_);
        // 输出物理场
        data_process item_output;
        double scale_dis = 1.0;
        item_output.export_dis_2_vtk(data_cae_, result_path, scale_dis);
        if (is_save_stiffness)
        {
            // move 原始节点坐标
            data_cae_.single_dis_vec_o_ = std::move(data_cae_.single_dis_vec_);
            data_cae_.resort_free_nodes_o_ = std::move(data_cae_.resort_free_nodes_);
            data_cae_.coords_o_ = std::move(data_cae_.coords_);
            data_cae_.single_load_vec_o_ = std::move(data_cae_.single_load_vec_);
        }
    }

    void CAE_process::CA_pre_process(string mesh_path, string node_now_path)
    {
        // 读取网格信息
        path_ = mesh_path;
        ReadInfo item_info(path_, node_now_path);

        // 读取删除单元拓扑信息
        item_info.CA_read_change(data_cae_);

        // 读取现有模型节点
        item_info.CA_read_now_node(data_cae_);
    }

    void CAE_process::CA_ReAnalysis(string result_path, int n_basis, bool is_Update)
    {
        // 重分析变量初始化
        clock_t start, end; // 定义clock_t变量
        vector<double> ca_solution;
        assamble_stiffness delt_K;
        assamble_stiffness current_K;

        // 判断是否有新增单元的场景
        bool incre_flag = false;
        if (incre_flag)
        {
            // TODO 当有单元增加的时候，修改参考刚度矩阵索引：row_idx_，col_idx_
        }
        // 当前刚度矩阵无增加的情况
        current_K.num_row_ = data_cae_.item_assam_implicit_.num_row_;
        current_K.num_col_ = data_cae_.item_assam_implicit_.num_col_;
        current_K.num_nz_val_ = data_cae_.item_assam_implicit_.num_nz_val_;
        
        current_K.row_idx_.assign(data_cae_.item_assam_implicit_.row_idx_.begin(), data_cae_.item_assam_implicit_.row_idx_.end());
        current_K.col_idx_.assign(data_cae_.item_assam_implicit_.col_idx_.begin(), data_cae_.item_assam_implicit_.col_idx_.end());
        current_K.nz_val_.resize(current_K.num_nz_val_);
        std::fill(current_K.nz_val_.begin(), current_K.nz_val_.end(), 0.);

        // 刚度矩阵变化量 计时
        start = clock();
        ca_get_delt_stiffness(data_cae_, delt_K, mat_); // 计算delt_K

        // 填充当前刚度矩阵的值
        for (int i = 0; i < current_K.num_nz_val_; i++)
        {
            current_K.nz_val_[i] = data_cae_.item_assam_implicit_.nz_val_[i] + delt_K.nz_val_[i];
        }
        end = clock();
        cout << "It took " << double(end - start) / CLOCKS_PER_SEC << " s to compute the amount of change in the stiffness matrix" << endl;
        
        // 构造组合近似降阶模型计时
        start = clock();
        ca_build_rom(data_cae_, delt_K, n_basis); // 计算组合近似降阶模型
        end = clock();
        cout << "It took " << double(end - start) / CLOCKS_PER_SEC << " s to compute the CA model" << endl;

        data_cae_.single_dis_vec_.clear();
        data_cae_.single_full_dis_vec_.clear();
        
        // 求解计时
        start = clock();
        ca_solve(data_cae_, current_K, ca_solution); // 求解降阶后的模型
        end = clock();
        cout << "It took " << double(end - start) / CLOCKS_PER_SEC << " s to solve the reduced model" << endl;

        // 提取节点位移, 此处位移场仍是参考模型上所有结点位移，在VTK输出中，只输出修改后模型的位移
        simulation_post post_item;
        post_item.reset_ca_displacement(data_cae_, ca_solution);

        // 输出物理场
        data_process item_output;
        double scale_dis = 1.0;
        item_output.CA_export_dis_2_vtk(data_cae_, scale_dis, result_path);

        if (is_Update)
        {
            // TODO:留个接口，将原模型更新为修改后的模型（节点、拓扑、刚度、位移等）
            // 可能以后CAD/CAE会用到
        }
    }

    // 执行结构光滑有限元分析
    void CAE_process::implict_SFEManalysis(string result_path, string path_abaqus)
    {
        //=====建立光滑有限元相关数据=====
        SFEM3D SFEMData(&data_cae_);
        SFEMData.build();

        set_BCs item_bcs;

        // 设置边界条件
        item_bcs.build_free_index(data_cae_);

        // 建立单载荷向量
        item_bcs.build_single_load(data_cae_);

        // 组装刚度矩阵
        assamble_stiffness item_assam;

        item_assam.SFEM_build_CSR(&SFEMData);

        item_assam.SFEM_fill_CSR_sparse_mat(&SFEMData, mat_);

        // 求解
        int num_free_nodes = data_cae_.nd_ - data_cae_.dis_bc_set_.size();
        data_cae_.single_dis_vec_.resize(3 * num_free_nodes);
        string type_solver = "SuperLU_func"; // "SuperLU"; "Pardiso_func"; "Pardiso_class"; "CA"
        solution_api(item_assam, data_cae_, type_solver);

        // 输出物理场
        simulation_post post_item;
        post_item.reset_displacement(data_cae_);
        // 输出物理场
        data_process item_output;
        double scale_dis = 1.0;
        item_output.export_dis_2_vtk(data_cae_, result_path, scale_dis);
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
        double time_scale_ = 0.9; // use to scale the timestep
        int output_gap_ = 10;
        // start explicit solve
        // 1.allocate and inite disp_tp1, disp_t0, vel_tphalf, vel_thalf, acc_t0, InFroce_, OutFroce_, Mass_, stress_, strain_, strain_p_, real_coords_;
        int nnode_ = data_cae_.nd_, nele_ = data_cae_.ne_;
        vector<double> disp_tp1(3 * nnode_), disp_t0(3 * nnode_), disp_d(3 * nnode_), vel_tphalf(3 * nnode_),
            vel_thalf(3 * nnode_), acc_t0(3 * nnode_), InFroce_(3 * nnode_), OutFroce_(3 * nnode_), Mass_(nnode_);
        vector<vector<double>> stress_(nele_), strain_(nele_), strain_p_(nele_), real_coords_(data_cae_.coords_);
        for (int i = 0; i < nele_; i++)
        {
            int idx = data_cae_.ele_map_list_[data_cae_.ele_list_idx_[i]];
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
        for (int node_idx_ : data_cae_.load_set_)
        {
            int node_idx = node_idx_ - 1;
            OutFroce_[3 * node_idx + force_dof] = force_value;
        }
        // 3.update timestep & calcul half time vel_thalf for solve loop
        // 3.1 update timestep(at this step, we can use time_step as parameter directly)
        if (data_cae_.time_step_ == 0.0)
        {
            auto_time_ = true;
            UpdateTimeStep(data_cae_.node_topos_, real_coords_, data_cae_.ele_list_, data_cae_.ele_list_idx_, data_cae_.ele_map_list_, time_step_);
            time_step_ *= time_scale_;
        }
        else
        {
            auto_time_ = false;
            time_step_ = data_cae_.time_step_;
        }
        // 3.2 calcul half time velocity(vel_thalf)
        for (int i = 0; i < nnode_; i++)
        {
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
        int output_count = 0;  // save times of output
        int step_count = 0;    // save times of step
        bool save_VTK = 0;     // parameter of save result
        bool finish_solve = 0; // parameter of finish solve
        // std::string save_file = result_path + "ExplicitCae_result_" + std::to_string(output_count) + ".vtk";
        std::string save_file = result_path + "_" + std::to_string(output_count) + ".vtk";
        // 4.2 save the origin struct
        // 输出物理场
        data_process item_output;
        double scale_dis = 1.0;
        data_cae_.single_full_dis_vec_ = disp_t0;
        item_output.export_dis_2_vtk(data_cae_, save_file, scale_dis);
        // 4.3 change time to first step, swap the result
        time_now_ += time_step_;
        SwapData(disp_tp1, disp_t0, vel_tphalf, vel_thalf, acc_t0, InFroce_, time_step_, time_step_old_, auto_time_);
        // 4.4 start loop
        while (time_now_ <= data_cae_.time_total_)
        {
            step_count++;
            // 4.4.1 update displacement by last step result
            // disp_d = vel_tphalf * time_step_old
            // disp_t0 = disp_tp1 + disp_d
            for (int i = 0; i < disp_t0.size(); i++)
            {
                disp_d[i] = vel_tphalf[i] * time_step_old_;
                disp_t0[i] = disp_tp1[i] + disp_d[i];
            }
            // 4.4.2 boundary condiction
            for (int node_idx_ : data_cae_.dis_bc_set_)
            {
                int node_idx = node_idx_ - 1;
                disp_d[3 * node_idx] = 0;
                disp_d[3 * node_idx + 1] = 0;
                disp_d[3 * node_idx + 2] = 0;
                disp_t0[3 * node_idx] = 0;
                disp_t0[3 * node_idx + 1] = 0;
                disp_t0[3 * node_idx + 2] = 0;
            }
            // 4.4.3 iterate over all elements to update Infroce
            for (int i = 0; i < nele_; i++)
            {
                int idx = data_cae_.ele_map_list_[data_cae_.ele_list_idx_[i]];
                data_cae_.ele_list_[idx]->cal_in_force(data_cae_.node_topos_[i], data_cae_.coords_, disp_d, stress_[i], strain_[i], InFroce_);
            }
            // 4.4.4 update timestep
            // update real_coords
            update_coords(data_cae_.coords_, disp_t0, real_coords_);
            if (auto_time_)
            {
                UpdateTimeStep(data_cae_.node_topos_, real_coords_, data_cae_.ele_list_, data_cae_.ele_list_idx_, data_cae_.ele_map_list_, time_step_);
                time_step_ *= time_scale_;
            }
            // check wether if time_now + timestep > time_total or next output time
            // next output time = (time_total / output_gap_) * (output_count + 1)
            if (CheckTime(time_step_, (data_cae_.time_total_ / output_gap_) * (output_count + 1) - time_now_))
            {
                output_count++;
            }
            CheckTime(time_step_, data_cae_.time_total_ - time_now_);
            // 4.4.5 calcul half time velocity(vel_thalf)
            for (int i = 0; i < nnode_; i++)
            {
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
            if (time_now_ == data_cae_.time_total_ || time_now_ == (data_cae_.time_total_ / output_gap_) * (output_count))
            {
                // save_file = result_path + "ExplicitCae_result_" + std::to_string(output_count) + ".vtk";
                save_file = result_path + "_" + std::to_string(output_count) + ".vtk";
                data_cae_.single_full_dis_vec_ = disp_t0;
                item_output.export_dis_2_vtk(data_cae_, save_file, scale_dis);
                if (time_now_ == data_cae_.time_total_)
                {
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
    void CAE_process::Save_stiffness(assamble_stiffness &item_assam)
    {
        // 存储刚度矩阵，用于重分析
        data_cae_.item_assam_implicit_.num_row_ = item_assam.num_row_;
        data_cae_.item_assam_implicit_.num_col_ = item_assam.num_col_;
        data_cae_.item_assam_implicit_.num_nz_val_ = item_assam.num_nz_val_;
        data_cae_.item_assam_implicit_.nz_val_ = std::move(item_assam.nz_val_);
        data_cae_.item_assam_implicit_.row_idx_ = std::move(item_assam.row_idx_);
        data_cae_.item_assam_implicit_.col_idx_ = std::move(item_assam.col_idx_);
    }
}