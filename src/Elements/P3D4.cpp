///**************************************************************************
//
//Copyright:  WH team
//
//Author: YinJichao <jichaoyinyjc@163.com>、ChenBinXiang <chen_mech99@163.com>
//
//Completion date:  XXX
//
//Description: XXX
//
//**************************************************************************/
//
//#pragma once
//
//#include "include/P3D4.h"
//
//namespace CAE
//{
//    REGISTER(ele_base, plate4_ele_elastic, "P3D4");
//    // 建立本构矩阵
//    void plate4_ele_elastic::build_cons_mat()
//    {
//
//        double factor = E / (1 - v * v);
//
//        std::vector<std::vector<double>> D(3, std::vector<double>(3));
//
//        D[0][0] = factor;
//        D[0][1] = factor * v;
//        D[0][2] = 0;
//
//        D[1][0] = factor * v;
//        D[1][1] = factor;
//        D[1][2] = 0;
//
//        D[2][0] = 0;
//        D[2][1] = 0;
//        D[2][2] = factor * (1 - v) / 2;
//
//        double em = matrial_struc_.young_modulus;
//        double nu = matrial_struc_.poisson_ratio;
//        double fac = em / (1 - nu * nu);
//        double fac_a = fac * nu;
//        double fac_b = fac * (1.0 - nu) / 2.0;
//        // 本构矩阵赋值
//        C_matrix_ << fac, fac_a, 0.,
//            fac_a, fac, 0.,
//            0., 0., fac_b;
//    }
//
//    // 建立应变矩阵
//    void plate4_ele_elastic::build_strain_mat(Eigen::Ref<Eigen::MatrixXd> node_coords, Eigen::Ref<Eigen::MatrixXd> strain_mat)
//    {
//        // TODO 后期加入关键词或者参数，控制 积分方式选择
//        det_jacobi_ = build_strain_mat_hammer(node_coords, strain_mat);
//    }
//
//    // 建立单元刚度矩阵
//    void plate4_ele_elastic::build_ele_stiff_mat(Eigen::Ref<Eigen::MatrixXd> node_coords, Eigen::Ref<Eigen::MatrixXd> stiffness_matrix)
//    {
//        det_jacobi_ = build_ele_stiff_mat_hammer(node_coords, stiffness_matrix, C_matrix_);
//        if (det_jacobi_ < 0) {
//            std::cout << "det_jacobi<0" << std::endl;
//            count_det_++;
//            std::cout << count_det_ << std::endl;
//        }
//    }
//
//    // 建立切线单元刚度矩阵（几何非线性）
//    void plate4_ele_elastic::build_ele_nl_stiff_mat(Eigen::Ref<Eigen::MatrixXd> node_coords, Eigen::Ref<Eigen::MatrixXd> node_dis,
//                                                   Eigen::Ref<Eigen::MatrixXd> stiffness_matrix, Eigen::Ref<Eigen::MatrixXd> inter_force)
//    {
//        build_ele_nl_stiff_mat_hammer(node_coords, node_dis, stiffness_matrix, inter_force, C_matrix_);
//    }
//
//    // 计算节点处位移应变转换矩阵
//    void plate4_ele_elastic::build_strain_node_mat(Eigen::Ref<Eigen::MatrixXd> node_coords, Eigen::Ref<Eigen::MatrixXd> strain_mat)
//    {
//        build_strain_node_mat_hammer(node_coords, strain_mat);
//    }
//
//    // 计算节点处应力
//    void plate4_ele_elastic::get_stress_node(Eigen::Ref<Eigen::MatrixXd> node_coords, Eigen::Ref<Eigen::MatrixXd> stress_mat, Eigen::Ref<Eigen::MatrixXd> dis_vec)
//    {
//        get_stress_node_hammer(node_coords, stress_mat, C_matrix_, dis_vec);
//    }
//
//    // 建立单元质量矩阵
//    void plate4_ele_elastic::build_ele_mass(const vector<int> &node_topos, const vector<vector<double>> &coords, vector<double> &Mass)
//    {
//        // define & initialize variables
//        double volume, mass;
//        Eigen::Matrix4d nodes_coor;
//        nodes_coor.setZero();
//        int item_node;
//        // trans nodes coor to a matrix
//        for (int i = 0; i < 4; i++)
//        {
//            item_node = node_topos[i] - 1;
//            nodes_coor(i, 0) = coords[item_node][0];
//            nodes_coor(i, 1) = coords[item_node][1];
//            nodes_coor(i, 2) = coords[item_node][2];
//            nodes_coor(i, 3) = 1.0;
//        }
//        // calculate the element volume & mass = volums * density
//        volume = std::abs(nodes_coor.determinant() / 6.0);
//        mass = volume * this->matrial_struc_.density;
//        // assemble the Global mass vector
//        for (int i = 0; i < 4; i++)
//        {
//            Mass[node_topos[i] - 1] = Mass[node_topos[i] - 1] + mass / 4.0;
//        }
//    }
//
//    // 计算单元内力
//    void plate4_ele_elastic::cal_in_force(const vector<int> &node_topos, const vector<vector<double>> &real_coords, const vector<double> &disp_d,
//                                         vector<double> &stress, vector<double> &strain, vector<double> &InFroce)
//    {
//        Matrix4d3 nodes_coor;
//        nodes_coor.setZero();
//        int item_node;
//        // trans nodes coor to a matrix
//        for (int i = 0; i < 4; i++)
//        {
//            item_node = node_topos[i] - 1;
//            nodes_coor(i, 0) = real_coords[item_node][0];
//            nodes_coor(i, 1) = real_coords[item_node][1];
//            nodes_coor(i, 2) = real_coords[item_node][2];
//        }
//        Eigen::MatrixXd stiffness_matrix, disp, dstrain, inforce;
//        stiffness_matrix.resize(12, 12);
//        disp.resize(12, 1);
//        int idx_temp;
//        for (int i = 0; i < 4; i++)
//        {
//            idx_temp = node_topos[i] - 1;
//            disp(3 * i, 0) = disp_d[3 * idx_temp];
//            disp(3 * i + 1, 0) = disp_d[3 * idx_temp + 1];
//            disp(3 * i + 2, 0) = disp_d[3 * idx_temp + 2];
//        }
//        build_ele_stiff_mat(nodes_coor, stiffness_matrix);
//        inforce = stiffness_matrix * disp; // 12*12 12*1 = 12*1
//        for (int i = 0; i < 4; i++)
//        {
//            idx_temp = node_topos[i] - 1;
//            // 内力为负，故减去
//            InFroce[3 * idx_temp] -= inforce(3 * i, 0);
//            InFroce[3 * idx_temp + 1] -= inforce(3 * i + 1, 0);
//            InFroce[3 * idx_temp + 2] -= inforce(3 * i + 2, 0);
//        }
//    }
//
//    // 计算单元时间步长
//    void plate4_ele_elastic::update_timestep(Eigen::Ref<Eigen::MatrixXd> node_coords, double &time_step)
//    {
//        Eigen::MatrixXd stiffness_matrix;
//        stiffness_matrix.resize(12, 12);
//        build_ele_stiff_mat(node_coords, stiffness_matrix);
//        double max_v = DBL_MIN;
//        double temp;
//        for (int i = 0; i < 12; i++)
//        {
//            temp = 0;
//            for (int j = 0; j < 12; j++)
//            {
//                temp = temp + abs(stiffness_matrix(i, j));
//            }
//            if (max_v < temp)
//            {
//                max_v = temp;
//            }
//        }
//        if (time_step > 2 / sqrt(max_v))
//        {
//            time_step = 2 / sqrt(max_v);
//        }
//    }
//
//    // 交界面积分点物理空间坐标（非协调）
//    void plate4_ele_elastic::gps_phy_coords(Eigen::Ref<Eigen::MatrixXd> nodes1,
//        Eigen::Ref<Eigen::MatrixXd> phy_gps, vector<double>& W_1,
//        vector<Eigen::Vector3d>& Normal)
//    {
//        // hammer积分，权重为1/2
//        double hm_values = 1. / 3;
//        Eigen::MatrixXd hps(1, 2); //积分点
//        hps << hm_values, hm_values;
//
//        Eigen::MatrixXd dNdxi_1(2, 3);
//        Eigen::MatrixXd N_1(1, 3);
//        for (int q = 0; q < 1; q++)
//        {
//            dNdxi_1(0, 0) = -1;
//            dNdxi_1(0, 1) = 1;
//            dNdxi_1(0, 2) = 0;
//            dNdxi_1(1, 0) = -1;
//            dNdxi_1(1, 1) = 0;
//            dNdxi_1(1, 2) = 1;
//
//            N_1(0, 0) = 1 - hps(q, 0) - hps(q, 1);
//            N_1(0, 1) = hps(q, 0);
//            N_1(0, 2) = hps(q, 1);
//
//            Eigen::MatrixXd Jac = dNdxi_1 * nodes1;
//            Eigen::Vector3d a1 = Jac.row(0);
//            Eigen::Vector3d a2 = Jac.row(1);
//            Eigen::Vector3d a3 = a1.cross(a2);
//            double norm_a3 = a3.norm();
//            Eigen::Vector3d unit_a3 = a3.normalized();
//
//            phy_gps.row(q) = N_1 * nodes1;
//            W_1[q] = norm_a3;
//            Normal[q] = unit_a3;
//                
//        }
//
//    }
//
//    
//    void plate4_ele_elastic::text_gps_phy_coords(Eigen::Ref<Eigen::MatrixXd> nodes1,
//        vector<vector<double>>& text_gps, vector<double>& text_W_1,
//        vector<Eigen::Vector3d>& text_Normal)
//    {
//        // hammer积分，权重为1/2
//        double hm_values = 1. / 3;
//        Eigen::MatrixXd hps(1, 2); //积分点
//        hps << hm_values, hm_values;
//
//        Eigen::MatrixXd dNdxi_1(2, 3);
//        Eigen::MatrixXd N_1(1, 3);
//        Eigen::MatrixXd result(1, 3);
//        for (int q = 0; q < 1; q++)
//        {
//            dNdxi_1(0, 0) = -1;
//            dNdxi_1(0, 1) = 1;
//            dNdxi_1(0, 2) = 0;
//            dNdxi_1(1, 0) = -1;
//            dNdxi_1(1, 1) = 0;
//            dNdxi_1(1, 2) = 1;
//
//            N_1(0, 0) = 1 - hps(q, 0) - hps(q, 1);
//            N_1(0, 1) = hps(q, 0);
//            N_1(0, 2) = hps(q, 1);
//
//            Eigen::MatrixXd Jac = dNdxi_1 * nodes1;
//            Eigen::Vector3d a1 = Jac.row(0);
//            Eigen::Vector3d a2 = Jac.row(1);
//            Eigen::Vector3d a3 = a1.cross(a2);
//            double norm_a3 = a3.norm();
//            Eigen::Vector3d unit_a3 = a3.normalized();
//            
//            result = N_1 * nodes1;
//            text_W_1.push_back( norm_a3);
//            text_Normal.push_back(unit_a3);
//            text_gps.emplace_back(vector<double>{result(0, 0), result(0, 1), result(0, 2)});
//
//            /*vector<double> result1(3);
//            result1[0] = result(0, 0);
//            result1[1] = result(0, 1);
//            result1[2] = result(0, 2);
//
//            text_gps.push_back(result1);*/
//            //text_gps.push_back( N_1 * nodes1);
//        }
//
//    }
//
//}