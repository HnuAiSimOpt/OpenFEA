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
//#include "ele_base.h"
//
//namespace CAE
//{
//    class plate4_ele_elastic : public ele_base
//    {
//    public:
//        double det_jacobi_; // 雅可比矩阵行列式(中心积分点)
//        elastic_mat matrial_struc_;
//        Matrix3d3 C_matrix_;
//        Matrix6d12 strain_mat;
//        int count_det_ = 0;
//
//    public:
//        // 构造函数，析构函数
//
//        plate4_ele_elastic() { type_ = "P3D4"; node_dof_ = 3; nnode_ = 4; ngps_ = 4; face_node = 4; face_gps = 4; };
//        plate4_ele_elastic(elastic_mat matrial_struc) : matrial_struc_(matrial_struc) { type_ = "P3D4"; node_dof_ = 3; nnode_ = 4; ngps_ = 4; face_node = 4; face_gps = 4; };
//
//
//        // 材料赋属性
//        void set_matrial(elastic_mat matrial_struc) override { matrial_struc_ = matrial_struc; };
//
//        // 建立本构矩阵
//        virtual void build_cons_mat();
//
//        // 建立位移应变转换矩阵
//        virtual void build_strain_mat(Eigen::Ref<Eigen::MatrixXd> node_coords, Eigen::Ref<Eigen::MatrixXd> strain_mat);
//
//        // 计算节点处位移应变转换矩阵
//        virtual void build_strain_node_mat(Eigen::Ref<Eigen::MatrixXd> node_coords, Eigen::Ref<Eigen::MatrixXd> strain_mat);
//
//        // 计算节点处应力
//        virtual void get_stress_node(Eigen::Ref<Eigen::MatrixXd> node_coords, Eigen::Ref<Eigen::MatrixXd> stress_mat, Eigen::Ref<Eigen::MatrixXd> dis_vec);
//
//        // 建立单元刚度矩阵
//        void build_ele_stiff_mat(Eigen::Ref<Eigen::MatrixXd> node_coords, Eigen::Ref<Eigen::MatrixXd> stiffness_matrix) override;
//
//        // 建立切线单元刚度矩阵（几何非线性）
//        void build_ele_nl_stiff_mat(Eigen::Ref<Eigen::MatrixXd> node_coords, Eigen::Ref<Eigen::MatrixXd> node_dis,
//                                    Eigen::Ref<Eigen::MatrixXd> stiffness_matrix, Eigen::Ref<Eigen::MatrixXd> inter_force) override;
//
//        // 建立单元质量矩阵
//        void build_ele_mass(const vector<int> &node_topos, const vector<vector<double>> &coords, vector<double> &Mass) override;
//
//        // 计算单元内力
//        void cal_in_force(const vector<int> &node_topos, const vector<vector<double>> &real_coords, const vector<double> &disp_d,
//                          vector<double> &stress, vector<double> &strain, vector<double> &InFroce) override;
//
//        // 计算单元时间步长
//
//        void update_timestep(Eigen::Ref<Eigen::MatrixXd> node_coords, double& time_step) override;
//
//        // 交界面积分点物理空间坐标（非协调）
//        void gps_phy_coords(Eigen::Ref<Eigen::MatrixXd> nodes1,
//            Eigen::Ref<Eigen::MatrixXd> phy_gps, vector<double>& W_1,
//            vector<Eigen::Vector3d>& Normal) override;
//
//         void text_gps_phy_coords(Eigen::Ref<Eigen::MatrixXd> nodes1,
//            vector<vector<double>>& text_gps, vector<double>& text_W_1,
//            vector<Eigen::Vector3d>& text_Normal) override;
//    };
//}