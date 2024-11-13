/**************************************************************************

Copyright:  WH team

Author: YinJichao <jichaoyinyjc@163.com>

Completion date:  XXX

Description: XXX

**************************************************************************/

#pragma once

#include "include/data2vtk.h"

namespace CAE
{
    // 读取Abaqus位移场
    void data_process::read_abaqus_dis(string path, vector<double> &abaqus_dis)
    {
        // read displacement file
        std::ifstream infile(path.c_str(), std::ios::in);
        string line;
        // record data
        string dis_x, dis_y, dis_z;
        int item_ele_id = 0;
        while (getline(infile, line))
        {
            std::istringstream iss_x(line);
            iss_x >> dis_x;
            abaqus_dis[3 * item_ele_id] = stod(dis_x);
            getline(infile, line);
            std::istringstream iss_y(line);
            iss_y >> dis_y;
            abaqus_dis[3 * item_ele_id + 1] = stod(dis_y);
            getline(infile, line);
            std::istringstream iss_z(line);
            iss_z >> dis_z;
            abaqus_dis[3 * item_ele_id + 2] = stod(dis_z);
            item_ele_id += 1;
        }
        cout << "the displacement of ABAQUS has been readed." << endl;
    }

    // 导出 位移场 VTK
    void data_process::export_dis_2_vtk(data_management &data_cae, string path, double scale_dis)
    {
        std::ofstream fout;
        fout.open(path, std::ios::out);
        if (!fout)
        {
            std::cerr << "Error: Cannot open " << path << std::endl;
            exit(EXIT_FAILURE);
        }
        fout << std::unitbuf; // 关闭缓冲
        fout << "# vtk DataFile Version 3.0\n";                 // Version Statement
        fout << "The density field of the optimized results\n"; // title
        fout << "ASCII\n";                                      // file format statement
        fout << "DATASET UNSTRUCTURED_GRID\n\n";                // data format: unstructured grid

        // 输入节点坐标
        int num_node = data_cae.nd_;
        fout << "POINTS\t" << num_node << "\tdouble\n";

        int id = 1;
        for (int i = 0; i < num_node; i++)
        {
            fout << data_cae.coords_[i][0] + scale_dis * data_cae.single_full_dis_vec_[3 * i] << "\t\t"
                 << data_cae.coords_[i][1] + scale_dis * data_cae.single_full_dis_vec_[3 * i + 1] << "\t\t"
                 << data_cae.coords_[i][2] + scale_dis * data_cae.single_full_dis_vec_[3 * i + 2] << "\n";
        }

        // 统计单元类型
        int num_ele = data_cae.ne_;
        int num_ele_C3D4 = 0;
        int num_ele_C3D8 = 0;
        int num_ele_C3D8R = 0;
        for (int i = 0; i < num_ele; i++)
        {
            int ele_type = data_cae.ele_list_idx_[i];
            int map_idx = data_cae.ele_map_list_[ele_type];
            string item_ele_type = data_cae.ele_list_[map_idx]->type_;
            switch (ELE_TYPES[item_ele_type])
            {
            case 1:
            {
                num_ele_C3D4 += 1;
                break;
            }
            case 2:
            {
                num_ele_C3D8 += 1;
                break;
            }
            case 3:
            {
                num_ele_C3D8R += 1;
                break;
            }
            default:
            {
                cout << "This type does not exist in the element library" << endl;
                break;
            }
            }
        }
        // 输入节点拓扑关系
        fout << "CELLS\t" << num_ele << "\t" << num_ele_C3D8 * (8 + 1) + num_ele_C3D8R * (8 + 1) + num_ele_C3D4 * (4 + 1) << "\n";
        for (int i = 0; i < num_ele; i++)
        {
            int ele_type = data_cae.ele_list_idx_[i];
            int map_idx = data_cae.ele_map_list_[ele_type];
            string item_ele_type = data_cae.ele_list_[map_idx]->type_;
            switch (ELE_TYPES[item_ele_type])
            {
            case 1:
            {
                fout << 4 << "\t" << data_cae.node_topos_[i][0] - 1 << "\t" << data_cae.node_topos_[i][2] - 1 << "\t"
                     << data_cae.node_topos_[i][1] - 1 << "\t" << data_cae.node_topos_[i][3] - 1 << "\n";

                break;
            }
            case 2:
            {
                fout << 8 << "\t" << data_cae.node_topos_[i][0] - 1 << "\t" << data_cae.node_topos_[i][1] - 1 << "\t"
                     << data_cae.node_topos_[i][3] - 1 << "\t" << data_cae.node_topos_[i][2] - 1 << "\t"
                     << data_cae.node_topos_[i][4] - 1 << "\t" << data_cae.node_topos_[i][5] - 1 << "\t"
                     << data_cae.node_topos_[i][7] - 1 << "\t" << data_cae.node_topos_[i][6] - 1 << "\n";

                break;
            }
            case 3:
            {
                fout << 8 << "\t" << data_cae.node_topos_[i][0] - 1 << "\t" << data_cae.node_topos_[i][1] - 1 << "\t"
                     << data_cae.node_topos_[i][3] - 1 << "\t" << data_cae.node_topos_[i][2] - 1 << "\t"
                     << data_cae.node_topos_[i][4] - 1 << "\t" << data_cae.node_topos_[i][5] - 1 << "\t"
                     << data_cae.node_topos_[i][7] - 1 << "\t" << data_cae.node_topos_[i][6] - 1 << "\n";

                break;
            }
            default:
            {
                cout << "This type does not exist in the element library" << endl;
                break;
            }
            }
        }
        // 写入单元类型
        fout << "CELL_TYPES\t\t" << num_ele << "\n";
        for (int i = 0; i < num_ele; i++)
        {
            int ele_type = data_cae.ele_list_idx_[i];
            int map_idx = data_cae.ele_map_list_[ele_type];
            string item_ele_type = data_cae.ele_list_[map_idx]->type_;
            switch (ELE_TYPES[item_ele_type])
            {
            case 1:
            {
                fout << 10 << "\n";

                break;
            }
            case 2:
            {
                fout << 11 << "\n";

                break;
            }
            case 3:
            {
                fout << 11 << "\n";

                break;
            }
            default:
            {
                cout << "This type does not exist in the element library" << endl;

                break;
            }
            }
        }

        // 写入单元位移
        // X
        fout << "\nPOINT_DATA\t" << num_node
             << "\nSCALARS u_x double 1\n"
             << "LOOKUP_TABLE  table1\n";
        for (int i = 0; i < num_node; i++)
        {
            fout << data_cae.single_full_dis_vec_[3 * i] << "\n";
        }
        // Y
        fout << "\nSCALARS u_y double 1\n"
             << "LOOKUP_TABLE  table2\n";
        for (int i = 0; i < num_node; i++)
        {
            fout << data_cae.single_full_dis_vec_[3 * i + 1] << "\n";
        }
        // Z
        fout << "\nSCALARS u_z double 1\n"
             << "LOOKUP_TABLE  table3\n";
        for (int i = 0; i < num_node; i++)
        {
            fout << data_cae.single_full_dis_vec_[3 * i + 2] << "\n";
        }
        // 合位移
        fout << "\nSCALARS u_magnitude double 1\n"
             << "LOOKUP_TABLE  table4\n";
        for (int i = 0; i < num_node; i++)
        {
            double u_ = sqrt(data_cae.single_full_dis_vec_[3 * i + 2] * data_cae.single_full_dis_vec_[3 * i + 2] +
                             data_cae.single_full_dis_vec_[3 * i + 1] * data_cae.single_full_dis_vec_[3 * i + 1] +
                             data_cae.single_full_dis_vec_[3 * i] * data_cae.single_full_dis_vec_[3 * i]);
            fout << u_ << "\n";
        }

        // 写入单元应力
        if (data_cae.stress_node_mat_.size() != 0)
        {
            // xx应力
            fout << "\nSCALARS S_xx double 1\n"
                 << "LOOKUP_TABLE  table5\n";
            for (int i = 0; i < num_node; i++)
            {
                fout << data_cae.stress_node_mat_[0][i] << "\n";
            }
            // yy应力
            fout << "\nSCALARS S_yy double 1\n"
                 << "LOOKUP_TABLE  table6\n";
            for (int i = 0; i < num_node; i++)
            {
                fout << data_cae.stress_node_mat_[1][i] << "\n";
            }
            // zz应力
            fout << "\nSCALARS S_zz double 1\n"
                 << "LOOKUP_TABLE  table7\n";
            for (int i = 0; i < num_node; i++)
            {
                fout << data_cae.stress_node_mat_[2][i] << "\n";
            }
            // xy应力
            fout << "\nSCALARS S_xy double 1\n"
                 << "LOOKUP_TABLE  table8\n";
            for (int i = 0; i < num_node; i++)
            {
                fout << data_cae.stress_node_mat_[3][i] << "\n";
            }
            // yz应力
            fout << "\nSCALARS S_yz double 1\n"
                 << "LOOKUP_TABLE  table9\n";
            for (int i = 0; i < num_node; i++)
            {
                fout << data_cae.stress_node_mat_[4][i] << "\n";
            }
            // xz应力
            fout << "\nSCALARS S_xz double 1\n"
                 << "LOOKUP_TABLE  table10\n";
            for (int i = 0; i < num_node; i++)
            {
                fout << data_cae.stress_node_mat_[5][i] << "\n";
            }
            // Von_Mises应力
            fout << "\nSCALARS S_Mises double 1\n"
                 << "LOOKUP_TABLE  table11\n";
            for (int i = 0; i < num_node; i++)
            {
                fout << data_cae.stress_node_mat_[6][i] << "\n";
            }
        }
        fout.close();
    };

    void CAE::data_process::CA_export_dis_2_vtk(data_management &data_cae, double scale_dis, string result_path)
    {
        std::ofstream fout;
        fout.open(result_path, std::ios::out);
        if (!fout)
        {
            std::cerr << "Error: Cannot open " << result_path << std::endl;
            exit(EXIT_FAILURE);
        }
        fout << std::unitbuf; // 关闭缓冲
        fout << "# vtk DataFile Version 3.0\n";                 // Version Statement
        fout << "The density field of the optimized results\n"; // title
        fout << "ASCII\n";                                      // file format statement
        fout << "DATASET UNSTRUCTURED_GRID\n\n";                // data format: unstructured grid

        // 输入节点坐标
        int num_node = data_cae.coords_mfull_.size();
        fout << "POINTS\t" << num_node << "\tdouble\n";
        int id = 1;
        for (int i = 0; i < num_node; i++)
        {
            double dis1 = data_cae.coords_mfull_[i][0] + scale_dis * data_cae.single_full_ca_dis_vec_[3 * i];
            double dis2 = data_cae.coords_mfull_[i][1] + scale_dis * data_cae.single_full_ca_dis_vec_[3 * i + 1];
            double dis3 = data_cae.coords_mfull_[i][2] + scale_dis * data_cae.single_full_ca_dis_vec_[3 * i + 2];
            fout << dis1 << "\t\t" << dis2 << "\t\t" << dis3 << "\n";
        }

        // 统计单元类型
        int num_ele = data_cae.node_topos_mfull_.size();
        int num_ele_C3D4 = 0;
        int num_ele_C3D8 = 0;
        int num_ele_C3D8R = 0;
        for (int i = 0; i < num_ele; i++)
        {
            int ele_type = data_cae.ele_list_idx_m_[i];
            int map_idx = data_cae.ele_map_list_[ele_type];
            string item_ele_type = data_cae.ele_list_[map_idx]->type_;
            switch (ELE_TYPES[item_ele_type])
            {
            case 1:
            {
                num_ele_C3D4 += 1;
                break;
            }
            case 2:
            {
                num_ele_C3D8 += 1;
                break;
            }
            case 3:
            {
                num_ele_C3D8R += 1;
                break;
            }
            default:
            {
                cout << "This type does not exist in the element library" << endl;
                break;
            }
            }
        }
        // 输入节点拓扑关系
        fout << "CELLS\t" << num_ele << "\t" << num_ele_C3D8 * (8 + 1) + num_ele_C3D8R * (8 + 1) + num_ele_C3D4 * (4 + 1) << "\n";
        for (int i = 0; i < num_ele; i++)
        {
            int ele_type = data_cae.ele_list_idx_m_[i];
            int map_idx = data_cae.ele_map_list_[ele_type];
            string item_ele_type = data_cae.ele_list_[map_idx]->type_;
            switch (ELE_TYPES[item_ele_type])
            {
            case 1:
            {
                fout << 4 << "\t" << data_cae.node_topos_mfull_[i][0] - 1 << "\t" << data_cae.node_topos_mfull_[i][2] - 1 << "\t"
                     << data_cae.node_topos_mfull_[i][1] - 1 << "\t" << data_cae.node_topos_mfull_[i][3] - 1 << "\n";
                break;
            }
            case 2:
            {
                fout << 8 << "\t" << data_cae.node_topos_mfull_[i][0] - 1 << "\t" << data_cae.node_topos_mfull_[i][1] - 1 << "\t"
                     << data_cae.node_topos_mfull_[i][3] - 1 << "\t" << data_cae.node_topos_mfull_[i][2] - 1 << "\t"
                     << data_cae.node_topos_mfull_[i][4] - 1 << "\t" << data_cae.node_topos_mfull_[i][5] - 1 << "\t"
                     << data_cae.node_topos_mfull_[i][7] - 1 << "\t" << data_cae.node_topos_mfull_[i][6] - 1 << "\n";
                break;
            }
            case 3:
            {
                fout << 8 << "\t" << data_cae.node_topos_mfull_[i][0] - 1 << "\t" << data_cae.node_topos_mfull_[i][1] - 1 << "\t"
                     << data_cae.node_topos_mfull_[i][3] - 1 << "\t" << data_cae.node_topos_mfull_[i][2] - 1 << "\t"
                     << data_cae.node_topos_mfull_[i][4] - 1 << "\t" << data_cae.node_topos_mfull_[i][5] - 1 << "\t"
                     << data_cae.node_topos_mfull_[i][7] - 1 << "\t" << data_cae.node_topos_mfull_[i][6] - 1 << "\n";
                break;
            }
            default:
            {
                cout << "This type does not exist in the element library" << endl;
                break;
            }
            }
        }
        // 写入单元类型
        fout << "CELL_TYPES\t\t" << num_ele << "\n";
        for (int i = 0; i < num_ele; i++)
        {
            int ele_type = data_cae.ele_list_idx_m_[i];
            int map_idx = data_cae.ele_map_list_[ele_type];
            string item_ele_type = data_cae.ele_list_[map_idx]->type_;
            switch (ELE_TYPES[item_ele_type])
            {
            case 1:
            {
                fout << 10 << "\n";
                break;
            }
            case 2:
            {
                fout << 11 << "\n";
                break;
            }
            case 3:
            {
                fout << 11 << "\n";
                break;
            }
            default:
            {
                cout << "This type does not exist in the element library" << endl;
                break;
            }
            }
        }

        // 写入单元位移
        // X
        fout << "\nPOINT_DATA\t" << num_node
             << "\nSCALARS u_x double 1\n"
             << "LOOKUP_TABLE  table1\n";
        for (int i = 0; i < num_node; i++)
        {
            fout << data_cae.single_full_ca_dis_vec_[3 * i] << "\n";
        }
        // Y
        fout << "\nSCALARS u_y double 1\n"
             << "LOOKUP_TABLE  table2\n";
        for (int i = 0; i < num_node; i++)
        {
            fout << data_cae.single_full_ca_dis_vec_[3 * i + 1] << "\n";
        }
        // Z
        fout << "\nSCALARS u_z double 1\n"
             << "LOOKUP_TABLE  table3\n";
        for (int i = 0; i < num_node; i++)
        {
            fout << data_cae.single_full_ca_dis_vec_[3 * i + 2] << "\n";
        }
        // 合位移
        fout << "\nSCALARS u_magnitude double 1\n"
             << "LOOKUP_TABLE  table4\n";
        for (int i = 0; i < num_node; i++)
        {
            double u_ = sqrt(data_cae.single_full_ca_dis_vec_[3 * i + 2] * data_cae.single_full_ca_dis_vec_[3 * i + 2] +
                             data_cae.single_full_ca_dis_vec_[3 * i + 1] * data_cae.single_full_ca_dis_vec_[3 * i + 1] +
                             data_cae.single_full_ca_dis_vec_[3 * i] * data_cae.single_full_ca_dis_vec_[3 * i]);
            fout << u_ << "\n";
        }
        fout.close();
    }
}
