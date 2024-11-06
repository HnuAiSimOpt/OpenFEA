/**************************************************************************

Copyright:  WH team

Author: YinJichao <jichaoyinyjc@163.com>

Completion date:  XXX

Description: XXX

**************************************************************************/

#pragma once
#include <iostream>
#include "include/cae.h"
#include "include/elastic_mat.h"
#include "include/sample_eigen_superlu_mkl_svd.h"

void code_test();

int main(int argc, char* argv[])
{
    bool test_ = false;//设置true为以前的启动模式
    if (test_) {
        code_test();
        return 0;
    }
    else {
        // 建立CAE分析对象
        CAE::CAE_process cae_item;
        //初始化
        cae_item.Init(argc, argv);
        //求解
        cae_item.Solve();
        return 0;
    }
}

void code_test()
{
    int case_num = 205;
    if (case_num == 1)
    {
        // 材料属性赋值
        CAE::elastic_mat mat_item{ 2.1e5, 0.3, 7800 };
        // 材料路径
        std::string path = "E:\\WH_CAE\\test_model\\cadcae-linear.inp";
        // 关键字
        string load_set_keyword = "Set-load";
        string load_value_keyword = "Cload";
        string dis_set_keyword = "Set-fix";
        // 建立CAE分析对象
        CAE::CAE_process cae_item(path, mat_item);
        // 读取计算文件
        cae_item.pre_info(load_set_keyword, load_value_keyword, dis_set_keyword);
        // 执行结构响应分析
        string result_path = "E:\\WH_CAE\\test_model\\output\\cadcae-linear.vtk";
        cae_item.implict_analysis(result_path);
    }
    else if (case_num == 2)
    {
        // 非协调计算
        // 材料属性赋值
        CAE::elastic_mat mat_item{ 2.1e5, 0.3, 7800 }; // 此材料参数用于非协调计算
        // 材料路径
        string path = "E:/project/CADCAE/bracket_nonconforming/inp/Job-c3d4-ele12k.inp"; // 非协调路径
        // 关键字
        string load_set_keyword = "Set-load";
        string load_value_keyword = "Cload";
        string dis_set_keyword = "Set-fix";
        // 建立CAE分析对象
        CAE::CAE_process cae_item(path, mat_item);
        // 读取计算文件
        cae_item.pre_info(load_set_keyword, load_value_keyword, dis_set_keyword);
        // 执行结构响应分析
        string result_path = "E:/project/CADCAE/bracket_nonconforming/output/dis_with_error.vtk";
        cae_item.implict_analysis(result_path);
    }
    else if (case_num == 101) // 隐式完整分析+重分析
    {
        // 开始隐式分析
        // 材料属性赋值
        // CAE::elastic_mat mat_item{2.1e5, 0.3, 7800};
        // 材料路径
        // std::string path = "E:\\WH_CAE\\test_model\\CA_cadcae\\mesh_o.inp";
        // 关键字
        // string load_set_keyword = "Set-load";
        // string load_value_keyword = "Cload";
        // string dis_set_keyword = "Set-fix";
        // 建立CAE分析对象
        // CAE::CAE_process cae_item(path, mat_item);
        // 读取计算文件
        // cae_item.pre_info(load_set_keyword, load_value_keyword, dis_set_keyword);
        // 执行结构响应分析
        // string result_path = "E:\\WH_CAE\\test_model\\output\\original.vtk";
        // string path_abaqus = "E:\\WH_CAE\\test_model\\original.txt";
        // cae_item.implict_analysis(result_path, path_abaqus);

        // 开始重分析
        // 修改网格路径(仅包含1个part的网格)
        // string mesh_path = "E:\\WH_CAE\\test_model\\original.inp"; // 此处使用original模型，里面包含删除单元
        // 关键字
        // string CA_del_set_keyword = "Set-del";
        // 建立重分析CAE分析对象
        // cae_item.path_ = mesh_path;
        // 处理重复节点网格
        // vector<int> del_topo;
        // cae_item.CA_pre_process(CA_del_set_keyword, del_topo);
        // 执行重分析
        // string CA_result_path = "E:\\WH_CAE\\test_model\\output\\CA_modefy.vtk";
        // string CA_path_abaqus = "E:\\WH_CAE\\test_model\\CA_modefy.txt";
        // bool Is_Update = false; // 是否将原模型更新为修改后的模型
        // cae_item.CA_ReAnalysis(CA_result_path, CA_path_abaqus, del_topo, Is_Update);
    }
    else if (case_num == 102) // 隐式完整分析+重分析
    {
        // ----------------------------------------------------------------------------------------------------------- 
        // 隐式 全分析
        // ----------------------------------------------------------------------------------------------------------- 
        // 材料属性赋值
        CAE::elastic_mat mat_item{ 2.1e5, 0.3, 7800 };
        // 材料路径
        std::string path = "E:\\WH_CAE\\test_model\\CA_cadcae\\mesh_o.inp";
        // 关键字
        string load_set_keyword = "Set-load";
        string load_value_keyword = "Cload";
        string dis_set_keyword = "Set-fix";
        // 建立CAE分析对象
        CAE::CAE_process cae_item(path, mat_item);
        // 读取计算文件
        cae_item.pre_info(load_set_keyword, load_value_keyword, dis_set_keyword);
        // 执行结构响应分析
        string result_path = "E:\\WH_CAE\\test_model\\output\\cadcae_o.vtk";
        bool is_save_stiffness = true;
        cae_item.implict_analysis(result_path, is_save_stiffness);

        // ----------------------------------------------------------------------------------------------------------- 
        // 重分析
        // ----------------------------------------------------------------------------------------------------------- 
        // 修改网格路径(仅包含1个part的网格)
        string mesh_path1 = "E:\\WH_CAE\\test_model\\CA_cadcae\\CA_info.txt";
        string mesh_path2 = "E:\\WH_CAE\\test_model\\CA_cadcae\\mesh_m.inp";
        cae_item.CA_pre_process(mesh_path1, mesh_path2);
        // 开始执行重分析
        string CA_result_path = "E:\\WH_CAE\\test_model\\output\\cadcae_ca.vtk";
        int n_basis = 4;
        cae_item.CA_ReAnalysis(CA_result_path, n_basis);
    }
    else if (case_num == 201)
    {
        // ----------------------------------------------------------------------------------------------------------- 
        // 显式动力学分析
        // ----------------------------------------------------------------------------------------------------------- 
            // 材料属性赋值
            CAE::elastic_mat mat_item{ 2.1e5, 0.3, 7800 };
            // 材料路径
            std::string path = "E:\\WH_CAE\\test_model\\CA_cadcae\\mesh_o.inp";
            // 关键字
            string load_set_keyword = "Set-load";
            string load_value_keyword = "Cload";
            string dis_set_keyword = "Set-fix";
            // 建立CAE分析对象
            CAE::CAE_process cae_item(path, mat_item);
            // 读取计算文件
            cae_item.pre_info(load_set_keyword, load_value_keyword, dis_set_keyword);
            // 执行动力学分析
             string result_path = "F:/OpenFEM/model/";
             string path_abaqus = "F:/OpenFEM/model/Abaqus_U.txt";
             // 可以手动指定计算步长，会以此处修改生效
             // cae_item.data_cae_.time_total_ = 0.1;       
             // cae_item.data_cae_.time_step_ = 0.0;
             cae_item.explicit_analysis(result_path, path_abaqus);
    }
    else if (case_num == 202)
    {
        // ----------------------------------------------------------------------------------------------------------- 
        // 隐式 全分析
        // ----------------------------------------------------------------------------------------------------------- 
        // 材料属性赋值
        CAE::elastic_mat mat_item{ 2.1e5, 0.3, 7800 };
        // 材料路径
        std::string path = "F:\\OpenFEM\\model\\CAppt\\full_analysis_fine.inp";
        // 关键字
        string load_set_keyword = "Set-load";
        string load_value_keyword = "Cload";
        string dis_set_keyword = "Set-fix";
        // 建立CAE分析对象
        CAE::CAE_process cae_item(path, mat_item);
        // 读取计算文件
        cae_item.pre_info(load_set_keyword, load_value_keyword, dis_set_keyword);
        // 执行结构响应分析
        string result_path = "F:\\OpenFEM\\model\\CAppt\\full_analysis_fine.vtk";
        bool is_save_stiffness = true;
        cae_item.implict_analysis(result_path, is_save_stiffness);

        // ----------------------------------------------------------------------------------------------------------- 
        // 重分析
        // ----------------------------------------------------------------------------------------------------------- 
        // 修改网格路径(仅包含1个part的网格)
        string mesh_path1 = "F:\\OpenFEM\\model\\CAppt\\fine_info_ca.txt";
        string mesh_path2 = "F:\\OpenFEM\\model\\CAppt\\fine_model_m.inp";
        cae_item.CA_pre_process(mesh_path1, mesh_path2);
        // 开始执行重分析
        string CA_result_path = "F:\\OpenFEM\\model\\CAppt\\fine_model_m.vtk";
        int n_basis = 4;
        cae_item.CA_ReAnalysis(CA_result_path, n_basis);
    }
    else if (case_num == 203)
    {
        // ----------------------------------------------------------------------------------------------------------- 
        // 隐式 全分析
        // ----------------------------------------------------------------------------------------------------------- 
        // 材料属性赋值
        CAE::elastic_mat mat_item{ 2.1e5, 0.3, 7800 };
        // 材料路径
        std::string path = "F:\\OpenFEM\\model\\bugfix-stress\\solid3_o.inp";
        // 关键字
        string load_set_keyword = "Set-load";
        string load_value_keyword = "Cload";
        string dis_set_keyword = "Set-fix";
        // 建立CAE分析对象
        CAE::CAE_process cae_item(path, mat_item);
        // 读取计算文件
        cae_item.pre_info(load_set_keyword, load_value_keyword, dis_set_keyword);
        // 执行结构响应分析
        string result_path = "F:\\OpenFEM\\model\\bugfix-stress\\solid3_o.vtk";
        bool is_save_stiffness = true;
        cae_item.implict_analysis(result_path, is_save_stiffness);
    }
    else if (case_num == 204)
    {
        // ----------------------------------------------------------------------------------------------------------- 
        // 非协调计算
        // ----------------------------------------------------------------------------------------------------------- 
        // 材料属性赋值
        CAE::elastic_mat mat_item{ 2.1e5, 0.3, 7800 }; // 此材料参数用于非协调计算
        // 材料路径
        // std::string path = "E:/project/CADCAE/bracket_nonconforming/inp/model.inp";//c3d8非协调路径
        // std::string path = "E:/project/CADCAE/bracket_nonconforming/inp/Job-c3d4-random.inp";//非协调路径
        // std::string path = "E:/project/CADCAE/bracket_nonconforming/inp/Job-c3d4-ncf.inp";//非协调路径
        string path = "F:\\OpenFEM\\model\\noncomfortable\\Job-c3d4-ele12k.inp"; // 非协调路径
        // std::string path = "E:/project/CADCAE/bracket_nonconforming/inp/Job-c3d4_110k.inp";//非协调路径
        // std::string path = "E:/project/CADCAE/bracket_nonconforming/inp/Job-c3d4-1000k.inp";//非协调路径
        //  关键字
        string load_set_keyword = "Set-load";
        string load_value_keyword = "Cload";
        string dis_set_keyword = "Set-fix";
        // 建立CAE分析对象
        CAE::CAE_process cae_item(path, mat_item);
        // 读取计算文件
        cae_item.pre_info(load_set_keyword, load_value_keyword, dis_set_keyword);
        // 执行结构响应分析
        string result_path = "F:\\OpenFEM\\model\\noncomfortable\\dis_with_error.vtk";
        cae_item.implict_analysis(result_path);
    }
    else if (case_num == 205)
    {
        // ----------------------------------------------------------------------------------------------------------- 
        // SFEM
        // ----------------------------------------------------------------------------------------------------------- 
        // 材料属性赋值
        CAE::elastic_mat mat_item{ 2.1e5, 0.3, 7800 };
        // 材料路径
        std::string path = "F:\\OpenFEM\\model\\CAppt\\full_analysis_fine.inp";
        // 关键字
        string load_set_keyword = "Set-load";
        string load_value_keyword = "Cload";
        string dis_set_keyword = "Set-fix";
        // 建立CAE分析对象 
        CAE::CAE_process cae_item(path, mat_item);
        // 读取计算文件
        cae_item.pre_info(load_set_keyword, load_value_keyword, dis_set_keyword);
        // 执行结构响应分析
        string result_path = "F:\\OpenFEM\\model\\CAppt\\full_analysis_fine_SFEM. vtk";
        cae_item.implict_SFEManalysis(result_path, "");
    }
    else
    {
        std::cout << "Please check your code\n";
    }
}










//--------------------------------------------------------------------------------------------//
///**************************************************************************
//
//Copyright:  WH team
//
//Author: YinJichao <jichaoyinyjc@163.com>
//
//Completion date:  XXX
//
//Description: XXX
//
//**************************************************************************/
//
//#pragma once
//#include <iostream>
//#include "include/cae.h"
//#include "include/elastic_mat.h"
//#include "include/sample_eigen_superlu_mkl_svd.h"
//
//const char* const prog_name = "OpenFEA.exe";
//const char* const prog_version = "version 0.1 (2024/06/25)";
//
//void code_test();
//void implicit_(string ifile_name, string ofile_name);
//void sfem_(string ifile_name, string ofile_name);
//void explicit_(string ifile_name, string ofile_name);
//
//inline void usage(int exit_value = 0) {
//    cerr << prog_name << " (-i) input_filename (-o output_filename)" << endl;
//    exit(exit_value);
//}
//
//int main(int argc, char* argv[])
//{
//    bool test_ = false;
//    if (test_) {
//        code_test();
//        return 0;
//    }
//    else {
//        // 建立CAE分析对象
//        CAE::CAE_process cae_item;
//
//        //初始化
//        cae_item.Init(argc, argv);
//
//        //求解
//        cae_item.Solve();
//
//
//
//
//        //声明用于记录用户指定选项的变量
//        bool ifile_on = false;
//        bool ofile_on = false;
//
//        string ifile_name;//记录出现的输入文件名
//        string ofile_name;//记录出现的输出文件名
//
//        // 记录求解类型
//        bool impl = true;// 隐式；默认使用隐式求解器
//        bool sfem = false;// sfem
//        bool exp = false;// 显式
//
//        //cout << "argc:" << argc << endl;
//        char* pchar;
//        for (int i = 1; i < argc; ++i) {//读取argv中的每个选项
//            //输出第i+1个参量,argv[0]是程序本身
//            cout << "argv[" << i << "]:" << argv[i] << endl;
//
//            pchar = argv[i];
//            if (strcmp(pchar, "-imp") == 0) {
//                option.analysis_type = 1;
//            }
//            else if (strcmp(pchar, "-exp") == 0) {
//                option.analysis_type = 2;
//            }
//            else if (strcmp(pchar, "-sfem") == 0) {
//                option.sfem_flag = true;
//            }
//            else if (strcmp(pchar, "-re") == 0) {
//                option.reanalysis_flag = true;
//            }
//
//
//            switch (pchar[0]) {//确定选项类型：-i,-o,-h,-v;或者其他
//            case '-': {
//                //cout << "case \'-\' found" << endl;
//                switch (pchar[1]) {//确定用户指定的选项：i,o,h,v
//                case 'i'://处理输入文件：
//                    //cout << "-i found:input file turned on!" << endl;
//                    ifile_on = true;
//                    break;
//                case 'o'://处理输出文件
//                    //cout << "-o found:output file!" << endl;
//                    ofile_on = true;
//                    break;
//                case 's'://处理输出文件
//                    //cout << "sfem" << endl;
//                    impl = false;
//                    sfem = true;
//                    break;
//                case 'e'://处理输出文件
//                    //cout << "explicit" << endl;
//                    impl = false;
//                    exp = true;
//                    break;
//                case 'h'://处理帮助
//                    //cout << "-h found:help info!" << endl;
//                    usage();
//                case 'v'://处理版本请求
//                    //cout << "-v found:version info displayed!" << endl;
//                    cout << prog_name << ":" << prog_version << endl;
//                    return 0;
//                default://无法识别的选项
//                    cerr << prog_name << ":error:unrecognition option -:" << pchar << endl;
//                    usage(-1);
//                }
//                break;
//            }
//            default://不以'-'开头，是文件名
//                if (ifile_on) {//输出文件名
//                    cout << "input file path:" << pchar << endl;
//                    ifile_name = pchar;
//                    ifile_on = false;//复位
//                }
//                else if (ofile_on) {//输出文件名
//                    cout << "set output file path:" << pchar << endl;
//                    ofile_name = pchar;
//                    ofile_on = false;//复位
//                }
//                else {//文件名
//                    cout << "input file path:" << pchar << endl;
//                    ifile_name = pchar;
//                }
//                break;
//            }
//        }
//
//        if (ifile_name.empty()) {
//            cout << "Please enter the correct path to the input file ending in .inp" << ofile_name << endl;
//            usage(-1);
//        }
//        else {
//            // 找到最后一个'/'或'\'的位置
//            size_t last_slash_pos = ifile_name.find_last_of("/\\");
//            // 找到最后一个'.'的位置
//            size_t last_dot_pos = ifile_name.find_last_of('.');
//            // 如果没有找到'/'或'\'，则文件名就是整个路径;如果没有找到'.'，则没有后缀名;或者后缀不是inp
//            if (last_slash_pos == std::string::npos || last_dot_pos == std::string::npos || ifile_name.substr(last_dot_pos + 1) != "inp") {
//                cout << "Please enter the correct path to the input file ending in .inp" << ofile_name << endl;
//                usage(-1);
//            }
//            if (!ofile_name.empty()) {
//                size_t last_dot_pos_out = ofile_name.find_last_of('.');
//                // 如果没有找到'.'，则没有后缀名;或者后缀不是vtk
//                if (last_dot_pos == std::string::npos || ofile_name.substr(last_dot_pos_out + 1) != "vtk") {
//                    cout << "Please enter the correct path to the output file ending in .vtk" << ofile_name << endl;
//                    usage(-1);
//                }
//            }
//            else {
//                ofile_name = ifile_name.substr(0, last_dot_pos) + ".vtk";
//            }
//        }
//
//        if (impl) {
//            implicit_(ifile_name, ofile_name);
//        }
//        else if (sfem) {
//            sfem_(ifile_name, ofile_name);
//        }
//        else if (exp) {
//            explicit_(ifile_name, ofile_name);
//        }
//
//        return 0;
//    }
//
//}
//
//void implicit_(string ifile_name, string ofile_name) {
//    // 建立CAE分析对象
//    CAE::CAE_process cae_item(ifile_name);
//    // 关键字
//    string load_set_keyword = "Set-load";
//    string load_value_keyword = "Cload";
//    string dis_set_keyword = "Set-fix";
//    // 读取计算文件
//    cae_item.pre_info(load_set_keyword, load_value_keyword, dis_set_keyword);
//    bool is_save_stiffness = true;
//    cae_item.implict_analysis(ofile_name, is_save_stiffness);
//}
//
//void sfem_(string ifile_name, string ofile_name) {
//    // ----------------------------------------------------------------------------------------------------------- 
//    // SFEM
//    // ----------------------------------------------------------------------------------------------------------- 
//    // 建立CAE分析对象
//    CAE::CAE_process cae_item(ifile_name);
//    // 关键字
//    string load_set_keyword = "Set-load";
//    string load_value_keyword = "Cload";
//    string dis_set_keyword = "Set-fix";
//    // 读取计算文件
//    cae_item.pre_info(load_set_keyword, load_value_keyword, dis_set_keyword);
//    // 执行结构响应分析
//    cae_item.implict_SFEManalysis(ofile_name, "");
//}
//
//void explicit_(string ifile_name, string ofile_name) {
//
//}
//
//void code_test()
//{
//    int case_num = 205;
//    if (case_num == 1)
//    {
//        // 材料属性赋值
//        CAE::elastic_mat mat_item{ 2.1e5, 0.3, 7800 };
//        // 材料路径
//        std::string path = "E:\\WH_CAE\\test_model\\cadcae-linear.inp";
//        // 关键字
//        string load_set_keyword = "Set-load";
//        string load_value_keyword = "Cload";
//        string dis_set_keyword = "Set-fix";
//        // 建立CAE分析对象
//        CAE::CAE_process cae_item(path, mat_item);
//        // 读取计算文件
//        cae_item.pre_info(load_set_keyword, load_value_keyword, dis_set_keyword);
//        // 执行结构响应分析
//        string result_path = "E:\\WH_CAE\\test_model\\output\\cadcae-linear.vtk";
//        cae_item.implict_analysis(result_path);
//    }
//    else if (case_num == 2)
//    {
//        // 非协调计算
//        // 材料属性赋值
//        CAE::elastic_mat mat_item{ 2.1e5, 0.3, 7800 }; // 此材料参数用于非协调计算
//        // 材料路径
//        string path = "E:/project/CADCAE/bracket_nonconforming/inp/Job-c3d4-ele12k.inp"; // 非协调路径
//        // 关键字
//        string load_set_keyword = "Set-load";
//        string load_value_keyword = "Cload";
//        string dis_set_keyword = "Set-fix";
//        // 建立CAE分析对象
//        CAE::CAE_process cae_item(path, mat_item);
//        // 读取计算文件
//        cae_item.pre_info(load_set_keyword, load_value_keyword, dis_set_keyword);
//        // 执行结构响应分析
//        string result_path = "E:/project/CADCAE/bracket_nonconforming/output/dis_with_error.vtk";
//        cae_item.implict_analysis(result_path);
//    }
//    else if (case_num == 101) // 隐式完整分析+重分析
//    {
//        // 开始隐式分析
//        // 材料属性赋值
//        // CAE::elastic_mat mat_item{2.1e5, 0.3, 7800};
//        // 材料路径
//        // std::string path = "E:\\WH_CAE\\test_model\\CA_cadcae\\mesh_o.inp";
//        // 关键字
//        // string load_set_keyword = "Set-load";
//        // string load_value_keyword = "Cload";
//        // string dis_set_keyword = "Set-fix";
//        // 建立CAE分析对象
//        // CAE::CAE_process cae_item(path, mat_item);
//        // 读取计算文件
//        // cae_item.pre_info(load_set_keyword, load_value_keyword, dis_set_keyword);
//        // 执行结构响应分析
//        // string result_path = "E:\\WH_CAE\\test_model\\output\\original.vtk";
//        // string path_abaqus = "E:\\WH_CAE\\test_model\\original.txt";
//        // cae_item.implict_analysis(result_path, path_abaqus);
//
//        // 开始重分析
//        // 修改网格路径(仅包含1个part的网格)
//        // string mesh_path = "E:\\WH_CAE\\test_model\\original.inp"; // 此处使用original模型，里面包含删除单元
//        // 关键字
//        // string CA_del_set_keyword = "Set-del";
//        // 建立重分析CAE分析对象
//        // cae_item.path_ = mesh_path;
//        // 处理重复节点网格
//        // vector<int> del_topo;
//        // cae_item.CA_pre_process(CA_del_set_keyword, del_topo);
//        // 执行重分析
//        // string CA_result_path = "E:\\WH_CAE\\test_model\\output\\CA_modefy.vtk";
//        // string CA_path_abaqus = "E:\\WH_CAE\\test_model\\CA_modefy.txt";
//        // bool Is_Update = false; // 是否将原模型更新为修改后的模型
//        // cae_item.CA_ReAnalysis(CA_result_path, CA_path_abaqus, del_topo, Is_Update);
//    }
//    else if (case_num == 102) // 隐式完整分析+重分析
//    {
//        // ----------------------------------------------------------------------------------------------------------- 
//        // 隐式 全分析
//        // ----------------------------------------------------------------------------------------------------------- 
//        // 材料属性赋值
//        CAE::elastic_mat mat_item{ 2.1e5, 0.3, 7800 };
//        // 材料路径
//        std::string path = "E:\\WH_CAE\\test_model\\CA_cadcae\\mesh_o.inp";
//        // 关键字
//        string load_set_keyword = "Set-load";
//        string load_value_keyword = "Cload";
//        string dis_set_keyword = "Set-fix";
//        // 建立CAE分析对象
//        CAE::CAE_process cae_item(path, mat_item);
//        // 读取计算文件
//        cae_item.pre_info(load_set_keyword, load_value_keyword, dis_set_keyword);
//        // 执行结构响应分析
//        string result_path = "E:\\WH_CAE\\test_model\\output\\cadcae_o.vtk";
//        bool is_save_stiffness = true;
//        cae_item.implict_analysis(result_path, is_save_stiffness);
//
//        // ----------------------------------------------------------------------------------------------------------- 
//        // 重分析
//        // ----------------------------------------------------------------------------------------------------------- 
//        // 修改网格路径(仅包含1个part的网格)
//        string mesh_path1 = "E:\\WH_CAE\\test_model\\CA_cadcae\\CA_info.txt";
//        string mesh_path2 = "E:\\WH_CAE\\test_model\\CA_cadcae\\mesh_m.inp";
//        cae_item.CA_pre_process(mesh_path1, mesh_path2);
//        // 开始执行重分析
//        string CA_result_path = "E:\\WH_CAE\\test_model\\output\\cadcae_ca.vtk";
//        int n_basis = 4;
//        cae_item.CA_ReAnalysis(CA_result_path, n_basis);
//    }
//    else if (case_num == 201)
//    {
//        // ----------------------------------------------------------------------------------------------------------- 
//        // 显式动力学分析
//        // ----------------------------------------------------------------------------------------------------------- 
//
//    }
//    else if (case_num == 202)
//    {
//        // ----------------------------------------------------------------------------------------------------------- 
//        // 隐式 全分析
//        // ----------------------------------------------------------------------------------------------------------- 
//        // 材料属性赋值
//        CAE::elastic_mat mat_item{ 2.1e5, 0.3, 7800 };
//        // 材料路径
//        std::string path = "F:\\OpenFEM\\model\\CAppt\\full_analysis_fine.inp";
//        // 关键字
//        string load_set_keyword = "Set-load";
//        string load_value_keyword = "Cload";
//        string dis_set_keyword = "Set-fix";
//        // 建立CAE分析对象
//        CAE::CAE_process cae_item(path, mat_item);
//        // 读取计算文件
//        cae_item.pre_info(load_set_keyword, load_value_keyword, dis_set_keyword);
//        // 执行结构响应分析
//        string result_path = "F:\\OpenFEM\\model\\CAppt\\full_analysis_fine.vtk";
//        bool is_save_stiffness = true;
//        cae_item.implict_analysis(result_path, is_save_stiffness);
//
//        // ----------------------------------------------------------------------------------------------------------- 
//        // 重分析
//        // ----------------------------------------------------------------------------------------------------------- 
//        // 修改网格路径(仅包含1个part的网格)
//        string mesh_path1 = "F:\\OpenFEM\\model\\CAppt\\fine_info_ca.txt";
//        string mesh_path2 = "F:\\OpenFEM\\model\\CAppt\\fine_model_m.inp";
//        cae_item.CA_pre_process(mesh_path1, mesh_path2);
//        // 开始执行重分析
//        string CA_result_path = "F:\\OpenFEM\\model\\CAppt\\fine_model_m.vtk";
//        int n_basis = 4;
//        cae_item.CA_ReAnalysis(CA_result_path, n_basis);
//    }
//    else if (case_num == 203)
//    {
//        // ----------------------------------------------------------------------------------------------------------- 
//        // 隐式 全分析
//        // ----------------------------------------------------------------------------------------------------------- 
//        // 材料属性赋值
//        CAE::elastic_mat mat_item{ 2.1e5, 0.3, 7800 };
//        // 材料路径
//        std::string path = "F:\\OpenFEM\\model\\bugfix-stress\\solid3_o.inp";
//        // 关键字
//        string load_set_keyword = "Set-load";
//        string load_value_keyword = "Cload";
//        string dis_set_keyword = "Set-fix";
//        // 建立CAE分析对象
//        CAE::CAE_process cae_item(path, mat_item);
//        // 读取计算文件
//        cae_item.pre_info(load_set_keyword, load_value_keyword, dis_set_keyword);
//        // 执行结构响应分析
//        string result_path = "F:\\OpenFEM\\model\\bugfix-stress\\solid3_o.vtk";
//        bool is_save_stiffness = true;
//        cae_item.implict_analysis(result_path, is_save_stiffness);
//    }
//    else if (case_num == 204)
//    {
//        // ----------------------------------------------------------------------------------------------------------- 
//        // 非协调计算
//        // ----------------------------------------------------------------------------------------------------------- 
//        // 材料属性赋值
//        CAE::elastic_mat mat_item{ 2.1e5, 0.3, 7800 }; // 此材料参数用于非协调计算
//        // 材料路径
//        // std::string path = "E:/project/CADCAE/bracket_nonconforming/inp/model.inp";//c3d8非协调路径
//        // std::string path = "E:/project/CADCAE/bracket_nonconforming/inp/Job-c3d4-random.inp";//非协调路径
//        // std::string path = "E:/project/CADCAE/bracket_nonconforming/inp/Job-c3d4-ncf.inp";//非协调路径
//        string path = "F:\\OpenFEM\\model\\noncomfortable\\Job-c3d4-ele12k.inp"; // 非协调路径
//        // std::string path = "E:/project/CADCAE/bracket_nonconforming/inp/Job-c3d4_110k.inp";//非协调路径
//        // std::string path = "E:/project/CADCAE/bracket_nonconforming/inp/Job-c3d4-1000k.inp";//非协调路径
//        //  关键字
//        string load_set_keyword = "Set-load";
//        string load_value_keyword = "Cload";
//        string dis_set_keyword = "Set-fix";
//        // 建立CAE分析对象
//        CAE::CAE_process cae_item(path, mat_item);
//        // 读取计算文件
//        cae_item.pre_info(load_set_keyword, load_value_keyword, dis_set_keyword);
//        // 执行结构响应分析
//        string result_path = "F:\\OpenFEM\\model\\noncomfortable\\dis_with_error.vtk";
//        cae_item.implict_analysis(result_path);
//    }
//    else if (case_num == 205)
//    {
//        // ----------------------------------------------------------------------------------------------------------- 
//        // SFEM
//        // ----------------------------------------------------------------------------------------------------------- 
//        // 材料属性赋值
//        CAE::elastic_mat mat_item{ 2.1e5, 0.3, 7800 };
//        // 材料路径
//        std::string path = "F:\\OpenFEM\\model\\CAppt\\full_analysis_fine.inp";
//        // 关键字
//        string load_set_keyword = "Set-load";
//        string load_value_keyword = "Cload";
//        string dis_set_keyword = "Set-fix";
//        // 建立CAE分析对象 
//        CAE::CAE_process cae_item(path, mat_item);
//        // 读取计算文件
//        cae_item.pre_info(load_set_keyword, load_value_keyword, dis_set_keyword);
//        // 执行结构响应分析
//        string result_path = "F:\\OpenFEM\\model\\CAppt\\full_analysis_fine_SFEM. vtk";
//        cae_item.implict_SFEManalysis(result_path, "");
//    }
//    else
//    {
//        std::cout << "Please check your code\n";
//    }
//}
