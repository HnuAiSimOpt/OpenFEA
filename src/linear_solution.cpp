/**************************************************************************

Copyright:  WH team

Author: YinJichao <jichaoyinyjc@163.com>

Completion date:  XXX

Description: XXX

**************************************************************************/

#include "include/linear_solution.h"

namespace CAE
{
    bool solution_api(assamble_stiffness &item_ass, data_management &item_cae_data, string solver_type)
    {
        if (solver_type == "SuperLU_func")
        {
            superlu_solver_func(item_ass.nz_val_, item_ass.row_idx_, item_ass.col_idx_, item_cae_data.single_load_vec_, item_cae_data.single_dis_vec_);
        }
        else if (solver_type == "SuperLU_class")
        {
            item_cae_data.item_superlu = SuperLUSolution();
            // ----------------------------------------------------------------------------------------
            // vector<double> F_copy(item_cae_data.single_load_vec_.size(), 0.);
            // for(int i=0;i<item_cae_data.single_load_vec_.size();i++)
            //     F_copy[i] = 10. * item_cae_data.single_load_vec_[i];
            // ----------------------------------------------------------------------------------------
            // 初次求解
            bool a1 = item_cae_data.item_superlu.superlu_init(item_ass.nz_val_, item_ass.row_idx_, item_ass.col_idx_);
            // 求解
            vector<double> F_copy(item_cae_data.single_load_vec_.size(), 0.);
            for (int i = 0; i < item_cae_data.single_load_vec_.size(); i++)
                F_copy[i] = item_cae_data.single_load_vec_[i];
            bool a2 = item_cae_data.item_superlu.superlu_solution_1st(F_copy, item_cae_data.single_dis_vec_, false);
            // bool a3 = item_cae_data.item_superlu.superlu_solution_next(F_copy, item_cae_data.single_dis_vec_, true);
        }
        else if (solver_type == "Pardiso_func")
        {
            pardiso_solver_func(item_ass.nz_val_, item_ass.row_idx_, item_ass.col_idx_,
                                item_cae_data.single_load_vec_, item_cae_data.single_dis_vec_);
        }
        else if (solver_type == "Pardiso_class")
        {
            item_cae_data.item_pardiso = PardisoSolution();
            int n = int(item_cae_data.single_load_vec_.size());
            item_cae_data.item_pardiso.phase_00 = item_cae_data.item_pardiso.pardiso_init(item_ass.nz_val_, item_ass.row_idx_,
                                                                                          item_ass.col_idx_, n);
            item_cae_data.item_pardiso.phase_1122 = item_cae_data.item_pardiso.pardiso_decomposition();
            item_cae_data.item_pardiso.phase_33 = item_cae_data.item_pardiso.pardiso_solution(item_cae_data.single_load_vec_,
                                                                                              item_cae_data.single_dis_vec_);
        }
        else if (solver_type == "CA")
        {
        }
        else
        {
            std::cout << "Please check your flag of solver !!!";
            return false;
        }
        return true;
    }

    bool solution_nl_api(assamble_nl_stiffness &item_ass, data_management &item_cae_data, vector<double> &b, vector<double> &x, string solver_type)
    {
        if (solver_type == "SuperLU_func")
        {
            vector<double> F_copy(b.size(), 0.);
            for (int i = 0; i < b.size(); i++)
                F_copy[i] = b[i];
            superlu_solver_func(item_ass.nz_val, item_ass.row_idx, item_ass.col_idx, F_copy, x);
        }
        else if (solver_type == "SuperLU_class")
        {
            item_cae_data.item_superlu = SuperLUSolution();
            // 初次求解
            bool a1 = item_cae_data.item_superlu.superlu_init(item_ass.nz_val, item_ass.row_idx, item_ass.col_idx);
            // 求解
            vector<double> F_copy(b.size(), 0.);
            for (int i = 0; i < b.size(); i++)
                F_copy[i] = b[i];
            bool a2 = item_cae_data.item_superlu.superlu_solution_1st(F_copy, x, false);
            // bool a3 = item_cae_data.item_superlu.superlu_solution_next(F_copy, item_cae_data.single_dis_vec_, true);
        }
        else if (solver_type == "Pardiso_func")
        {
            pardiso_solver_func(item_ass.nz_val, item_ass.row_idx, item_ass.col_idx, b, x);
        }
        else if (solver_type == "Pardiso_class")
        {
            item_cae_data.item_pardiso = PardisoSolution();
            item_cae_data.item_pardiso.phase_00 = item_cae_data.item_pardiso.pardiso_init(item_ass.nz_val, item_ass.row_idx,
                                                                                          item_ass.col_idx, b, x);
            item_cae_data.item_pardiso.phase_1122 = item_cae_data.item_pardiso.pardiso_decomposition();
            item_cae_data.item_pardiso.phase_33 = item_cae_data.item_pardiso.pardiso_solution();
        }
        else if (solver_type == "CA")
        {
        }
        else
        {
            std::cout << "Please check your flag of solver !!!";
            return false;
        }
        return true;
    }
}
