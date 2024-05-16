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
        if (solver_type == "SuperLU")
        {
            superlu_solver(item_ass.nz_val, item_ass.row_idx, item_ass.col_idx,
                           item_cae_data.single_load_vec_, item_cae_data.single_dis_vec_);
        }
        else if (solver_type == "Pardiso_func")
        {
            pardiso_solver_func(item_ass.nz_val, item_ass.row_idx, item_ass.col_idx,
                                item_cae_data.single_load_vec_, item_cae_data.single_dis_vec_);
        }
        else if (solver_type == "Pardiso_class")
        {
            item_cae_data.item_pardiso = PardisoSolution();
            int n = int(item_cae_data.single_load_vec_.size());
            item_cae_data.item_pardiso.phase_00 = item_cae_data.item_pardiso.pardiso_init(item_ass.nz_val, item_ass.row_idx,
                                                                                          item_ass.col_idx, n);
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
}
