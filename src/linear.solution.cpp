/**************************************************************************

Copyright:  WH team

Author: YinJichao <jichaoyinyjc@163.com>

Completion date:  XXX

Description: XXX

**************************************************************************/

#include "include/linear_solution.h"

namespace CAE
{
    bool solution_api(vector<double> &nz_val, vector<int> &row_idx, vector<int> &col_idx, vector<double> &b,
                      vector<double> &x, string solver_type)
    {
        if (solver_type == "SuperLU")
        {
            superlu_solver(nz_val, row_idx, col_idx, b, x);
        }
        else if (solver_type == "Pardiso_func")
        {
            pardiso_solver_func(nz_val, row_idx, col_idx, b, x);
        }
        else if (solver_type == "Pardiso_class")
        {
            PardisoSolution item_pardiso;
            item_pardiso.phase_00 = item_pardiso.pardiso_init(nz_val, row_idx, col_idx, b, x);
            item_pardiso.phase_1122 = item_pardiso.pardiso_decomposition();
            item_pardiso.phase_33 = item_pardiso.pardiso_solution();
        }
        else
        {
            std::cout << "Please check your flag of solver !!!";
            return false;
        }
        return true;
    }
}
