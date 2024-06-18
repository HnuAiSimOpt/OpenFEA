/**************************************************************************

Copyright:  WH team

Author: YinJichao <jichaoyinyjc@163.com>

Completion date:  XXX

Description: XXX

**************************************************************************/

#pragma once

#include <vector>
#include "slu_ddefs.h"

using std::vector;

namespace CAE
{
    class SuperLUSolution
    {
    public:
        bool phase_ = 0;
        int *ia_;
        int *ja_;
        double *a_;
        int neqs_;
        int nz_;
        SuperMatrix A_, L_, U_;
        int *perm_r_; /* row permutations from partial pivoting */
        int *perm_c_; /* column permutation vector */
        superlu_options_t options_;
        SuperLUStat_t stat_;

    public:
        // 构造函数，析构函数
        SuperLUSolution() {};
        // 初次求解
        bool superlu_init(vector<double> &nz_val, vector<int> &row_idx, vector<int> &col_idx);
        // 求解
        bool superlu_solution_1st(vector<double> &b, vector<double> &x, bool free_mem);
        bool superlu_solution_next(vector<double> &b, vector<double> &x, bool free_mem);

    private:
    };
    void superlu_solver_func(vector<double> &nz_val, vector<int> &row_idx, vector<int> &col_idx, vector<double> &b, vector<double> &x);
}