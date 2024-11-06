/**************************************************************************

Copyright:  WH team

Author: YinJichao <jichaoyinyjc@163.com>

Completion date:  XXX

Description: XXX

**************************************************************************/

#pragma once

#include "include/solver_superlu.h"

namespace CAE
{
    // 初始化
    bool SuperLUSolution::superlu_init(vector<double> &nz_val, vector<int> &row_idx, vector<int> &col_idx)
    {
        /* Initialize matrix A. */
        ia_ = col_idx.data();
        ja_ = row_idx.data();
        a_ = nz_val.data();
        nz_ = int(nz_val.size());
        return true;
    }

    // 首次求解
    bool SuperLUSolution::superlu_solution_1st(vector<double> &b, vector<double> &x, bool free_mem)
    {
        neqs_ = int(b.size());
        double *rhs = b.data();
        SuperMatrix B;

        /*Create matrix A in the format expected by SuperLU.*/
        dCreate_CompCol_Matrix(&A_, neqs_, neqs_, nz_, a_, ja_, ia_, SLU_NC, SLU_D, SLU_GE);

        /* Create right-hand side matrix B. */
        int nrhs = 1;
        dCreate_Dense_Matrix(&B, neqs_, nrhs, rhs, neqs_, SLU_DN, SLU_D, SLU_GE);
        if (!(perm_r_ = intMalloc(neqs_)))
            ABORT("Malloc fails for perm_r[].");
        if (!(perm_c_ = intMalloc(neqs_)))
            ABORT("Malloc fails for perm_c[].");

        /* Set the default input options. */
        set_default_options(&options_);
        options_.ColPerm = COLAMD;
        /* Initialize the statistics variables. */
        StatInit(&stat_);

        /* Solve the linear system. */
        int info;
        dgssv(&options_, &A_, perm_c_, perm_r_, &L_, &U_, &B, &stat_, &info);
        // 传递解
        x = b;
        phase_ = 1;
        /* De-allocate storage */
        if (free_mem)
        {
            // SUPERLU_FREE(rhs);
            SUPERLU_FREE(perm_r_);
            SUPERLU_FREE(perm_c_);
            // Destroy_CompCol_Matrix(&A);
            // Destroy_SuperMatrix_Store(&B);
            Destroy_SuperNode_Matrix(&L_);
            Destroy_CompCol_Matrix(&U_);
            StatFree(&stat_);
        }
        return true;
    }

    // 利用首次求解的分解矩阵信息计算
    bool SuperLUSolution::superlu_solution_next(vector<double> &b, vector<double> &x, bool free_mem)
    {
        if (phase_)
        {
            double *rhs = b.data();
            trans_t trans = NOTRANS;
            int info1 = 0;
            SuperMatrix B;
            /* Create right-hand side matrix B. */
            int nrhs = 1;
            dCreate_Dense_Matrix(&B, neqs_, nrhs, rhs, neqs_, SLU_DN, SLU_D, SLU_GE);
            // 回代求解（无矩阵分解）
            dgstrs(trans, &L_, &U_, perm_c_, perm_r_, &B, &stat_, &info1);
            // 传递解
            x = b;

            /* De-allocate storage */
            if (free_mem)
            {
                // SUPERLU_FREE(rhs);
                SUPERLU_FREE(perm_r_);
                SUPERLU_FREE(perm_c_);
                // Destroy_CompCol_Matrix(&A);
                // Destroy_SuperMatrix_Store(&B);
                Destroy_SuperNode_Matrix(&L_);
                Destroy_CompCol_Matrix(&U_);
                StatFree(&stat_);
            }
            return true;
        }
        else
        {
            return false;
        }
    }

    // -----------------------------------------------------------------------------------------------------------------------------------------------------
    // 函数形式
    // -----------------------------------------------------------------------------------------------------------------------------------------------------
    void superlu_solver_func(vector<double> &nz_val, vector<int> &row_idx, vector<int> &col_idx, vector<double> &F, vector<double> &dis)
    {

        if (F.size() == 0)
        {
            printf("有无穷多解\n");
        }
        else
        {
            SuperMatrix A, L, U, B;
            int *perm_r; /* row permutations from partial pivoting */
            int *perm_c; /* column permutation vector */
            superlu_options_t options;
            SuperLUStat_t stat;

            /* Initialize matrix A. */
            int *ia = col_idx.data();
            int *ja = row_idx.data();
            double *a = nz_val.data();
            double *rhs = F.data();
            int neqs = int(F.size());
            int nz = int(nz_val.size());

            /*Create matrix A in the format expected by SuperLU.*/
            dCreate_CompCol_Matrix(&A, neqs, neqs, nz, a, ja, ia, SLU_NC, SLU_D, SLU_GE);

            /* Create right-hand side matrix B. */
            int nrhs = 1;
            dCreate_Dense_Matrix(&B, neqs, nrhs, rhs, neqs, SLU_DN, SLU_D, SLU_GE);
            if (!(perm_r = intMalloc(neqs)))
                ABORT("Malloc fails for perm_r[].");
            if (!(perm_c = intMalloc(neqs)))
                ABORT("Malloc fails for perm_c[].");

            /* Set the default input options. */
            set_default_options(&options);
            options.ColPerm = COLAMD;
            //  有限元方法求解时将options.ColPerm设置为COLAMD速度可以提高好几倍(SuperLU默认的ColPerm值就是COLAMD)，
            //	ColPer参数在ug中的描述是：Specifies how to permute the columns of the matrix for sparsity preservation.
            //	主要涉及A矩阵稀疏储存的变换方式，这个设置主要根据A来定。

            /* Initialize the statistics variables. */
            StatInit(&stat);

            /* Solve the linear system. */
            int info;
            dgssv(&options, &A, perm_c, perm_r, &L, &U, &B, &stat, &info);
            // 传递解
            dis = F;

            /* De-allocate storage */
            // SUPERLU_FREE(rhs);
            SUPERLU_FREE(perm_r);
            SUPERLU_FREE(perm_c);
            // Destroy_CompCol_Matrix(&A);
            // Destroy_SuperMatrix_Store(&B);
            Destroy_SuperNode_Matrix(&L);
            Destroy_CompCol_Matrix(&U);
            StatFree(&stat);
        }
    }
}