/**************************************************************************

Copyright:  WH team

Author: YinJichao <jichaoyinyjc@163.com>

Completion date:  XXX

Description: XXX

**************************************************************************/

#pragma once

#include "include/solver_linear.h"

namespace CAE
{
    void superlu_solver(assamble_stiffness &K, vector<double> &F, vector<double> &dis)
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
            int *ia = K.col_idx.data();
            int *ja = K.row_idx.data();
            double *a = K.nz_val.data();
            double *rhs = F.data();
            int neqs = F.size();
            int nz = K.num_nz_val;

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
            //	主要涉及A矩阵稀疏储存的变换方式，这个设置主要根据A来定，具体的选择涉及理论知识，小弟不才，实在不懂。

            /* Initialize the statistics variables. */
            StatInit(&stat);

            /* Solve the linear system. */
            int info;
            dgssv(&options, &A, perm_c, perm_r, &L, &U, &B, &stat, &info);

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