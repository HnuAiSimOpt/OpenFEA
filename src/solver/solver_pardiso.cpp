/**************************************************************************

Copyright:  WH team

Author: YinJichao <jichaoyinyjc@163.com>

Completion date:  XXX

Description: XXX

**************************************************************************/

#pragma once

#include "include/solver_pardiso.h"

namespace CAE
{
    void pardiso_solver(vector<double> &nz_val, vector<int> &row_idx, vector<int> &col_idx, vector<double> &b, vector<double> &x)
    {
        std::cout << "The pardiso solver is perfomed ......\n";
        /*
        pardiso 初始化
        */
        int *ia = col_idx.data();
        int *ja = row_idx.data();
        double *a = nz_val.data();
        double *rhs = b.data();
        double *solution = x.data();
        //
        void *pt[64];
        pt[1] = 0;
        int mtype = 1; /* 矩阵类型  实数对称正定时为1   实数非对称时为11   复数对称时为3   复数非对称时为13 */
        int m_iparm[64];
        for (int i = 0; i < 64; i++)
        {
            m_iparm[i] = 0;
        }
        m_iparm[59] = 1;
        pardisoinit(pt, &mtype, m_iparm);
        std::cout << "Initialization has been finished ......\n";
        /*
        分析
        */
        int n = int(b.size());
        int nnz = int(nz_val.size());
        int nrhs = 1;
        int maxfct = 1; /* Maximum number of numerical factorizations */
        int mnum = 1;   /* Which factorization to use */
        int msglvl = 0; /* 0 Suppress printing, 1 Print statistical information */
        int phase = 11;
        int error = 0;
        pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ja, ia, NULL, &nrhs, m_iparm, &msglvl, NULL, NULL, &error);
        std::cout << "Phase 11 has been finished ......\n";
        /*
        分解
        */
        phase = 22;
        bool m_iparm3 = false;
        m_iparm[3] = (m_iparm3 ? 61 : 0);
        error = 0;
        pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, NULL, &nrhs, m_iparm, &msglvl, NULL, NULL, &error);
        std::cout << "Phase 22 has been finished ......\n";
        /*
        回代
        */
        phase = 33;
        m_iparm[7] = 1; /* Maximum number of iterative refinement steps */
        error = 0;
        pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, NULL, &nrhs, m_iparm, &msglvl, rhs, solution, &error);
        std::cout << "Phase 33 has been finished ......\n";
    }
}