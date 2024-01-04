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
    bool PardisoSolution::pardiso_init(vector<double> &nz_val, vector<int> &row_idx, vector<int> &col_idx, int n_size)
    {
        /* -------------------------- */
        n = n_size;
        mtype = 11; /* Real and nonsymmetric matrix */
        nrhs = 1;   /* Number of right hand sides. */
        ia = col_idx.data();
        ja = row_idx.data();
        a = nz_val.data();
        /* PARDISO control parameters */
        for (int i = 0; i < 64; i++)
        {
            iparm[i] = 0;
        }
        iparm[0] = 1;  /* No solver default */
        iparm[1] = 3;  /* Fill-in reordering from METIS, Numbers of processors, value of OMP_NUM_THREADS */
        iparm[7] = 2;  /* Max numbers of iterative refinement steps */
        iparm[34] = 1; /* Zero-based indexing */
        for (int i = 0; i < 64; i++)
        {
            pt[i] = 0;
        }
        /* -------------------------- */
        maxfct = 1; /* Maximum number of numerical factorizations. */
        mnum = 1;   /* Which factorization to use. */
        msglvl = 0; /* Print statistical information: 0(not show) or 1(show)*/
        return true;
    }

    bool PardisoSolution::pardiso_decomposition()
    {
        /* ---------------------------------------------------------------------------------------------------------------
        分析
        ---------------------------------------------------------------------------------------------------------------- */
        MKL_INT error = 0; /* Initialize error flag */
        MKL_INT phase = 11;
        PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
                &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
        if (error != 0)
        {
            printf("\nERROR during symbolic factorization: " IFORMAT, error);
            exit(1);
        }
        std::cout << "Phase 11 has been finished ......\n";
        /* ---------------------------------------------------------------------------------------------------------------
        分解
        ---------------------------------------------------------------------------------------------------------------- */
        phase = 22;
        PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
                &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
        if (error != 0)
        {
            printf("\nERROR during numerical factorization: " IFORMAT, error);
            exit(2);
        }
        std::cout << "Phase 22 has been finished ......\n";
        return true;
    }

    bool PardisoSolution::pardiso_solution(vector<double> &b, vector<double> &x)
    {
        /* ---------------------------------------------------------------------------------------------------------------
        回代
        ---------------------------------------------------------------------------------------------------------------- */
        rhs = b.data();
        solution = x.data();
        MKL_INT error = 0; /* Initialize error flag */
        MKL_INT phase = 33;
        iparm[7] = 2; /* Max numbers of iterative refinement steps. */
        PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
                &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, rhs, solution, &error);
        if (error != 0)
        {
            printf("\nERROR during solution: " IFORMAT, error);
            exit(3);
        }
        std::cout << "Phase 33 has been finished ......\n";
    }

    void pardiso_solver_func(vector<double> &nz_val, vector<int> &row_idx, vector<int> &col_idx, vector<double> &b, vector<double> &x)
    {
        /*
        pardiso 初始化
        */
        MKL_INT n = int(b.size());
        int *ia = col_idx.data();
        int *ja = row_idx.data();
        double *a = nz_val.data();
        double *rhs = b.data();
        double *solution = x.data();
        //
        MKL_INT mtype = 11; /* Real and nonsymmetric matrix */
        MKL_INT nrhs = 1;   /* Number of right hand sides. */

        /* PARDISO control parameters. */
        MKL_INT iparm[64];
        for (int i = 0; i < 64; i++)
        {
            iparm[i] = 0;
        }
        iparm[0] = 1;  /* No solver default */
        iparm[1] = 3;  /* Fill-in reordering from METIS, Numbers of processors, value of OMP_NUM_THREADS */
        iparm[7] = 2;  /* Max numbers of iterative refinement steps */
        iparm[34] = 1; /* Zero-based indexing */
        //
        MKL_INT maxfct = 1; /* Maximum number of numerical factorizations. */
        MKL_INT mnum = 1;   /* Which factorization to use. */
        MKL_INT msglvl = 0; /* Print statistical information in file */
        MKL_INT error = 0;  /* Initialize error flag */
        MKL_INT idum;       /* Integer dummy. */
        double ddum;        /* Double dummy */
        MKL_INT phase;
        // Initialize the internal solver memory pointer. This is only necessary for the FIRST call of the PARDISO solver.
        void *pt[64];
        for (int i = 0; i < 64; i++)
        {
            pt[i] = 0;
        }
        std::cout << "Initialization has been finished ......\n";
        /* ---------------------------------------------------------------------------------------------------------------
        分析
        ---------------------------------------------------------------------------------------------------------------- */
        phase = 11;
        PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
                &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
        if (error != 0)
        {
            printf("\nERROR during symbolic factorization: " IFORMAT, error);
            exit(1);
        }
        printf("\nNumber of nonzeros in factors = " IFORMAT, iparm[17]);
        printf("\nNumber of factorization MFLOPS = " IFORMAT, iparm[18]);
        std::cout << "Phase 11 has been finished ......\n";
        /* ---------------------------------------------------------------------------------------------------------------
        分解
        ---------------------------------------------------------------------------------------------------------------- */
        phase = 22;
        PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
                &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
        if (error != 0)
        {
            printf("\nERROR during numerical factorization: " IFORMAT, error);
            exit(2);
        }
        std::cout << "Phase 22 has been finished ......\n";
        /* ---------------------------------------------------------------------------------------------------------------
        回代
        ---------------------------------------------------------------------------------------------------------------- */
        phase = 33;
        iparm[7] = 2; /* Max numbers of iterative refinement steps. */
        PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
                &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, rhs, solution, &error);
        if (error != 0)
        {
            printf("\nERROR during solution: " IFORMAT, error);
            exit(3);
        }
        std::cout << "Phase 33 has been finished ......\n";
    }
}