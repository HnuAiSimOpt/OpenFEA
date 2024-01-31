// /**************************************************************************

// Copyright:  WH team

// Author: YinJichao <jichaoyinyjc@163.com>

// Completion date:  XXX

// Description: XXX

// **************************************************************************/

// #pragma once

// #include "include/sample_eigen_superlu_mkl.h"

// bool sample_eigen()
// {
//     MatrixXd d_(3, 3);
//     MatrixXd m_(3, 3);

//     d_ = MatrixXd::Identity(3, 3);
//     d_(0, 1) = 3;

//     m_ = MatrixXd::Ones(3, 3);
//     m_(0, 1) = 3;

//     cout << "d_:" << endl
//          << d_ << endl;
//     cout << endl;
//     cout << "m_:" << endl
//          << m_ << endl;
//     cout << endl;
//     cout << "d_+m_:" << endl
//          << d_ + m_ << endl;
//     return true;
// }

// bool sample_superlu()
// {
//     /*
//      * Purpose
//      * =======
//      *
//      * This is the small 5x5 example used in the Sections 2 and 3 of the
//      * Users' Guide to illustrate how to call a SuperLU routine, and the
//      * matrix data structures used by SuperLU.
//      *
//      */
//     SuperMatrix A, L, U, B;
//     double *a, *rhs;
//     double s, u, p, e, r, l;
//     int *asub, *xa;
//     int *perm_r; /* row permutations from partial pivoting */
//     int *perm_c; /* column permutation vector */
//     int nrhs, info, i, m, n, nnz, permc_spec;
//     superlu_options_t options;
//     SuperLUStat_t stat;

//     /* Initialize matrix A. */
//     m = n = 5;
//     nnz = 12;
//     if (!(a = doubleMalloc(nnz)))
//         ABORT("Malloc fails for a[].");
//     if (!(asub = intMalloc(nnz)))
//         ABORT("Malloc fails for asub[].");
//     if (!(xa = intMalloc(n + 1)))
//         ABORT("Malloc fails for xa[].");
//     s = 19.0;
//     u = 21.0;
//     p = 16.0;
//     e = 5.0;
//     r = 18.0;
//     l = 12.0;
//     a[0] = s;
//     a[1] = l;
//     a[2] = l;
//     a[3] = u;
//     a[4] = l;
//     a[5] = l;
//     a[6] = u;
//     a[7] = p;
//     a[8] = u;
//     a[9] = e;
//     a[10] = u;
//     a[11] = r;
//     asub[0] = 0;
//     asub[1] = 1;
//     asub[2] = 4;
//     asub[3] = 1;
//     asub[4] = 2;
//     asub[5] = 4;
//     asub[6] = 0;
//     asub[7] = 2;
//     asub[8] = 0;
//     asub[9] = 3;
//     asub[10] = 3;
//     asub[11] = 4;
//     xa[0] = 0;
//     xa[1] = 3;
//     xa[2] = 6;
//     xa[3] = 8;
//     xa[4] = 10;
//     xa[5] = 12;

//     /* Create matrix A in the format expected by SuperLU. */
//     dCreate_CompCol_Matrix(&A, m, n, nnz, a, asub, xa, SLU_NC, SLU_D, SLU_GE);

//     /* Create right-hand side matrix B. */
//     nrhs = 1;
//     if (!(rhs = doubleMalloc(m * nrhs)))
//         ABORT("Malloc fails for rhs[].");
//     for (i = 0; i < m; ++i)
//         rhs[i] = 1.0;
//     dCreate_Dense_Matrix(&B, m, nrhs, rhs, m, SLU_DN, SLU_D, SLU_GE);

//     if (!(perm_r = intMalloc(m)))
//         ABORT("Malloc fails for perm_r[].");
//     if (!(perm_c = intMalloc(n)))
//         ABORT("Malloc fails for perm_c[].");

//     /* Set the default input options. */
//     set_default_options(&options);
//     options.ColPerm = NATURAL;

//     /* Initialize the statistics variables. */
//     StatInit(&stat);

//     /* Solve the linear system. */
//     dgssv(&options, &A, perm_c, perm_r, &L, &U, &B, &stat, &info);

//     dPrint_CompCol_Matrix("A", &A);
//     dPrint_CompCol_Matrix("U", &U);
//     dPrint_SuperNode_Matrix("L", &L);
//     print_int_vec("\nperm_r", m, perm_r);
//     dPrint_Dense_Matrix("B", &B);

//     /* De-allocate storage */
//     SUPERLU_FREE(rhs);
//     SUPERLU_FREE(perm_r);
//     SUPERLU_FREE(perm_c);
//     Destroy_CompCol_Matrix(&A);
//     Destroy_SuperMatrix_Store(&B);
//     Destroy_SuperNode_Matrix(&L);
//     Destroy_CompCol_Matrix(&U);
//     StatFree(&stat);

//     getchar();

//     return true;
// }

// bool sample_mkl()
// {
//     double *A, *B, *C;
//     int m, n, k, i, j;
//     double alpha, beta;

//     printf("\n This example computes real matrix C=alpha*A*B+beta*C using \n"
//            " Intel(R) MKL function dgemm, where A, B, and  C are matrices and \n"
//            " alpha and beta are double precision scalars\n\n");

//     m = 2000, k = 200, n = 1000;
//     printf(" Initializing data for matrix multiplication C=A*B for matrix \n"
//            " A(%ix%i) and matrix B(%ix%i)\n\n",
//            m, k, k, n);
//     alpha = 1.0;
//     beta = 0.0;

//     printf(" Allocating memory for matrices aligned on 64-byte boundary for better \n"
//            " performance \n\n");
//     A = (double *)mkl_malloc(m * k * sizeof(double), 64);
//     B = (double *)mkl_malloc(k * n * sizeof(double), 64);
//     C = (double *)mkl_malloc(m * n * sizeof(double), 64);
//     if (A == NULL || B == NULL || C == NULL)
//     {
//         printf("\n ERROR: Can't allocate memory for matrices. Aborting... \n\n");
//         mkl_free(A);
//         mkl_free(B);
//         mkl_free(C);
//         return 1;
//     }

//     printf(" Intializing matrix data \n\n");
//     for (i = 0; i < (m * k); i++)
//     {
//         A[i] = (double)(i + 1);
//     }

//     for (i = 0; i < (k * n); i++)
//     {
//         B[i] = (double)(-i - 1);
//     }

//     for (i = 0; i < (m * n); i++)
//     {
//         C[i] = 0.0;
//     }

//     printf(" Computing matrix product using Intel(R) MKL dgemm function via CBLAS interface \n\n");
//     cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
//                 m, n, k, alpha, A, k, B, n, beta, C, n);
//     printf("\n Computations completed.\n\n");

//     printf(" Top left corner of matrix A: \n");
//     for (i = 0; i < min(m, 6); i++)
//     {
//         for (j = 0; j < min(k, 6); j++)
//         {
//             printf("%12.0f", A[j + i * k]);
//         }
//         printf("\n");
//     }

//     printf("\n Top left corner of matrix B: \n");
//     for (i = 0; i < min(k, 6); i++)
//     {
//         for (j = 0; j < min(n, 6); j++)
//         {
//             printf("%12.0f", B[j + i * n]);
//         }
//         printf("\n");
//     }

//     printf("\n Top left corner of matrix C: \n");
//     for (i = 0; i < min(m, 6); i++)
//     {
//         for (j = 0; j < min(n, 6); j++)
//         {
//             printf("%12.5G", C[j + i * n]);
//         }
//         printf("\n");
//     }

//     printf("\n Deallocating memory \n\n");
//     mkl_free(A);
//     mkl_free(B);
//     mkl_free(C);

//     printf(" Example completed. \n\n");
//     return true;
// }

// bool sample_pardiso()
// {
//     /* Matrix data. */
//     MKL_INT n = 8;
//     MKL_INT ia[9] = {1, 5, 8, 10, 12, 15, 17, 18, 19};
//     MKL_INT ja[18] =
//         {1, 3, 6, 7,
//          2, 3, 5,
//          3, 8,
//          4, 7,
//          5, 6, 7,
//          6, 8,
//          7,
//          8};
//     double a[18] =
//         {7.0, 1.0, 2.0, 7.0,
//          -4.0, 8.0, 2.0,
//          1.0, 5.0,
//          7.0, 9.0,
//          5.0, 1.0, 5.0,
//          -1.0, 5.0,
//          11.0,
//          5.0};
//     MKL_INT mtype = -2; /* Real symmetric matrix */
//     /* RHS and solution vectors. */
//     double b[8], x[8];
//     MKL_INT nrhs = 1; /* Number of right hand sides. */
//     /* Internal solver memory pointer pt, */
//     /* 32-bit: int pt[64]; 64-bit: long int pt[64] */
//     /* or void *pt[64] should be OK on both architectures */
//     void *pt[64];
//     /* PARDISO control parameters. */
//     MKL_INT iparm[64];
//     MKL_INT maxfct, mnum, phase, error, msglvl;
//     /* Auxiliary variables. */
//     MKL_INT i;
//     double ddum;  /* Double dummy */
//     MKL_INT idum; /* Integer dummy. */
//                   /* -------------------------------------------------------------------- */
//                   /* .. Setup PARDISO control parameters. */
//                   /* -------------------------------------------------------------------- */
//     for (i = 0; i < 64; i++)
//     {
//         iparm[i] = 0;
//     }
//     iparm[0] = 1;   /* No solver default */
//     iparm[1] = 2;   /* Fill-in reordering from METIS */
//     iparm[3] = 0;   /* No iterative-direct algorithm */
//     iparm[4] = 0;   /* No user fill-in reducing permutation */
//     iparm[5] = 0;   /* Write solution into x */
//     iparm[6] = 0;   /* Not in use */
//     iparm[7] = 2;   /* Max numbers of iterative refinement steps */
//     iparm[8] = 0;   /* Not in use */
//     iparm[9] = 13;  /* Perturb the pivot elements with 1E-13 */
//     iparm[10] = 1;  /* Use nonsymmetric permutation and scaling MPS */
//     iparm[11] = 0;  /* Not in use */
//     iparm[12] = 0;  /* Maximum weighted matching algorithm is switched-off (default for symmetric). Try iparm[12] = 1 in case of inappropriate accuracy */
//     iparm[13] = 0;  /* Output: Number of perturbed pivots */
//     iparm[14] = 0;  /* Not in use */
//     iparm[15] = 0;  /* Not in use */
//     iparm[16] = 0;  /* Not in use */
//     iparm[17] = -1; /* Output: Number of nonzeros in the factor LU */
//     iparm[18] = -1; /* Output: Mflops for LU factorization */
//     iparm[19] = 0;  /* Output: Numbers of CG Iterations */
//     maxfct = 1;     /* Maximum number of numerical factorizations. */
//     mnum = 1;       /* Which factorization to use. */
//     msglvl = 1;     /* Print statistical information in file */
//     error = 0;      /* Initialize error flag */
//                     /* -------------------------------------------------------------------- */
//                     /* .. Initialize the internal solver memory pointer. This is only */
//                     /* necessary for the FIRST call of the PARDISO solver. */
//                     /* -------------------------------------------------------------------- */
//     for (i = 0; i < 64; i++)
//     {
//         pt[i] = 0;
//     }
//     /* -------------------------------------------------------------------- */
//     /* .. Reordering and Symbolic Factorization. This step also allocates */
//     /* all memory that is necessary for the factorization. */
//     /* -------------------------------------------------------------------- */
//     phase = 11;
//     PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
//             &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
//     if (error != 0)
//     {
//         printf("\nERROR during symbolic factorization: " IFORMAT, error);
//         exit(1);
//     }
//     printf("\nReordering completed ... ");
//     printf("\nNumber of nonzeros in factors = " IFORMAT, iparm[17]);
//     printf("\nNumber of factorization MFLOPS = " IFORMAT, iparm[18]);
//     /* -------------------------------------------------------------------- */
//     /* .. Numerical factorization. */
//     /* -------------------------------------------------------------------- */
//     phase = 22;
//     PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
//             &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
//     if (error != 0)
//     {
//         printf("\nERROR during numerical factorization: " IFORMAT, error);
//         exit(2);
//     }
//     printf("\nFactorization completed ... ");
//     /* -------------------------------------------------------------------- */
//     /* .. Back substitution and iterative refinement. */
//     /* -------------------------------------------------------------------- */
//     phase = 33;
//     iparm[7] = 2; /* Max numbers of iterative refinement steps. */
//     /* Set right hand side to one. */
//     for (i = 0; i < n; i++)
//     {
//         b[i] = 1;
//     }
//     PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
//             &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, b, x, &error);
//     if (error != 0)
//     {
//         printf("\nERROR during solution: " IFORMAT, error);
//         exit(3);
//     }
//     printf("\nSolve completed ... ");
//     printf("\nThe solution of the system is: ");
//     for (i = 0; i < n; i++)
//     {
//         printf("\n x [" IFORMAT "] = % f", i, x[i]);
//     }
//     printf("\n");
//     /* -------------------------------------------------------------------- */
//     /* .. Termination and release of memory. */
//     /* -------------------------------------------------------------------- */
//     phase = -1; /* Release internal memory. */
//     PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
//             &n, &ddum, ia, ja, &idum, &nrhs,
//             iparm, &msglvl, &ddum, &ddum, &error);
// }