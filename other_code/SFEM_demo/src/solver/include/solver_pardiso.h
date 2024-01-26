/**************************************************************************

Copyright:  WH team

Author: YinJichao <jichaoyinyjc@163.com>

Completion date:  XXX

Description: XXX

**************************************************************************/

#pragma once

#if !defined(MKL_ILP64)
#define IFORMAT "%i"
#else
#define IFORMAT "%lli"
#endif

#include <iostream>
#include <vector>
#include <mkl.h>
#include <mkl_pardiso.h>

using std::vector;

namespace CAE
{
    class PardisoSolution
    {
    public:
        bool phase_00;
        bool phase_1122;
        bool phase_33;

    public:
        // 构造函数，析构函数
        PardisoSolution(){};
        // 数据预处理
        bool pardiso_init(vector<double> &nz_val, vector<int> &row_idx, vector<int> &col_idx, vector<double> &b, vector<double> &x);
        // 分解矩阵
        bool pardiso_decomposition();
        // 求解
        bool pardiso_solution();

    private:
        /* -------------------------- */
        MKL_INT n;
        MKL_INT mtype;
        MKL_INT nrhs;
        int *ia;
        int *ja;
        double *a;
        double *rhs;
        double *solution;
        /* PARDISO control parameters */
        MKL_INT iparm[64];
        void *pt[64];
        /* -------------------------- */
        MKL_INT maxfct;
        MKL_INT mnum;
        MKL_INT msglvl;
        MKL_INT idum;
        double ddum;
    };

    void pardiso_solver_func(vector<double> &nz_val, vector<int> &row_idx, vector<int> &col_idx, vector<double> &b, vector<double> &x);
}