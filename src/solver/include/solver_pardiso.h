/**************************************************************************

Copyright:  WH team

Author: YinJichao <jichaoyinyjc@163.com>

Completion date:  XXX

Description: XXX

**************************************************************************/

#pragma once

#include<iostream>
#include <vector>
#include <mkl.h>
#include <mkl_pardiso.h>

using std::vector;

namespace CAE
{
    void pardiso_solver(vector<double> &nz_val, vector<int> &row_idx, vector<int> &col_idx, vector<double> &b, vector<double> &x);
}