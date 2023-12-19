/**************************************************************************

Copyright:  WH team

Author: YinJichao <jichaoyinyjc@163.com>

Completion date:  XXX

Description: XXX

**************************************************************************/

#include <iostream>
#include <vector>
#include "solver/include/solver_superlu.h"
#include "solver/include/solver_pardiso.h"

using std::string;
using std::vector;

namespace CAE
{
    bool solution_api(vector<double> &nz_val, vector<int> &row_idx,
                      vector<int> &col_idx, vector<double> &b, vector<double> &x, string solver_type);
}