/**************************************************************************

Copyright:  WH team

Author: YinJichao <jichaoyinyjc@163.com>

Completion date:  XXX

Description: XXX

**************************************************************************/

#pragma once

#include <vector>
#include "slu_ddefs.h"
#include "include/assemble.h"
#include "include/set_bcs.h"

using std::vector;

namespace CAE
{
    void superlu_solver(assamble_stiffness &A, vector<double> &b, vector<double> &x);
}