/**************************************************************************

Copyright:  WH team

Author: YinJichao <jichaoyinyjc@163.com>

Completion date:  XXX

Description: XXX

**************************************************************************/
#pragma once

#include <iostream>
#include <vector>
#include "include/data_management.h"
#include "include/assemble.h"
#include "include/nl_assemble.h"
#include "solver/include/solver_superlu.h"
#include "solver/include/solver_pardiso.h"

using std::string;
using std::vector;

namespace CAE
{
    bool solution_api(assamble_stiffness &item_ass, data_management &item_cae_data, string solver_type);

    bool solution_nl_api(assamble_nl_stiffness &item_ass, data_management &item_cae_data, vector<double> &b, vector<double> &x, string solver_type);
}