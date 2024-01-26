/**************************************************************************

Copyright:  WH team

Author: YinJichao <jichaoyinyjc@163.com>

Completion date:  XXX

Description: XXX

**************************************************************************/

#pragma once

#define min(x,y) (((x) < (y)) ? (x) : (y))

// Define the format to printf MKL_INT values
#if !defined(MKL_ILP64)
#define IFORMAT "%i"
#else
#define IFORMAT "%lli"
#endif

#include <iostream>
#include <Eigen/Dense>
#include "slu_ddefs.h"
#include <mkl.h>
#include <mkl_pardiso.h>

using namespace std;
using namespace Eigen;

bool sample_eigen();

bool sample_superlu();

bool sample_mkl();

bool sample_pardiso();