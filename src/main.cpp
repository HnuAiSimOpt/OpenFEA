#pragma once

#include <iostream>
#include <Eigen/Dense>
#include "slu_ddefs.h"
#include "include/sample_eigen_superlu.h"

using namespace std;
using namespace Eigen;

void main()
{
    bool tt = sample_superlu();
}