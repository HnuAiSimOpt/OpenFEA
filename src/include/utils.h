/**************************************************************************

Copyright:  WH team

Author: YinJichao <jichaoyinyjc@163.com>

Completion date:  XXX

Description: XXX

**************************************************************************/

#pragma once
#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>

using namespace std;
namespace CAE
{
    // 删除负数
    void delete_negative(vector<int> &vec);

    void write_spmat(vector<int> & row_idx, vector<int> &col_idx, vector<double> &nz_val, string path_row_value, string path_col);

    void write_mat(vector<vector<double>> & density_mat, string path);

    void write_vec(vector<double> & vec, string path);

    void write_vec(vector<int> & vec, string path);
}