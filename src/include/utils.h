/**************************************************************************

Copyright:  WH team

Author: YinJichao <jichaoyinyjc@163.com>

Completion date:  XXX

Description: XXX

**************************************************************************/

#pragma once
#include <set>
#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>

using namespace std;
using std::set;
namespace CAE
{
    // 删除负数
    void delete_negative(vector<int> &vec);

    void write_spmat(vector<int> &row_idx, vector<int> &col_idx, vector<double> &nz_val, string path_row_value, string path_col);

    void write_mat(vector<vector<double>> &density_mat, string path);

    void write_vec(vector<double> &vec, string path);

    void write_vec(vector<int> &vec, string path);

    // --------------------------------------------------------------------------------------------------------------------------------
    class SpMatrix
    {
    protected:
        int nRows_;
        int nColumns_;
        int nz_;               // 非零元素的总个数
        vector<double> nzval_; // 非零元素的值
        vector<int> rowind_;   // rowind记录每一个元素的行号,数量等于非零元素的个数(nz)
        vector<int> colptr_;   // colptr记录每一列第一个非零元素的行号【默认size = 总自由度+1】

    public:
        // 构造函数
        SpMatrix() : nRows_(0), nColumns_(0) {}
        SpMatrix(int m, int n) : nRows_(m), nColumns_(n) {}
        ~SpMatrix() {}

        // 构造稀疏函数
        void build_sp(int n, vector<set<int>> &columns, vector<double> &values);

        // 与向量做矩阵乘法
        void this_dot_vector(vector<double> &vec, vector<double> &answer);
    };
}