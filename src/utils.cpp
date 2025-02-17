/**************************************************************************

Copyright:  WH team

Author: YinJichao <jichaoyinyjc@163.com>

Completion date:  XXX

Description: XXX

**************************************************************************/

#include "include/utils.h"

using namespace std;
namespace CAE
{
    // 删除负数
    void delete_negative(vector<int> &vec)
    {
        for (auto it = vec.begin(); it != vec.end();)
        {
            if ((*it) < 0)
            {
                it = vec.erase(it);
            }
            else
            {
                ++it;
            }
        }
    };

    void write_spmat(vector<int> &row_idx, vector<int> &col_idx, vector<double> &nz_val, string path_row_value, string path_col)
    {
        std::ofstream outFile;
        outFile.open(path_row_value);
        if (!outFile.is_open())
        {
            std::cerr << "无法打开row_value文件: " << path_row_value << std::endl;
        }
        for (int i = 0; i < row_idx.size(); i++) // 写入一些文本到文件
        {
            outFile << row_idx[i] << ",    " << setw(12) << setiosflags(ios::fixed) << setprecision(8) << nz_val[i] << "\n";
        }
        outFile.close(); // 关闭文件
        //
        std::ofstream outFile2;
        outFile2.open(path_col);
        if (!outFile2.is_open())
        {
            std::cerr << "无法打开col文件: " << path_col << std::endl;
        }
        for (int i = 0; i < col_idx.size(); i++)
        {
            outFile2 << col_idx[i] << "\n";
        }
        outFile2.close();
    };

    void write_mat(vector<vector<double>> &density_mat, string path)
    {
        std::ofstream outFile;
        outFile.open(path);
        if (!outFile.is_open())
        {
            std::cerr << "无法打开carom文件: " << path << std::endl;
        }
        for (int i = 0; i < density_mat[0].size(); i++)
        {
            for (int j = 0; j < density_mat.size(); j++)
            {
                outFile << density_mat[j][i] << ",  ";
            }
            outFile << "\n";
        }
        outFile.close();
    };

    void write_vec(vector<double> &vec, string path)
    {
        std::ofstream outFile;
        outFile.open(path);
        if (!outFile.is_open())
        {
            std::cerr << "无法打开vec文件: " << path << std::endl;
        }
        for (int j = 0; j < vec.size(); j++)
        {
            outFile << vec[j] << "\n";
        }
        outFile.close();
    };

    void write_vec(vector<int> & vec, string path)
    {
        std::ofstream outFile;
        outFile.open(path);
        if (!outFile.is_open())
        {
            std::cerr << "无法打开vec文件: " << path << std::endl;
        }
        for (int j = 0; j < vec.size(); j++)
        {
            outFile << vec[j] << "\n";
        }
        outFile.close();
    };
}