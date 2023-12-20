/**************************************************************************

Copyright:  WH team

Author: YinJichao <jichaoyinyjc@163.com>、ChenBinXiang <chen_mech99@163.com>

Completion date:  XXX

Description: XXX---------------------

**************************************************************************/
#include <vector>
#pragma once

#include "include/P4.h"

namespace CAE
{
    REGISTER(ele_base, tetra_ele_elastic, "C3D3");
    // 建立本构矩阵



    std::vector<std::vector<double>> hooke(double E, double v) {
    double factor = E / (1 - v * v);

    std::vector<std::vector<double>> D(3, std::vector<double>(3));

    D[0][0] = factor;
    D[0][1] = factor * v;
    D[0][2] = 0;

    D[1][0] = factor * v;
    D[1][1] = factor;
    D[1][2] = 0;

    D[2][0] = 0;
    D[2][1] = 0;
    D[2][2] = factor * (1 - v) / 2;

    return D;
}
    // 生成单元刚度矩阵
vector<vector<double>> platere( std::vector<double> ex, std::vector<double> ey, double ep,  std::vector<std::vector<double>> D) {
    double Lx = ex[2] - ex[0];
    double Ly = ey[2] - ey[0];
    double t = ep;

    // Calculate D matrix
    double D_val = (t * t * t) / 12.0;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            D[i][j] *= D_val;
        }
    }

    // Define coefficients A1 to A9
    double A1 = Ly / (Lx * Lx * Lx);
    double A2 = Lx / (Ly * Ly * Ly);
    double A3 = 1.0 / (Lx * Ly);
    double A4 = Ly / (Lx * Lx);
    double A5 = Lx / (Ly * Ly);
    double A6 = 1.0 / Lx;
    double A7 = 1.0 / Ly;
    double A8 = Ly / Lx;
    double A9 = Lx / Ly;

    // Define coefficients C1 to C21
    double C1 = 4.0 * A1 * D[0][0] + 4.0 * A2 * D[1][1] + 2.0 * A3 * D[0][1] + 5.6 * A3 * D[2][2];
    double C2 = -4.0 * A1 * D[0][0] + 2.0 * A2 * D[1][1] - 2.0 * A3 * D[0][1] - 5.6 * A3 * D[2][2];
    double C3 = 2.0 * A1 * D[0][0] - 4.0 * A2 * D[1][1] - 2.0 * A3 * D[0][1] - 5.6 * A3 * D[2][2];
    double C4 = -2.0 * A1 * D[0][0] - 2.0 * A2 * D[1][1] + 2.0 * A3 * D[0][1] + 5.6 * A3 * D[2][2];
    double C5 = 2.0 * A5 * D[1][1] + A6 * D[0][1] + 0.4 * A6 * D[2][2];
    double C6 = 2.0 * A4 * D[0][0] + A7 * D[0][1] + 0.4 * A7 * D[2][2];
    double C7 = 2.0 * A5 * D[1][1] + 0.4 * A6 * D[2][2];
    double C8 = 2.0 * A4 * D[0][0] + 0.4 * A7 * D[2][2];
    double C9 = A5 * D[1][1] - A6 * D[0][1] - 0.4 * A6 * D[2][2];
    double C10 = A4 * D[0][0] - A7 * D[0][1] - 0.4 * A7 * D[2][2];
    double C11 = A5 * D[1][1] - 0.4 * A6 * D[2][2];
    double C12 = A4 * D[0][0] - 0.4 * A7 * D[2][2];
    double C13 = (4.0 / 3.0) * A9 * D[1][1] + (8.0 / 15.0) * A8 * D[2][2];
    double C14 = (4.0 / 3.0) * A8 * D[0][0] + (8.0 / 15.0) * A9 * D[2][2];
    double C15 = (2.0 / 3.0) * A9 * D[1][1] - (8.0 / 15.0) * A8 * D[2][2];
    double C16 = (2.0 / 3.0) * A8 * D[0][0] - (8.0 / 15.0) * A9 * D[2][2];
    double C17 = (2.0 / 3.0) * A9 * D[1][1] - (2.0 / 15.0) * A8 * D[2][2];
    double C18 = (2.0 / 3.0) * A8 * D[0][0] - (2.0 / 15.0) * A9 * D[2][2];
    double C19 = (1.0 / 3.0) * A9 * D[1][1] + (2.0 / 15.0) * A8 * D[2][2];
    double C20 = (1.0 / 3.0) * A8 * D[0][0] + (2.0 / 15.0) * A9 * D[2][2];
    double C21 = D[0][1];

    // Create the 12x12 stiffness matrix Keq
    vector<vector<double>> Keq(12, vector<double>(12, 0.0));

    Keq[0][0] = C1;
    Keq[0][1] = C5;
    Keq[0][2] = -C6;
    Keq[0][3] = C2;
    Keq[0][4] = C9;
    Keq[0][5] = -C8;
    Keq[0][6] = C4;
    Keq[0][7] = C11;
    Keq[0][8] = -C12;
    Keq[0][9] = C3;
    Keq[0][10] = C7;
    Keq[0][11] = -C10;

    Keq[1][1] = C13;
    Keq[1][2] = -C21;
    Keq[1][3] = C9;
    Keq[1][4] = C15;
    Keq[1][5] = 0; 
    Keq[1][6] = -C11; 
    Keq[1][7] =C19;
    Keq[1][8] = 0; 
    Keq[1][9] = -C7;
    Keq[1][10] = C17;
    Keq[1][11] = 0;


    Keq[2][2] = C14;
    Keq[2][3] = C8;
    Keq[2][4] = 0;
    Keq[2][5] = C18;
    Keq[2][6] = C12;
    Keq[2][7] = 0;
    Keq[2][8] = C20;
    Keq[2][9] =-C10;
    Keq[2][10] = 0;
    Keq[2][11] = C16;


    Keq[3][3] = C2;
    Keq[3][4] = -C9;
    Keq[3][5] = C8;
    Keq[3][6] = -C4;
    Keq[3][7] = C11;
    Keq[3][8] = -C7;
    Keq[3][9] = C2;
    Keq[3][10] = C2;
    Keq[3][11] = C2;


    Keq[3][3] = C1;
    Keq[3][4] = C5;
    Keq[3][5] = C6;
    Keq[3][6] = C3;
    Keq[3][7] = C7;
    Keq[3][8] = C10;
    Keq[3][9] = C4;
    Keq[3][10] = C11;
    Keq[3][11] = C12;



    Keq[4][4] = C16;
    Keq[4][5] = -C21;
    Keq[4][6] = -C21;
    Keq[4][7] = C12;
    Keq[4][8] = -C21;
    Keq[4][9] = -C19;
    Keq[4][10] = -C21;
    Keq[4][11] = -C21;



    Keq[4][4] = C13;
    Keq[4][5] = C21;
    Keq[4][6] = -C7;
    Keq[4][7] = C17;
    Keq[4][8] = 0;
    Keq[4][9] = -C11;
    Keq[4][10] =C19;
    Keq[4][11] =0;


    Keq[5][5] = C14;
    Keq[5][6] = C10;
    Keq[5][7] = 0;
    Keq[5][8] = C16;
    Keq[5][9] = -C12;
    Keq[5][10] = 0;
    Keq[5][11] = C20;


    Keq[6][6] = C3;
    Keq[6][7] = -C11;
    Keq[6][8] = C10;
    Keq[6][8] = C10;
    Keq[6][8] = C10;
    Keq[6][8] = C10;

    
    Keq[6][6] = C1;
    Keq[6][7] = -C5;
    Keq[6][8] = C6;
    Keq[6][9] = C2;
    Keq[6][10] = -C9;
    Keq[6][11] = C8;



    Keq[7][7] = C13;
    Keq[7][8] = -C21;
    Keq[7][9] = -C9;
    Keq[7][10] = C15;
    Keq[7][11] = 0;



    Keq[8][8] = C14;
    Keq[8][9] = -C8;
    Keq[8][10] = 0;
    Keq[8][11] = C18;



    Keq[9][9] = C1;
    Keq[9][10] = -C5;
    Keq[9][11] = -C6;




    Keq[10][10] = C13;
    Keq[10][11] = C21;

    Keq[11][11] = C14;

    // Fill the lower triangle of the symmetric matrix
    for (int i = 1; i < 12; i++) {
        for (int j = 0; j < i; j++) {
            Keq[i][j] = Keq[j][i];
        }
    }

    return Keq;
}

}


