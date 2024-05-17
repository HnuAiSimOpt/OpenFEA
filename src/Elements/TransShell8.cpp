#include "TransShell8.h"



mat TransShell8(mat v1i,mat v2i,mat v3i) {


   mat T = zeros<mat>(40, 40);
   mat theta;

   mat ri = zeros<mat>(5, 5);
   mat temp;

    for (int i = 0; i < 8; ++i) {

        theta = join_cols(v1i.row(i), join_cols(v2i.row(i), v3i.row(i)));
        ri(span(0, 2), span(0, 2)) = theta.t();

        temp = {{v1i(i,0), v1i(i,1)},{v2i(i,0), v2i(i,1)}};
         ri(span(3, 4), span(3, 4)) = temp;
        T(span(i * 5, 4 + i * 5), span(i * 5, 4 + i * 5)) = ri;


 }









 return T;

}
