#include "BlocalShell.h"



BlocalShellRes BlocalShell(shapeFunction sf,mat H3i_Ele_Node,mat elemCoordinates,
    mat elemThickness,mat gaosi) {

    // 公式21，dudx=J-1*duds// 公式20，du'dx'=DCMT*dudx*DCM

    mat v3i = H3i_Ele_Node;
    mat tempi = {1, 0, 0} ;

    mat xv2i(8, 3);

    mat xv1i(8, 3);
    mat xv3i(8, 3);

    mat temp1;
    for (int i = 0; i < 8; ++i) {
        temp1 = cross(v3i.row(i).t(), tempi.t());
        xv2i.row(i) = temp1.t() / norm(temp1);
        temp1 = cross(xv2i.row(i), v3i.row(i));
        xv1i.row(i) = temp1 /norm(temp1);
        temp1 = v3i.row(i);

        xv3i.row(i) = temp1 / norm(temp1);


    }

    tempi = {1, 0, 0} ;

    mat v3 = sf.fun(0) * H3i_Ele_Node.row(0).t() + sf.fun(1) * H3i_Ele_Node.row(1).t() + 
        sf.fun(2) * H3i_Ele_Node.row(2).t() + sf.fun(3) * H3i_Ele_Node.row(3).t() + sf.fun(4) * H3i_Ele_Node.row(4).t() + sf.fun(5) 
        * H3i_Ele_Node.row(5).t() + sf.fun(6) * H3i_Ele_Node.row(6).t() + sf.fun(7) * H3i_Ele_Node.row(7).t();  //公式17，z'




    temp1 = cross(v3, tempi.t());
    double tn2 = sqrt(temp1(0) * temp1(0) + temp1(1) * temp1(1) + temp1(2) * temp1(2));
    mat xv2 = temp1 /tn2;
    temp1 = cross(xv2, v3.t());
    double tn1 = sqrt(temp1(0) * temp1(0) + temp1(1) * temp1(1) + temp1(2) * temp1(2));
    mat xv1 = temp1 /tn1;

    double tn3 = sqrt(v3(0) * v3(0) + v3(1) * v3(1) + v3(2) * v3(2));
    mat xv3 = v3 /tn3;

    mat theta = join_horiz(xv1, join_horiz(xv2, xv3));

    mat zmtemp = elemCoordinates.t();
    mat pni = sf.dfun;
    mat ni = sf.fun;

    mat JAC = zeros<mat>(3, 3);





    for (int I = 0; I < 8; ++I) {
        JAC(0, 0) = (zmtemp(0, I) + gaosi(2) * xv3(0) * 0.5 * elemThickness(I)) * pni(I, 0) + JAC(0, 0);
        JAC(0, 1) = (zmtemp(1, I) + gaosi(2) * xv3(1) * 0.5 * elemThickness(I)) * pni(I, 0) + JAC(0, 1);
        JAC(0, 2) = (zmtemp(2, I) + gaosi(2) * xv3(2) * 0.5 * elemThickness(I)) * pni(I, 0) + JAC(0, 2);
        JAC(1, 0) = (zmtemp(0, I) + gaosi(2) * xv3(0) * 0.5 * elemThickness(I)) * pni(I, 1) + JAC(1, 0);

        JAC(1, 1) = (zmtemp(1, I) + gaosi(2) * xv3(1) * 0.5 * elemThickness(I)) * pni(I, 1) + JAC(1, 1);
        JAC(1, 2) = (zmtemp(2, I) + gaosi(2) * xv3(2) * 0.5 * elemThickness(I)) * pni(I, 1) + JAC(1, 2);
        JAC(2, 0) = ni(I) * 0.5 * elemThickness(I) * xv3(0) + JAC(2, 0);
        JAC(2, 1) = ni(I) * 0.5 * elemThickness(I) * xv3(1) + JAC(2, 1);

        JAC(2, 2) = ni(I) * 0.5 * elemThickness(I) * xv3(2) + JAC(2, 2);

  }

    mat invj = inv(JAC);

    mat t1 = {pni(0, 0), pni(1, 0), pni(2, 0), pni(3, 0), pni(4, 0), pni(5, 0), pni(6, 0), pni(7, 0)};
    mat t2 = {pni(0, 1), pni(1, 1), pni(2, 1), pni(3, 1), pni(4, 1), pni(5, 1), pni(6, 1), pni(7, 1)};

    mat t3 = {0, 0, 0, 0, 0, 0, 0, 0};
    mat dni = invj * join_cols(t1, t2, t3);


    t1 = {gaosi(2) * pni(0, 0), gaosi(2) * pni(1, 0), gaosi(2) * pni(2, 0), gaosi(2) * pni(3, 0), 
      gaosi(2) * pni(4, 0), gaosi(2) * pni(5, 0), gaosi(2) * pni(6, 0), gaosi(2) * pni(7, 0)};
   t2 = {gaosi(2) * pni(0, 1), gaosi(2) * pni(1, 1), gaosi(2) * pni(2, 1), gaosi(2) * pni(3, 1), gaosi(2) * 
        pni(4, 1), gaosi(2) * pni(5, 1), gaosi(2) * pni(6, 1), gaosi(2) * pni(7, 1)};
    t3 = {ni(0), ni(1), ni(2), ni(3), ni(4), ni(5), ni(6), ni(7)};
    mat dmi = invj * join_cols(t1, t2, t3);



    mat btemp = zeros<mat>(5, 5);
    mat B = zeros<mat>(5, 40);
    for (int bi = 0; bi < 8; ++bi) 
    {
        mat alpha = {theta(0, 0) * dni(0, bi) + theta(1, 0) * dni(1, bi) + theta(2, 0) * dni(2, bi),
                     theta(0, 1) * dni(0, bi) + theta(1, 1) * dni(1, bi) + theta(2, 1) * dni(2, bi),
                     theta(0, 2) * dni(0, bi) + theta(1, 2) * dni(1, bi) + theta(2, 2) * dni(2, bi)};
        mat beta = {(theta(0, 0) * dmi(0, bi) + theta(1, 0) * dmi(1, bi) + theta(2, 0) * dmi(2, bi)) * elemThickness(bi) / 2,
                    (theta(0, 1) * dmi(0, bi) + theta(1, 1) * dmi(1, bi) + theta(2, 1) * dmi(2, bi)) * elemThickness(bi) / 2,
                    (theta(0, 2) * dmi(0, bi) + theta(1, 2) * dmi(1, bi) + theta(2, 2) * dmi(2, bi)) * elemThickness(bi) / 2};
        mat gama = {theta(0, 0) * xv1i(bi, 0) + theta(1, 0) * xv1i(bi, 1) + theta(2, 0) * xv1i(bi, 2),
                    theta(0, 1) * xv1i(bi, 0) + theta(1, 1) * xv1i(bi, 1) + theta(2, 1) * xv1i(bi, 2),
                    theta(0, 2) * xv1i(bi, 0) + theta(1, 2) * xv1i(bi, 1) + theta(2, 2) * xv1i(bi, 2)};
        mat lmd = {theta(0, 0) * xv2i(bi, 0) + theta(1, 0) * xv2i(bi, 1) + theta(2, 0) * xv2i(bi, 2),
                   theta(0, 1) * xv2i(bi, 0) + theta(1, 1) * xv2i(bi, 1) + theta(2, 1) * xv2i(bi, 2),
                   theta(0, 2) * xv2i(bi, 0) + theta(1, 2) * xv2i(bi, 1) + theta(2, 2) * xv2i(bi, 2)};

        btemp = {theta(0, 0) * alpha(0), theta(1, 0) * alpha(0), theta(2, 0) * alpha(0), beta(0) * gama(0), beta(0) * lmd(0),
                     theta(0, 1) * alpha(1), theta(1, 1) * alpha(1), theta(2, 1) * alpha(1), beta(1) * gama(1), beta(1) * lmd(1),
                     theta(0, 0) * alpha(1) + theta(0, 1) * alpha(0), theta(1, 0) * alpha(1) + theta(1, 1) * alpha(0), theta(2, 0) 
            * alpha(1) + theta(2, 1) * alpha(0), beta(0) * gama(1) + beta(1) * gama(0), beta(0) * lmd(1) + beta(1) * lmd(0),
                     theta(0, 1) * alpha(2) + theta(0, 2) * alpha(1), theta(1, 1) * alpha(2) + theta(1, 2) * alpha(1), theta(2, 1) * alpha(2) 
            + theta(2, 2) * alpha(1), beta(1) * gama(2) + beta(2) * gama(1), beta(1) * lmd(2) + beta(2) * lmd(1),
                     theta(0, 2) * alpha(0) + theta(0, 0) * alpha(2), theta(1, 2) * alpha(0) + theta(1, 0) * alpha(2), theta(2, 2) * alpha(0) 
            + theta(2, 0) * alpha(2), beta(2) * gama(0) + beta(0) * gama(2), beta(2) * lmd(0) + beta(0) * lmd(2)};
       
        
        
        
        // reshape to 5*5
        btemp = reshape(btemp, 5, 5).t();
        B(span(0, 4), span(bi * 5, bi * 5 + 4)) = btemp;

    
    
    
  }


    // djt=det(JAC);
    double djt = det(JAC);
    // djt_2d=det(JAC(1:2,1:2));
    double djt_2d = det(JAC(span(0, 1), span(0, 1)));  //




     BlocalShellRes res = {B, djt, djt_2d};






     return res;










}

