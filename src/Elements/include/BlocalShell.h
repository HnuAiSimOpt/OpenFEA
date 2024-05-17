#ifndef BLOCALSHELL_H
#define BLOCALSHELL_H
#include "ShapeFunction.h"



struct BlocalShellRes {


    mat B;
    double djt;
    double djt_2d;


};


BlocalShellRes BlocalShell(shapeFunction sf,mat H3i_Ele_Node,mat elemCoordinates,
    mat elemThickness,mat gaosi);








#endif //BLOCALSHELL_H
