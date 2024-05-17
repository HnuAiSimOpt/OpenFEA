#ifndef SHAPEFUNCTION_H
#define SHAPEFUNCTION_H
#include "armadillo"
using namespace arma;
struct shapeFunction {




	mat fun;
	mat dfun;


};


shapeFunction ShapeFunction(double s, double t);



#endif //SHAPEFUNCTION_H
