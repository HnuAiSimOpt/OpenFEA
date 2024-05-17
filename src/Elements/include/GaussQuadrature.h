#ifndef GAUSSQUADRATURE_H
#define GAUSSQUADRATURE_H
#include <string>
#include <armadillo>
using namespace arma;
using namespace std;





struct Quadrature {

	mat points;
	vec weights;


};



Quadrature GaussQuadrature(const string &GaussOption);







#endif //GAUSSQUADRATURE_H
