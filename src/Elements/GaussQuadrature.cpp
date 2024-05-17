#include "GaussQuadrature.h"



Quadrature GaussQuadrature(const string &GaussOption) {


	Quadrature quadrature;
	
	if (GaussOption == "gauss1") {


		quadrature.points = zeros<mat>(1, 3);

		quadrature.weights = ones<vec>(1, 1) * 8;
	
	} 
	else if (GaussOption == "gauss2") {

		quadrature.points = zeros<mat>(8, 3);
		quadrature.points.col(0) = {1, 1, -1, -1, 1, 1, -1, -1};
		quadrature.points.col(1) = {1, -1, 1, -1, 1, -1, 1, -1};

		quadrature.points.col(2) = {1, 1, 1, 1, -1, -1, -1, -1};
		quadrature.points /= sqrt(3);
		quadrature.weights = ones<vec>(8, 1);
}
	
	return quadrature;
}
