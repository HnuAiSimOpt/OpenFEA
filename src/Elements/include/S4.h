#ifndef STIFFNESSMATRIX_H
#define STIFFNESSMATRIX_H
#include "GaussQuadrature.h"
#include "MeshGenerator.h"
#include "Material.h"



mat StiffnessMatrix(Material material, Mesh mesh, Quadrature quadrature, string ShapeOption);






#endif //STIFFNESSMATRIX_H
