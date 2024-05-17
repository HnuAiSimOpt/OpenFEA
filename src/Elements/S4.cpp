#include "StiffnessMatrix.h"
#include "ShapeFunction.h"
#include "BlocalShell.h"
#include "TransShell8.h"





mat StiffnessMatrix(Material material, Mesh mesh, Quadrature quadrature, string ShapeOption) {


	int nn = mesh.Nid.n_rows;
	int ne = mesh.Eid.n_rows;
	int shape_order;
	
	if (ShapeOption == "Q4") {
		shape_order = 4;

   } else if (ShapeOption == "Q8") {
		shape_order = 8;


	} else if (ShapeOption == "Q9") {
		shape_order = 9;
	}
	
	mat K = zeros<mat>(nn * 5, nn * 5);
	mat M = zeros<mat>(nn * 5, nn * 5);


	double k = 1.2;//系数

	mat D = material.E / (1 - material.v * material.v) * mat({{1,          material.v,  0,                    0,     0},
	                                                          {material.v, 1,          0,                    0,     0},
	                                                          {0,          0,          (1 - material.v) / 2, 0,     0},
	                                                          {0,          0,          0,                    (1 -
	                                                                                                          material.v) /
	                                                                                                         2 / k, 0},
	                                                          {0,          0,          0,                    0,     (1 -
	                                                                                                                 material.v) /
	                                                                                                                2 /
	                                                                                                                k}}); //D

	
	
	mat KElem_loc = zeros<mat>(shape_order * 5, shape_order * 5);
	mat KElem = zeros<mat>(shape_order * 5, shape_order * 5);


	for (int iel = 0; iel < mesh.ne; iel++) 
	{
		KElem_loc = zeros<mat>(shape_order * 5, shape_order * 5);
		for (int ig = 0; ig < quadrature.points.n_rows; ig++) 
		{
			shapeFunction sf = ShapeFunction(quadrature.points(ig, 0), quadrature.points(ig, 1));
			mat elemCoordinates(mesh.Eid.n_cols -1, 3);
			mat elemThickness(mesh.Eid.n_cols -1, 1);
			mat H3i_Ele_Node(mesh.Eid.n_cols -1, mesh.H3i.n_cols);
			mat H2i_Ele_Node(mesh.Eid.n_cols -1, mesh.H2i.n_cols);
			mat H1i_Ele_Node(mesh.Eid.n_cols -1, mesh.H1i.n_cols);
			mat gaosi(1,quadrature.points.n_cols);


			for (int i = 1; i <= mesh.Eid.n_cols -1 ; ++i) 
			{
				elemCoordinates(i-1, span(0,2)) = mesh.Nid( mesh.Eid(iel, i)-1, span(1,3) );
				elemThickness(i-1, 0) = mesh.t(mesh.Eid(iel, i) - 1);
				H3i_Ele_Node(i-1,span(0,mesh.H3i.n_cols-1))=mesh.H3i(mesh.Eid(iel,i)-1,span(0,mesh.H3i.n_cols-1));
				H2i_Ele_Node(i-1,span(0,mesh.H2i.n_cols-1))=mesh.H2i(mesh.Eid(iel,i)-1,span(0,mesh.H2i.n_cols-1));
				H1i_Ele_Node(i-1,span(0,mesh.H1i.n_cols-1))=mesh.H1i(mesh.Eid(iel,i)-1,span(0,mesh.H1i.n_cols-1));
				gaosi = quadrature.points(ig,span(0,quadrature.points.n_cols-1));


			}

			BlocalShellRes bsr = BlocalShell(sf,H3i_Ele_Node,elemCoordinates,elemThickness,gaosi);

			// KElem_loc=KElem_loc+quadrature.weights(ig)*B'*D*B*detJ;% % 论文 公式28
			KElem_loc = KElem_loc + quadrature.weights(ig) * bsr.B.t() * D * bsr.B * bsr.djt;
			mat T = TransShell8(H1i_Ele_Node,H2i_Ele_Node,H3i_Ele_Node);

			// KElem=T'*KElem_loc*T;
			KElem = T.t() * KElem_loc * T;
	}


		imat ntot = mesh.Eid(iel, span(1,mesh.Eid.n_cols-1))*5-4; //index of DOF

		for (int j = 0; j < ntot.n_elem; ++j) {         //组装总刚
			for (int k = 0; k < ntot.n_elem; ++k) {
				int j_i = ntot(j) - 1;
				int k_i = ntot(k) - 1;
				K(span(j_i, 4 + j_i), span(k_i, 4 + k_i)) = K(span(j_i, 4 + j_i), span(k_i, 4 + k_i)) + KElem(span(j*5, 4 + j*5), span(k*5, 4 + k*5));
			}
	}

	}







	return K;
}
