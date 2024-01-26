/**************************************************************************

Copyright:  WH team

Author: ChenBinXiang <chen_mech99@163.com>

Completion date:  XXX

Description: XXX

**************************************************************************/


#include "include/SFEM_C3D4.h"
#include "../include/SFEM3D.h"
#include "../include/SupportForSFEM.h"
#include "../include/data_management.h"
#define M_PI  3.14159265358979323846;

namespace CAE
{
	//REGISTER(ele_base, SFEM_C3D4, "SFEM_C3D4");
	double SFEM_C3D4::H[4][4] = {	{ 0.62200846792814624, 0.16666666666666663, 0.044658198738520435, 0.16666666666666663 },
									{ 0.16666666666666663, 0.62200846792814624, 0.16666666666666663, 0.044658198738520435 },
									{ 0.044658198738520435, 0.16666666666666663, 0.62200846792814624, 0.16666666666666663 },
									{ 0.16666666666666663, 0.044658198738520435, 0.16666666666666663, 0.62200846792814624 } };

	double SFEM_C3D4::Gr[4][4] = {	{ -0.39433756729740643, 0.39433756729740643, 0.10566243270259354, -0.10566243270259354 },
									{ -0.39433756729740643, 0.39433756729740643, 0.10566243270259354, -0.10566243270259354 },
									{ -0.10566243270259354,  0.10566243270259354, 0.39433756729740643, -0.39433756729740643 },
									{ -0.10566243270259354, 0.10566243270259354,  0.39433756729740643, -0.39433756729740643 }};

	double SFEM_C3D4::Gs[4][4] = {  { -0.39433756729740643, -0.10566243270259354, 0.10566243270259354, 0.39433756729740643 },
									{ -0.10566243270259354, -0.39433756729740643,  0.39433756729740643,  0.10566243270259354 },
									{ -0.10566243270259354, -0.39433756729740643,  0.39433756729740643, 0.10566243270259354 },
									{ -0.39433756729740643, -0.10566243270259354, 0.10566243270259354, 0.39433756729740643 } };



	void SFEM_C3D4::build_cons_mat(double D[6][6])
	{
		double em = matrial_struc_.young_modulus;
		double nu = matrial_struc_.poisson_ratio;
		double fac = (em * (1.0 - nu)) / ((1.0 + nu) * (1.0 - 2.0 * nu));
		double fac_a = fac * nu / (1.0 - nu);
		double fac_b = fac * (1.0 - 2.0 * nu) / (2.0 * (1.0 - nu));
		D[0][0] = fac; D[0][1] = fac_a; D[0][2] = fac_a;
		D[1][0] = fac_a; D[1][1] = fac; D[1][2] = fac_a;
		D[2][0] = fac_a; D[2][1] = fac_a; D[2][2] = fac;
		D[3][3] = fac_b;
		D[4][4] = fac_b;
		D[5][5] = fac_b;
	}

    // 建立单元刚度矩阵
    void SFEM_C3D4::build_ele_stiff_mat(vector<vector<double>>& ke)
    {
		//计算形函数偏导
		Point G[50];	Mat3 GG[50];
		shapeGradient(G, GG);

		//材料矩阵
		double D[6][6] = { 0.0 };
		build_cons_mat(D);

		//计算刚度矩阵
		double dv = sfemData->nodeVolume[n] / 6.;
		nnode_ = sfemData->nodeArNode[n].size();
		ke.resize(nnode_ * 3);
		for (auto& i : ke) { i.assign(nnode_ * 3, 0.0); }
		for(int m = 0; m < 6; m++){
			double DB[6][3] = {0.};
			double sign = m <= 2 ? -1.: 1.;
			double in = m <= 2 ? m : m - 3;
			for(int i = 0, i3 = 0; i < nnode_; ++i, i3 += 3){
				double Gxi = G[i].x + sign * GG[i](0, in);
				double Gyi = G[i].y + sign * GG[i](1, in);
				double Gzi = G[i].z + sign * GG[i](2, in);
				for(int j = 0, j3 = 0; j < nnode_; ++j, j3 += 3){
					double Gxj = G[j].x + sign * GG[j](0, in);
					double Gyj = G[j].y + sign * GG[j](1, in);
					double Gzj = G[j].z + sign * GG[j](2, in);

					DB[0][0] = (D[0][0] * Gxj + D[0][4] * Gzj + D[0][5] * Gyj);
					DB[0][1] = (D[0][1] * Gyj + D[0][3] * Gzj + D[0][5] * Gxj);
					DB[0][2] = (D[0][2] * Gzj + D[0][3] * Gyj + D[0][4] * Gxj);

					DB[1][0] = (D[1][0] * Gxj + D[1][4] * Gzj + D[1][5] * Gyj);
					DB[1][1] = (D[1][1] * Gyj + D[1][3] * Gzj + D[1][5] * Gxj);
					DB[1][2] = (D[1][2] * Gzj + D[1][3] * Gyj + D[1][4] * Gxj);

					DB[2][0] = (D[2][0] * Gxj + D[2][4] * Gzj + D[2][5] * Gyj);
					DB[2][1] = (D[2][1] * Gyj + D[2][3] * Gzj + D[2][5] * Gxj);
					DB[2][2] = (D[2][2] * Gzj + D[2][3] * Gyj + D[2][4] * Gxj);

					DB[3][0] = (D[3][0] * Gxj + D[3][4] * Gzj + D[3][5] * Gyj);
					DB[3][1] = (D[3][1] * Gyj + D[3][3] * Gzj + D[3][5] * Gxj);
					DB[3][2] = (D[3][2] * Gzj + D[3][3] * Gyj + D[3][4] * Gxj);

					DB[4][0] = (D[4][0] * Gxj + D[4][4] * Gzj + D[4][5] * Gyj);
					DB[4][1] = (D[4][1] * Gyj + D[4][3] * Gzj + D[4][5] * Gxj);
					DB[4][2] = (D[4][2] * Gzj + D[4][3] * Gyj + D[4][4] * Gxj);

					DB[5][0] = (D[5][0] * Gxj + D[5][4] * Gzj + D[5][5] * Gyj);
					DB[5][1] = (D[5][1] * Gyj + D[5][3] * Gzj + D[5][5] * Gxj);
					DB[5][2] = (D[5][2] * Gzj + D[5][3] * Gyj + D[5][4] * Gxj);

					ke[i3 + 0][j3 + 0] += (Gxi * DB[0][0] + Gzi * DB[4][0] + Gyi * DB[5][0]) * dv;
					ke[i3 + 0][j3 + 1] += (Gxi * DB[0][1] + Gzi * DB[4][1] + Gyi * DB[5][1]) * dv;
					ke[i3 + 0][j3 + 2] += (Gxi * DB[0][2] + Gzi * DB[4][2] + Gyi * DB[5][2]) * dv;
					ke[i3 + 1][j3 + 0] += (Gyi * DB[1][0] + Gzi * DB[3][0] + Gxi * DB[5][0]) * dv;
					ke[i3 + 1][j3 + 1] += (Gyi * DB[1][1] + Gzi * DB[3][1] + Gxi * DB[5][1]) * dv;
					ke[i3 + 1][j3 + 2] += (Gyi * DB[1][2] + Gzi * DB[3][2] + Gxi * DB[5][2]) * dv;
					ke[i3 + 2][j3 + 0] += (Gzi * DB[2][0] + Gyi * DB[3][0] + Gxi * DB[4][0]) * dv;
					ke[i3 + 2][j3 + 1] += (Gzi * DB[2][1] + Gyi * DB[3][1] + Gxi * DB[4][1]) * dv;
					ke[i3 + 2][j3 + 2] += (Gzi * DB[2][2] + Gyi * DB[3][2] + Gxi * DB[4][2]) * dv;
				}
			}
		}
    }

	void SFEM_C3D4::shapeGradient(Point* G, Mat3* GG)
	{
		Point t4g[4];
		bool boundary = sfemData->surfaceNode.find(n) == sfemData->surfaceNode.end()?false : true;
	    for (int i = 0; i < sfemData->nodeElement[n].size(); i++) {//遍历节点所在所有单元
	    	int element = sfemData->nodeElement[n][i];//单元编号
	    	vector<int> nodes = sfemData->data_cae->node_topos_[element];//单元节点
			int nd[4];//旋转后的单元节点编号
	    	//基于节点在单元的位置旋转单元
	    	if (n == nodes[0]) {
				nd[0] = nodes[0];
				nd[1] = nodes[1];
				nd[2] = nodes[2];
				nd[3] = nodes[3];
	    	}
	    	else if (n == nodes[1]) {
				nd[0] = nodes[1];
				nd[1] = nodes[2];
				nd[2] = nodes[3];
				nd[3] = nodes[0];
	    	}
	    	else if (n == nodes[2]) {
				nd[0] = nodes[2];
				nd[1] = nodes[3];
				nd[2] = nodes[0];
				nd[3] = nodes[1];
	    	}
   			 else if (n == nodes[3]) {
				nd[0] = nodes[3];
				nd[1] = nodes[0];
				nd[2] = nodes[1];
				nd[3] = nodes[2];
	    	}
			//获取节点坐标
			Point node1 = *sfemData->nodeCoord[nd[0]];
			Point node2 = *sfemData->nodeCoord[nd[1]];
			Point node3 = *sfemData->nodeCoord[nd[2]];
			Point node4 = *sfemData->nodeCoord[nd[3]];

			Point xm1 = (node1 + node2) / 2.;//边中点坐标
			Point xm2 = (node1 + node3) / 2.;//边中点坐标
			Point xm3 = (node1 + node4) / 2.;//边中点坐标
			Point xm4 = (node2 + node2) / 2.;//边中点坐标
			Point xm5 = (node3 + node2) / 2.;//边中点坐标
			Point xm6 = (node2 + node2) / 2.;//边中点坐标

			Point xc1 = (node1 + node2 + node3) / 3.;//计算单元面心坐标
			Point xc2 = (node1 + node3 + node4) / 3.;//计算单元面心坐标
			Point xc3 = (node1 + node2 + node4) / 3.;//计算单元面心坐标
			Point xc4 = (node2 + node3 + node4) / 3.;//计算单元面心坐标
			Point xcc = (node1 + node2 + node3 + node4) / 4.;//计算单元中心坐标

			if (GG) {
				this->T4G(t4g, node1, node2, node3, node4);
			}

			vector<Point> lCoord; lCoord.resize(4);
			double phi_c[4][4] = { 0. };
			for(int isurf = 1; isurf <= 3; isurf++){//单元边界面循环
				if(isurf == 1){
					lCoord[0] = xc3;
					lCoord[1] = xm1;
					lCoord[2] = xc1;
					lCoord[3] = xcc;
					phi_c[0][0] = 1. /3.;		phi_c[1][0] = 1. /2.;		phi_c[2][0] = 1. /3.;		phi_c[3][0] = 1. /4.;
					phi_c[0][1] = 1. /3.;		phi_c[1][1] =  1. /2.;		phi_c[2][1] = 1. /3.;		phi_c[3][1] = 1. /4.;
					phi_c[0][2] = 0.;			phi_c[1][2] = 0.;			phi_c[2][2] = 1. /3.;		phi_c[3][2] = 1. /4.;
					phi_c[0][3] = 1. /3.;		phi_c[1][3] = 0.;			phi_c[2][3] = 0.;			phi_c[3][3] = 1. /4.;
				}
				else if(isurf == 2){
					lCoord[0] = xc1;
					lCoord[1] = xm2;
					lCoord[2] = xc2;
					lCoord[3] = xcc;
					phi_c[0][0] = 1. /3.;		phi_c[1][0] = 1. /2.;		phi_c[2][0] = 1. /3.;		phi_c[3][0] = 1. /4.;
					phi_c[0][1] = 1. /3.;		phi_c[1][1] = 0.;			phi_c[2][1] = 0.;			phi_c[3][1] = 1. /4.;
					phi_c[0][2] = 1. /3.;		phi_c[1][2] = 1. /2.;		phi_c[2][2] = 1. /3.;		phi_c[3][2] = 1. /4.;
					phi_c[0][3] = 0.;			phi_c[1][3] = 0.;			phi_c[2][3] = 1. /3.;		phi_c[3][3] = 1. /4.;
				}
				else if(isurf == 3){
					lCoord[0] = xc2;
					lCoord[1] = xm3;
					lCoord[2] = xc3;
					lCoord[3] = xcc;
					phi_c[0][0] = 1. /3.;		phi_c[1][0] = 1. /2.;		phi_c[2][0] = 1. /3.;		phi_c[3][0] = 1. /4.;
					phi_c[0][1] = 0.;			phi_c[1][1] = 0.;			phi_c[2][1] = 1. /3.;		phi_c[3][1] = 1. /4.;
					phi_c[0][2] = 1. /3.;		phi_c[1][2] = 0.;			phi_c[2][2] = 0.;			phi_c[3][2] = 1. /4.;
					phi_c[0][3] = 1. /3.;		phi_c[1][3] = 1. /2.;		phi_c[2][3] = 1. /3.;		phi_c[3][3] = 1. /4.;
				}
				
				shapeAtSurf(lCoord, node1, phi_c, sfemData->nodeArNode[n], nd, sfemData->nodeVolume[n], t4g, G, GG);
			}

			if(!boundary){
				continue;
			}
			int bnodes[3];
			int boundaryNumBer;
			for(int isurf = 1; isurf <= 3; isurf++){
				if(isurf == 1){
					bnodes[0] = nd[0];
					bnodes[1] = nd[2];
					bnodes[2] = nd[3];

					boundaryNumBer = findFirstIndexOf(nodes, nd[1]);
				}
				else if(isurf == 2){
					bnodes[0] = nd[0];
					bnodes[1] = nd[1];
					bnodes[2] = nd[2];

					boundaryNumBer = findFirstIndexOf(nodes, nd[3]);
				}
				else if(isurf == 3){
					bnodes[0] = nd[0];
					bnodes[1] = nd[1];
					bnodes[2] = nd[3];

					boundaryNumBer = findFirstIndexOf(nodes, nd[2]);
				}

				if(boundaryNumBer == 0)			boundaryNumBer = 2;
				else if(boundaryNumBer == 1)	boundaryNumBer = 3;
				else if(boundaryNumBer == 2)	boundaryNumBer = 1;
				else 							boundaryNumBer = 0;

				if(sfemData->elementFaceFlag[element][boundaryNumBer]){
					if(isurf == 1){
						lCoord[0] = xm2;
						lCoord[1] = xc2;
						lCoord[2] = xm3;
						lCoord[3] = node1;
						phi_c[0][0] = 1. /2.;		phi_c[1][0] = 1. /3.;		phi_c[2][0] = 1. /2.;		phi_c[3][0] = 1.;
						phi_c[0][1] = 0.;			phi_c[1][1] = 0.;			phi_c[2][1] = 0.;			phi_c[3][1] = 0.;
						phi_c[0][2] = 1. /2.;		phi_c[1][2] = 1. /3.;		phi_c[2][2] = 0.;			phi_c[3][2] = 0.;
						phi_c[0][3] = 0.;			phi_c[1][3] = 1. /3.;		phi_c[2][3] = 1. /2.;		phi_c[3][3] = 0.;
					}
					else if(isurf == 2){
						lCoord[0] = xm1;
						lCoord[1] = xc1;
						lCoord[2] = xm2;
						lCoord[3] = node1;
						phi_c[0][0] = 1. /2.;		phi_c[1][0] = 1. /3.;		phi_c[2][0] = 1. /2.;		phi_c[3][0] = 1.;
						phi_c[0][1] = 1. /2.;		phi_c[1][1] = 1. /3.;		phi_c[2][1] = 0.;			phi_c[3][1] = 0.;
						phi_c[0][2] = 0.;			phi_c[1][2] = 1. /3.;		phi_c[2][2] = 1. /2.;		phi_c[3][2] = 0.;
						phi_c[0][3] = 0.;			phi_c[1][3] = 0.;			phi_c[2][3] = 0.;			phi_c[3][3] = 0.;
					}
					else if(isurf == 3){
						lCoord[0] = xm3;
						lCoord[1] = xc3;
						lCoord[2] = xm1;
						lCoord[3] = node1;
						phi_c[0][0] = 1. /2.;		phi_c[1][0] = 1. /3.;		phi_c[2][0] = 1. /2.;		phi_c[3][0] = 1.;
						phi_c[0][1] = 0.;			phi_c[1][1] = 1. /3.;		phi_c[2][1] = 1. /2.;		phi_c[3][1] = 0.;
						phi_c[0][2] = 0.;			phi_c[1][2] = 0.;			phi_c[2][2] = 0.;			phi_c[3][2] = 0.;
						phi_c[0][3] = 1. /2.;		phi_c[1][3] = 1. /3.;		phi_c[2][3] = 0.;			phi_c[3][3] = 0.;	
					}
					shapeAtSurf(lCoord, xcc, phi_c, sfemData->nodeArNode[n], nd, sfemData->nodeVolume[n], t4g, G, GG);
				}
			}

		}

		if(GG){
			int node = (int) sfemData->nodeArNode[n].size();

			double a = 3 * sfemData->nodeVolume[n] / 4. / M_PI;

			double l = pow(a, 1. / 3.);

			for(int i = 0;i < node; i++){
				GG[i] *= l;
			}
		}	
	}

	void SFEM_C3D4::T4G(Point* t4g, Point& x1, Point& x2, Point& x3, Point& x4)
	{
		double J[3][3] = { 0. };
		double Ji[3][3] = { 0. };

		vector<double> Gr = { -1.0, 1.0, 0.0, 0.0 };
		vector<double> Gs = { -1.0, 0.0, 1.0, 0.0 };
		vector<double> Gt = { -1.0, 0.0, 0.0, 1.0 };
		//计算雅克比矩阵
		J[0][0] = Gr[0] * x1.x + Gr[1] * x2.x + Gr[2] * x3.x + Gr[3] * x4.x;
		J[0][1] = Gs[0] * x1.x + Gs[1] * x2.x + Gs[2] * x3.x + Gs[3] * x4.x;
		J[0][2] = Gt[0] * x1.x + Gt[1] * x2.x + Gt[2] * x3.x + Gt[3] * x4.x;
		J[1][0] = Gr[0] * x1.y + Gr[1] * x2.y + Gr[2] * x3.y + Gr[3] * x4.y;
		J[1][1] = Gs[0] * x1.y + Gs[1] * x2.y + Gs[2] * x3.y + Gs[3] * x4.y;
		J[1][2] = Gt[0] * x1.y + Gt[1] * x2.y + Gt[2] * x3.y + Gt[3] * x4.y;
		J[2][0] = Gr[0] * x1.z + Gr[1] * x2.z + Gr[2] * x3.z + Gr[3] * x4.z;
		J[2][1] = Gs[0] * x1.z + Gs[1] * x2.z + Gs[2] * x3.z + Gs[3] * x4.z;
		J[2][2] = Gt[0] * x1.z + Gt[1] * x2.z + Gt[2] * x3.z + Gt[3] * x4.z;

		double detJ = J[0][0] * (J[1][1] * J[2][2] - J[1][2] * J[2][1])
					+ J[0][1] * (J[1][2] * J[2][0] - J[2][2] * J[1][0])
					+ J[0][2] * (J[1][0] * J[2][1] - J[1][1] * J[2][0]);
		
		double deti = 1.0 / detJ;

		Ji[0][0] = deti * (J[1][1] * J[2][2] - J[1][2] * J[2][1]);
		Ji[1][0] = deti * (J[1][2] * J[2][0] - J[1][0] * J[2][2]);
		Ji[2][0] = deti * (J[1][0] * J[2][1] - J[1][1] * J[2][0]);

		Ji[0][1] = deti * (J[0][2] * J[2][1] - J[0][1] * J[2][2]);
		Ji[1][1] = deti * (J[0][0] * J[2][2] - J[0][2] * J[2][0]);
		Ji[2][1] = deti * (J[0][1] * J[2][0] - J[0][0] * J[2][1]);

		Ji[0][2] = deti * (J[0][1] * J[1][2] - J[1][1] * J[0][2]);
		Ji[1][2] = deti * (J[0][2] * J[1][0] - J[0][0] * J[1][2]);
		Ji[2][2] = deti * (J[0][0] * J[1][1] - J[0][1] * J[1][0]);

		for (int i = 0; i < 4; i++) {
			t4g[i].x = Ji[0][0] * Gr[i] + Ji[1][0] * Gs[i] + Ji[2][0] * Gt[i];
			t4g[i].y = Ji[0][1] * Gr[i] + Ji[1][1] * Gs[i] + Ji[2][1] * Gt[i];
			t4g[i].z = Ji[0][2] * Gr[i] + Ji[1][2] * Gs[i] + Ji[2][2] * Gt[i];
		}


	}
	
	void SFEM_C3D4::shapeAtSurf(vector<Point>& lCoord, Point& xcc, double phi_c[4][4],
	 							const vector<int>& nodeArray, int n[4], double vol, 
								const Point* t4g, Point* sg, Mat3* gg)
	{
		//确定法向矢量
		Point v12 = lCoord[1] - lCoord[0];
		Point v13 = lCoord[2] - lCoord[0];
		Point v14 = lCoord[3] - lCoord[0];

		Point norm = v12 ^ v13;
		norm.unit();

		Point vec = (lCoord[0] + lCoord[1] + lCoord[2] + lCoord[3]) / 4 - xcc;
		if(vec * norm< 0){
			norm.negate();
		}

		//求解该表面的局部坐标
		double d12 = v12.norm();
		double d13 = v13.norm();
		double d14 = v14.norm();
		double a13 = acos(v12 * v13 / d12 / d13);
		double a14 = acos(v12 * v14 / d12 / d14);

		double llcoords[4][2] = { 0. };
		llcoords[0][0] = 0.;
		llcoords[0][1] = 0.;
		llcoords[1][0] = d12;
		llcoords[1][1] = 0.;
		llcoords[2][0] = d13 * cos(a13);
		llcoords[2][1] = d13 * sin(a13);
		llcoords[3][0] = d14 * cos(a14);
		llcoords[3][1] = d14 * sin(a14);

		double J[2][2];
		double detJ;
		double phi[4];
		for(int i = 0; i < 4; ++i){
	
	
	

			J[0][0] = J[0][1] = J[1][0] = J[1][1] = 0.;
			for(int j = 0; j < 4; j++){
				J[0][0] += Gr[i][j] * llcoords[j][0];
				J[0][1] += Gr[i][j] * llcoords[j][1];
				J[1][0] += Gs[i][j] * llcoords[j][0];
				J[1][1] += Gs[i][j] * llcoords[j][1];

				phi[j] = H[i][0] * phi_c[0][j] + H[i][1] * phi_c[1][j] + H[i][2] * phi_c[2][j] + H[i][3] * phi_c[3][j];
			}
			detJ = J[0][0] * J[1][1] - J[0][1] * J[1][0];

			for(int j = 0; j < 4; j++){
				int index = findFirstIndexOf(nodeArray, n[j]);

				sg[index].x += phi[j] * norm.x * detJ / vol;
				sg[index].y += phi[j] * norm.y * detJ / vol;
				sg[index].z += phi[j] * norm.z * detJ / vol;
				if(gg){
					gg[index].d[0][0] += t4g[j].x * norm.x * detJ / vol;
					gg[index].d[0][1] += t4g[j].x * norm.y * detJ / vol;
					gg[index].d[0][2] += t4g[j].x * norm.z * detJ / vol;
					gg[index].d[1][0] += t4g[j].y * norm.x * detJ / vol;
					gg[index].d[1][1] += t4g[j].y * norm.y * detJ / vol;
					gg[index].d[1][2] += t4g[j].y * norm.z * detJ / vol;
					gg[index].d[2][0] += t4g[j].z * norm.x * detJ / vol;
					gg[index].d[2][1] += t4g[j].z * norm.y * detJ / vol;
					gg[index].d[2][2] += t4g[j].z * norm.z * detJ / vol;
				}
			}
		}
	}

	int SFEM_C3D4::findFirstIndexOf(const vector<int>& nodeArray, int n)
	{
		auto it = std::find(nodeArray.begin(), nodeArray.end(), n);
		if(it == nodeArray.end()){
			return -1;
		}
		else{
			return (int)(it - nodeArray.begin());
		}		
	}
}