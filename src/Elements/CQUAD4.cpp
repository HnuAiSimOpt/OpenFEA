/**************************************************************************

Copyright:  WH team

Author: YinJichao <jichaoyinyjc@163.com>、ChenBinXiang <chen_mech99@163.com>

Completion date:  XXX

Description: XXX

**************************************************************************/

#pragma once

#include "include/CQUAD4.h"
namespace CAE
{
	REGISTER(ele_base, CQUAD4, "CQUAD4");

	int CQUAD4::IDI[4][9] = { { 0,1,2,3,4,5,12,13,14}, {3,4,5,6,7,8,12,13,14 }, {6,7,8,9,10,11,12,13,14}, {9,10,11,0,1,2,12,13,14} };

	int CQUAD4::IDM[12] = { 2,3,4,8,9,10,14,15,16,20,21,22 };

	int CQUAD4::IDV[9] = { 0,3,6,1,4,7,2,5,8 };

	int CQUAD4::kmIndex[8] = { 0,1,6,7,12,13,18,19 }; //HBAR = 0时膜在总刚矩阵中的索引

	// 建立本构矩阵
	void CQUAD4::build_cons_mat(double D[8][8])
	{
		double e = matrial_struc_.young_modulus;
		double nu = matrial_struc_.poisson_ratio;
		double temp = e / (1. - nu * nu);

		double G1[3][3] = { 0. };
		G1[0][0] = temp;
		G1[0][1] = nu * temp;
		G1[1][0] = nu * temp;
		G1[1][1] = temp;
		G1[2][2] = (1. - nu) / 2. * temp;

		double G3[2][2] = { 0. };
		G3[0][0] = 0.5 * e / (1 + nu);
		G3[1][1] = 0.5 * e / (1 + nu);

		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				D[i][j] = G1[i][j] * t;
				D[i + 3][j + 3] = G1[i][j] * t * t * t / 12.;
				if (i == 2 || j == 2) {
					continue;
				}
				D[i + 6][j + 6] = G3[i][j] * t * 5. / 6.;
			}
		}
	}

	//建立单元局部坐标系
	void CQUAD4::computeLocalSystems(Eigen::Ref<Eigen::MatrixXd> node_coords)
	{
		Point n1 = Point(node_coords(0, 0), node_coords(0, 1), node_coords(0, 2));
		Point n2 = Point(node_coords(1, 0), node_coords(1, 1), node_coords(1, 2));
		Point n3 = Point(node_coords(2, 0), node_coords(2, 1), node_coords(2, 2));
		Point n4 = Point(node_coords(3, 0), node_coords(3, 1), node_coords(3, 2));

		Point V12 = n2 - n1;
		Point V13 = n3 - n1;
		Point V24 = n4 - n2;

		Point KVEC;

		KVEC.x = V13.y * V24.z - V13.z * V24.y;
		KVEC.y = V13.z * V24.x - V13.x * V24.z;
		KVEC.z = V13.x * V24.y - V13.y * V24.x;
		KVEC = KVEC.normalized();

		HBAR = ((V12.x * KVEC.x + V12.y * KVEC.y + V12.z * KVEC.z));
		Point IVEC = V12 - KVEC * 2.0 * HBAR;
		IVEC = IVEC.normalized();

		Point JVEC;
		JVEC.x = KVEC.y * IVEC.z - KVEC.z * IVEC.y;
		JVEC.y = KVEC.z * IVEC.x - KVEC.x * IVEC.z;
		JVEC.z = KVEC.x * IVEC.y - KVEC.y * IVEC.x;
		JVEC = JVEC.normalized();

		double TE_12[3][3];
		TE_12[0][0] = IVEC.x;      TE_12[0][1] = IVEC.y;      TE_12[0][2] = IVEC.z;
		TE_12[1][0] = JVEC.x;      TE_12[1][1] = JVEC.y;      TE_12[1][2] = JVEC.z;
		TE_12[2][0] = KVEC.x;      TE_12[2][1] = KVEC.y;      TE_12[2][2] = KVEC.z;

		double THETA = acos((V13.x * IVEC.x + V13.y * IVEC.y + V13.z * IVEC.z) / V13.norm());
		double GAMMA = acos((-V24.x * IVEC.x - V24.y * IVEC.y - V24.z * IVEC.z) / V24.norm());
		double DELTA = (THETA - GAMMA) / 2.0;

		double CT_QD[3][3] = { 0. };
		CT_QD[0][0] = cos(DELTA);      CT_QD[0][1] = sin(DELTA);
		CT_QD[1][0] = -sin(DELTA);     CT_QD[1][1] = cos(DELTA);
		CT_QD[2][2] = 1.0;

		
		TEG[0][0] = CT_QD[0][0] * TE_12[0][0] + CT_QD[0][1] * TE_12[1][0] + CT_QD[0][2] * TE_12[2][0];
		TEG[0][1] = CT_QD[0][0] * TE_12[0][1] + CT_QD[0][1] * TE_12[1][1] + CT_QD[0][2] * TE_12[2][1];
		TEG[0][2] = CT_QD[0][0] * TE_12[0][2] + CT_QD[0][1] * TE_12[1][2] + CT_QD[0][2] * TE_12[2][2];
		TEG[1][0] = CT_QD[1][0] * TE_12[0][0] + CT_QD[1][1] * TE_12[1][0] + CT_QD[1][2] * TE_12[2][0];
		TEG[1][1] = CT_QD[1][0] * TE_12[0][1] + CT_QD[1][1] * TE_12[1][1] + CT_QD[1][2] * TE_12[2][1];
		TEG[1][2] = CT_QD[1][0] * TE_12[0][2] + CT_QD[1][1] * TE_12[1][2] + CT_QD[1][2] * TE_12[2][2];
		TEG[2][0] = CT_QD[2][0] * TE_12[0][0] + CT_QD[2][1] * TE_12[1][0] + CT_QD[2][2] * TE_12[2][0];
		TEG[2][1] = CT_QD[2][0] * TE_12[0][1] + CT_QD[2][1] * TE_12[1][1] + CT_QD[2][2] * TE_12[2][1];
		TEG[2][2] = CT_QD[2][0] * TE_12[0][2] + CT_QD[2][1] * TE_12[1][2] + CT_QD[2][2] * TE_12[2][2];


		lcoords.push_back(TEG * (n1 - n1));
		lcoords.push_back(TEG * (n2 - n1));
		lcoords.push_back(TEG * (n3 - n1));
		lcoords.push_back(TEG * (n4 - n1));

		XSD[0] = lcoords[0].x - lcoords[1].x;
		XSD[1] = lcoords[1].x - lcoords[2].x;
		XSD[2] = lcoords[2].x - lcoords[3].x;
		XSD[3] = lcoords[3].x - lcoords[0].x;

		YSD[0] = lcoords[0].y - lcoords[1].y;
		YSD[1] = lcoords[1].y - lcoords[2].y;
		YSD[2] = lcoords[2].y - lcoords[3].y;
		YSD[3] = lcoords[3].y - lcoords[0].y;

		double NUM = lcoords[3].y * (lcoords[1].x - lcoords[3].x) - lcoords[3].x * (lcoords[1].y - lcoords[3].y);
		double DEN = lcoords[2].y * (lcoords[1].x - lcoords[3].x) - lcoords[2].x * (lcoords[1].y - lcoords[3].y);
		double RAT = NUM / DEN;
		double x = lcoords[2].x * RAT;
		double y = lcoords[2].y * RAT;

		lcoords.push_back(Point(x, y, 0.0));//新增中心虚拟节点

		double area = 0.5 * ((lcoords[1].x - lcoords[3].x) * (lcoords[2].y - lcoords[0].y) -
			(lcoords[1].y - lcoords[3].y) * (lcoords[2].x - lcoords[0].x));
	}

	// 建立单元刚度矩阵
	void CQUAD4::build_ele_stiff_mat(Eigen::Ref<Eigen::MatrixXd> node_coords, Eigen::Ref<Eigen::MatrixXd> stiffness_matrix)
	{

		double matD[8][8] = { 0. };
		build_cons_mat(matD);

		computeLocalSystems(node_coords);

		double ke[24][24] = { 0. };
		KmPart(node_coords, matD, ke);//膜部分：四高斯点积分

		KpPart(matD, ke);//板部分：四子三角形分别计算

		TEG.trans(ke);
	}

	void CQUAD4::KmPart(Eigen::Ref<Eigen::MatrixXd> node_coords, double D[8][8], double ke[24][24])
	{
		Point n1 = Point(node_coords(0, 0), node_coords(0, 1), node_coords(0, 2));
		Point n2 = Point(node_coords(1, 0), node_coords(1, 1), node_coords(1, 2));
		Point n3 = Point(node_coords(2, 0), node_coords(2, 1), node_coords(2, 2));
		Point n4 = Point(node_coords(3, 0), node_coords(3, 1), node_coords(3, 2));

		double  DETJ[4] = { 0. };
		double Bm[4][3][8] = { 0. };

		for (int n = 0; n < 4; ++n)
		{

			double _xi = 0.0, _eta = 0.0, _w = 1 / sqrt(3);
			if (n == 0) { _xi = -_w;	_eta = -_w; }
			if (n == 1) { _xi = -_w;	_eta = _w; }
			if (n == 2) { _xi = _w;		_eta = -_w; }
			if (n == 3) { _xi = _w;		_eta = _w; }

			//计算形函数
			//double N[4]; 
			//N[0] = 0.25 * (1 - _xi) * (1 - _eta);
			//N[1] = 0.25 * (1 + _xi) * (1 - _eta);
			//N[2] = 0.25 * (1 + _xi) * (1 + _eta);
			//N[3] = 0.25 * (1 - _xi) * (1 + _eta);
			double dNdxi[4], dNdeta[4];
			dNdxi[0] = -0.25 * (1 - _eta);
			dNdxi[1] = 0.25 * (1 - _eta);
			dNdxi[2] = 0.25 * (1 + _eta);
			dNdxi[3] = -0.25 * (1 + _eta);

			dNdeta[0] = -0.25 * (1 - _xi);
			dNdeta[1] = -0.25 * (1 + _xi);
			dNdeta[2] = 0.25 * (1 + _xi);
			dNdeta[3] = 0.25 * (1 - _xi);

			//计算雅克比矩阵
			double J[2][2], InvJ[2][2];
			J[0][0] = (-(1.0 - _eta) * XSD[0] + (1.0 + _eta) * XSD[2]) / 4.0;
			J[0][1] = (-(1.0 - _eta) * YSD[0] + (1.0 + _eta) * YSD[2]) / 4.0;
			J[1][0] = ((1.0 - _xi) * XSD[3] - (1.0 + _xi) * XSD[1]) / 4.0;
			J[1][1] = ((1.0 - _xi) * YSD[3] - (1.0 + _xi) * YSD[1]) / 4.0;

			DETJ[n] = J[0][0] * J[1][1] - J[0][1] * J[1][0];

			InvJ[0][0] = J[1][1] / DETJ[n];		InvJ[0][1] = -J[0][1] / DETJ[n];
			InvJ[1][0] = -J[1][0] / DETJ[n];	InvJ[1][1] = J[0][0] / DETJ[n];

			Bm[n][0][0] = InvJ[0][0] * dNdxi[0] + InvJ[0][1] * dNdeta[0];
			Bm[n][0][2] = InvJ[0][0] * dNdxi[1] + InvJ[0][1] * dNdeta[1];
			Bm[n][0][4] = InvJ[0][0] * dNdxi[2] + InvJ[0][1] * dNdeta[2];
			Bm[n][0][6] = InvJ[0][0] * dNdxi[3] + InvJ[0][1] * dNdeta[3];
			Bm[n][1][1] = InvJ[1][0] * dNdxi[0] + InvJ[1][1] * dNdeta[0];
			Bm[n][1][3] = InvJ[1][0] * dNdxi[1] + InvJ[1][1] * dNdeta[1];
			Bm[n][1][5] = InvJ[1][0] * dNdxi[2] + InvJ[1][1] * dNdeta[2];
			Bm[n][1][7] = InvJ[1][0] * dNdxi[3] + InvJ[1][1] * dNdeta[3];
			Bm[n][2][0] = InvJ[1][0] * dNdxi[0] + InvJ[1][1] * dNdeta[0];
			Bm[n][2][1] = InvJ[0][0] * dNdxi[0] + InvJ[0][1] * dNdeta[0];
			Bm[n][2][2] = InvJ[1][0] * dNdxi[1] + InvJ[1][1] * dNdeta[1];
			Bm[n][2][3] = InvJ[0][0] * dNdxi[1] + InvJ[0][1] * dNdeta[1];
			Bm[n][2][4] = InvJ[1][0] * dNdxi[2] + InvJ[1][1] * dNdeta[2];
			Bm[n][2][5] = InvJ[0][0] * dNdxi[2] + InvJ[0][1] * dNdeta[2];
			Bm[n][2][6] = InvJ[1][0] * dNdxi[3] + InvJ[1][1] * dNdeta[3];
			Bm[n][2][7] = InvJ[0][0] * dNdxi[3] + InvJ[0][1] * dNdeta[3];
		}
		double SUMD = DETJ[0] + DETJ[1] + DETJ[2] + DETJ[3];
		for (int i = 0; i < 8; i++) {
			double SUMB = 0.0;
			for (int K = 0; K < 4; K++) {
				SUMB += DETJ[K] * Bm[K][2][i];
			}
			for (int K = 0; K < 4; K++) {
				Bm[K][2][i] = SUMB / SUMD;
			}
		}

		if (abs(HBAR) > 0.0) {

			Point V12B = n2 - n1;
			Point V13B = n3 - n1;
			Point V24B = n4 - n2;

			Point IVEC, JVEC, KVEC;
			KVEC.x = V13B.y * V24B.z - V13B.z * V24B.y;
			KVEC.y = V13B.z * V24B.x - V13B.x * V24B.z;
			KVEC.z = V13B.x * V24B.y - V13B.y * V24B.x;
			KVEC = KVEC.normalized();
			IVEC = V12B - KVEC * 2.0 * HBAR;
			IVEC = IVEC.normalized();
			JVEC.x = KVEC.y * IVEC.z - KVEC.z * IVEC.y;
			JVEC.y = KVEC.z * IVEC.x - KVEC.x * IVEC.z;
			JVEC.z = KVEC.x * IVEC.y - KVEC.y * IVEC.x;
			JVEC = JVEC.normalized();

			double X12 = -(V12B.x * IVEC.x + V12B.y * IVEC.y + V12B.z * IVEC.z);
			double X13 = -(V13B.x * IVEC.x + V13B.y * IVEC.y + V13B.z * IVEC.z);
			double X24 = -(V24B.x * IVEC.x + V24B.y * IVEC.y + V24B.z * IVEC.z);
			double X14 = X12 + X24;
			double X23 = X13 - X12;
			double X34 = X14 - X13;
			double X3 = -X13;
			double X4 = -X14;
			double Y3 = (V13B.x * JVEC.x + V13B.y * JVEC.y + V13B.z * JVEC.z);
			double Y4 = (V24B.x * JVEC.x + V24B.y * JVEC.y + V24B.z * JVEC.z);
			double Y34 = Y3 - Y4;
			double L12 = abs(X12);
			double L23 = sqrt(X23 * X23 + Y3 * Y3);
			double L34 = sqrt(X34 * X34 + Y34 * Y34);
			double L41 = sqrt(X14 * X14 + Y4 * Y4);

			double SIN_TH1 = Y4 / L41;
			double COS_TH1 = -X14 / L41;

			double SIN_TH2 = Y3 / L23;
			double COS_TH2 = X23 / L23;

			double SIN_GAM = (Y4 - Y3) / L34;
			double COS_GAM = (X3 - X4) / L34;

			double CTN_TH1 = COS_TH1 / SIN_TH1;
			double CTN_TH2 = COS_TH2 / SIN_TH2;

			double DELTA1 = SIN_TH2 * COS_GAM - COS_TH2 * SIN_GAM;
			double DELTA2 = SIN_TH1 * COS_GAM + COS_TH1 * SIN_GAM;

			double BMEANT[8][12] = { 0. };
			BMEANT[0][0] = 1.0;
			BMEANT[1][1] = 1.0;
			BMEANT[2][3] = 1.0;
			BMEANT[3][4] = 1.0;
			BMEANT[4][6] = 1.0;
			BMEANT[5][7] = 1.0;
			BMEANT[6][9] = 1.0;
			BMEANT[7][1] = 1.0;
			BMEANT[0][2] = HBAR / L12;
			BMEANT[1][2] = -HBAR * (CTN_TH1 / L12 - 1.0 / (L41 * SIN_TH1));
			BMEANT[2][2] = -BMEANT[0][2];
			BMEANT[3][2] = -HBAR * CTN_TH2 / L12;
			BMEANT[6][2] = -HBAR * SIN_GAM / (L41 * DELTA2);
			BMEANT[7][2] = -HBAR * COS_GAM / (L41 * DELTA2);
			BMEANT[0][5] = -BMEANT[0][2];
			BMEANT[1][5] = HBAR * CTN_TH1 / L12;
			BMEANT[2][5] = BMEANT[0][2];
			BMEANT[3][5] = HBAR * (CTN_TH2 / L12 - 1.0 / (L23 * SIN_TH2));
			BMEANT[4][5] = HBAR * SIN_GAM / (L23 * DELTA1);
			BMEANT[5][5] = HBAR * COS_GAM / (L23 * DELTA1);
			BMEANT[3][8] = HBAR / (L23 * SIN_TH2);
			BMEANT[4][8] = -HBAR * (SIN_GAM / L23 + SIN_TH2 / L34) / DELTA1;
			BMEANT[5][8] = -HBAR * (COS_GAM / L23 + COS_TH2 / L34) / DELTA1;
			BMEANT[6][8] = HBAR * SIN_TH1 / (L34 * DELTA2);
			BMEANT[7][8] = -HBAR * COS_TH1 / (L34 * DELTA2);
			BMEANT[1][11] = -HBAR / (L41 * SIN_TH1);
			BMEANT[4][11] = HBAR * SIN_TH2 / (L34 * DELTA1);
			BMEANT[5][11] = HBAR * COS_TH2 / (L34 * DELTA1);
			BMEANT[6][11] = -HBAR * (SIN_TH1 / L34 - SIN_GAM / L41) / DELTA2;
			BMEANT[7][11] = HBAR * (COS_TH1 / L34 + COS_GAM / L41) / DELTA2;

			for (int K = 0; K < 4; K++) {
				double B[3][12] = { 0. };
				B[0][0] = Bm[K][0][0];
				B[0][2] = Bm[K][0][0] * BMEANT[0][2] + Bm[K][0][2] * BMEANT[2][2] + Bm[K][0][6] * BMEANT[6][2];
				B[0][3] = Bm[K][0][2];
				B[0][5] = Bm[K][0][0] * BMEANT[0][5] + Bm[K][0][2] * BMEANT[2][5] + Bm[K][0][4] * BMEANT[4][5];
				B[0][6] = Bm[K][0][4];
				B[0][8] = Bm[K][0][4] * BMEANT[4][8] + Bm[K][0][6] * BMEANT[6][8];
				B[0][9] = Bm[K][0][6];
				B[0][11] = Bm[K][0][4] * BMEANT[4][11] + Bm[K][0][6] * BMEANT[6][11];

				B[1][1] = Bm[K][1][1];
				B[1][2] = Bm[K][1][1] * BMEANT[1][2] + Bm[K][1][3] * BMEANT[3][2] + Bm[K][1][7] * BMEANT[7][2];
				B[1][4] = Bm[K][1][3];
				B[1][5] = Bm[K][1][1] * BMEANT[1][5] + Bm[K][1][3] * BMEANT[3][5] + Bm[K][1][5] * BMEANT[5][5];
				B[1][7] = Bm[K][1][5];
				B[1][8] = Bm[K][1][3] * BMEANT[3][8] + Bm[K][1][5] * BMEANT[5][8] + Bm[K][1][7] * BMEANT[7][8];
				B[1][10] = Bm[K][1][7];
				B[1][11] = Bm[K][1][1] * BMEANT[1][11] + Bm[K][1][5] * BMEANT[5][11] + Bm[K][1][7] * BMEANT[7][11];

				B[2][0] = Bm[K][2][0];
				B[2][1] = Bm[K][2][1];
				B[2][2] = Bm[K][2][0] * BMEANT[0][2] + Bm[K][2][1] * BMEANT[1][2] + Bm[K][2][2] * BMEANT[2][2] + Bm[K][2][3] * BMEANT[3][2] + Bm[K][2][6] * BMEANT[6][2] + Bm[K][2][7] * BMEANT[7][2];
				B[2][3] = Bm[K][2][2];
				B[2][4] = Bm[K][2][3];
				B[2][5] = Bm[K][2][0] * BMEANT[0][5] + Bm[K][2][1] * BMEANT[1][5] + Bm[K][2][2] * BMEANT[2][5] + Bm[K][2][3] * BMEANT[3][5] + Bm[K][2][4] * BMEANT[4][5] + Bm[K][2][5] * BMEANT[5][5];
				B[2][6] = Bm[K][2][4];
				B[2][7] = Bm[K][2][5];
				B[2][8] = Bm[K][2][3] * BMEANT[3][8] + Bm[K][2][4] * BMEANT[4][8] + Bm[K][2][5] * BMEANT[5][8] + Bm[K][2][6] * BMEANT[6][8] + Bm[K][2][7] * BMEANT[7][8];
				B[2][9] = Bm[K][2][6];
				B[2][10] = Bm[K][2][7];
				B[2][11] = Bm[K][2][1] * BMEANT[1][11] + Bm[K][2][4] * BMEANT[4][11] + Bm[K][2][5] * BMEANT[5][11] + Bm[K][2][6] * BMEANT[6][11] + Bm[K][2][7] * BMEANT[7][11];

				double DB[3][12] = { 0. };
				DB[0][0] = D[0][0] * B[0][0] + D[0][1] * B[1][0] + D[0][2] * B[2][0];
				DB[0][1] = D[0][0] * B[0][1] + D[0][1] * B[1][1] + D[0][2] * B[2][1];
				DB[0][2] = D[0][0] * B[0][2] + D[0][1] * B[1][2] + D[0][2] * B[2][2];
				DB[0][3] = D[0][0] * B[0][3] + D[0][1] * B[1][3] + D[0][2] * B[2][3];
				DB[0][4] = D[0][0] * B[0][4] + D[0][1] * B[1][4] + D[0][2] * B[2][4];
				DB[0][5] = D[0][0] * B[0][5] + D[0][1] * B[1][5] + D[0][2] * B[2][5];
				DB[0][6] = D[0][0] * B[0][6] + D[0][1] * B[1][6] + D[0][2] * B[2][6];
				DB[0][7] = D[0][0] * B[0][7] + D[0][1] * B[1][7] + D[0][2] * B[2][7];
				DB[0][8] = D[0][0] * B[0][8] + D[0][1] * B[1][8] + D[0][2] * B[2][8];
				DB[0][9] = D[0][0] * B[0][9] + D[0][1] * B[1][9] + D[0][2] * B[2][9];
				DB[0][10] = D[0][0] * B[0][10] + D[0][1] * B[1][10] + D[0][2] * B[2][10];
				DB[0][11] = D[0][0] * B[0][11] + D[0][1] * B[1][11] + D[0][2] * B[2][11];

				DB[1][0] = D[1][0] * B[0][0] + D[1][1] * B[1][0] + D[1][2] * B[2][0];
				DB[1][1] = D[1][0] * B[0][1] + D[1][1] * B[1][1] + D[1][2] * B[2][1];
				DB[1][2] = D[1][0] * B[0][2] + D[1][1] * B[1][2] + D[1][2] * B[2][2];
				DB[1][3] = D[1][0] * B[0][3] + D[1][1] * B[1][3] + D[1][2] * B[2][3];
				DB[1][4] = D[1][0] * B[0][4] + D[1][1] * B[1][4] + D[1][2] * B[2][4];
				DB[1][5] = D[1][0] * B[0][5] + D[1][1] * B[1][5] + D[1][2] * B[2][5];
				DB[1][6] = D[1][0] * B[0][6] + D[1][1] * B[1][6] + D[1][2] * B[2][6];
				DB[1][7] = D[1][0] * B[0][7] + D[1][1] * B[1][7] + D[1][2] * B[2][7];
				DB[1][8] = D[1][0] * B[0][8] + D[1][1] * B[1][8] + D[1][2] * B[2][8];
				DB[1][9] = D[1][0] * B[0][9] + D[1][1] * B[1][9] + D[1][2] * B[2][9];
				DB[1][10] = D[1][0] * B[0][10] + D[1][1] * B[1][10] + D[1][2] * B[2][10];
				DB[1][11] = D[1][0] * B[0][11] + D[1][1] * B[1][11] + D[1][2] * B[2][11];

				DB[2][0] = D[2][0] * B[0][0] + D[2][1] * B[1][0] + D[2][2] * B[2][0];
				DB[2][1] = D[2][0] * B[0][1] + D[2][1] * B[1][1] + D[2][2] * B[2][1];
				DB[2][2] = D[2][0] * B[0][2] + D[2][1] * B[1][2] + D[2][2] * B[2][2];
				DB[2][3] = D[2][0] * B[0][3] + D[2][1] * B[1][3] + D[2][2] * B[2][3];
				DB[2][4] = D[2][0] * B[0][4] + D[2][1] * B[1][4] + D[2][2] * B[2][4];
				DB[2][5] = D[2][0] * B[0][5] + D[2][1] * B[1][5] + D[2][2] * B[2][5];
				DB[2][6] = D[2][0] * B[0][6] + D[2][1] * B[1][6] + D[2][2] * B[2][6];
				DB[2][7] = D[2][0] * B[0][7] + D[2][1] * B[1][7] + D[2][2] * B[2][7];
				DB[2][8] = D[2][0] * B[0][8] + D[2][1] * B[1][8] + D[2][2] * B[2][8];
				DB[2][9] = D[2][0] * B[0][9] + D[2][1] * B[1][9] + D[2][2] * B[2][9];
				DB[2][10] = D[2][0] * B[0][10] + D[2][1] * B[1][10] + D[2][2] * B[2][10];
				DB[2][11] = D[2][0] * B[0][11] + D[2][1] * B[1][11] + D[2][2] * B[2][11];

				for (int i = 0, i6 = 0; i < 4; ++i, i6 += 3) {
					for (int j = 0, j6 = 0; j < 4; ++j, j6 += 3) {
						ke[i6 * 2 + 0][j6 * 2 + 0] += (B[0][i6] * DB[0][j6] + B[1][i6] * DB[1][j6] + B[2][i6] * DB[2][j6]) * DETJ[K];
						ke[i6 * 2 + 0][j6 * 2 + 1] += (B[0][i6] * DB[0][j6 + 1] + B[1][i6] * DB[1][j6 + 1] + B[2][i6] * DB[2][j6 + 1]) * DETJ[K];
						ke[i6 * 2 + 0][j6 * 2 + 2] += (B[0][i6] * DB[0][j6 + 2] + B[1][i6] * DB[1][j6 + 2] + B[2][i6] * DB[2][j6 + 2]) * DETJ[K];

						ke[i6 * 2 + 1][j6 * 2 + 0] += (B[0][i6 + 1] * DB[0][j6] + B[1][i6 + 1] * DB[1][j6] + B[2][i6 + 1] * DB[2][j6]) * DETJ[K];
						ke[i6 * 2 + 1][j6 * 2 + 1] += (B[0][i6 + 1] * DB[0][j6 + 1] + B[1][i6 + 1] * DB[1][j6 + 1] + B[2][i6 + 1] * DB[2][j6 + 1]) * DETJ[K];
						ke[i6 * 2 + 1][j6 * 2 + 2] += (B[0][i6 + 1] * DB[0][j6 + 2] + B[1][i6 + 1] * DB[1][j6 + 2] + B[2][i6 + 1] * DB[2][j6 + 2]) * DETJ[K];

						ke[i6 * 2 + 2][j6 * 2 + 0] += (B[0][i6 + 2] * DB[0][j6] + B[1][i6 + 2] * DB[1][j6] + B[2][i6 + 2] * DB[2][j6]) * DETJ[K];
						ke[i6 * 2 + 2][j6 * 2 + 1] += (B[0][i6 + 2] * DB[0][j6 + 1] + B[1][i6 + 2] * DB[1][j6 + 1] + B[2][i6 + 2] * DB[2][j6 + 1]) * DETJ[K];
						ke[i6 * 2 + 2][j6 * 2 + 2] += (B[0][i6 + 2] * DB[0][j6 + 2] + B[1][i6 + 2] * DB[1][j6 + 2] + B[2][i6 + 2] * DB[2][j6 + 2]) * DETJ[K];
					}
				}
			}
		}
		else {
			for (int K = 0; K < 4; K++) {
				double DB[3][12] = { 0. };
				DB[0][0] = D[0][0] * Bm[K][0][0] + D[0][1] * Bm[K][1][0] + D[0][2] * Bm[K][2][0];
				DB[0][1] = D[0][0] * Bm[K][0][1] + D[0][1] * Bm[K][1][1] + D[0][2] * Bm[K][2][1];
				DB[0][2] = D[0][0] * Bm[K][0][2] + D[0][1] * Bm[K][1][2] + D[0][2] * Bm[K][2][2];
				DB[0][3] = D[0][0] * Bm[K][0][3] + D[0][1] * Bm[K][1][3] + D[0][2] * Bm[K][2][3];
				DB[0][4] = D[0][0] * Bm[K][0][4] + D[0][1] * Bm[K][1][4] + D[0][2] * Bm[K][2][4];
				DB[0][5] = D[0][0] * Bm[K][0][5] + D[0][1] * Bm[K][1][5] + D[0][2] * Bm[K][2][5];
				DB[0][6] = D[0][0] * Bm[K][0][6] + D[0][1] * Bm[K][1][6] + D[0][2] * Bm[K][2][6];
				DB[0][7] = D[0][0] * Bm[K][0][7] + D[0][1] * Bm[K][1][7] + D[0][2] * Bm[K][2][7];
				DB[1][0] = D[1][0] * Bm[K][0][0] + D[1][1] * Bm[K][1][0] + D[1][2] * Bm[K][2][0];
				DB[1][1] = D[1][0] * Bm[K][0][1] + D[1][1] * Bm[K][1][1] + D[1][2] * Bm[K][2][1];
				DB[1][2] = D[1][0] * Bm[K][0][2] + D[1][1] * Bm[K][1][2] + D[1][2] * Bm[K][2][2];
				DB[1][3] = D[1][0] * Bm[K][0][3] + D[1][1] * Bm[K][1][3] + D[1][2] * Bm[K][2][3];
				DB[1][4] = D[1][0] * Bm[K][0][4] + D[1][1] * Bm[K][1][4] + D[1][2] * Bm[K][2][4];
				DB[1][5] = D[1][0] * Bm[K][0][5] + D[1][1] * Bm[K][1][5] + D[1][2] * Bm[K][2][5];
				DB[1][6] = D[1][0] * Bm[K][0][6] + D[1][1] * Bm[K][1][6] + D[1][2] * Bm[K][2][6];
				DB[1][7] = D[1][0] * Bm[K][0][7] + D[1][1] * Bm[K][1][7] + D[1][2] * Bm[K][2][7];
				DB[2][0] = D[2][0] * Bm[K][0][0] + D[2][1] * Bm[K][1][0] + D[2][2] * Bm[K][2][0];
				DB[2][1] = D[2][0] * Bm[K][0][1] + D[2][1] * Bm[K][1][1] + D[2][2] * Bm[K][2][1];
				DB[2][2] = D[2][0] * Bm[K][0][2] + D[2][1] * Bm[K][1][2] + D[2][2] * Bm[K][2][2];
				DB[2][3] = D[2][0] * Bm[K][0][3] + D[2][1] * Bm[K][1][3] + D[2][2] * Bm[K][2][3];
				DB[2][4] = D[2][0] * Bm[K][0][4] + D[2][1] * Bm[K][1][4] + D[2][2] * Bm[K][2][4];
				DB[2][5] = D[2][0] * Bm[K][0][5] + D[2][1] * Bm[K][1][5] + D[2][2] * Bm[K][2][5];
				DB[2][6] = D[2][0] * Bm[K][0][6] + D[2][1] * Bm[K][1][6] + D[2][2] * Bm[K][2][6];
				DB[2][7] = D[2][0] * Bm[K][0][7] + D[2][1] * Bm[K][1][7] + D[2][2] * Bm[K][2][7];

				for (int i = 0; i < 8; i++) {
					for (int j = 0; j < 8; j++) {
						ke[kmIndex[i]][kmIndex[j]] += (Bm[K][0][i] * DB[0][j] + Bm[K][1][i] * DB[1][j] + Bm[K][2][i] * DB[2][j]) * DETJ[K];
					}
				}
			}
		}

		//BE1
		{
			double shp[2][4];
			shp[0][0] = -0.25;	shp[0][1] = 0.25;		shp[0][2] = 0.25;		shp[0][3] = -0.25;
			shp[1][0] = -0.25;	shp[1][1] = -0.25;		shp[1][2] = 0.25;		shp[1][3] = 0.25;

			double J[2][2];
			J[0][0] = (XSD[2] - XSD[0]) / 4.0;		J[0][1] = (YSD[2] - YSD[0]) / 4.0;
			J[1][0] = (XSD[3] - XSD[1]) / 4.0;		J[1][1] = (YSD[3] - YSD[1]) / 4.0;

			double DETJ = J[0][0] * J[1][1] - J[0][1] * J[1][0];

			double InvJ[2][2];
			InvJ[0][0] = J[1][1] / DETJ;		InvJ[0][1] = -J[0][1] / DETJ;
			InvJ[1][0] = -J[1][0] / DETJ;	InvJ[1][1] = J[0][0] / DETJ;

			this->BE1[0][0] = InvJ[0][0] * shp[0][0] + InvJ[0][1] * shp[1][0];
			this->BE1[0][1] = InvJ[0][0] * shp[0][1] + InvJ[0][1] * shp[1][1];
			this->BE1[0][2] = InvJ[0][0] * shp[0][2] + InvJ[0][1] * shp[1][2];
			this->BE1[0][3] = InvJ[0][0] * shp[0][3] + InvJ[0][1] * shp[1][3];
			this->BE1[1][0] = InvJ[1][0] * shp[0][0] + InvJ[1][1] * shp[1][0];
			this->BE1[1][1] = InvJ[1][0] * shp[0][1] + InvJ[1][1] * shp[1][1];
			this->BE1[1][2] = InvJ[1][0] * shp[0][2] + InvJ[1][1] * shp[1][2];
			this->BE1[1][3] = InvJ[1][0] * shp[0][3] + InvJ[1][1] * shp[1][3];
		}
	}

	void CQUAD4::KpPart(double matD[8][8], double k0[24][24])
	{
		double KM_QQ_5[15][15] = { 0. }, B2M_QQ_5[3][15] = { 0. }, B3M_QQ_5[3][15] = { 0. };

		double ALPHA[4], BETA[4], GAMMA[4];//子三角形几何数据

		for (int K = 0; K < 4; K++) {

			//****************************************************************************
			//			ELMGM_TRIA      函数开始: 此函数进行子三角形的几何处理
			//****************************************************************************
			Point V12, V13, V23;
			if (K == 0) {
				V12 = lcoords[1] - lcoords[0];	V13 = lcoords[4] - lcoords[0];	V23 = lcoords[4] - lcoords[1];
			}
			else if (K == 1) {
				V12 = lcoords[2] - lcoords[1];	V13 = lcoords[4] - lcoords[1];	V23 = lcoords[4] - lcoords[2];
			}
			else if (K == 2) {
				V12 = lcoords[3] - lcoords[2];	V13 = lcoords[4] - lcoords[2];	V23 = lcoords[4] - lcoords[3];
			}
			else if (K == 3) {
				V12 = lcoords[0] - lcoords[3];	V13 = lcoords[4] - lcoords[3];	V23 = lcoords[4] - lcoords[0];
			}

			double L12 = sqrt(V12.x * V12.x + V12.y * V12.y);
			double L13 = sqrt(V13.x * V13.x + V13.y * V13.y);
			double L23 = sqrt(V23.x * V23.x + V23.y * V23.y);

			double CTHETA = V12.x / L12;		//cosθ
			double STHETA = V12.y / L12;		//sinθ
			double THETA = atan2(V12.y, V12.x);	//θ

			//这里的D矩阵不是材料矩阵，这么命名是方便和MYSTRAN对标Debug。(此函数中材料矩阵用matD表示)
			double D[9][9] = { 0. };
			D[0][0] = 1;
			D[1][1] = CTHETA;   D[1][2] = STHETA;
			D[2][1] = -STHETA;  D[2][2] = CTHETA;
			D[3][3] = 1;
			D[4][4] = CTHETA;   D[4][5] = STHETA;
			D[5][4] = -STHETA;  D[5][5] = CTHETA;
			D[6][6] = 1;
			D[7][7] = CTHETA;   D[7][8] = STHETA;
			D[8][7] = -STHETA;  D[8][8] = CTHETA;

			double CALPHA = (L13 * L13 + L12 * L12 - L23 * L23) / (2.0 * L13 * L12);//cosα
			double CBETA = (L23 * L23 + L12 * L12 - L23 * L23) / (2.0 * L23 * L12);	//cosβ

			ALPHA[K] = acos(CALPHA);			//α
			BETA[K] = acos(CBETA);				//β
			double X3TL = L13 * CALPHA;
			double Y3TL = L13 * sin(ALPHA[K]);
			double X2TL = L12;
			//****************************************************************************
			//			ELMGM_TRIA      函数结束: 此函数进行子三角形的几何处理
			//****************************************************************************

			AREA_TRIA[K] = X2TL * Y3TL / 2;//三角形面积

			//****************************************************************************
			//			TPLT2			函数开始
			//****************************************************************************
			double A[3] = { X3TL - X2TL, -X3TL, X2TL };
			double B[3] = { -Y3TL, Y3TL, 0. };
			double A4 = 4.0 * AREA_TRIA[K];
			double A42 = A4 * AREA_TRIA[K];




			//计算KB,KS
			double KB[9][9] = { 0. }, KS[9][9] = { 0. };

			caculateKB(A, B, AREA_TRIA[K], A4, A42, matD, KB, KS, PHI_SQ_TRIA[K]);

			//KM_TT = KB + PHI_SQ_TRIA * KS
			double KM_TT[9][9] = { 0. };
			for (int i = 0; i < 9; i++) {
				for (int j = 0; j < 9; j++) {
					KM_TT[IDV[i]][IDV[j]] = KB[i][j] + PHI_SQ_TRIA[K] * KS[i][j];
				}
			}

			//KM_TTD = KM_TT * D
			double KM_TTD[9][9] = { 0. };
			for (int i = 0; i < 9; i++) {
				for (int j = 0; j < 9; j++) {
					KM_TTD[i][j] = KM_TT[i][0] * D[0][j] + KM_TT[i][1] * D[1][j] + KM_TT[i][2] * D[2][j]
						+ KM_TT[i][3] * D[3][j] + KM_TT[i][4] * D[4][j] + KM_TT[i][5] * D[5][j]
						+ KM_TT[i][6] * D[6][j] + KM_TT[i][7] * D[7][j] + KM_TT[i][8] * D[8][j];
				}
			}
			//KM_TQ = D.transpose() * (KM_TT * D)	因为KM_TQ会组装到KM_QQ_5中，这里不创建KM_TQ直接组装了
			for (int i = 0; i < 9; i++) {
				for (int j = 0; j < 9; j++) {
					KM_QQ_5[IDI[K][i]][IDI[K][j]] += D[0][i] * KM_TTD[0][j] + D[1][i] * KM_TTD[1][j] + D[2][i] * KM_TTD[2][j]
						+ D[3][i] * KM_TTD[3][j] + D[4][i] * KM_TTD[4][j] + D[5][i] * KM_TTD[5][j]
						+ D[6][i] * KM_TTD[6][j] + D[7][i] * KM_TTD[7][j] + D[8][i] * KM_TTD[8][j];
				}
			}

			//==========================================开始计算B矩阵==========================================
			double CTH2 = cos(2.0 * THETA);
			double STH2 = sin(2.0 * THETA);

			double TE_STRESS[3][3];
			TE_STRESS[0][0] = 0.5 * (1.0 + CTH2);   TE_STRESS[0][1] = 0.5 * (1.0 - CTH2);   TE_STRESS[0][2] = STH2;
			TE_STRESS[1][0] = 0.5 * (1.0 - CTH2);   TE_STRESS[1][1] = 0.5 * (1.0 + CTH2);   TE_STRESS[1][2] = -STH2;
			TE_STRESS[2][0] = 0.5 * STH2;			TE_STRESS[2][1] = -0.5 * STH2;			TE_STRESS[2][2] = CTH2;

			double TE_SHEAR[3][3];
			TE_SHEAR[0][0] = CTHETA;	TE_SHEAR[0][1] = -STHETA;	TE_SHEAR[0][2] = 0.0;
			TE_SHEAR[1][0] = STHETA;	TE_SHEAR[1][1] = CTHETA;	TE_SHEAR[01][2] = 0.0;
			TE_SHEAR[2][0] = 0.0;		TE_SHEAR[2][1] = 0.0;		TE_SHEAR[2][2] = 1.0;

			double XI = 1.0 / 3.0, _A4 = 4.0 * AREA_TRIA[K];

			double B2V[3][9] = { 0. }, B3V[3][9] = { 0. };
			B2V[0][6] = B[0] / (2.0 * AREA_TRIA[K]);
			B2V[0][7] = B[1] / (2.0 * AREA_TRIA[K]);
			B2V[0][8] = B[2] / (2.0 * AREA_TRIA[K]);
			B2V[1][3] = -A[0] / (2.0 * AREA_TRIA[K]);
			B2V[1][4] = -A[1] / (2.0 * AREA_TRIA[K]);
			B2V[1][5] = -A[2] / (2.0 * AREA_TRIA[K]);
			B2V[2][3] = -B2V[0][6];
			B2V[2][4] = -B2V[0][7];
			B2V[2][5] = -B2V[0][8];
			B2V[2][6] = -B2V[1][3];
			B2V[2][7] = -B2V[1][4];
			B2V[2][8] = -B2V[1][5];

			B3V[0][0] = B[0] / (2.0 * AREA_TRIA[K]);
			B3V[1][0] = A[0] / (2.0 * AREA_TRIA[K]);
			B3V[0][1] = B[1] / (2.0 * AREA_TRIA[K]);
			B3V[1][1] = A[1] / (2.0 * AREA_TRIA[K]);
			B3V[0][2] = B[2] / (2.0 * AREA_TRIA[K]);
			B3V[1][2] = A[2] / (2.0 * AREA_TRIA[K]);
			B3V[0][3] = -(B[0] * (XI * B[2] - XI * B[1]) / A4);
			B3V[0][4] = -(B[1] * (-XI * B[2] + XI * B[0]) / A4);
			B3V[0][5] = -(B[2] * (XI * B[1] - XI * B[0]) / A4);
			B3V[0][6] = (XI * (A4 - B[1] * A[2] + B[2] * A[1]) + B[0] * (-XI * A[2] + XI * A[1])) / A4;
			B3V[0][7] = (XI * (A4 + B[0] * A[2] - B[2] * A[0]) + B[1] * (XI * A[2] - XI * A[0])) / A4;
			B3V[0][8] = (XI * (A4 + B[1] * A[0] - B[0] * A[1]) + B[2] * (-XI * A[1] + XI * A[0])) / A4;
			B3V[1][3] = -(XI * (A4 + A[1] * B[2] - A[2] * B[1]) + A[0] * (XI * B[2] - XI * B[1])) / A4;
			B3V[1][4] = -(XI * (A4 - A[0] * B[2] + A[2] * B[0]) + A[1] * (-XI * B[2] + XI * B[0])) / A4;
			B3V[1][5] = -(XI * (A4 - A[1] * B[0] + A[0] * B[1]) + A[2] * (XI * B[1] - XI * B[0])) / A4;
			B3V[1][6] = A[0] * (-XI * A[2] + XI * A[1]) / A4;
			B3V[1][7] = A[1] * (XI * A[2] - XI * A[0]) / A4;
			B3V[1][8] = A[2] * (-XI * A[1] + XI * A[0]) / A4;

			//============计算当前子三角形的B2M_QQk，并根据组装至B2M_QQ_5
			double B2M_TT[3][9];
			for (int i = 0; i < 9; i++) {
				B2M_TT[0][IDV[i]] = B2V[0][i];
				B2M_TT[1][IDV[i]] = B2V[1][i];
				B2M_TT[2][IDV[i]] = B2V[2][i];
			}
			double B2M_TQ[3][9];
			//B2M_TQ = B2M_TT * D
			for (int i = 0; i < 9; i++) {//这里还有一定优化空间，D矩阵有部分值固定为0，可以省略计算，但是需要展开循环，会导致代码篇幅比较大
				B2M_TQ[0][i] = B2M_TT[0][0] * D[0][i] + B2M_TT[0][1] * D[1][i] + B2M_TT[0][2] * D[2][i]
					+ B2M_TT[0][3] * D[3][i] + B2M_TT[0][4] * D[4][i] + B2M_TT[0][5] * D[5][i]
					+ B2M_TT[0][6] * D[6][i] + B2M_TT[0][7] * D[7][i] + B2M_TT[0][8] * D[8][i];
				B2M_TQ[1][i] = B2M_TT[1][0] * D[0][i] + B2M_TT[1][1] * D[1][i] + B2M_TT[1][2] * D[2][i]
					+ B2M_TT[1][3] * D[3][i] + B2M_TT[1][4] * D[4][i] + B2M_TT[1][5] * D[5][i]
					+ B2M_TT[1][6] * D[6][i] + B2M_TT[1][7] * D[7][i] + B2M_TT[1][8] * D[8][i];
				B2M_TQ[2][i] = B2M_TT[2][0] * D[0][i] + B2M_TT[2][1] * D[1][i] + B2M_TT[2][2] * D[2][i]
					+ B2M_TT[2][3] * D[3][i] + B2M_TT[2][4] * D[4][i] + B2M_TT[2][5] * D[5][i]
					+ B2M_TT[2][6] * D[6][i] + B2M_TT[2][7] * D[7][i] + B2M_TT[2][8] * D[8][i];
			}

			//B2M_QQk = TE_STRESS * B2M_TQ			因为B2M_QQk会组装到B2M_QQ_5中，这里不创建B2M_QQk直接组装了
			for (int i = 0; i < 9; i++) {//不计算铺层所需BIG_BB，这里就直接累加了，且是在加完后循环外再除以4
				B2M_QQ_5[0][IDI[K][i]] += TE_STRESS[0][0] * B2M_TQ[0][i] + TE_STRESS[0][1] * B2M_TQ[1][i] + TE_STRESS[0][2] * B2M_TQ[2][i];
				B2M_QQ_5[1][IDI[K][i]] += TE_STRESS[1][0] * B2M_TQ[0][i] + TE_STRESS[1][1] * B2M_TQ[1][i] + TE_STRESS[1][2] * B2M_TQ[2][i];
				B2M_QQ_5[2][IDI[K][i]] += TE_STRESS[2][0] * B2M_TQ[0][i] + TE_STRESS[2][1] * B2M_TQ[1][i] + TE_STRESS[2][2] * B2M_TQ[2][i];
			}


			//============计算当前子三角形的B3M_QQk，并根据组装至B3M_QQ_5
			double B3M_TT[3][9];
			for (int i = 0; i < 9; i++) {
				B3M_TT[0][IDV[i]] = B3V[0][i];
				B3M_TT[1][IDV[i]] = B3V[1][i];
				B3M_TT[2][IDV[i]] = B3V[2][i];
			}
			double B3M_TQ[3][9];
			//B3M_TQ = B3M_TT * D
			for (int i = 0; i < 9; i++) {//这里还有一定优化空间，D矩阵有部分值固定为0，可以省略计算，但是需要展开循环，会导致代码篇幅比较大
				B3M_TQ[0][i] = B3M_TT[0][0] * D[0][i] + B3M_TT[0][1] * D[1][i] + B3M_TT[0][2] * D[2][i]
					+ B3M_TT[0][3] * D[3][i] + B3M_TT[0][4] * D[4][i] + B3M_TT[0][5] * D[5][i]
					+ B3M_TT[0][6] * D[6][i] + B3M_TT[0][7] * D[7][i] + B3M_TT[0][8] * D[8][i];
				B3M_TQ[1][i] = B3M_TT[1][0] * D[0][i] + B3M_TT[1][1] * D[1][i] + B3M_TT[1][2] * D[2][i]
					+ B3M_TT[1][3] * D[3][i] + B3M_TT[1][4] * D[4][i] + B3M_TT[1][5] * D[5][i]
					+ B3M_TT[1][6] * D[6][i] + B3M_TT[1][7] * D[7][i] + B3M_TT[1][8] * D[8][i];
				B3M_TQ[2][i] = 0.0;/*B3M_TT[2][0] * D[0][i] + B3M_TT[2][1] * D[1][i] + B3M_TT[2][2] * D[2][i]
							 + B3M_TT[2][3] * D[3][i] + B3M_TT[2][4] * D[4][i] + B3M_TT[2][5] * D[5][i]
							 + B3M_TT[2][6] * D[6][i] + B3M_TT[2][7] * D[7][i] + B3M_TT[2][8] * D[8][i]*/
			}
			//B3M_QQk = TE_SHEAR * B3M_TQ			因为B3M_QQk会组装到B3M_QQ_5中，这里不创建B3M_QQk直接组装了
			for (int i = 0; i < 9; i++) {
				B3M_QQ_5[0][IDI[K][i]] += TE_SHEAR[0][0] * B3M_TQ[0][i] + TE_SHEAR[0][1] * B3M_TQ[1][i] + TE_SHEAR[0][2] * B3M_TQ[2][i];
				B3M_QQ_5[1][IDI[K][i]] += TE_SHEAR[1][0] * B3M_TQ[0][i] + TE_SHEAR[1][1] * B3M_TQ[1][i] + TE_SHEAR[1][2] * B3M_TQ[2][i];
				/*B3M_QQ_5[2][IDI[K][i]] += TE_SHEAR[2][0] * B3M_TQ[0][i] + TE_SHEAR[2][1] * B3M_TQ[1][i] + TE_SHEAR[2][2] * B3M_TQ[2][i];*/
			}
		}

		//四边形的PHI_SQ为四个子三角形的平均、四边形的面积AREA_QUAD为四个子三角形的和
		for (int i = 0; i < 4; i++) { PHI_SQ += PHI_SQ_TRIA[i]; }	PHI_SQ *= 0.25;
		for (int i = 0; i < 4; i++) { AREA_QUAD += AREA_TRIA[i]; }

		//来自四个三角形的组装，需要除以4，组装完成后执行此操作效率更高
		for (int i = 0; i < 15; i++) {
			B2M_QQ_5[0][i] *= 0.25;			B3M_QQ_5[0][i] *= 0.25;
			B2M_QQ_5[1][i] *= 0.25;			B3M_QQ_5[1][i] *= 0.25;
			B2M_QQ_5[2][i] *= 0.25;			/*B3M_QQ_5[2][i] *= 0.25;*/
		}

		STATIC_CONDENSATION(KM_QQ_5, B2M_QQ_5, B3M_QQ_5, k0);
	}

	void CQUAD4::caculateKB(double* A, 
							double* B, 
							double AREA, 
							double A4, 
							double A42, 
							double D[8][8], 
							double KB[9][9], 
							double KS[9][9], 
							double& PHI_SQ)
	{
		double FXX[3][3] = { 0. }, FYY[3][3] = { 0. }, FXY[3][3] = { 0. };

		FXX[0][0] = B[0] * B[0] / A42;
		FXX[0][1] = B[0] * B[1] / A42;
		FXX[1][0] = B[1] * B[0] / A42;
		FXX[1][1] = B[1] * B[1] / A42;

		FYY[0][0] = A[0] * A[0] / A42;
		FYY[0][1] = A[0] * A[1] / A42;
		FYY[0][2] = A[0] * A[2] / A42;
		FYY[1][0] = A[1] * A[0] / A42;
		FYY[1][1] = A[1] * A[1] / A42;
		FYY[1][2] = A[1] * A[2] / A42;
		FYY[2][0] = A[2] * A[0] / A42;
		FYY[2][1] = A[2] * A[1] / A42;
		FYY[2][2] = A[2] * A[2] / A42;

		FXY[0][0] = B[0] * A[0] / A42;
		FXY[0][1] = B[0] * A[1] / A42;
		FXY[0][2] = B[0] * A[2] / A42;
		FXY[1][0] = B[1] * A[0] / A42;
		FXY[1][1] = B[1] * A[1] / A42;
		FXY[1][2] = B[1] * A[2] / A42;

		KB[3][3] = AREA * (D[4][4] * FYY[0][0] + D[4][5] * (FXY[0][0] + FXY[0][0]) + D[5][5] * FXX[0][0]);
		KB[3][4] = AREA * (D[4][4] * FYY[0][1] + D[4][5] * (FXY[0][1] + FXY[1][0]) + D[5][5] * FXX[0][1]);
		KB[3][5] = AREA * (D[4][4] * FYY[0][2] + D[4][5] * (FXY[0][2]/*+FXY[2][0]*/)/*+D[5][5] * FXX[0][2]*/);				//注释部分因为FXX，FXY固定为0
		KB[4][3] = KB[3][4];
		KB[4][4] = AREA * (D[4][4] * FYY[1][1] + D[4][5] * (FXY[1][1] + FXY[1][1]) + D[5][5] * FXX[1][1]);
		KB[4][5] = AREA * (D[4][4] * FYY[1][2] + D[4][5] * (FXY[1][2]/*+FXY[2][1]*/)/*+D[5][5] * FXX[1][2]*/);				//注释部分因为FXX，FXY固定为0
		KB[5][3] = KB[3][5];
		KB[5][4] = KB[4][5];
		KB[5][5] = AREA * (D[4][4] * FYY[2][2] + D[4][5] * (FXY[2][2]/*+FXY[2][2]*/)/*+D[5][5] * FXX[2][2]*/);				//注释部分因为FXX，FXY固定为0


		KB[3][6] = -AREA * (D[3][4] * FXY[0][0] + D[4][5] * FYY[0][0] + D[3][5] * FXX[0][0] + D[5][5] * FXY[0][0]);
		KB[3][7] = -AREA * (D[3][4] * FXY[1][0] + D[4][5] * FYY[0][1] + D[3][5] * FXX[0][1] + D[5][5] * FXY[0][1]);
		KB[3][8] = -AREA * (/*D[3][4]*FXY[2][0]*/+D[4][5] * FYY[0][2] +/*D[3][5]*FXX[0][2]+*/ D[5][5] * FXY[0][2]);			//注释部分因为FXX，FXY固定为0
		KB[4][6] = -AREA * (D[3][4] * FXY[0][1] + D[4][5] * FYY[1][0] + D[3][5] * FXX[1][0] + D[5][5] * FXY[1][0]);
		KB[4][7] = -AREA * (D[3][4] * FXY[1][1] + D[4][5] * FYY[1][1] + D[3][5] * FXX[1][1] + D[5][5] * FXY[1][1]);
		KB[4][8] = -AREA * (/*D[3][4]*FXY[2][1]*/+D[4][5] * FYY[1][2] /*+D[3][5]*FXX[1][2]*/ + D[5][5] * FXY[1][2]);		//注释部分因为FXX，FXY固定为0
		KB[5][6] = -AREA * (D[3][4] * FXY[0][2] + D[4][5] * FYY[2][0] /*+D[3][5]*FXX[2][0]  + D[5][5] * FXY[2][0]*/);		//注释部分因为FXX，FXY固定为0
		KB[5][7] = -AREA * (D[3][4] * FXY[1][2] + D[4][5] * FYY[2][1] /*+D[3][5]*FXX[2][1]  + D[5][5] * FXY[2][1]*/);		//注释部分因为FXX，FXY固定为0
		KB[5][8] = -AREA * (/*D[3][4]*FXY[2][2]*/+D[4][5] * FYY[2][2] /*+D[3][5]*FXX[2][2]  + D[5][5] * FXY[2][2]*/);		//注释部分因为FXX，FXY固定为0
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				KB[j + 6][i + 3] = KB[i + 3][j + 6];
			}
		}

		KB[6][6] = AREA * (D[3][3] * FXX[0][0] + D[3][5] * (FXY[0][0] + FXY[0][0]) + D[5][5] * FYY[0][0]);
		KB[6][7] = AREA * (D[3][3] * FXX[0][1] + D[3][5] * (FXY[0][1] + FXY[1][0]) + D[5][5] * FYY[0][1]);
		KB[6][8] = AREA * (/*D[3][3]*FXX[0][2]+*/D[3][5] * (FXY[0][2] + FXY[2][0]) + D[5][5] * FYY[0][2]);//注释部分因为FXX为0
		KB[7][6] = KB[6][7];
		KB[7][7] = AREA * (D[3][3] * FXX[1][1] + D[3][5] * (FXY[1][1] + FXY[1][1]) + D[5][5] * FYY[1][1]);
		KB[7][8] = AREA * (/*D[3][3]* XX[1][2]+*/D[3][5] * (FXY[1][2] + FXY[2][1]) + D[5][5] * FYY[1][2]);//注释部分因为FXX为0
		KB[8][6] = KB[6][8];
		KB[8][7] = KB[7][8];
		KB[8][8] = AREA * (/*D[3][3]*FXX[2][2]+*/D[3][5] * (FXY[2][2] + FXY[2][2]) + D[5][5] * FYY[2][2]);//注释部分因为FXX为0

		double A1[3][3] = { 0. }, B2[3][3] = { 0. }, S1[3][3] = { 0. }, S2[3][3] = { 0. };
		A1[0][1] = A[1] * A[2] / A4;
		A1[0][2] = -A1[0][1];
		A1[1][0] = -A[0] * A[2] / A4;
		A1[1][2] = -A1[1][0];
		A1[2][0] = A[0] * A[1] / A4;
		A1[2][1] = -A1[2][0];

		B2[2][0] = -B[0] * B[1] / A4;
		B2[2][1] = -B2[2][0];

		S1[0][0] = (A[1] * B[2] - A[2] * B[1]) / A4 + 1.0;
		S1[0][2] = A[2] * B[1] / A4;
		S1[1][1] = (A[2] * B[0] - A[0] * B[2]) / A4 + 1.0;
		S1[1][2] = -A[2] * B[0] / A4;
		S1[2][0] = -A[0] * B[1] / A4;
		S1[2][1] = A[1] * B[0] / A4;
		S1[2][2] = (A[0] * B[1] - A[1] * B[0]) / A4 + 1.0;

		S2[0][0] = (A[1] * B[2] - A[2] * B[1]) / A4 + 1.0;
		S2[0][1] = B[1] * A[2] / A4;
		S2[1][0] = -B[0] * A[2] / A4;
		S2[1][1] = (A[2] * B[0] - A[0] * B[2]) / A4 + 1.0;
		S2[2][0] = B[0] * A[1] / A4;
		S2[2][1] = -B[1] * A[0] / A4;
		S2[2][2] = (A[0] * B[1] - A[1] * B[0]) / A4 + 1.0;

		double I00[3][3] = { 0. }, IX0[3][3] = { 0. }, IY0[3][3] = { 0. };
		I00[0][0] = AREA / 6.0;
		I00[0][1] = AREA / 12.0;
		I00[0][2] = I00[0][1];
		I00[1][0] = I00[0][1];
		I00[1][1] = I00[0][0];
		I00[1][2] = I00[0][1];
		I00[2][0] = I00[0][1];
		I00[2][1] = I00[0][1];
		I00[2][2] = I00[0][0];

		IX0[0][0] = B[0] / 6.0;
		IX0[0][1] = IX0[0][0];
		IX0[0][2] = IX0[0][0];
		IX0[1][0] = B[1] / 6.0;
		IX0[1][1] = IX0[1][0];
		IX0[1][2] = IX0[1][0];

		IY0[0][0] = A[0] / 6.0;
		IY0[0][1] = IY0[0][0];
		IY0[0][2] = IY0[0][0];
		IY0[1][0] = A[1] / 6.0;
		IY0[1][1] = IY0[1][0];
		IY0[1][2] = IY0[1][0];
		IY0[2][0] = A[2] / 6.0;
		IY0[2][1] = IY0[2][0];
		IY0[2][2] = IY0[2][0];

		double T1[3][3] = { 0. }, T2[3][3] = { 0. }, T3[3][3] = { 0. }, T4[3][3] = { 0. };
		T1[0][0] =/*D[6][6]*B2[0][0]+*/ D[6][7] * S1[0][0];
		T1[0][1] =/*D[6][6]*B2[0][1]  + D[6][7] * S1[0][1]*/0.0;
		T1[0][2] =/*D[6][6]*B2[0][2]+*/ D[6][7] * S1[0][2];
		T1[1][0] =/*D[6][6]*B2[1][0]  + D[6][7] * S1[1][0]*/0.0;
		T1[1][1] =/*D[6][6]*B2[1][1]+*/ D[6][7] * S1[1][1];
		T1[1][2] =/*D[6][6]*B2[1][2]+*/ D[6][7] * S1[1][2];
		T1[2][0] = D[6][6] * B2[2][0] + D[6][7] * S1[2][0];
		T1[2][1] = D[6][6] * B2[2][1] + D[6][7] * S1[2][1];
		T1[2][2] =/*D[6][6]*B2[2][2]+*/ D[6][7] * S1[2][2];

		T2[0][0] = /*D[6][7]*B2[0][0]+*/D[7][7] * S1[0][0];
		T2[0][1] = /*D[6][7]*B2[0][1] + D[7][7] * S1[0][1]*/0.0;
		T2[0][2] = /*D[6][7]*B2[0][2]+*/D[7][7] * S1[0][2];
		T2[1][0] = /*D[6][7]*B2[1][0] + D[7][7] * S1[1][0]*/0.0;
		T2[1][1] = /*D[6][7]*B2[1][1]+*/D[7][7] * S1[1][1];
		T2[1][2] = /*D[6][7]*B2[1][2]+*/D[7][7] * S1[1][2];
		T2[2][0] = D[6][7] * B2[2][0] + D[7][7] * S1[2][0];
		T2[2][1] = D[6][7] * B2[2][1] + D[7][7] * S1[2][1];
		T2[2][2] = D[6][7] * B2[2][2] + D[7][7] * S1[2][2];

		T3[0][0] = D[6][6] * S2[0][0]/*+D[6][7] * A1[0][0]*/;
		T3[0][1] = D[6][6] * S2[0][1] + D[6][7] * A1[0][1];
		T3[0][2] =/*D[6][6]*S2[0][2]+*/ D[6][7] * A1[0][2];
		T3[1][0] = D[6][6] * S2[1][0] + D[6][7] * A1[1][0];
		T3[1][1] = D[6][6] * S2[1][1]/*+D[6][7] * A1[1][1]*/;
		T3[1][2] =/*D[6][6]*S2[1][2]+*/ D[6][7] * A1[1][2];
		T3[2][0] = D[6][6] * S2[2][0] + D[6][7] * A1[2][0];
		T3[2][1] = D[6][6] * S2[2][1] + D[6][7] * A1[2][1];
		T3[2][2] = D[6][6] * S2[2][2]/*+D[6][7] * A1[2][2]*/;

		T4[0][0] = D[6][7] * S2[0][0] /*+D[7][7] * A1[0][0]*/;
		T4[0][1] = D[6][7] * S2[0][1] + D[7][7] * A1[0][1];
		T4[0][2] =/*D[6][7] * S2[0][2]+*/D[7][7] * A1[0][2];
		T4[1][0] = D[6][7] * S2[1][0] + D[7][7] * A1[1][0];
		T4[1][1] = D[6][7] * S2[1][1] /*+D[7][7] * A1[1][1]*/;
		T4[1][2] =/*D[6][7] * S2[1][2]+*/D[7][7] * A1[1][2];
		T4[2][0] = D[6][7] * S2[2][0] + D[7][7] * A1[2][0];
		T4[2][1] = D[6][7] * S2[2][1] + D[7][7] * A1[2][1];
		T4[2][2] = D[6][7] * S2[2][2] /*+D[7][7] * A1[2][2]*/;

		KS[0][0] = AREA * (D[6][6] * FXX[0][0] + D[6][7] * (FXY[0][0] + FXY[0][0]) + D[7][7] * FYY[0][0]);
		KS[0][1] = AREA * (D[6][6] * FXX[0][1] + D[6][7] * (FXY[0][1] + FXY[1][0]) + D[7][7] * FYY[0][1]);
		KS[0][2] = AREA * (/*D[6][6] * FXX[0][2]+*/ D[6][7] * (FXY[0][2] /*+ FXY[2][0]*/) + D[7][7] * FYY[0][2]);
		KS[1][0] = AREA * (D[6][6] * FXX[1][0] + D[6][7] * (FXY[1][0] + FXY[0][1]) + D[7][7] * FYY[1][0]);
		KS[1][1] = AREA * (D[6][6] * FXX[1][1] + D[6][7] * (FXY[1][1] + FXY[1][1]) + D[7][7] * FYY[1][1]);
		KS[1][2] = AREA * (/*D[6][6] * FXX[1][2]+*/ D[6][7] * (FXY[1][2] /*+ FXY[2][1]*/) + D[7][7] * FYY[1][2]);
		KS[2][0] = AREA * (/*D[6][6] * FXX[2][0]+*/ D[6][7] * (/*FXY[2][0]+*/FXY[0][2]) + D[7][7] * FYY[2][0]);
		KS[2][1] = AREA * (/*D[6][6] * FXX[2][1]+*/ D[6][7] * (/*FXY[2][1]+*/FXY[1][2]) + D[7][7] * FYY[2][1]);
		KS[2][2] = AREA * (/*D[6][6] * FXX[2][2]+*/ D[6][7] * (/*FXY[2][2]+*/FXY[2][2]) + D[7][7] * FYY[2][2]);

		double DUM1[3][3] = { 0. }, DUM2[3][3] = { 0. };
		//! DUM1 = IX0 * T1;
		DUM1[0][0] = IX0[0][0] * T1[0][0] + /*IX0[0][1] * T1[1][0] +*/ IX0[0][2] * T1[2][0];
		DUM1[0][1] = /*IX0[0][0] * T1[0][1] +*/ IX0[0][1] * T1[1][1] + IX0[0][2] * T1[2][1];
		DUM1[0][2] = IX0[0][0] * T1[0][2] + IX0[0][1] * T1[1][2] + IX0[0][2] * T1[2][2];
		DUM1[1][0] = IX0[1][0] * T1[0][0] + /*IX0[1][1] * T1[1][0] +*/ IX0[1][2] * T1[2][0];
		DUM1[1][1] = /*IX0[1][0] * T1[0][1] +*/ IX0[1][1] * T1[1][1] + IX0[1][2] * T1[2][1];
		DUM1[1][2] = IX0[1][0] * T1[0][2] + IX0[1][1] * T1[1][2] + IX0[1][2] * T1[2][2];
		DUM1[2][0] = /*IX0[2][0] * T1[0][0] + IX0[2][1] * T1[1][0] + IX0[2][2] * T1[2][0]*/ 0.0;
		DUM1[2][1] = /*IX0[2][0] * T1[0][1] + IX0[2][1] * T1[1][1] + IX0[2][2] * T1[2][1]*/ 0.0;
		DUM1[2][2] = /*IX0[2][0] * T1[0][2] + IX0[2][1] * T1[1][2] + IX0[2][2] * T1[2][2]*/ 0.0;
		//! DUM2 = IY0 * T2;
		DUM2[0][0] = IY0[0][0] * T2[0][0] +/*IY0[0][1]*T2[1][0]+*/  IY0[0][2] * T2[2][0];
		DUM2[0][1] = /*IY0[0][0]*T2[0][1]+*/IY0[0][1] * T2[1][1] + IY0[0][2] * T2[2][1];
		DUM2[0][2] = IY0[0][0] * T2[0][2] + IY0[0][1] * T2[1][2] + IY0[0][2] * T2[2][2];
		DUM2[1][0] = IY0[1][0] * T2[0][0] +/*IY0[1][1]*T2[1][0]+*/  IY0[1][2] * T2[2][0];
		DUM2[1][1] = /*IY0[1][0]*T2[0][1]+*/IY0[1][1] * T2[1][1] + IY0[1][2] * T2[2][1];
		DUM2[1][2] = IY0[1][0] * T2[0][2] + IY0[1][1] * T2[1][2] + IY0[1][2] * T2[2][2];
		DUM2[2][0] = IY0[2][0] * T2[0][0] +/*IY0[2][1]*T2[1][0]+*/  IY0[2][2] * T2[2][0];
		DUM2[2][1] = /*IY0[2][0]*T2[0][1]+*/IY0[2][1] * T2[1][1] + IY0[2][2] * T2[2][1];
		DUM2[2][2] = IY0[2][0] * T2[0][2] + IY0[2][1] * T2[1][2] + IY0[2][2] * T2[2][2];
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				KS[i][j + 3] = -DUM1[i][j] - DUM2[i][j];
				KS[j + 3][i] = KS[i][j + 3];
			}
		}

		//! DUM1 = IX0 * T3;
		DUM1[0][0] = IX0[0][0] * T3[0][0] + IX0[0][1] * T3[1][0] + IX0[0][2] * T3[2][0];
		DUM1[0][1] = IX0[0][0] * T3[0][1] + IX0[0][1] * T3[1][1] + IX0[0][2] * T3[2][1];
		DUM1[0][2] = IX0[0][0] * T3[0][2] + IX0[0][1] * T3[1][2] + IX0[0][2] * T3[2][2];
		DUM1[1][0] = IX0[1][0] * T3[0][0] + IX0[1][1] * T3[1][0] + IX0[1][2] * T3[2][0];
		DUM1[1][1] = IX0[1][0] * T3[0][1] + IX0[1][1] * T3[1][1] + IX0[1][2] * T3[2][1];
		DUM1[1][2] = IX0[1][0] * T3[0][2] + IX0[1][1] * T3[1][2] + IX0[1][2] * T3[2][2];
		DUM1[2][0] = /*IX0[2][0] * T3[0][0] + IX0[2][1] * T3[1][0] + IX0[2][2] * T3[2][0]*/0.0;
		DUM1[2][1] = /*IX0[2][0] * T3[0][1] + IX0[2][1] * T3[1][1] + IX0[2][2] * T3[2][1]*/0.0;
		DUM1[2][2] = /*IX0[2][0] * T3[0][2] + IX0[2][1] * T3[1][2] + IX0[2][2] * T3[2][2]*/0.0;

		//! DUM2 = IY0 * T4;
		DUM2[0][0] = IY0[0][0] * T4[0][0] + IY0[0][1] * T4[1][0] + IY0[0][2] * T4[2][0];
		DUM2[0][1] = IY0[0][0] * T4[0][1] + IY0[0][1] * T4[1][1] + IY0[0][2] * T4[2][1];
		DUM2[0][2] = IY0[0][0] * T4[0][2] + IY0[0][1] * T4[1][2] + IY0[0][2] * T4[2][2];
		DUM2[1][0] = IY0[1][0] * T4[0][0] + IY0[1][1] * T4[1][0] + IY0[1][2] * T4[2][0];
		DUM2[1][1] = IY0[1][0] * T4[0][1] + IY0[1][1] * T4[1][1] + IY0[1][2] * T4[2][1];
		DUM2[1][2] = IY0[1][0] * T4[0][2] + IY0[1][1] * T4[1][2] + IY0[1][2] * T4[2][2];
		DUM2[2][0] = IY0[2][0] * T4[0][0] + IY0[2][1] * T4[1][0] + IY0[2][2] * T4[2][0];
		DUM2[2][1] = IY0[2][0] * T4[0][1] + IY0[2][1] * T4[1][1] + IY0[2][2] * T4[2][1];
		DUM2[2][2] = IY0[2][0] * T4[0][2] + IY0[2][1] * T4[1][2] + IY0[2][2] * T4[2][2];
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				KS[i][j + 6] = DUM1[i][j] + DUM2[i][j];
				KS[j + 6][i] = KS[i][j + 6];
			}
		}

		double I00T1[3][3] = { 0. }, I00T2[3][3] = { 0. }, I00T3[3][3] = { 0. }, I00T4[3][3] = { 0. };
		//I00T1 = I00 * T1;
		I00T1[0][0] = I00[0][0] * T1[0][0] +   /*I00[0][1] * T1[1][0] +*/ I00[0][2] * T1[2][0];
		I00T1[0][1] = /*I00[0][0] * T1[0][1]+*/  I00[0][1] * T1[1][1] + I00[0][2] * T1[2][1];
		I00T1[0][2] = I00[0][0] * T1[0][2] + I00[0][1] * T1[1][2] + I00[0][2] * T1[2][2];
		I00T1[1][0] = I00[1][0] * T1[0][0] +   /*I00[1][1] * T1[1][0] +*/ I00[1][2] * T1[2][0];
		I00T1[1][1] = /*I00[1][0] * T1[0][1] +*/ I00[1][1] * T1[1][1] + I00[1][2] * T1[2][1];
		I00T1[1][2] = I00[1][0] * T1[0][2] + I00[1][1] * T1[1][2] + I00[1][2] * T1[2][2];
		I00T1[2][0] = I00[2][0] * T1[0][0] +   /*I00[2][1] * T1[1][0] +*/ I00[2][2] * T1[2][0];
		I00T1[2][1] = /*I00[2][0] * T1[0][1] +*/ I00[2][1] * T1[1][1] + I00[2][2] * T1[2][1];
		I00T1[2][2] = I00[2][0] * T1[0][2] + I00[2][1] * T1[1][2] + I00[2][2] * T1[2][2];
		//I00T2 = I00 * T2;
		I00T2[0][0] = I00[0][0] * T2[0][0] +  /*I00[0][1] * T2[1][0] +*/ I00[0][2] * T2[2][0];
		I00T2[0][1] = /*I00[0][0] * T2[0][1]+*/ I00[0][1] * T2[1][1] + I00[0][2] * T2[2][1];
		I00T2[0][2] = I00[0][0] * T2[0][2] + I00[0][1] * T2[1][2] + I00[0][2] * T2[2][2];
		I00T2[1][0] = I00[1][0] * T2[0][0] +  /*I00[1][1] * T2[1][0] +*/ I00[1][2] * T2[2][0];
		I00T2[1][1] = /*I00[1][0] * T2[0][1]+*/ I00[1][1] * T2[1][1] + I00[1][2] * T2[2][1];
		I00T2[1][2] = I00[1][0] * T2[0][2] + I00[1][1] * T2[1][2] + I00[1][2] * T2[2][2];
		I00T2[2][0] = I00[2][0] * T2[0][0] +  /*I00[2][1] * T2[1][0] +*/ I00[2][2] * T2[2][0];
		I00T2[2][1] = /*I00[2][0] * T2[0][1]+*/ I00[2][1] * T2[1][1] + I00[2][2] * T2[2][1];
		I00T2[2][2] = I00[2][0] * T2[0][2] + I00[2][1] * T2[1][2] + I00[2][2] * T2[2][2];
		//I00T3 = I00 * T3;
		I00T3[0][0] = I00[0][0] * T3[0][0] + I00[0][1] * T3[1][0] + I00[0][2] * T3[2][0];
		I00T3[0][1] = I00[0][0] * T3[0][1] + I00[0][1] * T3[1][1] + I00[0][2] * T3[2][1];
		I00T3[0][2] = I00[0][0] * T3[0][2] + I00[0][1] * T3[1][2] + I00[0][2] * T3[2][2];
		I00T3[1][0] = I00[1][0] * T3[0][0] + I00[1][1] * T3[1][0] + I00[1][2] * T3[2][0];
		I00T3[1][1] = I00[1][0] * T3[0][1] + I00[1][1] * T3[1][1] + I00[1][2] * T3[2][1];
		I00T3[1][2] = I00[1][0] * T3[0][2] + I00[1][1] * T3[1][2] + I00[1][2] * T3[2][2];
		I00T3[2][0] = I00[2][0] * T3[0][0] + I00[2][1] * T3[1][0] + I00[2][2] * T3[2][0];
		I00T3[2][1] = I00[2][0] * T3[0][1] + I00[2][1] * T3[1][1] + I00[2][2] * T3[2][1];
		I00T3[2][2] = I00[2][0] * T3[0][2] + I00[2][1] * T3[1][2] + I00[2][2] * T3[2][2];
		//I00T4 = I00 * T4;
		I00T4[0][0] = I00[0][0] * T4[0][0] + I00[0][1] * T4[1][0] + I00[0][2] * T4[2][0];
		I00T4[0][1] = I00[0][0] * T4[0][1] + I00[0][1] * T4[1][1] + I00[0][2] * T4[2][1];
		I00T4[0][2] = I00[0][0] * T4[0][2] + I00[0][1] * T4[1][2] + I00[0][2] * T4[2][2];
		I00T4[1][0] = I00[1][0] * T4[0][0] + I00[1][1] * T4[1][0] + I00[1][2] * T4[2][0];
		I00T4[1][1] = I00[1][0] * T4[0][1] + I00[1][1] * T4[1][1] + I00[1][2] * T4[2][1];
		I00T4[1][2] = I00[1][0] * T4[0][2] + I00[1][1] * T4[1][2] + I00[1][2] * T4[2][2];
		I00T4[2][0] = I00[2][0] * T4[0][0] + I00[2][1] * T4[1][0] + I00[2][2] * T4[2][0];
		I00T4[2][1] = I00[2][0] * T4[0][1] + I00[2][1] * T4[1][1] + I00[2][2] * T4[2][1];
		I00T4[2][2] = I00[2][0] * T4[0][2] + I00[2][1] * T4[1][2] + I00[2][2] * T4[2][2];

		//DUM1 = B2.transpose() * I00 * T1;
		DUM1[0][0] = /*B2[0][0] * I00T1[0][0] + B2[1][0] * I00T1[1][0] +*/ B2[2][0] * I00T1[2][0];
		DUM1[0][1] = /*B2[0][0] * I00T1[0][1] + B2[1][0] * I00T1[1][1] +*/ B2[2][0] * I00T1[2][1];
		DUM1[0][2] = /*B2[0][0] * I00T1[0][2] + B2[1][0] * I00T1[1][2] +*/ B2[2][0] * I00T1[2][2];
		DUM1[1][0] = /*B2[0][1] * I00T1[0][0] + B2[1][1] * I00T1[1][0] +*/ B2[2][1] * I00T1[2][0];
		DUM1[1][1] = /*B2[0][1] * I00T1[0][1] + B2[1][1] * I00T1[1][1] +*/ B2[2][1] * I00T1[2][1];
		DUM1[1][2] = /*B2[0][1] * I00T1[0][2] + B2[1][1] * I00T1[1][2] +*/ B2[2][1] * I00T1[2][2];
		DUM1[2][0] = /*B2[0][2] * I00T1[0][0] + B2[1][2] * I00T1[1][0] +   B2[2][2] * I00T1[2][0]*/0.0;
		DUM1[2][1] = /*B2[0][2] * I00T1[0][1] + B2[1][2] * I00T1[1][1] +   B2[2][2] * I00T1[2][1]*/0.0;
		DUM1[2][2] = /*B2[0][2] * I00T1[0][2] + B2[1][2] * I00T1[1][2] +   B2[2][2] * I00T1[2][2]*/0.0;
		//DUM2 = S1.transpose() * I00 * T2;
		DUM2[0][0] = S1[0][0] * I00T2[0][0] + /*S1[1][0] * I00T2[1][0]+*/ S1[2][0] * I00T2[2][0];
		DUM2[0][1] = S1[0][0] * I00T2[0][1] + /*S1[1][0] * I00T2[1][1]+*/ S1[2][0] * I00T2[2][1];
		DUM2[0][2] = S1[0][0] * I00T2[0][2] + /*S1[1][0] * I00T2[1][2]+*/ S1[2][0] * I00T2[2][2];
		DUM2[1][0] = /*S1[0][1] * I00T2[0][0]+*/  S1[1][1] * I00T2[1][0] + S1[2][1] * I00T2[2][0];
		DUM2[1][1] = /*S1[0][1] * I00T2[0][1]+*/  S1[1][1] * I00T2[1][1] + S1[2][1] * I00T2[2][1];
		DUM2[1][2] = /*S1[0][1] * I00T2[0][2]+*/  S1[1][1] * I00T2[1][2] + S1[2][1] * I00T2[2][2];
		DUM2[2][0] = S1[0][2] * I00T2[0][0] + S1[1][2] * I00T2[1][0] + S1[2][2] * I00T2[2][0];
		DUM2[2][1] = S1[0][2] * I00T2[0][1] + S1[1][2] * I00T2[1][1] + S1[2][2] * I00T2[2][1];
		DUM2[2][2] = S1[0][2] * I00T2[0][2] + S1[1][2] * I00T2[1][2] + S1[2][2] * I00T2[2][2];
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				KS[i + 3][j + 3] = DUM1[i][j] + DUM2[i][j];
			}
		}

		//DUM1 = B2.transpose()* I00 * T3;
		DUM1[0][0] = /*B2[0][0] * I00T3[0][0] + B2T[1][0] * I00T3[1][0] +*/ B2[2][0] * I00T3[2][0];
		DUM1[0][1] = /*B2[0][0] * I00T3[0][1] + B2T[1][0] * I00T3[1][1] +*/ B2[2][0] * I00T3[2][1];
		DUM1[0][2] = /*B2[0][0] * I00T3[0][2] + B2T[1][0] * I00T3[1][2] +*/ B2[2][0] * I00T3[2][2];
		DUM1[1][0] = /*B2[0][1] * I00T3[0][0] + B2T[1][1] * I00T3[1][0] +*/ B2[2][1] * I00T3[2][0];
		DUM1[1][1] = /*B2[0][1] * I00T3[0][1] + B2T[1][1] * I00T3[1][1] +*/ B2[2][1] * I00T3[2][1];
		DUM1[1][2] = /*B2[0][1] * I00T3[0][2] + B2T[1][1] * I00T3[1][2] +*/ B2[2][1] * I00T3[2][2];
		DUM1[2][0] = /*B2[0][2] * I00T3[0][0] + B2T[1][2] * I00T3[1][0] +   B2[2][2] * I00T3[2][0]*/0.0;
		DUM1[2][1] = /*B2[0][2] * I00T3[0][1] + B2T[1][2] * I00T3[1][1] +   B2[2][2] * I00T3[2][1]*/0.0;
		DUM1[2][2] = /*B2[0][2] * I00T3[0][2] + B2T[1][2] * I00T3[1][2] +   B2[2][2] * I00T3[2][2]*/0.0;
		//DUM2 = S1.transpose() * I00 * T4;
		DUM2[0][0] = S1[0][0] * I00T4[0][0] + /*S1[1][0] * I00T4[1][0] +*/ S1[2][0] * I00T4[2][0];
		DUM2[0][1] = S1[0][0] * I00T4[0][1] + /*S1[1][0] * I00T4[1][1] +*/ S1[2][0] * I00T4[2][1];
		DUM2[0][2] = S1[0][0] * I00T4[0][2] + /*S1[1][0] * I00T4[1][2] +*/ S1[2][0] * I00T4[2][2];
		DUM2[1][0] = /*S1[0][1] * I00T4[0][0] +*/ S1[1][1] * I00T4[1][0] + S1[2][1] * I00T4[2][0];
		DUM2[1][1] = /*S1[0][1] * I00T4[0][1] +*/ S1[1][1] * I00T4[1][1] + S1[2][1] * I00T4[2][1];
		DUM2[1][2] = /*S1[0][1] * I00T4[0][2] +*/ S1[1][1] * I00T4[1][2] + S1[2][1] * I00T4[2][2];
		DUM2[2][0] = S1[0][2] * I00T4[0][0] + S1[1][2] * I00T4[1][0] + S1[2][2] * I00T4[2][0];
		DUM2[2][1] = S1[0][2] * I00T4[0][1] + S1[1][2] * I00T4[1][1] + S1[2][2] * I00T4[2][1];
		DUM2[2][2] = S1[0][2] * I00T4[0][2] + S1[1][2] * I00T4[1][2] + S1[2][2] * I00T4[2][2];
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				KS[i + 3][j + 6] = -DUM1[i][j] - DUM2[i][j];
				KS[j + 6][i + 3] = KS[i + 3][j + 6];
			}
		}

		//DUM1 = S2.transpose() * I00 * T3;
		DUM1[0][0] = S2[0][0] * I00T3[0][0] + S2[1][0] * I00T3[1][0] + S2[2][0] * I00T3[2][0];
		DUM1[0][1] = S2[0][0] * I00T3[0][1] + S2[1][0] * I00T3[1][1] + S2[2][0] * I00T3[2][1];
		DUM1[0][2] = S2[0][0] * I00T3[0][2] + S2[1][0] * I00T3[1][2] + S2[2][0] * I00T3[2][2];
		DUM1[1][0] = S2[0][1] * I00T3[0][0] + S2[1][1] * I00T3[1][0] + S2[2][1] * I00T3[2][0];
		DUM1[1][1] = S2[0][1] * I00T3[0][1] + S2[1][1] * I00T3[1][1] + S2[2][1] * I00T3[2][1];
		DUM1[1][2] = S2[0][1] * I00T3[0][2] + S2[1][1] * I00T3[1][2] + S2[2][1] * I00T3[2][2];
		DUM1[2][0] = /*S2[0][2] * I00T3[0][0] + S2[1][2] * I00T3[1][0]*/ +S2[2][2] * I00T3[2][0];
		DUM1[2][1] = /*S2[0][2] * I00T3[0][1] + S2[1][2] * I00T3[1][1]*/ +S2[2][2] * I00T3[2][1];
		DUM1[2][2] = /*S2[0][2] * I00T3[0][2] + S2[1][2] * I00T3[1][2]*/ +S2[2][2] * I00T3[2][2];
		//DUM2 = A1.transpose() * I00 * T4;
		DUM2[0][0] = /*A1[0][0] * I00T4[0][0] +*/ A1[1][0] * I00T4[1][0] + A1[2][0] * I00T4[2][0];
		DUM2[0][1] = /*A1[0][0] * I00T4[0][1] +*/ A1[1][0] * I00T4[1][1] + A1[2][0] * I00T4[2][1];
		DUM2[0][2] = /*A1[0][0] * I00T4[0][2] +*/ A1[1][0] * I00T4[1][2] + A1[2][0] * I00T4[2][2];
		DUM2[1][0] = A1[0][1] * I00T4[0][0] + /*A1[1][1] * I00T4[1][0] +*/ A1[2][1] * I00T4[2][0];
		DUM2[1][1] = A1[0][1] * I00T4[0][1] + /*A1[1][1] * I00T4[1][1] +*/ A1[2][1] * I00T4[2][1];
		DUM2[1][2] = A1[0][1] * I00T4[0][2] + /*A1[1][1] * I00T4[1][2] +*/ A1[2][1] * I00T4[2][2];
		DUM2[2][0] = A1[0][2] * I00T4[0][0] + A1[1][2] * I00T4[1][0] /*+ A1[2][2] * I00T4[2][0]*/;
		DUM2[2][1] = A1[0][2] * I00T4[0][1] + A1[1][2] * I00T4[1][1] /*+ A1[2][2] * I00T4[2][1]*/;
		DUM2[2][2] = A1[0][2] * I00T4[0][2] + A1[1][2] * I00T4[1][2] /*+ A1[2][2] * I00T4[2][2]*/;
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				KS[i + 6][j + 6] = DUM1[i][j] + DUM2[i][j];
			}
		}

		double BENSUM = 0.0, SHRSUM = 0.0;;
		for (int i = 3; i < 9; i++) { BENSUM += KB[i][i];		SHRSUM += KS[i][i]; }

		double PSI_HAT = BENSUM / SHRSUM;
		PHI_SQ = 2.0 * PSI_HAT / (1.0 + 2.0 * PSI_HAT);
	}

	void CQUAD4::STATIC_CONDENSATION(double KM_QQ_5[15][15], 
									 double B2M_QQ_5[3][15], 
									 double B3M_QQ_5[3][15], 
									 double k0[24][24])
	{
		//此函数将B矩阵由5节点转换为4节点，消除中心虚拟节点自由度
		//GOA = KM_QQ_5.inv()(13,14,15自由度) * KM_QQ_5(1`12自由度);			   //这里是刚度矩阵后三行，虚拟节点与其它节点的耦合项(3, 15)
		//KM_QQ_4 = KM_QQ_5(1`12自由度) - GOA.transpose * KM_QQ_5(13,14,15自由度)
		//B2M_QQ_4 = B2M_QQ_5(1`12自由度) - B2M_QQ_5(13,14,15自由度) * GOA
		//B3M_QQ_4 = B3M_QQ_5(1`12自由度) - B3M_QQ_5(13,14,15自由度) * GOA

		//GOA = KOO.inv() * KM_QQ_5(1`12自由度)
		Mat3 KOO, KOOi;		double GOA[3][12];
		KOO[0][0] = KM_QQ_5[12][12];	KOO[1][0] = KM_QQ_5[13][12];   KOO[2][0] = KM_QQ_5[14][12];
		KOO[0][1] = KM_QQ_5[12][13];	KOO[1][1] = KM_QQ_5[13][13];   KOO[2][1] = KM_QQ_5[14][13];
		KOO[0][2] = KM_QQ_5[12][14];	KOO[1][2] = KM_QQ_5[13][14];   KOO[2][2] = KM_QQ_5[14][14];
		KOOi = KOO.inv();

		GOA[0][0] = KOOi[0][0] * KM_QQ_5[12][0] + KOOi[0][1] * KM_QQ_5[13][0] + KOOi[0][2] * KM_QQ_5[14][0];
		GOA[0][1] = KOOi[0][0] * KM_QQ_5[12][1] + KOOi[0][1] * KM_QQ_5[13][1] + KOOi[0][2] * KM_QQ_5[14][1];
		GOA[0][2] = KOOi[0][0] * KM_QQ_5[12][2] + KOOi[0][1] * KM_QQ_5[13][2] + KOOi[0][2] * KM_QQ_5[14][2];
		GOA[0][3] = KOOi[0][0] * KM_QQ_5[12][3] + KOOi[0][1] * KM_QQ_5[13][3] + KOOi[0][2] * KM_QQ_5[14][3];
		GOA[0][4] = KOOi[0][0] * KM_QQ_5[12][4] + KOOi[0][1] * KM_QQ_5[13][4] + KOOi[0][2] * KM_QQ_5[14][4];
		GOA[0][5] = KOOi[0][0] * KM_QQ_5[12][5] + KOOi[0][1] * KM_QQ_5[13][5] + KOOi[0][2] * KM_QQ_5[14][5];
		GOA[0][6] = KOOi[0][0] * KM_QQ_5[12][6] + KOOi[0][1] * KM_QQ_5[13][6] + KOOi[0][2] * KM_QQ_5[14][6];
		GOA[0][7] = KOOi[0][0] * KM_QQ_5[12][7] + KOOi[0][1] * KM_QQ_5[13][7] + KOOi[0][2] * KM_QQ_5[14][7];
		GOA[0][8] = KOOi[0][0] * KM_QQ_5[12][8] + KOOi[0][1] * KM_QQ_5[13][8] + KOOi[0][2] * KM_QQ_5[14][8];
		GOA[0][9] = KOOi[0][0] * KM_QQ_5[12][9] + KOOi[0][1] * KM_QQ_5[13][9] + KOOi[0][2] * KM_QQ_5[14][9];
		GOA[0][10] = KOOi[0][0] * KM_QQ_5[12][10] + KOOi[0][1] * KM_QQ_5[13][10] + KOOi[0][2] * KM_QQ_5[14][10];
		GOA[0][11] = KOOi[0][0] * KM_QQ_5[12][11] + KOOi[0][1] * KM_QQ_5[13][11] + KOOi[0][2] * KM_QQ_5[14][11];

		GOA[1][0] = KOOi[1][0] * KM_QQ_5[12][0] + KOOi[1][1] * KM_QQ_5[13][0] + KOOi[1][2] * KM_QQ_5[14][0];
		GOA[1][1] = KOOi[1][0] * KM_QQ_5[12][1] + KOOi[1][1] * KM_QQ_5[13][1] + KOOi[1][2] * KM_QQ_5[14][1];
		GOA[1][2] = KOOi[1][0] * KM_QQ_5[12][2] + KOOi[1][1] * KM_QQ_5[13][2] + KOOi[1][2] * KM_QQ_5[14][2];
		GOA[1][3] = KOOi[1][0] * KM_QQ_5[12][3] + KOOi[1][1] * KM_QQ_5[13][3] + KOOi[1][2] * KM_QQ_5[14][3];
		GOA[1][4] = KOOi[1][0] * KM_QQ_5[12][4] + KOOi[1][1] * KM_QQ_5[13][4] + KOOi[1][2] * KM_QQ_5[14][4];
		GOA[1][5] = KOOi[1][0] * KM_QQ_5[12][5] + KOOi[1][1] * KM_QQ_5[13][5] + KOOi[1][2] * KM_QQ_5[14][5];
		GOA[1][6] = KOOi[1][0] * KM_QQ_5[12][6] + KOOi[1][1] * KM_QQ_5[13][6] + KOOi[1][2] * KM_QQ_5[14][6];
		GOA[1][7] = KOOi[1][0] * KM_QQ_5[12][7] + KOOi[1][1] * KM_QQ_5[13][7] + KOOi[1][2] * KM_QQ_5[14][7];
		GOA[1][8] = KOOi[1][0] * KM_QQ_5[12][8] + KOOi[1][1] * KM_QQ_5[13][8] + KOOi[1][2] * KM_QQ_5[14][8];
		GOA[1][9] = KOOi[1][0] * KM_QQ_5[12][9] + KOOi[1][1] * KM_QQ_5[13][9] + KOOi[1][2] * KM_QQ_5[14][9];
		GOA[1][10] = KOOi[1][0] * KM_QQ_5[12][10] + KOOi[1][1] * KM_QQ_5[13][10] + KOOi[1][2] * KM_QQ_5[14][10];
		GOA[1][11] = KOOi[1][0] * KM_QQ_5[12][11] + KOOi[1][1] * KM_QQ_5[13][11] + KOOi[1][2] * KM_QQ_5[14][11];

		GOA[2][0] = KOOi[2][0] * KM_QQ_5[12][0] + KOOi[2][1] * KM_QQ_5[13][0] + KOOi[2][2] * KM_QQ_5[14][0];
		GOA[2][1] = KOOi[2][0] * KM_QQ_5[12][1] + KOOi[2][1] * KM_QQ_5[13][1] + KOOi[2][2] * KM_QQ_5[14][1];
		GOA[2][2] = KOOi[2][0] * KM_QQ_5[12][2] + KOOi[2][1] * KM_QQ_5[13][2] + KOOi[2][2] * KM_QQ_5[14][2];
		GOA[2][3] = KOOi[2][0] * KM_QQ_5[12][3] + KOOi[2][1] * KM_QQ_5[13][3] + KOOi[2][2] * KM_QQ_5[14][3];
		GOA[2][4] = KOOi[2][0] * KM_QQ_5[12][4] + KOOi[2][1] * KM_QQ_5[13][4] + KOOi[2][2] * KM_QQ_5[14][4];
		GOA[2][5] = KOOi[2][0] * KM_QQ_5[12][5] + KOOi[2][1] * KM_QQ_5[13][5] + KOOi[2][2] * KM_QQ_5[14][5];
		GOA[2][6] = KOOi[2][0] * KM_QQ_5[12][6] + KOOi[2][1] * KM_QQ_5[13][6] + KOOi[2][2] * KM_QQ_5[14][6];
		GOA[2][7] = KOOi[2][0] * KM_QQ_5[12][7] + KOOi[2][1] * KM_QQ_5[13][7] + KOOi[2][2] * KM_QQ_5[14][7];
		GOA[2][8] = KOOi[2][0] * KM_QQ_5[12][8] + KOOi[2][1] * KM_QQ_5[13][8] + KOOi[2][2] * KM_QQ_5[14][8];
		GOA[2][9] = KOOi[2][0] * KM_QQ_5[12][9] + KOOi[2][1] * KM_QQ_5[13][9] + KOOi[2][2] * KM_QQ_5[14][9];
		GOA[2][10] = KOOi[2][0] * KM_QQ_5[12][10] + KOOi[2][1] * KM_QQ_5[13][10] + KOOi[2][2] * KM_QQ_5[14][10];
		GOA[2][11] = KOOi[2][0] * KM_QQ_5[12][11] + KOOi[2][1] * KM_QQ_5[13][11] + KOOi[2][2] * KM_QQ_5[14][11];

		//KM_QQ_4 = KM_QQ_5(1`12自由度) - GOA.transpose * KM_QQ_5(13,14,15自由度)	 因为KM_QQ_4会组装到k0中，这里不创建KM_QQ_4直接组装了
		for (int i = 0; i < 12; i++) {
			for (int j = 0; j < 12; j++) {
				k0[IDM[i]][IDM[j]] += KM_QQ_5[i][j] - (GOA[0][i] * KM_QQ_5[12][j] + GOA[1][i] * KM_QQ_5[13][j] + GOA[2][i] * KM_QQ_5[14][j]);
			}
		}

		//B2M_QQ_4 = B2M_QQ_5(1`12自由度,B2O) - B2M_QQ_5(13,14,15自由度,B2AB) * GOA
		BE2[0][0] = B2M_QQ_5[0][0] - (B2M_QQ_5[0][12] * GOA[0][0] + B2M_QQ_5[0][13] * GOA[1][0] + B2M_QQ_5[0][14] * GOA[2][0]);
		BE2[0][1] = B2M_QQ_5[0][1] - (B2M_QQ_5[0][12] * GOA[0][1] + B2M_QQ_5[0][13] * GOA[1][1] + B2M_QQ_5[0][14] * GOA[2][1]);
		BE2[0][2] = B2M_QQ_5[0][2] - (B2M_QQ_5[0][12] * GOA[0][2] + B2M_QQ_5[0][13] * GOA[1][2] + B2M_QQ_5[0][14] * GOA[2][2]);
		BE2[0][3] = B2M_QQ_5[0][3] - (B2M_QQ_5[0][12] * GOA[0][3] + B2M_QQ_5[0][13] * GOA[1][3] + B2M_QQ_5[0][14] * GOA[2][3]);
		BE2[0][4] = B2M_QQ_5[0][4] - (B2M_QQ_5[0][12] * GOA[0][4] + B2M_QQ_5[0][13] * GOA[1][4] + B2M_QQ_5[0][14] * GOA[2][4]);
		BE2[0][5] = B2M_QQ_5[0][5] - (B2M_QQ_5[0][12] * GOA[0][5] + B2M_QQ_5[0][13] * GOA[1][5] + B2M_QQ_5[0][14] * GOA[2][5]);
		BE2[0][6] = B2M_QQ_5[0][6] - (B2M_QQ_5[0][12] * GOA[0][6] + B2M_QQ_5[0][13] * GOA[1][6] + B2M_QQ_5[0][14] * GOA[2][6]);
		BE2[0][7] = B2M_QQ_5[0][7] - (B2M_QQ_5[0][12] * GOA[0][7] + B2M_QQ_5[0][13] * GOA[1][7] + B2M_QQ_5[0][14] * GOA[2][7]);
		BE2[0][8] = B2M_QQ_5[0][8] - (B2M_QQ_5[0][12] * GOA[0][8] + B2M_QQ_5[0][13] * GOA[1][8] + B2M_QQ_5[0][14] * GOA[2][8]);
		BE2[0][9] = B2M_QQ_5[0][9] - (B2M_QQ_5[0][12] * GOA[0][9] + B2M_QQ_5[0][13] * GOA[1][9] + B2M_QQ_5[0][14] * GOA[2][9]);
		BE2[0][10] = B2M_QQ_5[0][10] - (B2M_QQ_5[0][12] * GOA[0][10] + B2M_QQ_5[0][13] * GOA[1][10] + B2M_QQ_5[0][14] * GOA[2][10]);
		BE2[0][11] = B2M_QQ_5[0][11] - (B2M_QQ_5[0][12] * GOA[0][11] + B2M_QQ_5[0][13] * GOA[1][11] + B2M_QQ_5[0][14] * GOA[2][11]);

		BE2[1][0] = B2M_QQ_5[1][0] - (B2M_QQ_5[1][12] * GOA[0][0] + B2M_QQ_5[1][13] * GOA[1][0] + B2M_QQ_5[1][14] * GOA[2][0]);
		BE2[1][1] = B2M_QQ_5[1][1] - (B2M_QQ_5[1][12] * GOA[0][1] + B2M_QQ_5[1][13] * GOA[1][1] + B2M_QQ_5[1][14] * GOA[2][1]);
		BE2[1][2] = B2M_QQ_5[1][2] - (B2M_QQ_5[1][12] * GOA[0][2] + B2M_QQ_5[1][13] * GOA[1][2] + B2M_QQ_5[1][14] * GOA[2][2]);
		BE2[1][3] = B2M_QQ_5[1][3] - (B2M_QQ_5[1][12] * GOA[0][3] + B2M_QQ_5[1][13] * GOA[1][3] + B2M_QQ_5[1][14] * GOA[2][3]);
		BE2[1][4] = B2M_QQ_5[1][4] - (B2M_QQ_5[1][12] * GOA[0][4] + B2M_QQ_5[1][13] * GOA[1][4] + B2M_QQ_5[1][14] * GOA[2][4]);
		BE2[1][5] = B2M_QQ_5[1][5] - (B2M_QQ_5[1][12] * GOA[0][5] + B2M_QQ_5[1][13] * GOA[1][5] + B2M_QQ_5[1][14] * GOA[2][5]);
		BE2[1][6] = B2M_QQ_5[1][6] - (B2M_QQ_5[1][12] * GOA[0][6] + B2M_QQ_5[1][13] * GOA[1][6] + B2M_QQ_5[1][14] * GOA[2][6]);
		BE2[1][7] = B2M_QQ_5[1][7] - (B2M_QQ_5[1][12] * GOA[0][7] + B2M_QQ_5[1][13] * GOA[1][7] + B2M_QQ_5[1][14] * GOA[2][7]);
		BE2[1][8] = B2M_QQ_5[1][8] - (B2M_QQ_5[1][12] * GOA[0][8] + B2M_QQ_5[1][13] * GOA[1][8] + B2M_QQ_5[1][14] * GOA[2][8]);
		BE2[1][9] = B2M_QQ_5[1][9] - (B2M_QQ_5[1][12] * GOA[0][9] + B2M_QQ_5[1][13] * GOA[1][9] + B2M_QQ_5[1][14] * GOA[2][9]);
		BE2[1][10] = B2M_QQ_5[1][10] - (B2M_QQ_5[1][12] * GOA[0][10] + B2M_QQ_5[1][13] * GOA[1][10] + B2M_QQ_5[1][14] * GOA[2][10]);
		BE2[1][11] = B2M_QQ_5[1][11] - (B2M_QQ_5[1][12] * GOA[0][11] + B2M_QQ_5[1][13] * GOA[1][11] + B2M_QQ_5[1][14] * GOA[2][11]);

		BE2[2][0] = B2M_QQ_5[2][0] - (B2M_QQ_5[2][12] * GOA[0][0] + B2M_QQ_5[2][13] * GOA[1][0] + B2M_QQ_5[2][14] * GOA[2][0]);
		BE2[2][1] = B2M_QQ_5[2][1] - (B2M_QQ_5[2][12] * GOA[0][1] + B2M_QQ_5[2][13] * GOA[1][1] + B2M_QQ_5[2][14] * GOA[2][1]);
		BE2[2][2] = B2M_QQ_5[2][2] - (B2M_QQ_5[2][12] * GOA[0][2] + B2M_QQ_5[2][13] * GOA[1][2] + B2M_QQ_5[2][14] * GOA[2][2]);
		BE2[2][3] = B2M_QQ_5[2][3] - (B2M_QQ_5[2][12] * GOA[0][3] + B2M_QQ_5[2][13] * GOA[1][3] + B2M_QQ_5[2][14] * GOA[2][3]);
		BE2[2][4] = B2M_QQ_5[2][4] - (B2M_QQ_5[2][12] * GOA[0][4] + B2M_QQ_5[2][13] * GOA[1][4] + B2M_QQ_5[2][14] * GOA[2][4]);
		BE2[2][5] = B2M_QQ_5[2][5] - (B2M_QQ_5[2][12] * GOA[0][5] + B2M_QQ_5[2][13] * GOA[1][5] + B2M_QQ_5[2][14] * GOA[2][5]);
		BE2[2][6] = B2M_QQ_5[2][6] - (B2M_QQ_5[2][12] * GOA[0][6] + B2M_QQ_5[2][13] * GOA[1][6] + B2M_QQ_5[2][14] * GOA[2][6]);
		BE2[2][7] = B2M_QQ_5[2][7] - (B2M_QQ_5[2][12] * GOA[0][7] + B2M_QQ_5[2][13] * GOA[1][7] + B2M_QQ_5[2][14] * GOA[2][7]);
		BE2[2][8] = B2M_QQ_5[2][8] - (B2M_QQ_5[2][12] * GOA[0][8] + B2M_QQ_5[2][13] * GOA[1][8] + B2M_QQ_5[2][14] * GOA[2][8]);
		BE2[2][9] = B2M_QQ_5[2][9] - (B2M_QQ_5[2][12] * GOA[0][9] + B2M_QQ_5[2][13] * GOA[1][9] + B2M_QQ_5[2][14] * GOA[2][9]);
		BE2[2][10] = B2M_QQ_5[2][10] - (B2M_QQ_5[2][12] * GOA[0][10] + B2M_QQ_5[2][13] * GOA[1][10] + B2M_QQ_5[2][14] * GOA[2][10]);
		BE2[2][11] = B2M_QQ_5[2][11] - (B2M_QQ_5[2][12] * GOA[0][11] + B2M_QQ_5[2][13] * GOA[1][11] + B2M_QQ_5[2][14] * GOA[2][11]);

		//B3M_QQ_4 = B3M_QQ_5(1`12自由度) - B3M_QQ_5(13,14,15自由度) * GOA
		BE3[0][0] = B3M_QQ_5[0][0] - (B3M_QQ_5[0][12] * GOA[0][0] + B3M_QQ_5[0][13] * GOA[1][0] + B3M_QQ_5[0][14] * GOA[2][0]);
		BE3[0][1] = B3M_QQ_5[0][1] - (B3M_QQ_5[0][12] * GOA[0][1] + B3M_QQ_5[0][13] * GOA[1][1] + B3M_QQ_5[0][14] * GOA[2][1]);
		BE3[0][2] = B3M_QQ_5[0][2] - (B3M_QQ_5[0][12] * GOA[0][2] + B3M_QQ_5[0][13] * GOA[1][2] + B3M_QQ_5[0][14] * GOA[2][2]);
		BE3[0][3] = B3M_QQ_5[0][3] - (B3M_QQ_5[0][12] * GOA[0][3] + B3M_QQ_5[0][13] * GOA[1][3] + B3M_QQ_5[0][14] * GOA[2][3]);
		BE3[0][4] = B3M_QQ_5[0][4] - (B3M_QQ_5[0][12] * GOA[0][4] + B3M_QQ_5[0][13] * GOA[1][4] + B3M_QQ_5[0][14] * GOA[2][4]);
		BE3[0][5] = B3M_QQ_5[0][5] - (B3M_QQ_5[0][12] * GOA[0][5] + B3M_QQ_5[0][13] * GOA[1][5] + B3M_QQ_5[0][14] * GOA[2][5]);
		BE3[0][6] = B3M_QQ_5[0][6] - (B3M_QQ_5[0][12] * GOA[0][6] + B3M_QQ_5[0][13] * GOA[1][6] + B3M_QQ_5[0][14] * GOA[2][6]);
		BE3[0][7] = B3M_QQ_5[0][7] - (B3M_QQ_5[0][12] * GOA[0][7] + B3M_QQ_5[0][13] * GOA[1][7] + B3M_QQ_5[0][14] * GOA[2][7]);
		BE3[0][8] = B3M_QQ_5[0][8] - (B3M_QQ_5[0][12] * GOA[0][8] + B3M_QQ_5[0][13] * GOA[1][8] + B3M_QQ_5[0][14] * GOA[2][8]);
		BE3[0][9] = B3M_QQ_5[0][9] - (B3M_QQ_5[0][12] * GOA[0][9] + B3M_QQ_5[0][13] * GOA[1][9] + B3M_QQ_5[0][14] * GOA[2][9]);
		BE3[0][10] = B3M_QQ_5[0][10] - (B3M_QQ_5[0][12] * GOA[0][10] + B3M_QQ_5[0][13] * GOA[1][10] + B3M_QQ_5[0][14] * GOA[2][10]);
		BE3[0][11] = B3M_QQ_5[0][11] - (B3M_QQ_5[0][12] * GOA[0][11] + B3M_QQ_5[0][13] * GOA[1][11] + B3M_QQ_5[0][14] * GOA[2][11]);

		BE3[1][0] = B3M_QQ_5[1][0] - (B3M_QQ_5[1][12] * GOA[0][0] + B3M_QQ_5[1][13] * GOA[1][0] + B3M_QQ_5[1][14] * GOA[2][0]);
		BE3[1][1] = B3M_QQ_5[1][1] - (B3M_QQ_5[1][12] * GOA[0][1] + B3M_QQ_5[1][13] * GOA[1][1] + B3M_QQ_5[1][14] * GOA[2][1]);
		BE3[1][2] = B3M_QQ_5[1][2] - (B3M_QQ_5[1][12] * GOA[0][2] + B3M_QQ_5[1][13] * GOA[1][2] + B3M_QQ_5[1][14] * GOA[2][2]);
		BE3[1][3] = B3M_QQ_5[1][3] - (B3M_QQ_5[1][12] * GOA[0][3] + B3M_QQ_5[1][13] * GOA[1][3] + B3M_QQ_5[1][14] * GOA[2][3]);
		BE3[1][4] = B3M_QQ_5[1][4] - (B3M_QQ_5[1][12] * GOA[0][4] + B3M_QQ_5[1][13] * GOA[1][4] + B3M_QQ_5[1][14] * GOA[2][4]);
		BE3[1][5] = B3M_QQ_5[1][5] - (B3M_QQ_5[1][12] * GOA[0][5] + B3M_QQ_5[1][13] * GOA[1][5] + B3M_QQ_5[1][14] * GOA[2][5]);
		BE3[1][6] = B3M_QQ_5[1][6] - (B3M_QQ_5[1][12] * GOA[0][6] + B3M_QQ_5[1][13] * GOA[1][6] + B3M_QQ_5[1][14] * GOA[2][6]);
		BE3[1][7] = B3M_QQ_5[1][7] - (B3M_QQ_5[1][12] * GOA[0][7] + B3M_QQ_5[1][13] * GOA[1][7] + B3M_QQ_5[1][14] * GOA[2][7]);
		BE3[1][8] = B3M_QQ_5[1][8] - (B3M_QQ_5[1][12] * GOA[0][8] + B3M_QQ_5[1][13] * GOA[1][8] + B3M_QQ_5[1][14] * GOA[2][8]);
		BE3[1][9] = B3M_QQ_5[1][9] - (B3M_QQ_5[1][12] * GOA[0][9] + B3M_QQ_5[1][13] * GOA[1][9] + B3M_QQ_5[1][14] * GOA[2][9]);
		BE3[1][10] = B3M_QQ_5[1][10] - (B3M_QQ_5[1][12] * GOA[0][10] + B3M_QQ_5[1][13] * GOA[1][10] + B3M_QQ_5[1][14] * GOA[2][10]);
		BE3[1][11] = B3M_QQ_5[1][11] - (B3M_QQ_5[1][12] * GOA[0][11] + B3M_QQ_5[1][13] * GOA[1][11] + B3M_QQ_5[1][14] * GOA[2][11]);

		BE3[2][0] = 0.0;//B3M_QQ_5[2][0] - (B3M_QQ_5[2][12] * GOA[0][0] + B3M_QQ_5[2][13] * GOA[1][0] + B3M_QQ_5[2][14] * GOA[2][0]);
		BE3[2][1] = 0.0;//B3M_QQ_5[2][1] - (B3M_QQ_5[2][12] * GOA[0][1] + B3M_QQ_5[2][13] * GOA[1][1] + B3M_QQ_5[2][14] * GOA[2][1]);
		BE3[2][2] = 0.0;//B3M_QQ_5[2][2] - (B3M_QQ_5[2][12] * GOA[0][2] + B3M_QQ_5[2][13] * GOA[1][2] + B3M_QQ_5[2][14] * GOA[2][2]);
		BE3[2][3] = 0.0;//B3M_QQ_5[2][3] - (B3M_QQ_5[2][12] * GOA[0][3] + B3M_QQ_5[2][13] * GOA[1][3] + B3M_QQ_5[2][14] * GOA[2][3]);
		BE3[2][4] = 0.0;//B3M_QQ_5[2][4] - (B3M_QQ_5[2][12] * GOA[0][4] + B3M_QQ_5[2][13] * GOA[1][4] + B3M_QQ_5[2][14] * GOA[2][4]);
		BE3[2][5] = 0.0;//B3M_QQ_5[2][5] - (B3M_QQ_5[2][12] * GOA[0][5] + B3M_QQ_5[2][13] * GOA[1][5] + B3M_QQ_5[2][14] * GOA[2][5]);
		BE3[2][6] = 0.0;//B3M_QQ_5[2][6] - (B3M_QQ_5[2][12] * GOA[0][6] + B3M_QQ_5[2][13] * GOA[1][6] + B3M_QQ_5[2][14] * GOA[2][6]);
		BE3[2][7] = 0.0;//B3M_QQ_5[2][7] - (B3M_QQ_5[2][12] * GOA[0][7] + B3M_QQ_5[2][13] * GOA[1][7] + B3M_QQ_5[2][14] * GOA[2][7]);
		BE3[2][8] = 0.0;//B3M_QQ_5[2][8] - (B3M_QQ_5[2][12] * GOA[0][8] + B3M_QQ_5[2][13] * GOA[1][8] + B3M_QQ_5[2][14] * GOA[2][8]);
		BE3[2][9] = 0.0;//B3M_QQ_5[2][9] - (B3M_QQ_5[2][12] * GOA[0][9] + B3M_QQ_5[2][13] * GOA[1][9] + B3M_QQ_5[2][14] * GOA[2][9]);
		BE3[2][10] = 0.0;//B3M_QQ_5[2][10] - (B3M_QQ_5[2][12] * GOA[0][10] + B3M_QQ_5[2][13] * GOA[1][10] + B3M_QQ_5[2][14] * GOA[2][10]);
		BE3[2][11] = 0.0;//B3M_QQ_5[2][11] - (B3M_QQ_5[2][12] * GOA[0][11] + B3M_QQ_5[2][13] * GOA[1][11] + B3M_QQ_5[2][14] * GOA[2][11]);

	}

	void CQUAD4::Stress(vector<double>& UEL)
	{

		double matD[8][8] = { 0. };
		build_cons_mat(matD);

		double STRESS[9] = { 0. };
		double STRAIN[3][3] = { 0. };

		STRAIN[0][0] = BE1[0][0] * UEL[0] + BE1[0][1] * UEL[6] + BE1[0][2] * UEL[12] + BE1[0][3] * UEL[18];
		STRAIN[0][1] = BE1[1][0] * UEL[1] + BE1[1][1] * UEL[7] + BE1[1][2] * UEL[13] + BE1[1][3] * UEL[19];
		STRAIN[0][2] = BE1[1][0] * UEL[0] + BE1[0][0] * UEL[1] + BE1[1][1] * UEL[6] + BE1[0][1] * UEL[7]
			+ BE1[1][2] * UEL[12] + BE1[0][2] * UEL[13] + BE1[1][3] * UEL[18] + BE1[0][3] * UEL[19];

		STRAIN[1][0] = BE2[0][0] * UEL[2] + BE2[0][1] * UEL[3] + BE2[0][2] * UEL[4] + BE2[0][3] * UEL[8]
			+ BE2[0][4] * UEL[9] + BE2[0][5] * UEL[10] + BE2[0][6] * UEL[14] + BE2[0][7] * UEL[15]
			+ BE2[0][8] * UEL[16] + BE2[0][9] * UEL[20] + BE2[0][10] * UEL[21] + BE2[0][11] * UEL[22];
		STRAIN[1][1] = BE2[1][0] * UEL[2] + BE2[1][1] * UEL[3] + BE2[1][2] * UEL[4] + BE2[1][3] * UEL[8]
			+ BE2[1][4] * UEL[9] + BE2[1][5] * UEL[10] + BE2[1][6] * UEL[14] + BE2[1][7] * UEL[15]
			+ BE2[1][8] * UEL[16] + BE2[1][9] * UEL[20] + BE2[1][10] * UEL[21] + BE2[1][11] * UEL[22];
		STRAIN[1][2] = BE2[2][0] * UEL[2] + BE2[2][1] * UEL[3] + BE2[2][2] * UEL[4] + BE2[2][3] * UEL[8]
			+ BE2[2][4] * UEL[9] + BE2[2][5] * UEL[10] + BE2[2][6] * UEL[14] + BE2[2][7] * UEL[15]
			+ BE2[2][8] * UEL[16] + BE2[2][9] * UEL[20] + BE2[2][10] * UEL[21] + BE2[2][11] * UEL[22];

		STRAIN[2][0] = BE3[0][0] * UEL[2] + BE3[0][1] * UEL[3] + BE3[0][2] * UEL[4] + BE3[0][3] * UEL[8]
			+ BE3[0][4] * UEL[9] + BE3[0][5] * UEL[10] + BE3[0][6] * UEL[14] + BE3[0][7] * UEL[15]
			+ BE3[0][8] * UEL[16] + BE3[0][9] * UEL[20] + BE3[0][10] * UEL[21] + BE3[0][11] * UEL[22];
		STRAIN[2][1] = BE3[1][0] * UEL[2] + BE3[1][1] * UEL[3] + BE3[1][2] * UEL[4] + BE3[1][3] * UEL[8]
			+ BE3[1][4] * UEL[9] + BE3[1][5] * UEL[10] + BE3[1][6] * UEL[14] + BE3[1][7] * UEL[15]
			+ BE3[1][8] * UEL[16] + BE3[1][9] * UEL[20] + BE3[1][10] * UEL[21] + BE3[1][11] * UEL[22];
		STRAIN[2][2] = BE3[2][0] * UEL[2] + BE3[2][1] * UEL[3] + BE3[2][2] * UEL[4] + BE3[2][3] * UEL[8]
			+ BE3[2][4] * UEL[9] + BE3[2][5] * UEL[10] + BE3[2][6] * UEL[14] + BE3[2][7] * UEL[15]
			+ BE3[2][8] * UEL[16] + BE3[2][9] * UEL[20] + BE3[2][10] * UEL[21] + BE3[2][11] * UEL[22];

		double EM[3][3];
		EM[0][0] = matD[0][0] / t;			EM[0][1] = matD[0][1] / t;			EM[0][2] = matD[0][2] / t;
		EM[1][0] = matD[1][0] / t;			EM[1][1] = matD[1][1] / t;			EM[1][2] = matD[1][2] / t;
		EM[2][0] = matD[2][0] / t;			EM[2][1] = matD[2][1] / t;			EM[2][2] = matD[2][2] / t;

		STRESS[0] = EM[0][0] * STRAIN[0][0] + EM[0][1] * STRAIN[0][1] + EM[0][2] * STRAIN[0][2];
		STRESS[1] = EM[1][0] * STRAIN[0][0] + EM[1][1] * STRAIN[0][1] + EM[1][2] * STRAIN[0][2];
		STRESS[2] = EM[2][0] * STRAIN[0][0] + EM[2][1] * STRAIN[0][1] + EM[2][2] * STRAIN[0][2];

		STRESS[3] = EM[0][0] * STRAIN[1][0] + EM[0][1] * STRAIN[1][1] + EM[0][2] * STRAIN[1][2];
		STRESS[4] = EM[1][0] * STRAIN[1][0] + EM[1][1] * STRAIN[1][1] + EM[1][2] * STRAIN[1][2];
		STRESS[5] = EM[2][0] * STRAIN[1][0] + EM[2][1] * STRAIN[1][1] + EM[2][2] * STRAIN[1][2];

		STRESS[6] = EM[2][2] * STRAIN[2][0] * PHI_SQ;
		STRESS[7] = EM[2][2] * STRAIN[2][1] * PHI_SQ;
		STRESS[8] = 0.0;
	}

}