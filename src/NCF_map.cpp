/**************************************************************************

Copyright:  WH team

Author: Zhumingjun <1765380405@qq.com>

Completion date:  XXX

Description: XXX

**************************************************************************/

#pragma once
#include "include/NCF_map.h"



namespace CAE
{
	

	void NCF_map::PhySpaceGPs(data_management& data_cae, elastic_mat& data_mat)
	{
		int nF_bmesh = data_cae.BndMesh_F.size();
		MatrixXd nodes1;//交界面(三角形、四边形）节点编号
		//MatrixXd pts1, pts2;//交界处粗、细单元节点坐标
		
		//*******点在平面*********
		Point A(-1.1059770001, 0, 0.25);
		Point B(-1.1059770002, 0.1875000, 0.375000000);
		Point C(-1.1059770005, 0, 0.50);

		Point D(-1.105977, 0.062500000, 0.375000000);

		double s = CalculateArea(A, B, C);

		double s1 = CalculateArea(D, A, B);
		double s2 = CalculateArea(D, B, C);
		double s3 = CalculateArea(D, C, A);
		
		double b = (s1+s2+s3) / s;
		
		
		double c=s;
		



		

		//交界面单元循环
		for (int e = 0; e < nF_bmesh; e++) //交界面处细网格个数
		{   
			int id_ele_F = data_cae.BndMesh_F[e] - 1;//索引-1
			//int id_ele_C = data_cae.BndMesh_C[e] - 1;

			int face_nodes_ = data_cae.ele_list_[data_cae.ele_list_idx_[id_ele_F]]->face_node;
			GetIntF_face_Inform(data_cae, nodes1, face_nodes_, e);
			
			//交界处粗、细单元节点个数
			//int n_node_F = data_cae.ele_list_[data_cae.ele_list_idx_[id_ele_F]]->nnode_;
			//int n_node_C = data_cae.ele_list_[data_cae.ele_list_idx_[id_ele_C]]->nnode_;
			//pts1.resize(n_node_F, 3);
			//pts2.resize(n_node_C, 3);
			//GetIntF_ele_Inform(data_cae, pts1, pts2, e);

			//一个面上的高斯积分点个数
			//int face_gps_ = data_cae.ele_list_[data_cae.ele_list_idx_[id_ele_F]]->face_gps;
			//储存一个面上的高斯积分点物理坐标
			//MatrixXd phy_gps(face_gps_,3);
			//MatrixXd X(face_gps_,4);//权重和法向量
			
			 // 计算交界面积分点物理空间坐标
			data_cae.ele_list_[data_cae.ele_list_idx_[id_ele_F]]->text_gps_phy_coords(nodes1, text_gps,text_W_1,text_Normal);
				
			//高斯积分点循环 4组  
			
			//for (int q = 0; q < face_gps_; q++)
			//{
			//	
			//	//积分点物理空间坐标
			//	MatrixXd xx = phy_gps.row(q);

			//	//积分点分别向父空间映射
			//	
			//	MatrixXd X1, X2;
			//	//-----------------------------------测试
			//	std::cout << "F_mesh number :" << data_cae.BndMesh_F[e] << std::endl;
			//	//X1 = GlobalMap3D(xx, pts1);//细网格
			//	std::cout << "C_mesh number :" << data_cae.BndMesh_C[e] << std::endl;
			//	//X2 = GlobalMap3D(xx, pts2);//粗网格
			//	//-----------------------------------
			//	
			//	//P_GP1.push_back(X1);
			//	//P_GP2.push_back(X2);
			//	

			//}
		}
		//int a = 10;

		
		//需要找到积分点对应的粗网格单元编号

		int n_gps=text_gps.size();

		for (int i = 0; i < n_gps; i++)
		{
			Point G(text_gps[i][0], text_gps[i][1], text_gps[i][2]);








		}



	}
	

	double NCF_map::CalculateArea(const Point& A, const Point& B, const Point& C)
	{
		// Calculate vectors AB and AC
		double ABx = B.x - A.x;
		double ABy = B.y - A.y;
		double ABz = B.z - A.z;

		double ACx = C.x - A.x;
		double ACy = C.y - A.y;
		double ACz = C.z - A.z;

		// Calculate cross product of vectors AB and AC
		double crossX = ABy * ACz - ABz * ACy;
		double crossY = ABz * ACx - ABx * ACz;
		double crossZ = ABx * ACy - ABy * ACx;

		// Calculate magnitude of cross product and divide by 2 to get triangle area
		double area = 0.5 * sqrt(crossX * crossX + crossY * crossY + crossZ * crossZ);
		return area;
	}


	//获取交界面处（四边形、三角形）面节点信息
	void NCF_map::GetIntF_face_Inform(data_management& data_cae, MatrixXd& nodes1,
		int& face_nodes,int& e )
	{   
		nodes1.resize(face_nodes, 3);
		for (int i = 0; i < face_nodes; ++i)
		{
			int node_index = data_cae.bndFace_finemesh[e][i] - 1;
			for (int j = 0; j < 3; ++j)
			{
				nodes1(i, j) = data_cae.coords_[node_index][j];
			}

		}
		

	}

	//void NCF_map::GetIntF_ele_Inform(data_management& data_cae,MatrixXd& pts1,
	//	MatrixXd& pts2, int& e)
	//{
	//	const auto& tps_F = data_cae.node_topos_[data_cae.BndMesh_F[e]-1];//索引要减1
	//	const auto& tps_C = data_cae.node_topos_[data_cae.BndMesh_C[e]-1];

	//	for (int i = 0; i < 8; ++i)
	//	{
	//		const auto& coord_sctr1 = data_cae.coords_[tps_F[i] - 1];
	//		const auto& coord_sctr2 = data_cae.coords_[tps_C[i] - 1];

	//		for (int j = 0; j < 3; ++j)
	//		{
	//			pts1(i, j) = coord_sctr1[j];
	//			pts2(i, j) = coord_sctr2[j];
	//		}
	//	}
	//	
	//}

	/*
	void NCF_map::GetIntF_ele_Inform(data_management& data_cae, vector<int>& sctr1, vector<int>& sctr2,
		MatrixXd& pts1, MatrixXd& pts2, int& e)
	{
		sctr1[0] = data_cae.node_topos_[data_cae.BndMesh_F[e]][0]; 
		sctr1[1] = data_cae.node_topos_[data_cae.BndMesh_F[e]][1];
		sctr1[2] = data_cae.node_topos_[data_cae.BndMesh_F[e]][2];
		sctr1[3] = data_cae.node_topos_[data_cae.BndMesh_F[e]][3];
		sctr1[4] = data_cae.node_topos_[data_cae.BndMesh_F[e]][4];
		sctr1[5] = data_cae.node_topos_[data_cae.BndMesh_F[e]][5];
		sctr1[6] = data_cae.node_topos_[data_cae.BndMesh_F[e]][6];
		sctr1[7] = data_cae.node_topos_[data_cae.BndMesh_F[e]][7];
		sctr1[8] = data_cae.node_topos_[data_cae.BndMesh_F[e]][8];
		
		sctr2[0] = data_cae.node_topos_[data_cae.BndMesh_C[e]][0];
		sctr2[1] = data_cae.node_topos_[data_cae.BndMesh_C[e]][1];
		sctr2[2] = data_cae.node_topos_[data_cae.BndMesh_C[e]][2];
		sctr2[3] = data_cae.node_topos_[data_cae.BndMesh_C[e]][3];
		sctr2[4] = data_cae.node_topos_[data_cae.BndMesh_C[e]][4];
		sctr2[5] = data_cae.node_topos_[data_cae.BndMesh_C[e]][5];
		sctr2[6] = data_cae.node_topos_[data_cae.BndMesh_C[e]][6];
		sctr2[7] = data_cae.node_topos_[data_cae.BndMesh_C[e]][7];
		

		for (int i = 0; i < 8; i++) //交界处六面体单元节点坐标
		{
			pts1(i, 0) = data_cae.coords_[sctr1[i]][0];
			pts1(i, 1) = data_cae.coords_[sctr1[i]][1];
			pts1(i, 2) = data_cae.coords_[sctr1[i]][2];
			
			pts2(i, 0) = data_cae.coords_[sctr2[i]][0];
			pts2(i, 1) = data_cae.coords_[sctr2[i]][1];
			pts2(i, 2) = data_cae.coords_[sctr2[i]][2];

			
		}
	}
	*/
	


	void NCF_map::InterfacialStifMatrix(data_management& data_cae, elastic_mat& data_mat,
		vector<double>& nz_val, vector<int>& row_idx, vector<int>& col_idx)
	{

		
		MatrixXd Ce;
		Get_Ce(data_mat,Ce);
		MatrixXd pts1, pts2;//交界处细、粗六面体单元节点坐标
		vector<int> F_eper_dof, C_eper_dof; // 粗细单元节点自由度
	    //初始化B,Nm,耦合矩阵,界面总刚
		MatrixXd B1, B2, Nm1, Nm2, Kp11, Kp12, Kp22,Kd11, Kd12, Kd21, Kd22,K11, K12, K21, K22;
		MatrixXd nodes1;//细网格交界面(三角形、四边形）节点编号
		MatrixXd phy_gps;//储存一个面上的高斯积分点物理坐标
		vector<double> W_1;//权重
		vector<Eigen::Vector3d> Normal;//法向量
		//法向量矩阵
	    MatrixXd n(3, 6);
		int nF_bmesh = data_cae.BndMesh_F.size();//交界面处细网格个数
		for (int e = 0; e < nF_bmesh; e++) //每一组"粗细网格单元的界面刚度矩阵"计算
		{  
			 int id_ele_F = data_cae.BndMesh_F[e] - 1;//索引-1
			 int id_ele_C = data_cae.BndMesh_C[e] - 1;
			 int face_nodes_ = data_cae.ele_list_[data_cae.ele_list_idx_[id_ele_F]]->face_node;
			 GetIntF_face_Inform(data_cae, nodes1, face_nodes_, e);

			 //交界处粗、细单元节点个数
			 int n_node_F = data_cae.ele_list_[data_cae.ele_list_idx_[id_ele_F]]->nnode_;
			 int n_node_C = data_cae.ele_list_[data_cae.ele_list_idx_[id_ele_C]]->nnode_;
			 pts1.resize(n_node_F, 3);
			 pts2.resize(n_node_C, 3);
			 //获取"粗细单元"节点信息
			 GetIntF_ele_Inform(pts1, pts2, F_eper_dof, C_eper_dof, data_cae,
				 e, n_node_F, n_node_C);
			 //初始化每一组“粗细单元”的矩阵
			 Kp11.resize(n_node_F * 3, n_node_F * 3); Kp12.resize(n_node_F * 3, n_node_C * 3);
			 Kp22.resize(n_node_C * 3, n_node_C * 3);
			 Kd11.resize(n_node_F * 3, n_node_F * 3);Kd12.resize(n_node_F * 3, n_node_C * 3);
			 Kd21.resize(n_node_C * 3, n_node_F * 3);Kd22.resize(n_node_C * 3, n_node_C * 3);
			 B1.resize(6, n_node_F * 3); Nm1.resize(3, n_node_F * 3);
			 B2.resize(6, n_node_C * 3); Nm2.resize(3, n_node_C * 3);
			 Kp11.setZero(); Kp12.setZero(); Kp22.setZero();
			 Kd11.setZero(); Kd12.setZero(); Kd21.setZero(); Kd22.setZero();
			 
			 //一个面上的高斯积分点个数
			 int face_gps_ = data_cae.ele_list_[data_cae.ele_list_idx_[id_ele_F]]->face_gps;
			 phy_gps.resize(face_gps_, 3);
			 W_1.resize(face_gps_);//权重
			 Normal.resize(face_gps_);//法向量
			 // 计算交界面积分点物理空间坐标、权重、法向量
			 data_cae.ele_list_[data_cae.ele_list_idx_[id_ele_F]]->gps_phy_coords(nodes1, phy_gps, W_1, Normal);
			 
			 
			 //text****************
			 vector<double> p(3);
			 vector<vector<double>> np(3, vector<double>(3));
			
			 p[0] = phy_gps(0, 0);
			 p[1] = phy_gps(0, 1);
			 p[2] = phy_gps(0, 2);
			 for (int i = 0; i < 3; i++)
			 {
				 for (int j = 0; j < 3; j++)
				 {
					 np[i][j] = nodes1(i, j);
				 }

			 }
			 //*************


			 //积分点循环 
			 for (int q = 0; q < face_gps_; q++)
			 {
				 //积分点物理空间坐标
				 MatrixXd xx = phy_gps.row(q);

				 //判断积分点位于哪一个粗网格面上


				 //积分点分别向"粗细网格"父空间映射
				 MatrixXd X1, X2; //"粗细网格"积分点父空间坐标
				 //-----------------------------------测试
				 std::cout << "F_mesh number :" << data_cae.BndMesh_F[e] << std::endl;
				 X1 = GlobalMap3D(xx, pts1, n_node_F);//细网格
				 std::cout << "C_mesh number :" << data_cae.BndMesh_C[e] << std::endl;
				 X2 = GlobalMap3D(xx, pts2, n_node_C);//粗网格
				 //-----------------------------------

				 Calculate_B_Nm(X1, pts1, B1, Nm1, n_node_F);
				 Calculate_B_Nm(X2, pts2, B2, Nm2, n_node_C);

				 n.setZero();//交界面细网格法向量方向
				 n(0, 0) = Normal[q](0);  n(0, 3) = Normal[q](1); n(0, 5) = Normal[q](2);
				 n(1, 1) = Normal[q](1);  n(1, 3) = Normal[q](0); n(1, 4) = Normal[q](2);
				 n(2, 2) = Normal[q](2);  n(2, 4) = Normal[q](1); n(2, 5) = Normal[q](0);
				 //权重
				 double wt1 = W_1[q];
				 //计算罚参数
				 double alpha = Get_Alpha(data_cae, data_mat, pts1);
				 //计算Kp,Kd
				 Calculate_Kp_Kd(Kp11, Kp12, Kp22, Kd11, Kd12, Kd21, Kd22,
					 Nm1, Nm2, B1, B2, n, Ce, alpha, wt1);

			 }

			 //计算界面总刚
			 Calculate_InterFMatrix(K11, K12, K21, K22, Kd11, Kd12, Kd21, Kd22,
				 Kp11, Kp12, Kp22, n_node_F, n_node_C);
			 //界面刚度矩阵储存为CSR
			 Fill_InterFMatrix(F_eper_dof, F_eper_dof, nz_val, row_idx, col_idx, K11);
			 Fill_InterFMatrix(F_eper_dof, C_eper_dof, nz_val, row_idx, col_idx, K12);
			 Fill_InterFMatrix(C_eper_dof, F_eper_dof, nz_val, row_idx, col_idx, K21);
			 Fill_InterFMatrix(C_eper_dof, C_eper_dof, nz_val, row_idx, col_idx, K22);
			 //text*******************
			//int n = 0;
			//std::cout << std::fixed << std::setprecision(6);
			//for (int i = 0; i < 12; i++)
			//{
			   // for (int j = 0; j < 12; j++)
			   // {
			   //	 std::cout <<n<<":"<< K12(i, j) << endl; // 输出矩阵元素值
			   //	 n++;
			   // }
			//}
		}

	}


	//Ce计算
	void NCF_map::Get_Ce(elastic_mat& data_mat, MatrixXd& Ce)
	{
		double em = data_mat.young_modulus;
		double nu = data_mat.poisson_ratio;
		double fac = (em * (1.0 - nu)) / ((1.0 + nu) * (1.0 - 2.0 * nu));
		double fac_a = fac * nu / (1.0 - nu);
		double fac_b = fac * (1.0 - 2.0 * nu) / (2.0 * (1.0 - nu));
		// 本构矩阵赋值
		Ce.resize(6, 6);
		Ce << fac, fac_a, fac_a, 0., 0., 0.,
			fac_a, fac, fac_a, 0., 0., 0.,
			fac_a, fac_a, fac, 0., 0., 0.,
			0., 0., 0., fac_b, 0., 0.,
			0., 0., 0., 0., fac_b, 0.,
			0., 0., 0., 0., 0., fac_b;

	}


	//计算B、Nm
	void NCF_map::Calculate_B_Nm(MatrixXd& p_gps, MatrixXd& nodes, MatrixXd& B,
		MatrixXd& Nm,int& n_node_mesh)
	{
		LagrangeBR r_out = lagrange_basis(p_gps, n_node_mesh);
		MatrixXd J, inv_J, dNdx, T_dNdx;
		J = nodes.transpose() * r_out.dNdxi_out;
		inv_J = J.inverse();
		dNdx = r_out.dNdxi_out * inv_J;
		T_dNdx = dNdx.transpose();
		B.setZero();
		Nm.setZero();
		for (int i = 0; i < n_node_mesh; i++)
		{
			B(0, 3 * i)     = T_dNdx(0, i);
			B(1, 3 * i + 1) = T_dNdx(1, i);
			B(2, 3 * i + 2) = T_dNdx(2, i);
			B(3, 3 * i)     = T_dNdx(1, i);   B(3, 3 * i + 1) = T_dNdx(0, i);
			B(4, 3 * i + 1) = T_dNdx(2, i);   B(4, 3 * i + 2) = T_dNdx(1, i);
			B(5, 3 * i)     = T_dNdx(2, i);   B(5, 3 * i + 2) = T_dNdx(0, i);

			Nm(0, 3 * i)     = r_out.N_out(0, i);
			Nm(1, 3 * i + 1) = r_out.N_out(0, i);
			Nm(2, 3 * i + 2) = r_out.N_out(0, i);

		}

	}



	//罚参数计算
	double NCF_map::Get_Alpha(data_management& data_cae, elastic_mat& data_mat, MatrixXd& pts1)
	{
		double x1, y1, z1, x2, y2, z2;
		x1 = pts1(0, 0); y1 = pts1(0, 1); z1 = pts1(0, 2);
		x2 = pts1(1, 0); y2 = pts1(1, 1); z2 = pts1(1, 2);
		double h = sqrt(pow(x1 - x2, 2) + pow(y1 - y2, 2) + pow(z1 - z2, 2));//单元尺寸大小
		double E = data_mat.young_modulus, nu = data_mat.poisson_ratio;
		double lambda = E * nu / (1 - 2 * nu) * (1 + nu);
		double mu = E / 2 * (1 + nu);
		double theta = 12;
		double alpha = theta * (lambda + mu) / (2 * h);
		return alpha;
	}


	//计算Kp，Kd
	void NCF_map::Calculate_Kp_Kd(MatrixXd& Kp11, MatrixXd& Kp12, MatrixXd& Kp22,
		MatrixXd& Kd11, MatrixXd& Kd12, MatrixXd& Kd21, MatrixXd& Kd22,
		MatrixXd& Nm1, MatrixXd& Nm2, MatrixXd& B1, MatrixXd& B2,
		MatrixXd& n, MatrixXd& Ce, double& alpha, double& wt1)
	{
		MatrixXd Tp_Nm11 = Nm1.transpose() * Nm1 * (alpha * wt1);
		MatrixXd Tp_Nm12 = Nm1.transpose() * Nm2 * (alpha * wt1);
		MatrixXd Tp_Nm22 = Nm2.transpose() * Nm2 * (alpha * wt1);

		Kp11 = Kp11+Tp_Nm11;
		Kp12 = Kp12+Tp_Nm12;
		Kp22 = Kp22+Tp_Nm22;

		MatrixXd Tp_Nm1_n = Nm1.transpose() * n * 0.5;
		MatrixXd Tp_Nm1_n_c = Tp_Nm1_n * Ce;
		MatrixXd Tp_Nm1_n_c_B1 = Tp_Nm1_n_c * B1 * wt1;
		MatrixXd Tp_Nm1_n_c_B2 = Tp_Nm1_n_c * B2 * wt1;

		MatrixXd Tp_Nm2_n = Nm2.transpose() * n * 0.5;
		MatrixXd Tp_Nm2_n_c = Tp_Nm2_n * Ce;
		MatrixXd Tp_Nm2_n_c_B1 = Tp_Nm2_n_c * B1 * wt1;
		MatrixXd Tp_Nm2_n_c_B2 = Tp_Nm2_n_c * B2 * wt1;

		Kd11 = Kd11 + Tp_Nm1_n_c_B1;
		Kd12 = Kd12 + Tp_Nm1_n_c_B2;
		Kd21 = Kd21 + Tp_Nm2_n_c_B1;
		Kd22 = Kd22 + Tp_Nm2_n_c_B2;
	}

	void NCF_map::GetIntF_ele_Inform(MatrixXd& pts1,MatrixXd& pts2,vector<int>& F_eper_dof,
		vector<int>& C_eper_dof,data_management& data_cae, int& e, int& n_node_F, int& n_node_C)
	{
		const auto& tps_F = data_cae.node_topos_[data_cae.BndMesh_F[e] - 1];//索引要减1
		const auto& tps_C = data_cae.node_topos_[data_cae.BndMesh_C[e] - 1];
		for (int i = 0; i < n_node_F; ++i)
		{
			const auto& coord_sctr1 = data_cae.coords_[tps_F[i] - 1];
			for (int j = 0; j < 3; ++j)
			{
				pts1(i, j) = coord_sctr1[j];
			}
		}
		for (int i = 0; i < n_node_C; ++i)
		{
			const auto& coord_sctr2 = data_cae.coords_[tps_C[i] - 1];
			for (int j = 0; j < 3; ++j)
			{
				pts2(i, j) = coord_sctr2[j];
			}
		}

		F_eper_dof.resize(n_node_F * 3);
		C_eper_dof.resize(n_node_C * 3);

		for (int i = 0; i < n_node_F; i++)
		{  
			// 细网格自由度
			int ncf_F = data_cae.BndMesh_F[e] - 1;//细网格单元编号从0开始
			int item_dof_F = data_cae.resort_free_nodes_[data_cae.node_topos_[ncf_F][i] - 1];
			F_eper_dof[3 * i] = 3 * item_dof_F;
			F_eper_dof[3 * i + 1] = 3 * item_dof_F + 1;
			F_eper_dof[3 * i + 2] = 3 * item_dof_F + 2;
			
		}
		for (int i = 0; i < n_node_C; i++)
		{  //粗网格自由度
			int ncf_C = data_cae.BndMesh_C[e] - 1;//粗网格单元编号从0开始
			int item_dof_C = data_cae.resort_free_nodes_[data_cae.node_topos_[ncf_C][i] - 1];
			C_eper_dof[3 * i] = 3 * item_dof_C;
			C_eper_dof[3 * i + 1] = 3 * item_dof_C + 1;
			C_eper_dof[3 * i + 2] = 3 * item_dof_C + 2;

		}
	}

	//计算总界面刚度矩阵
	void NCF_map::Calculate_InterFMatrix(MatrixXd& K11, MatrixXd& K12,MatrixXd& K21, MatrixXd& K22,
		MatrixXd& Kd11, MatrixXd& Kd12, MatrixXd& Kd21, MatrixXd& Kd22,
		MatrixXd& Kp11, MatrixXd& Kp12, MatrixXd& Kp22, int& n_node_F, int& n_node_C)
	{
		K11.resize(n_node_F * 3, n_node_F * 3);
		K12.resize(n_node_F * 3, n_node_C * 3);
		K21.resize(n_node_C * 3, n_node_F * 3);
		K22.resize(n_node_C * 3, n_node_C * 3);

		MatrixXd T_Kd11 = Kd11.transpose();
		MatrixXd T_Kd21 = Kd21.transpose();
		MatrixXd T_Kd12 = Kd12.transpose();
		MatrixXd T_Kd22 = Kd22.transpose();
		MatrixXd T_Kp12 = Kp12.transpose();
		
	
		for (int i = 0; i < n_node_F * 3; i++)
		{
			for (int j = 0; j < n_node_F * 3; j++)
			{
				K11(i, j) = -Kd11(i, j) - T_Kd11(i, j) + Kp11(i, j);
				/*K12(i, j) = -Kd12(i, j) + T_Kd21(i, j) - Kp12(i, j);
				K21(i, j) = Kd21(i, j) - T_Kd12(i, j) - T_Kp12(i, j);
				K22(i, j) = Kd22(i, j) + T_Kd22(i, j) + Kp22(i, j);*/
			}
		}

		for (int i = 0; i < n_node_F * 3; i++)
		{
			for (int j = 0; j < n_node_C * 3; j++)
			{
				K12(i, j) = -Kd12(i, j) + T_Kd21(i, j) - Kp12(i, j);
			}
		}

		for (int i = 0; i < n_node_C * 3; i++)
		{
			for (int j = 0; j < n_node_F * 3; j++)
			{
				K21(i, j) = Kd21(i, j) - T_Kd12(i, j) - T_Kp12(i, j);
			}
		}

		for (int i = 0; i < n_node_C* 3; i++)
		{
			for (int j = 0; j < n_node_C * 3; j++)
			{
				K22(i, j) = Kd22(i, j) + T_Kd22(i, j) + Kp22(i, j);
			}
		}


	}

	void NCF_map::Fill_InterFMatrix(vector<int>& j_eper_dof, vector<int>& i_eper_dof, 
		vector<double>& nz_val, vector<int>& row_idx, vector<int>& col_idx,
		MatrixXd& K_interface)
	{
		// 组装
		int ii_dof, jj_dof, loop_size = j_eper_dof.size();
		for (int mm = 0; mm < loop_size; mm++)
		{
			jj_dof = j_eper_dof[mm];
			if (jj_dof >= 0)
			{
				int start = col_idx[jj_dof];
				for (int nn = 0; nn < loop_size; nn++)
				{
					int t = start;
					ii_dof = i_eper_dof[nn];
					if (ii_dof >= 0)
					{
						for (; row_idx[t] < ii_dof; t++)
						{
						} // 使用上三角矩阵
						nz_val[t] = nz_val[t] + K_interface(mm, nn);
					}
				}
			}
		}


	}

	
	//积分点网格单元向父空间映射
	MatrixXd NCF_map::GlobalMap3D(MatrixXd gpoint, MatrixXd nodes, int& n_node_mesh)
	{
		const int nMax = 10;//d迭代次数改为15  10.11
		const double tol = 1e-14;
		double tolSquared = tol * tol;
		//double tolSquared = tol ;
		vector<double> xm(3);
		vector<double> columnSums(3);
		
		/*int rows = nodes.rows();
		int cols = nodes.cols();*/
		//移动坐标到中心
		for (int j = 0; j < 3; j++)
		{
			for (int i = 0; i < n_node_mesh; i++)
			{
				columnSums[j] += nodes(i, j);
			}
			xm[j] = columnSums[j] / n_node_mesh;
		}

		for (int j = 0; j < 3; j++)
		{
			for (int i = 0; i < n_node_mesh; i++)
			{
				nodes(i, j) -= xm[j];
			}
		}

		for (int j = 0; j < 3; j++)
		{
			gpoint(0, j) -= xm[j];
		}

		//Newton-Raphson iterations
		double dSi = 1.0;
		int n = 1;

		MatrixXd Xi(1, 3);
		Xi.setZero();
		
		while (dSi > tolSquared && n < nMax)
		{

			//调函数求形函数及偏导
			LagrangeBR result_out = lagrange_basis(Xi, n_node_mesh);
			//传出形函数及偏导
			MatrixXd N, dNdxi;
			dNdxi = result_out.dNdxi_out;
			N = result_out.N_out;

			MatrixXd x, hessian, inv_hessian, f, dxi;
			x = N * nodes;
			MatrixXd T_nodes = nodes.transpose();
			hessian = T_nodes * dNdxi;
			f = x - gpoint;
			inv_hessian = hessian.inverse();
			MatrixXd T_f = f.transpose();
			dxi = inv_hessian * T_f;
			MatrixXd T_dxi = dxi.transpose();

			Xi = Xi - T_dxi;
			n = n + 1;
			dSi = dxi(0, 0) * dxi(0, 0) + dxi(1, 0) * dxi(1, 0) + dxi(2, 0) * dxi(2, 0);

			if (n == nMax && dSi > tolSquared)//9.5
			{
				std::cout << "Warning: Mapping Gauss points; residual, dX = " << std::sqrt(dSi) << std::endl;
			}

		}
		
		return Xi;
	}


	NCF_map::LagrangeBR NCF_map::lagrange_basis(MatrixXd& coord, int& n_node_mesh)
	{
		double xi = coord(0, 0), eta = coord(0, 1), zeta = coord(0, 2);
		LagrangeBR result;
		if (n_node_mesh==8)
		{
			vector<double>  I1(3), I2(3);
			I1[0] = 0.5 - 0.5 * xi; I1[1] = 0.5 - 0.5 * eta; I1[2] = 0.5 - 0.5 * zeta;
			I2[0] = 0.5 + 0.5 * xi; I2[1] = 0.5 + 0.5 * eta; I2[2] = 0.5 + 0.5 * zeta;

			MatrixXd N(1, 8);//转置后的N
			N(0, 0) = I1[0] * I1[1] * I1[2];
			N(0, 1) = I2[0] * I1[1] * I1[2];
			N(0, 2) = I2[0] * I2[1] * I1[2];
			N(0, 3) = I1[0] * I2[1] * I1[2];
			N(0, 4) = I1[0] * I1[1] * I2[2];
			N(0, 5) = I2[0] * I1[1] * I2[2];
			N(0, 6) = I2[0] * I2[1] * I2[2];
			N(0, 7) = I1[0] * I2[1] * I2[2];

			MatrixXd dNdxi(8, 3);

			dNdxi(0, 0) = 0.125 * (-1 + eta + zeta - eta * zeta);
			dNdxi(1, 0) = 0.125 * (1 - eta - zeta + eta * zeta);
			dNdxi(2, 0) = 0.125 * (1 + eta - zeta - eta * zeta);
			dNdxi(3, 0) = 0.125 * (-1 - eta + zeta + eta * zeta);
			dNdxi(4, 0) = 0.125 * (-1 + eta - zeta + eta * zeta);
			dNdxi(5, 0) = 0.125 * (1 - eta + zeta - eta * zeta);
			dNdxi(6, 0) = 0.125 * (1 + eta + zeta + eta * zeta);
			dNdxi(7, 0) = 0.125 * (-1 - eta - zeta - eta * zeta);

			dNdxi(0, 1) = 0.125 * (-1 + xi + zeta - xi * zeta);
			dNdxi(1, 1) = 0.125 * (-1 - xi + zeta + xi * zeta);
			dNdxi(2, 1) = 0.125 * (1 + xi - zeta - xi * zeta);
			dNdxi(3, 1) = 0.125 * (1 - xi - zeta + xi * zeta);
			dNdxi(4, 1) = 0.125 * (-1 + xi - zeta + xi * zeta);
			dNdxi(5, 1) = 0.125 * (-1 - xi - zeta - xi * zeta);
			dNdxi(6, 1) = 0.125 * (1 + xi + zeta + xi * zeta);
			dNdxi(7, 1) = 0.125 * (1 - xi + zeta - xi * zeta);

			dNdxi(0, 2) = 0.125 * (-1 + xi + eta - xi * eta);
			dNdxi(1, 2) = 0.125 * (-1 - xi + eta + xi * eta);
			dNdxi(2, 2) = 0.125 * (-1 - xi - eta - xi * eta);
			dNdxi(3, 2) = 0.125 * (-1 + xi - eta + xi * eta);
			dNdxi(4, 2) = 0.125 * (1 - xi - eta + xi * eta);
			dNdxi(5, 2) = 0.125 * (1 + xi - eta - xi * eta);
			dNdxi(6, 2) = 0.125 * (1 + xi + eta + xi * eta);
			dNdxi(7, 2) = 0.125 * (1 - xi + eta - xi * eta);

			result.N_out = N;
			result.dNdxi_out = dNdxi;
		}
		else if(n_node_mesh==4)
		{
			MatrixXd N(1, 4);//转置后的N
			//N.resize(1, 8);
			N(0, 0) = 1 - xi - eta - zeta;
			N(0, 1) = xi;
			N(0, 2) = eta;
			N(0, 3) = zeta;

			MatrixXd dNdxi(4, 3);
			dNdxi(0, 0) = -1; dNdxi(0, 1) = -1; dNdxi(0, 2) = -1;
			dNdxi(1, 0) = 1;  dNdxi(1, 1) = 0;  dNdxi(1, 2) = 0;
			dNdxi(2, 0) = 0;  dNdxi(2, 1) = 1;  dNdxi(2, 2) = 0;
			dNdxi(3, 0) = 0;  dNdxi(3, 1) = 0;  dNdxi(3, 2) = 1;

			result.N_out = N;
			result.dNdxi_out = dNdxi;
		}
		
		return result;
	}

	;

}