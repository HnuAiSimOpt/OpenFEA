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
		// employing two - point Gaussian integral权重为1
		double gp_values = 1. / sqrt(3.);
		MatrixXd gps(4, 2);//高斯积分点
		gps << gp_values, gp_values,
			gp_values, -gp_values,
			-gp_values, gp_values,
			-gp_values, -gp_values;
		
		int n_interF = data_cae.BndMesh_F.size();
		int n_gps = gps.rows();
		GP1.resize(n_interF * n_gps, 7);
		GP2.resize(n_interF * n_gps, 3);//初始化。行数：积分点组数*交界处六面体单元个数
		int t = 0;         //记录交界处细网格面上高斯点总个数
		//交界面单元循环
		for (int e = 0; e < n_interF; e++) //交界面处细网格个数
		{
			MatrixXd nodes1(4, 3);//四边形交界面节点编号
			MatrixXd pts1(8, 3), pts2(8, 3);//交界处粗、细六面体单元节点坐标
			vector<int>  sctr1(8), sctr2(8);//交界处粗、细六面体单元节点编号
			int b_node1 = data_cae.bndFace_finemesh[e][0];
			int b_node2 = data_cae.bndFace_finemesh[e][1];
			int b_node3 = data_cae.bndFace_finemesh[e][2];
			int b_node4 = data_cae.bndFace_finemesh[e][3];

			nodes1 << data_cae.coords_[b_node1][0], data_cae.coords_[b_node1][1], data_cae.coords_[b_node1][2],
				data_cae.coords_[b_node2][0], data_cae.coords_[b_node2][1], data_cae.coords_[b_node2][2],
				data_cae.coords_[b_node3][0], data_cae.coords_[b_node3][1], data_cae.coords_[b_node3][2],
				data_cae.coords_[b_node4][0], data_cae.coords_[b_node4][1], data_cae.coords_[b_node4][2];
				
			GetIntF_ele_Inform(data_cae, sctr1, sctr2, pts1, pts2, e);
				
			//高斯积分点循环 4组  
			for (int q = 0; q < 4; q++)
			{
				
				MatrixXd dNdxi_1(2, 4);
				dNdxi_1(0, 0) = -0.25 * (1 - gps(q, 1));
				dNdxi_1(0, 1) =  0.25 * (1 - gps(q, 1));
				dNdxi_1(0, 2) =  0.25 * (1 + gps(q, 1));
				dNdxi_1(0, 3) = -0.25 * (1 + gps(q, 1));
				dNdxi_1(1, 0) = -0.25 * (1 - gps(q, 0));
				dNdxi_1(1, 1) = -0.25 * (1 + gps(q, 0));
				dNdxi_1(1, 2) =  0.25 * (1 + gps(q, 0));
				dNdxi_1(1, 3) =  0.25 * (1 - gps(q, 0));

				MatrixXd Jac = dNdxi_1* nodes1;
				Eigen::Vector3d a1 = Jac.row(0);
				Eigen::Vector3d a2 = Jac.row(1);
				Eigen::Vector3d a3 = a1.cross(a2);
				double norm_a3 = a3.norm();
				Eigen::Vector3d unit_a3 = a3.normalized();

				
				/*
				vector<double> a1(3), a2(3), a3(3), a3_normlized(3);
				a1[0] = Jac(0, 0); a1[1] = Jac(0, 1); a1[2] = Jac(0, 2);
				a2[0] = Jac(1, 0); a2[1] = Jac(1, 1); a2[2] = Jac(1, 2);
				*/
				//交界处单元高斯积分点物理坐标 
				MatrixXd N_1(1, 4);//形函数
				N_1(0, 0) = 0.25 * (1 - gps(q, 0)) * (1 - gps(q, 1));
				N_1(0, 1) = 0.25 * (1 + gps(q, 0)) * (1 - gps(q, 1));
				N_1(0, 2) = 0.25 * (1 + gps(q, 0)) * (1 + gps(q, 1));
				N_1(0, 3) = 0.25 * (1 - gps(q, 0)) * (1 + gps(q, 1));
				MatrixXd X = N_1 * nodes1;//高斯积分点物理坐标
				

				//积分点网格单元向父空间映射
				MatrixXd X1, X2;
				//-----------------------------------10.11测试
				std::cout << "F_mesh number :" << data_cae.BndMesh_F[e] << std::endl;
				X1 = GlobalMap3D(X, pts1);//细网格9.5
				std::cout << "C_mesh number :" << data_cae.BndMesh_C[e] << std::endl;
				X2 = GlobalMap3D(X, pts2);//粗网格9.5
				//-----------------------------------


				GP1(t, 0) = X1(0, 0);
				GP1(t, 1) = X1(0, 1);
				GP1(t, 2) = X1(0, 2);
				GP1(t, 3) = norm_a3;
				GP1(t, 4) = unit_a3[0];
				GP1(t, 5) = unit_a3[1];
				GP1(t, 4) = unit_a3[2];

				GP2(t, 0) = X2(0, 0);
				GP2(t, 1) = X2(0, 1);
				GP2(t, 2) = X2(0, 2);
				t++;

			}
		}

	}

	void NCF_map::GetIntF_ele_Inform(data_management& data_cae, vector<int>& sctr1, 
		vector<int>& sctr2,MatrixXd& pts1, MatrixXd& pts2, int& e)
	{
		const auto& node_topos_BndMesh_F = data_cae.node_topos_[data_cae.BndMesh_F[e]];
		const auto& node_topos_BndMesh_C = data_cae.node_topos_[data_cae.BndMesh_C[e]];

		const auto& coords = data_cae.coords_;

		for (int i = 0; i < 8; ++i)
		{
			sctr1[i] = node_topos_BndMesh_F[i];
			sctr2[i] = node_topos_BndMesh_C[i];

			const auto& coord_sctr1_i = coords[sctr1[i]];
			const auto& coord_sctr2_i = coords[sctr2[i]];

			for (int j = 0; j < 3; ++j)
			{
				pts1(i, j) = coord_sctr1_i[j];
				pts2(i, j) = coord_sctr2_i[j];
			}
		}
	}
	/*
	void NCF_map::GetIntF_ele_Inform(data_management& data_cae, vector<int>& sctr1, vector<int>& sctr2,
		MatrixXd& pts1, MatrixXd& pts2, int& e)
	{//交界处细网格
		sctr1[0] = data_cae.node_topos_[data_cae.BndMesh_F[e]][0]; 
		sctr1[1] = data_cae.node_topos_[data_cae.BndMesh_F[e]][1];
		sctr1[2] = data_cae.node_topos_[data_cae.BndMesh_F[e]][2];
		sctr1[3] = data_cae.node_topos_[data_cae.BndMesh_F[e]][3];
		sctr1[4] = data_cae.node_topos_[data_cae.BndMesh_F[e]][4];
		sctr1[5] = data_cae.node_topos_[data_cae.BndMesh_F[e]][5];
		sctr1[6] = data_cae.node_topos_[data_cae.BndMesh_F[e]][6];
		sctr1[7] = data_cae.node_topos_[data_cae.BndMesh_F[e]][7];
		sctr1[8] = data_cae.node_topos_[data_cae.BndMesh_F[e]][8];
		//交界处粗网格
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
		double em = data_mat.young_modulus;
		double nu = data_mat.poisson_ratio;
		double fac = (em * (1.0 - nu)) / ((1.0 + nu) * (1.0 - 2.0 * nu));
		double fac_a = fac * nu / (1.0 - nu);
		double fac_b = fac * (1.0 - 2.0 * nu) / (2.0 * (1.0 - nu));
		// 本构矩阵赋值
		MatrixXd Ce(6,6);
		 Ce << fac, fac_a, fac_a, 0., 0., 0.,
			fac_a, fac, fac_a, 0., 0., 0.,
			fac_a, fac_a, fac, 0., 0., 0.,
			0., 0., 0., fac_b, 0., 0.,
			0., 0., 0., 0., fac_b, 0.,
			0., 0., 0., 0., 0., fac_b;

		 int nF_bmesh = data_cae.BndMesh_F.size();
		 for (int e = 0; e < nF_bmesh; e++) //交界面处细网格个数
		 {
			 MatrixXd pts1(8, 3), pts2(8, 3);//交界处细、粗六面体单元节点坐标
			 vector<int>  sctr1(8), sctr2(8);//交界处细、粗六面体单元节点编号
			 GetIntF_ele_Inform(data_cae, sctr1, sctr2, pts1, pts2, e);
			
			 //初始化耦合矩阵
			 MatrixXd  Kp11(24, 24), Kp12(24, 24), Kp22(24, 24),
				 Kd11(24, 24), Kd12(24, 24), Kd21(24, 24), Kd22(24, 24);
			 Kp11.setZero();
			 Kp12.setZero();
			 Kp22.setZero();
			 Kd11.setZero();
			 Kd12.setZero();
			 Kd21.setZero();
			 Kd22.setZero();

			 //【***一个粗网格接受16个积分点，当粗细网格对应情况不同时，此处需要修改***】
			 for (int q = 4 * e ; q < 4 * e  + 4; q++)
			 {

				 MatrixXd pt1(1, 3), pt2(1, 3), n(3, 6);
				 pt1.row(0) = GP1.row(q).head(3);
				 pt2.row(0) = GP2.row(q).head(3);
				 double wt1 = GP1(q, 3);

				 vector<double> normal(3);
				 normal[0] = GP1(q, 4); normal[1] = GP1(q, 5); normal[2] = GP1(q, 6);
				 //法向量方向
				 n.setZero();
				 n(0, 0) = normal[0]; n(0, 3) = normal[1]; n(0, 5) = normal[2];
				 n(1, 1) = normal[1]; n(1, 3) = normal[0]; n(1, 4) = normal[2];
				 n(2, 2) = normal[2]; n(2, 4) = normal[1]; n(2, 5) = normal[0];

				 LagrangeBR r_out1 = lagrange_basis(pt1);
				 LagrangeBR r_out2 = lagrange_basis(pt2);
				 
				 MatrixXd J1, J2, inv_J1, inv_J2;
				 J1 = pts1.transpose() * r_out1.dNdxi_out;
				 J2 = pts2.transpose() * r_out2.dNdxi_out;
				 inv_J1 = J1.inverse();
				 inv_J2 = J2.inverse();

				 MatrixXd dN1dx, dN2dx, T_dN1dx, T_dN2dx;
				 dN1dx = r_out1.dNdxi_out * inv_J1;
				 dN2dx = r_out2.dNdxi_out * inv_J2;

				 T_dN1dx = dN1dx.transpose();
				 T_dN2dx = dN2dx.transpose();

				 MatrixXd B1(6, 24), B2(6, 24), Nm1(3, 24), Nm2(3, 24);
				 B1.setZero();
				 B2.setZero();
				 Nm1.setZero();
				 Nm2.setZero();

				 for (int i = 0; i < 8; i++) 
				 {
					 B1(0, 3 * i)     = T_dN1dx(0, i);
					 B1(1, 3 * i + 1) = T_dN1dx(1, i);
					 B1(2, 3 * i + 2) = T_dN1dx(2, i);
					 B1(3, 3 * i)     = T_dN1dx(1, i);   B1(3, 3 * i + 1) = T_dN1dx(0, i);
					 B1(4, 3 * i + 1) = T_dN1dx(2, i);   B1(4, 3 * i + 2) = T_dN1dx(1, i); 
					 B1(5, 3 * i)     = T_dN1dx(2, i);   B1(5, 3 * i + 2) = T_dN1dx(0, i);
					 
					 B2(0, 3 * i)     = T_dN2dx(0, i);
					 B2(1, 3 * i + 1) = T_dN2dx(1, i);
					 B2(2, 3 * i + 2) = T_dN2dx(2, i);
					 B2(3, 3 * i)     = T_dN2dx(1, i); B2(3, 3 * i + 1) = T_dN2dx(0, i);
					 B2(4, 3 * i + 1) = T_dN2dx(2, i); B2(4, 3 * i + 2) = T_dN2dx(1, i);
					 B2(5, 3 * i)     = T_dN2dx(2, i); B2(5, 3 * i + 2) = T_dN2dx(0, i);
					
					 Nm1(0, 3 * i)     = r_out1.N_out(0, i);
					 Nm1(1, 3 * i + 1) = r_out1.N_out(0, i);
					 Nm1(2, 3 * i + 2) = r_out1.N_out(0, i);

					 Nm2(0, 3 * i)     = r_out2.N_out(0, i);
					 Nm2(1, 3 * i + 1) = r_out2.N_out(0, i);
					 Nm2(2, 3 * i + 2) = r_out2.N_out(0, i);
				 }
				 
				 //计算罚参数
				 double alpha = Get_Alpha(data_cae, data_mat, pts1);
				 //计算Kp,Kd
				 Calculate_Kp_Kd(Kp11, Kp12, Kp22, Kd11, Kd12, Kd21, Kd22,
					 Nm1, Nm2, B1, B2, n, Ce, alpha, wt1);

			 }
			 //计算界面总刚
			 MatrixXd K11(24, 24), K12(24, 24), K21(24, 24), K22(24, 24);
			 vector<int> F_eper_dof(24), C_eper_dof(24);
			 Calculate_InterFMatrix(K11, K12, K21, K22, F_eper_dof, C_eper_dof,
				 Kd11, Kd12, Kd21, Kd22, Kp11, Kp12, Kp22, data_cae, e);

			 //界面刚度矩阵储存为CSR
			 Fill_InterFMatrix(F_eper_dof, F_eper_dof, nz_val, row_idx, col_idx, K11);
			 Fill_InterFMatrix(F_eper_dof, C_eper_dof, nz_val, row_idx, col_idx, K12);
			 Fill_InterFMatrix(C_eper_dof, F_eper_dof, nz_val, row_idx, col_idx, K21);
			 Fill_InterFMatrix(C_eper_dof, C_eper_dof, nz_val, row_idx, col_idx, K22);

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

		Kp11 += Tp_Nm11;
		Kp12 += Tp_Nm12;
		Kp22 += Tp_Nm22;

		MatrixXd Tp_Nm1_n = Nm1.transpose() * n * 0.5;
		MatrixXd Tp_Nm1_n_c = Tp_Nm1_n * Ce;
		MatrixXd Tp_Nm1_n_c_B1 = Tp_Nm1_n_c * B1 * wt1;
		MatrixXd Tp_Nm1_n_c_B2 = Tp_Nm1_n_c * B2 * wt1;

		MatrixXd Tp_Nm2_n = Nm2.transpose() * n * 0.5;
		MatrixXd Tp_Nm2_n_c = Tp_Nm2_n * Ce;
		MatrixXd Tp_Nm2_n_c_B1 = Tp_Nm2_n_c * B1 * wt1;
		MatrixXd Tp_Nm2_n_c_B2 = Tp_Nm2_n_c * B2 * wt1;

		Kd11 += Tp_Nm1_n_c_B1;
		Kd12 += Tp_Nm1_n_c_B2;
		Kd21 += Tp_Nm2_n_c_B1;
		Kd22 += Tp_Nm2_n_c_B2;
	}

	//计算总界面刚度矩阵
	void NCF_map::Calculate_InterFMatrix(MatrixXd& K11, MatrixXd& K12,
		MatrixXd& K21, MatrixXd& K22, vector<int> F_eper_dof, vector<int> C_eper_dof,
		MatrixXd& Kd11, MatrixXd& Kd12, MatrixXd& Kd21, MatrixXd& Kd22,
		MatrixXd& Kp11, MatrixXd& Kp12, MatrixXd& Kp22,
		data_management& data_cae, int& e)
	{
		MatrixXd T_Kd11 = Kd11.transpose();
		MatrixXd T_Kd21 = Kd21.transpose();
		MatrixXd T_Kd12 = Kd12.transpose();
		MatrixXd T_Kd22 = Kd22.transpose();
		MatrixXd T_Kp12 = Kp12.transpose();
		
		for (int i = 0; i < 8; i++)
		{
			// 自由度
			int item_dof_F = data_cae.resort_free_nodes_[data_cae.BndMesh_F[e] - 1];
			F_eper_dof[3 * i] = 3 * item_dof_F;
			F_eper_dof[3 * i + 1] = 3 * item_dof_F + 1;
			F_eper_dof[3 * i + 2] = 3 * item_dof_F + 2;

			int item_dof_C = data_cae.resort_free_nodes_[data_cae.BndMesh_C[e] - 1];
			C_eper_dof[3 * i] = 3 * item_dof_C;
			C_eper_dof[3 * i + 1] = 3 * item_dof_C + 1;
			C_eper_dof[3 * i + 2] = 3 * item_dof_C + 2;
		}
	
		for (int i = 0; i < 24; i++)
		{
			for (int j = 0; j < 24; j++)
			{
				K11(i, j) = -Kd11(i, j) - T_Kd11(i, j) + Kp11(i, j);
				K12(i, j) = -Kd12(i, j) + T_Kd21(i, j) - Kp12(i, j);
				K21(i, j) = Kd21(i, j) - T_Kd12(i, j) - T_Kp12(i, j);
				K22(i, j) = Kd22(i, j) + T_Kd22(i, j) + Kp22(i, j);
			}
		}
	}

	void NCF_map::Fill_InterFMatrix(vector<int> j_eper_dof, vector<int> i_eper_dof, 
		vector<double>& nz_val, vector<int> row_idx, vector<int>& col_idx, MatrixXd& K_interface)
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
	MatrixXd NCF_map::GlobalMap3D(MatrixXd gpoint, MatrixXd nodes)
	{
		const int nMax = 10;//d迭代次数改为15  10.11
		const double tol = 1e-14;
		double tolSquared = tol * tol;
		//double tolSquared = tol ;
		vector<double> xm(3);
		vector<double> columnSums(3);
		int rows = 8;
		int cols = 3;
		//移动坐标到中心
		for (int j = 0; j < cols; j++)
		{
			for (int i = 0; i < rows; i++)
			{
				columnSums[j] += nodes(i, j);
			}
			xm[j] = columnSums[j] / 8;//8节点
		}

		for (int j = 0; j < cols; j++)
		{
			for (int i = 0; i < rows; i++)
			{
				nodes(i, j) -= xm[j];
			}
		}

		for (int j = 0; j < cols; j++)
		{
			gpoint(0, j) -= xm[j];
		}

		//Newton-Raphson iterations
		double dSi = 1.0;
		int n = 1;

		MatrixXd Xi(1, 3);
		while (dSi > tolSquared && n < nMax)
		{

			//调函数求形函数及偏导
			LagrangeBR result_out = lagrange_basis(Xi);
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


	NCF_map::LagrangeBR NCF_map::lagrange_basis(MatrixXd& coord)
	{
		double xi = coord(0, 0), eta = coord(0, 1), zeta = coord(0, 2);
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

		LagrangeBR result;
		result.N_out = N;
		result.dNdxi_out = dNdxi;

		return result;
	}


















	;

}