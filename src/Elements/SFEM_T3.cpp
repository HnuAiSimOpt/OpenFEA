/**************************************************************************

Copyright:  WH team

Author: ChenBinXiang <chen_mech99@163.com>

Completion date:  XXX

Description: XXX

**************************************************************************/


#include "include/SFEM_T3.h"
#include "../include/SupportForSFEM.h"
#include "../include/SFEM2D.h"
#include "../include/data_management.h"
namespace CAE
{
    //REGISTER(ele_base, SFEM_T3, "SFEM_T3");

    // 建立单元刚度矩阵
    void SFEM_T3::build_ele_stiff_mat()
    {
		double area = sfemData->nodeArea[n];

        int ind = 0;//节点域的相关节点数量
	    vector<int> nv;//节点域的相关节点
	    double NS_cartd[2][20] = { 0. };
	    vector<double> NS_shape;

		vector <vector<double>> Bmax(3);
	    for (int i = 0; i < sfemData->nodeElement[n].size(); i++) {//遍历节点所在所有单元
	    	int element = sfemData->nodeElement[n][i];
	    	vector<int> lnv = sfemData->femData->node_topos_[element];
    
    
	    	int n1 = lnv[0], n2 = lnv[1], n3 = lnv[2];
    
	    	Point node1 = *sfemData->nodeCoord[n1];
	    	Point node2 = *sfemData->nodeCoord[n2];
	    	Point node3 = *sfemData->nodeCoord[n3];
    
	    	Point edge1 = (node1 + node2) / 2.;
	    	Point edge2 = (node2 + node3) / 2.;
	    	Point edge3 = (node1 + node3) / 2.;
	    	Point center = (node1 + node2 + node3) / 3.;//TODO 计算边中点和单元中心坐标
    
	    	//基于节点在单元的位置，构建子积分域的数据
	    	vector<Point> xgg(3);
	    	double phi_c[4][3] = { 0. };
	    	int k1, k2, k3;
	    	if (n == n1) {
	    		k1 = n1;			k2 = n2;				k3 = n3;
	    		xgg[0] = edge1;		xgg[1] = center;		xgg[2] = edge3;
    
	    		phi_c[0][0] = 0.5;		phi_c[0][1] = 0.5;
	    		phi_c[1][0] = 1.0 / 3;	phi_c[1][1] = 1.0 / 3;	phi_c[1][2] = 1.0 / 3;
	    		phi_c[2][0] = 0.5;		phi_c[2][2] = 0.5;
	    		phi_c[3][0] = 1.0;
	    	}
	    	else if (n == n2) {
	    		k1 = n2;			k2 = n3;				k3 = n1;
	    		xgg[0] = edge2;		xgg[1] = center;		xgg[2] = edge1;
    
	    		phi_c[0][1] = 0.5;		phi_c[0][2] = 0.5;
	    		phi_c[1][0] = 1.0 / 3;	phi_c[1][1] = 1.0 / 3;	phi_c[1][2] = 1.0 / 3;
	    		phi_c[2][0] = 0.5;		phi_c[2][1] = 0.5;
	    		phi_c[3][1] = 1.0;
	    	}
	    	else if (n == n3) {
	    		k1 = n3;			k2 = n1;				k3 = n2;
	    		xgg[0] = edge3;		xgg[1] = center;		xgg[2] = edge2;
    
	    		phi_c[0][0] = 0.5;		phi_c[0][2] = 0.5;
	    		phi_c[1][0] = 1.0 / 3;	phi_c[1][1] = 1.0 / 3;	phi_c[1][2] = 1.0 / 3;
	    		phi_c[2][1] = 0.5;		phi_c[2][2] = 0.5;
	    		phi_c[3][2] = 1.0;
	    	}
    
	    	//子积分域的两部分
	    	for (int i2 = 0; i2 < 2; i2++) {
	    		Point p1 = xgg[i2], p2 = xgg[i2 + 1];
    
	    		double xL = p1.distance(p2);
	    		double nx = (p2.y - p1.y) / xL;
	    		double ny = (p1.x - p2.x) / xL;
	    		double xg = (p1.x + p2.x) / 2;
	    		double yg = (p1.y + p2.y) / 2;
    
	    		double phi[3] = { 0. };
	    		phi[0] = (phi_c[i2][0] + phi_c[i2 + 1][0]) * 0.5;
	    		phi[1] = (phi_c[i2][1] + phi_c[i2 + 1][1]) * 0.5;
	    		phi[2] = (phi_c[i2][2] + phi_c[i2 + 1][2]) * 0.5;
    
	    		if (i + i2 == 0) {//如果是第一个单元的第一个子域（第一块）
	    			nv.push_back(lnv[0]);
	    			nv.push_back(lnv[1]);
	    			nv.push_back(lnv[2]);
	    			ind = 3;//节点域涉及的节点数量
	    			for (int i3 = 0; i3 < 3; i3++) {
	    				NS_cartd[0][i3] = phi[i3] * nx * xL / area;
	    				NS_cartd[1][i3] = phi[i3] * ny * xL / area;
	    			}
	    		}
	    		else {
	    			for (int i3 = 0; i3 < 3; i3++) {
	    				int flag = 0;
	    				for (int j = 0; j < ind; j++) {
	    					if (lnv[i3] == nv[j]) {
	    						NS_cartd[0][j] += phi[i3] * nx * xL / area;
	    						NS_cartd[1][j] += phi[i3] * ny * xL / area;
	    						flag = 1;
	    						break;
	    					}
	    				}
	    				if (flag == 0) {
	    					ind++;
	    					nv.push_back(lnv[i3]);
	    					NS_cartd[0][ind - 1] = phi[i3] * nx * xL / area;
	    					NS_cartd[1][ind - 1] = phi[i3] * ny * xL / area;
	    				}
	    			}
	    		}
	    	}
    
	    	int kk1 = sfemData->nodeEdgeFlag[k1] ? 1 : 0;
	    	int kk2 = sfemData->nodeEdgeFlag[k2] ? 1 : 0;
	    	int kk3 = sfemData->nodeEdgeFlag[k3] ? 1 : 0;
	    	int kk4 = kk1 + kk2 + kk3;
    
	    	if (kk1 == 1 & kk4 >= 2) {//如果第一个节点是边上且还有其他节点也在边上
	    		if (kk2 == 1) {//如果节点2也在边上
	    			int kkp = k2;
	    			Point xp1 = *sfemData->nodeCoord[k1];
	    			Point xp2 = (*sfemData->nodeCoord[k1] + *sfemData->nodeCoord[k2]) / 2;
	    			double xL = xp1.distance(xp2);
	    			double nx = (xp2.y - xp1.y) / xL;
	    			double ny = (xp1.x - xp2.x) / xL;
	    			double xg = (xp1.x + xp2.x) / 2;
	    			double yg = (xp1.y + xp2.y) / 2;
    
	    			double phi[3] = { 0. };
	    			phi[0] = (phi_c[3][0] + phi_c[0][0]) * 0.5;
	    			phi[1] = (phi_c[3][1] + phi_c[0][1]) * 0.5;
	    			phi[2] = (phi_c[3][2] + phi_c[0][2]) * 0.5;
	    			for (int i3 = 0; i3 < 3; i3++) {
	    				int flag = 0;
	    				for (int j = 0; j < ind; j++) {
	    					if (lnv[i3] == nv[j]) {
	    						NS_cartd[0][j] += phi[i3] * nx * xL / area;
	    						NS_cartd[1][j] += phi[i3] * ny * xL / area;
	    						flag = 1;
	    						break;
	    					}
	    				}
	    				if (flag == 0) {
	    					ind++;
	    					nv.push_back(lnv[i3]);
	    					NS_cartd[0][ind - 1] = phi[i3] * nx * xL / area;
	    					NS_cartd[1][ind - 1] = phi[i3] * ny * xL / area;
	    				}
	    			}
	    		}
    
	    		if (kk3 == 1) {//如果节点2也在边上
	    			int kkp = k3;
	    			Point xp1 = (*sfemData->nodeCoord[k1] + *sfemData->nodeCoord[kkp]) / 2;
	    			Point xp2 = *sfemData->nodeCoord[k1];
	    			double xL = xp1.distance(xp2);
	    			double nx = (xp2.y - xp1.y) / xL;
	    			double ny = (xp1.x - xp2.x) / xL;
	    			double xg = (xp1.x + xp2.x) / 2;
	    			double yg = (xp1.y + xp2.y) / 2;
    
	    			double phi[3] = { 0. };
	    			phi[0] = (phi_c[3][0] + phi_c[2][0]) * 0.5;
	    			phi[1] = (phi_c[3][1] + phi_c[2][1]) * 0.5;
	    			phi[2] = (phi_c[3][2] + phi_c[2][2]) * 0.5;
	    			for (int i3 = 0; i3 < 3; i3++) {
	    				int flag = 0;
	    				for (int j = 0; j < ind; j++) {
	    					if (lnv[i3] == nv[j]) {
	    						NS_cartd[0][j] += phi[i3] * nx * xL / area;
	    						NS_cartd[1][j] += phi[i3] * ny * xL / area;
	    						flag = 1;
	    						break;
	    					}
	    				}
	    				if (flag == 0) {
	    					ind++;
	    					nv.push_back(lnv[i3]);
	    					NS_cartd[0][ind - 1] = phi[i3] * nx * xL / area;
	    					NS_cartd[1][ind - 1] = phi[i3] * ny * xL / area;
	    				}
	    			}
	    		}
	    	}
	    	for (int i = 0; i < ind; i++) {
	    		if (nv[i] == n) {
	    			NS_shape.push_back(1.0);
	    		}
	    		else {
	    			NS_shape.push_back(0.0);
	    		}
	    	}
	    }
    
	    for (int i = 0; i < ind; i++) {
	    	Bmax[0].push_back(NS_cartd[0][i]);
	    	Bmax[1].push_back(0.0);
	    	Bmax[2].push_back(NS_cartd[1][i]);
    
	    	Bmax[0].push_back(0.0);
	    	Bmax[1].push_back(NS_cartd[1][i]);
	    	Bmax[2].push_back(NS_cartd[0][i]);
	    }
    }

}