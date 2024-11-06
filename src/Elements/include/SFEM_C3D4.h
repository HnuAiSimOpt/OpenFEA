/**************************************************************************

Copyright:  WH team

Author: GuoDaozhen <397908710@qq.com>

Completion date:  XXX

Description: XXX

**************************************************************************/

#pragma once

#include "ele_base.h"
namespace CAE
{
    class SFEM3D;
    class Point;
    class Mat3;
    class SFEM_C3D4 : public ele_base
    {
    public:
        elastic_mat matrial_struc_;
     

        SFEM3D* sfemData = nullptr;//光滑有限元数据

	 
	    int n = 0;//节点积分域编号

        static double H[4][4];
        static double Gr[4][4];
        static double Gs[4][4];


    public:
        // 构造函数，析构函数
        SFEM_C3D4() { type_ = "SFEM_C3D4"; nnode_ = 0; node_dof_ = 3; ngps_ = 0; };
        SFEM_C3D4(SFEM3D* sfemData, elastic_mat matrial_struc, int n ) : sfemData(sfemData), matrial_struc_(matrial_struc),n(n) { };
        
    
        void build_cons_mat(double D[6][6]);
        // 建立单元刚度矩阵
        void build_ele_stiff_mat(vector<vector<double>>& ke) ;

    private:
        void shapeGradient(Point* G, Mat3* GG);
        void T4G(Point* t4g, Point& node1, Point& node2, Point& node3, Point& node4);
        void shapeAtSurf(vector<Point>& lCoord, Point& xcc, double phi_c[4][4],
                           const vector<int>& nodeArray, int n[4], double vol,
                           const Point* t4g, Point* sg, Mat3* gg);
        int findFirstIndexOf(const vector<int>& nodeArray, int n);
        
    };
}