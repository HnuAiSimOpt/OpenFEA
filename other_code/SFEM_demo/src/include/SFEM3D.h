/**************************************************************************

Copyright:  WH team

Author: GuoDaozhen <397908710@qq.com>

Completion date:  XXX

Description: XXX

**************************************************************************/

#pragma once

#include <iostream>
#include <string>
#include <map>
#include <set>
#include <vector>
#include <algorithm>

using namespace std;
namespace CAE
{
    class data_management;
    class Point;
    class Mat3;
    class elastic_mat;
    class SFEM3D
    {
    public:
    data_management* data_cae;

    
    map<int, double> elementVolume;                 //四面体单元体积
    map<int, vector<bool>> elementFaceFlag;         //单元位于表面的标志

    int ns = 0;
    map<int, Point*> nodeCoord;                     //节点坐标
    map<int, vector<int>> nodeElement;              //节点所在的单元
    map<int, vector<int>> nodePosition;             //节点所在的单元位置
    map<int, vector<int>> nodeArNode;               //节点周围的节点
    map<int, double> nodeVolume;                    //节点体积
    set<int> surfaceNode;                           //节点位于表面的标志


    public:
        // 构造函数，析构函数
        SFEM3D(){};
        SFEM3D(data_management* femData) : data_cae(femData){};

        void build();

        void stiffness_Ns(vector<vector<double>>& ke, int Ns, elastic_mat data_mat);

        private:
            void eleData();
	        void nodeData();
            double C3D4Volume(const vector<double>& a, const vector<double>& b, const vector<double>& c, const vector<double>& d);
    };
}
