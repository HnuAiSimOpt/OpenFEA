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
    class SFEM2D
    {
    public:
    data_management* femData;
    map<int, double> elementArea;       //三角形单元面积
    map<int, vector<double>> elementAngle;   //三角形单元角度（各个节点）

    map<int, Point*> nodeCoord;         //节点坐标
    map<int, vector<int>> nodeElement;  //节点所在的单元
    map<int, vector<int>> nodePosition; //节点所在的单元位置
    map<int, vector<double>> nodeAngle; //节点角度
    map<int, bool> nodeEdgeFlag;        //节点位于边缘的标志
    map<int, double> nodeArea;          //节点面积

    public:
        // 构造函数，析构函数
        SFEM2D(){};
        SFEM2D(data_management* femData) : femData(femData){};

        void build();
        private:
            void eleData();
	        void nodeData();
            double T3Area(const vector<double>& a, const vector<double>& b, const vector<double>& c);
            vector<double> T3Angle(const vector<double>& a, const vector<double>& b, const vector<double>& c);
    };
}
