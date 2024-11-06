#include <iostream>
#include <string>
#include <map>
#include <set>
#include <vector>
#include <algorithm>

#pragma once

using namespace std;
namespace CAE
{
    //三维向量数据类
    class Point {
    public:
    	Point() {
            x = 0.0;
            y = 0.0;
            z = 0.0;
        };
    	Point(double x, double y, double z) { this->x = x; this->y = y; this->z = z; };
    	~Point() {

        };
    	double x, y, z;
    	double distance(Point& a) 
    	{
    		double dis = pow(x - a.x,2.0) + pow(y - a.y, 2.0) + pow(z - a.z, 2.0);
    		return pow(dis, 0.5);
    	}
    	Point operator +(const Point& a) { return Point(x + a.x, y + a.y, z + a.z); }
        Point operator -(const Point& a) { return Point(x - a.x, y - a.y, z - a.z); }

    	Point operator /(double a) { return Point(x / a, y / a, z / a); }
        Point operator ^(const Point& a) const {return Point(y * a.z - z * a.y, z * a.x - x * a.z, x * a.y - y * a.x);}
        Point operator *(double a) { return Point(x * a, y * a, z * a); }
        double operator *(const Point& a) { return (x * a.x + y * a.y + z * a.z); }

        double unit(){
            double d = sqrt(x * x + y * y + z * z);
            if(d != 0){
                x /= d;
                y /= d;
                z /= d;
            }
            return d;
        }
        double norm() { return sqrt(x * x + y * y + z * z); }
        Point normalized() {
            double d = sqrt(x * x + y * y + z * z);
            d = d == 0.0 ? 1.0 : 1.0 / d;
            return Point(x * d, y * d, z * d);

        }
        void negate() { x = -x; y = -y; z = -z; }


    };

    //3x3矩阵类
    class Mat3 {
    public:
    	Mat3() {};
    	~Mat3() {};
        double d[3][3]={0.};

        double& operator () (int i, int j)  {return d[i][j];}
        double* operator [] (int i){ return d[i]; }
        Point operator*(const Point& r) const
        {
            return Point(   d[0][0] * r.x + d[0][1] * r.y + d[0][2] * r.z,
                            d[1][0] * r.x + d[1][1] * r.y + d[1][2] * r.z,
                            d[2][0] * r.x + d[2][1] * r.y + d[2][2] * r.z);
        }


        Mat3& operator *= (double a) 
        {
            d[0][0] *= a;
            d[0][1] *= a;
            d[0][2] *= a;
            d[1][0] *= a;
            d[1][1] *= a;
            d[1][2] *= a;
            d[2][0] *= a;
            d[2][1] *= a;
            d[2][2] *= a;
            return(*this);
        }

        Mat3& inv() {
            double deti = 1 / (d[0][0] * (d[1][1] * d[2][2] - d[1][2] * d[2][1])
                + d[0][1] * (d[1][2] * d[2][0] - d[2][2] * d[1][0])
                + d[0][2] * (d[1][0] * d[2][1] - d[1][1] * d[2][0]));

            Mat3 ans;
            ans[0][0] = deti * (d[1][1] * d[2][2] - d[1][2] * d[2][1]);
            ans[1][0] = deti * (d[1][2] * d[2][0] - d[1][0] * d[2][2]);
            ans[2][0] = deti * (d[1][0] * d[2][1] - d[1][1] * d[2][0]);
            ans[0][1] = deti * (d[0][2] * d[2][1] - d[0][1] * d[2][2]);
            ans[1][1] = deti * (d[0][0] * d[2][2] - d[0][2] * d[2][0]);
            ans[2][1] = deti * (d[0][1] * d[2][0] - d[0][0] * d[2][1]);
            ans[0][2] = deti * (d[0][1] * d[1][2] - d[1][1] * d[0][2]);
            ans[1][2] = deti * (d[0][2] * d[1][0] - d[0][0] * d[1][2]);
            ans[2][2] = deti * (d[0][0] * d[1][1] - d[0][1] * d[1][0]);
            return ans;
        }

        bool trans(double ke[24][24]) 
        {
            double T[24][24] = { 0. };

            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    T[i][j]         = d[i][j];
                    T[i+3][j+3]     = d[i][j];
                    T[i+6][j+6]     = d[i][j];
                    T[i+9][j+9]     = d[i][j];
                    T[i+12][j+12]   = d[i][j];
                    T[i+15][j+15]   = d[i][j];
                    T[i+18][j+18]   = d[i][j];
                    T[i+21][j+21]   = d[i][j];
                }
            }


            double DUM1[24][24] = {};

            for (int i = 0; i < 24; ++i) {
                for (int j = 0; j < 24; ++j) {
                    for (int k = 0; k < 24; ++k) {
                        DUM1[i][j] += ke[i][k] * T[k][j];
                    }
                }
            }

            memset(ke, 0.0, sizeof(double) * 24 * 24);


            for (int i = 0; i < 24; ++i) {
                for (int j = 0; j < 24; ++j) {
                    for (int k = 0; k < 24; ++k) {
                        ke[i][j] += T[k][i] * DUM1[k][j];
                    }
                }
            }
            
            return true;
        }

        bool trans(double ke[18][18]) 
        {

            double T[18][18] = { 0. };

            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    T[i][j] = d[i][j];
                    T[i + 3][j + 3] = d[i][j];
                    T[i + 6][j + 6] = d[i][j];
                    T[i + 9][j + 9] = d[i][j];
                    T[i + 12][j + 12] = d[i][j];
                    T[i + 15][j + 15] = d[i][j];
                }
            }


            double DUM1[18][18] = {}; // 存放结果的二维数组

            for (int i = 0; i < 18; ++i) {
                for (int j = 0; j < 18; ++j) {
                    for (int k = 0; k < 18; ++k) {
                        DUM1[i][j] += ke[i][k] * T[k][j];
                    }
                }
            }

            memset(ke, 0.0, sizeof(double) * 18 * 18);


            for (int i = 0; i < 18; ++i) {
                for (int j = 0; j < 18; ++j) {
                    for (int k = 0; k < 18; ++k) {
                        ke[i][j] += T[k][i] * DUM1[k][j];
                    }
                }
            }

            return true;
        }

    };
}