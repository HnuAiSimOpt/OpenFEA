#include "include/SFEM3D.h"
#include "include/SupportForSFEM.h"
#include "include/SurfaceNode.h"
#include "Elements/include/SFEM_C3D4.h"
#include "include/data_management.h"
namespace CAE
{
    void SFEM3D::build()
    {
        if (!this->data_cae) {
            
        }
    	eleData();
        nodeData();



        //for (int i = 1; i <= data_cae->ne_; i++) {
        //    SFEM_C3D4* e = new SFEM_C3D4();
        //    e->n = i;
        //    e->sfemData = this;
        //    data_cae->ele_list_.push_back(e);
        //    data_cae->ele_list_idx_.push_back(i);
        //}
    }

    void SFEM3D::stiffness_Ns(vector<vector<double>>& ke, int Ns, elastic_mat data_mat)
    {
        // 材料属性赋值
        SFEM_C3D4 NS(this, data_mat, Ns);
        NS.build_ele_stiff_mat(ke);
    }

    void SFEM3D::eleData()
    {
        //检查是否所有单元都是四面体TODO
        for(int i =0; i < data_cae->ne_;i++){
             std::cout << "Only C3D4 ELement for 3D SFEM ......\n";
        }

        SurfaceNode surfaceJudge;
        //计算所有单元的体积
        for(int i = 0; i < data_cae->ne_;i++){
            int n1 = data_cae->node_topos_[i][0];
            int n2 = data_cae->node_topos_[i][1];
            int n3 = data_cae->node_topos_[i][2];
            int n4 = data_cae->node_topos_[i][3];

            const vector<double>& p1 = data_cae->coords_[n1];
            const vector<double>& p2 = data_cae->coords_[n2];
            const vector<double>& p3 = data_cae->coords_[n3];
            const vector<double>& p4 = data_cae->coords_[n4];

            elementVolume[i] = C3D4Volume(p1,p2,p3,p4);//四面体体积

            surfaceJudge.insertFace(n1, n2, n3);//单元面1
            surfaceJudge.insertFace(n1, n2, n4);//单元面2
            surfaceJudge.insertFace(n2, n3, n4);//单元面3
            surfaceJudge.insertFace(n1, n3, n4);//单元面4
        }

        for (int i = 0; i < data_cae->ne_; i++) {
            int n1 = data_cae->node_topos_[i][0];
            int n2 = data_cae->node_topos_[i][1];
            int n3 = data_cae->node_topos_[i][2];
            int n4 = data_cae->node_topos_[i][3];

            const vector<double>& p1 = data_cae->coords_[n1];
            const vector<double>& p2 = data_cae->coords_[n2];
            const vector<double>& p3 = data_cae->coords_[n3];
            const vector<double>& p4 = data_cae->coords_[n4];

            elementFaceFlag[i].push_back(surfaceJudge.isSurface(n1, n2, n3));//单元面1
            elementFaceFlag[i].push_back(surfaceJudge.isSurface(n1, n2, n4));//单元面2
            elementFaceFlag[i].push_back(surfaceJudge.isSurface(n2, n3, n4));//单元面3
            elementFaceFlag[i].push_back(surfaceJudge.isSurface(n1, n3, n4));//单元面4
        }



        this->surfaceNode = surfaceJudge.node();

        for (int i = 0; i < data_cae->ne_; i++) {
            for (auto node : data_cae->node_topos_[i]) {
                if (surfaceNode.find(node) != surfaceNode.end()) {
                    //elementEdgeFlag
                    break;
                }
            }
        }
        
    }

    void SFEM3D::nodeData()
    {
        //节点积分域数量
        this->ns = data_cae->nd_;

        //节点坐标
        for(int i = 0;i < data_cae->nd_;i++){
            double _x = data_cae->coords_[i][0];
            double _y = data_cae->coords_[i][1];
            double _z = data_cae->coords_[i][2];
            nodeCoord[i] = new Point(_x,_y,_z);
        }
        
        //节点所在的单元及在单元内部位置
        for(int i = 0; i < data_cae->ne_; i++){
            int _position = 0;
    		for (auto n : data_cae->node_topos_[i]) {
    			nodeElement[n].push_back(i);//记录节点所在单元
                nodePosition[n].push_back(_position);//记录节点在单元内部的位置(从0开始)
                _position++;
    		}
        }

         //节点周围的节点（节点积分域）
        for(auto node : nodeElement){
            set<int> _arNode;//临时使用set过滤重复节点
            for(int ie : node.second){//历遍此节点所在的所有单元
                for (auto n : data_cae->node_topos_[ie]){//此单元的节点
                    _arNode.insert(n);
                }
            }
            for(int n:_arNode){//将set中的数据存入vector
                nodeArNode[node.first].push_back(n);
            }
        }

        //节点体积
        for(auto node : nodeElement){
            for(int ie : node.second){//历遍此节点所在的所有单元
                nodeVolume[node.first] += elementVolume[ie] / 4;
            }
        }
    }
   
    //四面体面积计算
    double SFEM3D::C3D4Volume(const vector<double>& a, const vector<double>& b, const vector<double>& c, const vector<double>& d)
    {
        vector<double> v1(3);
        vector<double> v2(3);
        vector<double> v3(3);
        v1[0] = a[0] - d[0];
        v1[1] = a[1] - d[1];
        v1[2] = a[2] - d[2];
        v2[0] = b[0] - d[0];
        v2[1] = b[1] - d[1];
        v2[2] = b[2] - d[2];
        v3[0] = c[0] - d[0];
        v3[1] = c[1] - d[1];
        v3[2] = c[2] - d[2];

        double val = v1[0] * v2[1] * v3[2] + v1[1] * v2[2] * v3[0] + v1[2] * v2[0] * v3[1]
            - v1[0] * v2[2] * v3[1] - v1[1] * v2[0] * v3[2] - v1[2] * v2[1] * v3[0];
        return abs(val / 6.0);
    }
}