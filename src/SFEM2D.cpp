#include "include/SFEM2D.h"
#include "include/SupportForSFEM.h"
#include "include/data_management.h"
#include "Elements/include/SFEM_T3.h"
namespace CAE
{
    void SFEM2D::build()
    {
    	eleData();
        nodeData();
        for (int i = 0; i < femData->nd_; i++) {
            SFEM_T3* e = new SFEM_T3();
            e->n = i;
            e->sfemData = this;
            femData->ele_list_.push_back(e);
        }
    }

    void SFEM2D::eleData()
    {
        //检查是否所有单元都是三角形
        for(int i =0; i < femData->ne_;i++){
          /*  if(!data_cae->ele_list_->type_== "T3" ){
                std::cout << "Only T3 ELement for 2D SFEM ......\n";
            }*/
        }

        //计算所有单元的面积和角度
        for(int i = 0; i < femData->ne_;i++){
            int n1 = femData->node_topos_[i][0];
            int n2 = femData->node_topos_[i][1];
            int n3 = femData->node_topos_[i][2];

            const vector<double>& p1 = femData->coords_[n1];
            const vector<double>& p2 = femData->coords_[n2];
            const vector<double>& p3 = femData->coords_[n3];
              


            elementArea[i] = T3Area(p1,p2,p3);//三角形面积

            elementAngle[i] = T3Angle(p1,p2,p3);//三角形角度
        }
    }

    void SFEM2D::nodeData()
    {
        //节点坐标
        for(int i = 0;i < femData->nd_;i++){
            double x = femData->coords_[i][0];
            double y = femData->coords_[i][1];
            double z = femData->coords_[i][2];
            nodeCoord[i] = new Point(x,y,z);
        }
       
        //节点所在的单元及位置
        for(int i = 0; i < femData->ne_; i++){
            int position = 0;
    		for (auto n : femData->node_topos_[i]) {
    			nodeElement[n].push_back(i);//记录节点所在单元
                nodePosition[n].push_back(position);//记录节点在单元内部的位置(从0开始)
              
                position++;
    		}
        }

        //节点面积
        for(auto node : nodeElement){
            for(int i : node.second){//历遍此节点所在的所有单元
            nodeArea[node.first] += elementArea[i] / 3;
            }
        }

        //节点的角度
        for(auto node : nodeElement){
            vector<int> & ele = node.second;//节点所在的所有单元
            vector<int> & pos = nodePosition[node.first];//节点在单元内位置

            for(int i = 0; i < ele.size();i++){//历遍此节点所在的所有单元
                int e = ele[i];//单元编号
                int p = pos[i];//节点位置
                nodeAngle[node.first].push_back(elementAngle[e][p]);
            }
        }

        //节点位于边缘的判断
        for(auto node : nodeAngle){
            double totalAngle = 0.0;
            for(double angle : node.second){
                totalAngle += angle;
            }

            if(abs(totalAngle-360) < 0.00000001){
                nodeEdgeFlag[node.first] = false;
            }
            else{
                nodeEdgeFlag[node.first] = true;
            }
        }
    }

    //三角形单元面积计算
    double SFEM2D::T3Area(const vector<double>& a, const vector<double>& b, const vector<double>& c)
    {
        return 0.0;
    }

    //三角形单元角度计算
    vector<double> SFEM2D::T3Angle(const vector<double>& a, const vector<double>& b, const vector<double>& c)
    {
        vector<double> angle(3);

        return angle;
    }

}