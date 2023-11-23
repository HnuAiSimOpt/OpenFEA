#include "include/data2cae.h"

namespace CAE
{
    // 读取单元、节点总数
    void ReadInfo::read_ele_node_num(data_management &data_cae)
    {
        int nd, ne;
        std::ifstream infile(path_.c_str(), std::ios::in);
        if (!infile)
        {
            std::cerr << "Error: Cannot open " << path_ << std::endl;
            exit(EXIT_FAILURE);
        }
        string line;
        while (getline(infile, line))
        {
            if (line.find("*Nset") != string::npos)
            {
                getline(infile, line);
                std::istringstream iss(line);
                string t1, nd_temp, t2;
                iss >> t1 >> nd_temp >> t2;
                nd = atoi(nd_temp.c_str());
            }
            if (line.find("*Elset") != string::npos)
            {
                getline(infile, line);
                std::istringstream iss(line);
                string t1, ne_temp, t2;
                iss >> t1 >> ne_temp >> t2;
                ne = atoi(ne_temp.c_str());
                break;
            }
        }
        infile.close();
        data_cae.ne_ = ne;
        data_cae.nd_ = nd;
        cout << "the number of element is: " << data_cae.ne_ << endl
             << "the number of node is: " << data_cae.nd_ << endl;
    }

    // 基于指定符号划分字符串
    vector<string> ReadInfo::split_str(const string str, const string part)
    {
        vector<string> resVec;
        if (str == "")
        {
            return resVec;
        }
        string strs = str + part;
        size_t pos = strs.find(part);
        size_t size = strs.size();
        while (pos != string::npos)
        {
            string x = strs.substr(0, pos);
            resVec.push_back(x);
            strs = strs.substr(pos + 1, size);
            pos = strs.find(part);
        }
        return resVec;
    }

    // 读取节点坐标，节点拓扑关系，单元类型
    void ReadInfo::read_geo_mesh(data_management &data_cae)
    {
        // 初始化单元数据
        data_cae.coords_.resize(data_cae.nd_, vector<double>(3));
        data_cae.node_topos_.resize(data_cae.ne_, vector<int>(8));
        data_cae.ele_type_.resize(data_cae.ne_);
        // 读取计算文件
        std::ifstream infile(path_.c_str(), std::ios::in);
        string line;
        int id_node = 0, id_ele = 0;
        while (getline(infile, line))
        {
            // 读取节点坐标

            if (line.find("*Node") != string::npos)
            {
                while (getline(infile, line))
                {
                    if (line.find("*") != string::npos)
                        break;
                    else
                    {
                        string id, x, y, z;
                        double x_, y_, z_;
                        std::istringstream iss(line);
                        iss >> id >> x >> y >> z;
                        x.erase(x.end() - 1); // 删除字符串最后的符号
                        x_ = stod(x);         // 转换字符串为double
                        y.erase(y.end() - 1);
                        y_ = stod(y);
                        z_ = stod(z);
                        data_cae.coords_[id_node][0] = x_;
                        data_cae.coords_[id_node][1] = y_;
                        data_cae.coords_[id_node][2] = z_;
                        id_node = id_node + 1;
                    }
                }
            }

            // 读取单元类型及节点拓扑关系

            std::vector<string> type_temp;
            if (line.find("*Element") != string::npos)
            {
                type_temp = split_str(line, "=");
                while (getline(infile, line))
                {
                    bool flag = true;
                    if (line.find("*Element") != string::npos)
                    {
                        type_temp = split_str(line, "=");
                        flag = false;
                        getline(infile, line);
                    }
                    if (line.find("*") != string::npos && flag)
                    {
                        break;
                    }
                    if (type_temp[1] == "C3D4")
                    {
                        string id, node1, node2, node3, node4;
                        int id_, node1_, node2_, node3_, node4_;
                        std::istringstream iss(line);
                        iss >> id >> node1 >> node2 >> node3 >> node4;
                        id.erase(id.end() - 1);               // 删除字符串最后的符号
                        id_ = atoi(id.c_str());               // 转换字符串为int
                        data_cae.ele_type_[id_ - 1] = "C3D4"; // 记录单元类型
                        node1.erase(node1.end() - 1);
                        node1_ = atoi(node1.c_str());
                        node2.erase(node2.end() - 1);
                        node2_ = atoi(node2.c_str());
                        node3.erase(node3.end() - 1);
                        node3_ = atoi(node3.c_str());
                        node4_ = atoi(node4.c_str());
                        data_cae.node_topos_[id_ - 1][0] = node1_;
                        data_cae.node_topos_[id_ - 1][1] = node2_;
                        data_cae.node_topos_[id_ - 1][2] = node3_;
                        data_cae.node_topos_[id_ - 1][3] = node4_;
                    }
                    else if (type_temp[1] == "C3D8R")
                    {
                        string id, node1, node2, node3, node4, node5, node6, node7, node8;
                        int id_, node1_, node2_, node3_, node4_, node5_, node6_, node7_, node8_;
                        std::istringstream iss(line);
                        iss >> id >> node1 >> node2 >> node3 >> node4 >> node5 >> node6 >> node7 >> node8;
                        id.erase(id.end() - 1);
                        id_ = atoi(id.c_str());
                        data_cae.ele_type_[id_ - 1] = "C3D8R";
                        node1.erase(node1.end() - 1);
                        node1_ = atoi(node1.c_str());
                        node2.erase(node2.end() - 1);
                        node2_ = atoi(node2.c_str());
                        node3.erase(node3.end() - 1);
                        node3_ = atoi(node3.c_str());
                        node4.erase(node4.end() - 1);
                        node4_ = atoi(node4.c_str());
                        node5.erase(node5.end() - 1);
                        node5_ = atoi(node5.c_str());
                        node6.erase(node6.end() - 1);
                        node6_ = atoi(node6.c_str());
                        node7.erase(node7.end() - 1);
                        node7_ = atoi(node7.c_str());
                        node8_ = atoi(node8.c_str());
                        data_cae.node_topos_[id_ - 1][0] = node1_;
                        data_cae.node_topos_[id_ - 1][1] = node2_;
                        data_cae.node_topos_[id_ - 1][2] = node3_;
                        data_cae.node_topos_[id_ - 1][3] = node4_;
                        data_cae.node_topos_[id_ - 1][4] = node5_;
                        data_cae.node_topos_[id_ - 1][5] = node6_;
                        data_cae.node_topos_[id_ - 1][6] = node7_;
                        data_cae.node_topos_[id_ - 1][7] = node8_;
                        id_ele = id_;
                    }
                }
            }
        }
        cout << "a total of " << id_node << " nodal coordinates have been readed." << endl;
        cout << "a total of " << id_ele << " nodal connectivities have been readed." << endl;
    }
}
