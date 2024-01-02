/**************************************************************************

Copyright:  WH team

Author: YinJichao <jichaoyinyjc@163.com>

Completion date:  2023/11/23

Description: read a calculation file

**************************************************************************/

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
        data_cae.ele_list_idx_.resize(data_cae.ne_);
        // 读取计算文件
        std::ifstream infile(path_.c_str(), std::ios::in);
        string line;
        int id_node = 0, id_ele = 0, ele_type_idx = 0;
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
                int ele_nnode;
                del_blank(type_temp[1]);
                ele_nnode = data_cae.add_ele(type_temp[1]);
                while (getline(infile, line))
                {
                    bool flag = true;
                    if (line.find("*Element") != string::npos)
                    {
                        type_temp = split_str(line, "=");
                        del_blank(type_temp[1]);
                        ele_nnode = data_cae.add_ele(type_temp[1]);
                        ele_type_idx++;
                        flag = false;
                        getline(infile, line);
                    }
                    if (line.find("*") != string::npos && flag)
                    {
                        break;
                    }
                    std::istringstream iss(line);
                    vector<string> temp_a;
                    string temp;
                    while (getline(iss, temp, ',')) 
                    {
                        del_blank(temp);
                        temp_a.push_back(temp);
                    }
                    int id_ = atoi(temp_a[0].c_str());
                    data_cae.ele_list_idx_[id_ - 1] = ele_type_idx;
                    for (int i = 0; i < ele_nnode; i++)
                    {
                        data_cae.node_topos_[id_ - 1][i] = atoi(temp_a[i + 1].c_str());
                    }
                    id_ele++;
                }
            }
        }
        cout << "a total of " << id_node << " nodal coordinates have been readed." << endl;
        cout << "a total of " << id_ele << " nodal connectivities have been readed." << endl;
        cout << "a total of " << ele_type_idx + 1 << " types of elements have been readed." << endl;
    }

    //读取非协调信息
    void ReadInfo::readNconformingMessage(data_management& data_cae)
    {
        std::ifstream infile(path_.c_str(), std::ios::in);
        string line;
        while (getline(infile, line))
        {
            string node_id;
            int node_id_;
            //读取细网格单元8节点信息
            if (line.find("*BndMesh_F") != string::npos)
            {
                while (getline(infile, line))
                {
                    if (line.find("*") != string::npos)
                        break;
                    else
                    {
                        std::istringstream iss(line);
                        while (iss >> node_id)
                        {
                            if (node_id.find(",") != string::npos)
                            {
                                node_id.erase(node_id.end() - 1);   // delete the last symbol in the string
                            }
                            node_id_ = atoi(node_id.c_str());       // convert to int type
                            data_cae.BndMesh_F.push_back(node_id_);
                        }
                    }
                }
            }
            //读取粗网格单元8节点信息
            if (line.find("*BndMesh_C") != string::npos)
            {
                while (getline(infile, line))
                {
                    if (line.find("*") != string::npos)
                        break;
                    else
                    {
                        std::istringstream iss(line);
                        while (iss >> node_id)
                        {
                            if (node_id.find(",") != string::npos)
                            {
                                node_id.erase(node_id.end() - 1);   // delete the last symbol in the string
                            }
                            node_id_ = atoi(node_id.c_str());       // convert to int type
                            data_cae.BndMesh_C.push_back(node_id_);
                        }
                    }
                }
            }
            //读取细网格 交界面4节点信息
            if (line.find("*BndFace_finemesh") != string::npos)
            {  
                int i;
                i = data_cae.BndMesh_F.size();
                data_cae.bndFace_finemesh.resize(i, vector<int>(4));
                int id_node = 0;
                while (getline(infile, line))
                {
                    if (line.find("*") != string::npos)
                        break;
                    else
                    {
                        string  x1, x2, x3, x4;
                        double x1_, x2_, x3_, x4_;
                        std::istringstream iss(line);
                        iss >> x1 >> x2 >> x3 >> x4;
                        x1.erase(x1.end() - 1);  
                        x1_ = stod(x1);
                        x2.erase(x2.end() - 1);
                        x2_ = stod(x2);
                        x3.erase(x3.end() - 1);
                        x3_ = stod(x3);
                        x4_ = stod(x4);
                        data_cae.bndFace_finemesh[id_node][0] = x1_;
                        data_cae.bndFace_finemesh[id_node][1] = x2_;
                        data_cae.bndFace_finemesh[id_node][2] = x3_;
                        data_cae.bndFace_finemesh[id_node][3] = x4_;
                        id_node = id_node + 1;
                    }
                }
            }
            //读取粗网格 交界面4节点信息
            if (line.find("*BndFace_coarsemesh") != string::npos)
            {    
                int j;
                j = data_cae.BndMesh_C.size();
                data_cae.bndFace_coarsemesh.resize(j, vector<int>(4));
                int id_node = 0;
                while (getline(infile, line))
                {
                    if (line.find("*") != string::npos)
                        break;
                    else
                    {
                        string  x1, x2, x3, x4;
                        double x1_, x2_, x3_, x4_;
                        std::istringstream iss(line);
                        iss >> x1 >> x2 >> x3 >> x4;
                        x1.erase(x1.end() - 1);  
                        x1_ = stod(x1);
                        x2.erase(x2.end() - 1);
                        x2_ = stod(x2);
                        x3.erase(x3.end() - 1);
                        x3_ = stod(x3);
                        x4_ = stod(x4);
                        data_cae.bndFace_coarsemesh[id_node][0] = x1_;
                        data_cae.bndFace_coarsemesh[id_node][1] = x2_;
                        data_cae.bndFace_coarsemesh[id_node][2] = x3_;
                        data_cae.bndFace_coarsemesh[id_node][3] = x4_;
                        id_node = id_node + 1;
                    }
                }
            }
            if (line.find("*End Nconforming") != string::npos)
            {
                break;
            }

        }
        infile.close();

    }

    // 读取载荷信息
    void ReadInfo::read_load_bcs(string load_set_keyword, string load_value_keyword, data_management &data_cae)
    {
        // 读取计算文件
        std::ifstream infile(path_.c_str(), std::ios::in);
        string line;
        // 读取载荷节点集合
        while (getline(infile, line))
        {
            string node_id;
            int node_id_;
            if (line.find(load_set_keyword) != string::npos)
            {
                while (getline(infile, line))
                {
                    if (line.find("*") != string::npos)
                        break;
                    else
                    {
                        std::istringstream iss(line);
                        while (iss >> node_id)
                        {
                            if (node_id.find(",") != string::npos)
                            {
                                node_id.erase(node_id.end() - 1); // 删除字符串最后的符号
                            }
                            node_id_ = atoi(node_id.c_str()); // 转换字符串为int
                            data_cae.load_set_.push_back(node_id_);
                        }
                    }
                }
            }
            if (line.find("*End Assembly") != string::npos)
            {
                break;
            }
        }
        // 读取载荷自由度和幅值
        while (getline(infile, line))
        {
            string load_name, dof_id, value_str;
            int dof_id_;
            double value_;
            if (line.find(load_value_keyword) != string::npos)
            {
                getline(infile, line);
                std::istringstream iss(line);
                iss >> load_name >> dof_id >> value_str;
                dof_id.erase(dof_id.end() - 1); // 删除字符串最后的符号
                dof_id_ = atoi(dof_id.c_str()); // 转换字符串为int
                value_ = stod(value_str);       // 转换字符串为double
                data_cae.load_dof_ = dof_id_;
                data_cae.load_value_ = value_;
                break;
            }
        }
        infile.close();
        cout << "the information of load boundary (" << data_cae.load_set_.size() << ") have been readed." << endl;
        cout << "the fixed DOF is " << data_cae.load_dof_ << ", and the value is " << data_cae.load_value_ << endl;
    }

    // 读取约束信息
    void ReadInfo::read_dis_bcs(string dis_set_keyword, data_management &data_cae)
    {
        // 读取计算文件
        std::ifstream infile(path_.c_str(), std::ios::in);
        string line;
        // 读取载荷节点集合
        while (getline(infile, line))
        {
            string node_id;
            int node_id_;
            if (line.find(dis_set_keyword) != string::npos)
            {
                while (getline(infile, line))
                {
                    if (line.find("*") != string::npos)
                        break;
                    else
                    {
                        std::istringstream iss(line);
                        while (iss >> node_id)
                        {
                            if (node_id.find(",") != string::npos)
                            {
                                node_id.erase(node_id.end() - 1); // 删除字符串最后的符号
                            }
                            node_id_ = atoi(node_id.c_str()); // 转换字符串为int
                            data_cae.dis_bc_set_.push_back(node_id_);
                        }
                    }
                }
            }
            if (line.find("*End Assembly") != string::npos)
            {
                break;
            }
        }
        infile.close();
        cout << "the information of displacement boundary (" << data_cae.dis_bc_set_.size() << ")  have been readed." << endl;
    }

    void ReadInfo::del_blank(string& str)
    {
        int index = 0;
        if (!str.empty())
        {
            while ((index = str.find(' ', index)) != string::npos)
            {
                str.erase(index, 1);
            }
        }
    }
}
