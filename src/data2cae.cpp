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
        string line, temp_node_topo;
        std::vector<string> type_temp;
        int id_node = 0, id_ele = 0, ele_type_idx = 0, I_ele_type, node_per_ele;
        bool node_record = false, node_topo_record = false;
        while (getline(infile, line))
        {
            /* --------------------- 读取节点坐标 --------------------- */
            if (line.find("*Node") != string::npos)
            {
                node_record = true;
                getline(infile, line);
            }
            if (node_record & line.find("*") != string::npos)
            {
                node_record = false;
            }
            if (node_record)
            {
                string id, x, y, z;
                double id_, x_, y_, z_;
                std::istringstream iss(line);
                iss >> id >> x >> y >> z;
                x.erase(x.end() - 1); // 删除字符串最后的符号
                x_ = stod(x);         // 转换字符串为double
                y.erase(y.end() - 1);
                y_ = stod(y);
                z_ = stod(z);
                id.erase(id.end() - 1);
                id_ = stod(id);
                data_cae.coords_[id_node][0] = x_;
                data_cae.coords_[id_node][1] = y_;
                data_cae.coords_[id_node][2] = z_;
                data_cae.coords_idx_[id_] = id_node;
                id_node = id_node + 1;
            }
            /* --------------------- 读取单元类型及节点拓扑关系 --------------------- */
            if (line.find("*Element") != string::npos)
            {
                // 读取单元类型
                type_temp = split_str(line, "=");
                I_ele_type = del_blank(type_temp[1], data_cae.ELE_TYPES);
                data_cae.ele_class_.insert(I_ele_type);
                //
                switch (I_ele_type)
                {
                case 0:
                    node_per_ele = 4;
                    break;
                case 1:
                    node_per_ele = 8;
                    break;
                default:
                    cout << "Error in element type !!!\n";
                    break;
                }
                node_topo_record = true;
                getline(infile, line);
            }
            if (node_topo_record & line.find("*") != string::npos)
            {
                node_topo_record = false;
            }
            if (node_topo_record)
            {
                std::istringstream iss(line);
                // 储存单元类型
                data_cae.ele_list_idx_[id_ele] = I_ele_type;
                int temp_id = -1;
                while (getline(iss, temp_node_topo, ','))
                {
                    del_blank(temp_node_topo);
                    if (temp_id >= 0)
                        data_cae.node_topos_[id_ele][temp_id] = atoi(temp_node_topo.c_str());
                    temp_id += 1;
                }
                id_ele += 1;
            }
            if (line.find("*End Assembly") != string::npos)
            {
                break;
            }
        }
        // 单元类型在set自动去重，初始化单元
        for (auto temp_ele_type : data_cae.ele_class_)
        {
            string temp_ele_type_str;
            for (std::map<string, int>::iterator it = data_cae.ELE_TYPES.begin(); it != data_cae.ELE_TYPES.end(); it++)
            {
                if (it->second == temp_ele_type)
                    temp_ele_type_str = it->first;
            }
            data_cae.add_ele(temp_ele_type_str);
            cout << "The element " << temp_ele_type_str << " has been built !!!" << endl;
        }
        cout << "a total of " << id_node << " nodal coordinates have been readed." << endl;
        cout << "a total of " << id_ele << " nodal connectivities have been readed." << endl;
        cout << "a total of " << ele_type_idx + 1 << " types of elements have been readed." << endl;
    }

    // 读取非协调信息
    void ReadInfo::readNconformingMessage(data_management &data_cae)
    {
        std::ifstream infile(path_.c_str(), std::ios::in);
        string line;
        while (getline(infile, line))
        {
            string node_id;
            int node_id_;
            // 读取细网格单元8节点信息
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
                                node_id.erase(node_id.end() - 1); // delete the last symbol in the string
                            }
                            node_id_ = atoi(node_id.c_str()); // convert to int type
                            data_cae.BndMesh_F.push_back(node_id_);
                        }
                    }
                }
            }
            // 读取粗网格单元8节点信息
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
                                node_id.erase(node_id.end() - 1); // delete the last symbol in the string
                            }
                            node_id_ = atoi(node_id.c_str()); // convert to int type
                            data_cae.BndMesh_C.push_back(node_id_);
                        }
                    }
                }
            }

            //读取细网格 交界面（四边形、三角形）节点信息
              std::vector<string> facenode_temp;
            if (line.find("*BndFace_finemesh") != string::npos)
            {  
                facenode_temp = split_str(line, "=");
                int face_nnode;
                del_blank(facenode_temp[1]);
                face_nnode = std::stoi(facenode_temp[1]);
                
                 int r = data_cae.BndMesh_F.size();
                 data_cae.bndFace_finemesh.resize(r, vector<int>(face_nnode));

                int id_node = 0;
                while (getline(infile, line))
                {
                    if (line.find("*") != string::npos)
                    {
                        break;
                    }
                        
                    std::istringstream iss(line);
                    vector<string> tempN_a;
                    string tempN;
                    while (getline(iss, tempN, ','))
                    {
                        del_blank(tempN);
                        tempN_a.push_back(tempN);
                    }
                    
                    for (int i = 0; i < face_nnode; i++)
                    {
                        data_cae.bndFace_finemesh[id_node][i] = atoi(tempN_a[i].c_str());
                    }
                    id_node++;
                }
            }
            if (line.find("*BndFace_coarsemesh") != string::npos)
            {
                facenode_temp = split_str(line, "=");
                int face_nnode;
                del_blank(facenode_temp[1]);
                face_nnode = std::stoi(facenode_temp[1]);

                int r = data_cae.BndMesh_C.size();
                data_cae.bndFace_coarsemesh.resize(r, vector<int>(face_nnode));

                int id_node = 0;
                while (getline(infile, line))
                {
                    if (line.find("*") != string::npos)
                    {
                        break;
                    }

                    std::istringstream iss(line);
                    vector<string> tempN_a;
                    string tempN;
                    while (getline(iss, tempN, ','))
                    {
                        del_blank(tempN);
                        tempN_a.push_back(tempN);
                    }

                    for (int i = 0; i < face_nnode; i++)
                    {
                        data_cae.bndFace_coarsemesh[id_node][i] = atoi(tempN_a[i].c_str());
                    }
                    id_node++;
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

    void ReadInfo::del_blank(string &str)
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

    int ReadInfo::del_blank(string &str, map<string, int> &ele_map)
    {
        int index = 0;
        if (!str.empty())
        {
            while ((index = str.find(' ', index)) != string::npos)
            {
                str.erase(index, 1);
            }
        }
        int I_ele_type = ele_map[str];
        return I_ele_type;
    }

    void CAE::ReadInfo::CA_data_convert(data_management &data_cae_Im, data_management &data_cae)
    {
        // 只读不改data_cae_Im，构建一个新的data_cae
        /* 
        预留修改坐标处理
        修改模型的节点坐标我读进来了，因为我只知道删除，所以下面直接把旧模型的节点坐标复制到新的，新的相当于没用上
        如果是修改坐标的情况你看下要怎么处理。
        data_managment中加入了这个成员，存储修改节点信息。
        map<int, int> coords_idx_;// 节点编号索引(文件编号：coords_的行索引）
        */
        data_cae.coords_ = data_cae_Im.coords_;
        data_cae.nd_ = data_cae.coords_.size();
        // 自由度索引
        data_cae.resort_free_nodes_ = data_cae_Im.resort_free_nodes_;
        // 荷载列阵
        data_cae.single_load_vec_ = data_cae_Im.single_load_vec_;
        // 无约束自由度的位移向量
        data_cae.single_dis_vec_ = data_cae_Im.single_dis_vec_;
        // paidiso求解器对象
        data_cae.item_pardiso = data_cae_Im.item_pardiso;
        // 查看CA代码仅用到上面几个量，其它量并未拷贝，如后续有需求可加上
    }

}// namespace