//
// Created by Yanjia on 2024/5/9.
//


#include "hvector.h"


arma::vec crossProduct(const arma::vec &a, const arma::vec &b) {
    return arma::vectorise(arma::cross(a.subvec(0, 2), b.subvec(0, 2)));
}

// 辅助函数，计算四边形单元的面法向量
arma::vec computeFaceNormal(const arma::mat &Nid, int id1, int id2, int id3, int id4) {
    // 创建一个映射，将节点ID映射到它们在 Nid 矩阵中的行号
    std::unordered_map<int, int> nodeIdToRowMap;

    for (arma::uword i = 0; i < Nid.n_rows; ++i) {
        int nodeId = Nid(i, 0);
        nodeIdToRowMap[nodeId] = i; //
    }

// 现在 nodeIdToRowMap 已经被填充，可以用于查找节点ID对应的行号

    int rowId1 = nodeIdToRowMap[id1];
    int rowId2 = nodeIdToRowMap[id2];
    int rowId3 = nodeIdToRowMap[id3];
    int rowId4 = nodeIdToRowMap[id4];

//submat提取的是一个行向量，而vec接收的是一个列向量。因此用t()将行向量转换为列向量。
    arma::vec v1 = Nid.submat(rowId2, 1, rowId2, 3).row(0).t() - Nid.submat(rowId1, 1, rowId1, 3).row(0).t();
    arma::vec v2 = Nid.submat(rowId3, 1, rowId3, 3).row(0).t() - Nid.submat(rowId1, 1, rowId1, 3).row(0).t();
    return arma::normalise(crossProduct(v1, v2));
}

// 计算给定节点的法向量
arma::vec computeNodeNormal(const arma::imat &Eid, const arma::mat &Nid, int nodeId) {
    arma::vec normalSum(3, arma::fill::zeros);
    int numContributions = 0;

    // 遍历所有单元，检查是否包含给定节点
    for (int i = 0; i < Eid.n_rows; ++i) {//0~13
        bool found = false;
        for (int j = 1; j < Eid.n_cols; ++j) {//1~8
            if (Eid(i, j) == nodeId) {
                // 根据i行数就可以找到对应的单元号，以及该单元的所有节点
                found = true;
                break;
            }
        }
        if (!found) continue;

        // 找到包含当前节点的单元的所有节点编号
        int id1 = Eid(i, 1);
        int id2 = Eid(i, 2);
        int id3 = Eid(i, 3);
        int id4 = Eid(i, 4); // 假设Eid只有4个节点编号，不考虑中间的节点

        // 计算单元法向量并累加到总和中
        normalSum += computeFaceNormal(Nid, id1, id2, id3, id4);
        ++numContributions;
    }

    // 标准化法向量并返回
    return arma::normalise(normalSum / numContributions);
}

arma::mat nodeNormal(imat Eid, mat Nid, double nn) {

//    mat H3i;
    mat H3i(nn, 3, fill::zeros);
    for (int i = 0; i < Nid.n_rows; ++i) {

        H3i.row(i) = computeNodeNormal(Eid, Nid, Nid(i, 0)).t(); // 假设节点编号在Nid的第一列
    }

    return H3i;

}

HVector hVector(imat Eid, mat Nid, double nn) {


    HVector hv;
    mat H3i = nodeNormal(Eid, Nid, nn);
    vec ii = vec({1, 0, 0}); // x-axis
    mat H2i(nn, 3, fill::zeros);

    for (int i = 0; i < nn; ++i) {
        vec temp = trans(cross(H3i.row(i), ii));
        H2i.row(i) = temp.t() / norm(temp);
    }
    mat H1i(nn, 3, fill::zeros);
    for (int i = 0; i < nn; ++i) {
        vec temp = trans(cross(H3i.row(i), H2i.row(i)));
        H1i.row(i) = temp.t() / norm(temp);
    }

    hv.H3i = H3i;
    hv.H2i = H2i;
    hv.H1i = H1i;

    return hv;
}