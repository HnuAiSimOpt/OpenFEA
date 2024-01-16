# -------------------------------------------------------------------------
# Copyright:  WH team
# Author: YinJichao <jichaoyinyjc@163.com>
# Completion date:  XXX
# Description: XXX
# -------------------------------------------------------------------------
import numpy as np
import math
import sys
from scipy.sparse import coo_matrix


def get_dof(e_node_id):
    e_dof = np.zeros(8, dtype=int)
    for j in range(4):
        node_id = e_node_id[j] - 1  # 从0开始计数
        e_dof[2 * j] = 2 * node_id  # X 方向
        e_dof[2 * j + 1] = 2 * node_id + 1  # Y 方向
    return e_dof


class GeoNL3D:
    def __init__(self):
        self.ne_ = None
        self.nd_ = None
        self.em_ = None
        self.nu_ = None
        self.D0_ = None
        self.Kt_row = None
        self.Kt_col = None
        self.dN_ = {}  # 高斯积分点的形函数偏导
        self.BL0_ = {}  # 高斯积分点的线性项（BL0）
        self.BL1_ = {}  # 高斯积分点的线性项（BL1）
        self.BNL_ = {}  # 高斯积分点的二次项（BNL）
        self.pk2s_ = {}  # 高斯积分点的PK-II 应力
        self.det_jacobi = {}  # 高斯积分点的雅可比矩阵行列式

    def geo_nonlinear_analysis(self, node, element, load, fixed, em, nu, load_axis, load_value):
        # 设置材料
        self.set_mat(em, nu)
        # 建立本构矩阵
        self.get_constitutive_model()
        # 初始化参数
        self.ne_ = len(element)
        self.nd_ = len(node)
        dis_nl = np.zeros(2 * self.nd_)
        # 计算切线刚度矩阵
        f_int, KT = self.get_global_stiff_f_int(node, element, dis_nl)
        #

    def build_sp_idx(self, element):
        n_dof_ele, n_value = 8, 64
        self.Kt_row = np.zeros(self.ne_ * n_value, dtype=int)
        self.Kt_col = np.zeros(self.ne_ * n_value, dtype=int)
        for i in range(self.ne_):
            e_node_id = element[i]
            e_dof = get_dof(e_node_id)
            row_temp = np.zeros(n_value)
            col_temp = np.zeros(n_value)
            for m in range(n_dof_ele):
                for n in range(n_dof_ele):
                    id_ = m * n_dof_ele + n
                    row_temp[id_] = e_dof[m]
                    col_temp[id_] = e_dof[n]
            self.Kt_row[n_value * i: n_value * (i + 1)] = row_temp
            self.Kt_col[n_value * i: n_value * (i + 1)] = col_temp

    def get_global_stiff_f_int(self, node, element, pre_dis_nl):
        self.build_sp_idx(element)
        #
        n_k_value = 64
        value_kt = np.zeros(self.ne_ * n_k_value)
        f_int = np.zeros((3 * self.nd_, 1))
        for i in range(self.ne_):
            # 点集坐标
            points_e = node[element[i] - 1]
            pre_dis_e = np.vstack((pre_dis_nl[2 * (element[i] - 1)], pre_dis_nl[2 * (element[i] - 1) + 1]))
            # 计算性函数的偏导，高斯积分点的雅可比矩阵行列式，应力转换矩阵
            self.get_dN(points_e)
            self.get_trans_mat(pre_dis_e)
            # 计算小位移下的单元刚度矩阵
            KL0_e = self.get_linear_stiff()
            # 计算位移引起的单元刚度矩阵
            KL1_e = self.get_dis_stiff(pre_dis_e)
            # 计算应力引起的单元刚度矩阵
            KNL_e = self.get_stress_stiff(pre_dis_e)
            #
            KT_e = KL0_e + KL1_e + KNL_e
            value_kt[n_k_value * i: n_k_value * (i + 1)] = KT_e.reshape([n_k_value, 1], order='F').flatten()
            # 计算内力
            f_int_e = self.get_f_int()
            e_dof = get_dof(element[i])
            f_int[e_dof] = f_int[e_dof] + f_int_e
        KT = coo_matrix((value_kt, (self.Kt_row, self.Kt_col)), shape=(3 * self.nd_, 3 * self.nd_)).tocsc()
        return f_int, KT

    def get_linear_stiff(self):
        KL0_e = np.zeros((8, 8))
        for i in range(4):
            name = 'point_' + str(i + 1)
            det_jacobi = self.det_jacobi[name]
            BL0 = self.BL0_[name]
            KL0_e_i = (BL0.T @ self.D0_) @ BL0
            KL0_e = KL0_e + det_jacobi * KL0_e_i
        return KL0_e

    def get_dis_stiff(self, dis):
        KL1_e = np.zeros((8, 8))
        for i in range(4):
            # 形函数偏导，雅可比行列式
            name = 'point_' + str(i + 1)
            det_jacobi = self.det_jacobi[name]
            #
            BL0, BL1 = self.BL0_[name], self.BL1_[name]
            KL1_e_i_1 = (BL0.T @ self.D0_) @ BL1
            KL1_e_i_2 = (BL1.T @ self.D0_) @ BL0
            KL1_e_i_3 = (BL1.T @ self.D0_) @ BL1
            KL1_e = KL1_e + det_jacobi * (KL1_e_i_1 + KL1_e_i_2 + KL1_e_i_3)
        return KL1_e

    def get_stress_stiff(self, dis):
        Ks_e = np.zeros((8, 8))
        dis_vec = np.array([dis[0, 0], dis[1, 0], dis[0, 1], dis[1, 1],
                            dis[0, 2], dis[1, 2], dis[0, 3], dis[1, 3]])
        pk2s_mat = np.zeros((4, 4))
        for i in range(4):
            # 形函数偏导，雅可比行列式
            name = 'point_' + str(i + 1)
            det_jacobi = self.det_jacobi[name]
            #
            BL = self.BL0_[name] + self.BL0_[name]
            # PK-II 应力
            self.pk2s_[name] = BL @ dis_vec
            pk2s_mat[0, 0], pk2s_mat[0, 1] = self.pk2s_[name][0], self.pk2s_[name][2]
            pk2s_mat[1, 0], pk2s_mat[1, 1] = self.pk2s_[name][2], self.pk2s_[name][1]
            pk2s_mat[2, 2], pk2s_mat[2, 3] = self.pk2s_[name][0], self.pk2s_[name][2]
            pk2s_mat[3, 2], pk2s_mat[3, 3] = self.pk2s_[name][2], self.pk2s_[name][1]
            #
            BNL = self.BNL_[name]
            KNL_e_i = (BNL.T @ pk2s_mat) @ BNL
            Ks_e = Ks_e + det_jacobi * KNL_e_i
        return Ks_e

    def get_f_int(self):
        f_int_e = np.zeros((8, 1))
        for i in range(4):
            # 形函数偏导，雅可比行列式
            name = 'point_' + str(i + 1)
            det_jacobi = self.det_jacobi[name]
            #
            BL = self.BL0_[name] + self.BL0_[name]
            pk2s = self.pk2s_[name]
            f_int_e_i = BL.T @ pk2s
            f_int_e = f_int_e + det_jacobi * f_int_e_i.reshape((8, 1))
        return f_int_e

    def set_mat(self, em, nu):
        self.em_ = em
        self.nu_ = nu

    def get_constitutive_model(self):
        # 判断材料是否存在
        if self.em_ is None:
            print("em is None")
            sys.exit(1)
        if self.nu_ is None:
            print("nu is None")
            sys.exit(1)
        # 计算本构矩阵
        temp_1 = self.em_ / (1. - self.nu_ ** 2)
        temp_2 = (1. - self.nu_) / 2.
        self.D0_ = temp_1 * np.array([[1., self.nu_, 0.],
                                      [self.nu_, 1., 0.],
                                      [0., 0., temp_2]])

    def get_dN(self, point):
        gauss_points = 1. / math.sqrt(3.) * np.array([[-1, -1], [1, -1], [1, 1], [-1, 1]])
        for i in range(4):
            s, t = gauss_points[i, 0], gauss_points[i, 1]
            dNdr = 0.25 * np.array([[-(1 - t), 1 - t, 1 + t, -(1 + t)],
                                    [-(1 - s), -(1 + s), 1 + s, 1 - s]])
            jacobi = dNdr @ point
            inv_jacbi = np.linalg.inv(jacobi)
            dN_i = inv_jacbi @ dNdr
            name = 'point_' + str(i + 1)
            self.dN_[name] = dN_i
            det_jacobi = np.linalg.det(jacobi)
            self.det_jacobi[name] = det_jacobi

    def get_trans_mat(self, dis):
        for i in range(4):
            # 形函数偏导，雅可比行列式
            name = 'point_' + str(i + 1)
            # 应变转换矩阵预处理
            dN_i = self.dN_[name]
            dNdx, dNdy = dN_i[0, :].reshape((4, 1)), dN_i[1, :].reshape((4, 1))
            #
            ux, uy = dis[0, :].reshape((4, 1)), dis[1, :].reshape((4, 1))
            Lxx, Lxy, Lyx, Lyy = dNdx.T @ ux, dNdx.T @ uy, dNdy.T @ ux, dNdy.T @ uy
            #
            BL0, BL1, BNL = np.zeros((3, 8)), np.zeros((3, 8)), np.zeros((4, 8))
            for j in range(4):
                #
                BL0[0, 2 * j] = dN_i[0, j]
                BL0[1, 2 * j + 1] = dN_i[1, j]
                BL0[2, 2 * j], BL0[2, 2 * j + 1] = dN_i[1, j], dN_i[0, j]
                #
                BL1[0, 2 * j], BL1[0, 2 * j + 1] = Lxx * dN_i[0, j], Lyx * dN_i[0, j]
                BL1[1, 2 * j], BL1[1, 2 * j + 1] = Lxy * dN_i[1, j], Lyy * dN_i[1, j]
                BL1[2, 2 * j] = Lxx * dN_i[1, j] + Lxy * dN_i[1, j]
                BL1[2, 2 * j + 1] = Lyx * dN_i[0, j] + Lyy * dN_i[0, j]
                #
                BNL[0, 2 * j], BNL[1, 2 * j] = dN_i[0, j], dN_i[1, j]
                BNL[2, 2 * j + 1], BNL[2, 2 * j + 1] = dN_i[0, j], dN_i[1, j]
            self.BL0_[name] = BL0
            self.BL1_[name] = BL1
            self.BNL_[name] = BNL
