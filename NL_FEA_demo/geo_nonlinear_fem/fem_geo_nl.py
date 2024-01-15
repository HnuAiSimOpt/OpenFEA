# -------------------------------------------------------------------------
# Copyright:  WH team
# Author: YinJichao <jichaoyinyjc@163.com>
# Completion date:  XXX
# Description: XXX
# -------------------------------------------------------------------------
import numpy as np
import math
import sys


class GeoNL3D:
    def __init__(self):
        self.ne_ = None
        self.nd_ = None
        self.em_ = None
        self.nu_ = None
        self.D0_ = None
        self.dN_ = {}
        self.det_jacobi = {}

    def set_mat(self, em, nu):
        self.em_ = em
        self.nu_ = nu

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
        KT = self.get_global_stiff(node, element, dis_nl)

    def get_global_stiff(self, node, element, pre_dis_nl):
        for i in range(self.ne_):
            # 点集坐标
            points_e = node[element[i] - 1]
            pre_dis_e = np.vstack((pre_dis_nl[2 * (element[i] - 1)], pre_dis_nl[2 * (element[i] - 1) + 1]))
            # 计算性函数的偏导，高斯积分点的雅可比矩阵行列式
            self.get_dN(points_e)
            # 计算小位移下的单元刚度矩阵
            KL0_e = self.get_linear_stiff()
            # 计算位移引起的单元刚度矩阵
            KL1_e = self.get_dis_stiff(pre_dis_e)
            # 计算应力引起的单元刚度矩阵
            KNL_e = self.get_stress_stiff(dis)
            #
            KT_e = KL0_e + KL1_e
        return KT_e

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

    def get_linear_stiff(self):
        KL0_e = np.zeros((8, 8))
        for i in range(4):
            name = 'point_' + str(i + 1)
            dN_i = self.dN_[name]
            det_jacobi = self.det_jacobi[name]
            BL0 = np.zeros((3, 8))
            for j in range(4):
                BL0[0, 2 * j] = dN_i[0, j]
                BL0[1, 2 * j + 1] = dN_i[1, j]
                BL0[2, 2 * j], BL0[2, 2 * j + 1] = dN_i[1, j], dN_i[0, j]
            KL0_e_i = (BL0.T @ self.D0_) @ BL0
            KL0_e = KL0_e + det_jacobi * KL0_e_i
        return KL0_e

    def get_dis_stiff(self, dis):
        KL1_e = np.zeros((8, 8))
        for i in range(4):
            # 形函数偏导，雅可比行列式
            name = 'point_' + str(i + 1)
            dN_i = self.dN_[name]
            det_jacobi = self.det_jacobi[name]
            # 应变转换矩阵预处理
            dNdx, dNdy = dN_i[0, :].reshape((4, 1)), dN_i[1, :].reshape((4, 1))
            ux, uy = dis[0, :].reshape((4, 1)), dis[1, :].reshape((4, 1))
            Lxx, Lxy, Lyx, Lyy = dNdx.T @ ux, dNdx.T @ uy, dNdy.T @ ux, dNdy.T @ uy
            #
            BL0, BL1 = np.zeros((3, 8)), np.zeros((3, 8))
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
            KL1_e_i_1 = (BL0.T @ self.D0_) @ BL1
            KL1_e_i_2 = (BL1.T @ self.D0_) @ BL0
            KL1_e_i_3 = (BL1.T @ self.D0_) @ BL1
            KL1_e = KL1_e + det_jacobi * (KL1_e_i_1 + KL1_e_i_2 + KL1_e_i_3)
        return KL1_e

    def get_stress_stiff(self, dis):
        return 0.
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
