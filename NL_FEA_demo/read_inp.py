# -------------------------------------------------------------------------
# Copyright:  WH team
# Author: YinJichao <jichaoyinyjc@163.com>
# Completion date:  XXX
# Description: XXX
# -------------------------------------------------------------------------
import numpy as np


def read_inp(path):
    # read the inp file
    inp_file = open(path)
    lines = inp_file.readlines()  # 按行读取
    inp_file.close()
    # 定义变量
    node, element, load_, fixed_ = [], [], [], []
    node_flag, ele_flag, load_flag, fixed_flag = False, False, False, False
    set_generate = False
    for _, line in enumerate(lines):
        if node_flag:
            if '*' in line:
                break
            temp = line.split(',')
            list_ = []
            for id_, value in enumerate(temp):
                list_.append(float(value.split()[0]))
            node.append(np.array(list_))
        if '*Node' in line:
            node_flag = True
    node = np.asarray(node)
    # 读取连接关系
    for _, line in enumerate(lines):
        if ele_flag:
            if '*' in line:
                break
            temp = line.split(',')
            list_ = []
            for id_, value in enumerate(temp):
                list_.append(int(value.split()[0]))
            element.append(np.array(list_))
        if '*Element' in line:
            ele_flag = True
    element = np.asarray(element)
    # 读取载荷节点
    for _, line in enumerate(lines):
        if load_flag:
            if '*' in line:
                break
            temp = line.split(',')
            for id_, value in enumerate(temp):
                load_.append(int(value.split()[0]))
        if 'Set-load' in line:
            load_flag = True
        if 'generate' in line:
            set_generate = True
    load = None
    if set_generate:
        load = np.arange(load_[0], load_[1] + 1, load_[2], dtype=int)
    else:
        load = load_.copy()
    set_generate = False
    load = np.asarray(load)
    # 读取被固定节点
    for _, line in enumerate(lines):
        if fixed_flag:
            if '*' in line:
                break
            temp = line.split(',')
            for id_, value in enumerate(temp):
                fixed_.append(int(value.split()[0]))
        if 'Set-fix' in line:
            fixed_flag = True
        if 'generate' in line:
            set_generate = True
    fixed = None
    if set_generate:
        fixed = np.arange(fixed_[0], fixed_[1] + 1, fixed_[2], dtype=int)
    else:
        fixed = fixed_.copy()
    set_generate = False
    fixed = np.asarray(fixed)
    # 读取材料
    mat_flag = False
    em = 0.
    nu = 0.
    for _, line in enumerate(lines):
        if mat_flag:
            temp = line.split(',')
            em = float(temp[0])
            nu = float(temp[1])
            mat_flag = False
        if '*Elastic' in line:
            mat_flag = True
    # 读取载荷方向及大小
    value_flag = False
    load_axis = 0
    load_value = 0.
    for _, line in enumerate(lines):
        if value_flag:
            temp = line.split(',')
            load_axis = int(temp[1])
            load_value = float(temp[2])
            value_flag = False
        if '*Cload' in line:
            value_flag = True
    return node, element, load, fixed, em, nu, load_axis, load_value
