import pandas as pd
import numpy as np

# Vms = np.zeros((19393, 1))
# vms_f = np.asarray(pd.read_fwf("./try_1/VMS.rpt", encoding='GB2312'))
# flag = False
# n_vms = len(vms_f)
# id_ = 0
# for i in range(n_vms):
#     temp = vms_f[i, 0].split()
#     if flag and len(temp) > 1:
#         Vms[id_, 0] = float(temp[4])
#         id_ = id_+1
#     if ('Mises' in temp) and ('积分点' in temp):
#         flag = True
# with open('Abaqus_VMS.txt', 'w') as f:
#     for i in range(len(Vms)):
#         f.write(str(Vms[i, 0]) + '\n')

U = []
node_id = []
u_f = np.asarray(pd.read_fwf("E:/CADCAE_project/OpenFEA/model/mix_ele_model/shaft_bracket/abaqus_U123.rpt", encoding='GB2312'))
flag = False
n_U = len(u_f)
id_ = 0
for i in range(n_U):
    temp = u_f[i, 0].split()
    if flag and len(temp) > 1:
        if temp[2] not in node_id:
            U.append(float(temp[4]))
            U.append(float(temp[5]))
            U.append(float(temp[6]))
            node_id.append(temp[2])
        else:
            pass
    if ('U1' in temp) and ('U2' in temp) and ('U3' in temp):
        flag = True


with open('Abaqus_U.txt', 'w') as f:
    for i in range(len(U)):
        f.write(str(U[i]) + '\n')
