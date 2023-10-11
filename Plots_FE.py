import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import hFE_functions as hFEf


def read_RFnodeFile_opt(file_):
    # read data from text file
    df_ = np.loadtxt(file_, delimiter=',')
    rfy_ = np.array(df_[:, 1])
    uy_ = np.array(df_[:, 2])
    return uy_, rfy_


def read_icotec_experiment(sheet_):
    file_ = '/home/biomech/Documents/01_Icotec/01_Experiments/99_Others/Screw_Bending_test/F1798Cantilever_UniBern.xls'
    f_ = pd.read_excel(file_, sheet_name=sheet_, header=[1, 2], )
    return f_


# %%
peek_samples = [2, 5, 7, 8, 10, 13, 15, 16, 18, 21, 23, 24, 26, 29, 31, 32]  # with 24
ti_samples = [3, 4, 6, 9, 11, 12, 14, 17, 19, 20, 22, 25, 27, 28, 30, 33]  # with 25
FE_points = [10, 31, 52, 73, 94, 115, 136]


# %% PEEK material optimisation
path = '/home/biomech/Documents/01_Icotec/02_FEA/00_Model/'
# '' = 18000 GPa
#  1 = 20000 GPa
#  2 = 25000 GPa (orig)
#  3 = 25000 GPa (44 --> 50 mm lever)
#  4 = 30000 GPa (44 --> 50 mm lever)
#  5 = 28000 GPa (44 --> 50 mm lever)
#  6 = 28000 GPa, plasticity (0/280, 0.03/300), (44 --> 50 mm lever)
#  7 = 28000 GPa, plasticity (0/200, 0.03/220), (44 --> 50 mm lever)
#  8 = 30000 GPa, plasticity (0/200, 0.03/220), (44 --> 50 mm lever)
#  9 = 30000 GPa, plasticity (0/150, 0.03/180), (44 --> 50 mm lever)
# 10 = 35000 GPa, plasticity (0/180, 0.03/200), (44 --> 50 mm lever)
ext = ['', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10']

plt.figure()
for i in [3, 10]:
    file = '96_screw_Osteoporosis_new_Bending_RFnode_' + ext[i] + '.txt'
    file2 = '97_screw_Osteoporosis_new_Bending_screw_no-nlgeom_RFnode_' + ext[i] + '.txt'
    [uy, rfy] = read_RFnodeFile_opt(path + file)
    [uy2, rfy2] = read_RFnodeFile_opt(path + file2)
    if i < len(ext) - 1:
        n = np.where(-rfy > 5)[0][0]
        plt.plot((-uy - -uy[n])*1.5, -rfy*1.5, label='FE: linear', color='#2a52be')
        n2 = np.where(-rfy2 > 5)[0][0]
        # plt.plot(-uy2 - -uy2[n2], -rfy2, label='FE: Screw only', color='#e77471', alpha=0)
    else:
        n = np.where(-rfy > 5)[0][0]
        plt.plot(-uy - -uy[n], -rfy, label='FE: optimised', color='#045f5f')
        n2 = np.where(-rfy2 > 5)[0][0]
        # plt.plot(-uy2 - -uy2[n2], -rfy2, label='FE material model', color='#045f5f')
data = {}
plt.xlabel('Displacement / mm')
plt.ylabel('Force / N')
plt.legend()
tests = ['1.1', '1.2', '1.3', '1.4', '1.5', '1.6', '2.1', '2.2', '2.3']
for i in range(len(tests)):
    data[tests[i]] = read_icotec_experiment(tests[i])
    # if i < 6:
    #     plt.plot(data[tests[i]]['Dehnung'], data[tests[i]]['Standardkraft'], color='#045f5f')
    # else:
    #     plt.plot(data[tests[i]]['Dehnung'], data[tests[i]]['Standardkraft'], color='#e77471')
mean_x = np.mean((data[tests[0]]['Dehnung'][:12000],
                  data[tests[1]]['Dehnung'][:12000],
                  data[tests[2]]['Dehnung'][:12000],
                  data[tests[3]]['Dehnung'][:12000],
                  data[tests[4]]['Dehnung'][:12000],
                  data[tests[5]]['Dehnung'][:12000],
                  data[tests[6]]['Dehnung'][:12000],
                  data[tests[7]]['Dehnung'][:12000],
                  data[tests[8]]['Dehnung'][:12000],), axis=0)
mean_y = np.mean((data[tests[0]]['Standardkraft'][:12000],
                  data[tests[1]]['Standardkraft'][:12000],
                  data[tests[2]]['Standardkraft'][:12000],
                  data[tests[3]]['Standardkraft'][:12000],
                  data[tests[4]]['Standardkraft'][:12000],
                  data[tests[5]]['Standardkraft'][:12000],
                  data[tests[6]]['Standardkraft'][:12000],
                  data[tests[7]]['Standardkraft'][:12000],
                  data[tests[8]]['Standardkraft'][:12000],), axis=0)
plt.plot(mean_x, mean_y, color='r', label='Mean experiments icotec')
plt.legend()
# %% uFE result and energy
sample_list = open('/home/biomech/Documents/01_Icotec/Specimens.txt', 'r').read().splitlines()
no = 7
specimen = sample_list[no]
print(specimen)
number = '16'
fe_path = '/home/biomech/DATA/01_Icotec/02_FEA/02_uFE/01_MainStudy/' + specimen + '/'
fe = number + '_' + specimen + '_0121399_'

file_exp = '/home/biomech/Documents/01_Icotec/01_Experiments/00_Data/01_MainStudy/' + specimen + '_resample.csv'
U, _ = hFEf.read_RFnodeFile(fe_path + fe + 'RFnode.txt')
# U, _ = read_RFnodeFile('/home/biomech/DATA/01_Icotec/02_FEA/02_uFE/Tests/test_B_0T_' + 'RFnode.txt')
_, F = hFEf.read_RFnodeFile(fe_path + fe + 'RFnodeFix.txt')
# _, F = read_RFnodeFile('/home/biomech/DATA/01_Icotec/02_FEA/02_uFE/Tests/test_B_0T_' + 'RFnodeFix.txt')
[_, A_y, _, _, _, _, a_y, a_f, _] = hFEf.read_resample(file_exp)
energy = 0
if energy:
    t_ke, ke = hFEf.read_energy(fe_path + fe + 'ke.txt')
    t_ie, ie = hFEf.read_energy(fe_path + fe + 'ie.txt')
    plt.figure()
    # plt.plot(t_ke, ke)
    # plt.plot(t_ie, ie)
    plt.plot(t_ke, ke/ie)
    plt.plot([t_ke[0], t_ke[-1]], [0.1, 0.1])
plt.figure()
if no in [2, 5, 7, 8, 10, 13, 15, 16, 18, 21, 23, 24, 26, 29, 31, 32]:  # without 0
    plt.title(number + '_' + specimen + ' (PEEK)')
    print('PEEK')
elif no in [3, 4, 6, 9, 11, 12, 14, 17, 19, 20, 22, 25, 27, 28, 30, 33]:  # without 1
    plt.title(number + '_' + specimen + ' (Ti)')
    print('Ti')

plt.plot(a_y, a_f-a_f[0])
plt.plot(U, -F)
