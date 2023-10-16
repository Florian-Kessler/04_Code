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

samples = [5, 7, 8, 10, 15, 16, 18, 21]  # 16 no hfe
F_extrem = np.zeros((6, 22))
#                    [0, 1, 2, 3, 4,     5, 6,      7,     8, 9,    10, 11, 12, 13, 14,     15,    16, 17,    18, 19, 20,    21]
peaks_4mm = np.array([0, 0, 0, 0, 0, -92.6, 0, -101.5, -76.1, 0, -65.4,  0,  0,  0,  0, -114.0, -91.9,  0, -76.0,  0,  0, -80.2])
for no in samples:
    specimen = sample_list[no]
    print(specimen)
    number = '17'
    fe_path = '/home/biomech/DATA/01_Icotec/02_FEA/02_uFE/01_MainStudy/' + specimen + '/'
    fe = number + '_' + specimen + '_0121399_'
    hfe_path = '/home/biomech/Documents/01_Icotec/02_FEA/01_MainStudy/' + specimen + '/60_L50_S50_D45/60_L50_S50_D45_d1_05_P_'

    file_exp = '/home/biomech/Documents/01_Icotec/01_Experiments/00_Data/01_MainStudy/' + specimen + '_resample.csv'
    U_uFE, _ = hFEf.read_RFnodeFile(fe_path + fe + 'RFnode.txt')
    U_hFE, _ = hFEf.read_RFnodeFile(hfe_path + 'RFnode.txt')
    _, F_uFE = hFEf.read_RFnodeFile(fe_path + fe + 'RFnodeFix.txt')
    _, F_hFE = hFEf.read_RFnodeFile(hfe_path + 'RFnodeFix.txt')
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

    plt.plot(a_y, a_f-a_f[0], label='Experiment')
    plt.plot(U_uFE, -F_uFE, label='uFE')
    plt.plot(U_hFE, -F_hFE, label='hFE')
    # plt.title(specimen + ' (No. ' + str(no) + ')')
    plt.xlabel('Displacement / mm')
    plt.ylabel('Force / N')
    plt.legend()
    plt.savefig('/home/biomech/Documents/01_Icotec/02_FEA/91_Pictures/02_uFE/' + number + '_'
                + specimen + '_ForceDisplacement.png')
    plt.close()
    if no != 16:
        F_extrem[4, no] = F_hFE[73]
        F_extrem[5, no] = F_hFE[94]
        F_extrem[0, no] = hFEf.Peak_exp(3, no)
        F_extrem[1, no] = hFEf.Peak_exp(4, no)
        F_extrem[2, no] = F_uFE[45]
        F_extrem[3, no] = F_uFE[-1]

plt.figure()
plt.scatter(F_extrem[0, :], F_extrem[2, :], label='2 mm Amplitude')
plt.scatter(F_extrem[1, :], F_extrem[3, :], label='4 mm Amplitude')
plt.plot([0, 120], [0, 120], 'k--')
plt.xlabel('Force Experiment / N')
plt.ylabel('Force uFE / N')
plt.legend()
plt.savefig('/home/biomech/Documents/01_Icotec/02_FEA/91_Pictures/02_uFE/scatter_uFE_exp')

plt.figure()
plt.scatter(F_extrem[4, :], F_extrem[2, :], label='2 mm Amplitude')
plt.scatter(F_extrem[5, :], F_extrem[3, :], label='4 mm Amplitude')
plt.plot([0, 120], [0, 120], 'k--')
plt.xlabel('Force hFE / N')
plt.ylabel('Force uFE / N')
plt.legend()
plt.savefig('/home/biomech/Documents/01_Icotec/02_FEA/91_Pictures/02_uFE/scatter_hFE_uFE')

plt.figure()
plt.scatter(F_extrem[0, :], F_extrem[4, :], label='2 mm Amplitude')
plt.scatter(F_extrem[1, :], F_extrem[5, :], label='4 mm Amplitude')
plt.plot([0, 120], [0, 120], 'k--')
plt.xlabel('Force Experiment / N')
plt.ylabel('Force hFE / N')
plt.legend()
plt.savefig('/home/biomech/Documents/01_Icotec/02_FEA/91_Pictures/02_uFE/scatter_hFE_exp')

plt.figure()
plt.scatter(F_extrem[0, :], F_extrem[4, :], label='2 mm Amplitude hFE')
plt.scatter(F_extrem[1, :], F_extrem[5, :], label='4 mm Amplitude hFE')
plt.scatter(F_extrem[0, :], F_extrem[2, :], label='2 mm Amplitude uFE', color='#1f77b4', marker='x')
plt.scatter(F_extrem[1, :], F_extrem[3, :], label='4 mm Amplitude uFE', color='#ff7f0e', marker='x')
plt.plot([0, 120], [0, 120], 'k--')
plt.xlabel('Force Experiment / N')
plt.ylabel('Force FE / N')
plt.legend()
plt.savefig('/home/biomech/Documents/01_Icotec/02_FEA/91_Pictures/02_uFE/scatter_hFE_uFE_exp')

plt.close('all')