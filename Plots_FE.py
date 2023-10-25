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


def read_icotec_experiment(sheet_, mat_):
    if mat_ == 'peek':
        file_ = '/home/biomech/Documents/01_Icotec/01_Experiments/99_Others/Screw_Bending_test/' \
                'F1798Cantilever_UniBern.xls'
    elif mat_ == 'ti':
        file_ = '/home/biomech/Documents/01_Icotec/01_Experiments/99_Others/Screw_Bending_DPS/' \
                'F1798Cantilever_UniBern_DPS.xls'
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
    data[tests[i]] = read_icotec_experiment(tests[i], 'peek')
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

#%% Ti material optimisation
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
def read_RFnodeFile(file_):
    # read data from text file
    df_ = np.loadtxt(file_, delimiter=',')
    uy_ = np.array(df_[:, 2])
    fy_ = np.array(df_[:, 1])
    return uy_, fy_


#  0: YM = 100'000, no plastic
#  1: YM = 110'000, no plastic
#  2: YM = 113'000, no plastic
#  3: YM = 113'000, pl (880/0)


plt.figure()
for i in [0, 2, 3]:
    uy, fy = read_RFnodeFile('/home/biomech/Documents/01_Icotec/02_FEA/00_Model/99_screw_DPS_Bending_job_RFnode_'
                             + str(i) + '.txt')
    plt.plot(-uy, -fy)
data = read_icotec_experiment('1.3', 'ti')
plt.plot(data['Dehnung'], data['Standardkraft'], color='k')
# %% uFE result
sample_list = open('/home/biomech/Documents/01_Icotec/Specimens.txt', 'r').read().splitlines()
# plt.close('all')
saving = 0
fs = 13.5
numbers = ['17']
Fe0_1 = []
Fe0 = []
Fe1 = []
Fe2_1 = []
Fe2 = []
Fe3 = []
Fe4_1 = []
Fe4 = []
Fe5 = []
# numbers = ['29', '30', '31', '32', '33', '34', '35']  # alpha all
# damp = ['0.01', '0.03', '0.05', '0.07', '0.1', '0.15', '0.2']  # alpha all
# numbers = ['17', '29', '33', '35']  # alpha some
# damp = ['No damping', '0.01', '0.1', '0.2']  # alpha some
# numbers = ['27', '28', '20', '24', '25', '19', '21', '22', ]  # beta all
# damp = ['0.0001', '0.0005', '0.005', '0.01', '0.02', '0.05', '0.1', '1', '10']  # beta all
# numbers = ['27', '24', '21', '22', '23']  # beta some
# damp = ['0.0001', '0.01', '0.1', '1', '10']  # beta some
# samples = [3, 4, 14, 17, 19, 20]
samples = [5, 7, 8, 10, 15, 18, 21, 26, 29, 31]  # 13 no uFE, 16, 23 no hfe. 24 exp bad. 2 and 32 still computing
# samples = [10]  # for damping
# samples = [21]  # for entire cycle
F_extrem = np.zeros((6, 32))
exp_hfe_plots = 0
plt.figure()
for n in range(len(numbers)):
    number = numbers[n]
    for no in samples:
        specimen = sample_list[no]
        print(specimen)
        fe_path = '/home/biomech/DATA/01_Icotec/02_FEA/02_uFE/01_MainStudy/' + \
                  specimen + '/'
        fe = number + '_' + specimen + '_0121399_'
        if no in [2, 5, 7, 8, 10, 13, 15, 16, 18, 21, 23, 24, 26, 29, 31, 32]:  # without 0
            hfe_path = '/home/biomech/Documents/01_Icotec/02_FEA/01_MainStudy/' + \
                      specimen + '/60_L50_S50_D45/60_L50_S50_D45_d1_05_P_'
        elif no in [3, 4, 6, 9, 11, 12, 14, 17, 19, 20, 22, 25, 27, 28, 30, 33]:  # without 1
            hfe_path = '/home/biomech/Documents/01_Icotec/02_FEA/01_MainStudy/' + \
                       specimen + '/63_L50_S50_D45/63_L50_S50_D45_d1_02_T_'

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
            plt.plot(t_ke, ke/ie*100)
            plt.plot([t_ke[0], t_ke[-1]], [10, 10])
            plt.xlabel('Time / s')
            plt.ylabel('Energy ratio / %')
            plt.ylim([0, 15])
            plt.figure()
        if no in [2, 5, 7, 8, 10, 13, 15, 16, 18, 21, 23, 24, 26, 29, 31, 32]:  # without 0
            # plt.title(number + '_' + specimen + ' (PEEK)')
            print('PEEK')
            scr = 'icotec'
        elif no in [3, 4, 6, 9, 11, 12, 14, 17, 19, 20, 22, 25, 27, 28, 30, 33]:  # without 1
            # plt.title(number + '_' + specimen + ' (Ti)')
            print('Ti')
            scr = 'DPS'
        if not exp_hfe_plots:
            plt.plot(a_y, a_f-a_f[0], label='Experiment (' + scr + ')')
            plt.plot(U_hFE, -F_hFE, label='hFE')
        if len(numbers) > 1:
            plt.plot(U_uFE, -F_uFE, label=damp[n])
            plt.xlim([-4.5, 0])
            plt.ylim([-150, 50])
        else:
            plt.plot(U_uFE, -F_uFE, label='uFE, fast')
        # plt.title(specimen + ' (No. ' + str(no) + ')')
        plt.xlabel('Displacement / mm', fontsize=fs)
        plt.ylabel('Force / N', fontsize=fs)
        plt.legend(fontsize=fs)
        plt.tick_params(axis='both', which='major', labelsize=fs)
        plt.tick_params(axis='both', which='minor', labelsize=fs - 2)
        plt.subplots_adjust(left=0.15)
        # exp_hfe_plots = 1
        if saving:
            plt.savefig('/home/biomech/Documents/GitHub/05_Report/03_Pictures_Res/' + number + '_'
                        + specimen + '_ForceDisplacement.eps')
        if len(samples) > 1:
            plt.close()
        if no != 16:
            Fe0_1.append(hFEf.Peak_exp(2, no))  # Experiment 1 mm
            F_extrem[0, no] = hFEf.Peak_exp(3, no)  # Experiment 2 mm
            Fe0.append(hFEf.Peak_exp(3, no))
            F_extrem[1, no] = hFEf.Peak_exp(4, no)  # Experiment 4 mm
            Fe1.append(hFEf.Peak_exp(4, no))

            Fe2_1.append((F_uFE[32] + F_uFE[33]) / 2)  # uFE 1 mm, average (no datapoint at 1 mm)
            F_extrem[2, no] = F_uFE[45]  # uFE 2 mm
            Fe2.append(F_uFE[45])
            F_extrem[3, no] = F_uFE[90]  # uFE 4 mm
            Fe3.append(F_uFE[90])
            try:
                Fe4_1.append(F_hFE[52])  # hFE 1 mm
                F_extrem[4, no] = F_hFE[73]  # hFE 2 mm
                Fe4.append(F_hFE[73])
                F_extrem[5, no] = F_hFE[94]  # hFE 4 mm
                Fe5.append(F_hFE[94])
            except:
                print('no hFE data')
# %%
F_range = np.array([0, 140])
plt.figure()
# plt.scatter(F_extrem[0, :], F_extrem[2, :], label='2 mm Amplitude')
# plt.scatter(F_extrem[1, :], F_extrem[3, :], label='4 mm Amplitude')
plt.scatter(Fe0_1, Fe2_1, label='1 mm Amplitude')
plt.scatter(Fe0, Fe2, label='2 mm Amplitude')
# plt.scatter(Fe1, Fe3, label='4 mm Amplitude')
plt.plot(F_range, F_range, 'k')
plt.xlabel('Force Experiment / N', fontsize=fs)
plt.ylabel('Force uFE / N', fontsize=fs)
plt.xlim(F_range)
plt.ylim(F_range)
plt.tick_params(axis='both', which='major', labelsize=fs)
plt.tick_params(axis='both', which='minor', labelsize=fs-2)
# regression_P, xx_P, yy_P = hFEf.lin_reg(np.array(np.append(np.append(Fe0_1, Fe0), Fe1)),
#                                         np.array(np.append(np.append(Fe2_1, Fe2), Fe3)))
regression_P, xx_P, yy_P = hFEf.lin_reg(np.array(np.append(Fe0_1, Fe0)),
                                        np.array(np.append(Fe2_1, Fe2)))
if regression_P.params[0] < 0:
    lab_regr = str(round(regression_P.params[1], 2)) + 'x - ' + str(round(-regression_P.params[0], 2))
else:
    lab_regr = str(round(regression_P.params[1], 2)) + 'x + ' + str(round(regression_P.params[0], 2))
plt.plot(F_range, F_range * regression_P.params[1] + regression_P.params[0], color='k', linestyle='dashdot',
         label=lab_regr)
if regression_P.pvalues[1] >= 0.05:
    lab_pvalue_P = 'p = ' + str(np.round(regression_P.pvalues[1], 2))
else:
    lab_pvalue_P = 'p < 0.05'
plt.plot([-1, 0], [-1, 0], color='w', label='R$^2$ = {:0.2f}'.format(np.round(regression_P.rsquared, 2)) +
                                            ', ' + lab_pvalue_P)
plt.legend(fontsize=fs)
if saving:
    if no in [2, 5, 7, 8, 10, 13, 15, 16, 18, 21, 23, 24, 26, 29, 31, 32]:  # without 0
        plt.savefig('/home/biomech/Documents/GitHub/05_Report/03_Pictures_Res/' + number + '_scatter_uFE_exp_PEEK12.eps')
    elif no in [3, 4, 6, 9, 11, 12, 14, 17, 19, 20, 22, 25, 27, 28, 30, 33]:  # without 1
        plt.savefig('/home/biomech/Documents/GitHub/05_Report/03_Pictures_Res/' + number + '_scatter_uFE_exp_Ti12.eps')

plt.figure()
plt.scatter(Fe4_1, Fe2_1, label='1 mm Amplitude')
plt.scatter(Fe4, Fe2, label='2 mm Amplitude')
# plt.scatter(Fe5, Fe3, label='4 mm Amplitude')
plt.plot(F_range, F_range, 'k')
plt.xlabel('Force hFE / N', fontsize=fs)
plt.ylabel('Force uFE / N', fontsize=fs)
plt.xlim(F_range)
plt.ylim(F_range)
plt.tick_params(axis='both', which='major', labelsize=fs)
plt.tick_params(axis='both', which='minor', labelsize=fs-2)
# regression_P, xx_P, yy_P = hFEf.lin_reg(np.array(np.append(np.append(Fe4_1, Fe4), Fe5)),
#                                         np.array(np.append(np.append(Fe2_1, Fe2), Fe3)))
regression_P, xx_P, yy_P = hFEf.lin_reg(np.array(np.append(Fe4_1, Fe4)),
                                        np.array(np.append(Fe2_1, Fe2)))
if regression_P.params[0] < 0:
    lab_regr = str(round(regression_P.params[1], 2)) + 'x - ' + str(round(-regression_P.params[0], 2))
else:
    lab_regr = str(round(regression_P.params[1], 2)) + 'x + ' + str(round(regression_P.params[0], 2))
plt.plot(F_range, F_range * regression_P.params[1] + regression_P.params[0], color='k', linestyle='dashdot',
         label=lab_regr)
if regression_P.pvalues[1] >= 0.05:
    lab_pvalue_P = 'p = ' + str(np.round(regression_P.pvalues[1], 2))
else:
    lab_pvalue_P = 'p < 0.05'
plt.plot([-1, 0], [-1, 0], color='w', label='R$^2$ = {:0.2f}'.format(np.round(regression_P.rsquared, 2)) +
                                            ', ' + lab_pvalue_P)
plt.legend(fontsize=fs)
if saving:
    if no in [2, 5, 7, 8, 10, 13, 15, 16, 18, 21, 23, 24, 26, 29, 31, 32]:  # without 0
        plt.savefig('/home/biomech/Documents/GitHub/05_Report/03_Pictures_Res/' + number + '_scatter_hFE_uFE_PEEK12.eps')
    elif no in [3, 4, 6, 9, 11, 12, 14, 17, 19, 20, 22, 25, 27, 28, 30, 33]:  # without 1
        plt.savefig('/home/biomech/Documents/GitHub/05_Report/03_Pictures_Res/' + number + '_scatter_hFE_uFE_Ti12.eps')

plt.figure()
plt.scatter(Fe0_1, Fe4_1, label='1 mm Amplitude')
plt.scatter(Fe0, Fe4, label='2 mm Amplitude')
# plt.scatter(Fe1, Fe5, label='4 mm Amplitude')
plt.plot(F_range, F_range, 'k')
plt.xlabel('Force Experiment / N', fontsize=fs)
plt.ylabel('Force hFE / N', fontsize=fs)
plt.xlim(F_range)
plt.ylim(F_range)
plt.tick_params(axis='both', which='major', labelsize=fs)
plt.tick_params(axis='both', which='minor', labelsize=fs-2)
# regression_P, xx_P, yy_P = hFEf.lin_reg(np.array(np.append(np.append(Fe0_1, Fe0), Fe1)),
#                                         np.array(np.append(np.append(Fe4_1, Fe4), Fe5)))
regression_P, xx_P, yy_P = hFEf.lin_reg(np.array(np.append(Fe0_1, Fe0)),
                                        np.array(np.append(Fe4_1, Fe4)))
if regression_P.params[0] < 0:
    lab_regr = str(round(regression_P.params[1], 2)) + 'x - ' + str(round(-regression_P.params[0], 2))
else:
    lab_regr = str(round(regression_P.params[1], 2)) + 'x + ' + str(round(regression_P.params[0], 2))
plt.plot(F_range, F_range * regression_P.params[1] + regression_P.params[0], color='k', linestyle='dashdot',
         label=lab_regr)
if regression_P.pvalues[1] >= 0.05:
    lab_pvalue_P = 'p = ' + str(np.round(regression_P.pvalues[1], 2))
else:
    lab_pvalue_P = 'p < 0.05'
plt.plot([-1, 0], [-1, 0], color='w', label='R$^2$ = {:0.2f}'.format(np.round(regression_P.rsquared, 2)) +
                                            ', ' + lab_pvalue_P)
plt.legend(fontsize=fs)
if saving:
    if no in [2, 5, 7, 8, 10, 13, 15, 16, 18, 21, 23, 24, 26, 29, 31, 32]:  # without 0
        plt.savefig('/home/biomech/Documents/GitHub/05_Report/03_Pictures_Res/' + number + '_scatter_hFE_exp_PEEK124.eps')
    elif no in [3, 4, 6, 9, 11, 12, 14, 17, 19, 20, 22, 25, 27, 28, 30, 33]:  # without 1
        plt.savefig('/home/biomech/Documents/GitHub/05_Report/03_Pictures_Res/' + number + '_scatter_hFE_exp_Ti124.eps')

plt.figure()
plt.scatter(Fe0_1, Fe4_1, label='1 mm Amplitude')
plt.scatter(Fe0, Fe4, label='2 mm Amplitude')
# plt.scatter(Fe1, Fe5, label='4 mm Amplitude')
plt.scatter(Fe0_1, Fe2_1, label='_nolegend_', color='#1f77b4', marker='x')
plt.scatter(Fe0, Fe2, label='_nolegend_', color='#ff7f0e', marker='x')
# plt.scatter(Fe1, Fe3, label='_nolegend_', color='#2ca02c', marker='x')
plt.scatter(-1e9, -1e9, label='hFE', color='k')
plt.scatter(-1e9, -1e9, label='uFE', color='k', marker='x')
plt.plot(F_range, F_range, 'k')
plt.xlabel('Force Experiment / N', fontsize=fs)
plt.ylabel('Force FE / N', fontsize=fs)
plt.xlim(F_range)
plt.ylim(F_range)
plt.tick_params(axis='both', which='major', labelsize=fs)
plt.tick_params(axis='both', which='minor', labelsize=fs-2)
plt.legend(fontsize=fs)
if saving:
    if no in [2, 5, 7, 8, 10, 13, 15, 16, 18, 21, 23, 24, 26, 29, 31, 32]:  # without 0
        plt.savefig('/home/biomech/Documents/GitHub/05_Report/03_Pictures_Res/' + number + '_scatter_hFE_uFE_exp_PEEK12.eps')
    elif no in [3, 4, 6, 9, 11, 12, 14, 17, 19, 20, 22, 25, 27, 28, 30, 33]:  # without 1
        plt.savefig('/home/biomech/Documents/GitHub/05_Report/03_Pictures_Res/' + number + '_scatter_hFE_uFE_exp_Ti12.eps')

#plt.close('all')

# %% Amplitude pattern

plt.figure()
ampl = np.array([0, 0.25, 0, 0.5, 0, 1, 0, 2, 0, 4, 0, 8, 0, 16, 0])
plt.plot(-ampl)
fs = 13.5
plt.xlabel('Time', fontsize=fs)
plt.ylabel('Displacement / mm', fontsize=fs)
plt.tick_params(axis='both', which='major', labelsize=fs)
plt.tick_params(axis='both', which='minor', labelsize=fs-2)
plt.yticks([-16, -8, -4, -2, -1, 0])
plt.xticks([])