import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pickle

from hFE_functions import read_RFnodeFile, read_resample, Peak_exp, lin_reg, Peak_exp_d


def read_FE_(number, model_code, plot, fric_):
    # Locations
    model_code1 = []
    model_code2 = []
    specimens = open('/home/biomech/Documents/01_Icotec/Specimens.txt', 'r').read().splitlines()  # Read specimens
    loc_Exp = '/home/biomech/Documents/01_Icotec/01_Experiments/00_Data/01_MainStudy/'  # location experimental results
    loc_FEA = '/home/biomech/Documents/01_Icotec/02_FEA/01_MainStudy/'  # location of fea results
    if number in [0, 2, 5, 7, 8, 10, 13, 15, 16, 18, 21, 23, 24, 26, 29, 31, 32]:  # PEEK
        model_code1 = str(int(model_code[:2]) - 0) + model_code[2:19] + fric_.split('.')[-1] + '_P'
        model_code2 = str(int(model_code[:2]) - 0) + model_code[2:19] + fric_.split('.')[-1] + '_P'  # HERE -0 --> -2

    elif number in [1, 3, 4, 6, 9, 11, 12, 14, 17, 19, 20, 22, 25, 27, 28, 30, 33]:  # Ti
        fric_ = '0.2'  # HERE friction for Ti changed
        model_code1 = str(int(model_code[:2]) + 3) + model_code[2:19] + fric_.split('.')[-1] + '_T'  # HERE +1 --> -1
        model_code2 = str(int(model_code[:2]) + 3) + model_code[2:19] + fric_.split('.')[-1] + '_T'  # HERE -1 --> -3
    else:
        print('Invalid model code!')
    # print(model_code1)
    # print(model_code2)
    specimen = specimens[number]
    file = [loc_Exp + specimen + '_resample.csv',
            loc_FEA + specimen + '/' + model_code1[:14] + '/' + model_code1 + '_RFnode.txt',
            loc_FEA + specimen + '/' + model_code1[:14] + '/' + model_code1 + '_RFnodeFix.txt',
            loc_FEA + specimen + '/' + model_code2[:14] + '/' + model_code2 + '_RFnode.txt',
            loc_FEA + specimen + '/' + model_code2[:14] + '/' + model_code2 + '_RFnodeFix.txt'
            ]
    # Load data
    # [A_x, A_y, A_z, A_rx, A_ry, A_rz, a_y, a_f, a_c]
    RFy = []
    Uy = []
    [_, A_y, _, _, _, _, a_y, a_f, _] = read_resample(file[0])# load experimental result file (csv)
    # print(file[0])
    [Uy1, _] = read_RFnodeFile(file[1])  # read y displacement of moving reference node
    [_, RFy1] = read_RFnodeFile(file[2])  # read y reaction force of fixed reference node
    [Uy2, _] = read_RFnodeFile(file[3])  # read y displacement of moving reference node
    [_, RFy2] = read_RFnodeFile(file[4])  # read y reaction force of fixed reference node
    if len(Uy1) >= len(Uy2):
        Uy = Uy1
        RFy = RFy1
    elif len(Uy1) < len(Uy2):
        Uy = Uy2
        RFy = RFy2
    if plot == 1:
        fig1, ax1 = plt.subplots(1, 1, figsize=(9, 6))  # set figure size
        plt.title('Experimental results ' + specimen + ' ' + model_code1.split('_')[-1])
        if model_code1.split('_')[-1] == 'P':
            ax1.plot(Uy, -RFy + RFy[0], label='FEA PEEK', color='#1f77b4')  # plot fea results
        else:
            ax1.plot(Uy, -RFy + RFy[0], label='FEA Ti', color='#1f77b4')  # plot fea results
        if number in range(0, 40):  # [0, 1, 2, 10, 14, 17]:
            # print('Using Acumen displacement.')
            ax1.plot(a_y, a_f - a_f[0], label='Experiment', color='#ff7f0e')  # plot experimental results
        # else:
        #     ax1.plot(A_y, a_f - a_f[0], label='Experiment', color='#ff7f0e')  # plot experimental results
        # ax1.plot(a_y, a_f - a_f[0], '--', label='_nolegend_', color='#ff7f0e')  # plot y data from acumen

        ax1.legend()
        ax1.set_xlabel('Displacement / mm')
        ax1.set_ylabel('Force / N')

        fig1.savefig('/home/biomech/Documents/01_Icotec/02_FEA/91_Pictures/00_Exp_FE/Exp_FE_' +
                     specimen + '_' + model_code1 + '.png')
    elif plot == 2:
        fig2, ax2 = plt.subplots(1, 1, figsize=(9, 6))  # set figure size
        plt.title('Experimental results ' + specimen + ' ' + model_code1.split('_')[-1] + ' (no. ' + str(number) + ')')
        # if model_code.split('_')[-1] == 'P':
        #     ax2.plot(Uy1, -RFy1 + RFy1[0], label='FEA PEEK 88', color='#1f77b4')  # plot fea results
        #     ax2.plot(Uy2, -RFy2 + RFy2[0], '--', label='FEA PEEK 86', color='#1f77b4')  # plot fea results
        # else:
        #     ax2.plot(Uy1, -RFy1 + RFy1[0], label='FEA Ti 87', color='#1f77b4')  # plot fea results
        #     ax2.plot(Uy2, -RFy2 + RFy2[0], label='FEA Ti 85', color='#1f77b4')  # plot fea results
        # ax2.plot(a_y, a_f - a_f[0], label='Experiment Acumen', color='#ff7f0e')  # plot experimental (acumen)
        # ax2.plot(a_y, a_f, label='Experiment Acumen', color='#ff7f0e')  # plot experimental (acumen)
        # ax2.plot(A_y, a_f - a_f[0], label='Experiment Aramis', color='#ff7f0e')  # plot experimental (Aramis)
        ax2.axhline(y=0, color='k', linestyle='--', lw=0.5)
        ax2.axvline(x=0, color='k', linestyle='--', lw=0.5)
        ax2.legend()
        ax2.set_xlabel('Displacement / mm')
        ax2.set_ylabel('Force / N')
    return Uy, RFy, A_y, a_y, a_f


# %% Linear regression
fs = 13.5

# plt.close('all')
# PEEK, without 0 (diff ampl), 24 (Exp. weird)
peek_samples = [2, 5, 7, 8, 10, 13, 15, 16, 18, 21, 23, 26, 29, 31, 32]  # without 24
# Titanium, without 1 (diff ampl)
ti_samples = [3, 4, 6, 9, 11, 12, 14, 17, 19, 20, 22, 27, 28, 30, 33]  # without 25

x = 0  # 0 = 0.25 mm, 1 = 0.5 mm, 2 = 1 mm, 3 = 2 mm, 4 = 4 mm, 5 = 8 mm, 6 = 16 mm
lab = ['0.25 mm', '0.5 mm', '1 mm', '2 mm', '4 mm', '8 mm', '16 mm']
# lab = ['_nolegend_', '_nolegend_', '_nolegend_', '_nolegend_', '_nolegend_', '_nolegend_', '_nolegend_', '_nolegend_']
x0 = 0
x1 = 7  # max 7
# model = '88_L50_S50_D45_d1_02_P'  # automatically switches to titanium for respective samples
model = '60_L50_S50_D45_d1_05_P'

# peak_FE
RFy_FE = np.zeros((x1, 34))
RFy_exp = np.zeros((x1, 34))
# col = ['#1f77b4', '#1f77b4', '#1f77b4', '#1f77b4', '#1f77b4', '#1f77b4', '#1f77b4', '#1f77b4']
# col = ['#ff7f0e', '#ff7f0e', '#ff7f0e', '#ff7f0e', '#ff7f0e', '#ff7f0e', '#ff7f0e', '#ff7f0e']
col = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f']
fig, axs = plt.subplots(1, 1)
fig.set_figheight(7)
fig.set_figwidth(7)
RFy_FE_P = []
RFy_FE_T = []
RFy_exp_P = []
RFy_exp_T = []

RFy_exp_all = np.zeros((x1, 34))
loglog = 0
alp = 0.3
if loglog:
    F_range = np.array([-0.5, 2.6])
else:
    F_range = np.array([-10, 410])
plt.scatter(-1e9, -1e9, color='k', marker='v', label='icotec')
plt.scatter(-1e9, -1e9, color='k', marker='s', label='DPS')
# plt.scatter(-1e9, -1e9, color=col[0], marker='v', label='Mat A')
friction = '0.5'
for x in range(x0, x1):
    for i in range(2, 32):  # 2-34 because 0, 1 not existing in data frame
    # for i in peek_samples:
        # print('x: ' + str(x) + ' , i: ' + str(i))
        RFy_exp_all[x, i] = Peak_exp(x, i)
        try:
            [_, RFy_, _, _, _] = read_FE_(i, model, 0, friction)
        except FileNotFoundError:
            continue
        try:
            if loglog:
                if RFy_[x * 21 + 10] < 1:
                    RFy_FE[x, i] = 1
                else:
                    RFy_FE[x, i] = np.log10(RFy_[x * 21 + 10])
            else:
                RFy_FE[x, i] = RFy_[x * 21 + 10]
        except IndexError:
            continue
        if loglog:
            if Peak_exp(x, i) < 1:
                RFy_exp[x, i] = 1
            else:
                RFy_exp[x, i] = np.log10(Peak_exp(x, i))
        else:
            RFy_exp[x, i] = Peak_exp(x, i)
        if i == 29:  # Displacement labels
            plt.scatter(RFy_exp[x, i], RFy_FE[x, i], color=col[x], label=lab[x], marker='v')
        if i in peek_samples:  # P
            plt.scatter(RFy_exp[x, i], RFy_FE[x, i], color=col[x], label='_nolegend_', marker='v')
            RFy_FE_P.append(RFy_FE[x, i])
            RFy_exp_P.append(RFy_exp[x, i])
        elif i in ti_samples:
            plt.scatter(RFy_exp[x, i], RFy_FE[x, i], color=col[x], label='_nolegend_', marker='s')
            RFy_FE_T.append(RFy_FE[x, i])
            RFy_exp_T.append(RFy_exp[x, i])
        elif i == 24:  # HERE exclude sample from regression, plot transparent
            print('Excluded ' + str(i))
            # plt.scatter(RFy_exp[x, i], RFy_FE[x, i], color=col[x], label='_nolegend_', marker='v', alpha=alp)
        elif i == 25:  # HERE exclude sample from regression, plot transparent
            print('Excluded ' + str(i))
            # plt.scatter(RFy_exp[x, i], RFy_FE[x, i], color=col[x], label='_nolegend_', marker='s', alpha=alp)
        del RFy_
axs.plot(F_range, F_range, 'k', label='1:1')
if loglog:
    axs.set_xlabel('log$_{10}$(Experiment / N)', fontsize=fs)
    axs.set_ylabel('log$_{10}$(FE / N)', fontsize=fs)
else:
    axs.set_xlabel('Experiment / N', fontsize=fs)
    axs.set_ylabel('FE / N', fontsize=fs)
# axs.set_title('Peak Forces at ' + str(2 ** (x0 - 2)) + ' mm - ' + str(2 ** (x1 - 3)) + ' mm amplitudes')
# ', $\mu$ = ' + friction)
axs.set_aspect('equal')
axs.set_xlim(F_range)
axs.set_ylim(F_range)

specimen_names = open('/home/biomech/Documents/01_Icotec/Specimens.txt', 'r').read().splitlines()  # Read specimens
pd.DataFrame(RFy_FE).to_csv('/home/biomech/Downloads/corr.csv', index_label='Amplitude',
                            header=specimen_names)
# ['0.25 mm', '0.5 mm', '1 mm', '2 mm', '4 mm', '8 mm', '16 mm', ])
pd.DataFrame(RFy_exp_all).to_csv('/home/biomech/Downloads/corr2.csv', index_label='Amplitude',
                                 header=specimen_names)
print('done1')
regression_P, xx_P, yy_P = lin_reg(np.array(RFy_exp_P), np.array(RFy_FE_P))
axs.plot(F_range, F_range * regression_P.params[1] + regression_P.params[0], color='k', linestyle='dashdot',
         label='icotec:')
# axs.plot(F_range, F_range * regression_P.params[1] + regression_P.params[0], color='k', linestyle='dashdot',
#          label='0.5:')
if regression_P.pvalues[1] >= 0.05:
    lab_pvalue_P = 'p = ' + str(np.round(regression_P.pvalues[1], 2))
else:
    lab_pvalue_P = 'p < 0.05'
axs.plot([-1, 0], [-1, 0], color='w', linestyle='dashdot',
         label='R$^2$ = {:0.2f}'.format(np.round(regression_P.rsquared, 2)))
axs.plot([-1, 0], [-1, 0], color='w', label=lab_pvalue_P + '\n')

regression_T, xx_T, yy_T = lin_reg(np.array(RFy_exp_T), np.array(RFy_FE_T))
axs.plot(F_range, F_range * regression_T.params[1] + regression_T.params[0], color='k', linestyle='dotted',
         label='DPS:')
if regression_T.pvalues[1] >= 0.05:
    lab_pvalue_T = 'p = ' + str(np.round(regression_T.pvalues[1], 2))
else:
    lab_pvalue_T = 'p < 0.05'
axs.plot([-1, 0], [-1, 0], color='w', linestyle='dashed',
         label='R$^2$ = {:0.2f}'.format(np.round(regression_T.rsquared, 2)))
axs.plot([-1, 0], [-1, 0], color='w', label=lab_pvalue_T)

if loglog:
    plt.legend(framealpha=1, loc='upper left', fontsize=fs)
else:
    plt.legend(framealpha=1, loc='lower right', fontsize=fs)

# axs.set_ylabel('Difference / %', fontsize=fs)
# ax4.legend(fontsize=fs)
# axs.set_ylim([-5, 5])
axs.tick_params(axis='both', which='major', labelsize=fs)
axs.tick_params(axis='both', which='minor', labelsize=fs-2)
# plt.subplots_adjust(left=0.2)
# ax4.tick_params(
#      axis='x',           # changes apply to the x-axis
#      which='both',       # both major and minor ticks are affected
#      bottom=False,       # ticks along the bottom edge are off
#      top=False,          # ticks along the top edge are off
#      labelbottom=False)  # labels along the bottom edge are off
# if loglog:
#     plt.savefig('/home/biomech/Documents/GitHub/05_Report/03_Pictures_Res/hFE_regression_log.eps')
# else:
#     plt.savefig('/home/biomech/Documents/GitHub/05_Report/03_Pictures_Res/hFE_regression.eps')
# %% Stiffness hFE
fs = 13.5
plt.close('all')
peek_samples = [2, 5, 7, 8, 10, 13, 15, 16, 18, 21, 23, 26, 29, 31, 32]  # PEEK, without 0, 24
ti_samples = [3, 4, 6, 9, 11, 12, 14, 17, 19, 20, 22, 27, 28, 30, 33]  # Ti, without 0, 25
x = 0  # 0 = 0.25 mm, 1 = 0.5 mm, 2 = 1 mm, 3 = 2 mm, 4 = 4 mm, 5 = 8 mm, 6 = 16 mm
lab = ['0.25 mm', '0.5 mm', '1 mm', '2 mm', '4 mm', '8 mm', '16 mm']
x0 = 4
x1 = 7  # max 7
model = '60_L50_S50_D45_d1_05_P'
RFy_FE = np.zeros((x1, 34))
RFy_exp = np.zeros((x1, 34))
col = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f']
fig, axs = plt.subplots(1, 1)
fig.set_figheight(7)
fig.set_figwidth(7)
RFy_FE_P = []
RFy_FE_T = []
RFy_exp_P = []
RFy_exp_T = []

stiffness_FE = np.zeros((x1, 34))
stiffness_exp = np.zeros((x1, 34))
st_FE = []
st_exp = []

RFy_exp_all = np.zeros((x1, 34))
loglog = 0
alp = 0.3
if loglog:
    F_range = np.array([-0.5, 2.6])
else:
    F_range = np.array([-10, 410])
friction = '0.5'
for x in range(x0, x1):
    for i in range(2, 32):  # 2-34 because 0, 1 not existing in data frame
        try:
            [_, RFy_, _, _, _] = read_FE_(i, model, 0, friction)
        except FileNotFoundError:
            continue
        try:
            st_FE.append(RFy_[x * 21 + 10] / (2 ** (x - 2)))
            st_exp.append(Peak_exp(x, i) / (2 ** (x - 2)))
        except IndexError:
            continue

axs.plot(F_range, F_range, 'k', label='1:1')
axs.set_xlabel('Stiffness Experiment / N/mm', fontsize=fs)
axs.set_ylabel('Stiffness FE / N/mm', fontsize=fs)
axs.set_aspect('equal')
axs.set_xlim([0, 60])
axs.set_ylim([0, 60])
axs.scatter(st_exp, st_FE, color='k')
specimen_names = open('/home/biomech/Documents/01_Icotec/Specimens.txt', 'r').read().splitlines()  # Read specimens
pd.DataFrame(RFy_FE).to_csv('/home/biomech/Downloads/corr.csv', index_label='Amplitude',
                            header=specimen_names)
pd.DataFrame(RFy_exp_all).to_csv('/home/biomech/Downloads/corr2.csv', index_label='Amplitude',
                                 header=specimen_names)
regression_P, xx_P, yy_P = lin_reg(np.array(st_exp), np.array(st_FE))
if regression_P.params[0] < 0:
    lab_regr = str(round(regression_P.params[1], 2)) + 'x - ' + str(round(-regression_P.params[0], 2))
else:
    lab_regr = str(round(regression_P.params[1], 2)) + 'x + ' + str(round(regression_P.params[0], 2))
axs.plot(F_range, F_range * regression_P.params[1] + regression_P.params[0], color='k', linestyle='dashdot',
         label=lab_regr)
if regression_P.pvalues[1] >= 0.05:
    lab_pvalue_P ='p = ' + str(np.round(regression_P.pvalues[1], 2))
else:
    lab_pvalue_P ='R$^2$ = {:0.2f}'.format(np.round(regression_P.rsquared, 2)) +  ', p < 0.05'
axs.plot([-1, 0], [-1, 0], color='w', label=lab_pvalue_P + '\n')
plt.legend(framealpha=1, loc='lower right', fontsize=fs)
axs.tick_params(axis='both', which='major', labelsize=fs)
axs.tick_params(axis='both', which='minor', labelsize=fs-2)
low = str(2**(x0-2)).replace('.', '')
high = str(2**(x-2)).replace('.', '')
print(low)
print(high)
plt.savefig('/home/biomech/Documents/GitHub/05_Report/03_Pictures_Res/hFE_stiffness_' + low + 'mm_' + high + 'mm.eps')
# plt.close()

# %% Each amplitude
fig5, axs5 = plt.subplots(1, 1)
fig5.set_figheight(7)
fig5.set_figwidth(7)
with open('/home/biomech/Documents/01_Icotec/01_Experiments/03_Analysis/mergedDf.pkl', 'rb') as f:
    merged = pickle.load(f)

for x in range(x0, x1):
    for i in range(2, 34):
        peakF = merged['MaxForce'][(i - 2) * 7 + x]
        if i == 8 and x == 0:
            axs5.scatter(x - 0.3, peakF, color='r', marker='v', label='Experiment PEEK')
        elif i == 9 and x == 0:
            axs5.scatter(x + 0.1, peakF, color='r', marker='s', label='Experiment Ti')
        elif i == 24:  # HERE exclude sample
            axs5.scatter(x - 0.3, peakF, color='r', marker='v', alpha=alp, label='_nolegend_')
        else:
            if i in peek_samples:  # P (missing 0)
                axs5.scatter(x - 0.3, peakF, color='r', marker='v', label='_nolegend_')
            elif i in ti_samples:  # T (missing 1)
                axs5.scatter(x + 0.1, peakF, color='r', marker='s', label='_nolegend_')
        try:
            [_, RFy_, _, _, _] = read_FE_(i, model, 0, '0.2')
        except FileNotFoundError:
            continue
        try:
            RFy_FE = RFy_[x * 21 + 10]

        except IndexError:
            continue
        if i == 8 and x == 0:
            axs5.scatter(x - 0.1, RFy_FE, color='b', marker='v', label='FE PEEK')
        elif i == 9 and x == 0:
            axs5.scatter(x + 0.3, RFy_FE, color='b', marker='s', label='FE Ti')
        elif i == 25:  # HERE exclude sample
            axs5.scatter(x + 0.3, RFy_FE, color='b', marker='s', alpha=alp, label='_nolegend_')
        else:
            if i in peek_samples:  # P
                axs5.scatter(x - 0.1, RFy_FE, color='b', marker='v', label='_nolegend_')
            elif i in ti_samples:
                axs5.scatter(x + 0.3, RFy_FE, color='b', marker='s', label='_nolegend_')
    plt.plot([-0.5, -0.5], [0, 400], 'k--')
    plt.plot([x + 0.5, x + 0.5], [0, 400], 'k--')

plt.legend(framealpha=1)
plt.xlabel('Amplitude')
plt.ylabel('Max. Force / N')
plt.xticks(np.arange(0, 7), lab)
plt.title('Peak Forces')
# %% Friction comparison
fig7, axs7 = plt.subplots(1, 1)

with open('/home/biomech/Documents/01_Icotec/01_Experiments/03_Analysis/mergedDf.pkl', 'rb') as f:
    merged = pickle.load(f)

for x in range(x0, x1):
    for i in range(2, 34):
        peakF = merged['MaxForce'][(i - 2) * 7 + x]
        if i == 8 and x == 0:
            axs7.scatter(x - 0.3, peakF, color='r', marker='v', label='Experiment PEEK')
        elif i == 9 and x == 0:
            axs7.scatter(x + 0.1, peakF, color='r', marker='s', label='Experiment Ti')
        elif i == 24:
            axs7.scatter(x - 0.3, peakF, color='r', marker='v', alpha=alp, label='_nolegend_')
        else:
            if i in peek_samples:  # P (missing 0)
                axs7.scatter(x - 0.3, peakF, color='r', marker='v', label='_nolegend_')
            elif i in ti_samples:  # T (missing 1)
                axs7.scatter(x + 0.1, peakF, color='r', marker='s', label='_nolegend_')
        try:
            [_, RFy_, _, _, _] = read_FE_(i, model, 0, '0.5')
        except FileNotFoundError:
            continue
        try:
            [_, RFy_2, _, _, _] = read_FE_(i, model, 0, '0.2')
        except FileNotFoundError:
            continue
        try:
            RFy_FE = RFy_[x * 21 + 10]
        except IndexError:
            continue
        try:
            RFy_FE2 = RFy_2[x * 21 + 10]
        except IndexError:
            continue
        if i == 8 and x == 0:
            axs7.scatter(x - 0.1, RFy_FE, color='b', marker='v', label='FE PEEK $\mu$ = 0.5')
            axs7.scatter(x - 0.0, RFy_FE2, color='k', marker='v', label='FE PEEK $\mu$ = 0.2')
        elif i == 9 and x == 0:
            axs7.scatter(x + 0.3, RFy_FE, color='b', marker='s', label='FE Ti $\mu$ = 0.5')
            axs7.scatter(x + 0.4, RFy_FE2, color='k', marker='s', label='FE Ti $\mu$ = 0.2')
        elif i == 25:
            axs7.scatter(x + 0.3, RFy_FE, color='b', marker='s', alpha=alp, label='_nolegend_')
            axs7.scatter(x + 0.4, RFy_FE2, color='k', marker='s', alpha=alp, label='_nolegend_')
        else:
            if i in peek_samples:  # P
                axs7.scatter(x - 0.1, RFy_FE, color='b', marker='v', label='_nolegend_')
                axs7.scatter(x - 0.0, RFy_FE2, color='k', marker='v', label='_nolegend_')
            elif i in ti_samples:
                axs7.scatter(x + 0.3, RFy_FE, color='b', marker='s', label='_nolegend_')
                axs7.scatter(x + 0.4, RFy_FE2, color='k', marker='s', label='_nolegend_')
    plt.plot([-0.5, -0.5], [0, 400], 'k--')
    plt.plot([x + 0.5, x + 0.5], [0, 400], 'k--')

plt.legend(framealpha=1)
plt.xlabel('Amplitude')
plt.ylabel('Max. Force / N')
plt.xticks(np.arange(0, 7), lab)
plt.title('Peak Forces')

# %% Residual Displacement
fig8, axs8 = plt.subplots(1, 1)
fig8.set_figheight(7)
fig8.set_figwidth(7)
with open('/home/biomech/Documents/01_Icotec/01_Experiments/03_Analysis/mergedDf.pkl', 'rb') as f:
    merged = pickle.load(f)
for x in range(x0, 1):  # x1):
    for i in [7]:
        # print('x: ' + str(x) + ', i: ' + str(i))
        # resExp = merged['ResidualDisplacement'][(i - 2) * 7 + x]
        try:
            [UY, Fy, _, _, _] = read_FE_(i, model, 0, '0.5')
        except FileNotFoundError:
            # print('no File')
            continue
        # try:
        # Uy_res = UY[x * 21 + 19]
        # except IndexError:
        # print('no Datapoints')
        # continue
        axs8.plot(UY, Fy)
        axs8.plot([-20, 2], [0, 0], 'k--')
        # axs8.scatter(resExp, resFE)
        # print(resFE)

# %% Linear regression incl friction


peek_samples = [2, 5, 7, 8, 10, 13, 15, 16, 18, 21, 23, 26, 29, 31, 32]
# Titanium, without 1 (diff ampl)
ti_samples = [3, 4, 6, 9, 11, 12, 14, 17, 19, 20, 22, 25, 27, 28, 30, 33]

x = 0  # 0 = 0.25 mm, 1 = 0.5 mm, 2 = 1 mm, 3 = 2 mm, 4 = 4 mm, 5 = 8 mm, 6 = 16 mm
lab = ['0.25 mm', '0.5 mm', '1 mm', '2 mm', '4 mm', '8 mm', '16 mm']
x0 = 0
x1 = 7  # max 7
model = '86_L50_S50_D45_d1_05_P'  # automatically switches to titanium for respective samples

# peak_FE
RFy_FE = np.zeros((x1, 34))
RFy_exp = np.zeros((x1, 34))
# col = ['k', 'k', 'k', 'k', 'k', 'k', 'k', 'k']
col = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f']
fig, axs = plt.subplots(1, 1)
fig.set_figheight(7)
fig.set_figwidth(7)
RFy_FE_P = []
RFy_FE_T = []
RFy_exp_P = []
RFy_exp_T = []

RFy_exp_all = np.zeros((x1, 34))
loglog = 0
alp = 0.3
if loglog:
    F_range = np.array([0, 2.6])
else:
    F_range = np.array([-10, 450])
plt.scatter(-1e9, -1e9, color='k', marker='v', label='PEEK, $\mu$ = 0.2')
plt.scatter(-1e9, -1e9, color='k', marker='^', label='PEEK, $\mu$ = 0.5')
plt.scatter(-1e9, -1e9, color='k', marker='s', label='Titanium, $\mu$ = 0.2')
plt.scatter(-1e9, -1e9, color='k', marker='o', label='Titanium, $\mu$ = 0.5')
mark = ['v', 's', '^', 'o']
friction = ['0.2', '0.5']
# friction = ['0.5']
for fric in range(len(friction)):
    for x in range(x0, x1):
        for i in range(2, 34):  # [2, 3, 4, 5, 10, 11]:  # 2-34 because 0, 1 not existing in data frame
            # print('x: ' + str(x) + ' , i: ' + str(i))
            RFy_exp_all[x, i] = Peak_exp(x, i)
            try:
                [_, RFy_, _, _, _] = read_FE_(i, model, 0, friction[fric])
            except FileNotFoundError:
                continue
            try:
                if loglog:
                    if RFy_[x * 21 + 10] < 1:
                        RFy_FE[x, i] = 1
                    else:
                        RFy_FE[x, i] = np.log10(RFy_[x * 21 + 10])
                else:
                    RFy_FE[x, i] = RFy_[x * 21 + 10]
            except IndexError:
                continue
            # print('Specimen: ' + str(i) + ', amplitude: ' + str(x) + ', Force FE: ' + str(RFy_FE))
            print(mark[2 * fric])
            print(mark[2 * fric + 1])
            if loglog:
                if Peak_exp(x, i) < 1:
                    RFy_exp[x, i] = 1
                else:
                    RFy_exp[x, i] = np.log10(Peak_exp(x, i))
            else:
                RFy_exp[x, i] = Peak_exp(x, i)
            if i == 8:
                plt.scatter(RFy_exp[x, i], RFy_FE[x, i], color=col[x], label=lab[x], marker=mark[2 * fric])
            if i in peek_samples:  # P
                plt.scatter(RFy_exp[x, i], RFy_FE[x, i], color=col[x], label='_nolegend_', marker=mark[2 * fric])
                RFy_FE_P.append(RFy_FE[x, i])
                RFy_exp_P.append(RFy_exp[x, i])
            elif i in ti_samples:
                plt.scatter(RFy_exp[x, i], RFy_FE[x, i], color=col[x], label='_nolegend_', marker=mark[2 * fric + 1])
                RFy_FE_T.append(RFy_FE[x, i])
                RFy_exp_T.append(RFy_exp[x, i])
            elif i == 24:  # HERE exclude sample
                plt.scatter(RFy_exp[x, i], RFy_FE[x, i], color=col[x], label='_nolegend_', marker=mark[2 * fric],
                            alpha=alp)
            # elif i == 25:  # HERE exclude sample
            #     plt.scatter(RFy_exp[x, i], RFy_FE[x, i], color=col[x], label='_nolegend_', marker='s', alpha=alp)
            del RFy_
axs.plot(F_range, F_range, 'k', label='1:1')
if loglog:
    axs.set_xlabel('log$_{10}$(Experiment / N)')
    axs.set_ylabel('log$_{10}$(FE / N)')
else:
    axs.set_xlabel('Experiment / N')
    axs.set_ylabel('FE / N')
axs.set_title('Peak Forces at ' + str(2 ** (x0 - 2)) + ' mm - ' + str(2 ** (x1 - 3)) + ' mm amplitudes')
axs.set_aspect('equal')
axs.set_xlim(F_range)
axs.set_ylim(F_range)

specimen_names = open('/home/biomech/Documents/01_Icotec/Specimens.txt', 'r').read().splitlines()  # Read specimens
pd.DataFrame(RFy_FE).to_csv('/home/biomech/Downloads/corr.csv', index_label='Amplitude',
                            header=specimen_names)
# ['0.25 mm', '0.5 mm', '1 mm', '2 mm', '4 mm', '8 mm', '16 mm', ])
pd.DataFrame(RFy_exp_all).to_csv('/home/biomech/Downloads/corr2.csv', index_label='Amplitude',
                                 header=specimen_names)
print('done1')
regression_P, xx_P, yy_P = lin_reg(np.array(RFy_exp_P), np.array(RFy_FE_P))
axs.plot(F_range, F_range * regression_P.params[1] + regression_P.params[0], color='k', linestyle='dashdot',
         label='PEEK:')
if regression_P.pvalues[1] >= 0.05:
    lab_pvalue_P = 'p = ' + str(np.round(regression_P.pvalues[1], 2))
else:
    lab_pvalue_P = 'p < 0.05'
axs.plot([-1, 0], [-1, 0], color='w', linestyle='dashdot',
         label='R$^2$ = {:0.2f}'.format(np.round(regression_P.rsquared, 2)))
axs.plot([-1, 0], [-1, 0], color='w', label=lab_pvalue_P + '\n')

regression_T, xx_T, yy_T = lin_reg(np.array(RFy_exp_T), np.array(RFy_FE_T))
axs.plot(F_range, F_range * regression_T.params[1] + regression_T.params[0], color='k', linestyle='dotted',
         label='Titanium:')
if regression_T.pvalues[1] >= 0.05:
    lab_pvalue_T = 'p = ' + str(np.round(regression_T.pvalues[1], 2))
else:
    lab_pvalue_T = 'p < 0.05'
axs.plot([-1, 0], [-1, 0], color='w', linestyle='dashed',
         label='R$^2$ = {:0.2f}'.format(np.round(regression_T.rsquared, 2)))
axs.plot([-1, 0], [-1, 0], color='w', label=lab_pvalue_T)

plt.legend(framealpha=1)
