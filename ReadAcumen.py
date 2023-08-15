import numpy as np
import matplotlib.pyplot as plt
import pickle
import statsmodels.api as sm
import pandas as pd


# %%


def read_RFnodeFile(file_):
    # read data from text file
    df_ = np.loadtxt(file_, delimiter=',')

    # Contains the frame number:
    # frame = np.array(df_[:, 0])

    # reaction moment of the implant:
    # rmx_ = np.array(df_[:, 1])
    # rmy_ = np.array(df_[:, 2])
    # rmz_ = np.array(df_[:, 3])

    # reaction force of the implant:
    # rfx_ = np.array(df_[:, 4])
    rfy_ = np.array(df_[:, 5])
    # rfz_ = np.array(df_[:, 6])

    # transverse disp. of the implant:
    # ux_ = np.array(df_[:, 7])
    uy_ = np.array(df_[:, 8])
    # uz_ = np.array(df_[:, 9])

    # rotational disp. of the implant:
    # urx_ = np.array(df_[:, 10])
    # ury_ = np.array(df_[:, 11])
    # urz_ = np.array(df_[:, 12])

    return uy_, rfy_


def read_acumen(file_a):
    df_ = pd.read_csv(file_a, delimiter=';', skiprows=[0, 2])
    t_ = pd.DataFrame(df_, columns=['Time ']).to_numpy()
    d_ = pd.DataFrame(df_, columns=['Axial Displacement ']).to_numpy()
    # d_ = d_ - d_[0]  # calibrate displacement to zero at beginning
    f_ = pd.DataFrame(df_, columns=['Axial Force ']).to_numpy()
    # f_set_ = pd.DataFrame(df_, columns=['Axial Force Command ']).to_numpy()
    cycle_ = pd.DataFrame(df_, columns=['Axial Count ']).to_numpy()
    # arr_ = 0
    peak_ = np.zeros(int(np.max(cycle_)))
    vall_ = np.zeros(int(np.max(cycle_)))
    for j_ in range(2, int(np.max(cycle_))):
        # del arr_
        arr_ = np.where((cycle_ == j_) | (cycle_ == j_ + .5))[0]
        peak_[j_] = arr_[int(np.argmin(f_[arr_]))]
        vall_[j_] = arr_[int(np.argmax(f_[arr_]))]

    peak_ = peak_.astype(int)
    vall_ = vall_.astype(int)

    return cycle_, d_, f_, peak_, vall_, t_


def read_ARAMIS(file_A):
    df_ = pd.read_csv(file_A, delimiter=';', skiprows=[0])
    stop = np.where(df_.isnull())[0][0]
    lx = pd.DataFrame(df_, columns=['RodCsys→BoneCsys.LX [mm]'])
    ly = pd.DataFrame(df_, columns=['RodCsys→BoneCsys.LY [mm]'])
    lz = pd.DataFrame(df_, columns=['RodCsys→BoneCsys.LZ [mm]'])
    phiX = pd.DataFrame(df_, columns=['RodCsys→BoneCsys.Phi(X) [°]'])
    thetaY = pd.DataFrame(df_, columns=['RodCsys→BoneCsys.Theta(Y) [°]'])
    psiZ = pd.DataFrame(df_, columns=['RodCsys→BoneCsys.Psi(Z) [°]'])
    t_ = pd.DataFrame(df_, columns=['Time UTC'])

    return lx[:stop], ly[:stop], lz[:stop], phiX[:stop], thetaY[:stop], psiZ[:stop], t_[:stop]


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx


def find_first(array, value):
    array = np.asarray(array)
    idx = next(xd for xd, val in enumerate(array)
               if val <= value)
    return idx


def read_resample(file_r):
    df_ = pd.read_csv(file_r, delimiter=',')
    A_x_ = pd.DataFrame(df_, columns=['Aramis X']).to_numpy()
    A_y_ = pd.DataFrame(df_, columns=['Aramis Y']).to_numpy()
    A_z_ = pd.DataFrame(df_, columns=['Aramis Z']).to_numpy()
    A_rx_ = pd.DataFrame(df_, columns=['Aramis rX']).to_numpy()
    A_ry_ = pd.DataFrame(df_, columns=['Aramis rY']).to_numpy()
    A_rz_ = pd.DataFrame(df_, columns=['Aramis rZ']).to_numpy()
    a_y_ = pd.DataFrame(df_, columns=['Acumen Y']).to_numpy()
    a_f_ = pd.DataFrame(df_, columns=['Acumen Fy']).to_numpy()
    a_c_ = pd.DataFrame(df_, columns=['Acumen C']).to_numpy()

    return A_x_, A_y_, A_z_, A_rx_, A_ry_, A_rz_, a_y_, a_f_, a_c_


def read_exp_peaks():
    with open('/home/biomech/Documents/01_Icotec/01_Experiments/03_Analysis/mergedDf.pkl', 'rb') as f_:
        data = pickle.load(f_)
    return data


def read_FE_(number, model_code, plot, fric_):
    # Locations
    specimens = open('/home/biomech/Documents/01_Icotec/Specimens.txt', 'r').read().splitlines()  # Read specimens
    loc_Exp = '/home/biomech/Documents/01_Icotec/01_Experiments/00_Data/01_MainStudy/'  # location experimental results
    loc_FEA = '/home/biomech/Documents/01_Icotec/02_FEA/01_MainStudy/'  # location of fea results
    if number in [0, 2, 5, 7, 8, 10, 13, 15, 16, 18, 21, 23, 24, 26, 29, 31, 32]:
        model_code = model_code[:19] + fric_.split('.')[-1] + '_P'
    elif number in [1, 3, 4, 6, 9, 11, 12, 14, 17, 19, 20, 22, 25, 27, 28, 30, 33]:
        model_code = '87' + model_code[2:19] + fric_.split('.')[-1] + '_T'
    else:
        print('Invalid model code!')
    specimen = specimens[number]
    file = [loc_Exp + specimen + '_resample.csv',
            loc_FEA + specimen + '/' + model_code[:14] + '/' + model_code + '_RFnode.txt',
            loc_FEA + specimen + '/' + model_code[:14] + '/' + model_code + '_RFnodeFix.txt']

    # Load data
    # [A_x, A_y, A_z, A_rx, A_ry, A_rz, a_y, a_f, a_c]
    [_, A_y, _, _, _, _, a_y, a_f, _] = read_resample(file[0])  # load experimental result file (csv)
    [Uy, _] = read_RFnodeFile(file[1])  # read y displacement of moving reference node
    [_, RFy] = read_RFnodeFile(file[2])  # read y reaction force of fixed reference node
    if plot:
        fig1, ax1 = plt.subplots(1, 1, figsize=(9, 6))  # set figure size
        plt.title('Experimental results ' + specimen + ' ' + model_code.split('_')[-1])
        if model_code.split('_')[-1] == 'P':
            ax1.plot(Uy, -RFy + RFy[0], label='FEA (YM PEEK = 25 GPa, u = 0.5)', color='#1f77b4')  # plot fea results
        else:
            ax1.plot(Uy, -RFy + RFy[0], label='FEA (YM Ti = 100 GPa, u = 0.5)', color='#1f77b4')  # plot fea results
        if number in [0, 1, 2, 10, 14, 17]:
            print('Using Acumen displacement.')
            ax1.plot(a_y, a_f - a_f[0], '--', label='Experiment', color='#ff7f0e')  # plot experimental results
        else:
            ax1.plot(A_y, a_f - a_f[0], label='Experiment', color='#ff7f0e')  # plot experimental results
            # ax1.plot(a_y, a_f - a_f[0], '--', label='_nolegend_', color='#ff7f0e')  # plot y data from acumen

        ax1.legend()
        ax1.set_xlabel('Displacement / mm')
        ax1.set_ylabel('Force / N')

        fig1.savefig('/home/biomech/Documents/01_Icotec/02_FEA/91_Pictures/00_Exp_FE/Exp_FE_' +
                     specimen + '_' + model_code + '.png')
    return Uy, RFy, A_y, a_y, a_f


def Peak_exp(ampl_, number_):
    d_ = read_exp_peaks()
    # level = 2 ** (ampl - 2)
    peakF_ = d_['MaxForce'][(number_ - 2) * 7 + ampl_]
    # print(d['DisplacementLevel'][(number - 2) * 7 + ampl])
    # print(d['Specimen'][(number - 2) * 7 + ampl])
    # print(peakF)
    return peakF_


def lin_reg(X, Y):
    X = X.flatten().ravel()
    Y = Y.flatten()
    X = X[X != 0]
    Y = Y[X != 0]
    X = X[Y != 0]
    Y = Y[Y != 0]
    X = sm.add_constant(X)  # Add a constant term to the independent variable array
    mod = sm.OLS(Y, X)  # y, X
    reg = mod.fit()
    return reg, X, Y


# t1 = time.time()
# print('Execution time: '+str(int((time.time()-t1)/60)) + ' min '+str(round(np.mod(time.time()-t1, 60), 1))+' sec.')

# %% Linear regression

# plt.close('all')
# PEEK, without 0 (diff ampl), 24 (Exp. weird)
peek_samples = [2, 5, 7, 8, 10, 13, 15, 16, 18, 21, 23, 26, 29, 31, 32]
# Titanium, without 1 (diff ampl)
ti_samples = [3, 4, 6, 9, 11, 12, 14, 17, 19, 20, 22, 25, 27, 28, 30, 33]

x = 0  # 0 = 0.25 mm, 1 = 0.5 mm, 2 = 1 mm, 3 = 2 mm, 4 = 4 mm, 5 = 8 mm, 6 = 16 mm
lab = ['0.25 mm', '0.5 mm', '1 mm', '2 mm', '4 mm', '8 mm', '16 mm']
x0 = 0
x1 = 7  # max 7
model = '86_L50_S50_D45_d1_02_P'  # automatically switches to titanium for respective samples

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
loglog = 1
alp = 0.3
if loglog:
    F_range = np.array([0, 2.6])
else:
    F_range = np.array([-10, 450])
plt.scatter(-1e9, -1e9, color='k', marker='v', label='PEEK')
plt.scatter(-1e9, -1e9, color='k', marker='s', label='Titanium')
friction = '0.2'
for x in range(x0, x1):
    for i in range(2, 34):  # [2, 3, 4, 5, 10, 11]:  # 2-34 because 0, 1 not existing in data frame
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
        if i == 8:
            plt.scatter(RFy_exp[x, i], RFy_FE[x, i], color=col[x], label=lab[x], marker='v')
        if i in peek_samples:  # P
            plt.scatter(RFy_exp[x, i], RFy_FE[x, i], color=col[x], label='_nolegend_', marker='v')
            RFy_FE_P.append(RFy_FE[x, i])
            RFy_exp_P.append(RFy_exp[x, i])
        elif i in ti_samples:
            plt.scatter(RFy_exp[x, i], RFy_FE[x, i], color=col[x], label='_nolegend_', marker='s')
            RFy_FE_T.append(RFy_FE[x, i])
            RFy_exp_T.append(RFy_exp[x, i])
        elif i == 24:  # HERE exclude sample
            plt.scatter(RFy_exp[x, i], RFy_FE[x, i], color=col[x], label='_nolegend_', marker='v', alpha=alp)
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
axs.set_title('Peak Forces at ' + str(2 ** (x0 - 2)) + ' mm - ' + str(2 ** (x1 - 3)) + ' mm amplitudes, $\mu$ = ' +
              friction)
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

# plt.close('all')
# PEEK, without 0 (diff ampl), 24 (Exp. weird)
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
