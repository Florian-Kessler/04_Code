import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import time
import pandas as pd
import scipy
import statsmodels.api as sm


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


def read_RFnodeFile_uFE(file_):
    # read data from text file
    df_ = np.loadtxt(file_, delimiter='\t')

    t_ = np.array(df_[:, 0])
    d_ = np.array(df_[:, 1])

    return t_, d_


def read_energy(file_):
    # read data from text file
    df_ = np.loadtxt(file_)  # , delimiter='   ')
    t_ = np.array(df_[:, 0])
    e_ = np.array(df_[:, 1])
    return t_, e_


def read_ARAMIS(file_A):
    stop_ = 0
    df_ = pd.read_csv(file_A, delimiter=';', skiprows=[0])  # , quoting=csv.QUOTE_NONE, error_bad_lines=False)
    if np.isnan(np.array(df_)[-1][-1]):
        stop_ = np.where(df_.isnull())[0][0]
    lx = pd.DataFrame(df_, columns=['RodCsys→BoneCsys.LX [mm]'])
    ly = pd.DataFrame(df_, columns=['RodCsys→BoneCsys.LY [mm]'])
    lz = pd.DataFrame(df_, columns=['RodCsys→BoneCsys.LZ [mm]'])
    phiX = pd.DataFrame(df_, columns=['RodCsys→BoneCsys.Phi(X) [°]'])
    thetaY = pd.DataFrame(df_, columns=['RodCsys→BoneCsys.Theta(Y) [°]'])
    psiZ = pd.DataFrame(df_, columns=['RodCsys→BoneCsys.Psi(Z) [°]'])
    t_ = pd.DataFrame(df_, columns=['Time UTC'])
    # print(len(lx))
    if stop_:
        return lx[:stop_], ly[:stop_], lz[:stop_], phiX[:stop_], thetaY[:stop_], psiZ[:stop_], t_[:stop_]
    else:
        return lx, ly, lz, phiX, thetaY, psiZ, t_


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
    # print(len(d_))
    peak_ = peak_.astype(int)
    vall_ = vall_.astype(int)

    return cycle_, d_, f_, peak_, vall_, t_


def read_resample(file_r):
    df_ = pd.read_csv(file_r)
    ArX_ = pd.DataFrame(df_, columns=['Aramis X'])
    ArY_ = pd.DataFrame(df_, columns=['Aramis Y'])
    ArZ_ = pd.DataFrame(df_, columns=['Aramis Z'])
    ArrX_ = pd.DataFrame(df_, columns=['Aramis rX'])
    ArrY_ = pd.DataFrame(df_, columns=['Aramis rY'])
    ArrZ_ = pd.DataFrame(df_, columns=['Aramis rZ'])
    AcY_ = pd.DataFrame(df_, columns=['Acumen Y'])
    AcFy_ = pd.DataFrame(df_, columns=['Acumen Fy'])
    AcC_ = pd.DataFrame(df_, columns=['Acumen C'])
    return ArX_, ArY_, ArZ_, ArrX_, ArrY_, ArrZ_, AcY_, AcFy_, AcC_


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx


def find_first(array, value):
    array = np.asarray(array)
    idx = next(xd for xd, val in enumerate(array)
               if val <= value)
    return idx


def smooth(y_, box_pts):
    box = np.ones(box_pts) / box_pts
    y_smooth = np.convolve(y_, box, mode='same')
    return y_smooth


def lin_reg(X, Y):
    X = X.flatten().ravel()
    Y = Y.flatten()
    # X = X[X != 0]
    # Y = Y[X != 0]
    # X = X[Y != 0]
    # Y = Y[Y != 0]
    X = sm.add_constant(X)  # Add a constant term to the independent variable array
    mod = sm.OLS(Y, X)  # y, X
    reg = mod.fit()
    return reg, X, Y


t1 = time.time()
# plt.close('all')
# mue = ['07', '05', '03', '02', '01', '00']

# # # # # INPUT # # # # #
loc = '/home/biomech/Documents/01_Icotec/01_Experiments/00_Data/01_MainStudy/'
specimen_names = open('/home/biomech/Documents/01_Icotec/Specimens.txt', 'r').read().splitlines()  # Read specimens
col = ['#0072BD', '#D95319', '#EDB120', '#7E2F8E', '#77AC30', '#4DBEEE', '#A2142F', '#A214CC', '#A2DD2F']
number = ['86']  # simulations


#%% Stiffness Exp vs FEA
plt.close('all')

peek_samples = [2, 5, 7, 8, 10, 13, 15, 16, 18, 21, 23, 26, 29, 31, 32]  # without 24
ti_samples = [3, 4, 6, 9, 11, 12, 14, 17, 19, 20, 22, 27, 28, 30, 33]  # without 25
all_samples = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16,
               17, 18, 19, 20, 21, 22, 23, 26, 27, 28, 29, 30, 31, 32, 33]
slope = np.zeros(34)
slope2 = np.zeros(34)
slopeFE02 = np.zeros(34)
slopeFE02_2 = np.zeros(34)
slopeFE05 = np.zeros(34)
slopeFE05_2 = np.zeros(34)
f_rel = np.zeros(34)
f_rel2 = np.zeros(34)
f_abs = np.zeros(34)
cycle = 2
plots = 0
noData02 = []
noData05 = []
samples = all_samples  # peek_samples

for i in samples:  # ti_samples:  # range(2, 34):
    specimen = specimen_names[i]  # 'S131318_L4_right'
    uy = 0
    rfy = 0
    file_path = []
    file_path1 = []
    file_path2 = []
    # # # # # Experiments # # # # #

    sample = loc + specimen + '_resample.csv'
    [ArX, ArY, ArZ, ArrX, ArrY, ArrZ, AcY, AcFy, AcC] = read_resample(sample)
    AcFy = AcFy[5:]
    AcY = AcY[5:]
    AcFy_smooth = smooth(np.array(AcFy).reshape(len(AcFy), ), 4)
    # ArY_smooth = smooth(np.array(ArY).reshape(len(ArY),), 2)
    AcY_smooth = smooth(np.array(AcY).reshape(len(AcY), ), 4)
    peaks = np.array(scipy.signal.argrelextrema(AcY_smooth, np.less))[0]
    valls = np.array(scipy.signal.argrelextrema(AcY_smooth, np.greater))[0]
    s = [peaks[cycle - 1], int((valls[cycle - 1] + peaks[cycle - 1]) / 2)]  # secant between extremum and halfway back
    s2 = [peaks[cycle - 1], valls[cycle - 1]]  # secant between both extrema
    if plots:
        # axs1.plot(np.array(ArY_smooth), AcFy_smooth, color=col[0])
        # axs1.scatter(np.array(ArY_smooth)[valls], np.array(AcFy_smooth)[valls], color='r')
        # axs1.scatter(np.array(ArY_smooth)[peaks], np.array(AcFy_smooth)[peaks], color='r')
        # axs1.plot(ArY_smooth, color=col[0])
        # axs1.plot(AcY, color=col[1])
        fig1, axs1 = plt.subplots(1, 1, figsize=(9, 6))
        plt.title(specimen)
        axs1.plot(AcY[:valls[2]], AcFy_smooth[:valls[2]], color=col[0])
        axs1.scatter(AcY_smooth[valls[:3]], AcFy_smooth[valls[:3]], color=col[1])
        axs1.scatter(AcY_smooth[peaks[:3]], AcFy_smooth[peaks[:3]], color=col[2])
        # axs1.scatter(AcY_smooth[peaks[cycle-1]], AcFy_smooth[peaks[cycle-1]], color=col[2])
        axs1.scatter(AcY_smooth[s[1]], AcFy_smooth[s[1]], color=col[3])
        axs1.plot(AcY_smooth[s], AcFy_smooth[s], 'r--')
        axs1.scatter(AcY_smooth[s2[1]], AcFy_smooth[s2[1]], color=col[3])
        axs1.plot(AcY_smooth[s2], AcFy_smooth[s2], 'g--')
        axs1.set_xlabel('Displacement / mm')
        axs1.set_ylabel('Force / N')
    slope[i] = (AcFy_smooth[s[1]] - AcFy_smooth[s[0]]) / (AcY_smooth[s[1]] - AcY_smooth[s[0]])
    slope2[i] = (AcFy_smooth[s2[1]] - AcFy_smooth[s2[0]]) / (AcY_smooth[s2[1]] - AcY_smooth[s2[0]])
    f_rel[i] = (AcFy_smooth[s[1]] - AcFy_smooth[s[0]])
    f_rel2[i] = (AcFy_smooth[s2[1]] - AcFy_smooth[s2[0]])
    f_abs[i] = AcFy_smooth[s[1]]
    sample = sample.split('_resample')[0].split('/')[-1]
    if i in peek_samples:
        file_path1 = '/home/biomech/Documents/01_Icotec/02_FEA/01_MainStudy/' + sample + '/86_L50_S50_D45' + \
                    '/86_L50_S50_D45_d1_02_P_RFnode.txt'
        file_path2 = '/home/biomech/Documents/01_Icotec/02_FEA/01_MainStudy/' + sample + '/88_L50_S50_D45' + \
                    '/88_L50_S50_D45_d1_02_P_RFnode.txt'

    elif i in ti_samples:
        file_path = '/home/biomech/Documents/01_Icotec/02_FEA/01_MainStudy/' + sample + '/87_L50_S50_D45' + \
                     '/87_L50_S50_D45_d1_02_T_RFnode.txt'
        # file_path2 = '/home/biomech/Documents/01_Icotec/02_FEA/01_MainStudy/' + sample + '/85_L50_S50_D45' + \
        #              '/85_L50_S50_D45_d1_02_T_RFnode.txt'

    if i in peek_samples and len(open(file_path1, 'r').readline()) > 1 and len(open(file_path2, 'r').readline()) > 1:
        [uy1, _] = read_RFnodeFile(file_path1)
        [_, rfy1] = read_RFnodeFile(file_path1.split('.txt')[0] + 'Fix.txt')
        [uy2, _] = read_RFnodeFile(file_path2)
        [_, rfy2] = read_RFnodeFile(file_path2.split('.txt')[0] + 'Fix.txt')
        if len(uy1) < len(uy2):
            uy = uy2
            rfy = rfy2
        else:
            uy = uy1
            rfy = rfy1
    elif i in ti_samples and len(open(file_path, 'r').readline()) > 1:
        [uy, _] = read_RFnodeFile(file_path)
        [_, rfy] = read_RFnodeFile(file_path.split('.txt')[0] + 'Fix.txt')
    if len(uy) >= 42:
        # slopeFE02[i] = -(rfy[-11] - rfy[-6]) / (uy[-11] - uy[-6])
        slopeFE02[i] = -(rfy[31] - rfy[36]) / (uy[31] - uy[36])
        # slopeFE02_2[i] = -(rfy[-11] - rfy[-1]) / (uy[-11] - uy[-1])
        slopeFE02_2[i] = -(rfy[31] - rfy[41]) / (uy[31] - uy[41])
        # print('Exp.: ' + str(np.round(slope[i], 1)) + ' N/mm')
        # print('Exp2: ' + str(np.round(slope2[i], 1)) + ' N/mm')
        # print('FE:   ' + str(np.round(slopeFE02[i], 1)) + ' N/mm')
        # print('FE2:  ' + str(np.round(slopeFE02_2[i], 1)) + ' N/mm\n')
    else:
        print('Data missing for ' + str(sample) + '.\n')
        noData02.append(i)
    if plots:
        plt.figure()
        plt.plot(uy, -rfy)
        plt.scatter(uy[36], -rfy[36], color=col[1])
        plt.scatter(uy[41], -rfy[41], color=col[1])
        plt.scatter(uy[31], -rfy[31], color=col[2])
        plt.plot(uy[[31, 36]], -rfy[[31, 36]], 'r--')
        plt.plot(uy[[31, 41]], -rfy[[31, 41]], 'g--')
        plt.xlabel('Displacement / mm')
        plt.ylabel('Force / N')
    '''
    file_path = '/home/biomech/Documents/01_Icotec/02_FEA/01_MainStudy/' + sample + '/86_L50_S50_D45' + \
                '/86_L50_S50_D45_d1_05_P_RFnode.txt'
    if len(open(file_path, 'r').readline()) > 0:
        [uy, rf_] = read_RFnodeFile(file_path)
        [u_, rfy] = read_RFnodeFile(file_path.split('.txt')[0] + 'Fix.txt')
        if len(uy) == 42:
            slopeFE05[i] = -(rfy[31] - rfy[36]) / (uy[31] - uy[36])
            slopeFE05_2[i] = -(rfy[31] - rfy[41]) / (uy[31] - uy[41])
            print('Exp.: ' + str(np.round(slope[i], 1)) + ' N/mm')
            print('Exp2: ' + str(np.round(slope2[i], 1)) + ' N/mm')
            print('FE:   ' + str(np.round(slopeFE05[i], 1)) + ' N/mm')
            print('FE:   ' + str(np.round(slopeFE05_2[i], 1)) + ' N/mm\n')
        else:
            print('Data missing for ' + str(sample) + '.\n')
            noData05.append(i)
    else:
        print('Data missing for ' + str(sample) + '.\n')
        noData05.append(i)
    if plots:
        plt.figure()
        plt.plot(uy, -rfy)
        plt.scatter(uy[-6], -rfy[-6], color=col[1])
        plt.scatter(uy[-1], -rfy[-1], color=col[1])
        plt.scatter(uy[-11], -rfy[-11], color=col[2])
        plt.plot(uy[[31, 36]], -rfy[[31, 36]], 'r--')
        plt.plot(uy[[31, 41]], -rfy[[31, 41]], 'g--')
        plt.xlabel('Displacement / mm')
        plt.ylabel('Force / N')
    '''

fig2, axs2 = plt.subplots(1, 1)

axs2.scatter(slope[samples], slopeFE02[samples], color=col[0], label='$\mu$ = 0.2')
axs2.scatter(slope2[samples], slopeFE02_2[samples], color=col[0], marker='x', label='$\mu$ = 0.2 (extrema)')
# axs2.scatter(slope[[3, 4, 6, 9, 11, 12, 14, 22, 27, 28, 30, 33]],
#             slopeFE05[[3, 4, 6, 9, 11, 12, 14, 22, 27, 28, 30, 33]], color=col[1], label='$\mu$ = 0.5', alpha=0.7)
# axs2.scatter(slope2[[3, 4, 6, 9, 11, 12, 14, 22, 27, 28, 30, 33]],
#             slopeFE05_2[[3, 4, 6, 9, 11, 12, 14, 22, 27, 28, 30, 33]],
#             color=col[1], marker='x', label='$\mu$ = 0.5 (extrema)', alpha=0.7)
axs2.plot([0, 100], [0, 100], 'k')
axs2.set_xlabel('Stiffness Experiment / N/mm')
axs2.set_ylabel('Stiffness FE / N/mm')
axs2.set_aspect('equal')
plt.legend()

offset = open('/home/biomech/Documents/01_Icotec/01_Experiments/02_Scans/Lever_4mmrad.txt').read().splitlines()
for i in range(len(offset)):
    offset[i] = float(offset[i])
offset = np.array(offset)
fig3, axs3 = plt.subplots(1, 1)
axs3.scatter(offset[samples], slopeFE02[samples], color=col[0])
# plt.scatter(offset[samples], slopeFE02_2[samples], color=col[0], marker='x')
axs3.set_xlabel('Offset / mm')
axs3.set_ylabel('Stiffness / N/mm')

# plt.plot([0, 0], [0, 0], color='w', label=lab_pvalue_T)
# plt.figure()
# plt.plot(slope[peek_samples] / slope[ti_samples])

# plt.figure()
# plt.scatter(slope[peek_samples], slope[ti_samples], color='k')
# plt.xlim([0, 80])
# plt.ylim([0, 80])

# # # # # Stiffness vs F_rel # # # # #
fig4, axs4 = plt.subplots(1, 1)
xdata = slope2[samples]
ydata = f_rel2[samples]
datarange = np.array([0, 90])
axs4.scatter(xdata, ydata, color=col[0])
axs4.set_xlabel('Stiffness Experiment / N/mm')
axs4.set_ylabel('F$_{rel}$ / N')
regression, xx, yy = lin_reg(np.array(xdata), np.array(ydata))
axs4.plot(datarange, datarange * regression.params[1] + regression.params[0], color='k', linestyle='dotted',
          label='R$^2$ = {:0.2f}'.format(np.round(regression.rsquared, 2)))
if regression.pvalues[1] >= 0.05:
    lab_pvalue = 'p = ' + str(np.round(regression.pvalues[1], 2))
else:
    lab_pvalue = 'p < 0.05'
axs4.plot([-1, 0], [-1, 0], color='w', label=lab_pvalue)
plt.legend()

# # # # # Stiffness normalised by F_rel # # # # #
fig41, axs41 = plt.subplots(1, 1)
datarange = np.array([0, 2500])
axs41.scatter(slope[ti_samples] * f_rel[ti_samples], slopeFE02[ti_samples] * f_rel[ti_samples],
              color=col[0], label='Ti')
# axs41.scatter(slope2[ti_samples] * f_rel2[ti_samples], slopeFE02_2[ti_samples] * f_rel2[ti_samples],
#               color=col[0], marker='x', label='Ti (extrema)')
axs41.scatter(slope[peek_samples] * f_rel[peek_samples], slopeFE02[peek_samples] * f_rel[peek_samples],
              color=col[1], label='PEEK')
# axs41.scatter(slope2[peek_samples] * f_rel2[peek_samples], slopeFE02_2[peek_samples] * f_rel2[peek_samples],
#               color=col[1], marker='x', label='PEEK (extrema)')

regression_T, xx_T, yy_T = lin_reg(np.array(slope[ti_samples] * f_rel[ti_samples]),
                                   np.array(slopeFE02[ti_samples] * f_rel[ti_samples]))
axs41.plot(datarange, datarange * regression_T.params[1] + regression_T.params[0], color='k', linestyle='dotted',
           label='R$^2$ = {:0.2f}'.format(np.round(regression_T.rsquared, 2)))
if regression_T.pvalues[1] >= 0.05:
    lab_pvalue_T = 'p = ' + str(np.round(regression_T.pvalues[1], 2))
else:
    lab_pvalue_T = 'p < 0.05'
axs41.plot([-1, 0], [-1, 0], color='w', label=lab_pvalue_T)
regression_P, xx_P, yy_P = lin_reg(np.array(slope[peek_samples] * f_rel[peek_samples]),
                                   np.array(slopeFE02[peek_samples] * f_rel[peek_samples]))
axs41.plot(datarange, datarange * regression_P.params[1] + regression_P.params[0], color='k', linestyle='dashed',
           label='R$^2$ = {:0.2f}'.format(np.round(regression_P.rsquared, 2)))
if regression_P.pvalues[1] >= 0.05:
    lab_pvalue_P = 'p = ' + str(np.round(regression_P.pvalues[1], 2))
else:
    lab_pvalue_P = 'p < 0.05'
axs41.plot([-1, 0], [-1, 0], color='w', label=lab_pvalue_P)
'''
regression_Tx, xx_Tx, yy_Tx = lin_reg(np.array(slope2[ti_samples] * f_rel2[ti_samples]),
                                      np.array(slopeFE02_2[ti_samples] * f_rel2[ti_samples]))
axs41.plot(datarange, datarange * regression_Tx.params[1] + regression_Tx.params[0], color='k', linestyle='dotted',
           label='R$^2$ = {:0.2f}'.format(np.round(regression_Tx.rsquared, 2)))
if regression_Tx.pvalues[1] >= 0.05:
    lab_pvalue_Tx = 'p = ' + str(np.round(regression_Tx.pvalues[1], 2))
else:
    lab_pvalue_Tx = 'p < 0.05'
axs41.plot([-1, 0], [-1, 0], color='w', label=lab_pvalue_Tx)
regression_Px, xx_P, yy_P = lin_reg(np.array(slope2[peek_samples] * f_rel2[peek_samples]),
                                    np.array(slopeFE02_2[peek_samples] * f_rel2[peek_samples]))
axs41.plot(datarange, datarange * regression_Px.params[1] + regression_Px.params[0], color='k', linestyle='dashed',
           label='R$^2$ = {:0.2f}'.format(np.round(regression_Px.rsquared, 2)))
if regression_Px.pvalues[1] >= 0.05:
    lab_pvalue_Px = 'p = ' + str(np.round(regression_Px.pvalues[1], 2))
else:
    lab_pvalue_Px = 'p < 0.05'
axs41.plot([-1, 0], [-1, 0], color='w', label=lab_pvalue_Px)
'''
plt.legend()
axs41.set_xlabel('Stiffness Experiment * F$_{rel}$')
axs41.set_ylabel('Stiffness FEA * F$_{rel}$')
axs41.set_aspect('equal')
axs41.set_xlim(datarange)
axs41.set_ylim(datarange)

# # # # # Stiffness vs F_abs # # # # #
fig5, axs5 = plt.subplots(1, 1)
axs5.scatter(slope[ti_samples], f_abs[ti_samples], color=col[0], label='Ti')
axs5.scatter(slope[peek_samples], f_abs[peek_samples], color=col[1], label='PEEK')
axs5.set_xlabel('Stiffness Experiment / N/mm')
axs5.set_ylabel('F$_{abs}$ / N')
datarange = np.array([0, 90])

regression_T, xx_T, yy_T = lin_reg(np.array(slope[ti_samples]), f_abs[ti_samples])
axs5.plot(datarange, datarange * regression_T.params[1] + regression_T.params[0], color='k', linestyle='dotted',
           label='R$^2$ = {:0.2f}'.format(np.round(regression_T.rsquared, 2)))
if regression_T.pvalues[1] >= 0.05:
    lab_pvalue_T = 'p = ' + str(np.round(regression_T.pvalues[1], 2))
else:
    lab_pvalue_T = 'p < 0.05'
axs5.plot([-1, 0], [-1, 0], color='w', label=lab_pvalue_T)

regression_P, xx_P, yy_P = lin_reg(np.array(slope[peek_samples]), f_abs[peek_samples])
axs5.plot(datarange, datarange * regression_P.params[1] + regression_P.params[0], color='k', linestyle='dashed',
           label='R$^2$ = {:0.2f}'.format(np.round(regression_P.rsquared, 2)))
if regression_P.pvalues[1] >= 0.05:
    lab_pvalue_P = 'p = ' + str(np.round(regression_P.pvalues[1], 2))
else:
    lab_pvalue_P = 'p < 0.05'
axs5.plot([-1, 0], [-1, 0], color='w', label=lab_pvalue_P)

plt.legend()
axs5.set_xlim(datarange)
#%% FEA catalogue
for i in ti_samples:  # range(2, 34):
    specimen = specimen_names[i]  # 'S131318_L4_right'
    uy = 0
    rfy = 0
    # # # # # Experiments # # # # #
    sample = loc + specimen + '_resample.csv'
    [ArX, ArY, ArZ, ArrX, ArrY, ArrZ, AcY, AcFy, AcC] = read_resample(sample)
    print(specimen)
    sample = sample.split('_resample')[0].split('/')[-1]
    if i in peek_samples:
        file_path1 = '/home/biomech/Documents/01_Icotec/02_FEA/01_MainStudy/' + sample + '/86_L50_S50_D45' + \
                    '/86_L50_S50_D45_d1_02_P_RFnode.txt'
        file_path2 = '/home/biomech/Documents/01_Icotec/02_FEA/01_MainStudy/' + sample + '/88_L50_S50_D45' + \
                    '/88_L50_S50_D45_d1_02_P_RFnode.txt'
    elif i in ti_samples:
        # file_path1 = '/home/biomech/Documents/01_Icotec/02_FEA/01_MainStudy/' + sample + '/85_L50_S50_D45' + \
        #             '/85_L50_S50_D45_d1_02_T_RFnode.txt'
        file_path2 = '/home/biomech/Documents/01_Icotec/02_FEA/01_MainStudy/' + sample + '/87_L50_S50_D45' + \
                    '/87_L50_S50_D45_d1_02_T_RFnode.txt'
        file_path1 = file_path2
    else:
        print('Specimen ' + specimen + ' not in list.')
        continue
    if len(open(file_path1, 'r').readline()) > 0 and len(open(file_path2, 'r').readline()) > 0:
        [uy1, _] = read_RFnodeFile(file_path1)
        [_, rfy1] = read_RFnodeFile(file_path1.split('.txt')[0] + 'Fix.txt')
        [uy2, _] = read_RFnodeFile(file_path2)
        [_, rfy2] = read_RFnodeFile(file_path2.split('.txt')[0] + 'Fix.txt')
        if len(uy1) < len(uy2):
            uy = uy2
            rfy = rfy2
        else:
            uy = uy1
            rfy = rfy1
    plt.figure()
    plt.plot(AcY, AcFy-AcFy['Acumen Fy'][0], label='Experiment')
    plt.plot(uy, -rfy, label='FEA')
    plt.xlabel('Displacement / mm')
    plt.ylabel('Force / N')
    plt.title(specimen)
    plt.legend()
    plt.savefig('/home/biomech/Documents/01_Icotec/02_FEA/91_Pictures/01_Catalogue_FEA_Exp/' + sample + '.png')
#%% uFE test files
# file = '/home/biomech/DATA/01_Icotec/02_FEA/02_uFE/Tests/test_14_RFnode.txt'
# file = '/home/biomech/DATA/01_Icotec/02_FEA/02_uFE/Tests/hFE_P3_RFnode.txt'
# file = '/home/biomech/DATA/01_Icotec/02_FEA/02_uFE/Tests/test_A_3T_RFnode.txt'
file = '/home/biomech/DATA/01_Icotec/02_FEA/02_uFE/Tests/test_A_2_RFnode.txt'
[_, uy] = read_RFnodeFile_uFE(file)
[_, rfy] = read_RFnodeFile_uFE(file.split('.txt')[0] + 'Fix.txt')

sample = loc + specimen_names[2] + '_resample.csv'
[ArX, ArY, ArZ, ArrX, ArrY, ArrZ, AcY, AcFy, AcC] = read_resample(sample)

# file2 = '/home/biomech/DATA/01_Icotec/02_FEA/02_uFE/Tests/test_15_RFnode.txt'
# file2 = '/home/biomech/DATA/01_Icotec/02_FEA/02_uFE/Tests/test_A_2T_RFnode.txt'
# [_, uy2] = read_RFnodeFile_uFE(file2)
# [_, rfy2] = read_RFnodeFile_uFE(file2.split('.txt')[0] + 'Fix.txt')

plt.figure()
plt.plot(AcY, AcFy-AcFy['Acumen Fy'][0], label='Experiment')
plt.plot(uy, rfy, label='FE')
# plt.plot(uy2, rfy2, label='Ti')
# plt.scatter(uy[-1], rfy[-1], color='r', marker='x')
# plt.scatter(uy2[-1], rfy2[-1], color='r', marker='x')
plt.title('μFE tests (Ti)')
plt.xlabel('Displacement / mm')
plt.ylabel('Force / N')
plt.legend()

#%% Energy plot
file = '/home/biomech/DATA/01_Icotec/02_FEA/02_uFE/Tests/test_7_ke.txt'
t, ke = read_energy(file)
t = t[1:]
ke = ke[1:]
file = '/home/biomech/DATA/01_Icotec/02_FEA/02_uFE/Tests/test_7_ie.txt'
_, ie = read_energy(file)
ie = ie[1:]
plt.figure()
# plt.plot(t, ke, label='$E_k$')
# plt.plot(t, ie, label='$E_i$')
plt.plot(t, ke/(ie + ke), label='$E_k$/$E_{total}$')
plt.plot([t[0], t[-1]], [0.1, 0.1], 'k--')
plt.legend()

file = '/home/biomech/DATA/01_Icotec/02_FEA/02_uFE/Tests/test_8_ke.txt'
t, ke = read_energy(file)
t = t[1:]
ke = ke[1:]
file = '/home/biomech/DATA/01_Icotec/02_FEA/02_uFE/Tests/test_8_ie.txt'
_, ie = read_energy(file)
ie = ie[1:]
plt.figure()
# plt.plot(t, ke)
# plt.plot(t, ie)
plt.plot(t, ke/(ie + ke), label='$E_k$/$E_{total}$')
plt.plot([t[0], t[-1]], [0.1, 0.1], 'k--')
plt.legend()
#%% BVTV Plot
radius = [6]
radius_txt = [6]
# fig, axs = plt.subplots(1, 1, figsize=(9, 6))
plt.figure()
for j in range(len(radius)):
    bvtvList = []
    specList = []
    indexList = []
    for i in range(len(specimen_names)):
        loc_ = '/home/biomech/DATA/01_Icotec/01_Experiments/02_Scans/BVTV/BVTV_along_load_'
        bvtv = np.load(loc_ + specimen_names[i] + '_' + str(radius[j]) + 'mm.npy')
        # axs.plot(bvtv)
        # test line plot seaborn
        bvtvList += bvtv.T.tolist()
        specList += len(bvtv)*[i]
        indexList += range(len(bvtv))
    dfPlot = pd.DataFrame()
    dfPlot['BV/TV'] = bvtvList
    dfPlot['specimen'] = specList
    dfPlot['Slice Number'] = indexList
    sns.lineplot(data=dfPlot, x='Slice Number', y='BV/TV', errorbar='sd', label='Radius: ' + str(radius_txt[j]) + ' mm')
#%% BVTV vs Exp

bvtv_range = np.array([0, 0.6])
radius = [4]
offset = 0
start = 0
RR = np.array([])
# stop = [235, 236, 237]  # 236 for along, mean(3500:5700)
# stop = [263, 264, 265]  # 264 for along, mean(3500:5700)
# stop = [90, 91, 92]  # 91 for along_load, mean(3500:5700)
# stop = [81, 82, 83]  # 82 for along_load, min(3500:5700)
stop = [82]
weight = 1
temp = np.array([])
col = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf',
       '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
mark = ['o', 's']
for j in range(len(stop)):
    plt.figure()
    xdata = []
    ydata = []
    for i in range(len(specimen_names)):
        loc_ = '/home/biomech/DATA/01_Icotec/01_Experiments/02_Scans/BVTV/BVTV_along_'
        bvtv = np.load(loc_ + specimen_names[i] + '_' + str(radius[0]) + 'mm.npy')
        offset = int((bvtv != 0).argmax(axis=0)/weight)
        sample = loc + specimen_names[i] + '_resample.csv'
        [ArX, ArY, ArZ, ArrX, ArrY, ArrZ, AcY, AcFy, AcC] = read_resample(sample)
        # plt.scatter(np.mean(bvtv[start:stop[j]], axis=0), np.max(-AcFy, axis=0))
        AcFy = AcFy - AcFy['Acumen Fy'][0]
        xdata = np.append(xdata, np.mean(bvtv[start+offset:stop[j]+offset], axis=0))
        ydata = np.append(ydata, np.min(AcFy[3500:5700], axis=0))
        print('Offset for ' + specimen_names[i] + ' = \t' + str(np.round(offset*0.0606995, 3)) + ' mm')
        temp = np.append(temp, offset*0.0606995)
        if i in ti_samples:
            plt.scatter(xdata[i], ydata[i], color=col[int(i/2)], marker=mark[0])
        elif i in peek_samples:
            plt.scatter(xdata[i], ydata[i], color=col[int(i/2)], marker=mark[1])
    plt.xlabel('BV/TV for slices ' + str(start) + ' to ' + str(stop[j]))
    plt.ylabel('Mean(force / N) of last cycle')
    plt.title('Radius: ' + str(radius[0]) + ', weighted offset: w = 1/' + str(weight))

    regression_T, xx_T, yy_T = lin_reg(np.array(xdata), np.array(ydata))
    plt.plot(bvtv_range, bvtv_range * regression_T.params[1] + regression_T.params[0], color='k', linestyle='dotted',
             label='Titanium:')
    if regression_T.pvalues[1] >= 0.05:
        lab_pvalue_T = 'p = ' + str(np.round(regression_T.pvalues[1], 2))
    else:
        lab_pvalue_T = 'p < 0.05'
    plt.plot([0, 0], [0, 0], color='w', linestyle='dashed',
             label='R$^2$ = {:0.2f}'.format(np.round(regression_T.rsquared, 2)))
    plt.plot([0, 0], [0, 0], color='w', label=lab_pvalue_T)
    plt.legend()
    RR = np.append(RR, regression_T.rsquared)
    # plt.close('all')
# plt.figure()
# plt.plot(temp)
# plt.plot([0, 33], [np.mean(temp), np.mean(temp)])
print('Mean offset:\t\t\t\t\t' + str(np.round(np.mean(temp), 3)) + ' mm')

#%% BVTV vs Exp MOMENT
plot = 0
bvtv_range = np.array([0, 0.4])
radius = [4]
offset = 0
start = 0
RR = np.array([])
# stop = [235, 236, 237]  # 236 for along, mean(3500:5700)
# stop = [263, 264, 265]  # 264 for along, mean(3500:5700)
# stop = [90, 91, 92]  # 91 for along_load, mean(3500:5700)
# stop = [81, 82, 83]  # 82 for along_load, min(3500:5700)
# stop = [290] for along, moment(min(3500:5700))
stop = [290]

temp = np.array([])
col = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf',
       '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
mark = ['o', 's']
for j in range(len(stop)):
    fig6, axs6 = plt.subplots(1, 1)
    xdata = []
    ydata = []
    cog = []
    moment = []
    for i in range(len(specimen_names)):
        loc_ = '/home/biomech/DATA/01_Icotec/01_Experiments/02_Scans/BVTV/BVTV_along_'
        bvtv = np.load(loc_ + specimen_names[i] + '_' + str(radius[0]) + 'mm.npy')
        offset = int((bvtv != 0).argmax(axis=0))
        sample = loc + specimen_names[i] + '_resample.csv'
        [ArX, ArY, ArZ, ArrX, ArrY, ArrZ, AcY, AcFy, AcC] = read_resample(sample)
        # plt.scatter(np.mean(bvtv[start:stop[j]], axis=0), np.max(-AcFy, axis=0))
        AcFy = AcFy - AcFy['Acumen Fy'][0]
        xdata = np.append(xdata, np.mean(bvtv[start+offset:stop[j]+offset], axis=0))
        #ydata = np.append(ydata, np.min(AcFy[3500:5700], axis=0))

        # Center of Gravity to calculate moment
        x = np.arange(start+offset, stop[j]+offset)
        m = np.array(bvtv[start+offset:stop[j]+offset])
        COG = np.sum(m*x) / np.sum(m)
        cog = np.append(cog, COG)
        # print('COG:\t' + str(np.round(COG, 2)) + ' mm')

        # Moment
        lever = COG + offset*0.0606995
        moment = np.append(moment, lever/1000 * np.min(AcFy[3500:5700], axis=0))
        # print('Moment:\t' + str(np.round(moment[i], 2)) + ' Nm\n')
        ydata = np.append(ydata, moment[i])

        # print('Offset for ' + specimen_names[i] + ' = \t' + str(np.round(offset*0.0606995, 3)) + ' mm')
        temp = np.append(temp, offset*0.0606995)
        if i in ti_samples:
            axs6.scatter(xdata[i], ydata[i], color=col[int(i/2)], marker=mark[0])
        elif i in peek_samples:
            axs6.scatter(xdata[i], ydata[i], color=col[int(i/2)], marker=mark[1])
    axs6.set_xlabel('BV/TV for slices ' + str(start) + ' to ' + str(stop[j]))
    axs6.set_ylabel('Mean(force / N) of last cycle')
    axs6.set_title('Radius: ' + str(radius[0]) + ' mm')

    regression_T, xx_T, yy_T = lin_reg(np.array(xdata), np.array(ydata))
    axs6.plot(bvtv_range, bvtv_range * regression_T.params[1] + regression_T.params[0], color='k', linestyle='dotted',
              label='_nolegend_')
    if regression_T.pvalues[1] >= 0.05:
        lab_pvalue_T = 'p = ' + str(np.round(regression_T.pvalues[1], 2))
    else:
        lab_pvalue_T = 'p < 0.05'
    axs6.plot([0, 0], [0, 0], color='w', linestyle='dashed',
              label='R$^2$ = {:0.2f}'.format(np.round(regression_T.rsquared, 2)))
    axs6.plot([0, 0], [0, 0], color='w', label=lab_pvalue_T)
    axs6.legend()
    RR = np.append(RR, regression_T.rsquared)
    # plt.close('all')
# plt.figure()
# plt.plot(temp)
# plt.plot([0, 33], [np.mean(temp), np.mean(temp)])
print('Mean offset:\t' + str(np.round(np.mean(temp), 3)) + ' mm')
