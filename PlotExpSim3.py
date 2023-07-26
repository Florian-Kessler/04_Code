import numpy as np
import matplotlib.pyplot as plt
import time
import pandas as pd
import os
import scipy


# from scipy.signal import argrelextrema


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


def read_ARAMIS(file_A):
    stop = 0
    df_ = pd.read_csv(file_A, delimiter=';', skiprows=[0])  # , quoting=csv.QUOTE_NONE, error_bad_lines=False)
    if np.isnan(np.array(df_)[-1][-1]):
        stop = np.where(df_.isnull())[0][0]
    lx = pd.DataFrame(df_, columns=['RodCsys→BoneCsys.LX [mm]'])
    ly = pd.DataFrame(df_, columns=['RodCsys→BoneCsys.LY [mm]'])
    lz = pd.DataFrame(df_, columns=['RodCsys→BoneCsys.LZ [mm]'])
    phiX = pd.DataFrame(df_, columns=['RodCsys→BoneCsys.Phi(X) [°]'])
    thetaY = pd.DataFrame(df_, columns=['RodCsys→BoneCsys.Theta(Y) [°]'])
    psiZ = pd.DataFrame(df_, columns=['RodCsys→BoneCsys.Psi(Z) [°]'])
    t_ = pd.DataFrame(df_, columns=['Time UTC'])
    # print(len(lx))
    if stop:
        return lx[:stop], ly[:stop], lz[:stop], phiX[:stop], thetaY[:stop], psiZ[:stop], t_[:stop]
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
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y_, box, mode='same')
    return y_smooth


t1 = time.time()
plt.close('all')
mue = ['07', '05', '03', '02', '01', '00']

# # # # # INPUT # # # # #
loc = '/home/biomech/Documents/01_Icotec/01_Experiments/00_Data/01_MainStudy/'
specimen = 'S131318_L1_left'
number = ['75']  # simulations


fig1, axs1 = plt.subplots(1, 1, figsize=(9, 6))
plt.title(specimen)

col = ['#0072BD', '#D95319', '#EDB120', '#7E2F8E', '#77AC30', '#4DBEEE', '#A2142F', '#A214CC', '#A2DD2F']

# # # # # Experiments # # # # #

sample = loc + specimen + '_resample.csv'
[ArX, ArY, ArZ, ArrX, ArrY, ArrZ, AcY, AcFy, AcC] = read_resample(sample)
AcFy_smooth = smooth(np.array(AcFy).reshape(len(AcFy),), 3)
peaks = np.array(scipy.signal.argrelextrema(AcFy_smooth, np.less))[0]
valls = np.array(scipy.signal.argrelextrema(AcFy_smooth, np.greater))[0]
axs1.plot(-np.array(ArY), -AcFy_smooth,
          color=col[0])#, alpha=0.2, linestyle=style[0], label='Experiment ' + l_s[0])



#%%


style = ['solid', 'dashed', 'dashed']

l_s = ['Icotec', 'Icotec2', 'DPS']
for i in range(len(samples)):
    if samples[i]:
        [ArX, ArY, ArZ, ArrX, ArrY, ArrZ, AcY, AcFy, AcC] = read_resample(loc + folder + samples[i])
        AcFy_smooth = smooth(np.array(AcFy).reshape(len(AcFy),), 3)
        peaks = np.array(scipy.signal.argrelextrema(AcFy_smooth, np.less))[0]
        valls = np.array(scipy.signal.argrelextrema(AcFy_smooth, np.greater))[0]
        if i == 0 or i == 1:
            # figP.scatter(-np.array(ArY)[peaks], -AcFy_smooth[peaks], color=col[i], s=1, label='Experiment ' + l_s[i])
            figP.plot(-np.array(ArY), -AcFy_smooth,
                      color=col[i], alpha=0.2, linestyle=style[i], label='Experiment ' + l_s[i])
        elif i == 2:
            # figT.scatter(-np.array(ArY)[peaks], -AcFy_smooth[peaks], color=col[0], s=1, label='Experiment ' + l_s[i])
            figT.plot(-np.array(ArY), -AcFy_smooth,
                      color=col[i], alpha=0.2, linestyle=style[0], label='Experiment ' + l_s[i])


# # # # # Simulations # # # # #
loc = '/home/biomech/Documents/01_Icotec/02_FEA/98_Pilots/' + specimen + '/'
for n in range(len(mue)):
    friction = mue[n]
    lab = '\u03BC = 0.' + str(int(friction))
    if friction == '00':
        lab += '5'
    for i in range(len(number)):
        folder = [filename for filename in os.listdir(loc) if filename.startswith(number[i])][0] + '/'
        samples = [filename for filename in sorted(os.listdir(loc + folder + '/'))
                   if filename.endswith('RFnode.txt') and '_' + friction + '_' in filename]
        [uy, rfy] = 2 * [0]
        screw_force = np.zeros([5, 21])
        ang = np.zeros([5, 21])
        Screw_mat = ['P', 'T']
        Sim_mat = ['P', 'T']
        Fpoint = (-np.array([0, 75, 150])) * -37.037 - 1851
        for j in range(len(samples)):
            file_path = loc + folder + samples[j]
            [uy, rf_] = read_RFnodeFile(file_path)
            [u_, rfy] = read_RFnodeFile(file_path.split('.txt')[0] + 'Fix.txt')
            # print('\n' + samples[j])
            # print('Displacement:')
            # print(uy)
            # print('Force:')
            # print(rfy)

            if 'P_P' in samples[j]:
                figP.plot(-uy, rfy, color=col[n], linestyle='solid', label=lab)
                if uy[-1] > 1.01:
                    figP.scatter(-uy[-1], rfy[-1], color='k', marker='x', label='_')
            elif 'T_T' in samples[j]:
                figT.plot(-uy, rfy, color=col[n], linestyle='solid', label=lab)
                if uy[-1] > 1.01:
                    figT.scatter(-uy[-1], rfy[-1], color='k', marker='x', label='_')
            elif 'P_T' in samples[j]:
                figT.plot(-uy, rfy, color=col[n], linestyle='dashed', label='_')
                if uy[-1] > 1.01:
                    figT.scatter(-uy[-1], rfy[-1], color='k', marker='x', label='_')
            elif 'T_P' in samples[j]:
                figP.plot(-uy, rfy, color=col[n], linestyle='solid', label='_')
                if uy[-1] > 1.01:
                    figP.scatter(-uy[-1], rfy[-1], color='k', marker='x', label='_')
            else:
                print('\n . . . Invalid file!\n')
figP.axis([0, 10.5, 0, 250])
figT.axis([0, 15.0, 0, 300])

figP.legend()
figP.set_xlabel('Displacement / mm')
figP.set_ylabel('Force / N')
figT.legend()
figT.set_xlabel('Displacement / mm')
figT.set_ylabel('Force / N')


tRun = time.time()-t1
if tRun >= 3600:
    print('Execution time: ' + str(int(tRun/3600)) + ' h ' + str(int(np.mod(tRun, 3600)/60)) + ' min ' +
          str(round(np.mod(tRun, 60), 1)) + ' sec.')
elif tRun >= 60:
    print('Execution time: ' + str(int(tRun/60)) + ' min ' + str(round(np.mod(tRun, 60), 1)) + ' sec.')
else:
    print('Execution time: ' + str(round(tRun, 1)) + ' sec.')
