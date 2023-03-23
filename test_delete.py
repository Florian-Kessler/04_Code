import numpy as np
import matplotlib.pyplot as plt
import time
import pandas as pd
import os


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
    print(len(lx))
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

    peak_ = peak_.astype(int)
    vall_ = vall_.astype(int)

    return cycle_, d_, f_, peak_, vall_, t_


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx


def find_first(array, value):
    array = np.asarray(array)
    idx = next(xd for xd, val in enumerate(array)
               if val <= value)
    return idx


t1 = time.time()
plt.close('all')

# # # # # INPUT # # # # #
loc = '/home/biomech/Documents/01_Icotec/01_Experiments/00_Data/'
specimen = '05_Pilot5'
number = ['13']  # , '13']  # simulations

fig1, figP = plt.subplots(1, 1, figsize=(9, 6))
plt.title('PEEK (YM = 15 GPa)')
fig2, figT = plt.subplots(1, 1, figsize=(9, 6))
plt.title('Ti (YM = 100 GPa)')
col = ['#0072BD', '#D95319', '#EDB120', '#7E2F8E', '#77AC30', '#4DBEEE', '#A2142F', '#A214CC', '#A2DD2F']

no = specimen.split('_')[0]
folder = [filename for filename in os.listdir(loc) if filename.startswith(no)][0] + '/'
sampleIPD = ([filename for filename in os.listdir(loc + folder)
              if 'aramis' in filename and 'icotec' in filename and 'kwire' not in filename] or [None])[0]
sampleIPF = ([filename for filename in os.listdir(loc + folder)
              if 'acumen' in filename and 'icotec' in filename and 'kwire' not in filename] or [None])[0]
sampleIKD = ([filename for filename in os.listdir(loc + folder)
              if 'aramis' in filename and 'kwire' in filename] or [None])[0]
sampleIKF = ([filename for filename in os.listdir(loc + folder)
              if 'acumen' in filename and 'kwire' in filename] or [None])[0]
sampleTiD = ([filename for filename in os.listdir(loc + folder)
              if 'aramis' in filename and 'DPS' in filename] or [None])[0]
sampleTiF = ([filename for filename in os.listdir(loc + folder)
              if 'acumen' in filename and 'DPS' in filename] or [None])[0]
samplesD = sampleIPD, sampleIKD, sampleTiD
samplesF = sampleIPF, sampleIKF, sampleTiF
label_screw = ['Icotec', 'Icotec2', 'DPS']
for i in range(len(samplesD)):
    if samplesD[i] and samplesF[i]:
        [x, y, z, rX, rY, rZ, t] = read_ARAMIS(loc + folder + samplesD[i])
        y = y.to_numpy().flatten()
        t = np.array(t).flatten()
        for j in range(len(t)):
            hhmmss = t[j].split(' ')[1]
            hh = hhmmss.split(':')[0]
            mm = hhmmss.split(':')[1]
            ss = hhmmss.split(':')[2].split(',')[0]
            fr = hhmmss.split(':')[2].split(',')[1]
            t[j] = int(hh) * 3600 + int(mm) * 60 + int(ss) + int(fr) / 1000
        peak = [0]
        vall = [0]
        for s in range(int(t[0]), int(np.max(t) - 1)):
            arr = np.where(t.astype(int) == s)[0]
            if y[arr[np.argmin(y[arr])]] < y[peak[-1]]:
                peak.append(arr[int(np.argmin(y[arr]))])
                vall.append(arr[int(np.argmax(y[arr]))])
            else:
                peak.append(peak[-1])
                vall.append(vall[-1])

        [C, D, F, _, _, T] = read_acumen(loc + folder + samplesF[i])
        peakAc = [0]
        vallAc = [0]
        for s in range(int(T[0]), int(np.max(T) - 1)):
            arrAc = np.where(T.astype(int) == s)[0]
            if D[arrAc[np.argmin(D[arrAc])]] < D[peakAc[-1]]:
                peakAc.append(arrAc[int(np.argmin(D[arrAc]))])
                vallAc.append(arrAc[int(np.argmax(D[arrAc]))])
            else:
                peakAc.append(peakAc[-1])
                vallAc.append(vallAc[-1])
        if len(peakAc) > len(peak):
            start = len(peakAc) - len(peak)
            if i == 0 or 1:
                figP.scatter(-y[peak], -F[peakAc[start:]], color=col[i], s=1, label='Experiment ' + label_screw[i])
            elif i == 2:
                figT.scatter(-y[peak], -F[peakAc[start:]], color=col[i], s=1, label='Experiment ' + label_screw[i])
        elif len(peak) > len(peakAc):
            start = len(peak) - len(peakAc)
            if i == 0 or 1:
                figP.scatter(-y[peak[start:]], -F[peakAc], color=col[i], s=1, label='Experiment ' + label_screw[i])
            elif i == 2:
                figT.scatter(-y[peak[start:]], -F[peakAc], color=col[i], s=1, label='Experiment ' + label_screw[i])
        else:
            if i == 0 or 1:
                figP.scatter(-y[peak], -F[peakAc], color=col[i], s=1, label='Experiment ' + label_screw[i])
            elif i == 2:
                figT.scatter(-y[peak], -F[peakAc], color=col[i], s=1, label='Experiment ' + label_screw[i])


loc = '/home/biomech/Documents/01_Icotec/02_FEA/98_Pilots/' + specimen + '/'
for i in range(len(number)):
    folder = [filename for filename in os.listdir(loc) if filename.startswith(number[i])][0] + '/'
    samples = [filename for filename in os.listdir(loc + folder + '/')
               if filename.endswith('RFnode.txt') and '_02_' in filename]
    [uy, rfy] = 2 * [0]
    screw_force = np.zeros([5, 21])
    ang = np.zeros([5, 21])
    Screw_mat = ['P', 'T']
    Sim_mat = ['P', 'T']
    Fpoint = (-np.array([0, 75, 150])) * -37.037 - 1851
    for j in range(len(samples)):
        lab = 'Screw excess = ' + folder.split('_S')[-1][0] + '.' + folder.split('_S')[-1][1] + ' mm, ' \
              # + 'diameter = ' + folder.split('_D')[-1][0] + '.' + folder.split('_D')[-1][1] + ' mm'
        file_path = loc + folder + samples[j]
        [uy, rf_] = read_RFnodeFile(file_path)
        [u_, rfy] = read_RFnodeFile(file_path.split('.txt')[0] + 'Fix.txt')
        if 'P_P' in samples[j]:
            figP.plot(-uy, rfy, color=col[i], linestyle='solid', label='FE Icotec')
            if rfy[-1] > 20:
                figP.scatter(-uy[-1], rfy[-1], color='k', marker='x', label='_nolegend_')
        elif 'T_T' in samples[j]:
            figT.plot(-uy, rfy, color=col[i], linestyle='solid', label='DPS')  # figT
            if rfy[-1] > 20:
                figT.scatter(-uy[-1], rfy[-1], color='k', marker='x', label='_nolegend_')  # figT
        elif 'P_T' in samples[j]:
            figT.plot(-uy, rfy, color=col[i], linestyle='dashed', label='DPS (C)')  # figT
            if rfy[-1] > 20:
                figT.scatter(-uy[-1], rfy[-1], color='k', marker='x', label='_nolegend_')  # figT
        elif 'T_P' in samples[j]:
            figP.plot(-uy, rfy, color=col[1], linestyle='solid', label='FE Icotec2')  # here:
            # changed dashed to solid and color i to 1
            if rfy[-1] > 20:
                figP.scatter(-uy[-1], rfy[-1], color='k', marker='x', label='_nolegend_')
        else:
            print('\n . . . Invalid file!\n')
figP.axis([0, 10.5, 0, 250])
figT.axis([0, 10.5, 0, 350])

figP.legend()


'''
fileAramis = '/home/biomech/Documents/01_Icotec/01_Experiments/00_Data/05_Pilot5/ICOTEC_S130684_L4_aramis.csv'
fileAcumen = '/home/biomech/Documents/01_Icotec/01_Experiments/00_Data/05_Pilot5/ICOTEC_S130684_L4_accumen.csv'
# df = read_ARAMIS(fileAramis)
[x, y, z, rX, rY, rZ, t] = read_ARAMIS(fileAramis)
[C, D, F, P, V, T] = read_acumen(fileAcumen)


y = y.to_numpy().flatten()

t = np.array(t).flatten()
for i in range(len(t)):
    hhmmss = t[i].split(' ')[1]
    hh = hhmmss.split(':')[0]
    mm = hhmmss.split(':')[1]
    ss = hhmmss.split(':')[2].split(',')[0]
    fr = hhmmss.split(':')[2].split(',')[1]
    t[i] = int(hh)*3600 + int(mm)*60 + int(ss) + int(fr)/1000

plt.figure()
peak = [0]
vall = [0]
for s in range(int(t[0]), int(np.max(t)-1)):
    arr = np.where(t.astype(int) == s)[0]
    if y[arr[np.argmin(y[arr])]] < y[peak[-1]]:  # and y[arr[np.argmax(y[arr])]] < y[vall[-1]]:
        peak.append(arr[int(np.argmin(y[arr]))])
        vall.append(arr[int(np.argmax(y[arr]))])
    else:
        peak.append(peak[-1])
        vall.append(vall[-1])


start = 18
peakAc = [0]
vallAc = [0]
for s in range(int(T[0]), int(np.max(T)-1)):
    arrAc = np.where(T.astype(int) == s)[0]
    if D[arrAc[np.argmin(D[arrAc])]] < D[peakAc[-1]]:  # and D[arrAc[np.argmax(D[arrAc])]] > D[vallAc[-1]]:
        peakAc.append(arrAc[int(np.argmin(D[arrAc]))])
        vallAc.append(arrAc[int(np.argmax(D[arrAc]))])
    else:
        peakAc.append(peakAc[-1])
        vallAc.append(vallAc[-1])

figP.scatter(-y[peak[start:]], -F[peakAc])
'''

tRun = time.time()-t1
if tRun >= 3600:
    print('Execution time: ' + str(int(tRun/3600)) + ' h ' + str(int(np.mod(tRun, 3600)/60)) + ' min ' +
          str(round(np.mod(tRun, 60), 1)) + ' sec.')
elif tRun >= 60:
    print('Execution time: ' + str(int(tRun/60)) + ' min ' + str(round(np.mod(tRun, 60), 1)) + ' sec.')
else:
    print('Execution time: ' + str(round(tRun, 1)) + ' sec.')
