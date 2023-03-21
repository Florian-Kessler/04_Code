import numpy as np
import matplotlib.pyplot as plt
import time
import pandas as pd
import os
from scipy.signal import find_peaks


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

specimen = '05_Pilot5'

'''
fig, ax1 = plt.subplots(1, 1, figsize=(9, 6))
ax2 = ax1.twinx()



# PEEK1 = '01_Pilot1/S131840_L4_S1_PEEK.csv'
# TITAN1 = '01_Pilot1/S191840_L4_S2_DPS.csv'
PEEK = ''
TITAN = ''
cor = 0
if '03' in specimen:
    PEEK = '03_Pilot3/ICOTEC_S130672_L5_icotec_accumen.csv'
    TITAN = '03_Pilot3/ICOTEC_S130672_L5_DPS_accumen.csv'
    cor = float(2.0)  # Correction for displacement offset for Ti
elif '04' in specimen:
    PEEK = '04_Pilot4/ICOTEC_S130672_L4_icotec_accumen.csv'
    TITAN = '04_Pilot4/ICOTEC_S130672_L4_icotec_kwire_accumen.csv'
elif '05' in specimen:
    PEEK = '05_Pilot5/ICOTEC_S130684_L4_accumen.csv'
    TITAN = '05_Pilot5/ICOTEC_S130684_L4_kwire_accumen.csv'

loc = '/home/biomech/Documents/01_Icotec/01_Experiments/00_Data/'

file = [loc + PEEK, loc + TITAN]
exp = {}
col = [['#1f6a06', '#0a1eb2'],
       ['#373737', '#b20a0a']]
labels = ['PEEK', 'Ti']
for i in range(0, len(file)):
    [exp['cycle' + labels[i]], exp['d' + labels[i]], exp['f' + labels[i]], exp['peak' + labels[i]],
     exp['vall' + labels[i]]] = read_acumen(str(file[i]))
    if 'Ti' in labels[i]:
        exp['dTi'] = exp['dTi'] + cor
    ax1.plot(exp['cycle' + labels[i]][exp['peak' + labels[i]]], exp['d' + labels[i]][exp['peak' + labels[i]]],
             color=col[i][0], label='_nolegend_')
    ax1.plot(exp['cycle' + labels[i]][exp['vall' + labels[i]]], exp['d' + labels[i]][exp['vall' + labels[i]]],
             color=col[i][0], alpha=0.75, label='_nolegend_')
    ax2.plot(exp['cycle' + labels[i]][exp['peak' + labels[i]]], exp['f' + labels[i]][exp['peak' + labels[i]]],
             color=col[i][1], label='Force ' + str(labels[i]))
    ax2.plot(exp['cycle' + labels[i]][exp['vall' + labels[i]]], exp['f' + labels[i]][exp['vall' + labels[i]]],
             color=col[i][1], alpha=0.75, label='_nolegend_')
    if 'Ti' in labels[i]:
        ax2.plot([-.2, .2], color=col[i][0], label='Displacement ' + str(labels[i]) + ', corr.: ' + str(cor) + ' mm')
    else:
        ax2.plot([-.2, .2], color=col[i][0], label='Displacement ' + str(labels[i]))

ax2.spines.right.set_visible(True)
ax1.set_xlabel('Cycle Number')
ax1.set_ylabel('Displacement / mm')
ax2.set_ylabel('Force / N')
ax1.axis([0, int(np.ceil(np.max(exp['cycle' + labels[1]])/1000)*1000),  -10, 0])
ax2.axis([0, int(np.ceil(np.max(exp['cycle' + labels[1]])/1000)*1000), -300, 0])
ax2.legend(['PEEK Force', 'PEEK Displacement', 'Ti Force', 'Ti Displacement'])

try:
    cyc1T = find_first(exp['fTi'][exp['peakTi']], -75)
    try:
        cyc2T = find_first(exp['fTi'][exp['peakTi']], -150)
    except StopIteration:
        cyc2T = 0
except StopIteration:
    cyc1T = 0
    cyc2T = 0

try:
    cyc1P = find_first(exp['fPEEK'][exp['peakPEEK']], -75)
    try:
        cyc2P = find_first(exp['fPEEK'][exp['peakPEEK']], -150)
    except StopIteration:
        cyc2P = 0
except StopIteration:
    cyc1P = 0
    cyc2P = 0
'''

fig1, figP = plt.subplots(1, 1, figsize=(9, 6))
plt.title('Inverse simulations PEEK (YM = 15 GPa), force controlled')
fig2, figT = plt.subplots(1, 1, figsize=(9, 6))
plt.title('Inverse simulations Ti (YM = 110 GPa), force controlled')
col = ['#0072BD', '#D95319', '#EDB120', '#7E2F8E', '#77AC30', '#4DBEEE', '#A2142F', '#A214CC', '#A2DD2F']

# number = ['00', '01', '02', '10', '11', '12']
number = ['03', '13']
for i in range(len(number)):
    loc = '/home/biomech/Documents/01_Icotec/02_FEA/98_Pilots/' + specimen + '/'
    folder = [filename for filename in os.listdir(loc) if filename.startswith(number[i])][0] + '/'
    samples = [filename for filename in os.listdir(loc + folder + '/') if filename.endswith('RFnode.txt')
               and '_02_' in filename]
    [uy, rfy] = 2 * [0]
    screw_force = np.zeros([5, 21])
    ang = np.zeros([5, 21])
    Screw_mat = ['P', 'T']
    Sim_mat = ['P', 'T']
    Fpoint = (-np.array([0, 75, 150])) * -37.037 - 1851
    for j in range(len(samples)):
        lab = 'Screw excess = ' + folder.split('_S')[-1][0] + '.' + folder.split('_S')[-1][1] + ' mm, ' \
              # + 'diameter = ' + folder.split('_D')[-1][0] + '.' + folder.split('_D')[-1][1] + ' mm'
        # lab = 'Unidirectional'
        file_path = loc + folder + samples[j]
        [uy, rf_] = read_RFnodeFile(file_path)  # switch for inverse loading
        [u_, rfy] = read_RFnodeFile(file_path.split('.txt')[0] + 'Fix.txt')  # switch for inverse loading
        if 'P_P' in samples[j]:
            figP.plot(-uy, rfy, color=col[i], linestyle='solid', label='_nolegend_')
            if rfy[-1] > 20:
                figP.scatter(-uy[-1], rfy[-1], color='k', marker='x', label='_nolegend_')
        elif 'T_T' in samples[j]:
            figT.plot(-uy, rfy, color=col[i], linestyle='solid', label='_nolegend_')  # figT
            if rfy[-1] > 20:
                figT.scatter(-uy[-1], rfy[-1], color='k', marker='x', label='_nolegend_')  # figT
        elif 'P_T' in samples[j]:
            figT.plot(-uy, rfy, color=col[i], linestyle='dashed', label='_nolegend_')  # figT
            if rfy[-1] > 20:
                figT.scatter(-uy[-1], rfy[-1], color='k', marker='x', label='_nolegend_')  # figT
        elif 'T_P' in samples[j]:
            figP.plot(-uy, rfy, color=col[i], linestyle='dashed', label='_nolegend_')
            if rfy[-1] > 20:
                figP.scatter(-uy[-1], rfy[-1], color='k', marker='x', label='_nolegend_')
        else:
            print('\n . . . Invalid file!\n')

# figP.plot(-1E5, -1E5, color='k', label='Tested side')
# figP.plot(-1E5, -1E5, color='k', linestyle='dashed', label='Collateral side')
# figP.plot(-1E5, -1E5, color=col[0], label='Screw-excess = 5 mm')
# figP.plot(-1E5, -1E5, color=col[1], label='Screw-excess = 0 mm')
# figT.plot(-1E5, -1E5, color='k', label='Tested side')
# figT.plot(-1E5, -1E5, color='k', linestyle='dashed', label='Collateral side')
# figT.plot(-1E5, -1E5, color=col[0], label='Screw-excess = 5 mm')
# figT.plot(-1E5, -1E5, color=col[1], label='Screw-excess = 0 mm')


#%%
fileAramis = '/home/biomech/Documents/01_Icotec/01_Experiments/00_Data/05_Pilot5/ICOTEC_S130684_L4_aramis.csv'
fileAcumen = '/home/biomech/Documents/01_Icotec/01_Experiments/00_Data/05_Pilot5/ICOTEC_S130684_L4_accumen.csv'
# df = read_ARAMIS(fileAramis)
[x, y, z, rX, rY, rZ, t] = read_ARAMIS(fileAramis)
[C, D, F, P, V, T] = read_acumen(fileAcumen)

'''plt.figure()
plt.title('Displacement')
plt.plot(x)
plt.plot(y)
plt.plot(z)
plt.legend(['X', 'Y', 'Z'])

plt.figure()
plt.title('Angles')
plt.plot(rX)
plt.plot(rY)
plt.plot(rZ)
plt.legend(['rX', 'rY', 'rZ'])'''

y = y.to_numpy().flatten()
valls, _ = find_peaks(y, distance=30)
peaks, _ = find_peaks(-y, distance=30)

'''plt.figure()
plt.plot(y[valls])
plt.plot(y[peaks])'''

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
for s in range(int(t[0]), int(t[0]+100)):  # int(np.max(t)-1)):
    arr = np.where(t.astype(int) == s)[0]
    print('\nNew:')
    print(y[np.argmin(y[arr])])  # here: problem is that argmin gives only relative index of value, but absolute needed!
    print('should be smaller than')
    print(y[peak[-1]])
    #print(y[np.argmax(y[arr])])
    #print(y[vall[-1]])
    if y[np.argmin(y[arr])] < y[peak[-1]]: #and y[np.argmax(y[arr])] < y[vall[-1]]:
        print('\n\n\nok at ' + str(s) + '\n\n\n')
        peak.append(arr[int(np.argmin(y[arr]))])
        vall.append(arr[int(np.argmax(y[arr]))])
    # else:
    #     peak.append(arr[int(np.argmin(y[arr]))])
    #     vall.append(arr[int(np.argmax(y[arr]))])

'''
start = 17
peakAc = [0]
vallAc = [0]
for s in range(int(T[0]), int(T[0]+5)):  # int(np.max(T))):
    arrAc = np.where(T.astype(int) == s)[0]
    print('\nNew:')
    print(np.argmin(D[arrAc]))
    print(D[peakAc[-1]])
    print(np.argmax(D[arrAc]))
    print(D[vallAc[-1]])
    if np.argmin(D[arrAc]) < D[peakAc[-1]] and np.argmax(D[arrAc]) > D[vallAc[-1]]:
        peakAc.append(arrAc[int(np.argmin(D[arrAc]))])
        # if np.argmax(D[arrAc]) > D[vallAc[-1]]:
        vallAc.append(arrAc[int(np.argmax(D[arrAc]))])
    # else:
    #     peakAc.append(arrAc[int(np.argmin(D[arrAc]))])
    #     vallAc.append(arrAc[int(np.argmax(D[arrAc]))])

plt.plot(y[peak[start:]], c='r', label='ARAMIS')
plt.plot(y[vall[start:]], c='r')
plt.plot(D[peakAc], c='b', label='ACUMEN')
plt.plot(D[vallAc], c='b')
plt.plot(F[peakAc], c='k', label='Force')
plt.plot(F[vallAc], c='k')
plt.legend()

plt.figure()
plt.plot(y[peak[start:]], F[peakAc])
plt.plot(uy, rfy)

figP.scatter(-y[peak[start:]], -F[peakAc])
'''