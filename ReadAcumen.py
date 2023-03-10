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


def read_acumen(file_a):
    df_ = pd.read_csv(file_a, delimiter=';', skiprows=[0, 2])
    # t_ = pd.DataFrame(df_, columns=['Time ']).to_numpy()
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

    return cycle_, d_, f_, peak_, vall_


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx


def find_first(array, value):
    array = np.asarray(array)
    idx = next(x for x, val in enumerate(array)
               if val <= value)
    return idx


t1 = time.time()
plt.close('all')

fig, ax1 = plt.subplots(1, 1, figsize=(9, 6))
ax2 = ax1.twinx()

specimen = '05_Pilot5'

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

fig1, figP = plt.subplots(1, 1, figsize=(9, 6))
plt.title('PEEK screw')
# fig2, figT = plt.subplots(1, 1, figsize=(9, 6))
# plt.title('Ti screw')
col = ['#0072BD', '#D95319', '#EDB120', '#7E2F8E', '#77AC30', '#4DBEEE', '#A2142F', '#A214CC', '#A2DD2F']

# number = ['00', '01', '02', '10', '11', '12']
number = ['50']  # , '55']
for i in range(len(number)):
    loc = '/home/biomech/Documents/01_Icotec/02_FEA/98_Pilots/' + specimen + '/'
    folder = [filename for filename in os.listdir(loc) if filename.startswith(number[i])][0] + '/'
    samples = [filename for filename in os.listdir(loc + folder + '/') if filename.endswith('RFnode.txt')
               and '_02_' in filename]
    print(samples)
    [uy, rfy] = 2*[0]
    screw_force = np.zeros([5, 21])
    ang = np.zeros([5, 21])
    Screw_mat = ['P', 'T']
    Sim_mat = ['P', 'T']
    Fpoint = (-np.array([0, 75, 150])) * -37.037 - 1851
    for j in range(len(samples)):
        # lab = 'Screw excess = ' + folder.split('_S')[-1][0] + '.' + folder.split('_S')[-1][1] + ' mm, ' + \
        #       'diameter = ' + folder.split('_D')[-1][0] + '.' + folder.split('_D')[-1][1] + ' mm'
        lab = 'Unidirectional'
        if '55' in number[i]:
            lab = 'Toggling'
        file_path = loc + folder + samples[j]
        [uy, r_] = read_RFnodeFile(file_path)
        [u_, rfy] = read_RFnodeFile(file_path.split('.txt')[0] + 'Fix.txt')
        if 'P_P' in samples[j]:
            figP.plot(-uy, rfy, color=col[i], linestyle='solid', label=lab)
            if rfy[-1] > 15:
                figP.scatter(-uy[-1], rfy[-1], color='k', marker='x', label='_nolegend_')
        elif 'T_T' in samples[j]:
            figP.plot(-uy, rfy, color=col[i], linestyle='solid', label=lab)  # figT
            if rfy[-1] > 15:
                figP.scatter(-uy[-1], rfy[-1], color='k', marker='x', label='_nolegend_')  # figT
        elif 'P_T' in samples[j]:
            figP.plot(-uy, rfy, color=col[i], linestyle='dashdot', label='_nolegend_')  # figT
            if rfy[-1] > 15:
                figP.scatter(-uy[-1], rfy[-1], color='k', marker='x', label='_nolegend_')  # figT
        elif 'T_P' in samples[j]:
            figP.plot(-uy, rfy, color=col[i], linestyle='dashdot', label='_nolegend_')
            if rfy[-1] > 15:
                figP.scatter(-uy[-1], rfy[-1], color='k', marker='x', label='_nolegend_')
        else:
            print('\n . . . Invalid file!\n')

# figP.plot(-1E5, -1E5, color='k', label='Tested side')
# figP.plot(-1E5, -1E5, color='k', linestyle='dashdot', label='Collateral side')
# figT.plot(-1E5, -1E5, color='k', label='Tested side')
# figT.plot(-1E5, -1E5, color='k', linestyle='dashdot', label='Collateral side')

if cyc2P:
    figP.scatter(-exp['dPEEK'][exp['peakPEEK'][:cyc2P+1000]],
                 -exp['fPEEK'][exp['peakPEEK'][:cyc2P+1000]], color='k', s=0.5, label='Experiment')
else:
    figP.scatter(-exp['dPEEK'][exp['peakPEEK']],
                 -exp['fPEEK'][exp['peakPEEK']], color='k', s=0.5, label='Experiment')
if cyc2T:
    figP.scatter(-exp['dTi'][exp['peakTi'][:cyc2T+1000]],  # figT, 'k'
                 -exp['fTi'][exp['peakTi'][:cyc2T+1000]], color='r', s=0.5, label='Experiment, corr.: '+str(cor)+' mm')
else:
    figP.scatter(-exp['dTi'][exp['peakTi']],  # figT, 'k'
                 -exp['fTi'][exp['peakTi']], color='r', s=0.5, label='Experiment, corr.: ' + str(cor) + ' mm')

figP.axis([0, 15, 0, 180])  # 30
figP.axis([0, 15, 0, 180])  # figP
figP.set_xlabel('Displacement / mm')
figP.set_ylabel('Force / N')
figP.set_xlabel('Displacement / mm')  # figP
figP.set_ylabel('Force / N')  # figP
figP.legend(loc='lower right')
figP.legend(loc='lower right')  # figP

'''
# temp #
fig5, ax5 = plt.subplots(1, 1, figsize=(9, 6))
file = '/home/biomech/Documents/01_Icotec/02_FEA/98_Pilots/03_Pilot3/50_L50_S00_D30/E_PEEK.txt'
df = np.loadtxt(file, delimiter='\t')
timeX = np.array(df[:, 0])
strain = np.array(df[:, 1])
plt.plot(timeX[40:]-1, strain[40:])
plt.plot([0, 1], [0.01, 0.01], linestyle='dashed', c='k')
plt.xlabel('Step time')
plt.ylabel('Logarithmic strain')

print('\nRuntime: ' + str(round(time.time() - t1, 2)) + ' seconds.')
'''
