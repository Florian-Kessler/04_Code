import numpy as np
import matplotlib.pyplot as plt
import time
import pandas as pd
import os
# from scipy.signal import argrelextrema


def read_RFnodeFile(file_):
    # read data from text file
    df = np.loadtxt(file_, delimiter=',')

    # Contains the frame number:
    # frame = np.array(df[:, 0])

    # reaction moment of the implant:
    # rmx_ = np.array(df[:, 1])
    # rmy_ = np.array(df[:, 2])
    # rmz_ = np.array(df[:, 3])

    # reaction force of the implant:
    # rfx_ = np.array(df[:, 4])
    rfy_ = np.array(df[:, 5])
    # rfz_ = np.array(df[:, 6])

    # transverse disp. of the implant:
    # ux_ = np.array(df[:, 7])
    uy_ = np.array(df[:, 8])
    # uz_ = np.array(df[:, 9])

    # rotational disp. of the implant:
    # urx_ = np.array(df[:, 10])
    # ury_ = np.array(df[:, 11])
    # urz_ = np.array(df[:, 12])

    return uy_, rfy_


def read_acumen(file_a):
    df = pd.read_csv(file_a, delimiter=';', skiprows=[0, 2])
    # t_ = pd.DataFrame(df, columns=['Time ']).to_numpy()
    d_ = pd.DataFrame(df, columns=['Axial Displacement ']).to_numpy()
    # d_ = d_ - d_[0]  # calibrate displacement to zero at beginning
    f_ = pd.DataFrame(df, columns=['Axial Force ']).to_numpy()
    # f_set_ = pd.DataFrame(df, columns=['Axial Force Command ']).to_numpy()
    cycle_ = pd.DataFrame(df, columns=['Axial Count ']).to_numpy()
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


# To do: 10_...TP RF node files
t1 = time.time()
plt.close('all')

fig, ax1 = plt.subplots(1, 1, figsize=(9, 6))
ax2 = ax1.twinx()

PEEK1 = 'Pilot1/S131840_L4_S1_PEEK.csv'
TITAN1 = 'Pilot1/S191840_L4_S2_DPS.csv'
PEEK3 = 'Pilot3/ICOTEC_S130672_L5_icotec_accumen.csv'
TITAN3 = 'Pilot3/ICOTEC_S130672_L5_DPS_accumen.csv'
loc = '/home/biomech/Documents/01_Icotec/01_Experiments/00_Data/'

file = [loc + PEEK3, loc + TITAN3]
exp = {}
col = [['#1f6a06', '#0a1eb2'],
       ['#373737', '#b20a0a']]
labels = ['PEEK', 'Ti']
for i in range(0, len(file)):
    [exp['cycle' + labels[i]], exp['d' + labels[i]], exp['f' + labels[i]], exp['peak' + labels[i]],
     exp['vall' + labels[i]]] = read_acumen(str(file[i]))
    ax1.plot(exp['cycle' + labels[i]][exp['peak' + labels[i]]], exp['d' + labels[i]][exp['peak' + labels[i]]],
             color=col[i][0], label='_nolegend_')
    ax1.plot(exp['cycle' + labels[i]][exp['vall' + labels[i]]], exp['d' + labels[i]][exp['vall' + labels[i]]],
             color=col[i][0], alpha=0.75, label='_nolegend_')
    ax2.plot(exp['cycle' + labels[i]][exp['peak' + labels[i]]], exp['f' + labels[i]][exp['peak' + labels[i]]],
             color=col[i][1], label='Force ' + str(labels[i]))
    ax2.plot(exp['cycle' + labels[i]][exp['vall' + labels[i]]], exp['f' + labels[i]][exp['vall' + labels[i]]],
             color=col[i][1], alpha=0.75, label='_nolegend_')
    ax2.plot([-.2, .2], color=col[i][0], label='Displacement ' + str(labels[i]))

ax2.spines.right.set_visible(True)
ax1.set_xlabel('Cycle Number')
ax1.set_ylabel('Displacement / mm')
ax2.set_ylabel('Force / N')
ax1.axis([0, int(np.ceil(np.max(exp['cycle' + labels[1]])/1000)*1000),  -25, 0])
ax2.axis([0, int(np.ceil(np.max(exp['cycle' + labels[1]])/1000)*1000), -750, 0])
ax2.legend(['PEEK Force', 'PEEK Displacement', 'Ti Force', 'Ti Displacement'])

cyc1T = find_first(exp['fTi'][exp['peakTi']], -75)
cyc2T = find_first(exp['fTi'][exp['peakTi']], -150)
cyc1P = find_first(exp['fPEEK'][exp['peakPEEK'][:1330]], -75)
cyc2P = find_first(exp['fPEEK'][exp['peakPEEK'][:1330]], -150)


fig1, figP = plt.subplots(1, 1, figsize=(9, 6))
plt.title('PEEK screw')
fig2, figT = plt.subplots(1, 1, figsize=(9, 6))
plt.title('Ti screw')
col = ['#0072BD', '#D95319', '#EDB120', '#7E2F8E', '#77AC30', '#4DBEEE', '#A2142F', '#A214CC', '#A2DD2F']

number = ['00', '01', '02', '10', '11', '12']
for i in range(len(number)):
    loc = '/home/biomech/Documents/01_Icotec/02_FEA/98_Pilots/03_Pilot3/'
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
        lab = 'Screw excess = ' + folder.split('_S')[-1][0] + '.' + folder.split('_S')[-1][1] + ' mm, ' + \
              'diameter = ' + folder.split('_D')[-1][0] + '.' + folder.split('_D')[-1][1] + ' mm'
        file_path = loc + folder + samples[j]
        [uy, r_] = read_RFnodeFile(file_path)
        [u_, rfy] = read_RFnodeFile(file_path.split('.txt')[0] + 'Fix.txt')
        if 'P_P' in samples[j]:
            figP.plot(-uy, rfy, color=col[i], linestyle='solid', label=lab)
            if rfy[-1] > 30:
                figP.scatter(-uy[-1], rfy[-1], color='k', marker='x', label='_nolegend_')
        elif 'T_T' in samples[j]:
            figT.plot(-uy, rfy, color=col[i], linestyle='solid', label=lab)
            if rfy[-1] > 30:
                figT.scatter(-uy[-1], rfy[-1], color='k', marker='x', label='_nolegend_')
        elif 'P_T' in samples[j]:
            figT.plot(-uy, rfy, color=col[i], linestyle='dashdot', label='_nolegend_')
            if rfy[-1] > 30:
                figT.scatter(-uy[-1], rfy[-1], color='k', marker='x', label='_nolegend_')
        elif 'T_P' in samples[j]:
            figP.plot(-uy, rfy, color=col[i], linestyle='dashdot', label='_nolegend_')
            if rfy[-1] > 30:
                figP.scatter(-uy[-1], rfy[-1], color='k', marker='x', label='_nolegend_')
        else:
            print('\n . . . Invalid file!\n')

'''number = ['00', '01', '02', '10']  # , '11', '12']
for i in range(len(number)):
    loc = '/home/biomech/Documents/01_Icotec/02_FEA/98_Pilots/03_Pilot3/'
    folder = [filename for filename in os.listdir(loc) if filename.startswith(number[i])][0] + '/'
    samples = [filename for filename in os.listdir(loc + folder + '/') if filename.endswith('RFnode.txt')
               and '_07_' in filename]
    [uy, rfy] = 2*[0]
    screw_force = np.zeros([5, 21])
    ang = np.zeros([5, 21])
    Screw_mat = ['P', 'T']
    Sim_mat = ['P', 'T']
    Fpoint = (-np.array([0, 75, 150])) * -37.037 - 1851
    for j in range(len(samples)):
        lab = 'Screw excess = ' + folder.split('_S')[-1][0] + '.' + folder.split('_S')[-1][1] + ' mm, ' + \
              'diameter = ' + folder.split('_D')[-1][0] + '.' + folder.split('_D')[-1][1] + ' mm'
        file_path = loc + folder + samples[j]
        [uy, r_] = read_RFnodeFile(file_path)
        [u_, rfy] = read_RFnodeFile(file_path.split('.txt')[0] + 'Fix.txt')
        if 'P_P' in samples[j]:
            figP.plot(-uy, rfy, color=col[i], linestyle='dotted', label='_nolegend_', alpha=0.5)
            if rfy[-1] > 30:
                figP.scatter(-uy[-1], rfy[-1], color='k', marker='x', label='_nolegend_', alpha=0.5)
        elif 'T_T' in samples[j]:
            figT.plot(-uy, rfy, color=col[i], linestyle='dotted', label='_nolegend_', alpha=0.5)
            if rfy[-1] > 30:
                figT.scatter(-uy[-1], rfy[-1], color='k', marker='x', label='_nolegend_', alpha=0.5)
        elif 'P_T' in samples[j]:
            figT.plot(-uy, rfy, color=col[i], linestyle='dashdot', label='_nolegend_', alpha=0.5)
            if rfy[-1] > 30:
                figT.scatter(-uy[-1], rfy[-1], color='k', marker='x', label='_nolegend_', alpha=0.5)
        elif 'T_P' in samples[j]:
            figP.plot(-uy, rfy, color=col[i], linestyle='dashdot', label='_nolegend_', alpha=0.5)
            if rfy[-1] > 30:
                figP.scatter(-uy[-1], rfy[-1], color='k', marker='x', label='_nolegend_', alpha=0.5)
        else:
            print('\n . . . Invalid file!\n')'''

figP.plot(-1E5, -1E5, color='k', label='Tested side')
figP.plot(-1E5, -1E5, color='k', linestyle='dashdot', label='Collateral side')
figT.plot(-1E5, -1E5, color='k', label='Tested side')
figT.plot(-1E5, -1E5, color='k', linestyle='dashdot', label='Collateral side')

# figP.plot(-exp['dPEEK'][exp['peakPEEK'][cyc1P]-39:exp['peakPEEK'][cyc1P]+39],
#           -exp['fPEEK'][exp['peakPEEK'][cyc1P]-39:exp['peakPEEK'][cyc1P]+39], color='k', label='_nolegend_')
# figP.plot(-exp['dPEEK'][exp['peakPEEK'][cyc2P]-39:exp['peakPEEK'][cyc2P]+39],
#           -exp['fPEEK'][exp['peakPEEK'][cyc2P]-39:exp['peakPEEK'][cyc2P]+39], color='k', label='_nolegend_')
# figT.plot(-exp['dTi'][exp['peakTi'][cyc1T]-39:exp['peakTi'][cyc1T]+39],
#           -exp['fTi'][exp['peakTi'][cyc1T]-39:exp['peakTi'][cyc1T]+39], color='k', label='_nolegend_')
# figT.plot(-exp['dTi'][exp['peakTi'][cyc2T]-39:exp['peakTi'][cyc2T]+39],
#           -exp['fTi'][exp['peakTi'][cyc2T]-39:exp['peakTi'][cyc2T]+39], color='k', label='_nolegend_')

figP.scatter(-exp['dPEEK'][exp['peakPEEK'][:cyc2P+1000]],
             -exp['fPEEK'][exp['peakPEEK'][:cyc2P+1000]], color='k', s=0.5, label='Experiment')
figT.scatter(-exp['dTi'][exp['peakTi'][:cyc2T+1000]],
             -exp['fTi'][exp['peakTi'][:cyc2T+1000]], color='k', s=0.5, label='Experiment')

figP.axis([0, 30, 0, 180])
figT.axis([0, 15, 0, 180])
figP.set_xlabel('Displacement / mm')
figP.set_ylabel('Force / N')
figT.set_xlabel('Displacement / mm')
figT.set_ylabel('Force / N')
figP.legend(loc='lower right')
figT.legend(loc='lower right')
print('\nRuntime: ' + str(round(time.time() - t1, 2)) + ' seconds.')
