import numpy as np
import matplotlib.pyplot as plt
import time
import pandas as pd
from scipy.signal import argrelextrema


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
    f_ = pd.DataFrame(df, columns=['Axial Force ']).to_numpy()
    # f_set_ = pd.DataFrame(df, columns=['Axial Force Command ']).to_numpy()
    cycle_ = pd.DataFrame(df, columns=['Axial Count ']).to_numpy()
    # arr_ = 0
    peak_ = np.zeros(int(np.max(cycle_)))
    vall_ = np.zeros(int(np.max(cycle_)))
    for j in range(2, int(np.max(cycle_))):
        # del arr_
        arr_ = np.where((cycle_ == j) | (cycle_ == j + .5))[0]
        peak_[j] = arr_[int(np.argmin(f_[arr_]))]
        vall_[j] = arr_[int(np.argmax(f_[arr_]))]

    peak_ = peak_.astype(int)
    vall_ = vall_.astype(int)

    return cycle_, d_, f_, peak_, vall_


t1 = time.time()
plt.close('all')


fig, ax1 = plt.subplots(1)
ax2 = ax1.twinx()

PEEK1 = 'Pilot1/S131840_L4_S1_PEEK.csv'
TITAN1 = 'Pilot1/S191840_L4_S2_DPS.csv'
PEEK3 = 'Pilot3/ICOTEC_S130672_L5_icotec_accumen.csv'
TITAN3 = 'Pilot3/ICOTEC_S130672_L5_DPS_accumen.csv'
loc = '/home/biomech/Documents/01_Icotec/01_Experiments/00_Data/'

file = [loc + PEEK3, loc + TITAN3]

col = [['#1f6a06', '#0a1eb2'],
       ['#373737', '#b20a0a']]
labels = ['ICOTEC', 'Titanium']
cycle = 0
for i in range(0, len(file)):
    [cycle, d, f, peak, vall] = read_acumen(str(file[i]))
    ax1.plot(cycle[peak], d[peak], color=col[i][0], label='_nolegend_')
    ax1.plot(cycle[vall], d[vall], color=col[i][0], alpha=0.75, label='_nolegend_')
    ax2.plot(cycle[peak], f[peak], color=col[i][1], label='Force ' + str(labels[i]))
    ax2.plot(cycle[vall], f[vall], color=col[i][1], alpha=0.75, label='_nolegend_')
    ax2.plot([-.2, .2], color=col[i][0], label='Displacement ' + str(labels[i]))

ax2.spines.right.set_visible(True)
ax1.set_xlabel('Cycle Number')
ax1.set_ylabel('Displacement / mm')
ax2.set_ylabel('Force / N')
ax1.axis([0, int(np.ceil(np.max(cycle)/1000)*1000),  -20, 0])
ax2.axis([0, int(np.ceil(np.max(cycle)/1000)*1000), -750, 0])

plt.legend(loc='lower right')


loc = '/home/biomech/Documents/01_Icotec/02_FEA/99_Tests/Pilot3/archive/'
sample = '00_Pilot3'
[uy, rfy] = 2*[0]
col = ['#0072BD', '#D95319', '#EDB120', '#7E2F8E', '#77AC30', '#4DBEEE', '#A2142F', '#A214CC', '#A2DD2F']
screw_force = np.zeros([5, 21])
ang = np.zeros([5, 21])
Screw_mat = ['T', 'P']
Sim_mat = ['T', 'P']
plt.figure()
Fpoint = (-np.array([0, 50, 75, 150, 185, 225, 300])) * (-4950/185)-50
for i in range(len(Screw_mat)):
    for j in range(len(Sim_mat)):
        file = sample + Screw_mat[i] + '_' + Sim_mat[j]
        del uy, rfy
        pathRF = loc + str(file) + '_RFnode.txt'
        pathRFfix = loc + str(file) + '_RFnodeFix.txt'
        [uy, r_] = read_RFnodeFile(pathRF)
        [u_, rfy] = read_RFnodeFile(pathRFfix)
        plt.plot(uy, rfy)

ax2.scatter(Fpoint, -np.array([0, 50, 75, 150, 185, 225, 300]))
'''
        Fpmax = np.append(argrelextrema(rfy, np.greater), len(rfy) - 1)
        Fpmin = np.append(0, argrelextrema(rfy, np.less))
        Fpmax_cycles = rfy[Fpmax] * 20 - 1000
        Fpmin_cycles = rfy[Fpmax] * 20 - 1000
        if '_P' in file:
            c = i
            ax1.scatter(F_point, uy[Fpmax], color=col[c], marker='^', alpha=0.5)
            #ax1.plot(F_point, uy[Fpmax], color=col[c], alpha=0.5)
            ax2.scatter(-1E9, 1E9, color=col[c], marker='^')
            ax1.scatter(F_point, uy[Fpmin], color=col[c], marker='^')
            ax1.plot(F_point, uy[Fpmin], color=col[c])
        elif '_T' in file:
            c = i
            ax1.scatter(F_point, uy[Fpmax], color=col[c], marker='o', alpha=0.5)
            #ax1.plot(F_point, uy[Fpmax], color=col[c], alpha=0.5, ls=':')
            ax2.scatter(-1E9, 1E9, color=col[c], marker='o')
            ax1.scatter(F_point, uy[Fpmin], color=col[c], marker='o')
            ax1.plot(F_point, uy[Fpmin], color=col[c], ls=':')
'''

print('\nRuntime: ' + str(round(time.time() - t1, 2)) + ' seconds.')
