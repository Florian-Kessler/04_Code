import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import time
import pandas as pd
import scipy


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


def read_energy(file_):
    # read data from text file
    df_ = np.loadtxt(file_)#, delimiter='   ')
    t_ = np.array(df_[:, 0])
    e_ = np.array(df_[:, 1])
    return t_, e_


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
    box = np.ones(box_pts) / box_pts
    y_smooth = np.convolve(y_, box, mode='same')
    return y_smooth


t1 = time.time()
# plt.close('all')
# mue = ['07', '05', '03', '02', '01', '00']

# # # # # INPUT # # # # #
loc = '/home/biomech/Documents/01_Icotec/01_Experiments/00_Data/01_MainStudy/'
specimen_names = open('/home/biomech/Documents/01_Icotec/Specimens.txt', 'r').read().splitlines()  # Read specimens
col = ['#0072BD', '#D95319', '#EDB120', '#7E2F8E', '#77AC30', '#4DBEEE', '#A2142F', '#A214CC', '#A2DD2F']
number = ['75']  # simulations


#%%
plt.close('all')

peek_samples = [2, 5, 7, 8, 10, 13, 15, 16, 18, 21, 23, 26, 29, 31, 32]  # without 24
ti_samples = [3, 4, 6, 9, 11, 12, 14, 17, 19, 20, 22, 27, 28, 30, 33]  # without 25
slope = np.zeros(34)
slope2 = np.zeros(34)
slopeFE02 = np.zeros(34)
slopeFE02_2 = np.zeros(34)
slopeFE05 = np.zeros(34)
slopeFE05_2 = np.zeros(34)
cycle = 2
plots = 1
noData02 = []
noData05 = []

for i in [9]:  # ti_samples:  # range(2, 34):
    specimen = specimen_names[i]  # 'S131318_L4_right'
    uy = 0
    rfy = 0
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

    sample = sample.split('_resample')[0].split('/')[-1]
    file_path = '/home/biomech/Documents/01_Icotec/02_FEA/01_MainStudy/' + sample + '/87_L50_S50_D45' + \
                '/87_L50_S50_D45_d1_02_T_RFnode.txt'
    if len(open(file_path, 'r').readline()) > 0:
        [uy, rf_] = read_RFnodeFile(file_path)
        [u_, rfy] = read_RFnodeFile(file_path.split('.txt')[0] + 'Fix.txt')
        if len(uy) == 42:
            slopeFE02[i] = -(rfy[-11] - rfy[-6]) / (uy[-11] - uy[-6])
            slopeFE02_2[i] = -(rfy[-11] - rfy[-1]) / (uy[-11] - uy[-1])
            print('Exp.: ' + str(np.round(slope[i], 1)) + ' N/mm')
            print('Exp2: ' + str(np.round(slope2[i], 1)) + ' N/mm')
            print('FE:   ' + str(np.round(slopeFE02[i], 1)) + ' N/mm')
            print('FE2:  ' + str(np.round(slopeFE02_2[i], 1)) + ' N/mm\n')
        else:
            print('Data missing for ' + str(sample) + '.\n')
            noData02.append(i)
    else:
        print('Data missing for ' + str(sample) + '.\n')
        noData02.append(i)
    if plots:
        plt.figure()
        plt.plot(uy, -rfy)
        plt.scatter(uy[-6], -rfy[-6], color=col[1])
        plt.scatter(uy[-1], -rfy[-1], color=col[1])
        plt.scatter(uy[-11], -rfy[-11], color=col[2])
        plt.plot(uy[[-11, -6]], -rfy[[-11, -6]], 'r--')
        plt.plot(uy[[-11, -1]], -rfy[[-11, -1]], 'g--')
        plt.xlabel('Displacement / mm')
        plt.ylabel('Force / N')
    file_path = '/home/biomech/Documents/01_Icotec/02_FEA/01_MainStudy/' + sample + '/87_L50_S50_D45' + \
                '/87_L50_S50_D45_d1_05_T_RFnode.txt'
    if len(open(file_path, 'r').readline()) > 0:
        [uy, rf_] = read_RFnodeFile(file_path)
        [u_, rfy] = read_RFnodeFile(file_path.split('.txt')[0] + 'Fix.txt')
        if len(uy) == 42:
            slopeFE05[i] = -(rfy[-11] - rfy[-6]) / (uy[-11] - uy[-6])
            slopeFE05_2[i] = -(rfy[-11] - rfy[-1]) / (uy[-11] - uy[-1])
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
        plt.plot(uy[[-11, -6]], -rfy[[-11, -6]], 'r--')
        plt.plot(uy[[-11, -1]], -rfy[[-11, -1]], 'g--')
        plt.xlabel('Displacement / mm')
        plt.ylabel('Force / N')



plt.figure()

plt.scatter(slope[ti_samples], slopeFE02[ti_samples], color=col[0], label='$\mu$ = 0.2')
plt.scatter(slope2[ti_samples], slopeFE02_2[ti_samples], color=col[0], marker='x', label='$\mu$ = 0.2 (extrema)')
plt.scatter(slope[[3, 4, 6, 9, 11, 12, 14, 22, 27, 28, 30, 33]],
            slopeFE05[[3, 4, 6, 9, 11, 12, 14, 22, 27, 28, 30, 33]], color=col[1], label='$\mu$ = 0.5', alpha=0.7)
plt.scatter(slope2[[3, 4, 6, 9, 11, 12, 14, 22, 27, 28, 30, 33]],
            slopeFE05_2[[3, 4, 6, 9, 11, 12, 14, 22, 27, 28, 30, 33]],
            color=col[1], marker='x', label='$\mu$ = 0.5 (extrema)', alpha=0.7)
plt.plot([0, 100], [0, 100], 'k')
plt.xlabel('Stiffness Experiment / N/mm')
plt.ylabel('Stiffness FE / N/mm')
plt.legend()

# plt.figure()
# plt.plot(slope[peek_samples] / slope[ti_samples])

# plt.figure()
# plt.scatter(slope[peek_samples], slope[ti_samples], color='k')
# plt.xlim([0, 80])
# plt.ylim([0, 80])


'''

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
    print('Execution time: ' + str(round(tRun, 1)) + ' sec.')'''

'''
fric = '02'
sample = sample.split('_resample')[0].split('/')[-1]
file_path = '/home/biomech/Documents/01_Icotec/02_FEA/01_MainStudy/' + sample + '/87_L50_S50_D45' + 
            '/87_L50_S50_D45_d1_' + fric + '_T_RFnode.txt'
[uy, rf_] = read_RFnodeFile(file_path)
[u_, rfy] = read_RFnodeFile(file_path.split('.txt')[0] + 'Fix.txt')
plt.figure()
plt.plot(uy, -rfy)
plt.plot(uy[-11:-5], -rfy[-11:-5])
slopeFE = -(rfy[-11]-rfy[-6])/(uy[-11]-uy[-6])
print(slopeFE)
'''
#%%
file = '/home/biomech/DATA/01_Icotec/02_FEA/02_uFE/Tests/test_7_RFnode.txt'
[uy, rf_] = read_RFnodeFile(file)
[u_, rfy] = read_RFnodeFile(file.split('.txt')[0] + 'Fix.txt')

sample = loc + specimen_names[8] + '_resample.csv'
[ArX, ArY, ArZ, ArrX, ArrY, ArrZ, AcY, AcFy, AcC] = read_resample(sample)

file = '/home/biomech/DATA/01_Icotec/02_FEA/02_uFE/Tests/test_8_RFnode.txt'
[uy2, rf_] = read_RFnodeFile(file)
[u_, rfy2] = read_RFnodeFile(file.split('.txt')[0] + 'Fix.txt')

plt.figure()
plt.plot(uy, -rfy, label='Ti')

plt.plot(AcY, AcFy-AcFy['Acumen Fy'][0], label='Experiment')
plt.plot(uy2, -rfy2, label='PEEK')
plt.scatter(uy[-1], -rfy[-1], color='r', marker='x')
plt.scatter(uy2[-1], -rfy2[-1], color='r', marker='x')
plt.title('uFE tests')
plt.xlabel('Displacement / mm')
plt.ylabel('Force / N')

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
radius = [4, 6]
# fig, axs = plt.subplots(1, 1, figsize=(9, 6))
for j in range(len(radius)):
    bvtvList = []
    specList = []
    indexList = []
    for i in range(len(specimen_names)):
        loc_ = '/home/biomech/DATA/01_Icotec/01_Experiments/02_Scans/BVTV/BVTV_along_'
        bvtv = np.load(loc_ + specimen_names[i] + '_' + str(radius[j]) + 'mm.npy')
        # axs.plot(bvtv)
        # test line plot seaborn
        bvtvList += bvtv.T.tolist()
        specList += len(bvtv)*[i]
        indexList += range(len(bvtv))
    dfPlot = pd.DataFrame()
    dfPlot['bvtv'] = bvtvList
    dfPlot['specimen'] = specList
    dfPlot['index'] = indexList
    sns.lineplot(data=dfPlot, x='index', y='bvtv', errorbar='sd', label='Radius: ' + str(radius[j]) + ' mm')
#%%
plt.figure()
radius = [4, 45, 5, 6]
for i in range(len(specimen_names)):
    loc_ = '/home/biomech/DATA/01_Icotec/01_Experiments/02_Scans/BVTV/BVTV_along_'
    bvtv = np.load(loc_ + specimen_names[i] + '_' + str(radius[0]) + 'mm.npy')

    sample = loc + specimen_names[i] + '_resample.csv'
    [ArX, ArY, ArZ, ArrX, ArrY, ArrZ, AcY, AcFy, AcC] = read_resample(sample)
    plt.scatter(np.max(AcFy, axis=0), np.mean(bvtv[:200], axis=0))
