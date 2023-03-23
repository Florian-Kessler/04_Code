import numpy as np
import matplotlib.pyplot as plt
import time
import pandas as pd
import os
from scipy.signal import resample


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

'''fig1, figP = plt.subplots(1, 1, figsize=(9, 6))
plt.title('PEEK (YM = 15 GPa)')
fig2, figT = plt.subplots(1, 1, figsize=(9, 6))
plt.title('Ti (YM = 100 GPa)')
col = ['#0072BD', '#D95319', '#EDB120', '#7E2F8E', '#77AC30', '#4DBEEE', '#A2142F', '#A214CC', '#A2DD2F']'''

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
plt.figure()

i = 1
cut = 553

print(samplesD[i])
print(samplesF[i])
[x, y, z, rX, rY, rZ, t] = read_ARAMIS(loc + folder + samplesD[i])
x = x.to_numpy().flatten()[cut:-(9320-252)]
y = y.to_numpy().flatten()[cut:-(9320-252)]
z = z.to_numpy().flatten()[cut:-(9320-252)]
rX = rX.to_numpy().flatten()[cut:-(9320-252)]
rY = rY.to_numpy().flatten()[cut:-(9320-252)]
rZ = rZ.to_numpy().flatten()[cut:-(9320-252)]
t = np.array(t).flatten()
plt.plot(y)
for j in range(len(t)):
    hhmmss = t[j].split(' ')[1]
    hh = hhmmss.split(':')[0]
    mm = hhmmss.split(':')[1]
    ss = hhmmss.split(':')[2].split(',')[0]
    fr = hhmmss.split(':')[2].split(',')[1]
    t[j] = int(hh) * 3600 + int(mm) * 60 + int(ss) + int(fr) / 1000

[C, D, F, _, _, T] = read_acumen(loc + folder + samplesF[i])
Ds = resample(D, int(len(D)/128*30))#[:len(y)]
Ds = Ds.reshape(len(Ds),)
Fs = resample(F, int(len(F)/128*30))#[:len(y)]
Fs = Fs.reshape(len(Fs),)
Cs = resample(C, int(len(D)/128*30))#[:len(y)]
Cs = Cs.reshape(len(Cs),)
plt.plot(Ds)

if not len(y) - len(Ds):
    print('\nEqual length.')
    df = pd.DataFrame({'Aramis X': x,
                       'Aramis Y': y,
                       'Aramis Z': z,
                       'Aramis rX': rX,
                       'Aramis rY': rY,
                       'Aramis rZ': rZ,
                       'Acumen Y': Ds,
                       'Acumen Fy': Fs,
                       'Acumen C': Cs})
    df.to_csv(loc + folder + samplesD[i].split('_a')[0] + '_resample.csv')
else:
    print(len(y))
    print(len(Ds))
