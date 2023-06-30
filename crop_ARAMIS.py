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
    df_ = pd.read_csv(file_a, delimiter='\t', skiprows=[0, 1, 2, 3, 5])
    t_ = pd.DataFrame(df_, columns=['Time ']).to_numpy()
    d_ = pd.DataFrame(df_, columns=['Axial Displacement ']).to_numpy()
    # d_ = d_ - d_[0]  # calibrate displacement to zero at beginning
    f_ = pd.DataFrame(df_, columns=['Axial Force ']).to_numpy()
    # f_set_ = pd.DataFrame(df_, columns=['Axial Force Command ']).to_numpy()
    cycle_ = pd.DataFrame(df_, columns=['Axial Count ']).to_numpy()
    # arr_ = 0
    # peak_ = np.zeros(int(np.max(cycle_)))
    # vall_ = np.zeros(int(np.max(cycle_)))
    # for j_ in range(2, int(np.max(cycle_))):
    #     # del arr_
    #     arr_ = np.where((cycle_ == j_) | (cycle_ == j_ + .5))[0]
    #     peak_[j_] = arr_[int(np.argmin(f_[arr_]))]
    #     vall_[j_] = arr_[int(np.argmax(f_[arr_]))]
    # print(len(d_))
    # peak_ = peak_.astype(int)
    # vall_ = vall_.astype(int)

    # return cycle_, d_, f_, peak_, vall_, t_
    return cycle_, d_, f_, t_

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
loc = '/home/biomech/Documents/01_Icotec/01_Experiments/00_Data/01_MainStudy/'
specimens = open('/home/biomech/Documents/01_Icotec/Specimens.txt', 'r').read().splitlines()
# # # # # INPUT # # # # #
specimen = specimens[12]

# for i in range(34):
#     read_file = pd.read_csv(r'/home/biomech/Documents/01_Icotec/01_Experiments/00_Data/01_MainStudy/' +
#     specimens[i] + '_acumen.txt')
#     read_file.to_csv(r'/home/biomech/Documents/01_Icotec/01_Experiments/00_Data/01_MainStudy/' +
#     specimens[i].split('_a')[0] + '_acumen.csv', index=None)

file_acumen = [filename for filename in os.listdir(loc) if filename.startswith(specimen)
               and filename.endswith('acumen.csv')][0]
file_aramis = [filename for filename in os.listdir(loc) if filename.startswith(specimen)
               and filename.endswith('aramis.csv')][0]
print('Aramis: ' + file_aramis)
[x, y, z, rX, rY, rZ, t] = read_ARAMIS(loc + file_aramis)


'''

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


print(samplesD[i])
print(samplesF[i])

[x, y, z, rX, rY, rZ, t] = read_ARAMIS(loc + folder + samplesD[i])
'''
# # # # # INPUT # # # # #
startBlue = 0
cutBlue = 1  # >0
a = []
# a = np.arange(21073, 21551)
# a = [len(y)-1-startBlue-cutBlue]

x = x.to_numpy().flatten()[startBlue:-cutBlue]
y = y.to_numpy().flatten()[startBlue:-cutBlue]
z = z.to_numpy().flatten()[startBlue:-cutBlue]
rX = rX.to_numpy().flatten()[startBlue:-cutBlue]
rY = rY.to_numpy().flatten()[startBlue:-cutBlue]
rZ = rZ.to_numpy().flatten()[startBlue:-cutBlue]

if a:
    x = np.delete(x, a)
    y = np.delete(y, a)
    z = np.delete(z, a)
    rX = np.delete(rX, a)
    rY = np.delete(rY, a)
    rZ = np.delete(rZ, a)
plt.plot(y, color='b', label='ARAMIS')
t = np.array(t).flatten()
for j in range(len(t)):
    hhmmss = t[j].split(' ')[1]
    hh = hhmmss.split(':')[0]
    mm = hhmmss.split(':')[1]
    ss = hhmmss.split(':')[2].split(',')[0]
    fr = hhmmss.split(':')[2].split(',')[1]
    t[j] = int(hh) * 3600 + int(mm) * 60 + int(ss) + int(fr) / 1000

[C, D, F, T] = read_acumen(loc + file_acumen)

# # # # # INPUT # # # # #
startRed = 30
cutRed = 1  # >0
corr = 0  # correct amplitude for overlay

res = 10

Ds = resample(D, int(len(D)/128*res))[startRed:-cutRed]
Ds = Ds.reshape(len(Ds),)
Fs = resample(F, int(len(F)/128*res))[startRed:-cutRed]
Fs = Fs.reshape(len(Fs),)
Cs = resample(C, int(len(D)/128*res))[startRed:-cutRed]
Cs = Cs.reshape(len(Cs),)
print(Ds)
plt.plot(Ds-corr, color='r', label='ACUMEN')
plt.legend()

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
    df.to_csv(loc + file_aramis.split('_a')[0] + '_resample.csv')
    print('resample.csv written.')
else:
    print('Length Aramis (blue): ' + str(len(y)))
    print('Length Acumen (red):  ' + str(len(Ds)))
    print('Cut blue: ' + str(cutBlue + len(y)-len(Ds)))
