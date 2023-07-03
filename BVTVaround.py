# import mapping_noRot as mappNR
import numpy as np
import matplotlib.pyplot as plt
# import matplotlib.cm as cm
import SimpleITK as sitk
import time
import os
import ReadRawMHD as rR
import pandas as pd
from scipy.signal import butter, filtfilt, find_peaks


t0 = time.time()


def butter_lowpass_filter(data_, cutoff_, order=9):
    fs = 10  # sample rate, Hz
    nyq = 0.5 * fs
    normal_cutoff = cutoff_ / nyq
    # Get the filter coefficients
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    y = filtfilt(b, a, data_)
    return y


def read_resample(file_r):
    df_ = pd.read_csv(file_r, delimiter=',')
    A_x_ = pd.DataFrame(df_, columns=['Aramis X']).to_numpy().ravel()
    A_y_ = pd.DataFrame(df_, columns=['Aramis Y']).to_numpy().ravel()
    A_z_ = pd.DataFrame(df_, columns=['Aramis Z']).to_numpy().ravel()
    A_rx_ = pd.DataFrame(df_, columns=['Aramis rX']).to_numpy().ravel()
    A_ry_ = pd.DataFrame(df_, columns=['Aramis rY']).to_numpy().ravel()
    A_rz_ = pd.DataFrame(df_, columns=['Aramis rZ']).to_numpy().ravel()
    a_y_ = pd.DataFrame(df_, columns=['Acumen Y']).to_numpy().ravel()
    a_f_ = pd.DataFrame(df_, columns=['Acumen Fy']).to_numpy().ravel()
    a_c_ = pd.DataFrame(df_, columns=['Acumen C']).to_numpy().ravel()

    return A_x_, A_y_, A_z_, A_rx_, A_ry_, A_rz_, a_y_, a_f_, a_c_


def eval_bvtv(sample, radius):
    t1 = time.time()
    check = 0
    sample_code = sample
    path_project = '/home/biomech/Documents/01_Icotec/'  # General project folder
    path_ct = path_project + '01_Experiments/02_Scans/' + sample_code + '/04_Registered/'  # Folder of CT dat
    file_bone = [filename for filename in os.listdir(path_ct + '/') if filename.endswith('image_corr.mhd')
                 and str(sample_code) in filename][0]
    file = path_ct + file_bone

    bone_grey = sitk.ReadImage(file)
    bone_img = np.transpose(sitk.GetArrayFromImage(bone_grey), [2, 1, 0])
    bone_bvtv = rR.zeros_and_ones(bone_img, 320)
    check_image = rR.zeros_and_ones(bone_img, 320)
    res = max(np.array(bone_grey.GetSpacing()))
    ori = abs((np.array(bone_grey.GetOrigin())/res).astype(int))

    # Area to evaluate
    r_mm = radius  # radius in mm
    r = int(np.rint(r_mm / res))
    length = np.rint(np.array([-45, 0]) / res).astype(int)
    drill = int(1.4/res)  # radius drill

    b = 0
    o = 0
    for z in range(min(length), max(length)):
        for y in range(-r, r):
            for x in range(-r, r):
                if r**2 >= x**2 + y**2 > drill**2:
                    check_image[x+ori[0], y+ori[1], z+ori[2]] = check_image[x+ori[0], y+ori[1], z+ori[2]] + 2
                    if bone_bvtv[x+ori[0], y+ori[1], z+ori[2]] == 1:
                        b = b+1
                    elif bone_bvtv[x+ori[0], y+ori[1], z+ori[2]] == 0:
                        o = o+1
                    else:
                        print('**ERROR**')

    bvtv = round(b/(b+o), 3)
    if check:
        print('BV/BV: ' + str(bvtv))
        plt.figure()
        plt.imshow(check_image[:, :, 800])
        plt.figure()
        plt.imshow(check_image[:, :, 600])
        plt.figure()
        plt.imshow(check_image[:, :, 400])
        plt.figure()
        plt.imshow(check_image[:, :, 200])
        plt.figure()
        plt.imshow(check_image[:, ori[1], :])
        plt.figure()
        plt.imshow(check_image[ori[0], :, :])
    tRun = time.time() - t1
    if tRun >= 3600:
        print('Execution time: ' + str(int(tRun / 3600)) + ' h ' + str(int(np.mod(tRun, 3600) / 60)) + ' min ' +
              str(round(np.mod(tRun, 60), 1)) + ' sec.')
    elif tRun >= 60:
        print('Execution time: ' + str(int(tRun / 60)) + ' min ' + str(round(np.mod(tRun, 60), 1)) + ' sec.')
    else:
        print('Execution time: ' + str(round(tRun, 1)) + ' sec.')
    return bvtv


def findPeaks(number_, co, plot_):
    sample_list_ = open('/home/biomech/Documents/01_Icotec/Specimens.txt', 'r').read().splitlines()
    data = {}
    [data['A_x'], data['A_y'], data['A_z'], data['A_rx'], data['A_ry'],
     data['A_rz'], data['a_y'], data['a_f'], data['a_c']] = \
        read_resample('/home/biomech/Documents/01_Icotec/01_Experiments/00_Data/01_MainStudy/' +
                      sample_list_[number_] + '_resample.csv')

    cutoff = co

    data_filtered = {}
    data_filtered['A_y'] = butter_lowpass_filter(data['A_y'], cutoff)
    data_filtered['a_y'] = butter_lowpass_filter(data['a_y'], cutoff)
    data_filtered['a_f'] = butter_lowpass_filter(data['a_f'], cutoff)
    peakdata = data_filtered['A_y']
    peakdata = data_filtered['a_y']
    peakf = data_filtered['a_f']
    [extAy, _] = find_peaks(-peakdata, width=20)
    n_peaks = len(extAy)
    if plot_:
        plt.close('all')
        plt.figure()
        plt.plot(peakdata, data_filtered['a_f'])
        #plt.plot(data_filtered['a_y'], data_filtered['a_f'])
        plt.figure()
        plt.plot(data_filtered['a_f'], label='force')
        plt.plot(data_filtered['A_y'], label='disp')
        plt.scatter(extAy, data_filtered['A_y'][extAy], label='disp_ext_peaks')
        plt.scatter(extAy, data_filtered['a_f'][extAy], label='f_peaks')
        plt.legend()
    return n_peaks, extAy, data_filtered['A_y'][extAy], data_filtered['a_f'][extAy]-data['a_f'][0]


sample_list = open('/home/biomech/Documents/01_Icotec/Specimens.txt', 'r').read().splitlines()

# try:
#     os.remove('/home/biomech/Documents/01_Icotec/01_Experiments/02_Scans/BVTV.txt')
#     print('Existing file has been deleted. Creating new file')
# except:
#     print('Creating new file')

# radius_mm = [3, 4, 5, 6]

radius_mm = [5]

for j in range(len(radius_mm)):
    try:
        os.remove('/home/biomech/Documents/01_Icotec/01_Experiments/02_Scans/BVTV_' + str(radius_mm[j]) +
                  '.txt')
        os.remove('/home/biomech/Documents/01_Icotec/01_Experiments/02_Scans/BVTV_' + str(radius_mm[j]) +
                  's.txt')
        print('Existing file has been deleted. Creating new file')
    except:
        print('Creating new file')

    for i in range(len(sample_list)):
        BVTV = eval_bvtv(sample_list[i], radius_mm[j])
        print('\n' + str(i) + '/' + str(len(sample_list)))
        print(BVTV)
        with open('/home/biomech/Documents/01_Icotec/01_Experiments/02_Scans/BVTV_' + str(radius_mm[j]) +
                  '.txt', 'a') as f:
            f.write(sample_list[i] + '\n')
            f.write(str(BVTV) + '\n')
        f.close()
        with open('/home/biomech/Documents/01_Icotec/01_Experiments/02_Scans/BVTV_' + str(radius_mm[j]) +
                  's.txt', 'a') as f:
            f.write(str(BVTV) + '\n')
        f.close()

tRunT = time.time() - t0
if tRunT >= 3600:
    print('Execution time (total): ' + str(int(tRunT / 3600)) + ' h ' + str(int(np.mod(tRunT, 3600) / 60)) + ' min ' +
          str(round(np.mod(tRunT, 60), 1)) + ' sec.')
elif tRunT >= 60:
    print('Execution time (total): ' + str(int(tRunT / 60)) + ' min ' + str(round(np.mod(tRunT, 60), 1)) + ' sec.')
else:
    print('Execution time (total): ' + str(round(tRunT, 1)) + ' sec.')


#%%
# Peaks Experiment
sample_list = open('/home/biomech/Documents/01_Icotec/Specimens.txt', 'r').read().splitlines()
peaks_A = np.zeros((34, 7))
peaks_F = np.zeros((34, 7))

for i in range(len(sample_list)):
    if i not in [0, 1, 2, 4, 14]:
        n_p, _, A, F = findPeaks(i, 4, 0)
        if n_p != 7:
            print(sample_list[i] + ': ' + str(n_p) + '/7 peaks detected.')
            A = [0, 0, 0, 0, 0, 0, 0]
            F = [0, 0, 0, 0, 0, 0, 0]
        elif n_p == 7:
            print(sample_list[i] + ': OK.')
            peaks_A[i, :] = A
            peaks_F[i, :] = F
    else:
        A = [0, 0, 0, 0, 0, 0, 0]
        F = [0, 0, 0, 0, 0, 0, 0]
    A_str = ' '.join(map(str, A))
    F_str = ' '.join(map(str, F))
    with open('/home/biomech/Documents/01_Icotec/01_Experiments/03_Analysis/peaks_Ax.txt', 'a') as f:
        f.write(A_str + '\n')
    f.close()
    with open('/home/biomech/Documents/01_Icotec/01_Experiments/03_Analysis/peaks_Fx.txt', 'a') as f:
        f.write(F_str + '\n')
    f.close()

#%%
# Plots to BVTV

sample_list = open('/home/biomech/Documents/01_Icotec/Specimens.txt', 'r').read().splitlines()
n = np.zeros((1, 34))

f = open('/home/biomech/Documents/01_Icotec/01_Experiments/02_Scans/BVTV_5s_corr.txt', 'r').read().splitlines()
n = np.array(f).astype(float)

fA = open('/home/biomech/Documents/01_Icotec/01_Experiments/03_Analysis/peaks_Ax.txt', 'r').read().splitlines()
A2 = np.array(fA)#.astype(float)

fF = open('/home/biomech/Documents/01_Icotec/01_Experiments/03_Analysis/peaks_Fx.txt', 'r').read().splitlines()
F2 = np.array(fF)#.astype(float)

peak = 2  # 0...6
#plt.figure()
#for i in range(34):
#    if i not in [0, 0, 2, 4, 14, 30]:
#        plt.scatter(n[i], )

