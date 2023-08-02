import numpy as np
import matplotlib.pyplot as plt
import SimpleITK as sitk
import time
import os
import ReadRawMHD as rR
import pandas as pd
from scipy.signal import butter, filtfilt, find_peaks
from ConcaveHull import ConcaveHull
from PIL import Image, ImageDraw
from skimage import morphology


def BoneEnvelope(img, res_, tolerance=3, plot=False, path='', name=''):
    '''
    Creates the concave mask from the envelope of a porouse structure. The input image needs to be binray
    :param img: 2d numpy array binary (segmented)
    :param res_: resolution of image in mm
    :param tolerance: tolernce for creating the concave envelope
                    (see also: https://gist.github.com/AndreLester/589ea1eddd3a28d00f3d7e47bd9f28fb)
    :param plot: Ture --> generates controle plots of the generated mask
    :param path: path for control plot, by default the plot will be stored in the project folder
    :param name: name of the image
    :return: binary mask --> 2d array
    '''

    # get the ponits of segmented area
    coordinates = np.array(np.where(img == 1)).T * res_

    # concave hull with ConcaveHull.py
    ch = ConcaveHull()
    ch.loadpoints(coordinates)
    ch.calculatehull(tol=tolerance)
    boundary_points = np.vstack(ch.boundary.exterior.coords.xy).T
    # # control plot of the point cloud generated form segmented image
    # plt.plot(coordinates[:, 0], coordinates[:, 1], 'o')
    # plt.plot(boundary_points[:, 0], boundary_points[:, 1], 'x--', color='r')
    # plt.show()

    # transforming the points back into the 2d array
    positionX = boundary_points[:, 1] / res_
    positionY = boundary_points[:, 0] / res_
    # # Control plot of selected points of the concave hull for the mask
    # plt.imshow(img)
    # plt.plot(positionX, positionY, '-r')
    # plt.imshow(img)
    # plt.show()

    # create mask form selected points forming a concave object
    polygon_coords = [(x, y) for x, y in zip(positionX, positionY)]
    im = np.zeros(img.shape).T

    # Create mask from coordinates
    img_ = Image.new('L', im.shape, 0)
    ImageDraw.Draw(img_).polygon(polygon_coords, outline=1, fill=1)
    mask = np.array(img_, dtype=float)
    # crate control plot of the masked region
    if plot == True:
        plt.imshow(mask + img)
        # plt.show()
        plt.savefig(path + name + 'ROIimpPosition.png', bbox_inches='tight')
        plt.close()

    return mask


def BoneMask(array_3dC, reso, tolerance, islandSize=80):
    """
    Generates a 3d array with by sampling trough ech slice and creates an envelope of the bone. The result is a 3D mask
    the bone's outer surface as boundary.
    :param array_3dC: segmented 3d array
    :return: 3d array with bone mask
    """
    # create bone mask
    # create empty 3d array for the bone mask
    mask_3d = np.zeros_like(array_3dC)

    for i in range(array_3dC.shape[1]):
        # get a slice of with bone
        slice_ = array_3dC[:, i, :]
        # clean image form noise --> delete islands smaller then 10 pixels
        slice_ = morphology.remove_small_objects(np.array(slice_, bool), islandSize)
        slice_ = slice_ * 1
        if np.sum(slice_) >= 10:
            mask_3d[:, i, :] = BoneEnvelope(slice_, reso, tolerance)
            # plt.imshow(mask_3d[:, i, :] + slice_)
            # plt.show()
    return mask_3d


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
    ori = abs((np.array(bone_grey.GetOrigin()) / res).astype(int))

    # Area to evaluate
    r_mm = radius  # radius in mm
    r = int(np.rint(r_mm / res))
    length = np.rint(np.array([-45, 0]) / res).astype(int)
    drill = int(1.4 / res)  # radius drill

    b = 0
    o = 0
    for z in range(min(length), max(length)):
        for y in range(-r, r):
            for x in range(-r, r):
                if r ** 2 >= x ** 2 + y ** 2 > drill ** 2:
                    check_image[x + ori[0], y + ori[1], z + ori[2]] = check_image[
                                                                          x + ori[0], y + ori[1], z + ori[2]] + 2
                    if bone_bvtv[x + ori[0], y + ori[1], z + ori[2]] == 1:
                        b = b + 1
                    elif bone_bvtv[x + ori[0], y + ori[1], z + ori[2]] == 0:
                        o = o + 1
                    else:
                        print('**ERROR**')

    bvtv = round(b / (b + o), 3)
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
        # plt.plot(data_filtered['a_y'], data_filtered['a_f'])
        plt.figure()
        plt.plot(data_filtered['a_f'], label='force')
        plt.plot(data_filtered['A_y'], label='disp')
        plt.scatter(extAy, data_filtered['A_y'][extAy], label='disp_ext_peaks')
        plt.scatter(extAy, data_filtered['a_f'][extAy], label='f_peaks')
        plt.legend()
    return n_peaks, extAy, data_filtered['A_y'][extAy], data_filtered['a_f'][extAy] - data['a_f'][0]


sample_list = open('/home/biomech/Documents/01_Icotec/Specimens.txt', 'r').read().splitlines()

for i in [5]:
    sample_code_ = sample_list[i]
    path_project_ = '/home/biomech/Documents/01_Icotec/'  # General project folder
    path_ct_ = path_project_ + '01_Experiments/02_Scans/' + sample_code_ + '/04_Registered/'  # Folder of CT dat
    file_bone_ = [filename for filename in os.listdir(path_ct_ + '/') if filename.endswith('image_corr.mhd')
                  and str(sample_code_) in filename][0]
    file_ = path_ct_ + file_bone_

    bone_grey_ = sitk.ReadImage(file_)
    resolution = bone_grey_.GetSpacing()[1]
    bone_img_ = np.transpose(sitk.GetArrayFromImage(bone_grey_), [2, 1, 0])
    bone_bvtv_ = rR.zeros_and_ones(bone_img_, 320)
    mask = BoneMask(bone_bvtv_, resolution, 2)

# %%
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

# %%
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

# %%
# Plots to BVTV

sample_list = open('/home/biomech/Documents/01_Icotec/Specimens.txt', 'r').read().splitlines()
n = np.zeros((1, 34))

f = open('/home/biomech/Documents/01_Icotec/01_Experiments/02_Scans/BVTV_5s_corr.txt', 'r').read().splitlines()
n = np.array(f).astype(float)

fA = pd.read_csv('/home/biomech/Documents/01_Icotec/01_Experiments/03_Analysis/peaks_Ax.txt', sep=' ')
A2 = np.array(fA)

fF = pd.read_csv('/home/biomech/Documents/01_Icotec/01_Experiments/03_Analysis/peaks_Fx.txt', sep=' ')
F2 = np.array(fF)

peak = 5  # 0...6
plt.close('all')
plt.figure()
for i in range(34):
    if i in [5, 7, 8, 10, 13, 15, 16, 18, 21, 23, 24, 26, 29, 32]:  # P
        plt.scatter(n[i], F2[i, peak], color='#1f77b4')
    elif i in [3, 6, 9, 11, 12, 17, 19, 20, 22, 25, 27, 28, 30, 33]:
        plt.scatter(n[i], F2[i, peak], color='#ff7f0e')

plt.figure()
for i in range(34):
    if not np.mod(i, 2):  # uneven
        diff = (n[i] - n[i + 1]) * 2 / (n[i] + n[i + 1])
        if i in [0, 1]:
            plt.scatter(i, diff, color='#1f77b4', label='S130684')
        elif i in [2]:
            plt.scatter(i, diff, color='#ff7f0e', label='S131318')
        elif i in [3, 4, 5, 6, 7, 8, 9, 10, 11]:
            plt.scatter(i, diff, color='#ff7f0e', label='_nolegend_')
        elif i in [12]:
            plt.scatter(i, diff, color='#2ca02c', label='S131788')
        elif i in [13, 14, 15, 16, 17, 18, 19]:
            plt.scatter(i, diff, color='#2ca02c', label='_nolegend_')
        elif i in [20]:
            plt.scatter(i, diff, color='#d62728', label='S131835')
        elif i in [21, 22, 23, 24, 25, 26, 27]:
            plt.scatter(i, diff, color='#d62728', label='_nolegend_')
        elif i in [28]:
            plt.scatter(i, diff, color='#9467bd', label='S131840')
        elif i in [29, 30, 31, 32, 33]:
            plt.scatter(i, diff, color='#9467bd', label='_nolegend_')
plt.xticks([0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32],
           ['L2', 'L1', 'L2', 'L3', 'L4', 'L5', 'L1', 'L2', 'L3', 'L4', 'L1', 'L2', 'L3', 'L4', 'L1', 'L2', 'L3'])
plt.ylim((-0.3, 0.3))
plt.title('Left/right comparison')
plt.ylabel('Difference in %')
plt.yticks([-0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3],
           ['+30% right', '+20% right', '+10% right', '0%', '+10% left', '+20% left', '+30% left', ])
plt.legend()
