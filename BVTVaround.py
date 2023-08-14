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
import statsmodels.api as sm


def BoneEnvelope(img, res_, tolerance=3, plot=False, path='', name=''):
    """
    Creates the concave mask from the envelope of a porous structure. The input image needs to be binary
    :param img: 2d numpy array binary (segmented)
    :param res_: resolution of image in mm
    :param tolerance: tolerance for creating the concave envelope
                    (see also: https://gist.github.com/AndreLester/589ea1eddd3a28d00f3d7e47bd9f28fb)
    :param plot: Ture --> generates control plots of the generated mask
    :param path: path for control plot, by default the plot will be stored in the project folder
    :param name: name of the image
    :return: binary mask --> 2d array
    """

    # get the points of segmented area
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
    if plot:
        plt.imshow(mask + img)
        # plt.show()
        plt.savefig(path + name + 'ROIimpPosition.png', bbox_inches='tight')
        plt.close()

    return mask


def BoneMask(array_3dC, reso, axis, tolerance, islandSize=2):
    """
    Generates a 3d array with by sampling trough ech slice and creates an envelope of the bone. The result is a 3D mask
    the bone's outer surface as boundary.
    :param array_3dC: segmented 3d array
    :param reso: resolution in mm
    :param axis: 0, 1, 2 axis
    :param tolerance: tolerance of mask shape
    :param islandSize: size of islands to delete, default=2
    :return: 3d array with bone mask
    """
    # create bone mask
    # create empty 3d array for the bone mask
    mask_3d = np.zeros_like(array_3dC)
    if axis == 0:
        for i_ in range(array_3dC.shape[0]):
            # get a slice of with bone
            slice_ = array_3dC[i_, :, :]
            # clean image form noise --> delete islands smaller then 10 pixels
            slice_ = morphology.remove_small_objects(np.array(slice_, bool), islandSize)
            slice_ = slice_ * 1
            if np.sum(slice_) >= 10:
                mask_3d[i_, :, :] = BoneEnvelope(slice_, reso, tolerance)
                # plt.imshow(mask_3d[i, :, :] + slice_)
                # plt.show()
    elif axis == 1:
        for i_ in range(array_3dC.shape[1]):
            # get a slice of with bone
            slice_ = array_3dC[:, i_, :]
            # clean image form noise --> delete islands smaller then 10 pixels
            slice_ = morphology.remove_small_objects(np.array(slice_, bool), islandSize)
            slice_ = slice_ * 1
            if np.sum(slice_) >= 10:
                mask_3d[:, i_, :] = BoneEnvelope(slice_, reso, tolerance)
                # plt.imshow(mask_3d[:, i, :] + slice_)
                # plt.show()
    elif axis == 2:
        for i_ in range(array_3dC.shape[2]):
            # get a slice of with bone
            slice_ = array_3dC[:, :, i_]
            # clean image form noise --> delete islands smaller then 10 pixels
            slice_ = morphology.remove_small_objects(np.array(slice_, bool), islandSize)
            slice_ = slice_ * 1
            if np.sum(slice_) >= 10:
                mask_3d[:, :, i_] = BoneEnvelope(slice_, reso, tolerance)
                # plt.imshow(mask_3d[:, :, i] + slice_)
                # plt.show()
    return mask_3d


def butter_lowpass_filter(data_, cutoff_, order=9):
    fs = 10  # sample rate, Hz
    nyq = 0.5 * fs
    normal_cutoff = cutoff_ / nyq
    # Get the filter coefficients
    b, a = butter(order, normal_cutoff, btype='low', analog=False, output='ba')
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
    file_bone = [filename for filename in os.listdir(path_ct + '/') if filename.endswith('image.mhd')
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

    bvtv_ = round(b / (b + o), 3)
    print('BV/TV: ' + str(bvtv_))
    if check:
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
    tRun_ = time.time() - t1
    if tRun_ >= 3600:
        print('Execution time: ' + str(int(tRun_ / 3600)) + ' h ' + str(int(np.mod(tRun_, 3600) / 60)) + ' min ' +
              str(round(np.mod(tRun_, 60), 1)) + ' sec.')
    elif tRun_ >= 60:
        print('Execution time: ' + str(int(tRun_ / 60)) + ' min ' + str(round(np.mod(tRun_, 60), 1)) + ' sec.')
    else:
        print('Execution time: ' + str(round(tRun_, 1)) + ' sec.')
    return bvtv_


def eval_bvtv_mask(sample, radius):
    t1 = time.time()
    check = 0
    sample_code = sample
    path_project = '/home/biomech/Documents/01_Icotec/'  # General project folder
    path_ct = path_project + '01_Experiments/02_Scans/' + sample_code + '/04_Registered/'  # Folder of CT dat
    path_mask = '/home/biomech/DATA/01_Icotec/01_Experiments/02_Scans/BVTV/' + sample_code
    mask_X_ = np.load(path_mask + '_mask_x.npy')
    mask_Y_ = np.load(path_mask + '_mask_y.npy')
    mask_Z_ = np.load(path_mask + '_mask_z.npy')
    mask = ((mask_X_ + mask_Y_ + mask_Z_) >= 1).astype(int)
    file_bone = [filename for filename in os.listdir(path_ct + '/') if filename.endswith('image.mhd')
                 and str(sample_code) in filename][0]
    file = path_ct + file_bone

    bone_grey = sitk.ReadImage(file)
    bone_img = np.transpose(sitk.GetArrayFromImage(bone_grey), [2, 1, 0])
    bone_bvtv = rR.zeros_and_ones(bone_img, 320) + mask
    check_image = rR.zeros_and_ones(bone_img, 320) + mask * 4
    res = max(np.array(bone_grey.GetSpacing()))
    ori = abs((np.array(bone_grey.GetOrigin()) / res).astype(int))

    # Area to evaluate
    r_mm = radius  # radius in mm
    r = int(np.rint(r_mm / res))
    length = np.rint(np.array([-45, 0]) / res).astype(int)
    drill = int(1.4 / res)  # radius drill

    b = 0
    o = 0
    e = 0
    for z in range(min(length), max(length)):
        for y in range(-r, r):
            for x in range(-r, r):
                if r ** 2 >= x ** 2 + y ** 2 > drill ** 2:
                    check_image[x + ori[0], y + ori[1], z + ori[2]] = check_image[
                                                                          x + ori[0], y + ori[1], z + ori[2]] + 2
                    if bone_bvtv[x + ori[0], y + ori[1], z + ori[2]] == 2:
                        b = b + 1
                    elif bone_bvtv[x + ori[0], y + ori[1], z + ori[2]] == 1:
                        o = o + 1
                    else:
                        e = e + 1

    bvtv_ = round(b / (b + o), 3)
    print('BV/TV: ' + str(bvtv_))
    if check:
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
    tRun_ = time.time() - t1
    if tRun_ >= 3600:
        print('Execution time: ' + str(int(tRun_ / 3600)) + ' h ' + str(int(np.mod(tRun_, 3600) / 60)) + ' min ' +
              str(round(np.mod(tRun_, 60), 1)) + ' sec.')
    elif tRun_ >= 60:
        print('Execution time: ' + str(int(tRun_ / 60)) + ' min ' + str(round(np.mod(tRun_, 60), 1)) + ' sec.')
    else:
        print('Execution time: ' + str(round(tRun_, 1)) + ' sec.')
    return bvtv_


def eval_bvtv_mask_along(sample, radius):
    t1 = time.time()
    check = 0
    sample_code = sample
    path_project = '/home/biomech/Documents/01_Icotec/'  # General project folder
    path_ct = path_project + '01_Experiments/02_Scans/' + sample_code + '/04_Registered/'  # Folder of CT dat
    path_mask = '/home/biomech/DATA/01_Icotec/01_Experiments/02_Scans/BVTV/' + sample_code
    mask_X_ = np.load(path_mask + '_mask_x.npy')
    mask_Y_ = np.load(path_mask + '_mask_y.npy')
    mask_Z_ = np.load(path_mask + '_mask_z.npy')
    mask = ((mask_X_ + mask_Y_ + mask_Z_) >= 1).astype(int)
    file_bone = [filename for filename in os.listdir(path_ct + '/') if filename.endswith('image.mhd')
                 and str(sample_code) in filename][0]
    file = path_ct + file_bone

    bone_grey = sitk.ReadImage(file)
    bone_img = np.transpose(sitk.GetArrayFromImage(bone_grey), [2, 1, 0])
    bone_bvtv = rR.zeros_and_ones(bone_img, 320) + mask * 3
    check_image = rR.zeros_and_ones(bone_img, 320) + mask * 4
    res = max(np.array(bone_grey.GetSpacing()))
    ori = abs((np.array(bone_grey.GetOrigin()) / res).astype(int))

    # Area to evaluate
    r_mm = radius  # radius in mm
    r = int(np.rint(r_mm / res))
    length = np.rint(np.array([-45, 0]) / res).astype(int)
    drill = int(1.4 / res)  # radius drill
    bvtv_along_ = np.zeros(len(range(min(length), max(length))))
    bv = 0
    ev = 0
    for z in range(min(length), max(length)):
        bv_z = 0
        ev_z = 0
        for y in range(-r, r):
            for x in range(-r, r):
                if r ** 2 >= x ** 2 + y ** 2 > drill ** 2:
                    check_image[x + ori[0], y + ori[1], z + ori[2]] = check_image[
                                                                          x + ori[0], y + ori[1], z + ori[2]] + 2
                    if bone_bvtv[x + ori[0], y + ori[1], z + ori[2]] == 4:
                        bv = bv + 1
                        bv_z = bv_z + 1
                    elif bone_bvtv[x + ori[0], y + ori[1], z + ori[2]] == 3:
                        ev = ev + 1
                        ev_z = ev_z + 1
        if (bv_z + ev_z) != 0:
            bvtv_along_[z] = bv_z / (bv_z + ev_z)
        else:
            bvtv_along_[z] = 0
    tv = bv + ev
    bvtv_ = round(bv / tv, 3)
    print('BV/TV: ' + str(bvtv_))
    if check:
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
    tRun_ = time.time() - t1
    if tRun_ >= 3600:
        print('Execution time: ' + str(int(tRun_ / 3600)) + ' h ' + str(int(np.mod(tRun_, 3600) / 60)) + ' min ' +
              str(round(np.mod(tRun_, 60), 1)) + ' sec.')
    elif tRun_ >= 60:
        print('Execution time: ' + str(int(tRun_ / 60)) + ' min ' + str(round(np.mod(tRun_, 60), 1)) + ' sec.')
    else:
        print('Execution time: ' + str(round(tRun_, 1)) + ' sec.')
    return bvtv_, bvtv_along_


def eval_bvtv_mask_along_load(sample, radius):
    t1 = time.time()
    check = 0
    sample_code = sample
    path_project = '/home/biomech/Documents/01_Icotec/'  # General project folder
    path_ct = path_project + '01_Experiments/02_Scans/' + sample_code + '/04_Registered/'  # Folder of CT dat
    path_mask = '/home/biomech/DATA/01_Icotec/01_Experiments/02_Scans/BVTV/' + sample_code
    mask_X_ = np.load(path_mask + '_mask_x.npy')
    mask_Y_ = np.load(path_mask + '_mask_y.npy')
    mask_Z_ = np.load(path_mask + '_mask_z.npy')
    mask = ((mask_X_ + mask_Y_ + mask_Z_) >= 1).astype(int)
    file_bone = [filename for filename in os.listdir(path_ct + '/') if filename.endswith('image.mhd')
                 and str(sample_code) in filename][0]
    file = path_ct + file_bone

    bone_grey = sitk.ReadImage(file)
    bone_img = np.transpose(sitk.GetArrayFromImage(bone_grey), [2, 1, 0])
    bone_bvtv = rR.zeros_and_ones(bone_img, 320) + mask * 3
    check_image = rR.zeros_and_ones(bone_img, 320) + mask * 4
    res = max(np.array(bone_grey.GetSpacing()))
    ori = abs((np.array(bone_grey.GetOrigin()) / res).astype(int))

    # Area to evaluate
    r_mm = radius  # radius in mm
    r = int(np.rint(r_mm / res))
    length = np.rint(np.array([-45, 0]) / res).astype(int)
    drill = int(1.4 / res)  # radius drill
    bvtv_along_ = np.zeros(len(range(min(length), max(length))))
    bv = 0
    ev = 0
    for z in range(min(length), max(length)):
        bv_z = 0
        ev_z = 0
        for y in range(-r, 0):  # r -> 0 to only include loaded side. maybe even reduce to -2r/3??
            for x in range(-r, r):
                if r ** 2 >= x ** 2 + y ** 2 > drill ** 2:
                    check_image[x + ori[0], y + ori[1], z + ori[2]] = check_image[
                                                                          x + ori[0], y + ori[1], z + ori[2]] + 2
                    if bone_bvtv[x + ori[0], y + ori[1], z + ori[2]] == 4:
                        bv = bv + 1
                        bv_z = bv_z + 1
                    elif bone_bvtv[x + ori[0], y + ori[1], z + ori[2]] == 3:
                        ev = ev + 1
                        ev_z = ev_z + 1
        if (bv_z + ev_z) != 0:
            bvtv_along_[z] = bv_z / (bv_z + ev_z)
        else:
            bvtv_along_[z] = 0
    tv = bv + ev
    bvtv_ = round(bv / tv, 3)
    print('BV/TV: ' + str(bvtv_))
    if check:
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
    tRun_ = time.time() - t1
    if tRun_ >= 3600:
        print('Execution time: ' + str(int(tRun_ / 3600)) + ' h ' + str(int(np.mod(tRun_, 3600) / 60)) + ' min ' +
              str(round(np.mod(tRun_, 60), 1)) + ' sec.')
    elif tRun_ >= 60:
        print('Execution time: ' + str(int(tRun_ / 60)) + ' min ' + str(round(np.mod(tRun_, 60), 1)) + ' sec.')
    else:
        print('Execution time: ' + str(round(tRun_, 1)) + ' sec.')
    return bvtv_, bvtv_along_


def findPeaks(number_, cutoff, plot_):
    sample_list_ = open('/home/biomech/Documents/01_Icotec/Specimens.txt', 'r').read().splitlines()
    data = {}
    [data['A_x'], data['A_y'], data['A_z'], data['A_rx'], data['A_ry'],
     data['A_rz'], data['a_y'], data['a_f'], data['a_c']] = \
        read_resample('/home/biomech/Documents/01_Icotec/01_Experiments/00_Data/01_MainStudy/' +
                      sample_list_[number_] + '_resample.csv')

    data_filtered = {}
    data_filtered['A_y'] = butter_lowpass_filter(data['A_y'], cutoff)
    data_filtered['a_y'] = butter_lowpass_filter(data['a_y'], cutoff)
    data_filtered['a_f'] = butter_lowpass_filter(data['a_f'], cutoff)
    # peakdata = data_filtered['A_y']
    peakdata = data_filtered['a_y']
    # peakf = data_filtered['a_f']
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


def lin_reg(X, Y):
    X = X.flatten().ravel()
    Y = Y.flatten()
    X = X[X != 0]
    Y = Y[X != 0]
    X = X[Y != 0]
    Y = Y[Y != 0]
    X = sm.add_constant(X)  # Add a constant term to the independent variable array
    mod = sm.OLS(Y, X)  # y, X
    reg = mod.fit()
    return reg, X, Y


class IndexTracker(object):
    def __init__(self, ax_, X):
        self.ax = ax_
        ax_.set_title('use scroll wheel to navigate images')

        self.X = X
        rows, self.slices, cols = X.shape
        self.ind = self.slices // 2

        self.im = ax_.imshow(self.X[:, self.ind, :])
        self.update()

    def onscroll(self, event):
        print("%s %s" % (event.button, event.step))
        if event.button == 'up':
            self.ind = (self.ind + 1) % self.slices
        else:
            self.ind = (self.ind - 1) % self.slices
        self.update()

    def update(self):
        self.im.set_data(self.X[:, self.ind, :])
        self.ax.set_ylabel('slice %s' % self.ind)
        self.im.axes.figure.canvas.draw()


t0 = time.time()
sample_list = open('/home/biomech/Documents/01_Icotec/Specimens.txt', 'r').read().splitlines()
path_bvtv = '/home/biomech/DATA/01_Icotec/01_Experiments/02_Scans/BVTV/'  # on DATA drive, not in Documents!!!
path_project_ = '/home/biomech/Documents/01_Icotec/'  # General project folder

'''
for i in [0, 1]:  # range(3, 34):
    t0 = time.time()
    sample_code_ = sample_list[i]
    path_ct_ = path_project_ + '01_Experiments/02_Scans/' + sample_code_ + '/04_Registered/'  # Folder of CT dat
    file_bone_ = [filename for filename in os.listdir(path_ct_ + '/') if filename.endswith('image.mhd')
                  and str(sample_code_) in filename][0]
    file_ = path_ct_ + file_bone_
    print('\nWorking on ' + sample_code_ + '...')
    bone_grey_ = sitk.ReadImage(file_)
    resolution = bone_grey_.GetSpacing()[1]
    bone_img_ = np.transpose(sitk.GetArrayFromImage(bone_grey_), [2, 1, 0])
    bone_bvtv_ = rR.zeros_and_ones(bone_img_, 320)
    mask_x = BoneMask(bone_bvtv_, resolution, 0, 2)
    mask_y = BoneMask(bone_bvtv_, resolution, 1, 2)
    mask_z = BoneMask(bone_bvtv_, resolution, 2, 2)

    np.save(path_bvtv + sample_code_ + '_mask_x.npy', mask_x)
    np.save(path_bvtv + sample_code_ + '_mask_y.npy', mask_y)
    np.save(path_bvtv + sample_code_ + '_mask_z.npy', mask_z)
    tRun = time.time() - t0
    if tRun >= 3600:
        print('Execution time: ' + str(int(tRun / 3600)) + ' h ' + str(int(np.mod(tRun, 3600) / 60)) + ' min ' +
              str(round(np.mod(tRun, 60), 1)) + ' sec.\n')
    elif tRun >= 60:
        print('Execution time: ' + str(int(tRun / 60)) + ' min ' + str(round(np.mod(tRun, 60), 1)) + ' sec.\n')
    else:
        print('Execution time: ' + str(round(tRun, 1)) + ' sec.\n')
'''

# %%
'''
ii = 5
sample_code_ = sample_list[ii]
maskX = np.load(path_bvtv + sample_code_ + '_mask_x.npy')
maskY = np.load(path_bvtv + sample_code_ + '_mask_y.npy')
maskZ = np.load(path_bvtv + sample_code_ + '_mask_z.npy')

#%%
mask_add = ((maskX + maskY + maskZ) >= 1).astype(int)
mask_mix = ((maskX + maskY + maskZ) >= 2).astype(int)
mask_mul = maskX * maskY * maskZ
'''
# %%
'''
samples = [5]  # np.arange(0, 34)
radius_mm = [6]
radius_mm_str = ['6']
bvtv_mask = np.zeros((4, 33))
bvtv = np.zeros((4, 33))

for jj in range(len(radius_mm)):
    print('Radius: ' + str(radius_mm[jj]))
    try:
        os.remove('/home/biomech/Documents/01_Icotec/01_Experiments/02_Scans/BVTV/x3BVTV_mask_'
                  + str(radius_mm_str[jj]) + '.txt')
        os.remove('/home/biomech/Documents/01_Icotec/01_Experiments/02_Scans/BVTV/x3BVTV_mask_'
                  + str(radius_mm_str[jj]) + 's.txt')
        print('Existing mask file has been deleted. Creating new mask file.')
    except FileNotFoundError:
        print('Creating new file.')
    try:
        os.remove('/home/biomech/Documents/01_Icotec/01_Experiments/02_Scans/BVTV/x3BVTV_'
                  + str(radius_mm_str[jj]) + '.txt')
        os.remove('/home/biomech/Documents/01_Icotec/01_Experiments/02_Scans/BVTV/x3BVTV_'
                  + str(radius_mm_str[jj]) + 's.txt')
        print('Existing file has been deleted. Creating new file.')
    except FileNotFoundError:
        print('Creating new mask file.')
    for ii in samples:
        t1 = time.time()
        BVTV_mask, BVTV_mask_along = eval_bvtv_mask_along(sample_list[ii], radius_mm[jj])
        print('\n-----\n-----\nBV/TV along:')
        print(BVTV_mask_along)
        BVTV = eval_bvtv(sample_list[ii], radius_mm[jj])
        print('\n' + str(ii) + '/' + str(len(sample_list)))
        print('Sample: ' + sample_list[ii])
        with open('/home/biomech/Documents/01_Icotec/01_Experiments/02_Scans/BVTV/x3BVTV_' + str(radius_mm_str[jj]) +
                  '.txt', 'a') as f:
            f.write(sample_list[ii] + '\n')
            f.write(str(BVTV) + '\n')
        f.close()
        with open('/home/biomech/Documents/01_Icotec/01_Experiments/02_Scans/BVTV/x3BVTV_' + str(radius_mm_str[jj]) +
                  's.txt', 'a') as f:
            f.write(str(BVTV) + '\n')
        f.close()
        with open('/home/biomech/Documents/01_Icotec/01_Experiments/02_Scans/BVTV/x3BVTV_mask_' + str(radius_mm_str[jj])
                  + '.txt', 'a') as f:
            f.write(sample_list[ii] + '\n')
            f.write(str(BVTV_mask) + '\n')
        f.close()
        with open('/home/biomech/Documents/01_Icotec/01_Experiments/02_Scans/BVTV/x3BVTV_mask_' + str(radius_mm_str[jj])
                  + 's.txt', 'a') as f:
            f.write(str(BVTV_mask) + '\n')
        f.close()
        tRun_ = time.time() - t1
        if tRun_ >= 3600:
            print('Execution time: ' + str(int(tRun_ / 3600)) + ' h ' + str(int(np.mod(tRun_, 3600) / 60)) + ' min ' +
                  str(round(np.mod(tRun_, 60), 1)) + ' sec.\n')
        elif tRun_ >= 60:
            print('Execution time: ' + str(int(tRun_ / 60)) + ' min ' + str(round(np.mod(tRun_, 60), 1)) + ' sec.\n')
        else:
            print('Execution time: ' + str(round(tRun_, 1)) + ' sec.\n')
    print(str(radius_mm[jj]) + ' mm ROI done.')
tRun = time.time() - t0
if tRun >= 3600:
    print('Execution time: ' + str(int(tRun / 3600)) + ' h ' + str(int(np.mod(tRun, 3600) / 60)) + ' min ' +
          str(round(np.mod(tRun, 60), 1)) + ' sec.\n')
elif tRun >= 60:
    print('Execution time: ' + str(int(tRun / 60)) + ' min ' + str(round(np.mod(tRun, 60), 1)) + ' sec.\n')
else:
    print('Execution time: ' + str(round(tRun, 1)) + ' sec.\n')
'''
# %%
'''
fig, ax = plt.subplots(1, 1)

# X = np.random.rand(20, 20, 40)

tracker = IndexTracker(ax, mask_add)


fig.canvas.mpl_connect('scroll_event', tracker.onscroll)
plt.show()'''
# %%
'''
plt.figure()
s = 350
plt.imshow(mask_y[:, s, :] + bone_bvtv_[:, s, :])
plt.figure()
plt.imshow(mask_x[200, :, :] + bone_bvtv_[200, :, :])
plt.figure()
plt.imshow(mask_z[:, :, 500] + bone_bvtv_[:, :, 500])

mask_add = mask_x + mask_y + mask_z
mask_add_b = (mask_add >= 1).astype(int)
mask_mul = mask_x * mask_y * mask_z
plt.figure()
plt.imshow(mask_add_b[:, s, :] + bone_bvtv_[:, s, :])
plt.figure()
plt.imshow(mask_add_b[200, :, :] + bone_bvtv_[200, :, :])
plt.figure()
plt.imshow(mask_add_b[:, :, 500] + bone_bvtv_[:, :, 500])
plt.figure()
plt.imshow(mask_mul[:, s, :] + bone_bvtv_[:, s, :])
plt.figure()
plt.imshow(mask_mul[200, :, :] + bone_bvtv_[200, :, :])
plt.figure()
plt.imshow(mask_mul[:, :, 500] + bone_bvtv_[:, :, 500])
'''

# %%
# try:
#     os.remove('/home/biomech/Documents/01_Icotec/01_Experiments/02_Scans/BVTV.txt')
#     print('Existing file has been deleted. Creating new file')
# except:
#     print('Creating new file')

# radius_mm = [3, 4, 5, 6]
'''
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
'''
# %%
'''
radius = 6
radius_file = 6
file = '/home/biomech/Documents/01_Icotec/01_Experiments/02_Scans/BVTV/xBVTV_' + str(radius_file) + 's.txt'
bvtv_wo = pd.read_csv(file)
file = '/home/biomech/Documents/01_Icotec/01_Experiments/02_Scans/BVTV/xBVTV_mask_' + str(radius_file) + 's.txt'
bvtv_wm = pd.read_csv(file)
plt.figure()
plt.scatter(bvtv_wo, bvtv_wm)

regression_T, xx_T, yy_T = lin_reg(np.array(bvtv_wo), np.array(bvtv_wm))
plt.plot(np.array([0, 0.225]), np.array([0, 0.225]) * regression_T.params[1] + regression_T.params[0],
         color='k', linestyle='dotted', label=str(np.round(regression_T.params[1], 3)) + 'x + ' +
                                              str(np.round(regression_T.params[0], 3)))
if regression_T.pvalues[1] >= 0.05:
    lab_pvalue_T = 'p = ' + str(np.round(regression_T.pvalues[1], 2))
else:
    lab_pvalue_T = 'p < 0.05'
plt.plot([-1, 0], [-1, 0], color='w', linestyle='dashed',
         label='R$^2$ = {:0.2f}'.format(np.round(regression_T.rsquared, 2)))
plt.plot([-1, 0], [-1, 0], color='w', label=lab_pvalue_T)

plt.title('Radius: ' + str(radius) + ' mm')
plt.plot([0, 0.225], [0, 0.225], 'k')
plt.xlabel('BVTV w/o mask')
plt.ylabel('BVTV with mask')
plt.legend()
plt.xlim([0, 0.225])
plt.ylim([0, 0.225])

# %%
sample_code_ = sample_list[8]
path_bvtv = '/home/biomech/DATA/01_Icotec/01_Experiments/02_Scans/BVTV/'  # on DATA drive, not in Documents!!!
maskX = np.load(path_bvtv + sample_code_ + '_mask_x.npy')
maskY = np.load(path_bvtv + sample_code_ + '_mask_y.npy')
maskZ = np.load(path_bvtv + sample_code_ + '_mask_z.npy')
mask = ((maskX + maskY + maskZ)>=1).astype(int)
#%%
plt.figure()
plt.imshow(mask[:, 300, :])
'''
#%%

radius_list = [6]  # working on: 5, done: 4, 4.5
radius_list_str = ['6']
# BVTV_mask, BVTV_mask_along = eval_bvtv_mask_along(sample_list[8], 4.5)
BVTV_mask = 0
BVTV_mask_along = 0
for i in range(len(sample_list)):
    print('\nStart ' + sample_list[i] + '...')
    for j in range(len(radius_list)):
        del BVTV_mask
        del BVTV_mask_along
        BVTV_mask, BVTV_mask_along = eval_bvtv_mask_along_load(sample_list[i], radius_list[j])
        plt.figure()
        plt.plot(BVTV_mask_along)
        plt.savefig(path_bvtv + 'BVTV_along_load_' + sample_list[i] + '_' + radius_list_str[j] + 'mm.png')
        plt.close('all')
        np.save(path_bvtv + 'BVTV_along_load_' + sample_list[i] + '_' + radius_list_str[j] + 'mm', BVTV_mask_along)
        print(sample_list[i] + ' on radius ' + str(radius_list[j]) + ' mm finished.')
    print(sample_list[i] + ' finished.\n')
tRun = time.time() - t0
if tRun >= 3600:
    print('Execution time: ' + str(int(tRun / 3600)) + ' h ' + str(int(np.mod(tRun, 3600) / 60)) + ' min ' +
          str(round(np.mod(tRun, 60), 1)) + ' sec.\n')
elif tRun >= 60:
    print('Execution time: ' + str(int(tRun / 60)) + ' min ' + str(round(np.mod(tRun, 60), 1)) + ' sec.\n')
else:
    print('Execution time: ' + str(round(tRun, 1)) + ' sec.\n')
