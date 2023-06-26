import mapping_noRot as mappNR
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import SimpleITK as sitk
import time
import os
import ReadRawMHD as rR


def eval_bvtv(sample):
    t1 = time.time()
    sample_code = sample
    path_project = '/home/biomech/Documents/01_Icotec/'  # General project folder
    path_ct = path_project + '01_Experiments/02_Scans/' + sample_code + '/04_Registered/'  # Folder of CT dat
    file_bone = [filename for filename in os.listdir(path_ct + '/') if filename.endswith('image.mhd') and str(sample_code) in filename][0]
    file = path_ct + file_bone

    bone_grey = sitk.ReadImage(file)
    bone_img = np.transpose(sitk.GetArrayFromImage(bone_grey), [2, 1, 0])
    bone_bvtv = rR.zeros_and_ones(bone_img, 320)
    check_image = rR.zeros_and_ones(bone_img, 320)
    res = max(np.array(bone_grey.GetSpacing()))
    ori = abs((np.array(bone_grey.GetOrigin())/res).astype(int))
    print(res)

    # Area to evaluate
    r_mm = 6  # radius in mm
    r = int(np.rint(r_mm / res))
    length = np.rint(np.array([-45, 0]) / res).astype(int)
    drill = int(1.4/res)  # radius drill

    b = 0
    o = 0
    ii = 0
    jj = 0
    for z in range(min(length), max(length)):
        for y in range(-r, r):
            for x in range(-r, r):
                if r**2 >= x**2 + y**2 > drill**2:
                    check_image[x+ori[0], y+ori[1], z+ori[2]] = check_image[x+ori[0], y+ori[1], z+ori[2]] + 3
                    if bone_bvtv[x+ori[0], y+ori[1], z+ori[2]] == 1:
                        b = b+1
                    elif bone_bvtv[x+ori[0], y+ori[1], z+ori[2]] == 0:
                        o = o+1
                    else:
                        print('**ERROR**')
    bvtv = round(b/(b+o), 3)
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


sample_list = ['S130684_L2_left', 'S130684_L2_right',
               'S131318_L1_left', 'S131318_L1_right', 'S131318_L2_left', 'S131318_L2_right',
               'S131318_L3_left', 'S131318_L3_right', 'S131318_L4_left', 'S131318_L4_right',
               'S131318_L5_left', 'S131318_L5_right',
               'S131788_L1_left', 'S131788_L1_right', 'S131788_L2_left', 'S131788_L2_right',
               'S131788_L3_left', 'S131788_L3_right', 'S131788_L4_left', 'S131788_L4_right',
               'S131835_L1_left', 'S131835_L1_right', 'S131835_L2_left', 'S131835_L2_right',
               'S131835_L3_left', 'S131835_L3_right', 'S131835_L4_left', 'S131835_L4_right',
               'S131840_L1_left', 'S131840_L1_right', 'S131840_L2_left', 'S131840_L2_right',
               'S131840_L3_left', 'S131840_L3_right']

i = 2
BVTV = eval_bvtv(sample_list[i])
print(sample_list[i])
print(BVTV)
