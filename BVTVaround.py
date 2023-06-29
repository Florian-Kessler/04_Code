# import mapping_noRot as mappNR
import numpy as np
import matplotlib.pyplot as plt
# import matplotlib.cm as cm
import SimpleITK as sitk
import time
import os
import ReadRawMHD as rR
import pandas as pd


t0 = time.time()
def eval_bvtv(sample, radius):
    t1 = time.time()
    check = 1
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


sample_list = open('/home/biomech/Documents/01_Icotec/Specimens.txt', 'r').read().splitlines()

# try:
#     os.remove('/home/biomech/Documents/01_Icotec/01_Experiments/02_Scans/BVTV.txt')
#     print('Existing file has been deleted. Creating new file')
# except:
#     print('Creating new file')

# radius_mm = [3, 4, 5, 6]
radius_mm = [5]

for k in range(1):
    k = 3
    for j in range(len(radius_mm)):
        for i in range(2, 3):#len(sample_list)):
            BVTV = eval_bvtv(sample_list[i], radius_mm[j])
            print(BVTV)
            if k == 0:
                with open('/home/biomech/Documents/01_Icotec/01_Experiments/02_Scans/BVTV_' + str(radius_mm[j]) +
                          '.txt', 'a') as f:
                    f.write(sample_list[i] + '\n')
                    f.write(str(BVTV) + '\n')
                f.close()
            elif k == 1:
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
# Plots to BVTV

radius_mm = [3, 4, 5, 6]
average = np.zeros((len(radius_mm), 1))
n = np.zeros((len(radius_mm), 34))
for i in range(len(radius_mm)):
    f = open('/home/biomech/Documents/01_Icotec/01_Experiments/02_Scans/BVTV/BVTV_' + str(radius_mm[i])
             + 's.txt', 'r').read().splitlines()
    n[i, :] = np.array(f).astype(float)

for i in range(len(radius_mm)):
    average[i] = np.mean(n[:, i])
    print('Radius: ' + str(radius_mm[i]) + ' mm: ' + str(average[i]))

# plt.figure()
# plt.bar(sample_list, n)

plt.plot(n[:12, 2], )
