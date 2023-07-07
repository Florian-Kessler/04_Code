import numpy as np
import SimpleITK as sitk
import matplotlib.pyplot as plt
import imageio.v3 as iio


plt.close('all')
sample_list = open('/home/biomech/Documents/01_Icotec/Specimens.txt', 'r').read().splitlines()

sample_code = sample_list[10]
# Drill weird: 8

path_project = '/home/biomech/Documents/01_Icotec/'  # General project folder
path_ct = path_project + '01_Experiments/02_Scans/' + sample_code + '/04_Registered/'  # Folder of CT data
filename = path_ct + sample_code + '_screw_image.mhd'
filename_corr = path_ct + sample_code + '_screw_image_corr.mhd'
im = []
del im
save = 0


im = sitk.ReadImage(filename)
bone_img = np.transpose(sitk.GetArrayFromImage(im), [2, 1, 0])
spacing = im.GetSpacing()[0]
origin = np.array(im.GetOrigin())
origin_v = -(np.rint(origin/spacing)).astype(int)
print(origin_v[2])

png = iio.imread(path_project + '02_FEA/01_MainStudy/' + sample_code + '/80_L50_S50_D45/noCorr/mappingControlPlot.png')
plt.imshow(png)


corr_mm = 0
corr = (np.rint(corr_mm/spacing)).astype(int)


plt.figure()
plt.scatter(origin_v[2], origin_v[0], c='r')
plt.scatter(origin_v[2]-(np.rint(45/spacing)).astype(int), origin_v[0], c='r')
if corr:
    plt.scatter(origin_v[2]+corr, origin_v[0], c='b')
    plt.scatter(origin_v[2]+corr - (np.rint(45 / spacing)).astype(int), origin_v[0], c='b')
plt.imshow(bone_img[:, origin_v[1], :])


plt.figure()
plt.scatter(origin_v[2], origin_v[1], c='r')
plt.scatter(origin_v[2]-(np.rint(45/spacing)).astype(int), origin_v[1], c='r')
if corr:
    plt.scatter(origin_v[2]+corr, origin_v[1], c='b')
    plt.scatter(origin_v[2]+corr - (np.rint(45 / spacing)).astype(int), origin_v[1], c='b')
plt.imshow(bone_img[origin_v[0], :, :])


if save:
    im_corr = im
    im_corr.SetOrigin((origin[0], origin[1], origin[2] - corr_mm))
    im_corr.SetSpacing(im.GetSpacing())
    sitk.WriteImage(im_corr, filename_corr)


# # # # # Control plot # # # # #
contr = 0
if contr:
    if save:
        im_cc = sitk.ReadImage(filename_corr)
        im_corr_contr = np.transpose(sitk.GetArrayFromImage(im_cc), [2, 1, 0])
        origin_corr_contr = np.array(im_cc.GetOrigin())
        origin_v_corr_contr = -(np.rint(origin_corr_contr / spacing)).astype(int)
        plt.figure()
        plt.scatter(origin_v_corr_contr[2], origin_v_corr_contr[1], c='r')
        plt.scatter(origin_v_corr_contr[2] - (np.rint(45 / spacing)).astype(int), origin_v_corr_contr[1], c='r')
        plt.imshow(bone_img[origin_v_corr_contr[0], :, :])
