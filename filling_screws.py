import numpy as np
import matplotlib.pyplot as plt
import SimpleITK as sitk


im_itk = sitk.ReadImage('/home/biomech/Documents/01_Icotec/02_FEA/02_uFE/99_Templates/LowStiffnessScrew_v0.mhd')
im0 = sitk.GetArrayFromImage(im_itk)
im0 = np.transpose(im0, (2, 1, 0))
plt.figure()
plt.imshow(im0[:, int(im0.shape[1]/2), :])
origin = np.array([-5.75, -8.75, -45])
spacing = np.array(list(reversed(im_itk.GetSpacing())))

im1 = im0[:, :, :]
# Fill cannulation
for i in range(len(im1[1])):
    for j in range(len(im1[1])):
        if (i-66)**2+(j-66)**2 <= 18**2:
            im1[:, i, j] = 1
# Fill torx key
for i in range(len(im1[1])):
    for j in range(len(im1[1])):
        if (i-66)**2+(j-66)**2 <= 40**2:
            im1[:90, i, j] = 1

im1 = np.flip(np.swapaxes(im1, 0, 2), 1)
im1 = im1[:, ::-1, :]
im1 = np.transpose(im1, (2, 1, 0))
plt.figure()
plt.imshow(im1[:, int(im1.shape[1]/2), :])
im_itk1 = sitk.GetImageFromArray(im1, isVector=None)
im_itk1.SetSpacing(spacing)
im_itk1.SetOrigin(origin)
sitk.WriteImage(im_itk1, '/home/biomech/Documents/01_Icotec/02_FEA/02_uFE/99_Templates/'
                         'LowStiffnessScrew_v0_closed.mhd')


im_itk = sitk.ReadImage('/home/biomech/Documents/01_Icotec/02_FEA/02_uFE/99_Templates/LowStiffnessScrew_v0.mhd')
im0 = sitk.GetArrayFromImage(im_itk)
im0 = np.transpose(im0, (2, 1, 0))
imc = np.flip(np.swapaxes(im0, 0, 2), 1)
imc = imc[:, ::-1, :]
imc = np.transpose(imc, (2, 1, 0))
plt.figure()
plt.imshow(imc[:, int(imc.shape[1]/2), :])
im_itkc = sitk.GetImageFromArray(imc, isVector=None)
im_itkc.SetSpacing(spacing)
im_itkc.SetOrigin(origin)
sitk.WriteImage(im_itkc, '/home/biomech/Documents/01_Icotec/02_FEA/02_uFE/99_Templates/'
                         'LowStiffnessScrew_v0_can.mhd')
#%% Cut head off from DPS screw
plt.close('all')
im = sitk.ReadImage('/home/biomech/Documents/01_Icotec/02_FEA/02_uFE/99_Templates/DPSScrew_v0.mhd')
img = sitk.GetArrayFromImage(im)
img = np.transpose(img, (2, 1, 0))
plt.figure()
plt.imshow(img[:, int(img.shape[1]/2), :])
plt.title('Original')
print(img.shape)
# Cut tulip axially
img1 = img[:, :, 8:832]

# Cut tulip radially
for i in range(img1.shape[0]):
    for j in range(img1.shape[1]):
        if (i-92)**2+(j-111)**2 >= 59**2:
            img1[i, j, :] = 0
for i in range(img1.shape[0]):
    for j in range(img1.shape[1]):
        if (i-92)**2+(j-111)**2 >= 54**2:
            img1[i, j, :744] = 0
# Add material to screw head
for i in range(img1.shape[0]):
    for j in range(img1.shape[1]):
        if 36**2 <= (i-92)**2+(j-113)**2 <= 59**2:
            img1[i, j, 770:] = 1

img1 = img1[25:-26, 45:-42, :]  # x1 32-> 31?
plt.figure()
plt.imshow(img1[:, int(img1.shape[1]/2), :])
plt.title('Cut tulip, for FE')
img1 = np.transpose(img1, (2, 1, 0))
im_itk_ = sitk.GetImageFromArray(img1, isVector=None)
im_itk_.SetSpacing(spacing)
im_itk_.SetOrigin(origin)
sitk.WriteImage(im_itk_, '/home/biomech/Documents/01_Icotec/02_FEA/02_uFE/99_Templates/DPSScrew_v0_can.mhd')

img1 = np.transpose(img, (2, 1, 0))
img1 = np.flip(np.swapaxes(img1, 0, 2), 1)
img1 = img1[:, ::-1, 8:832]


# Fill cannulation
for i in range(img1.shape[0]):
    for j in range(img1.shape[1]):
        if (i-92)**2+(j-113)**2 <= 18**2:
            img1[i, j, :] = 1
# Fill torx key
for i in range(img1.shape[0]):
    for j in range(img1.shape[1]):
        if (i-92)**2+(j-113)**2 <= 47**2:
            img1[i, j, 737:] = 1
# Fill fenestration
for i in range(img1.shape[0]):
    for j in range(img1.shape[1]):
        if (i-92)**2+(j-113)**2 <= 28**2:
            img1[i, j, 135:280] = 1
for i in range(img1.shape[0]):
    for j in range(img1.shape[1]):
        if (i-92)**2+(j-113)**2 <= 25**2:
            img1[i, j, 90:280] = 1
# Cut image to shape
img1 = img1[25:-26, 45:-42, :]  # x1 32-> 31?
print(img1.shape)

plt.figure()
plt.imshow(img1[:, int(img1.shape[1]/2), :])
plt.title('Y-cut')
plt.figure()
plt.imshow(img1[int(img1.shape[0]/2), :, :])
plt.title('X-cut')
# plt.figure()
# plt.imshow(img1[:, :, 777])


# img1 = np.flip(np.swapaxes(img1, 0, 2), 1)
# img1 = img1[:, :, ::-1]
img1 = np.transpose(img1, (2, 1, 0))
# plt.figure()
# plt.imshow(img1[:, int(img1.shape[1]/2), :])
im_itk_ = sitk.GetImageFromArray(img1, isVector=None)
im_itk_.SetSpacing(spacing)
im_itk_.SetOrigin(origin)
sitk.WriteImage(im_itk_, '/home/biomech/Documents/01_Icotec/02_FEA/02_uFE/99_Templates/DPSScrew_v0_closed.mhd')
