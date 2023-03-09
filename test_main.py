import ReadRawMHD as rR
import mapping as mapp
import numpy as np
import matplotlib.pyplot as plt
import SimpleITK as sitk
import time
import os
t1 = time.time()


# # # # # Input # # # # #
# Input imaging data
models = ['00_L50_S50_D30', '01_L50_S50_D35', '02_L50_S50_D40', '10_L50_S00_D30', '11_L50_S00_D35', '12_L50_S00_D40',
          '50_L50_S00_D30', '55_L50_S00_D30']
model_code = models[0]  # FEA model name
sample_code = '05_Pilot5'  # Experimental sample. contains number and sample name, e.g. 00_Pilot3.
ExpScrew = 'T'  # T or P, site from experiment. Two input files (SimScrew Ti and PEEK) will be generated for each site
# file_bone = 'XCT_Icotec_S130672_L5_intact_planned.mhd'  # gray value bone ct scan, gauss filter applied  pilot3
# file_bone = 'Icotec_S130672_L4_intact.mhd'  # pilot4
file_bone = 'Icotec_S130684_L4_intact.mhd'  # pilot5
# file_inst = 'ICOTEC_S130672_L5_implants_XCTres.mhd'  # segmented screw scan  pilot3
# file_inst = 'Icotec_S130672_L4_mask.mhd'  # pilot4
file_inst = 'Icotec_S130684_L4_mask.mhd'  # pilot5

# Specify file locations
# ---sample_no, sample = sample_code.split('_')  # sample name and number
path_project = '/home/biomech/Documents/01_Icotec/'  # General project folder
path_ct = path_project + '01_Experiments/02_Scans/' + sample_code + '/04_Registered/'  # Folder of CT data
path_fea = path_project + '02_FEA/98_Pilots/' + sample_code + '/' + model_code + '/'  # Folder of FEA files
# Include general path for mesh/template later
info = sample_code + ExpScrew + '_info.txt'  # .txt file containing info about landmarks. Specific structure

# # # # # Input # # # # #
# Input FEA data
Input = {}
Input['Project_Folder'] = path_project
Input['FEA_loc'] = path_fea  # path to FEA files
Input['Model_Code'] = model_code  # model code (FEA model name)
Input['Screw'] = ExpScrew  # which site of scan will be processed
Input['Resolution'] = 0.0607  # scan resolution, should be 0.0607 for HR-pQCT

Input['Load_mode'] = 'f'  # no effect 'd' or 'f', displacement or force controlled
Input['F_dir'] = '-'  # force direction (negative corresponds to experiment, positive = inverse). Input: '-' or '+'
Input['F_max'] = 150  # peak load
Input['F_min'] = 10  # valley load, not included yet!!!!!!!!!!!!!!!!!!!!!!
Input['Cycles'] = 2  # how many cycles to simulate, not included yet!!!!!!!!!!!!!!!!!!!!!!
Input['d_dir'] = '-'  # no effect
Input['d_max'] = '15'  # no effect
Input['Friction'] = 0.2  # friction between screw and bone, not included yet!!!!!!!!!!!!!!!!!!!!!!
Input['ElType'] = 'Quadratic'  # element type quadratic or linear C3D8 = linear, not included yet!!!!!!!!!!!!!!!!!!!!!!
Input['Mapping_Diameter'] = 2.5  # diameter of sphere for mapping, in mm. should be larger than element size

# parametrise screw length? --> different mesh sizes // ROI? // step size?

# Submit on cortex or ubelix?
# -> Other parameters as e-mail, runtime, memory etc. can be changed in the template file
Input['Submit'] = 'cortex'
mapp.write_submit(Input)

# Write output images? segmented image and mask, for visual check
write_output = 0

# Write mesh input file
mapp.write_mesh(Input)  # Original input file, path for mesh.inp

# Load data screw information from .txt file
# About file:
# Origin is where screw enters bone, will be origin of bone-mesh
# Rotation axis positive x axis: Point on rotation axis (where bone is fixed in testing machine), positive x-coordinate
# Rotation axis negative x axis: Point on rotation axis (where bone is fixed in testing machine), negative x-coordinate
imD = rR.load_itk(path_ct + file_inst)  # screw image
with open(path_ct + info, 'r') as f:
    content = f.read()
ori = content.split('Origin: ')[1].split('\n')[0]
ori = np.array([int(ori.split(' ')[0]), int(ori.split(' ')[1]), int(ori.split(' ')[2])])
p1P = content.split('positive x axis: ')[1].split('\n')[0]
p1P = np.array([int(p1P.split(' ')[0]), int(p1P.split(' ')[1]), int(p1P.split(' ')[2])])
p1N = content.split('negative x axis: ')[1].split('\n')[0]
p1N = np.array([int(p1N.split(' ')[0]), int(p1N.split(' ')[1]), int(p1N.split(' ')[2])])

# # # # # Input # # # # #
# Define screw vector
v3 = []
if ExpScrew == 'T':
    # lineT = rR.axis3D(imD[0], 670, 1100, 'x')  # if starting from screw tip towards head, add (-) to v3  pilot3
    # lineT = rR.axis3D(imD[0], 915, 1485, 'x')  # pilot4
    lineT = rR.axis3D(imD[0], 832, 1464, 'x')  # pilot5
    v3 = -lineT.vector  # z-axis = screw axis, found by function. Specify +/-!!
elif ExpScrew == 'P':
    # lineT = rR.axis3D(imD[0], 0, 600, 'x')  # if starting from screw tip towards head, add (-) to v3  pilot3
    # lineT = rR.axis3D(imD[0], 483, 772, 'x')  # pilot4
    lineT = rR.axis3D(imD[0], 418, 667, 'x')  # pilot5
    v3 = lineT.vector  # z-axis = screw axis, found by function. Specify +/-!!
del imD

# Compute matrix
v2 = np.cross(v3, p1P-p1N)/np.linalg.norm(np.cross(v3, p1P-p1N))  # y-axis
v1 = np.cross(v2, v3)  # x-axis
M = np.vstack((np.append(v1, ori[0]), np.append(v2, ori[1]), np.append(v3, ori[2]), np.array([0, 0, 0, 1])))
print(M)

# Check coordinate system
# z-axis (blue) should point from screw head (origin) towards screw tip (dot)
# y-axis (green) should be perpendicular to screw-rotAxis plane and point upwards (against loading direction)
# x-axis (red): right-handed coordinate system

#plt.figure()
#ax = plt.axes(projection='3d')
#fact = 100
#for i in range(0, 5):
#    ax.scatter3D(ori[0]+i*fact*v1[0], ori[1]+i*fact*v1[1], ori[2]+i*fact*v1[2], c='r', alpha=1)
#    ax.scatter3D(ori[0]+i*fact*v2[0], ori[1]+i*fact*v2[1], ori[2]+i*fact*v2[2], c='g', alpha=1)
#    ax.scatter3D(ori[0]+i*fact*v3[0], ori[1]+i*fact*v3[1], ori[2]+i*fact*v3[2], c='b', alpha=1)
#ax.scatter3D([0, 1500], [0, 1500], [0, 1500], alpha=0)
#ax.scatter3D(p1P[0], p1P[1], p1P[2], c='k')
#ax.scatter3D(p1N[0], p1N[1], p1N[2], c='k')
#ax.scatter3D(ori[0], ori[1], ori[2], c='c')
# ax.scatter3D(879, 486, 799, c='m')  # tip 621 476 806
# ax.scatter3D(879, 100, 799, c='r')  # force direction (negative!)
#ax.set_xlabel('x')
#ax.set_ylabel('y')
#ax.set_zlabel('z')
#plt.show()

bone = {}

bone = mapp.readInpBoneDummy(bone, Input)  # Read bone mesh from abaqus. Read elements, nodes

bone = mapp.boneMeshMask(bone, Input)  # Create mask from abaqus bone mesh

bone = mapp.load_BVTVdata(bone, path_ct + file_bone)

# Read mask
imMask = sitk.ReadImage(Input['FEA_loc'] + Input['Model_Code'] + Input['Screw'] + '_mask.mhd')
imMask_np = np.transpose(sitk.GetArrayFromImage(imMask), [2, 1, 0])
bone["BVTVscaled"].SetOrigin([0, 0, 0])
# Define rotation and translation
[theta1, theta2, theta3] = rR.rotation_angles_from_matrix(M[:3, :3], 'zyx')
theta2 = -theta2 + np.pi  # +pi because arc-tan has no unique solution. WHY NEGATIVE, because mask[:, ::, :] reshaped?
center = np.array([imMask_np.shape[0]/2, imMask_np.shape[1]/2, 0]) * imMask.GetSpacing()
trans = M[:3, 3] * imMask.GetSpacing() - center

# Write transformation file
f = open(path_fea + 'transformation_' + model_code + ExpScrew + '.tfm', "w")
f.write(
    "#Insight Transform File V1.0\n"
    "#Transform 0\n"
    "Transform: CompositeTransform_double_3_3\n"
    "#Transform 1\n"
    "Transform: Euler3DTransform_double_3_3\n"
    "Parameters:  " + f'{theta1}' + " " + f'{theta2}' + " " + f'{theta3}'
    + " " + f'{trans[0]}' + " " + f'{trans[1]}' + " " + f'{trans[2]}' + "\n"  # transformation
    "FixedParameters: " + f'{center[0]}' + " " + f'{center[1]}' + " " + f'{center[2]}' + " 0\n")  # CoR
f.close()

# Apply transformation to mask using the inverse transformation
bone['Transform'] = sitk.ReadTransform(path_fea + 'transformation_' + model_code + ExpScrew + '.tfm')
bone['Transform_inv'] = bone['Transform'].GetInverse()
imMask_trans = sitk.Resample(imMask, bone['BVTVscaled'], bone['Transform_inv'], sitk.sitkNearestNeighbor, 0.0,
                             bone['BVTVscaled'].GetPixelID())
# Delete some files / variables to save memory
os.remove(path_fea + model_code + ExpScrew + '_mask.mhd')
os.remove(path_fea + model_code + ExpScrew + '_mask.raw')
if write_output:
    img_seg = sitk.GetImageFromArray(np.transpose(bone['BVTVscaled'], [2, 1, 0]))
    img_seg.SetOrigin(bone["BVTVscaled"].GetOrigin())
    img_seg.SetSpacing(bone["BVTVscaled"].GetSpacing())
    sitk.WriteImage(img_seg, path_fea + sample_code + '_seg.mhd')
    imMask_trans.SetOrigin(bone["BVTVscaled"].GetOrigin())
    imMask_trans.SetSpacing(bone["BVTVscaled"].GetSpacing())
    sitk.WriteImage(imMask_trans, path_fea + model_code + '_' + ExpScrew + '_maskTrans.mhd')
    del img_seg
# del bone['GreyImage']

# BVTV segmentation / calibration HR-pQCT to uCT
bone['MASK_array_T'] = np.transpose(sitk.GetArrayFromImage(imMask_trans), [2, 1, 0])
# scaling factor/intercept from Schenk et al. 2022, has to be discussed w Ph

mapp.HFE_mapping_trans(bone, Input)

sliceNo = 460
plt.imshow(bone['BVTVscaled'][:, sliceNo, :] + bone['MASK_array_T'][:, sliceNo, :])

plt.show()

# Write final input file
mapp.HFE_inp_creator(Input)
if not write_output:
    os.remove(Input['FEA_loc'] + Input['Model_Code'] + Input['Screw'] + '_elsets.inp')
    os.remove(Input['FEA_loc'] + Input['Model_Code'] + Input['Screw'] + '_materials.inp')
    os.remove(Input['FEA_loc'] + Input['Model_Code'] + '_mesh.inp')
    os.remove(Input['FEA_loc'] + 'transformation_' + Input['Model_Code'] + Input['Screw'] + '.tfm')

print('Execution time: ' + str(int((time.time()-t1)/60)) + ' m ' + str(round(np.mod(time.time()-t1, 60), 1)) + ' s.')
