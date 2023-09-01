import mapping_noRot as mappNR
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import SimpleITK as sitk
import time
import os


def mapping(sample, mod, fric_):
    t1 = time.time()
    # # # # # Input # # # # #
    # FE version
    models = ['00_L50_S50_D30', '01_L50_S50_D35', '02_L50_S50_D40', '03_L50_S50_D45', '04_L50_S50_D27',  # 0, 1, 2, 3, 4
              '10_L50_S00_D30', '11_L50_S00_D35', '12_L50_S00_D40', '13_L50_S00_D45', '14_L50_S00_D27',  # 5, 6, 7, 8, 9
              '15_L50_S00_D44',  # 10
              '31_L50_S50_D35', '43_L50_S00_D45',  # 11, 12
              '50_L50_S00_D30', '55_L50_S00_D30',  # 13, 14
              '63_L50_S50_D45',  # 15
              '65_L50_S50_D45', '66_L50_S50_D45',  # 16, 17
              '74_L50_S50_D45', '75_L50_S50_D45', '76_L50_S50_D45', '77_L50_S50_D45_expl',  # 18, 19, 20, 21
              '80_L50_S50_D45', '81_L50_S50_D45', '82_L50_S50_D45', '83_L50_S50_D45',  # 22, 23, 24, 25
              '85_L50_S50_D45', '86_L50_S50_D45', '87_L50_S50_D45', '88_L50_S50_D45',  # 26, 27, 28, 29
              '94_OSTP']  # 30

    model_code = models[mod]  # FEA model name
    print('Model: ' + str(model_code))
    sample_code = sample  # Experimental sample. contains number and sample name

    # Specify file locations
    path_project = '/home/biomech/Documents/01_Icotec/'  # General project folder
    path_ct = path_project + '01_Experiments/02_Scans/' + sample_code + '/04_Registered/'  # Folder of CT data
    path_fea = path_project + '02_FEA/01_MainStudy/' + sample_code + '/' + model_code + '/'  # Folder of FEA files
    file_bone = [filename for filename in os.listdir(path_ct + '/') if
                 filename.endswith('image.mhd') and str(sample_code) in filename][0]
    print(file_bone)

    # # # # # Input # # # # #
    # Input FEA data
    Input = {}
    Input['Project_Folder'] = path_project
    Input['FEA_loc'] = path_fea  # path to FEA files
    Input['Model_Code'] = model_code  # model code (FEA model name)
    Input['Resolution'] = 0.0606995  # scan resolution, should be 0.0606995 for HR-pQCT II

    Input['d_dir'] = '-'  # displ direction (negative corresponds to experiment, positive = inverse). Input: '-' or '+'
    Input['d_max'] = 1  # peak displ

    Input['Friction'] = fric_  # friction between screw and bone
    Input['Mapping_Diameter'] = 2  # diameter of sphere for mapping, in mm. should be larger than element size

    Input['YM_peek'] = str(25000)  # young's modulus peek screw
    Input['v_peek'] = str(0.3)  # poisson ratio peek screw
    Input['YM_titan'] = str(100000)  # young's modulus titanium screw
    Input['v_titan'] = str(0.3)  # poisson ratio titanium screw

    # Check if folder exists, otherwise create it
    isExist = os.path.exists(Input['FEA_loc'])
    if not isExist:
        os.makedirs(Input['FEA_loc'])
        print('New directory created: ' + Input['FEA_loc'])

    # Submit on cortex or ubelix?
    # -> Other parameters as e-mail, runtime, memory etc. can be changed in the template file
    Input['Submit'] = 'ubelix'
    mappNR.write_submit(Input)

    # Write output images? segmented image and mask, for visual check
    write_output = 0

    # Write mesh input file
    mappNR.write_mesh(Input)  # Original input file, path for mesh.inp

    # Bone dictionary will contain all information about the CT image
    bone = {}
    bone = mappNR.readInpBoneDummy(bone, Input)  # Read bone mesh from abaqus. Read elements, nodes
    bone = mappNR.load_BVTVdata(bone, path_ct + file_bone)
    bone = mappNR.boneMeshMask(bone, Input)  # Create mask from abaqus bone mesh

    # Read mask
    imMask = sitk.ReadImage(Input['FEA_loc'] + Input['Model_Code'] + '_mask.mhd')
    imMask_np = np.transpose(sitk.GetArrayFromImage(imMask), [2, 1, 0])
    orM = np.array(imMask.GetOrigin())
    orB = np.array(bone['GreyImage'].GetOrigin())
    insBefore = np.rint(abs(orB - orM) / Input['Resolution']).astype(int)
    bone['insBefore'] = insBefore
    print('insBefore: ' + str(insBefore))
    dimMask = np.array(imMask_np.shape)
    dimBone = np.array(bone['GreyImage'].GetSize())
    insAfter = (dimBone - dimMask - insBefore).astype(int)
    bone['insAfter'] = insAfter
    print('insAfter: ' + str(insAfter))
    print('mask dimension: ' + str(dimMask))
    print('bone dimension: ' + str(dimBone))
    imMask_np_corr = imMask_np
    imMask_np_corr = np.append(np.insert(imMask_np_corr, 0, np.zeros((insBefore[0], imMask_np_corr.shape[1],
                                                                      imMask_np_corr.shape[2])), 0),
                               np.zeros((insAfter[0], imMask_np_corr.shape[1], imMask_np_corr.shape[2])), 0)

    imMask_np_corr = np.append(
        np.insert(imMask_np_corr, 0, np.zeros((insBefore[1], imMask_np_corr.shape[0], imMask_np_corr.shape[2])), 1),
        np.zeros((imMask_np_corr.shape[0], insAfter[1], imMask_np_corr.shape[2])), 1)

    imMask_np_corr = np.append(
        np.insert(imMask_np_corr, 0, np.zeros((insBefore[2], imMask_np_corr.shape[0], imMask_np_corr.shape[1])), 2),
        np.zeros((imMask_np_corr.shape[0], imMask_np_corr.shape[1], insAfter[2])), 2)
    print('new mask dimension: ' + str(np.array(imMask_np_corr.shape)))

    plt.figure()
    plt.imshow(imMask_np[int(dimMask[0] / 2), :, :], cmap=cm.RdBu_r)
    plt.show()
    plt.savefig(Input['FEA_loc'] + 'mappingControlPlot.png')
    plt.close()

    bone['MASK'] = imMask_np_corr
    shape_mask_x = np.load('/home/biomech/DATA/01_Icotec/01_Experiments/02_Scans/BVTV/' + sample_code + '_mask_x.npy')
    shape_mask_y = np.load('/home/biomech/DATA/01_Icotec/01_Experiments/02_Scans/BVTV/' + sample_code + '_mask_y.npy')
    shape_mask_z = np.load('/home/biomech/DATA/01_Icotec/01_Experiments/02_Scans/BVTV/' + sample_code + '_mask_z.npy')
    bone['shape_mask'] = ((shape_mask_x + shape_mask_y + shape_mask_z) > 0).astype(int)
    del shape_mask_x
    del shape_mask_y
    del shape_mask_z
    if write_output:
        img_screw = sitk.GetImageFromArray(np.transpose(bone['MASK'], [2, 1, 0]))
        img_screw.SetOrigin(bone['GreyImage'].GetOrigin())
        img_screw.SetSpacing(bone['GreyImage'].GetSpacing())
        sitk.WriteImage(img_screw, path_fea + sample_code + '_screw.mhd')
        print('Screw saved.')
        img_seg = sitk.GetImageFromArray(np.transpose(bone['BVTVscaled'], [2, 1, 0]))
        img_seg.SetOrigin(bone['GreyImage'].GetOrigin())
        img_seg.SetSpacing(bone['GreyImage'].GetSpacing())
        sitk.WriteImage(img_seg, path_fea + sample_code + '_seg.mhd')
        print('Segmented image saved.')
        del img_seg

    mappNR.HFE_mapping_trans(bone, Input)

    # Write final input file
    mappNR.HFE_inp_creator(Input)
    if not write_output:
        os.remove(Input['FEA_loc'] + Input['Model_Code'] + '_elsets.inp')
        os.remove(Input['FEA_loc'] + Input['Model_Code'] + '_materials.inp')
        os.remove(Input['FEA_loc'] + Input['Model_Code'] + '_mesh.inp')
        os.remove(path_fea + model_code + '_mask.mhd')
        os.remove(path_fea + model_code + '_mask.raw')

        tRun = time.time() - t1
        if tRun >= 3600:
            print('Execution time: ' + str(int(tRun / 3600)) + ' h ' + str(int(np.mod(tRun, 3600) / 60)) + ' min ' +
                  str(round(np.mod(tRun, 60), 1)) + ' sec.\n')
        elif tRun >= 60:
            print('Execution time: ' + str(int(tRun / 60)) + ' min ' + str(round(np.mod(tRun, 60), 1)) + ' sec.\n')
        else:
            print('Execution time: ' + str(round(tRun, 1)) + ' sec.\n')


sample_list = open('/home/biomech/Documents/01_Icotec/Specimens.txt', 'r').read().splitlines()

peek_samples = [2, 5, 7, 8, 10, 13, 15, 16, 18, 21, 23, 24, 26, 29, 31, 32]  # without 0
ti_samples = [3, 4, 6, 9, 11, 12, 14, 17, 19, 20, 22, 25, 27, 28, 30, 33]  # without 1

for i in ti_samples:  # ti_samples:  # range(2, len(sample_list)):  # range(12, 20):  # len(sample_list)):

    print(sample_list[i])
    mapping(sample_list[i], 16, 0.2)  # samples, model, friction
