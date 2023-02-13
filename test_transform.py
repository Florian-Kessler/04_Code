import SimpleITK as sitk
import numpy as np


def rotation_angles_from_matrix(matrix, order):
    """
    input
        matrix = 3x3 rotation matrix (numpy array)
        order(str) = rotation order of x, y, z : e.g, rotation XZY -- 'xzy'
    output
        theta1, theta2, theta3 = rotation angles in rotation order
    """
    r11, r12, r13 = matrix[0]
    r21, r22, r23 = matrix[1]
    r31, r32, r33 = matrix[2]
    theta1 = []
    theta2 = []
    theta3 = []

    if order == 'xzx':
        theta1 = np.arctan(r31 / r21)
        theta2 = np.arctan(r21 / (r11 * np.cos(theta1)))
        theta3 = np.arctan(-r13 / r12)

    elif order == 'xyx':
        theta1 = np.arctan(-r21 / r31)
        theta2 = np.arctan(-r31 / (r11 * np.cos(theta1)))
        theta3 = np.arctan(r12 / r13)

    elif order == 'yxy':
        theta1 = np.arctan(r12 / r32)
        theta2 = np.arctan(r32 / (r22 * np.cos(theta1)))
        theta3 = np.arctan(-r21 / r23)

    elif order == 'yzy':
        theta1 = np.arctan(-r32 / r12)
        theta2 = np.arctan(-r12 / (r22 * np.cos(theta1)))
        theta3 = np.arctan(r23 / r21)

    elif order == 'zyz':
        theta1 = np.arctan(r23 / r13)
        theta2 = np.arctan(r13 / (r33 * np.cos(theta1)))
        theta3 = np.arctan(-r32 / r31)

    elif order == 'zxz':
        theta1 = np.arctan(-r13 / r23)
        theta2 = np.arctan(-r23 / (r33 * np.cos(theta1)))
        theta3 = np.arctan(r31 / r32)

    elif order == 'xzy':
        theta1 = np.arctan(r32 / r22)
        theta2 = np.arctan(-r12 * np.cos(theta1) / r22)
        theta3 = np.arctan(r13 / r11)

    elif order == 'xyz':
        theta1 = np.arctan(-r23 / r33)
        theta2 = np.arctan(r13 * np.cos(theta1) / r33)
        theta3 = np.arctan(-r12 / r11)

    elif order == 'yxz':
        theta1 = np.arctan(r13 / r33)
        theta2 = np.arctan(-r23 * np.cos(theta1) / r33)
        theta3 = np.arctan(r21 / r22)

    elif order == 'yzx':
        theta1 = np.arctan(-r31 / r11)
        theta2 = np.arctan(r21 * np.cos(theta1) / r11)
        theta3 = np.arctan(-r23 / r22)

    elif order == 'zyx':
        theta1 = np.arctan(r21 / r11)
        theta2 = np.arctan(-r31 * np.cos(theta1) / r11)
        theta3 = np.arctan(r32 / r33)

    elif order == 'zxy':
        theta1 = np.arctan(-r12 / r22)
        theta2 = np.arctan(r32 * np.cos(theta1) / r22)
        theta3 = np.arctan(-r31 / r33)

    #theta1 = theta1 * 180 / np.pi
    #theta2 = theta2 * 180 / np.pi
    #theta3 = theta3 * 180 / np.pi

    return theta1, theta2, theta3


# test_transform
Test_trans_4_4 = np.array([[-7.91777893e-01,  5.34755369e-02, -6.08463749e-01, -1.33535091e+03],
                   [ 2.25184329e-02,  9.98038607e-01,  5.84111263e-02,  5.22183419e+02],
                   [ 6.10393879e-01,  3.25469884e-02, -7.91429091e-01,  9.98370437e+01],
                   [ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  1.00000000e+00]])

mask_path='/home/biomech/Downloads/'

img_mask = sitk.ReadImage(mask_path+'test.mhd')

img_mask_np=np.transpose(sitk.GetArrayFromImage(img_mask),[2,1,0])

greypath = '/home/biomech/Documents/01_Icotec/01_Experiments/02_Scans/Pilot3/04_Registered/'

img_grey = sitk.ReadImage(greypath+'XCT_Icotec_S130672_L5_intact_planned.mhd')
img_grey.SetOrigin([0,0,0])

[theta1, theta2, theta3] = rotation_angles_from_matrix(Test_trans_4_4[:3, :3], 'zyx')
print(theta1)

Center=np.array(img_mask_np.shape)/2*img_mask.GetSpacing()
print(Center)
winkel= [theta1, theta2, theta3]

trans= Test_trans_4_4[:3,3]*img_mask.GetSpacing()-Center

print(winkel)

f = open(mask_path + 'Test_transformation.tfm', "w")
f.write(
    "#Insight Transform File V1.0\n"
    "#Transform 0\n"
    "Transform: CompositeTransform_double_3_3\n"
    "#Transform 1\n"
    "Transform: Euler3DTransform_double_3_3\n"
    "Parameters:  "+f'{theta1}' + " " +f'{theta2}' + " " +f'{theta3}' + " 10 10 10\n"
    "FixedParameters: "+f'{Center[0]}' + " " +f'{Center[1]}' + " " +f'{Center[2]}' + " 0\n")
f.close()

transform = sitk.ReadTransform(mask_path + 'Test_transformation.tfm')
transform_inv=transform.GetInverse()
print(transform)
print(transform_inv)

img_mask_trans = sitk.Resample(img_mask, img_grey, transform_inv, sitk.sitkNearestNeighbor,
                                             0.0, img_grey.GetPixelID())
sitk.WriteImage(img_mask_trans,mask_path+'Test_mask_trans.mhd')

#"Parameters: " +f'{theta1}' + " " +f'{theta2}' + " " +f'{theta3}' + " " +f'{trans[0]}' + " " +f'{trans[1]}' + " " +f'{trans[2]}' + "\n"