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
    theta1_ = []
    theta2_ = []
    theta3_ = []

    if order == 'xzx':
        theta1_ = np.arctan(r31 / r21)
        theta2_ = np.arctan(r21 / (r11 * np.cos(theta1_)))
        theta3_ = np.arctan(-r13 / r12)

    elif order == 'xyx':
        theta1_ = np.arctan(-r21 / r31)
        theta2_ = np.arctan(-r31 / (r11 * np.cos(theta1_)))
        theta3_ = np.arctan(r12 / r13)

    elif order == 'yxy':
        theta1_ = np.arctan(r12 / r32)
        theta2_ = np.arctan(r32 / (r22 * np.cos(theta1_)))
        theta3_ = np.arctan(-r21 / r23)

    elif order == 'yzy':
        theta1_ = np.arctan(-r32 / r12)
        theta2_ = np.arctan(-r12 / (r22 * np.cos(theta1_)))
        theta3_ = np.arctan(r23 / r21)

    elif order == 'zyz':
        theta1_ = np.arctan(r23 / r13)
        theta2_ = np.arctan(r13 / (r33 * np.cos(theta1_)))
        theta3_ = np.arctan(-r32 / r31)

    elif order == 'zxz':
        theta1_ = np.arctan(-r13 / r23)
        theta2_ = np.arctan(-r23 / (r33 * np.cos(theta1_)))
        theta3_ = np.arctan(r31 / r32)

    elif order == 'xzy':
        theta1_ = np.arctan(r32 / r22)
        theta2_ = np.arctan(-r12 * np.cos(theta1_) / r22)
        theta3_ = np.arctan(r13 / r11)

    elif order == 'xyz':
        theta1_ = np.arctan(-r23 / r33)
        theta2_ = np.arctan(r13 * np.cos(theta1_) / r33)
        theta3_ = np.arctan(-r12 / r11)

    elif order == 'yxz':
        theta1_ = np.arctan(r13 / r33)
        theta2_ = np.arctan(-r23 * np.cos(theta1_) / r33)
        theta3_ = np.arctan(r21 / r22)

    elif order == 'yzx':
        theta1_ = np.arctan(-r31 / r11)
        theta2_ = np.arctan(r21 * np.cos(theta1_) / r11)
        theta3_ = np.arctan(-r23 / r22)

    elif order == 'zyx':
        theta1_ = np.arctan(r21 / r11)
        theta2_ = np.arctan(-r31 * np.cos(theta1_) / r11)
        theta3_ = np.arctan(r32 / r33)

    elif order == 'zxy':
        theta1_ = np.arctan(-r12 / r22)
        theta2_ = np.arctan(r32 * np.cos(theta1_) / r22)
        theta3_ = np.arctan(-r31 / r33)

    # theta1_ = theta1_ * 180 / np.pi
    # theta2_ = theta2_ * 180 / np.pi
    # theta3_ = theta3_ * 180 / np.pi

    return theta1_, theta2_, theta3_


# test_transform, original
rotTransl_matrix = np.array([[-0.79327354, -0.01318472,  0.60872264, 1130],
                             [ 0.02211240,  0.99848211,  0.05044315,  453],
                             [-0.60846375,  0.05347554, -0.79177789,  764],
                             [ 0.00000000,  0.00000000,  0.00000000,    1]])

Test_trans_4_4 = rotTransl_matrix
mask_path = '/home/biomech/Downloads/'

img_mask = sitk.ReadImage(mask_path+'test.mhd')

img_mask_np = np.transpose(sitk.GetArrayFromImage(img_mask), [2, 1, 0])

greypath = '/home/biomech/Documents/01_Icotec/01_Experiments/02_Scans/Pilot3/04_Registered/'

img_grey = sitk.ReadImage(greypath+'XCT_Icotec_S130672_L5_intact_planned.mhd')
img_grey.SetOrigin([0, 0, 0])

[theta1, theta2, theta3] = rotation_angles_from_matrix(Test_trans_4_4[:3, :3], 'zyx')
print(theta1)

Center = np.array([img_mask_np.shape[0]/2, img_mask_np.shape[1]/2, 0])*img_mask.GetSpacing()
print(Center)
angle = [theta1, theta2, theta3]

trans = Test_trans_4_4[:3, 3]*img_mask.GetSpacing()-Center

print(angle)

f = open(mask_path + 'Test_transformation.tfm', "w")
f.write(
    "#Insight Transform File V1.0\n"
    "#Transform 0\n"
    "Transform: CompositeTransform_double_3_3\n"
    "#Transform 1\n"
    "Transform: Euler3DTransform_double_3_3\n"
#    "Parameters:  " + f'{theta1}' + " " + f'{theta2}' + " " + f'{theta3}' + f'{trans[0]}' + " " + f'{trans[1]}' + " " + f'{trans[2]}' + "\n"
    "Parameters:  " + f'{theta1}' + " " + f'{theta2}' + " " + f'{theta3}' + " 0 0 0\n"
    "FixedParameters: " + f'{Center[0]}' + " " + f'{Center[1]}' + " " + f'{Center[2]}' + " 0\n")  # Center of rotation
f.close()

transform = sitk.ReadTransform(mask_path + 'Test_transformation.tfm')
transform_inv = transform.GetInverse()
print(transform)
print(transform_inv)

img_mask_trans = sitk.Resample(img_mask, img_grey, transform, sitk.sitkNearestNeighbor,
                                             0.0, img_grey.GetPixelID())
sitk.WriteImage(img_mask_trans, mask_path + 'Test_mask_trans.mhd')
sitk.WriteImage(img_grey, mask_path + 'test_CT.mhd')

# "Parameters: " +f'{theta1}' + " " +f'{theta2}' + " " +f'{theta3}'
# + " " +f'{trans[0]}' + " " +f'{trans[1]}' + " " +f'{trans[2]}' + "\n"
