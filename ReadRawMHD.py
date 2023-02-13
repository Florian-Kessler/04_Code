import SimpleITK as sitk
from skimage.measure import label
import numpy as np
import matplotlib.pyplot as plt
import time
from scipy.ndimage import rotate
from colorama import Fore, Style
from scipy import stats
from skspatial.objects import Line
import imutils


def load_itk(filename):
    # Load scan, filename = path of .mhd (.raw in same directory)

    # Reads the image using SimpleITK
    itkimage = sitk.ReadImage(filename)

    # Convert the image to a  numpy array first and then shuffle the dimensions to get axis in the order z,y,x
    ct_scan = sitk.GetArrayFromImage(itkimage)

    # Transform image from z,y,x to x,y,z
    ct_scan = np.transpose(ct_scan, [2, 1, 0])

    # Read the origin of the ct_scan, will be used to convert the coordinates from world to voxel and vice versa.
    origin = np.array(list(reversed(itkimage.GetOrigin())))

    # Read the spacing along each dimension
    spacing = np.array(list(reversed(itkimage.GetSpacing())))

    return ct_scan, origin, spacing


def plot_grey(image, axis, pos):
    plt.figure()
    if axis == 0:
        plt.imshow(image[pos, :, :], cmap='gray')
    elif axis == 1:
        plt.imshow(image[:, pos, :], cmap='gray')
    elif axis == 2:
        plt.imshow(image[:, :, pos], cmap='gray')
    else:
        print(Fore.CYAN + 'Error: Enter valid axis number (0, 1, 2).')
        print(Style.RESET_ALL)


def plot_HU(image, axis, pos, level, window):
    """Plot image in defined HU window
    Input: Image array (2d), level of window centre, window size"""
    ma = level + window / 2
    mi = level - window / 2
    slic = image.copy()
    slic[slic < mi] = mi
    slic[slic > ma] = ma

    plt.figure()
    if axis == 0:
        plt.imshow(slic[pos, :, :], cmap='gray')
    elif axis == 1:
        plt.imshow(slic[:, pos, :], cmap='gray')
    elif axis == 2:
        plt.imshow(slic[:, :, pos], cmap='gray')
    else:
        print(Fore.CYAN + 'Error: Enter valid axis number (0, 1, 2).')
        print(Style.RESET_ALL)
    return slice


def rot(image, centre_xy, angle):
    """Rotates an image around a rotation centre.
    Input: Image array (2d), rotation centre, angle in degrees
    Output: Rotated image array (2d), new rotation centre"""
    im_rot = rotate(image, angle)
    org_center = (np.array(image.shape[:2][::-1]) - 1) / 2.
    rot_center = (np.array(im_rot.shape[:2][::-1]) - 1) / 2.
    org = centre_xy - org_center
    a_ = np.deg2rad(angle)
    new = np.array([org[0] * np.cos(a_) + org[1] * np.sin(a_),
                    - org[0] * np.sin(a_) + org[1] * np.cos(a_)])
    return im_rot, new + rot_center


def xray(image, axis):
    plt.figure()
    plt.imshow(np.sum(image, axis=axis), cmap='gray')
    plt.show()


def getLargestCC(segm_im):
    labels = label(segm_im)
    assert (labels.max() != 0)  # assume at least 1 CC
    largestCC = labels == np.argmax(np.bincount(labels.flat)[1:]) + 1
    return largestCC


def countOccurrence(arr, x_):
    """Counts the occurrence of x in array arr"""
    res = 0
    arrF = arr.flatten()
    for j in range(len(arrF)):
        if x_ == arrF[j]:
            res += 1
    return res


def islandFinder(img):
    """Finds connected islands (connectivity = 1, 2, ..., size(img)).
    Returns an array containing the boolean array of each island and
    the size of each island."""
    [imglabeled, island_count] = label(img, background=0, return_num=True, connectivity=1)
    islands = np.zeros(np.append(island_count, np.append(1, imglabeled.shape)))
    island_size = np.zeros(island_count)
    for j in range(island_count):
        islands[j] = np.array(imglabeled == j + 1).astype(int)
        island_size[j] = countOccurrence(imglabeled, j + 1)
    return islands, island_size


def axis2D3D(img, start, end, test, plot):
    # centre of gravity for each slice of the screw
    xs = []
    ys = []
    zs = []
    for i_ in range(start, end):
        # Find slices where screw appears
        loc_ = np.array(np.where(img[:, :, i_] == 1))
        if loc_.shape[1]:
            # Find centre of gravity for slice i
            cog = [np.mean(loc_[0]), np.mean(loc_[1])]
            xs.append(cog[0])
            ys.append(cog[1])
            zs.append(i_)
            if test == 1:
                # Print each slice to visually check centre of gravity
                print('Screw found on layer ' + str(i_))
                plt.figure()
                plt.scatter(cog[1], cog[0], c='r')
                plt.imshow(img[:, :, i_], cmap='gray')
                plt.pause(2)
                plt.close()
        else:
            if test == 1:
                print('No screw on layer ' + str(i_))
    # Linear regression in every plane (correct?)
    [slope_xy, intercept_xy, _, _, _] = stats.linregress(xs, ys)
    [slope_xz, intercept_xz, _, _, _] = stats.linregress(xs, zs)
    [_, _, _, _, _] = stats.linregress(ys, zs)
    # Minima and maxima to define direction of screw axis vector
    x_ax = [np.min(xs), np.max(xs)]
    y_ax = np.multiply(x_ax, slope_xy) + intercept_xy
    z_ax = np.multiply(x_ax, slope_xz) + intercept_xz
    # Return start and end points
    points = np.zeros((2, 3))
    points[0, :] = [x_ax[0], y_ax[0], z_ax[0]]
    points[1, :] = [x_ax[-1], y_ax[-1], z_ax[-1]]
    if plot == 1:
        plt.figure()
        ax1 = plt.axes(projection='3d')
        ax1.scatter3D(xs, ys, zs, c='b')
        ax1.plot3D(x_ax, y_ax, z_ax, c='r')
    return points


def seg(img, th, th2):
    img[img < th] = 0
    img[img > th2] = 0
    img_seg = np.array(img > 0).astype(int)
    return img_seg


def axis3D(img, start, end, test, plot):
    # centre of gravity for each slice of the screw
    xs = []
    ys = []
    zs = []
    if plot == 1:
        plt.figure(20)
        ax1 = plt.axes(projection='3d')
        ax1.scatter3D([-img.shape[0], img.shape[0]], [-img.shape[1], img.shape[1]], [start, end], alpha=0)
    for i_ in range(start, end):
        # Find slices where screw appears
        loc_ = np.array(np.where(img[:, :, i_] >= 1))
        if loc_.shape[1]:
            # Find centre of gravity for slice i
            cog = [np.mean(loc_[0]), np.mean(loc_[1])]
            xs.append(cog[0])
            ys.append(cog[1])
            zs.append(i_)
            if test == 1:
                # Print each slice to visually check centre of gravity
                print('Screw found on layer ' + str(i_))
                plt.figure(19)
                plt.scatter(cog[1], cog[0], c='r')
                plt.imshow(img[:, :, i_], cmap='gray')
                plt.pause(2)
                plt.close()
            if plot == 1:
                # ax1 = plt.axes(projection='3d')
                # ax1.scatter3D(xs, ys, zs, c='b', alpha=0.4)
                ax1.scatter3D(loc_[0], loc_[1], i_, c='k', alpha=0.1)
                plt.pause(0.01)
        else:
            if test == 1:
                print('No screw on layer ' + str(i_))
    line_fit = Line.best_fit(np.array([xs, ys, zs]).T)
    if plot == 1:
        t_ = 100
        plt.figure(21)
        ax2 = plt.axes(projection='3d')
        ax2.plot3D(
            [line_fit.point[0] - line_fit.vector[0] * t_, line_fit.point[0] + line_fit.vector[0] * t_],
            [line_fit.point[1] - line_fit.vector[1] * t_, line_fit.point[1] + line_fit.vector[1] * t_],
            [line_fit.point[2] - line_fit.vector[2] * t_, line_fit.point[2] + line_fit.vector[2] * t_],
            c='r')
        # min_ = np.min([line_fit.point - line_fit.vector * t_, line_fit.point + line_fit.vector * t_])
        # max_ = np.max([line_fit.point - line_fit.vector * t_, line_fit.point + line_fit.vector * t_])
        ax2.scatter3D([-img.shape[0], img.shape[0]], [-img.shape[1], img.shape[1]], [start, end], alpha=0)
    return line_fit


def COGIcoTip(img, start, end):
    xs = []
    ys = []
    zs = []
    for i_ in range(start, end):
        # Find slices where screw tip appears
        loc_ = np.array(np.where(img[:, :, i_] == 1))
        if loc_.shape[1]:
            # Find centre of gravity for slice i
            print('Tip on slice ' + str(i_))
            cog = [np.mean(loc_[0]), np.mean(loc_[1])]
            xs.append(cog[0])
            ys.append(cog[1])
            zs.append(i_)
    S = [np.mean(xs), np.mean(ys), np.mean(zs)]
    return S, xs, ys, zs


def rotate3Daround(img, angle, axis):
    dim = img.shape
    img_rotated = np.zeros_like(img)
    if axis == 2:
        for i_ in range(dim[2]):
            img_rotated[:, :, i_] = imutils.rotate(img[:, :, i_], angle)
    if axis == 0:
        for i_ in range(dim[0]):
            img_rotated[i_, :, :] = imutils.rotate(img[i_, :, :], angle)
    if axis == 1:
        for i_ in range(dim[1]):
            img_rotated[:, i_, :] = imutils.rotate(img[:, i_, :], angle)
    return img_rotated


def normal_vector(p1_, p2_, p3_):
    # Find normal vectors to the plane of three points
    v1_ = p2_ - p1_
    v2_ = p3_ - p1_
    normal = np.cross(v1_, v2_)
    return normal / np.linalg.norm(normal)


def angle_between(v1_, v2_):
    # Angle between two vectors
    # Input form: np.array([1, 0, 0])
    dot = np.dot(v1_, v2_)
    mag_v1 = np.linalg.norm(v1_)
    mag_v2 = np.linalg.norm(v2_)
    return np.rad2deg(np.arccos(dot / (mag_v1 * mag_v2)))


def set_axes_equal(axy):
    """Make axes of 3D plot have equal scale so that spheres appear as spheres,
    cubes as cubes, etc.  This is one possible solution to Matplotlib 's
    ax.set_aspect('equal') and ax.axis('equal') not working for 3D.

    Input
      ax: a matplotlib axis, e.g., as output from plt.gca().
    """

    x_limits = axy.get_xlim3d()
    y_limits = axy.get_ylim3d()
    z_limits = axy.get_zlim3d()

    x_range = abs(x_limits[1] - x_limits[0])
    x_middle = np.mean(x_limits)
    y_range = abs(y_limits[1] - y_limits[0])
    y_middle = np.mean(y_limits)
    z_range = abs(z_limits[1] - z_limits[0])
    z_middle = np.mean(z_limits)

    # The plot bounding box is a sphere in the sense of the infinity
    # norm, hence I call half the max range the plot radius.
    plot_radius = 0.5*max([x_range, y_range, z_range])

    axy.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    axy.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    axy.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])


def fit_plane(points):
    # convert the points to a matrix format
    X_ = points[:, 0]
    Y_ = points[:, 1]
    Z_ = points[:, 2]
    A_ = np.column_stack((X_, Y_, np.ones(len(X_))))
    plane_, _, _, _ = np.linalg.lstsq(A_, Z_, rcond=None)
    return plane_


def rotation_matrix_from_vectors(vec1, vec2):
    """ Find the rotation matrix that aligns vec1 to vec2
    :param vec1: A 3d "source" vector
    :param vec2: A 3d "destination" vector
    :return mat: A transform matrix (3x3) which when applied to vec1, aligns it with vec2.
    """
    a_, b_ = (vec1 / np.linalg.norm(vec1)).reshape(3), (vec2 / np.linalg.norm(vec2)).reshape(3)
    v = np.cross(a_, b_)
    c = np.dot(a_, b_)
    s = np.linalg.norm(v)
    kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    rota_matrix = np.eye(3) + kmat + kmat.dot(kmat) * ((1 - c) / (s ** 2))
    return rota_matrix


def rotation_matrix_from_angles(theta1, theta2, theta3, order='xyz'):
    """
    input
        theta1, theta2, theta3 = rotation angles in rotation order (degrees)
        order = rotation order of x,y,z　e.g. XZY rotation -- 'xzy'
    output
        3x3 rotation matrix (numpy array)
    """
    c1 = np.cos(theta1 * np.pi / 180)
    s1 = np.sin(theta1 * np.pi / 180)
    c2 = np.cos(theta2 * np.pi / 180)
    s2 = np.sin(theta2 * np.pi / 180)
    c3 = np.cos(theta3 * np.pi / 180)
    s3 = np.sin(theta3 * np.pi / 180)
    matrix = []

    if order == 'xzx':
        matrix = np.array([[c2, -c3*s2, s2*s3],
                          [c1*s2, c1*c2*c3-s1*s3, -c3*s1-c1*c2*s3],
                          [s1*s2, c1*s3+c2*c3*s1, c1*c3-c2*s1*s3]])
    elif order == 'xyx':
        matrix = np.array([[c2, s2*s3, c3*s2],
                          [s1*s2, c1*c3-c2*s1*s3, -c1*s3-c2*c3*s1],
                          [-c1*s2, c3*s1+c1*c2*s3, c1*c2*c3-s1*s3]])
    elif order == 'yxy':
        matrix = np.array([[c1*c3-c2*s1*s3, s1*s2, c1*s3+c2*c3*s1],
                          [s2*s3, c2, -c3*s2],
                          [-c3*s1-c1*c2*s3, c1*s2, c1*c2*c3-s1*s3]])
    elif order == 'yzy':
        matrix = np.array([[c1*c2*c3-s1*s3, -c1*s2, c3*s1+c1*c2*s3],
                          [c3*s2, c2, s2*s3],
                          [-c1*s3-c2*c3*s1, s1*s2, c1*c3-c2*s1*s3]])
    elif order == 'zyz':
        matrix = np.array([[c1*c2*c3-s1*s3, -c3*s1-c1*c2*s3, c1*s2],
                          [c1*s3+c2*c3*s1, c1*c3-c2*s1*s3, s1*s2],
                          [-c3*s2, s2*s3, c2]])
    elif order == 'zxz':
        matrix = np.array([[c1*c3-c2*s1*s3, -c1*s3-c2*c3*s1, s1*s2],
                          [c3*s1+c1*c2*s3, c1*c2*c3-s1*s3, -c1*s2],
                          [s2*s3, c3*s2, c2]])
    elif order == 'xyz':
        matrix = np.array([[c2*c3, -c2*s3, s2],
                          [c1*s3+c3*s1*s2, c1*c3-s1*s2*s3, -c2*s1],
                          [s1*s3-c1*c3*s2, c3*s1+c1*s2*s3, c1*c2]])
    elif order == 'xzy':
        matrix = np.array([[c2*c3, -s2, c2*s3],
                          [s1*s3+c1*c3*s2, c1*c2, c1*s2*s3-c3*s1],
                          [c3*s1*s2-c1*s3, c2*s1, c1*c3+s1*s2*s3]])
    elif order == 'yxz':
        matrix = np.array([[c1*c3+s1*s2*s3, c3*s1*s2-c1*s3, c2*s1],
                          [c2*s3, c2*c3, -s2],
                          [c1*s2*s3-c3*s1, c1*c3*s2+s1*s3, c1*c2]])
    elif order == 'yzx':
        matrix = np.array([[c1*c2, s1*s3-c1*c3*s2, c3*s1+c1*s2*s3],
                          [s2, c2*c3, -c2*s3],
                          [-c2*s1, c1*s3+c3*s1*s2, c1*c3-s1*s2*s3]])
    elif order == 'zyx':
        matrix = np.array([[c1*c2, c1*s2*s3-c3*s1, s1*s3+c1*c3*s2],
                          [c2*s1, c1*c3+s1*s2*s3, c3*s1*s2-c1*s3],
                          [-s2, c2*s3, c2*c3]])
    elif order == 'zxy':
        matrix = np.array([[c1*c3-s1*s2*s3, -c2*s1, c1*s3+c3*s1*s2],
                          [c3*s1+c1*s2*s3, c1*c2, s1*s3-c1*c3*s2],
                          [-c2*s3, s2, c2*c3]])
    return matrix


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

    theta1 = theta1 * 180 / np.pi
    theta2 = theta2 * 180 / np.pi
    theta3 = theta3 * 180 / np.pi

    return theta1, theta2, theta3


def zeros_and_ones(img, th_):
    """Creates a new image with zeros (below threshold) and ones only"""
    img01 = np.array((img >= th_).astype(int))
    return img01



t1 = time.time()
plt.close('all')

loc = '/home/biomech/Documents/01_Icotec/01_Experiments/02_Scans/Pilot3/04_Registered/'

bone = 'XCT_Icotec_S130672_L5_intact_planned.mhd'
inst = 'ICOTEC_S130672_L5_implants_XCTres.mhd'

im0 = load_itk(loc + bone)
imD = load_itk(loc + inst)

# Convert to COS:
lineT = axis3D(imD[0], 670, 1100, 0, 0)
'''
lineI = np.array([187, 466, 451]) - np.array([820, 474, 54])
MT = rotation_matrix_from_vectors(lineT.vector, [0, 1, 0])
MI = rotation_matrix_from_vectors(lineI, [0, 1, 0])
[angTX, angTXY, angTXYZ] = rotation_angles_from_matrix(MT, 'xyz')
[angIX, angIXY, angIXYZ] = rotation_angles_from_matrix(MI, 'xyz')
w = int(np.ceil(11.5/im0[2][0]/2))
h = int(np.ceil(17.5/im0[2][0]/2))

im0T_rotX = rotate3Daround(im0[0], angTX, 0)
im0I_rotX = rotate3Daround(im0[0], angIX, 0)
del im0
imDT_rotX = rotate3Daround(imD[0], angTX, 0)
imDI_rotX = rotate3Daround(imD[0], angIX, 0)
del imD
print('RotX done.')

im0T_rotXY = rotate3Daround(im0T_rotX, angTXY, 1)
im0I_rotXY = rotate3Daround(im0I_rotX, angIXY, 1)
del im0T_rotX
del im0I_rotX
imDT_rotXY = rotate3Daround(imDT_rotX, angTXY, 1)
imDI_rotXY = rotate3Daround(imDI_rotX, angIXY, 1)
del imDT_rotX
del imDI_rotX
print('RotY done.')

im0T_rotXYZ = rotate3Daround(im0T_rotXY, angTXYZ, 2)
im0I_rotXYZ = rotate3Daround(im0I_rotXY, angIXYZ, 2)
del im0T_rotXY
del im0I_rotXY
imDT_rotXYZ = rotate3Daround(imDT_rotXY, angTXYZ, 2)
imDI_rotXYZ = rotate3Daround(imDI_rotXY, angIXYZ, 2)
del imDT_rotXY
del imDI_rotXY
print('RotZ done.')

im0T_fin = rotate3Daround(im0T_rotXYZ, -38.5, 1)
im0I_fin = rotate3Daround(im0I_rotXYZ, 0, 1)
del im0T_rotXYZ
del im0I_rotXYZ
imDT_fin = rotate3Daround(imDT_rotXYZ, -38.5, 1)
imDI_fin = rotate3Daround(imDI_rotXYZ, 0, 1)
del imDT_rotXYZ
del imDI_rotXYZ
imDT_fin = zeros_and_ones(imDT_fin, 0.5)
imDI_fin = zeros_and_ones(imDI_fin, 0.5)
xray(im0I_fin*0+imDI_fin*5000, 0)
xray(im0I_fin*0+imDI_fin*5000, 1)
xray(im0I_fin*0+imDI_fin*5000, 2)

boneT = im0T_fin[584-h:584+h, :, 947-w:947+w]
# boneI = im0I_fin[584-h:584+h, :, 947-w:947+w]
#xray(boneT, 0)
#xray(boneT, 1)
#xray(boneT, 2)
'''

'''
# angX = 49
im3_rotX = np.zeros((len(im0c[:, 0, 0]),
                     len(rot(im0c[0, :, :], [0, 0], angX)[0]),
                     len(rot(im0c[0, :, :], [0, 0], angX)[0][0])))
for i in range(len(im0c[:, 0, 0])):
    [im3_rotX[i, :, :], _] = rot(im0c[i, :, :], [0, 0], angX)
    if np.mod(i, 50) == 0:
        print(str(i) + '/' + str(len(im0c[:, 0, 0])))
im3_rotX = np.array(im3_rotX)
print('First rotation completed.')
#xray(im0c, 0)
#xray(im3_rotX, 0)
del im0c

im3_rotXY = np.zeros((len(rot(im3_rotX[:, 0, :], [0, 0], angXY)[0]),
                      len(im3_rotX[0, :, 0]),
                      len(rot(im3_rotX[:, 0, :], [0, 0], angXY)[0][0])))
for i in range(len(im3_rotX[0, :, 0])):
    [im3_rotXY[:, i, :], _] = rot(im3_rotX[:, i, :], [0, 0], angXY)
    if np.mod(i, 50) == 0:
        print(str(i) + '/' + str(len(im3_rotX[0, :, 0])))
im3_rotXY = np.array(im3_rotXY)
print('Second rotation completed.')
#xray(im3_rotX, 1)
#xray(im3_rotXY, 1)
del im3_rotX
im3_rotXYZ = np.zeros((len(rot(im3_rotXY[:, :, 0], [0, 0], angXYZ)[0]),
                       len(rot(im3_rotXY[:, :, 0], [0, 0], angXYZ)[0][0]),
                       len(im3_rotXY[0, 0, :])))
for i in range(len(im3_rotXY[0, 0, :])):
    [im3_rotXYZ[:, :, i], _] = rot(im3_rotXY[:, :, i], [0, 0], angXYZ)
    if np.mod(i, 50) == 0:
        print(str(i) + '/' + str(len(im3_rotXY[0, 0, :])))
im3_rotXYZ = np.array(im3_rotXYZ)
print('Third rotation completed.')
xray(im3_rotXYZ, 0)
xray(im3_rotXYZ, 1)
xray(im3_rotXYZ, 2)
del im3_rotXY
'''

'''
fig = plt.figure()
ax = fig.add_subplot(projection='3d')
vox = np.where(imD[0])
plt.scatter(vox[0] + imD[1][0], vox[1] + imD[1][1], vox[2] + imD[1][2])
vox = np.where(imI[0])
plt.scatter(vox[0] + imI[1][0], vox[1] + imI[1][1], vox[2] + imI[1][2])
'''

'''
loc = '/home/biomech/Documents/01_Icotec/01_Experiments/02_Scans/Pilot2/04_Registered/'
# file1 = 'Icotec_S130684_L5_intact.mhd'
# file2 = 'Icotec_S130684_L5_predrilled.mhd'
file3 = 'Icotec_S130684_L5_instrumented.mhd'
im3 = load_itk(loc + file3)[0]
im3_seg = seg(load_itk(loc + file3)[0], 18000, np.inf)
lfV = axis3D(im3_seg, 140, 230, 0, 0)
tipP = COGIcoTip(im3_seg, 100, 140)

direct = [0, 1, 0]
mRot = np.outer(lfV.vector, direct)
im3_PMMA = seg(load_itk(loc + file3)[0], 6900, 7320)

# Points from PMMA surface read from ITK-SNAP
PMMA_points = np.array([[  4, 211,  57],
                        [  4, 253, 140],
                        [ 33, 163, 160],
                        [  1, 200,  25],
                        [ 11, 200,  88],
                        [ 19, 200, 143],
                        [ 27, 200, 194],
                        [ 33, 200, 241],
                        [ 31, 125,  70],
                        [ 31, 166, 158],
                        [ 31, 202, 225],
                        [ 35, 136, 121],
                        [ 40, 150, 178],
                        [ 48, 158, 233]])
PMMA_centre = [np.mean(PMMA_points[:3, 0]), np.mean(PMMA_points[:3, 1]), np.mean(PMMA_points[:3, 2])]

# ICOTEC mounting
# [118, 176, 222]
# [ 65, 261,  28]
p1_PEEK = np.array([118, 176, 222])
p2_PEEK = np.array([ 65, 261,  28])
v_PEEK = (p1_PEEK - p2_PEEK) / np.linalg.norm(p1_PEEK - p2_PEEK)
# ICOTEC screw
# [ 91, 206, 104] tip
# [111,  96,  51]
tip_PEEK = np.array([ 87, 204, 104])
head_PEEK = np.array([111,  96,  51])
screw_PEEK = (tip_PEEK - head_PEEK) / np.linalg.norm(tip_PEEK - head_PEEK)

# Titanium mounting
# [100, 140, 50] ok
# [86, 283, 216]  # guessed!
p1_Ti = np.array([ 100, 140,  50])
p2_Ti = np.array([ 86, 290, 211])
v_Titanium = (p1_Ti - p2_Ti) / np.linalg.norm(p1_Ti - p2_Ti)
# Titanium screw
# [ 83, 214, 147] tip
# [123, 141, 228]
tip_Ti = np.array([ 83, 214, 143])
head_Ti = np.array([130, 129, 242])
screw_Ti = (tip_Ti - head_Ti) / np.linalg.norm(tip_Ti - head_Ti)

for i in range(len(PMMA_points)):
    im3_PMMA[PMMA_points[i, 0], PMMA_points[i, 1], PMMA_points[i, 2]] = 200

xray(im3_PMMA, 0)
xray(im3_PMMA, 1)
xray(im3_PMMA, 2)
'''

'''
n_PMMA = normal_vector(PMMA_points[0], PMMA_points[1], PMMA_points[2])
plt.figure()
ax = plt.axes(projection='3d')
# Plot PMMA points (lower surface, read in itk)
for i in range(len(PMMA_points)):
    ax.scatter3D(PMMA_points[i, 0], PMMA_points[i, 1], PMMA_points[i, 2], c='k')
    for j in range(i+1, len(PMMA_points)):
        ax.plot3D([PMMA_points[i, 0], PMMA_points[np.mod(j, 3), 0]],
                  [PMMA_points[i, 1], PMMA_points[np.mod(j, 3), 1]],
                  [PMMA_points[i, 2], PMMA_points[np.mod(j, 3), 2]],
                  c='k')

# Plot PMMA point connections
ax.plot3D([PMMA_points[0, 0], PMMA_points[1, 0]],
          [PMMA_points[0, 1], PMMA_points[1, 1]],
          [PMMA_points[0, 2], PMMA_points[1, 2]],
          c='k')
# Plot PMMA normal vector
t = 150
ax.plot3D([PMMA_centre[0] - n_PMMA[0] * 0, PMMA_centre[0] + n_PMMA[0] * t],
          [PMMA_centre[1] - n_PMMA[1] * 0, PMMA_centre[1] + n_PMMA[1] * t],
          [PMMA_centre[2] - n_PMMA[2] * 0, PMMA_centre[2] + n_PMMA[2] * t],
          c='k')
# PMMA surface
plane = fit_plane(PMMA_points)
X, Y = np.meshgrid(np.linspace(np.min(PMMA_points[:, 0]), np.max(PMMA_points[:, 0]), 20),
                   np.linspace(np.min(PMMA_points[:, 1]), np.max(PMMA_points[:, 1]), 20))
Z = np.array([[plane[0] * x + plane[1] * y + plane[2] for x, y in zip(X_row, Y_row)] for X_row, Y_row in zip(X, Y)])
ax.scatter3D(X, Y, Z, c='g', alpha=0.1)

# Plot Ti screw vector (found with script)
ax.plot3D([lfV.point[0] - lfV.vector[0] * t, lfV.point[0] + lfV.vector[0] * t],
          [lfV.point[1] - lfV.vector[1] * t, lfV.point[1] + lfV.vector[1] * t],
          [lfV.point[2] - lfV.vector[2] * t, lfV.point[2] + lfV.vector[2] * t],
          c='b')
# Plot Ti mounting
ax.scatter3D(p1_Ti[0], p1_Ti[1], p1_Ti[2], c='b')
ax.scatter3D(p2_Ti[0], p2_Ti[1], p2_Ti[2], c='b')
ax.plot3D([p1_Ti[0], p2_Ti[0]],
          [p1_Ti[1], p2_Ti[1]],
          [p1_Ti[2], p2_Ti[2]],
          c='b')
# Plot Ti screw (read with itk)
ax.scatter3D(tip_Ti[0], tip_Ti[1], tip_Ti[2], c='b')
ax.scatter3D(head_Ti[0], head_Ti[1], head_Ti[2], c='b')
ax.plot3D([tip_Ti[0], head_Ti[0]],
          [tip_Ti[1], head_Ti[1]],
          [tip_Ti[2], head_Ti[2]],
          c='b', linestyle='dashed')

# Plot PEEK mounting
ax.scatter3D(p1_PEEK[0], p1_PEEK[1], p1_PEEK[2], c='r')
ax.scatter3D(p2_PEEK[0], p2_PEEK[1], p2_PEEK[2], c='r')
ax.plot3D([p1_PEEK[0], p2_PEEK[0]],
          [p1_PEEK[1], p2_PEEK[1]],
          [p1_PEEK[2], p2_PEEK[2]],
          c='r')
# Plot PEEK screw (read with itk)
ax.scatter3D(tip_PEEK[0], tip_PEEK[1], tip_PEEK[2], c='r')
ax.scatter3D(head_PEEK[0], head_PEEK[1], head_PEEK[2], c='r')
ax.plot3D([tip_PEEK[0], head_PEEK[0]],
          [tip_PEEK[1], head_PEEK[1]],
          [tip_PEEK[2], head_PEEK[2]],
          c='r', linestyle='dashed')

# Angle information about Titanium screw
print('\nTitanium screw:')
print('Angle between PMMA-surface and screw: ' + str(angle_between(n_PMMA, lfV.vector)))
print('Angle between fixation and screw: ' + str(angle_between(v_Titanium, lfV.vector)))
print('Angle between PMMA-surface and fixation: ' + str(angle_between(n_PMMA, v_Titanium)))

# Angle information about PEEK screw
print('\nPEEK screw:')
print('Angle between PMMA-surface and screw: ' + str(angle_between(n_PMMA, screw_PEEK)))
print('Angle between fixation and screw: ' + str(angle_between(v_PEEK, screw_PEEK)))
print('Angle between PMMA-surface and fixation: ' + str(angle_between(n_PMMA, v_PEEK)))
'''


print('\nRuntime: ' + str(round(time.time() - t1, 2)) + ' seconds.')


p3 = np.array([1130, 453, 764])  # point on Ti screw, z-axis, origin of COS

p1 = np.array([54, 474, 820])  # point on rotation axis, x-axis
v3 = -lineT.vector  # z-axis, found by function
v2 = np.cross(p1-p3, v3)/np.linalg.norm(np.cross(p1-p3, v3))  # y-axis
v1 = np.cross(v2, v3)  # x-axis
rotation_matrix = np.vstack((np.append(v1, 0), np.append(v2, 0), np.append(v3, 0), np.array([0, 0, 0, 1])))
translation_matrix = np.array([[1, 0, 0, p3[0]],
                               [0, 1, 0, p3[1]],
                               [0, 0, 1, p3[2]],
                               [0, 0, 0,    1]])
COS_CT = np.dot(rotation_matrix, translation_matrix)
COS_FE = np.eye(4)
