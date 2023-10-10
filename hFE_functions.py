import SimpleITK as sitk
from skimage.measure import label
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import rotate
from colorama import Fore, Style
from scipy import stats
from skspatial.objects import Line
import imutils
import sys
import os
import pickle
import statsmodels.api as sm
import pandas as pd


def load_itk(filename):
    """
    Function to load a raw/mhd image. Includes transformation from zyx to xyz.
    :param filename: path to mhd file, raw path should be specified in mhd header
    :rtype filename: str
    :return: tuple(ct_scan, origin, spacing) ---
             ct_scan is the image array ---
             origin is the origin of the image read from the mhd-file ---
             spacing is the voxel size
    """
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


def plot_grey(image, axis=0, pos=0):
    """
    Grey plot of a 3d image slice
    :param image: 3d array
    :param axis: Optional; axis of cut (0, 1, 2)
    :rtype axis: int
    :param pos: Optional; cutting position
    :rtype pos: int
    :return:
    """
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


def plot_HU(image, axis=0, pos=0, level=0, window=1000):
    """
    Plot image in defined HU window
    :param image: 3d array of image in HU values
    :param axis: Optional; axis of cut (0, 1, 2)
    :rtype axis: int
    :param pos: Optional; position within image on axis
    :rtype pos: int
    :param level: Optional; level of centre of HU-window
    :rtype level: float
    :param window: Optional; size of window
    :rtype window: float
    :return: slice with HU-parameters
    """
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
    return slic


def rot(image, centre_xy, angle):
    """
    Rotates an image around a rotation centre.
    :param image: 2d array of image
    :param centre_xy: rotation centre
    :rtype centre_xy: tuple[int]
    :param angle: rotation angle in degrees
    :return: tuple(Rotated image 2d array, new rotation centre)
    """
    im_rot = rotate(image, angle)
    org_center = (np.array(image.shape[:2][::-1]) - 1) / 2.
    rot_center = (np.array(im_rot.shape[:2][::-1]) - 1) / 2.
    org = centre_xy - org_center
    a_ = np.deg2rad(angle)
    new = np.array([org[0] * np.cos(a_) + org[1] * np.sin(a_),
                    - org[0] * np.sin(a_) + org[1] * np.cos(a_)])
    return im_rot, new + rot_center


def xray(image, axis_=0):
    """
    Plot an x-ray of ad 3d image (=2d projection)
    :param image: 3d array of image
    :param axis_: Optional; axis (0, 1, 2), default = 0
    :rtype axis_: int
    :return:
    """
    plt.figure()
    plt.imshow(np.sum(image, axis=axis_), cmap='gray')
    plt.show()


def getLargestCC(segm_im):
    """
    Function to find the largest connected region
    :param segm_im: segmented image array consisting of 0 and 1
    :return: segmented image array with only largest connected region
    """
    labels = label(segm_im)
    assert (labels.max() != 0)  # assume at least 1 CC
    largestCC = labels == np.argmax(np.bincount(labels.flat)[1:]) + 1
    return largestCC


def countOccurrence(arr, x_):
    """
    Counts the occurrence of x in array arr
    :param arr: any nd-array
    :param x_: value to count occurrence in array
    :rtype x_: int
    :return: number of occurrences
    """
    res = 0
    arrF = arr.flatten()
    for j in range(len(arrF)):
        if x_ == arrF[j]:
            res += 1
    return res


def islandFinder(img):
    """
    Finds connected islands (connectivity = 1, 2, ..., size(img), ).
    Returns an array containing the boolean array of each island and
    the size of each island.
    :param img: 3d image array
    :returns: tuple (islands, island_size) ---
              island is array with islands ---
              island_size is array with size of islands, correspond to islands-array
    """
    [imglabeled, island_count] = label(img, background=0, return_num=True, connectivity=1)
    islands = np.zeros(np.append(island_count, np.append(1, imglabeled.shape)))
    island_size = np.zeros(island_count)
    for j in range(island_count):
        islands[j] = np.array(imglabeled == j + 1).astype(int)
        island_size[j] = countOccurrence(imglabeled, j + 1)
    return islands, island_size


def axis2D3D(img, start, end, test=0, plot=0):
    """
    Find screw axis with linear regressions on two 2d-projections along 3rd (z) axis
    :param img: 3d image array
    :param start: start of ROI
    :type
    :rtype start: int
    :param end: end of ROI
    :rtype end: int
    :param test: Optional; if 1 shows figures for visual check on each layer
    :rtype test: int
    :param plot: Optional; if 1 shows 3d plot at the end
    :rtype plot: int
    :return: two points on the regression axis in 3d
    """
    # centre of gravity for each slice of the screw
    xs = []
    ys = []
    zs = []
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


def axis3D(img, start, end, axis_):
    """
    Function to find screw axis in 3d image. Orthogonal fit
    :param img: 3d array of the segmented image. Background = 0, screw >= 1.
    :param start: Value <= ROI
    :param end: Value >= ROI
    :param axis_: Along which axis to evaluate ('x', 'y' or 'z')
    :return: point and vector
    """
    # centre of gravity for each slice of the screw
    xs = []
    ys = []
    zs = []
    if axis_ == 'x':
        for i_ in range(start, end):
            # Find slices where screw appears
            loc_ = np.array(np.where(img[i_, :, :] >= 2))
            if loc_.shape[1]:
                # Find centre of gravity for slice i
                cog = [np.mean(loc_[0]), np.mean(loc_[1])]
                xs.append(i_)
                ys.append(cog[0])
                zs.append(cog[1])
    elif axis_ == 'y':
        for i_ in range(start, end):
            # Find slices where screw appears
            loc_ = np.array(np.where(img[:, i_, :] >= 2))
            if loc_.shape[1]:
                # Find centre of gravity for slice i
                cog = [np.mean(loc_[0]), np.mean(loc_[1])]
                xs.append(cog[0])
                ys.append(i_)
                zs.append(cog[1])
    elif axis_ == 'z':
        for i_ in range(start, end):
            # Find slices where screw appears
            loc_ = np.array(np.where(img[:, :, i_] >= 2))
            if loc_.shape[1]:
                # Find centre of gravity for slice i
                cog = [np.mean(loc_[0]), np.mean(loc_[1])]
                xs.append(cog[0])
                ys.append(cog[1])
                zs.append(i_)
    line_fit = Line.best_fit(np.array([xs, ys, zs]).T)
    return line_fit


def COGIcoTip(img, start, end, axis):
    """
    Function to find the center of gravity (COG) in a segmented image
    :param img: segmented 3d image
    :param start: starting point for evaluation on specified axis
    :type start: int
    :param end: end point for evaluation on specified axis
    :type end: int
    :param axis:  along which the evaluation should follow (0/1/2)
    :type axis: int
    :return: center of gravity
    """
    xs = []
    ys = []
    zs = []
    S = 0
    if axis == 0:
        for i_ in range(start, end):
            # Find slices where screw tip appears
            loc_ = np.array(np.where(img[i_, :, :] >= 1))
            if loc_.shape[1]:
                # Find centre of gravity for slice i_
                # print('Tip on slice ' + str(i_))
                cog = [np.mean(loc_[0]), np.mean(loc_[1])]
                xs.append(i_)
                ys.append(cog[0])
                zs.append(cog[1])
        S = [np.mean(xs), np.mean(ys), np.mean(zs)]
    elif axis == 1:
        for i_ in range(start, end):
            # Find slices where screw tip appears
            loc_ = np.array(np.where(img[:, i_, :] >= 1))
            if loc_.shape[1]:
                # Find centre of gravity for slice i_
                # print('Tip on slice ' + str(i_))
                cog = [np.mean(loc_[0]), np.mean(loc_[1])]
                xs.append(cog[0])
                ys.append(i_)
                zs.append(cog[1])
        S = [np.mean(xs), np.mean(ys), np.mean(zs)]
    elif axis == 2:

        for i_ in range(start, end):
            # Find slices where screw tip appears
            loc_ = np.array(np.where(img[:, :, i_] >= 1))
            if loc_.shape[1]:
                # Find centre of gravity for slice i_
                # print('Tip on slice ' + str(i_))
                cog = [np.mean(loc_[0]), np.mean(loc_[1])]
                xs.append(cog[0])
                ys.append(cog[1])
                zs.append(i_)
        S = [np.mean(xs), np.mean(ys), np.mean(zs)]
    return S  # , xs, ys, zs


def rotate3Daround(img, angle, axis):
    """
    Function to rotate a 3d image around a coordinate axis by a specified angle
    :param img: 3d image array
    :param angle: angle in degrees (positive direction = anti-clockwise)
    :param axis: 0/1/2 (standing for x, y, z axis)
    :return: rotated image in same dimension as input (hence image data gets lost)
    """
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
    """
    Function to find the normal vector to a plane defined by three points
    :param p1_: Point1 as 3d np array
    :param p2_: Point2 as 3d np array
    :param p3_: Point3 as 3d np array
    :return: Normalised vector normal to the plane. (Note: direction depending on order of points)
    """
    v1_ = p2_ - p1_
    v2_ = p3_ - p1_
    normal = np.cross(v1_, v2_)
    return normal / np.linalg.norm(normal)


def angle_between(v1_, v2_):
    """
    Function to find the angle between two vectors
    :param v1_: First 3d-vector as np array
    :param v2_: Second 3d-vector as np array
    :return: angle in degrees
    """
    dot = np.dot(v1_, v2_)
    mag_v1 = np.linalg.norm(v1_)
    mag_v2 = np.linalg.norm(v2_)
    return np.rad2deg(np.arccos(dot / (mag_v1 * mag_v2)))


def set_axes_equal(axy):
    """
    Make axes of 3D plot have equal scale so that spheres appear as spheres,
    cubes as cubes, etc.  This is one possible solution to Matplotlib 's
    ax.set_aspect('equal') and ax.axis('equal') not working for 3D.
    :param axy: a matplotlib axis, e.g., as output from plt.gca().
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
    """
    Find the rotation matrix that aligns vec1 to vec2
    :param vec1: A 3d "source" vector
    :param vec2: A 3d "destination" vector
    :return: A transform matrix (3x3) which when applied to vec1, aligns it with vec2.
    """
    a_, b_ = (vec1 / np.linalg.norm(vec1)).reshape(3), (vec2 / np.linalg.norm(vec2)).reshape(3)
    v = np.cross(a_, b_)
    c = np.dot(a_, b_)
    s = np.linalg.norm(v)
    kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    rota_matrix = np.eye(3) + kmat + kmat.dot(kmat) * ((1 - c) / (s ** 2))
    return rota_matrix


def rotation_matrix_from_angles(theta1, theta2, theta3, order):
    """
    Function to create a rotation matrix from given angles in a specified order
    :param theta1: angle in rad, rotation around axis 1 (x-axis)
    :param theta2: angle in rad, rotation around axis 2 (y-axis)
    :param theta3: angle in rad, rotation around axis 3 (z-axis)
    :param order: rotation order of x,y,z　e.g. XZY rotation --> 'xzy'
    :return: 3x3 rotation matrix (numpy array)
    """
    c1 = np.cos(theta1)  # * np.pi / 180)
    s1 = np.sin(theta1)  # * np.pi / 180)
    c2 = np.cos(theta2)  # * np.pi / 180)
    s2 = np.sin(theta2)  # * np.pi / 180)
    c3 = np.cos(theta3)  # * np.pi / 180)
    s3 = np.sin(theta3)  # * np.pi / 180)
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
    This function finds rotation angles from a rotation matrix in a specified order of rotation.
    :param matrix: 3x3 rotation matrix (numpy array)
    :param order: (str) rotation order of x, y, z : e.g, rotation XZY -- 'xzy'
    :return: theta1, theta2, theta3 = rotation angles in rotation order
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

    # theta1 = theta1 * 180 / np.pi
    # theta2 = theta2 * 180 / np.pi
    # theta3 = theta3 * 180 / np.pi

    return theta1, theta2, theta3


def zeros_and_ones(img, th_):
    """
    This function creates a new image with zeros (below threshold) and ones only
    :param img: grey image array
    :param th_: threshold; below = 0, above = 1
    :return: np array from segmented image"""
    return np.array((img >= th_).astype(int))


def blockPrint():
    """
    Function to block print output
    :return:
    """
    sys.stdout = open(os.devnull, 'w')  # disable printing


def enablePrint():
    """
    Function to enable print output
    :return:
    """
    sys.stdout = sys.__stdout__   # enable printing


def read_RFnodeFile(file_):
    # read data from text file
    df_ = np.loadtxt(file_, delimiter=',')

    # Contains the frame number:
    # frame = np.array(df_[:, 0])

    # reaction moment of the implant:
    # rmx_ = np.array(df_[:, 1])
    # rmy_ = np.array(df_[:, 2])
    # rmz_ = np.array(df_[:, 3])

    # reaction force of the implant:
    # rfx_ = np.array(df_[:, 4])
    rfy_ = np.array(df_[:, 5])
    # rfz_ = np.array(df_[:, 6])

    # transverse disp. of the implant:
    # ux_ = np.array(df_[:, 7])
    uy_ = np.array(df_[:, 8])
    # uz_ = np.array(df_[:, 9])

    # rotational disp. of the implant:
    # urx_ = np.array(df_[:, 10])
    # ury_ = np.array(df_[:, 11])
    # urz_ = np.array(df_[:, 12])

    return uy_, rfy_


def read_RFnodeFile_uFE(file_):
    # read data from text file
    df_ = np.loadtxt(file_, delimiter='\t')

    t_ = np.array(df_[:, 0])
    d_ = np.array(df_[:, 1])

    return t_, d_


def read_energy(file_):
    # read data from text file
    df_ = np.loadtxt(file_)  # , delimiter='   ')
    t_ = np.array(df_[:, 0])
    e_ = np.array(df_[:, 1])
    return t_, e_


def read_acumen(file_a):
    df_ = pd.read_csv(file_a, delimiter=';', skiprows=[0, 2])
    t_ = pd.DataFrame(df_, columns=['Time ']).to_numpy()
    d_ = pd.DataFrame(df_, columns=['Axial Displacement ']).to_numpy()
    # d_ = d_ - d_[0]  # calibrate displacement to zero at beginning
    f_ = pd.DataFrame(df_, columns=['Axial Force ']).to_numpy()
    # f_set_ = pd.DataFrame(df_, columns=['Axial Force Command ']).to_numpy()
    cycle_ = pd.DataFrame(df_, columns=['Axial Count ']).to_numpy()
    # arr_ = 0
    peak_ = np.zeros(int(np.max(cycle_)))
    vall_ = np.zeros(int(np.max(cycle_)))
    for j_ in range(2, int(np.max(cycle_))):
        # del arr_
        arr_ = np.where((cycle_ == j_) | (cycle_ == j_ + .5))[0]
        peak_[j_] = arr_[int(np.argmin(f_[arr_]))]
        vall_[j_] = arr_[int(np.argmax(f_[arr_]))]

    peak_ = peak_.astype(int)
    vall_ = vall_.astype(int)

    return cycle_, d_, f_, peak_, vall_, t_


def peaks_raw(i_, file_a='/home/biomech/Documents/01_Icotec/01_Experiments/03_Analysis/peaks_raw.csv'):
    df_ = pd.read_csv(file_a, delimiter=',', skiprows=[1])
    Acumen_ = df_.iloc[i_]
    return Acumen_


def read_ARAMIS(file_A):
    df_ = pd.read_csv(file_A, delimiter=';', skiprows=[0])
    stop = np.where(df_.isnull())[0][0]
    lx = pd.DataFrame(df_, columns=['RodCsys→BoneCsys.LX [mm]'])
    ly = pd.DataFrame(df_, columns=['RodCsys→BoneCsys.LY [mm]'])
    lz = pd.DataFrame(df_, columns=['RodCsys→BoneCsys.LZ [mm]'])
    phiX = pd.DataFrame(df_, columns=['RodCsys→BoneCsys.Phi(X) [°]'])
    thetaY = pd.DataFrame(df_, columns=['RodCsys→BoneCsys.Theta(Y) [°]'])
    psiZ = pd.DataFrame(df_, columns=['RodCsys→BoneCsys.Psi(Z) [°]'])
    t_ = pd.DataFrame(df_, columns=['Time UTC'])

    return lx[:stop], ly[:stop], lz[:stop], phiX[:stop], thetaY[:stop], psiZ[:stop], t_[:stop]


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx


def find_first(array, value):
    array = np.asarray(array)
    idx = next(xd for xd, val in enumerate(array)
               if val <= value)
    return idx


def read_resample(file_r):
    df_ = pd.read_csv(file_r, delimiter=',')
    A_x_ = pd.DataFrame(df_, columns=['Aramis X']).to_numpy()
    A_y_ = pd.DataFrame(df_, columns=['Aramis Y']).to_numpy()
    A_z_ = pd.DataFrame(df_, columns=['Aramis Z']).to_numpy()
    A_rx_ = pd.DataFrame(df_, columns=['Aramis rX']).to_numpy()
    A_ry_ = pd.DataFrame(df_, columns=['Aramis rY']).to_numpy()
    A_rz_ = pd.DataFrame(df_, columns=['Aramis rZ']).to_numpy()
    a_y_ = pd.DataFrame(df_, columns=['Acumen Y']).to_numpy()
    a_f_ = pd.DataFrame(df_, columns=['Acumen Fy']).to_numpy()
    a_c_ = pd.DataFrame(df_, columns=['Acumen C']).to_numpy()

    return A_x_, A_y_, A_z_, A_rx_, A_ry_, A_rz_, a_y_, a_f_, a_c_


def read_exp_peaks():
    with open('/home/biomech/Documents/01_Icotec/01_Experiments/03_Analysis/mergedDf.pkl', 'rb') as f_:
        data = pickle.load(f_)
    return data


def Peak_exp(ampl_, number_):
    d_ = read_exp_peaks()
    # level = 2 ** (ampl - 2)
    peakF_ = d_['MaxForce'][(number_ - 2) * 7 + ampl_]
    # print(d['DisplacementLevel'][(number - 2) * 7 + ampl])
    # print(d['Specimen'][(number - 2) * 7 + ampl])
    # print(peakF)
    return peakF_


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


def smooth(y_, box_pts):
    box = np.ones(box_pts) / box_pts
    y_smooth = np.convolve(y_, box, mode='same')
    return y_smooth

