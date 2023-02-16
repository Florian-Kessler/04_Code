import SimpleITK as sitk
from skimage.measure import label
import numpy as np
import matplotlib.pyplot as plt
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


def xray(image, axis_):
    plt.figure()
    plt.imshow(np.sum(image, axis=axis_), cmap='gray')
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


def axis3D(img, start, end, axis_):
    """
    Function to find screw axis in 3d image.
    :param img: 3d array of the segmented image. Background = 0, screw >= 1.
    :param start: Value <= ROI
    :param end: Value >= ROI
    :param axis_: Along which axis to evaluate ('x', 'y' or 'z')
    :return: Line containing point and vector
    """
    # centre of gravity for each slice of the screw
    xs = []
    ys = []
    zs = []
    if axis_ == 'x':
        for i_ in range(start, end):
            # Find slices where screw appears
            loc_ = np.array(np.where(img[i_, :, :] >= 1))
            if loc_.shape[1]:
                # Find centre of gravity for slice i
                cog = [np.mean(loc_[0]), np.mean(loc_[1])]
                xs.append(i_)
                ys.append(cog[0])
                zs.append(cog[1])
    elif axis_ == 'y':
        for i_ in range(start, end):
            # Find slices where screw appears
            loc_ = np.array(np.where(img[:, i_, :] >= 1))
            if loc_.shape[1]:
                # Find centre of gravity for slice i
                cog = [np.mean(loc_[0]), np.mean(loc_[1])]
                xs.append(cog[0])
                ys.append(i_)
                zs.append(cog[1])
    elif axis_ == 'z':
        for i_ in range(start, end):
            # Find slices where screw appears
            loc_ = np.array(np.where(img[:, :, i_] >= 1))
            if loc_.shape[1]:
                # Find centre of gravity for slice i
                cog = [np.mean(loc_[0]), np.mean(loc_[1])]
                xs.append(cog[0])
                ys.append(cog[1])
                zs.append(i_)
    line_fit = Line.best_fit(np.array([xs, ys, zs]).T)
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

    # theta1 = theta1 * 180 / np.pi
    # theta2 = theta2 * 180 / np.pi
    # theta3 = theta3 * 180 / np.pi

    return theta1, theta2, theta3


def zeros_and_ones(img, th_):
    """Creates a new image with zeros (below threshold) and ones only"""
    img01 = np.array((img >= th_).astype(int))
    return img01
