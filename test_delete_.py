from MedtoolFunctions import medtool_functions as mf
from Imaging_Functions import CTdata as ct
import numpy as np
from scipy.spatial import ConvexHull, transform
from scipy import ndimage
import time
import os
import pandas as pd
import pyvista as pv
import SimpleITK as sitk
import skimage.morphology as morph
import utils_SA as utils
import io_utils_SA as io_utils
import sys
import scipy
import matplotlib.pyplot as plt
# import stltovoxel


# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Attention !!!!!!!!!!!!!!!!!!!!!
# discuss calibration curve for BVTV with Ph, currently the calibration curve of Schenk et al. 2022 is included

def HFE_mapping_trans(bone, filename1, filename2):
    """
    Material Mapping, including PSL ghost padding layers as copy of most distal and proximal layers
    """

    print('... start material mapping with copying boundary layers as ghost layers')

    BVTVscaled = bone["BVTVscaled"]
    MASK_array_T = bone["MASK_array_T"]  # COS vom CT gray-scale Bild
    FEelSize = bone["FEelSize"]
    Spacing = bone["Spacing"]
    elems = bone["elems"]
    nodes = bone["nodes"]
    elsets = bone["elsets"]
    evalue = np.array([1, 1, 1])  # bone["evalue"]
    evect = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])  # bone["evect"]
    maskX = bone['MaskX']  # max and min coordinates of the FE mesh
    maskY = bone['MaskY']
    maskZ = bone['MaskZ']
    offset2COS = np.array([np.min(maskX), np.min(maskY), np.min(maskZ)])

    elsets["BONE"] = []
    RHOb = {}
    RHOb_FE = {}
    RHOb_corrected = {}
    PHIb = {}
    mm = {}
    m = {}
    cogs = {}

    ROI_BVTV_size = 1.25  # config["roi_bvtv_size"]  # in mm

    print("FEelSize material mapping = " + str(FEelSize))

    for i, elem in enumerate(elems):

        # Find COG of element
        cog_real = np.mean(
            [np.asarray(nodes[node].get_coord() - offset2COS) for node in elems[elem].get_nodes()],
            # [np.asarray(nodes[node].get_coord()) for node in elems[elem].get_nodes()],
            axis=0,
        )

        elems[elem].set_cog(cog_real)

        cog = bone['Transform'].TransformPoint(cog_real)

        PHIbone = 1
        #
        if PHIbone > 0.0:
            # if elset_marker == "bone":
            elsets["BONE"].append(elem)
            elems[elem].set_elset('BONE')

            BVTVbone = utils.computeBVTV_onephase(
                cog, Spacing, ROI_BVTV_size, BVTVscaled, MASK_array_T, PHIbone
            )
            BVTVbone_FE = utils.computeBVTV_FEel(
                cog, Spacing, FEelSize, BVTVscaled, MASK_array_T
            )

            RHOb[elem] = BVTVbone
            RHOb_FE[elem] = BVTVbone_FE
            PHIb[elem] = PHIbone
            RHOb_corrected[elem] = BVTVbone * PHIbone

            m[elem], mm[elem] = evalue, evect  # fabric anisotropy, evaleu =Eigenvalue, evect= Eigen vector
            cogs[elem] = cog  # included just for debugging reasons

            # Update instance variables of element
            elems[elem].set_bvtv(BVTVbone)
            elems[elem].set_bvtv_fe(BVTVbone_FE)
            elems[elem].set_m(m)
            elems[elem].set_mm(mm)

        sys.stdout.write(
            "\r" + " ... material mapping element " + str(i) + "/" + str(len(elems))
        )
        sys.stdout.flush()

    # conversion to np array for calculation of BVTV and selection of only elements contained in ELSET BONE
    # -----------------------------------------------------------------------------------------------------
    # these arrays are only used for computations in the summary file, and not for the abaqus input file.
    # Therefore, only elements belonging to ELSET BONE are considered.

    PHIb_array = np.array([PHIb[k] for k in elsets["BONE"] if k in PHIb])
    RHOb_orig_array = np.array([RHOb[k] for k in elsets["BONE"] if k in RHOb])
    RHOb_FE_array = np.array([RHOb_FE[k] for k in elsets["BONE"] if k in RHOb_FE])
    cogs_array = np.array([cogs[k] for k in elsets["BONE"] if k in cogs])

    elems = {elem: elems[elem] for elem in elsets["BONE"]}
    elems_bone = elems

    INPname1 = filename1
    INPname2 = filename2
    # mf.writeAbaqus(INPname, None, nodes, None, elems, elsets, NscaResults=None)
    # *****************************************************************
    marray = np.real(np.mean([np.asarray(m[elem]) for elem in m.keys()], axis=0))
    mmarray1 = np.real(np.mean([np.asarray(mm[elem][0]) for elem in m.keys()], axis=0))
    mmarray2 = np.real(np.mean([np.asarray(mm[elem][1]) for elem in m.keys()], axis=0))
    mmarray3 = np.real(np.mean([np.asarray(mm[elem][2]) for elem in m.keys()], axis=0))

    # io_utils.print_mem_usage()

    print("... material mapping finished")

    # set bone values
    # bone["RHOb_array"] = RHOb_array
    bone["RHOb_orig_array"] = RHOb_orig_array
    bone["RHOb_FE_array"] = RHOb_FE_array
    bone["PHIb_array"] = PHIb_array
    bone["elems"] = elems
    bone["elems_BONE"] = elems_bone
    bone["nodes"] = nodes
    # bone["nodes_BONE"] = nodes_bone
    bone["elsets"] = elsets
    bone["marray"] = marray
    bone["mmarray1"] = mmarray1
    bone["mmarray2"] = mmarray2
    bone["mmarray3"] = mmarray3
    bone["cogs"] = cogs_array
    bone["CoarseFactor"] = bone["FEelSize"][0] / bone["Spacing"][0]

    # Write elements and material properties to input file
    print("\n ... update ABAQUS file       :  " + INPname1 + " and " + INPname2)
    outfile1 = open(INPname1, 'w')
    outfile2 = open(INPname2, 'w')
    outfile1.write("***********************************************************\n")
    outfile2.write("***********************************************************\n")
    outfile2.write("** MATERIALS\n")

    # Write node sets as elements with material properties
    for elem in elems:
        outfile1.write("*Elset, ELSET=Elset" + str(elem) + "\n")
        outfile1.write(str(elem) + ", \n")
        outfile1.write("**POSITION: X = " + str(cogs[elem][0]) + " Y = " + str(cogs[elem][1]) + " Z = " + str(
            cogs[elem][2]) + "\n")
        outfile1.write("*SOLID SECTION, ELSET=Elset" + str(elem) + ",  MATERIAL=Mat" + str(elem) + "\n")
        outfile1.write("***********************************************************\n")

        # Definition of the material parameters for the UMAT
        outfile2.write("*MATERIAL, NAME=Mat" + str(elem) + "\n")
        outfile2.write("*DEPVAR\n")
        outfile2.write("22,\n")
        outfile2.write("*USER MATERIAL, CONSTANTS=5\n")
        outfile2.write("**Probs1 (bone probs), BVTV of element, m1, m2, m2\n")
        outfile2.write("0.0 " + ", " + str(np.round(RHOb[elem], 5)) + ", " + "1., 1., 1. \n")
        outfile2.write("***********************************************************\n")

    outfile1.close()
    outfile2.close()

    return bone


def readInpBoneDummy(bone, inpFileBoneDummy):
    """
    This function will generate an Abaqus input file which defines for each element in the mesh file a separate
    element-set.
    :param meshInp: Abaqus input file which contains only the mesh
    :return: An Abaqus input file with the mesh information and element sets
    """

    # Read information from Dummy input file of bone
    inpDummy = mf.readAbaqus(inpFileBoneDummy)
    # title = inp[0]
    nodes = inpDummy[1]
    # nsets = inpDummy[2]
    elems = inpDummy[3]
    elsets = inpDummy[4]


    # callculate the volume of each element
    evol = []
    for i, eleid in enumerate(elems):
        nodeid = elems[eleid].get_nodes()
        node_coo = np.array([nodes[id].get_coord() for id in nodeid])
        evol.append(ConvexHull(node_coo).volume)

    # set bone values
    bone["elems"] = elems
    bone["nodes"] = nodes
    bone["elsets"] = elsets
    bone["FEelSize"] = evol

    # print('createEleSets')
    return bone


def list_txt_files(path):
    # list all text files in path:
    files = []
    for i in os.listdir(path):
        if i.endswith('.txt'):
            files.append(i)
    return files


def boneMeshMask(bone, path, filename, resolution, mask_name, controlplot=False, reshape=True, closing=True):
    """
    This function creates a mask form any stl file and returns a 3d array mask - and store the mask as mhd in the given
    path.
    :param path: path to store a mhd file of the mask
    :param filename: name of the stl file
    :param resolution: resolution (voxel size) of the mask
    :param mask_name: name of the mhd file
    :param controlplot: If true a control 3d image of the stl will pop up - close it to proceed
    :param reshape: sometimes the order of the slices for the 3d array does not match - activate it if mhd looks weird
    :param closing: sometimes there are some small holes in the mask - activate it if needed
    :return: 3d array mask
    """

    # read in the stl to generate the mask
    reader = pv.get_reader(filename)
    mesh = reader.read()

    if controlplot==True:
        mesh.plot()
        voxels = pv.voxelize(mesh, density=resolution)
        p = pv.Plotter()
        p.add_mesh(voxels, color=True, show_edges=True)
        p.add_mesh(mesh, color="red", opacity=0.5)
        p.show()

    x_min, x_max, y_min, y_max, z_min, z_max = mesh.bounds
    x = np.arange(x_min, x_max, resolution)
    y = np.arange(y_min, y_max, resolution)
    z = np.arange(z_min, z_max, resolution)
    x, y, z = np.meshgrid(x, y, z)

    # Create unstructured grid from the structured grid
    grid = pv.StructuredGrid(x, y, z)
    ugrid = pv.UnstructuredGrid(grid)

    # get part of the mesh within the mesh's bounding surface.
    selection = ugrid.select_enclosed_points(mesh.extract_surface(), tolerance=0.0, check_surface=False)
    mask_ = selection.point_data['SelectedPoints'].view(np.bool_)

    if reshape == True:
        # sometimes the order of the matrix gets changed
        mask = mask_.reshape([z.shape[2], z.shape[1], z.shape[0]])
        mask = mask[:, ::-1, :]
        mask = np.rot90(mask, k=-1, axes=(1, 2))
    else:
        mask = mask_.reshape(x.shape)

    mask = np.array(mask)

    if closing == True:
        # close small holes and gabs:
        for i in range(0, mask.shape[0]):
            mask[i, :, :] = morph.closing(mask[i, :, :], np.ones([3, 3]))
            # mask_copy[i, :, :] = morph.dilation(mask_copy[i, :, :], np.ones([3, 3]))
            # mask_copy[i, :, :] = morph.erosion(mask_copy[i, :, :], np.ones([2, 2]))


    origin = [0, 0, 0]
    spacing = np.array([1, 1, 1]) * resolution

    mask_trans = mask.astype(np.short)
    itkmask = sitk.GetImageFromArray(mask_trans, isVector=None)
    itkmask.SetSpacing(spacing)
    itkmask.SetOrigin(origin)

    path_to_local_folder = path
    sitk.WriteImage(itkmask, f'{path_to_local_folder}/{mask_name}')

    # set bone values
    bone['MASK_array'] = mask_trans.T
    # Wird benötigt, um das KOS an die richtige position in der Maske zu verschieben
    bone['MaskX'] = np.array([x_min, x_max])
    bone['MaskY'] = np.array([y_min, y_max])
    bone['MaskZ'] = np.array([z_min, z_max])

    print('BoneMeshMask')
    return bone


def load_BVTVdata(bone, filename):
    bone_img = ct.load_itk(filename)

    BVTVscaled = bone_img[0]*0.651+0.05646  # scaling factor/intercept from Schenk et al. 2022, has to be discussed w Ph

    # Reorientation of axes
    BVTVscaled = BVTVscaled

    # Flip image 180° to get same COS origin
    BVTVscaled = BVTVscaled[:, :, ::-1]

    bone["BVTVscaled"] = BVTVscaled
    Spacing = bone_img[2]
    bone["Spacing"] = Spacing

    return bone


def HFE_inp_creator(inpDummy, eleSets, material, inpName):
    outfile = open(inpName, 'w')

    f_inpDummy = open(inpDummy)
    f_eleSets = open(eleSets)
    f_material = open(material)

    for lines in f_inpDummy:
        outfile.write(lines)
        if '**ADD ELEMENT SETS' in lines:
            for lines_sets in f_eleSets:
                outfile.write(lines_sets)
        if '**ADD MATERIALS' in lines:
            for lines_mat in f_material:
                outfile.write(lines_mat)
    outfile.close()

    print("Ende HFE_inp_creator")


class IndexTracker(object):
    def __init__(self, ax, Xa):
        self.ax = ax
        ax.set_title('use scroll wheel to navigate images')

        self.Xa = Xa
        rows, cols, self.slices = Xa.shape
        self.ind = self.slices // 2

        self.ima = ax.imshow(self.Xa[:, :, self.ind], interpolation='nearest', cmap="bone")
        # self.cba = plt.colorbar(self.ima)
        # self.cba.set_label('[]')
        self.update()

    def onscroll(self, event):
        print("%s %s" % (event.button, event.step))
        if event.button == 'up':
            self.ind = (self.ind + 1) % self.slices
        else:
            self.ind = (self.ind - 1) % self.slices
        self.update()

    def update(self):
        self.ima.set_data(self.Xa[:, :, self.ind])
        self.ax.set_ylabel('slice %s' % self.ind)
        self.ima.axes.figure.canvas.draw()


def plot3d(imagea):
    fig, ax = plt.subplots(1, 1)
    tracker = IndexTracker(ax, imagea)
    fig.canvas.mpl_connect('scroll_event', tracker.onscroll)
    plt.show()


def controlDirection(bone):
    mask = bone['MASK_array']
    boneCT = bone["BVTVscaled"]

    control_img = mask * boneCT

    plot3d(control_img)

    print('Ende - control direction')

'''
t1 = time.time()
# create an empty bone dictionary to store all important information
bone = {}

# Project path
path_pro = f'/home/biomech/Documents/01_Icotec/'
# Path to mesh dummy and stl - is needed for mapping (node and element coordinates) and mask
path_mesh = f'{path_pro}02_FEA/99_Tests/Pilot3/'
# Path to the bone image - sample
# path_CT = f'C:/Users/kirta/Documents/2022_Fynmann_01_local/04_Segmented_Test_72um'
path_CT = f'{path_pro}01_Experiments/02_Scans/Pilot3/04_Registered/'


# File for mapping
bone_file = f'{path_CT}/XCT_Icotec_S130672_L5_intact_planned.mhd'
bone['SampleName'] = bone_file

# File to generat mask
stl_file = f'{path_mesh}/mesh.stl'

# File to get bone mesh information (nodes, elements and cog)
meshBonnyDummy = f'{path_mesh}mesh.inp'

# Abq input file to generate the file for simulation
inpDummy = f'{path_mesh}Pilot3_mesh.inp'


# Read the bone mesh created in Abaqus - to read out all elements, the connecting nodes and their coordinates
bone = readInpBoneDummy(bone, meshBonnyDummy)

print('a)')

# Create a mask form bone mesh created in Abaqus
bone = boneMeshMask(bone, path_mesh, stl_file, 0.0607, sample + '_mask.mhd')  # , controlplot=True)
print('b)')

mask = ct.load_itk(f'{path_mesh}/BoneTest.mhd')
print('c)')


# Read in the sample for BVTV mapping
# The sample needs to be oriented and should have the same size as the bone-mask!!!
bone = load_BVTVdata(bone, bone_file)
print('--> Files loaded')

bone['Bone_Mask'] = np.zeros(bone['BVTVscaled'].shape)
'''
'''
mask_pos = np.array(np.where(mask[0] == 1))
mask_pos_ = np.zeros_like(mask_pos)
for i in range(len(mask_pos[0])):
    point = np.array(np.append(mask_pos[:, i], 1)).reshape(4, 1)
    mask_pos_[:, i] = (np.round(np.dot(np.linalg.inv(COS_CT), point)[:3].ravel())).astype(int)
print(str(round(time.time()-t1, 2)))
im = sitk.GetImageFromArray(mask_pos_, isVector=False)
im.SetSpacing((0.0607, 0.0607, 0.0607))
im.SetOrigin((-27.9994, -157.164, -38.7857))

rotated_image = ndimage.affine_transform(bone['MASK_array'], COS_CT_inv[:3, :3],
                                         offset=COS_CT_inv[:3, 3], mode='nearest', output=bone['Bone_Mask'])


## Make CT-image to same size as mask
def sameSize(bone):
    mask = bone['MASK_array']
    boneCT = bone["BVTVscaled"]
    same_size = np.zeros_like(mask)

    if mask.shape != boneCT.shape:
        same_size[0:boneCT.shape[0], 0:boneCT.shape[1], 0:boneCT.shape[2]] = boneCT


        plot3d(same_size * mask)
        print('HALLO')

sameSize(bone)

print('1)')

## Check direction
controlDirection(bone)

print('2)')

## Do the mapping between CT-image in Abq-Mesh
## It will create tow separate files. One containing all elements sets and one the related material parameters
bone = EasyHFE_mapping(bone, 'Element_Set.inp', 'Material.inp')
print('3)')

HFE_inp_creator(inpDummy, 'Element_Set.inp', 'Material.inp', 'Test_inpFile.inp')
print('4)')
# BVTVscaled muss an Mesh-Maske angepasst werden
# KOS müssen noch richtig zueinander angeordnet werden!!!! -- offest in cog addieren!


# Maske und BVTV-scaeld müssen noch in die richtige richtung gedreht werden!!
# --> kan wahrscheindlich mit .T gemacht werden


path = f'C:/Users/kirta/OneDrive - Universitaet Bern/01_MBArtorg/2021_Projects/'\
       f'2021_Schroedinger/04_SampleInfos/05_IMG_cropped_orientation'
files = list_txt_files(path)

startpoints = np.array([pd.read_csv(f'{path}/{txt}')['startPoint'] for txt in files])
print(startpoints.max())


print('main')

    # for the mask --> the max startpoint is at position 202 so the mask need to be about XX elements long.
'''
