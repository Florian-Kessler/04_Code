from MedtoolFunctions import medtool_functions as mf
import numpy as np
from scipy.spatial import ConvexHull
import os
import pyvista as pv
import SimpleITK as sitk
import skimage.morphology as morph
import utils_SA as utils
import sys
import matplotlib.pyplot as plt
import ReadRawMHD as rR
# import time
# import scipy


def HFE_mapping_trans(bone, inp):
    """
    Material Mapping, including PSL ghost padding layers as copy of most distal and proximal layers
    :param bone: bone dictionary
    :param inp: inputfile properties
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

    ROI_BVTV_size = inp['Mapping_Diameter']  # sphere diameter in mm

    # print("FEelSize material mapping = " + str(FEelSize))

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

            RHOb[elem] = BVTVbone*0.651+0.05646  # Schenk et al. 2022
            RHOb_FE[elem] = BVTVbone_FE*0.651+0.05646
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
            "\r" + "... material mapping element " + str(i+1) + "/" + str(len(elems))
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

    # mf.writeAbaqus(INPname, None, nodes, None, elems, elsets, NscaResults=None)
    # *****************************************************************
    marray = np.real(np.mean([np.asarray(m[elem]) for elem in m.keys()], axis=0))
    mmarray1 = np.real(np.mean([np.asarray(mm[elem][0]) for elem in m.keys()], axis=0))
    mmarray2 = np.real(np.mean([np.asarray(mm[elem][1]) for elem in m.keys()], axis=0))
    mmarray3 = np.real(np.mean([np.asarray(mm[elem][2]) for elem in m.keys()], axis=0))

    # io_utils.print_mem_usage()

    print("\n... material mapping finished")

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
    print("\n ... update ABAQUS files:\n" + inp['FEA_loc'] + inp['Model_Code'] + inp['Screw'] + '_elsets.inp' +
          "\nand\n" + inp['FEA_loc'] + inp['Model_Code'] + inp['Screw'] + '_materials.inp')
    outfile1 = open(inp['FEA_loc'] + inp['Model_Code'] + inp['Screw'] + '_elsets.inp', 'w')
    outfile2 = open(inp['FEA_loc'] + inp['Model_Code'] + inp['Screw'] + '_materials.inp', 'w')
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


def readInpBoneDummy(bone, inp):
    """
    This function will generate an Abaqus input file which defines for each element in the mesh file a separate
    element-set.
    :param bone: bone dictionary
    :param inp: inp dictionary with abaqus input file details
    :return: An Abaqus input file with the mesh information and element sets
    """

    # Read information from Dummy input file of bone
    inpDummy = mf.readAbaqus(inp['FEA_loc'] + inp['Model_Code'] + '_mesh.inp')
    # title = inp[0]
    nodes = inpDummy[1]
    # nsets = inpDummy[2]
    elems = inpDummy[3]
    elsets = inpDummy[4]

    # calculate the volume of each element
    evol = []
    for i, eleid in enumerate(elems):
        nodeid = elems[eleid].get_nodes()
        node_coo = np.array([nodes[ids].get_coord() for ids in nodeid])
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


def boneMeshMask(bone, inp, controlplot=False, reshape=True, closing=True):
    """
    This function creates a mask form any stl file and returns a 3d array mask - and store the mask as mhd in the given
    path.
    :param bone: bone dictionary
    :param inp: input file properties
    :param controlplot: If true a control 3d image of the stl will pop up - close it to proceed
    :param reshape: sometimes the order of the slices for the 3d array does not match - activate it if mhd looks weird
    :param closing: sometimes there are some small holes in the mask - activate it if needed
    :return: 3d array mask
    """
    # read in the stl to generate the mask
    reader = pv.get_reader(inp['FEA_loc'] + inp['Model_Code'] + '_mesh.stl')
    mesh = reader.read()

    if controlplot:
        mesh.plot()
        voxels = pv.voxelize(mesh, density=inp['Resolution'])
        p = pv.Plotter()
        p.add_mesh(voxels, color=True, show_edges=True)
        p.add_mesh(mesh, color="red", opacity=0.5)
        p.show()

    x_min, x_max, y_min, y_max, z_min, z_max = mesh.bounds
    x = np.arange(x_min, x_max, inp['Resolution'])
    y = np.arange(y_min, y_max, inp['Resolution'])
    z = np.arange(z_min, z_max, inp['Resolution'])
    x, y, z = np.meshgrid(x, y, z)

    # Create unstructured grid from the structured grid
    grid = pv.StructuredGrid(x, y, z)
    ugrid = pv.UnstructuredGrid(grid)

    # get part of the mesh within the mesh's bounding surface.
    selection = ugrid.select_enclosed_points(mesh.extract_surface(), tolerance=0.0, check_surface=False)
    mask_ = selection.point_data['SelectedPoints'].view(np.bool_)

    if reshape:
        # sometimes the order of the matrix gets changed
        mask = mask_.reshape([z.shape[2], z.shape[1], z.shape[0]])
        mask = mask[:, ::-1, :]  # mirroring
        mask = np.rot90(mask, k=-1, axes=(1, 2))
    else:
        mask = mask_.reshape(x.shape)

    mask = np.array(mask)

    if closing:
        # close small holes and gabs:
        for i in range(0, mask.shape[0]):
            mask[i, :, :] = morph.closing(mask[i, :, :], np.ones([3, 3]))
            # mask_copy[i, :, :] = morph.dilation(mask_copy[i, :, :], np.ones([3, 3]))
            # mask_copy[i, :, :] = morph.erosion(mask_copy[i, :, :], np.ones([2, 2]))

    origin = [0, 0, 0]
    spacing = np.array([1, 1, 1]) * inp['Resolution']

    mask_trans = mask.astype(np.short)
    itkmask = sitk.GetImageFromArray(mask_trans, isVector=None)
    itkmask.SetSpacing(spacing)
    itkmask.SetOrigin(origin)

    sitk.WriteImage(itkmask, inp['FEA_loc'] + inp['Model_Code'] + inp['Screw'] + '_mask.mhd')

    # set bone values
    bone['MASK_array'] = mask_trans.T
    # To move COS to right place in image
    bone['MaskX'] = np.array([x_min, x_max])
    bone['MaskY'] = np.array([y_min, y_max])
    bone['MaskZ'] = np.array([z_min, z_max])

    print('BoneMeshMask')
    return bone


def load_BVTVdata(bone, filename):
    """
    Load and segment bone CT
    :param bone: bone dictionary
    :param filename: path to grey image of bone CT
    :return: segmented BVTV image in bone dictionary
    """
    bone["GreyImage"] = sitk.ReadImage(filename)
    print('Image read.')
    # bone["GreyImage"] = scipy.ndimage.gaussian_filter(bone["GreyImage"], sigma=0.8, truncate=1.25)  # Schenk 2022
    # print('Start Gauss filtering')
    # tG = time.time()
    # bone["GreyImage"] = scipy.ndimage.gaussian_filter(bone["GreyImage"], sigma=0.8, truncate=1.25)  # Schenk 2022
    # print('Gauss filtering time: ' + str(int((time.time() - tG) / 60)) + ' min ' + str(
    #    round(np.mod(time.time() - tG, 60), 1)) + ' sec.')

    # Convert the image to a  numpy array first and then shuffle the dimensions to get axis in the order z,y,x
    # Transform image from z,y,x to x,y,z
    bone_img = np.transpose(sitk.GetArrayFromImage(bone["GreyImage"]), [2, 1, 0])
    print('Transpose.')
    # Transform image from z,y,x to x,y,z
    # bone_img = np.transpose(bone_img, [2, 1, 0])

    # Read the origin of the ct_scan, will be used to convert the coordinates from world to voxel and vice versa
    # origin = np.array(list(reversed(itkimage.GetOrigin())))

    # Read the spacing along each dimension
    bone["Spacing"] = np.array(list(reversed(bone["GreyImage"].GetSpacing())))
    print('Spacing.')
    # scaling factor/intercept from Schenk et al. 2022. Threshold of 320 was found in paper for trabecular (cort: 450)
    # BVTVscaled = rR.zeros_and_ones(bone_img, 320)

    # Flip image 180?? to get same COS origin
    # bone["BVTVscaled"] = BVTVscaled  # [:, :, ::-1]

    bone["BVTVscaled"] = rR.zeros_and_ones(bone_img, 320)
    print('Segmented.')
    return bone


def HFE_inp_creator(inp):
    """
    Creates inputfile with specified details
    :param inp: inputfile properties
    :return: writes input file
    """
    SimMat = ['T', 'P']  # simulated materials
    for i in range(len(SimMat)):
        step = 0
        f_inpDummy = open(inp['FEA_loc'] + inp['Model_Code'] + '_model.inp')
        f_eleSets = open(inp['FEA_loc'] + inp['Model_Code'] + inp['Screw'] + '_elsets.inp')
        f_material = open(inp['FEA_loc'] + inp['Model_Code'] + inp['Screw'] + '_materials.inp')
        # Included in input file-name:
        # Model code (geometry), force maximum, friction, experiment screw material, simulation screw material
        if inp['F_dir'] == '-':
            outfile = open(inp['FEA_loc'] + inp['Model_Code'] + '_F' + str(inp['F_max']) + '_' +
                           str(inp['Friction']).replace('.', '') + '_' + inp['Screw'] + '_' + SimMat[i] + '.inp', 'w')
        elif inp['F_dir'] == '+':
            outfile = open(inp['FEA_loc'] + inp['Model_Code'] + '_Finv' + str(inp['F_max']) + '_' +
                           str(inp['Friction']).replace('.', '') + '_' + inp['Screw'] + '_' + SimMat[i] + '.inp', 'w')
        else:
            print('Wrong F_dir input')
            exit()
        for lines in f_inpDummy:

            # Define step 2
            if '*Step, name=Step-2,' in lines:
                step = 2

            # Set material of implant
            if '*Solid Section, elset=Set-Impl, material=PEEK' in lines:
                # replace material section
                if SimMat[i] == 'T':
                    # outfile.truncate()
                    outfile.write('*Solid Section, elset=Set-Impl, material=Ti\n')
                    print('Section set to Ti.')
                elif SimMat[i] == 'P':
                    outfile.write('*Solid Section, elset=Set-Impl, material=PEEK\n')
                    print('Section set to PEEK.')

            # Add bone element sets
            elif '*Solid Section, elset=Set-Bone, material=Bone' in lines:
                outfile.write(lines)
                # add elsets
                for lines_sets in f_eleSets:
                    outfile.write(lines_sets)

            # Add bone element set material properties for UMAT
            elif '0., 0.06,   1.,   1.,   1.' in lines:
                outfile.write(lines)
                # add bone mat
                for lines_mat in f_material:
                    outfile.write(lines_mat)

            # Set force amplitude and direction
            elif 'Set-RP, 2,' in lines:
                if step == 5:  # HERE EDITED FOR TESTING DISPLACEMENT CONTROLLED !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    outfile.write('Set-RP, 2, ' + str(inp['F_dir']) + str(inp['F_max']) + '\n')
                else:
                    outfile.write(lines)

            # Amplitude settings (number of cycles, peak and valley load)
            # here proceed with amplitude etc.

            # lines.replace(old string, new string) could be used

            # elif '*Amplitude, name=Amp-1' in lines:
                # outfile.write('0., 0.,\n')
                # for j in range(inp['Cycles']):
                #    amp = (j+1)/inp['Cycles']

                # outfile.write('
            else:
                outfile.write(lines)

        outfile.close()
        f_inpDummy.close()
        f_eleSets.close()
        f_material.close()
    print("End HFE_inp_creator")


def write_submit(inp):
    """
    Function to write a submit file, either for cortex or for ubelix
    :param inp: inputfile properties
    :return: writes .sh files
    """
    SimMat = ['P', 'T']
    for i in range(len(SimMat)):
        inputfile = ''
        template = open(inp['Project_Folder'] + '02_FEA/97_Templates/abq_submit_template_' + inp['Submit'] + '.sh')
        abq_file = open(inp['FEA_loc'] + 'abq_submit_' + inp['Screw'] + SimMat[i] + '.sh', 'w')
        if inp['F_dir'] == '-':
            inputfile = inp['Model_Code'] + '_F' + str(inp['F_max']) + '_' + str(inp['Friction']).replace('.', '') \
                        + '_' + inp['Screw'] + '_' + SimMat[i]
        elif inp['F_dir'] == '+':
            inputfile = inp['Model_Code'] + '_Finv' + str(inp['F_max']) + '_' + str(inp['Friction']).replace('.', '') \
                        + '_' + inp['Screw'] + '_' + SimMat[i]
        material_code = inp['Screw'] + SimMat[i]
        for lines in template:
            if 'inpname' in lines:
                new_line = lines.replace('inpname', inputfile)
                abq_file.write(new_line)
            elif 'materials' in lines:
                new_line = lines.replace('materials', material_code)
                abq_file.write(new_line)
            else:
                abq_file.write(lines)
        abq_file.close()
        os.chmod(inp['FEA_loc'] + 'abq_submit_' + inp['Screw'] + SimMat[i] + '.sh', 0o744)


def write_mesh(inp):
    """
    Extracts the bone mesh from the input file. Part must be called name=Bone
    :param inp: input file dictionary
    :return: no return variable, writes an input file containing the bone mesh only
    """
    orig = open(inp['FEA_loc'] + inp['Model_Code'] + '_model.inp')
    mesh = open(inp['FEA_loc'] + inp['Model_Code'] + '_mesh.inp', 'w')
    start = 0
    for lines in orig:
        if '*Part, name=Bone' in lines:
            start = 1
        if start == 1:
            if '*Nset' in lines:
                print('Finished extracting mesh file.')
                break
            mesh.write(lines)
    mesh.close()


class IndexTracker(object):
    def __init__(self, ax, xa):
        self.ax = ax
        ax.set_title('use scroll wheel to navigate images')

        self.xa = xa
        rows, cols, self.slices = xa.shape
        self.ind = self.slices // 2

        self.ima = ax.imshow(self.xa[:, :, self.ind], interpolation='nearest', cmap="bone")
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
        self.ima.set_data(self.xa[:, :, self.ind])
        self.ax.set_ylabel('slice %s' % self.ind)
        self.ima.axes.figure.canvas.draw()


def plot3d(image_a):
    fig, ax = plt.subplots(1, 1)
    tracker = IndexTracker(ax, image_a)
    fig.canvas.mpl_connect('scroll_event', tracker.onscroll)
    plt.show()


def controlDirection(bone):
    mask = bone['MASK_array']
    boneCT = bone["BVTVscaled"]

    control_img = mask * boneCT

    plot3d(control_img)

    print('End - control direction')
