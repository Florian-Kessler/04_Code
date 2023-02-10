"""
This script contains input/output functions.

Author: Jarunan Panyasantisuk, ETH Scientific IT Services
Date: 15 November 2018
"""

from __future__ import print_function

import os
import yaml
import socket

import numpy as np

import utils_SA as utils

debug = 1


def read_config_file(filename):
    """ read intialization values """

    print(" ... reading initialization file", filename)
    with open(filename, "r") as f:
        config = yaml.load(f, Loader=yaml.FullLoader)

    return config


def set_umat_parameters(config):
    """
    compute umat parameters
    (is not influencing UMAT parameters in fortran file at the moment, so obsolete)
    """

    parameters = config["umat_parameters"]

    parameters["E0"] = 9759.0 * (12000.0 / 9759.0) * config["sca1"]
    parameters["MU0"] = 3117.0 * (12000.0 / 9759.0) * config["sca1"]
    parameters["SIGD0P"] = 57.69 * config["sca2"]
    parameters["SIGD0N"] = 73.1 * config["sca2"]
    parameters["TAUD0"] = 29.61 * config["sca2"]

    return parameters


def set_filenames(config, sample, pipeline="fast"):
    """
    set filenames for each grayscale file
    filenames depend on pipeline (fast/accurate)
    Always:
        - native image for header
        - BMD or Native image croppped to ROI
    additional for fast pipeline:
        - periosteal mask
    additional for accurate pipeline:
        - trabecular mask
        - cortical mask
        - two-phase segmentation
    """

    # Always, not depending on nphase
    current_version = config["current_version"]
    folder_id = config["folder_id"]
    folder = folder_id[sample]
    aimdir = os.path.join(config["aimdir"], folder)
    origaimdir = os.path.join(config["origaimdir"], folder)

    filename_postfix_bmd = config["filename_postfix_bmd"]

    filename = {}
    filename["FILEBMD"] = "{}{}".format(sample, filename_postfix_bmd)
    filename["FILEGRAY"] = "{}{}".format(sample, ".AIM")
    filename["RAWname"] = os.path.join(origaimdir, filename["FILEGRAY"])
    filename["BMDname"] = os.path.join(origaimdir, filename["FILEBMD"])
    filename["boundary_conditions"] = config["boundary_conditions"]

    # additional for fast pipeline
    if pipeline == "fast":
        filename_postfix_mask = config["filename_postfix_mask"]
        filename["FILEMASK"] = "{}{}".format(sample, filename_postfix_mask)
        filename["MASKname"] = os.path.join(origaimdir, filename["FILEMASK"])

    # Additionl for accurate pipeline
    if pipeline == "accurate":
        filename_postfix_trab_mask = config["filename_postfix_trab_mask"]
        filename_postfix_cort_mask = config["filename_postfix_cort_mask"]
        filename_postfix_seg = config["filename_postfix_seg"]

        filename["FILEMASKCORT"] = "{}{}".format(sample, filename_postfix_cort_mask)
        filename["FILEMASKTRAB"] = "{}{}".format(sample, filename_postfix_trab_mask)
        filename["FILESEG"] = "{}{}".format(sample, filename_postfix_seg)

        filename["CORTMASKname"] = os.path.join(origaimdir, filename["FILEMASKCORT"])
        filename["TRABMASKname"] = os.path.join(origaimdir, filename["FILEMASKTRAB"])
        filename["SEGname"] = os.path.join(origaimdir, filename["FILESEG"])

    # General filenames
    print(filename["BMDname"])
    new_filename = "{}_V_{}.inp".format(sample, current_version)
    filename["INPname"] = os.path.join(aimdir, new_filename)

    new_filename = "{}_V_{}_summary.txt".format(sample, current_version)
    filename["SUMname"] = os.path.join(aimdir, new_filename)

    new_filename = "{}_V_{}_BPVb".format(sample, current_version)
    filename["VER_BPVname"] = os.path.join(aimdir, new_filename)

    return filename


def read_aim_old(filenames):
    """
    read AIM files
    --------------
    All necessary AIM files are imported.
    The scaling is taken from an original AIM file
    for Data processed with medtool, as the header is missing.
    For later hFE simulations (NODARATIS),
    the header can be read directly from the BMD image
    depending on the CASE, MASK, COR, TRA and SEG are imported
    """

    print("\n ... read AIM files")

    """
    read image informations from raw image
    (only for files processed in medtool
    """
    XXX, Spacing, Scaling, Slope, Intercept, Header = utils.AIMreader(
        filenames["RAWname"], 0
    )
    if debug:
        print("\n ... ... reading original AIM completed\n")

    BMD_vtk = utils.AIMreader(filenames["BMDname"], Spacing)[0]
    if debug:
        print("\n ... ... reading BMD image completed\n")  # read BMD image

    MASK_vtk = utils.AIMreader(filenames["MASKname"], Spacing)[0]
    if debug:
        print("\n ... ... reading mask completed\n")


    """
    convert AIM files to numpy arrays
    ---------------------------------
    All vtk files are converted to numpy arrays and
    vtk arrays are deleted to save memory
    """

    print("\n ... AIM data to numpy arrays")
    BMD_array = utils.vtk2numpy(BMD_vtk)
    MASK_array = utils.vtk2numpy(MASK_vtk)


    bone = {}
    bone["Spacing"] = Spacing
    bone["Scaling"] = Scaling
    bone["Slope"] = Slope
    bone["Intercept"] = Intercept
    bone["BMD_array"] = BMD_array
    bone["MASK_array"] = MASK_array

    print("Spacing = "+str(bone["Spacing"]))
    print("Scaling = "+str(bone["Scaling"]))
    print("Slope = "+str(bone["Slope"]))
    print("Intercept = "+str(bone["Intercept"]))

    return bone


def read_aim(name, filenames, bone):
    """
    read AIM image
    --------------
    All necessary AIM files are imported and stored in bone dict
    Input: name specifier, filenames dict, bone dict
    Output: bone dict
    - numpy array containing AIM image
    """

    print("\n ... read file :" + name)

    Spacing = bone["Spacing"]
    # Read image as vtk
    IMG_vtk = utils.AIMreader(filenames[name + 'name'], Spacing)[0]
    # convert AIM files to numpy arrays
    IMG_array = utils.vtk2numpy(IMG_vtk)
    bone[name + "_array"] = IMG_array

    return bone


def log_summary(bone, config, filenames, var):

    # bone
    BMD_array = bone["BMD_array"]
    Slope = bone["Slope"]
    Intercept = bone["Intercept"]
    Spacing = bone["Spacing"]
    Scaling = bone["Scaling"]
    elems = bone["elems"]
    elems_bone = bone["elems_BONE"]
    FEelSize = bone["FEelSize"]
    marray = bone["marray"]
    mmarray1 = bone["mmarray1"]
    mmarray2 = bone["mmarray2"]
    mmarray3 = bone["mmarray3"]

    # config
    SCA1 = config["sca1"]
    SCA2 = config["sca2"]
    KMAX = config["kmax"]

    # PSL
    mode_ghost_layer = config["mode_ghost_layer"]
    n_ghost_proximal = config["padding_elements_proximal"]
    n_ghost_distal = config["padding_elements_distal"]

    # filenames
    BMDname = filenames["BMDname"]
    SUMname = filenames["SUMname"]

    # Variables
    mean_BMD = var["mean_BMD"]
    BVTV_tissue = var["BVTV_tissue"]
    mask_volume_MASK = var["Mask_Volume_MASK"]
    mask_volume_FE = var["Mask_Volume_FE"]
    mask_volume_quality = var["Volume_ratio"]
    BV_tissue = var["BV_tissue"]
    BVTV_FE_tissue_ROI = var["simulation_BVTV_FE_tissue_ROI"]
    BVTV_FE_tissue_ELEM = var["simulation_BVTV_FE_tissue_ELEM"]

    BVTV_FE_elements_ROI = var["simulation_BVTV_FE_elements_ROI"]
    BVTV_FE_elements_ELEM = var["simulation_BVTV_FE_elements_ELEM"]

    BVTV_quality_ROI = var["BVTV_ratio_ROI"]
    BVTV_quality_ELEM = var["BVTV_ratio_ELEM"]

    BMC_tissue = var["BMC_tissue"]
    BMC_FE_tissue_ROI = var["simulation_BMC_FE_tissue_ROI"]
    BMC_FE_tissue_ROI_orig = var["simulation_BMC_FE_tissue_ROI_orig"]
    BMC_FE_tissue_ELEM = var["simulation_BMC_FE_tissue_ELEM"]
    BMC_quality_ROI = var["BMC_ratio_ROI"]
    BMC_quality_ELEM = var["BMC_ratio_ELEM"]
    BMC_compensation = BMC_tissue/BMC_FE_tissue_ROI_orig

    summary = "\n".join(
        [
            """
******************************************************************
**                         SUMMARY FILE                         **
**                hFE pipeline Denis Schenk 2018                **
******************************************************************""",
            "File                 : {}".format(BMDname),
            "System computed on   : {}".format(socket.gethostname()),
            "Simulation Type      : Fast model (one phase / isotropic)",
            "Fitting Variables    : {:.2f} | {:.2f} | {:.3f} |".format(
                SCA1, SCA2, KMAX
            ),
            "*****************************************************************",
            "Patient specific loading",
            "-------------------------------------------------------------------",
            "Ghost layer mode     : {}".format(mode_ghost_layer),
            "Ghost layers prox.   : {}".format(n_ghost_proximal),
            "Ghost layers dist.   : {}".format(n_ghost_distal),
            "*****************************************************************",
            "Image Dimension      : {}".format(np.shape(BMD_array)),
            "Scaling              : {}".format(Scaling),
            "Slope                : {:.3f}".format(Slope),
            "Intercept            : {:.3f}".format(Intercept),
            "Spacing              : {:.3f}, {:.3f}, {:.3f} mm".format(*Spacing),
            "FE element size      : {:.3f}, {:.3f}, {:.3f} mm".format(*FEelSize),
            "Number of elements   : {:d} + {:d}".format(len(elems_bone), len(elems)-len(elems_bone)),
            "******************************************************************",
            "Apparent variables computed from original BMD, BVTV images and mask",
            "-------------------------------------------------------------------",
            "Apparent BMD grays.  : {:.2f} mgHA/ccm".format(mean_BMD),
            "Apparent BVTV grays. : {:.3f} %".format(BVTV_tissue * 100.0),
            "Apparent BMC grays.  : {:.1f} mgHA".format(BMC_tissue),
            "Apparent BV grays.   : {:.1f} mm^3".format(BV_tissue),
            "******************************************************************",
            "Volumes of mask",
            "-------------------------------------------------------------------",
            "Mask Volume image    : {:.1f} mm^3".format(mask_volume_MASK),
            "Mask Volume FE mesh  : {:.1f} mm^3".format(mask_volume_FE),
            "Mask Volume ratio    : {:.3f} ".format(mask_volume_quality),
            "******************************************************************",
            "BVTV",
            "-------------------------------------------------------------------",
            "BVTV tissue           : {:.2f} %".format(BVTV_tissue * 100.0),
            "BVTV FE tissue ROI    : {:.2f} %".format(BVTV_FE_tissue_ROI * 100.0),
            "BVTV FE tissue ELEM   : {:.2f} %".format(BVTV_FE_tissue_ELEM * 100.0),
            "BVTV FE elements ROI  : {:.2f} %".format(BVTV_FE_elements_ROI * 100.0),
            "BVTV FE elements ELEM : {:.2f} %".format(BVTV_FE_elements_ELEM * 100.0),
            "BVTV ratio ROI        : {:.3f} ".format(BVTV_quality_ROI),
            "BVTV ratio ELEM       : {:.3f} ".format(BVTV_quality_ELEM),
            "******************************************************************",
            "BMC",
            "-------------------------------------------------------------------",
            "BMC tissue            : {:.1f} mgHA".format(BMC_tissue),
            "BMC FE tissue ELEM    : {:.1f} mgHA".format(BMC_FE_tissue_ELEM),
            "BMC FE tissue ROI orig: {:.1f} mgHA".format(BMC_FE_tissue_ROI_orig),
            "Compensated for BMC conservation:",
            "BMC FE tissue ROI     : {:.1f} mgHA".format(BMC_FE_tissue_ROI),
            "BMC ratio ROI         : {:.3f} ".format(BMC_quality_ROI),
            "BMC ratio ELEM        : {:.3f} ".format(BMC_quality_ELEM),
            "BMC comp reco/sim     : {:.3f} ".format(BMC_compensation),
            "******************************************************************",
            "Average anisotropy",
            "-------------------------------------------------------------------",
            "Eigenvalues (max, mid, min)    : {0} {1} {2}".format(*marray),
            "Degree of anisotropy (max/min) : {:.3f}".format(marray[0] / marray[2]),
            "Eigenvector max (x,y,z)        : {}".format(mmarray1),
            "Eigenvector mid (x,y,z)        : {}".format(mmarray2),
            "Eigenvector min (x,y,z)        : {}".format(mmarray3),
            "******************************************************************",
            "Summary Benchmark Tests",
            "-------------------------------------------------------------------",
            "Benchmark tests compare variables computed from the full volume",
            "against variables computed from the FE elements",
            "Volume Image/Subvolumes        : {:.3f}".format(mask_volume_quality),
            "BVTV ROI Image/Subvolumes      : {:.3f}".format(BVTV_quality_ROI),
            "BVTV ELEM Image/Subvolumes     : {:.3f}".format(BVTV_quality_ELEM),
            "BMC ROI Image/Subvolumes       : {:.3f}".format(BMC_quality_ROI),
            "BMC ELEM Image/Subvolumes      : {:.3f}".format(BMC_quality_ELEM),
            "******************************************************************",
        ]
    )

    print(summary)

    with open(SUMname, "w") as sumFile:
        sumFile.write(summary)


def log_append_processingtime(filename, time):
    SUMname = filename
    time_summary = "\n".join(
        ["Summary Processing Time",
         "Full processing time           : {:.3f} [s]".format(time),
         "******************************************************************",])
    print(time_summary)

    with open(SUMname, "a") as sumUpdate:
        sumUpdate.write("\n")
        sumUpdate.write(time_summary)
    print("... added processing time to summary file")


def datfilereader(infilename):
    outfilename = infilename.replace(".dat", ".txt")

    # Read the .dat file
    with open(infilename, "r") as infile:
        lines = infile.readlines()

    # Lists used
    ref_nodedata = []
    # Node: RF3 & U3
    U3 = " "
    RF3 = " "
    UR1 = " "
    RM1 = " "

    j = 0
    for i in range(0, len(lines)):
        if lines[i].find("U3") > -1:
            j = j + 1  # value of increment
            member1 = lines[i + 3].split(" ")  # value of U3
            member2 = lines[i + 12].split(" ")  # value of RF3
            for k in range(0, len(member1)):
                if member1[k].find(".") != -1:
                    U3 = member1[k]
            for k in range(0, len(member2)):
                if member2[k].find(".") != -1:
                    if member2[k].find("+") != -1:
                        if member2[k].find("E+") != -1:
                            RF3 = member2[k]
                        else:
                            RF3 = member2[k].replace("+", "E+")
                    else:
                        RF3 = member2[k]
            ref_nodedata.append(
                str(j) + " " + U3.replace("\n", "") + " " + RF3.replace("\n", "")
            )  # absolute values are taken
    # Create the output files
    with open(outfilename, "w") as out:
        out.write(
            """
*********************************************************
* Displacement and reaction force of the reference node *
*                 Ghislain MAQUER 2011                  *
*********************************************************
inc U3 RF3
0 0 0
"""
        )
        for i in range(0, len(ref_nodedata)):
            out.write(ref_nodedata[i] + "\n")


def read_img_param(filenames, bone):
    """
    Read image pararmeters from AIM image header.
    Input: AIM image (Scanco Medical)
    Output: (bone dictionary)
    - Spacing
    - Scaling
    - Slope
    - Intercept
    """
    print("\n ... read AIM files")

    """
    read image informations from raw image
    (only for files processed in medtool
    """
    XXX, Spacing, Scaling, Slope, Intercept, Header = utils.AIMreader(
        filenames["RAWname"], 0
    )

    bone["Spacing"] = Spacing
    bone["Scaling"] = Scaling
    bone["Slope"] = Slope
    bone["Intercept"] = Intercept

    return bone
