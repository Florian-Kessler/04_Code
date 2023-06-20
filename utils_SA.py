from __future__ import print_function

import sys
import struct
import numpy
import numpy as np
import scipy
import vtk
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
import itertools as it
from itertools import chain
from MedtoolFunctions import medtool_functions as mf
import os
import json
from csv import DictWriter, writer
import matplotlib.pyplot as plt

# *****************************************************************
# I. Functions
# *****************************************************************


def sphere_array(shape, radius, position):

    semisizes = (float(radius),) * 3
    grid = [slice(-x0, dim - x0) for x0, dim in zip(position, shape)]
    position = numpy.ogrid[grid]
    arr = numpy.zeros(numpy.asarray(shape).astype(int), dtype=float)
    for x_i, semisize in zip(position, semisizes):
        arr += numpy.abs(x_i / semisize) ** 2
    return (arr <= 1.0).astype("int")


def vtk2numpy(imvtk):
    """turns a vtk image data into a numpy array"""
    dim = imvtk.GetDimensions()
    data = imvtk.GetPointData().GetScalars()
    imnp = vtk_to_numpy(data)
    # vtk and numpy have different array conventions
    imnp = imnp.reshape(dim[2], dim[1], dim[0])
    imnp = imnp.transpose(2, 1, 0)
    return imnp


def numpy2vtk(imnp, spacing):
    """turns a numpy array into a vtk image data"""
    # vtk and numpy have different array conventions
    imnp_flat = imnp.transpose(2, 1, 0).flatten()
    if imnp.dtype == "int8":
        arraytype = vtk.VTK_CHAR
    elif imnp.dtype == "int16":
        arraytype = vtk.VTK_SHORT
    else:
        arraytype = vtk.VTK_FLOAT
    imvtk = numpy_to_vtk(num_array=imnp_flat, deep=True, array_type=arraytype)
    image = vtk.vtkImageData()
    image.SetDimensions(imnp.shape)
    image.SetSpacing(spacing)
    points = image.GetPointData()
    points.SetScalars(imvtk)
    return image


def numpy2mhd(imnp, spacing, filename, header=None):
    """writes a numpy array in metaImage (mhd+raw)"""
    # turns the numpy array to vtk array
    writer = vtk.vtkMetaImageWriter()
    try:
        writer.SetInputData(numpy2vtk(imnp, spacing))
    except:
        writer.SetInput(numpy2vtk(imnp, spacing))
    # writes it as a mhd+raw format
    writer.SetFileName(filename)
    writer.Write()
    # writer AIM header if provided
    if header is not None:
        with open(filename, "a") as f:
            f.write(
                """
!-------------------------------------------------------------------------------
                               AIM Log  
!-------------------------------------------------------------------------------"""
            )
            for line in header:
                f.write(line + "\n")


def numpy2mhdNoWrite(imnp, spacing, filename, header=None):
    """writes a numpy array in metaImage (mhd+raw)"""
    # turns the numpy array to vtk array
    writer = vtk.vtkMetaImageWriter()
    try:
        writer.SetInputData(numpy2vtk(imnp, spacing))
    except:
        writer.SetInput(numpy2vtk(imnp, spacing))
    # writes it as a mhd+raw format
    writer.SetFileName(filename)
    # writer.Write()
    # writer AIM header if provided
    # if header is not None:
    #     f = open(filename, 'a')
    #     f.write("\n!-------------------------------------------------------------------------------\n")
    #     f.write("                                   AIM Log                                       ")
    #     f.write("\n!-------------------------------------------------------------------------------\n")
    #     for line in header: f.write(line + "\n")
    #     f.close()


def ext(filename, new_ext):
    """changes the file extension"""
    return filename.replace("." + filename.split(".")[-1], new_ext)


def get_AIM_ints(f):
    """Function by Glen L. Niebur, University of Notre Dame (2010)
    reads the integer data of an AIM file to find its header length"""
    nheaderints = 32
    nheaderfloats = 8
    f.seek(0)
    binints = f.read(nheaderints * 4)
    header_int = struct.unpack("=32i", binints)
    return header_int


def AIMreader(fileINname, Spacing):
    """reads an AIM file and provides the corresponding vtk image with spacing, calibration data and header"""
    # read header
    print("     " + fileINname)
    with open(fileINname, 'rb') as f:
        AIM_ints = get_AIM_ints(f)
        # check AIM version
        if int(AIM_ints[5]) == 16:
            print("     -> version 020")
            if int(AIM_ints[10]) == 131074:
                format = "short"
                print("     -> format " + format)
            elif int(AIM_ints[10]) == 65537:
                format = "char"
                print("     -> format " + format)
            elif int(AIM_ints[10]) == 1376257:
                format = "bin compressed"
                print("     -> format " + format + " not supported! Exiting!")
                exit(1)
            else:
                format = "unknown"
                print("     -> format " + format + "! Exiting!")
                exit(1)
            header = f.read(AIM_ints[2])
            header_len = len(header) + 160
            extents = (0, AIM_ints[14] - 1, 0, AIM_ints[15] - 1, 0, AIM_ints[16] - 1)
        else:
            print("     -> version 030")
            if int(AIM_ints[17]) == 131074:
                format = "short"
                print("     -> format " + format)
            elif int(AIM_ints[17]) == 65537:
                format = "char"
                print("     -> format " + format)
            elif int(AIM_ints[17]) == 1376257:
                format = "bin compressed"
                print("     -> format " + format + " not supported! Exiting!")
                exit(1)
            else:
                format = "unknown"
                print("     -> format " + format + "! Exiting!")
                exit(1)
            header = f.read(AIM_ints[8])
            header_len = len(header) + 280
            extents = (0, AIM_ints[24] - 1, 0, AIM_ints[26] - 1, 0, AIM_ints[28] - 1)

    # collect data from header if existing
    # header = re.sub('(?i) +', ' ', header)
    header = header.split('\n'.encode())
    header.pop(0)
    header.pop(0)
    header.pop(0)
    header.pop(0)
    Scaling = None
    Slope = None
    Intercept = None
    IPLPostScanScaling = 1
    for line in header:
        if line.find(b"Orig-ISQ-Dim-p") > -1:
            origdimp = ([int(s) for s in line.split(b" ") if s.isdigit()])

        if line.find("Orig-ISQ-Dim-um".encode()) > -1:
            origdimum = ([int(s) for s in line.split(b" ") if s.isdigit()])

        if line.find("Orig-GOBJ-Dim-p".encode()) > -1:
            origdimp = ([int(s) for s in line.split(b" ") if s.isdigit()])

        if line.find("Orig-GOBJ-Dim-um".encode()) > -1:
            origdimum = ([int(s) for s in line.split(b" ") if s.isdigit()])

        if line.find("Scaled by factor".encode()) > -1:
            Scaling = float(line.split(" ".encode())[-1])
        if line.find("Density: intercept".encode()) > -1:
            Intercept = float(line.split(" ".encode())[-1])
        if line.find("Density: slope".encode()) > -1:
            Slope = float(line.split(" ".encode())[-1])
            # if el_size scale was applied, the above still takes the original voxel size. This function works
            # only if an isotropic scaling was applied!!!!
        if line.find("scale".encode()) > -1:
            IPLPostScanScaling = float(line.split(" ".encode())[-1])
    # Spacing is calculated from Original Dimensions. This is wrong, when the images were coarsened and
    # the voxel size is not anymore corresponding to the original scanning resolution!

    try:
        Spacing = IPLPostScanScaling * (
            numpy.around(numpy.asarray(origdimum) / numpy.asarray(origdimp) / 1000, 5)
        )
    except:
        pass
    # read AIM
    reader = vtk.vtkImageReader2()
    reader.SetFileName(fileINname)
    reader.SetDataByteOrderToLittleEndian()
    reader.SetFileDimensionality(3)
    reader.SetDataExtent(extents)
    reader.SetHeaderSize(header_len)
    if format == "short":
        reader.SetDataScalarTypeToShort()
    elif format == "char":
        reader.SetDataScalarTypeToChar()
    reader.SetDataSpacing(Spacing)
    reader.Update()
    imvtk = reader.GetOutput()
    return imvtk, Spacing, Scaling, Slope, Intercept, header


def computeBVTV_FEel(cog, Spacing, FE_elsize_mm, imarray, maskarray):
    """computes BVTV from a numpy array containing the BVTV values for a region of size 'ROIsize' centered in the
    center of gravity of the element"""

    # Cut out ROI from image array
    x, y, z = cog / Spacing
    FE_elsize = FE_elsize_mm[0] / Spacing[0]

    # ROImask_sphere = sphere_array(numpy.shape(imarray), FE_elsize / 2, [x, y, z])

    X = [max(x - FE_elsize / 2, 0), min(x + FE_elsize / 2, imarray.shape[0])]
    Y = [max(y - FE_elsize / 2, 0), min(y + FE_elsize / 2, imarray.shape[1])]
    Z = [max(z - FE_elsize / 2, 0), min(z + FE_elsize / 2, imarray.shape[2])]
    ROI = imarray[
        int(numpy.rint(X[0])) : int(numpy.rint(X[1])),
        int(numpy.rint(Y[0])) : int(numpy.rint(Y[1])),
        int(numpy.rint(Z[0])) : int(numpy.rint(Z[1])),
    ]

    # The ROI for BVTV computation corresponds to the homogenized element
    ROImask = maskarray[
        int(numpy.rint(X[0])) : int(numpy.rint(X[1])),
        int(numpy.rint(Y[0])) : int(numpy.rint(Y[1])),
        int(numpy.rint(Z[0])) : int(numpy.rint(Z[1])),
    ]

    BVTV_FE = numpy.mean(ROI[ROImask != 0])

    if numpy.isnan(BVTV_FE):
        BVTV_FE = 0.0

    return BVTV_FE
    # return BVTVinMask


def computeBVTV_twophase(cog, Spacing, ROIsize_cort_mm, ROIsize_trab_mm, imarray, cortmask, trabmask, phicort, phitrab):
    """computes BVTV from a numpy array containing the BVTV values for a region of size 'ROIsize' centered in the
    center of gravity of the element"""
    # ROI size: If PHItrab = 0, ROI size = sphere with equal volume as FE element
    #           If PHItrab =! 0, ROI size = ROIsize_mm

    # for mixed elements: BVTV_trab computes BVTV from trab and cort, while BVTV_cort only accounts values in cort mask
    # Cut out ROI from image array
    x, y, z = cog / Spacing

    ROIsize_trab = ROIsize_trab_mm / Spacing[0]

    X = [max(x - ROIsize_trab / 2, 0), min(x + ROIsize_trab / 2, imarray.shape[0])]
    Y = [max(y - ROIsize_trab / 2, 0), min(y + ROIsize_trab / 2, imarray.shape[1])]
    Z = [max(z - ROIsize_trab / 2, 0), min(z + ROIsize_trab / 2, imarray.shape[2])]
    ROI = imarray[
        int(numpy.rint(X[0])) : int(numpy.rint(X[1])),
        int(numpy.rint(Y[0])) : int(numpy.rint(Y[1])),
        int(numpy.rint(Z[0])) : int(numpy.rint(Z[1])),
    ]

    ROI_cort_mask = cortmask[
        int(numpy.rint(X[0])) : int(numpy.rint(X[1])),
        int(numpy.rint(Y[0])) : int(numpy.rint(Y[1])),
        int(numpy.rint(Z[0])) : int(numpy.rint(Z[1])),
    ]
    ROI_trab_mask = trabmask[
        int(numpy.rint(X[0])) : int(numpy.rint(X[1])),
        int(numpy.rint(Y[0])) : int(numpy.rint(Y[1])),
        int(numpy.rint(Z[0])) : int(numpy.rint(Z[1])),
    ]

    mean_BVTV_trab = 0.0
    mean_BVTV_cort = 0.0

    # calculate center of sphere in new image
    xc = x - X[0]
    yc = y - Y[0]
    zc = z - Z[0]

    if phitrab > 0.0:
        # Compute trabecular BVTV
        # ------------------------

        # create masking array
        ROImask_sphere_trab = sphere_array(
            numpy.shape(ROI), ROIsize_trab / 2, [xc, yc, zc]
        )

        BVTVimg_trab = ROI[ROImask_sphere_trab + ROI_cort_mask + ROI_trab_mask == 2]

        mean_BVTV_trab = numpy.mean(BVTVimg_trab)

        # check for meaningfull output
        if numpy.isnan(mean_BVTV_trab):
            mean_BVTV_trab = 0.0
        if mean_BVTV_trab > 1:
            mean_BVTV_trab = 1

    # BVTVimg = imarray[numpy.logical_and(ROImask_sphere_trab == 1, (cortmask + trabmask) >= 1)]
    # imarray[cortmask + trabmask == 0] = 0
    # mean_BVTV_trab = numpy.mean(imarray[imarray != 0])

    if phicort > 0.0:
        # Compute cortical BVTV
        # ------------------------
        ROIsize_cort = ROIsize_cort_mm / Spacing[0]
        ROImask_sphere_cort = sphere_array(
            numpy.shape(ROI), ROIsize_cort / 2, [xc, yc, zc]
        )
        BVTVimg_cort = ROI[ROImask_sphere_cort + ROI_cort_mask == 2]

        mean_BVTV_cort = numpy.mean(BVTVimg_cort)

        # check for meaningfull output
        if numpy.isnan(mean_BVTV_cort):
            mean_BVTV_cort = 0.0
        if mean_BVTV_cort > 1:
            mean_BVTV_cort = 1

        # imarray[ROImask_sphere_cort == 0] = 0
        # imarray[cortmask == 0] = 0
        # mean_BVTV_cort = numpy.mean(imarray[imarray != 0])

    # if numpy.isnan(mean_BVTV_trab):
    #     mean_BVTV_trab = 0
    #
    # if numpy.isnan(mean_BVTV_cort):
    #         mean_BVTV_cort = 0

    return mean_BVTV_cort, mean_BVTV_trab
    #
    # x, y, z = cog / Spacing
    # ROIsize = ROIsize_mm / Spacing[0]
    # X = [max(x - ROIsize / 2, 0), min(x + ROIsize / 2, imarray.shape[0])]
    # Y = [max(y - ROIsize / 2, 0), min(y + ROIsize / 2, imarray.shape[1])]
    # Z = [max(z - ROIsize / 2, 0), min(z + ROIsize / 2, imarray.shape[2])]
    # ROI = imarray[int(numpy.rint(X[0])):int(numpy.rint(X[1])), int(numpy.rint(Y[0])):int(numpy.rint(Y[1])),
    #       int(numpy.rint(Z[0])):int(numpy.rint(Z[1]))]
    #
    #
    # # The ROI for BVTV computation corresponds to the homogenized element
    # ROImask = maskarray[int(numpy.rint(X[0])):int(numpy.rint(X[1])), int(numpy.rint(Y[0])):int(numpy.rint(Y[1])),
    #           int(numpy.rint(Z[0])):int(numpy.rint(Z[1]))]
    # ROImask_mean = numpy.mean(ROImask)
    # # In the following, several types of BVTV calculation are implemented
    #
    # # - all voxels are considered, inclusive background (not set to zero!)
    # #   This makes the most sense, as lower resolution is changing this BVTV the least
    # BVTVall = numpy.mean(ROI)
    #
    # # - only voxels inside the mask are considered
    # #   Must be corrected with partial volume
    # BVTVinMask = numpy.mean(ROI[ROImask != 0])
    #
    # # - all voxels, but background (outside mask) voxels outside mask are set to zero
    # #   Be aware, that this changes the histogram and introduce a bias because of image noise (positive and negative)
    # ROInoBG = numpy.copy(ROI)
    # ROInoBG[ROImask == 0] = 0
    # BVTVallnoBG = numpy.mean(ROInoBG)
    #
    # # - all voxels, but negative values are set to zero before AVG (in and outside of the mask)
    # #   This only makes sense, if a lot of air bubbles are visible in the image
    # ROInoNeg = numpy.copy(ROI)
    # ROInoNeg[ROInoNeg < 0] = 0
    # BVTVnoNeg = numpy.mean(ROInoNeg)
    #
    # if numpy.isnan(BVTVall):
    #     BVTVall = 0.0
    # if numpy.isnan(BVTVinMask):
    #     BVTVinMask = 0.0
    # if numpy.isnan(BVTVallnoBG):
    #     BVTVallnoBG = 0.0
    # if numpy.isnan(BVTVnoNeg):
    #     BVTVnoNeg = 0.0
    #
    # del ROInoNeg
    # del ROInoBG
    #
    # return BVTVinMask, BVTVall, BVTVallnoBG, BVTVnoNeg
    # # return BVTVinMask


def computeBVTV_onephase(cog, Spacing, ROIsize_mm, imarray, mask, phi=1.0):
    """computes BVTV from a numpy array containing the BVTV values for a region of size 'ROIsize' centered in the
    center of gravity of the element

    The sphere will only take elements which are in the mask. Zero elements in the mask will be ignored for the BVTV
    calculation.
    """
    # ROI size = sphere with 6.6mm diameter

    # for mixed elements: BVTV_trab computes BVTV from trab and cort, while BVTV_cort only accounts values in cort mask
    # Cut out ROI from image array
    x, y, z = cog / Spacing

    # Number of elements representing the diameter of the sphere
    ROIsize = ROIsize_mm / Spacing[0]
    # Define the position of the element in the image
    X = [max(x - ROIsize / 2, 0), min(x + ROIsize / 2, imarray.shape[0])]
    Y = [max(y - ROIsize / 2, 0), min(y + ROIsize / 2, imarray.shape[1])]
    Z = [max(z - ROIsize / 2, 0), min(z + ROIsize / 2, imarray.shape[2])]

    # Separated region of interest form the bone image
    ROI = imarray[
        int(numpy.rint(X[0])) : int(numpy.rint(X[1])),
        int(numpy.rint(Y[0])) : int(numpy.rint(Y[1])),
        int(numpy.rint(Z[0])) : int(numpy.rint(Z[1])),
    ]

    # Separate the same region of interest from the mask
    ROI_mask = mask[
        int(numpy.rint(X[0])) : int(numpy.rint(X[1])),
        int(numpy.rint(Y[0])) : int(numpy.rint(Y[1])),
        int(numpy.rint(Z[0])) : int(numpy.rint(Z[1])),
    ]
    print('\nROI, ROI_mask: \n')
    print(X)
    print(Y)
    print(Z)
    print('\nmask')
    print(len(mask))
    print(len(mask[0]))
    print(len(mask[0][0]))
    print('\nimarray.shape')
    print(imarray.shape)

    ROI_mask[ROI_mask > 0] = 1
    mean_BVTV = 0.01

    # calculate center of sphere in new image
    xc = x - X[0]
    yc = y - Y[0]
    zc = z - Z[0]
    print('\nx: ' + str(x))
    print(', X: ' + str(X[0]))
    print(', xc: ' + str(xc))
    print('\ny: ' + str(y))
    print(', Y: ' + str(Y[0]))
    print(', yc: ' + str(yc))
    print('\nz: ' + str(z))
    print(', Z: ' + str(Z[0]))
    print(', zc: ' + str(zc))


    if phi > 0.0:
        # Compute BVTV
        # ------------------------
        # create masking array with the shape of the sphere

        #print('\nnumpy.shape(ROI), ROIsize, xc, yc, zc\n')

        #print(numpy.shape(ROI))
        #print(ROIsize)
        #print(xc)
        #print(yc)
        #print(zc)

        ROImask_sphere = sphere_array(numpy.shape(ROI), ROIsize / 2, [xc, yc, zc])

        # Overly the ROImask_sphere with the ROI_mask
        # All voxel containing a 2 are inside the sphere and will be used to calculate the BV/TV
        # take alle voxel entries inside the intersection between the two masks for BVTV calculation

        BVTVimg = ROI[ROImask_sphere + ROI_mask == 2]

        # The mean over whole sphere gives the BVTV
        mean_BVTV = numpy.round(numpy.mean(BVTVimg), 2)

        # check for meaningfull output
        if numpy.isnan(mean_BVTV):
            mean_BVTV = 0.01
        if mean_BVTV > 1:
            mean_BVTV = 1
        if mean_BVTV == 0.00:
            mean_BVTV = 0.01

    return mean_BVTV


def computePHI(cog, Spacing, ROIsize, imarray):
    """computes bone partial volume from a numpy array containing the MASK values for a region of size 'ROIsize' centered in the center of gravity of the element"""
    x, y, z = cog / Spacing
    ROIsize = ROIsize / Spacing[0]
    X = [max(x - ROIsize / 2, 0), min(x + ROIsize / 2, imarray.shape[0])]
    Y = [max(y - ROIsize / 2, 0), min(y + ROIsize / 2, imarray.shape[1])]
    Z = [max(z - ROIsize / 2, 0), min(z + ROIsize / 2, imarray.shape[2])]

    ROI = imarray[
        int(numpy.rint(X[0])) : int(numpy.rint(X[1])),
        int(numpy.rint(Y[0])) : int(numpy.rint(Y[1])),
        int(numpy.rint(Z[0])) : int(numpy.rint(Z[1])),
    ]
    try:
        PHI = float(numpy.count_nonzero(ROI)) / ROI.size
        # print(PHI)
    except:
        PHI = 0

    # check for meaningfull output
    if numpy.isnan(PHI):
        PHI = 0.0
    if PHI > 1:
        PHI = 1
    return PHI, X, Y, Z


def writeoutImageElement(cog, Spacing, ROIsize, imarray, i, filename):
    x, y, z = cog / Spacing
    ROIsize = ROIsize / Spacing[0] * 3
    X = [max(x - ROIsize / 2, 0), min(x + ROIsize / 2, imarray.shape[0])]
    Y = [max(y - ROIsize / 2, 0), min(y + ROIsize / 2, imarray.shape[1])]
    Z = [max(z - ROIsize / 2, 0), min(z + ROIsize / 2, imarray.shape[2])]
    ROI = imarray[
        int(numpy.rint(X[0])) : int(numpy.rint(X[1])),
        int(numpy.rint(Y[0])) : int(numpy.rint(Y[1])),
        int(numpy.rint(Z[0])) : int(numpy.rint(Z[1])),
    ]
    filename2 = str(filename + "_" + str(i) + ".mhd")
    numpy2mhd(ROI, Spacing, filename2)
    # return ROI, filename2


# Function for computing global fabric (Denis)
# Computed fabric on a global scale (in contrast to ROI for local fabric evaluation)
# Spacing = real element size in mm, imarray = numpy array containing SEG values


def compute_globalFAB(Spacing, imarray):
    # computes fabric from a numpy array containing the SEG values
    # create STL from isosurface
    BV = float(numpy.sum(imarray))
    VTK = numpy2vtk(imarray, Spacing)
    STL = vtk.vtkDiscreteMarchingCubes()
    try:
        STL.SetInputData(VTK)
    except:
        STL.SetInput(VTK)
    STL.GenerateValues(1, 1, 1)
    STL.Update()

    # decimate STL
    STLdeci = vtk.vtkDecimatePro()
    STLdeci.SetInputConnection(STL.GetOutputPort())
    STLdeci.SetTargetReduction(0.9)
    STLdeci.PreserveTopologyOn()
    STLdeci.SplittingOn()
    STLdeci.BoundaryVertexDeletionOn()
    STLdeci.Update()
    if VERIF == "Yes" or VERIF == "yes":
        STLname = ext(BMDname, ".stl")
        writer = vtk.vtkSTLWriter()
        writer.SetInputConnection(STLdeci.GetOutputPort())
        writer.SetFileName(STLname)
        writer.Write()
    if debug:
        print("\n ... ... decimate STL completed\n")

    # compute MSL
    vtkSTL = STLdeci.GetOutput()
    nfacet = numpy.arange(vtkSTL.GetNumberOfCells())
    facet = {
        i: [vtkSTL.GetCell(i).GetPoints().GetPoint(j) for j in range(3)] for i in nfacet
    }

    vtkNOR = vtk.vtkPolyDataNormals()
    # vtkNOR.SetInputData(vtkSTL)
    # IMPORTANT NOTE:
    #  SetInputData did not work!
    # Probably because locally I'm running VTK 5. This needs to be checked when running on the servers
    try:
        vtkNOR.SetInputData(vtkSTL)
    except:
        vtkNOR.SetInput(vtkSTL)

    vtkNOR.ComputeCellNormalsOn()
    vtkNOR.ComputePointNormalsOff()
    vtkNOR.Update()
    ndata = vtkNOR.GetOutput().GetCellData().GetNormals()
    dyad = {
        i: numpy.outer(numpy.array(ndata.GetTuple(i)), numpy.array(ndata.GetTuple(i)))
        for i in nfacet
    }

    vtkARE = vtk.vtkTriangle()
    area = {
        i: vtkARE.TriangleArea(facet[i][0], facet[i][1], facet[i][2]) for i in nfacet
    }

    # Computation of MSL fabric according to Hosseini et al. 2017
    A = sum(area[i] * dyad[i] for i in nfacet)  # Eq. 1
    H = scipy.linalg.inv(A) * (2.0 * BV)  # Eq. 2
    MSL = 3.0 * H / numpy.trace(H)  # Eq. 3
    evalue, evect = scipy.linalg.eig(MSL)

    if debug == 1:
        print("\n... ... compute MSL completed\n")

    # order eigenvalues 0=max, 1=mid, 2=min
    idx = evalue.argsort()[::-1]
    evalue = evalue[idx]
    evect = evect[:, idx]
    evalue = [e.real for e in evalue]
    evect = [evect[:, i] for i in [0, 1, 2]]
    return evalue, evect


# Function for computing local fabric (Denis)
# In contrast to global fabric, for each element, the fabric is calculated in ROI only.
# cog = center of gravity, Spacing = real element size in mm, ROIsize = size of region of interesst
# imarray = numpy array containing SEG values
# def compute_localFAB_Hadi(Spacing, ROIsize, ROI, BV, name):
#     """computes fabric from a numpy array containing the SEG values for a region of size rad centered in the center of gravity of the element"""
#
#     threshold = 1.0  # gray value of the bone voxels in the segmented image.
#     reduction = 0.9  # percentage of triangle decimation between 0 (no decimation) and 1 (100% decimation)
#     topol = 1  # 0: does not preserve the actual topology; 1: preserves the actual topology
#     # NB: the feature angle is 15 degres by default in Paraview
#
#     # Compute dimensions of ROI in mm
#     # why /2?
#     dimX = ROI.shape[0] * Spacing[0] / 2
#     dimY = ROI.shape[1] * Spacing[1] / 2
#     dimZ = ROI.shape[2] * Spacing[2] / 2
#
#     diameter = ROIsize
#
#     image_name = name[:-4]
#     stl = image_name + '.stl'
#
#     if BV > 0.1:
#         try:
#             mhd = MetaFileSeriesReader(FileNames=[image_name + '.mhd'])  # Paraview 4.3
#         except NameError:
#             mhd = MetaImagereader(FileName=image_name + '.mhd')  # Paraview 3.8
#
#         # Create contour with specified parameters
#         Cont = Contour(PointMergeMethod="Uniform Binning")
#         # specifies an incremental point locator for merging dublicate/coincident points
#         Cont.PointMergeMethod = "Uniform Binning"
#         #  specifies the name of the scalar array from which the contour filter will compute isolines and/or isosurfaces
#         Cont.ContourBy = ['POINTS', 'MetaImage']
#         # specifies the values at which to compute isosurfaces/isolines and also the number of such values
#         Cont.Isosurfaces = [threshold]
#
#         # The Decimate filter reduces the number of triangles in a polygonal data set. Because this filter only operates on
#         # triangles, first run the Triangulate filter on a dataset that contains polygons other than triangles.
#         # https://www.paraview.org/ParaView/Doc/Nightly/www/py-doc/paraview.simple.Decimate.html
#         Decim = Decimate()
#
#         # desired reduction in the total number of polygons in the output dataset. For example, if the TargetReduction
#         # value is 0.9, the Decimate filter will attempt to produce an output dataset that is 10% the size of the input
#         Decim.TargetReduction = reduction
#
#         # f this property is set to 1, decimation will not split the dataset or produce holes, but it may keep the filter
#         # from reaching the reduction target. If it is set to 0, better reduction can occur (reaching the reduction target),
#         # but holes in the model may be produced.
#         Decim.PreserveTopology = topol
#
#         # writes stereo lithography (.stl) files in either ASCII or binary form. Stereo lithography files only
#         # contain triangles.
#         writer = PSTLWriter(FileName=image_name + '.stl', Input=Decim, FileType='Ascii')
#         writer.UpdatePipeline()
#
#         # print "Read stl mesh"
#         stlFile = open(stl, 'r')
#         linesstl = stlFile.readlines()
#         stlFile.close()
#
#         # Evaluate MSL
#         elemno, Adpsum = MSL_f.msl_function(BV, stl, dimX, dimY, dimZ, diameter)  # equation 1 in Hosseini et al. 2017
#         H = LA.inv(Adpsum)  # equation 2 in Hosseini et al. 2017 (2*BV was already done in msl_function)
#         MSL = 3.0 * H / numpy.trace(H)  # equation 3 in Hosseini et al. 2017
#
#         eigvals, eigvects = LA.eig(MSL)
#
#         eigval1 = numpy.round(numpy.real(eigvals[0]), 4)
#         eigval2 = numpy.round(numpy.real(eigvals[1]), 4)
#         eigval3 = numpy.round(numpy.real(eigvals[2]), 4)
#
#         eigvec11 = numpy.round(numpy.real(eigvects[0][0]), 4)
#         eigvec12 = numpy.round(numpy.real(eigvects[1][0]), 4)
#         eigvec13 = numpy.round(numpy.real(eigvects[2][0]), 4)
#
#         eigvec21 = numpy.round(numpy.real(eigvects[0][1]), 4)
#         eigvec22 = numpy.round(numpy.real(eigvects[1][1]), 4)
#         eigvec23 = numpy.round(numpy.real(eigvects[2][1]), 4)
#
#         eigvec31 = numpy.round(numpy.real(eigvects[0][2]), 4)
#         eigvec32 = numpy.round(numpy.real(eigvects[1][2]), 4)
#         eigvec33 = numpy.round(numpy.real(eigvects[2][2]), 4)
#
#         return [eigval1, eigval2, eigval3], [[eigvec11, eigvec12, eigvec13], [eigvec21, eigvec22, eigvec23],
#                                              [eigvec31, eigvec32, eigvec33]]
#
#         del stl
#
#         # print "hello"
#
#     else:
#         return [1.0, 1.0, 1.0], [[0.0, 0.0, 1.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0]]
#         # set dimensions
#         # # IMPORTANT: This value depends on the segmented image. We should calculate fabric only in trabecular structure.
#         # # If seg image has cort = 1; trab = 2; background = 0, the following value must be ROI>1 to delete everything
#         # # which is not trabecular bone.
#         # ROI[ROI > 1] = 1
#         #
#         # # create square mask qith shape 6mm x 6mm x 6mm
#         # r2 = numpy.arange(-ROIsize / 2, ROIsize / 2) ** 2
#         # dist2 = r2[:, None, None] + r2[:, None] + r2
#         # SPH = numpy.zeros_like(dist2)
#         # # all values within R^2 are set to 1 (sphere mask)
#         # SPH[dist2 <= (ROIsize / 2) ** 2] = 1
#         #
#         # # mask ROI
#         # # NOTE: This is depending on how the ROI should be set. If the ROI for fabric evaluation is only trabecular structure,
#         # # ROI[ROI == 1] = 1, and ROI[ROI == 2] = 0.
#         # # If the BV for fabric evaluation should include cortex as well: ROI[ROI == 1] = 1 and ROI[ROI == 2] = 1
#         # ROI = numpy.resize(ROI, SPH.shape)
#         # ROI = numpy.add(SPH, ROI)
#         # ROI[ROI == 1] = 1
#         # ROI[ROI == 2] = 1
#         # BVnumber = numpy.sum(ROI)
#         # print "\nBV = " + str(BV) + "\n"
#         #
#         # if BV > 0.1:
#         #     # create STL from isosurface
#         #     ROI = numpy2vtk(ROI, Spacing)
#         #     STL = vtk.vtkDiscreteMarchingCubes()
#         #     try:
#         #         STL.SetInputData(ROI)
#         #     except:
#         #         STL.SetInput(ROI)
#         #     STL.GenerateValues(1, 1, 1)
#         #     STL.Update()
#         #
#         #     # decimate STL
#         #     STLdeci = vtk.vtkDecimatePro()
#         #     STLdeci.SetInputConnection(STL.GetOutputPort())
#         #     STLdeci.SetTargetReduction(0.9)
#         #     STLdeci.PreserveTopologyOn()
#         #     STLdeci.Update()
#         #
#         #
#         #     writer = vtk.vtkSTLWriter()
#         #     writer.SetFileName('/home/schenk/Documents/PhD/06_hFE/11_Fabric/01_Benjamin_Code/ROI.stl')
#         #     writer.SetInputConnection(STLdeci.GetOutputPort())
#         #     writer.Write()
#         #
#         #
#         #     elemno, Adpsum = MSL_f.msl_function(BV, '/home/schenk/Documents/PhD/06_hFE/11_Fabric/01_Benjamin_Code/ROI.stl', dimX, dimY, dimZ, diameter)  # equation 1 in Hosseini et al. 2017
#         #     try:
#         #         H = scipy.linalg.inv(Adpsum) * (2.0 * BV)  # Eq. 2
#         #         MSL = 3.0 * H / numpy.trace(H)  # Eq. 3
#         #         evalue, evect = scipy.linalg.eig(MSL)
#         #
#         #         # order eigenvalues 0=max, 1=mid, 2=min
#         #         # idx = evalue.argsort()[::-1]
#         #         # evalue = evalue[idx]
#         #         # evect = evect[:, idx]
#         #         # evalue = [e.real for e in evalue]
#         #         # evect = [evect[:, i] for i in [0, 1, 2]]
#         #         return evalue, evect
#         #     except:
#         #         print "Expect function was called --> isotropic fabric!"
#         #         return [1.0, 1.0, 1.0], [[0.0, 0.0, 1.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0]]
#         # else:
#         #     return [1.0, 1.0, 1.0], [[0.0, 0.0, 1.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0]]


def compute_localFAB(Spacing, ROIsize, ROI):
    """computes fabric from a numpy array containing the SEG values for a region of size rad centered in the center of gravity of the element"""
    # set dimensions
    # IMPORTANT: This value depends on the segmented image. We should calculate fabric only in trabecular structure.
    # If seg image has cort = 1; trab = 2; background = 0, the following value must be ROI>1 to delete everything
    # which is not trabecular bone.
    ROI[ROI > 1] = 1

    # create square mask with shape 6mm x 6mm x 6mm
    r2 = numpy.arange(-ROIsize / 2, ROIsize / 2) ** 2
    dist2 = r2[:, None, None] + r2[:, None] + r2
    SPH = numpy.zeros_like(dist2)
    # all values within R^2 are set to 1 (sphere mask)
    SPH[dist2 <= (ROIsize / 2) ** 2] = 1

    # mask ROI
    # NOTE: This is depending on how the ROI should be set. If the ROI for fabric evaluation is only trabecular structure,
    # ROI[ROI == 1] = 1, and ROI[ROI == 2] = 0.
    # If the BV for fabric evaluation should include cortex as well: ROI[ROI == 1] = 1 and ROI[ROI == 2] = 1
    ROI = numpy.resize(ROI, SPH.shape)
    ROI = numpy.add(SPH, ROI)
    ROI[ROI == 1] = 1
    ROI[ROI == 2] = 1
    BV = numpy.sum(ROI)
    print("\nBV = " + str(BV) + "\n")

    if BV > 0.1:
        # create STL from isosurface
        ROI = numpy2vtk(ROI, Spacing)
        STL = vtk.vtkDiscreteMarchingCubes()
        try:
            STL.SetInputData(ROI)
        except:
            STL.SetInput(ROI)
        STL.GenerateValues(1, 1, 1)
        STL.Update()

        # decimate STL
        STLdeci = vtk.vtkDecimatePro()
        STLdeci.SetInputConnection(STL.GetOutputPort())
        STLdeci.SetTargetReduction(0.9)
        STLdeci.PreserveTopologyOn()
        STLdeci.Update()

        # compute MSL
        vtkSTL = STLdeci.GetOutput()
        nfacet = numpy.arange(vtkSTL.GetNumberOfCells())
        facet = {
            i: [vtkSTL.GetCell(i).GetPoints().GetPoint(j) for j in range(3)]
            for i in nfacet
        }

        vtkNOR = vtk.vtkPolyDataNormals()
        # vtkNOR.SetInputData(vtkSTL)
        # IMPORTANT NOTE:
        #  SetInputData did not work! Probably because locally I'm running VTK 5. This needs to be checked when running on the servers
        try:
            vtkNOR.SetInputData(vtkSTL)
        except:
            vtkNOR.SetInput(vtkSTL)
        vtkNOR.ComputeCellNormalsOn()
        vtkNOR.ComputePointNormalsOff()
        vtkNOR.Update()
        ndata = vtkNOR.GetOutput().GetCellData().GetNormals()
        dyad = {
            i: numpy.outer(
                numpy.array(ndata.GetTuple(i)), numpy.array(ndata.GetTuple(i))
            )
            for i in nfacet
        }

        vtkARE = vtk.vtkTriangle()
        area = {
            i: vtkARE.TriangleArea(facet[i][0], facet[i][1], facet[i][2])
            for i in nfacet
        }
        A = sum(area[i] * dyad[i] for i in nfacet)  # Eq. 1

        # # chech if square matrix
        # rows = len(A)
        # is_square = True
        # for row in A:
        #     if len(row) != rows:
        #         is_square = False
        # if is_square:
        #     H = scipy.linalg.inv(A) * (2.0 * BV)  # Eq. 2
        #     MSL = 3.0 * H / numpy.trace(H)  # Eq. 3
        #     evalue, evect = scipy.linalg.eig(MSL)
        #
        #     # order eigenvalues 0=max, 1=mid, 2=min
        #     idx = evalue.argsort()[::-1]
        #     evalue = evalue[idx]
        #     evect = evect[:, idx]
        #     evalue = [e.real for e in evalue]
        #     evect = [evect[:, i] for i in [0, 1, 2]]
        #     return evalue, evect
        # else:
        #     return [1.0, 1.0, 1.0], [[0.0, 0.0, 1.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0]]
        try:
            H = scipy.linalg.inv(A) * (2.0 * BV)  # Eq. 2
            MSL = 3.0 * H / numpy.trace(H)  # Eq. 3
            evalue, evect = scipy.linalg.eig(MSL)

            # order eigenvalues 0=max, 1=mid, 2=min
            idx = evalue.argsort()[::-1]
            evalue = evalue[idx]
            evect = evect[:, idx]
            evalue = [e.real for e in evalue]
            evect = [evect[:, i] for i in [0, 1, 2]]
            print("Fabric_worked")
            return evalue, evect
        except:
            print("Expect function was called --> isotropic fabric!")
            return [1.0, 1.0, 1.0], [[0.0, 0.0, 1.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0]]

    else:
        print("BV less than 0.1, isotropic fabric!")
        return [1.0, 1.0, 1.0], [[0.0, 0.0, 1.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0]]


def compute_isoFAB():
    """
    Returns isotropic fabric
    return [eval1, eval2, eval3], [evecxx, evecxy, evecxz], [evecyx, evecyy, evecyz], [eveczx, eveczy, eveczz]]
    """
    return [1.0, 1.0, 1.0], [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]

# write ensight functions for writing .case files
# def inp2case2phase(INPname, RHOc, RHOt, PHIc, PHIt):
#     """writes case files for checking the material mapping"""
#     # read abaqus input file
#     # myfec = fec.fec()
#     # inp = myfec.readAbaqus(INPname)
#     inp = mf.readAbaqus(INPname)
#     title = inp[0]
#     nodes = inp[1]
#     nsets = inp[2]
#     elems = inp[3]
#     elsets = inp[4]
#     # write case files
#     # myfec.writeEnsight(
#     mf.writeEnsight(
#         ext(INPname, "_BVTVcort.case"),
#         None,
#         nodes,
#         None,
#         elems,
#         elsets,
#         NscaResults=None,
#         EscaResults=[RHOc],
#         vecResults=None,
#         EvecResults=None,
#     )
#     # myfec.writeEnsight(
#     mf.writeEnsight(
#         ext(INPname, "_BVTVtrab.case"),
#         None,
#         nodes,
#         None,
#         elems,
#         elsets,
#         NscaResults=None,
#         EscaResults=[RHOt],
#         vecResults=None,
#         EvecResults=None,
#     )
#     # myfec.writeEnsight(
#     mf.writeEnsight(
#         ext(INPname, "_PHIcort.case"),
#         None,
#         nodes,
#         None,
#         elems,
#         elsets,
#         NscaResults=None,
#         EscaResults=[PHIc],
#         vecResults=None,
#         EvecResults=None,
#     )
#     # myfec.writeEnsight(
#     mf.writeEnsight(
#         ext(INPname, "_PHItrab.case"),
#         None,
#         nodes,
#         None,
#         elems,
#         elsets,
#         NscaResults=None,
#         EscaResults=[PHIt],
#         vecResults=None,
#         EvecResults=None,
#     )


# def inp2case1phase(INPname, RHOb, PHIb):
#     """writes case files for checking the material mapping"""
#     # read abaqus input file
#     # myfec = fec.fec()
#     # inp = myfec.readAbaqus(INPname)
#     inp = mf.readAbaqus(INPname)
#     title = inp[0]
#     nodes = inp[1]
#     nsets = inp[2]
#     elems = inp[3]
#     elsets = inp[4]
#     # write case files
#     # myfec.writeEnsight(
#     mf.writeEnsight(
#         ext(INPname, "_BVTVbone.case"),
#         None,
#         nodes,
#         None,
#         elems,
#         elsets,
#         NscaResults=None,
#         EscaResults=[RHOb],
#         vecResults=None,
#         EvecResults=None,
#     )
#     # myfec.writeEnsight(
#     mf.writeEnsight(
#         ext(INPname, "_PHIbone.case"),
#         None,
#         nodes,
#         None,
#         elems,
#         elsets,
#         NscaResults=None,
#         EscaResults=[PHIb],
#         vecResults=None,
#         EvecResults=None,
#     )


def fab2vtk(INPname, m, mm):
    """writes vtk files displaying the eigenvectors computed for a mesh (only for elements with fabric)"""
    # read abaqus input file
    # myfec = fec.fec()
    # inp = myfec.readAbaqus(INPname)
    inp = mf.readAbaqus(INPname)
    title = inp[0]
    nodes = inp[1]
    nsets = inp[2]
    elems = inp[3]
    elsets = inp[4]
    # calculate center of gravity of each element with fabric
    cog = {
        elem: numpy.mean(
            [
                numpy.asarray(nodes[node].get_coord())
                for node in elems[elem].get_nodes()
            ],
            axis=0,
        )
        for elem in m.keys()
    }
    # write vtk files for each eigenvector
    for i in [0, 1, 2]:
        if i == 0:
            vtkname = ext(INPname, "_FABmax.vtk")
        elif i == 1:
            vtkname = ext(INPname, "_FABmid.vtk")
        elif i == 2:
            vtkname = ext(INPname, "_FABmin.vtk")
        print(" ... write vtk file: " + vtkname)
        with open(vtkname, "w") as vtkFile:
            vtkFile.write("# vtk DataFile Version 2.0\n")
            vtkFile.write("Reconstructed Lagrangian Field Data\n")
            vtkFile.write("ASCII\n")
            vtkFile.write("DATASET UNSTRUCTURED_GRID\n")
            vtkFile.write("\nPOINTS " + str(2 * len(m.keys())) + " float\n")
            for elem in m.keys():
                vtkFile.write(
                    str(cog[elem][0] - m[elem][i] * mm[elem][i][0])
                    + " "
                    + str(cog[elem][1] - m[elem][i] * mm[elem][i][1])
                    + " "
                    + str(cog[elem][2] - m[elem][i] * mm[elem][i][2])
                    + "\n"
                )
            for elem in m.keys():
                vtkFile.write(
                    str(cog[elem][0] + m[elem][i] * mm[elem][i][0])
                    + " "
                    + str(cog[elem][1] + m[elem][i] * mm[elem][i][1])
                    + " "
                    + str(cog[elem][2] + m[elem][i] * mm[elem][i][2])
                    + "\n"
                )
            vtkFile.write(
                "\nCELLS " + str(len(m.keys())) + " " + str(3 * len(m.keys())) + "\n"
            )
            count = -1
            for elem in m.keys():
                count = count + 1
                vtkFile.write(
                    "2 " + str(count) + " " + str(count + len(m.keys())) + "\n"
                )
            vtkFile.write("\nCELL_TYPES " + str(len(m.keys())) + "\n")
            for elem in m.keys():
                vtkFile.write("3\n")
            if i == 0:
                vtkFile.write("\nCELL_DATA " + str(len(m.keys())) + "\n")
                vtkFile.write("scalars DOA_max float\n")
                vtkFile.write("LOOKUP_TABLE default\n")
                for elem in m.keys():
                    vtkFile.write(str(m[elem][0] / m[elem][2]) + "\n")


def assign_MSL(SEGim_vtk):
    from vtk.numpy_interface import dataset_adapter as dsa

    # Create STL file from segmented image
    STL = vtk.vtkDiscreteMarchingCubes()
    STL.SetInputData(SEGim_vtk)  # TODO try and except
    STL.GenerateValues(1, 1, 1)
    STL.Update()
    print("STL file creation finished")

    # decimate STL
    STLdeci = vtk.vtkDecimatePro()
    STLdeci.SetInputConnection(STL.GetOutputPort())
    STLdeci.SetTargetReduction(0.9)
    STLdeci.PreserveTopologyOn()
    STLdeci.Update()
    print("Decimation finished")

    # Calculate number of cells in triangulated mesh
    vtkSTL = STLdeci.GetOutput()
    nfacet = numpy.arange(vtkSTL.GetNumberOfCells())
    print("Number of cells calculated")

    triangle_points = []
    cog_points = []
    # calculate center of gravity for each triangle (xc = (x1+x2+x3)/3...)
    # only keep cogs which are not at the border (z-direction)
    indizes = []
    a = dsa.WrapDataObject(vtkSTL)
    for i in nfacet:
        triangle_points.append(
            [
                a.GetCell(i).GetPoints().GetPoint(0),
                a.GetCell(i).GetPoints().GetPoint(1),
                a.GetCell(i).GetPoints().GetPoint(2),
            ]
        )
        cog_points.append(
            [
                (
                    triangle_points[i][0][0]
                    + triangle_points[i][1][0]
                    + triangle_points[i][2][0]
                )
                / 3,
                (
                    triangle_points[i][0][1]
                    + triangle_points[i][1][1]
                    + triangle_points[i][2][1]
                )
                / 3,
                (
                    triangle_points[i][0][2]
                    + triangle_points[i][1][2]
                    + triangle_points[i][2][2]
                )
                / 3,
            ]
        )
    print("Computation COG finished")

    # calc cell normals and dyadic product
    vtkNormals = vtk.vtkPolyDataNormals()
    vtkNormals.SetInputConnection(STLdeci.GetOutputPort())
    vtkNormals.ComputeCellNormalsOn()
    vtkNormals.ComputePointNormalsOff()
    vtkNormals.ConsistencyOn()
    vtkNormals.AutoOrientNormalsOn()  # Only works with closed surface. All Normals will point outward.
    # vtkNormals.SetNonManifoldTraversal()  # TODO is this function helpful?
    vtkNormals.Update()
    PointNormalArray = vtkNormals.GetOutput().GetCellData().GetNormals()

    dyad = numpy.zeros([3, 3])
    dyadic = []
    for i in nfacet:
        vtk.vtkMath.Outer(
            PointNormalArray.GetTuple(i), PointNormalArray.GetTuple(i), dyad
        )
        dyadic.append(dyad.tolist())
    print("Computation dyadic products finished")

    # Get Cell Area https://www.vtk.org/Wiki/VTK/Examples/Python/MeshLabelImage
    triangleCellAN = vtk.vtkMeshQuality()
    triangleCellAN.SetInputConnection(vtkNormals.GetOutputPort())
    triangleCellAN.SetTriangleQualityMeasureToArea()
    triangleCellAN.SaveCellQualityOn()  # default
    triangleCellAN.Update()  # creates vtkDataSet
    qualityArray = triangleCellAN.GetOutput().GetCellData().GetArray("Quality")
    np_cog_points = numpy.array(cog_points)
    minZ = np_cog_points[:, 2].min()
    maxZ = np_cog_points[:, 2].max()

    # if cog is at zmax or zmin, area is set to zero
    # so that triangles on the cut surfaces are ignored for fabric
    tolerance = 0.5
    area = []
    for i in nfacet:
        if np_cog_points[i, 2] <= tolerance or np_cog_points[i, 2] >= maxZ - tolerance:
            area.append(0.0)
        else:
            area.append(qualityArray.GetValue(i))
    print("Computation areas finished")
    #
    # for i in nfacet:
    #     area.append(qualityArray.GetValue(i))
    # print "Computation areas finished"

    # area dyadic represents the multiplication of the area with the cross-product of the normals of each triangle
    # these values can now be assigned to the elements according to the center of gravity of the triangle
    # all lists are sorted according to index (position in list is index identifier for triangles)
    areadyadic = [numpy.multiply(a, b) for a, b in zip(area, dyadic)]

    # filename = '/home/schenk/Documents/PhD/06_hFE/08_Pipeline_homogenization/01_AIM2FE_homogenized/1_AIM/C00000X_Isosurface.stl'
    #
    # stlWriter = vtk.vtkSTLWriter()
    # stlWriter.SetFileName(filename)
    # stlWriter.SetInputConnection(STLdeci.GetOutputPort())
    # stlWriter.Write()

    return STLdeci, cog_points, areadyadic, nfacet, area


def assign_MSL_v2(SEGim_vtk, image_dimensions, tolerance):
    from vtk.numpy_interface import dataset_adapter as dsa

    dimX = image_dimensions[0]
    dimY = image_dimensions[1]
    dimZ = image_dimensions[2]

    # Create STL file from segmented image
    STL = vtk.vtkDiscreteMarchingCubes()
    STL.SetInputData(SEGim_vtk)  # TODO try and except
    STL.GenerateValues(1, 1, 1)
    STL.Update()
    print("1/7 STL file creation finished")

    # decimate STL
    STLdeci = vtk.vtkDecimatePro()
    STLdeci.SetInputConnection(STL.GetOutputPort())
    STLdeci.SetTargetReduction(0.9)
    STLdeci.PreserveTopologyOn()
    STLdeci.Update()
    print("2/7 Decimation finished")

    # Calculate number of cells in triangulated mesh
    vtkSTL = STLdeci.GetOutput()
    nfacet = numpy.arange(vtkSTL.GetNumberOfCells())
    print("3/7 Number of cells calculated")

    triangle_points = []
    cog_points = []
    # calculate center of gravity for each triangle (xc = (x1+x2+x3)/3...)
    # only keep cogs which are not at the border (z-direction)
    indizes = []
    a = dsa.WrapDataObject(vtkSTL)
    for i in nfacet:
        triangle_points.append(
            [
                a.GetCell(i).GetPoints().GetPoint(0),
                a.GetCell(i).GetPoints().GetPoint(1),
                a.GetCell(i).GetPoints().GetPoint(2),
            ]
        )
        cog_points_temp = [
            (
                triangle_points[i][0][0]
                + triangle_points[i][1][0]
                + triangle_points[i][2][0]
            )
            / 3,
            (
                triangle_points[i][0][1]
                + triangle_points[i][1][1]
                + triangle_points[i][2][1]
            )
            / 3,
            (
                triangle_points[i][0][2]
                + triangle_points[i][1][2]
                + triangle_points[i][2][2]
            )
            / 3,
        ]
        if 0 + tolerance <= cog_points_temp[2] <= dimZ - tolerance:
            cog_points.append(cog_points_temp)
            indizes.append(i)
    print("4/7 Computation COG finished")

    # calc cell normals and dyadic product
    vtkNormals = vtk.vtkPolyDataNormals()
    vtkNormals.SetInputConnection(STLdeci.GetOutputPort())
    vtkNormals.ComputeCellNormalsOn()
    vtkNormals.ComputePointNormalsOff()
    vtkNormals.ConsistencyOn()
    vtkNormals.AutoOrientNormalsOn()  # Only works with closed surface. All Normals will point outward.
    # vtkNormals.SetNonManifoldTraversal()  # TODO is this function helpful?
    vtkNormals.Update()
    PointNormalArray = vtkNormals.GetOutput().GetCellData().GetNormals()

    dyad = numpy.zeros([3, 3])
    dyadic = []
    for i in indizes:
        vtk.vtkMath.Outer(
            PointNormalArray.GetTuple(i), PointNormalArray.GetTuple(i), dyad
        )
        dyadic.append(dyad.tolist())
    print("5/7 Computation dyadic products finished")

    # Get Cell Area https://www.vtk.org/Wiki/VTK/Examples/Python/MeshLabelImage
    triangleCellAN = vtk.vtkMeshQuality()
    triangleCellAN.SetInputConnection(vtkNormals.GetOutputPort())
    triangleCellAN.SetTriangleQualityMeasureToArea()
    triangleCellAN.SaveCellQualityOn()  # default
    triangleCellAN.Update()  # creates vtkDataSet
    qualityArray = triangleCellAN.GetOutput().GetCellData().GetArray("Quality")

    # if cog is at zmax or zmin, area is set to zero
    # so that triangles on the cut surfaces are ignored for fabric
    area = []
    for i in indizes:
        area.append(qualityArray.GetValue(i))
    print("6/7 Computation areas finished")
    #
    # for i in nfacet:
    #     area.append(qualityArray.GetValue(i))
    # print "Computation areas finished"

    # area dyadic represents the multiplication of the area with the cross-product of the normals of each triangle
    # these values can now be assigned to the elements according to the center of gravity of the triangle
    # all lists are sorted according to index (position in list is index identifier for triangles)
    areadyadic = [numpy.multiply(a, b) for a, b in zip(area, dyadic)]

    # filename = '/home/schenk/Documents/PhD/06_hFE/08_Pipeline_homogenization/01_AIM2FE_homogenized/1_AIM/C00000X_Isosurface.stl'
    #
    # stlWriter = vtk.vtkSTLWriter()
    # stlWriter.SetFileName(filename)
    # stlWriter.SetInputConnection(STLdeci.GetOutputPort())
    # stlWriter.Write()
    print("7/7 Computation areas finished")

    return STLdeci, cog_points, areadyadic, nfacet, indizes, area


def factorize(n):
    l = []
    for i in chain([2], range(3, n // 2 + 1, 2)):
        while n % i == 0:
            l.append(i)
            n = n // i
        if i > n:
            break
    return l


def closest_divisor(numerator, divisor):
    assert isinstance(numerator, int)
    assert isinstance(divisor, int)
    assert numerator > divisor

    if numerator % divisor == 0:
        return divisor

    distance = 1

    while distance < divisor:
        if numerator % (divisor - distance) == 0:
            return divisor - distance
        elif numerator % (divisor + distance) == 0:
            return divisor + distance
        else:
            distance += 1

    raise ArithmeticError("Logical error minimal divisor should at least be 1")


def cap_permutations(s):
    lu_sequence = ((c.lower(), c.upper()) for c in s)
    return ["".join(x) for x in it.product(*lu_sequence)]


def find_simulation_type(ftype, phase):
    str_iso = str(cap_permutations("iso"))
    str_global = str(cap_permutations("global"))
    str_test = str(cap_permutations("test"))
    str_local = str(cap_permutations("local"))
    str_one = str(cap_permutations("one"))
    str_two = str(cap_permutations("two"))

    # test for fabric evaluation
    if str_iso.find(ftype) > -1 or str_test.find(ftype) > -1:
        fabric_eval = 0
        # print "found iso/test"
    elif str_local.find(ftype) > -1 or str_global.find(ftype) > -1:
        fabric_eval = 1
        # print "found local/global"
    else:
        raise NameError("Could not interpret fabric evaluation type")

    # test for number of phases
    if str_one.find(phase) > -1:
        phase_eval = 1
        # print "found 1 phase"
    elif str_two.find(phase) > -1:
        phase_eval = 2
        # print "found 2 phase"
    else:
        raise NameError("Could not interpret phase evaluation type")

    try:
        return fabric_eval, phase_eval
    except:
        return -1, -1


def vectoronplane(evect_max, evect_mid, evect_min, direction):
    # evect_max_projected is computed by projection of direction (usually [0,0,1]) into the plane evect_max, evect_mid.
    # evect_min_projected is orthogonal to max and mid, so computed from their cross product
    # evect_mid_projected is orthogonal to max and min, so computed from their cross product
    normal = numpy.cross(evect_max, evect_mid)
    normalized = normal / numpy.linalg.norm(normal)
    evect_max_projected_normalized = (
        direction - numpy.dot(direction, normalized) * normalized
    )
    evect_max_projected = evect_max_projected_normalized * numpy.linalg.norm(evect_max)

    evect_min_projected = normalized
    evect_mid_projected = numpy.cross(evect_max_projected, normalized)

    return evect_max_projected, evect_mid_projected, evect_min_projected


def checkkernel(number):
    # even number
    if number % 2 == 0:
        new_number = number + 1
        print(
            "Warning: Kernel size can't be even! Kernel changed to next higher integer: "
            + str(new_number)
        )
    # uneven number
    elif number % 2 == 1:
        new_number = number
    # zero
    elif number == 0:
        new_number = 1
        print("Warning: Kernel size can not be equal 0. Kernel size was changed to 1.")
    # no integer
    elif number % int(number) > 0:
        if int(number) % 2 == 0:
            new_number = int(number) + 1
            print(
                "Warning: Kernel size must be integer! Kernel changed to next higher uneven integer: "
                + str(new_number)
            )
        elif int(number) % 2 == 1:
            new_number = int(number)
            print(
                "Warning: Kernel size must be integer! Kernel changed to next higher uneven integer: "
                + str(new_number)
            )

    return new_number


def cort_ssss(mat_param_cort, phic, rhoc, mm3, mm2, mm1):

    if phic > 0:

        # rhoc = 0.8

        # Cortical material
        CCCC = numpy.zeros((6, 6))

        # Compliance tensor according to UMAT_Debug.f (TIRBCT, JJS, 2013)
        CCCC[0, 0] = 1.0 / mat_param_cort["E0"] / rhoc ** mat_param_cort["KS"]
        CCCC[1, 1] = CCCC[0, 0]
        CCCC[2, 2] = 1.0 / mat_param_cort["EAA"] / rhoc ** mat_param_cort["KS"]
        CCCC[3, 3] = (
            0.5
            / (mat_param_cort["E0"] / 2.0 / (1.0 + mat_param_cort["V0"]))
            / rhoc ** mat_param_cort["KS"]
        )
        CCCC[4, 4] = 0.5 / mat_param_cort["MUA0"] / rhoc ** mat_param_cort["KS"]
        CCCC[5, 5] = CCCC[4, 4]

        CCCC[1, 0] = (
            -1.0
            * mat_param_cort["V0"]
            / mat_param_cort["E0"]
            / rhoc ** mat_param_cort["KS"]
        )
        CCCC[2, 0] = (
            -1.0
            * mat_param_cort["VA0"]
            / mat_param_cort["EAA"]
            / rhoc ** mat_param_cort["KS"]
        )
        CCCC[2, 1] = CCCC[2, 0]

        CCCC[0, 1] = CCCC[1, 0]
        CCCC[0, 2] = CCCC[2, 0]
        CCCC[1, 2] = CCCC[2, 1]

        SSSS = numpy.linalg.inv(CCCC)

        # Transversly isotropic quadratic fourth order tensor FFFF
        S0 = (
            (mat_param_cort["SIGD0P"] + mat_param_cort["SIGD0N"])
            / 2.0
            / mat_param_cort["SIGD0P"]
            / mat_param_cort["SIGD0N"]
        )
        SA = (
            (mat_param_cort["SIGDAP"] + mat_param_cort["SIGDAN"])
            / 2.0
            / mat_param_cort["SIGDAP"]
            / mat_param_cort["SIGDAN"]
        )
        TAUD0 = numpy.sqrt(0.5 / S0 ** 2.0 / (1.0 + mat_param_cort["ZETA0"]))

        FFFF = numpy.zeros((6, 6))
        FFFF[0, 0] = S0 ** 2.0  # FFFF 11
        FFFF[1, 1] = S0 ** 2.0  # FFFF 22
        FFFF[2, 2] = SA ** 2.0  # FFFF 33
        FFFF[3, 3] = 0.5 / TAUD0 ** 2.0  # FFFF 44
        FFFF[4, 4] = (
            0.5 / mat_param_cort["TAUDA0"] ** 2.0
        )  # FFFF 55 #TODO not the same in UMAT and mathematica file
        FFFF[5, 5] = (
            0.5 / mat_param_cort["TAUDA0"] ** 2.0
        )  # FFFF 66 #TODO not the same in UMAT and mathematica file

        FFFF[0, 1] = (-1.0 * mat_param_cort["ZETA0"]) * S0 ** 2.0
        FFFF[0, 2] = (-1.0 * mat_param_cort["ZETAA0"]) * SA ** 2.0

        FFFF[1, 2] = FFFF[0, 2]
        FFFF[1, 0] = FFFF[0, 1]
        FFFF[2, 0] = FFFF[0, 2]
        FFFF[2, 1] = FFFF[1, 2]

        FFFF = FFFF / rhoc ** mat_param_cort["PP"] / rhoc ** mat_param_cort["PP"]

        # Transversly isotropic quadratic second order tensor FF
        FF = numpy.zeros((3, 3))
        FF[0, 0] = (
            (-mat_param_cort["SIGD0P"] + mat_param_cort["SIGD0N"])
            / 2.0
            / mat_param_cort["SIGD0P"]
            / mat_param_cort["SIGD0N"]
        )
        FF[1, 1] = FF[0, 0]
        FF[2, 2] = (
            (-mat_param_cort["SIGDAP"] + mat_param_cort["SIGDAN"])
            / 2.0
            / mat_param_cort["SIGDAP"]
            / mat_param_cort["SIGDAN"]
        )
        FF = FF / rhoc ** mat_param_cort["PP"]

        return SSSS, CCCC, FFFF, FF
    else:
        empty = numpy.zeros(
            (6, 6)
        )  # TODO maybe this needs to be changed, that phic == 0 doesn't come to the function at all
        return empty, empty, empty, empty


def trab_ssss(mat_param_trab, phit, rhot, mm3, mm2, mm1):
    if phit > 0:
        #
        # rhot = 0.25
        # mm1 = 1.0
        # mm2 = 1.0
        # mm3 = 1.0

        lambda_plus_two_mu = (mat_param_trab["E0"] * (1.0 - mat_param_trab["V0"])) / (
            (1.0 + mat_param_trab["V0"]) * (1.0 - 2.0 * mat_param_trab["V0"])
        )
        lambda_zero = (
            mat_param_trab["E0"]
            * mat_param_trab["V0"]
            / ((1.0 + mat_param_trab["V0"]) * (1.0 - 2.0 * mat_param_trab["V0"]))
        )
        SSSS = numpy.zeros((6, 6))

        # Stiffness tensor acc. to UMAT_Debug.f (FABRGWTB PMUBC, PHZ, 2013 / FABRGWTB, PHZ, 2013)
        SSSS[0, 0] = (
            mm1 ** (2.0 * mat_param_trab["LS"])
            * rhot ** mat_param_trab["KS"]
            * lambda_plus_two_mu
        )
        SSSS[1, 1] = (
            mm2 ** (2.0 * mat_param_trab["LS"])
            * rhot ** mat_param_trab["KS"]
            * lambda_plus_two_mu
        )
        SSSS[2, 2] = (
            mm3 ** (2.0 * mat_param_trab["LS"])
            * rhot ** mat_param_trab["KS"]
            * lambda_plus_two_mu
        )
        SSSS[3, 3] = (
            2.0
            * mat_param_trab["MU0"]
            * (mm1 ** mat_param_trab["LS"])
            * (mm2 ** mat_param_trab["LS"])
            * (rhot ** mat_param_trab["KS"])
        )
        SSSS[4, 4] = (
            2.0
            * mat_param_trab["MU0"]
            * (mm3 ** mat_param_trab["LS"])
            * (mm1 ** mat_param_trab["LS"])
            * (rhot ** mat_param_trab["KS"])
        )
        SSSS[5, 5] = (
            2.0
            * mat_param_trab["MU0"]
            * (mm2 ** mat_param_trab["LS"])
            * (mm3 ** mat_param_trab["LS"])
            * (rhot ** mat_param_trab["KS"])
        )

        SSSS[1, 0] = (
            lambda_zero
            * (mm1 ** mat_param_trab["LS"])
            * (mm2 ** mat_param_trab["LS"])
            * (rhot ** mat_param_trab["KS"])
        )
        SSSS[2, 0] = (
            lambda_zero
            * (mm3 ** mat_param_trab["LS"])
            * (mm1 ** mat_param_trab["LS"])
            * (rhot ** mat_param_trab["KS"])
        )
        SSSS[2, 1] = (
            lambda_zero
            * (mm2 ** mat_param_trab["LS"])
            * (mm3 ** mat_param_trab["LS"])
            * (rhot ** mat_param_trab["KS"])
        )

        SSSS[0, 1] = SSSS[1, 0]
        SSSS[0, 2] = SSSS[2, 0]
        SSSS[1, 2] = SSSS[2, 1]

        CCCC = numpy.linalg.inv(SSSS)

        # Quadratic fourth order tensor FFFF
        FFFF = numpy.zeros((6, 6))
        S0 = (
            (mat_param_trab["SIGD0P"] + mat_param_trab["SIGD0N"])
            / 2.0
            / mat_param_trab["SIGD0P"]
            / mat_param_trab["SIGD0N"]
        )
        FFFF[0, 0] = (
            S0 ** 2
            / ((rhot ** mat_param_trab["PP"]) * mm1 ** (2.0 * mat_param_trab["QQ"]))
            ** 2
        )
        FFFF[1, 1] = (
            S0 ** 2
            / ((rhot ** mat_param_trab["PP"]) * mm2 ** (2.0 * mat_param_trab["QQ"]))
            ** 2
        )
        FFFF[2, 2] = (
            S0 ** 2
            / ((rhot ** mat_param_trab["PP"]) * mm3 ** (2.0 * mat_param_trab["QQ"]))
            ** 2
        )
        FFFF[3, 3] = 0.5 / (
            (
                mat_param_trab["TAUD0"]
                * (rhot ** mat_param_trab["PP"])
                * (mm1 * mm2) ** mat_param_trab["QQ"]
            )
            ** 2
        )
        FFFF[4, 4] = 0.5 / (
            (
                mat_param_trab["TAUD0"]
                * (rhot ** mat_param_trab["PP"])
                * (mm1 * mm3) ** mat_param_trab["QQ"]
            )
            ** 2
        )
        FFFF[5, 5] = 0.5 / (
            (
                mat_param_trab["TAUD0"]
                * (rhot ** mat_param_trab["PP"])
                * (mm2 * mm3) ** mat_param_trab["QQ"]
            )
            ** 2
        )

        FFFF[0, 1] = (
            -(mat_param_trab["ZETA0"] * (mm1 / mm2) ** (2.0 * mat_param_trab["QQ"]))
            * S0 ** 2
            / ((rhot ** mat_param_trab["PP"]) * mm1 ** (2.0 * mat_param_trab["QQ"]))
            ** 2
        )
        FFFF[0, 2] = (
            -(mat_param_trab["ZETA0"] * (mm1 / mm3) ** (2.0 * mat_param_trab["QQ"]))
            * S0 ** 2
            / ((rhot ** mat_param_trab["PP"]) * mm1 ** (2.0 * mat_param_trab["QQ"]))
            ** 2
        )
        FFFF[1, 2] = (
            -(mat_param_trab["ZETA0"] * (mm2 / mm3) ** (2.0 * mat_param_trab["QQ"]))
            * S0 ** 2
            / ((rhot ** mat_param_trab["PP"]) * mm2 ** (2.0 * mat_param_trab["QQ"]))
            ** 2
        )

        FFFF[1, 0] = FFFF[0, 1]
        FFFF[2, 0] = FFFF[0, 2]
        FFFF[2, 1] = FFFF[1, 2]

        # Quadratic second order tensor FF
        FF = numpy.zeros((3, 3))
        FF[0, 0] = (
            (-mat_param_trab["SIGD0P"] + mat_param_trab["SIGD0N"])
            / 2.0
            / mat_param_trab["SIGD0P"]
            / mat_param_trab["SIGD0N"]
            / ((rhot ** mat_param_trab["PP"]) * mm1 ** (2.0 * mat_param_trab["QQ"]))
        )
        FF[1, 1] = (
            (-mat_param_trab["SIGD0P"] + mat_param_trab["SIGD0N"])
            / 2.0
            / mat_param_trab["SIGD0P"]
            / mat_param_trab["SIGD0N"]
            / ((rhot ** mat_param_trab["PP"]) * mm2 ** (2.0 * mat_param_trab["QQ"]))
        )
        FF[2, 2] = (
            (-mat_param_trab["SIGD0P"] + mat_param_trab["SIGD0N"])
            / 2.0
            / mat_param_trab["SIGD0P"]
            / mat_param_trab["SIGD0N"]
            / ((rhot ** mat_param_trab["PP"]) * mm3 ** (2.0 * mat_param_trab["QQ"]))
        )

        # Quadratic fourth order tensor FFFF
        FFFF = numpy.zeros((6, 6))
        S0 = (
            (mat_param_trab["SIGD0P"] + mat_param_trab["SIGD0N"])
            / 2.0
            / mat_param_trab["SIGD0P"]
            / mat_param_trab["SIGD0N"]
        )
        FFFF[0, 0] = (
            S0 ** 2
            / ((rhot ** mat_param_trab["PP"]) * mm1 ** (2.0 * mat_param_trab["QQ"]))
            ** 2
        )
        FFFF[1, 1] = (
            S0 ** 2
            / ((rhot ** mat_param_trab["PP"]) * mm2 ** (2.0 * mat_param_trab["QQ"]))
            ** 2
        )
        FFFF[2, 2] = (
            S0 ** 2
            / ((rhot ** mat_param_trab["PP"]) * mm3 ** (2.0 * mat_param_trab["QQ"]))
            ** 2
        )
        FFFF[3, 3] = 0.5 / (
            (
                mat_param_trab["TAUD0"]
                * (rhot ** mat_param_trab["PP"])
                * (mm1 * mm2) ** mat_param_trab["QQ"]
            )
            ** 2
        )
        FFFF[4, 4] = 0.5 / (
            (
                mat_param_trab["TAUD0"]
                * (rhot ** mat_param_trab["PP"])
                * (mm1 * mm3) ** mat_param_trab["QQ"]
            )
            ** 2
        )
        FFFF[5, 5] = 0.5 / (
            (
                mat_param_trab["TAUD0"]
                * (rhot ** mat_param_trab["PP"])
                * (mm2 * mm3) ** mat_param_trab["QQ"]
            )
            ** 2
        )

        FFFF[0, 1] = (
            -(mat_param_trab["ZETA0"] * (mm1 / mm2) ** (2.0 * mat_param_trab["QQ"]))
            * S0 ** 2
            / ((rhot ** mat_param_trab["PP"]) * mm1 ** (2.0 * mat_param_trab["QQ"]))
            ** 2
        )
        FFFF[0, 2] = (
            -(mat_param_trab["ZETA0"] * (mm1 / mm3) ** (2.0 * mat_param_trab["QQ"]))
            * S0 ** 2
            / ((rhot ** mat_param_trab["PP"]) * mm1 ** (2.0 * mat_param_trab["QQ"]))
            ** 2
        )
        FFFF[1, 2] = (
            -(mat_param_trab["ZETA0"] * (mm2 / mm3) ** (2.0 * mat_param_trab["QQ"]))
            * S0 ** 2
            / ((rhot ** mat_param_trab["PP"]) * mm2 ** (2.0 * mat_param_trab["QQ"]))
            ** 2
        )

        FFFF[1, 0] = FFFF[0, 1]
        FFFF[2, 0] = FFFF[0, 2]
        FFFF[2, 1] = FFFF[1, 2]

        # Quadratic second order tensor FF
        FF = numpy.zeros((3, 3))
        FF[0, 0] = (
            (-mat_param_trab["SIGD0P"] + mat_param_trab["SIGD0N"])
            / 2.0
            / mat_param_trab["SIGD0P"]
            / mat_param_trab["SIGD0N"]
            / ((rhot ** mat_param_trab["PP"]) * mm1 ** (2.0 * mat_param_trab["QQ"]))
        )
        FF[1, 1] = (
            (-mat_param_trab["SIGD0P"] + mat_param_trab["SIGD0N"])
            / 2.0
            / mat_param_trab["SIGD0P"]
            / mat_param_trab["SIGD0N"]
            / ((rhot ** mat_param_trab["PP"]) * mm2 ** (2.0 * mat_param_trab["QQ"]))
        )
        FF[2, 2] = (
            (-mat_param_trab["SIGD0P"] + mat_param_trab["SIGD0N"])
            / 2.0
            / mat_param_trab["SIGD0P"]
            / mat_param_trab["SIGD0N"]
            / ((rhot ** mat_param_trab["PP"]) * mm3 ** (2.0 * mat_param_trab["QQ"]))
        )

        return SSSS, CCCC, FFFF, FF

    else:
        empty = numpy.zeros((6, 6))
        return empty, empty, empty, empty


def material_superposition(SSSSc, SSSSt, FFFFc, FFFFt, FFc, FFt, phic, phit):
    # Material superposition stiffness
    if phit > 0 and phic > 0:
        try:
            SSSSct = phic * SSSSc + phit * SSSSt
        except:
            print(
                "WARNING: Material superposition did not work! Stiffness tensor of trabecular bone is used."
            )
            SSSSct = SSSSt
        try:
            CCCCct = numpy.linalg.inv(SSSSct)
        except:
            print(
                "WARNING: Material superposition did not work! Compliance tensor of trabecular bone is used."
            )
            CCCCct = numpy.linalg.inv(SSSSt)

        # Material superposition FFFF
        try:
            FFFFct = numpy.linalg.inv(
                phic * numpy.linalg.inv(FFFFc) + phit * numpy.linalg.inv(FFFFt)
            )
        except:
            try:
                print(
                    "WARNING: Material superposition did not work! FFFF tensor of trabecular bone is used."
                )
                FFFFct = FFFFt
            except:
                print(
                    "WARNING: Material superposition did not work! FFFF tensor of compact bone is used."
                )
                FFFFct = FFFFc

        try:
            FFct = numpy.linalg.inv(
                phic * numpy.linalg.inv(FFc) + phit * numpy.linalg.inv(FFt)
            )
        except:
            try:
                print(
                    "WARNING: Material superposition did not work! FF tensor of trabecular bone is used."
                )
                FFct = FFt
            except:
                print(
                    "WARNING: Material superposition did not work! FF tensor of compact bone is used."
                )
                FFct = FFc

    elif phit == 0:
        SSSSct = SSSSt
        CCCCct = CCCCt
        FFFFct = FFFFt
        FFct = FFt

    elif phic == 0:
        SSSSct = SSSSc
        CCCCct = CCCCc
        FFFFct = FFFFc
        FFct = FFc

    return SSSSct, CCCCct, FFFFct, FFct


def material_superposition_old(mat_param_cort, mat_param_trab, phic, phit, rhoc, rhot, mm1, mm2, mm3):

    import time

    start_time = time.time()
    # phic = 0.2
    # phit = 0.75
    # rhoc = 0.8
    # rhot = 0.25
    # mm1 = 1.0
    # mm2 = 1.0
    # mm3 = 1.0

    # CORTICAL compliance tensor
    CCCC = numpy.zeros((6, 6))
    CCCC[0, 0] = 1 / mat_param_cort["E0"] / rhoc ** mat_param_cort["KS"]
    CCCC[1, 1] = CCCC[0, 0]
    CCCC[2, 2] = 1 / mat_param_cort["EAA"] / rhoc ** mat_param_cort["KS"]
    CCCC[3, 3] = (
        0.5
        / mat_param_cort["E0"]
        / 2
        / (1 + mat_param_cort["V0"])
        / rhoc ** mat_param_cort["KS"]
    )
    CCCC[4, 4] = 0.5 / mat_param_cort["MUA0"] / rhoc ** mat_param_cort["KS"]
    CCCC[5, 5] = CCCC[4, 4]

    CCCC[1, 0] = (
        -1 * mat_param_cort["V0"] / mat_param_cort["E0"] / rhoc ** mat_param_cort["KS"]
    )
    CCCC[2, 0] = (
        -1
        * mat_param_cort["VA0"]
        / mat_param_cort["EAA"]
        / rhoc ** mat_param_cort["KS"]
    )
    CCCC[2, 1] = CCCC[2, 0]

    CCCC[0, 1] = CCCC[1, 0]
    CCCC[0, 2] = CCCC[2, 0]
    CCCC[1, 2] = CCCC[2, 1]

    # CORTICAL stiffness tensor
    SSSSc = numpy.linalg.inv(CCCC)

    # TRABECULAR stiffness tensor
    lambda_plus_two_mu = (mat_param_trab["E0"] * (1.0 - mat_param_trab["V0"])) / (
        (1.0 + mat_param_trab["V0"]) * (1.0 - 2.0 * mat_param_trab["V0"])
    )
    lambda_zero = (
        mat_param_trab["E0"]
        * mat_param_trab["V0"]
        / ((1.0 + mat_param_trab["V0"]) * (1.0 - 2.0 * mat_param_trab["V0"]))
    )
    SSSSt = numpy.zeros((6, 6))

    SSSSt[0, 0] = (
        mm1 ** (2.0 * mat_param_trab["LS"])
        * rhot ** mat_param_trab["KS"]
        * lambda_plus_two_mu
    )
    SSSSt[1, 1] = (
        mm2 ** (2.0 * mat_param_trab["LS"])
        * rhot ** mat_param_trab["KS"]
        * lambda_plus_two_mu
    )
    SSSSt[2, 2] = (
        mm3 ** (2.0 * mat_param_trab["LS"])
        * rhot ** mat_param_trab["KS"]
        * lambda_plus_two_mu
    )
    SSSSt[3, 3] = (
        2.0
        * mat_param_trab["MU0"]
        * (mm1 ** mat_param_trab["LS"])
        * (mm2 ** mat_param_trab["LS"])
        * (rhot ** mat_param_trab["LS"])
    )
    SSSSt[4, 4] = (
        2.0
        * mat_param_trab["MU0"]
        * (mm3 ** mat_param_trab["LS"])
        * (mm1 ** mat_param_trab["LS"])
        * (rhot ** mat_param_trab["LS"])
    )
    SSSSt[5, 5] = (
        2.0
        * mat_param_trab["MU0"]
        * (mm2 ** mat_param_trab["LS"])
        * (mm3 ** mat_param_trab["LS"])
        * (rhot ** mat_param_trab["LS"])
    )

    SSSSt[1, 0] = (
        lambda_zero
        * (mm1 ** mat_param_trab["LS"])
        * (mm2 * mat_param_trab["LS"])
        * (rhot ** mat_param_trab["LS"])
    )
    SSSSt[2, 0] = (
        lambda_zero
        * (mm3 ** mat_param_trab["LS"])
        * (mm1 * mat_param_trab["LS"])
        * (rhot ** mat_param_trab["LS"])
    )
    SSSSt[2, 1] = (
        lambda_zero
        * (mm2 ** mat_param_trab["LS"])
        * (mm3 * mat_param_trab["LS"])
        * (rhot ** mat_param_trab["LS"])
    )

    SSSSt[0, 1] = SSSSt[1, 0]
    SSSSt[0, 2] = SSSSt[2, 0]
    SSSSt[1, 2] = SSSSt[2, 1]

    # CORTICAL FFFF & FF

    # Material superpositon
    try:
        SSSSct = phic * SSSSc + phit * SSSSt
    except:
        print(
            "WARNING: Material superposition did not work! Stiffness tensor of trabecular bone is used."
        )
        SSSSct = SSSSt

    CCCC = numpy.linalg.inv(SSSSct)

    print("--- %s seconds ---" % (time.time() - start_time))
    return SSSSct, CCCC


def adjust_image_size(image, coarsefactor, crop_z):

    # measure image shape
    IMDimX = numpy.shape(image)[0]
    IMDimY = numpy.shape(image)[1]
    IMDimZ = numpy.shape(image)[2]

    AddDimX = coarsefactor - (IMDimX % coarsefactor)
    AddDimY = coarsefactor - (IMDimY % coarsefactor)

    # adjust in x and y direction
    shape_diff = [AddDimX, AddDimY]
    xy_image_adjusted = numpy.lib.pad(
        image,
        ((0, shape_diff[0]), (0, shape_diff[1]), (0, 0)),
        "constant",
        constant_values=(0),
    )

    if crop_z == CropType.crop:
        image_adjusted = xy_image_adjusted

    if crop_z == CropType.expand:
        AddDimZ = coarsefactor - (IMDimZ % coarsefactor)
        shape_diff = [AddDimX, AddDimY, AddDimZ]
        image_adjusted = numpy.lib.pad(
            xy_image_adjusted, ((0, 0), (0, 0), (0, shape_diff[2])), "edge"
        )

    if crop_z == CropType.variable:
        limit = coarsefactor / 2.0
        if IMDimZ % coarsefactor > limit:
            AddDimZ = coarsefactor - (IMDimZ % coarsefactor)
            shape_diff = [AddDimX, AddDimY, AddDimZ]
            image_adjusted = numpy.lib.pad(
                xy_image_adjusted, ((0, 0), (0, 0), (0, shape_diff[2])), "edge"
            )
        if IMDimZ % coarsefactor < limit:
            image_adjusted = xy_image_adjusted

    return image_adjusted


def create_folder(dir):
    try:
        os.mkdir(dir)
    except OSError:
        if os.path.exists(dir):
            print("%s already exists" % dir)
        else:
            print("Creation of the directory %s failed" % dir)
    else:
        print("Successfully created the directory %s " % dir)


class CropType:
    expand = 0
    crop = 1
    variable = 2


def load_path(filepath):
    """Given a path like /path/to/my_module.pyc (or .py) imports the
    module and returns it
    """

    path, fname = os.path.split(filepath)
    modulename, _ = os.path.splitext(fname)

    if path not in sys.path:
        sys.path.insert(0, path)

    return __import__(modulename)


def pad_array_top_bottom(array, n_pad_layer_proximal, n_pad_layers_distal):
    """returns 1-padded array
    array: numpy array
    n_pad_layers_proximal/distal: number of padding layers in FE elements(int)
    """
    npad = ((0, 0), (0, 0), (n_pad_layer_proximal, n_pad_layers_distal))
    padded_array = numpy.pad(array, pad_width=npad, mode='constant', constant_values=1)
    return padded_array


def save_optim_to_json(optim, config, sample):
    print("save json")
    optim['sample'] = sample

    # optim['OF_value_map'] = numpy.where(optim['OF_value_map'] == 'nan', '', optim['OF_value_map'])
    # optim['OF_value_map_FX'] = numpy.where(optim['OF_value_map_FX'] == 'nan', '', optim['OF_value_map_FX'])
    # optim['OF_value_map_FY'] = numpy.where(optim['OF_value_map_FY'] == 'nan', '', optim['OF_value_map_FY'])
    # optim['OF_value_map_FZ'] = numpy.where(optim['OF_value_map_FZ'] == 'nan', '', optim['OF_value_map_FZ'])
    # optim['OF_value_map_MX'] = numpy.where(optim['OF_value_map_MX'] == 'nan', '', optim['OF_value_map_MX'])
    # optim['OF_value_map_MY'] = numpy.where(optim['OF_value_map_MY'] == 'nan', '', optim['OF_value_map_MY'])
    # optim['OF_value_map_MZ'] = numpy.where(optim['OF_value_map_MZ'] == 'nan', '', optim['OF_value_map_MZ'])


    class NumpyEncoder(json.JSONEncoder):
        def default(self, obj):
            if isinstance(obj, numpy.ndarray):
                return obj.tolist()
            return json.JSONEncoder.default(self, obj)

    key_list = [# General variables
                'sample', 'processing_time',
                # optimized alpha values for 6DOF and single
                'alpha',
                'alpha_single_FX', 'alpha_single_FY', 'alpha_single_FZ',
                'alpha_single_MX', 'alpha_single_MY', 'alpha_single_MZ',
                # objective function values
                'OF_value_optim',
                'OF_value_FX', 'OF_value_FY', 'OF_value_FZ',
                'OF_value_MX', 'OF_value_MY', 'OF_value_MZ',
                # Objective function derivatives
                'OF_derivative_sum_OPT',
                'OF_derivative_sum_FX', 'OF_derivative_sum_MX',
                'OF_derivative_sum_FY', 'OF_derivative_sum_MY',
                'OF_derivative_sum_FZ', 'OF_derivative_sum_MZ',
                # Stiffness and compliance
                'stiffness6D', 'compliance6D',
                'stiff_1D_FZ', 'stiff_projected_FZ', 'stiff_projected_OPT',
                # linear load cases
                'energy_single_FX', 'energy_single_FY', 'energy_single_FZ',
                'energy_single_MX', 'energy_single_MY', 'energy_single_MZ',
                'energy_single_OPT', 'energy_vector_single_OPT',
                # Maximum load cases
                'max_torsor_norm_fitted_OPT_MAX', 'disp_norm_at_max_torsor_norm_fitted_OPT_MAX',
                'nonlinear_energy_at_max_torsor_norm_OPT_MAX','sum_nonlinear_energy_at_max_torsor_norm_OPT_MAX',
                'reached_max_in_norms_OPT_MAX',
                'max_torsor_norm_fitted_FZ_MAX', 'disp_norm_at_max_torsor_norm_fitted_FZ_MAX',
                'nonlinear_energy_at_max_torsor_norm_FZ_MAX', 'sum_nonlinear_energy_at_max_torsor_norm_FZ_MAX',
                'reached_max_in_norms_FZ_MAX',
                'max_torsor_norm_fitted_FX_MAX', 'disp_norm_at_max_torsor_norm_fitted_FX_MAX',
                'nonlinear_energy_at_max_torsor_norm_FX_MAX', 'sum_nonlinear_energy_at_max_torsor_norm_FX_MAX',
                'reached_max_in_norms_FX_MAX',
                'max_torsor_norm_fitted_FY_MAX', 'disp_norm_at_max_torsor_norm_fitted_FY_MAX',
                'nonlinear_energy_at_max_torsor_norm_FY_MAX', 'sum_nonlinear_energy_at_max_torsor_norm_FY_MAX',
                'reached_max_in_norms_FY_MAX',
                'max_torsor_norm_fitted_MX_MAX', 'disp_norm_at_max_torsor_norm_fitted_MX_MAX',
                'nonlinear_energy_at_max_torsor_norm_MX_MAX', 'sum_nonlinear_energy_at_max_torsor_norm_MX_MAX',
                'reached_max_in_norms_MX_MAX',
                'max_torsor_norm_fitted_MY_MAX', 'disp_norm_at_max_torsor_norm_fitted_MY_MAX',
                'nonlinear_energy_at_max_torsor_norm_MY_MAX', 'sum_nonlinear_energy_at_max_torsor_norm_MY_MAX',
                'reached_max_in_norms_MY_MAX',
                'max_torsor_norm_fitted_MZ_MAX', 'disp_norm_at_max_torsor_norm_fitted_MZ_MAX',
                'nonlinear_energy_at_max_torsor_norm_MZ_MAX', 'sum_nonlinear_energy_at_max_torsor_norm_MZ_MAX',
                'reached_max_in_norms_MZ_MAX',
                # mean strain fields linear load cases
                'mean_E_field_raw_FX', 'mean_E_field_optim_FX', 'std_E_field_raw_FX', 'std_E_field_optim_FX',
                'mean_E_field_raw_FY', 'mean_E_field_optim_FY', 'std_E_field_raw_FY', 'std_E_field_optim_FY',
                'mean_E_field_raw_FZ', 'mean_E_field_optim_FZ', 'std_E_field_raw_FZ', 'std_E_field_optim_FZ',
                'mean_E_field_raw_MX', 'mean_E_field_optim_MX', 'std_E_field_raw_MX', 'std_E_field_optim_MX',
                'mean_E_field_raw_MY', 'mean_E_field_optim_MY', 'std_E_field_raw_MY', 'std_E_field_optim_MY',
                'mean_E_field_raw_MZ', 'mean_E_field_optim_MZ', 'std_E_field_raw_MZ', 'std_E_field_optim_MZ',
                'mean_E_field_raw_OPT', 'mean_E_field_optim_OPT', 'std_E_field_raw_OPT', 'std_E_field_optim_OPT',
                'mean_E_field_raw_OPT_theoretical', 'std_E_field_raw_OPT_theoretical',
                ]

    dict_reduced  = {}

    for key in key_list:
        dict_reduced[key] = optim[key]

    try:
        os.mkdir(config['feadir'] + "/optimdicts")
    except:
        pass

    dict_json_file = config["feadir"] + "/" + config["folder_id"][sample] + "/" + "Optimization/" + sample + '_' + \
                     config['current_version'][0:2] + "_optimdict.json"
    dict_json_file2 = config['feadir'] + "/optimdicts/" + sample + '_' + \
                      config['current_version'][0:2] + "_optimdict.json"
    file = open(dict_json_file, 'w')
    json.dump(dict_reduced, file, cls=NumpyEncoder)
    file.close()

    file2 = open(dict_json_file2, 'w')
    json.dump(dict_reduced, file2, cls=NumpyEncoder)
    file2.close()


def write_data_summary(config: dict, optim: dict, bone:dict, sample: str)-> None:
    """
    Function that writes a summary csv file with all OPTIM parameters, for easy postprocessing in R studio
    File is stored in summary folder.
    @param config:
    @param optim:
    @return:
    """

    # Make summary file folder if not yet existing
    summary_path = config['feadir'] + '/summaries'
    try:
        os.mkdir(summary_path)
        print(" Directory '%s' created" %summary_path)
    except:
        pass

    # Create summary file if not yet existing
    filename = summary_path + '/' + config['current_version'] + '_data_summary.csv'

    field_names_dict = [sample,
                        optim['alpha'][0], optim['alpha'][1], optim['alpha'][2],
                        optim['alpha'][3], optim['alpha'][4], optim['alpha'][5],
                        optim['alpha_single_FX'][0], optim['alpha_single_FY'][1], optim['alpha_single_FZ'][2],
                        optim['alpha_single_MX'][3], optim['alpha_single_MY'][4], optim['alpha_single_MZ'][5],
                        optim['OF_value_optim'],
                        optim['OF_value_FX'], optim['OF_value_FY'], optim['OF_value_FZ'],
                        optim['OF_value_MX'], optim['OF_value_MY'], optim['OF_value_MZ'],
                        optim['energy_single_OPT'],
                        optim['energy_single_FX'][0], optim['energy_single_FY'][0], optim['energy_single_FZ'][0],
                        optim['energy_single_MX'][0], optim['energy_single_MY'][0], optim['energy_single_MZ'][0],
                        optim['CV_linear_model_FX'], optim['CV_linear_model_FY'], optim['CV_linear_model_FZ'],
                        optim['CV_linear_model_MX'], optim['CV_linear_model_MY'], optim['CV_linear_model_MZ'],
                        optim['CV_linear_model_OPT'],
                        optim['max_torsor_norm_fitted_OPT_MAX'],
                        optim['max_torsor_norm_fitted_FX_MAX'], optim['max_torsor_norm_fitted_FY_MAX'],
                        optim['max_torsor_norm_fitted_FZ_MAX'], optim['max_torsor_norm_fitted_MX_MAX'],
                        optim['max_torsor_norm_fitted_MY_MAX'], optim['max_torsor_norm_fitted_MZ_MAX'],
                        optim['disp_norm_at_max_torsor_norm_fitted_OPT_MAX'],
                        optim['disp_norm_at_max_torsor_norm_fitted_FX_MAX'], optim['disp_norm_at_max_torsor_norm_fitted_FY_MAX'],
                        optim['disp_norm_at_max_torsor_norm_fitted_FZ_MAX'], optim['disp_norm_at_max_torsor_norm_fitted_MX_MAX'],
                        optim['disp_norm_at_max_torsor_norm_fitted_MY_MAX'], optim['disp_norm_at_max_torsor_norm_fitted_MZ_MAX'],
                        optim['sum_nonlinear_energy_at_max_torsor_norm_OPT_MAX'],
                        optim['sum_nonlinear_energy_at_max_torsor_norm_FX_MAX'], optim['sum_nonlinear_energy_at_max_torsor_norm_FY_MAX'],
                        optim['sum_nonlinear_energy_at_max_torsor_norm_FZ_MAX'], optim['sum_nonlinear_energy_at_max_torsor_norm_MX_MAX'],
                        optim['sum_nonlinear_energy_at_max_torsor_norm_MY_MAX'], optim['sum_nonlinear_energy_at_max_torsor_norm_MZ_MAX'],
                        numpy.mean(optim['mean_E_field_optim_OPT'][0:3]),
                        numpy.mean(optim['mean_E_field_optim_OPT'][3:6]),
                        numpy.mean(optim['mean_E_field_optim_FX'][0:3]),
                        numpy.mean(optim['mean_E_field_optim_FX'][3:6]),
                        numpy.mean(optim['mean_E_field_optim_FY'][0:3]),
                        numpy.mean(optim['mean_E_field_optim_FY'][3:6]),
                        numpy.mean(optim['mean_E_field_optim_FZ'][0:3]),
                        numpy.mean(optim['mean_E_field_optim_FZ'][3:6]),
                        numpy.mean(optim['mean_E_field_optim_MX'][0:3]),
                        numpy.mean(optim['mean_E_field_optim_MX'][3:6]),
                        numpy.mean(optim['mean_E_field_optim_MY'][0:3]),
                        numpy.mean(optim['mean_E_field_optim_MY'][3:6]),
                        numpy.mean(optim['mean_E_field_optim_MZ'][0:3]),
                        numpy.mean(optim['mean_E_field_optim_MZ'][3:6]),
                        numpy.mean(optim['mean_E_field_raw_OPT'][0:3]),
                        numpy.mean(optim['mean_E_field_raw_OPT'][3:6]),
                        numpy.mean(optim['mean_E_field_raw_FX'][0:3]),
                        numpy.mean(optim['mean_E_field_raw_FX'][3:6]),
                        numpy.mean(optim['mean_E_field_raw_FY'][0:3]),
                        numpy.mean(optim['mean_E_field_raw_FY'][3:6]),
                        numpy.mean(optim['mean_E_field_raw_FZ'][0:3]),
                        numpy.mean(optim['mean_E_field_raw_FZ'][3:6]),
                        numpy.mean(optim['mean_E_field_raw_MX'][0:3]),
                        numpy.mean(optim['mean_E_field_raw_MX'][3:6]),
                        numpy.mean(optim['mean_E_field_raw_MY'][0:3]),
                        numpy.mean(optim['mean_E_field_raw_MY'][3:6]),
                        numpy.mean(optim['mean_E_field_raw_MZ'][0:3]),
                        numpy.mean(optim['mean_E_field_raw_MZ'][3:6]),
                        optim['max_force_disp_FZ_MAX'][0][2], optim['max_force_disp_FZ_MAX'][1],
                        optim['force_FZ_MAX'][0][2] / optim['disp_FZ_MAX'][0][2],
                        bone["BMC_tissue"], bone["simulation_BMC_FE_tissue_ROI_orig"],
                        bone["simulation_BMC_FE_tissue_ROI"], bone["BMC_ratio_ROI"],
                        ]

    field_names_titles = ['Sample',
                          'alpha1', 'alpha2', 'alpha3', 'alpha4', 'alpha5', 'alpha6',
                          'alpha_single_FX', 'alpha_single_FY', 'alpha_single_FZ',
                          'alpha_single_MX', 'alpha_single_MY', 'alpha_single_MZ',
                          'OF_value_OPT',
                          'OF_value_FX', 'OF_value_FY', 'OF_value_FZ',
                          'OF_value_MX', 'OF_value_MY', 'OF_value_MZ',
                          'energy_single_OPT',
                          'energy_single_FX', 'energy_single_FY', 'energy_single_FZ',
                          'energy_single_MX', 'energy_single_MY', 'energy_single_MZ',
                          'CV_model_FX', 'CV_model_FY', 'CV_model_FZ',
                          'CV_model_MX', 'CV_model_MY', 'CV_model_MZ',
                          'CV_model_OPT',
                          'max_t_norm_fit_OPT_MAX',
                          'max_t_norm_fit_FX_MAX', 'max_t_norm_fit_FY_MAX', 'max_t_norm_fit_FZ_MAX',
                          'max_t_norm_fit_MX_MAX', 'max_t_norm_fit_MY_MAX', 'max_t_norm_fit_MZ_MAX',
                          'd_norm_at_max_t_norm_OPT_MAX',
                          'd_norm_at_max_t_norm_FX_MAX', 'd_norm_at_max_t_norm_FY_MAX', 'd_norm_at_max_t_norm_FZ_MAX',
                          'd_norm_at_max_t_norm_MX_MAX', 'd_norm_at_max_t_norm_MY_MAX', 'd_norm_at_max_t_norm_MZ_MAX',
                          'energy_at_max_torsor_norm_OPT_MAX',
                          'energy_at_max_torsor_norm_FX_MAX', 'energy_at_max_torsor_norm_FY_MAX',
                          'energy_at_max_torsor_norm_FZ_MAX', 'energy_at_max_torsor_norm_MX_MAX',
                          'energy_at_max_torsor_norm_MY_MAX', 'energy_at_max_torsor_norm_MZ_MAX',
                          'mean_normal_strains_optim_OPT', 'mean_shear_strains_optim_OPT',
                          'mean_normal_strains_optim_FX', 'mean_shear_strains_optim_FX',
                          'mean_normal_strains_optim_FY', 'mean_shear_strains_optim_FY',
                          'mean_normal_strains_optim_FZ', 'mean_shear_strains_optim_FZ',
                          'mean_normal_strains_optim_MX', 'mean_shear_strains_optim_MX',
                          'mean_normal_strains_optim_MY', 'mean_shear_strains_optim_MY',
                          'mean_normal_strains_optim_MZ', 'mean_shear_strains_optim_MZ',
                          'mean_normal_strains_raw_OPT',
                          'mean_shear_strains_raw_OPT',
                          'mean_normal_strains_raw_FX',
                          'mean_shear_strains_raw_FX',
                          'mean_normal_strains_raw_FY',
                          'mean_shear_strains_raw_FY',
                          'mean_normal_strains_raw_FZ',
                          'mean_shear_strains_raw_FZ',
                          'mean_normal_strains_raw_MX',
                          'mean_shear_strains_raw_MX',
                          'mean_normal_strains_raw_MY',
                          'mean_shear_strains_raw_MY',
                          'mean_normal_strains_raw_MZ',
                          'mean_shear_strains_raw_MZ',
                          'max_force_FZ_MAX', 'disp_at_max_force_FZ_MAX',
                           'stiffness_1D_FZ_MAX',
                          'BMC_tissue', 'BMC_FE_tissue_ROI_orig', 'BMC_FE_tissue_ROI', 'BMC_ratio_ROI'
                          ]

    def append_list_as_row(filename: str, list_of_elem: list) -> None:
        with open(filename, 'a+', newline='') as write_obj:
            csv_writer = writer(write_obj)
            csv_writer.writerow(list_of_elem)

    if os.path.exists(filename):
        append_list_as_row(filename, field_names_dict)
    else:
        file = open(filename, 'a+', newline='')
        csv_writer = writer(file)
        csv_writer.writerow(field_names_titles)
        file.close()
        append_list_as_row(filename, field_names_dict)