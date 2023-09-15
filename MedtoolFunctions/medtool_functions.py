from MedtoolFunctions.med_classes import *
import sys
from sys import stdout
import time
import numpy as np
import re
import os
from sys import exit
import dpUtils


def __init__(self, modelName='default'):
    self.modelName = modelName

def throwError(string):
    sys.stdout.write('\n **ERROR** %s\n\n' % string)
    sys.stdout.write('\n E N D E D  with ERRORS \n\n')
    sys.stdout.flush()
    exit(1)


def progressStart(text):
    global curProgress
    sys.stdout.write(text + '|')
    curProgress = 0
    sys.stdout.flush()


def progressNext(progress):
    global curProgress
    if progress > curProgress:
        curProgress += 1
        sys.stdout.write('=')
        sys.stdout.flush()


def progressEnd():
    sys.stdout.write('|\n')
    sys.stdout.flush()


def get_Shape(voxelModel):
    """
    Returns the shape of a voxel model.
    NOT CHANGED!!!
    Changed by Schenk, May 2021 to shapes[0], shapes[1], shapes[2] for X, Y, Z.
    Before it was Z, Y, X. Function is only used for writeAbaqusGeneral.
    """
    shapes = voxelModel.shape
    return (
        shapes[2], shapes[1], shapes[0])


def computeNumpyMinMax(curVoxelModel, echo=0):
    return (
        np.min(curVoxelModel), np.max(curVoxelModel))


def getFilenameAndExtension(fileName):
    """ Function returns file extension and file name. """
    parts = fileName.split('.')
    nparts = len(parts)
    ext = parts[nparts - 1]
    filename = ''
    for part in range(nparts - 1):
        if part < nparts - 2:
            filename = filename + parts[part] + '.'
        else:
            filename = filename + parts[part]

    return (
        filename, ext)


def getShortFilenameAndExtension(fileName):
    filename, ext = getFilenameAndExtension(fileName)
    filename = filename.split('/')
    return (
        filename[len(filename) - 1], ext)


def getAbaqusArgument(string, argName):
    argument = ''
    string = string.replace('\n', '')
    STRING = string.upper()
    ARGNAME = argName.upper()
    if STRING.find(ARGNAME) > 0:
        string1 = string.split(',')
        # string1 = split(string, ',')
        for stringpart in string1:
            stringpart = stringpart.replace(' ', '')
            if stringpart.upper().find(ARGNAME) == 0:
                command, argument = stringpart.split('=')
                # command, argument = split(stringpart, '=')

    else:
        print(" **ERROR** getAbaqusArgument(): Argument '%s' not found!" % argName)
        sys.stdout.write('\n E N D E D  with ERRORS \n\n')
        sys.stdout.flush()
        sys.stdout.flush()
        exit(1)
    return argument


def readAbaqus(inFileName, props=False):
    # print(' ... read Abaqus file       : ', inFileName)
    sys.stdout.flush()
    try:
        inStream = open(inFileName, 'r')
    except IOError:
        sys.stdout.write("\n **ERROR**: '-in' intput file '%s' not found!\n\n" % inFileName)
        sys.stdout.flush()
        sys.stdout.write('\n E N D E D  with ERRORS \n\n')
        sys.stdout.flush()
        exit(1)

    lines = []
    line = ' '
    while line != '':
        line = inStream.readline()
        if line != '':
            lines.append(line)
        LINE = line.upper()
        if LINE.find('*INCLUDE') == 0:
            inputfile = getAbaqusArgument(line, 'INPUT')
            inStream2 = open(inputfile, 'r')
            line2 = ' '
            while line2 != '':
                line2 = inStream2.readline()
                if line2 != '':
                    lines.append(line2)
    read = False
    title = None
    nodes = {}
    elems = {}
    nsets = {}
    elsets = {}
    properties = {}
    unknownElems = []
    lineNo = 0
    while lineNo < len(lines):
        lineNo = lineNo + 1
        line = lines[lineNo - 1]
        LINE = line.upper()
        if LINE.find('*HEADING') == 0:
            lineNo = lineNo + 1
            line = lines[lineNo - 1]
            LINE = line.upper()
            title = lines[lineNo - 1]
            title = title.replace('\n', '')
        if LINE.find('*NODE') == 0 and LINE.upper().find('*NODE PRINT') == -1 and LINE.upper().find(
                '*NODE FILE') == -1 and LINE.upper().find('*NODE OUTPUT') == -1:
            nsetName = None
            if LINE.upper().find('NSET') > 0:
                nsetName = getAbaqusArgument(line, 'NSET')
            while lineNo < len(lines):
                lineNo = lineNo + 1
                line = lines[lineNo - 1]
                if line.find('*') == 0 and line.find('**') == -1:
                    lineNo = lineNo - 1
                    break
                if line.find('**') == 0 or line.find('\n') == 0:
                    pass
                else:
                    vList = line.split(',')
                    nNo = int(vList[0])
                    curNode = node_class(nNo, float(vList[1]))
                    if len(vList) > 2:
                        curNode.set_y(float(vList[2]))
                    if len(vList) > 3:
                        curNode.set_z(float(vList[3]))
                    nodes[nNo] = curNode
                    if nsetName != None:
                        # if nsets.has_key(nsetName):
                        if nsetName in nsets:
                            nsets[nsetName].append(nNo)
                        else:
                            nsets[nsetName] = [nNo]

            continue
        LINE = line.upper()
        if LINE.find('*NSET') == 0:
            print('  -> found *NSET    at Line %s' % repr(lineNo))
            nsetName = getAbaqusArgument(line, 'NSET')
            while lineNo < len(lines):
                lineNo = lineNo + 1
                line = lines[lineNo - 1]
                if line.find('*') == 0 and line.find('**') == -1:
                    lineNo = lineNo - 1
                    break
                if line.find('**') == 0 or line.find('\n') == 0:
                    pass
                else:
                    line = line.replace('\n', '')
                    line = line.replace(' ', '')
                    vList = line.split(',')
                    for Id in vList:
                        if len(Id) > 0:
                            # if nsets.has_key(nsetName):
                            if nsetName in nsets:
                                nsets[nsetName].append(int(Id))
                            else:
                                nsets[nsetName] = [int(Id)]

            continue
        LINE = line.upper()
        if LINE.find('*ELEMENT') == 0 and LINE.upper().find('*ELEMENT OUTPUT') == -1:
            elType = ''
            aElType = getAbaqusArgument(line, 'TYPE')
            aElType = aElType.upper()
            nExpNo = 0
            if aElType.find('B32') == 0:
                elType = 'bar3'
                noExpNo = 3
            elif aElType.find('B3') == 0 or aElType.find('T3') == 0:
                elType = 'bar2'
                noExpNo = 2
            elif aElType.find('CPS3') == 0 or aElType.find('CPE3') == 0 or aElType.find('S3') == 0 or aElType.find(
                    'STRI3') == 0:
                elType = 'tria3'
                noExpNo = 3
            elif aElType.find('STRI65') == 0:
                elType = 'tria6'
                noExpNo = 6
            elif aElType.find('CPS4') == 0 or aElType.find('CPE4') == 0 or aElType.find('S4') == 0:
                elType = 'quad4'
                noExpNo = 4
            elif aElType.find('CPS8') == 0 or aElType.find('CPE8') == 0 or aElType.find('S8') == 0:
                elType = 'quad8'
                noExpNo = 8
            elif aElType.find('C3D4') == 0:
                elType = 'tetra4'
                noExpNo = 4
            elif aElType.find('C3D5') == 0:
                elType = 'pyra5'
                noExpNo = 5
            elif aElType.find('C3D8') == 0 or aElType.find('SC8') == 0:
                elType = 'hexa8'
                noExpNo = 8
            elif aElType.find('C3D6') == 0 or aElType.find('SC6') == 0:
                elType = 'penta6'
                noExpNo = 6
            elif aElType.find('C3D10') == 0:
                elType = 'tetra10'
                noExpNo = 10
            elif aElType.find('C3D15') == 0:
                elType = 'penta15'
                noExpNo = 15
            elif aElType.find('C3D20') == 0:
                elType = 'hexa20'
                noExpNo = 20
            else:
                if aElType not in unknownElems:
                    unknownElems.append(aElType)
                continue
            elsetName = ''
            if LINE.find('ELSET') > 0:
                elsetName = getAbaqusArgument(line, 'ELSET')
            while lineNo < len(lines):
                lineNo += 1
                line = lines[lineNo - 1]
                vList = []
                if line.find('*') == 0 and line.find('**') == -1:
                    lineNo = lineNo - 1
                    break
                if line.find('**') == 0 or line.find('\n') == 0:
                    pass
                else:
                    line = line.replace('\n', '')
                    line = line.replace(',', ' ')
                    vList1 = line.split()
                    # vList1 = split(line)
                    if len(vList1) - 1 != noExpNo:
                        lineNo += 1
                        line = lines[lineNo - 1]
                        line = line.replace('\n', '')
                        line = line.replace(',', ' ')
                        vList2 = line.split()
                        # vList2 = split(line)
                        if len(vList1) + len(vList2) - 1 != noExpNo:
                            lineNo += 1
                            line = lines[lineNo - 1]
                            line = line.replace('\n', '')
                            line = line.replace(',', ' ')
                            vList3 = line.split()
                            # vList3 = split(line)
                            if len(vList1) + len(vList2) + len(vList3) - 1 != noExpNo:
                                sys.stdout.write(
                                    '\n **ERROR**: fec.readAbaqus(): Line %i ff: Number of nodes for this' % (
                                            lineNo - 2))
                                sys.stdout.write('\n            element and expected nodes to not coincide !\n\n')
                                sys.stdout.write('\n E N D E D  with ERRORS \n\n')
                                sys.stdout.flush()
                                exit(1)
                            else:
                                vList = vList1 + vList2 + vList3
                        else:
                            vList = vList1 + vList2
                    else:
                        vList = vList1
                    eNo = int(vList[0])
                    nList = []
                    for nNo in range(1, len(vList)):
                        nList.append(int(vList[nNo]))

                    curElem = element(eNo, nList, elType)
                    elems[eNo] = curElem
                    if elsetName in elsets:
                        # if elsets.has_key(elsetName) > 0:
                        elsets[elsetName].append(eNo)
                    else:
                        elsets[elsetName] = [eNo]

            continue
        if LINE.find('*ELSET') == 0:
            print('\n ** WARNING ** :  *ELSET keyword not supported\n ')
        if LINE.find('*BEAM SECTION') == 0:
            elsetName = getAbaqusArgument(line, 'ELSET')
            sectName = getAbaqusArgument(line, 'SECTION')
            matName = getAbaqusArgument(line, 'MATERIAL')
            if sectName.find('CIRC') == 0:
                lineNo += 1
                line = lines[lineNo - 1]
                data = line.split(',')
                if len(data) == 1:
                    radius = [float(data[0])]
                else:
                    radius = [float(data[0]), float(data[1])]
            else:
                sys.stdout.write('\n ** WARNING ** :  *BEAM SECTION, SECTION=%s not implemented\n ' % sectName)
            properties[elsetName] = {'type': 'BEAM',
                                     'material': matName,
                                     'section': sectName,
                                     'geometry': radius}
            continue
        if LINE.find('*SHELL SECTION') == 0:
            elsetName = getAbaqusArgument(line, 'ELSET')
            matName = getAbaqusArgument(line, 'MATERIAL')
            lineNo += 1
            line = lines[lineNo - 1]
            thickness = float(line)
            properties[elsetName] = {'type': 'SHELL',
                                     'material': matName,
                                     'thickness': thickness}
            continue
        if LINE.find('*SOLID SECTION') == 0:
            elsetName = getAbaqusArgument(line, 'ELSET')
            matName = getAbaqusArgument(line, 'MATERIAL')
            properties[elsetName] = {'type': 'SOLID',
                                     'material': matName}
            lineNo += 1
            continue

    if len(unknownElems) > 0:
        sys.stdout.write("\n **WARNING**: fec.readAbaqus() Element Types '%s' not implemented!\n" % str(unknownElems))
        sys.stdout.flush()
    if props == True:
        return (title,
                nodes,
                nsets,
                elems,
                elsets,
                properties)
    else:
        return (title,
                nodes,
                nsets,
                elems,
                elsets)


def castType(curVoxelModel, format):
    numpyVersion = float(np.__version__[0:3])
    minVox = 10000000
    maxVox = 10000000
    if format == 'B' or format == 'H' or format == 'h':
        maxVox = curVoxelModel.max()
        minVox = curVoxelModel.min()
    if format == 'B':
        if int(minVox) < 0 or int(maxVox) > 255:
            sys.stdout.write(
                '\n **ERROR** castType(). min=%s, max=%s, format=%s!\n' % (repr(minVox), repr(maxVox), format))
            sys.stdout.flush()
            sys.stdout.write(' *********** Use "-scale" option to scale your data from 0..255 first. \n')
            sys.stdout.flush()
            sys.stdout.write('\n E N D E D  with ERRORS \n\n')
            sys.stdout.flush()
            exit(1)
        elif numpyVersion > 1.6:
            curVoxelModel = curVoxelModel.astype(np.uint8, order='F')
        else:
            curVoxelModel = curVoxelModel.astype(np.uint8)
    elif format == 'H':
        if int(minVox) < 0 or int(maxVox) > 65535:
            sys.stdout.write(
                '\n **ERROR** castType(). min=%s, max=%s, format=%s!\n' % (repr(minVox), repr(maxVox), format))
            sys.stdout.flush()
            sys.stdout.write(' *********** Use "-scale" option to scale your data from 0..65535 first. \n')
            sys.stdout.flush()
            sys.stdout.write('\n E N D E D  with ERRORS \n\n')
            sys.stdout.flush()
            exit(1)
        elif numpyVersion > 1.6:
            curVoxelModel = curVoxelModel.astype(np.uint16, order='F')
        else:
            curVoxelModel = curVoxelModel.astype(np.uint16)
    elif format == 'h':
        if int(minVox) < -32768 or int(maxVox) > 32767:
            sys.stdout.write(
                '\n **ERROR** castType(). min=%s, max=%s, format=%s!\n' % (repr(minVox), repr(maxVox), format))
            sys.stdout.flush()
            sys.stdout.write(' *********** Use "-scale" option to scale your data from -32768..+32767 first. \n')
            sys.stdout.flush()
            sys.stdout.write('\n E N D E D  with ERRORS \n\n')
            sys.stdout.flush()
            exit(1)
        elif numpyVersion > 1.6:
            curVoxelModel = curVoxelModel.astype(np.int16, order='F')
        else:
            curVoxelModel = curVoxelModel.astype(np.int16)
    elif format == 'i':
        if numpyVersion > 1.6:
            curVoxelModel = curVoxelModel.astype(np.int32, order='F')
        else:
            curVoxelModel = curVoxelModel.astype(np.int32)
    elif format == 'f':
        if numpyVersion > 1.6:
            curVoxelModel = curVoxelModel.astype(np.float32, order='F')
        else:
            curVoxelModel = curVoxelModel.astype(np.float32)
    else:
        sys.stdout.write('\n **ERROR** castType(). format=%s! not implemented\n' % format)
        sys.stdout.flush()
        sys.stdout.write('\n E N D E D  with ERRORS \n\n')
        sys.stdout.flush()
        exit(1)
    return curVoxelModel


def userSplit(oldString):
    if oldString.find(':') > -1 and oldString.find(';') > -1:
        throwError("Option value '%s' shows a not allowed mixture of  ':' and ';' delimiters!" % oldString)
    newString = oldString.replace(':', ';')
    findList = re.findall('[A-Z];\\\\', newString)
    for val in findList:
        newString = newString.replace(val[0] + ';', val[0] + ':')

    findList = re.findall('[A-Z];/', newString)
    for val in findList:
        newString = newString.replace(val[0] + ';', val[0] + ':')

    return newString.split(';')


def writeAbaqusGeneral(outFileName, curVoxelModel, dimList, templateFile, smooth):
    """
    General Abaqus *.inp file writer. For these materials a default material will be
    applied. Supported commands:
      *USER NODE
      *USER ELEMENT
      *USER NSET, type=point, location=arbitrary
        generate NSET: ARB_NODE_S, ARB_NODE_N, ARB_NODE_E, ARB_NODE_W, ARB_NODE_T, ARB_NODE_B
      *USER NSET, type=point, location=addcorner
        generate NSET: ACOR_NODE_SWB, ACOR_NODE_SEB, ACOR_NODE_NEB, ACOR_NODE_NWB,
                       ACOR_NODE_SWT, ACOR_NODE_SET, ACOR_NODE_NET, ACOR_NODE_NWT
      *USER NSET, type=face
        generate NSET: ALL_NODE_S, ALL_NODE_N, ALL_NODE_E, ALL_NODE_W, ALL_NODE_T, ALL_NODE_B
      *USER ELSET, type=face
        generate ELSET: ALL_S, ALL_ELEM_N, ALL_ELEM_E, ALL_ELEM_W, ALL_ELEM_T, ALL_ELEM_B
      *USER PROPERTY, file=property_temp.inp, range=5:367
        generate multiple material cards, internal variables are "SetName, CardName, GrayValue"
        for the given example: GrayValues > 5 and GrayValues <= 367 are written
        This card can be used multiple times
        If range=... is not given, material cards for all GrayValues are written
    Elements are only written for the given ranges in *USER PROPERTY

    @param outFileName: name of the output file
    @param curVoxelModel: voxel model of the RVE
        - TYPE: np.array[iX, jY, kZ] = grayValue
        - int iX, jY, kZ ... voxels number ID in x,y,z start a 0, x fastest.
        - int grayValue  ... value of voxel
    @param  dimList: list of voxel dimension's
        - TYPE: list[0] = lenX, list[0] = lenY, list[0] = lenZ
        - float lenX, lenY, lenZ ... physical voxel dimensions in x,y,z
    @param  templateFile: name of the template file
    @param  smoothParam: taubin voxel model smoothing parameter
        - TYPE: list[0] = iter, list[1] = lambda, list[2] = kPB,
                list[3] = nearIntf, list[4] = bcid, list[5] = shrink
        - int iter, float lambda, float kPB, int nearIntf

    @return:
      no return value
    """
    sys.stdout.write(' ... setup ABAQUS *.inp file from template\n')
    sys.stdout.write("     -> recast model from '%s' to 'i'\n" % curVoxelModel.dtype.char)
    sys.stdout.flush()
    curVoxelModel = castType(curVoxelModel, 'i')
    time1 = time.time()
    if dimList.all() == None:  # 12.01.01 change: if dimList == None:     to if dimList.all() == None
        print('\n **ERROR** writeAbaqusGeneral(): Voxel size not optional for this function!\n')
        sys.stdout.write('\n E N D E D  with ERRORS \n\n')
        sys.stdout.flush()
        exit(1)
    xvox = dimList[0]
    yvox = dimList[1]
    zvox = dimList[2]
    nx, ny, nz = get_Shape(curVoxelModel)
    minVox, maxVox = computeNumpyMinMax(curVoxelModel, 0)
    minVox = int(minVox + 0.5)
    maxVox = int(maxVox + 0.5)
    noid = 0
    activeNodes = {}
    nodeSets = {}
    nodeSets['ALL_NODE_S'] = []
    nodeSets['ALL_NODE_N'] = []
    nodeSets['ALL_NODE_E'] = []
    nodeSets['ALL_NODE_W'] = []
    nodeSets['ALL_NODE_T'] = []
    nodeSets['ALL_NODE_B'] = []
    elemSets = {}
    elemSets['ALL_ELEM_S'] = []
    elemSets['ALL_ELEM_N'] = []
    elemSets['ALL_ELEM_E'] = []
    elemSets['ALL_ELEM_W'] = []
    elemSets['ALL_ELEM_T'] = []
    elemSets['ALL_ELEM_B'] = []
    tempflag = False
    if templateFile == None:
        tempflag = True
        OS = open('../Imaging_Meshing/temp.inp', 'w')
        OS.write('*USER NODE\n*USER ELEMENT\n*USER PROPERTY, file=prop.inp, range=1:255\n')
        templateFile = 'temp.inp'
        OS.close()
        OS = open('../Imaging_Meshing/prop.inp', 'w')
        OS.write('*SOLID SECTION, ELSET=SetName, MATERIAL=CardName\n1.\n')
        OS.write('*MATERIAL,NAME=CardName\n')
        OS.write('*ELASTIC\n')
        OS.write('20000., 0.3\n')
        OS.close()
        templateFile = 'temp.inp'
    OS = open(outFileName, 'w')
    try:
        osTempFile = open(templateFile, 'r')
    except IOError:
        sys.stdout.write(
            "\n **ERROR** mic.writeAbaqusGeneral(): Abaqus Template file '%s' not found!\n\n" % templateFile)
        sys.stdout.flush()
        sys.stdout.write('\n E N D E D  with ERRORS \n\n')
        sys.stdout.flush()
        exit(1)

    lines = osTempFile.readlines()
    elsetNodes = {}
    ranges = {}
    thresList = []
    rangeMin = 0
    rangeMax = 255
    outFlag = False
    overlap = np.zeros(rangeMax + 1, np.int)
    for line in lines:
        line = line.replace('\n', '')
        if line.upper().find('*USER PROPERTY') == 0:
            line = line.replace(' ', '')
            args = line.split(',')
            matTemplateFilename = ''
            for arg in args:
                if arg.upper().find('RANGE') == 0:
                    dummy, rangeStr = arg.split('=')
                    rangeMin, rangeMax = userSplit(rangeStr)
                    rangeMin = int(rangeMin)
                    rangeMax = int(rangeMax)
                    if rangeMin < 1:
                        sys.stdout.write('\n **ERROR** mic.writeAbaqusGeneral(): Minimum Range < 1!\n\n')
                        sys.stdout.write('\n E N D E D  with ERRORS \n\n')
                        sys.stdout.flush()
                        exit(1)
                    if rangeMax > maxVox:
                        outFlag = True
                    for ii in range(rangeMax - rangeMin + 1):
                        overlap[rangeMin + ii] += 1

                if arg.upper().find('FILE') == 0:
                    dummy, matTemplateFilename = arg.split('=')

            ranges[matTemplateFilename] = (
                rangeMin, rangeMax)

    if len(ranges) == 0:
        sys.stdout.write('\n **ERROR** mic.writeAbaqusGeneral(): *USER PROPERTY: keyword missing!\n\n')
        sys.stdout.write('\n E N D E D  with ERRORS \n\n')
        sys.stdout.flush()
        exit(1)
    if rangeMax > maxVox:
        sys.stdout.write(
            '\n **WARNING** mic.writeAbaqusGeneral(): *USER PROPERTY: Max GV Range (%i) > Max Image GV (%i)!\n\n' % (
                rangeMax, maxVox))
    if np.sum(np.greater(overlap, 1)) > 0:
        sys.stdout.write(
            '\n **ERROR** mic.writeAbaqusGeneral(): *USER PROPERTY: Ranges in property template overlap!\n\n')
        sys.stdout.write('\n E N D E D  with ERRORS \n\n')
        sys.stdout.flush()
        exit(1)
    for crange in ranges:
        for matId in range(ranges[crange][0], ranges[crange][1] + 1):
            # print('matID', matId)
            elsetNodes[repr(matId)] = []
            thresList.append(matId)

    elid = 0
    nx1 = nx + 1
    nxy1 = (ny + 1) * (nx + 1)
    sum = 0
    progressStart('     -> setup Element Data  : ')
    for k in range(nz):
        sum += 1
        progress = float(sum) / float(nz) * 10.0
        for j in range(ny):
            for i in range(nx):
                grayValue = curVoxelModel[k, j, i]
                # print('\ni: %s, j: %s, k:%s' % (i, j, k))
                # print(grayValue)
                # print(elsetNodes['1'])
                # print(elsetNodes[repr(grayValue)])
                if repr(grayValue) in elsetNodes:
                    # if elsetNodes.has_key(repr(grayValue)):
                    elid = elid + 1
                    elnds = [nxy1 * k + nx1 * j + (i + 1),
                             nxy1 * k + nx1 * j + (i + 2),
                             nxy1 * k + nx1 * (j + 1) + (i + 2),
                             nxy1 * k + nx1 * (j + 1) + (i + 1),
                             nxy1 * (k + 1) + nx1 * j + (i + 1),
                             nxy1 * (k + 1) + nx1 * j + (i + 2),
                             nxy1 * (k + 1) + nx1 * (j + 1) + (i + 2),
                             nxy1 * (k + 1) + nx1 * (j + 1) + (i + 1)]
                    elsetNodes[repr(grayValue)].append((elid, elnds))
                    if k == 0:
                        elemSets['ALL_ELEM_B'].append(elid)
                    if k == nz - 1:
                        elemSets['ALL_ELEM_T'].append(elid)
                    if j == 0:
                        elemSets['ALL_ELEM_S'].append(elid)
                    if j == ny - 1:
                        elemSets['ALL_ELEM_N'].append(elid)
                    if i == 0:
                        elemSets['ALL_ELEM_W'].append(elid)
                    if i == nx - 1:
                        elemSets['ALL_ELEM_E'].append(elid)

        progressNext(progress)

    progressEnd()
    sys.stdout.write('     -> setup Node Data     :')
    for matid in thresList:
        if len(elsetNodes[repr(matid)]) > 0:
            matidStr = 'SET' + repr(matid)
            for elnds in elsetNodes[repr(matid)]:
                elid = elnds[0]
                for elnd in elnds[1]:
                    activeNodes[elnd] = 1

    noid = 0
    for k in range(nz + 1):
        for j in range(ny + 1):
            for i in range(nx + 1):
                noid = noid + 1
                if noid in activeNodes:
                    # if activeNodes.has_key(noid):
                    if k == 0:
                        nodeSets['ALL_NODE_B'].append(noid)
                    if k == nz:
                        nodeSets['ALL_NODE_T'].append(noid)
                    if j == 0:
                        nodeSets['ALL_NODE_S'].append(noid)
                    if j == ny:
                        nodeSets['ALL_NODE_N'].append(noid)
                    if i == 0:
                        nodeSets['ALL_NODE_W'].append(noid)
                    if i == nx:
                        nodeSets['ALL_NODE_E'].append(noid)

    sys.stdout.write(' Done\n')
    sys.stdout.flush()
    nodeCoord = {}
    nodeCoordOrig = {}

    if smooth != None:
        print('error: dpMesherf77 needed...')
    #     activeNodes2, nElem, nNode, nIntElem, nIntFaces, nIntNode = dpMesherf77.check_voxmesh2d(curVoxelModel,
    #                                                                                             smooth[4], 2)
    #     nodeCoordF77, nodeCoordInt, noidF77PY, noidIntVoxF77 = dpMesherf77.get_voxmesh2d_nodes(curVoxelModel,
    #                                                                                            dimList, smooth,
    #                                                                                            activeNodes2, nNode,
    #                                                                                            nIntElem, nIntNode,
    #                                                                                            2)
    #     for noidF77 in range(len(noidF77PY)):
    #         noid = noidF77PY[noidF77]
    #         nodeCoord[noid] = (nodeCoordF77[noidF77][0], nodeCoordF77[noidF77][1], nodeCoordF77[noidF77][2])
    else:
        noid = 0
        for k in range(nz + 1):
            for j in range(ny + 1):
                for i in range(nx + 1):
                    noid = noid + 1
                    if noid in activeNodes:
                        # if activeNodes.has_key(noid):
                        nodeCoord[noid] = (
                            float(xvox * i), float(yvox * j), float(zvox * k))

    curPathFilename, ext = getFilenameAndExtension(outFileName)
    curFilename, ext = getShortFilenameAndExtension(outFileName)
    sys.stdout.write(' ... write ABAQUS *.inp file from template\n')
    for line in lines:
        line = line.replace('\n', '')
        line = line.replace('$filename', curFilename)
        line = line.replace('$pathfilename', curPathFilename)
        if line.upper().find('*USER NODE') > -1:
            OS.write('*NODE\n')
            noid2 = 0
            noid = 0
            progressStart('     -> process Node IDs    : ')
            for k in range(nz + 1):
                progress = float(k + 1) / float(nz + 1) * 10.0
                for j in range(ny + 1):
                    for i in range(nx + 1):
                        noid = noid + 1
                        if noid in activeNodes:
                            # if activeNodes.has_key(noid):
                            noid2 = noid2 + 1
                            OS.write('%12i,%13.6g,%13.6g,%13.6g\n' % (
                                noid, nodeCoord[noid][0], nodeCoord[noid][1], nodeCoord[noid][2]))

                progressNext(progress)

            progressEnd()
            sys.stdout.write('     -> write Nodes         : %10i \n' % noid2)
            sys.stdout.flush()
        elif line.upper().find('*USER ELEMENT') > -1:
            count = 0
            progressStart('     -> process Elements    : ')
            for matid in thresList:
                count += 1
                progress = count / float(len(thresList)) * 10.0
                if len(elsetNodes[repr(matid)]) > 0:
                    matidStr = 'SET' + repr(matid)
                    OS.write('*ELEMENT, TYPE=C3D8, ELSET=%s\n' % matidStr)
                    for elnds in elsetNodes[repr(matid)]:
                        elid = elnds[0]
                        OS.write('%s,%s,%s,%s,%s,%s,%s,%s,%s\n' % (
                            elid, elnds[1][0], elnds[1][1], elnds[1][2], elnds[1][3], elnds[1][4], elnds[1][5],
                            elnds[1][6], elnds[1][7]))
                        for elnd in elnds[1]:
                            activeNodes[elnd] = 1

                progressNext(progress)

            progressEnd()
            sys.stdout.write('     -> write Elements      : %10i             \n' % elid)
            sys.stdout.flush()
        elif line.upper().find('*USER NSET') > -1:
            if line.upper().find('TYPE=FACE') > -1:
                sys.stdout.write('     -> write BCs Node Sets     \n')
                sys.stdout.flush()
                for nsetName in nodeSets:
                    OS.write('*NSET, NSET=%s\n' % nsetName)
                    entry = 0
                    for noid in nodeSets[nsetName]:
                        entry = entry + 1
                        if entry == 16:
                            OS.write('%s' % repr(noid))
                            entry = 0
                            OS.write('\n')
                        else:
                            OS.write('%s,' % repr(noid))

                    OS.write('\n')

            if line.upper().find('TYPE=POINT') > -1:
                if line.upper().find('LOCATION=ARBITRARY') > -1:
                    for nsetName in nodeSets:
                        if len(nodeSets[nsetName]) > 0:
                            nid = nodeSets[nsetName][0]
                            name = nsetName.replace('ALL_NODE_', 'ARB_NODE_')
                            OS.write('*NSET, NSET=%s\n' % name)
                            OS.write('%s\n' % repr(nid))

                if line.upper().find('LOCATION=ADDCORNER') > -1:
                    nid = (nx + 1) * (ny + 1) * (nz + 1)
                    OS.write('*NODE, NSET=ACOR_NODE_SWB\n')
                    OS.write('%i, %13.6g,  %13.6g, %13.6g\n' % (nid + 1, 0.0, 0.0, 0.0))
                    OS.write('*NODE, NSET=ACOR_NODE_SEB\n')
                    OS.write('%i, %13.6g,  %13.6g, %13.6g\n' % (nid + 2, nx * xvox, 0.0, 0.0))
                    OS.write('*NODE, NSET=ACOR_NODE_NEB\n')
                    OS.write('%i, %13.6g,  %13.6g, %13.6g\n' % (nid + 3, nx * xvox, ny * yvox, 0.0))
                    OS.write('*NODE, NSET=ACOR_NODE_NWB\n')
                    OS.write('%i, %13.6g,  %13.6g, %13.6g\n' % (nid + 4, 0.0, ny * yvox, 0.0))
                    OS.write('*NODE, NSET=ACOR_NODE_SWT\n')
                    OS.write('%i, %13.6g,  %13.6g, %13.6g\n' % (nid + 5, 0.0, 0.0, nz * zvox))
                    OS.write('*NODE, NSET=ACOR_NODE_SET\n')
                    OS.write('%i, %13.6g,  %13.6g, %13.6g\n' % (nid + 6, nx * xvox, 0.0, nz * zvox))
                    OS.write('*NODE, NSET=ACOR_NODE_NET\n')
                    OS.write('%i, %13.6g,  %13.6g, %13.6g\n' % (nid + 7, nx * xvox, ny * yvox, nz * zvox))
                    OS.write('*NODE, NSET=ACOR_NODE_NWT\n')
                    OS.write('%i, %13.6g,  %13.6g, %13.6g\n' % (nid + 8, 0.0, ny * yvox, nz * zvox))
        elif line.upper().find('*USER ELSET') > -1:
            if line.upper().find('TYPE=FACE') > -1:
                sys.stdout.write('     -> Write BCs Elem Sets          \n')
                sys.stdout.flush()
                for elsetName in elemSets:
                    OS.write('*ELSET, ELSET=%s\n' % elsetName)
                    entry = 0
                    for elid in elemSets[elsetName]:
                        entry = entry + 1
                        if entry == 16:
                            OS.write('%s' % repr(elid))
                            entry = 0
                            OS.write('\n')
                        else:
                            OS.write('%s,' % repr(elid))

                    OS.write('\n')

        elif line.upper().find('*USER PROPERTY') > -1:
            line = line.replace(' ', '')
            args = line.split(',')
            rangeMin = minVox
            rangeMax = maxVox
            matTemplateFilename = ''
            for arg in args:
                if arg.upper().find('RANGE') == 0:
                    dummy, rangeStr = arg.split('=')
                    rangeMin, rangeMax = userSplit(rangeStr)
                    rangeMin = int(rangeMin)
                    rangeMax = int(rangeMax)
                if arg.upper().find('FILE') == 0:
                    dummy, matTemplateFilename = arg.split('=')

            sys.stdout.write('     -> Write Property      : %s \n' % matTemplateFilename)
            try:
                osMatCard = open(matTemplateFilename, 'r')
            except IOError:
                sys.stdout.write(
                    "\n **ERROR** writeAbaqusGeneral(): Material template file '%s' not found!\n\n" % matTemplateFilename)
                sys.stdout.flush()
                sys.stdout.write('\n E N D E D  with ERRORS \n\n')
                sys.stdout.flush()
                exit(1)

            lines = osMatCard.readlines()
            for matid in thresList:
                GrayValue = matid
                if len(elsetNodes[repr(matid)]) > 0:
                    if matid >= rangeMin and matid <= rangeMax:
                        matidStr = 'SET' + repr(matid)
                        GrayValue = matid
                        for line in lines:
                            line = line.replace('\n', '')
                            if line.find('SetName') > -1:
                                line = line.replace('SetName', matidStr)
                            if line.find('CardName') > -1:
                                line = line.replace('CardName', 'MAT' + matidStr)
                            if line.find('GrayValue') > -1:
                                exprList = line.split(',')
                                count = 0
                                for expr in exprList:
                                    if expr.find('GrayValue') > -1:
                                        compValue = eval(expr)
                                        OS.write('%s' % repr(compValue))
                                    else:
                                        OS.write('%s' % expr)
                                    if count < len(exprList) - 1:
                                        OS.write(',')
                                    count += 1

                                OS.write('\n')
                            else:
                                OS.write('%s\n' % line)

            osMatCard.close()
        else:
            OS.write('%s\n' % line)

    osTempFile.close()
    if tempflag:
        os.remove('../Imaging_Meshing/temp.inp')
        os.remove('../Imaging_Meshing/prop.inp')
    time2 = time.time()
    sys.stdout.write('     -> Write finished in   :   %8.1f sec  \n' % (time2 - time1))
    sys.stdout.flush()
    OS.close
    return


def writeAbaqus(outFileName, title, nodes, nsets, elems, elsets, NscaResults=None):
    time1 = time.time()
    print(' ... write ABAQUS file       : ', outFileName)
    sys.stdout.flush()
    keys = list(nodes.keys())
    nkey1 = keys[1]
    del keys
    noSpatDim = nodes[nkey1].get_dimension()
    os = open(outFileName, 'w')
    if not title == None:
        os.write('*HEADING\n')
        os.write('%s\n' % title)
    os.write('***********************************************************\n')
    os.write('*NODE\n')
    for nodeId in nodes:
        os.write('%s, ' % repr(nodeId))
        os.write('%13.7e, ' % nodes[nodeId].get_x())
        if noSpatDim > 1:
            os.write('%13.7e, ' % nodes[nodeId].get_y())
        else:
            os.write(', ')
        if noSpatDim > 2:
            os.write('%13.7e ' % nodes[nodeId].get_z())
        os.write('\n')

    if NscaResults != None:
        os.write('***********************************************************\n')
        os.write('*NODAL THICKNESS\n')
        nodeThick = NscaResults[0]
        for nodeId in nodes:
            os.write('%s, ' % repr(nodeId))
            os.write('%13.7e\n' % nodeThick[nodeId])

    os.write('***********************************************************\n')
    if nsets != None:
        if len(nsets) > 0:
            for setName in nsets:
                if setName != '':
                    os.write('*NSET, NSET=%s\n' % setName)
                    count = 0
                    for nodeId in nsets[setName]:
                        count += 1
                        if count == 16:
                            os.write('%s' % nodeId)
                            os.write('\n')
                            count = 0
                        else:
                            os.write('%s, ' % nodeId)

                    if count != 0:
                        os.write('\n')

    else:
        os.write('** no NSET written\n')
    os.write('***********************************************************\n')
    if elsets != None:
        if len(elsets) > 0:
            for setName in elsets:
                aElType = ''
                elType = elems[elsets[setName][0]].get_type()
                if elType == 'bar2':
                    aElType = 'B3'
                elif elType == 'tria3':
                    aElType = 'S3'
                elif elType == 'quad4':
                    aElType = 'S4'
                elif elType == 'penta6':
                    aElType = 'C3D6'
                elif elType == 'hexa8':
                    aElType = 'C3D8'
                elif elType == 'tetra4':
                    aElType = 'C3D4'
                elif elType == 'pyra5':
                    aElType = 'C3D5'
                elif elType == 'bar3':
                    aElType = 'B32'
                elif elType == 'tria6':
                    aElType = 'STRI65'
                elif elType == 'quad8':
                    aElType = 'S8'
                elif elType == 'penta15':
                    aElType = 'C3D15'
                elif elType == 'hexa20':
                    aElType = 'C3D20'
                elif elType == 'tetra10':
                    aElType = 'C3D10'
                else:
                    sys.stdout.write(
                        "\n **ERROR** writeAbaqus() : Element Type '%s' not implemented!\n\n" % repr(elType))
                    sys.stdout.flush()
                    sys.stdout.write('\n E N D E D  with ERRORS \n\n')
                    sys.stdout.flush()
                    exit(1)
                os.write('*ELEMENT, TYPE=%s, ELSET=%s\n' % (aElType, setName))
                for elId in elsets[setName]:
                    os.write('%s' % repr(elId))
                    count = 1
                    for node in elems[elId].get_nodes():
                        count += 1
                        if count == 8:
                            os.write(', %s,\n' % repr(node))
                            count = 0
                        elif count == 1:
                            os.write('%s' % repr(node))
                        else:
                            os.write(', %s' % repr(node))

                    if count != 0:
                        os.write('\n')

    else:
        os.write('** no ELEMENTS and ELSET written\n')
    time2 = time.time()
    print('     -> write finished in    :   %8.1f sec' % (time2 - time1))
    os.close
    return


def writeEnsight(self, outFileName, title, nodes, nsets, elems, elsets, NscaResults=None, EscaResults=None,
                 vecResults=None, EvecResults=None):
    AllEResults = {}
    if elsets and len(elsets) > 0:
        for elset in elsets:
            AllEResults[elems[elsets[elset][0]].get_type()] = {}

        matid = {}
        i = 1
        for setName in elsets:
            matid[setName] = i
            i += 1

        for setName in elsets:
            for elid in elsets[setName]:
                AllEResults[elems[elid].get_type()][elid] = float(matid[setName])

        for elid in AllEResults:
            if EscaResults:
                EscaResults.append(AllEResults[elid])
            else:
                EscaResults = [
                    AllEResults[elid]]

    fileList = outFileName.split('/')
    ffilename = fileList.pop()
    filename, ext = self.getFilenameAndExtension(ffilename)
    pathFilename, ext = self.getFilenameAndExtension(outFileName)
    cfilename = filename + '.case'
    print(' ... write Ensight case file :', cfilename)
    caseOS = open(pathFilename + '.case', 'w')
    caseOS.write('FORMAT\n')
    caseOS.write('type: ensight gold\n\n')
    caseOS.write('GEOMETRY\n')
    caseOS.write('model:                           %s\n\n' % (filename + '.geo'))
    if NscaResults != None or EscaResults != None or vecResults != None or EvecResults != None:
        caseOS.write('VARIABLE\n')
        if vecResults != None:
            vres = 0
            for vecRes in vecResults:
                vres += 1
                caseOS.write('vector per node:  %10i Vector  %s\n' % (1, filename + '.vec' + repr(vres)))

        if EvecResults != None:
            vres = 0
            for vecRes in EvecResults:
                vres += 1
                caseOS.write('vector per element:  %10i Vector  %s\n' % (1, filename + '.Evec' + repr(vres)))

        if NscaResults != None:
            sres = 0
            for scaRes in NscaResults:
                sres += 1
                caseOS.write('scalar per node:  %10i Scalar  %s\n' % (1, filename + '.sca' + repr(sres)))

        if EscaResults != None:
            sres = 0
            for scaRes in EscaResults:
                sres += 1
                caseOS.write(
                    'scalar per element:  ESca%s%s  %s\n' % (repr(sres), filename, filename + '.Esca' + repr(sres)))

    caseOS.close()
    gfilename = filename + '.geo'
    print(' ... write Ensight Data File :', gfilename)
    geoOS = open(pathFilename + '.geo', 'w')
    geoOS.write('Title: %s\n' % title)
    geoOS.write('Description 2\n')
    geoOS.write('node id given\n')
    geoOS.write('element id given\n')
    geoOS.write('part \n')
    geoOS.write('         1\n')
    geoOS.write('Description PART \n')
    geoOS.write('coordinates \n')
    geoOS.write('%10i\n' % len(nodes))
    nnode = 0
    ensightNodeDict = {}
    for key in nodes.keys():
        nnode += 1
        geoOS.write('%10i\n' % key)
        ensightNodeDict[key] = nnode

    for key in nodes.keys():
        geoOS.write('%12.5e\n' % nodes[key].get_x())

    for key in nodes.keys():
        geoOS.write('%12.5e\n' % nodes[key].get_y())

    for key in nodes.keys():
        geoOS.write('%12.5e\n' % nodes[key].get_z())

    geoOS.write('\n')
    elTypes = {}
    for key in elems.keys():
        elType = elems[key].get_type()
        # if not elTypes.has_key(elType):
        if elType not in elTypes:
            elTypes[elType] = 1
        else:
            elTypes[elType] += 1

    for elType in elTypes:
        geoOS.write('%s\n' % elType)
        geoOS.write('%10i\n' % elTypes[elType])
        for key in elems.keys():
            if elems[key].get_type() == elType:
                geoOS.write('%10i\n' % key)

        for key in elems.keys():
            if elems[key].get_type() == elType:
                nodeList = elems[key].get_nodes()
                for noid in nodeList:
                    geoOS.write('%10i' % ensightNodeDict[noid])

                geoOS.write('\n')

    geoOS.write('\n')
    geoOS.write('\n')
    geoOS.close()
    if NscaResults != None or vecResults != None or EscaResults != None or EvecResults != None:
        if vecResults != None:
            vres = 0
            for vecRes in vecResults:
                vres += 1
                filename2 = pathFilename + '.vec' + repr(vres)
                disOS = open(filename2, 'w')
                disOS.write('Disp\n')
                disOS.write('part \n')
                disOS.write('         1\n')
                disOS.write('coordinates\n')
                for dir in range(3):
                    for resId in vecRes:
                        disOS.write('%12.5e\n' % vecRes[resId][dir])

                disOS.close()

        if EvecResults != None:
            vres = 0
            for vecRes in EvecResults:
                vres += 1
                filename2 = pathFilename + '.Evec' + repr(vres)
                disOS = open(filename2, 'w')
                disOS.write('Disp\n')
                disOS.write('part \n')
                disOS.write('         1\n')
                disOS.write('%s\n' % elType)
                for dir in range(3):
                    for resId in vecRes:
                        disOS.write('%12.5e\n' % vecRes[resId][dir])

                disOS.close()

        if NscaResults != None:
            sres = 0
            for scaRes in NscaResults:
                sres += 1
                filename2 = pathFilename + '.sca' + repr(sres)
                strOS = open(filename2, 'w')
                strOS.write('%s \n' % ('Scalar' + repr(sres)))
                strOS.write('part \n')
                strOS.write('         1\n')
                strOS.write('coordinates\n')
                for resId in scaRes:
                    strOS.write('%12.5e\n' % scaRes[resId])

                strOS.close()

        if EscaResults != None:
            sres = 0
            for scaRes in EscaResults:
                sres += 1
                filename2 = pathFilename + '.Esca' + repr(sres)
                strOS = open(filename2, 'w')
                strOS.write('%s \n' % ('ESca' + repr(sres) + filename))
                strOS.write('part \n')
                strOS.write('         1\n')
                for resId in scaRes:
                    elType = elems[resId].get_type()
                    break

                strOS.write('%s\n' % elType)
                for resId in scaRes:
                    strOS.write('%12.5e\n' % scaRes[resId])

                strOS.close()

    return

def clean(self, curVoxelModel, thresList, cleanList, echo=False):
    """
    Functions cleans the voxel model - removes islands and not proper connected regions.
    @param curVoxelModel: voxel model of the RVE
          - TYPE: numpy.array[iZ, jY, kX] = grayValue
          - int iX, jY, kZ ... voxels number ID in x,y,z start a 0, x fastest.
          - int grayValue  ... gray value of voxel,
    @param  thresList:  list of thresholds - only one threshold implemented!
          - TYPE: list[thresId] = thres
          - int thresID ... threshold id, sorted ascending
          - int thres   ... threshold value 0..25
    @param  cleanList: list of clean parameters
          - TYPE: list[0] = nodeShared, list[1] = islandSize
          - int nodeShared ... number of shared nodes of a bone to bone connection
          - int islandSize ... number of Voxel isolated bone/marrow voxel region
              which should be removed.
    @param echo: Flag if extended echo should be written on stdout.
          - TYPE: bool

    @return:
       cleanVoxelModel: cleaned voxel model of the RVE
          - TYPE: numpy.array[iZ, jY, kX] = grayValue
          - int iX, jY, kZ ... voxels number ID in x,y,z start a 0, x fastest.
          - int grayValue  ... gray value of voxel,
    """
    print(' ... clean model')
    nodeShared = cleanList[0]
    islandSize = cleanList[1]
    time1 = time.clock()
    ok = False
    nx, ny, nz = self.get_Shape(curVoxelModel)
    cleanVoxelModel = self.createVoxelModel(nx, ny, nz, 'f')
    cleanVoxelModel = curVoxelModel
    step = 0
    while not ok:
        step = step + 1
        stdout.write('     Analyses Step          : %10i       \n' % step)
        stdout.flush()
        sumid = 0
        if thresList == None or len(thresList) != 1:
            stdout.write('\n **ERROR** find_parts(): Exactly one threshold needed for this function!\n\n')
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        partMatNo = 0
        checkList = []
        partIdElemDict = {}
        partIdsVolume = np.zeros((nz, ny, nx))
        nbPartIdList = []
        matid = thresList[0]
        dpUtils.progressStart('     -> look for neighbours : ')
        for k in range(nz):
            progress = float(sumid) / float(nx * ny * nz) * 10.0
            for j in range(ny):
                for i in range(nx):
                    sumid = sumid + 1
                    if i - 1 >= 0:
                        checkList.append((k, j, i - 1))
                    if j - 1 >= 0:
                        checkList.append((k, j - 1, i))
                    if k - 1 >= 0:
                        checkList.append((k - 1, j, i))
                    if nodeShared <= 2:
                        if i - 1 >= 0 and j - 1 >= 0:
                            checkList.append((k, j - 1, i - 1))
                        if i - 1 >= 0 and k - 1 >= 0:
                            checkList.append((k - 1, j, i - 1))
                        if k - 1 >= 0 and j - 1 >= 0:
                            checkList.append((k - 1, j - 1, i))
                    if nodeShared <= 1:
                        if i - 1 >= 0 and j - 1 >= 0 and k - 1 >= 0:
                            checkList.append((k - 1, j - 1, i - 1))
                    for checkElem in checkList:
                        iN = checkElem[2]
                        jN = checkElem[1]
                        kN = checkElem[0]
                        if cleanVoxelModel[k, j, i] == cleanVoxelModel[kN, jN, iN]:
                            partId = partIdsVolume[kN, jN, iN]
                            if partId not in nbPartIdList:
                                nbPartIdList.append(partId)

                    partIdNo = len(nbPartIdList)
                    if partIdNo == 0:
                        partMatNo = partMatNo + 1
                        partIdsVolume[k, j, i] = partMatNo
                        partIdElemDict[partMatNo] = [(k, j, i)]
                    elif partIdNo == 1:
                        partIdsVolume[k, j, i] = nbPartIdList[0]
                        partIdElemDict[nbPartIdList[0]].append((k, j, i))
                    elif partIdNo > 1:
                        maxElem = 0
                        newPartId = 0
                        for partId in nbPartIdList:
                            if len(partIdElemDict[partId]) > maxElem:
                                maxElem = len(partIdElemDict[partId])
                                newPartId = partId

                        nbPartIdList.remove(newPartId)
                        for partId in nbPartIdList:
                            for element in partIdElemDict[partId]:
                                partIdElemDict[newPartId].append(element)

                            for elemTuple in partIdElemDict[partId]:
                                partIdsVolume[elemTuple[0], elemTuple[1], elemTuple[2]] = newPartId

                            del partIdElemDict[partId]

                        partIdsVolume[k, j, i] = newPartId
                        partIdElemDict[newPartId].append((k, j, i))
                    del nbPartIdList[0:len(nbPartIdList)]
                    del checkList[0:len(checkList)]

            dpUtils.progressNext(progress)

        dpUtils.progressEnd()
        ok = self.updateModel(partIdElemDict, cleanVoxelModel, thresList, islandSize)
        if echo == True or ok == True:
            for partId in partIdElemDict.keys():
                ck = partIdElemDict[partId][0][0]
                cj = partIdElemDict[partId][0][1]
                ci = partIdElemDict[partId][0][2]
                stdout.write('     -> Part %4i;  Elements in Part = %8i;  Material = %3i\n' % (
                partId, len(partIdElemDict[partId]), cleanVoxelModel[ck, cj, ci]))
                stdout.flush()

    time2 = time.clock()
    stdout.write('     -> clean finished in   :   %8.1f sec         \n' % (time2 - time1))
    stdout.flush()
    return cleanVoxelModel