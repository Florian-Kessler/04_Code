#   Takes input of odb name, nodeset in odb (of a driver o rpilot node),
#   and saves RM into a .txt file for all frames
#
#   Usage
#
#   abaqus python pilot-node-data-extraction.py -o <odb name with .odb extension> -n <pilot driver nodeset name>
#
#   Patrik Wili
#   Juni 2020


def main(ODBn):
    # from __future__ import division
    from abaqus import *
    from abaqusConstants import *
    from odbAccess import *
    from sys import argv, exit
    from math import fabs
    import sys
    from sys import path
    import os
    import platform
    def rightTrim(input, suffix):
        if (input.find(suffix) == -1):
            input = input + suffix
        return input

    ODBname = ODBn + '.odb'
    OUTname = ODBn + '_RFnode_3.txt'
    OUTname2 = ODBn
    OUTname3 = ODBn + '_BDI.txt'
    try:
        os.remove(OUTname)
        print('Old file deleted.')
    except:
        print('Creating new file.')
    odb = openOdb(ODBname)
    a = odb.rootAssembly
    nsetName2 = ['M_SET-3']
    opFile = OUTname

    if os.path.isfile(opFile):  # Append to file if it already exists
        try:
            opFileU = open(opFile, 'a')
        except IOError:
            print('cannot open ', opFile)
            exit(0)
    else:  # Otherwise create the file
        try:
            opFileU = open(opFile, 'w')
        except IOError:
            print('cannot open ', opFile)
            exit(0)
    steps = ['Step-1']  # ['Step-2', 'Step-3', 'Step-4', 'Step-5', 'Step-6', 'Step-7', 'Step-8']
    for s in range(len(steps)):
        for set_ in range(len(nsetName2)):
            step = odb.steps[steps[s]]
            frames = step.frames
            numFrames = len(frames)
            # print(numFrames)
            for i in range(numFrames):  # here start at 1, should be 0 to include first frame

                frame = step.frames[i]
                #field3 = frame.fieldOutputs['RM']
                field4 = frame.fieldOutputs['RF']
                field5 = frame.fieldOutputs['U']
                #field6 = frame.fieldOutputs['UR']
                try:
                    nodeSet2 = a.nodeSets[nsetName2[set_]]
                except:
                    print('Nodeset {0:s} not found in ODB assembly'.format(nsetName2[set_]))
                    sys.exit

                #subField3 = field3.getSubset(region=nodeSet2)
                subField4 = field4.getSubset(region=nodeSet2)
                subField5 = field5.getSubset(region=nodeSet2)
                #subField6 = field6.getSubset(region=nodeSet2)

                opFileU.write("%10d,%10f,%10f\n"  # ,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f\n" \
                              % (i, \
                                 #subField3.values[0].data[0], \
                                 #subField3.values[0].data[1], \
                                 #subField3.values[0].data[2], \
                                 #subField4.values[0].data[0], \
                                 subField4.values[0].data[1], \
                                 #subField4.values[0].data[2], \
                                 #subField5.values[0].data[0], \
                                 subField5.values[0].data[1], \
                                 #subField5.values[0].data[2], \
                                 #subField6.values[0].data[0], \
                                 #subField6.values[0].data[1], \
                                 #subField6.values[0].data[2]
                                 ))

    odb.close()
    opFileU.close()
    print('Done.')
if __name__ == '__main__':
    main('96_screw_Osteoporosis_new_Bending')
    main('97_screw_Osteoporosis_new_Bending_screw_no-nlgeom')
    print('Script has finished.')
