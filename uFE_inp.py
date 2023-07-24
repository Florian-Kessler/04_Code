import os


path = '/home/biomech/Downloads/FE_test/'
inp_bone = 'Mesh_test_0121399icotest0.inp'
inp_screw = 'Mesh_test_0121399icoscrew.inp'
inp_template = 'inp_temp.inp'

bone_mesh = open(path + inp_bone, 'r')
screw_mesh = open(path + inp_screw, 'r')
template = open(path + inp_template, 'r')
try:
    os.remove(path + 'test_0.inp')
except FileNotFoundError:
    print('New file.')

outfile = open(path + 'test_0.inp', 'w')

for lines in template:
    if '** Import Bone' in lines:
        outfile.write('** Bone mesh inserted\n')
        for line in bone_mesh:
            if '*ELEMENT, TYPE=C3D8, ELSET=SET1' in line:
                outfile.write('*ELEMENT, TYPE=C3D8, ELSET=SET-BONE\n')
            elif '*HEADING' in line:
                print('*HEADING removed.')
            elif 'CP3: CUBE' in line:
                print('CP3: CUBE removed.')
            elif '*NSET, NSET=ALL_NODE_T' in line:
                outfile.write('*Nset, nset=ALL_NODE_TB\n')
            elif '*NSET, NSET=ALL_NODE_B' in line:
                print('Combining node sets.')
            else:
                outfile.write(line)
        outfile.write('\n**\n')
    elif '** Import Impl' in lines:
        outfile.write('** Impl mesh inserted\n')
        for line in screw_mesh:
            if '*ELEMENT, TYPE=C3D8, ELSET=SET1' in line:
                outfile.write('*ELEMENT, TYPE=C3D8, ELSET=SET-IMPL\n')
            elif '*HEADING' in line:
                print('*HEADING removed.')
            elif 'CP3: CUBE' in line:
                print('CP3: CUBE removed.')
            else:
                outfile.write(line)
        outfile.write('\n**\n')
    else:
        outfile.write(lines)
outfile.close()

'''
def HFE_inp_creator(inp):
    """
    Creates inputfile with specified details
    :param inp: inputfile properties
    :return: writes input file
    """
    SimMat = ['T', 'P']  # simulated materials
    for i in range(len(SimMat)):
        friction = 0
        f_inpDummy = open(inp['Project_Folder'] + '02_FEA/00_Model/' + inp['Model_Code'] + '_model.inp')
        f_eleSets = open(inp['FEA_loc'] + inp['Model_Code'] + '_elsets.inp')
        f_material = open(inp['FEA_loc'] + inp['Model_Code'] + '_materials.inp')
        # Included in input file-name:
        # Model code (geometry), force maximum, friction, experiment screw material, simulation screw material
        if inp['d_dir'] == '-':
            outfile = open(inp['FEA_loc'] + inp['Model_Code'] + '_d' + str(inp['d_max']) + '_' +
                           str(inp['Friction']).replace('.', '') + '_' + SimMat[i] + '.inp', 'w')
        elif inp['d_dir'] == '+':
            outfile = open(inp['FEA_loc'] + inp['Model_Code'] + '_dinv' + str(inp['d_max']) + '_' +
                           str(inp['Friction']).replace('.', '') + '_' + SimMat[i] + '.inp', 'w')
        else:
            print('Wrong d_dir input')
            exit()
        for lines in f_inpDummy:

            # Friction value between screw and bone
            if friction == 1:
                lines = str(inp['Friction']) + '\n'
            friction = 0  # set friction again to 0
            # Find friction line, set it to 1 to change value during next iteration
            if '*Friction, slip' in lines:
                friction = 1

            # Set material of implant
            if ', material=PEEK' in lines:
                # replace material section
                if SimMat[i] == 'T':
                    # outfile.truncate()
                    outfile.write('*Solid Section, elset=Set-Impl, material=Ti\n')
                    print('Section set to Ti.')
                elif SimMat[i] == 'P':
                    outfile.write('*Solid Section, elset=Set-Impl, material=PEEK\n')
                    print('Section set to PEEK.')
            # Change material properties of implant
            elif '15000., 0.3' in lines:
                outfile.write(inp['YM_peek'] + '., ' + inp['v_peek'] + '\n')
            elif '100000., 0.3' in lines:
                outfile.write(inp['YM_titan'] + '., ' + inp['v_titan'] + '\n')

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

            else:
                outfile.write(lines)
            if '*Heading' in lines:
                outfile.write('** Sample name: ' + inp['FEA_loc'].split('/')[-3] + '\n')

        outfile.close()
        f_inpDummy.close()
        f_eleSets.close()
        f_material.close()
    print("End HFE_inp_creator")
    '''
