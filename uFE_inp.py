import os


path = '/home/biomech/Downloads/FE_test/'
inp_bone = 'Mesh_test_A_0121399icotest_A01.inp'
inp_screw = 'Mesh_test_A_0121399icoscrew_A01.inp'
inp_template = 'inp_temp2.inp'

bone_mesh = open(path + inp_bone, 'r')
screw_mesh = open(path + inp_screw, 'r')
template = open(path + inp_template, 'r')
input_file = 'test_A_0.inp'

try:
    os.remove(path + input_file)
except FileNotFoundError:
    print('New file.')

outfile = open(path + input_file, 'w')

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

