import os
import numpy as np


def copy(spec, mod, doc, dir):
    specimens = open('/home/biomech/Documents/01_Icotec/Specimens.txt', 'r').read().splitlines()  # Read specimens

    user = 'fk21b515'
    server = 'submit03.unibe.ch:'
    remote = '/storage/workspaces/artorg_msb/hpc_main/Florian_Simulation/01_Icotec/01_MainStudy/' \
             + specimens[spec] + '/' + mod + '/'
    local = '/home/biomech/Documents/01_Icotec/02_FEA/01_MainStudy/' + specimens[spec] + '/' + mod + '/'
    pw = open('/home/biomech/Documents/99_General/P_fk.txt', 'r').read()
    if dir == 0:
        os.system('sshpass -p ' + pw + ' scp ' + local + doc + ' ' + user + '@' + server + remote)
        print(doc + ' successfully copied from ' + local + ' \nto \n' + user + '@' + server + remote + '.')
    elif dir == 1:
        os.system('sshpass -p ' + pw + ' scp ' + user + '@' + server + remote + doc + ' ' + local)
        print(doc + ' successfully copied from ' + user + '@' + server + remote + ' \nto \n' + local + '.')
    else:
        print('Invalid direction.')


specimen = np.arange(0, 34)
document = '*02*T.inp'
direction = 0  # 0 = from local to remote, 1 = from remote to local
model = '81_L50_S50_D45'

# [20, 22, 25, 27, 28, 30, 33]:  # range(12, 19):  # range(len(specimen)):
for i in [12, 14, 17, 19, 20, 22, 25, 27, 28, 30, 33]:
    copy(specimen[i], model, document, direction)
