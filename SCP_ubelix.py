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


specimen = np.arange(20, 29)
document = '82*05*.txt'
direction = 1  # 0 = to remote, 1 = from remote
model = '82_L50_S50_D45'
for i in range(len(specimen)):
    copy(specimen[i], model, document, direction)
