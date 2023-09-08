import os


def copy(spec, mod, doc, dir_):
    specimens = open('/home/biomech/Documents/01_Icotec/Specimens.txt', 'r').read().splitlines()  # Read specimens

    user = 'fk21b515'
    server = 'submit03.unibe.ch:'
    remote = '/storage/workspaces/artorg_msb/hpc_main/Florian_Simulation/01_Icotec/01_MainStudy/' \
             + specimens[spec] + '/' + mod + '/'
    local = '/home/biomech/Documents/01_Icotec/02_FEA/01_MainStudy/' + specimens[spec] + '/' + mod + '/'
    pw = open('/home/biomech/Documents/99_General/.P_fk.txt', 'r').read()
    if dir_ == 'to':
        os.system('sshpass -p ' + pw + ' scp ' + local + doc + ' ' + user + '@' + server + remote)
        print(doc + ' successfully copied from ' + local + ' \nto \n' + user + '@' + server + remote + '.')
    elif dir_ == 'from':
        os.system('sshpass -p ' + pw + ' scp ' + user + '@' + server + remote + doc + ' ' + local)
        print(doc + ' successfully copied from ' + user + '@' + server + remote + ' \nto \n' + local + '.\n')
    else:
        print('Invalid direction.')


document = '*.txt'
direction = 'from'  # from local TO REMOTE or FROM REMOTE to local
model = '66_L50_S50_D45'

peek_samples = [2, 5, 7, 8, 10, 13, 15, 16, 18, 21, 23, 24, 26, 29, 31, 32]
ti_samples = [3, 4, 6, 9, 11, 12, 14, 17, 19, 20, 22, 25, 27, 28, 30, 33]
for i in peek_samples:
    copy(i, model, document, direction)
