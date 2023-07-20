import matplotlib.pyplot as plt
import numpy as np


def read_RFnodeFile(file_):
    # read data from text file
    df_ = np.loadtxt(file_, delimiter=',')
    uy_ = np.array(df_[:, 1])

    return uy_


labels = ['Original', 'Original filled', 'Simple filled', 'Simple cannulated']
# labels = ['1.5 mm', '0.5 mm', '0.2 mm', '0.2 mm, opt']

c = ['r', 'g', 'b', 'k']
fig1, ax1 = plt.subplots()
for i in [1, 2, 3]:
    uy = read_RFnodeFile('/home/biomech/Documents/01_Icotec/02_FEA/00_Model/94_screw_Osteoporosis_new_RFnode'
                         + str(i) + '.txt')

    for samples in [0, 3]:
        if i == 1:
            ax1.plot([4, 3, 2, 1, 0], np.append(uy[samples*4:samples*4+4], 0),
                     label=labels[samples], ls='-', marker='o', color=c[samples])
        else:
            ax1.plot([4, 3, 2, 1, 0], np.append(uy[samples * 4:samples * 4 + 4], 0),
                     label='_nolegend_', ls='--', marker='o', color=c[samples])
    ax1.set_xlabel('RP')
    ax1.set_ylabel('Displacement / mm')
    ax1.legend()
