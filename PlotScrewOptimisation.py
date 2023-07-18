import matplotlib.pyplot as plt
import numpy as np


def read_RFnodeFile(file_):
    # read data from text file
    df_ = np.loadtxt(file_, delimiter=',')
    uy_ = np.array(df_[:, 1])

    return uy_


# labels = ['Original', 'Original filled', 'Simple filled', 'Simple cannulated']
labels = ['1.5 mm', '0.5 mm', '0.2 mm', '0.2 mm, opt']
c = ['r', 'g', 'b', 'k']
plt.figure()
for i in [0, 1, 2, 3]:
    uy = read_RFnodeFile('/home/biomech/Documents/01_Icotec/02_FEA/00_Model/94_screw_Osteoporosis_new_RFnode'
                         + str(i) + '.txt')

    for samples in [0, 3]:
        if samples == 0:
            plt.plot([3, 2, 1, 0], uy[samples*4:samples*4+4], label=labels[i], ls='-', marker='o', color=c[i])
        else:
            plt.plot([3, 2, 1, 0], uy[samples * 4:samples * 4 + 4], label='_nolegend_', ls='--', marker='o', color=c[i])
    plt.xlabel('RP')
    plt.ylabel('Displacement / mm')
    plt.legend()
