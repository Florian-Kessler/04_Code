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
fig2, ax2 = plt.subplots()
for i in [3, 6, 7]:
    uy = read_RFnodeFile('/home/biomech/Documents/01_Icotec/02_FEA/00_Model/94_screw_Osteoporosis_new_RFnode'
                         + str(i) + '.txt')

    for samples in [0, 3]:
        if i == 3:
            ax1.plot([4, 3, 2, 1, 0], np.append(uy[samples*4:samples*4+4], 0),
                     label=labels[samples], ls='-', marker='o', color=c[samples])
            ax2.plot([4, 3, 2, 1, 0], np.append(uy[samples*4:samples*4+4], 1) / -np.append(uy[0:4], 1)*100+100,
                     label=labels[samples], ls='-', marker='o', color=c[samples])
        elif i == 4:
            ax1.plot([4, 3, 2, 1, 0], np.append(uy[samples * 4:samples * 4 + 4], 0),
                     label='_nolegend_', ls='--', marker='o', color=c[samples])
            ax2.plot([4, 3, 2, 1, 0], np.append(uy[samples * 4:samples * 4 + 4], 1) / -np.append(uy[0:4], 1)*100+100,
                     label='_nolegend_', ls='--', marker='o', color=c[samples])
        elif i == 5:
            ax1.plot([4, 3, 2, 1, 0], np.append(uy[samples * 4:samples * 4 + 4], 0),
                     label='_nolegend_', ls=':', marker='o', color=c[samples])
            ax2.plot([4, 3, 2, 1, 0],
                     np.append(uy[samples * 4:samples * 4 + 4], 1) / -np.append(uy[0:4], 1) * 100 + 100,
                     label='_nolegend_', ls=':', marker='o', color=c[samples])
        elif i == 6:
            ax1.plot([4, 3, 2, 1, 0], np.append(uy[samples * 4:samples * 4 + 4], 0),
                     label='_nolegend_', ls='dashdot', marker='o', color=c[samples])
            ax2.plot([4, 3, 2, 1, 0],
                     np.append(uy[samples * 4:samples * 4 + 4], 1) / -np.append(uy[0:4], 1) * 100 + 100,
                     label='_nolegend_', ls='dashdot', marker='o', color=c[samples])
        elif i == 7:
            ax1.plot([4, 3, 2, 1, 0], np.append(uy[samples * 4:samples * 4 + 4], 0),
                     label='_nolegend_', ls='--', marker='o', color=c[samples])
            ax2.plot([4, 3, 2, 1, 0], np.append(uy[samples * 4:samples * 4 + 4], 1) / -np.append(uy[0:4], 1)*100+100,
                     label='_nolegend_', ls='--', marker='o', color=c[samples])
ax1.set_xlabel('RP')
plt.xticks([0, 1, 2, 3, 4])
ax1.set_ylabel('Displacement / mm')
ax1.legend()
ax2.set_xlabel('RP')
plt.xticks([0, 1, 2, 3, 4])
ax2.set_ylabel('Difference / %')
ax2.legend()
