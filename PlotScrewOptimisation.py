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
for i in [10]:
    uy = read_RFnodeFile('/home/biomech/Documents/01_Icotec/02_FEA/00_Model/95_screw_DPS_RFnode'
                         + str(i) + '.txt')/4
    if i == 10:
        sam = [0, 2, 3]
    else:
        sam = [3]
    for samples in sam:
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
                     label=labels[samples], ls='--', marker='o', color=c[samples])
            ax2.plot([4, 3, 2, 1, 0], np.append(uy[samples * 4:samples * 4 + 4], 1) / -np.append(uy[0:4], 1)*100+100,
                     label=labels[samples], ls='--', marker='o', color=c[samples])
        elif i == 10:
            ax1.plot([4, 3, 2, 1, 0], np.append(uy[samples * 4:samples * 4 + 4], 0),
                     label=labels[samples], ls='--', marker='o', color=c[samples])
            ax2.plot([4, 3, 2, 1, 0], np.append(uy[samples * 4:samples * 4 + 4], 1) / -np.append(uy[0:4], 1)*100+100,
                     label=labels[samples], ls='--', marker='o', color=c[samples])
ax1.set_xlabel('RP')
plt.xticks([0, 1, 2, 3, 4])
ax1.set_ylabel('Displacement / mm')
ax1.legend()
ax2.set_xlabel('RP')
ax2.set_xticks([0, 1, 2, 3, 4])
plt.xticks([0, 1, 2, 3, 4])
ax2.set_ylabel('Difference / %')
ax2.legend()
ax2.set_ylim([-5, 5])
