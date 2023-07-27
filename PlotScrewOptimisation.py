import matplotlib.pyplot as plt
import numpy as np


def read_RFnodeFile(file_):
    # read data from text file
    df_ = np.loadtxt(file_, delimiter=',')
    uy_ = np.array(df_[:, 1])

    return uy_


plt.close('all')
labels = ['Original', 'Simple cannulated', 'Simple filled']

c = ['r', 'g', 'b', 'k']
fig1, ax1 = plt.subplots()
fig2, ax2 = plt.subplots()
rp = [50, 38.5, 26.5, 21.5, 9.5, 0]
for i in [0, 4]:
    uy = read_RFnodeFile('/home/biomech/Documents/01_Icotec/02_FEA/00_Model/95_screw_DPS_RFnode'
                         + str(i) + '.txt')
    if i == 0:
        sam = [0]
    else:
        sam = [1, 2]
    for samples in sam:
        print(uy[samples*5:samples*5+5])
        if i == 0:
            ax1.plot(rp, np.append(uy[samples*5:samples*5+5], 0),
                     label=labels[samples], ls='-', marker='o', color=c[samples])
            ax2.plot(rp, np.append(uy[samples*5:samples*5+5], 1) / -np.append(uy[0:5], 1)*100+100,
                     label=labels[samples], ls='-', marker='o', color=c[samples])
        else:
            ax1.plot(rp, np.append(uy[samples*5:samples*5+5], 0),
                     label=labels[samples], ls='--', marker='o', color=c[samples])
            ax2.plot(rp, np.append(uy[samples*5:samples*5+5], 1) / -np.append(uy[0:5], 1)*100+100,
                     label=labels[samples], ls='--', marker='o', color=c[samples])
ax1.set_xlabel('RP')
ax1.set_xticks(rp)
ax1.set_ylabel('Displacement / mm')
ax1.legend()
ax2.set_xlabel('RP')
ax2.set_xticks(rp)
ax2.set_ylabel('Difference / %')
ax2.legend()
ax2.set_ylim([-10, 10])
