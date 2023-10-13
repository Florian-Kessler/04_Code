# import matplotlib.pyplot as plt
# import numpy as np
#
#
# def read_RFnodeFile_opt(file_):
##     read data from text file
    # df_ = np.loadtxt(file_, delimiter=',')
    # uy_ = np.array(df_[:, 1])
    #
    # return uy_
#
#
# plt.close('all')
# labels = ['Original', 'Simple cannulated', 'Simple filled']
#
# mat = 'peek'
# c = ['r', 'b', 'g', 'k']
# fig1, ax1 = plt.subplots()
# fig2, ax2 = plt.subplots()
# rp = [50, 38.5, 26.5, 21.5, 9.5, 0]  # DPS
# rp = [50.0, 36.4, 23.5, 14.4, 0]  # PEEK
# for i in [0, 4]:  # DPS
# for i in [3, 6, 7]:  # PEEK
#     uy = read_RFnodeFile_opt('/home/biomech/Documents/01_Icotec/02_FEA/00_Model/95_screw_DPS_RFnode'
#                              + str(i) + '.txt')  # DPS
    # uy = read_RFnodeFile_opt('/home/biomech/Documents/01_Icotec/02_FEA/00_Model/'
    #                          '94_screw_Osteoporosis_new_RFnode'
    #                          + str(i) + '.txt')  # PEEK
    # if mat == 'ti':
    #     if i == 0:
    #         sam = [0]
    #     else:
    #         sam = [1]
    #     for samples in sam:
    #         print(uy[samples * 5:samples * 5 + 5])
    #         if i == 0:
    #             ax1.plot(rp, np.append(uy[samples * 5:samples * 5 + 5], 0),
    #                      label=labels[samples], ls='-', marker='o', color=c[samples])
    #             ax2.plot(rp, np.append(uy[samples * 5:samples * 5 + 5], 1) / -np.append(uy[0:5], 1) * 100 + 100,
    #                      label=labels[samples], ls='-', marker='o', color=c[samples])
    #         else:
    #             ax1.plot(rp, np.append(uy[samples * 5:samples * 5 + 5], 0),
    #                      label=labels[samples], ls='--', marker='o', color=c[samples])
    #             ax2.plot(rp, np.append(uy[samples * 5:samples * 5 + 5], 1) / -np.append(uy[0:5], 1) * 100 + 100,
    #                      label=labels[samples], ls='--', marker='o', color=c[samples])
    # if mat == 'peek':
    #     for samples in sam:
    #         if i == 3:
    #             ax1.plot([4, 3, 2, 1, 0], np.append(uy[samples * 4:samples * 4 + 4], 0),
    #                      label=labels[samples], ls='-', marker='o', color=c[samples])
    #             ax2.plot([4, 3, 2, 1, 0], np.append(uy[samples * 4:samples * 4 + 4], 1) / -np.append(uy[0:4], 1) * 100 + 100,
    #                      label=labels[samples], ls='-', marker='o', color=c[samples])
    #         else:
    #             ax1.plot([4, 3, 2, 1, 0], np.append(uy[samples * 4:samples * 4 + 4], 0),
    #                      label='_nolegend_', ls='--', marker='o', color=c[samples])
    #             ax2.plot([4, 3, 2, 1, 0], np.append(uy[samples * 4:samples * 4 + 4], 1) / -np.append(uy[0:4], 1)*100+100,
    #                      label='_nolegend_', ls='--', marker='o', color=c[samples])
#
# ax1.set_xlabel('RP')
# ax1.set_xticks(rp)
# ax1.set_ylabel('Displacement / mm')
# ax1.legend()
# ax2.set_xlabel('RP')
# ax2.set_xticks(rp)
# ax2.set_ylabel('Difference / %')
# ax2.legend()
# ax2.set_ylim([-5, 5])

#%%

import matplotlib.pyplot as plt
import matplotlib
import numpy as np
fs = 22

def read_RFnodeFile(file_):
    # read data from text file
    df_ = np.loadtxt(file_, delimiter=',')
    uy_ = np.array(df_[:, 1])

    return uy_


plt.close('all')
labels = ['Original', 'Simple cannulated', 'Simple filled']

c = ['r', 'b', 'b', 'k']
fig1, ax1 = plt.subplots()
fig1.set_figheight(5)
fig1.set_figwidth(7)
fig2, ax2 = plt.subplots()
fig2.set_figheight(5)
fig2.set_figwidth(7)

rp = [50, 38.5, 26.5, 21.5, 9.5, 0]
for i in [0, 4]:
    uy = read_RFnodeFile('/home/biomech/Documents/01_Icotec/02_FEA/00_Model/95_screw_DPS_RFnode'
                         + str(i) + '.txt')
    if i == 0:
        sam = [0]
    else:
        sam = [1]
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
# ax1.set_xlabel('RP', fontsize=fs)
# ax1.set_xticks(rp)
ax1.set_ylabel('Displacement / mm', fontsize=fs)
ax1.legend(fontsize=fs)
ax1.tick_params(axis='both', which='major', labelsize=fs)
ax1.tick_params(axis='both', which='minor', labelsize=fs-2)
fig1.subplots_adjust(left=0.2)
ax1.tick_params(
     axis='x',           # changes apply to the x-axis
     which='both',       # both major and minor ticks are affected
     bottom=False,       # ticks along the bottom edge are off
     top=False,          # ticks along the top edge are off
     labelbottom=False)  # labels along the bottom edge are off
fig1.savefig('/home/biomech/Documents/GitHub/05_Report/02_Pictures_MM/Opt_DPS_abs.eps')
# ax2.set_xlabel('RP', fontsize=fs)
# ax2.set_xticks(rp)
ax2.set_ylabel('Difference / %', fontsize=fs)
# ax2.legend(fontsize=fs)
ax2.set_ylim([-5, 5])
ax2.tick_params(axis='both', which='major', labelsize=fs)
ax2.tick_params(axis='both', which='minor', labelsize=fs-2)
fig2.subplots_adjust(left=0.2)
ax2.tick_params(
     axis='x',           # changes apply to the x-axis
     which='both',       # both major and minor ticks are affected
     bottom=False,       # ticks along the bottom edge are off
     top=False,          # ticks along the top edge are off
     labelbottom=False)  # labels along the bottom edge are off
fig2.savefig('/home/biomech/Documents/GitHub/05_Report/02_Pictures_MM/Opt_DPS_rel.eps')
#%%

import matplotlib.pyplot as plt
import matplotlib
import numpy as np
def read_RFnodeFile(file_):
    # read data from text file
    df_ = np.loadtxt(file_, delimiter=',')
    uy_ = np.array(df_[:, 1])
    return uy_
labels = ['Original', 'Original filled', 'Simple filled', 'Simple cannulated']
# labels = ['1.5 mm', '0.5 mm', '0.2 mm', '0.2 mm, opt']
c = ['r', 'b', 'b', 'b']
fig3, ax3 = plt.subplots()
fig3.set_figheight(5)
fig3.set_figwidth(7)
fig4, ax4 = plt.subplots()
fig4.set_figheight(5)
fig4.set_figwidth(7)
for i in [3, 7]:
    uy = read_RFnodeFile('/home/biomech/Documents/01_Icotec/02_FEA/00_Model/94_screw_Osteoporosis_new_RFnode'
                         + str(i) + '.txt')
    if i == 3:
        sam = [0]
    else:
        sam = [3]
    for samples in sam:
        if i == 3:
            ax3.plot([4, 3, 2, 1, 0], np.append(uy[samples * 4:samples * 4 + 4], 0),
                     label=labels[samples], ls='-', marker='o', color=c[samples])
            ax4.plot([4, 3, 2, 1, 0],
                     np.append(uy[samples * 4:samples * 4 + 4], 1) / -np.append(uy[0:4], 1) * 100 + 100,
                     label=labels[samples], ls='-', marker='o', color=c[samples])
        elif i == 4:
            ax3.plot([4, 3, 2, 1, 0], np.append(uy[samples * 4:samples * 4 + 4], 0),
                     label='_nolegend_', ls='--', marker='o', color=c[samples])
            ax4.plot([4, 3, 2, 1, 0],
                     np.append(uy[samples * 4:samples * 4 + 4], 1) / -np.append(uy[0:4], 1) * 100 + 100,
                     label='_nolegend_', ls='--', marker='o', color=c[samples])
        elif i == 5:
            ax3.plot([4, 3, 2, 1, 0], np.append(uy[samples * 4:samples * 4 + 4], 0),
                     label='_nolegend_', ls=':', marker='o', color=c[samples])
            ax4.plot([4, 3, 2, 1, 0],
                     np.append(uy[samples * 4:samples * 4 + 4], 1) / -np.append(uy[0:4], 1) * 100 + 100,
                     label='_nolegend_', ls=':', marker='o', color=c[samples])
        elif i == 6:
            ax3.plot([4, 3, 2, 1, 0], np.append(uy[samples * 4:samples * 4 + 4], 0),
                     label='_nolegend_', ls='dashdot', marker='o', color=c[samples])
            ax4.plot([4, 3, 2, 1, 0],
                     np.append(uy[samples * 4:samples * 4 + 4], 1) / -np.append(uy[0:4], 1) * 100 + 100,
                     label='_nolegend_', ls='dashdot', marker='o', color=c[samples])
        elif i == 7:
            ax3.plot([4, 3, 2, 1, 0], np.append(uy[samples * 4:samples * 4 + 4], 0),
                     label=labels[samples], ls='--', marker='o', color=c[samples])
            ax4.plot([4, 3, 2, 1, 0], np.append(uy[samples * 4:samples * 4 + 4], 1) / -np.append(uy[0:4], 1)*100+100,
                     label=labels[samples], ls='--', marker='o', color=c[samples])
# ax3.set_xlabel('RP', fontsize=fs)
# ax3.set_xticks(rp)
ax3.set_ylabel('Displacement / mm', fontsize=fs)
ax3.legend(fontsize=fs)
ax3.tick_params(axis='both', which='major', labelsize=fs)
ax3.tick_params(axis='both', which='minor', labelsize=fs-2)
fig3.subplots_adjust(left=0.2)
ax3.tick_params(
     axis='x',           # changes apply to the x-axis
     which='both',       # both major and minor ticks are affected
     bottom=False,       # ticks along the bottom edge are off
     top=False,          # ticks along the top edge are off
     labelbottom=False)  # labels along the bottom edge are off
fig3.savefig('/home/biomech/Documents/GitHub/05_Report/02_Pictures_MM/Opt_PEEK_abs.eps')
# ax4.set_xlabel('RP', fontsize=fs)
# ax4.set_xticks(rp)
ax4.set_ylabel('Difference / %', fontsize=fs)
# ax4.legend(fontsize=fs)
ax4.set_ylim([-5, 5])
ax4.tick_params(axis='both', which='major', labelsize=fs)
ax4.tick_params(axis='both', which='minor', labelsize=fs-2)
fig4.subplots_adjust(left=0.2)
ax4.tick_params(
     axis='x',           # changes apply to the x-axis
     which='both',       # both major and minor ticks are affected
     bottom=False,       # ticks along the bottom edge are off
     top=False,          # ticks along the top edge are off
     labelbottom=False)  # labels along the bottom edge are off
fig4.savefig('/home/biomech/Documents/GitHub/05_Report/02_Pictures_MM/Opt_PEEK_rel.eps')
