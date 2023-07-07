import numpy as np
import matplotlib.pyplot as plt
import time
import pandas as pd
import pickle
#%%


def read_RFnodeFile(file_):
    # read data from text file
    df_ = np.loadtxt(file_, delimiter=',')

    # Contains the frame number:
    # frame = np.array(df_[:, 0])

    # reaction moment of the implant:
    # rmx_ = np.array(df_[:, 1])
    # rmy_ = np.array(df_[:, 2])
    # rmz_ = np.array(df_[:, 3])

    # reaction force of the implant:
    # rfx_ = np.array(df_[:, 4])
    rfy_ = np.array(df_[:, 5])
    # rfz_ = np.array(df_[:, 6])

    # transverse disp. of the implant:
    # ux_ = np.array(df_[:, 7])
    uy_ = np.array(df_[:, 8])
    # uz_ = np.array(df_[:, 9])

    # rotational disp. of the implant:
    # urx_ = np.array(df_[:, 10])
    # ury_ = np.array(df_[:, 11])
    # urz_ = np.array(df_[:, 12])

    return uy_, rfy_


def read_acumen(file_a):
    df_ = pd.read_csv(file_a, delimiter=';', skiprows=[0, 2])
    t_ = pd.DataFrame(df_, columns=['Time ']).to_numpy()
    d_ = pd.DataFrame(df_, columns=['Axial Displacement ']).to_numpy()
    # d_ = d_ - d_[0]  # calibrate displacement to zero at beginning
    f_ = pd.DataFrame(df_, columns=['Axial Force ']).to_numpy()
    # f_set_ = pd.DataFrame(df_, columns=['Axial Force Command ']).to_numpy()
    cycle_ = pd.DataFrame(df_, columns=['Axial Count ']).to_numpy()
    # arr_ = 0
    peak_ = np.zeros(int(np.max(cycle_)))
    vall_ = np.zeros(int(np.max(cycle_)))
    for j_ in range(2, int(np.max(cycle_))):
        # del arr_
        arr_ = np.where((cycle_ == j_) | (cycle_ == j_ + .5))[0]
        peak_[j_] = arr_[int(np.argmin(f_[arr_]))]
        vall_[j_] = arr_[int(np.argmax(f_[arr_]))]

    peak_ = peak_.astype(int)
    vall_ = vall_.astype(int)

    return cycle_, d_, f_, peak_, vall_, t_


def read_ARAMIS(file_A):
    df_ = pd.read_csv(file_A, delimiter=';', skiprows=[0])
    stop = np.where(df_.isnull())[0][0]
    lx = pd.DataFrame(df_, columns=['RodCsys→BoneCsys.LX [mm]'])
    ly = pd.DataFrame(df_, columns=['RodCsys→BoneCsys.LY [mm]'])
    lz = pd.DataFrame(df_, columns=['RodCsys→BoneCsys.LZ [mm]'])
    phiX = pd.DataFrame(df_, columns=['RodCsys→BoneCsys.Phi(X) [°]'])
    thetaY = pd.DataFrame(df_, columns=['RodCsys→BoneCsys.Theta(Y) [°]'])
    psiZ = pd.DataFrame(df_, columns=['RodCsys→BoneCsys.Psi(Z) [°]'])
    t_ = pd.DataFrame(df_, columns=['Time UTC'])

    return lx[:stop], ly[:stop], lz[:stop], phiX[:stop], thetaY[:stop], psiZ[:stop], t_[:stop]


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx


def find_first(array, value):
    array = np.asarray(array)
    idx = next(xd for xd, val in enumerate(array)
               if val <= value)
    return idx


def read_resample(file_r):
    df_ = pd.read_csv(file_r, delimiter=',')
    A_x_ = pd.DataFrame(df_, columns=['Aramis X']).to_numpy()
    A_y_ = pd.DataFrame(df_, columns=['Aramis Y']).to_numpy()
    A_z_ = pd.DataFrame(df_, columns=['Aramis Z']).to_numpy()
    A_rx_ = pd.DataFrame(df_, columns=['Aramis rX']).to_numpy()
    A_ry_ = pd.DataFrame(df_, columns=['Aramis rY']).to_numpy()
    A_rz_ = pd.DataFrame(df_, columns=['Aramis rZ']).to_numpy()
    a_y_ = pd.DataFrame(df_, columns=['Acumen Y']).to_numpy()
    a_f_ = pd.DataFrame(df_, columns=['Acumen Fy']).to_numpy()
    a_c_ = pd.DataFrame(df_, columns=['Acumen C']).to_numpy()

    return A_x_, A_y_, A_z_, A_rx_, A_ry_, A_rz_, a_y_, a_f_, a_c_


def read_exp_peaks():
    with open('/home/biomech/Downloads/mergedDf.pkl', 'rb') as f:
        data = pickle.load(f)
    return data


t1 = time.time()
# plt.close('all')


def read_FE_exp(number, model_code, plot):
    # Locations
    specimens = open('/home/biomech/Documents/01_Icotec/Specimens.txt', 'r').read().splitlines()  # Read specimens
    loc_Exp = '/home/biomech/Documents/01_Icotec/01_Experiments/00_Data/01_MainStudy/'  # location experimental results
    loc_FEA = '/home/biomech/Documents/01_Icotec/02_FEA/01_MainStudy/'  # location of fea results
    if number in [0, 2, 5, 7, 8, 10, 13, 15, 16, 18, 21, 23, 24, 26, 29, 31, 32]:
        model_code = model_code[:21] + 'P'
    elif number in [1, 3, 4, 6, 9, 11, 12, 14, 17, 19, 20, 22, 25, 27, 28, 30, 33]:
        model_code = model_code[:21] + 'T'
    else:
        print('Invalid model code!')
    specimen = specimens[number]

    file = [loc_Exp + specimen + '_resample.csv',
            loc_FEA + specimen + '/' + model_code[:14] + '/' + model_code + '_RFnode.txt',
            loc_FEA + specimen + '/' + model_code[:14] + '/' + model_code + '_RFnodeFix.txt']  # here

    # Load data
    [A_x, A_y, A_z, A_rx, A_ry, A_rz, a_y, a_f, a_c] = read_resample(file[0])  # load experimental result file (csv)
    [Uy, _] = read_RFnodeFile(file[1])  # read y displacement of moving reference node
    [_, RFy] = read_RFnodeFile(file[2])  # read y reaction force of fixed reference node
    if plot:
        D_F_plot(specimen, number, model_code, Uy, RFy, A_y, a_y, a_f)
    return Uy, RFy, A_y, a_y, a_f


def D_F_plot(specimen, number, model_code, Uy, RFy, A_y, a_y, a_f):
    # Figure
    fig, ax1 = plt.subplots(1, 1, figsize=(9, 6))  # set figure size
    plt.title('Experimental results ' + specimen + ' ' + model_code.split('_')[-1])
    if model_code.split('_')[-1] == 'P':
        ax1.plot(Uy, -RFy+RFy[0], label='FEA (YM PEEK = 25 GPa)', color='#1f77b4')  # plot fea results
    else:
        ax1.plot(Uy, -RFy+RFy[0], label='FEA (YM Ti = 100 GPa)', color='#1f77b4')  # plot fea results
    if number in [0, 1, 2, 10, 14, 17]:
        print('Using Acumen displacement.')
        ax1.plot(a_y, a_f - a_f[0], '--', label='Experiment', color='#ff7f0e')  # plot experimental results
    else:
        ax1.plot(A_y, a_f - a_f[0], label='Experiment', color='#ff7f0e')  # plot experimental results
        # ax1.plot(a_y, a_f - a_f[0], '--', label='_nolegend_', color='#ff7f0e')  # plot y data from acumen

    ax1.legend()
    ax1.set_xlabel('Displacement / mm')
    ax1.set_ylabel('Force / N')


def Peak_exp(number, ampl):
    d = read_exp_peaks()
    level = 2 ** (ampl - 2)
    peakF = d['MaxForce'][(number - 2) * 7 + ampl]
    print(d['DisplacementLevel'][(number - 2) * 7 + ampl])
    print(d['Specimen'][(number - 2) * 7 + ampl])
    print(peakF)
    return peakF


print('Execution time: ' + str(int((time.time()-t1)/60)) + ' min '+str(round(np.mod(time.time()-t1, 60), 1)) + ' sec.')



#%%


x = 0  # 0 = 0.25 mm, 1 = 0.5 mm, 2 = 1 mm, 3 = 2 mm, 4 = 4 mm, 5 = 8 mm, 6 = 16 mm
model = '82_L50_S50_D45_d1_05_T'

# peak_FE

fig, axs = plt.subplots(1, 1)
for x in range(2, 6):
    for i in range(34):
        try:
            [_, RFy_FE, _, _, _] = read_FE_exp(i, model, 0)
        except:
            continue
        try:
            RFy_FE = RFy_FE[x*21+10]
        except:
            continue
        RFy_exp = Peak_exp(i, x)
        plt.scatter(RFy_exp, RFy_FE, color='k')
        del RFy_FE
plt.plot([0, 285], [0, 285], 'k--')
axs.set_xlabel('Experiment')
axs.set_ylabel('FE')
# axs.set_title('Peaks at ' + str(2**(x-2)))
axs.set_aspect('equal')
