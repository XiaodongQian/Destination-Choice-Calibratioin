import csv
import numpy as np
import math
import codecs
import os
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

def f_Ttij_Dij(Tt_ij, D_ij,Station_index):
    den_1 = float(0)
    for origin in range(len(Station_index)):
        for dest in range(len(Station_index)):
            if origin != dest:
                den_1 += Tt_ij[origin][dest] * np.log10(D_ij[origin][dest])
    return den_1


def g_Ttij_Aij(Tt_ij, A_ij,Station_index):
    den_1 = float(0)
    for origin in range(len(Station_index)):
        for dest in range(len(Station_index)):
            if origin != dest:
                den_1 += Tt_ij[origin][dest] * np.log10(A_ij[origin][dest])
    return den_1


def f_value_0(T_ij, Tt_ij, D_ij, Station_index):
    den_1 = 0
    for origin in range(len(Station_index)):
        for dest in range(len(Station_index)):
            if origin != dest:
                den_1 += T_ij[origin][dest]*np.log10(D_ij[origin][dest])

    den_2 = 0
    for origin in range(len(Station_index)):
        for dest in range(len(Station_index)):
            if origin != dest:
                den_2 += Tt_ij[origin][dest]*np.log10(D_ij[origin][dest])


    f_value = den_1-den_2

    return f_value


def g_value_0(T_ij, Tt_ij, A_ij,Station_index):
    den_1 = float(0)
    for origin in range(len(Station_index)):
        for dest in range(len(Station_index)):
            if origin != dest:
                den_1 += T_ij[origin][dest] * np.log10(A_ij[origin][dest])

    den_2 = float(0)
    for origin in range(len(Station_index)):
        for dest in range(len(Station_index)):
            if origin != dest:
                den_2 += Tt_ij[origin][dest] * np.log10(A_ij[origin][dest])


    g_value = den_1 - den_2

    return g_value


# read the station data
#path = '/Volumes/GoogleDrive/My Drive/Dissertation/Data/Chicago/OD'
path = 'G:\My Drive\Dissertation\Data\Chicago\OD'
## build a staiton index dictionary
Station_index = {}
with open(os.path.join(path, 'Chi_BS_Station_2015.csv'), 'r') as f:
    reader = csv.reader(f)
    reader.next()
    dim = 0
    for row in reader:
        time = row[5]
        year = int(time.split('/')[2])
        if year < 2016:
            Station_index[int(row[0])] = dim
            dim += 1
del f
del reader


# read OD matrix
## calculate the sum(OD_ij)
with codecs.open(os.path.join(path, 'Chi_Bikeshare_OD_2016.csv'),'r',encoding="utf-8-sig") as f:
    reader = csv.reader(f)
    x = list(reader)
    Tt_ij = np.array(x).astype("float")
for origin in range(len(Station_index)):
    Tt_ij[origin][origin] = 0
## calculate the Ot_i and Dt_j
Ot_i = np.sum(Tt_ij, axis=1)
Dt_j = np.sum(Tt_ij, axis=0)

# read distance matrix
with codecs.open(os.path.join(path,'Chi_BS_Dist_Matrix_0_calibration_simple.csv'),'r',encoding="utf-8-sig") as f:
    reader = csv.reader(f)
    x = list(reader)
    D_ij = np.array(x).astype("float")
# read price matrix and we compare the parameter

# calcualte the destination within reach matrix
## we build a connection matrix A_jk show if k is within 500 meters of j
A_jk = np.zeros((len(Station_index), len(Station_index)))
for i in range(len(Station_index)):
    for j in range(len(Station_index)):
        if i!= j:
            if D_ij[i][j]<= 3000:
                A_jk[i][j] = 1

## calculate A_ij
rho_0 = -1
A_ij = np.zeros((len(Station_index), len(Station_index)))
for origin in range(len(Station_index)):
    for dest in range(len(Station_index)):
        if origin != dest:
            a_ij = 0
            for mid in range(len(Station_index)):
                if mid != origin and mid != dest:
                    a_ij += A_jk[dest][mid] * Dt_j[mid] * (D_ij[dest][mid])**(rho_0)
            A_ij[origin][dest] = a_ij

f_Ttij_Dij_value = f_Ttij_Dij(Tt_ij, D_ij,Station_index)
g_Ttij_Aij_value = g_Ttij_Aij(Tt_ij, A_ij,Station_index)

Z_i = np.zeros((len(Station_index),1))
T_ij = np.zeros((len(Station_index),len(Station_index)))



Beta = np.arange(-5,0,0.1)
Sigma = np.arange(-5,5,0.2)
Sigma,Beta = np.meshgrid(Sigma,Beta)
f_value = np.zeros((np.shape(Beta)[0], np.shape(Beta)[1]))
g_value = np.zeros((np.shape(Sigma)[0], np.shape(Sigma)[1]))
for i in range(np.shape(Beta)[0]):
    for j in range(np.shape(Sigma)[1]):
        print i,j
        beta_0 = np.asarray([Beta[i][j]] * len(Station_index), dtype=np.float64)
        sigma_0 = np.asarray([Sigma[i][j]] * len(Station_index), dtype=np.float64)
        # calculate the Z_i
        for origin in range(len(Station_index)):
            z_i = 0
            for dest in range(len(Station_index)):
                if dest != origin:
                    z_i += Dt_j[dest] * math.pow(A_ij[origin][dest], sigma_0[origin]) * math.pow(D_ij[origin][dest], beta_0[origin])
            Z_i[origin] = math.pow(z_i, -1)

        # calculate T_ij
        for origin in range(len(Station_index)):
            for dest in range(len(Station_index)):
                if origin != dest:
                    T_ij[origin][dest] = Z_i[origin] * Ot_i[origin] * Dt_j[dest] * math.pow(A_ij[origin][dest],sigma_0[origin]) * math.pow(D_ij[origin][dest], beta_0[origin])

        # calculate f and g
        f_1 = f_value_0(T_ij, Tt_ij, D_ij, Station_index)
        f_value[i][j] = f_1
        g_1 = g_value_0(T_ij, Tt_ij, A_ij, Station_index)
        g_value[i][j] = g_1


fig = plt.figure(figsize = (9,6))
ax = plt.subplot(1,1,1,projection = '3d')
surf = ax.plot_surface(Sigma, Beta, g_value, cmap=cm.coolwarm, linewidth=0, antialiased=False)
ax.set_xlabel('Sigma')
ax.set_ylabel('Beta')
ax.set_zlabel('G_value')
ax.view_init(20,100)
plt.colorbar(surf, shrink=0.5, aspect=5)#标注
plt.show()

# subplot 2 data line
ax2 = plt.subplot(1,2,2,projection = '3d')
ax2.plot(Sigma[25,:],Beta[25,:],f_value[25,:]) # fix beta
ax2.plot(Sigma[:,25],Beta[:,25],f_value[:,25]) # fix sigma
plt.show()

fig = plt.figure(figsize = (9,6))
ax1 = plt.subplot(1,2,1,projection = '3d')
ax1.plot(Sigma[25,:],Beta[25,:],g_value[25,:]) # fix beta
ax1.view_init(20,100)
ax2 = plt.subplot(1,2,2,projection = '3d')
ax2.plot(Sigma[:,30],Beta[:,30],g_value[:,30]) # fix sigma
ax2.view_init(20,100)
plt.show()

Data_Point = []
for i in range(np.shape(Beta)[0]):
    for j in range(np.shape(Sigma)[1]):
        Data_Point.append({'Beta':Beta[i][j],'Sigma':Sigma[i][j],'G_value':g_value[i][j]})
F_file = open(os.path.join(path,'Chi_CDmodel_Gvalue_matrix.csv'),'wb')
with F_file:
    name = ['Beta','Sigma','G_value']
    writer = csv.DictWriter(F_file, fieldnames=name)
    writer.writeheader()
    writer.writerows(Data_Point)