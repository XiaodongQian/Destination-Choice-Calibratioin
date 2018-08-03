import csv
import numpy as np
import math
import codecs
import os
from numpy import matrix
from numpy import linalg as LA

def FG_value(T_ij, Tt_ij, D_ij, A_ij, origin_station, dest_station):
    den_1 = 0
    for origin in range(len(origin_station)):
        for dest in range(len(dest_station)):
            if origin_station[origin] != dest_station[dest]:
                den_1 += T_ij[origin][dest] * np.log10(D_ij[origin][dest])

    den_2 = sum(sum(T_ij))

    den_3 = 0
    for origin in range(len(origin_station)):
        for dest in range(len(dest_station)):
            if origin_station[origin] != dest_station[dest]:
                den_3 += Tt_ij[origin][dest] * np.log10(D_ij[origin][dest])

    den_4 = sum(sum(Tt_ij))

    f_value = den_1/den_4-den_3/den_4

    den_5 = 0
    for origin in range(len(origin_station)):
        for dest in range(len(dest_station)):
            if origin_station[origin] != dest_station[dest]:
                den_5 += T_ij[origin][dest] * np.log10(A_ij[origin][dest])

    den_6 = 0
    for origin in range(len(origin_station)):
        for dest in range(len(dest_station)):
            if origin_station[origin] != dest_station[dest]:
                den_6 += Tt_ij[origin][dest] * np.log10(A_ij[origin][dest])

    g_value = den_5/den_4-den_6/den_4

    return [f_value,g_value]


def Z_i_value(Dt_j,A_ij,sigma,D_ij,beta,origin_station,dest_station):
    Z_i = np.zeros((len(origin_station), 1))
    for origin in range(len(origin_station)):
        z_i = 0
        for dest in range(len(dest_station)):
            if origin_station[origin] != dest_station[dest]:
                z_i += Dt_j[dest] * math.pow(A_ij[origin][dest], sigma) * math.pow(D_ij[origin][dest], beta)
        Z_i[origin] = math.pow(z_i, -1)
    return Z_i


def T_ij_value(Z_i,Ot_i,Dt_j,A_ij,sigma,D_ij,beta,origin_station, dest_station):
    T_ij = np.zeros((len(origin_station), len(dest_station)))
    for origin in range(len(origin_station)):
        for dest in range(len(dest_station)):
            if origin_station[origin] != dest_station[dest]:
                T_ij[origin][dest] = (Z_i[origin] * Ot_i[origin] * Dt_j[dest] * math.pow(A_ij[origin][dest], sigma) * math.pow(D_ij[origin][dest], beta))
    return T_ij


def Jacobi_Matrix(T_ij,Tt_ij,D_ij,A_ij,origin_station, dest_station):
    jacobi_matrix = np.zeros((2,2), dtype = np.float32)
    # the partial deviation
    den_1 = 0
    for origin in range(len(origin_station)):
        for dest in range(len(dest_station)):
            if origin_station[origin] != dest_station[dest]:
                den_1 += T_ij[origin][dest] * (np.log10(D_ij[origin][dest]))**2
    '''
    den_2 = 0
    for origin in range(len(origin_station)):
        for dest in range(len(dest_station)):
            if origin_station[origin] != dest_station[dest]:
                den_2 += T_ij[origin][dest] * np.log10(D_ij[origin][dest])
    '''
    den_3 = 0
    for origin in range(len(origin_station)):
        for dest in range(len(dest_station)):
            if origin_station[origin] != dest_station[dest]:
                den_3 += T_ij[origin][dest] * np.log10(D_ij[origin][dest]) * np.log10(A_ij[origin][dest])
    '''
    den_4 = 0
    for origin in range(len(origin_station)):
        for dest in range(len(dest_station)):
            if origin_station[origin] != dest_station[dest]:
                den_4 += T_ij[origin][dest] * np.log10(A_ij[origin][dest])
    '''
    den_5 = 0
    for origin in range(len(origin_station)):
        for dest in range(len(dest_station)):
            if origin_station[origin] != dest_station[dest]:
                den_5 += T_ij[origin][dest] * (np.log10(A_ij[origin][dest])) ** 2

    den_6 = sum(sum(Tt_ij))

    jacobi_matrix[0][0] = den_1/den_6
    jacobi_matrix[0][1] = den_3/den_6
    jacobi_matrix[1][0] = jacobi_matrix[0][1]
    jacobi_matrix[1][1] = den_5/den_6

    return jacobi_matrix


def A_ij_matrix(Dt_j,D_ij,Station_index,distance,rho):
    # calcualte the destination within reach matrix
    ## we build a connection matrix A_jk show if k is within distance meters of j
    iter = len(Station_index)
    A_jk = np.zeros((iter, iter))
    for i in range(iter):
        for j in range(iter):
            if i != j:
                if D_ij[i][j] <= distance:
                    A_jk[i][j] = 1

    ## calculate A_ij
    A_ij = np.zeros((iter,iter))
    for origin in range(iter):
        for dest in range(iter):
            if origin != dest:
                a_ij = 0
                for mid in range(iter):
                    if mid != origin and mid != dest:
                        a_ij += A_jk[dest][mid] * Dt_j[mid] * (D_ij[dest][mid]) ** (rho)
                A_ij[origin][dest] = a_ij
    return A_ij


def OD_Diff(T_ij,Tt_ij,origin_station,dest_station):
    Diff = np.zeros((len(origin_station), len(dest_station)))
    diff = 0
    for origin in range(len(origin_station)):
        for dest in range(len(dest_station)):
            if origin_station[origin] != dest_station[dest]:
                if Tt_ij[origin][dest] != 0:
                    diff = float(int(T_ij[origin][dest]) - int(Tt_ij[origin][dest])) / Tt_ij[origin][dest]
                else:
                    if T_ij[origin][dest] != 0:
                        diff = np.inf
            Diff[origin][dest] = diff
    i = 0
    diff_value = 0
    for origin in range(len(origin_station)):
        for dest in range(len(dest_station)):
            if Diff[origin][dest] != np.inf:
                diff_value += Diff[origin][dest]
                i += 1
    return diff_value / i


#### Part 1 Data compile
## read the station data
path = '/Volumes/GoogleDrive/My Drive/Dissertation/Data/Chicago'
#path = 'G:\My Drive\Dissertation\Data\Chicago\OD'
## build a staiton index dictionary
Station_index = {}
with open(os.path.join(path, 'OD', 'Chi_BS_Station_2015.csv'), 'r') as f:
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

## read the stations in disadvantaged communities and other areas
dis_station = []
other_station = []
with open(os.path.join(path, 'Station', 'Chi_Bikeshare_Ridership_Regression.csv'), 'r') as f:
    reader = csv.reader(f)
    reader.next()
    for row in reader:
        type = row[58]
        station_id = int(row[2])
        if type == 'Disadvantaged areas':
            dis_station.append(station_id)
        else:
            other_station.append(station_id)
del f
del reader

dis_station.sort()
other_station.sort()

origin_station = other_station
dest_station = other_station

## read OD matrix
with codecs.open(os.path.join(path, 'OD', 'Chi_Bikeshare_OD_2016.csv'),'r',encoding="utf-8-sig") as f:
    reader = csv.reader(f)
    x = list(reader)
    Tt_ij = np.array(x).astype("float")
for origin in range(len(Station_index)):
    Tt_ij[origin][origin] = 0

Ot_i = np.sum(Tt_ij, axis=1)
Dt_j = np.sum(Tt_ij, axis=0)
## create a new matrix based on categories
Tt_ij_origin_dest = np.zeros((len(origin_station),len(dest_station)))
for origin in origin_station:
    for dest in dest_station:
        Tt_ij_origin_dest[origin_station.index(origin)][dest_station.index(dest)] = Tt_ij[Station_index[origin]][Station_index[dest]]
Ot_i_origin = np.sum(Tt_ij_origin_dest, axis=1)
Dt_j_dest = np.sum(Tt_ij_origin_dest, axis=0)

## read distance matrix
with codecs.open(os.path.join(path, 'OD', 'Chi_BS_Dist_Matrix_0_calibration_simple.csv'),'r',encoding="utf-8-sig") as f:
    reader = csv.reader(f)
    x = list(reader)
    D_ij = np.array(x).astype("float")
## create a new distance based on categories
D_ij_origin_dest = np.zeros((len(origin_station),len(dest_station)))
for origin in origin_station:
    for dest in dest_station:
        D_ij_origin_dest[origin_station.index(origin)][dest_station.index(dest)] = D_ij[Station_index[origin]][Station_index[dest]]
# read price matrix and we compare the parameter

## Calculate the A_ij Matrix
rho = -1
distance = 3000
A_ij = A_ij_matrix(Dt_j,D_ij,Station_index,distance,rho)
## create a new A_ij based on categories
A_ij_origin_dest = np.zeros((len(origin_station),len(dest_station)))
for origin in origin_station:
    for dest in dest_station:
        A_ij_origin_dest[origin_station.index(origin)][dest_station.index(dest)] = A_ij[Station_index[origin]][Station_index[dest]]


##############
#  calibration process
## initial value
beta_0 = -4
sigma_0 = 1
stepsize = 0.1 # (0,1)
iteration = 1
threshold = 0.00000001

## keep the value change for every iteration
fg_value_store = []
beta_value_store = [beta_0]
sigma_value_store = [sigma_0]
while iteration < 10000:
    print '############'
    print 'This is the ' + str(iteration) + 'th iteration'
    Z_i = Z_i_value(Dt_j_dest, A_ij_origin_dest, sigma_0, D_ij_origin_dest, beta_0, origin_station, dest_station)
    T_ij = T_ij_value(Z_i, Ot_i_origin, Dt_j_dest, A_ij_origin_dest, sigma_0, D_ij_origin_dest, beta_0, origin_station, dest_station)
    jacobi_matrix = Jacobi_Matrix(T_ij, D_ij_origin_dest, A_ij_origin_dest, origin_station, dest_station)
    jacobi_matrix = matrix(jacobi_matrix)
    var_0 = matrix([[beta_0],[sigma_0]])
    fg_value = FG_value(T_ij, Tt_ij_origin_dest, D_ij_origin_dest, A_ij_origin_dest, origin_station, dest_station)
    print fg_value
    fg_value_store.append(fg_value)
    LA_fg = LA.norm(fg_value)
    if LA_fg > threshold:
        direction = - jacobi_matrix.I * matrix(fg_value).T
        var_1 = var_0 + stepsize * direction
        beta_0 = float(var_1[0][0])
        sigma_0 = float(var_1[1][0])
        beta_value_store.append(beta_0)
        sigma_value_store.append(sigma_0)
        iteration += 1
    else:
        print 'done'
        break

# calculate the overall goodness of fit
Ave_diff = OD_Diff(T_ij,Tt_ij_origin_dest,origin_station,dest_station)


        
