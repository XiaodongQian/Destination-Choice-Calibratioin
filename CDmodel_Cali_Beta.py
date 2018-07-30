import csv
import numpy as np
import math
import codecs
import os

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


def f_pd_0(T_ij,D_ij,A_ij,origin):
    den_1 = 0
    for dest in range(len(Station_index)):
        if dest != origin:
            den_1 += T_ij[origin][dest] * (np.log10(D_ij[origin][dest]))**2

    f_partial_deviation_beta = den_1

    den_2 = 0
    for dest in range(len(Station_index)):
        if dest != origin:
            den_2 += T_ij[origin][dest] * (np.log10(D_ij[origin][dest])) * (np.log10(A_ij[origin][dest]))

    f_partial_deviation_sigma = den_2

    return [f_partial_deviation_beta,f_partial_deviation_sigma]

def g_pd_0(T_ij,A_ij,D_ij,origin):
    den_1 = 0
    for dest in range(len(Station_index)):
        if dest != origin:
            den_1 += T_ij[origin][dest] * (np.log10(A_ij[origin][dest])) ** 2

    g_partial_deviation_sigma = den_1

    den_2 = 0
    for dest in range(len(Station_index)):
        if dest != origin:
            den_2 += T_ij[origin][dest] * (np.log10(D_ij[origin][dest])) * (np.log10(A_ij[origin][dest]))

    g_partial_deviation_beta = den_2

    return [g_partial_deviation_beta, g_partial_deviation_sigma]


# read the station data
path = '/Volumes/GoogleDrive/My Drive/Dissertation/Data/Chicago/OD'
#path = 'G:\My Drive\Dissertation\Data\Chicago\OD'
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

# calibration process
## initial value
beta_0 = np.asarray(np.random.normal(-1.5, 0.5, len(Station_index)))
sigma_0 = np.asarray(np.random.normal(2.5, 1, len(Station_index)))
beta_1 = np.asarray([0.0] *len(Station_index), dtype=np.float64)
sigma_1 = np.asarray([0.0] *len(Station_index),dtype=np.float64)

f_0 = 20
g_0 = 20
stepsize_beta = 0.01
stepsize_sigma = 0.01
stepsize_beta_0 =0.01
stepsize_sigma_0 = 0.01
iteration = 1
threshold = 500
Z_i = np.zeros((len(Station_index),1))
T_ij = np.zeros((len(Station_index),len(Station_index)))

## keep the value change for every iteration
f_value_store = [f_0]
g_value_store = [g_0]
beta_value_store = beta_0
sigma_value_store = sigma_0
while iteration < 100:
    print '############'
    print 'This is the ' + str(iteration) + 'th iteration'
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
    ## we have interation through the i j because i!=j
    #f_1 = g_value_0(T_ij, Tt_ij, D_ij, Station_index)
    g_1 = g_value_0(T_ij, Tt_ij, A_ij, Station_index)
    print 'g_1 is: ' + str(g_1)
    if 30000 <= g_1 < np.inf:
        stepsize_beta = 1
        stepsize_sigma = 1
    elif 10000 <= g_1 < 30000:
        stepsize_beta = 0.1
        stepsize_sigma = 0.1
    else:
        stepsize_beta = 0.01
        stepsize_sigma = 0.01

    # decide the stepsize
    #if abs(f_0-f_1)/f_0

    #if abs(f_0 - f_1) > threshold and abs(g_0 - g_1) > threshold:
    if abs(g_1) > threshold:
        print 'we need to find new beta and sigma'
        # refine the beta and sigma
        for origin in range(len(Station_index)):
            beta_1[origin] = beta_0[origin] - g_pd_0(T_ij,D_ij,A_ij,origin)[0]/g_1 * stepsize_beta
            sigma_1[origin] = sigma_0[origin] - g_pd_0(T_ij,D_ij,A_ij,origin)[1]/g_1 * stepsize_sigma
        beta_0 = beta_1
        beta_value_store = np.vstack([beta_value_store, beta_0])
        #beta_value_store.append(beta_0)
        sigma_0 = sigma_1
        sigma_value_store = np.vstack([sigma_value_store, sigma_0])
        #sigma_value_store.append(sigma_0)
        g_0 = g_1
        g_value_store.append(g_0)
        iteration += 1
    else:
        print 'we need to reset the beta and sigma'
        break

# calculate the overall goodness of fit
## calculate the O_i and D_j

Diff = np.zeros((len(Station_index),len(Station_index)))
diff = 0
for origin in range(len(Station_index)):
    for dest in range(len(Station_index)):
        if origin != dest:
            if Tt_ij[origin][dest] != 0:
                diff = (T_ij[origin][dest] - Tt_ij[origin][dest])/Tt_ij[origin][dest]
            else:
                if T_ij[origin][dest] !=0:
                    diff = np.inf
        Diff[origin][dest] = diff
i = 0
diff_value = 0
for origin in range(len(Station_index)):
    for dest in range(len(Station_index)):
        if Diff[origin][dest] != np.inf:
            diff_value += Diff[origin][dest]
            i += 1
Ave_diff = diff_value/i