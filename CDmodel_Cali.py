import csv
import numpy as np
import math
import codecs
import os
from numpy import matrix

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

def f_value_1(T_ij, Tt_ij, D_ij, Station_index):
    den_1 = 0
    for origin in range(len(Station_index)):
        for dest in range(len(Station_index)):
            if origin != dest:
                den_1 += T_ij[origin][dest]*np.log10(D_ij[origin][dest])

    den_2 = 0
    for origin in range(len(Station_index)):
        for dest in range(len(Station_index)):
            if origin != dest:
                den_2 += T_ij[origin][dest]

    mol_1 = 0
    for origin in range(len(Station_index)):
        for dest in range(len(Station_index)):
            if origin != dest:
                mol_1 += Tt_ij[origin][dest]*np.log10(D_ij[origin][dest])

    mol_2 = 0
    for origin in range(len(Station_index)):
        for dest in range(len(Station_index)):
            if origin != dest:
                mol_2 += Tt_ij[origin][dest]

    f_value = den_1/mol_1-den_2/mol_2

    return f_value

def g_value_1(T_ij, Tt_ij, A_ij,Station_index):
    den_1 = float(0)
    for origin in range(len(Station_index)):
        for dest in range(len(Station_index)):
            if origin != dest:
                den_1 += T_ij[origin][dest] * np.log10(A_ij[origin][dest])

    den_2 = float(0)
    for origin in range(len(Station_index)):
        for dest in range(len(Station_index)):
            if origin != dest:
                den_2 += T_ij[origin][dest]

    mol_1 = float(0)
    for origin in range(len(Station_index)):
        for dest in range(len(Station_index)):
            if origin != dest:
                mol_1 += Tt_ij[origin][dest] * np.log10(A_ij[origin][dest])

    mol_2 = float(0)
    for origin in range(len(Station_index)):
        for dest in range(len(Station_index)):
            if origin != dest:
                mol_2 += Tt_ij[origin][dest]

    g_value = den_1 / mol_1 - den_2 / mol_2

    return g_value

def f_value(T_ij, Tt_ij, D_ij, Station_index):
    den_1 = 0
    for origin in range(len(Station_index)):
        for dest in range(len(Station_index)):
            if origin != dest:
                den_1 += T_ij[origin][dest]*np.log10(D_ij[origin][dest])

    mol_1 = 0
    for origin in range(len(Station_index)):
        for dest in range(len(Station_index)):
            if origin != dest:
                mol_1 += T_ij[origin][dest]

    den_2 = 0
    for origin in range(len(Station_index)):
        for dest in range(len(Station_index)):
            if origin != dest:
                den_2 += Tt_ij[origin][dest]*np.log10(D_ij[origin][dest])

    mol_2 = 0
    for origin in range(len(Station_index)):
        for dest in range(len(Station_index)):
            if origin != dest:
                mol_2 += Tt_ij[origin][dest]

    f_value = den_1/mol_1-den_2/mol_2

    return f_value

def g_value(T_ij, Tt_ij, A_ij,Station_index):
    den_1 = float(0)
    for origin in range(len(Station_index)):
        for dest in range(len(Station_index)):
            if origin != dest:
                den_1 += T_ij[origin][dest] * np.log10(A_ij[origin][dest])

    mol_1 = float(0)
    for origin in range(len(Station_index)):
        for dest in range(len(Station_index)):
            if origin != dest:
                mol_1 += T_ij[origin][dest]

    den_2 = float(0)
    for origin in range(len(Station_index)):
        for dest in range(len(Station_index)):
            if origin != dest:
                den_2 += Tt_ij[origin][dest] * np.log10(A_ij[origin][dest])

    mol_2 = float(0)
    for origin in range(len(Station_index)):
        for dest in range(len(Station_index)):
            if origin != dest:
                mol_2 += Tt_ij[origin][dest]

    g_value = den_1 / mol_1 - den_2 / mol_2

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

def g_pd_0(T_ij,A_ij,origin):
    den_1 = 0
    for dest in range(len(Station_index)):
        if dest != origin:
            den_1 += T_ij[origin][dest] * (np.log10(A_ij[origin][dest])) ** 2

    g_partial_deviation_sigma = den_1

    return g_partial_deviation_sigma

def f_pd(T_ij,D_ij,origin):
    den_1 = 0
    for dest in range(len(Station_index)):
        if dest != origin:
            den_1 += T_ij[origin][dest] * (np.log10(D_ij[origin][dest]))**2

    den_2 = 0
    for dest in range(len(Station_index)):
        if dest != origin:
            den_2 += T_ij[origin][dest]

    den_3 = 0
    for dest in range(len(Station_index)):
        if dest != origin:
            den_3 += T_ij[origin][dest] * np.log10(D_ij[origin][dest])

    mol = den_2**2

    f_partial_deviation = (den_1 * den_2 - den_3**2)/mol

    return f_partial_deviation

def g_pd(T_ij,A_ij,origin):
    den_1 = 0
    for dest in range(len(Station_index)):
        if dest != origin:
            den_1 += T_ij[origin][dest] * (np.log10(A_ij[origin][dest])) ** 2

    den_2 = 0
    for dest in range(len(Station_index)):
        if dest != origin:
            den_2 += T_ij[origin][dest]

    den_3 = 0
    for dest in range(len(Station_index)):
        if dest != origin:
            den_3 += T_ij[origin][dest] * np.log10(A_ij[origin][dest])

    mol = den_2**2

    g_partial_deviation = (den_1*den_2-den_3**2)/mol

    return g_partial_deviation

def Z_i_value(Dt_j,A_ij,sigma_0,D_ij,beta_0,Station_index):
    Z_i = np.zeros((len(Station_index), 1))
    for origin in range(len(Station_index)):
        z_i = 0
        for dest in range(len(Station_index)):
            if dest != origin:
                z_i += Dt_j[dest] * math.pow(A_ij[origin][dest], sigma_0[origin]) * math.pow(D_ij[origin][dest], beta_0[origin])
        Z_i[origin] = math.pow(z_i, -1)

    return Z_i


def T_ij_value(Z_i,Ot_i,Dt_j,A_ij,sigma_0,D_ij,beta_0,Station_index):
    T_ij = np.zeros((len(Station_index), len(Station_index)))
    for origin in range(len(Station_index)):
        for dest in range(len(Station_index)):
            if origin != dest:
                T_ij[origin][dest] = Z_i[origin] * Ot_i[origin] * Dt_j[dest] * math.pow(A_ij[origin][dest],sigma_0[origin]) * math.pow(D_ij[origin][dest], beta_0[origin])
    return T_ij

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
#beta_1 = np.asarray([0.0] *len(Station_index), dtype=np.float64)
#sigma_1 = np.asarray([0.0] *len(Station_index),dtype=np.float64)

f_0 = 20
f_0_mid = np.inf
g_0 = 20
g_0_mid = np.inf
stepsize = 0.01
iteration = 1
threshold = 100
Z_i = np.zeros((len(Station_index),1))
T_ij = np.zeros((len(Station_index),len(Station_index)))

## keep the value change for every iteration
f_value_store = [f_0]
g_value_store = [g_0]
beta_value_store = beta_0
sigma_value_store = sigma_0
while iteration < 1000:
    print '############'
    print 'This is the ' + str(iteration) + 'th iteration'
    # calculate the Z_i
    Z_i = Z_i_value(Dt_j,A_ij,sigma_0,D_ij,beta_0,Station_index)
    # calculate T_ij
    T_ij = T_ij_value(Z_i,Ot_i,Dt_j,A_ij,sigma_0,D_ij,beta_0,Station_index)
    # calculate f and g
    ## we have interation through the i j because i!=j
    f_1 = f_value_0(T_ij, Tt_ij, D_ij, Station_index)
    g_1 = g_value_0(T_ij, Tt_ij, A_ij, Station_index)
    print 'f_1 and g_1 are: '
    print f_1, g_1
    print 'f_0 - f_1 = ' + str(f_0 - f_1)
    print 'g_0 - g_1 = ' + str(g_0 - g_1)

    # decide the stepsize
    #if abs(f_0-f_1)/f_0

    #if abs(f_0 - f_1) > threshold and abs(g_0 - g_1) > threshold:
    if abs(f_1) > threshold or abs(g_1) > threshold:
        print 'we need to find new beta and sigma'
        # refine the beta and sigma
        var_0 = np.concatenate((beta_0,sigma_0),axis = 0)
        Jacobi_Matrix = np.zeros((2, len(beta_0)+len(sigma_0)))
        # the partial deviation
        f_pd_beta = np.zeros(len(beta_0))
        f_pd_sigma = np.zeros(len(beta_0))
        for origin in range(len(Station_index)):
            f_pd_beta[origin] = f_pd_0(T_ij,D_ij,A_ij,origin)[0]
            f_pd_sigma[origin] = f_pd_0(T_ij,D_ij,A_ij,origin)[1]
        g_pd_beta = f_pd_sigma
        g_pd_sigma = np.zeros(len(sigma_0))
        for origin in range(len(Station_index)):
            g_pd_sigma[origin] = g_pd_0(T_ij,A_ij,origin)
        # calculate the jacobi matrix
        Jacobi_Matrix[0,:] = np.concatenate((f_pd_beta,f_pd_sigma),axis=0)
        Jacobi_Matrix[1,:] = np.concatenate((g_pd_beta,g_pd_sigma),axis=0)
        Jacobi_Matrix = matrix(Jacobi_Matrix)
        var_0 = matrix(var_0)
        direction = Jacobi_Matrix.I*matrix([[f_1],[g_1]])
        stepsize = np.zeros(len(var_0))
        for i in range(10):
            var_1 = var_0.T - Jacobi_Matrix.I*matrix([[f_1],[g_1]])*stepsize*(10**(-i))
            beta_0 = np.asarray(var_1.T)[0][0:len(beta_0)]
            sigma_0 = np.asarray(var_1.T)[0][len(beta_0):len(var_1)]
            Z_i = Z_i_value(Dt_j, A_ij, sigma_0, D_ij, beta_0, Station_index)
            T_ij = T_ij_value(Z_i, Ot_i, Dt_j, A_ij, sigma_0, D_ij, beta_0, Station_index)
            f_2 = f_value_0(T_ij, Tt_ij, D_ij, Station_index)
            g_2 = g_value_0(T_ij, Tt_ij, A_ij, Station_index)
            if f_2 <f_0_mid and g_2<g_0_mid:
                break
            else:
                continue
        # store the value
        #beta_0 = np.asarray(var_1.T)[0][0:len(beta_0)]
        #sigma_0 = np.asarray(var_1.T)[0][len(beta_0):len(var_1)]
        beta_value_store = np.vstack([beta_value_store, beta_0])
        sigma_value_store = np.vstack([sigma_value_store, sigma_0])
        f_0_mid = f_0
        f_0 = f_1
        f_value_store.append(f_0)
        g_0_mid = g_0
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