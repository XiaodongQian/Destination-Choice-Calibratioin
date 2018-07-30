import csv
import numpy as np
import math
import codecs
import os
from numpy import matrix
from numpy import linalg as LA

def FG_value(T_ij, Tt_ij, D_ij, A_ij, Station_index):
    fg_value = np.zeros(len(Station_index)*2)
    iter = len(Station_index)
    # add f value
    for origin in range(iter):
        f_value = 0
        for dest in range(iter):
            if origin != dest:
                f_value += (T_ij[origin][dest]- Tt_ij[origin][dest])*np.log10(D_ij[origin][dest])
        fg_value[origin] = f_value
    # add g value
    for origin in range(iter):
        g_value = 0
        for dest in range(iter):
            if origin != dest:
                g_value += (T_ij[origin][dest]- Tt_ij[origin][dest])*np.log10(A_ij[origin][dest])
        fg_value[origin + iter] = g_value


    return fg_value


def Z_i_value(Dt_j,A_ij,sigma_0,D_ij,beta_0,Station_index):
    Z_i = np.zeros((len(Station_index), 1))
    for origin in range(len(Station_index)):
        z_i = 0
        for dest in range(len(Station_index)):
            if dest != origin:
                #z_i += Dt_j[dest] * math.pow(A_ij[origin][dest], sigma_0[origin]) * math.pow(D_ij[origin][dest],beta_0[origin])
                try:
                    z_i += Dt_j[dest] * math.pow(A_ij[origin][dest], sigma_0[origin]) * math.pow(D_ij[origin][dest], beta_0[origin])
                except OverflowError:
                    z_i = np.inf
                    break

        Z_i[origin] = math.pow(z_i, -1)
    return Z_i


def T_ij_value(Z_i,Ot_i,Dt_j,A_ij,sigma_0,D_ij,beta_0,Station_index):
    T_ij = np.zeros((len(Station_index), len(Station_index)))
    for origin in range(len(Station_index)):
        for dest in range(len(Station_index)):
            if origin != dest:
                T_ij[origin][dest] = (Z_i[origin][0] * Ot_i[origin] * Dt_j[dest] * math.pow(A_ij[origin][dest], sigma_0[origin]) * math.pow(D_ij[origin][dest], beta_0[origin]))
                '''
                if Z_i[origin][0] == 0:
                    T_ij[origin][dest] = 0
                else:
                    print "math.pow(A_ij[origin][dest],sigma_0[origin])"
                    print origin, dest
                    print math.pow(A_ij[origin][dest],sigma_0[origin])
                    print 'math.pow(D_ij[origin][dest], beta_0[origin])'
                    print math.pow(D_ij[origin][dest], beta_0[origin])
                    T_ij[origin][dest] = (Z_i[origin][0] * Ot_i[origin] * Dt_j[dest] * math.pow(A_ij[origin][dest],sigma_0[origin]) * math.pow(D_ij[origin][dest], beta_0[origin]))
                '''
    return T_ij


def Jacobi_Matrix(T_ij,D_ij,A_ij,Station_index):
    iter = len(Station_index)
    jacobi_matrix = np.zeros((iter*2,iter*2), dtype = np.float32)
    # the partial deviation
    for origin in range(iter):
        jacobi_matrix[origin][origin] = F_beta(T_ij,D_ij,origin,Station_index)
        jacobi_matrix[origin][origin+iter] = F_sigma(T_ij, D_ij,A_ij, origin, Station_index)
        jacobi_matrix[iter+origin][origin] = G_beta(T_ij, D_ij, A_ij, origin, Station_index)
        jacobi_matrix[iter + origin][iter+ origin] = G_sigma(T_ij, A_ij, origin, Station_index)

    return jacobi_matrix


def F_beta(T_ij,D_ij,origin,Station_index):
    iter = len(Station_index)
    f_beta = 0
    for dest in range(iter):
        if dest != origin:
            f_beta += T_ij[origin][dest] * (np.log10(D_ij[origin][dest]))**2

    return f_beta


def F_sigma(T_ij,D_ij,A_ij,origin,Station_index):
    iter = len(Station_index)
    f_sigma = 0
    for dest in range(iter):
        if dest != origin:
            f_sigma += T_ij[origin][dest] * np.log10(D_ij[origin][dest]) * np.log10(A_ij[origin][dest])

    return f_sigma


def G_beta(T_ij,D_ij,A_ij,origin,Station_index):
    iter = len(Station_index)
    g_beta = 0
    for dest in range(iter):
        if dest != origin:
            g_beta += T_ij[origin][dest] * np.log10(D_ij[origin][dest]) * np.log10(A_ij[origin][dest])

    return g_beta


def G_sigma(T_ij,A_ij,origin,Station_index):
    iter = len(Station_index)
    g_sigma = 0
    for dest in range(iter):
        if dest != origin:
            g_sigma += T_ij[origin][dest] * (np.log10(A_ij[origin][dest]))**2

    return g_sigma


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

# Calculate the A_ij Matrix
rho = -1
distance = 3000
A_ij = A_ij_matrix(Dt_j,D_ij,Station_index,distance,rho)

##############
#  calibration process
## initial value
#beta_0 = np.asarray(np.random.normal(-1.5, 0.5, len(Station_index)))
#sigma_0 = np.asarray(np.random.normal(2.5, 1, len(Station_index)))
beta_0 = np.asarray([-1.5] *len(Station_index), dtype=np.float64)
sigma_0 = np.asarray([2.5] *len(Station_index),dtype=np.float64)
stepsize_1 = 2 # (0,1)
stepsize_2 = 0.25 # rho and (0, 0.5)
iteration = 1
threshold = 100
Z_i = np.zeros((len(Station_index),1))
T_ij = np.zeros((len(Station_index),len(Station_index)))

## keep the value change for every iteration
fg_value_store = []
beta_value_store = beta_0
sigma_value_store = sigma_0
while iteration < 1000:
    print '############'
    print 'This is the ' + str(iteration) + 'th iteration'
    Z_i = Z_i_value(Dt_j, A_ij, sigma_0, D_ij, beta_0, Station_index)
    T_ij = T_ij_value(Z_i, Ot_i, Dt_j, A_ij, sigma_0, D_ij, beta_0, Station_index)
    jacobi_matrix = Jacobi_Matrix(T_ij,D_ij,A_ij,Station_index)
    jacobi_matrix = matrix(jacobi_matrix)
    var_0 = np.concatenate((beta_0, sigma_0), axis=0)
    var_0 = matrix(var_0)
    fg_value = FG_value(T_ij, Tt_ij, D_ij, A_ij, Station_index)
    LA_diff = LA.norm(fg_value)
    fg_value_store.append(LA_diff)
    print LA_diff
    if LA_diff > threshold:
        try:
            direction = - 0.01 * jacobi_matrix.I * matrix(fg_value).T
        except:
            [m, n] = jacobi_matrix.shape;
            jacobi_matrix_new = jacobi_matrix
            for index in range(m):
                jacobi_matrix_new[index,index] += 0.1
            direction = - 0.01 * jacobi_matrix_new.I * matrix(fg_value).T

        # decide the stepsize
        m = 0
        fg_value_mid = LA_diff
        while m < 8:
            var_mid = var_0.T + (stepsize_1 ** m) * direction
            beta_mid = np.asarray(var_mid.T)[0][0:len(beta_0)]
            sigma_mid = np.asarray(var_mid.T)[0][len(beta_0):len(var_mid)]
            try:
                Z_i_mid = Z_i_value(Dt_j, A_ij, sigma_mid, D_ij, beta_mid, Station_index)
                T_ij_mid = T_ij_value(Z_i_mid, Ot_i, Dt_j, A_ij, sigma_mid, D_ij, beta_mid, Station_index)
                fg_value_left = LA.norm(FG_value(T_ij_mid, Tt_ij, D_ij, A_ij, Station_index))
                print fg_value_left
                if fg_value_left < fg_value_mid:
                    m += 1
                    fg_value_mid = fg_value_left
                    print m
                else:
                    m -= 1
                    break
            except:
                m -= 1
                stepsize = stepsize_1 ** m
                var_1 = var_0.T + stepsize * direction
                beta_0 = np.asarray(var_1.T)[0][0:len(beta_0)]
                sigma_0 = np.asarray(var_1.T)[0][len(beta_0):len(var_1)]
                try:
                    Z_i_value(Dt_j, A_ij, sigma_0, D_ij, beta_0, Station_index)
                    T_ij_value(Z_i, Ot_i, Dt_j, A_ij, sigma_0, D_ij, beta_0, Station_index)
                    break
                except:
                    m = -3
                    break
        '''
        while m <20:
            var_mid = var_0.T + (stepsize_1**m)*direction
            beta_mid = np.asarray(var_mid.T)[0][0:len(beta_0)]
            sigma_mid = np.asarray(var_mid.T)[0][len(beta_0):len(var_mid)]
            Z_i_mid = Z_i_value(Dt_j, A_ij, sigma_mid, D_ij, beta_mid, Station_index)
            T_ij_mid = T_ij_value(Z_i_mid, Ot_i, Dt_j, A_ij, sigma_mid, D_ij, beta_mid, Station_index)
            fg_value_left = FG_value(T_ij_mid, Tt_ij, D_ij, A_ij, Station_index)
            add = LA.norm((stepsize_1**m) * stepsize_2 * matrix(fg_value) * direction)
            fg_value_right = [i + add for i in fg_value]
            if LA.norm(fg_value_left) <= LA.norm(fg_value_right):
                print m
                break
            else:
                m += 1
        '''
        stepsize = stepsize_1**m
        print stepsize
        var_1 = var_0.T + stepsize * direction
        beta_0 = np.asarray(var_1.T)[0][0:len(beta_0)]
        sigma_0 = np.asarray(var_1.T)[0][len(beta_0):len(var_1)]
        beta_value_store = np.vstack([beta_value_store, beta_0])
        sigma_value_store = np.vstack([sigma_value_store, sigma_0])
        iteration += 1
    else:
        print 'done'
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