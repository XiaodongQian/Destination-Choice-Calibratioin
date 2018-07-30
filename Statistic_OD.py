import csv
import numpy as np
import math
import codecs
import os

# read OD matrix
## calculate the sum(OD_ij)
path = '/Volumes/GoogleDrive/My Drive/Dissertation/Data/Chicago/OD'
with codecs.open(os.path.join(path, 'Chi_Annual_Bike_OD_20151231_20161231_AnnualMember.csv'),'r',encoding="utf-8-sig") as f:
    reader = csv.reader(f)
    x = list(reader)
    Tt_ij = np.array(x).astype("float")
## calculate the Ot_i and Dt_j
Ot_i = np.sum(Tt_ij, axis=1)
Dt_j = np.sum(Tt_ij, axis=0)

OD_diag = np.zeros((len(Tt_ij[0]), 1))
O_diag_per = np.zeros((len(Tt_ij[0]), 1))
D_diag_per = np.zeros((len(Tt_ij[0]), 1))
for origin in range(len(Tt_ij[0])):
    OD_diag[origin] = Tt_ij[origin][origin]
    O_diag_per[origin] = Tt_ij[origin][origin]/Ot_i[origin]
    D_diag_per[origin] = Tt_ij[origin][origin] /Dt_j[origin]

# plot the histgram
import matplotlib.pyplot as plt
plt.hist(D_diag_per)  # arguments are passed to np.histogram
plt.show()

OD_diff = np.zeros((len(Tt_ij[0]),len(Tt_ij[0])))
od_diff = []
for origin in range(len(Tt_ij[0])):
    for dest in range(len(Tt_ij[0])):
        if origin != dest:
            if Tt_ij[origin][dest]+Tt_ij[dest][origin] != 0:
                OD_diff[origin][dest] = (Tt_ij[origin][dest]-Tt_ij[dest][origin])*2/(Tt_ij[origin][dest]+Tt_ij[dest][origin])
                od_diff.append(OD_diff[origin][dest])

plt.hist(od_diff)  # arguments are passed to np.histogram
plt.show()
