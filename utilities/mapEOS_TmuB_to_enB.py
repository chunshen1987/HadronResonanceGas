#/usr/bin/env python3
# Copyright Chun Shen @ 2018

import numpy as np
import sys
from os import path
from scipy import interpolate


try:
    EOS_table_file = str(sys.argv[1])
    mode = int(sys.argv[2])
except:
    print("Usage: python3 {} EOS_table_file mode".format(sys.argv[0]))
    print("       mode = 0: NEoS-B")
    print("       mode = 1: NEoS-BS")
    print("       mode = 2: NEoS-BQS")
    exit(1)

ACCURACY = 1e-6
MAXITER  = 100
hbarC    = 0.19733

eos_table = np.loadtxt(EOS_table_file)

n_muB = 801
n_T   = 181
muBarr = np.linspace(0, 0.8, n_muB)
Tarr = np.linspace(0.01, 0.19, n_T)

T_table   = eos_table[:, 0].reshape(n_T, n_muB)
muB_table = eos_table[:, 1].reshape(n_T, n_muB)
ed_table  = eos_table[:, 2].reshape(n_T, n_muB)*(T_table**4.)/(hbarC**3.)        # GeV/fm^3
nB_table  = eos_table[:, 3].reshape(n_T, n_muB)*(T_table**3.)/(hbarC**3.)        # 1/fm^3
s_table   = eos_table[:, 4].reshape(n_T, n_muB)*(T_table**3.)/(hbarC**3.)         # 1/fm^3
P_table   = eos_table[:, 5].reshape(n_T, n_muB)*(T_table**4.)/(hbarC**3.)        # GeV/fm^3
muS_table = np.zeros([n_T, n_muB])
muQ_table = np.zeros([n_T, n_muB])

if mode == 1:
    muS_table = eos_table[:, 6].reshape(n_T, n_muB)
elif mode == 2:
    muS_table = eos_table[:, 6].reshape(n_T, n_muB)
    muQ_table = eos_table[:, 7].reshape(n_T, n_muB)

f_p   = interpolate.RectBivariateSpline(Tarr, muBarr, P_table )
f_e   = interpolate.RectBivariateSpline(Tarr, muBarr, ed_table)
f_nB  = interpolate.RectBivariateSpline(Tarr, muBarr, nB_table)
f_muS = interpolate.RectBivariateSpline(Tarr, muBarr, muS_table)
f_muQ = interpolate.RectBivariateSpline(Tarr, muBarr, muQ_table)


def binary_search_1d(ed_local, muB_local):
    iteration = 0
    T_min = 0.01; T_max = 0.20
    e_low = f_e(T_min, muB_local)
    e_up  = f_e(T_max, muB_local)
    if (ed_local < e_low):
        return(T_min)
    elif (ed_local > e_up):
        return(T_max)
    else:
        T_mid = (T_max + T_min)/2.
        e_mid = f_e(T_mid, muB_local)
        abs_err = abs(e_mid - ed_local)
        rel_err = abs_err/abs(e_mid + ed_local + 1e-15)
        while (rel_err > ACCURACY and abs_err > ACCURACY*1e-2
                and iteration < MAXITER):
            if (ed_local < e_mid):
                T_max = T_mid
            else:
                T_min = T_mid
            T_mid = (T_max + T_min)/2.
            e_mid = f_e(T_mid, muB_local)
            abs_err = abs(e_mid - ed_local)
            rel_err = abs_err/abs(e_mid + ed_local + 1e-15)
            iteration += 1
        return(T_mid)


def binary_search_2d(ed_local, nB_local):
    iteration = 0
    muB_min = 0.0; muB_max = 0.8
    T_min = 0.01; T_max = 0.20
    T_max = binary_search_1d(ed_local, muB_min)
    nB_min = f_nB(T_max, muB_min)
    T_min = binary_search_1d(ed_local, muB_max)
    nB_max = f_nB(T_min, muB_max)
    if (nB_local < nB_min):
        return(T_max, muB_min)
    elif (nB_local > nB_max):
        return(T_min, muB_max)
    else:
        muB_mid = (muB_min + muB_max)/2.
        T_mid = binary_search_1d(ed_local, muB_mid)
        nB_mid = f_nB(T_mid, muB_mid)
        abs_err = abs(nB_mid - nB_local)
        rel_err = abs_err/abs(nB_mid + nB_local + 1e-15)
        while (rel_err > ACCURACY and abs_err > ACCURACY*1e-2
                and iteration < MAXITER):
            if (nB_local < nB_mid):
                muB_max = muB_mid
            else:
                muB_min = muB_mid
            muB_mid = (muB_max + muB_min)/2.
            T_mid = binary_search_1d(ed_local, muB_mid)
            nB_mid = f_nB(T_mid, muB_mid)
            abs_err = abs(nB_mid - nB_local)
            rel_err = abs_err/abs(nB_mid + nB_local)
            iteration += 1
        return(T_mid, muB_mid)


def invert_EOS_tables(ed_local, nB_local):
    T_local, muB_local = binary_search_2d(ed_local, nB_local)
    P_local            = f_p(T_local, muB_local)[0]
    muS_local          = f_muS(T_local, muB_local)[0]
    muQ_local          = f_muQ(T_local, muB_local)[0]
    return(ed_local, nB_local ,P_local, T_local,
           muB_local, muS_local, muQ_local)

#T_local, muB_local = binary_search_2d(1.0, 0.02)
#print(T_local, muB_local, f_e(T_local, muB_local), f_nB(T_local, muB_local))

Ne = 200
ed_list = np.linspace(1e-2, 2, Ne)
nBmax_list = np.interp(ed_list, ed_table[:, 800], nB_table[:, 800])
NnB = 200

# generate tables
print("Generating table ... ")

# generate the grid arrays
output = []
for i, e_i in enumerate(ed_list):
    nB_list = np.linspace(0.0, nBmax_list[i], NnB)
    for j, nB_j in enumerate(nB_list):
        output.append(invert_EOS_tables(e_i, nB_j))

# save to files
np.savetxt("NEOS_converted.dat", output, fmt='%.6e', delimiter="  ",
           header=("e[GeV/fm^3]  nB[1/fm^3]  P[GeV/fm^3]  T[GeV]  "
                   + "muB[GeV]  muS[GeV]  muQ[GeV]"))
