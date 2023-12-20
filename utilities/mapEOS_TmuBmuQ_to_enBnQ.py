#/usr/bin/env python3
# Copyright Chun Shen @ 2018

import numpy as np
import sys
from os import path
from scipy import interpolate


try:
    EOS_table_file = str(sys.argv[1])
except:
    print("Usage: python3 {} EOS_table_file".format(sys.argv[0]))
    exit(1)

ACCURACY = 1e-6
MAXITER  = 100
hbarC    = 0.19733

eos_table = np.loadtxt(EOS_table_file)

n_muQ = 275
n_muB = 801
n_T   = 181
muQarr = np.linspace(-0.137, 0.137, n_muQ)
muBarr = np.linspace(0, 0.8, n_muB)
Tarr = np.linspace(0.01, 0.19, n_T)

T_table   = eos_table[:, 0].reshape(n_T, n_muB, n_muQ)
muB_table = eos_table[:, 1].reshape(n_T, n_muB, n_muQ)
muQ_table = eos_table[:, 2].reshape(n_T, n_muB, n_muQ)
ed_table  = eos_table[:, 3].reshape(n_T, n_muB, n_muQ)*(T_table**4.)/(hbarC**3.)        # GeV/fm^3
nB_table  = eos_table[:, 4].reshape(n_T, n_muB, n_muQ)*(T_table**3.)/(hbarC**3.)        # 1/fm^3
s_table   = eos_table[:, 5].reshape(n_T, n_muB, n_muQ)*(T_table**3.)/(hbarC**3.)         # 1/fm^3
P_table   = eos_table[:, 6].reshape(n_T, n_muB, n_muQ)*(T_table**4.)/(hbarC**3.)        # GeV/fm^3
muS_table = eos_table[:, 7].reshape(n_T, n_muB, n_muQ)
nQ_table  = eos_table[:, -1].reshape(n_T, n_muB, n_muQ)*(T_table**3.)/(hbarC**3.)        # 1/fm^3

f_p   = interpolate.RegularGridInterpolator((Tarr, muBarr, muQ_arr), P_table )
f_e   = interpolate.RegularGridInterpolator((Tarr, muBarr, muQ_arr), ed_table)
f_nB  = interpolate.RegularGridInterpolator((Tarr, muBarr, muQ_arr), nB_table)
f_muS = interpolate.RegularGridInterpolator((Tarr, muBarr, muQ_arr), muS_table)
f_nQ  = interpolate.RegularGridInterpolator((Tarr, muBarr, muQ_arr), nQ_table)


def binary_search_1d(ed_local, muB_local, muQ_local):
    iteration = 0
    T_min = 0.01; T_max = 0.18
    e_low = f_e(T_min, muB_local, muQ_local)
    e_up  = f_e(T_max, muB_local, muQ_local)
    if (ed_local < e_low):
        return(T_min)
    elif (ed_local > e_up):
        return(T_max)
    else:
        T_mid = (T_max + T_min)/2.
        e_mid = f_e(T_mid, muB_local, muQ_local)
        abs_err = abs(e_mid - ed_local)
        rel_err = abs_err/abs(e_mid + ed_local + 1e-15)
        while (rel_err > ACCURACY and abs_err > ACCURACY*1e-2
                and iteration < MAXITER):
            if (ed_local < e_mid):
                T_max = T_mid
            else:
                T_min = T_mid
            T_mid = (T_max + T_min)/2.
            e_mid = f_e(T_mid, muB_local, muQ_local)
            abs_err = abs(e_mid - ed_local)
            rel_err = abs_err/abs(e_mid + ed_local + 1e-15)
            iteration += 1
        return(T_mid)


def binary_search_2d(ed_local, nB_local, muQ_local):
    iteration = 0
    muB_min = 0.0; muB_max = 0.8
    T_min = 0.01; T_max = 0.18
    T_max = binary_search_1d(ed_local, muB_min, muQ_local)
    nB_min = f_nB(T_max, muB_min, muQ_local)
    T_min = binary_search_1d(ed_local, muB_max, muQ_local)
    nB_max = f_nB(T_min, muB_max, muQ_local)
    if (nB_local < nB_min):
        return(T_max, muB_min)
    elif (nB_local > nB_max):
        return(T_min, muB_max)
    else:
        muB_mid = (muB_min + muB_max)/2.
        T_mid = binary_search_1d(ed_local, muB_mid)
        nB_mid = f_nB(T_mid, muB_mid, muQ_local)
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
            nB_mid = f_nB(T_mid, muB_mid, muQ_local)
            abs_err = abs(nB_mid - nB_local)
            rel_err = abs_err/abs(nB_mid + nB_local)
            iteration += 1
        return(T_mid, muB_mid)


def binary_search_3d(ed_local, nB_local, nQ_local):
    iteration = 0
    muQ_min = -0.137; muB_max = 0.137
    muB_min = 0.0; muB_max = 0.8
    T_min = 0.01; T_max = 0.18
    T_max, muB_min = binary_search_2d(ed_local, nB_local, muQ_min)
    nQ_min = f_nQ(T_max, muB_min, muQ_min)
    T_min, muB_max = binary_search_2d(ed_local, nB_local, muQ_max)
    nQ_max = f_nQ(T_min, muB_max, muQ_max)
    if (nQ_local < nQ_min):
        return(T_max, muB_min, muQ_min)
    elif (nQ_local > nQ_max):
        return(T_min, muB_max, muQ_max)
    else:
        muQ_mid = (muQ_min + muQ_max)/2.
        T_mid, muB_mid = binary_search_2d(ed_local, nB_local, muQ_mid)
        nQ_mid = f_nQ(T_mid, muB_mid, muQ_mid)
        abs_err = abs(nQ_mid - nQ_local)
        rel_err = abs_err/abs(nQ_mid + nQ_local + 1e-15)
        while (rel_err > ACCURACY and abs_err > ACCURACY*1e-2
                and iteration < MAXITER):
            if (nQ_local < nQ_mid):
                muQ_max = muQ_mid
            else:
                muQ_min = muQ_mid
            muQ_mid = (muQ_max + muQ_min)/2.
            T_mid, muB_mid = binary_search_2d(ed_local, nB_local, muQ_mid)
            nQ_mid = f_nQ(T_mid, muB_mid, muQ_mid)
            abs_err = abs(nQ_mid - nQ_local)
            rel_err = abs_err/abs(nQ_mid + nQ_local)
            iteration += 1
        return(T_mid, muB_mid, muQ_mid)


def invert_EOS_tables(ed_local, nB_local, nQ_local):
    T_local, muB_local, muQ_local = binary_search_3d(ed_local, nB_local, nQ_local)
    P_local            = f_p(T_local, muB_local, muQ_local)[0]
    muS_local          = f_muS(T_local, muB_local, muQ_local)[0]
    return(ed_local, nB_local, nQ_local, P_local, T_local,
           muB_local, muS_local, muQ_local)

#T_local, muB_local, muQ_local = binary_search_3d(1.0, 0.02, 0.01)
#print(T_local, muB_local, muQ_local, f_e(T_local, muB_local, muQ_local),
#      f_nB(T_local, muB_local, muQ_local), f_nQ(T_local, muB_local, muQ_local))

Ne = 200
ed_list = np.linspace(1e-2, 2, Ne)
nBmax_list = np.interp(ed_list, ed_table[:, 800, :], nB_table[:, 800, :])
nQmin_list = np.interp(ed_list, ed_table[:, :, 0], nQ_table[:, :, 0])
nQmax_list = np.interp(ed_list, ed_table[:, :, 274], nQ_table[:, :, 274])
NnB = 200
NnQ = 100

# generate tables
print("Generating table ... ")

# generate the grid arrays
output = []
for i, e_i in enumerate(ed_list):
    nB_list = np.linspace(0.0, nBmax_list[i], NnB)
    for j, nB_j in enumerate(nB_list):
        nQ_list = np.linspace(nQmin_list[j], nQmax_list[j], NnQ)
        for k, nQ_k in enumerate(nQ_list):
            output.append(invert_EOS_tables(e_i, nB_j, nQ_k))

# save to files
np.savetxt("NEOS_converted.dat", output, fmt='%.6e', delimiter="  ",
           header=("e[GeV/fm^3]  nB[1/fm^3]  nQ[1/fm^3]  P[GeV/fm^3]  T[GeV]  "
                   + "muB[GeV]  muS[GeV]  muQ[GeV]"))
