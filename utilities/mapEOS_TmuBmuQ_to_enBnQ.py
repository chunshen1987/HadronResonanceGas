#/usr/bin/env python3
# Copyright Chun Shen @ 2018

import numpy as np
import sys
from os import path
from scipy import interpolate
from joblib import Parallel, delayed


try:
    EOS_table_file = str(sys.argv[1])
except:
    print("Usage: python3 {} EOS_table_file".format(sys.argv[0]))
    exit(1)

ACCURACY = 1e-6
MAXITER  = 100
hbarC    = 0.19733

eos_table = np.loadtxt(EOS_table_file)

n_muQ = 53; muQMIN = -0.13; muQMAX = 0.13
n_muB = 161; muBMAX = 0.8
n_T   = 33; TMIN = 0.01; TMAX = 0.17
muQarr = np.linspace(muQMIN, muQMAX, n_muQ)
muBarr = np.linspace(0, muBMAX, n_muB)
Tarr = np.linspace(TMIN, TMAX, n_T)

T_table   = eos_table[:, 0].reshape(n_T, n_muB, n_muQ)
muB_table = eos_table[:, 1].reshape(n_T, n_muB, n_muQ)
muQ_table = eos_table[:, 2].reshape(n_T, n_muB, n_muQ)
ed_table  = eos_table[:, 3].reshape(n_T, n_muB, n_muQ)*(T_table**4.)/(hbarC**3.)        # GeV/fm^3
nB_table  = eos_table[:, 4].reshape(n_T, n_muB, n_muQ)*(T_table**3.)/(hbarC**3.)        # 1/fm^3
s_table   = eos_table[:, 5].reshape(n_T, n_muB, n_muQ)*(T_table**3.)/(hbarC**3.)        # 1/fm^3
P_table   = eos_table[:, 6].reshape(n_T, n_muB, n_muQ)*(T_table**4.)/(hbarC**3.)        # GeV/fm^3
muS_table = eos_table[:, 7].reshape(n_T, n_muB, n_muQ)
nQ_table  = eos_table[:, -1].reshape(n_T, n_muB, n_muQ)*(T_table**3.)/(hbarC**3.)       # 1/fm^3

f_p   = interpolate.RegularGridInterpolator((Tarr, muBarr, muQarr), P_table )
f_e   = interpolate.RegularGridInterpolator((Tarr, muBarr, muQarr), ed_table)
f_nB  = interpolate.RegularGridInterpolator((Tarr, muBarr, muQarr), nB_table)
f_muS = interpolate.RegularGridInterpolator((Tarr, muBarr, muQarr), muS_table)
f_nQ  = interpolate.RegularGridInterpolator((Tarr, muBarr, muQarr), nQ_table)


def binary_search_1d(ed_local, muB_local, muQ_local):
    """
        return temperature at a given (e, muB, muQ)
    """
    iteration = 0
    T_min = TMIN; T_max = TMAX
    e_low = f_e(np.array([T_min, muB_local, muQ_local]))[0]
    e_up  = f_e(np.array([T_max, muB_local, muQ_local]))[0]
    if (ed_local < e_low):
        return(T_min)
    elif (ed_local > e_up):
        return(T_max)
    else:
        T_mid = (T_max + T_min)/2.
        e_mid = f_e(np.array([T_mid, muB_local, muQ_local]))[0]
        abs_err = abs(e_mid - ed_local)
        rel_err = abs_err/abs(e_mid + ed_local + 1e-15)
        while (rel_err > ACCURACY and abs_err > ACCURACY*1e-2
                and iteration < MAXITER):
            if (ed_local < e_mid):
                T_max = T_mid
            else:
                T_min = T_mid
            T_mid = (T_max + T_min)/2.
            e_mid = f_e(np.array([T_mid, muB_local, muQ_local]))[0]
            abs_err = abs(e_mid - ed_local)
            rel_err = abs_err/abs(e_mid + ed_local + 1e-15)
            iteration += 1
        return T_mid


def binary_search_2d(ed_local, nB_local, muQ_local):
    """
        return temperature and muB at a given (e, nB, muQ)
    """
    iteration = 0
    muB_min = 0.0; muB_max = muBMAX
    T_max = binary_search_1d(ed_local, muB_min, muQ_local)
    nB_min = f_nB(np.array([T_max, muB_min, muQ_local]))[0]
    T_min = binary_search_1d(ed_local, muB_max, muQ_local)
    nB_max = f_nB(np.array([T_min, muB_max, muQ_local]))[0]
    if (nB_local < nB_min):
        return T_max, muB_min
    elif (nB_local > nB_max):
        return T_min, muB_max
    else:
        muB_mid = (muB_min + muB_max)/2.
        T_mid = binary_search_1d(ed_local, muB_mid, muQ_local)
        nB_mid = f_nB(np.array([T_mid, muB_mid, muQ_local]))[0]
        abs_err = abs(nB_mid - nB_local)
        rel_err = abs_err/max(1e-15, abs(nB_mid + nB_local))
        while (rel_err > ACCURACY and abs_err > ACCURACY*1e-2
                and iteration < MAXITER):
            if (nB_local < nB_mid):
                muB_max = muB_mid
            else:
                muB_min = muB_mid
            muB_mid = (muB_max + muB_min)/2.
            T_mid = binary_search_1d(ed_local, muB_mid, muQ_local)
            nB_mid = f_nB(np.array([T_mid, muB_mid, muQ_local]))[0]
            abs_err = abs(nB_mid - nB_local)
            rel_err = abs_err/max(1e-15, abs(nB_mid + nB_local))
            iteration += 1
        return T_mid, muB_mid


def binary_search_3d(ed_local, nB_local, nQ_local):
    """
        return (T, muB, muQ) at a given (e, nB, nQ)
    """
    iteration = 0
    muQ_min = muQMIN; muQ_max = muQMAX
    T1, muB1 = binary_search_2d(ed_local, nB_local, muQ_min)
    nQ_min = f_nQ(np.array([T1, muB1, muQ_min]))
    T2, muB2 = binary_search_2d(ed_local, nB_local, muQ_max)
    nQ_max = f_nQ(np.array([T2, muB2, muQ_max]))
    if (nQ_local < nQ_min):
        return(T1, muB1, muQ_min)
    elif (nQ_local > nQ_max):
        return(T2, muB2, muQ_max)
    else:
        muQ_mid = (muQ_min + muQ_max)/2.
        T_mid, muB_mid = binary_search_2d(ed_local, nB_local, muQ_mid)
        nQ_mid = f_nQ(np.array([T_mid, muB_mid, muQ_mid]))
        abs_err = abs(nQ_mid - nQ_local)
        rel_err = abs_err/max(1e-15, abs(nQ_mid + nQ_local))
        while (rel_err > ACCURACY and abs_err > ACCURACY*1e-2
                and iteration < MAXITER):
            if (nQ_local < nQ_mid):
                muQ_max = muQ_mid
            else:
                muQ_min = muQ_mid
            muQ_mid = (muQ_max + muQ_min)/2.
            T_mid, muB_mid = binary_search_2d(ed_local, nB_local, muQ_mid)
            nQ_mid = f_nQ(np.array([T_mid, muB_mid, muQ_mid]))
            abs_err = abs(nQ_mid - nQ_local)
            rel_err = abs_err/max(1e-15, abs(nQ_mid + nQ_local))
            iteration += 1
        return T_mid, muB_mid, muQ_mid


def invert_EOS_tables(ed_local, nB_local, nQ_local):
    T_local, muB_local, muQ_local = binary_search_3d(ed_local, nB_local,
                                                     nQ_local)
    P_local = f_p(np.array([T_local, muB_local, muQ_local]))[0]
    muS_local = f_muS(np.array([T_local, muB_local, muQ_local]))[0]
    return [ed_local, nB_local, nQ_local, P_local, T_local,
            muB_local, muS_local, muQ_local]

e_check = 0.6; nB_check = 0.0; nQ_check = 0.0
T_local, muB_local, muQ_local = binary_search_3d(e_check, nB_check, nQ_check)
print(f"check e = {e_check}, nB = {nB_check}, nQ = {nQ_check}")
print(f"T = {T_local:.3e}, muB = {muB_local:.3e}, muQ = {muQ_local:.3e},"
      + f"e0 = {f_e(np.array([T_local, muB_local, muQ_local]))[0]:.3f},"
      + f"nB = {f_nB(np.array([T_local, muB_local, muQ_local]))[0]:.3f},"
      + f"nQ = {f_nQ(np.array([T_local, muB_local, muQ_local]))[0]:.3f}")
exit(0)

Ne = 60
NnB = 50
NnQ = 30
ed_list = np.linspace(1e-2, 0.6, Ne)

# generate tables
print(f"Generating table ... ")

# generate the grid arrays
output = []
for i, e_i in enumerate(ed_list):
    print(f"e = {e_i:.3f} GeV/fm^3 ...")
    nBmax = 0.0
    for k in range(n_muQ):
        nBmax_tmp = np.interp(e_i, ed_table[:, n_muB-1, k],
                              nB_table[:, n_muB-1, k])
        if nBmax < nBmax_tmp:
            nBmax = nBmax_tmp
    nB_list = np.linspace(0.0, nBmax, NnB)
    for j, nB_j in enumerate(nB_list):
        print(f"nB = {nB_j:.3e} fm^{-3} ...")
        T1, muB1 = binary_search_2d(e_i, nB_j, muQMIN)
        nQmin = f_nQ(np.array([T1, muB1, muQMIN]))[0]
        T2, muB2 = binary_search_2d(e_i, nB_j, muQMAX)
        nQmax = f_nQ(np.array([T2, muB2, muQMAX]))[0]
        nQ_list = np.linspace(nQmin, nQmax, NnQ)
        print(f"  nQmin = {nQmin:.3e} fm^{-3}, nQmax = {nQmax:.3e} fm^{-3}")
        #for k, nQ_k in enumerate(nQ_list):
        #    output.append(invert_EOS_tables(e_i, nB_j, nQ_k))
        num_jobs = -1  # Set to -1 to use all available cores
        results = Parallel(n_jobs=num_jobs)(
            delayed(invert_EOS_tables)(e_i, nB_j, nQ_k) for nQ_k in nQ_list)
        results = np.array(results)
        sorted_indices = np.argsort(results[:, 2])
        sorted_results = results[sorted_indices]
        #print(sorted_results)
        output.append(sorted_results)

# save to files
np.savetxt("NEOS_converted.dat", output, fmt='%.6e', delimiter="  ",
           header=("e[GeV/fm^3]  nB[1/fm^3]  nQ[1/fm^3]  P[GeV/fm^3]  T[GeV]  "
                   + "muB[GeV]  muS[GeV]  muQ[GeV]"))
