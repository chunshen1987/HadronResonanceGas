#/usr/bin/env python3
# Copyright Chun Shen @ 2022

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

f_p   = interpolate.interp1d(eos_table[:, 0], eos_table[:, 4], kind='cubic')
f_e   = interpolate.interp1d(eos_table[:, 0], eos_table[:, 1], kind='cubic')
f_cs2 = interpolate.interp1d(eos_table[:, 0], eos_table[:, 5], kind='cubic')

def binary_search_1d(ed_local):
    iteration = 0
    T_min = 0.01; T_max = 0.19
    e_low = f_e(T_min)
    e_up  = f_e(T_max)
    if (ed_local < e_low):
        return(T_min)
    elif (ed_local > e_up):
        return(T_max)
    else:
        T_mid = (T_max + T_min)/2.
        e_mid = f_e(T_mid)
        abs_err = abs(e_mid - ed_local)
        rel_err = abs_err/abs(e_mid + ed_local + 1e-15)
        while (rel_err > ACCURACY and abs_err > ACCURACY*1e-2
                and iteration < MAXITER):
            if (ed_local < e_mid):
                T_max = T_mid
            else:
                T_min = T_mid
            T_mid = (T_max + T_min)/2.
            e_mid = f_e(T_mid)
            abs_err = abs(e_mid - ed_local)
            rel_err = abs_err/abs(e_mid + ed_local + 1e-15)
            iteration += 1
        return(T_mid)

#T_local = binary_search_1d(0.234)
#print(T_local, f_e(T_local))

ed_list = np.linspace(1e-3, 2, 2000)
de = ed_list[1] - ed_list[0]
p_list = np.zeros(len(ed_list))
T_list = np.zeros(len(ed_list))
cs2_list = np.zeros(len(ed_list))

# generate tables
print("Generating table ... ")

# fill into the tables
for idx in range(len(ed_list)):
    T_local = binary_search_1d(ed_list[idx])
    T_list[idx] = T_local
    p_list[idx] = f_p(T_local)
    cs2_list[idx] = f_cs2(T_local)

s_list = (ed_list + p_list)/T_list
# save to files
output = np.array([ed_list, p_list, s_list, T_list, cs2_list])
np.savetxt("EOSTable_PST.dat", output.transpose(), fmt="%.6e", delimiter="  ",
           header="e [GeV/fm^3]  P [GeV/fm^3]  s [1/fm^3]  T [GeV]  cs^2")
