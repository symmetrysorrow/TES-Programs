# -*- coding: utf-8 -*-

# --------last updated 2018/12/01 by kurume-------------------

import math
from pprint import pprint
import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg as LA
import scipy.linalg
import scipy.sparse.linalg
import matplotlib.cm as cm
import matplotlib.ticker as ticker
import pandas as pd
import shutil
import json
import tqdm
import plt_config
import os
import natsort
import glob
import getpara as gp

# --------------------------------------------------------------
k_b = 1.381 * 1.0e-23  # Boltzmann's constant
ptfn_Flink = 0.5
e = 1.602e-19


def model(para):
    C_abs= para["C_abs"] / para['n_abs']
    G_abs_abs = para["G_abs-abs"] * (para['n_abs'] - 1) 
    energy = para['E']*1e3

    I = np.sqrt((para['G_tes-bath'] * para['T_c'] * (1 - ((para['T_bath'] / para['T_c']) ** para['n']))) / (para['n'] * para['R']))  # I_tes

    t_el = para['L']/ (para['R_l'] + para['R'] * (1 + para['beta']))  # tau_electron
    L_I = (para['alpha'] * (I**2) * para['R']) / (para['G_tes-bath'] * para['T_c'])  # Loop gain
    t_I = para['C_tes'] / ((1 - L_I) * para['G_tes-bath'])  # tau_?

    pixel = []
    cnt = 1
    for i in para['position']:
        if i > 0:
            pixel.append(cnt)
        cnt+=1
   
    def matrix_X(n_abs, posi):
        X = np.zeros((n_abs + 2, n_abs + 4))  # initialize matrix

        for i in range(n_abs + 2):
            if i == 0 or i == n_abs + 1:
                for j in range(n_abs + 4):
                    if j == i + 1:
                        X[i, j] = energy * e / para['C_tes']
            elif i in posi:
                for j in range(n_abs + 4):
                    if j == i + 1:
                        X[i, j] = energy * e / C_abs

            """
            else:
                for j in range(n_abs + 4):
                    if j == i + 1:
                        X[i, j] = E * e / C_abs
            """
        return X

    # matrix M without omega
    def matrix_M(n_abs):
        X = np.zeros((n_abs + 4, n_abs + 4))  # initialize matrix
        for i in range(n_abs + 4):
            if i == 0:
                for j in range(n_abs + 4):
                    if j == 0:
                        X[i, j] = 1 / t_el
                    elif j == 1:
                        X[i, j] = L_I * para['G_tes-bath'] / (I * para['L'])
            elif i == 1:
                for j in range(n_abs + 4):
                    if j == 0:
                        X[i, j] = -I * para['R'] * (2 + para['beta']) / para['C_tes']
                    elif j == 1:
                        X[i, j] = 1 / t_I + (para['G_abs-tes'] / para['C_tes'])
                    elif j == 2:
                        X[i, j] = -para['G_abs-tes'] / para['C_tes']
            elif i == 2:
                for j in range(n_abs + 4):
                    if j == 1:
                        X[i, j] = -para['G_abs-tes'] / C_abs
                    elif j == 2:
                        X[i, j] = para['G_abs-tes'] / C_abs + G_abs_abs / C_abs
                    elif j == 3:
                        X[i, j] = -G_abs_abs / C_abs
            elif i == n_abs + 1:
                for j in range(n_abs + 4):
                    if j == n_abs:
                        X[i, j] = -G_abs_abs / C_abs
                    elif j == n_abs + 1:
                        X[i, j] = para['G_abs-tes'] / C_abs + G_abs_abs / C_abs
                    elif j == n_abs + 2:
                        X[i, j] = -para['G_abs-tes'] / C_abs
            elif i == n_abs + 2:
                for j in range(n_abs + 4):
                    if j == n_abs + 1:
                        X[i, j] = -para['G_abs-tes'] / para['C_tes']
                    elif j == n_abs + 2:
                        X[i, j] = 1 / t_I + (para['G_abs-tes'] / para['C_tes'])
                    elif j == n_abs + 3:
                        X[i, j] = -I * para['R'] * (2 + para['beta']) / para['C_tes']
            elif i == n_abs + 3:
                for j in range(n_abs + 4):
                    if j == n_abs + 2:
                        X[i, j] = L_I * para['G_tes-bath'] / (I * para['L'])
                    elif j == n_abs + 3:
                        X[i, j] = 1 / t_el
            else:
                for j in range(n_abs + 4):
                    if j == i - 1:
                        X[i, j] = -G_abs_abs / C_abs
                    elif j == i:
                        X[i, j] = 2 * G_abs_abs / C_abs
                    elif j == i + 1:
                        X[i, j] = -G_abs_abs / C_abs
        return X

    M = matrix_M(para['n_abs']) * -1
    #print(pd.DataFrame(M))

    X = matrix_X(para['n_abs'], para["position"])

    val, vec = scipy.linalg.eig(
        M, left=False, right=True, overwrite_a=False, check_finite=True
    )

    arb = [np.linalg.solve(vec, i) for i in X]

    time = np.linspace(0, para['samples'] / para['rate'], int(para['samples']))

    I_t_0 = []
    I_t_sub_0 = []
    I_t_1 = []
    I_t_sub_1 = []
    
    for i in range(len(pixel)):
        I_t_sub_0.append([])
        I_t_sub_1.append([])

    cnt=0
    for i in para["position"]:
        for j in range(para['n_abs'] + 4):
            I_t_sub_0[cnt].append(arb[i][j] * vec[0][j] * np.exp(val[j] * time))
        pulse = (sum(I_t_sub_0[cnt])).real
        I_t_0.append(pulse)
        cnt+=1


    

    cnt=0
    for i in para["position"]:
        for j in range(para['n_abs'] + 4):
            I_t_sub_1[cnt].append(arb[i][j] * vec[para['n_abs'] + 3][j] * np.exp(val[j] * time))
        pulse = (sum(I_t_sub_1[cnt])).real
        I_t_1.append(pulse)
        cnt+=1
    
    
    return np.array(I_t_0)*-1,np.array(I_t_1)*-1
