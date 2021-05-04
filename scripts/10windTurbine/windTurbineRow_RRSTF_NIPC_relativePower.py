#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 20 23:48:12 2021

@author: p285464
"""

# <codecell> Import Required Modules and env vars
from __future__ import print_function

import os, sys
import numpy as np
import matplotlib.pyplot as plt

SCRIPTS = os.environ['SCRIPTS']
DATA    = os.environ['windTurbineData']

sys.path.append(SCRIPTS+'/10windTurbine')
sys.path.append(SCRIPTS+'/pythonScripts')

cwd = os.getcwd() + "/"

import myUQlib

# <codecell> Case details
caseName   = cwd.split('/')[-2]
casePngDir = '/2RRSTF/NIPC/' + caseName+'/png'

if (os.path.exists(DATA+casePngDir))==False:
    print('Making new case data directory...')
    os.makedirs(DATA+casePngDir, exist_ok=True)

# <codecell> Parametrs
nWT = 8#10

# AD powers HornRev_WD270_WS8_TI7.7_Ti2019.csv
# ADpow_EXP = np.loadtxt(DATA+"/0EXP/HornRev_WD270_WS8_TI7.7_Wu2015.csv")[:,1]
ADpow_EXP = np.loadtxt(DATA+"/0EXP/HornRev_WD270_WS8_TI7.7_Ti2019.csv")[:,1]
ADpow_kEpStd_7_7_x_up_D_2_0 = np.loadtxt("baseCase_kEpStd_TI_7.7_x_up_D_2.0_8WTs/power")
ADpow_LRRRSM_7_7_x_up_D_2_0 = np.loadtxt("baseCase_LRRRSM_TI_7.7_x_up_D_2.0_8WTs/power")
ADpow_LRRRSM_7_7_x_up_D_0_1_2_0 = np.loadtxt("baseCase_LRRRSM_TI_7.7_x_up_D_0.1_2.0_8WTs/power")
ADpow_LRRRSM_7_7_x_up_D_0_5_2_0 = np.loadtxt("baseCase_LRRRSM_TI_7.7_x_up_D_0.5_2.0_8WTs/power")
ADpow_LRRRSM_10_x_up_D_2_0 = np.loadtxt("baseCase_LRRRSM_TI_10._x_up_D_2.0_8WTs/power")
ADpow_LRRRSM_10_x_up_D_0_1_2_0 = np.loadtxt("baseCase_LRRRSM_TI_10._x_up_D_0.1_2.0_8WTs/power")
ADpow_LRRRSM_15_x_up_D_2_0 = np.loadtxt("baseCase_LRRRSM_TI_15._x_up_D_2.0_8WTs/power")

myUQlib.rcParamsSettings(15)
fig, ax = plt.subplots(ncols=1, nrows=1, constrained_layout=True, figsize=(6,4))
WTlist = np.arange(1,nWT+1)

ax.plot(WTlist, ADpow_EXP[:nWT],'ks', label='EXP')
ax.plot(WTlist, ADpow_kEpStd_7_7_x_up_D_2_0[:nWT]/ADpow_kEpStd_7_7_x_up_D_2_0[0],'rs', label='kEp TI 7.7 xbyDup 2.0')
# ax.plot(WTlist, ADpow_LRRRSM_7_7_x_up_D_2_0[:nWT]/ADpow_LRRRSM_7_7_x_up_D_2_0[0],'--bo', label='LRR 7.7 2.0')
# ax.plot(WTlist, ADpow_LRRRSM_7_7_x_up_D_0_1_2_0[:nWT]/ADpow_LRRRSM_7_7_x_up_D_0_1_2_0[0],'--bs', label='LRR 7.7 0.1 2.0')
# ax.plot(WTlist, ADpow_LRRRSM_7_7_x_up_D_0_5_2_0[:nWT]/ADpow_LRRRSM_7_7_x_up_D_0_5_2_0[0],'--bx', label='LRR 7.7 0.5 2.0')
# ax.plot(WTlist, ADpow_LRRRSM_10_x_up_D_2_0[:nWT]/ADpow_LRRRSM_10_x_up_D_2_0[0],'--go', label='LRR 10 2.0')
# ax.plot(WTlist, ADpow_LRRRSM_10_x_up_D_0_1_2_0[:nWT]/ADpow_LRRRSM_10_x_up_D_0_1_2_0[0],'--gx', label='LRR 10 0.1 2.0')
ax.plot(WTlist, ADpow_LRRRSM_15_x_up_D_2_0[:nWT]/ADpow_LRRRSM_15_x_up_D_2_0[0],'bs', label='LRR TI 15 xbyDup 2.0')

ax.set_ylim(0.2,1.1)
ax.set_xlim(0.8,8.2)
ax.legend(ncol=1)
ax.set_xticks(WTlist)
ax.set_yticks(np.linspace(0.2,1,9))
ax.grid()

fig.savefig(DATA+casePngDir+'/power.png', dpi=300, bbox_inches='tight')