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

# %% Case details
caseName   = cwd.split('/')[-2]
casePngDir = '/2RRSTF/NIPC/uqPaperCases/' + caseName+'/png'

if (os.path.exists(DATA+casePngDir))==False:
    print('Making new case data directory...')
    os.makedirs(DATA+casePngDir, exist_ok=True)

# <codecell> Parametrs
nWT = 8#10

# AD powers HornRev_WD270_WS8_TI7.7_Ti2019.csv
# ADpow_EXP = np.loadtxt(DATA+"/0EXP/HornRev_WD270_WS8_TI7.7_Wu2015.csv")[:,1]
ADpow_EXP = np.loadtxt(DATA+"/0EXP/HornRev_WD270_WS8_TI7.7_Ti2019.csv")[:,1]
ADpow_LES = np.array([1, 0.45, 0.53, 0.525, 0.525, 0.551])
ADpow_realKE_5_8_x_up_D_0_1 = np.loadtxt("baseCase_realKE_TI_5.8_x_up_D_0.1_6WTs_topInlet/power")
ADpow_realKE_5_8_x_up_D_1_0 = np.loadtxt("baseCase_realKE_TI_5.8_x_up_D_1.0_6WTs_topInlet/power")
ADpow_LRRRSM_5_8_x_up_D_1_0 = np.loadtxt("baseCase_LRRRSM_TI_5.8_x_up_D_1.0_6WTs_topInlet/power")

# ADpow_kEpStd_7_7_x_up_D_2_0 = np.loadtxt("baseCase_kEpStd_TI_7.7_x_up_D_2.0_8WTs/power")
# ADpow_LRRRSM_7_7_x_up_D_2_0 = np.loadtxt("baseCase_LRRRSM_TI_7.7_x_up_D_2.0_8WTs/power")
# ADpow_LRRRSM_7_7_x_up_D_0_1_2_0 = np.loadtxt("baseCase_LRRRSM_TI_7.7_x_up_D_0.1_2.0_8WTs/power")
# ADpow_LRRRSM_7_7_x_up_D_0_5_2_0 = np.loadtxt("baseCase_LRRRSM_TI_7.7_x_up_D_0.5_2.0_8WTs/power")
# ADpow_LRRRSM_10_x_up_D_2_0 = np.loadtxt("baseCase_LRRRSM_TI_10._x_up_D_2.0_8WTs/power")
# ADpow_LRRRSM_10_x_up_D_0_1_2_0 = np.loadtxt("baseCase_LRRRSM_TI_10._x_up_D_0.1_2.0_8WTs/power")
# ADpow_LRRRSM_15_x_up_D_2_0 = np.loadtxt("baseCase_LRRRSM_TI_15._x_up_D_2.0_8WTs/power")

myUQlib.rcParamsSettings(15)
fig, ax = plt.subplots(ncols=1, nrows=1, constrained_layout=True,
                       figsize=(8,3),dpi=150)
WTlist = np.arange(1,nWT+1)
DETClr  = (0.6350, 0.0780, 0.1840)

# ax.plot(WTlist[:-2], ADpow_EXP[:nWT-2],'ks', ls='-', label='EXP')
ax.scatter(WTlist[:-2], ADpow_LES, ec='k', fc='none', label='LES')
ax.scatter(WTlist[:-2], ADpow_LRRRSM_5_8_x_up_D_1_0[:nWT]/ADpow_LRRRSM_5_8_x_up_D_1_0[0],
        ec='k', fc='k',label='LRR TI 5.8 xbyDup 1.0')

ax.scatter(WTlist[:-2], ADpow_realKE_5_8_x_up_D_0_1[:nWT]/ADpow_realKE_5_8_x_up_D_0_1[0],
        ec='b', fc='b', label='Realizable$\ k-\epsilon$ TI 5.8 xbyDup 0.1')
ax.scatter(WTlist[:-2], ADpow_realKE_5_8_x_up_D_1_0[:nWT]/ADpow_realKE_5_8_x_up_D_1_0[0],
        ec=DETClr, fc=DETClr, label='Realizable\ $k-\epsilon$ TI 5.8 xbyDup 1.0')


# ax.plot(WTlist[:-2], ADpow_kEpStd_7_7_x_up_D_2_0[:nWT-2]/ADpow_kEpStd_7_7_x_up_D_2_0[0],
#         'gs', ls='--', label='kEp TI 7.7 xbyDup 2.0')
# ax.plot(WTlist, ADpow_LRRRSM_7_7_x_up_D_2_0[:nWT]/ADpow_LRRRSM_7_7_x_up_D_2_0[0],'--bo', label='LRR 7.7 2.0')
# ax.plot(WTlist, ADpow_LRRRSM_7_7_x_up_D_0_1_2_0[:nWT]/ADpow_LRRRSM_7_7_x_up_D_0_1_2_0[0],'--bs', label='LRR 7.7 0.1 2.0')
# ax.plot(WTlist, ADpow_LRRRSM_7_7_x_up_D_0_5_2_0[:nWT]/ADpow_LRRRSM_7_7_x_up_D_0_5_2_0[0],'--bx', label='LRR 7.7 0.5 2.0')
# ax.plot(WTlist, ADpow_LRRRSM_10_x_up_D_2_0[:nWT]/ADpow_LRRRSM_10_x_up_D_2_0[0],'--go', label='LRR 10 2.0')
# ax.plot(WTlist, ADpow_LRRRSM_10_x_up_D_0_1_2_0[:nWT]/ADpow_LRRRSM_10_x_up_D_0_1_2_0[0],'--gx', label='LRR 10 0.1 2.0')
# ax.plot(WTlist[:-2], ADpow_LRRRSM_15_x_up_D_2_0[:nWT-2]/ADpow_LRRRSM_15_x_up_D_2_0[0],
#         'bs', ls='--',label='LRR 15 xbyDup 2.0')

ax.set_ylim(0.2,1.2)
ax.set_xlim(0.5,6.5)
ax.legend(ncol=1)
ax.set_xticks(WTlist[:-2])
ax.set_yticks(np.linspace(0.2,1.2,6))
ax.set_xlabel('Turbine number')
ax.set_ylabel('Normalized power')
# ax.grid()

fig.savefig(DATA+casePngDir+'/power_baseCases.png', dpi=150, bbox_inches='tight')


# <codecell> CoCpmparison of distributions
nWT = 8#10

ADpow_EXP = np.loadtxt(DATA+"/0EXP/HornRev_WD270_WS8_TI7.7_Ti2019.csv")[:,1]
ADpow_LES = np.array([1, 0.45, 0.53, 0.525, 0.525, 0.551])
ADpow_realKE_5_8_x_up_D_0_1 = np.loadtxt("baseCase_realKE_TI_5.8_x_up_D_0.1_6WTs_topInlet/power")
ADpow_realKE_5_8_x_up_D_1_0 = np.loadtxt("baseCase_realKE_TI_5.8_x_up_D_1.0_6WTs_topInlet/power")
ADpow_LRRRSM_5_8_x_up_D_1_0 = np.loadtxt("baseCase_LRRRSM_TI_5.8_x_up_D_1.0_6WTs_topInlet/power")

myUQlib.rcParamsSettings(15)
fig = plt.figure(figsize=(8,3),dpi=150)
ax = fig.add_axes([0,0,1,1])

WTlist = np.arange(1,nWT+1)
DETClr  = (0.6350, 0.0780, 0.1840)

ax.scatter(WTlist[:-2], ADpow_LES, ec='k', fc='none', label='LES')
ax.scatter(WTlist[:-2], ADpow_LRRRSM_5_8_x_up_D_1_0[:nWT]/ADpow_LRRRSM_5_8_x_up_D_1_0[0],
        ec='k', fc='k',label='LRR')

ax.scatter(WTlist[:-2], ADpow_realKE_5_8_x_up_D_0_1[:nWT]/ADpow_realKE_5_8_x_up_D_0_1[0],
        ec='b', fc='b', label='Stochatic RANS Solver')
ax.scatter(WTlist[:-2], ADpow_realKE_5_8_x_up_D_1_0[:nWT]/ADpow_realKE_5_8_x_up_D_1_0[0],
        ec=DETClr, fc=DETClr, label='3D U-Net + Wake-superposition')

ax.set_ylim(0.0,1.2)
ax.set_xlim(0.5,6.5)
leg = ax.legend(ncol=2, loc='best', bbox_to_anchor=(0.5, 0.532, 0.5, 0.5),
                framealpha=1, edgecolor='none')
ax.set_xticks(WTlist[:-2])
ax.set_yticks(np.linspace(0.0,1.2,7))
ax.set_xlabel('Turbine number')
ax.set_ylabel('Normalized power')
# ax.spines[['right', 'top']].set_visible(False)

##%% Stats

pow_samples = []
for i in range(90):
    pow_samples.append(np.loadtxt(f"sample_{i}/power"))
pow_samples = np.array(pow_samples)
pow_samples_std = pow_samples.std(axis=0)

pow_samples_mean_normed = ADpow_realKE_5_8_x_up_D_0_1/ADpow_realKE_5_8_x_up_D_0_1[0]
pow_samples_std_normed = pow_samples_std/ADpow_realKE_5_8_x_up_D_0_1[0]


##%% Distributions

from scipy.stats import norm

i = 1
for i in range(1,6):
    N, scale = 10, 10
    x_axis = np.linspace(pow_samples_mean_normed[i]+N*pow_samples_std_normed[i], 
                         pow_samples_mean_normed[i]-N*pow_samples_std_normed[i], 100)
    pdf = norm.pdf(x_axis, pow_samples_mean_normed[i], 1*pow_samples_std_normed[i])
    ax.fill_betweenx(x_axis, pdf/scale + (i+1), (i+1), color='b', alpha=0.2)


pow_samples_mean_normed = ADpow_realKE_5_8_x_up_D_1_0/ADpow_realKE_5_8_x_up_D_1_0[0]
for i in range(1,6):
    N, scale = 10, 10
    x_axis = np.linspace(pow_samples_mean_normed[i]+N*pow_samples_std_normed[i], 
                         pow_samples_mean_normed[i]-N*pow_samples_std_normed[i], 100)
    pdf = norm.pdf(x_axis, pow_samples_mean_normed[i], 1.2*pow_samples_std_normed[i])
    ax.fill_betweenx(x_axis, pdf/scale + (i+1), (i+1), color=DETClr, alpha=0.2)


fig.savefig(DATA+casePngDir+'/power_MC_PCE_dist.png', dpi=150, bbox_inches='tight')

