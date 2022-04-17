#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  3 11:09:02 2021

@author: jigar
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
    
myUQlib.rcParamsSettings(usetex=True)


# %% Data
fillUpto = 8

# defU Sobol Idx
S_T_u_h_defU = np.array([0.5, 0.6, 0.8, 0.2, 0.2, 0.3])/100
S_T_I_h_defU = np.array([32, 30, 26, 22, 17, 16])/100

S_T_RRSTF_defU = np.zeros((6, 18))
S_T_RRSTF_defU[0,:fillUpto] = np.array([13, 21, 7, 5, 6, 9, 0.5, 0.8])/100
S_T_RRSTF_defU[1,:fillUpto] = np.array([18, 20, 6, 4, 8, 7, 1.1, 0.2])/100
S_T_RRSTF_defU[2,:fillUpto] = np.array([21, 23, 7, 3, 5, 7, 1.5, 0.8])/100
S_T_RRSTF_defU[3,:fillUpto] = np.array([19, 25, 5, 8, 4, 10, 0.5, 0.8])/100
S_T_RRSTF_defU[4,:fillUpto] = np.array([22, 27, 7, 5, 8, 8, 0.5, 0.8])/100
S_T_RRSTF_defU[5,:fillUpto] = np.array([21, 28, 9, 4.5, 7, 8, 0.5, 0.8])/100

# TI Sobol Idx
S_T_u_h_TI = np.array([0.6, 1.0, 0.5, 0.3, 0.1, 1])/100
S_T_I_h_TI = np.array([38, 25, 32, 31, 33, 32])/100

S_T_RRSTF_TI = np.zeros((6, 18))
S_T_RRSTF_TI[0,:fillUpto] = np.array([10, 18, 5, 8, 3, 11, 0.7, 0.5])/100
S_T_RRSTF_TI[1,:fillUpto] = np.array([21, 23, 7, 2, 6, 7, 1.5, 0.8])/100
S_T_RRSTF_TI[2,:fillUpto] = np.array([20, 24, 5, 3, 2, 8, 0.4, 0.3])/100
S_T_RRSTF_TI[3,:fillUpto] = np.array([18, 25, 3.5, 3, 5, 8, 0.5, 0.8])/100

S_T_RRSTF_TI[4,:fillUpto] = np.array([16, 20, 6, 3, 8, 7, 1.1, 0.2])/100
S_T_RRSTF_TI[5,:fillUpto] = np.array([17, 23, 5, 2, 4, 10, 0.5, 0.8])/100


# %% Figure 
myUQlib.rcParamsSettings(22.5)

DETClr = (0.6350, 0.0780, 0.1840)
uhClr = 'g'
IhClr = DETClr 
RRSTFClr  = 'b' 
meanClr = 'b'
barWidth = 0.7
alpha = 0.4
figsize = (15,6)

fig, ax = plt.subplots(ncols=6, nrows=2, constrained_layout=True, 
                       sharex=True, sharey=True, figsize=figsize)

spline = lambda x=False: ax[axIdx].spines[['right', 'top']].set_visible(x)
xticks = lambda x=False: ax[axIdx].tick_params(bottom=x)

bar_RRSTF = np.arange(fillUpto) + 2 + 1

for WT in range(1,7):
    ## defU ANOVA
    # u_h, I_h
    axIdx = (0,WT-1)
    ax[axIdx].bar(1, S_T_u_h_defU[WT-1], color=uhClr, width = barWidth, alpha=alpha, edgecolor='grey' )
    ax[axIdx].bar(2, S_T_I_h_defU[WT-1], color=IhClr, width = barWidth, alpha=alpha, edgecolor='grey' )
    
    # RRSTF
    ax[axIdx].bar(bar_RRSTF, S_T_RRSTF_defU[WT-1,:fillUpto], color = RRSTFClr, 
                  width = barWidth, alpha=alpha, edgecolor='grey')
    
    ax[axIdx].set_title('WT ' + str(WT), pad=20)
    if WT==1: 
        ax[axIdx].set_ylabel('Sobol Index \n for $\Delta u / u_h$')
        ax[axIdx].set_ylim([0, 0.4])
    spline(False)
    xticks(False)
    
    ## TI ANOVA
    # u_h, I_h
    axIdx = (1,WT-1)
    ax[axIdx].bar(1, S_T_u_h_TI[WT-1], color=uhClr, width = barWidth, alpha=alpha, edgecolor='grey' )
    ax[axIdx].bar(2, S_T_I_h_TI[WT-1], color=IhClr, width = barWidth, alpha=alpha, edgecolor='grey' )
    
    # RRSTF
    ax[axIdx].bar(bar_RRSTF, S_T_RRSTF_TI[WT-1,:fillUpto], color = RRSTFClr, 
                  width = barWidth, alpha=alpha, edgecolor='grey')
    
    ax[axIdx].set_xlabel('Mode')
    if WT==1:
        ax[axIdx].set_ylabel('Sobol Index \n for $I$')
        ax[axIdx].set_ylim([0, 0.4])
        ax[axIdx].set_xticks(np.arange(len(bar_RRSTF)+2)+1)
    xticks(False)
    spline(False)

leg = ['$S_T|_{u_{h}}$', '$S_T|_{I_{h}}$', '$S_{T_i}|_{RRSTF}$']
fig.legend(leg, ncol=3, frameon=False, columnspacing=2, 
           loc='lower center', bbox_to_anchor=(0.52, 1.0))

fig.savefig(DATA+casePngDir+'/SobolIndices_defU_TI.png', dpi=150, bbox_inches='tight')
