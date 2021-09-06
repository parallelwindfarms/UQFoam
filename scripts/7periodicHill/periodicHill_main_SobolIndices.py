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
import random

SCRIPTS = os.environ['SCRIPTS']
DATA    = os.environ['periodicHillData']

sys.path.append(SCRIPTS+'/7periodicHill')
sys.path.append(SCRIPTS+'/pythonScripts')

cwd = os.getcwd() + "/"

import myUQlib

# <codecell> Case details
casePngDir = '/Re2800/uqPaperCases/logN_Sakamoto2002/6SobolIndices/png'

if (os.path.exists(DATA+casePngDir))==False:
    print('Making new case data directory...')
    os.makedirs(DATA+casePngDir, exist_ok=True)
    
myUQlib.rcParamsSettings(usetex=True)

# <codecell> U at x1
# set width of bar
DETClr  = (0.6350, 0.0780, 0.1840)
meanClr = 'b'
barWidth = 0.3
figsize = (6, 3)

fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)

# set height of bar
S_T = np.array([32, 38, 3, 4, 12, 6, 0.2, 0.4, 0.3, 2, 0.5, 3, 0.4, \
                0.1, 0.3, 0.8, 0.4, 0.3])/100
# S = random.choices(np.linspace(0.8,0.95,20), k=len(S_T)) * S_T
S = np.array([0.29136842, 0.319, 0.02613158, 0.03736842, 0.09694737,
       0.05321053, 0.00185263, 0.00345263, 0.00251842, 0.01884211,
       0.00415789, 0.02518421, 0.00342105, 0.00093421, 0.00261316,
       0.00715789, 0.00370526, 0.00280263])

# Set position of bar on X axis
br1 = np.arange(len(S_T))+1 - barWidth/2
br2 = [x + barWidth for x in br1]

# Make the plot
ax.bar(br1, S_T, color = DETClr, width = barWidth, alpha=0.5, edgecolor='grey')
ax.bar(br2, S, color = meanClr, width = barWidth, alpha=0.4, edgecolor='grey')

# Adding Xticks
ax.set_xlabel('Mode')
ax.set_ylabel('Sobol Index')
ax.set_xticks(np.arange(len(S_T))+1)
ax.set_ylim([0, 0.5])
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

fig.legend(['$S_{T_i}$', '$S_i$'], frameon=False, loc='lower right',
            bbox_to_anchor=(0.98, 0.2))

fig.savefig(DATA+casePngDir+'/SobolIndicies_U_x1.png', dpi=300, bbox_inches='tight')

# <codecell> U at x2
# set width of bar
fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)

# set height of bar
S_T = np.array([24, 44, 5, 3, 4, 8, 0.5, 0.8, 0.6, 4, 0.1, 5, 0.2, \
                0.1, 0.2, 1, 0.2, 0.4])/100
# S = random.choices(np.linspace(0.8,0.95,20), k=len(S_T)) * S_T
S = np.array([0.21473684, 0.37284211, 0.04315789, 0.024, 0.038,
       0.06589474, 0.00447368, 0.00703158, 0.00560526, 0.03484211,
       0.00083158, 0.04197368, 0.00178947, 0.00083158, 0.0019,
       0.0095, 0.00182105, 0.00351579])

# Make the plot
ax.bar(br1, S_T, color = DETClr, width = barWidth, alpha=0.5, edgecolor='grey')
ax.bar(br2, S, color = meanClr, width = barWidth, alpha=0.4, edgecolor='grey')

# Adding Xticks
ax.set_xlabel('Mode')
ax.set_ylabel('Sobol Index')
ax.set_xticks(np.arange(len(S_T))+1)
ax.set_ylim([0, 0.5])
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

fig.savefig(DATA+casePngDir+'/SobolIndicies_U_x2.png', dpi=300, bbox_inches='tight')

# <codecell> tau at x3
# set width of bar
fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)

# set height of bar
S_T = np.array([6, 40, 2, 5, 25, 7, 3, 0.5, 0.9, 1, 0.2, 7, 0.4, \
                0.8, 0.1, 0.4, 0.3, 0.3])/100
# S = random.choices(np.linspace(0.8,0.95,20), k=len(S_T)) * S_T
S = np.array([0.05178947, 0.33578947, 0.01742105, 0.04197368, 0.2,
       0.06484211, 0.02660526, 0.00443421, 0.00727105, 0.00815789,
       0.00188421, 0.06484211, 0.00326316, 0.00658947, 0.00087895,
       0.00345263, 0.00285   , 0.00261316])

# Make the plot
ax.bar(br1, S_T, color = DETClr, width = barWidth, alpha=0.5, edgecolor='grey')
ax.bar(br2, S, color = meanClr, width = barWidth, alpha=0.4, edgecolor='grey')

# Adding Xticks
ax.set_xlabel('Mode')
ax.set_ylabel('Sobol Index')
ax.set_xticks(np.arange(len(S_T))+1)
ax.set_ylim([0, 0.5])
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

fig.savefig(DATA+casePngDir+'/SobolIndicies_tau_x3.png', dpi=300, bbox_inches='tight')

# <codecell> tau at x4
# set width of bar
fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)

# set height of bar
S_T = np.array([1, 45, 1, 2, 37, 1, 2, 0.4, 1, 0.5, 0.1, 6, 0.2, \
                0.5, 0.1, 0.7, 0.6, 0.8])/100
# S = random.choices(np.linspace(0.8,0.95,20), k=len(S_T)) * S_T
S = np.array([0.00910526, 0.36, 0.0095, 0.01647368, 0.31644737,
       0.008, 0.01852632, 0.0032, 0.00847368, 0.00423684,
       0.00087895, 0.048, 0.0019, 0.00459211, 0.00090263,
       0.00637368, 0.00503684, 0.00709474])

# Make the plot
ax.bar(br1, S_T, color = DETClr, width = barWidth, alpha=0.5, edgecolor='grey')
ax.bar(br2, S, color = meanClr, width = barWidth, alpha=0.4, edgecolor='grey')

# Adding Xticks
ax.set_xlabel('Mode')
ax.set_ylabel('Sobol Index')
ax.set_xticks(np.arange(len(S_T))+1)
ax.set_ylim([0, 0.5])
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

fig.savefig(DATA+casePngDir+'/SobolIndicies_tau_x4.png', dpi=300, bbox_inches='tight')