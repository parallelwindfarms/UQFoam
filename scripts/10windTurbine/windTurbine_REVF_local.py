#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 29 15:21:21 2021

@author: jigar
"""

# <codecell> Import Required Modules and env vars
from __future__ import print_function

import os, sys
import numpy as np
import pyvista as pv
import matplotlib.pyplot as plt

SCRIPTS = os.environ['SCRIPTS']
DATA    = os.environ['windTurbineData']

sys.path.append(SCRIPTS+'/10windTurbine')

cwd = os.getcwd() + "/"

# <codecell> Plot settings
# Font, Size
from matplotlib import rcParams
allFontSize = 15
params = {
   'axes.labelsize': allFontSize,
   'legend.fontsize': allFontSize,
   'xtick.labelsize': allFontSize,
   'ytick.labelsize': allFontSize,
   'text.usetex': True,
   'xtick.direction':'in',
   'ytick.direction':'in',
   'font.family':'serif',
   'font.serif': ['Computer Modern Roman']
   }
rcParams.update(params)

# <codecell> Case details
caseName   = '4_Vestas2MWV80_refineLocal_kOmSST_25Dx8Dx8D_WdxD_080_WdyzD_20_varScl'
casePngDir = '/1REVF/' + caseName+'/png'

if (os.path.exists(DATA+casePngDir))==False:
    print('Making new case data directory...')
    os.makedirs(DATA+casePngDir, exist_ok=True)
    
# <codecell> Parametrs
D = 80
h = 70
Uref = 8

# center of AD
ADloc = (0,0,h)

# mesh dependent
nx = 200
Wdx = D/20

# <codecell> Read OpenFOAM case mesh access to grid and cell values
vtkFile   = 'VTK/'+caseName+'_1000.vtk'
mesh      = pv.UnstructuredGrid(cwd + vtkFile)
nCells    = mesh.n_cells
cCenter   = np.array(mesh.cell_centers().points)

shiftedcCenter = cCenter[:,1:3] - ADloc[1:3]
cIdxInDisk     = np.where(np.linalg.norm(shiftedcCenter, axis=1) <= D/2)[0]
cCenterInDisk  = cCenter[cIdxInDisk]

x_D = cCenter[0:nx,0]/D

nCellsInDisk = int(cIdxInDisk.shape[0]/nx)

# <codecell> Importing LES, RANS, UQRANS data
U      = mesh.cell_arrays['U'][cIdxInDisk]
U0     = mesh.cell_arrays['U0'][cIdxInDisk]
USigma = mesh.cell_arrays['USigma'][cIdxInDisk]

nut0     = mesh.cell_arrays['nut0'][cIdxInDisk]
nutSigma = mesh.cell_arrays['nutSigma'][cIdxInDisk]

k      = mesh.cell_arrays['k'][cIdxInDisk]
k0     = mesh.cell_arrays['k0'][cIdxInDisk]
kSigma = k*nutSigma/nut0

# <codecell> Computing averaged U fields
Umag      = np.linalg.norm(U, axis=1)
U0mag     = np.linalg.norm(U0, axis=1)
USigmaMag = np.linalg.norm(USigma, axis=1)
defU      = (Uref-Umag)/Uref
defU0     = (Uref-U0mag)/Uref

defUAvgd  = np.average(np.reshape(defU, (nCellsInDisk, nx)), axis=0)
defU0Avgd = np.average(np.reshape(defU0, (nCellsInDisk, nx)), axis=0)
USigmaMagAvgd = np.average(np.reshape(USigmaMag, (nCellsInDisk, nx)), axis=0)

# <codecell> Computing averaged TI fields
TI  = np.sqrt(k*2/3)/Uref *100
TI0 = np.sqrt(k0*2/3)/Uref *100

N = 2
TIplsNTISigma = np.sqrt((k0+N*kSigma)*2/3)/Uref *100
TIminNTISigma = np.sqrt((k0-N*kSigma)*2/3)/Uref *100

TIAvgd  = np.average(np.reshape(TI, (nCellsInDisk, nx)), axis=0)
TI0Avgd = np.average(np.reshape(TI0, (nCellsInDisk, nx)), axis=0)
TIplsNTISigmaAvgd = np.average(np.reshape(TIplsNTISigma, (nCellsInDisk, nx)), axis=0)
TIminNTISigmaAvgd = np.average(np.reshape(TIminNTISigma, (nCellsInDisk, nx)), axis=0)

# <codecell> Figures settings
rcParams.update({'figure.figsize': [15,7]})

DETClr  = (0.6350, 0.0780, 0.1840)
meanClr = 'b'
fillClr = 'b'
DNSClr  = 'k'
diskClr = 'grey'

# Plotting averaged fields

fig = plt.figure()
ax  = (fig.add_subplot(211), fig.add_subplot(212))

axIdx = 0
ax[axIdx].plot(x_D, defUAvgd, color=DETClr)
ax[axIdx].plot(x_D, defU0Avgd, color=meanClr)
ax[axIdx].fill_between(x_D, defU0Avgd + N*USigmaMagAvgd/Uref, \
                            defU0Avgd - N*USigmaMagAvgd/Uref, \
                            alpha=0.2, linewidth=0, color='b')
ax[axIdx].axvspan(ADloc[0]/D-Wdx/D, ADloc[0]/D+Wdx/D, linewidth=0, color=diskClr)
ax[axIdx].set_ylim(bottom=-0.1,top=0.7)
ax[axIdx].set_xlim(left=-2,right=23)
ax[axIdx].set_xlabel('$x/D$')
ax[axIdx].set_ylabel('$\Delta u / u_h$')
ax[axIdx].spines['right'].set_visible(False)
ax[axIdx].spines['top'].set_visible(False)

axIdx = 1
ax[axIdx].plot(x_D, TIAvgd, color=DETClr)
ax[axIdx].plot(x_D, TI0Avgd, color=meanClr)
ax[axIdx].fill_between(x_D, TIplsNTISigmaAvgd, TIminNTISigmaAvgd, \
                       alpha=0.2, linewidth=0, color='b')
ax[axIdx].axvspan(ADloc[0]/D-Wdx/D, ADloc[0]/D+Wdx/D, linewidth=0, color=diskClr)

ax[axIdx].set_ylim(bottom=2.5,top=17.5)
ax[axIdx].set_xlim(left=-2,right=23)
ax[axIdx].set_xlabel('$x/D$')
ax[axIdx].set_ylabel('$I \\%$')
ax[axIdx].spines['right'].set_visible(False)
ax[axIdx].spines['top'].set_visible(False)

# fig.savefig(DATA+casePngDir+'/diskAvgdFields.png', dpi=300)

# <codecell>
