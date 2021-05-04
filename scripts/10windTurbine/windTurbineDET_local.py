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
casePngDir = '/' + caseName+'/png'

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
vtkFile   = 'VTK/'+caseName+'_2000.vtk'
mesh      = pv.UnstructuredGrid(cwd + vtkFile)
nCells    = mesh.n_cells
cCenter   = np.array(mesh.cell_centers().points)

shiftedcCenter = cCenter[:,1:3] - ADloc[1:3]
cIdxInDisk     = np.where(np.linalg.norm(shiftedcCenter, axis=1) <= D/2)[0]
cCenterInDisk  = cCenter[cIdxInDisk]

x_D = cCenter[0:nx,0]/D

# <codecell> Importing LES and RANS
U = mesh.cell_arrays['U'][cIdxInDisk]
k = mesh.cell_arrays['k'][cIdxInDisk]

# <codecell> Computing averaged fields
Umag = np.linalg.norm(U, axis=1)
defU = (Uref-Umag)/Uref
TI   = np.sqrt(k*2/3)/Uref *100

defUAvgd = np.zeros(nx)
TIAvgd   = np.zeros(nx)

# For structured / unstructured mesh
# for i in range(nx): # useful for unstructured mesh
#     cellIdxAtXPlane = np.where(cCenterInDisk[:,0] == cCenterInDisk[i,0])
#     defUAvgd[i] = np.average(defU[cellIdxAtXPlane])
#     TIAvgd[i]   = np.average(TI[cellIdxAtXPlane])

# For structured mesh
defUAvgd = np.average(np.reshape(defU, (int(cIdxInDisk.shape[0]/nx), nx)), axis=0)
TIAvgd   = np.average(np.reshape(TI,   (int(cIdxInDisk.shape[0]/nx), nx)), axis=0)

# <codecell> Plotting averaged fields
rcParams.update({'figure.figsize': [15,7]})

fig = plt.figure()
ax  = (fig.add_subplot(211), fig.add_subplot(212))

axIdx = 0
ax[axIdx].plot(x_D, defUAvgd, 'k')
ax[axIdx].axvline(ADloc[0]/D-Wdx/D/2, color='grey')
ax[axIdx].axvline(ADloc[0]/D+Wdx/D/2, color='grey')
ax[axIdx].set_ylim(bottom=-0.1,top=0.7)
ax[axIdx].set_xlim(left=-4,right=25)
ax[axIdx].set_xlabel('$x/D$')
ax[axIdx].set_ylabel('$\Delta u / u_h$')
ax[axIdx].spines['right'].set_visible(False)
ax[axIdx].spines['top'].set_visible(False)

axIdx = 1
ax[axIdx].plot(x_D, TIAvgd, 'k')
ax[axIdx].axvline(ADloc[0]/D-Wdx/D/2, color='grey')
ax[axIdx].axvline(ADloc[0]/D+Wdx/D/2, color='grey')
ax[axIdx].set_ylim(bottom=2.5,top=17.5)
ax[axIdx].set_xlim(left=-4,right=25)
ax[axIdx].set_xlabel('$x/D$')
ax[axIdx].set_ylabel('$I \\%$')
ax[axIdx].spines['right'].set_visible(False)
ax[axIdx].spines['top'].set_visible(False)

fig.savefig(DATA+casePngDir+'/diskAvgdFields.png', dpi=300)

# <codecell>
