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
import chaospy as cp
import pyvista as pv
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
D = 80
h = 70
Uref = 8

# center of AD
ADloc = (0,0,h)

# mesh dependent
nx = 150
Wdx = D/20

# <codecell> Read OpenFOAM case mesh access to grid and cell values
vtkFile = 'baseCase/VTK/baseCase_1000.vtk'
mesh    = pv.UnstructuredGrid(cwd + vtkFile)
nCells  = mesh.n_cells
cCenter = np.array(mesh.cell_centers().points)

cIdxInDuct = np.where((abs(cCenter[:,1]-ADloc[1]) <= D*1.30) & (cCenter[:,2] <= D*2.6))[0]
cCenterInDuct = cCenter[cIdxInDuct]

lines = (2,4,7,12)
numLines = len(lines)
cIdxInPlanes = []
for l in lines:
    tmpIdxInPlane = np.where((cCenter[:,0] >= l*D) & (cCenter[:,0] <= l*D+l*Wdx))[0]
    cIdxInPlanes.append(np.intersect1d(cIdxInDuct, np.array(tmpIdxInPlane)))
cIdxInPlanes = np.array(cIdxInPlanes)
cCenterInPlanes = cCenter[cIdxInPlanes]
nCells = len(cIdxInPlanes[0])
y, z = cCenterInPlanes[0,:,1]/D, cCenterInPlanes[0,:,2]/D

RDet = mesh.cell_arrays['turbulenceProperties:R'][cIdxInPlanes]
RDet = np.array([myUQlib.symmTensorToTensorv2012(RDet[l], nCells) for l in range(numLines)])
tkeDet = mesh.cell_arrays['k'][cIdxInPlanes]
a_ij_Det = np.array([myUQlib.anisotropyTensor(RDet[l], tkeDet[l], nCells) for l in range(numLines)])

# Random sample
s = 25
vtkFile = 'sample_'+str(s)+'/VTK/'+'sample_'+str(s)+'_1000.vtk'
mesh = pv.UnstructuredGrid(cwd + vtkFile)
RSample = mesh.cell_arrays['turbulenceProperties:R'][cIdxInPlanes]
RSample = np.array([myUQlib.symmTensorToTensorv2012(RSample[l], nCells) for l in range(numLines)])
tkeSample = mesh.cell_arrays['k'][cIdxInPlanes]
a_ij_Sample = np.array([myUQlib.anisotropyTensor(RSample[l], tkeSample[l], nCells) for l in range(numLines)])

# LES
s = 85
vtkFile = 'sample_'+str(s)+'/VTK/'+'sample_'+str(s)+'_1000.vtk'
mesh = pv.UnstructuredGrid(cwd + vtkFile)
RLES = mesh.cell_arrays['turbulenceProperties:R'][cIdxInPlanes]
RLES = np.array([myUQlib.symmTensorToTensorv2012(RLES[l], nCells) for l in range(numLines)])
tkeLES = mesh.cell_arrays['k'][cIdxInPlanes]
a_ij_LES = np.array([myUQlib.anisotropyTensor(RLES[l], tkeLES[l], nCells) for l in range(numLines)])


# <codecell> Eigen decomposition of a_ij
eValDet = np.zeros((numLines,nCells,3))
eVecDet = np.zeros((numLines,nCells,3,3))
for l in range(numLines):
    eValDet[l], eVecDet[l] = myUQlib.eigenDecomposition(a_ij_Det[l], nCells)

eValSample = np.zeros((numLines,nCells,3))
eVecSample = np.zeros((numLines,nCells,3,3))
for l in range(numLines):
    eValSample[l], eVecSample[l] = myUQlib.eigenDecomposition(a_ij_Sample[l], nCells)

eValLES = np.zeros((numLines,nCells,3))
eVecLES = np.zeros((numLines,nCells,3,3))
for l in range(numLines):
    eValLES[l], eVecLES[l] = myUQlib.eigenDecomposition(a_ij_LES[l], nCells)

# <codecell> Compute Barycentric Coordinates and its Viridis color values
C_Det = np.array([myUQlib.baryCentricCoordinates(eValDet[l]) for l in range(numLines)])
clr_C_Det = myUQlib.bccViridisCmap(C_Det)

C_Sample = np.array([myUQlib.baryCentricCoordinates(eValSample[l]) for l in range(numLines)])
clr_C_Sample = myUQlib.bccViridisCmap(C_Sample)

C_LES = np.array([myUQlib.baryCentricCoordinates(eValLES[l]) for l in range(numLines)])
clr_C_LES = myUQlib.bccViridisCmap(C_LES)

# <codecell> plotBaryCentricCoordinateSystemWithCmap
fig, ax = plt.subplots(ncols=1, nrows=1)
myUQlib.plotBaryCentricCoordinateSystemWithCmap(ax)
fig.savefig(DATA+casePngDir+'/baryCentricCoordinateSystemWithCmap.png', dpi=300, bbox_inches='tight')

# <codecell> Figures settings
myUQlib.rcParamsSettings(usetex=True)

nrows = 2#3

fig, ax = plt.subplots(ncols=numLines, nrows=nrows, constrained_layout=True,
                       figsize=(10,nrows*2.5), sharex=True, sharey=True)

yNr = 37 # Need some generic expression!!

row = 0
for l in range(numLines):
    axIdx = (row,l)

    c = clr_C_Det[l].reshape(yNr,int(C_Det[l].shape[0]/yNr),-1)
    ax[axIdx].imshow(c, interpolation='bilinear', origin='lower',
                      extent=[y.min(),y.max(),z.min(),z.max()])

    myUQlib.drawCircle(ax[axIdx], ADloc[1]/D, ADloc[2]/D, 1/2)
    ax[axIdx].set_xticks([-1.0, 0.0, 1.0])
    ax[axIdx].set_yticks([0.5, 1.5, 2.5])
    if l==0: ax[axIdx].set_ylabel('$z/D$')
    myUQlib.setSpines(ax[axIdx])
    myUQlib.setTiks(ax[axIdx])

    ax[axIdx].set_title('$x/D=$ ' + str(lines[l]), pad=20)

row = 1
for l in range(numLines):
    axIdx = (row,l)

    c=clr_C_Sample[l].reshape(yNr,int(C_Sample[l].shape[0]/yNr),-1)
    ax[axIdx].imshow(c, interpolation='bilinear', origin='lower',
                      extent=[y.min(),y.max(),z.min(),z.max()])

    myUQlib.drawCircle(ax[axIdx], ADloc[1]/D, ADloc[2]/D, 1/2)
    if l==0: ax[axIdx].set_ylabel('$z/D$')
    myUQlib.setSpines(ax[axIdx])
    myUQlib.setTiks(ax[axIdx])


# row = 2
# for l in range(numLines):
#     axIdx = (row,l)

#     c=clr_C_LES[l].reshape(yNr,int(C_LES[l].shape[0]/yNr),-1)
#     ax[axIdx].imshow(c, interpolation='bilinear', origin='lower',
#                       extent=[y.min(),y.max(),z.min(),z.max()])

#     myUQlib.drawCircle(ax[axIdx], ADloc[1]/D, ADloc[2]/D, 1/2)
#     ax[axIdx].set_xlabel('$y/D$')
#     if l==0: ax[axIdx].set_ylabel('$z/D$')
#     myUQlib.setSpines(ax[axIdx])
#     myUQlib.setTiks(ax[axIdx])

fig.savefig(DATA+casePngDir+'/anisoTensorContour.png', dpi=300, bbox_inches='tight')






















