#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 10 16:32:43 2020

@author: jigar
"""


# -*- coding: utf-8 -*-
"""
Purpose of this python script -
(1) Create a KL Expansion of a Gaussian field/process
(2) Compute the PC Expansion of a Random SPD Matrix field
(3) Compute inner products of polynomials
"""

# <codecell> Import Required Modules and env vars
#from __future__ import print_function

import numpy as np
import pyvista as pv
import os, sys
import chaospy as cp
import openturns as ot

SCRIPTS = os.environ['SCRIPTS']
sys.path.append(SCRIPTS+'/pythonScripts')

import myUQlib   # Import myUQlib
cwd = os.getcwd() + "/"


# <codecell> Read OpenFOAM case mesh acces to grid and cell values
phyDim  = 2

vtkFile = 'VTK/1nutKLEgPC_lx1.5H_ly0.5H_var0.50_M100x80_iter1_NstdFromEnKF_5000.vtk'
mesh    = pv.UnstructuredGrid(cwd + vtkFile)
nCells  = mesh.n_cells
cCenter = mesh.cell_centers().points
cBounds = mesh.cell_centers().bounds
nPoints = mesh.n_points
mBounds = mesh.bounds
mPoints = mesh.points
sized   = mesh.compute_cell_sizes()
cVols   = sized.cell_arrays["Volume"]
mVol    = mesh.volume
volAvg  = cVols/mVol

if (phyDim==2):

    Nx = 100; Ny = 80;

    myIndices = [Nx-1, Ny-1]
    mesher    = ot.IntervalMesher(myIndices)
    kleBounds = ot.Interval([cBounds[0], cBounds[2]], 
                            [cBounds[1], cBounds[3]])
    kleMesh = mesher.build(kleBounds)
    kleMesh.setVertices(cCenter[:,[0,1]])

if (phyDim==3):
        
    Nx = 100; Ny = 4; Nz = 5;

    myIndices = [Nx-1, Ny-1, Nz-1]
    mesher    = ot.IntervalMesher(myIndices)
    kleBounds = ot.Interval([cBounds[0], cBounds[2], cBounds[4]], 
                            [cBounds[1], cBounds[3], cBounds[5]])
    kleMesh = mesher.build(kleBounds)
    kleMesh.setVertices(cCenter[:,[0,1,2]])
    

# <codecell> Observation locations and matrix
h  = 0.028
cCenter_h = cCenter/h

      #0.05,0.5,1, 2, 3, 4, 5, 6, 7, 8
xByH_idx = [0,5,12,24,34,44,55,65,76,87]
numLines = len(xByH_idx);
yByH_at_xByH = np.zeros([Ny,numLines])
H_map = np.zeros([Ny*numLines,nCells])
H = np.zeros([Ny*numLines,nCells])
tol = 5e-2
yTurn_1 = 1.30
yTurn_2 = 1.30

for j in range(numLines):
    for i in range(Ny):
        yByH_at_xByH[i,j] = cCenter_h[i*Nx+xByH_idx[j],1]
        H_map[j+i*numLines,i*Nx+xByH_idx[j]] = 1
        # Setting Observation MAtrix at selected locaitons only
        if (j==2 or j==3 or j==8 or j==9) and (i==9 or i==19 or i==29):
            H[j+i*numLines,i*Nx+xByH_idx[j]] = 1
        if (j>=2 and j<=7) and (abs(yByH_at_xByH[i,j]-yTurn_1)<tol or abs(yByH_at_xByH[i,j]-yTurn_2)<tol):
            H[j+i*numLines,i*Nx+xByH_idx[j]] = 1
            
# <codecell> Load data
        
# U_DNS        
pathToDNS = os.environ['periodicHillData'] + '/Re2800/0DNS/'

U_DNS = []
for i in range(1,numLines):
    U_DNS.append(np.loadtxt(pathToDNS + 'DNS0'+str(i)+'.dat')[:,:2])
U_DNS.append(np.loadtxt(pathToDNS + 'DNS'+str(numLines)+'.dat')[:,:2])

# U_UQ
numModes = 55
U_UQ = np.zeros([numModes, nCells])
Std_UQ = np.zeros(nCells)
for i in range(numModes):
    U_UQ[i] = mesh.cell_arrays['U'+str(i)][:,0]
    if i>0:
        Std_UQ += U_UQ[i]**2
Mean_UQ = U_UQ[0]
Std_UQ = Std_UQ**0.5


# <codecell> Map DNS data to RANS mesh

Y_U_lineVal = np.zeros([Ny, numLines])
for i in range(0,numLines):
    x_arr = U_DNS[i][:,0]
    y_arr = U_DNS[i][:,1]
    x_val = yByH_at_xByH[:,i]
    Y_U_lineVal[:,i] = np.interp(x_val, x_arr[::-1], y_arr[::-1])

Y_U = np.zeros(Ny*numLines)
for i in range(Ny):
    Y_U[i*numLines:(i*numLines+numLines)] = Y_U_lineVal[i]


# <codecell> Sampling from U_UQ

d = 18
n = 3
q = 0.1

sampleSize = 1000

dist      = cp.Iid(cp.Normal(0,1), d)
xiSamples = (dist.sample(size=sampleSize, rule="L")).T

polyExpsn = cp.orth_ttr(n, dist, retall=False, normed=False, cross_truncation=1/q)

USamples    = np.zeros([sampleSize, nCells])
for i in range(sampleSize):
    USamples[i]    = np.matmul(U_UQ.T, np.array(polyExpsn(*(xiSamples[i])))[0:numModes])
USampleMean    = np.mean(USamples,0)
HSt = np.matmul(H, USamples.T)

print('Relative Error in means from UQ and Sampling = ', np.max(abs(USampleMean-Mean_UQ)/Mean_UQ))
print('Max of xiSamples = ', np.max(xiSamples))
print('Min of xiSamples = ', np.min(xiSamples))
print('Mean of xiSamples = ', xiSamples.mean(axis=0))
print('Std of xiSamples = ', xiSamples.std(axis=0))


# <codecell> Analysis step - Kalman gain matrix (K)
DNSnoise = abs(0.05*Y_U)
R = np.diag(DNSnoise*DNSnoise)

xi_p  = xiSamples - xiSamples.mean(axis=0,keepdims=1)
HSt_p = HSt - HSt.mean(axis=1,keepdims=1)

PHt  = 1/(1-sampleSize) * np.matmul(HSt_p,xi_p)
HPHt = 1/(1-sampleSize) * np.matmul(HSt_p,HSt_p.T)

K = PHt.T.dot(np.linalg.inv(HPHt+R))

print('Max of K = ', np.max(K))
print('Min of K = ', np.min(K))

# <codecell> Analysis step - updating ensemble

xiSamples_Updated = 0*xiSamples
for i in range(sampleSize):
    diff = Y_U - HSt.T[i]
    xiSamples_Updated[i] = xiSamples[i] + np.matmul(K,diff)

print('Max of xiSamples_Updated = ', np.max(xiSamples_Updated))
print('Max of xiSamples_Updated = ', np.min(xiSamples_Updated))
xiSamples_Updated_mean = xiSamples_Updated.mean(axis=0)
xiSamples_Updated_std = xiSamples_Updated.std(axis=0)
print('Mean of xiSamples_Updated = ', xiSamples_Updated_mean)
print('Std of xiSamples_Updated = ', xiSamples_Updated_std)


# <codecell> Updated REVF (run expGp_gPCKLE_logNProc.py (without dist) after this cell)

dist = cp.Normal(xiSamples_Updated_mean[0], xiSamples_Updated_std[0])
# dist = cp.Normal(0, xiSamples_Updated_std[0])

for i in range(1,d):
    ith_dist = cp.Normal(xiSamples_Updated_mean[i], xiSamples_Updated_std[i])
    # ith_dist = cp.Normal(0, xiSamples_Updated_std[i])
    dist = cp.J(dist, ith_dist)


# <codecell> Setting UQ info. for OpenFOAM simulation
#myUQlib.uqInfo(n-1, d, q, "REVF", dist)

