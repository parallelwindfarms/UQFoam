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
yTurn_2 = 2.00

for j in range(numLines):
    for i in range(Ny):
        yByH_at_xByH[i,j] = cCenter_h[i*Nx+xByH_idx[j],1]
        H_map[j+i*numLines,i*Nx+xByH_idx[j]] = 1
        # Setting Observation Matrix at selected locaitons only
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
numModes = 55#55
U_UQ = np.zeros([numModes, nCells])
Std_UQ = np.zeros(nCells)
for i in range(numModes):
    U_UQ[i] = mesh.cell_arrays['U'+str(i)][:,0]
    if i>0:
        Std_UQ += U_UQ[i]**2
Std_UQ = Std_UQ**0.5
Mean_UQ = U_UQ[0]

# tauW_DNS
tauW_DNS_x = np.loadtxt(pathToDNS + 'tauW.dat')[:,0]
tauW_DNS   = np.loadtxt(pathToDNS + 'tauW.dat')[:,1]

# tau_UQ
tauW_dataLocation = 'postProcessing/sampleTauW/surface/1/'
tauW_UQ_x = (np.loadtxt(tauW_dataLocation+'tauW0SpAvg_patch_hills.raw')[:,0])/h
tauW_UQ = np.zeros([numModes, Nx])
tauWStd_UQ = np.zeros(Nx)
for i in range(numModes):
    tauW_UQ[i] = -np.loadtxt(tauW_dataLocation+'tauW'+str(i)+'_patch_hills.raw')[:,3]
    if i>0:
        tauWStd_UQ += tauW_UQ[i]**2
tauWStd_UQ = tauWStd_UQ**0.5
tauWMean_UQ = tauW_UQ[0]

H_tauW = np.zeros([Nx,Nx])
H_map_tauW = np.zeros([Nx,Nx])
for i in range(Nx):
    H_map_tauW[i,i] = 1
    if (abs(tauW_UQ_x[i]-1)<tol or abs(tauW_UQ_x[i]-3)<tol or abs(tauW_UQ_x[i]-6)<tol):
    # if (abs(tauW_UQ_x[i]-1)<tol or abs(tauW_UQ_x[i]-7)<tol):
        H_tauW[i,i] = 1 
    

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

Y_tauW = np.interp(tauW_UQ_x,tauW_DNS_x,tauW_DNS)

# <codecell> Sampling from U_UQ

d = 18
n = 3
q = 0.1

sampleSize = 2000#00

dist      = cp.Iid(cp.Normal(0,1), d)
xiSamples = dist.sample(size=sampleSize, rule="L")

polyExpsn = cp.orth_ttr(n, dist, retall=False, normed=False, cross_truncation=1/q)

USamples    = np.zeros([sampleSize, nCells])
tauWSamples = np.zeros([sampleSize, Nx])
for i in range(sampleSize):
    USamples[i]    = np.matmul(U_UQ.T, np.array(polyExpsn(*(xiSamples.T[i])))[0:numModes])
    tauWSamples[i] = np.matmul(tauW_UQ.T, np.array(polyExpsn(*(xiSamples.T[i])))[0:numModes])
USampleMean    = np.mean(USamples,0)
tauWSampleMean = np.mean(tauWSamples,0)

print('Relative Error in means from UQ and Sampling = ', np.max(abs(USampleMean-Mean_UQ)/Mean_UQ))
print('Max of USamples = ', np.max(USamples))
print('Min of USamples = ', np.min(USamples))


# <codecell> Analysis step - Kalman gain matrix (K)

from time import time
t = time()

P = np.zeros([nCells, nCells])

for i in range(1,numModes):
    print('Mode : ', i)
    P += np.outer(U_UQ[i], U_UQ[i])

# tmpP = np.outer(USampleMean, USampleMean)
# for i in range(sampleSize):
#     print('Sample : ', i)
#     P += np.outer(USamples[i],USamples[i])
# P = (P - sampleSize*tmpP)/(sampleSize-1)


# from numba import njit, prange

# @njit(parallel=True)
# def prangeOterPdt(A):
#     s = np.zeros((A.shape[1], A.shape[1]))
#     for i in prange(A.shape[0]):
#         print('Sample : ', i, '\n')
#         s += np.outer(A[i], A[i])
#     return s

# P = (prangeOterPdt(USamples) - sampleSize*tmpP)/(sampleSize-1)


print('Time = ', time()-t, ' s')

# <codecell> U - Gain matrix
DNSnoise = abs(0.10*Y_U)
R = np.diag(DNSnoise*DNSnoise)
K = np.matmul(np.matmul(P,H.T), np.linalg.inv(np.matmul(H,np.matmul(P,H.T)) + R))

print('Max of K = ', np.max(K))
print('Min of K = ', np.min(K))

# <codecell> U - Analysis step - updating ensemble

USamples_Updated = 0*USamples
for i in range(sampleSize):
    diff = Y_U - np.matmul(H,USamples[i])
    USamples_Updated[i] = USamples[i] + np.matmul(K,diff)

print('Max of USamples_Updated = ', np.max(USamples_Updated))
print('Max of USamples_UPdated = ', np.min(USamples_Updated))

# <codecell> U - Evaluate the exansion coefficients

U_UQ_Updated = 0*U_UQ

for j in range(sampleSize):
    tmp = polyExpsn(*(xiSamples.T[j]))
    for i in range(numModes):
        U_UQ_Updated[i] += USamples_Updated[j] * tmp[i]
U_UQ_Updated = U_UQ_Updated/sampleSize 

Std_UQ_Updated = np.zeros(nCells)
for i in range(1,numModes):
    Std_UQ_Updated += U_UQ_Updated[i]**2
Mean_UQ_Updated = U_UQ_Updated[0]
Std_UQ_Updated  = Std_UQ_Updated**0.5


# <codecell> U - Plot
import matplotlib.pyplot as plt

y_val_UQ = np.zeros([Ny,numLines])
Std_val_UQ = np.zeros([Ny,numLines])
y_val_UQ_Updated = np.zeros([Ny,numLines])
Std_val_UQ_Updated = np.zeros([Ny,numLines])

tmp = np.matmul(H_map,Mean_UQ.T)
tmp_Std = np.matmul(H_map,Std_UQ.T)
tmp_Updated = np.matmul(H_map,Mean_UQ_Updated.T)
tmp_Std_Updated = np.matmul(H_map,Std_UQ_Updated.T)

scale = 2

for i in range(2,numLines):
   
    plt.plot(scale*Y_U_lineVal[:,i]+1*i-1,yByH_at_xByH[:,i],'k-')
    for j in range(Ny):
        y_val_UQ[j,i] = tmp[i+j*numLines] 
        Std_val_UQ[j,i] = tmp_Std[i+j*numLines] 
        y_val_UQ_Updated[j,i] = tmp_Updated[i+j*numLines] 
        Std_val_UQ_Updated[j,i] = tmp_Std_Updated[i+j*numLines] 
    
    plt.plot(scale*y_val_UQ[:,i]+1*i-1,yByH_at_xByH[:,i],'b')
    plt.plot(scale*(y_val_UQ[:,i]+2*Std_val_UQ[:,i])+1*i-1,yByH_at_xByH[:,i],'b--')
    plt.plot(scale*(y_val_UQ[:,i]-2*Std_val_UQ[:,i])+1*i-1,yByH_at_xByH[:,i],'b--')

    plt.plot(scale*y_val_UQ_Updated[:,i]+1*i-1,yByH_at_xByH[:,i],'r')
    plt.plot(scale*(y_val_UQ_Updated[:,i]+2*Std_val_UQ_Updated[:,i])+1*i-1,yByH_at_xByH[:,i],'r--')
    plt.plot(scale*(y_val_UQ_Updated[:,i]-2*Std_val_UQ_Updated[:,i])+1*i-1,yByH_at_xByH[:,i],'r--')
        

# <codecell> tauW - Analysis step - Kalman gain matrix (K)

from time import time
t = time()

P_tauW = np.zeros([Nx, Nx])
tmpP = np.outer(tauWSampleMean, tauWSampleMean)
for i in range(sampleSize):
    print('Sample : ', i)
    P_tauW += np.outer(tauWSamples[i],tauWSamples[i])
P_tauW = (P_tauW - sampleSize*tmpP)/(sampleSize-1)

print('Time = ', time()-t, ' s')

# <codecell> tauW - Gain matrix tauW
DNSnoise = abs(0.10*Y_tauW)
R = np.diag(DNSnoise*DNSnoise)
K = np.matmul(np.matmul(P_tauW,H_tauW.T), np.linalg.inv(np.matmul(H_tauW,np.matmul(P_tauW,H_tauW.T)) + R))

print('Max of K = ', np.max(K))
print('Min of K = ', np.min(K))

# <codecell> tauW - Analysis step - updating ensemble

tauWSamples_Updated = 0*tauWSamples
for i in range(sampleSize):
    diff = Y_tauW - np.matmul(H_tauW,tauWSamples[i])
    tauWSamples_Updated[i] = tauWSamples[i] + np.matmul(K,diff)

print('Max of tauWSamples_Updated = ', np.max(tauWSamples_Updated))
print('Max of tauWSamples_Updated = ', np.min(tauWSamples_Updated))

# <codecell> tauW - Evaluate the exansion coefficients

tauW_UQ_Updated = 0*tauW_UQ

for j in range(sampleSize):
    tmp = polyExpsn(*(xiSamples.T[j]))
    for i in range(numModes):
        tauW_UQ_Updated[i] += tauWSamples_Updated[j] * tmp[i]
tauW_UQ_Updated = tauW_UQ_Updated/sampleSize 

tauWStd_UQ_Updated = np.zeros(Nx)
for i in range(1,numModes):
    tauWStd_UQ_Updated += tauW_UQ_Updated[i]**2
tauWMean_UQ_Updated = tauW_UQ_Updated[0]
tauWStd_UQ_Updated  = tauWStd_UQ_Updated**0.5


# <codecell> tauW - Plot
import matplotlib.pyplot as plt

y_val_UQ = np.zeros([Ny,numLines])
Std_val_UQ = np.zeros([Ny,numLines])
y_val_UQ_Updated = np.zeros([Ny,numLines])
Std_val_UQ_Updated = np.zeros([Ny,numLines])

tmp = np.matmul(H_map_tauW,tauWMean_UQ.T)
tmp_Std = np.matmul(H_map_tauW,tauWStd_UQ.T)
tmp_Updated = np.matmul(H_map_tauW,tauWMean_UQ_Updated.T)
tmp_Std_Updated = np.matmul(H_map_tauW,tauWStd_UQ_Updated.T)

scale = 1
   
plt.plot(tauW_UQ_x,scale*Y_tauW,'k-')

plt.plot(tauW_UQ_x,scale*tmp,'b-')
plt.plot(tauW_UQ_x,scale*tmp+2*tmp_Std,'b--')
plt.plot(tauW_UQ_x,scale*tmp-2*tmp_Std,'b--')

plt.plot(tauW_UQ_x,scale*tmp_Updated,'r-')
plt.plot(tauW_UQ_x,scale*tmp_Updated+2*tmp_Std_Updated,'r--')
plt.plot(tauW_UQ_x,scale*tmp_Updated-2*tmp_Std_Updated,'r--')

plt.plot(tauW_UQ_x,tauW_UQ_x*0,'k--')
plt.ylim(-0.02,0.04)
        
# <codecell> Update U, tauW for the currect Kalman Filter iteration
    
# for i in range(numModes):
#     myUQlib.writeEnKF_U_UQ_Updated_Modes_PerHill(i,U_UQ_Updated[i].reshape(nCells,1), \
#                               nCells, "1/", phyDim) 