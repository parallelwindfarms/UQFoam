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
# vtkFile = 'VTK/1nutKLEgPC_lx1.5H_ly0.5H_var0.50_M100x80_iter1_NstdFromEnKF_10000.vtk'
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
    

# <codecell> KLE algorithm
H  = 0.028

lx = 1.50 * H #0.6#1.5#2.0#3.0
ly = 0.50 * H #0.2#0.5#1.0
lz = 0.75 * H #0.75

l  = [lx, ly]
if phyDim==3 :
    l.append(lz)

variance  = 0.5
mean      = -0.5*variance
covModel  = ot.SquaredExponential(l, [np.sqrt(variance)])
eigValK_0 = 1e-2


# <codecell> Set KLE Algorithm Karhunen Loeve P1/SVD Algorithm
# kleAlgo = ot.KarhunenLoeveP1Algorithm(kleMesh, covModel, eigValK_0)

f       = ot.SymbolicFunction(['x','y'],[str(mean)])
fTrend  = ot.TrendTransform(f, kleMesh)
sample  = ot.GaussianProcess(fTrend, covModel, kleMesh).getSample(10000)
# sample  = ot.GaussianProcess(covModel, kleMesh).getSample(10000)
kleAlgo = ot.KarhunenLoeveSVDAlgorithm(sample, eigValK_0)


# <codecell> Run KLE Algorithm    
kleAlgo.run()
kleResult = kleAlgo.getResult()


# <codecell> KLE results
eigVals = np.array(kleResult.getEigenValues())
modes   = kleResult.getModes()
g       = kleResult.getScaledModesAsProcessSample()
g       = np.array(g).reshape(len(g), nCells)


# <codecell> Number of dimension, Poly order, Hyperbolic truncation factor
d = len(g)
n = 3    
q = 0.1


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
        # Setting Observation MAtrix at selected locaitons only
        if (j==2 or j==3 or j==8 or j==9) and (i==9 or i==19 or i==29):
            H[j+i*numLines,i*Nx+xByH_idx[j]] = 1
        if (j>=2 and j<=7) and (abs(yByH_at_xByH[i,j]-yTurn_1)<tol):# or abs(yByH_at_xByH[i,j]-yTurn_2)<tol):
            H[j+i*numLines,i*Nx+xByH_idx[j]] = 1
            
            
# <codecell> Load data
        
# U_DNS        
pathToDNS = os.environ['periodicHillData'] + '/Re2800/0DNS/'

U_DNS = []
for i in range(1,numLines):
    U_DNS.append(np.loadtxt(pathToDNS + 'DNS0'+str(i)+'.dat')[:,:2])
U_DNS.append(np.loadtxt(pathToDNS + 'DNS'+str(numLines)+'.dat')[:,:2])


numModes = 55#55

# nut_UQ
nut_UQ = np.zeros([numModes, nCells])
nut_Std_UQ = np.zeros(nCells)
for i in range(numModes):
    nut_UQ[i] = mesh.cell_arrays['nut'+str(i)]
    if i>0:
        nut_Std_UQ += nut_UQ[i]**2
nut_Mean_UQ = nut_UQ[0]
nut_Std_UQ = nut_Std_UQ**0.5

# U_UQ
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

sampleSize = 1000#5000

dist      = cp.Iid(cp.Normal(0,1/2), d)
xiSamples = dist.sample(size=sampleSize, rule="L")

# polyExpsn = cp.orth_ttr(n, dist, retall=False, normed=True, cross_truncation=1/q)
polyExpsn = cp.orth_ttr(n, dist, retall=False, normed=False, cross_truncation=1/q)


nutSamples = np.zeros([sampleSize, nCells])
for i in range(sampleSize):
    # nutSamples[i]    = np.matmul(nut_UQ.T, np.array(polyExpsn(*(xiSamples.T[i])))[0:numModes])
    nutSamples[i]    = nut_UQ[0]*np.exp(np.dot(g.T,xiSamples[:,i]))
nutSampleMean    = np.mean(nutSamples,0)

USamples = np.zeros([sampleSize, nCells])
for i in range(sampleSize):
    USamples[i]    = np.matmul(U_UQ.T, np.array(polyExpsn(*(xiSamples.T[i])))[0:numModes])
USampleMean    = np.mean(USamples,0)
HSt = np.matmul(H, USamples.T)

print('Relative Error in means from UQ and Sampling = ', np.max(abs(nutSampleMean-nut_Mean_UQ)/nut_Mean_UQ))
print('Max of nutSamples = ', np.max(nutSamples))
print('Min of nutSamples = ', np.min(nutSamples))


# <codecell> Analysis step - Kalman gain matrix (K)
DNSnoise = (0.1*Y_U)
R = np.diag(DNSnoise*DNSnoise)

xi_p  = nutSamples - nutSamples.mean(axis=0,keepdims=1)
HSt_p = HSt - HSt.mean(axis=1,keepdims=1)

PHt  = 1/(1-sampleSize) * np.matmul(HSt_p,xi_p)
HPHt = 1/(1-sampleSize) * np.matmul(HSt_p,HSt_p.T)

K = PHt.T.dot(np.linalg.inv(HPHt+R))

print('Max of K = ', np.max(K))
print('Min of K = ', np.min(K))


# <codecell> Analysis step - updating ensemble

relaxFac = 1
nutSamples_Updated = nutSamples*0.0
for i in range(sampleSize):
    diff = Y_U - HSt.T[i]
    # dnut = (np.matmul(K,diff))
    dnut = (np.matmul(K,diff))
    nutSamples_Updated[i] = nutSamples[i] + dnut*relaxFac

print('Max of nutSamples_Updated = ', np.max(nutSamples_Updated))
print('Min of nutSamples_UPdated = ', np.min(nutSamples_Updated))

# <codecell> Evaluate the expansion coefficients

nut_UQ_Updated = nut_UQ*0.0

for j in range(sampleSize):
    tmp = polyExpsn(*(xiSamples.T[j]))
    for i in range(numModes):
        nut_UQ_Updated[i] += nutSamples_Updated[j] * tmp[i]
nut_UQ_Updated = nut_UQ_Updated/sampleSize 
# nut_UQ_Updated[0] = 1/3 * nut_UQ_Updated[0] ######## !TRICK! ########

nut_Std_UQ_Updated = np.zeros(nCells)
for i in range(1,numModes):
    nut_Std_UQ_Updated += nut_UQ_Updated[i]**2

nut_Mean_UQ_Updated = nut_UQ_Updated[0]
nut_Std_UQ_Updated  = nut_Std_UQ_Updated**0.5


# <codecell> Plot
import matplotlib.pyplot as plt

y_val_UQ = np.zeros([Ny,numLines])
Std_val_UQ = np.zeros([Ny,numLines])
y_val_UQ_Updated = np.zeros([Ny,numLines])
Std_val_UQ_Updated = np.zeros([Ny,numLines])

tmp = np.matmul(H_map,nut_Mean_UQ.T)
tmp_Std = np.matmul(H_map,nut_Std_UQ.T)
tmp_Updated = np.matmul(H_map,nut_Mean_UQ_Updated.T)
tmp_Std_Updated = np.matmul(H_map,nut_Std_UQ_Updated.T)

scale = 100

for i in range(2,numLines):
   
    # plt.plot(scale*Y_U_lineVal[:,i]+1*i-1,yByH_at_xByH[:,i],'k-')
    for j in range(Ny):
        y_val_UQ[j,i] = tmp[i+j*numLines] 
        Std_val_UQ[j,i] = tmp_Std[i+j*numLines] 
        y_val_UQ_Updated[j,i] = tmp_Updated[i+j*numLines] 
        Std_val_UQ_Updated[j,i] = tmp_Std_Updated[i+j*numLines] 
    
    plt.plot(scale*y_val_UQ[:,i]+1*i-1,yByH_at_xByH[:,i],'b')
    # plt.plot(scale*(y_val_UQ[:,i]+2*Std_val_UQ[:,i])+1*i-1,yByH_at_xByH[:,i],'b--')
    # plt.plot(scale*(y_val_UQ[:,i]-2*Std_val_UQ[:,i])+1*i-1,yByH_at_xByH[:,i],'b--')

    plt.plot(scale*y_val_UQ_Updated[:,i]+1*i-1,yByH_at_xByH[:,i],'r')
    # plt.plot(scale*(y_val_UQ_Updated[:,i]+2*Std_val_UQ_Updated[:,i])+1*i-1,yByH_at_xByH[:,i],'r--')
    # plt.plot(scale*(y_val_UQ_Updated[:,i]-2*Std_val_UQ_Updated[:,i])+1*i-1,yByH_at_xByH[:,i],'r--')
        
# <codecell> Update U for the currect Kalman Filter iteration
    
timeDir = "1/"
# timeDir = "2/"
    
for i in range(numModes):
    myUQlib.writeEnKF_nut_UQ_Updated_Modes_PerHill(i,nut_UQ_Updated[i].reshape(nCells,1), \
                              nCells, timeDir , phyDim) 