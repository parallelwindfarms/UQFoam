#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 29 15:21:21 2021

@author: jigar
"""

# <codecell> Import Required Modules and env vars
import os, sys

import numpy as np
import chaospy as cp
import matplotlib.pyplot as plt

from multiprocessing import Pool
from functools import partial
from timeit import default_timer as timer

SCRIPTS = os.environ['SCRIPTS']
DATA    = os.environ['windTurbineData']

sys.path.append(SCRIPTS+'/10windTurbine')
sys.path.append(SCRIPTS+'/pythonScripts')

cwd = os.getcwd() + "/"

import myUQlib

# <codecell> Plot settings
myUQlib.rcParamsSettings(15)

# <codecell> Case details
caseName   = cwd.split('/')[-2]
casePngDir = '/2RRSTF/NIPC/uqPaperCases/' + caseName+'/png'

if (os.path.exists(DATA+casePngDir))==False:
    print('Making new case data directory...')
    os.makedirs(DATA+casePngDir, exist_ok=True)
    
# <codecell> Read OpenFOAM case mesh access to grid and cell values
vtkFile   = 'baseCase/VTK/baseCase_2500.vtk'
mesh, nCells, cCenter = myUQlib.getMeshInfo(cwd+vtkFile)

# <codecell> Importing k, RDet and creating the anisotropic tensor a_ij_Det
RDet     = myUQlib.symmTensorToTensorv2012(mesh.cell_arrays['turbulenceProperties:R'], nCells)
tkeDet   = mesh.cell_arrays['k']
a_ij_Det = myUQlib.anisotropyTensor(RDet, tkeDet, nCells)

# <codecell> Loading lines data
# lines = (2,4,7,12)
# numLines = len(lines)
# yPts = 100
# nCells = numLines*yPts
# R_symmTensor = np.loadtxt("baseCase/postProcessing/sample/1000/X_by_D_2_turbulenceProperties:R.xy")[:,1:]
# for l in range(1,numLines):
#     fileName = 'baseCase/postProcessing/sample/1000/X_by_D_'+str(lines[l])+'_turbulenceProperties:R.xy'
#     R_symmTensor = np.vstack((R_symmTensor, np.loadtxt(fileName)[:,1:]))

# RDet    = myUQlib.symmTensorToTensorv1806(R_symmTensor, nCells)
# tkeDet = np.array([sum(RDet[c].diagonal())/2 for c in range(nCells)])
# a_ij_Det = myUQlib.anisotropyTensor(RDet, tkeDet, nCells)
    
# <codecell> Eigen decomposition of a_ij_Det
eValDet, eVecDet = myUQlib.eigenDecomposition(a_ij_Det, nCells)

# <codecell> Sampling and Perturbing
# Sampling 3 r.v.s s.t. eVal1* + eVal2* + eVal3* = 0
nSamples = 30

# delB = cp.Normal(0.5, 0.1)
# delB = cp.Uniform(0.25, 0.75)
# delB = cp.Uniform(0,1) # earlier used
delB = cp.Gamma(2,0.15,0)
# delB = cp.Iid(delB, 3)

np.random.seed(1)
delBSamples = delB.sample(nSamples, rule='L')
delBSamples.sort()
print(f"delBSamples min = {delBSamples.min()} , max = {delBSamples.max()}")

eVal_c = np.array([[2/3,-1/3,-1/3], [1/6,1/6,-1/3], [0,0,0]])

# Perturbed anisotropy tensor eigenvalues
eValPertSamples = np.array(
    [eValDet + delBSamples[s]*(eVal_c[0,:] - eValDet) 
                                for s in range(0, nSamples)])
eValPertSamples = np.vstack((eValPertSamples, np.array(
    [eValDet + delBSamples[s]*(eVal_c[1,:] - eValDet)
                                for s in range(0, nSamples)]) ))
eValPertSamples = np.vstack((eValPertSamples, np.array(
    [eValDet + delBSamples[s]*(eVal_c[2,:] - eValDet)
                                for s in range(0, nSamples)]) ))

eValPertSamplesMean = eValPertSamples.mean(axis=0)
                                           
nSamples = nSamples*3
                             
# <codecell> Computing the Bary Centric coordinates (for testing)
sampleNr = np.arange(nSamples)
cellNr   = [215]#np.random.randint(0, nCells, 1)

fig = plt.figure()
ax = fig.add_subplot(111)

myUQlib.plotBaryCentricCoordinateSystem(ax)

C_DET                 = myUQlib.baryCentricCoordinates(eValDet[cellNr])
C_eValPertSamples     = myUQlib.baryCentricCoordinates(eValPertSamples[sampleNr, cellNr])
C_eValPertSamplesMean = myUQlib.baryCentricCoordinates(eValPertSamplesMean[cellNr])

ax.scatter(C_DET[:,0], C_DET[:,1],  color='r', marker='x', alpha=1)
ax.scatter(C_eValPertSamples[:,0], C_eValPertSamples[:,1], color='k', alpha=0.2)
ax.scatter(C_eValPertSamplesMean[:,0], C_eValPertSamplesMean[:,1], color='k', alpha=1)

fig.savefig(DATA+casePngDir+'/RSamplesBaryCentric_Gamma.png', dpi=300)

# <codecell> Perturbed samples of anisotropy tensor a_ij_PertSamples
def get_a_ij_R_PertSample(eVec,eValPertSamples,nCells,t,nProcs,s):
    print(t*nProcs+s)
    start = timer()
    a_ij_PertSample = np.array([np.matmul(eVec[c], np.matmul(np.diag(eValPertSamples[s,c]),eVec[c].T)) 
                                        for c in range(nCells)])    
    R_PertSample = 2*tkeDet.reshape((nCells,1,1)) * (np.eye(3)/3 + a_ij_PertSample)
    # myUQlib.writeSymmTensorField(t*nProcs+s, R_PertSample, nCells, "R", "R_Samples/", "[0 2 -2 0 0 0 0]")
    myUQlib.writeSymmTensorField(t*nProcs+s, R_PertSample-RDet, nCells, "deltaR", "deltaR_Samples/", "[0 2 -2 0 0 0 0]")
    print(timer()-start, 's')

nProcs = 10
start = timer()
for t in range(9):
    pool = Pool(processes=nProcs)
    np.array(pool.map(partial(get_a_ij_R_PertSample, eVecDet,eValPertSamples[t*nProcs:(t+1)*nProcs],nCells,t,nProcs), range(nProcs)))
    pool.close()
    pool.join()
print(timer()-start, 's')

# <codecell> ###################################################################################
# The below code is for constructing a PCE of R and plotting its the barycentric coordinates
################################################################################################


# <codecell> Perturbed samples of anisotropy tensor a_ij_PertSamples
def get_a_ij_PertSample(eVec,eValPertSamples,nCells,s):
    print(s)
    a_ij_PertSample = np.array([np.matmul(eVec[c], np.matmul(np.diag(eValPertSamples[s,c]),eVec[c].T)) 
                                        for c in range(nCells)])
    return a_ij_PertSample

pool = Pool(processes=8)
start = timer()
a_ij_PertSamples = np.array(pool.map(partial(get_a_ij_PertSample, eVecDet,eValPertSamples,nCells), range(nSamples)))
print(timer()-start, 's')
pool.close()
pool.join()

a_ij_PertSamplesMean = a_ij_PertSamples.mean(axis=0)

# <codecell> Perturbed samples of R
R_PertSamples = np.zeros((nSamples, nCells, 3, 3))
for s in range(nSamples):
    R_PertSamples[s] = 2*tkeDet.reshape((nCells,1,1)) * (np.eye(3)/3 + a_ij_PertSamples[s])

R_PertSamplesMean = R_PertSamples.mean(axis=0)
R_PertSamplesStd  = R_PertSamples.std(axis=0)

# <codecell> PCE of R
# d = 5; n = 2; q = 1 # works okayish
d = 25; n = 1; q = 1 # perfect fit

dist       = cp.Iid(cp.Normal(0,1.25), d)
# dist       = cp.Iid(cp.Uniform(-3,3), d)
phi,norms  = cp.orth_ttr(n,dist,retall=True,normed=True,cross_truncation=1/q)    
Pplus1     = len(phi)

np.random.seed(1)
xiSamples  = dist.sample(nSamples, rule="L")

R_PCEModes = np.zeros((Pplus1, nCells, 3, 3))
 

start = timer()
for i in range(3):
    for j in range(3):
        if i>=j :
            R_PCEApp = cp.fit_regression(phi, xiSamples, R_PertSamples[:,:,i,j])
            R_PCEModes[:,:,i,j] = np.array(R_PCEApp.coefficients).T
            if i>j:
                R_PCEModes[:,:,j,i] = R_PCEModes[:,:,i,j]

print(timer()-start, 's')
R_PCEStd = np.sum(R_PCEModes[1:]**2, axis=0)**0.5

tke_R_PCEModes0 = np.array([sum(R_PCEModes[0,c].diagonal())/2 for c in range(nCells)])
a_ij_R_PCEModes0 = myUQlib.anisotropyTensor(R_PCEModes[0], tke_R_PCEModes0, nCells)
eVala_ij_R_PCEModes0, eVeca_ij_R_PCEModes0 = myUQlib.eigenDecomposition(a_ij_R_PCEModes0, nCells)

# <codecell> Sampling R from PCE of R (for testing)
nSamples_R_PCE = 100#400
np.random.seed(1)
xiSample_R_PCE = dist.sample(nSamples_R_PCE, rule="L")
eVala_ij_R_PCESamples = np.zeros((nSamples_R_PCE, nCells, 3))

def get_eVala_ij_R_PCESample(phi,xiSample_R_PCE,R_PCEModes,Pplus1,nCells,s):
    R_PCESamples = np.sum(np.array([(phi(*(xiSample_R_PCE[:,s]))[m] * R_PCEModes[m]) for m in range(Pplus1) ]), axis=0)
    # R_PCESamples = np.sum(np.array([(phi(xiSample_R_PCE[s])[m] * R_PCEModes[m]) for m in range(Pplus1) ]), axis=0)
    tke_R_PCESamples = np.array([sum(R_PCESamples[c].diagonal())/2 for c in range(nCells)])
    a_ij_R_PCESample = myUQlib.anisotropyTensor(R_PCESamples, tke_R_PCESamples, nCells)
    eVala_ij_R_PCESample, _ = myUQlib.eigenDecomposition(a_ij_R_PCESample, nCells)
    return eVala_ij_R_PCESample

pool = Pool(processes=4)
start = timer()
eVala_ij_R_PCESamples = np.array(pool.map(partial(get_eVala_ij_R_PCESample, 
                            phi,xiSample_R_PCE,R_PCEModes,Pplus1,nCells), range(nSamples_R_PCE)))
print(timer()-start, 's')
pool.close()
pool.join()

# <codecell> Computing the Bary Centric coordinates (for testing)
cellNr   = [215]
# cellNr   = np.random.randint(0, nCells, 1)

sampleNr = np.arange(nSamples)
sampleNr_PCESamples = np.arange(nSamples_R_PCE)

fig = plt.figure()
ax = fig.add_subplot(111)

myUQlib.plotBaryCentricCoordinateSystem(ax)

C_DET                 = myUQlib.baryCentricCoordinates(eValDet[cellNr])
C_eValPertSamples     = myUQlib.baryCentricCoordinates(eValPertSamples[sampleNr, cellNr])
C_eValPertSamplesMean = myUQlib.baryCentricCoordinates(eValPertSamplesMean[cellNr])
C_R_PCESample         = myUQlib.baryCentricCoordinates(eVala_ij_R_PCESamples[sampleNr_PCESamples, cellNr])
C_R_PCEModes0         = myUQlib.baryCentricCoordinates(eVala_ij_R_PCEModes0[cellNr])

ax.scatter(C_DET[:,0], C_DET[:,1],  color='r', marker='x', alpha=1)
ax.scatter(C_eValPertSamples[:,0], C_eValPertSamples[:,1], color='k', alpha=0.2)
ax.scatter(C_eValPertSamplesMean[:,0], C_eValPertSamplesMean[:,1], color='k', alpha=1)
ax.scatter(C_R_PCESample[:,0],  C_R_PCESample[:,1],  color='b', alpha=0.1)
ax.scatter(C_R_PCEModes0[:,0],  C_R_PCEModes0[:,1],  color='b', alpha=1)

# fig.savefig(DATA+casePngDir+'/RSamplesBaryCentric_Normal_d_25_n_1.png', dpi=300)
# fig.savefig(DATA+casePngDir+'/RSamplesBaryCentric_Uniform_d_25_n_1.png', dpi=300)

# <codecell> Write the R_PCEModes
# for s in range(Pplus1):
#     print(s)
#     myUQlib.writeSymmTensorField(s, R_PCEModes[s], nCells, "R", "R_PCEModes/", "[0 2 -2 0 0 0 0]")

# <codecell> Setting UQ info. for OpenFOAM simulation
# myUQlib.uqInfo(n, d, q, "eigPertRRSTF", dist)  
