#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 29 15:21:21 2021

@author: jigar
"""

# %% Import Required Modules and env vars
import os, sys

import numpy as np
import chaospy as cp
import matplotlib.pyplot as plt
import pickle

from multiprocessing import Pool
from functools import partial
from timeit import default_timer as timer

SCRIPTS = os.environ['SCRIPTS']
DATA    = os.environ['windTurbineData']

sys.path.append(SCRIPTS+'/10windTurbine')
sys.path.append(SCRIPTS+'/pythonScripts')

cwd = os.getcwd() + "/"

import myUQlib

# %% Plot settings
myUQlib.rcParamsSettings(15)

# %% Case details
caseName   = cwd.split('/')[-2]
casePngDir = '/2RRSTF/NIPC/uqPaperCases/' + caseName+'/png'

if (os.path.exists(DATA+casePngDir))==False:
    print('Making new case data directory...')
    os.makedirs(DATA+casePngDir, exist_ok=True)

# %% Read OpenFOAM case mesh access to grid and cell values
vtkFile   = 'baseCase/VTK/baseCase_10000.vtk'
mesh, nCells, cCenter = myUQlib.getMeshInfo(cwd+vtkFile)
idxCellsLow = np.arange(nCells)

# %% Importing k, RDet and creating the anisotropic tensor a_ij_Det
RDet     = myUQlib.symmTensorToTensorv2012(
            mesh.cell_data['turbulenceProperties:R'][idxCellsLow], nCells)
tkeDet   = mesh.cell_data['k'][idxCellsLow]
a_ij_Det = myUQlib.anisotropyTensor(RDet, tkeDet, nCells)

# %% Eigen decomposition of a_ij_Det
eValDet, eVecDet = myUQlib.eigenDecomposition(a_ij_Det, nCells)

# %% Sampling and Perturbing
# Sampling 3 r.v.s s.t. eVal1* + eVal2* + eVal3* = 0
# For 1 WT
# nSamples = 30
# delB = cp.Gamma(2,0.15,0)

# For >1 WTs
nSamples = 5
delB = cp.Gamma(3,0.15,0) 

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

# %% Computing the Bary Centric coordinates (for testing)
sampleNr = np.arange(nSamples)
cellNr   = [215]#np.random.randint(0, nCells, 1)

fig = plt.figure()
ax = fig.add_subplot(111)

myUQlib.plotBaryCentricCoordinateSystem(ax)

C_DET = myUQlib.baryCentricCoordinates(eValDet[cellNr])
C_eValPertSamples = myUQlib.baryCentricCoordinates(
    eValPertSamples[sampleNr, cellNr])
C_eValPertSamplesMean = myUQlib.baryCentricCoordinates(
    eValPertSamplesMean[cellNr])

ax.scatter(C_DET[:,0], C_DET[:,1],  
           color='r', marker='x', alpha=1)
ax.scatter(C_eValPertSamples[:,0], C_eValPertSamples[:,1], 
           color='k', alpha=0.2)
ax.scatter(C_eValPertSamplesMean[:,0], C_eValPertSamplesMean[:,1], 
           color='k', alpha=1)

fig.savefig(DATA+casePngDir+'/RSamplesBaryCentric_Gamma.png', dpi=300)

# %% Perturbed gamma samples of anisotropy tensor a_ij_PertSamples
# def get_a_ij_R_PertSample(eVec,eValPertSamples,nCells,t,nProcs,s):
#     print(t*nProcs+s)
#     start = timer()
#     a_ij_PertSample = np.array([np.matmul(eVec[c], np.matmul(
#           np.diag(eValPertSamples[s,c]),eVec[c].T)) for c in range(nCells)])
#     R_PertSample = 2*tkeDet.reshape((nCells,1,1)) * (np.eye(3)/3 + a_ij_PertSample)
#     myUQlib.writeSymmTensorField(t*nProcs+s, R_PertSample, nCells, "R",\
#                                       "Samples_R/", "[0 2 -2 0 0 0 0]")
#     myUQlib.writeSymmTensorField(t*nProcs+s, R_PertSample-RDet, nCells, \
#                         "deltaR", "Samples_deltaR/", "[0 2 -2 0 0 0 0]")
#     print(timer()-start, 's')

# nProcs = 10
# start = timer()
# for t in range(9):
#     pool = Pool(processes=nProcs)
#     np.array(pool.map(partial(get_a_ij_R_PertSample, eVecDet,\
#       eValPertSamples[t*nProcs:(t+1)*nProcs],nCells,t,nProcs), range(nProcs)))
#     pool.close()
#     pool.join()
# print(timer()-start, 's')

# %% ##################################################################
# The below code is for -
# 1. Getting perturbed a_ij, R
# 2. Constructing a PCE of R
# 3. Plotting the barycentric coordinates
#######################################################################

# %% Perturbed samples of anisotropy tensor a_ij_PertSamples
def get_a_ij_PertSample(eVec,eValPertSamples,nCells,s):
    print(s)
    a_ij_PertSample = np.array([np.matmul(eVec[c], np.matmul(
        np.diag(eValPertSamples[s,c]),eVec[c].T)) for c in range(nCells)])
    return a_ij_PertSample

pool = Pool(processes=32)
start = timer()
a_ij_PertSamples = np.array(pool.map(partial(get_a_ij_PertSample, \
                            eVecDet,eValPertSamples,nCells), range(nSamples)))
print(timer()-start, 's')
pool.close()
pool.join()

a_ij_PertSamplesMean = a_ij_PertSamples.mean(axis=0)

# %% Perturbed samples of R
R_PertSamples = np.zeros((nSamples, nCells, 3, 3))
for s in range(nSamples):
    R_PertSamples[s] = 2*tkeDet.reshape((nCells,1,1)) * \
                        (np.eye(3)/3 + a_ij_PertSamples[s])

R_PertSamplesMean = R_PertSamples.mean(axis=0)
R_PertSamplesStd  = R_PertSamples.std(axis=0)

# %% PCE of R
# For 1 WT
# d = 20; n = 1; q = 1
# dist = cp.Normal(0,10)
# for i in range(d-1):
#     dist = cp.J(dist, cp.Normal(0,10))
    
# For >1 WTs
d = 5; n = 2; q = 1 
dist = cp.Normal(0,1.25)
for i in range(d-1):
    dist = cp.J(dist, cp.Normal(0,1.25))
    
phi,norms = cp.orth_ttr(n,dist,retall=True,normed=True,cross_truncation=1/q)
Pplus1 = len(phi)

np.random.seed(1)
xiSamples = dist.sample(nSamples, rule="L")

# %% Computing PCE of R
# R_PCEModes = np.zeros((Pplus1, nCells, 3, 3))

# start = timer()
# for i in range(3):
#     for j in range(3):
#         if i>=j :
            # R_PCEApp = cp.fit_regression(
            #     phi, xiSamples, R_PertSamples[:,:,i,j]
            # )
#             R_PCEModes[:,:,i,j] = np.array(R_PCEApp.coefficients).T
#             if i>j:
#                 R_PCEModes[:,:,j,i] = R_PCEModes[:,:,i,j]

# print(timer()-start, 's')
# R_PCEStd = np.sum(R_PCEModes[1:]**2, axis=0)**0.5

# tke_R_PCEModes0 = np.array(
#     [sum(R_PCEModes[0,c].diagonal())/2 for c in range(nCells)])
# a_ij_R_PCEModes0 = myUQlib.anisotropyTensor(
#     R_PCEModes[0], tke_R_PCEModes0, nCells)
# eVala_ij_R_PCEModes0, eVeca_ij_R_PCEModes0 = myUQlib.eigenDecomposition(
#     a_ij_R_PCEModes0, nCells)

# Dumping PCE of R
# with open('R_PCEModes.pkl', 'wb') as f:
#     pickle.dump(R_PCEModes, f, protocol=-1)

# Importing PCE of R
with open('R_PCEModes.pkl', 'rb') as f:
    R_PCEModes = pickle.load(f)
    
    
tke_R_PCEModes0 = np.array(
    [sum(R_PCEModes[0,c].diagonal())/2 for c in range(nCells)])
a_ij_R_PCEModes0 = myUQlib.anisotropyTensor(
    R_PCEModes[0], tke_R_PCEModes0, nCells)
eVala_ij_R_PCEModes0, eVeca_ij_R_PCEModes0 = myUQlib.eigenDecomposition(
    a_ij_R_PCEModes0, nCells)

# %% Sampling R from PCE of R (for testing) via jointDist
UHub  = cp.Uniform(7,9)
TIHub = cp.Uniform(5,7)

jointDist = cp.J(dist, UHub, TIHub)
nSamples_joint = 100
np.random.seed(1)
xiSample = jointDist.sample(nSamples_joint, rule="L")

xiSample_R_PCE = xiSample[0:d]

#%% Get eVals of a_ij_R_PCESample
def get_eVala_ij_R_PCESample(phi,xiSample_R_PCE,R_PCEModes,Pplus1,nCells,s):
    print(s)

    R_PCESample = np.sum(np.array([(phi(*(xiSample_R_PCE[:,s]))[m] *\
                            R_PCEModes[m]) for m in range(Pplus1) ]), axis=0)
    tke_R_PCESamples = np.array([sum(R_PCESample[c].diagonal())/2 \
                                 for c in range(nCells)])
    a_ij_R_PCESample = myUQlib.anisotropyTensor(R_PCESample, \
                                                tke_R_PCESamples, nCells)
    # Uncomment for writing
    # myUQlib.writeSymmTensorField(s, R_PCESample-RDet, \
    #             nCells, 'deltaR', 'Samples_deltaR/', '[0 2 -2 0 0 0 0]')
    # myUQlib.writeSymmTensorField(s, a_ij_R_PCESample-a_ij_Det, \
    #             nCells, 'deltaAij', 'Samples_deltaAij/', '[0 0 0 0 0 0 0]')

    eVala_ij_R_PCESample, _ = myUQlib.eigenDecomposition(a_ij_R_PCESample, nCells)
    return eVala_ij_R_PCESample


pool = Pool(processes=32)
start = timer()
eVala_ij_R_PCESamples = np.array(pool.map(partial(get_eVala_ij_R_PCESample,
        phi,xiSample_R_PCE,R_PCEModes,Pplus1,nCells), range(nSamples_joint)))

print(timer()-start, 's')
pool.close()
pool.join()

# %% Plotting the Bary Centric coordinates (for testing)
cellNr   = [215]#np.random.randint(0, nCells, 1)

sampleNr = np.arange(nSamples)
sampleNr_PCESamples = np.arange(nSamples_joint)

fig = plt.figure()
ax = fig.add_subplot(111)

myUQlib.plotBaryCentricCoordinateSystem(ax)

C_DET                 = myUQlib.baryCentricCoordinates(eValDet[cellNr])
C_eValPertSamples     = myUQlib.baryCentricCoordinates(eValPertSamples[sampleNr, cellNr])
C_eValPertSamplesMean = myUQlib.baryCentricCoordinates(eValPertSamplesMean[cellNr])
C_R_PCESample         = myUQlib.baryCentricCoordinates(eVala_ij_R_PCESamples[sampleNr_PCESamples, cellNr])
C_R_PCEModes0         = myUQlib.baryCentricCoordinates(eVala_ij_R_PCEModes0[cellNr])

ax.scatter(C_DET[:,0], C_DET[:,1], color='k', marker='o', alpha=1)
# ax.scatter(C_eValPertSamples[:,0], C_eValPertSamples[:,1], color='k', alpha=0.2)
# ax.scatter(C_eValPertSamplesMean[:,0], C_eValPertSamplesMean[:,1], color='k', alpha=1)
ax.scatter(C_R_PCESample[:,0], C_R_PCESample[:,1], color='b', alpha=0.1)
ax.scatter(C_R_PCEModes0[:,0], C_R_PCEModes0[:,1], color='b', alpha=1)

fig.savefig(DATA+casePngDir+'/RSamplesBaryCentric.png', dpi=150)

# %% Write the Samples
kappa, Cmu = 0.4, 0.09
Uref, Zref  = xiSample[d], 70
num = kappa * (2/3)**0.5
den = xiSample[d+1]/100* Cmu**0.25
Z0 = Zref/np.exp(num/den)

for s in range(nSamples_joint):
    print(s)
    # Uncomment for writing
    # myUQlib.writeScalarValue(s, 'Samples_Uref/', 'Uref', Uref[s])
    # myUQlib.writeScalarValue(s, 'Samples_z0/', 'z0', Z0[s])
    
# %% Setting UQ info. for OpenFOAM simulation
# myUQlib.uqInfo(n, d, q, "eigPertRRSTF", dist)  
