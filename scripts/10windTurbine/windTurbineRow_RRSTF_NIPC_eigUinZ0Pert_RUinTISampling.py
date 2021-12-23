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
casePngDir = '/2RRSTF/NIPC/' + caseName+'/png'

if (os.path.exists(DATA+casePngDir))==False:
    print('Making new case data directory...')
    os.makedirs(DATA+casePngDir, exist_ok=True)

# <codecell> Read OpenFOAM case mesh access to grid and cell values
vtkFile   = 'baseCase/VTK/baseCase_10000.vtk'
mesh, nCells, cCenter = myUQlib.getMeshInfo(cwd+vtkFile)
idxCellsLow = np.arange(nCells)

# nCellsLow = 10000
# idxCellsLow = np.random.choice(nCells, size=nCellsLow, replace=False)
# cCenter = cCenter[idxCellsLow]
# nCells = nCellsLow

# <codecell> Importing k, RDet and creating the anisotropic tensor a_ij_Det
RDet     = myUQlib.symmTensorToTensorv2012(
            mesh.cell_arrays['turbulenceProperties:R'][idxCellsLow], nCells)
tkeDet   = mesh.cell_arrays['k'][idxCellsLow]
a_ij_Det = myUQlib.anisotropyTensor(RDet, tkeDet, nCells)

# <codecell> Loading lines data
# lines = (2,4,6)
# numLines = len(lines)
# yPts = 100
# nCells = numLines*yPts
# R_symmTensor = np.loadtxt("baseCase/postProcessing/sample/"+\
#                     "10000/AD1_X_by_D_2_turbulenceProperties:R.xy")[:,1:]
# for l in range(1,numLines):
#     fileName = 'baseCase/postProcessing/sample/10000/'+\
#                 'AD1_X_by_D_'+str(lines[l])+'_turbulenceProperties:R.xy'
#     R_symmTensor = np.vstack((R_symmTensor, np.loadtxt(fileName)[:,1:]))

# RDet    = myUQlib.symmTensorToTensorv1806(R_symmTensor, nCells)
# tkeDet = np.array([sum(RDet[c].diagonal())/2 for c in range(nCells)])
# a_ij_Det = myUQlib.anisotropyTensor(RDet, tkeDet, nCells)

# <codecell> Eigen decomposition of a_ij_Det
eValDet, eVecDet = myUQlib.eigenDecomposition(a_ij_Det, nCells)

# <codecell> Sampling and Perturbing
# Sampling 3 r.v.s s.t. eVal1* + eVal2* + eVal3* = 0
nSamples = 5#30

# delB = cp.Normal(0.5, 0.1)
# delB = cp.Uniform(0.25, 0.75)
# delB = cp.Uniform(0,1) # earlier used
delB = cp.Gamma(3,0.15,0)#cp.Gamma(2,0.15,0)
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

C_DET = myUQlib.baryCentricCoordinates(eValDet[cellNr])
C_eValPertSamples = myUQlib.baryCentricCoordinates(eValPertSamples[sampleNr, cellNr])
C_eValPertSamplesMean = myUQlib.baryCentricCoordinates(eValPertSamplesMean[cellNr])

ax.scatter(C_DET[:,0], C_DET[:,1],  color='r', marker='x', alpha=1)
ax.scatter(C_eValPertSamples[:,0], C_eValPertSamples[:,1], color='k', alpha=0.2)
ax.scatter(C_eValPertSamplesMean[:,0], C_eValPertSamplesMean[:,1], color='k', alpha=1)

fig.savefig(DATA+casePngDir+'/RSamplesBaryCentric_Gamma.png', dpi=300)

# <codecell> Perturbed samples of anisotropy tensor a_ij_PertSamples
# def get_a_ij_R_PertSample(eVec,eValPertSamples,nCells,t,nProcs,s):
#     print(t*nProcs+s)
#     start = timer()
#     a_ij_PertSample = np.array([np.matmul(eVec[c], np.matmul(
#           np.diag(eValPertSamples[s,c]),eVec[c].T)) for c in range(nCells)])
#     R_PertSample = 2*tkeDet.reshape((nCells,1,1)) * (np.eye(3)/3 + a_ij_PertSample)
#     myUQlib.writeSymmTensorField(t*nProcs+s, R_PertSample, nCells, "R",\
#                                       "Samples_R/", "[0 2 -2 0 0 0 0]")
#     myUQlib.writeSymmTensorField(t*nProcs+s, R_PertSample-RDet, nCells, \
#                        "deltaR", "Samples_deltaR/", "[0 2 -2 0 0 0 0]")
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


# <codecell> ##################################################################
# The below code is for constructing a PCE of R and
# plotting its the barycentric coordinates
###############################################################################


# <codecell> Perturbed samples of anisotropy tensor a_ij_PertSamples
def get_a_ij_PertSample(eVec,eValPertSamples,nCells,s):
    print(s)
    a_ij_PertSample = np.array([np.matmul(eVec[c], np.matmul(
        np.diag(eValPertSamples[s,c]),eVec[c].T)) for c in range(nCells)])
    return a_ij_PertSample

pool = Pool(processes=nSamples)
start = timer()
a_ij_PertSamples = np.array(pool.map(partial(get_a_ij_PertSample, \
                            eVecDet,eValPertSamples,nCells), range(nSamples)))
print(timer()-start, 's')
pool.close()
pool.join()

a_ij_PertSamplesMean = a_ij_PertSamples.mean(axis=0)

# <codecell> Perturbed samples of R
R_PertSamples = np.zeros((nSamples, nCells, 3, 3))
for s in range(nSamples):
    R_PertSamples[s] = 2*tkeDet.reshape((nCells,1,1)) * \
                        (np.eye(3)/3 + a_ij_PertSamples[s])

R_PertSamplesMean = R_PertSamples.mean(axis=0)
R_PertSamplesStd  = R_PertSamples.std(axis=0)

# <codecell> PCE of R
d = 5; n = 2; q = 1 # works okayish
# d = 25; n = 1; q = 1 # perfect fit

dist = cp.Normal(0,2)
for i in range(d-1):
    dist = cp.J(dist, cp.Normal(0,2))

# dist       = cp.Iid(cp.Normal(0,2), d)
# dist       = cp.Iid(cp.Uniform(-3,3), d)
phi,norms  = cp.orth_ttr(n,dist,retall=True,normed=True,cross_truncation=1/q)
Pplus1     = len(phi)

np.random.seed(1)
xiSamples  = dist.sample(nSamples, rule="L")

#%%
# R_PCEModes = np.zeros((Pplus1, nCells, 3, 3))

# start = timer()
# for i in range(3):
#     for j in range(3):
#         if i>=j :
#             R_PCEApp = cp.fit_regression(phi, xiSamples, R_PertSamples[:,:,i,j])
#             R_PCEModes[:,:,i,j] = np.array(R_PCEApp.coefficients).T
#             if i>j:
#                 R_PCEModes[:,:,j,i] = R_PCEModes[:,:,i,j]

# print(timer()-start, 's')
# R_PCEStd = np.sum(R_PCEModes[1:]**2, axis=0)**0.5

# tke_R_PCEModes0 = np.array([sum(R_PCEModes[0,c].diagonal())/2 for c in range(nCells)])
# a_ij_R_PCEModes0 = myUQlib.anisotropyTensor(R_PCEModes[0], tke_R_PCEModes0, nCells)
# eVala_ij_R_PCEModes0, eVeca_ij_R_PCEModes0 = myUQlib.eigenDecomposition(a_ij_R_PCEModes0, nCells)

#%%
import pickle

# with open('R_PCEModes.pkl', 'wb') as f:
#     pickle.dump(R_PCEModes, f, protocol=-1)

with open('R_PCEModes.pkl', 'rb') as f:
    R_PCEModes = pickle.load(f)
tke_R_PCEModes0 = np.array([sum(R_PCEModes[0,c].diagonal())/2 for c in range(nCells)])
a_ij_R_PCEModes0 = myUQlib.anisotropyTensor(R_PCEModes[0], tke_R_PCEModes0, nCells)
eVala_ij_R_PCEModes0, eVeca_ij_R_PCEModes0 = myUQlib.eigenDecomposition(a_ij_R_PCEModes0, nCells)

# <codecell> Sampling R from PCE of R (for testing) vis jointDist
UHub  = cp.Uniform(7,9)
TIHub = cp.Uniform(5,7)

jointDist = cp.J(dist, UHub, TIHub)
nSamples_joint = 90
np.random.seed(1)
xiSample  = jointDist.sample(nSamples_joint, rule="L")

xiSample_R_PCE = xiSample[0:d]

#%%
def get_eVala_ij_R_PCESample(phi,xiSample_R_PCE,R_PCEModes,Pplus1,nCells,s):
    print(s)

    R_PCESample = np.sum(np.array([(phi(*(xiSample_R_PCE[:,s]))[m] *\
                            R_PCEModes[m]) for m in range(Pplus1) ]), axis=0)
    tke_R_PCESamples = np.array([sum(R_PCESample[c].diagonal())/2 \
                                 for c in range(nCells)])
    a_ij_R_PCESample = myUQlib.anisotropyTensor(R_PCESample, \
                                                tke_R_PCESamples, nCells)

    myUQlib.writeSymmTensorField(s, R_PCESample-RDet, \
                nCells, 'deltaR', 'Samples_deltaR/', '[0 2 -2 0 0 0 0]')
    myUQlib.writeSymmTensorField(s, a_ij_R_PCESample-a_ij_Det, \
                nCells, 'deltaAij', 'Samples_deltaAij/', '[0 0 0 0 0 0 0]')

    eVala_ij_R_PCESample, _ = myUQlib.eigenDecomposition(a_ij_R_PCESample, nCells)
    return eVala_ij_R_PCESample


# pool = Pool(processes=4)
# start = timer()
# # pool.map(partial(get_eVala_ij_R_PCESample,phi,xiSample_R_PCE,R_PCEModes,\
# #                   Pplus1,nCells), range(nSamples_joint))

# eVala_ij_R_PCESamples = np.array(pool.map(partial(get_eVala_ij_R_PCESample,
#         phi,xiSample_R_PCE,R_PCEModes,Pplus1,nCells), range(nSamples_joint)))

# print(timer()-start, 's')
# pool.close()
# pool.join()

eVala_ij_R_PCESamples = np.array([get_eVala_ij_R_PCESample(
                        phi,xiSample_R_PCE,R_PCEModes,Pplus1,nCells,s) \
                        for s in range(nSamples_joint)])

# <codecell> Computing the Bary Centric coordinates (for testing)
cellNr   = [215]
# cellNr   = np.random.randint(0, nCells, 1)

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

ax.scatter(C_DET[:,0], C_DET[:,1],  color='r', marker='x', alpha=1)
ax.scatter(C_eValPertSamples[:,0], C_eValPertSamples[:,1], color='k', alpha=0.2)
ax.scatter(C_eValPertSamplesMean[:,0], C_eValPertSamplesMean[:,1], color='k', alpha=1)
ax.scatter(C_R_PCESample[:,0],  C_R_PCESample[:,1],  color='b', alpha=0.5)
ax.scatter(C_R_PCEModes0[:,0],  C_R_PCEModes0[:,1],  color='b', alpha=1)

fig.savefig(DATA+casePngDir+'/RSamplesBaryCentric_Normal_d_5_n_2_mu_0_sig_2.png', dpi=300)

# <codecell> Write the Samples
kappa, Cmu = 0.4, 0.09
Uref, Zref  = xiSample[d], 70
num = kappa * (2/3)**0.5
den = xiSample[d+1]/100* Cmu**0.25
Z0 = Zref/np.exp(num/den)

for s in range(nSamples_joint):
    print(s)
    myUQlib.writeScalarValue(s, 'Samples_Uref/', 'Uref', Uref[s])
    myUQlib.writeScalarValue(s, 'Samples_z0/', 'z0', Z0[s])

# <codecell> Setting UQ info. for OpenFOAM simulation
# myUQlib.uqInfo(n, d, q, "eigPertRRSTF", dist)

















































