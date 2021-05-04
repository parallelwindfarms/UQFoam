#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 20 18:05:20 2020

@author: jigar
"""

from __future__ import print_function

import numpy as np

import chaospy as cp
import pyvista as pv

import itertools as itr
import time

from joblib import Parallel, delayed

# Local modules
from getFields   import *
from writeFields import *
from plotSettings import *

# <codecell> Read OpenFOAM case mesh acces to grid and cell values
def getMeshInfo(vtkFile):
  
    mesh    = pv.UnstructuredGrid(vtkFile)
    nCells  = mesh.n_cells
    cCenter = mesh.cell_centers().points
    # sized   = mesh.compute_cell_sizes()
    # cVols   = sized.cell_arrays["Volume"]
    # mVol    = mesh.volume
    # volAvg  = cVols/mVol
    
    return mesh, nCells, cCenter

# <codecell> Function for setting UQ info. for OpenFOAM simulation
def triplePdt1(i, phi3, dist):
    start_time = time.time()
    tmpCijk = cp.E(phi3, dist)
    print(i, "\t t = ", time.time() - start_time)
    return tmpCijk

def triplePdt2(i, polyExpsn, dist, N):
    start_time = time.time()
    tmpCijk = np.zeros((N,N))
    for j in range(N):
        tmpCijk[j] = cp.E(polyExpsn[i]*polyExpsn[j]*polyExpsn, dist)
    print(i, "\t t = ", time.time() - start_time)
    return tmpCijk

def triplePdt3(i, polyExpsn, dist, N):
    start_time = time.time()
    tmpCijk = np.zeros((N,N))
    for j in range(N):
        if i>=j:
            for k in range(N):
                if j>=k:
                    tmpCijk[j,k] = cp.E(polyExpsn[i]*polyExpsn[j]*polyExpsn[k], dist)
    print(i, "\t t = ", time.time() - start_time)
    return tmpCijk


def uqInfo(n, d, q, RF, dist):
    if dist=='IidNormal':
        dist = cp.Iid(cp.Normal(0,1), d)
    
    phi,norms = cp.orth_ttr(n, dist, retall=True, normed=True, cross_truncation=1/q)
    polyExpsn = cp.orth_ttr(n, dist, retall=False, normed=True, cross_truncation=1/q)
    # phi,norms = cp.orth_ttr(n, dist, retall=True, normed=False, cross_truncation=1/q)
    # polyExpsn = cp.orth_ttr(n, dist, retall=False, normed=False, cross_truncation=1/q)
    
    Pplus1    = len(phi)
    print("Pplus1 = ", Pplus1)
    N = Pplus1
    SMALL = 1e-5#1e-12
    
    # Nomalized outer and inner products of the polynomials
    Ckk  = cp.E(phi*phi, dist)

    if N<=51:
        
        phi3 = cp.outer(phi, phi, phi)
        Cijk = cp.E(phi3, dist)
          
    else: 
        
        # Looping over i,j,k
        start_time_loop = time.time()
        Cijk = np.array(Parallel(n_jobs=4, backend='multiprocessing')\
            (delayed(triplePdt3)(i, polyExpsn, dist, N) for i in range(N)))
        print("Loop time = ", time.time() - start_time_loop)
        
        for i in range(N):
            for j in range(N):
                if i>=j:
                    for k in range(N):
                        if j>=k:
                            if abs(Cijk[i,j,k])>=SMALL:
                                tmpCijk = Cijk[i,j,k]
                                if j!=k:
                                    Cijk[i,k,j] = tmpCijk
                                Cijk[j,i,k] = tmpCijk
                                if k!=i:
                                    Cijk[j,k,i] = tmpCijk                        
                                Cijk[k,i,j] = tmpCijk
                                if i!=j:
                                    Cijk[k,j,i] = tmpCijk
                                           
        print("Loop time = ", time.time() - start_time_loop)
    
    
        # # Double outer pdt - parallel
        # start_time_loop = time.time()
        # phi2 = cp.outer(polyExpsn, polyExpsn)
        # Cijk_par = np.array(Parallel(n_jobs=4, backend='multiprocessing')\
        #     (delayed(triplePdt1)(i, polyExpsn[i]*phi2, dist) for i in range(N)))
        # print("Loop time = ", time.time() - start_time_loop)


        # # Triple outer pdt - parallel
        # start_time_loop = time.time()
        # Cijk_par = np.array(Parallel(n_jobs=4, backend='multiprocessing')\
        #     (delayed(triplePdt2)(i, polyExpsn, dist, N) for i in range(N)))
        # print("Loop time = ", time.time() - start_time_loop)

    # Write out Cikj, Cikjl in OpenFOAM format
    print ("Writing Mijk, Mijkl into constant/gPCcoeffs ...\n")
        
    UQProperties = open("constant/uqInfo","w+")

    UQProperties.write("\
/*--------------------------------*- C++ -*----------------------------------*\\\n\
| =========                 |                                                 |\n\
| \\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n\
|  \\\    /   O peration     | Version:  v1806                                 |\n\
|   \\\  /    A nd           | Web:      www.OpenFOAM.com                      |\n\
|    \\\/     M anipulation  |                                                 |\n\
\*---------------------------------------------------------------------------*/\n\
FoamFile\n\
{\n\
    version     2.0;\n\
    format      ascii;\n\
    class       dictionary;\n\
    location    \"constant\";\n\
    object      uqInfo;\n\
}\n\
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n")

    UQProperties.write("// Type of random field\n\n")
    UQProperties.write("randField\t"+RF+";\n\n")
    
    UQProperties.write("// Number of dimensions\n\n")
    UQProperties.write("d\t"+str(d)+";\n\n")

    UQProperties.write("// Number of modes\n\n")
    UQProperties.write("P\t"+str(Pplus1-1)+";\n\n")
       
    UQProperties.write("\n\n// Quad product of basis function k with itself \n\n")
    for i in range(Pplus1):
        UQProperties.write("M["+str(i)+"]\t" + str(Ckk[i])+";\n")
    
    UQProperties.write("\n\n// Triple product of basis functions\n\n")
    for i,j,k in itr.product(*[range(Pplus1)]*3):
        if (np.abs(Cijk[i,j,k])>=SMALL):
            UQProperties.write("M["+str(i)+"]["+str(j)+"]["+str(k)+"]\t" \
                                  +str(Cijk[i,j,k])+";\n")       
        
    UQProperties.close()
    
    print ("Done.\n")

# <codecell> Convert symmTensor to Tensor
def symmTensorToTensorv2012(symmTensor, nCells):
    tensor = np.zeros((nCells, 3, 3))

    tensor[:,0,0],tensor[:,0,1],tensor[:,0,2],\
    tensor[:,1,0],tensor[:,1,1],tensor[:,1,2],\
    tensor[:,2,0],tensor[:,2,1],tensor[:,2,2] = \
    symmTensor[:,0],symmTensor[:,3],symmTensor[:,5],\
    symmTensor[:,3],symmTensor[:,1],symmTensor[:,4],\
    symmTensor[:,5],symmTensor[:,4],symmTensor[:,2]    
    
    return tensor

def symmTensorToTensorv1806(symmTensor, nCells):
    tensor = np.zeros((nCells, 3, 3))
    
    tensor[:,0,0],tensor[:,0,1],tensor[:,0,2],\
    tensor[:,1,0],tensor[:,1,1],tensor[:,1,2],\
    tensor[:,2,0],tensor[:,2,1],tensor[:,2,2] = \
    symmTensor[:,0],symmTensor[:,1],symmTensor[:,2],\
    symmTensor[:,1],symmTensor[:,3],symmTensor[:,4],\
    symmTensor[:,2],symmTensor[:,4],symmTensor[:,5] 
    
    return tensor

# <codecell> Bary Centric Coordinates
def baryCentricCoordinates(eigValsArray):
    
    if eigValsArray.ndim==2 :
        numElem = eigValsArray.shape[0]
        
        C_vec = np.zeros((numElem, 2))
    
        for n in range(numElem):
            eigVals = eigValsArray[n]
            
            # eigVals.sort()
            # eigVals = eigVals[::-1]
            # print('eigVals =', eigVals)
            C = [ eigVals[0] - eigVals[1] ]
            C.append(2*(eigVals[1] - eigVals[2]))
            C.append(3*eigVals[2] + 1)
            # print('C = ', C)
            e1 = np.array([1,0]); e2 = np.array([-1,0]); e3 = np.array([0,np.sqrt(3)])
            
            C_vec[n] = C[0]*e1 + C[1]*e2 + C[2]*e3
        
        return C_vec
    
    else:
        eigVals = eigValsArray
        
        # eigVals.sort()
        # eigVals = eigVals[::-1]
        # print('eigVals =', eigVals)
        C = [ eigVals[0] - eigVals[1] ]
        C.append(2*(eigVals[1] - eigVals[2]))
        C.append(3*eigVals[2] + 1)
        # print('C = ', C)
        e1 = np.array([1,0]); e2 = np.array([-1,0]); e3 = np.array([0,np.sqrt(3)])
        
        C_vec = C[0]*e1 + C[1]*e2 + C[2]*e3
        
        return C_vec    

# <codecell> Return Anisotropy Tensor A
def anisotropyTensor(R, tke, nCells):
    return np.array(R/(2*tke.reshape(nCells,1,1)) - np.eye(3)/3)

# <codecell> Compute eigen decompistion of a tensor field
def eigenDecomposition(tensorField, nCells):
        
    eVals, eVecs = np.linalg.eig(tensorField)
    idx = eVals.argsort(axis=-1)[:,::-1]
    for c in range(nCells):
        eVals[c]=eVals[c,idx[c]]
        eVecs[c]=eVecs[c,:,idx[c]].T
    
    return eVals, eVecs

# <codecell> Perturbed samples of anisotropy tensor a_ij_PertSamples
def anisotropyTensorPert(eVec,eVal,nCells):
    a_ij = np.array([np.matmul(eVec[c], np.matmul(np.diag(eVal[c]),eVec[c].T)) 
                                        for c in range(nCells)])
    
    return a_ij

# <codecell> Sampling three r.v.s with zero mean 
def sampleThreeNormalRVsWithZeroMean(nSamples, sigma, seedValue):
    
    eValPert = np.zeros((nSamples, 3))
    
    dist = cp.Iid(cp.Normal(0,sigma), 2)
    
    np.random.seed(seedValue)
    eValPert[:, (0,1)] = dist.sample(nSamples, rule='L').T
    
    eValPert[:, 2] = -eValPert[:,0] -eValPert[:,1]
    
    # print(f'{eValPert}')
    # print(f'sum of eVal*s   = {eValPert.sum(axis=1)}')
    print(f'mean of samples = {eValPert.mean(axis=0)}')
    
    return eValPert
    
# <codecell> Add sampled perturbation to mean field of eigenvalues
def addPerturbationToEigenvalues(eVal, zeroMeanRVSamples, nSamples, nCells):
    eValPertSmaples = np.zeros((nSamples, nCells, 3))
    for s in range(nSamples):
        eValPertSmaples[s] = eVal  + zeroMeanRVSamples[s]
        # eValPertSmaples[s] = eVal * (1 + zeroMeanRVSamples[s])
    
    return eValPertSmaples