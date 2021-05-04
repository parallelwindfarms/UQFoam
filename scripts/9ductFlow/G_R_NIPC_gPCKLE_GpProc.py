# -*- coding: utf-8 -*-
"""
Purpose of this python script -
(1) Create a KL Expansion of a Gaussian field/process
(2) Compute the PC Expansion of a Random SPD Matrix field
(3) Compute inner products of polynomials
"""

# <codecell> Import Required Modules and env vars
from __future__ import print_function

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

vtkFile = 'VTK/0_4quarters_dp0.1_0.vtk'
mesh    = pv.UnstructuredGrid(cwd + vtkFile)
nCells  = mesh.n_cells
cCenter = mesh.cell_centers().points
cBounds = mesh.cell_centers().bounds

if (phyDim==2):

    Nx = 64; Ny =Nx;

    myIndices = [Nx-1, Ny-1]
    mesher    = ot.IntervalMesher(myIndices)
    kleBounds = ot.Interval([cBounds[0], cBounds[2]], 
                            [cBounds[1], cBounds[3]])
    kleMesh = mesher.build(kleBounds)
    kleMesh.setVertices(cCenter[:,[0,1]])

if (phyDim==3):
        
    Nx = 5; Ny = Nx; Nz = 40;

    myIndices = [Nx-1, Ny-1, Nz-1]
    mesher    = ot.IntervalMesher(myIndices)
    kleBounds = ot.Interval([cBounds[0], cBounds[2], cBounds[4]], 
                            [cBounds[1], cBounds[3], cBounds[5]])
    kleMesh = mesher.build(kleBounds)
    kleMesh.setVertices(cCenter[:,[0,1,2]])
    

# <codecell> Loading R_det and finding Cholesky Decomposition
R_det = mesh.cell_arrays['R'][:,:]
R_det_mat = np.zeros((nCells,3,3))
R_det_mat[:,0,0] = R_det[:,0]
R_det_mat[:,0,1] = R_det_mat[:,1,0] = R_det[:,1]
R_det_mat[:,0,2] = R_det_mat[:,2,0] = R_det[:,2]

R_det_mat[:,1,1] = R_det[:,3]
R_det_mat[:,1,2] = R_det_mat[:,2,1] = R_det[:,4]

R_det_mat[:,2,2] = R_det[:,5]

L_R = np.zeros(R_det_mat.shape)
L_R[:] = np.linalg.cholesky(R_det_mat[:])
for i in range(nCells):
    L_R[i] = L_R[i].T

# <codecell> KLE algorithm
D  = 2
h  = D/2

lx = 2 * h 
ly = 2 * h 
lz = 1 * h 

l  = [lx, ly]
if phyDim==3 :
    l.append(lz)

covModel  = ot.SquaredExponential(l, [1])
eigValK_0 = 1e-2

# <codecell> Set KLE Algorithm Karhunen Loeve P1/SVD Algorithm
kleAlgo = ot.KarhunenLoeveP1Algorithm(kleMesh, covModel, eigValK_0)
# sample  = ot.GaussianProcess(covModel, kleMesh).getSample(1000)
# kleAlgo = ot.KarhunenLoeveSVDAlgorithm(sample, eigValK_0)


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
q = 1

phyDim  = 3 # for 2D

# <codecell> Dispersion parameter 
dp     = 0.2 * np.ones(nCells)
stD_d  = dp/(np.sqrt(phyDim+1))
stD_ij = stD_d
kU_ii  = np.zeros((phyDim,nCells))
stD_ii = np.zeros((phyDim,nCells))
for phyDim_i in range(phyDim):
    kU_ii[phyDim_i]  = 1/2 * stD_d**-2 - phyDim_i/2
    stD_ii[phyDim_i] = np.sqrt(4*stD_d**4 * kU_ii[phyDim_i])


# <codecell> Joint dist of L_ij and L_ii
dist = cp.Normal(0,1)
for d_i in range(d-1):
    dist = cp.J(dist, cp.Normal(0,1))
# dist = cp.Iid(cp.Normal(0,1), d)

phi,norms = cp.orth_ttr(n,dist,retall=True,normed=True,cross_truncation=1/q)    
Pplus1    = len(phi)


# <codecell> Coeffs of PCE of Gamma dist wrt Normal dist
n_u = 5
U   = np.zeros((phyDim, nCells, n_u+1))

phi_u,norms_u = cp.orth_ttr(n_u,cp.Normal(0,1),retall=True,normed=True)    

for phyDim_i in range(phyDim):
    for cellNr in range(1):#nCells):
        U[phyDim_i,:] = myUQlib.getPCEGamma(cp.Gamma(kU_ii[phyDim_i,cellNr]), n_u)
        
        
# <codecell> LHS sampling experiment method - Sampling fields from Gaussian process
sampleSize = 1000
xiSamples  = dist.sample(size=sampleSize, rule="L")
gpSamples  = [myUQlib.getGpSample(g[0:d], node, d, nCells, 1) \
                  for node in xiSamples[0:d].T]


# <codecell> NIPC qudratue method - Compute fields from Gaussian process
# xiSamples,wts = cp.generate_quadrature(order=n,dist=dist,rule="gaussian",sparse=True)
# sampleSize    = len(wts)
# gpSamples     = [myUQlib.getGpSample(g[0:d], node, d, nCells, 1) \
#                     for node in xiSamples[0:d].T]
    
    
# <codecell> Sampling L_ij, L_ii and L fields
L_ijSamples = [myUQlib.getGpSample(g[0:d], node, d, nCells, stD_d) \
                for node in xiSamples[0:d].T]
    
L_iiSamples = [[stD_d * np.sqrt(2*np.sum(U[phyDim_i] * \
                 phi_u(gpSamples[s]).T, 1)) for phyDim_i in \
                 range(phyDim)] for s in range(sampleSize)]
    
# <codecell> Sobol indices of modes of the random field
# L_ijPCEApp   = cp.fit_regression(phi, xiSamples, L_ijSamples)
# L_ijPCEModes = L_ijPCEApp.coefficients
    
# sblIdx_i       = cp.Sens_m(L_ijPCEApp, dist)
# sblIdx_iVolAvg = sblIdx_i * volAvg
# #sblIdx_iTotVolAvg = cp.Sens_t(nutPCEApp, dist) * volAvg

# for i in range(d):
#     print(i+1,'\t',sum(sblIdx_iVolAvg[i])*100, '%')
# print("\n")
# for i in range(d):
#     print(i+1,'\t',(sblIdx_i[i][int(nCells/2)])*100, '%')
    

# <codecell> Sampling G and contructing PCE of G
L_Samples = np.zeros((((sampleSize,nCells,3,3))))
for idx_i in range(phyDim):
    for idx_j in range(phyDim):
        if (idx_i < idx_j):
        # if (idx_i < idx_j) and (idx_i < 2) and (idx_j < 2): # for 2D
            L_Samples[:,:,idx_i,idx_j] = np.array(L_ijSamples)[:,:]
        if idx_i == idx_j:
            L_Samples[:,:,idx_i,idx_i] = np.array(L_iiSamples)[:,idx_i,:]

G_Samples = np.zeros((((sampleSize,nCells,3,3))))
for s in range(sampleSize):
    for c in range(nCells):
        G_Samples[s,c] = np.matmul(L_Samples[s,c].T, L_Samples[s,c])

G_PCEModes = np.zeros((((nCells,Pplus1,3,3))))
for idx_i in range(phyDim):
    for idx_j in range(phyDim):
        G_PCEApp = cp.fit_regression(phi, xiSamples, G_Samples[:,:,idx_i,idx_j])
        # G_PCEApp = cp.fit_quadrature(phi, xiSamples, wts, G_Samples[:,:,idx_i,idx_j])
        G_PCEModes[:,:,idx_i,idx_j] = G_PCEApp.coefficients        

G_PCEModes_mean = G_PCEModes[:,0]
G_PCEModes_stD = np.sum(G_PCEModes[:,1:]**2, 1)**0.5

G_Samples_mean = np.mean(G_Samples,axis=0)
G_Samples_std  = np.std(G_Samples,axis=0)


# <codecell> Sampling R using G
R_Samples = np.zeros(G_Samples.shape)
for s in range(sampleSize):
    for c in range(nCells):
        R_Samples[s,c] = np.matmul(np.matmul(L_R[c].T, G_Samples[s,c]), L_R[c])

R_Samples_mean = np.mean(R_Samples, axis=0)
R_Samples_std  = np.std(R_Samples, axis=0)


# <codecell> Writing samples of R
phyDim  = 2
for s in range(sampleSize):
    myUQlib.writeRSamplesDuct(s, R_Samples[s], nCells, "RSamples/", phyDim)


# <codecell> Write the G_PCEModes
# phyDim  = 2
# for idx_k in range(Pplus1):
#     myUQlib.writeGmatPCEModesDuct(idx_k, G_PCEModes[:,idx_k], nCells, "0/", phyDim)
    
    
# <codecell> Setting UQ info. for OpenFOAM simulation
# myUQlib.uqInfo(n, d, q, "RRSTF",dist)

