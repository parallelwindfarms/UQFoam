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

vtkFile = 'G_PCEModes/VTK/G_PCEModes_0.vtk'
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

    Nx = 32; Ny = 12;

    myIndices = [Nx-1, Ny-1]
    mesher    = ot.IntervalMesher(myIndices)
    kleBounds = ot.Interval([cBounds[0], cBounds[2]], 
                            [cBounds[1], cBounds[3]])
    kleMesh = mesher.build(kleBounds)
    kleMesh.setVertices(cCenter[:,[0,1]])

if (phyDim==3):
        
    Nx = 0; Ny = 0; Nz = 0;

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

covModel  = ot.SquaredExponential(l, [1])
eigValK_0 = 1e-2#1e-2#1e-4


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
q = 0.1

phyDim  = 3 # for 2D

# <codecell> Dispersion parameter 
dp     = 0.3 * np.ones(nCells)
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
# sampleSize = 1000
# xiSamples  = dist.sample(size=sampleSize, rule="L")
# gpSamples  = [myUQlib.getGpSample(g[0:d], node, d, nCells, 1) \
#                   for node in xiSamples[0:d].T]


# <codecell> NIPC qudratue method - Compute fields from Gaussian process
xiSamples,wts = cp.generate_quadrature(order=n,dist=dist,rule="gaussian",sparse=True)
sampleSize    = len(wts)
gpSamples     = [myUQlib.getGpSample(g[0:d], node, d, nCells, 1) \
                    for node in xiSamples[0:d].T]
    
    
# <codecell> Sampling L_ij, L_ii and L fields and finding solbol indices
L_ijSamples = [myUQlib.getGpSample(g[0:d], node, d, nCells, stD_d) \
                for node in xiSamples[0:d].T]
    
L_iiSamples = [[stD_d * np.sqrt(2*np.sum(U[phyDim_i] * \
                 phi_u(gpSamples[sampleNr]).T, 1)) for phyDim_i in \
                 range(phyDim)] for sampleNr in range(sampleSize)]
    
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
    

# <codecell> Sampling G field and contructing PCE of G
L_Samples = np.zeros((((sampleSize,nCells,3,3))))
for idx_i in range(phyDim):
    for idx_j in range(phyDim):
        # if (idx_i < idx_j):
        if (idx_i < idx_j) and (idx_i < 2) and (idx_j < 2): # for 2D
            L_Samples[:,:,idx_i,idx_j] = np.array(L_ijSamples)[:,:]
        if idx_i == idx_j:
            L_Samples[:,:,idx_i,idx_i] = np.array(L_iiSamples)[:,idx_i,:]

G_Samples = np.zeros((((sampleSize,nCells,3,3))))
for sample in range(sampleSize):
    for cellNr in range(nCells):
        G_Samples[sample,cellNr] = np.matmul(L_Samples[sample,cellNr].T, \
                                            L_Samples[sample,cellNr])

G_PCEModes = np.zeros((((nCells,Pplus1,3,3))))
for idx_i in range(phyDim):
    for idx_j in range(phyDim):
        # G_PCEApp = cp.fit_regression(phi, xiSamples, G_Samples[:,:,idx_i,idx_j])
        G_PCEApp = cp.fit_quadrature(phi, xiSamples, wts, G_Samples[:,:,idx_i,idx_j])
        G_PCEModes[:,:,idx_i,idx_j] = G_PCEApp.coefficients        

stD_G = np.sum(G_PCEModes[:,1:]**2, 1)**0.5


# <codecell> Write the G_PCEModes
phyDim  = 2 # for 2D

for idx_k in range(Pplus1):
    myUQlib.writeGmatPCEModes(idx_k, G_PCEModes[:,idx_k], nCells, "G_PCEModes/0/", phyDim)


# <codecell> Setting UQ info. for OpenFOAM simulation
# myUQlib.uqInfo(n, d, q, "RRSTF")
























# <codecell> ALT-Sampling and contructing PCE for L_ij, L_ii and G fields

# phi3 = cp.outer(phi, phi, phi)
# Cijk = cp.E(phi3, dist)

# L_iiPCEApp   = cp.fit_regression(phi, xiSamples, L_iiSamples)
# L_iiPCEModes = L_iiPCEApp.coefficients

# L_PCEModes = np.zeros((((nCells,Pplus1,3,3))))
# for idx_i in range(phyDim):
#     for idx_j in range(phyDim):
#         L_PCEApp = cp.fit_regression(phi, xiSamples, L_Samples[:,:,idx_i,idx_j])
#         L_PCEModes[:,:,idx_i,idx_j] = L_PCEApp.coefficients 

# G_PCEModes_alt = np.zeros((((nCells,Pplus1,3,3))))

# for cellNr in range(nCells):
#     for idx_k,idx_l,idx_m in itr.product(*[range(Pplus1)]*3):
#         tmpCijk = Cijk[idx_k][idx_l][idx_m]
#         if abs(tmpCijk) >= 1e-12:
#             G_PCEModes_alt[cellNr,idx_k] += \
#                 np.matmul(L_PCEModes[cellNr,idx_l].T, L_PCEModes[cellNr,idx_m]) * tmpCijk

# for idx_k in range(Pplus1):
#     G_PCEModes_alt[:,idx_k] = G_PCEModes_alt[:,idx_k] / norms[idx_k]

# G_std_alt = np.sum(G_PCEModes_alt[:,1:]**2, 1)**0.5