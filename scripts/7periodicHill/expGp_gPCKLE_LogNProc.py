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

vtkFile = 'expGp_PCEModes/VTK/expGp_PCEModes_0.vtk'
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

    Nx = 50; Ny = 20;

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

# <codecell> d=nKLEModes , g=KLEmodes
stD       = np.sqrt(variance) * np.ones(nCells) 
dist      = cp.Iid(cp.Normal(0,1), d)
polyExpsn = cp.orth_ttr(n, dist, retall=False, normed=True, cross_truncation=1/q)
# polyExpsn = cp.orth_ttr(n, dist, retall=False, normed=False, cross_truncation=1/q)
expoMat   = cp.bertran.bindex(start=0, stop=n, dim=d, cross_truncation=1/q)
P         = len(expoMat)-1
Pplus1    = P+1

# <codecell> PCE of logN field based on Ghanem 2002
expGpPCEModes, var = myUQlib.getlogNFieldPCE(mean,stD,nCells,d,g,volAvg,phyDim,expoMat,P,dist)
print(var)

# <codecell>  PCE of logN field using NIPC (for lower number of dim)
# nodes, wts = cp.generate_quadrature\
#                     (order=n, dist=dist, rule="gaussian", sparse=True)
# numQuadPts = len(wts)
# expGpEvals   = [myUQlib.getExpGpSample(g[0:d].reshape(d,nCells), node, d, nCells, stD) \
#                 for node in nodes.T]
# expGpPCEApp  = cp.fit_quadrature(polyExpsn, nodes, wts, expGpEvals)

# <codecell> PCE of logN field using LHS experiment (robust for any number of dim)
# sampleSize   = 5000
# xiSamples    = dist.sample(size=sampleSize, rule="L")
# expGpSamples = [myUQlib.getExpGpSample(g[0:d].reshape(d,nCells), node, d, nCells, stD) \
#               for node in xiSamples.T]
# expGpPCEApp  = cp.fit_regression(polyExpsn, xiSamples, expGpSamples)    

# <codecell>  Finding sobol indices
# sblIdx_i       = cp.Sens_m(expGpPCEApp, dist)
# sblIdx_iVolAvg = sblIdx_i * volAvg
# #sblIdx_iTotVolAvg = cp.Sens_t(nutPCEApp, dist) * volAvg

# for i in range(d):
#     print(i+1,'\t',sum(sblIdx_iVolAvg[i])*100, '%')
# print("\n")
# for i in range(d):
#     print(i+1,'\t',(sblIdx_i[i][int(nCells/2)])*100, '%')

# <codecell>  Finding PCE
# expGpPCEModes = np.array(expGpPCEApp.coefficients.T)
# var = np.sum(expGpPCEModes[1:]**2, axis=0)
# print(var)

# <codecell> Write the nut_PCEModes
for i in range(Pplus1):
    myUQlib.writeExpGpPCECoeffs(i,expGpPCEModes[i].reshape(nCells,1), \
                              nCells, "expGp_PCEModes/0/", phyDim)        

# for i in range(d):
#     myUQlib.writeExpGpPCECoeffs(i,g[i].reshape(nCells,1), \
#                               nCells, "expGp_PCEModes/0/", phyDim) 
        
# <codecell> Setting UQ info. for OpenFOAM simulation
# myUQlib.uqInfo(n, d, q, "REVF", dist)
