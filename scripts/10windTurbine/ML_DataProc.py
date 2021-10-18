#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# <codecell> Import Required Modules and env vars
import os, sys
import numpy as np
import pyvista as pv
import pickle
import tqdm

from multiprocessing import Pool
from functools import partial
from timeit import default_timer as timer

SCRIPTS = os.environ['SCRIPTS']
DATA    = os.environ['windTurbineData']

sys.path.append(SCRIPTS+'/10windTurbine')
sys.path.append(SCRIPTS+'/pythonScripts')

cwd = os.getcwd() + "/"

import myUQlib

# <codecell> Case details
caseName   = cwd.split('/')[-2]
casePngDir = '/2RRSTF/ML/' + caseName+'/png'

if (os.path.exists(DATA+casePngDir))==False:
    print('Making new case data directory...')
    os.makedirs(DATA+casePngDir, exist_ok=True)

samplesDir = cwd+'DET_samples/'
mlDir = cwd+'ML_training/'
mlDataDir = mlDir+'data/'

# <codecell>
# Parametrs
D, h, Uref = 80, 70, 8
Wdx, d_D, nWT = D/60., 7, 1
ADloc = (0,0,h)

# # Original mesh params
# nx, ny, nz = int(5*30), int(4*1.25*15), int(1.6*(1.25+7/8.)*15)

# Projection mesh params
# mlMeshName, nx, ny, nz = 'M128/', 549, 128, 128 # M128
mlMeshName, nx, ny, nz = 'M64/', 275, 64, 64

# <codecell> Read OpenFOAM case mesh access to grid and cell values
# vtkFile = samplesDir+'baseCase/project2MLMesh/VTK/project2MLMesh_0.vtk'
vtkFile = samplesDir+'sample_57/project2MLMesh/VTK/project2MLMesh_0.vtk'

mesh    = pv.UnstructuredGrid(vtkFile)
nCells  = mesh.n_cells
cCenter = np.array(mesh.cell_centers().points)

x_up_D_WT = 1
xStart_WT, xEnd_WT = ADloc[0]-x_up_D_WT*D, ADloc[0]+(7.-x_up_D_WT)*D
cCenterWT_idx = np.where(
    (cCenter[:,0] >= xStart_WT) & (cCenter[:,0] <= xEnd_WT))[0]
cCenter_WT = cCenter[cCenterWT_idx]
nCells_WT = cCenter_WT.shape[0]
nx_WT = int(nCells_WT/ny/nz)
assert (nx_WT==nz==ny)
ML_meshShape = tuple([nz, ny, nx_WT])

startPlane_WT_idx = (cCenter_WT[:,0]==cCenter_WT[:,0].min())
endPlane_WT_idx   = (cCenter_WT[:,0]==cCenter_WT[:,0].max())

# <codecell> Read WT cell values for baseCase
UMagDet = np.linalg.norm(mesh.cell_data['U'][cCenterWT_idx], axis=1)
UHubDet = UMagDet[startPlane_WT_idx].mean()
defUDet = (UHubDet-UMagDet)/UHubDet

tkeDet   = np.array(mesh.cell_data['k'][cCenterWT_idx])
TIAddDet = np.sqrt(tkeDet*2/3)/Uref *100
TIHubDet = TIAddDet[startPlane_WT_idx].mean()
TIAddDet = TIAddDet - TIHubDet

RDet = myUQlib.symmTensorToTensorv2021(
    mesh.cell_data['turbulenceProperties:R'][cCenterWT_idx], nCells_WT)
ADet = myUQlib.anisotropyTensor(RDet, tkeDet, nCells_WT)
ADet = myUQlib.getSymmTensorValues(ADet, nCells_WT)

del UMagDet, tkeDet, RDet
del mesh, cCenter

# <codecell> Preprocess and dump samples
def getDataInWT(cCenterWT_idx, s):
    extn = '.pkl'

    vtkFile = samplesDir+\
        'sample_'+str(s)+'/project2MLMesh/VTK/project2MLMesh_0.vtk'
    mesh = pv.UnstructuredGrid(vtkFile)

    UMag = np.linalg.norm(mesh.cell_data['U'][cCenterWT_idx], axis=1)
    UHub = UMag[startPlane_WT_idx].mean()
    defU = (UHub-UMag)/UHub
    pickle.dump(UHub,
                open(mlDataDir+mlMeshName+'sample_'+str(s)+'/UHub'+extn,'wb'))
    pickle.dump(defU.reshape(*ML_meshShape, 1),
                open(mlDataDir+mlMeshName+'sample_'+str(s)+'/defU'+extn,'wb'))

    tke   = np.array(mesh.cell_data['k'][cCenterWT_idx])
    TIAdd = np.sqrt(tke*2/3)/Uref *100
    TIHub = TIAdd[startPlane_WT_idx].mean()
    TIAdd = TIAdd - TIHub
    pickle.dump(TIHub,
                open(mlDataDir+mlMeshName+'sample_'+str(s)+'/TIHub'+extn,'wb'))
    pickle.dump(TIAdd.reshape(*ML_meshShape, 1),
                open(mlDataDir+mlMeshName+'sample_'+str(s)+'/TIAdd'+extn,'wb'))

    R = myUQlib.symmTensorToTensorv2021(
        mesh.cell_data['turbulenceProperties:R'][cCenterWT_idx], nCells_WT)
    A = myUQlib.anisotropyTensor(R, tke, nCells_WT)
    A = myUQlib.getSymmTensorValues(A, nCells_WT)
    pickle.dump(A.reshape(*ML_meshShape, 6),
                open(mlDataDir+mlMeshName+'sample_'+str(s)+'/A'+extn,'wb'))

samplesRange = range(0, 1000)

pool = Pool(processes=10)
start = timer()
list(tqdm.tqdm(
        pool.imap(partial(getDataInWT, cCenterWT_idx), samplesRange),
        total=len(samplesRange)))
pool.close()
pool.join()
print(timer()-start, 's')

# <codecell> Load a batch [Testing]
batchSize = 100
loadSampleRange = range(0,100)#np.random.randint(50, 60, batchSize)

UHub = np.zeros((batchSize))
defU = np.zeros((batchSize, *ML_meshShape, 1))
TIHub = np.zeros((batchSize))
TIAdd = np.zeros((batchSize, *ML_meshShape, 1))
A = np.zeros((batchSize, *ML_meshShape, 6))

start = timer()
extn = '.pkl'
for i in tqdm.tqdm(range(batchSize)):
    s = loadSampleRange[i]
    print(s)
    UHub[i] = pickle.load(open(mlDataDir+mlMeshName+'sample_'+str(s)+'/UHub'+extn,'rb'))
    defU[i] = pickle.load(open(mlDataDir+mlMeshName+'sample_'+str(s)+'/defU'+extn,'rb'))
    TIHub[i] = pickle.load(open(mlDataDir+mlMeshName+'sample_'+str(s)+'/TIHub'+extn,'rb'))
    TIAdd[i] = pickle.load(open(mlDataDir+mlMeshName+'sample_'+str(s)+'/TIAdd'+extn,'rb'))
    A[i] = pickle.load(open(mlDataDir+mlMeshName+'sample_'+str(s)+'/A'+extn,'rb'))

print(timer()-start, 's')


# <codecell>


# <codecell>


# <codecell>


# <codecell>


# <codecell>
