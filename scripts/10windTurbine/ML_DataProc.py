#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# <codecell> Import Required Modules and env vars
import os, sys
import numpy as np
import pyvista as pv
import pickle
from tqdm import tqdm

from multiprocessing import Pool
from functools import partial
from timeit import default_timer as timer
from ML_Utils import getMeshData, getCaseData

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
x_up_D_WT = 1.
ADloc = (0,0,h)

# # Original mesh params
# nx, ny, nz = int(5*30), int(4*1.25*15), int(1.6*(1.25+7/8.)*15)

# Projection mesh params
# mlMeshName, nx, ny, nz = 'M128/', 549, 128, 128
mlMeshName, nx, ny, nz = 'M64/', 275, 64, 64

# <codecell> Read OpenFOAM case mesh access to grid and cell values
# case = samplesDir+'baseCase/'
case = samplesDir+'sample_0/'

vtkFile = case+'project2MLMesh_'+mlMeshName+\
    'VTK/project2MLMesh_'+mlMeshName[:-1]+'_0.vtk'
mesh, nCells, mlMeshShape, nCells_WT, cCenter_WT, \
    cCenterWT_idx, startPlane_WT_idx, endPlane_WT_idx, \
    y0Plane_WT_idx, zhPlane_WT_idx, cellsInDiskAtHubHeight = \
    getMeshData(vtkFile, D, h, ADloc, ny, nz, d_D, x_up_D_WT)
    
UMagDet, UHubDet, defUDet, tkeDet, TIDet, TIHubDet, RDet, ADet, nutDet = \
    getCaseData(myUQlib, mesh, nCells_WT, cCenterWT_idx, 
                    cellsInDiskAtHubHeight, Uref)

# <codecell> Preprocess and dump samples
def getDataInWT(cCenterWT_idx, s):
    extn = '.pkl'

    vtkFile = samplesDir+\
        'sample_'+str(s)+'/project2MLMesh_'+mlMeshName+\
            'VTK/project2MLMesh_'+mlMeshName[:-1]+'_0.vtk'
    mesh = pv.UnstructuredGrid(vtkFile)
    
    UMag, UHub, defU, tke, TI, TIHub, R, A, nut = \
        getCaseData(myUQlib, mesh, nCells_WT, cCenterWT_idx, 
                        cellsInDiskAtHubHeight, Uref)

    pickle.dump(UHub,
                open(mlDataDir+mlMeshName+'sample_'+str(s)+'/UHub'+extn,'wb'))
    pickle.dump(UMag.reshape(*mlMeshShape, 1),
                open(mlDataDir+mlMeshName+'sample_'+str(s)+'/UMag'+extn,'wb'))    
    pickle.dump(defU.reshape(*mlMeshShape, 1),
                open(mlDataDir+mlMeshName+'sample_'+str(s)+'/defU'+extn,'wb'))

    pickle.dump(TIHub,
                open(mlDataDir+mlMeshName+'sample_'+str(s)+'/TIHub'+extn,'wb'))
    pickle.dump(tke.reshape(*mlMeshShape, 1),
                open(mlDataDir+mlMeshName+'sample_'+str(s)+'/k'+extn,'wb'))
    pickle.dump(TI.reshape(*mlMeshShape, 1),
                open(mlDataDir+mlMeshName+'sample_'+str(s)+'/TI'+extn,'wb'))

    pickle.dump(nut.reshape(*mlMeshShape, 1),
                open(mlDataDir+mlMeshName+'sample_'+str(s)+'/nut'+extn,'wb'))
    pickle.dump(R.reshape(*mlMeshShape, 6),
                open(mlDataDir+mlMeshName+'sample_'+str(s)+'/R'+extn,'wb'))
    pickle.dump(A.reshape(*mlMeshShape, 6),
                open(mlDataDir+mlMeshName+'sample_'+str(s)+'/A'+extn,'wb'))      

samplesRange = range(0, 1000)

pool = Pool(processes=32)
start = timer()
list(tqdm(
        pool.imap(partial(getDataInWT, cCenterWT_idx), samplesRange),
        total=len(samplesRange))
)
pool.close()
pool.join()
print(timer()-start, 's')

# <codecell> Load a batch [Testing]
batchSize = 1000
loadSampleRange = range(0,1000)#np.random.randint(50, 60, batchSize)

UHub  = np.zeros((batchSize))
TIHub = np.zeros((batchSize))
defU = np.zeros((batchSize, *mlMeshShape, 1))
UMag = np.zeros((batchSize, *mlMeshShape, 1))
TI  = np.zeros((batchSize, *mlMeshShape, 1))
tke = np.zeros((batchSize, *mlMeshShape, 1))
nut = np.zeros((batchSize, *mlMeshShape, 1))
R = np.zeros((batchSize, *mlMeshShape, 6))
A = np.zeros((batchSize, *mlMeshShape, 6))

start = timer()
extn = '.pkl'
for i in tqdm(range(batchSize)):
    s = loadSampleRange[i]
    UHub[i] = pickle.load(open(mlDataDir+mlMeshName+'sample_'+str(s)+'/UHub'+extn,'rb'))
    # UMag[i] = pickle.load(open(mlDataDir+mlMeshName+'sample_'+str(s)+'/UMag'+extn,'rb'))
    # defU[i] = pickle.load(open(mlDataDir+mlMeshName+'sample_'+str(s)+'/defU'+extn,'rb'))
    TIHub[i] = pickle.load(open(mlDataDir+mlMeshName+'sample_'+str(s)+'/TIHub'+extn,'rb'))
    # tke[i] = pickle.load(open(mlDataDir+mlMeshName+'sample_'+str(s)+'/k'+extn,'rb'))
    # TI[i] = pickle.load(open(mlDataDir+mlMeshName+'sample_'+str(s)+'/TI'+extn,'rb'))
    # nut[i] = pickle.load(open(mlDataDir+mlMeshName+'sample_'+str(s)+'/nut'+extn,'rb'))
    # R[i] = pickle.load(open(mlDataDir+mlMeshName+'sample_'+str(s)+'/R'+extn,'rb'))
    # A[i] = pickle.load(open(mlDataDir+mlMeshName+'sample_'+str(s)+'/A'+extn,'rb'))

print(timer()-start, 's')

