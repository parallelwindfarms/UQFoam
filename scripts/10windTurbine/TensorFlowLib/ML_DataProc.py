#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# %% Import Required Modules and env vars
import os, sys
import numpy as np
import pyvista as pv
import pickle
from tqdm import tqdm

from multiprocessing import Pool
from functools import partial
from timeit import default_timer as timer

cwd = os.getcwd() + "/"

SCRIPTS = os.environ['SCRIPTS']
DATA    = os.environ['windTurbineData']

sys.path.append(SCRIPTS+'/10windTurbine')
sys.path.append(SCRIPTS+'/pythonScripts')
sys.path.append(SCRIPTS+'/10windTurbine'+'/TensorFlowLib')

import myUQlib
from ML_GetFuncs import getMeshData, getCaseData

# %% Hyper-parameters
# mlMeshName = 'M128/'
mlMeshName = 'M64/'

# %% Case details
caseName   = cwd.split('/')[-2]
casePngDir = '/2RRSTF/ML/' + caseName+'/png'

if (os.path.exists(DATA+casePngDir))==False:
    print('Making new case data directory...')
    os.makedirs(DATA+casePngDir, exist_ok=True)

samplesDir = cwd+'DET_samples/'
mlDir = cwd+'ML_training/'
mlDataDir = mlDir+'data/'

# %% Case Parameters
D, h, Uref = 80, 70, 8
# Influence until next 3 WTs
Wdx, d_D = D/60., 7.*3.
x_up_D_WT = 1.0
ADloc = lambda WT_num: (0 + (WT_num-1) *d_D*D, 0, h)

# Projection mesh params
if mlMeshName == 'M128/': nx, ny, nz = 549, 128, 128
elif mlMeshName == 'M64/': nx, ny, nz = 275, 64, 64

# %% Read OpenFOAM baseCase mesh access to grid and cell values
case = samplesDir+'baseCase/'
# case = samplesDir+'sample_0/'

vtkFile = case+'project2MLMesh_'+mlMeshName+'VTK/project2MLMesh_'+mlMeshName[:-1]+'_0.vtk'

mesh, nCells, mlMeshShape, nCells_WT, cCenter_WT, \
    cCenterWT_idx, startPlane_WT_idx, endPlane_WT_idx, \
    y0Plane_WT_idx, zhPlane_WT_idx, cellsInDiskAtHubHeight = \
    getMeshData(vtkFile, D, h, ADloc(1), ny, nz, d_D, x_up_D_WT)
    
UMagDet, UHubDet, defUDet, tkeDet, TIDet, TIHubDet, RDet, ADet, nutDet, C_ve_DET = \
    getCaseData(myUQlib, mesh, nCells_WT, cCenterWT_idx, 
                    cellsInDiskAtHubHeight)
    
print('UHubDet TIHubDet:', UHubDet, TIHubDet)

# %% Preprocess and dump samples
def getDataInWT(cCenterWT_idx, s):
    extn = '.pkl'

    vtkFile = samplesDir+\
        'sample_'+str(s)+'/project2MLMesh_'+mlMeshName+\
            'VTK/project2MLMesh_'+mlMeshName[:-1]+'_0.vtk'
    mesh = pv.UnstructuredGrid(vtkFile)
    
    UMag, UHub, defU, tke, TI, TIHub, R, A, nut, C_vec = \
        getCaseData(myUQlib, mesh, nCells_WT, cCenterWT_idx, 
                        cellsInDiskAtHubHeight)
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

    pickle.dump(C_vec.reshape(*mlMeshShape, 2),
                open(mlDataDir+mlMeshName+'sample_'+str(s)+'/C_vec'+extn,'wb'))    

samplesRange = range(0, 1000)

pool = Pool(processes=64)
start = timer()
list(tqdm(
        pool.imap(partial(getDataInWT, cCenterWT_idx), samplesRange),
        total=len(samplesRange))
)
pool.close()
pool.join()
print(timer()-start, 's')

# %% Load a batch and dump it
batchSize = 1000
loadSampleRange = range(0,1000)

UHub  = np.zeros((batchSize))
TIHub = np.zeros((batchSize))
defU = np.zeros((batchSize, *mlMeshShape, 1))
UMag = np.zeros((batchSize, *mlMeshShape, 1))
TI  = np.zeros((batchSize, *mlMeshShape, 1))
tke = np.zeros((batchSize, *mlMeshShape, 1))
nut = np.zeros((batchSize, *mlMeshShape, 1))
R = np.zeros((batchSize, *mlMeshShape, 6))
A = np.zeros((batchSize, *mlMeshShape, 6))
C_vec = np.zeros((batchSize, *mlMeshShape, 2))

start = timer()
extn = '.pkl'
for i in tqdm(range(batchSize)):
    s = loadSampleRange[i]
    UHub[i] = pickle.load(open(mlDataDir+mlMeshName+'sample_'+str(s)+'/UHub'+extn,'rb'))
    UMag[i] = pickle.load(open(mlDataDir+mlMeshName+'sample_'+str(s)+'/UMag'+extn,'rb'))
    defU[i] = pickle.load(open(mlDataDir+mlMeshName+'sample_'+str(s)+'/defU'+extn,'rb'))
    TIHub[i] = pickle.load(open(mlDataDir+mlMeshName+'sample_'+str(s)+'/TIHub'+extn,'rb'))
    tke[i] = pickle.load(open(mlDataDir+mlMeshName+'sample_'+str(s)+'/k'+extn,'rb'))
    TI[i] = pickle.load(open(mlDataDir+mlMeshName+'sample_'+str(s)+'/TI'+extn,'rb'))
    nut[i] = pickle.load(open(mlDataDir+mlMeshName+'sample_'+str(s)+'/nut'+extn,'rb'))
    R[i] = pickle.load(open(mlDataDir+mlMeshName+'sample_'+str(s)+'/R'+extn,'rb'))
    A[i] = pickle.load(open(mlDataDir+mlMeshName+'sample_'+str(s)+'/A'+extn,'rb'))
    C_vec[i] = pickle.load(open(mlDataDir+mlMeshName+'sample_'+str(s)+'/C_vec'+extn,'rb'))

pickle.dump(UHub, open(mlDataDir+mlMeshName+'UHub'+extn,'wb'))
pickle.dump(UMag, open(mlDataDir+mlMeshName+'UMag'+extn,'wb'))
pickle.dump(defU, open(mlDataDir+mlMeshName+'defU'+extn,'wb'))
pickle.dump(TIHub, open(mlDataDir+mlMeshName+'TIHub'+extn,'wb'))
pickle.dump(tke, open(mlDataDir+mlMeshName+'k'+extn,'wb'))
pickle.dump(TI, open(mlDataDir+mlMeshName+'TI'+extn,'wb'))
pickle.dump(nut, open(mlDataDir+mlMeshName+'nut'+extn,'wb'))
pickle.dump(R, open(mlDataDir+mlMeshName+'R'+extn,'wb'))
pickle.dump(A, open(mlDataDir+mlMeshName+'A'+extn,'wb'))
pickle.dump(C_vec, open(mlDataDir+mlMeshName+'C_vec'+extn,'wb'))
    
print(timer()-start, 's')

# %% Train-Test Split
dataSize = len(A)
trainSize = int(dataSize * 0.8)
valSize = trainSize + (dataSize-trainSize) // 2

ATrain = A[:trainSize]

UMagTrain, tkeTrain, TITrain, nutTrain, RTrain, ATrain, UHubTrain, TIHubTrain, C_vecTrain = \
    UMag[:trainSize], tke[:trainSize], TI[:trainSize], nut[:trainSize], \
    R[:trainSize], A[:trainSize], UHub[:trainSize], TIHub[:trainSize], C_vec[:trainSize]
    
UMagVal, tkeVal, TIVal, nutVal, RVal, AVal, UHubVal, TIHubVal, C_vecVal = \
    UMag[trainSize:valSize], tke[trainSize:valSize], TI[trainSize:valSize], nut[trainSize:valSize], \
    R[trainSize:valSize], A[trainSize:valSize], UHub[trainSize:valSize], TIHub[trainSize:valSize], C_vec[trainSize:valSize]    
    
UMagTest, tkeTest, TITest, nutTest, RTest, ATest, UHubTest, TIHubTest, C_vecTest = \
    UMag[valSize:], tke[valSize:], TI[valSize:], nut[valSize:], \
    R[valSize:], A[valSize:], UHub[valSize:], TIHub[valSize:], C_vec[valSize:]
    
print(UMagTrain.shape, UMagVal.shape, UMagTest.shape)

del nut, R, A, UMag, tke, UHub, TIHub, C_vec

# %% Data Statistics
UMagMin, tkeMin, TIMin= UMagTrain.min(), tkeTrain.min(), TITrain.min()
UMagMax, tkeMax, TIMax = UMagTrain.max(), tkeTrain.max(), TITrain.max()
nutMean, nutStd = nutTrain.mean(axis=0), nutTrain.std(axis=0)
RMean, RStd = RTrain.mean(axis=0), RTrain.std(axis=0)
AMean, AStd = ATrain.mean(axis=0), ATrain.std(axis=0)
C_vecMean, C_vecStd = C_vecTrain.mean(axis=0), C_vecTrain.std(axis=0)
UHubMean, UHubStd = UHubTrain.mean(), UHubTrain.std()
TIHubMean, TIHubStd = TIHubTrain.mean(), TIHubTrain.std()

pickle.dump(UMagMax, open(mlDataDir+mlMeshName+'UMagMax'+extn,'wb'))
pickle.dump(tkeMax, open(mlDataDir+mlMeshName+'tkeMax'+extn,'wb'))
pickle.dump(TIMax, open(mlDataDir+mlMeshName+'TIMax'+extn,'wb'))

pickle.dump(UMagMin, open(mlDataDir+mlMeshName+'UMagMin'+extn,'wb'))
pickle.dump(tkeMin, open(mlDataDir+mlMeshName+'tkeMin'+extn,'wb'))
pickle.dump(TIMin, open(mlDataDir+mlMeshName+'TIMin'+extn,'wb'))

pickle.dump(nutMean, open(mlDataDir+mlMeshName+'nutMean'+extn,'wb'))
pickle.dump(nutStd, open(mlDataDir+mlMeshName+'nutStd'+extn,'wb'))

pickle.dump(RMean, open(mlDataDir+mlMeshName+'RMean'+extn,'wb'))
pickle.dump(RStd, open(mlDataDir+mlMeshName+'RStd'+extn,'wb'))

pickle.dump(AMean, open(mlDataDir+mlMeshName+'AMean'+extn,'wb'))
pickle.dump(AStd, open(mlDataDir+mlMeshName+'AStd'+extn,'wb'))

pickle.dump(C_vecMean, open(mlDataDir+mlMeshName+'C_vecMean'+extn,'wb'))
pickle.dump(C_vecStd, open(mlDataDir+mlMeshName+'C_vecStd'+extn,'wb'))

pickle.dump(UHubMean, open(mlDataDir+mlMeshName+'UHubMean'+extn,'wb'))
pickle.dump(UHubStd, open(mlDataDir+mlMeshName+'UHubStd'+extn,'wb'))

pickle.dump(TIHubMean, open(mlDataDir+mlMeshName+'TIHubMean'+extn,'wb'))
pickle.dump(TIHubStd, open(mlDataDir+mlMeshName+'TIHubStd'+extn,'wb'))
