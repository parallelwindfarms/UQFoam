#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# %% Import Required Modules and Paths to Modules
from __future__ import print_function

import os, sys
import numpy as np
import pickle 

import tensorflow as tf

cwd = os.getcwd() + "/"

SCRIPTS = os.environ['SCRIPTS']
DATA    = os.environ['windTurbineData']

sys.path.append(SCRIPTS+'/pythonScripts')
sys.path.append(SCRIPTS+'/10windTurbine'+'/TensorFlowLib')

import myUQlib
from ML_GetFuncs import getMeshData, getCaseData
from ML_Utils import enableGPUMemGro, set_global_determinism
from ML_Utils import L1, L2, makePlots

# %% Set Env Vars and Global Settings
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '1'
os.environ['CUDA_VISIBLE_DEVICES'] = '2'
set_global_determinism()
enableGPUMemGro()

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
mlDir = '../4_Vestas2MWV_refineLocal_realKE_30Dx10Dx4.3D_'+\
    'WdxD_060_WdyzD_15_varScl_eigUinTIPert_1WT/ML_training/'
mlDataDir = mlDir+'data/'

# %% Model
modelName = mlDir+'models/UNet_transposed_Aij_d_D_21_x_up_D_WT_0.01_batch_10_M64.h5'
dependencies = {'L1': L1, 'L2': L2}
model = tf.keras.models.load_model(modelName, custom_objects=dependencies)
model.trainable = False
model.summary()

# %% Case Parameters
D, h = 80.0, 70.0
# Influence until next 3 WTs
d_D, x_up_D_WT = 7*3, 0.01 #1.0
ADloc = lambda WT_num: (0 + (WT_num-1) * 7.0*D, 0, h)

# Projection mesh params
if mlMeshName == 'M128/': ny, nz = 128, 128
elif mlMeshName == 'M64/': ny, nz = 64, 64

# %% From the Data Processing Step
UHubMean = pickle.load(open(mlDataDir+mlMeshName+'UHubMean.pkl', 'rb'))
UHubStd = pickle.load(open(mlDataDir+mlMeshName+'UHubStd.pkl', 'rb'))
TIHubMean = pickle.load(open(mlDataDir+mlMeshName+'TIHubMean.pkl', 'rb'))
TIHubStd = pickle.load(open(mlDataDir+mlMeshName+'TIHubStd.pkl', 'rb'))

# %% Read OpenFOAM 1 WT baseCase for Aij
case = '/projects/0/uqPint/10windTurbine/2RRSTF/ML/'+\
    '4_Vestas2MWV_refineLocal_realKE_30Dx10Dx4.3D_WdxD_060_WdyzD_15_'+\
    'varScl_eigUinTIPert_1WT'+'/DET_samples/baseCase/'

vtkFile = case+'project2MLMesh_'+mlMeshName+'VTK/project2MLMesh_'+\
    mlMeshName[:-1]+'_0.vtk'

for WT_num in [1]:
# for WT_num in range(1,7):
    print(f'{WT_num = }')
    mesh, nCells, mlMeshShape, nCells_WT, cCenter_WT, \
        cCenterWT_idx, startPlane_WT_idx, endPlane_WT_idx, \
        y0Plane_WT_idx, zhPlane_WT_idx, cellsInDiskAtHubHeight = \
        getMeshData(vtkFile, D, h, ADloc(WT_num), ny, nz, d_D, x_up_D_WT)
        
    _,_,_,_,_,_, _, ADet, _,_ = getCaseData(
            myUQlib, mesh, nCells_WT, cCenterWT_idx, 
            cellsInDiskAtHubHeight
    )
        
# %% Read OpenFOAM 6 WT sample access mesh and cell values
case = samplesDir+'baseCase/'
# case = samplesDir+'sample_0/'

vtkFile = case+'project2MLMesh_'+mlMeshName+'VTK/project2MLMesh_'+\
    mlMeshName[:-1]+'_0.vtk'

## %% Initialize Solution
print('## Initializing the solution ##')
n = 6
globalShape = (mlMeshShape[0],mlMeshShape[1],mlMeshShape[2]*n)
UMagTrue = np.zeros(globalShape)
UMagPred = np.zeros(globalShape)
TITrue = np.zeros(globalShape)
TIPred = np.zeros(globalShape)
defUTrue = np.zeros(globalShape)
defUPred = np.zeros(globalShape)
TITrue = np.zeros(globalShape)
TIPred = np.zeros(globalShape)

## %% Superposition Model  
for WT_num in range(1,7):
    print(f'\n{WT_num = }')
    mesh, nCells, mlMeshShape, nCells_WT, cCenter_WT, \
        cCenterWT_idx, startPlane_WT_idx, endPlane_WT_idx, \
        y0Plane_WT_idx, zhPlane_WT_idx, cellsInDiskAtHubHeight = \
        getMeshData(vtkFile, D, h, ADloc(WT_num), ny, nz, d_D, x_up_D_WT)
        
    UMagDet, UHubDet, defUDet, tkeDet, TIDet, TIHubDet, \
        RDet, _, nutDet, C_ve_DET = getCaseData(
            myUQlib, mesh, nCells_WT, cCenterWT_idx, 
            cellsInDiskAtHubHeight
    )
    print(f'Original: {UHubDet = }, {TIHubDet = }')
        
    if WT_num==1:
        Uin, TIin = UHubDet, TIHubDet
    # else:
    #     UHubDet = UMagPred[:,:,mlMeshShape[0]*(WT_num-2):mlMeshShape[0]*(WT_num-2)+mlMeshShape[2]]\
    #         [:,:,::-1].reshape(-1)[cellsInDiskAtHubHeight].mean() + 0.1
    #     TIHubDet = TIPred[:,:,mlMeshShape[0]*(WT_num-2):mlMeshShape[0]*(WT_num-2)+mlMeshShape[2]]\
    #         [:,:,::-1].reshape(-1)[cellsInDiskAtHubHeight].mean()
    #     print(f'Updated:  {UHubDet = }, {TIHubDet = }')

    UHubStand = (UHubDet-UHubMean)/UHubStd
    TIHubStand = (TIHubDet-TIHubMean)/TIHubStd

    ## %% Tensorflow inputs --> y_pred
    inputFields = np.concatenate(
        (
         ADet,
         np.ones_like(ADet[:,:1])*UHubStand,
         np.ones_like(ADet[:,:1])*TIHubStand
        ),
        axis=-1
    )
    inputFields =  tf.convert_to_tensor(
        inputFields.reshape(1,*mlMeshShape,-1), dtype=tf.float64
    )
    y_pred = model.predict(inputFields)
    
    ## %% Superposition Model
    UMagTestTrue = UMagDet
    UMagTestPred = y_pred[0,:,:,:,0].reshape(-1)
    UMagTestDiff = np.abs(UMagTestTrue-UMagTestPred).reshape(-1)
    
    TITestTrue = TIDet
    TITestPred = y_pred[0,:,:,:,1].reshape(-1)
    TITestDiff = np.abs(TITestTrue-TITestPred).reshape(-1)
    
    defUTestTrue = 1 - UMagTestTrue/Uin
    # defUTestPred = 1 - UMagTestPred/Uin
    defUTestPred = 1 - UMagTestPred/UHubDet
    defUTestDiff = np.abs(defUTestTrue-defUTestPred).reshape(-1)
 
    startIdx = int(mlMeshShape[2]*(WT_num-1)/3)
    endIdx = startIdx + mlMeshShape[2]
    
    if WT_num in [1,4]: 
        UMagTrue[:,:,startIdx:endIdx] = UMagTestTrue.reshape(mlMeshShape)
        TITrue[:,:,startIdx:endIdx] = TITestTrue.reshape(mlMeshShape)
        defUTrue[:,:,startIdx:endIdx] = defUTestTrue.reshape(mlMeshShape)
    
    offset = 0 if WT_num>1 else 0
    defUPred[:,:,startIdx+offset:endIdx] = \
        defUTestPred.reshape(mlMeshShape)[:,:,offset:]
    # defUPred[:,:,startIdx+offset:endIdx] += \
    #     defUTestPred.reshape(mlMeshShape)[:,:,offset:]
    defUDiff = np.abs(defUTrue-defUPred)

    # UMagPred = (1-defUPred)*Uin
    UMagPred[:,:,startIdx+offset:endIdx] = \
        UMagTestPred.reshape(mlMeshShape)[:,:,offset:]    
    UMagDiff = np.abs(UMagTrue-UMagPred)
    
    TIPred[:,:,startIdx+offset:endIdx] = \
        TITestPred.reshape(mlMeshShape)[:,:,offset:]    
    TIDiff = np.abs(TITrue-TIPred)

# %% Check few cases in test_data
_,_, mlPlotShape, _,_, _,_,_, y0Plane, zhPlane, _ = \
    getMeshData(vtkFile, D, h, ADloc(1), ny, nz, 7.0*6, x_up_D_WT)
    
makePlots(
    'test', mlPlotShape, y0Plane, zhPlane, 
    defUTrue[:,:,:mlPlotShape[2]].reshape(-1), 
    defUPred[:,:,:mlPlotShape[2]].reshape(-1),
    defUDiff[:,:,:mlPlotShape[2]].reshape(-1),
    
    TITrue[:,:,:mlPlotShape[2]].reshape(-1), 
    TIPred[:,:,:mlPlotShape[2]].reshape(-1),
    TIDiff[:,:,:mlPlotShape[2]].reshape(-1),
)
    
makePlots(
    'test', mlPlotShape, y0Plane, zhPlane, 
    UMagTrue[:,:,:mlPlotShape[2]].reshape(-1), 
    UMagPred[:,:,:mlPlotShape[2]].reshape(-1),
    UMagDiff[:,:,:mlPlotShape[2]].reshape(-1),
    
    TITrue[:,:,:mlPlotShape[2]].reshape(-1), 
    TIPred[:,:,:mlPlotShape[2]].reshape(-1),
    TIDiff[:,:,:mlPlotShape[2]].reshape(-1),
)

# %% End
print('Program ran successfully!\n')