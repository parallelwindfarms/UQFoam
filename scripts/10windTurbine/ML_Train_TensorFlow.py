#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# %% Import Required Modules and env vars
from __future__ import print_function

import os, sys
import numpy as np

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '1'
import tensorflow as tf
import tensorflow.keras.callbacks as callbacks

cwd = os.getcwd() + "/"

SCRIPTS = os.environ['SCRIPTS']
DATA    = os.environ['windTurbineData']

sys.path.append(SCRIPTS+'/10windTurbine')
sys.path.append(SCRIPTS+'/pythonScripts')
sys.path.append(SCRIPTS+'/10windTurbine'+'/TensorFlowLib')

import myUQlib
from ML_GetFuncs import getMeshData, getCaseData
from ML_Models import UNet, UNetAug
from ML_Utils import enableGPUMemGro, dataGenerator, relErrL1, relErrL2

enableGPUMemGro()

randomSeed = 42
np.random.seed = 42

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
Wdx, d_D = D/60., 7.
x_up_D_WT = 1.
ADloc = lambda WT_num: (0 + (WT_num-1) *d_D*D, 0, h)

# Projection mesh params
# mlMeshName, nx, ny, nz = 'M128/', 549, 128, 128
mlMeshName, nx, ny, nz = 'M64/', 275, 64, 64

# %% Read OpenFOAM case mesh access to grid and cell values
case = samplesDir+'baseCase/'
# case = samplesDir+'sample_0/'

# case = '/projects/0/uqPint/10windTurbine/2RRSTF/ML/'+\
#     '5_HornRevRow_refineLocal_realKE_65Dx10Dx4.3D_WdxD_060_WdyzD_15_varScl_eigUinTIPert_6WT'+\
#     '/baseCase/'

vtkFile = case+'project2MLMesh_'+mlMeshName+'VTK/project2MLMesh_'+mlMeshName[:-1]+'_0.vtk'

for WT_num in range(1,2):
# for WT_num in range(1,7):
    print(WT_num)
    mesh, nCells, mlMeshShape, nCells_WT, cCenter_WT, \
        cCenterWT_idx, startPlane_WT_idx, endPlane_WT_idx, \
        y0Plane_WT_idx, zhPlane_WT_idx, cellsInDiskAtHubHeight = \
        getMeshData(vtkFile, D, h, ADloc(WT_num), ny, nz, d_D, x_up_D_WT)
        
    UMagDet, UHubDet, defUDet, tkeDet, TIDet, TIHubDet, RDet, ADet, nutDet = \
        getCaseData(myUQlib, mesh, nCells_WT, cCenterWT_idx, 
                        cellsInDiskAtHubHeight, Uref)
        
    print(UHubDet, TIHubDet)

# %% Data generator
sampleRegx = r'sample_*'
generator = dataGenerator(mlDataDir, mlMeshName, sampleRegx,
                          mlMeshShape, trainFrac=0.8, batchSize=2)
fileList = generator.fileList

# %% Model
isUNet, isUNetAug, isAutoEnc_NN = 1, 0, 0

if isUNet:
    model = UNet(mlMeshShape, nRandFieldsChannels=8, dropFrac=0.25)
    modelName = mlDir+'models/UNet_'+mlMeshName[:-1]+'.h5'
    ioData = generator.UNetIOData
    train_data, valid_data, test_data = generator.UNetIOBatchedSplitData    
elif isUNetAug:
    model = UNetAug(mlMeshShape, dropFrac=0.)
    modelName = mlDir+'models/UNetAug_'+mlMeshName[:-1]+'.h5'
    ioData = generator.UNetAugIOData
    train_data, valid_data, test_data = generator.UNetAugIOBatchedSplitData
elif isAutoEnc_NN:
    pass
    
model.summary()

# %% Load and test the model
# dependencies = {'relErrL1': relErrL1, 'relErrL2': relErrL2}
# model = tf.keras.models.load_model(modelName, custom_objects=dependencies)

# %% Trian the model
epochs = 100

s = len(train_data)//5
lr = 1e-3
# lrS = tf.keras.optimizers.schedules.ExponentialDecay(lr, s, 0.9)
opt = tf.keras.optimizers.Adam(lr, beta_1=0.9, beta_2=0.999)

cbs = [callbacks.ModelCheckpoint(modelName, save_best_only=True),
        callbacks.EarlyStopping(patience=10, monitor='relErrL1')]
# cbs = None

model.compile(optimizer=opt, loss='mae', metrics=[relErrL1, relErrL2])

results = model.fit(train_data.shuffle(len(train_data)),
                    validation_data=valid_data,
                    callbacks=cbs, epochs=epochs)

# %% Prediction
sampleRegx_pred = r'sample_900'
generator_pred = dataGenerator(mlDataDir, mlMeshName, sampleRegx_pred,
                          mlMeshShape, trainFrac=1., batchSize=1)
fileList_pred = generator_pred.fileList
data_pred = generator_pred.UNetIOBatchedSplitData

U_true = [ i[1][0,:,:,:].numpy() for i in data_pred ][0][:,:,:,0]
TI_true = [ i[1][0,:,:,:].numpy() for i in data_pred ][0][:,:,:,1]

U_pred = model.predict(data_pred)[0][:,:,:,0]
TI_pred = model.predict(data_pred)[0][:,:,:,1]
















