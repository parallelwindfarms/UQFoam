#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# %% Import Required Modules and env vars
from __future__ import print_function

import os, sys
import numpy as np
import pyvista as pv

# Tensorflow and other local imports
import tensorflow as tf
from ML_UNetAug import UNetAug
from ML_Utils import enableGPUMemGro, dataGenerator, relErrL1, relErrL2
import tensorflow.keras.callbacks as callbacks

SCRIPTS = os.environ['SCRIPTS']
DATA    = os.environ['windTurbineData']

sys.path.append(SCRIPTS+'/10windTurbine')
sys.path.append(SCRIPTS+'/pythonScripts')

cwd = os.getcwd() + "/"
import myUQlib

enableGPUMemGro()
randomSeed = 42

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
Wdx, d_D, nWT = D/60., 7, 1
ADloc = (0,0,h)

# Projection mesh params
# mlMeshName, nx, ny, nz = 'M128/', 549, 128, 128 # M128
mlMeshName, nx, ny, nz = 'M64/', 275, 64, 64

# %% Read OpenFOAM case mesh access to grid and cell values
# case = 'baseCase/'
case = 'sample_0/'

vtkFile = samplesDir+case+'project2MLMesh_'+mlMeshName+'VTK/project2MLMesh_0.vtk'
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
mlMeshShape = tuple([nz, ny, nx_WT])

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

# %% Data generator
sampleRegx = r'sample_*'
generator = dataGenerator(mlDataDir, mlMeshName, sampleRegx,
                          mlMeshShape, trainFrac=0.8, batchSize=2)
fileList = generator.fileList
ioData = generator.ioData
train_data, valid_data, test_data = generator.ioBatchedSplitData

# %% Load the model
model = UNetAug(mlMeshShape)
model.summary()
modelName = mlDir+'models/UNetAug_'+mlMeshName[:-1]+'.h5'

# %% Load and test the model
# dependencies = {'scaledRes': scaledRes}
# model = tf.keras.models.load_model(modelName, custom_objects=dependencies)

# %% Trian the model
epochs = 100

s = len(train_data)//5
lr = 1e-3
lrS = tf.keras.optimizers.schedules.ExponentialDecay(lr, s, 0.9)
opt = tf.keras.optimizers.Adam(lrS, beta_1=0.5, beta_2=0.999)

cbs = [callbacks.ModelCheckpoint(modelName, save_best_only=True),
       callbacks.EarlyStopping(patience=10, monitor='relErrL2')]

model.compile(optimizer=opt, loss='mae', metrics=[relErrL1, relErrL2])

results = model.fit(train_data.shuffle(len(train_data)),
                    validation_data=valid_data,
                    callbacks=cbs, epochs=epochs)

