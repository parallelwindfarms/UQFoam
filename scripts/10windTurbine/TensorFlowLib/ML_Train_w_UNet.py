#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# %% Import Required Modules and Paths to Modules
from __future__ import print_function

import os, sys
import numpy as np
import pandas as pd
import pickle 
import matplotlib.pyplot as plt

import tensorflow as tf
import tensorflow.keras.callbacks as callbacks

cwd = os.getcwd() + "/"

SCRIPTS = os.environ['SCRIPTS']
DATA    = os.environ['windTurbineData']

sys.path.append(SCRIPTS+'/pythonScripts')
sys.path.append(SCRIPTS+'/10windTurbine'+'/TensorFlowLib')

import myUQlib
from ML_GetFuncs import getMeshData, getCaseData
from ML_Model_UNet import UNet
from ML_Utils import enableGPUMemGro, set_global_determinism
from ML_Utils import dataGenerator, L1, L2, makePlots

# %% Set Env Vars and Global Settings
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '1'
os.environ['CUDA_VISIBLE_DEVICES'] = '2'
set_global_determinism()
enableGPUMemGro()

# %% Hyper-parameters
# mlMeshName, BATCH_SIZE = 'M128/', 2
mlMeshName, BATCH_SIZE = 'M64/', 10

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
D, h = 80.0, 70.0
# Influence until next 3 WTs
d_D, x_up_D_WT = 7*3, 0.01 #1.0
ADloc = lambda WT_num: (0 + (WT_num-1) * 7.0*D, 0, h)

# Projection mesh params
if mlMeshName == 'M128/': ny, nz = 128, 128
elif mlMeshName == 'M64/': ny, nz = 64, 64

# %% Read OpenFOAM baseCase mesh access to grid and cell values
# case = samplesDir+'baseCase/'
case = samplesDir+'sample_0/'

# case = '/projects/0/uqPint/10windTurbine/2RRSTF/ML/'+\
#     '5_HornRevRow_refineLocal_realKE_65Dx10Dx4.3D_WdxD_060_WdyzD_15_'+\
#     'varScl_eigUinTIPert_6WT/DET_samples'+\
#     '/baseCase/'

vtkFile = case+'project2MLMesh_'+mlMeshName+'VTK/project2MLMesh_'+\
    mlMeshName[:-1]+'_0.vtk'

for WT_num in [1]:
# for WT_num in range(1,7):
    print(f'{WT_num = }')
    mesh, nCells, mlMeshShape, nCells_WT, cCenter_WT, \
        cCenterWT_idx, startPlane_WT_idx, endPlane_WT_idx, \
        y0Plane_WT_idx, zhPlane_WT_idx, cellsInDiskAtHubHeight = \
        getMeshData(vtkFile, D, h, ADloc(WT_num), ny, nz, d_D, x_up_D_WT)
        
    UMagDet, UHubDet, defUDet, tkeDet, TIDet, TIHubDet, \
        RDet, ADet, nutDet, C_ve_DET = getCaseData(
            myUQlib, mesh, nCells_WT, cCenterWT_idx, 
            cellsInDiskAtHubHeight
    )
    
    print('UHubDet TIHubDet:', UHubDet, TIHubDet, '\n')

# %% Data generator
fileNames = [mlDataDir+mlMeshName+'sample_'+str(i) for i in range(1000)] #IMP!
generator = dataGenerator(
    mlDataDir, mlMeshName, fileNames, mlMeshShape, batchSize=BATCH_SIZE
)
fileList = generator.fileList

# %% Model
transposed, channels, l1_lambda, dropFrac = 1, 64, 1e-3, 0.1
convType = 'transposed' if transposed else 'upsampled'

model = UNet(
    mlMeshShape, dropFrac=dropFrac, channels=channels,
    l1_lambda=l1_lambda, convType=convType
) 
modelName = mlDir+'models/UNet_'+convType+'_Aij_d_D_'+str(d_D)+'_x_up_D_WT_'+str(x_up_D_WT)+\
    '_batch_'+str(BATCH_SIZE)+'_'+mlMeshName[:-1]+'.h5'

ioData = generator.UNetIOData
train_data, valid_data, test_data = generator.UNetIOBatchedSplitData
    
model.summary()

# %% Model Parameters
epochs = 2000
s = len(train_data) * 20
lr = 1e-3
lrS = tf.keras.optimizers.schedules.ExponentialDecay(lr, s, 0.9)
opt = tf.keras.optimizers.Adam(lrS, beta_1=0.9, beta_2=0.999)
cbs = [callbacks.ModelCheckpoint(modelName, save_best_only=True),
       callbacks.EarlyStopping(patience=1000, monitor='L1')]
model.compile(optimizer=opt, loss='mae', metrics=[L1, L2])

# %% Train the model
model.fit(train_data.shuffle(len(train_data)), 
          validation_data=valid_data,
          callbacks=cbs, epochs=epochs
)

# %% Get latest lr and plot losses and errors
print('Latest lr =',opt._decayed_lr(tf.float32))

fig, ax = plt.subplots(ncols=3, nrows=1, figsize=(15,4), dpi=150,
                       sharex=True)

try: df
except NameError: df = pd.DataFrame(model.history.history)
 
df[['loss', 'val_loss']].plot(ax=ax[0])
df[['L1', 'val_L1']].plot(ax=ax[1])
df[['L2', 'val_L2']].plot(ax=ax[2])

ax[0].set_ylabel('MAE Loss')
ax[0].set_xlabel('Epoch')
ax[1].set_ylabel('Relative Error')
ax[1].set_xlabel('Epoch')
ax[2].set_ylabel('Relative Error')
ax[2].set_xlabel('Epoch')

ax[0].set_ylim(0,2)
ax[1].set_ylim(0,1)
ax[2].set_ylim(0,1)

pickle.dump(df, open(modelName[:-3]+'_history_df', 'wb'))

# %% Check few cases in test_data
load_model = 0 ###
data = test_data

# From the Data Processing Step
UHubMean, UHubStd = generator.UHubMean, generator.UHubStd
TIHubMean, TIHubStd = generator.TIHubMean, generator.TIHubStd

# Load the model
if load_model:
    dependencies = {'L1': L1, 'L2': L2}
    loaded_model_name = modelName
    loaded_model = tf.keras.models.load_model(
        loaded_model_name, custom_objects=dependencies
    )
    y_pred = loaded_model.predict(data)
else:
    y_pred = model.predict(data)
    
check_idx = np.random.randint(0, len(data), 50)
# check_idx = [49,  3,  1,  5, 53, 25, 88, 59, 40, 28]

for s, test_case in enumerate(data):
    if s in check_idx:
        print('#'*100+'\n')
        print('Smaple #',s,'\n')
        y_true = test_case[1][0]
        UMagTestTrue = y_true[:,:,:,0].numpy().reshape(-1)
        TITestTrue = y_true[:,:,:,1].numpy().reshape(-1)
        UMagTestPred = y_pred[s,:,:,:,0].reshape(-1)
        TITestPred = y_pred[s,:,:,:,1].reshape(-1)
        UMagDiff = np.abs(UMagTestTrue-UMagTestPred).reshape(-1)
        TIDiff = np.abs(TITestTrue-TITestPred).reshape(-1)
        
        UHubTest = test_case[0][0][0,0,0,6].numpy()
        TIHubTest = test_case[0][0][0,0,0,7].numpy()

        print(' UHub =',(UHubTest*UHubStd+UHubMean).item()*100//10/10, \
              'TIHub =',(TIHubTest*TIHubStd+TIHubMean).item()*100//10/10, 
              '\n'
        )
        
        print(' U: ',
              'L1 Error =',f'{L1(UMagTestTrue,UMagTestPred)*100:.1f} %',
              'L2 Error =',f'{L2(UMagTestTrue,UMagTestPred)*100:.1f} %'
        )
        print(' TI:',
              'L1 Error =', f'{L1(TITestTrue,TITestPred)*100:.1f} %',
              'L2 Error =', f'{L2(TITestTrue,TITestPred)*100:.1f} %'
        )
        
        random_idx = np.random.randint(0, len(UMagTestPred), 6)
        print('\n', 'U:  True', UMagTestTrue[random_idx],\
              '\n', 'U:  Pred', UMagTestPred[random_idx])
        print('\n', 'TI: True', TITestTrue[random_idx],\
              '\n', 'TI: Pred', TITestPred[random_idx], '\n')

        makePlots(
            s, mlMeshShape, y0Plane_WT_idx, zhPlane_WT_idx, 
            UMagTestTrue, UMagTestPred, UMagDiff/UMagTestTrue,
            TITestTrue, TITestPred, TIDiff/TITestTrue
        )
        
# %% End
print('Program ran successfully!\n')