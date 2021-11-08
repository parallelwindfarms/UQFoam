#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# %% Import Required Modules and env vars
from __future__ import print_function

import os, sys
import numpy as np
import matplotlib.pyplot as plt

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
from ML_Utils import enableGPUMemGro, dataGenerator, L1, L2

enableGPUMemGro()

randomSeed = 42
np.random.seed = 42

# %% Hyper-parameters
mlMeshName = 'M128/'
# mlMeshName = 'M64/'

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
if mlMeshName == 'M128/': nx, ny, nz = 549, 128, 128
elif mlMeshName == 'M64/': nx, ny, nz = 275, 64, 64

# %% Read OpenFOAM baseCase mesh access to grid and cell values
case = samplesDir+'baseCase/'
# case = samplesDir+'sample_0/'

# case = '/projects/0/uqPint/10windTurbine/2RRSTF/ML/'+\
#     '5_HornRevRow_refineLocal_realKE_65Dx10Dx4.3D_WdxD_060_WdyzD_15_'+\
#     'varScl_eigUinTIPert_6WT'+\
#     '/baseCase/'

vtkFile = case+'project2MLMesh_'+mlMeshName+'VTK/project2MLMesh_'+\
    mlMeshName[:-1]+'_0.vtk'

for WT_num in range(1,2):
# for WT_num in range(1,7):
    print('WT_num:', WT_num)
    mesh, nCells, mlMeshShape, nCells_WT, cCenter_WT, \
        cCenterWT_idx, startPlane_WT_idx, endPlane_WT_idx, \
        y0Plane_WT_idx, zhPlane_WT_idx, cellsInDiskAtHubHeight = \
        getMeshData(vtkFile, D, h, ADloc(WT_num), ny, nz, d_D, x_up_D_WT)
        
    UMagDet, UHubDet, defUDet, tkeDet, TIDet, TIHubDet, RDet, ADet, nutDet = \
        getCaseData(myUQlib, mesh, nCells_WT, cCenterWT_idx, 
                        cellsInDiskAtHubHeight, Uref)
        
    print('UHubDet TIHubDet:', UHubDet, TIHubDet)

# %% Data generator
fileNames = [mlDataDir+mlMeshName+'sample_'+str(i) for i in range(1000)] #IMP!
generator = dataGenerator(fileNames, mlMeshShape, batchSize=2)
fileList = generator.fileList

# %% Model
isUNet, isUNetAug = 1, 0
transposed = 1

if isUNet:
    if transposed:
        model = UNet(
            mlMeshShape, dropFrac=0.1, l1_lambda=1e-3, convType='transposed'
        )    
        modelName = mlDir+'models/UNet_DataGen_transposed_'+\
            mlMeshName[:-1]+'.h5'
    else:
        model = UNet(
            mlMeshShape, dropFrac=0.1, l1_lambda=1e-3, convType='upsampled'
        )    
        modelName = mlDir+'models/UNet_DataGen_upsampled_'+\
            mlMeshName[:-1]+'.h5'
    
    ioData = generator.UNetIOData
    train_data, valid_data, test_data = generator.UNetIOBatchedSplitData    
elif isUNetAug:
    if transposed:
        model = UNetAug(
            mlMeshShape, dropFrac=0.1, l1_lambda=1e-3, convType='transposed'
        )
        modelName = mlDir+'models/UNetAug_DataGen_transposed_'+\
            mlMeshName[:-1]+'.h5'
    else:
        model = UNetAug(
            mlMeshShape, dropFrac=0.1, l1_lambda=1e-3, convType='upsampled'
        )
        modelName = mlDir+'models/UNetAug_DataGen_upsampled_'+\
            mlMeshName[:-1]+'.h5'

    ioData = generator.UNetAugIOData
    train_data, valid_data, test_data = generator.UNetAugIOBatchedSplitData
    
model.summary()

s = len(train_data) * 10
lr = 1e-3
lrS = tf.keras.optimizers.schedules.ExponentialDecay(lr, s, 0.9)
opt = tf.keras.optimizers.Adam(lrS, beta_1=0.9, beta_2=0.999)

cbs = [callbacks.ModelCheckpoint(modelName, save_best_only=True),
       callbacks.EarlyStopping(patience=100, monitor='L1')]

model.compile(optimizer=opt, loss='mae', metrics=[L1, L2])

# %% Train the model
epochs = 10000

history = model.fit(
    train_data.shuffle(len(train_data)),
    validation_data=valid_data,
    callbacks=cbs, epochs=epochs
)

print('Latest lr =',opt._decayed_lr(tf.float32))

# %% Plot losses and errors
fig, ax = plt.subplots(ncols=2, nrows=1, figsize=(8,3), dpi=150)

history[['loss', 'val_loss']].plot(ax=ax[0])
history[['L1', 'L2', 'val_L1', 'val_L2']].plot(ax=ax[1])

ax[0].set_ylabel('MAE Loss')
ax[0].set_xlabel('Epoch')
ax[1].set_ylabel('Relative Error')
ax[1].set_xlabel('Epoch')

# %% Check few cases in test_data
load_model = 0

# From the Data Processing Step
UHub_mean, UHub_std = (6.279, 1.967)
TIHub_mean, TIHub_std = (12.969, 4.438)
UMagMax, TIMax, tkeMax = (14.616, 21.795, 8.815)

data = test_data
# data = generator.UNetIOData.batch(1).cache().prefetch(1)

# Load the model
if load_model:
    dependencies = {'L1': L1, 'L2': L2}
    loaded_model_name = mlDir+'models/UNet_DataGen_transposed_M64_best.h5'
    loaded_model = tf.keras.models.load_model(
        loaded_model_name, custom_objects=dependencies
    )
    y_pred = loaded_model.predict(data)
else:
    y_pred = model.predict(data)
    
check_idx = np.random.randint(0, len(data), 5)
for s, test_case in enumerate(data):
    if s in check_idx:
        print('##############################################\n')
        print('Smaple #',s,'\n')
        y_true = test_case[1][0]
        UMagTestTrue = y_true[:,:,:,0].numpy().reshape(-1)
        TITestTrue = y_true[:,:,:,1].numpy().reshape(-1)
        UMagTestPred = y_pred[s,:,:,:,0].reshape(-1)
        TITestPred = y_pred[s,:,:,:,1].reshape(-1)
        UMagDiff = np.abs(UMagTestTrue-UMagTestPred).reshape(-1)
        TIDiff = np.abs(TITestTrue-TITestPred).reshape(-1)
        
        if isUNetAug:
            UHubTest = test_case[0][1].numpy()[0,0]
            TIHubTest = test_case[0][1].numpy()[0,1]     
        else:
            UHubTest = test_case[0][0][0,0,0,6].numpy()
            TIHubTest = test_case[0][0][0,0,0,7].numpy()

        print(' UHub =', (UHubTest*UHub_std+UHub_mean).item()*100//10/10, \
              'TIHub =', (TIHubTest*TIHub_std+TIHub_mean).item()*100//10/10, 
              '\n'
        )
        
        print(' U:',
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

# %% Plotting
plot_soln = lambda ax, soln, plane, norm=None, cmap=None: ax.imshow(
    soln[plane].reshape(-1,mlMeshShape[0])[::-1], aspect='auto', 
    norm=norm, cmap=cmap, interpolation=None#'bicubic'
)
err_max_pct = 0.1

# Contour Plots
fig, ax = plt.subplots(ncols=3, nrows=4, constrained_layout=True, 
                       sharex=True, figsize=(16,6))
ax, CS = ax.flat, [0]*12

CS[0] = plot_soln(ax[0], UMagTestTrue, y0Plane_WT_idx)
CS[1] = plot_soln(ax[1], UMagTestPred, y0Plane_WT_idx)
norm = plt.Normalize(0, UMagTestTrue[y0Plane_WT_idx].max()*err_max_pct)
CS[2] = plot_soln(ax[2], UMagDiff, y0Plane_WT_idx, norm, cmap='gray')

CS[3] = plot_soln(ax[3], TITestTrue, y0Plane_WT_idx)
CS[4] = plot_soln(ax[4], TITestPred, y0Plane_WT_idx)
norm = plt.Normalize(0, TITestTrue[y0Plane_WT_idx].max()*err_max_pct)
CS[5] = plot_soln(ax[5], TIDiff, y0Plane_WT_idx, norm, cmap='gray')

CS[6] = plot_soln(ax[6], UMagTestTrue, zhPlane_WT_idx)
CS[7] = plot_soln(ax[7], UMagTestPred, zhPlane_WT_idx)
norm = plt.Normalize(0, UMagTestTrue[zhPlane_WT_idx].max()*err_max_pct)
CS[8] = plot_soln(ax[8], UMagDiff, zhPlane_WT_idx, norm, cmap='gray')

CS[9] = plot_soln(ax[9], TITestTrue, zhPlane_WT_idx)
CS[10] = plot_soln(ax[10], TITestPred, zhPlane_WT_idx)
norm = plt.Normalize(0, TITestTrue[zhPlane_WT_idx].max()*err_max_pct)
CS[11] = plot_soln(ax[11], TIDiff, zhPlane_WT_idx, norm, cmap='gray')

for i in range(12):
    ax[i].set_xticks([])
    ax[i].set_yticks([])
    if i in [1,4,7,10]:
        fig.colorbar(CS[i-1], ax=ax[i], aspect=50)
    else:
        fig.colorbar(CS[i], ax=ax[i], aspect=50)

ax[0].set_title('OpenFOAM (True)')
ax[1].set_title('U-Net (Pred)')
ax[2].set_title('|Error|')
ax[0].set_ylabel('UMag')
ax[3].set_ylabel('TI')
ax[6].set_ylabel('UMag')
ax[9].set_ylabel('TI')



