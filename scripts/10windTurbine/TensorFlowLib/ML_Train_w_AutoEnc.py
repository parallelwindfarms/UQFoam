#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# %% Import Required Modules and Env Vars
from __future__ import print_function

import os, sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import tensorflow as tf
import tensorflow.keras.callbacks as callbacks

cwd = os.getcwd() + "/"

SCRIPTS = os.environ['SCRIPTS']
DATA    = os.environ['windTurbineData']

sys.path.append(SCRIPTS+'/pythonScripts')
sys.path.append(SCRIPTS+'/10windTurbine'+'/TensorFlowLib')

import myUQlib
from ML_GetFuncs import getMeshData, getCaseData
from ML_Model_AutoEnc import AutoEnc, LatentSpaceRepr, TransposedCNN
from ML_Utils import enableGPUMemGro, set_global_determinism
from ML_Utils import dataGenerator, L1, L2

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '1'
os.environ['CUDA_VISIBLE_DEVICES'] = '1'
# set_global_determinism(27)
enableGPUMemGro()

# %% Hyper-parameters
# mlMeshName, BATCH_SIZE = 'M128/', 2
mlMeshName, BATCH_SIZE = 'M64/', 2

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
        getMeshData(
            vtkFile, D, h, ADloc(WT_num), 
            ny, nz, d_D, x_up_D_WT
        )
        
    UMagDet, UHubDet, defUDet, tkeDet, TIDet, \
        TIHubDet, RDet, ADet, nutDet, C_vecDet = \
        getCaseData(
            myUQlib, mesh, nCells_WT, cCenterWT_idx, 
            cellsInDiskAtHubHeight, Uref
        )
        
    print('UHubDet TIHubDet:', UHubDet, TIHubDet)

# %% Data generator
fileNames = [mlDataDir+mlMeshName+'sample_'+str(i) for i in range(1000)] #IMP!
generator = dataGenerator(fileNames, mlMeshShape, batchSize=BATCH_SIZE)
fileList = generator.fileList

# %% Autoencoder Model for Anisotropy
autoencoder = AutoEnc(
    mlMeshShape, nRandFieldsChannels=2, nOutChannels=2, nLatentChannels=10,
    dropFrac=0.1, l1_lambda=1e-3, convType='transposed', latent=True
)    
autoencoder.summary()
modelName = mlDir+'models/AutoEnc_transposed_C_vec_'+ mlMeshName[:-1]+'.h5'
ioData = generator.UNetIOData
train_data, valid_data, test_data = generator.AutoEncIOBatchedSplitData 
    
# %% Autoencoder Model Parameters
epochs = 200
s = len(train_data) * 20
lr = 1e-3
lrS = tf.keras.optimizers.schedules.ExponentialDecay(lr, s, 0.9)
opt = tf.keras.optimizers.Adam(lrS, beta_1=0.9, beta_2=0.999)
cbs = [callbacks.ModelCheckpoint(modelName, save_best_only=True),
       callbacks.EarlyStopping(patience=100, monitor='L1')]
autoencoder.compile(optimizer=opt, loss='mae', metrics=[L1, L2])

# %% Train the Autoencoder Model
autoencoder.fit(train_data.shuffle(len(train_data)), 
          validation_data=valid_data,
          callbacks=cbs, epochs=epochs
)

# %% Get latest lr and plot Losses and Errors
print('Latest lr =',opt._decayed_lr(tf.float32))

fig, ax = plt.subplots(ncols=2, nrows=1, figsize=(10,4), dpi=150)

try: autoencoder_df
except NameError: autoencoder_df = pd.DataFrame(autoencoder.history.history)
 
autoencoder_df[['loss', 'val_loss']].plot(ax=ax[0])
autoencoder_df[['L1', 'L2', 'val_L1', 'val_L2']].plot(ax=ax[1])

ax[0].set_ylabel('MAE Loss')
ax[0].set_xlabel('Epoch')
ax[1].set_ylabel('Relative Error')
ax[1].set_xlabel('Epoch')

# %% TransposedCNN Model
autoencoder.trainable = False
model = TransposedCNN(
    mlMeshShape, autoencoder, nHubData=2, nOutChannels=2, dropFrac=0.1, 
    channels=64, l1_lambda=1e-3, convType='transposed'
)    
model.summary()
modelName = mlDir+'models/TransposedCNN_C_vec_'+ mlMeshName[:-1]+'.h5'
ioData = generator.UNetIOData
train_data, valid_data, test_data = generator.TransposedCNNIOBatchedSplitData 

# %% TransposedCNN Model Parameters
epochs = 1000
s = len(train_data) * 20
lr = 1e-3
lrS = tf.keras.optimizers.schedules.ExponentialDecay(lr, s, 0.9)
opt = tf.keras.optimizers.Adam(lrS, beta_1=0.9, beta_2=0.999)
cbs = [callbacks.ModelCheckpoint(modelName, save_best_only=True),
       callbacks.EarlyStopping(patience=100, monitor='L1')]
model.compile(optimizer=opt, loss='mae', metrics=[L1, L2])

# %% Train the TransposedCNN Model
model.fit(train_data.shuffle(len(train_data)), 
          validation_data=valid_data,
          callbacks=cbs, epochs=epochs
)

# %% Get latest lr and plot Losses and Errors
print('Latest lr =',opt._decayed_lr(tf.float32))

fig, ax = plt.subplots(ncols=2, nrows=1, figsize=(10,4), dpi=150)

try: df
except NameError: df = pd.DataFrame(model.history.history)
 
df[['loss', 'val_loss']].plot(ax=ax[0])
df[['L1', 'L2', 'val_L1', 'val_L2']].plot(ax=ax[1])

ax[0].set_ylabel('MAE Loss')
ax[0].set_xlabel('Epoch')
ax[1].set_ylabel('Relative Error')
ax[1].set_xlabel('Epoch')

# %% Check few cases in test_data
# From the Data Processing Step
UHub_mean, UHub_std = (6.279, 1.967)
TIHub_mean, TIHub_std = (12.969, 4.438)
UMagMax, TIMax, tkeMax = (14.616, 21.795, 8.815)

data = train_data#test_data
# y_pred = autoencoder.predict(data)
y_pred = model.predict(data)
    
# check_idx = np.random.randint(0, len(data), 5)
check_idx = [49,  3,  1,  5, 53, 25, 88, 59, 40, 28]

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
        
        UHubTest = test_case[0][1][0,0].numpy()
        TIHubTest = test_case[0][1][0,1].numpy()

        print(' UHub =',(UHubTest*UHub_std+UHub_mean).item()*100//10/10, \
              'TIHub =',(TIHubTest*TIHub_std+TIHub_mean).item()*100//10/10, 
              '\n'
        )
        # latentVars = LatentSpaceRepr(autoencoder)(test_case[0])    
        latentVars = LatentSpaceRepr(autoencoder)(test_case[0][0])    
        print(' Latent Space Representation of Anisotropy\n', 
              latentVars.numpy().reshape(-1)
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

        ## %% Plotting
        plot_soln = lambda ax, soln, plane, norm=None: ax.imshow(
            soln[plane].reshape(-1,mlMeshShape[0])[::-1], 
            aspect='auto', 
            norm=norm
        )
        err_pct = 0.1
        
        # Contour Plots
        fig, ax = plt.subplots(ncols=3, nrows=4, constrained_layout=True, 
                               sharex=True, figsize=(16,6))
        ax, CS = ax.flat, [0]*12
        
        CS[0] = plot_soln(ax[0], UMagTestTrue, y0Plane_WT_idx)
        CS[1] = plot_soln(ax[1], UMagTestPred, y0Plane_WT_idx)
        norm = plt.Normalize(0, UMagTestTrue[y0Plane_WT_idx].max()*err_pct)
        CS[2] = plot_soln(ax[2], UMagDiff, y0Plane_WT_idx, norm)
        
        CS[3] = plot_soln(ax[3], TITestTrue, y0Plane_WT_idx)
        CS[4] = plot_soln(ax[4], TITestPred, y0Plane_WT_idx)
        norm = plt.Normalize(0, TITestTrue[y0Plane_WT_idx].max()*err_pct)
        CS[5] = plot_soln(ax[5], TIDiff, y0Plane_WT_idx, norm)
        
        CS[6] = plot_soln(ax[6], UMagTestTrue, zhPlane_WT_idx)
        CS[7] = plot_soln(ax[7], UMagTestPred, zhPlane_WT_idx)
        norm = plt.Normalize(0, UMagTestTrue[zhPlane_WT_idx].max()*err_pct)
        CS[8] = plot_soln(ax[8], UMagDiff, zhPlane_WT_idx, norm)
        
        CS[9] = plot_soln(ax[9], TITestTrue, zhPlane_WT_idx)
        CS[10] = plot_soln(ax[10], TITestPred, zhPlane_WT_idx)
        norm = plt.Normalize(0, TITestTrue[zhPlane_WT_idx].max()*err_pct)
        CS[11] = plot_soln(ax[11], TIDiff, zhPlane_WT_idx, norm)
        
        for i in range(12):
            ax[i].set_xticks([])
            ax[i].set_yticks([])
            if i in [1,4,7,10]:
                fig.colorbar(CS[i-1], ax=ax[i], aspect=50)
            else:
                fig.colorbar(CS[i], ax=ax[i], aspect=50)
        
        ax[0].set_title('OpenFOAM (True)')
        ax[1].set_title('U-Net (Pred)')
        ax[2].set_title('|True - Pred|')
        ax[0].set_ylabel('UMag')
        ax[3].set_ylabel('TI')
        ax[6].set_ylabel('UMag')
        ax[9].set_ylabel('TI')
        
# %% End
print('Program ran successfully!\n')