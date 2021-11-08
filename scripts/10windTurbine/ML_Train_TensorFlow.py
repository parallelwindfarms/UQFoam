#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# %% Import Required Modules and env vars
from __future__ import print_function

import os, sys
import numpy as np
from timeit import default_timer as timer
import pickle

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
from ML_Utils import enableGPUMemGro, relErrL1, relErrL2

enableGPUMemGro()

randomSeed = 42
np.random.seed = 42

# %% Hyper-parameters
tdtype = tf.float32
# mlMeshName = 'M128/'
mlMeshName = 'M64/'
REVF = 1

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
#     '5_HornRevRow_refineLocal_realKE_65Dx10Dx4.3D_WdxD_060_WdyzD_15_varScl_eigUinTIPert_6WT'+\
#     '/baseCase/'

vtkFile = case+'project2MLMesh_'+mlMeshName+'VTK/project2MLMesh_'+mlMeshName[:-1]+'_0.vtk'

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

# %% Load Data
shape = (-1,*mlMeshShape,1)
extn = '.pkl'

start = timer()

UHub  = pickle.load(open(mlDataDir+mlMeshName+'UHub'+extn,'rb'))
TIHub = pickle.load(open(mlDataDir+mlMeshName+'TIHub'+extn,'rb'))
defU = pickle.load(open(mlDataDir+mlMeshName+'defU'+extn,'rb')).reshape(shape)
UMag = pickle.load(open(mlDataDir+mlMeshName+'UMag'+extn,'rb')).reshape(shape)
TI  = pickle.load(open(mlDataDir+mlMeshName+'TI'+extn,'rb')).reshape(shape)
tke = pickle.load(open(mlDataDir+mlMeshName+'tke'+extn,'rb')).reshape(shape)
nut = pickle.load(open(mlDataDir+mlMeshName+'nut'+extn,'rb')).reshape(shape)
R = pickle.load(open(mlDataDir+mlMeshName+'R'+extn,'rb')).reshape(*(shape[:-1]),6)
A = pickle.load(open(mlDataDir+mlMeshName+'A'+extn,'rb')).reshape(*(shape[:-1]),6)

print(timer()-start, 's')

# %% Train-Test Split
dataSize = len(UHub)
testSize = dataSize // 5
valSize = testSize // 2

UMagTrain, tkeTrain, nutTrain, RTrain, ATrain, UHubTrain, TIHubTrain = \
    UMag[testSize:], tke[testSize:], nut[testSize:], \
    R[testSize:], A[testSize:], UHub[testSize:], TIHub[testSize:]
UMagVal, tkeVal, nutVal, RVal, AVal, UHubVal, TIHubVal = \
    UMag[:valSize], tke[:valSize], nut[:valSize], \
    R[:valSize], A[:valSize], UHub[:valSize], TIHub[:valSize]      
UMagTest, tkeTest, nutTest, RTest, ATest, UHubTest, TIHubTest = \
    UMag[valSize:testSize], tke[valSize:testSize], nut[valSize:testSize], \
    R[valSize:testSize], A[valSize:testSize], UHub[valSize:testSize], TIHub[valSize:testSize]  
    
print(UMagTrain.shape, UMagVal.shape, UMagTest.shape)

# del UMag, tke, nut, R, A, UHub, TIHub, 

# %% Scale Data
normalize = lambda x, m, s: (x-m)/s

# UMagMax, tkeMax = UMagTrain.max(), tkeTrain.max()
# UMagTrain, tkeTrain = UMagTrain/UMagMax, tkeTrain/tkeMax
# UMagVal, tkeVal = UMagVal/UMagMax, tkeVal/tkeMax
# UMagTest, tkeTest = UMagTest/UMagMax, tkeTest/tkeMax

nutMean, nutStd = nutTrain.mean(axis=0), nutTrain.std(axis=0)
nutTrain = normalize(nutTrain, nutMean, nutStd)
nutVal = normalize(nutVal, nutMean, nutStd)
nutTest = normalize(nutTest, nutMean, nutStd)

# RMean, RStd = RTrain.mean(axis=0), RTrain.std(axis=0)
# RTrain = normalize(RTrain, RMean, RStd)
# RVal = normalize(RVal, RMean, RStd)
# RTest = normalize(RTest, RMean, RStd)

AMean, AStd = ATrain.mean(axis=0), ATrain.std(axis=0)
ATrain = normalize(ATrain, AMean, AStd)
AVal = normalize(AVal, AMean, AStd)
ATest = normalize(ATest, AMean, AStd)

UHubMean, UHubStd = UHubTrain.mean(), UHubTrain.std()
UHubTrain = normalize(UHubTrain, UHubMean, UHubStd)
UHubVal = normalize(UHubVal, UHubMean, UHubStd)
UHubTest = normalize(UHubTest, UHubMean, UHubStd)

TIHubMean, TIHubStd = TIHubTrain.mean(), TIHubTrain.std()
TIHubTrain = normalize(TIHubTrain, TIHubMean, TIHubStd)
TIHubVal = normalize(TIHubVal, TIHubMean, TIHubStd)
TIHubTest = normalize(TIHubTest, TIHubMean, TIHubStd)

# %% Convert to TensorFlow tensors
UHubTrain = tf.convert_to_tensor(UHubTrain, dtype=tdtype)
TIHubTrain = tf.convert_to_tensor(TIHubTrain, dtype=tdtype)
UMagTrain = tf.convert_to_tensor(UMagTrain, dtype=tdtype)
tkeTrain = tf.convert_to_tensor(tkeTrain, dtype=tdtype)
outputTrain = tf.concat((UMagTrain, tkeTrain), axis=-1)

UHubVal = tf.convert_to_tensor(UHubVal, dtype=tdtype)
TIHubVal = tf.convert_to_tensor(TIHubVal, dtype=tdtype)
UMagVal = tf.convert_to_tensor(UMagVal, dtype=tdtype)
tkeVal = tf.convert_to_tensor(tkeVal, dtype=tdtype)
outputVal = tf.concat((UMagVal, tkeVal), axis=-1)

UHubTest = tf.convert_to_tensor(UHubTest, dtype=tdtype)
TIHubTest = tf.convert_to_tensor(TIHubTest, dtype=tdtype)
UMagTest = tf.convert_to_tensor(UMagTest, dtype=tdtype)
tkeTest = tf.convert_to_tensor(tkeTest, dtype=tdtype)
outputTest = tf.concat((UMagTest, tkeTest), axis=-1)

nutTrain = tf.convert_to_tensor(nutTrain, dtype=tdtype)
nutVal = tf.convert_to_tensor(nutVal, dtype=tdtype)
nutTest = tf.convert_to_tensor(nutTest, dtype=tdtype)
    
# RTrain= tf.convert_to_tensor(RTrain, dtype=tdtype)
# RVal= tf.convert_to_tensor(RVal, dtype=tdtype)
# RTest = tf.convert_to_tensor(RTest, dtype=tdtype)
    
ATrain = tf.convert_to_tensor(ATrain, dtype=tdtype)
AVal = tf.convert_to_tensor(AVal, dtype=tdtype)    
ATest = tf.convert_to_tensor(ATest, dtype=tdtype) 
                   
# %% Data collection in a batches
if REVF:
    trainBatch = tf.concat(
        ((tf.zeros_like(nutTrain)+tf.reshape(UHubTrain, [-1,*([1]*(len(nutTrain.shape)-1))])),
         (tf.zeros_like(nutTrain)+tf.reshape(TIHubTrain, [-1,*([1]*(len(nutTrain.shape)-1))])),
         nutTrain
         ), axis=-1
    )
    valBatch = tf.concat(
        ((tf.zeros_like(nutVal)+tf.reshape(UHubVal, [-1,*([1]*(len(nutVal.shape)-1))])),
         (tf.zeros_like(nutVal)+tf.reshape(TIHubVal, [-1,*([1]*(len(nutVal.shape)-1))])),
         nutVal
         ), axis=-1
    )  
    testBatch = tf.concat(
        ((tf.zeros_like(nutTest)+tf.reshape(UHubTest, [-1,*([1]*(len(nutTest.shape)-1))])),
         (tf.zeros_like(nutTest)+tf.reshape(TIHubTest, [-1,*([1]*(len(nutTest.shape)-1))])),
         nutTest
         ), axis=-1
    )       
else:
    trainBatch = tf.concat(
        ((tf.zeros_like(ATrain[:,:,:,:,:1])+tf.reshape(UHubTrain, [-1,*([1]*(len(ATrain.shape)-1))])),
         (tf.zeros_like(ATrain[:,:,:,:,:1])+tf.reshape(TIHubTrain, [-1,*([1]*(len(ATrain.shape)-1))])),
         # RTrain
          ATrain
         ), axis=-1
    )
    valBatch = tf.concat(
        ((tf.zeros_like(AVal[:,:,:,:,:1])+tf.reshape(UHubVal, [-1,*([1]*(len(AVal.shape)-1))])),
         (tf.zeros_like(AVal[:,:,:,:,:1])+tf.reshape(TIHubVal, [-1,*([1]*(len(AVal.shape)-1))])),
         # RVal
          AVal
         ), axis=-1
    )
    testBatch = tf.concat(
        ((tf.zeros_like(ATest[:,:,:,:,:1])+tf.reshape(UHubTest, [-1,*([1]*(len(ATest.shape)-1))])),
         (tf.zeros_like(ATest[:,:,:,:,:1])+tf.reshape(TIHubTest, [-1,*([1]*(len(ATest.shape)-1))])),
         # RTest
          ATest
         ), axis=-1
    )    
       
# %% Model
isUNet, isUNetAug, isAutoEnc_NN, isNN = 1, 0, 0, 0

if isUNet:
    # model = UNet(mlMeshShape, nRandFieldsChannels=2, dropFrac=0.)
    model = UNet(mlMeshShape, nRandFieldsChannels=3, dropFrac=0.)
    # model = UNet(mlMeshShape, nRandFieldsChannels=8, dropFrac=0., l1_lambda=0.)
    modelName = mlDir+'models/UNet_'+mlMeshName[:-1]+'.h5'
elif isUNetAug:
    model = UNetAug(mlMeshShape, dropFrac=0.)
    modelName = mlDir+'models/UNetAug_'+mlMeshName[:-1]+'.h5'
elif isAutoEnc_NN:
    model = UNet(mlMeshShape, nRandFieldsChannels=6, nOutChannels=6, dropFrac=0.)
    # model = UNet(mlMeshShape, nRandFieldsChannels=1, nOutChannels=1, dropFrac=0.)
    modelName = mlDir+'models/UNet_AutoEnc'+mlMeshName[:-1]+'.h5'
elif isNN:
    pass
    
model.summary()

# %% Load and test the model
# dependencies = {'relErrL1': relErrL1, 'relErrL2': relErrL2}
# model = tf.keras.models.load_model(modelName, custom_objects=dependencies)

# %% Trian the model
epochs = 10000

lr = 1e-4
lrS = tf.keras.optimizers.schedules.ExponentialDecay(
    lr, len(UHubTrain)*10, 0.9
)
opt = tf.keras.optimizers.Adam(lr, beta_1=0.9, beta_2=0.999)

# cbs = [callbacks.ModelCheckpoint(modelName, save_best_only=True),
#         callbacks.EarlyStopping(patience=10, monitor='relErrL1')]
cbs = None

model.compile(optimizer=opt, loss='mae', metrics=[relErrL1, relErrL2])

results = model.fit(trainBatch, outputTrain, batch_size=2,
                    validation_data=(valBatch, outputVal),
                    validation_batch_size=2,
                    callbacks=cbs, epochs=epochs)

# results = model.fit(ATrain, ATrain, batch_size=2,
#                     validation_data=(AVal, AVal),
#                     validation_batch_size=2,
#                     callbacks=cbs, epochs=epochs)

# results = model.fit(nutTrain, nutTrain, batch_size=2,
#                     validation_data=(nutVal, nutVal),
#                     validation_batch_size=2,
#                     callbacks=cbs, epochs=epochs)

print(opt._decayed_lr(tf.float32))

# %% Prediction
s = np.random.randint(0,len(testBatch))
print('Sample: ', s)
sampleRegx_pred = r'sample_'+str(s)

U_true = outputTest[s,:,:,:,0].numpy()
k_true = outputTest[s,:,:,:,1].numpy()

U_pred = model.predict(testBatch[s-1:s])[0][:,:,:,0]
k_pred = model.predict(testBatch[s-1:s])[0][:,:,:,1]

print('U_true: ',U_true.reshape(-1))
print('U_pred: ',U_pred.reshape(-1))
print('k_true: ',k_true.reshape(-1))
print('k_pred: ',k_pred.reshape(-1))
















