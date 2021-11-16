#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# %% Import packages, classes and funcitons
import os, sys

import numpy as np
import torch
import torch.nn.functional as F
import matplotlib.pyplot as plt

from timeit import default_timer as timer
from tqdm import tqdm
import pickle

torch.set_default_dtype(torch.float64)

cwd = os.getcwd() + "/"
os.environ['CUDA_LAUNCH_BLOCKING'] = "1"

SCRIPTS = os.environ['SCRIPTS']
DATA    = os.environ['windTurbineData']

sys.path.append(SCRIPTS+'/10windTurbine')
sys.path.append(SCRIPTS+'/pythonScripts')
sys.path.append(SCRIPTS+'/10windTurbine'+'/PyTorchLib')

import myUQlib
import nets, utilities
from ML_GetFuncs import getMeshData, getCaseData

device = utilities.getDefaultDevice(3)
print(f'Using device: {device}')

randomSeed = 42
np.random.seed = 42

# %% Hyper-parameters
tdtype=torch.float64
# mlMeshName = 'M128/'
mlMeshName = 'M64/'
N_IN_CELLS = 1000
N_LAYERS, N_NEURONS = 5, 256
HIDDEN_LAYERS = N_LAYERS*[N_NEURONS]
REVF = 0

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

# %% Load Data and Data Statistics
shape = (-1,mlMeshShape[0]**3)
extn = '.pkl'

start = timer()

UHub  = pickle.load(open(mlDataDir+mlMeshName+'UHub'+extn,'rb'))
TIHub = pickle.load(open(mlDataDir+mlMeshName+'TIHub'+extn,'rb'))
defU = pickle.load(open(mlDataDir+mlMeshName+'defU'+extn,'rb')).reshape(shape)
UMag = pickle.load(open(mlDataDir+mlMeshName+'UMag'+extn,'rb')).reshape(shape)
TI  = pickle.load(open(mlDataDir+mlMeshName+'TI'+extn,'rb')).reshape(shape)
tke = pickle.load(open(mlDataDir+mlMeshName+'k'+extn,'rb')).reshape(shape)
nut = pickle.load(open(mlDataDir+mlMeshName+'nut'+extn,'rb')).reshape(shape)
R = pickle.load(open(mlDataDir+mlMeshName+'R'+extn,'rb')).reshape(*shape,6)
A = pickle.load(open(mlDataDir+mlMeshName+'A'+extn,'rb')).reshape(*shape,6)

print(timer()-start, 's')

# %% Train-Test Split
dataSize = len(UHub)
trainSize = int(dataSize * 0.8)
valSize = trainSize + (dataSize-trainSize) // 2

UMagTrain, tkeTrain, TITrain, nutTrain, RTrain, ATrain, UHubTrain, TIHubTrain = \
    UMag[:trainSize], tke[:trainSize], TI[:trainSize], nut[:trainSize], \
    R[:trainSize], A[:trainSize], UHub[:trainSize], TIHub[:trainSize]
    
UMagVal, tkeVal, TIVal, nutVal, RVal, AVal, UHubVal, TIHubVal = \
    UMag[trainSize:valSize], tke[trainSize:valSize], TI[trainSize:valSize], nut[trainSize:valSize], \
    R[trainSize:valSize], A[trainSize:valSize], UHub[trainSize:valSize], TIHub[trainSize:valSize]     
    
UMagTest, tkeTest, TITest, nutTest, RTest, ATest, UHubTest, TIHubTest = \
    UMag[valSize:], tke[valSize:], TI[valSize:], nut[valSize:], \
    R[valSize:], A[valSize:], UHub[valSize:], TIHub[valSize:]  
    
print(UMagTrain.shape, UMagVal.shape, UMagTest.shape)

# del UMag, tke, nut, R, A, UHub, TIHub, 

# %% Scale Data
normalize = lambda x, m, s: (x-m)/s

UMagMax, tkeMax = UMagTrain.max(), tkeTrain.max()
UMagTrain, tkeTrain = UMagTrain/UMagMax, tkeTrain/tkeMax
UMagVal, tkeVal = UMagVal/UMagMax, tkeVal/tkeMax
UMagTest, tkeTest = UMagTest/UMagMax, tkeTest/tkeMax

# nutMean, nutStd = nutTrain.mean(axis=0), nutTrain.std(axis=0)
# nutTrain = normalize(nutTrain, nutMean, nutStd)
# nutVal = normalize(nutVal, nutMean, nutStd)
# nutTest = normalize(nutTest, nutMean, nutStd)

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

# %% Convert to torch tensors
UHubTrain = torch.tensor(UHubTrain, requires_grad=True)
TIHubTrain = torch.tensor(TIHubTrain, requires_grad=True)
UMagTrain = torch.tensor(UMagTrain, dtype=tdtype)
tkeTrain = torch.tensor(tkeTrain, dtype=tdtype)
outputTrain = torch.concat((UMagTrain, tkeTrain), dim=-1).to(device)

UHubVal = torch.tensor(UHubVal, requires_grad=True)
TIHubVal = torch.tensor(TIHubVal, requires_grad=True)
UMagVal = torch.tensor(UMagVal, dtype=tdtype)
tkeVal = torch.tensor(tkeVal, dtype=tdtype)
outputVal = torch.concat((UMagVal, tkeVal), dim=-1).to(device)

UHubTest = torch.tensor(UHubTest, requires_grad=True)
TIHubTest = torch.tensor(TIHubTest, requires_grad=True)
UMagTest = torch.tensor(UMagTest, dtype=tdtype)
tkeTest = torch.tensor(tkeTest, dtype=tdtype)
outputTest = torch.concat((UMagTest, tkeTest), dim=-1).to(device)
                          
# %% Data collection in batches
nInCells, nOutCells = N_IN_CELLS, nCells_WT
inCells = np.random.choice(np.arange(nOutCells),size=nInCells,replace=False)

nutTrain = torch.tensor(nutTrain[:,inCells], requires_grad=True)
nutVal = torch.tensor(nutVal[:,inCells], requires_grad=True)
nutTest = torch.tensor(nutTest[:,inCells], requires_grad=True)

# RTrain= torch.tensor(RTrain[:,inCells], requires_grad=True)
# RVal= torch.tensor(RVal[:,inCells], requires_grad=True)
# RTest = torch.tensor(RTest[:,inCells], requires_grad=True)
ATrain = torch.tensor(ATrain[:,inCells], requires_grad=True)
AVal = torch.tensor(AVal[:,inCells], requires_grad=True)    
ATest = torch.tensor(ATest[:,inCells], requires_grad=True) 
    
# UHub, TIHub, nut # UMag, tke
if REVF: nRFieldElems, nOutFieldsElems = 1+1+1 , 2 
# UHub, TIHub, A # UMag, tke
else: nRFieldElems, nOutFieldsElems = 6+1+1, 2 
# else: nRFieldElems, nOutFieldsElems = 1+1, 2 

if REVF:
    trainBatch = torch.concat(
        ((torch.zeros_like(nutTrain)+UHubTrain.reshape(-1,1)),
         (torch.zeros_like(nutTrain)+TIHubTrain.reshape(-1,1)), 
         nutTrain
         ), axis=1
    )
    valBatch = torch.concat(
        ((torch.zeros_like(nutVal)+UHubVal.reshape(-1,1)),
         (torch.zeros_like(nutVal)+TIHubVal.reshape(-1,1)),
         nutVal
         ), axis=1
    )  
    testBatch = torch.concat(
        ((torch.zeros_like(nutTest)+UHubTest.reshape(-1,1)),
         (torch.zeros_like(nutTest)+TIHubTest.reshape(-1,1)),
         nutTest
         ), axis=1
    )       
else:
    trainBatch = torch.concat(
        ((torch.zeros_like(ATrain[:,:,:1])+UHubTrain.reshape(-1,1,1)),
         (torch.zeros_like(ATrain[:,:,:1])+TIHubTrain.reshape(-1,1,1)),
         # RTrain
          ATrain
         ), axis=-1
    )
    valBatch = torch.concat(
        ((torch.zeros_like(AVal[:,:,:1])+UHubVal.reshape(-1,1,1)),
         (torch.zeros_like(AVal[:,:,:1])+TIHubVal.reshape(-1,1,1)),
         # RVal
          AVal
         ), axis=-1
    )
    testBatch = torch.concat(
        ((torch.zeros_like(ATest[:,:,:1])+UHubTest.reshape(-1,1,1)),
         (torch.zeros_like(ATest[:,:,:1])+TIHubTest.reshape(-1,1,1)),
         # RTest
          ATest
         ), axis=-1
    )    

# trainBatch = torch.concat((UHubTrain.reshape(-1,1), TIHubTrain.reshape(-1,1)), dim=-1)
# valBatch = torch.concat((UHubVal.reshape(-1,1), TIHubVal.reshape(-1,1)), dim=-1)
# testBatch = torch.concat((UHubTest.reshape(-1,1), TIHubTest.reshape(-1,1)), dim=-1)

# %% NN Model
layersizes = [nInCells*nRFieldElems] + HIDDEN_LAYERS + [nOutCells*nOutFieldsElems]

model = nets.NeuralNetwork(layersizes=layersizes, activation=torch.relu)
model.to(device)
utilities.totalParams(model)

# %% Training loop
epochs = 10000

opt = torch.optim.Adam(
    model.parameters(),
    # lr=1e-3
    lr=1e-4
    # lr=1e-5
)

lmbda = 1e-3 #1e-5
bestStatDict = model.state_dict()
bestValLoss  = np.inf
bestL1Err = np.inf
lossFunc = F.mse_loss

startTime = timer()
for epoch in tqdm(range(1,epochs+1)):

    # Training
    opt.zero_grad()
    soln =  model(trainBatch.reshape(-1,nInCells*nRFieldElems).to(device))
    trainLoss = lossFunc(soln, outputTrain) + \
        lmbda*utilities.L2Loss(model, device)
    trainLoss.backward()
    opt.step()

    L1ErrTrain = utilities.L1Err(soln, outputTrain)*100
    L2ErrTrain = utilities.L2Err(soln, outputTrain)*100

    # Validation
    if (epoch)%100 == 0:
        with torch.no_grad():
            soln =  model(valBatch.reshape(-1,nInCells*nRFieldElems).to(device))
            valLoss = lossFunc(soln, outputVal.to(device))
            L1ErrVal = utilities.L1Err(soln, outputVal)*100
            L2ErrVal = utilities.L2Err(soln, outputVal)*100

        print(
            # f' trainLoss:{trainLoss.item()*1e3:.3f}', \
            # f' valLoss:{valLoss.item()*1e3:.3f}', \
            f' L1ErrTrain:{L1ErrTrain:.1f}%', f' L2ErrTrain:{L2ErrTrain:.1f}%', \
            f' L1ErrVal:{L1ErrVal:.1f}%', f' L2ErrVal:{L2ErrVal:.1f}%', \
        )
        # if valLoss.item() < bestValLoss:
        #     bestStateDict = model.state_dict()
        #     bestValLoss = valLoss.item()
        if L1ErrVal < bestL1Err:
            bestStateDict = model.state_dict()
            bestL1Err = L1ErrVal
            
model.load_state_dict(bestStatDict)

endTime = timer()
print(f'Training time = {endTime-startTime} s')

# %% Test
soln = model(testBatch.reshape(-1,nInCells*nRFieldElems).to(device))
testLoss = lossFunc(soln, outputTest)
print(
    f'trainLoss:{trainLoss.item()*1e3:.3f} ', \
    f'valLoss:{valLoss.item()*1e3:.3f} ', \
    f'testLoss:{testLoss.item()*1e3:.3f} ', \
    f'L1Err:{utilities.L1Err(soln, outputTest)*100:.1f}% ', \
    f'L2Err:{utilities.L2Err(soln, outputTest)*100:.1f}% \n', \
)
    
# for s in [97,6,53]:
for s in np.random.randint(0, len(testBatch), 5):
    print('##############################################\n')
    print('Smaple #', s, '\n UHub =', \
          (UHubTest[s]*UHubStd+UHubMean).item()* 100 // 10 / 10, 'TIHub =', \
          (TIHubTest[s]*TIHubStd+TIHubMean).item()* 100 // 10 / 10
    )
    testCaseInp =  testBatch[s].reshape(-1,nInCells*nRFieldElems).to(device)
    UMagTestTrue, tkeTestTrue = UMagTest[s]*UMagMax, tkeTest[s]*tkeMax
    UMagTestPred, tkeTestPred = utilities.Net(model, testCaseInp, nOutCells)
    UMagTestPred, tkeTestPred = UMagTestPred[0]*UMagMax, tkeTestPred[0]*tkeMax

    UMagTestTrue = utilities.torchToNumpy(UMagTestTrue)
    UMagTestPred = utilities.torchToNumpy(UMagTestPred)
    tkeTestTrue = utilities.torchToNumpy(tkeTestTrue)
    tkeTestPred = utilities.torchToNumpy(tkeTestPred)
    UMagDiff = np.abs(UMagTestTrue-UMagTestPred)
    tkeDiff = np.abs(tkeTestTrue-tkeTestPred)
    
    print(' L1 Error in U =', f'{utilities.L1Err(UMagTestTrue, UMagTestPred)*100:.1f} %')
    print(' L2 Error in U =', f'{utilities.L2Err(UMagTestTrue, UMagTestPred)*100:.1f} %')
    print(' L1 Error in k =', f'{utilities.L1Err(tkeTestTrue, tkeTestPred)*100:.1f} %')
    print(' L2 Error in k =', f'{utilities.L2Err(tkeTestTrue, tkeTestPred)*100:.1f} %')
    
    print('\n', 'True', UMagTestTrue,'\n', 'Pred', UMagTestPred)
    print('\n', 'True', tkeTestTrue,'\n', 'Pred', tkeTestPred, '\n')

# %% Contour yPlane
yPlane = cCenter_WT[y0Plane_WT_idx]/D

fig, ax = plt.subplots(ncols=2, nrows=3, constrained_layout=True,
                        sharey=True, sharex=True, figsize=(10,6))
ax = ax.flat
ms = 25

norm = plt.Normalize(UMagTestTrue[y0Plane_WT_idx].min(), UMagTestTrue[y0Plane_WT_idx].max())
CS0 = ax[0].scatter(yPlane[:,0], yPlane[:,2], c=UMagTestPred[y0Plane_WT_idx], s=ms, norm=norm)
CS2 = ax[2].scatter(yPlane[:,0], yPlane[:,2], c=UMagTestTrue[y0Plane_WT_idx], s=ms, norm=norm)

norm = plt.Normalize(0, UMagTestTrue[y0Plane_WT_idx].max()*0.1)
CS4 = ax[4].scatter(yPlane[:,0], yPlane[:,2], c=UMagDiff[y0Plane_WT_idx], s=ms, norm=norm)

norm = plt.Normalize(tkeTestTrue[y0Plane_WT_idx].min(), tkeTestTrue[y0Plane_WT_idx].max())
CS1 = ax[1].scatter(yPlane[:,0], yPlane[:,2], c=tkeTestPred[y0Plane_WT_idx], s=ms, norm=norm)
CS3 = ax[3].scatter(yPlane[:,0], yPlane[:,2], c=tkeTestTrue[y0Plane_WT_idx], s=ms, norm=norm)

norm = plt.Normalize(0, tkeTestTrue[y0Plane_WT_idx].max()*0.1)
CS5 = ax[5].scatter(yPlane[:,0], yPlane[:,2], c=tkeDiff[y0Plane_WT_idx], s=ms, norm=norm)

ax[0].set_title('DNN: UMag'), ax[1].set_title('DNN: k')
ax[2].set_title('OpenFOAM: UMag'), ax[3].set_title('OpenFOAM: k')
ax[4].set_title('|DNN-OpenFOAM|: UMag'), ax[5].set_title('|DNN-OpenFOAM|: k')

fig.colorbar(CS2, ax=ax[0], aspect=50)
fig.colorbar(CS2, ax=ax[2], aspect=50)
fig.colorbar(CS4, ax=ax[4], aspect=50)
fig.colorbar(CS3, ax=ax[1], aspect=50)
fig.colorbar(CS3, ax=ax[3], aspect=50)
fig.colorbar(CS5, ax=ax[5], aspect=50)

## %% Contour zPlane

zPlane = cCenter_WT[zhPlane_WT_idx]/D

fig, ax = plt.subplots(ncols=2, nrows=3, constrained_layout=True,
                        sharey=True, sharex=True, figsize=(10,6))
ax = ax.flat

norm = plt.Normalize(UMagTestTrue[zhPlane_WT_idx].min(), UMagTestTrue[zhPlane_WT_idx].max())
CS0 = ax[0].scatter(zPlane[:,0], zPlane[:,1], c=UMagTestPred[zhPlane_WT_idx], s=ms, norm=norm)
CS2 = ax[2].scatter(zPlane[:,0], zPlane[:,1], c=UMagTestTrue[zhPlane_WT_idx], s=ms, norm=norm)

norm = plt.Normalize(0, UMagTestTrue[zhPlane_WT_idx].max()*0.1)
CS4 = ax[4].scatter(zPlane[:,0], zPlane[:,1], c=UMagDiff[zhPlane_WT_idx], s=ms, norm=norm)

norm = plt.Normalize(tkeTestTrue[zhPlane_WT_idx].min(), tkeTestTrue[zhPlane_WT_idx].max())
CS1 = ax[1].scatter(zPlane[:,0], zPlane[:,1], c=tkeTestPred[zhPlane_WT_idx], s=ms, norm=norm)
CS3 = ax[3].scatter(zPlane[:,0], zPlane[:,1], c=tkeTestTrue[zhPlane_WT_idx], s=ms, norm=norm)

norm = plt.Normalize(0, tkeTestTrue[zhPlane_WT_idx].max()*0.1)
CS5 = ax[5].scatter(zPlane[:,0], zPlane[:,1], c=tkeDiff[zhPlane_WT_idx], s=ms, norm=norm)

ax[0].set_title('DNN: UMag'), ax[1].set_title('DNN: k')
ax[2].set_title('OpenFOAM: UMag'), ax[3].set_title('OpenFOAM: k')
ax[4].set_title('|DNN-OpenFOAM|: UMag'), ax[5].set_title('|DNN-OpenFOAM|: k')

fig.colorbar(CS2, ax=ax[0], aspect=50)
fig.colorbar(CS2, ax=ax[2], aspect=50)
fig.colorbar(CS4, ax=ax[4], aspect=50)
fig.colorbar(CS3, ax=ax[1], aspect=50)
fig.colorbar(CS3, ax=ax[3], aspect=50)
fig.colorbar(CS5, ax=ax[5], aspect=50)
