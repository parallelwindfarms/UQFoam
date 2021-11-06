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

device = utilities.getDefaultDevice()
print(f'Using device: {device}')

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
# mlMeshName, nx, ny, nz = 'M128/', 549, 128, 128 # M128
mlMeshName, nx, ny, nz = 'M64/', 275, 64, 64

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

# %% Load a batch
batchSize = 1000
testSize = batchSize // 10
loadSampleRange = range(0,1000)

UHub  = np.zeros((batchSize))
TIHub = np.zeros((batchSize))
defU = np.zeros((batchSize, mlMeshShape[0]**3))
UMag = np.zeros((batchSize, mlMeshShape[0]**3))
TI  = np.zeros((batchSize, mlMeshShape[0]**3))
tke = np.zeros((batchSize, mlMeshShape[0]**3))
nut = np.zeros((batchSize, mlMeshShape[0]**3))
R = np.zeros((batchSize, mlMeshShape[0]**3, 6))
A = np.zeros((batchSize, mlMeshShape[0]**3, 6))

start = timer()
extn = '.pkl'
for i in tqdm(range(batchSize)):
    s = loadSampleRange[i]
    UHub[i] = pickle.load(open(mlDataDir+mlMeshName+'sample_'+str(s)+'/UHub'+extn,'rb'))
    UMag[i] = pickle.load(open(mlDataDir+mlMeshName+'sample_'+str(s)+'/UMag'+extn,'rb')).reshape(-1)
    defU[i] = pickle.load(open(mlDataDir+mlMeshName+'sample_'+str(s)+'/defU'+extn,'rb')).reshape(-1)
    TIHub[i] = pickle.load(open(mlDataDir+mlMeshName+'sample_'+str(s)+'/TIHub'+extn,'rb'))
    tke[i] = pickle.load(open(mlDataDir+mlMeshName+'sample_'+str(s)+'/k'+extn,'rb')).reshape(-1)
    TI[i] = pickle.load(open(mlDataDir+mlMeshName+'sample_'+str(s)+'/TI'+extn,'rb')).reshape(-1)
    nut[i] = pickle.load(open(mlDataDir+mlMeshName+'sample_'+str(s)+'/nut'+extn,'rb')).reshape(-1)
    R[i] = pickle.load(open(mlDataDir+mlMeshName+'sample_'+str(s)+'/R'+extn,'rb')).reshape(-1,6)
    A[i] = pickle.load(open(mlDataDir+mlMeshName+'sample_'+str(s)+'/A'+extn,'rb')).reshape(-1,6)

print(timer()-start, 's')

# %% Train-Test Split
UMagTrain, tkeTrain, nutTrain, RTrain, ATrain, UHubTrain, TIHubTrain = \
    UMag[testSize:], tke[testSize:], nut[testSize:], \
    R[testSize:], A[testSize:], UHub[testSize:], TIHub[testSize:]
UMagTest, tkeTest, nutTest, RTest, ATest, UHubTest, TIHubTest = \
    UMag[:testSize], tke[:testSize], nut[:testSize], \
    R[:testSize], A[:testSize], UHub[:testSize], TIHub[:testSize]  
    
# %% Normalize
normalize = lambda x, m, s: (x-m)/s

UMagMax, tkeMax = UMagTrain.max(), tkeTrain.max()
UHubMean, UHubStd = UHubTrain.mean(), UHubTrain.std()
TIHubMean, TIHubStd = TIHubTrain.mean(), TIHubTrain.std()
nutMean, nutStd = nutTrain.mean(), nutTrain.std()
RMean, RStd = RTrain.mean(axis=0), RTrain.std(axis=0)
AMean, AStd = ATrain.mean(axis=0), ATrain.std(axis=0)

UMagTrain, tkeTrain = UMagTrain/UMagMax, tkeTrain/tkeMax
nutTrain = normalize(nutTrain, nutMean, nutStd)
RTrain = normalize(RTrain, RMean, RStd)
ATrain = normalize(ATrain, AMean, AStd)
UHubTrain = normalize(UHubTrain, UHubMean, UHubStd)
TIHubTrain = normalize(TIHubTrain, TIHubMean, TIHubStd)

UMagTest, tkeTest = UMagTest/UMagMax, tkeTest/tkeMax
nutTest = normalize(nutTest, nutMean, nutStd)
RTest = normalize(RTest, RMean, RStd)
ATest = normalize(ATest, AMean, AStd)
UHubTest = normalize(UHubTest, UHubMean, UHubStd)
TIHubTest = normalize(TIHubTest, TIHubMean, TIHubStd)

# %% Convert to torch tensors
UMagTrain= torch.tensor(UMagTrain, dtype=torch.float64).to(device)
tkeTrain= torch.tensor(tkeTrain, dtype=torch.float64).to(device)

UMagTest= torch.tensor(UMagTest, dtype=torch.float64).to(device)
tkeTest= torch.tensor(tkeTest, dtype=torch.float64).to(device)
outputTest = torch.concat((UMagTest, tkeTest), dim=-1)

# %% Data collection in a single batch(s)
nInCells, nOutCells = 10000, nCells_WT
inCells = np.random.choice(np.arange(nOutCells),size=nInCells,replace=False)

REVF = 0
if REVF: nRFieldElems, nOutFieldsElems = 1+1+1 , 2 # UHub, TIHub nut # UMag, tke
else: nRFieldElems, nOutFieldsElems = 6+1+1, 2 # UHub, TIHub, A # UMag, tke

Nf = len(UHubTrain) //10 * 9

UinTrainBatch = torch.tensor(UHubTrain[:Nf], requires_grad=True).to(device)
UinValBatch = torch.tensor(UHubTrain[Nf:], requires_grad=True).to(device)
UinTestBatch = torch.tensor(UHubTest, requires_grad=True).to(device)
    
TIinTrainBatch = torch.tensor(TIHubTrain[:Nf], requires_grad=True).to(device)
TIinValBatch = torch.tensor(TIHubTrain[Nf:], requires_grad=True).to(device)
TIinTestBatch = torch.tensor(TIHubTest, requires_grad=True).to(device)

outputTrain = torch.concat((UMagTrain[:Nf], tkeTrain[:Nf]), dim=-1)
outputVal = torch.concat((UMagTrain[Nf:], tkeTrain[Nf:]), dim=-1)

if REVF:
    nutTrainBatch = torch.tensor(nutTrain[:Nf][:,inCells], requires_grad=True).to(device)
    nutValBatch = torch.tensor(nutTrain[Nf:][:,inCells], requires_grad=True).to(device)
    nutTestBatch = torch.tensor(nutTest[:,inCells], requires_grad=True).to(device)
    trainBatch = torch.concat(
        ((torch.zeros_like(nutTrainBatch)+UinTrainBatch.reshape(-1,1)),
         (torch.zeros_like(nutTrainBatch)+TIinTrainBatch.reshape(-1,1)), 
         nutTrainBatch
         ), axis=1
    )
    valBatch = torch.concat(
        ((torch.zeros_like(nutValBatch)+UinValBatch.reshape(-1,1)),
         (torch.zeros_like(nutValBatch)+TIinValBatch.reshape(-1,1)),
         nutValBatch
         ), axis=1
    )  
    testBatch = torch.concat(
        ((torch.zeros_like(nutTestBatch)+UinTestBatch.reshape(-1,1)),
         (torch.zeros_like(nutTestBatch)+TIinTestBatch.reshape(-1,1)),
         nutTestBatch
         ), axis=1
    )       
else:
    # RTrainBatch = torch.tensor(RTrain[:Nf][:,inCells], requires_grad=True).to(device)
    # RValBatch = torch.tensor(RTrain[Nf:][:,inCells], requires_grad=True).to(device)
    ATrainBatch = torch.tensor(ATrain[:Nf][:,inCells], requires_grad=True).to(device)
    AValBatch = torch.tensor(ATrain[Nf:][:,inCells], requires_grad=True).to(device)    
    ATestBatch = torch.tensor(ATest[:,inCells], requires_grad=True).to(device)    
    trainBatch = torch.concat(
        ((torch.zeros_like(ATrainBatch[:,:,0]).unsqueeze(-1)+UinTrainBatch.reshape(-1,1,1)),
         (torch.zeros_like(ATrainBatch[:,:,0]).unsqueeze(-1)+TIinTrainBatch.reshape(-1,1,1)),
         # RTrainBatch
         ATrainBatch
         ), axis=-1
    )
    valBatch = torch.concat(
        ((torch.zeros_like(AValBatch[:,:,0]).unsqueeze(-1)+UinValBatch.reshape(-1,1,1)),
         (torch.zeros_like(AValBatch[:,:,0]).unsqueeze(-1)+TIinValBatch.reshape(-1,1,1)),
         # RValBatch
         AValBatch
         ), axis=-1
    )
    testBatch = torch.concat(
        ((torch.zeros_like(ATestBatch[:,:,0]).unsqueeze(-1)+UinTestBatch.reshape(-1,1,1)),
         (torch.zeros_like(ATestBatch[:,:,0]).unsqueeze(-1)+TIinTestBatch.reshape(-1,1,1)),
         # RValBatch
         ATestBatch
         ), axis=-1
    )    

# %% NN parameters
layersizes = [nInCells*nRFieldElems] + [64]*2 + [nOutCells*nOutFieldsElems]

model = nets.NeuralNetwork(layersizes=layersizes, activation=torch.relu)
model.to(device)
utilities.totalParams(model)

# %% Training loop
epochs = 10000

# opt = torch.optim.Adam(model.parameters(), lr=1e-3)
opt = torch.optim.Adam(model.parameters(), lr=1e-4)

lmbda = 1e-5
bestStatDict = model.state_dict()
bestValLoss  = np.inf
lossFunc = F.mse_loss
torchZero = torch.tensor([0.]).to(device)

startTime = timer()
for epoch in tqdm(range(1,epochs+1)):

    # Training
    opt.zero_grad()
    soln =  model(trainBatch.reshape(-1,nInCells*nRFieldElems))
    trainLoss = lossFunc(soln, outputTrain) + lmbda*utilities.L2Loss(model, device)
    trainLoss.backward()
    opt.step()

    # Validation
    if (epoch)%100 == 0:
        soln =  model(valBatch.reshape(-1,nInCells*nRFieldElems))
        valLoss = lossFunc(soln, outputVal)

        print(
            # f' trainLoss:{trainLoss.item()*1e3:.3f}', \
            # f' valLoss:{valLoss.item()*1e3:.3f}', \
            f' L1Err:{utilities.L1Err(soln, outputVal)*100:.1f}%', \
            f' L2Err:{utilities.L2Err(soln, outputVal)*100:.1f}%', \
        )
        if valLoss.item() < bestValLoss:
            bestStateDict = model.state_dict()
            bestValLoss = valLoss.item()
            
model.load_state_dict(bestStatDict)

endTime = timer()
print(f'Training time = {endTime-startTime} s')

# %% Test
soln = model(testBatch.reshape(-1,nInCells*nRFieldElems))
testLoss = lossFunc(soln, outputTest)
print(
    f'trainLoss:{trainLoss.item()*1e3:.3f} ', \
    f'valLoss:{valLoss.item()*1e3:.3f} ', \
    f'testLoss:{testLoss.item()*1e3:.3f} ', \
    f'L1Err:{utilities.L1Err(soln, outputTest)*100:.1f}% ', \
    f'L2Err:{utilities.L2Err(soln, outputTest)*100:.1f}% \n', \
)
    
for s in np.random.randint(0, len(testBatch), 5):
    print('UHub = ', UHubTest[s]*UHubStd+UHubMean, 'TIHub = ', TIHubTest[s]*TIHubStd+TIHubMean)
    testCaseInp =  testBatch[s].reshape(-1,nInCells*nRFieldElems)
    UMagTestTrue, tkeTestTrue = UMagTest[s]*UMagMax, tkeTest[s]*tkeMax
    UMagTestPred, tkeTestPred = utilities.Net(model, testCaseInp, nOutCells)
    UMagTestPred, tkeTestPred = UMagTestPred[0]*UMagMax, tkeTestPred[0]*tkeMax
    
    print('L1 Error in U =', f'{utilities.L1Err(UMagTestTrue, UMagTestPred)*100:.1f} %')
    print('L2 Error in U =', f'{utilities.L2Err(UMagTestTrue, UMagTestPred)*100:.1f} %')
    print('L1 Error in k =', f'{utilities.L1Err(tkeTestTrue, tkeTestPred)*100:.1f} %')
    print('L2 Error in k =', f'{utilities.L2Err(tkeTestTrue, tkeTestPred)*100:.1f} %')
    
    UMagTestTrue = utilities.torchToNumpy(UMagTestTrue)
    UMagTestPred = utilities.torchToNumpy(UMagTestPred)
    tkeTestTrue = utilities.torchToNumpy(tkeTestTrue)
    tkeTestPred = utilities.torchToNumpy(tkeTestPred)
    UMagDiff = np.abs(UMagTestTrue-UMagTestPred)
    tkeDiff = np.abs(tkeTestTrue-tkeTestPred)
    
    # print('\n', 'True', UMagTestTrue,'\n', 'Pred', UMagTestPred)
    # print('\n', 'True', tkeTestTrue,'\n', 'Pred', tkeTestPred)

# %% Contour yPlane
yPlane = cCenter_WT[y0Plane_WT_idx]/D

fig, ax = plt.subplots(ncols=2, nrows=3, constrained_layout=True,
                        sharey=True, sharex=True, figsize=(10,6))
ax = ax.flat
ms = 15

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
ax[4].set_title('| PINN-OpenFOAM |: UMag'), ax[5].set_title('| DNN-OpenFOAM |: k')

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
ms = 15

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
ax[4].set_title('| PINN-OpenFOAM |: : UMag'), ax[5].set_title('| DNN-OpenFOAM |: k')

fig.colorbar(CS2, ax=ax[0], aspect=50)
fig.colorbar(CS2, ax=ax[2], aspect=50)
fig.colorbar(CS4, ax=ax[4], aspect=50)
fig.colorbar(CS3, ax=ax[1], aspect=50)
fig.colorbar(CS3, ax=ax[3], aspect=50)
fig.colorbar(CS5, ax=ax[5], aspect=50)
