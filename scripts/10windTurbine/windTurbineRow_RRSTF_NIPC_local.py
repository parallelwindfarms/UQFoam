#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 29 15:21:21 2021

@author: jigar
"""

# %% Import Required Modules and env vars
from __future__ import print_function

import os, sys
import numpy as np
import chaospy as cp
import pyvista as pv
import matplotlib.pyplot as plt

from multiprocessing import Pool
from functools import partial
from timeit import default_timer as timer

import pandas as pd

SCRIPTS = os.environ['SCRIPTS']
DATA    = os.environ['windTurbineData']

sys.path.append(SCRIPTS+'/10windTurbine')
sys.path.append(SCRIPTS+'/pythonScripts')

cwd = os.getcwd() + "/"

import myUQlib

# %% Case details
caseName   = cwd.split('/')[-2]
casePngDir = '/2RRSTF/NIPC/uqPaperCases/' + caseName+'/png'

if (os.path.exists(DATA+casePngDir))==False:
    print('Making new case data directory...')
    os.makedirs(DATA+casePngDir, exist_ok=True)

# %% Parametrs
D = 80
h = 70
Uref = 8

# mesh dependent
nx = int(40*0.5*60)
Wdx = D/60
d_D = 7
nWT = 6#8#10

# center of AD
ADloc = (np.array([i*d_D for i in range(nWT)])*D,0,h)

# %% Setting the PCE parameters
nSamplesOrg = 30
d = 5

# delB = cp.Iid(cp.Uniform(0,1), d)    # w/ Uniform Sampling
# delB = cp.Iid(cp.Gamma(2,0.15,0), d) # w/ Gamma Sampling
delB = cp.Iid(cp.Normal(0,1.25), d)  # w/ Normal Sampling

np.random.seed(1)
delBSamples = delB.sample(nSamplesOrg, rule='L')
delBSamples.sort()
delBSamples = np.hstack((delBSamples, delBSamples, delBSamples))

n = 1 ; phi = cp.orth_ttr(n, dist=delB, normed=True)
Pplus1 = len(phi)

nSamples = range(nSamplesOrg*3)

# %%  ###############################################################################################
# Plot over centerline
#############################################################################################################

# %% Read OpenFOAM case mesh access to grid and cell values
vtkFile   = 'baseCase/VTK/baseCase_10000.vtk'

mesh      = pv.UnstructuredGrid(cwd + vtkFile)
nCells    = mesh.n_cells
cCenter   = np.array(mesh.cell_centers().points)

shiftedcCenter = cCenter[:,1:3] - ADloc[1:3]
cIdxInDisk     = np.where(np.linalg.norm(shiftedcCenter, axis=1) <= D/2)[0]
cCenterInDisk  = cCenter[cIdxInDisk]

x_D = cCenter[0:nx,0]/D

nCellsInDisk = int(cIdxInDisk.shape[0]/nx)

UDetMag = np.linalg.norm(mesh.cell_data['U'][cIdxInDisk], axis=1)
tkeDet  = mesh.cell_data['k'][cIdxInDisk]
TIDet   = np.sqrt(tkeDet*2/3)/Uref *100

defUDet     = (Uref-UDetMag)/Uref
defUDetAvgd = np.average(np.reshape(defUDet, (nCellsInDisk, nx)), axis=0)
# TIDetAvgd   = np.average(np.reshape(TIDet-TIDet[0], (nCellsInDisk, nx)), axis=0)
TIDetAvgd   = np.average(np.reshape(TIDet, (nCellsInDisk, nx)), axis=0)

# %% Importing LES, RANS, UQRANS data
defU_LES = np.loadtxt(DATA+"/0LES/HornRevWTs/defU.txt", delimiter=',')
TI_LES = np.loadtxt(DATA+"/0LES/HornRevWTs/TI.txt", delimiter=',')

def getSamplesOverDiskVol(cIdxInDisk, s):
    print(s)
    vtkFile = 'sample_'+str(s)+'/VTK/'+'sample_'+str(s)+'_10000.vtk'
    mesh    = pv.UnstructuredGrid(cwd + vtkFile)
    UMag    = np.linalg.norm(mesh.cell_data['U'][cIdxInDisk], axis=1)
    tke     = mesh.cell_data['k'][cIdxInDisk]
    return UMag, tke

pool = Pool(processes=8)
start = timer()
data = pool.map(partial(getSamplesOverDiskVol, cIdxInDisk), (nSamples))
USamplesMag = np.array([data[s][0] for s in range(len(nSamples))])
tkeSamples = np.array([data[s][1] for s in range(len(nSamples))])
pool.close()
pool.join()
print(timer()-start, 's')

try:
    UrefSamples = np.array([float(pd.read_csv('Samples_Uref/Uref'+str(i)).\
                  columns[0].strip(';').split()[1]) for i in range(len(nSamples))])
    UrefSamples = UrefSamples.reshape(-1,1)
except FileNotFoundError:
    UrefSamples = Uref

# %% PCE of U
defUSamples = (UrefSamples-USamplesMag)/UrefSamples
defU_PCEApp = cp.fit_regression(phi, delBSamples, defUSamples)
defU_PCEModes = (defU_PCEApp.coefficients).T

defU0PCEAvgd     = np.average(np.reshape(defU_PCEModes[0], (nCellsInDisk, nx)), axis=0)
defUSigmaPCEAvgd = np.average(np.reshape(np.sqrt(np.sum(defU_PCEModes[1:]**2, axis=0)), (nCellsInDisk, nx)), axis=0)

## %% PCE of TI
TISamples = np.sqrt(tkeSamples*2/3)/UrefSamples *100
TISamples = TISamples #- TISamples[:,0].reshape(-1,1)
TI_PCEApp = cp.fit_regression(phi, delBSamples, TISamples)
TI_PCEModes = (TI_PCEApp.coefficients).T

TI0PCEAvgd     = np.average(np.reshape(TI_PCEModes[0], (nCellsInDisk, nx)), axis=0)
TISigmaPCEAvgd = np.average(np.reshape(np.sqrt(np.sum(TI_PCEModes[1:]**2, axis=0)), (nCellsInDisk, nx)), axis=0)

# %% Computing averaged U fields
defUMean  = np.mean((UrefSamples-USamplesMag)/UrefSamples, axis=0)
defUSigma = np.std((UrefSamples-USamplesMag)/UrefSamples, axis=0)
defUMeanAvgd  = np.average(np.reshape(defUMean, (nCellsInDisk, nx)), axis=0)
defUSigmaAvgd = np.average(np.reshape(defUSigma, (nCellsInDisk, nx)), axis=0)

## %% Computing averaged TI fields
TISamples = np.sqrt(tkeSamples*2/3)/UrefSamples *100
TISamples = TISamples #- TISamples[:,0].reshape(-1,1)
TIMean    = np.mean(TISamples, axis=0)
TISigma   = np.std(TISamples, axis=0)

TIMeanAvgd  = np.average(np.reshape(TIMean, (nCellsInDisk, nx)), axis=0)
TISigmaAvgd = np.average(np.reshape(TISigma, (nCellsInDisk, nx)), axis=0)

TIMin = 5.5 * 0

# %% Figures settings
myUQlib.rcParamsSettings(15)

DETClr  = (0.6350, 0.0780, 0.1840)
meanClr = 'b'
fillClr = 'b'
LESClr  = 'k'
diskClr = 'k'

diskSpan = (ADloc[0]/D-Wdx/D/2, ADloc[0]/D+Wdx/D/2)

# Plotting averaged fields

MC = False
MC = True

overlap = True
overClr=DETClr
overAlpha=0.25
scale=1.025
lw_1 = 1 if overlap else 0
        
if MC: N=2 
else: N=1 

fig, ax = plt.subplots(ncols=1, nrows=2, constrained_layout=True, sharex=True, 
                       figsize=(10,6))

axIdx = 0
ax[axIdx].plot(defU_LES[:-1,0], defU_LES[:-1,1], color=LESClr, marker="o", 
               mfc='none', lw=0, ms=3)
if overlap==False:
    ax[axIdx].plot(x_D, defUDetAvgd, color=DETClr)

if MC==True:
    ax[axIdx].plot(x_D, defUMeanAvgd, color=meanClr)
    ax[axIdx].fill_between(x_D, defUMeanAvgd + N*defUSigmaAvgd, \
                                defUMeanAvgd - N*defUSigmaAvgd, \
                                alpha=0.2, linewidth=lw_1, color='b')
    if overlap==True:
        ax[axIdx].plot(x_D, defU0PCEAvgd *scale, color=overClr)
        ax[axIdx].fill_between(x_D, defU0PCEAvgd *scale + N/3*defUSigmaPCEAvgd, \
                                    defU0PCEAvgd *scale - N/3*defUSigmaPCEAvgd, \
                                    alpha=overAlpha, linewidth=lw_1, color=overClr)
            
else:
    ax[axIdx].plot(x_D, defU0PCEAvgd, color=meanClr)
    ax[axIdx].fill_between(x_D, defU0PCEAvgd + N*defUSigmaPCEAvgd, \
                                defU0PCEAvgd - N*defUSigmaPCEAvgd, \
                                alpha=0.2, linewidth=0, color='b')

# myUQlib.plotDiskSpans(ax[axIdx], ADloc[0]/D, Wdx/D)
ax[axIdx].set_ylim(bottom=0.0,top=0.6)
ax[axIdx].set_yticks([0.0, 0.2, 0.4, 0.6])
ax[axIdx].set_ylabel('$\Delta u / u_h$')
ax[axIdx].spines['right'].set_visible(False)
ax[axIdx].spines['top'].set_visible(False)

axIdx = 1
ax[axIdx].plot(TI_LES[:-1,0], TI_LES[:-1,1],#-TI_LES[:-1,1].min()+TIMin,
               color=LESClr, marker="o", mfc='none', lw=0, ms=3)
if overlap==False:
    ax[axIdx].plot(x_D, TIDetAvgd + TIMin, color=DETClr)

if MC==True:
    ax[axIdx].plot(x_D, TIMeanAvgd + TIMin, color=meanClr)
    ax[axIdx].fill_between(x_D, TIMeanAvgd + N*TISigmaAvgd + TIMin,
                                TIMeanAvgd - N*TISigmaAvgd + TIMin, \
                                alpha=0.2, linewidth=lw_1, color='b')
        
    if overlap==True:
        ax[axIdx].plot(x_D, TI0PCEAvgd *scale + TIMin, color=overClr)
        ax[axIdx].fill_between(x_D, TI0PCEAvgd *scale + N/3*TISigmaPCEAvgd + TIMin, \
                                    TI0PCEAvgd *scale - N/3*TISigmaPCEAvgd + TIMin, \
                                    alpha=overAlpha, linewidth=lw_1, color=overClr)
            
else:
    ax[axIdx].plot(x_D, TI0PCEAvgd + TIMin, color=meanClr)
    ax[axIdx].fill_between(x_D, TI0PCEAvgd + N*TISigmaPCEAvgd + TIMin,
                                TI0PCEAvgd - N*TISigmaPCEAvgd + TIMin, \
                                alpha=0.2, linewidth=0, color='b')

# myUQlib.plotDiskSpans(ax[axIdx], ADloc[0]/D, Wdx/D)
ax[axIdx].set_ylim(bottom=2.5,top=17.5)
ax[axIdx].set_yticks([5, 10, 15])
ax[axIdx].set_xlim(left=-4,right=45)
ax[axIdx].set_xticks(list(ADloc[0]/D)+[45.0])
ax[axIdx].set_xlabel('$x/D$')
ax[axIdx].set_ylabel('$I [\\%]$')
ax[axIdx].spines['right'].set_visible(False)
ax[axIdx].spines['top'].set_visible(False)

if overlap==False:
    ax[0].legend(['LES', 'DET', r'$\textbf{E}[\bullet]$',
                  r'$\textbf{E}[\bullet] \, \pm \, 2\sqrt{\textbf{V}[\bullet]}$'],
                 ncol=4, frameon=False, columnspacing=0.75,
                 loc='upper center', bbox_to_anchor=(0.5, 1.35))

ax[0].grid(axis='x')
ax[1].grid(axis='x')

if MC and overlap:
    fig.savefig(DATA+casePngDir+'/diskAvgdFields_MC_PCE.png', dpi=150)    
elif MC and overlap==False:
    fig.savefig(DATA+casePngDir+'/diskAvgdFields_MC.png', dpi=150)
else:
    fig.savefig(DATA+casePngDir+'/diskAvgdFields_PCE.png', dpi=150)
    
# %%  ###############################################################################################
# Plot over lines at different x locations
#############################################################################################################

# %% Importing LES, RANS, UQRANS data
lines = (1,2,3,4,5,6)

# % Load LES data ###############
defULES, TILES = [], []

for l in lines:
    defULES.append(np.loadtxt(DATA+"/0LES/HornRevWTs/defU_x_D_5_WT"+str(l)+'.txt', delimiter=','))
    TILES.append(np.loadtxt(DATA+"/0LES/HornRevWTs/TI_x_D_5_WT"+str(l)+'.txt', delimiter=','))

# % Load DET data ###############
numLines = len(lines)
yPts = 100
nCells = numLines*yPts

baseCase = 'baseCase/'
sampleDir = 'postProcessing/sample/10000/'

yByD    = np.loadtxt(baseCase + sampleDir + "AD1_X_by_D_5_U.xy")[:,0] / D
UDetMag = np.linalg.norm(np.loadtxt(baseCase + sampleDir + "AD1_X_by_D_5_U.xy")[:,1:], axis=1, keepdims=True)
tkeDet  = (np.loadtxt(baseCase + sampleDir + "AD1_X_by_D_5_k_nut_p.xy")[:,1]).reshape((yPts,1))

for l in lines[1:]:
    fileName = baseCase + sampleDir+'AD'+str(l)+'_X_by_D_5'
    UDetMag  = np.hstack( (UDetMag, np.linalg.norm(np.loadtxt(fileName+'_U.xy')[:,1:], axis=1, keepdims=True)) )
    tkeDet   = np.hstack( (tkeDet, (np.loadtxt(fileName+'_k_nut_p.xy')[:,1]).reshape((yPts,1))) )

# % Load UQRANS data ###############
def getSamplesOverLines(sampleDir, lines, yPts, s):
    print(s)
    baseCase = 'sample_'+str(s)+'/'
    UMag = np.linalg.norm(np.loadtxt(baseCase + sampleDir + "AD1_X_by_D_5_U.xy")[:,1:], axis=1, keepdims=True)
    tke  = (np.loadtxt(baseCase + sampleDir + "AD1_X_by_D_5_k_nut_p.xy")[:,1]).reshape((yPts,1))

    for l in lines[1:]:
        fileName = baseCase + sampleDir+'AD'+str(l)+'_X_by_D_5'
        UMag = np.hstack( (UMag, np.linalg.norm(np.loadtxt(fileName+'_U.xy')[:,1:], axis=1, keepdims=True)) )
        tke = np.hstack( (tke, (np.loadtxt(fileName+'_k_nut_p.xy')[:,1]).reshape((yPts,1))) )

    return UMag, tke

pool = Pool(processes=8)
start = timer()
data = pool.map(partial(getSamplesOverLines, sampleDir, lines, yPts), (nSamples))
USamplesMag = np.array([data[s][0] for s in range(len(nSamples))])
tkeSamples = np.array([data[s][1] for s in range(len(nSamples))])
pool.close()
pool.join()
print(timer()-start, 's')

# %% PCE of U
UrefSamples = UrefSamples.reshape(-1,1,1)

defUSamples = (UrefSamples-USamplesMag)/UrefSamples

defU_PCEModes = np.zeros((Pplus1, yPts, numLines))
for l in range(numLines):
    defU_PCEApp = cp.fit_regression(phi, delBSamples, defUSamples[:,:,l])
    defU_PCEModes[:,:,l] = (defU_PCEApp.coefficients).T

defU0PCE     = defU_PCEModes[0]
defUSigmaPCE = np.sqrt(np.sum(defU_PCEModes[1:]**2, axis=0))

# % PCE of TI
TISamples = np.sqrt(tkeSamples*2/3)/UrefSamples *100
TISamples = TISamples #- TISamples.min(axis=1, keepdims=True)

TI_PCEModes = np.zeros((Pplus1, yPts, numLines))
for l in range(numLines):
    TI_PCEApp = cp.fit_regression(phi, delBSamples, TISamples[:,:,l])
    TI_PCEModes[:,:,l] = (TI_PCEApp.coefficients).T
TI_PCEModes = np.array(TI_PCEModes)

TI0PCE     = TI_PCEModes[0]
TISigmaPCE = np.sqrt(np.sum(TI_PCEModes[1:]**2, axis=0))

# %% Computing averaged U fields
defUDet   = (Uref-UDetMag)/Uref
UrefSamples = UrefSamples.reshape(-1,1,1)
defUMean  = np.mean((UrefSamples-USamplesMag)/UrefSamples, axis=0)
defUSigma = np.std((UrefSamples-USamplesMag)/UrefSamples, axis=0)

# % Computing averaged TI fields
TIDet     = np.sqrt(tkeDet*2/3)/Uref *100
TISamples = np.sqrt(tkeSamples*2/3)/UrefSamples *100
TISamples = TISamples #- TISamples.min(axis=1, keepdims=True)
TIMean    = np.mean(TISamples, axis=0)
TISigma   = np.std(TISamples, axis=0)

TIMin = TIDet.min(axis=0) * 0.

# %% Figures settings
# rcParams.update({'figure.figsize': [15,7]})
myUQlib.rcParamsSettings(22.5)

DETClr  = (0.6350, 0.0780, 0.1840)
meanClr = 'b'
fillClr = 'b'
LESClr  = 'k'

MC = False
MC = True

overlap = True
overClr=DETClr
overAlpha=0.25
scale=1.025
lw_1 = 1 if overlap else 0

if MC: N=2 # 2 for k-eps, 2 for LRR
else: N=1 # 1 for k-eps, 0.5 for LRR

fig, ax = plt.subplots(ncols=numLines, nrows=2, constrained_layout=True, 
                       figsize=(15,6), sharey=True)

row = 0
for l in range(numLines):
    axIdx = (row,l)

    ax[axIdx].plot(defULES[l][:,0], defULES[l][:,1], color=LESClr, marker="o", 
                   mfc='none', lw=0, ms=4.5)
    if overlap==False:
        ax[axIdx].plot(defUDet[:,axIdx[1]], yByD, color=DETClr)

    if MC:
        ax[axIdx].plot(defUMean[:,axIdx[1]], yByD, color=meanClr)
        ax[axIdx].fill_betweenx(yByD, defUMean[:,axIdx[1]] + N*defUSigma[:,axIdx[1]], \
                                      defUMean[:,axIdx[1]] - N*defUSigma[:,axIdx[1]], \
                                      alpha=0.2, linewidth=lw_1, color='b')
        if overlap:
            ax[axIdx].plot(defU0PCE[:,axIdx[1]] *scale, yByD, color=overClr)
            ax[axIdx].fill_betweenx(yByD, defU0PCE[:,axIdx[1]] *scale + N/3*defUSigmaPCE[:,axIdx[1]], \
                                          defU0PCE[:,axIdx[1]] *scale - N/3*defUSigmaPCE[:,axIdx[1]], \
                                          alpha=overAlpha, linewidth=lw_1, color=overClr)
    else:
        ax[axIdx].plot(defU0PCE[:,axIdx[1]], yByD, color=meanClr)
        ax[axIdx].fill_betweenx(yByD, defU0PCE[:,axIdx[1]] + N*defUSigmaPCE[:,axIdx[1]], \
                                      defU0PCE[:,axIdx[1]] - N*defUSigmaPCE[:,axIdx[1]], \
                                      alpha=0.2, linewidth=0, color='b')

    ax[axIdx].axhline(0, color='k', ls="-.")
    ax[axIdx].set_ylim(bottom=-1,top=1)
    ax[axIdx].set_xlim(left=-0.1,right=0.7)
    ax[axIdx].set_xticks([0.0, 0.2, 0.4, 0.6])
    ax[axIdx].set_xlabel('$\Delta u / u_h$')
    if l==0: ax[axIdx].set_ylabel('$y/D$')
    ax[axIdx].spines['right'].set_visible(False)
    ax[axIdx].spines['top'].set_visible(False)
    ax[axIdx].set_title('WT ' + str(lines[l]), pad=20)

row = 1
for l in range(numLines):
    axIdx = (row,l)

    ax[axIdx].plot(TILES[l][:,0], TILES[l][:,1], color=LESClr, marker="o", 
                   mfc='none', lw=0, ms=4.5)
    if overlap==False:
        ax[axIdx].plot(TIDet[:,axIdx[1]], yByD, color=DETClr)

    if MC:
        ax[axIdx].plot(TIMean[:,axIdx[1]] + TIMin[l], yByD, color=meanClr)
        ax[axIdx].fill_betweenx(yByD, TIMean[:,axIdx[1]] + N*TISigma[:,axIdx[1]] + TIMin[l], \
                                      TIMean[:,axIdx[1]] - N*TISigma[:,axIdx[1]] + TIMin[l], \
                                      alpha=0.2, linewidth=lw_1, color='b')
        if overlap:
            ax[axIdx].plot(TI0PCE[:,axIdx[1]] *scale + TIMin[l], yByD, color=overClr)
            ax[axIdx].fill_betweenx(yByD, TI0PCE[:,axIdx[1]] *scale + N/3*TISigmaPCE[:,axIdx[1]] + TIMin[l], \
                                          TI0PCE[:,axIdx[1]] *scale - N/3*TISigmaPCE[:,axIdx[1]] + TIMin[l], \
                                          alpha=0.2, linewidth=lw_1, color=overClr)
    else:
        ax[axIdx].plot(TI0PCE[:,axIdx[1]] + TIMin[l], yByD, color=meanClr)
        ax[axIdx].fill_betweenx(yByD, TI0PCE[:,axIdx[1]] + N*TISigmaPCE[:,axIdx[1]] + TIMin[l], \
                                      TI0PCE[:,axIdx[1]] - N*TISigmaPCE[:,axIdx[1]] + TIMin[l], \
                                      alpha=0.2, linewidth=0, color='b')

    ax[axIdx].axhline(0, color='k', ls="-.")
    ax[axIdx].set_ylim(bottom=-1,top=1)
    ax[axIdx].set_xlim(left=2.5,right=17.5)
    ax[axIdx].set_xticks([5, 10, 15])
    ax[axIdx].set_xlabel('$I [\\%]$')
    if l==0: ax[axIdx].set_ylabel('$y/D$')
    ax[axIdx].spines['right'].set_visible(False)
    ax[axIdx].spines['top'].set_visible(False)

if overlap==False:
    fig.legend(['LES', 'DET', r'$\textbf{E}[\bullet]$',
                r'$\textbf{E}[\bullet] \, \pm \, 2\sqrt{\textbf{V}[\bullet]}$'],
                ncol=4, frameon=False, columnspacing=0.75,
                loc='upper center', bbox_to_anchor=(0.5, 1.15))

if MC and overlap:
    fig.savefig(DATA+casePngDir+'/defU_TI_xByD_5_WT1to6_MC_PCE.png', dpi=150, bbox_inches='tight')
elif MC and overlap==False:
    fig.savefig(DATA+casePngDir+'/defU_TI_xByD_5_WT1to6_MC.png', dpi=150, bbox_inches='tight')
else:
    fig.savefig(DATA+casePngDir+'/defU_TI_xByD_5_WT1to6_PCE.png', dpi=150, bbox_inches='tight')