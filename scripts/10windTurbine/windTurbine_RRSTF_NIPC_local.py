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
import pandas as pd

from multiprocessing import Pool
from functools import partial
from timeit import default_timer as timer

SCRIPTS = os.environ['SCRIPTS']
DATA    = os.environ['windTurbineData']

sys.path.append(SCRIPTS+'/10windTurbine')
sys.path.append(SCRIPTS+'/pythonScripts')

cwd = os.getcwd() + "/"

import myUQlib

# %% Case details
caseName   = cwd.split('/')[-2]
# casePngDir = '/2RRSTF/NIPC/uqPaperCases/' + caseName+'/png'
casePngDir = '/2RRSTF/NIPC/' + caseName+'/png'

if (os.path.exists(DATA+casePngDir))==False:
    print('Making new case data directory...')
    os.makedirs(DATA+casePngDir, exist_ok=True)

# %% Parametrs
D = 80
h = 70
Uref = 8

# center of AD
ADloc = (0,0,h)

# mesh dependent
nx = 150
Wdx = D/20

# %% Setting the PCE parameters
nSamplesOrg = 30

# delB = cp.Uniform(0,1) # w/ Uniform Sampling
# delB = cp.Gamma(2,0.15,0) # w/ Gamma Sampling
d = 5
delB = cp.Iid(cp.Normal(0,1.25), d)

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
vtkFile   = 'baseCase/VTK/baseCase_1000.vtk'
# s = 80
# vtkFile = 'sample_'+str(s)+'/VTK/'+'sample_'+str(s)+'_1000.vtk'

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
TIDetAvgd   = np.average(np.reshape(TIDet-TIDet[0], (nCellsInDisk, nx)), axis=0)

# %% Importing LES, RANS, UQRANS data
U_LES  = np.loadtxt(DATA+"/0LES/Vestas2MWV80/LES_defU_centerline.csv")
TI_LES = np.loadtxt(DATA+"/0LES/Vestas2MWV80/LES_TI_centerline.csv")

def getSamplesOverDiskVol(cIdxInDisk, s):
    print(s)
    vtkFile = 'sample_'+str(s)+'/VTK/'+'sample_'+str(s)+'_1000.vtk'
    mesh    = pv.UnstructuredGrid(cwd + vtkFile)
    UMag    = np.linalg.norm(mesh.cell_data['U'][cIdxInDisk], axis=1)
    tke     = mesh.cell_data['k'][cIdxInDisk]
    return UMag, tke

pool = Pool(processes=32)
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
TISamples = TISamples - TISamples[:,0].reshape(-1,1)
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
TISamples = TISamples - TISamples[:,0].reshape(-1,1)
TIMean    = np.mean(TISamples, axis=0)
TISigma   = np.std(TISamples, axis=0)

TIMeanAvgd  = np.average(np.reshape(TIMean, (nCellsInDisk, nx)), axis=0)
TISigmaAvgd = np.average(np.reshape(TISigma, (nCellsInDisk, nx)), axis=0)

# %% Figures settings
# myUQlib.rcParamsSettings(15)

DETClr  = (0.6350, 0.0780, 0.1840)
meanClr = 'b'
fillClr = 'b'
LESClr  = 'k'
diskClr = 'k'

diskSpan = (ADloc[0]/D-Wdx/D/2, ADloc[0]/D+Wdx/D/2)

# Plotting averaged fields

MC = False
MC = True


if MC: N=2
else: N=2

fig, ax = plt.subplots(ncols=1, nrows=2, constrained_layout=True, sharex=True, figsize=(10,6))

axIdx = 0
ax[axIdx].plot(U_LES[:-1,0], U_LES[:-1,1], color=LESClr, marker="o", mfc='none', lw=0, ms=3)
ax[axIdx].plot(x_D, defUDetAvgd, color=DETClr)

if MC==True:
    ax[axIdx].plot(x_D, defUMeanAvgd, color=meanClr)
    ax[axIdx].fill_between(x_D, defUMeanAvgd + N*defUSigmaAvgd, \
                                defUMeanAvgd - N*defUSigmaAvgd, \
                                alpha=0.2, linewidth=0, color='b')
else:
    ax[axIdx].plot(x_D, defU0PCEAvgd, color=meanClr)
    ax[axIdx].fill_between(x_D, defU0PCEAvgd + N*defUSigmaPCEAvgd, \
                                defU0PCEAvgd - N*defUSigmaPCEAvgd, \
                                alpha=0.2, linewidth=0, color='b')

ax[axIdx].axvspan(*diskSpan, linewidth=0, color=diskClr)
ax[axIdx].set_ylim(bottom=0.0,top=0.6)
ax[axIdx].set_yticks([0.0, 0.2, 0.4, 0.6])
# ax[axIdx].set_xlim(left=-4,right=25)
# ax[axIdx].set_xlabel('$x/D$')
ax[axIdx].set_ylabel('$\Delta u / u_h$')
ax[axIdx].spines['right'].set_visible(False)
ax[axIdx].spines['top'].set_visible(False)

axIdx = 1
ax[axIdx].plot(TI_LES[:-1,0], TI_LES[:-1,1]-TI_LES[:-1,1].min(), 
               color=LESClr, marker="o", mfc='none', lw=0, ms=3)
ax[axIdx].plot(x_D, TIDetAvgd, color=DETClr)

if MC==True:
    ax[axIdx].plot(x_D, TIMeanAvgd, color=meanClr)
    ax[axIdx].fill_between(x_D, TIMeanAvgd + N*TISigmaAvgd,
                                TIMeanAvgd - N*TISigmaAvgd, \
                                alpha=0.2, linewidth=0, color='b')
else:
    ax[axIdx].plot(x_D, TI0PCEAvgd, color=meanClr)
    ax[axIdx].fill_between(x_D, TI0PCEAvgd + N*TISigmaPCEAvgd,
                                TI0PCEAvgd - N*TISigmaPCEAvgd, \
                                alpha=0.2, linewidth=0, color='b')

ax[axIdx].axvspan(*diskSpan, linewidth=0, color=diskClr)
ax[axIdx].set_ylim(bottom=0,top=12)
ax[axIdx].set_yticks([0, 4, 8, 12])
ax[axIdx].set_xlim(left=-4,right=25)
ax[axIdx].set_xticks([-4, 0, 5, 10, 15, 20, 25])
ax[axIdx].set_xlabel('$x/D$')
ax[axIdx].set_ylabel('$I_{add} [\\%]$')
ax[axIdx].spines['right'].set_visible(False)
ax[axIdx].spines['top'].set_visible(False)

# ax[0].legend(['LES', 'DET', r'$\textbf{E}[\bullet]$', '_nolegend_',
#               r'$\textbf{E}[\bullet] \, \pm \, 2\sqrt{\textbf{V}[\bullet]}$'],
#               ncol=2, frameon=False)
ax[0].legend(['LES', 'DET', r'${E}[\bullet]$', '_nolegend_',
              r'${E}[\bullet] \, \pm \, 2\sqrt{{V}[\bullet]}$'],
              ncol=4, frameon=False, columnspacing=0.75,
              loc='upper center', bbox_to_anchor=(0.5, 1.15))

if MC==True:
    fig.savefig(DATA+casePngDir+'/diskAvgdFields_MC.png', dpi=300)
else:
    fig.savefig(DATA+casePngDir+'/diskAvgdFields_PCE.png', dpi=300)






# %%  ###############################################################################################
# Plot over lines at different x locations
#############################################################################################################

# %% Importing LES, RANS, UQRANS data
lines = (2,4,7,12)
numLines = len(lines)
yPts = 100
nCells = numLines*yPts

# Load LES data ###############
yByDLES = np.loadtxt(DATA+"/0LES/Vestas2MWV80/LES_xByD_2_defU.csv")[:,1]
defULES = np.loadtxt(DATA+"/0LES/Vestas2MWV80/LES_xByD_2_defU.csv")[:,0].reshape((len(yByDLES),1))
TILES   = np.loadtxt(DATA+"/0LES/Vestas2MWV80/LES_xByD_2_TI.csv")[:,0].reshape((len(yByDLES),1))

for l in range(1,numLines):
    fileName = DATA+"/0LES/Vestas2MWV80/LES_xByD_"+str(lines[l])
    defULES  = np.hstack( (defULES, np.loadtxt(fileName+"_defU.csv")[:,0].reshape((len(yByDLES),1))) )
    TILES    = np.hstack( (TILES, np.loadtxt(fileName+"_TI.csv")[:,0].reshape((len(yByDLES),1))) )

# Load DET data ###############
baseCase = 'baseCase/'
sampleDir = 'postProcessing/sample/1000/'

yByD    = np.loadtxt(baseCase + sampleDir + "X_by_D_2_U.xy")[:,0] / D
UDetMag = np.linalg.norm(np.loadtxt(baseCase + sampleDir + "X_by_D_2_U.xy")[:,1:], axis=1, keepdims=True)
tkeDet  = (np.loadtxt(baseCase + sampleDir + "X_by_D_2_k_nut_p.xy")[:,1]).reshape((yPts,1))

for l in range(1,numLines):
    fileName = baseCase + sampleDir+'/X_by_D_'+str(lines[l])
    UDetMag  = np.hstack( (UDetMag, np.linalg.norm(np.loadtxt(fileName+'_U.xy')[:,1:], axis=1, keepdims=True)) )
    tkeDet   = np.hstack( (tkeDet, (np.loadtxt(fileName+'_k_nut_p.xy')[:,1]).reshape((yPts,1))) )

# Load UQRANS data ###############
def getSamplesOverLines(sampleDir, numLines, lines, yPts, s):
    print(s)
    baseCase = 'sample_'+str(s)+'/'
    UMag = np.linalg.norm(np.loadtxt(baseCase + sampleDir + "X_by_D_2_U.xy")[:,1:], axis=1, keepdims=True)
    tke  = (np.loadtxt(baseCase + sampleDir + "X_by_D_2_k_nut_p.xy")[:,1]).reshape((yPts,1))

    for l in range(1,numLines):
        fileName = baseCase + sampleDir+'/X_by_D_'+str(lines[l])
        UMag  = np.hstack( (UMag, np.linalg.norm(np.loadtxt(fileName+'_U.xy')[:,1:], axis=1, keepdims=True)) )
        tke   = np.hstack( (tke, (np.loadtxt(fileName+'_k_nut_p.xy')[:,1]).reshape((yPts,1))) )

    return UMag, tke

pool = Pool(processes=4)
start = timer()
data = pool.map(partial(getSamplesOverLines, sampleDir, numLines, lines, yPts), (nSamples))
USamplesMag = np.array([data[s][0] for s in range(len(nSamples))])
tkeSamples = np.array([data[s][1] for s in range(len(nSamples))])
pool.close()
pool.join()
print(timer()-start, 's')

# %% PCE of U
defUSamples   = (Uref-USamplesMag)/Uref

defU_PCEModes = np.zeros((Pplus1, yPts, numLines))
for l in range(numLines):
    defU_PCEApp = cp.fit_regression(phi, delBSamples, defUSamples[:,:,l])
    defU_PCEModes[:,:,l] = (defU_PCEApp.coefficients).T

defU0PCE     = defU_PCEModes[0]
defUSigmaPCE = np.sqrt(np.sum(defU_PCEModes[1:]**2, axis=0))

## %% PCE of TI
TISamples   = np.sqrt(tkeSamples*2/3)/Uref *100

TI_PCEModes = np.zeros((Pplus1, yPts, numLines))
for l in range(numLines):
    TI_PCEApp = cp.fit_regression(phi, delBSamples, TISamples[:,:,l])
    TI_PCEModes[:,:,l] = (TI_PCEApp.coefficients).T
TI_PCEModes = np.array(TI_PCEModes)

TI0PCE     = TI_PCEModes[0]
TISigmaPCE = np.sqrt(np.sum(TI_PCEModes[1:]**2, axis=0))

# %% Computing averaged U fields
USamplesMagMean  = np.mean(USamplesMag, axis=0)
USamplesMagSigma = np.std(USamplesMag, axis=0)

defUDet   = (Uref-UDetMag)/Uref
defUMean  = (Uref-USamplesMagMean)/Uref
defUSigma = USamplesMagSigma/Uref

# %% Computing averaged TI fields
TIDet     = np.sqrt(tkeDet*2/3)/Uref *100
TISamples = np.sqrt(tkeSamples*2/3)/Uref *100
TIMean  = np.mean(TISamples, axis=0)
TISigma = np.std(TISamples, axis=0)

# %% Figures settings
myUQlib.rcParamsSettings(15)

DETClr  = (0.6350, 0.0780, 0.1840)
meanClr = 'b'
fillClr = 'b'
LESClr  = 'k'

MC = False
MC = True

if MC: N=2
else:  N=2

fig, ax = plt.subplots(ncols=numLines, nrows=2, constrained_layout=True, figsize=(10,6), sharey=True)

row = 0
for l in range(numLines):
    axIdx = (row,l)

    ax[axIdx].plot(defULES[:,axIdx[1]], yByDLES, color=LESClr, marker="o", mfc='none', lw=0, ms=3)
    ax[axIdx].plot(defUDet[:,axIdx[1]], yByD, color=DETClr)

    if MC:
        ax[axIdx].plot(defUMean[:,axIdx[1]], yByD, color=meanClr)
        ax[axIdx].fill_betweenx(yByD, defUMean[:,axIdx[1]] + N*defUSigma[:,axIdx[1]], \
                                      defUMean[:,axIdx[1]] - N*defUSigma[:,axIdx[1]], \
                                      alpha=0.2, linewidth=0, color='b')
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
    ax[axIdx].set_title('$x/D=$ ' + str(lines[l]), pad=20)

row = 1
for l in range(numLines):
    axIdx = (row,l)

    ax[axIdx].plot(TILES[:,axIdx[1]], yByDLES, color=LESClr, marker="o", mfc='none', lw=0, ms=3)
    ax[axIdx].plot(TIDet[:,axIdx[1]], yByD, color=DETClr)

    if MC:
        ax[axIdx].plot(TIMean[:,axIdx[1]], yByD, color=meanClr)
        ax[axIdx].fill_betweenx(yByD, TIMean[:,axIdx[1]] + N*TISigma[:,axIdx[1]], \
                                      TIMean[:,axIdx[1]] - N*TISigma[:,axIdx[1]], \
                                      alpha=0.2, linewidth=0, color='b')
    else:
        ax[axIdx].plot(TI0PCE[:,axIdx[1]], yByD, color=meanClr)
        ax[axIdx].fill_betweenx(yByD, TI0PCE[:,axIdx[1]] + N*TISigmaPCE[:,axIdx[1]], \
                                      TI0PCE[:,axIdx[1]] - N*TISigmaPCE[:,axIdx[1]], \
                                      alpha=0.2, linewidth=0, color='b')

    ax[axIdx].axhline(0, color='k', ls="-.")
    ax[axIdx].set_ylim(bottom=-1,top=1)
    ax[axIdx].set_xlim(left=0,right=15)
    ax[axIdx].set_xticks([0, 5, 10, 15])
    ax[axIdx].set_xlabel('$I [\\%]$')
    if l==0: ax[axIdx].set_ylabel('$y/D$')
    ax[axIdx].spines['right'].set_visible(False)
    ax[axIdx].spines['top'].set_visible(False)
    # ax[axIdx].set_title('$x/D=$ ' + str(lines[l]))

fig.legend(['LES', 'DET', r'$\textbf{E}[\bullet]$', '_nolegend_',
            r'$\textbf{E}[\bullet] \, \pm \, 2\sqrt{\textbf{V}[\bullet]}$'],
            ncol=4, frameon=False, columnspacing=0.75,
            loc='upper center', bbox_to_anchor=(0.5, 1.15))

if MC:
    fig.savefig(DATA+casePngDir+'/defU_TI_xByD_2_4_7_12_MC.png', dpi=300, bbox_inches='tight')
else:
    fig.savefig(DATA+casePngDir+'/defU_TI_xByD_2_4_7_12_PCE.png', dpi=300, bbox_inches='tight')

