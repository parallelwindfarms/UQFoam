#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 29 15:21:21 2021

@author: jigar
"""

# <codecell> Import Required Modules and env vars
from __future__ import print_function

import os, sys
import numpy as np
import pyvista as pv
import matplotlib.pyplot as plt
import chaospy as cp

SCRIPTS = os.environ['SCRIPTS']
DATA    = os.environ['ductFlowData']

sys.path.append(SCRIPTS+'/9ductFlow')

cwd = os.getcwd() + "/"

# <codecell> Plot settings
# Font, Size
from matplotlib import rcParams
allFontSize = 15
params = {
   'axes.labelsize': allFontSize,
   'legend.fontsize': allFontSize,
   'xtick.labelsize': allFontSize,
   'ytick.labelsize': allFontSize,
   'text.usetex': True,
   'xtick.direction':'in',
   'ytick.direction':'in',
   'font.family':'serif',
   'font.serif': ['Computer Modern Roman']
   }
rcParams.update(params)

# <codecell> Case details
caseName   = '1URANS_RKLEgPC_dp0.2_stst'
casePngDir = '/uqPaperCases/'+caseName+'/png'

if (os.path.exists(DATA+casePngDir))==False:
    print('Making new case data directory...')
    os.makedirs(DATA+casePngDir, exist_ok=True)
    
# <codecell> Read OpenFOAM case mesh access to grid and cell values
vtkFile = 'VTK/1URANS_RKLEgPC_dp0.1_minimalOutput_500.vtk'
# vtkFile = 'VTK/'+caseName+'_500.vtk'
mesh    = pv.UnstructuredGrid(cwd + vtkFile)
nCells  = mesh.n_cells
cCenter = mesh.cell_centers().points

# <codecell> Importing DNS and RANS vectorsData
DNSdata    = np.loadtxt(DATA+'/0DNS/DATA_Re_1100/Re_1100_data.txt')
UQMeanData = mesh.cell_arrays['U0']
UQStdData  = mesh.cell_arrays['USigma']
DETRASdata = mesh.cell_arrays['U']

# <codecell> Making quiver plots
from scipy import interpolate

rcParams.update({'figure.figsize': [10.5,5]})

fig = plt.figure()
ticks = [0.0,0.2,0.4,0.6,0.8,1.0]
yLine = np.linspace(0,1,50)
x1, x2, x3 = 0.25, 0.50, 0.75

# UQ_RAS
cell_idx = ((cCenter[:,0]>1) & (cCenter[:,1]>1))
X_RAS, Y_RAS = np.array(cCenter[cell_idx,0]-1),np.array(cCenter[cell_idx,1]-1)
U, V = UQMeanData[cell_idx,0], UQMeanData[cell_idx,1]
ax2 = fig.add_subplot(121)
Q = ax2.quiver(X_RAS, Y_RAS, -U, -V, scale=0.06)
ax2.set_xticks(ticks)
ax2.set_yticks(ticks)
ax2.set_ylim(bottom=0,top=1)
ax2.set_xlim(left=0,right=1)

# Interpolating for plotting over line (next codecell)
f = interpolate.interp2d(X_RAS, Y_RAS, U, kind='linear')
U_UQ_lines = np.hstack([f(x1, yLine), f(x2, yLine), f(x3, yLine)])
f = interpolate.interp2d(X_RAS, Y_RAS, UQStdData[cell_idx,0], kind='linear')
UStd_UQ_lines = np.hstack([f(x1, yLine), f(x2, yLine), f(x3, yLine)])
f = interpolate.interp2d(X_RAS, Y_RAS, DETRASdata[cell_idx,0], kind='linear')
U_DET_lines = np.hstack([f(x1, yLine), f(x2, yLine), f(x3, yLine)])

# DNS
X_DNS, Y_DNS = DNSdata[:,0], DNSdata[:,1]
U, V = DNSdata[:,3], DNSdata[:,4]
ax1 = fig.add_subplot(122)
Q = ax1.quiver(X_DNS, Y_DNS, U, V, scale=0.25)
ax1.set_xticks(ticks)
ax1.set_yticks(ticks)
ax1.set_ylim(bottom=0,top=1)
ax1.set_xlim(left=0,right=1)

# Interpolating for plotting over line (next codecell)
f = interpolate.interp2d(Y_DNS, X_DNS, V, kind='linear')
U_DNS_lines = np.hstack([f(x1, yLine), f(x2, yLine), f(x3, yLine)])


# fig.savefig(DATA+casePngDir+'/vectorPlot.png', dpi=300)

# <codecell> Making line plots
rcParams.update({'figure.figsize': [15,4.5]})
DETClr  = (0.6350, 0.0780, 0.1840)
meanClr = 'b'
fillClr = 'b'
DNSClr  = 'k'

fig = plt.figure()
ax  = (fig.add_subplot(131), fig.add_subplot(132), fig.add_subplot(133))

numLines = 3
uTicks = [-0.04,-0.02,0.0,0.02,0.04]

stdN = 1.25
uN   = -2
ub = U_UQ_lines*uN + stdN*UStd_UQ_lines
lb = U_UQ_lines*uN - stdN*UStd_UQ_lines

for line in range(numLines):
    ax[line].plot(U_DNS_lines[:,line],    yLine, color=DNSClr)    # DNS
    ax[line].plot(U_DET_lines[:,line],    yLine, color=DETClr)    # DET
    ax[line].plot(U_UQ_lines[:,line]*uN, yLine, color=meanClr)    # UQ_Mean
    ax[line].fill_betweenx(yLine, ub[:,line], lb[:,line], \
                           alpha=0.2, linewidth=0, color='b')     # UQ Bounds
    ax[line].set_xticks(uTicks)
    ax[line].set_yticks(ticks)
    ax[line].set_xlabel('$u$')
    ax[line].set_ylabel('$y/H$')
    ax[line].spines['right'].set_visible(False)
    ax[line].spines['top'].set_visible(False)
    ax[line].set_ylim(bottom=0,top=1)
    ax[line].set_xlim(left=-0.05,right=0.05)
    
fig.legend(['DNS', 'DET', r'$\textbf{E}[u]$',
            r'$\textbf{E}[u] \, \pm \, 2\sqrt{\textbf{V}[u]}$'], 
            loc='upper center', ncol=4, frameon=False, 
            borderaxespad=-0.4, columnspacing=0.75)


# fig.savefig(DATA+casePngDir+'/U_xByD_DNS.png', dpi=300)

# <codecell> Useful functions
def symmTensorToTensorv2012(A):
    return np.array([ [A[0],A[3],A[5]],[A[3],A[1],A[4]],[A[5],A[4],A[2]] ])

def symmTensorToTensorv1806(A):
    return np.array([ [A[0],A[1],A[2]],[A[1],A[3],A[4]],[A[2],A[4],A[5]] ])

def bayCentricCoordinates(eigVals):
    eigVals.sort()
    eigVals = eigVals[::-1]
    # print('eigVals =', eigVals)
    
    C = [ eigVals[0] - eigVals[1] ]
    C.append(2*(eigVals[1] - eigVals[2]))
    C.append(3*eigVals[2] + 1)
    # print('C = ', C)
    
    e1 = np.array([1,0]); e2 = np.array([-1,0]); e3 = np.array([0,np.sqrt(3)])
    C_vec = C[0]*e1 + C[1]*e2 + C[2]*e3
    
    return C_vec

# <codecell> Loading Reynolds-stress tensor fields
selectedCellIdx = 10*32 + 10 # X_RAS = Y_RAS ~ 0.5
# for selectedCellIdx in range(10,100):
print('selectedCellIdx, x/H, y/H = ', selectedCellIdx, X_RAS[selectedCellIdx], Y_RAS[selectedCellIdx])

# R_DNS    = (np.loadtxt(DATA+'/0DNS/DATA_Re_1100/Re_1100_data.txt'))

R_UQMean = symmTensorToTensorv1806(mesh.cell_arrays['R0'][cell_idx][selectedCellIdx])
R_UQStd  = symmTensorToTensorv1806(mesh.cell_arrays['RSigma'][cell_idx][selectedCellIdx])
R_DET    = symmTensorToTensorv2012(pv.UnstructuredGrid\
            (cwd + 'VTK/'+caseName+'_500.vtk').cell_arrays['R'][cell_idx][selectedCellIdx])
        
# vtkFile_new  = 'VTK/'+caseName+'_500.vtk'
# mesh_new     = pv.UnstructuredGrid(cwd + vtkFile_new)
# R_DET    = symmTensorToTensorv2012(mesh_new.cell_arrays['R'][cell_idx][selectedCellIdx])
# R_UQMean = symmTensorToTensorv2012(mesh_new.cell_arrays['R0'][cell_idx][selectedCellIdx])
# R_UQStd  = symmTensorToTensorv2012(mesh_new.cell_arrays['RSigma'][cell_idx][selectedCellIdx])

## <codecell> Computing eigenvalues
tke_UQ  = 1.75*sum(R_UQMean.diagonal())/2
tke_DET = sum(R_DET.diagonal())/2
R_UQMean_eig,_ = np.linalg.eig(R_UQMean/(2*tke_UQ)-np.eye(3)/3)
R_DET_eig,_    = np.linalg.eig(R_DET/(2*tke_DET)-np.eye(3)/3)

## <codecell> Bary-centric triangle
C_DET    = bayCentricCoordinates(R_DET_eig)
C_UQMean = bayCentricCoordinates(R_UQMean_eig)

## <codecell> Plotting
rcParams.update({'figure.figsize': [5,5]})

fig = plt.figure()
ax = fig.add_subplot(111)

ax.plot([-1, 1],[0, 0],'k')
ax.plot([-1, 0],[0, np.sqrt(3)],'k')
ax.plot([1, 0],[0, np.sqrt(3)],'k')
ax.plot([-0.33,0],[0, np.sqrt(3)],'--k')

ax.scatter(C_DET[0], C_DET[1], color=DETClr, label='DET')
ax.scatter(C_UQMean[0], C_UQMean[1], color=meanClr, 
           label=r'$\textbf{E}[\textbf{R}]$')

stdN = 4
sampleSize = 500
RSamples = np.zeros((sampleSize,3,3))
for i in range(3):
    for j in range(3):
        if i>=j:
            dist = cp.Normal(R_UQMean[i,j], stdN * R_UQStd[i,j])
            # dist = cp.Uniform(R_UQMean[i,j]-stdN*R_UQStd[i,j],R_UQMean[i,j]+stdN*R_UQStd[i,j])
            RSamples[:,i,j] = RSamples[:,j,i] = dist.sample(sampleSize)

lbl = ''
C = np.zeros((sampleSize,2))
for s in range(sampleSize):
    tke_Sample = 1.3*sum(RSamples[s].diagonal())/2
    RSample_eig,_ = np.linalg.eig(RSamples[s]/(2*tke_Sample)-np.eye(3)/3)
    C[s] = bayCentricCoordinates(RSample_eig)
    if s==(sampleSize-1):
        lbl = 'Samples'
    if C[s,1]>0.5:
        ax.scatter(C[s,0], C[s,1], color='b', alpha=0.05, label=lbl)
            

# ax.tricontour(x, y, z, levels=8, linewidths=0.5, colors='k')
# cntr2 = ax2.tricontourf(x, y, z, levels=8, cmap="RdBu_r")

ax.legend(frameon=False)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.set_xticks([])
ax.set_yticks([])

# fig.savefig(DATA+casePngDir+'/barycentricScatter.png', dpi=300)

# <codecell>
import random 

R = np.linalg.norm(C_DET-C_UQMean)/3000
num_samples = 100
num_contours = 10
idx = np.linspace(0,num_samples,num_samples+1, dtype=int)[0:-1]
idx = random.choices(idx, k=10)

fig = plt.figure()
ax = fig.add_subplot(111)

lbl = ''
x_con = np.zeros((num_contours, len(idx)))
y_con = np.zeros((num_contours, len(idx)))
z_con = np.zeros((num_contours, len(idx)))
# theta = np.random.rand((num_samples)) * (2 * np.pi)
for c in range(num_contours):
    theta = np.random.rand((num_samples)) * (2 * np.pi)
    x = (np.exp(c))*R * np.cos(theta) + C_UQMean[0]
    y = (np.exp(c))*R * np.sin(theta) + C_UQMean[1]
    z = np.sqrt((x[idx]-C_UQMean[0])**2+(y[idx]-C_UQMean[1])**2)
    x_con[c]=x[idx]; y_con[c]=y[idx]; z_con[c]=z
    
# cmap = plt.cm.get_cmap("winter")
# cmap.set_under("yellow")
# cmap.set_over("magenta")
cs = ax.tricontourf(x_con.ravel(),y_con.ravel(),z_con.ravel(), num_contours)
cs.cmap.set_under('cyan')
cs.cmap.set_over('yellow')

ax.plot([-1, 1],[0, 0],'k')
ax.plot([-1, 0],[0, np.sqrt(3)],'k')
ax.plot([1, 0],[0, np.sqrt(3)],'k')
ax.plot([-0.33,0],[0, np.sqrt(3)],'--k')
ax.scatter(C_DET[0], C_DET[1], color=DETClr, label='DET')
ax.scatter(C_UQMean[0], C_UQMean[1], color=meanClr, 
            label=r'$\textbf{E}[\textbf{R}]$')

ax.legend(frameon=False)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.set_xticks([])
ax.set_yticks([])

# fig.savefig(DATA+casePngDir+'/barycentricContour.png', dpi=300)

# <codecell>

# <codecell>

# <codecell>


    
    
    
    
    
    