#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 25 14:03:18 2021

@author: jigar
"""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rcParams

# <codecell> Basic plot settings
DETClr  = (0.6350, 0.0780, 0.1840)
meanClr = 'b'
fillClr = 'b'
LESClr  = 'k'

# <codecell> rcParamsSettings
def rcParamsSettings(allFontSize=15, titlepad=20, usetex=True):
    params = {
       'axes.labelsize': allFontSize,
       'axes.titlesize': allFontSize,
       'axes.titlepad': titlepad,
       'legend.fontsize': allFontSize,
       'xtick.labelsize': allFontSize,
       'ytick.labelsize': allFontSize,
       'text.usetex': usetex,
       'xtick.direction':'in',
       'ytick.direction':'in',
       'font.family':'serif',
       'font.serif': ['Computer Modern Roman']
       }
    rcParams.update(params)
    
# <codecell> setSpines
def setSpines(axes, l=False,r=False,t=False,b=False):
    axes.spines['right'].set_visible(r)
    axes.spines['left'].set_visible(l)
    axes.spines['top'].set_visible(t)
    axes.spines['bottom'].set_visible(b)
 
# <codecell> setTiks
def setTiks(axes, x='none', y='none'):
    axes.xaxis.set_ticks_position(x) 
    axes.yaxis.set_ticks_position(y)
    
# <codecell> plotBaryCentricCoordinateSystem
def plotBaryCentricCoordinateSystem(ax, lw=1.5):
    ax.plot([-0.999, 0.999],[0, 0],'k', lw=lw)
    ax.plot([-1, 0],[0, np.sqrt(3)],'k', lw=lw)
    ax.plot([1, 0],[0, np.sqrt(3)],'k', lw=lw)
    ax.plot([-0.33,0],[0, np.sqrt(3)],'--k', lw=lw)
    ax.axis('off')

# <codecell> plotBaryCentricCoordinateSystemWithCmap(ax)
def plotBaryCentricCoordinateSystemWithCmap(ax):
    e1 = np.array([1,0])
    e2 = np.array([-1,0])
    e3 = np.array([0,np.sqrt(3)])
    
    nPts = 250
    grid_X = np.linspace(-1, 1, nPts)
    grid_Y = np.linspace(0, np.sqrt(3), nPts)
    
    XX, YY = np.meshgrid(grid_X, grid_Y)
    grid_XY = np.vstack((XX.ravel(), YY.ravel())).T
    
    from matplotlib.path import Path
    p = Path([e1, e2, e3])
    idxIn = p.contains_points(grid_XY)
    grid_XY = grid_XY[idxIn]
    
    clr_grid_XY = bccViridisCmap(grid_XY.reshape(1,grid_XY.shape[0],2))
    ax.scatter(grid_XY[:,0], grid_XY[:,1], c=clr_grid_XY[0], marker='.', s=20, lw=0)
    ax.axis('off')
    
    plotBaryCentricCoordinateSystem(ax, lw=2)

# <codecell> bccViridisCmap
def bccViridisCmap(C):
    cmap = plt.cm.get_cmap('viridis',3) # b, g, y
    
    v1 = np.array([1,0]);           clr_v1 = np.array(cmap(0))#np.array([[[1,0,0]]]); # b, r
    v2 = np.array([-1,0]);          clr_v2 = np.array(cmap(2))#np.array([[[0,1,0]]]); # y, g
    v3 = np.array([0,np.sqrt(3)]);  clr_v3 = np.array(cmap(1))#np.array([[[0,0,1]]])  # g, b
    
    Pt = C.reshape(-1,2)
    
    w_Pt_v1 = ((v2[1]-v3[1])*(Pt[:,0]-v3[0]) + (v3[0]-v2[0])*(Pt[:,1]-v3[1]))/ \
              ((v2[1]-v3[1])*(v1[0]-v3[0]) + (v3[0]-v2[0])*(v1[1]-v3[1]))
    w_Pt_v2 = ((v3[1]-v1[1])*(Pt[:,0]-v3[0]) + (v1[0]-v3[0])*(Pt[:,1]-v3[1]))/ \
              ((v2[1]-v3[1])*(v1[0]-v3[0]) + (v3[0]-v2[0])*(v1[1]-v3[1]))
    w_Pt_v3 = 1 - w_Pt_v1 - w_Pt_v2
    
    w_Pt_v1 = w_Pt_v1.reshape(C.shape[0], C.shape[1], 1)
    w_Pt_v2 = w_Pt_v2.reshape(C.shape[0], C.shape[1], 1)
    w_Pt_v3 = w_Pt_v3.reshape(C.shape[0], C.shape[1], 1)
    
    return (w_Pt_v1*clr_v1 + w_Pt_v2*clr_v2 + w_Pt_v3*clr_v3)/(w_Pt_v1 + w_Pt_v2 + w_Pt_v3)
    
# <codecell> drawCircle
def drawCircle(axes, x, y, r, fc='none', ec='k', ls='-.', lw='1'):
    axes.add_patch(plt.Circle((x, y), r, fc=fc, ec=ec, ls=ls, lw=lw))

# <codecell> plotDiskSpans
def plotDiskSpans(axes, ADloc, Wd, lw=1, diskClr='k'):
    for span in range(len(ADloc)):
        axes.axvspan(ADloc[span]-Wd, ADloc[span]+Wd, lw=lw, color=diskClr)

    
    