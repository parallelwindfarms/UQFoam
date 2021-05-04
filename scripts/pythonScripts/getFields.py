#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 25 13:36:59 2021

@author: jigar
"""

import numpy as np
import chaospy as cp

# <codecell> 
"""
getSmaple funtion for Gausaain process
"""
def getGpSample(kleModes, quadNode, dim, numCells, sigma):
    tmpSample = np.zeros(numCells)
    for j in range(dim):
        tmpSample += kleModes[j] * quadNode[j]
    return  sigma*tmpSample


# <codecell> 
"""
getSmaple funtion for a field (any) - cell value
"""
def getGpSampleCellVal(kleModes, quadNode, dim, sigma):
    tmpSample = 0
    for j in range(dim):
        tmpSample += kleModes[j] * quadNode[j]
    return  sigma*tmpSample


# <codecell> 
"""
getSmaple funtion for logN field
"""
def getExpGpSample(kleModes, quadNode, dim, numCells, sigma):
    tmpSample = np.zeros(numCells)
    for j in range(dim):
        tmpSample += (kleModes[j] * quadNode[j])
    return  np.exp(sigma*tmpSample)

    
# <codecell> Coeffs of PCE of Gamma dist wrt Normal dist
def getPCEGamma(u, n_u):
    
    xi = cp.Normal(0,1)
    
    nrDiv = 1000
    min_w = -10
    max_w = 10
    dw    = (max_w-min_w)/nrDiv
    w     = np.linspace(min_w, max_w, nrDiv+1)
    w_c   = (w[1:] + w[0:-1])/2
    dw    = w[1:] - w[0:-1]
    
    herm,nrms = cp.orth_ttr(n_u,xi,retall=True,normed=True)
    Np = len(herm)
    
    Ub = np.zeros(Np)
    
    for Np_i in range(Np):
        Ub[Np_i] = sum(u.inv(xi.cdf(w_c)) * herm[Np_i](w_c) * xi.pdf(w_c) * dw)

    return Ub


# <codecell> 
"""
PCE of logN field (Ghanem 2002) with vol avgd solbol indices
"""
def getlogNFieldPCE(mean,stD,nCells,d,g,volAvg,phyDim,expoMat,P,dist):
    
    Pplus1 = P+1
    
    
    xiSamples = (dist.sample(size=10000, rule="L")).T
    xiSamplesMean = xiSamples.mean(axis=0)
    
    temp = g[0]*0.0
    for i in range(d):
        temp += g[i]*xiSamplesMean[i]
    # print(temp)
    mode_0 = np.exp(temp)
    
    
    # mode_0 = np.exp(mean + 0.5*np.square(stD))
    modes  = np.zeros((Pplus1, nCells))        
    var    = np.zeros((1, nCells))
    
    for i in range(Pplus1):
        
        gPdt_i = 1.0
        p_i    = sum(expoMat[i])
        
        for j in range(d):
            gPdt_i *= np.power(g[j], expoMat[i][j])
            
        modes[i] = mode_0 * np.power(stD,p_i) * gPdt_i
        
        if i>0:
            var += np.square(modes[i]).T
        
    return modes, var
