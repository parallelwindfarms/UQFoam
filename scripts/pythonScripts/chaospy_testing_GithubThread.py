#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 16:33:15 2020

@author: jigar
"""


import chaospy as cp
import numpy as np
import time
from joblib import Parallel, delayed

n = 3
dn = 12
q = 1
# dg = 1

dist_ij   = cp.Iid(cp.Normal(0,1), dn)
# dist_ii   = cp.Iid(cp.Gamma(1), dg)

# dist      = cp.J(*dist_ij)#, *dist_ii)
dist      = dist_ij
polyExpsn = cp.orth_ttr(n, dist, retall=False, normed=True, cross_truncation=1/q)

N = len(polyExpsn)


# <codecell> Triple outer pdt
start_time_loop = time.time()
phi3 = cp.outer(polyExpsn, polyExpsn, polyExpsn)
Cijk_org  = cp.E(phi3, dist)
print("Loop time = ", time.time() - start_time_loop)




# <codecell> Looping ijk
Cijk = np.zeros(((N,N,N)))
start_time_loop = time.time()

for i in range(N):
    start_time = time.time()
    for j in range(N):
        if i>=j:
            for k in range(N):
                if j>=k:
                    
                    if Cijk[i,j,k]==0.0:
                        tmpCijk = cp.E(polyExpsn[i]*polyExpsn[j]*polyExpsn[k], dist) 
                        
                        if abs(tmpCijk)>=1e-12:
                            Cijk[i,j,k] = tmpCijk
                            if j!=k:
                                Cijk[i,k,j] = tmpCijk
                            
                            Cijk[j,i,k] = tmpCijk
                            if k!=i:
                                Cijk[j,k,i] = tmpCijk
                            
                            Cijk[k,i,j] = tmpCijk
                            if i!=j:
                                Cijk[k,j,i] = tmpCijk
                                
    print(i, "\t t = ", time.time() - start_time)
        
print("Loop time = ", time.time() - start_time_loop)

# <codecell> Looping ijk - parallel
def triplePdt3(i, polyExpsn, dist, N):
    start_time = time.time()
    tmpCijk = np.zeros((N,N))
    for j in range(N):
        if i>=j:
            for k in range(N):
                if j>=k:
                    tmpCijk[j,k] = cp.E(polyExpsn[i]*polyExpsn[j]*polyExpsn[k], dist)
    print(i, "\t t = ", time.time() - start_time)
    return tmpCijk

start_time_loop = time.time()
Cijk = np.array(Parallel(n_jobs=4, backend='multiprocessing')\
    (delayed(triplePdt3)(i, polyExpsn, dist, N) for i in range(N)))
print("Loop time = ", time.time() - start_time_loop)

for i in range(N):
    start_time = time.time()
    for j in range(N):
        if i>=j:
            for k in range(N):
                if j>=k:
                    if abs(Cijk[i,j,k])>=1e-12:
                        tmpCijk = Cijk[i,j,k]
                        if j!=k:
                            Cijk[i,k,j] = tmpCijk
                        Cijk[j,i,k] = tmpCijk
                        if k!=i:
                            Cijk[j,k,i] = tmpCijk                        
                        Cijk[k,i,j] = tmpCijk
                        if i!=j:
                            Cijk[k,j,i] = tmpCijk
                                        
print("Loop time = ", time.time() - start_time_loop)
    
# <codecell>
print("Max diff = ",np.max(Cijk_org-Cijk))
print("Min diff = ",np.min(Cijk_org-Cijk))



    
# <codecell> Double outer pdt
Cijk_new = np.zeros(((N,N,N)))
phi2 = cp.outer(polyExpsn, polyExpsn)

start_time_loop = time.time()
for i in range(N):
    start_time = time.time()
    phi3 = polyExpsn[i] * phi2
    Cijk_new[i] = cp.E(phi3, dist)
    print(i, "\t t = ", time.time() - start_time)
print("Loop time = ", time.time() - start_time_loop)

# <codecell>
print("Max diff = ",np.max(Cijk_org-Cijk_new))
print("Min diff = ",np.min(Cijk_org-Cijk_new))

# <codecell> Double outer pdt - parallel
def triplePdt1a(i, phi3, dist):
    start_time = time.time()
    tmpCijk = cp.E(phi3, dist)
    print(i, "\t t = ", time.time() - start_time)
    return tmpCijk


def triplePdt1b(i, polyExpsn, dist, N, stp):
    start_time = time.time()
    tmpCijk = np.zeros((N,N))
    for j in range(int(N/stp)):
        phi2 = cp.outer(polyExpsn[stp*j:stp*(j+1)], polyExpsn)
        tmpCijk[stp*j:stp*(j+1)] = cp.E(polyExpsn[i]*phi2, dist)
        if j==(int(N/stp)-1) and N>stp*(j+1):
            phi2 = cp.outer(polyExpsn[stp*(j+1):N], polyExpsn)
            tmpCijk[stp*(j+1):N] = cp.E(polyExpsn[i]*phi2, dist)
        print("j = ",j)
    print(i, "\t t = ", time.time() - start_time)
    return tmpCijk


def triplePdt2(i, polyExpsn, dist, N):
    start_time = time.time()
    tmpCijk = np.zeros((N,N))
    for j in range(N):
        tmpCijk[j] = cp.E(polyExpsn[i]*polyExpsn[j]*polyExpsn, dist)
    print(i, "\t t = ", time.time() - start_time)
    return tmpCijk



start_time_loop = time.time()

# phi2 = cp.outer(polyExpsn, polyExpsn)
# Cijk_par = np.array(Parallel(n_jobs=4, backend='multiprocessing')\
#     (delayed(triplePdt1a)(i, polyExpsn[i]*phi2, dist) for i in range(N)))
 
# Cijk_par = np.array(Parallel(n_jobs=4, backend='multiprocessing')\
#     (delayed(triplePdt1b)(i, polyExpsn, dist, N, 5) for i in range(N)))
   
Cijk_par = np.array(Parallel(n_jobs=4, backend='multiprocessing')\
    (delayed(triplePdt2)(i, polyExpsn, dist, N) for i in range(N)))
    
print("Loop time = ", time.time() - start_time_loop)

# <codecell>
print("Max diff = ",np.max(Cijk_org-Cijk_par))
print("Min diff = ",np.min(Cijk_org-Cijk_par))








    
    
    