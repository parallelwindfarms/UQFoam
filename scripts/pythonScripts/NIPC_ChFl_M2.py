# -*- coding: utf-8 -*-
"""

Non-intrusive polynomail chaos
for Channel Flow with mesh M1

"""

import numpy as np
import chaospy as cp
import matplotlib.pyplot as plt
import pandas as pd
import os
cwd = os.getcwd() + "/"

# <codecell>
a = 0.075           # Lower bound
b = 0.125           # Upper bound

Cs_mean  = 0.1      # mean
Cs_sigma = 0.0144   # std - by 2 factor for std of Normal !

p_n = 3             # Poly Order
q_n = 5             # Quad pts

nu = 2e-5   # Laminar viscosity

# <codecell>
Uniform = 1
Normal  = 0

if (Uniform):
    g_dist      = cp.Uniform(-1,1)                          # Germ distribution
    Cs_dist     = cp.Uniform(a,b)                           # Cs distribution
    orthpoly    = cp.orth_ttr(p_n, Cs_dist, normed=True)    # Orthonormal polys
#    orthpoly    = cp.orth_gs(p_n, Cs_dist, normed=True)    # Orthonormal polys
    nodes, wts  = cp.generate_quadrature(q_n, Cs_dist,'C')  # Nodes and Weights

if (Normal):
    g_dist      = cp.Normal(0,1)                            # Germ distribution
    Cs_dist     = cp.Normal(Cs_mean, Cs_sigma)              # Cs distribution
    orthpoly    = cp.orth_ttr(p_n, Cs_dist, normed=True)    # Orthonormal polys
    nodes, wts  = cp.generate_quadrature(q_n, Cs_dist,'G')  # Nodes and Weights

num_nodes = len(nodes[0])
   
# <codecell>
data = cwd  + "Samples/NIPC/M2/dpdx_5.00e-5/"
# data = cwd  + "Samples/NIPC/M2/dpdx_6.24e-5/"

case = data + "1/postProcessing/collapsedFields/latesTimeDir/"
y    = np.array(pd.read_csv(case + "UMean_X.xy",  header=None, sep='\s+', usecols=[0]))

yPl_evals   = np.zeros((num_nodes,1))
U_evals     = np.zeros((num_nodes,y.size))
R_evals     = np.zeros(((4,num_nodes,y.size)))
nut_evals   = np.zeros((num_nodes,y.size))
p_evals   = np.zeros((num_nodes,y.size))

for i in range(num_nodes):
    print("Loading data from Case =",i+1)
    case = data + str(i+1) + "/postProcessing/collapsedFields/latesTimeDir/"
    yPl_evals[i][:]  = pd.read_csv(case + "yPlusMean.xy",       header=None, sep='\s+', usecols=[1]).iloc[0,0]
    U_evals[i][:]    = pd.read_csv(case + "UMean_X.xy",         header=None, sep='\s+', usecols=[1]).iloc[:,0]
    R_evals[0][i][:] = pd.read_csv(case + "UPrime2Mean_XX.xy",  header=None, sep='\s+', usecols=[1]).iloc[:,0]
    R_evals[1][i][:] = pd.read_csv(case + "UPrime2Mean_YY.xy",  header=None, sep='\s+', usecols=[1]).iloc[:,0]
    R_evals[2][i][:] = pd.read_csv(case + "UPrime2Mean_ZZ.xy",  header=None, sep='\s+', usecols=[1]).iloc[:,0]
    R_evals[3][i][:] = pd.read_csv(case + "UPrime2Mean_XY.xy",  header=None, sep='\s+', usecols=[1]).iloc[:,0]
    nut_evals[i][:]  = pd.read_csv(case + "nutMean.xy",         header=None, sep='\s+', usecols=[1]).iloc[:,0]
    p_evals[i][:]    = pd.read_csv(case + "pMean.xy",           header=None, sep='\s+', usecols=[1]).iloc[:,0]

# <codecell>
utTmp       = yPl_evals*nu/y[1]

#U_app       = cp.fit_quadrature(orthpoly, nodes, wts, np.true_divide(U_evals, utTmp))
U_app       = cp.fit_quadrature(orthpoly, nodes, wts, U_evals)
UMean0      = cp.E(U_app, Cs_dist)
UMeanSigma  = cp.Std(U_app, Cs_dist)

RMean0      = np.zeros((4,y.size))
RMeanSigma  = np.zeros((4,y.size))
for i in range(4):
    #R_app           = cp.fit_quadrature(orthpoly, nodes, wts, np.true_divide(R_evals[i], utTmp**2))
    R_app           = cp.fit_quadrature(orthpoly, nodes, wts, R_evals[i])    
    RMean0[i]       = cp.E(R_app, Cs_dist)
    RMeanSigma[i]   = cp.Std(R_app, Cs_dist)

nut_app       = cp.fit_quadrature(orthpoly, nodes, wts, nut_evals)
nutMean0      = cp.E(nut_app, Cs_dist)
nutMeanSigma  = cp.Std(nut_app, Cs_dist)

p_app       = cp.fit_quadrature(orthpoly, nodes, wts, nut_evals)
pMean0      = cp.E(nut_app, Cs_dist)
pMeanSigma  = cp.Std(nut_app, Cs_dist)

# <codecell>
plt.figure(figsize=(8, 6))

for k in range(1,p_n+1):
    plt.plot(y,np.abs(cp.E(U_app*orthpoly[k], Cs_dist)), label='#node = %s' % k)
plt.legend()
plt.grid(True)
plt.savefig((data + "modesU.png"), dpi=600)

# <codecell>
Sum   = np.zeros(p_n+1)
cells = len(y)-1
delta_y = y[1:cells] - y[0:cells-1]
for k in range(p_n+1):
    for i in range(cells-1):
        Sum[k] += delta_y[i] * cp.E(U_app*orthpoly[k], Cs_dist)[i+1] 
    Sum[k] = Sum[k]/np.max(y)
    print("U_avg[",k,"] = ",Sum[k])

# <codecell>
Sum   = np.zeros(p_n+1)
cells = len(y)-1
delta_y = y[1:cells] - y[0:cells-1]
for k in range(p_n+1):
    for i in range(cells-1):
        Sum[k] += delta_y[i] * cp.E(p_app*orthpoly[k], Cs_dist)[i+1] 
    Sum[k] = Sum[k]/np.max(y)
    print("p_avg[",k,"] = ",Sum[k])
    
# <codecell>
N     = 2
    
Mean  = nutMean0*10**5
Sigma = nutMeanSigma*10**5

plt.figure(figsize=(8, 6))

plt.plot(y,Mean,linewidth=0.5)
plt.plot(y,Mean+N*Sigma,'k',linewidth=0.5)
plt.plot(y,Mean-N*Sigma,'k',linewidth=0.5)

plt.ylabel(r'$\nu_t \ x \ 1e5$')
plt.xlabel('$y/\delta$')
plt.axis([0, 1, 0, 2])
plt.grid(True)
plt.savefig((data + "nutMean_vs_y.png"), dpi=600)

# <codecell>
Mean  = pMean0*10**5
Sigma = pMeanSigma*10**5

plt.figure(figsize=(8, 6))

plt.plot(y,Mean,linewidth=0.5)
plt.plot(y,Mean+N*Sigma,'k',linewidth=0.5)
plt.plot(y,Mean-N*Sigma,'k',linewidth=0.5)

plt.ylabel(r'$p \ x \ 1e5$')
plt.xlabel('$y/\delta$')
#plt.axis([0, 1, 0, 5])
plt.grid(True)
plt.savefig((data + "pMean_vs_y.png"), dpi=600)

# <codecell>
utTmp = yPl_evals*nu/y[1]
ut    = np.mean(utTmp)
# ut    = 0.0072#0.9*nu/y[1]

Mean  = UMean0/ut
Sigma = UMeanSigma/ut

plt.figure(figsize=(8, 6))

plt.plot(y,Mean,linewidth=0.5)
plt.plot(y,Mean+N*Sigma,'k',linewidth=0.5)
plt.plot(y,Mean-N*Sigma,'k',linewidth=0.5)

plt.ylabel(r'$u/u_{\tau}$')
plt.xlabel('$y/\delta$')
plt.axis([0, 1, 0, 25])
plt.grid(True)
plt.savefig((data + "UMean_vs_y.png"), dpi=600)

# <codecell>
# Histogram
if (0):
    n_pts = 200
    r = np.random.randn(n_pts)
    
    #plt.figure(figsize=(8, 6))
    #plt.ylabel('$PDF$')
    #plt.xlabel('$U$')
    #plt.axis([0, 20, 0, 40])
    #plt.grid(True)
    #for i in range(25):
    #    yloc = i
    #    x = Mean[yloc] + Sigma[yloc] * r
    #    n, bins, patches = plt.hist(x, n_pts, density=1, alpha=0.75)
    #plt.savefig((data + "hist_UMean_vs_y.png"), dpi=600)
    
    plt.figure(figsize=(8, 6))
    plt.ylabel('$PDF$')
    plt.xlabel('$U$')
    plt.axis([14, 15, 0, 50])
    plt.grid(True)
    yloc = 8
    x = Mean[yloc] + Sigma[yloc] * r
    n, bins, patches = plt.hist(x, n_pts, density=1, facecolor='b', alpha=0.75)
    plt.savefig((data + "hist1_UMean_vs_y.png"), dpi=600)
    
    
    plt.figure(figsize=(8, 6))
    plt.ylabel('$PDF$')
    plt.xlabel('$U$')
    plt.axis([18.5, 19.5, 0, 50])
    plt.grid(True)
    yloc = 19
    x = Mean[yloc] + Sigma[yloc] * r
    n, bins, patches = plt.hist(x, n_pts, density=1, facecolor='b', alpha=0.75)
    plt.savefig((data + "hist2_UMean_vs_y.png"), dpi=600)
    
    plt.figure(figsize=(8, 6))
    plt.ylabel('$PDF$')
    plt.xlabel('$U$')
    plt.axis([19.5, 20.5, 0, 50])
    plt.grid(True)
    yloc = 25
    x = Mean[yloc] + Sigma[yloc] * r
    n, bins, patches = plt.hist(x, n_pts, density=1, facecolor='b', alpha=0.75)
    plt.savefig((data + "hist3_UMean_vs_y.png"), dpi=600)


# <codecell>
utTmp = yPl_evals*nu/y[1]
#utTmp = 1.8*nu/y[1]
ut    = np.mean(utTmp)
ut2   = ut**2

Mean  = RMean0/ut2
Sigma = RMeanSigma/ut2

plt.figure(figsize=(8, 6))

offset = 1
i = 1
plt.plot(y,np.sqrt(Mean[i]),linewidth=0.5)
plt.plot(y,np.sqrt(Mean[i]+N*Sigma[i]),'k',linewidth=0.5)
plt.plot(y,np.sqrt(Mean[i]-N*Sigma[i]),'k',linewidth=0.5)

i = 2
plt.plot(y,np.sqrt(Mean[i])+offset*1,linewidth=0.5)
plt.plot(y,np.sqrt(Mean[i]+N*Sigma[i])+offset*1,'k',linewidth=0.5)
plt.plot(y,np.sqrt(Mean[i]-N*Sigma[i])+offset*1,'k',linewidth=0.5)

i = 0
plt.plot(y,np.sqrt(Mean[i])+offset*2,linewidth=0.5)
plt.plot(y,np.sqrt(Mean[i]+N*Sigma[i])+offset*2,'k',linewidth=0.5)
plt.plot(y,np.sqrt(Mean[i]-N*Sigma[i])+offset*2,'k',linewidth=0.5)

plt.ylabel('$<u_{rms}>$')
plt.xlabel('$y/\delta$')
plt.axis([0, 1, 0, 6])
plt.grid(True)
plt.savefig((data + "RMean_vs_y.png"), dpi=600)


# <codecell>
plt.figure(figsize=(8, 6))

i = 3
plt.plot(y,(-Mean[i]),linewidth=0.5)
plt.plot(y,(-Mean[i]+N*Sigma[i]),'k',linewidth=0.5)
plt.plot(y,(-Mean[i]-N*Sigma[i]),'k',linewidth=0.5)

plt.ylabel('$-<uv>$')
plt.xlabel('$y/\delta$')
plt.axis([0, 1, 0, 1])
plt.grid(True)
plt.savefig((data + "uvMean_vs_y.png"), dpi=600)