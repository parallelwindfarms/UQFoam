#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import chaospy as cp
import pandas as pd
import matplotlib as plt

from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

df 	= pd.read_csv("Updf.csv")
cells  	= 100
lines	= 3
tol 	= 1e-3*np.ones((lines,cells))
U0 = np.zeros((lines,cells))
U1 = np.zeros((lines,cells))

for i in range(1,lines+1):
	U0[i-1] = df['U0_0'+str(i)]
	U1[i-1] = df['U1_0'+str(i)]
U0 = U0 + tol
U1 = np.abs(U1) + tol

u_evals = 500
u_min	= 0
u_max 	= 2
D		= 0.2
PDF 	= np.zeros(((lines,cells,u_evals)))
u_range	= np.linspace(u_min, u_max, u_evals)
u		= np.zeros((cells,u_evals))
r		= np.linspace(0.0,D/2.0,cells)/(D/2)
for i in range(lines):
	for j in range(cells):
		#u[j]	  = np.linspace(U0[i][j]-5*U1[i][j], U0[i][j]+5*U1[i][j], u_evals)
		#PDF[i][j] = cp.Normal(U0[i][j],U1[i][j]).pdf(u[j])
		PDF[i][j] = cp.Normal(U0[i][j],U1[i][j]).pdf(u_range)

#print (PDF[1][1])

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm


l 		= 2
Y		= r
X		= u_range
X, Y 	= np.meshgrid(X, Y)
Z		= PDF[l]/np.max(PDF[l])

fsize = 15
fig = plt.figure()
ax 	= plt.axes(projection='3d')
ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.viridis,linewidth=0, antialiased=False)#, rstride=8, cstride=8, alpha=0.3)
ax.set_xlabel("$u(y)$",fontsize=fsize)
ax.set_ylabel("$y/\delta$",fontsize=fsize)
ax.set_zlabel("$PDF(u(y))$",fontsize=fsize);

ax.set_xlim(2,0)
ax.set_ylim(0,1)
ax.set_zlim(0,1)

ax.set_xticks([0, .4, .8, 1.2, 1.6, 2],minor=False)
ax.set_yticks([0, .2, .4, .6, .8, 1],minor=False)
ax.set_zticks([0, .2, .4, .6, .8, 1],minor=False)

ax.view_init(40, -30)

plt.savefig('PDF_xByd_3.png', format='png', dpi=1200)
plt.savefig('PDF_xByd_3.eps', format='eps')

plt.show()

print ("Done")
