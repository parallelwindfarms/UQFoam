#!/usr/bin/env python3

import numpy as np
import chaospy as cp

nu = 0.1
sigma = 0.01
G = cp.Normal(nu, sigma)
g = cp.Normal(0,sigma)

def model_solver(q):
	return [G.pdf(q)*(np.sqrt(4*np.pi)*sigma)]
	#return [np.sqrt(2)*G.pdf(q)/G.pdf(0.1)]

N = 4
K = 10

orthpoly   = cp.orth_ttr(N, g, normed =False)
nodes, wts = cp.generate_quadrature(K, G, rule='G')
Cs = model_solver(nodes)
polyvals = orthpoly(nodes)
C = np.zeros(N+1)

for i in range(N+1):
	C[i] = sum(wts*polyvals[i][0]*Cs[0][0])
	print(C[i])
	
#print(nodes[0])
#print(wts)
#print(polyvals[1][0])
#print(Cs[0][0])

