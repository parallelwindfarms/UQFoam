#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Mapping a uniform distribution U(a,b)
to a standard uniform distribution U(-1,1)

"""

import chaospy
import numpy as np

# <codecell>
a = 0.075
b = 0.125
mu    = a+(b-a)/2
sigma = (b-a)/2/np.sqrt(3)

if(1):
    real_distribution = chaospy.Uniform(a, b)
    proxy_distribution = chaospy.Uniform(-1, 1)

if(0):
	real_distribution = chaospy.Normal(mu, sigma)
	proxy_distribution = chaospy.Normal(0, 1)

orthoNormal = True
#orthoNormal = False
polynomial_expansion = chaospy.orth_ttr(2, proxy_distribution, normed=orthoNormal)

proxy_X, proxy_W = chaospy.generate_quadrature(4, proxy_distribution)
real_X = real_distribution.inv(proxy_distribution.fwd(proxy_X))
real_W = proxy_W * real_distribution.pdf(real_X) / proxy_distribution.pdf(proxy_X)

evals = [(sample) for sample in real_X.T]

foo_approx, fourier_coeffs = chaospy.fit_quadrature(
    polynomial_expansion, proxy_X, real_W, evals, retall=True)

#print(foo_approx)    
print(fourier_coeffs)

