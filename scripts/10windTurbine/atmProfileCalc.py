#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 18:37:35 2021

@author: jigar
"""

import numpy as np


# <codecell> Constants
kappa = 0.4
Cmu   = 0.09
Uref  = 8
Zref  = 70
d     = 0


# <codecell> Turbulent quantities
I = 15/100

num = kappa * (2/3)**0.5
den = I* Cmu**0.25

Z0 = Zref/np.exp(num/den)
uf = kappa*Uref/np.log((Zref+Z0)/Z0)
k  = uf**2 / Cmu**0.5
e  = uf**3 / (kappa*(Zref+Z0))
w  = uf / (kappa* Cmu**0.5) / (Zref+Z0)
nut = Cmu * (k**2) / e
R11 = 2*k/3

print(f"I  {I} \nZ0 {Z0}\nuf {uf}\nk  {k} \
      \ne  {e}\nw  {w}\nnut  {nut}\nR11  {R11}")