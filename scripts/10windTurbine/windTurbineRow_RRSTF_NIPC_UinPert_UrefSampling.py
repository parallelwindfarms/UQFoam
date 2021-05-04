#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 29 15:21:21 2021

@author: jigar
"""

# <codecell> Import Required Modules and env vars
import os, sys

import numpy as np
import chaospy as cp
import matplotlib.pyplot as plt

from multiprocessing import Pool
from functools import partial
from timeit import default_timer as timer

SCRIPTS = os.environ['SCRIPTS']
DATA    = os.environ['windTurbineData']

sys.path.append(SCRIPTS+'/10windTurbine')
sys.path.append(SCRIPTS+'/pythonScripts')

cwd = os.getcwd() + "/"

import myUQlib

# <codecell> Sampling and Perturbing
# Sampling 3 r.v.s s.t. eVal1* + eVal2* + eVal3* = 0
nSamples = 30

UinDist = cp.Uniform(7.5, 8.5)

np.random.seed(1)
UinSamples = UinDist.sample(nSamples, rule='L')
UinSamples.sort()
print(f"UinSamples min = {UinSamples.min()} , max = {UinSamples.max()}")

UinSamples = np.vstack((UinSamples,UinSamples,UinSamples))
UinSamples = UinSamples.ravel()

nSamples = nSamples*3

# <codecell> Perturbed samples of inlet reference velocity Uin_Smaples
def get_Uin_PertSample(Uin_Sample, s):
    print(s)
    myUQlib.writeScalarValue(s, 'Uin_Samples/', 'Uref', Uin_Sample)

start = timer()
for s in range(nSamples):
    get_Uin_PertSample(UinSamples[s], s)
print(timer()-start, 's')

