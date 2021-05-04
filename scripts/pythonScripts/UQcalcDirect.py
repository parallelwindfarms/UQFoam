#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Computing the coefficients Cijk ...

"""
# <codecell>

# import libraries
import numpy as np
import itertools as itr
import chaospy as cp
import re

# <codecell>
### Read in UQ parameters and find Mijk, Mijkl products
print ("\nReading UQ paramenters from constant/transportProperties ...\n")
### Read in - order, dim , distType(uniform, normal ... )
UQtransportProperties = open("constant/transportProperties","r")
transportProperties  = UQtransportProperties.read()
UQtransportProperties.close()

print ("\nReading UQ paramenters from constant/turbulenceProperties ...\n")
### Read in - order, dim , distType(uniform, normal ... )
UQturbulenceProperties = open("constant/turbulenceProperties","r")
turbulenceProperties  = UQturbulenceProperties.read()
UQturbulenceProperties.close()

order 		= re.search(r'order.*?;', transportProperties , re.DOTALL).group()
order 		= int(order.strip('order').strip(";").strip())
  
dim   		= re.search(r'dim.*?;', transportProperties , re.DOTALL).group()
dim   		= int(dim.strip('dim').strip(";").strip())

print("Order = ",order,"\nDim = ",dim,"\n")

### Calculate number of terms in PC Expansion (P+1)
P 			= int(np.math.factorial(order+dim) \
  			/ (np.math.factorial(order)*np.math.factorial(dim)) - 1)
	
#nuDist = cp.Normal(2e-5,2e-7)
#CsDist = cp.Uniform(0.075,0.125) # Channel Flow
#CsDist = cp.Uniform(0.050,0.080) # Flow around a cylinder CsDist
CsDist = cp.Uniform(0.050,0.100) # Flow around a cylinder CsDist
#CsDist = cp.Uniform(0.540,0.640) # Flow around a cylinder UinDist

#dist = cp.J(nuDist, CsDist)
#dist = nuDist
dist = CsDist
print(dist)

# <codecell>
print ("Calculating the triple, quad products of basis functions ...\n")

# Determining ortho-gonal/normal polynomials w.r.t. distType
orthoNormal = True; # False;
orthpoly, norms 	= cp.orth_ttr(order, dist, retall=True, normed = orthoNormal)

# Nomalized inner products of the polynomials
phi  = orthpoly    
phi2 = cp.outer(phi, phi)
phi3 = cp.outer(phi, phi, phi)
CsDist, nuDist = cp.variable(2)
	
Ckk  = cp.E(phi2, dist)
Cijk = cp.E(phi3, dist)
Djk  = cp.E(nuDist * phi2, dist)
tjk  = cp.E(CsDist*CsDist * phi2, dist)
Tijk  = cp.E(CsDist*CsDist * phi3, dist)

# <codecell>
### Write out Cikj, Cikjl on OpenFOAM format
print ("Writing Mijk, Mijkl into constant/gPCcoeffs ...\n")

Pplus1 = P + 1
SMALL  = 1e-12

UQProperties = open("constant/gPCcoeffs","w+")

UQProperties.write("\
/*--------------------------------*- C++ -*----------------------------------*\\\n\
| =========                 |                                                 |\n\
| \\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n\
|  \\\    /   O peration     | Version:  v1806                                 |\n\
|   \\\  /    A nd           | Web:      www.OpenFOAM.com                      |\n\
|    \\\/     M anipulation  |                                                 |\n\
\*---------------------------------------------------------------------------*/\n\
FoamFile\n\
{\n\
    version     2.0;\n\
    format      ascii;\n\
    class       dictionary;\n\
    location    \"constant\";\n\
    object      gPCcoeffs;\n\
}\n\
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n")

UQProperties.write("// Quad product of basis function k with itself \n\n")
for i,j in itr.product(*[range(Pplus1)]*2):
    if (Ckk[i,j]!=0.0 and i==j):
        UQProperties.write("M["+str(i)+"]\t" + str(Ckk[i,j])+";\n")
        	        
UQProperties.write("\n\n// Triple product of basis functions\n\n")
for i,j,k in itr.product(*[range(Pplus1)]*3):
    if (np.abs(Cijk[i,j,k])>=SMALL):
        UQProperties.write("M["+str(i)+"]["+str(j)+"]["+str(k)+"]\t" \
                              +str(Cijk[i,j,k])+";\n")
	
UQProperties.write("\n\n// <Cs*Cs*phi*phi> product of basis functions\n\n")
for j,k in itr.product(*[range(Pplus1)]*2):
    if (np.abs(tjk[j,k])>=SMALL):
        UQProperties.write("tjk["+str(j)+"]["+str(k)+"]\t" \
                              +str(tjk[j,k])+";\n")
                                            	
UQProperties.write("\n\n// <Cs*Cs*phi*phi*phi> product of basis functions\n\n")
for i,j,k in itr.product(*[range(Pplus1)]*3):
    if (np.abs(Tijk[i,j,k])>=SMALL):
        UQProperties.write("Tijk["+str(i)+"]["+str(j)+"]["+str(k)+"]\t" \
                              +str(Tijk[i,j,k])+";\n")

UQProperties.close()

print ("Done.\n")