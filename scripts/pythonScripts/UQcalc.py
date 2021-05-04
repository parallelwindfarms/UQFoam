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

# Function to calculate the triple product of basis functins
def calcM(P, poly, distr):
    """
    Calculates product Mkk, Mijk, Mijkl, Mijklm of the 
    orthonomal polynomials based on the distribution 'distr'
    """
    
    phi    = poly    
    phi2   = cp.outer(phi, phi)
    phi3   = cp.outer(phi, phi, phi)
    phi4   = cp.outer(phi, phi, phi, phi)
    phi5   = cp.outer(phi, phi, phi, phi, phi)
    	
    Mkk    = cp.E(phi2, distr)
    Mijk   = cp.E(phi3, distr)
    Mijkl  = cp.E(phi4, distr)
    Mijklm = cp.E(phi5, distr)
       
    return Mkk, Mijk, Mijkl, Mijklm

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

nuDistName 	= re.search(r'distType.*?;', transportProperties , re.DOTALL).group()
nuDistName 	= nuDistName.strip('distType').strip(";").strip()

CsDistName 	= re.search(r'distType.*?;', turbulenceProperties , re.DOTALL).group()
CsDistName 	= CsDistName.strip('distType').strip(";").strip()

print("Order = ",order,"\nDim = ",dim,"\nDistribution type : nu -> ",nuDistName,"\nDistribution type : Cs -> ",CsDistName,"\n")

### Calculate number of terms in PC Expansion (P+1)
P 			= int(np.math.factorial(order+dim) \
  			/ (np.math.factorial(order)*np.math.factorial(dim)) - 1)

nu			= np.zeros(P+1)
Cs			= np.zeros(P+1)

### Finding the joint distribution
for i in range(P+1):
	start 	= "nu"+str(i)+" [0 2 -1 0 0 0 0]"
	if re.search(start, transportProperties):
		tmp		= re.search(r'nu'+str(i)+'.*?;', transportProperties, re.DOTALL).group()
		nu[i]	= tmp[tmp.find(start)+len(start):tmp.rfind(";")]
	start 	= "Cs"+str(i)
	if re.search(start, turbulenceProperties):
		tmp		= re.search(r'Cs'+str(i)+'.*?;', turbulenceProperties, re.DOTALL).group()
		Cs[i]	= tmp[tmp.find(start)+len(start):tmp.rfind(";")]
	
if nuDistName == 'normal':
	nuDistType = cp.Normal(0,1)
if nuDistName == 'uniform':
	nuDistType = cp.Uniform(-1,1)
	
if CsDistName == 'uniform':
	CsDistType = cp.Uniform(-1,1)
	#CsDistType = cp.Uniform(0.075,0.125)
if CsDistName == 'normal':
	CsDistType = cp.Normal(0,1)
	
#distType = cp.J(nuDistType, CsDistType)
distType = nuDistType
# distType = CsDistType
print(distType)

# <codecell>
print ("Calculating the triple, quad products of basis functions ...\n")

# Determining ortho-gonal/normal polynomials w.r.t. distType
orthoNormal = True; # False;
orthpoly, norms 	= cp.orth_ttr(order, distType, retall=True, normed = orthoNormal)

# Nomalized inner products of the polynomials
Ckk, Cikj, Cijkl, Cijklm = calcM(P, orthpoly, distType)

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
    if (np.abs(Cikj[i,j,k])>=SMALL):
        UQProperties.write("M["+str(i)+"]["+str(j)+"]["+str(k)+"]\t" \
                              +str(Cikj[i,j,k])+";\n")
	
UQProperties.write("\n\n// Bi-quad product of basis functions\n\n")
for i,j,k,l in itr.product(*[range(Pplus1)]*4):
    if (np.abs(Cijkl[i,j,k,l])>=SMALL):
        UQProperties.write("M["+str(i)+"]["+str(j)+"]["+str(k)+"]["+str(l)+"]\t" \
                              +str(Cijkl[i,j,k,l])+";\n")
                                            	
UQProperties.write("\n\n// Order 5 product of basis functions\n\n")
for i,j,k,l,m in itr.product(*[range(Pplus1)]*5):
    if (np.abs(Cijklm[i,j,k,l,m])>=SMALL):
        UQProperties.write("M["+str(i)+"]["+str(j)+"]["+str(k)+"]["+str(l)+"]["+str(m)+"]\t" \
                              +str(Cijklm[i,j,k,l,m])+";\n")

UQProperties.close()

print ("Done.\n")


