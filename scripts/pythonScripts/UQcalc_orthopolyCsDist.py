#!/usr/bin/env python3
# -*- coding: utf-8 -*-

###############################################################################

# import libraries
import numpy as np
import chaospy as cp
import re

###############################################################################

# Function to calculate the triple product of basis functins
def calcM(P, poly, distr):
    """
    Calculates triple product Mijk of the orthonomal polynomials
    based on the distribution 'distr'
    """
	
    SMALL	= 1e-12
    Pplus1 	= P + 1
    Mkk 	= np.zeros(Pplus1)
    Mijk 	= np.zeros(((Pplus1,Pplus1,Pplus1)))
    Mijkl 	= np.zeros((((Pplus1,Pplus1,Pplus1,Pplus1))))
    Mijklm	= np.zeros(((((Pplus1,Pplus1,Pplus1,Pplus1,Pplus1)))))
    Mlmpqk	= np.zeros(((((Pplus1,Pplus1,Pplus1,Pplus1,Pplus1)))))

    for k in range(Pplus1):
    	Mkk[k] = cp.E(poly[k]*poly[k], distr)

    for i in range(Pplus1):
        for j in range(Pplus1):
            for k in range(Pplus1):
                #if ((i+j+k)%2 != 0):
                #    Mijk[i,j,k] = 0.0
                #else :
                Mijk[i,j,k] = cp.E(poly[i]*poly[j]*poly[k],distr)# + 1e-12
                if (np.abs(Mijk[i,j,k])<SMALL):
                    Mijk[i,j,k] = 0.0
                
                for l in range(Pplus1):      
                        #if ((i+j+k+l)%2 != 0):
                        #	Mijkl[i,j,k,l] = 0.0
                        #else :
                        Mijkl[i,j,k,l] = cp.E(poly[i]*poly[j]*poly[k]*poly[l],distr)# + 1e-12
                        if (np.abs(Mijkl[i,j,k,l])<SMALL):
                        	Mijkl[i,j,k,l] = 0.0
                        
                        for m in range(Pplus1):
                        	#if ((i+j+k+l+m)%2 != 0):
                        	#	Mijklm[i,j,k,l,m] = 0.0
                        	#else :
                        	Mijklm[i,j,k,l,m] = cp.E(poly[i]*poly[j]*poly[k]*poly[l]*poly[m],distr)# + 1e-12
                        	if (np.abs(Mijklm[i,j,k,l,m])<SMALL):
                        		Mijklm[i,j,k,l,m] = 0.0
       
    return Mkk, Mijk, Mijkl, Mijklm, Mlmpqk

###############################################################################

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

# Calculate number of terms in PC Expansion (P+1)
P 			= int(np.math.factorial(order+dim) \
  			/ (np.math.factorial(order)*np.math.factorial(dim)) - 1)

nu			= np.zeros(P+1)
Cs			= np.zeros(P+1)

for i in range(P+1):
	start 	= "nu"+str(i)+" [0 2 -1 0 0 0 0]"
	if re.search(start, transportProperties):
		tmp		= re.search(r'nu'+str(i)+'.*?;', transportProperties, re.DOTALL).group()
		nu[i]	= tmp[tmp.find(start)+len(start):tmp.rfind(";")]
	start 	= "Cs"+str(i)
	if re.search(start, turbulenceProperties):
		tmp		= re.search(r'Cs'+str(i)+'.*?;', turbulenceProperties, re.DOTALL).group()
		Cs[i]	= tmp[tmp.find(start)+len(start):tmp.rfind(";")]

#tmp		= re.search(r'nu_s.*?;', transportProperties , re.DOTALL).group()
#nu_s	  	= float(tmp.strip('nu_s').strip("nu_s [0 2 -1 0 0 0 0]").strip(";"))
	
if nuDistName == 'normal':
	nuDistType = cp.Normal(0,1)
if CsDistName == 'normal':
	CsDistType = cp.Uniform(0.075,0.125)	
	#CsDistType = cp.Uniform(-1,1)
	#CsDistType = cp.Normal(0,1)

#distType = cp.J(nuDistType, CsDistType)
distType = CsDistType
#distType = nuDistType
#print(distType)

###############################################################################

print ("Calculating the triple, quad products of basis functions ...\n")

# Determining ortho-gonal/normal polynomials w.r.t. distType
orthoNormal = True;
#orthoNormal = False;
orthpoly 	= cp.orth_ttr(order, distType, normed = orthoNormal)
#print(orthpoly)

# Nomalized Inner product of the triple product of polynomials
Ckk, Cikj, Cijkl, Cijklm, Clmpqk	= calcM(P, orthpoly, distType)

#print ("Evaluating the polynomials at nu_s ...\n")

###############################################################################

### Write out Cikj, Cikjl on OpenFOAM format
print ("Writing Mijk, Mijkl into constant/gPCcoeffs ...\n")

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

Pplus1 = P + 1

UQProperties.write("// Quad product of basis function k with itself \n\n")
for k in range(Pplus1):
	if (orthoNormal==True):
		UQProperties.write("M["+str(k)+"]\t" + str(1.0)+";\n")
	else :
		UQProperties.write("M["+str(k)+"]\t" + str(Ckk[k])+";\n")
        	        
UQProperties.write("\n\n// Triple product of basis functions\n\n")
for i in range(Pplus1):
	for j in range(Pplus1):
		for k in range(Pplus1):
			if (Cikj[i,j,k]!=0.0):
	        	    UQProperties.write("M["+str(i)+"]["+str(j)+"]["+str(k)+"]\t" + str(Cikj[i,j,k])+";\n")
	
UQProperties.write("\n\n// Bi-quad product of basis functions\n\n")
for i in range(Pplus1):
	for j in range(Pplus1):
		for k in range(Pplus1):
			for l in range(Pplus1):
				if (Cijkl[i,j,k,l]!=0.0):
                                	UQProperties.write("M["+str(i)+"]["+str(j)+"]["+str(k)+"]["+str(l)+"]\t" \
                                        	    	+str(Cijkl[i,j,k,l])+";\n")
                                            	
UQProperties.write("\n\n// Order 5 product of basis functions\n\n")
for i in range(Pplus1):
	for j in range(Pplus1):
		for k in range(Pplus1):
			for l in range(Pplus1):
				for m in range(Pplus1):
					if (Cijklm[i,j,k,l,m]!=0.0):
		                        	UQProperties.write("M["+str(i)+"]["+str(j)+"]["+str(k)+"]["+str(l)+"]["+str(m)+"]\t" \
		                        	            	+str(Cijklm[i,j,k,l,m])+";\n")
                                            	                                          
#UQProperties.write("\n\n// Polynomials at nu_s for surrogate evaluation\n\n")

#for i in range(Pplus1):
#	UQProperties.write("psi["+str(i)+"]\t"+str(orthpoly[i](nu_s))+";\n")

UQProperties.close()

print ("Done.\n")


