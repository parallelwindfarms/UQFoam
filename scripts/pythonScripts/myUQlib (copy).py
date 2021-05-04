#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 20 18:05:20 2020

@author: jigar
"""

import numpy as np
import itertools as itr
import chaospy as cp
import pyvista as pv

# <codecell> 
"""
Read OpenFOAM case mesh acces to grid and cell values
"""
def getMeshInfo(cwd, vtkFile):
  
    mesh    = pv.UnstructuredGrid(cwd + vtkFile)
    nCells  = mesh.n_cells
    sized   = mesh.compute_cell_sizes()
    cVols   = sized.cell_arrays["Volume"]
    mVol    = mesh.volume
    volAvg  = cVols/mVol
    
    return mesh, nCells, volAvg


# <codecell> 
"""
getSmaple funtion for Gausaain process
"""
def getGpSample(kleModes, quadNode, dim, numCells, sigma):
    tmpSample = np.zeros(numCells)
    for j in range(dim):
        tmpSample += kleModes[j] * quadNode[j]
    return  sigma*tmpSample


# <codecell> 
"""
getSmaple funtion for a field (any) - cell value
"""
def getGpSampleCellVal(kleModes, quadNode, dim, sigma):
    tmpSample = 0
    for j in range(dim):
        tmpSample += kleModes[j] * quadNode[j]
    return  sigma*tmpSample


# <codecell> 
"""
getSmaple funtion for logN field
"""
def getExpGpSample(kleModes, quadNode, dim, numCells, sigma):
    tmpSample = np.zeros(numCells)
    for j in range(dim):
        tmpSample += (kleModes[j] * quadNode[j])
    return  np.exp(sigma*tmpSample)

    
# <codecell> Coeffs of PCE of Gamma dist wrt Normal dist
def getPCEGamma(u=cp.Gamma(1), n_u=5):
    
    xi = cp.Normal(0,1)
    
    nrDiv = 1000
    min_w = -10
    max_w = 10
    dw    = (max_w-min_w)/nrDiv
    w     = np.linspace(min_w, max_w, nrDiv+1)
    w_c   = (w[1:] + w[0:-1])/2
    dw    = w[1:] - w[0:-1]
    
    herm,nrms = cp.orth_ttr(n_u,xi,retall=True,normed=True)
    Np = len(herm)
    
    Ub = np.zeros(Np)
    
    for Np_i in range(Np):
        Ub[Np_i] = sum(u.inv(xi.cdf(w_c)) * herm[Np_i](w_c) * xi.pdf(w_c) * dw)

    return Ub


# <codecell> 
"""
Function for setting UQ info. for OpenFOAM simulation
"""
def uqInfo(n, d, q, RF):
              
    dist      = cp.Iid(cp.Normal(0,1), d)
    phi,norms = cp.orth_ttr(n, dist, retall=True, normed=True, cross_truncation=1/q)
    expoMat   = cp.bertran.bindex(start=0, stop=n, dim=d, cross_truncation=1/q)
    P         = len(expoMat)-1
    Pplus1    = P+1
    
    # Nomalized outer and inner products of the polynomials
    phi3 = cp.outer(phi, phi, phi)
    	
    Ckk  = cp.E(phi*phi, dist)
    Cijk = cp.E(phi3, dist)
                    
    # Write out Cikj, Cikjl in OpenFOAM format
    print ("Writing Mijk, Mijkl into constant/gPCcoeffs ...\n")
    
    SMALL  = 1e-12
    
    UQProperties = open("constant/uqInfo","w+")

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
    object      uqInfo;\n\
}\n\
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n")

    UQProperties.write("// Type of random field\n\n")
    UQProperties.write("randField\t"+RF+";\n\n")
    
    UQProperties.write("// Number of dimensions\n\n")
    UQProperties.write("d\t"+str(d)+";\n\n")

    UQProperties.write("// Number of modes\n\n")
    UQProperties.write("P\t"+str(P)+";\n\n")
       
    expoMat = np.array(expoMat)
    UQProperties.write("\n\n// Exponent matrix basis functions\n\n")
    for i in range(Pplus1):
        for j in range(d):
            if (expoMat[i,j]!=0):
                UQProperties.write("expoMat["+str(i)+"]["+str(j)+"]\t" \
                                       +str(expoMat[i,j])+";\n")
                    
    UQProperties.write("\n\n// Quad product of basis function k with itself \n\n")
    for i in range(Pplus1):
        UQProperties.write("M["+str(i)+"]\t" + str(Ckk[i])+";\n")
    
    UQProperties.write("\n\n// Triple product of basis functions\n\n")
    for i,j,k in itr.product(*[range(Pplus1)]*3):
        if (np.abs(Cijk[i,j,k])>=SMALL):
            UQProperties.write("M["+str(i)+"]["+str(j)+"]["+str(k)+"]\t" \
                                  +str(Cijk[i,j,k])+";\n")
    
    UQProperties.close()
    
    print ("Done.\n")


# <codecell> 
"""
PCE of logN field (Ghanem 2002) with vol avgd solbol indices
"""
def getlogNFieldPCE(nutBl,stD,nCells,n,d,q,g,volAvg,retained,phyDim):
    
    expoMat      = cp.bertran.bindex(start=0,stop=n,dim=d,cross_truncation=1/q)
    P            = len(expoMat)-1
    Pplus1       = P+1
    
    print('\nPplus1 = ', Pplus1)
    
    nut0         = np.exp(np.log(nutBl) + (1/2)*np.square(stD))
    
    nutVar       = np.zeros((1, nCells))
    nutVar_i     = np.zeros(((d, 1, nCells)))
    nutVar_i_tot = np.zeros(((d, 1, nCells)))
    
    for i in range(Pplus1):
        
        gPdt_i = 1
        p_i    = sum(expoMat[i])
        
        for j in range(d):
            gPdt_i *= np.power(g[j], expoMat[i][j])
            
        c_i = np.power(stD,i) * gPdt_i
    
        nut_i = c_i*nut0 
        
        if i>0:
            nutVar += np.square(nut_i).T
            
            for j in range(d):
                if expoMat[i][j]>0:
                    nutVar_i_tot[j] += np.square(nut_i).T
                    if expoMat[i][j]==p_i:
                        nutVar_i[j] += np.square(nut_i).T
        
        # if retained==True:
            # writeNutGPCModes(i,nut_i,nCells, "./nut_gPCKLE_LogNProc/", phyDim)
        
    sblIdx_i = (nutVar_i/nutVar * volAvg).sum(2) * 100
    
    # Printing volume averaged sobol indices
    if retained==False:
        gIdx = []
    for i in range(d):
        if sblIdx_i[i]>1:
            print(i+1,' ', sblIdx_i[i], '%')
            if retained==False:
                gIdx.append(i)
    
    if retained==False:
        return gIdx
        
# <codecell> 
"""
writing gPC-KLE modes of lognormal nut
"""
def writeNutPCEModes(i, nutMode, nCells, dirName, phyDim):
    outputFile = open(dirName+"nut"+str(i),"w+")
    
    outputFile.write('''\
/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\    /   O peration     | Version:  v1806                                 |
|   \\\  /    A nd           | Web:      www.OpenFOAM.com                      |
|    \\\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
version     2.0;
format      ascii;
class       volScalarField;
location    "0";
object      nut'''+str(i)+''';
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
dimensions      [0 2 -1 0 0 0 0];

internalField   nonuniform List<scalar> \n'''+

    str(nCells)+'''\n(\n''')

    for j in range(nCells):
        outputFile.write(str(nutMode[j][0]) + "\n")
    
    outputFile.write('''\
)
;

boundaryField
{
    "(inlet|outlet)"
    {
        type            cyclic;
    }''')
    
    if phyDim==2 :
        outputFile.write('''
        "(front|back)"
        {
            type            empty;
        }''')
        
    if phyDim==3 :
        outputFile.write('''
        "(front|back)"
        {
            type            cyclic;
        }''')
    
    outputFile.write('''
    "(top|hills)"
    {
        type            nutkWallFunction;
        value           uniform 0;
    }
}

''')
    outputFile.close()
    
# <codecell> 
"""
writing gPC-KLE modes of lognormal k
"""
def writeKGPCModes(i, kMode, nCells, dirName, phyDim):
    outputFile = open(dirName+"k"+str(i),"w+")
    
    outputFile.write('''\
/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\    /   O peration     | Version:  v1806                                 |
|   \\\  /    A nd           | Web:      www.OpenFOAM.com                      |
|    \\\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
version     2.0;
format      ascii;
class       volScalarField;
location    "0";
object      k'''+str(i)+''';
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
dimensions      [0 2 -2 0 0 0 0];

internalField   nonuniform List<scalar> \n'''+

    str(nCells)+'''\n(\n''')

    for j in range(nCells):
        outputFile.write(str(kMode[j][0]) + "\n")
    
    outputFile.write('''\
)
;

boundaryField
{
    "(inlet|outlet)"
    {
        type            cyclic;
    }''')
    
    if phyDim==2 :
        outputFile.write('''
        "(front|back)"
        {
            type            empty;
        }''')
        
    if phyDim==3 :
        outputFile.write('''
        "(front|back)"
        {
            type            cyclic;
        }''')
    
    outputFile.write('''
    "(top|hills)"
    {
        type            zeroGradient;
        value           uniform 0;
    }
}

''')
    outputFile.close()

# <codecell> 
"""
writing gPC coeffients 
"""
def writePCECoeffs(i, nutCoeffs, nCells, dirName, phyDim):
    outputFile = open(dirName+"nutCoeffs"+str(i),"w+")
    
    outputFile.write('''\
/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\    /   O peration     | Version:  v1806                                 |
|   \\\  /    A nd           | Web:      www.OpenFOAM.com                      |
|    \\\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
version     2.0;
format      ascii;
class       volScalarField;
location    "0";
object      nutCoeffs'''+str(i)+''';
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
dimensions      [0 0 0 0 0 0 0];

internalField   nonuniform List<scalar> \n'''+

    str(nCells)+'''\n(\n''')

    for j in range(nCells):
        outputFile.write(str(nutCoeffs[j][0]) + "\n")
    
    outputFile.write('''\
)
;

boundaryField
{
    "(inlet|outlet)"
    {
        type            cyclic;
    }''')
    
    if phyDim==2 :
        outputFile.write('''
        "(front|back)"
        {
            type            empty;
        }''')
        
    if phyDim==3 :
        outputFile.write('''
        "(front|back)"
        {
            type            cyclic;
        }''')
    
    outputFile.write('''
    "(top|hills)"
    {
        type            fixedValue;
        value           uniform 0;
    }
}

''')
    outputFile.close()
    
# <codecell> 
"""
writing gPC-KLE modes of lognormal nut
"""
def writeGmatPCEModes(i, GMode, nCells, dirName, phyDim):
    outputFile = open(dirName+"GCoeffs"+str(i),"w+")
    
    outputFile.write('''\
/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\    /   O peration     | Version:  v1806                                 |
|   \\\  /    A nd           | Web:      www.OpenFOAM.com                      |
|    \\\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
version     2.0;
format      ascii;
class       volTensorField;
location    "0";
object      GCoeffs'''+str(i)+''';
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
dimensions      [0 0 0 0 0 0 0];

internalField   nonuniform List<tensor> \n'''+

    str(nCells)+'''\n(\n''')

    for j in range(nCells):
        outputFile.write('''( ''')
        for p in range(3):
            for q in range(3):
                outputFile.write(str(GMode[j,p,q]) + " ")
        outputFile.write(''')\n''')
    
    outputFile.write('''\
)
;

boundaryField
{
    "(inlet|outlet)"
    {
        type            zeroGradient;
    }''')
    
    if phyDim==2 :
        outputFile.write('''
        "(front|back)"
        {
            type            empty;
        }''')
        
    if phyDim==3 :
        outputFile.write('''
        "(front|back)"
        {
            type            cyclic;
        }''')
    
    outputFile.write('''
    "(top|hills)"
    {
        type            zeroGradient;
    }
}

''')
    outputFile.close()