# -*- coding: utf-8 -*-
"""

KLE of Gaussian Process
on a coarse mesh to be
mapped on to the GPC 
mesh

"""
# <codecell> Import Required Modules
from __future__ import print_function

import numpy as np
import pyvista as pv
import openturns as ot
import os

cwd = os.getcwd() + "/"

phyDim = 2#2#3

# <codecell> Read OpenFOAM 2D case mesh
# acces to grid and cell values
if (phyDim==2):
    mesh    = pv.UnstructuredGrid(cwd + 'VTK/0KLE_GaussianProcess_0.vtk')
    nCells  = mesh.n_cells
    cCenter = mesh.cell_centers().points
    cBounds = mesh.cell_centers().bounds
    nPoints  = mesh.n_points
    mBounds = mesh.bounds
    mPoints = mesh.points
    # nArrays = mesh.n_arrays
    # mCentre = mesh.center
    
    # Nx = 32; Ny = 24;
    # Nx = 64; Ny = 24; 
    Nx = 128; Ny = 24; 
    # Nx = 256; Ny = 24; 
    
    myIndices = [Nx-1, Ny-1]
    mesher    = ot.IntervalMesher(myIndices)
    kleBounds = ot.Interval([cBounds[0], cBounds[2]], 
                            [cBounds[1], cBounds[3]])
    kleMesh = mesher.build(kleBounds)
    kleMesh.setVertices(cCenter[:,[0,1]])

# <codecell> Read OpenFOAM 3D case mesh
# acces to grid and cell values
if (phyDim==3):
    mesh    = pv.UnstructuredGrid(cwd + 'VTK/0KLE_GaussianProcess_0.vtk')
    nCells  = mesh.n_cells
    cCenter = mesh.cell_centers().points
    cBounds = mesh.cell_centers().bounds
    nPoints  = mesh.n_points
    mBounds = mesh.bounds
    mPoints = mesh.points
    # nArrays = mesh.n_arrays
    # mCentre = mesh.center
    
    # Nx = 16; Ny = 24; Nz = 8;
    Nx = 64; Ny = 24; Nz = 4;
    
    myIndices = [Nx-1, Ny-1, Nz-1]
    mesher    = ot.IntervalMesher(myIndices)
    kleBounds = ot.Interval([cBounds[0], cBounds[2], cBounds[4]], 
                            [cBounds[1], cBounds[3], cBounds[5]])
    kleMesh = mesher.build(kleBounds)
    kleMesh.setVertices(cCenter[:,[0,1,2]])

# <codecell> KLE algorithm
H  = 0.028

lx = 1.5 * H #0.6#1.5#3.0
ly = 0.5 * H #0.2#0.5#1.0
lz = 1.5 * H #1.5

l  = [lx, ly]
if phyDim==3 :
    l.append(lz)

var       = 0.5#0.5#1.0
sigma     = [np.sqrt(var)]
covModel  = ot.SquaredExponential(l, sigma)
eigValK_0 = 1e-2#1e-2

# <codecell> Set KLE Algorithm KarhunenLoeveP1Algorithm
kleAlgo = ot.KarhunenLoeveP1Algorithm(kleMesh, covModel, eigValK_0)

# <codecell> Set KLE Algorithm KarhunenLoeveSVDAlgorithm
# sample = ot.GaussianProcess(covModel, kleMesh).getSample(10)
# kleAlgo = ot.KarhunenLoeveSVDAlgorithm(sample, eigValK_0)

# <codecell> Run KLE Algorithm    
kleAlgo.run()
kleResult = kleAlgo.getResult()  

# <codecell> KLE results
eigVals = kleResult.getEigenValues()    
modes   = kleResult.getModes()
g       = kleResult.getScaledModesAsProcessSample()
xi      = kleResult.project(g)
g       = np.power(g,1)

tmpEvals = np.power(eigVals,1)

# <codecell> Writing scaled KLE modes
for i in range(len(eigVals)):

    if phyDim==2 :
        gFile = open("2D_KLE_M"+str(Nx)+"x"+str(Ny)+"_lx"+str(lx/H)+"H_ly"+ \
                     str(ly/H)+"H_var"+str(var)+"/g"+str(i+1),"w+")    
    if phyDim==3 :
        gFile = open("3D_KLE_M"+str(Nx)+"x"+str(Ny)+"x"+str(Nz)+"_lx"+str(lx/H)+"H_ly"+ \
                     str(ly/H)+"H_lz"+str(lz/H)+"H_var"+str(var)+"/g"+str(i+1),"w+")

    gFile.write('''\
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
object      g'''+str(i+1)+''';
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
dimensions      [0 2 -1 0 0 0 0];

internalField   nonuniform List<scalar> \n'''+

    str(nCells)+'''\n(\n''')

    for j in range(nCells):
        gFile.write(str(g[i][j][0]) + "\n")

    gFile.write('''\
)
;

boundaryField
{
    "(inlet|outlet)"
    {
        type            zeroGradient;
    }''')
    
    if phyDim==2 :
        gFile.write('''
        "(front|back)"
        {
            type            empty;
        }''')
        
    if phyDim==3 :
        gFile.write('''
        "(front|back)"
        {
            type            zeroGradient;
        }''')
    
    gFile.write('''
    "(top|hills)"
    {
        type            zeroGradient;
    }
}
''')
    gFile.close()


# <codecell> Writing scaled KLE modes
for i in range(len(eigVals)):

    gFile = open("2D_mapdKLE_M?x?_lx1.5H_ly0.5H_var0.5_from_M?x24/0/g"+str(i+1),"w+")    

    gFile.write('''\
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
object      g'''+str(i+1)+''';
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
dimensions      [0 2 -1 0 0 0 0];

internalField   uniform 0;\n

boundaryField
{
    "(inlet|outlet)"
    {
        type            cyclic;
    }''')
    
    if phyDim==2 :
        gFile.write('''
        "(front|back)"
        {
            type            empty;
        }''')
        
    if phyDim==3 :
        gFile.write('''
        "(front|back)"
        {
            type            zeroGradient;
        }''')
    
    gFile.write('''
    "(top|hills)"
    {
        type            zeroGradient;
    }
}
''')
    gFile.close()
