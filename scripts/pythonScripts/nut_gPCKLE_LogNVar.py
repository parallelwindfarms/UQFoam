# -*- coding: utf-8 -*-
"""

Non-intrusive polynomail chaos
for Periodic hill for uncertain
sigma parameter

"""

from __future__ import print_function

import numpy as np
import pyvista as pv
import os

cwd = os.getcwd() + "/"

# <codecell> Read OpenFOAM case mesh

# acces to grid and cell values
mesh    = pv.UnstructuredGrid(cwd + 'VTK/1baseline_kEpsilon_nutLogNVarialbe_nutgPCKLE_M128x192_10000.vtk')
nCells  = mesh.n_cells
cCenter = mesh.cell_centers().points
cBounds = mesh.cell_centers().bounds
nPoints  = mesh.n_points
mBounds = mesh.bounds
mPoints = mesh.points

nut0   = (mesh.cell_arrays['nut'])
k0     = (mesh.cell_arrays['k'])
sigma  = np.sqrt(0.10)
nutl0  = np.exp(np.log(nut0)+0.5*sigma*sigma)
kl0    = np.exp(np.log(k0)+0.5*sigma*sigma)
Nmodes = 4 # mean + 3 modes

# <codecell> writing KLE modes of lognormal nut
for i in range(Nmodes):

    nutFile = open("nut_gPCKLE_LogNVar/nut"+str(i),"w+")
    
    nutFile.write("\
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
    class       volScalarField;\n\
    location    \"0\";\n\
    object      nutMode;\n\
}\n\
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n\
dimensions      [0 2 -1 0 0 0 0];\n\n\
\
internalField   nonuniform List<scalar> \n")

    nutFile.write(str(nCells)+"\n(\n")

    c = np.power(sigma,i)/np.math.factorial(i)
    for j in range(len(nut0)):
        nutFile.write(str(nutl0[j]*c) + "\n")
    
    nutFile.write("\
)\n\
;\n\n\
\
boundaryField\n\
{\n\
    \"(inlet|outlet)\"\n\
    {\n\
        type            cyclic;\n\
    }\n\
\
    \"(front|back)\"\n\
    {\n\
        type            empty;\n\
    }\n\
\
    \"(top|hills)\"\n\
    {\n\
        type            nutkWallFunction;\n\
        value           uniform 0;\n\
    }\n\
}")
    
    nutFile.close()
    

# <codecell> writing KLE modes of lognormal k
for i in range(Nmodes):

    kFile = open("nut_gPCKLE_LogNVar/k"+str(i),"w+")
    
    kFile.write("\
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
    class       volScalarField;\n\
    location    \"0\";\n\
    object      kMoode;\n\
}\n\
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n\
dimensions      [0 2 -2 0 0 0 0];\n\n\
\
internalField   nonuniform List<scalar> \n")

    kFile.write(str(nCells)+"\n(\n")


    c = np.power(sigma,i)/np.math.factorial(i)
    for j in range(len(k0)):
        kFile.write(str(kl0[j]*c) + "\n")
    
    kFile.write("\
)\n\
;\n\n\
\
boundaryField\n\
{\n\
    \"(inlet|outlet)\"\n\
    {\n\
        type            cyclic;\n\
    }\n\
\
    \"(front|back)\"\n\
    {\n\
        type            empty;\n\
    }\n\
\
    \"(top|hills)\"\n\
    {\n\
        type            kqRWallFunction;\n\
        value           uniform 0;\n\
    }\n\
}")

    kFile.close()









