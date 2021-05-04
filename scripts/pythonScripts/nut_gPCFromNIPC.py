# -*- coding: utf-8 -*-
"""

Non-intrusive polynomail chaos
for Periodic hill for uncertain
sigma parameter

"""
import math
import numpy as np
import chaospy as cp
import pyvista as pv
import os
cwd = os.getcwd() + "/"

# <codecell>
a = 1.0           # Lower bound
b = 1.0           # Upper bound

sMean  = 8.35e-2  # mean
sSigma = 1.04e-2  # sigma

d   = 1
Np = 3           # Poly Order
Nq = 5           # Quad pts in each dim

P = int(math.factorial(Np + d)/math.factorial(Np)/math.factorial(d)-1)

nu = 2.643e-6     # Laminar viscosity

# <codecell>
# polynomial generation
Uniform = 0
Normal  = 1

if (Uniform):
    gDist      = cp.Uniform(-1,1)                          # Germ distribution
    sDist      = cp.Uniform(a,b)                           # Cs distribution
    orthpoly    = cp.orth_ttr(Np, sDist, normed=True)    # Orthonormal polys
    nodes, wts  = cp.generate_quadrature(Nq, sDist,'C')  # Nodes and Weights

if (Normal):
    gDist      = cp.Normal(0,1)                            # Germ distribution
    sDist      = cp.Normal(sMean, sSigma)              # Cs distribution
    orthpoly    = cp.orth_ttr(Np, sDist, normed=True)    # Orthonormal polys
    nodes, wts  = cp.generate_quadrature(Nq, sDist,'G')  # Nodes and Weights

Nn = len(nodes[0]) # Total number of quad nodes

# <codecell>
# acces to mesh and cell values
mesh    = pv.UnstructuredGrid(cwd + '1/VTK/1_20000.vtk')
nCells  = mesh.n_cells
cellCtr = mesh.cell_centers().points
#mbounds = mesh.bounds
#nArrays = mesh.n_arrays
#mCentre = mesh.center
#mPoints = mesh.points

# <codecell>
# initialization for samples/ evaluations
nutEvals = np.zeros((Nn,nCells))
pEvals   = np.zeros((Nn,nCells))
UEvals   = np.zeros(((Nn,nCells,3)))
kEvals   = np.zeros((Nn,nCells))

# initialization for PCE modes
nutModes = np.zeros((P+1,nCells))
pModes   = np.zeros((P+1,nCells))
UModes   = np.zeros(((P+1,nCells,3)))
kModes   = np.zeros((P+1,nCells))

# <codecell>
# read the fields
for i in range(Nn):
    vtkFile = cwd + str(i+1) + '/VTK/'+str(i+1)+'_20000.vtk'
    mesh = pv.UnstructuredGrid(vtkFile)
    nutEvals[i] = mesh.cell_arrays['nut'] 
    pEvals[i]   = mesh.cell_arrays['p'] 
    UEvals[i]   = mesh.cell_arrays['U']
    kEvals[i]   = mesh.cell_arrays['k'] 

# <codecell>        
# approximate the PCE of the fields
nutApp = cp.fit_quadrature(orthpoly, nodes, wts, nutEvals)
pApp   = cp.fit_quadrature(orthpoly, nodes, wts, pEvals)
UApp   = cp.fit_quadrature(orthpoly, nodes, wts, UEvals)
kApp   = cp.fit_quadrature(orthpoly, nodes, wts, kEvals)

nutSigma = np.sqrt(cp.Var(nutApp, sDist))
pSigma   = np.sqrt(cp.Var(pApp,   sDist))
USigma   = np.sqrt(cp.Var(UApp,   sDist))
kSigma   = np.sqrt(cp.Var(kApp,   sDist))

# <codecell>
# compute the modes of fields
for i in range(P+1):
    nutModes[i] = cp.E(nutApp*orthpoly[i], sDist)
    pModes[i]   = cp.E(pApp*orthpoly[i],   sDist)
    UModes[i]   = cp.E(UApp*orthpoly[i],   sDist)
    kModes[i]   = cp.E(kApp*orthpoly[i],   sDist)

# add modes to the mesh
for i in range(P+1):
    mesh._add_cell_array(nutModes[i], 'nut'+str(i))
    mesh._add_cell_array(pModes[i],   'p'+str(i))
    mesh._add_cell_array(UModes[i],   'U'+str(i))
    mesh._add_cell_array(kModes[i],   'k'+str(i))

mesh._add_cell_array(nutSigma, 'nutSigma')
mesh._add_cell_array(pSigma,   'pSigma')
mesh._add_cell_array(USigma,   'USigma')
mesh._add_cell_array(kSigma,   'nutSigma')

# save the modes to vtk format
mesh._add_cell_array(cellCtr, 'CVs')
mesh.save(cwd + 'baseline_StSt_nutPCE_10000_PCEs.vtk')   

# <codecell> writing KLE modes of lognormal nut
for i in range(len(nutModes)):

    nutFile = open("PCEs_PerHl_kEps_Cmu/nut"+str(i),"w+")
    
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
internalField   nonuniform List<scalar> \n\
8000\n\
(\n")

    for j in range(nCells):
        nutFile.write(str(nutModes[i][j]) + "\n")
    
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
for i in range(len(kModes)):

    kFile = open("PCEs_PerHl_kEps_Cmu/k"+str(i),"w+")
    
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
internalField   nonuniform List<scalar> \n\
8000\n\
(\n")

    for j in range(nCells):
        kFile.write(str(kModes[i][j]) + "\n")
    
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