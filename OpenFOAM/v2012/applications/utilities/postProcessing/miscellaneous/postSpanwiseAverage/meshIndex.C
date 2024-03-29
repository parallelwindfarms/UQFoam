/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "meshIndex.H"
#include "boolList.H"
#include "syncTools.H"
#include "OFstream.H"
#include "meshTools.H"
#include "Time.H"
#include "SortableList.H"

// * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * * //
/*
namespace Foam
{
    template<>
    const char* Foam::NamedEnum
    <
        Foam::vector::components,
        3
    >::names[] =
    {
        "x",
        "y",
        "z"
    };
}

const Foam::NamedEnum<Foam::vector::components, 3>
    Foam::meshIndex::vectorComponentsNames_;
*/

/*
const Foam::Enum
<
    Foam::vector::components
>
Foam::meshIndex::vectorComponentsNames_
(
    Foam::vector::components::X, { "x", "y", "z" }
);
*/

const Foam::Enum
<
    Foam::vector::components
>
Foam::meshIndex::vectorComponentsNames_
({
    { vector::components::X, "x" },
    { vector::components::Y, "y" },
    { vector::components::Z, "z" },
});

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Determines face blocking
void Foam::meshIndex::walkOppositeFaces
(
    const polyMesh& mesh,
    const labelList& startFaces,
    boolList& blockedFace
)
{
    const cellList& cells = mesh.cells();
    const faceList& faces = mesh.faces();
    label nBnd = mesh.nFaces() - mesh.nInternalFaces();

    forAll(faces, i)
    {
        blockedFace[i] = true;
    }

    DynamicList<label> frontFaces(startFaces);
    
    forAll(frontFaces, i)
    {
        label faceI = frontFaces[i];
        blockedFace[faceI] = true;
    }
    

    while (returnReduce(frontFaces.size(), sumOp<label>()) > 0)
    {
        // Transfer across.
        boolList isFrontBndFace(nBnd, false);
        forAll(frontFaces, i)
        {
            label faceI = frontFaces[i];

            if (!mesh.isInternalFace(faceI))
            {
                isFrontBndFace[faceI-mesh.nInternalFaces()] = true;
            }
        }

        syncTools::swapBoundaryFaceList(mesh, isFrontBndFace);

        // Add
        forAll(isFrontBndFace, i)
        {
            label faceI = mesh.nInternalFaces()+i;

            if (isFrontBndFace[i] && blockedFace[faceI])
            {
                blockedFace[faceI] = true;
                frontFaces.append(faceI);
            }
            
        }

        // Transfer across cells
        DynamicList<label> newFrontFaces(frontFaces.size());

        forAll(frontFaces, i)
        {
            label faceI = frontFaces[i];

            {
                const cell& ownCell = cells[mesh.faceOwner()[faceI]];

                label oppositeFaceI = ownCell.opposingFaceLabel(faceI, faces);

                if (oppositeFaceI == -1)
                {
                    FatalErrorIn("meshIndex::walkOppositeFaces(..)")
                        << "Face:" << faceI << " owner cell:" << ownCell
                        << " is not a hex?" << abort(FatalError);
                }
                else
                {
                    if (blockedFace[oppositeFaceI])
                    {
                        blockedFace[oppositeFaceI] = false;
                        newFrontFaces.append(oppositeFaceI);
                    }

                }
            }

            if (mesh.isInternalFace(faceI))
            {
                const cell& neiCell = mesh.cells()[mesh.faceNeighbour()[faceI]];

                label oppositeFaceI = neiCell.opposingFaceLabel(faceI, faces);

                if (oppositeFaceI == -1)
                {
                    FatalErrorIn("meshIndex::walkOppositeFaces(..)")
                        << "Face:" << faceI << " neighbour cell:" << neiCell
                        << " is not a hex?" << abort(FatalError);
                }
                else
                {
                    if (blockedFace[oppositeFaceI])
                    {
                        blockedFace[oppositeFaceI] = false;
                        newFrontFaces.append(oppositeFaceI);
                    }
                    
                }
            }
        }

        frontFaces.transfer(newFrontFaces);
    }
}


// Calculates regions.
void Foam::meshIndex::calcLayeredRegions
(
    const polyMesh& mesh,
    const labelList& startFaces
)
{
    boolList blockedFace(mesh.nFaces(), false);
    walkOppositeFaces
    (
        mesh,
        startFaces,
        blockedFace
    );


    if (false)
    {
        OFstream str(mesh.time().path()/"blockedFaces.obj");
        label vertI = 0;
        forAll(blockedFace, faceI)
        {
            if (blockedFace[faceI])
            {
                const face& f = mesh.faces()[faceI];
                forAll(f, fp)
                {
                    meshTools::writeOBJ(str, mesh.points()[f[fp]]);
                }
                str<< 'f';
                forAll(f, fp)
                {
                    str << ' ' << vertI+fp+1;
                }
                str << nl;
                vertI += f.size();
            }
        }
    }


    // Do analysis for connected regions
    cellRegion_.reset(new regionSplit(mesh, blockedFace));

    Info<< "Detected " << cellRegion_().nRegions() << " layers." << nl << endl;

    // Sum total number of cells
    cellCount_ = mesh.nCells();
    Info<< "Number of cells: " << cellCount_ << endl;

    // Sum number of entries per region
    regionCount_ = regionSum(scalarField(mesh.nCells(), 1.0));
//     Info<< "Number of cell regions: " << regionCount_ << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::meshIndex::meshIndex
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    dir_(vectorComponentsNames_.read(dict.lookup("component"))),
    dirAlt1_(vectorComponentsNames_.read(dict.lookup("componentAlt1"))),
    dirAlt2_(vectorComponentsNames_.read(dict.lookup("componentAlt2")))
{
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    const wordList patchNames(dict.lookup("patches"));

    label nFaces = 0;

    forAll(patchNames, i)
    {
        const label patchI = patches.findPatchID(patchNames[i]);

        if (patchI == -1)
        {
            FatalErrorIn("meshIndex::meshIndex(const polyMesh&)")
                << "Illegal patch " << patchNames[i]
                << ". Valid patches are " << patches.name()
                << exit(FatalError);
        }

        nFaces += patches[patchI].size();
    }

    labelList startFaces(nFaces);
    nFaces = 0;

    forAll(patchNames, i)
    {
        const polyPatch& pp = patches[patchNames[i]];

        forAll(pp, j)
        {
            startFaces[nFaces++] = pp.start()+j;
        }
    }

    // Calculate regions.
    calcLayeredRegions(mesh, startFaces);
}


Foam::meshIndex::meshIndex
(
    const polyMesh& mesh,
    const labelList& startFaces,
    const direction dir,
    const direction dirAlt1,
    const direction dirAlt2
)
:
    dir_(dir),
    dirAlt1_(dirAlt1),
    dirAlt2_(dirAlt2)
{
    // Calculate regions.
    calcLayeredRegions(mesh, startFaces);
}


// ************************************************************************* //
