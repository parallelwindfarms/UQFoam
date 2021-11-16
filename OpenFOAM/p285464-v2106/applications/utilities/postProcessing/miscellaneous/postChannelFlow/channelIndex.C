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

#include "channelIndex.H"
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
    Foam::channelIndex::vectorComponentsNames_;
*/

/*
const Foam::Enum
<
    Foam::vector::components
>
Foam::channelIndex::vectorComponentsNames_
(
    Foam::vector::components::X, { "x", "y", "z" }
);
*/

const Foam::Enum
<
    Foam::vector::components
>
Foam::channelIndex::vectorComponentsNames_
({
    { vector::components::X, "x" },
    { vector::components::Y, "y" },
    { vector::components::Z, "z" },
});


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Determines face blocking
void Foam::channelIndex::walkOppositeFaces
(
    const polyMesh& mesh,
    const labelList& startFaces,
    boolList& blockedFace
)
{
    const cellList& cells = mesh.cells();
    const faceList& faces = mesh.faces();

    //number of boundary faces 
    label nBnd = mesh.nFaces() - mesh.nInternalFaces();

    Info << "Number of boundary faces: " << nBnd << endl;

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
            if (isFrontBndFace[i] && !blockedFace[faceI])
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
                    FatalErrorIn("channelIndex::walkOppositeFaces(..)")
                        << "Face:" << faceI << " owner cell:" << ownCell
                        << " is not a hex?" << abort(FatalError);
                }
                else
                {
                    if (!blockedFace[oppositeFaceI])
                    {
                        blockedFace[oppositeFaceI] = true;
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
                    FatalErrorIn("channelIndex::walkOppositeFaces(..)")
                        << "Face:" << faceI << " neighbour cell:" << neiCell
                        << " is not a hex?" << abort(FatalError);
                }
                else
                {
                    if (!blockedFace[oppositeFaceI])
                    {
                        blockedFace[oppositeFaceI] = true;
                        newFrontFaces.append(oppositeFaceI);
                    }
                }
            }
        }

        frontFaces.transfer(newFrontFaces);
    }
}


// Calculate regions.
void Foam::channelIndex::calcLayeredRegions
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

    // Print out the blockedFaces for debugging puproses
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


    // regionSplit splits the mesh into regions.
    // Each region should be a slice of the mesh with the same value along dir_
    cellRegion_.reset(new regionSplit(mesh, blockedFace));

    Info<< "Detected " << cellRegion_().nRegions() << " layers." << nl << endl;

    // Number of cells per region
    regionCount_ = regionSum(scalarField(mesh.nCells(), 1.0));

    // Average cell centres to determine ordering.
    // Note: pointField = vectorField
    pointField regionCc
    (
        regionSum(mesh.cellCentres())
      / regionCount_
    );

    // Sort the cell-centers after dir_
    // Note: gets sorted on construction
    SortableList<scalar> sortComponent(regionCc.component(dir_));

    //Save the original indices
    sortMap_ = sortComponent.indices();

    // So y_ is a list of cell-centers sorted after dir_
    y_ = sortComponent;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::channelIndex::channelIndex
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    dir_(vectorComponentsNames_.read(dict.lookup("component")))
{
    //all the patches defined on the mesh
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    //patches to seed from
    const wordList patchNames(dict.lookup("patches"));

    label nFaces = 0;

    //calculate total amount of face on the patches in patchNames
    forAll(patchNames, i)
    {
        const label patchI = patches.findPatchID(patchNames[i]);

        if (patchI == -1)
        {
            FatalErrorIn("channelIndex::channelIndex(const polyMesh&)")
                << "Illegal patch " << patchNames[i]
                << ". Valid patches are " << patches.name()
                << exit(FatalError);
        }

        nFaces += patches[patchI].size();
    }

    //will hold the ids of the faces of the seed-patches
    labelList startFaces(nFaces);
    nFaces = 0;

    forAll(patchNames, i)
    {
        const polyPatch& pp = patches[patchNames[i]];

        //loop through the faces of the patch and store their id
        forAll(pp, j)
        {
            startFaces[nFaces++] = pp.start()+j;
        }
    }

    // Calculate regions.
    calcLayeredRegions(mesh, startFaces);
}


Foam::channelIndex::channelIndex
(
    const polyMesh& mesh,
    const labelList& startFaces,
    const direction dir
)
:
    dir_(dir)
{
    // Calculate regions.
    calcLayeredRegions(mesh, startFaces);
}


// ************************************************************************* //
