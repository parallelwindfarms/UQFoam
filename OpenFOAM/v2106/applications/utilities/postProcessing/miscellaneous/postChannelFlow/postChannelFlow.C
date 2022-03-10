/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
------------------------------------------------------------------------------- License This file is part of OpenFOAM.

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

Application
    postChannel

Description
    Post-processes data from channel flow calculations.

    For each time: calculate: txx, txy,tyy, txy,
    eps, prod, vorticity, enstrophy and helicity. Assuming that the mesh
    is periodic in the x and z directions, collapse Umeanx, Umeany, txx,
    txy and tyy to a line and print them as standard output.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "channelIndex.H"
#include "makeGraph.H"
#include "OSspecific.H"

#include <iostream> 
#include <regex>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//

template<class FieldType, class ReturnType>
ReturnType patchAverage
(
    const fvMesh& mesh,
    const IOobject& fieldHeader,
    const scalar area,
    const label patchI
)
{
    FieldType field(fieldHeader, mesh);

    typename FieldType::value_type sumField =
        pTraits<typename FieldType::value_type>::zero;

    if (area > 0)
    {
        sumField = gSum
        (
            mesh.magSf().boundaryField()[patchI]
            * field.boundaryField()[patchI]
        ) / area;
    }
        
    return sumField;
}
int main(int argc, char *argv[])
{
    argList::noParallel();
    timeSelector::addOptions();

#   include "setRootCase.H"
#   include "createTime.H"

    // Get times list
    instantList timeDirs = timeSelector::select0(runTime, args);

#   include "createMesh.H"
#   include "readTransportProperties.H"

    const word& gFormat = runTime.graphFormat();

    // Setup channel indexing for averaging over channel down to a line

    IOdictionary channelDict
    (
        IOobject
        (
            "postChannelDict",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );
    channelIndex channelIndexing(mesh, channelDict);

    wordList fieldNames
    (
        channelDict.lookup("fields")
    );

    word botPatchName(channelDict.lookup("botPatch"));
    word topPatchName(channelDict.lookup("topPatch"));

    const label patchBot = mesh.boundaryMesh().findPatchID(botPatchName);
    const label patchTop = mesh.boundaryMesh().findPatchID(topPatchName);

    scalar areaBot = gSum(mesh.magSf().boundaryField()[patchBot]);
    scalar areaTop = gSum(mesh.magSf().boundaryField()[patchTop]);

    scalar yBot = mesh.boundaryMesh()[patchBot].faceCentres()[0].y();
    scalar yTop = mesh.boundaryMesh()[patchTop].faceCentres()[0].y();

    // For each time step read all fields
    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Collapsing fields for time " << runTime.timeName() << endl;

        // Average fields over channel down to a line
#       include "collapse.H"
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
