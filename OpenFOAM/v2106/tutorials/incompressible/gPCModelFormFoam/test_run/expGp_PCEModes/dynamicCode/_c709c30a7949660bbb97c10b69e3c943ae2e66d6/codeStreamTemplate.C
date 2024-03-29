/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) YEAR AUTHOR, AFFILIATION
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

Description
    Template for use with codeStream.

\*---------------------------------------------------------------------------*/

#include "dictionary.H"
#include "Ostream.H"
#include "Pstream.H"
#include "pointField.H"
#include "tensor.H"
#include "unitConversion.H"

//{{{ begin codeInclude
#line 43 "/home/p285464/OpenFOAM/p285464-v2106/tutorials/incompressible/gPCModelFormSimpleFoam/1nutKLEgPC_lx1.5H_ly0.5H_var0.50_M100x80/expGp_PCEModes/system/blockMeshDict.#codeStream"
#include "pointField.H"
        #include "mathematicalConstants.H"
//}}} end codeInclude

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//{{{ begin localCode

//}}} end localCode


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

extern "C" void codeStream_c709c30a7949660bbb97c10b69e3c943ae2e66d6(Foam::Ostream& os, const Foam::dictionary& dict)
{
//{{{ begin code
    #line 49 "/home/p285464/OpenFOAM/p285464-v2106/tutorials/incompressible/gPCModelFormSimpleFoam/1nutKLEgPC_lx1.5H_ly0.5H_var0.50_M100x80/expGp_PCEModes/system/blockMeshDict.#codeStream"
const scalar xMin = 0;
        const scalar xMax = 252;
        const label nPoints = 1000;
        const scalar dx = (xMax - xMin)/scalar(nPoints - 1);

        os  << "(" << nl << "spline 0 1" << nl;
        pointField profile(nPoints, Zero);

        for (label i = 0; i < nPoints; ++i)
        {
            scalar x = xMin + i*dx;
            profile[i].x() = x;
            if (x > 198) x = 252 - x;

            if (x >= 0 && x < 9)
            {
                profile[i].y() =
                    28
                  + 6.775070969851E-03*x*x
                  - 2.124527775800E-03*x*x*x;
            }
            else if (x >= 9 && x < 14)
            {
                profile[i].y() =
                    25.07355893131
                  + 0.9754803562315*x
                  - 1.016116352781E-01*x*x
                  + 1.889794677828E-03*x*x*x;
            }
            else if (x >= 14 && x < 20)
            {
                profile[i].y() =
                    2.579601052357E+01
                  + 8.206693007457E-01*x
                  - 9.055370274339E-02*x*x
                  + 1.626510569859E-03*x*x*x;
            }
            else if (x >= 20 && x < 30)
            {
                profile[i].y() =
                    4.046435022819E+01
                  - 1.379581654948E+00*x
                  + 1.945884504128E-02*x*x
                  - 2.070318932190E-04*x*x*x;
            }
            else if (x >= 30 && x < 40)
            {
                profile[i].y() =
                    1.792461334664E+01
                  + 8.743920332081E-01*x
                  - 5.567361123058E-02*x*x
                  + 6.277731764683E-04*x*x*x;
            }
            else if (x >= 40 && x < 54)
            {
                profile[i].y() =
                    max
                    (
                        0,
                        5.639011190988E+01
                      - 2.010520359035E+00*x
                      + 1.644919857549E-02*x*x
                      + 2.674976141766E-05*x*x*x
                    );
            }
            profile[i].z() = 0;
        }
        os << profile << nl;

        os << "spline 4 5" << nl;
        profile.replace(2, 126);
        os << profile << nl;

        os  << ");" << nl;
//}}} end code
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

