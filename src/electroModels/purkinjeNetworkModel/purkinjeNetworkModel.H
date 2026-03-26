/*---------------------------------------------------------------------------*\
License
    This file is part of cardiacFoam.

    cardiacFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    cardiacFoam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with cardiacFoam.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::purkinjeNetworkModel

Description
    Concrete Purkinje network model (stub — full implementation in Task 2).

Author
    Simão Nieto de Castro, UCD.
    Philip Cardiff, UCD.

SourceFiles
    purkinjeNetworkModel.C
\*---------------------------------------------------------------------------*/

#ifndef purkinjeNetworkModel_H
#define purkinjeNetworkModel_H

#include "purkinjeModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class purkinjeNetworkModel Declaration
\*---------------------------------------------------------------------------*/

class purkinjeNetworkModel
:
    public purkinjeModel
{
public:

    TypeName("purkinjeNetworkModel");

    // Constructors

        purkinjeNetworkModel
        (
            const volScalarField& Vm,
            const dictionary&     dict,
            scalar                initialDeltaT
        );

    // Member Functions

        void evolve(scalar, scalar, volScalarField&) override {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
