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

\*---------------------------------------------------------------------------*/

#include "transmissionManufacturedSolution.H"
#include "obliqueSlabTransmission.H"
#include "sphericalShellTransmission.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::transmissionManufacturedSolution>
Foam::transmissionManufacturedSolution::New(const dictionary& verificationDict)
{
    const word geometry(verificationDict.get<word>("geometry"));

    if (geometry == "obliqueSlab")
    {
        return autoPtr<transmissionManufacturedSolution>
        (
            new obliqueSlabTransmission(verificationDict)
        );
    }
    if (geometry == "sphericalShell")
    {
        return autoPtr<transmissionManufacturedSolution>
        (
            new sphericalShellTransmission(verificationDict)
        );
    }

    FatalIOErrorInFunction(verificationDict)
        << "Unknown geometry '" << geometry << "'. Valid geometries: "
        << "(obliqueSlab sphericalShell)."
        << exit(FatalIOError);

    return nullptr;
}


// ************************************************************************* //
