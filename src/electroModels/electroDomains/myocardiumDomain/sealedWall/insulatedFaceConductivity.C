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

Function
    insulatedFaceConductivity

Description
    Face conductivity with zero value on sealed heart patches.

Author
    Simao Nieto de Castro, UCD.
\*---------------------------------------------------------------------------*/

#include "insulatedFaceConductivity.H"
#include "linear.H"
#include "zeroGradientFvPatchFields.H"
#include "conormalZeroFluxFvPatchScalarField.H"

namespace Foam
{

tmp<surfaceTensorField> insulatedFaceConductivity
(
    const volTensorField& G,
    const volScalarField& psi,
    const bool seal
)
{
    tmp<surfaceTensorField> tGf(linearInterpolate(G));
    if (!seal)
    {
        return tGf;
    }
    surfaceTensorField::Boundary& Gfb = tGf.ref().boundaryFieldRef();

    forAll(Gfb, patchI)
    {
        const fvPatchScalarField& psiP = psi.boundaryField()[patchI];

        if
        (
            !psiP.coupled()
         && (
                isA<zeroGradientFvPatchScalarField>(psiP)
             || isA<conormalZeroFluxFvPatchScalarField>(psiP)
            )
        )
        {
            Gfb[patchI] = Zero;
        }
    }

    return tGf;
}

}

// ************************************************************************* //
