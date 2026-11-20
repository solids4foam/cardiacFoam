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

Author
    Simao Nieto de Castro, UCD.
\*---------------------------------------------------------------------------*/

#include "compatibilityImbalance.H"

// * * * * * * * * * * * * * * * * Functions * * * * * * * * * * * * * * * //

Foam::scalar Foam::reportCompatibilityImbalance(const fvScalarMatrix& eqn)
{
    const fvMesh& mesh = eqn.psi().mesh();

    scalar imbalance = sum(eqn.source());
    scalar scale = sum(mag(eqn.source()));

    forAll(eqn.boundaryCoeffs(), patchI)
    {
        if (!mesh.boundary()[patchI].coupled())
        {
            imbalance += sum(eqn.boundaryCoeffs()[patchI]);
            scale += sum(mag(eqn.boundaryCoeffs()[patchI]));
        }
    }

    reduce(imbalance, sumOp<scalar>());
    reduce(scale, sumOp<scalar>());

    Info<< eqn.psi().name() << " pure-Neumann imbalance " << imbalance
        << " (relative " << mag(imbalance)/max(scale, VSMALL) << ')' << nl;

    return imbalance;
}


// ************************************************************************* //
