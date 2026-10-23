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

#include "compatibilityProjection.H"

// * * * * * * * * * * * * * * * * Functions * * * * * * * * * * * * * * * //

void Foam::projectCompatibility(fvScalarMatrix& eqn)
{
    const fvMesh& mesh = eqn.psi().mesh();

    scalar imbalance = sum(eqn.source());
    forAll(eqn.boundaryCoeffs(), patchI)
    {
        if (!mesh.boundary()[patchI].coupled())
        {
            imbalance += sum(eqn.boundaryCoeffs()[patchI]);
        }
    }
    reduce(imbalance, sumOp<scalar>());

    const scalarField& V = mesh.V();
    eqn.source() -= V*(imbalance/gSum(V));

    Info<< eqn.psi().name() << " compatibility projection: imbalance "
        << imbalance << nl;
}


// ************************************************************************* //
