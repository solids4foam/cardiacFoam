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

#include "facePeclet.H"
#include "PstreamReduceOps.H"

// * * * * * * * * * * * * * * * * Functions * * * * * * * * * * * * * * * //

Foam::tmp<Foam::surfaceScalarField> Foam::facePeclet
(
    const surfaceScalarField& phi,
    const surfaceTensorField& Gf
)
{
    const fvMesh& mesh = phi.mesh();

    tmp<surfaceScalarField> tPe
    (
        new surfaceScalarField
        (
            IOobject
            (
                "Peclet(" + phi.name() + ')',
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar(dimless, Zero)
        )
    );
    scalarField& Pe = tPe.ref().primitiveFieldRef();

    const vectorField& Sf = mesh.Sf();
    const scalarField& magSf = mesh.magSf();
    const scalarField& deltaCoeffs = mesh.nonOrthDeltaCoeffs();

    forAll(Pe, faceI)
    {
        const vector n(Sf[faceI]/magSf[faceI]);
        const scalar diffusion =
            (n & Gf[faceI] & n)*magSf[faceI]*deltaCoeffs[faceI];

        Pe[faceI] = mag(phi[faceI])/max(diffusion, VSMALL);
    }

    return tPe;
}


void Foam::reportFacePeclet(const surfaceScalarField& Pe, const word& fieldName)
{
    const scalarField& values = Pe.primitiveField();

    scalar maxPe = 0;
    scalar sumPe = 0;
    label nAbove = 0;
    forAll(values, faceI)
    {
        maxPe = max(maxPe, values[faceI]);
        sumPe += values[faceI];
        nAbove += (values[faceI] > 2);
    }
    label nFaces = values.size();

    reduce(maxPe, maxOp<scalar>());
    reduce(sumPe, sumOp<scalar>());
    reduce(nAbove, sumOp<label>());
    reduce(nFaces, sumOp<label>());

    Info<< "Face Peclet number of " << fieldName << ": max " << maxPe
        << ", mean " << sumPe/max(nFaces, 1)
        << ", internal faces above 2: "
        << 100.0*nAbove/max(nFaces, 1) << "%" << endl;
}


// ************************************************************************* //
