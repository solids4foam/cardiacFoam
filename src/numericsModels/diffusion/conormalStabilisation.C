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

#include "conormalStabilisation.H"
#include "fvm.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * Functions * * * * * * * * * * * * * * * //

Foam::tmp<Foam::fvScalarMatrix> Foam::conormalStabilisation
(
    const surfaceTensorField& Gf,
    const volScalarField& psi
)
{
    const fvMesh& mesh = psi.mesh();
    const word name("laplacian(" + Gf.name() + ',' + psi.name() + ')');

    const surfaceVectorField GSn((mesh.Sf()/mesh.magSf()) & Gf);
    const surfaceScalarField dk
    (
        "conormalStabilisation(" + Gf.name() + ')',
        mag(GSn) - (GSn & mesh.Sf())/mesh.magSf()
    );

    return fvm::laplacian(dk, psi, name) - fvc::laplacian(dk, psi, name);
}


// ************************************************************************* //
