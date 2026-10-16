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

#include "sealedWallLeastSquaresGrad.H"
#include "sealedWallLeastSquaresVectors.H"
#include "sealedWallPolicy.H"
#include "centredCPCCellToCellStencilObject.H"
#include "leastSquaresGrad.H"
#include "gaussGrad.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(sealedWallLeastSquaresGrad, 0);

    gradScheme<scalar>::addIstreamConstructorToTable
        <sealedWallLeastSquaresGrad>
        addSealedWallLeastSquaresGradScalarIstreamConstructorToTable_;
}
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volVectorField>
Foam::fv::sealedWallLeastSquaresGrad::calcGrad
(
    const volScalarField& vsf,
    const word& name
) const
{
    const fvMesh& mesh = vsf.mesh();
    const word key(sealedWallPolicy::fieldKey(name));
    const labelList& states = sealedWallPolicy::states(mesh, key);

    tmp<volVectorField> tGrad
    (
        leastSquaresGrad<scalar>(mesh).calcGrad(vsf, name)
    );

    if (!sealedWallPolicy::anyWall(states))
    {
        return tGrad;
    }

    const sealedWallLeastSquaresVectors& lsv =
        sealedWallLeastSquaresVectors::New
        (
            sealedWallLeastSquaresVectors::typeName + '_' + key,
            mesh,
            states
        );
    if (lsv.states() != states)
    {
        FatalErrorInFunction
            << "Wall-gradient states of " << key << " on mesh " << mesh.name()
            << " changed from " << lsv.states() << " to " << states
            << exit(FatalError);
    }

    const extendedCentredCellToCellStencil& stencil =
        centredCPCCellToCellStencilObject::New(mesh);

    // Field values on the stencil, as LeastSquaresGrad builds them
    List<scalar> flatVsf(stencil.map().constructSize(), Zero);
    forAll(vsf, celli)
    {
        flatVsf[celli] = vsf[celli];
    }
    forAll(vsf.boundaryField(), patchi)
    {
        const fvPatchScalarField& psf = vsf.boundaryField()[patchi];
        label nCompact = psf.patch().start() - mesh.nInternalFaces() + mesh.nCells();
        forAll(psf, i)
        {
            flatVsf[nCompact++] = psf[i];
        }
    }
    stencil.map().distribute(flatVsf);

    volVectorField& grad = tGrad.ref();
    vectorField& gradIf = grad.primitiveFieldRef();
    const labelListList& stencilAddr = stencil.stencil();
    const labelList& wallCells = lsv.wallCells();
    const List<List<vector>>& lsvs = lsv.vectors();

    forAll(wallCells, i)
    {
        const label celli = wallCells[i];
        const labelList& compactCells = stencilAddr[celli];
        const List<vector>& lsvc = lsvs[i];

        vector g(Zero);
        forAll(compactCells, j)
        {
            g += lsvc[j]*flatVsf[compactCells[j]];
        }
        gradIf[celli] = g;
    }

    grad.correctBoundaryConditions();
    gaussGrad<scalar>::correctBoundaryConditions(vsf, grad);

    return tGrad;
}


// ************************************************************************* //
