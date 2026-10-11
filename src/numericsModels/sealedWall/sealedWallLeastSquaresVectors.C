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

#include "sealedWallLeastSquaresVectors.H"
#include "sealedWallPolicy.H"
#include "centredCPCCellToCellStencilObject.H"
#include "volFields.H"
#include "bitSet.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(sealedWallLeastSquaresVectors, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sealedWallLeastSquaresVectors::sealedWallLeastSquaresVectors
(
    const word& name,
    const fvMesh& mesh,
    const labelUList& states
)
:
    MeshObject_type(name, mesh),
    states_(states)
{
    calcVectors();
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::sealedWallLeastSquaresVectors::calcVectors()
{
    const fvMesh& mesh = mesh_;
    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    const extendedCentredCellToCellStencil& stencil =
        centredCPCCellToCellStencilObject::New(mesh);

    // Sealed-face flag carried through the stencil to every entry
    volScalarField sealedFlag
    (
        IOobject
        (
            "sealedWallFlag",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        mesh,
        dimensionedScalar(dimless, Zero),
        fvPatchFieldBase::calculatedType()
    );
    bitSet isWall(mesh.nCells());
    forAll(states_, patchi)
    {
        if (states_[patchi] == sealedWallPolicy::sealed)
        {
            sealedFlag.boundaryFieldRef()[patchi] = 1;
        }
        if (states_[patchi] != sealedWallPolicy::stock)
        {
            isWall.set(patches[patchi].faceCells());
        }
    }
    wallCells_ = isWall.sortedToc();

    List<List<scalar>> flagData;
    stencil.collectData(sealedFlag, flagData);
    List<List<point>> posData;
    stencil.collectData(mesh.C(), posData);

    // Unit weight in the empty directions of a 1-D or 2-D mesh
    const symmTensor dd0(sqr((Vector<label>::one - mesh.geometricD())/2));

    vectors_.setSize(wallCells_.size());
    scalar minRatio = GREAT;

    forAll(wallCells_, i)
    {
        const label celli = wallCells_[i];
        const List<point>& pos = posData[celli];
        const List<scalar>& flag = flagData[celli];
        List<vector>& lsv = vectors_[i];
        lsv.setSize(pos.size(), Zero);

        symmTensor dd(dd0);
        for (label j = 1; j < pos.size(); ++j)
        {
            if (flag[j] > 0.5)
            {
                continue;
            }
            const vector d(pos[j] - pos[0]);
            const scalar w = 1/magSqr(d);
            dd += w*sqr(d);
            lsv[j] = w*d;
        }

        const vector ev(eigenValues(dd));
        minRatio = min(minRatio, cmptMin(ev)/max(cmptMax(ev), VSMALL));

        const symmTensor ddInv(inv(dd) - dd0);
        for (label j = 1; j < lsv.size(); ++j)
        {
            lsv[j] = ddInv & lsv[j];
            lsv[0] -= lsv[j];
        }
    }

    reduce(minRatio, minOp<scalar>());

    Info<< "sealedWallLeastSquares: " << name() << " on mesh " << mesh.name()
        << ", " << returnReduce(wallCells_.size(), sumOp<label>())
        << " wall cells, smallest eigenvalue ratio " << minRatio << endl;

    // Smallest-to-largest eigenvalue ratio of a usable wall-cell stencil
    if (minRatio < 1e-3)
    {
        FatalErrorInFunction
            << "Rank-deficient wall-cell stencil for " << name()
            << " on mesh " << mesh.name() << ": eigenvalue ratio "
            << minRatio << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::sealedWallLeastSquaresVectors::movePoints()
{
    calcVectors();
    return true;
}


// ************************************************************************* //
