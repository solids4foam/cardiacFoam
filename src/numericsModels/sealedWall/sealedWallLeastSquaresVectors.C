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

#include <algorithm>
#include <cmath>
#include <limits>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(sealedWallLeastSquaresVectors, 0);

    // Eigenvalues of a symmetric tensor in ascending order by cyclic Jacobi
    // rotations (Golub and Van Loan, Matrix Computations, 4th ed., Alg. 8.5.2
    // and 8.5.3), repeated until every off-diagonal entry is zero. A rotation
    // changes the two diagonal entries by -+t a_pq; an entry whose rotation
    // would leave both unchanged in floating point is set to zero instead.
    static vector symmEigenValues(const symmTensor& T)
    {
        scalar a[3][3] =
        {
            {T.xx(), T.xy(), T.xz()},
            {T.xy(), T.yy(), T.yz()},
            {T.xz(), T.yz(), T.zz()}
        };

        while (a[0][1] != 0 || a[0][2] != 0 || a[1][2] != 0)
        {
            for (label p = 0; p < 2; ++p)
            {
                for (label q = p + 1; q < 3; ++q)
                {
                    if (a[p][q] == 0)
                    {
                        continue;
                    }

                    const scalar tau = (a[q][q] - a[p][p])/(2*a[p][q]);
                    const scalar t =
                        (tau >= 0 ? 1 : -1)
                       /(mag(tau) + std::hypot(scalar(1), tau));
                    const scalar c = 1/std::hypot(scalar(1), t);
                    const scalar s = t*c;

                    const scalar app = a[p][p] - t*a[p][q];
                    const scalar aqq = a[q][q] + t*a[p][q];

                    if (app == a[p][p] && aqq == a[q][q])
                    {
                        a[p][q] = a[q][p] = 0;
                        continue;
                    }

                    a[p][p] = app;
                    a[q][q] = aqq;
                    a[p][q] = a[q][p] = 0;

                    const label r = 3 - p - q;
                    const scalar arp = a[r][p];
                    const scalar arq = a[r][q];
                    a[r][p] = a[p][r] = c*arp - s*arq;
                    a[r][q] = a[q][r] = s*arp + c*arq;
                }
            }
        }

        vector eVals(a[0][0], a[1][1], a[2][2]);
        std::sort(eVals.begin(), eVals.end());

        return eVals;
    }
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

    // A wall-cell normal matrix is numerically rank deficient when its
    // smallest eigenvalue is within the rounding error of its largest
    const scalar rankTol = 3*std::numeric_limits<scalar>::epsilon();

    scalar minRatio = 1;
    point minRatioCentre(Zero);

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

        const vector eVals(symmEigenValues(dd));

        if (eVals.x() <= rankTol*eVals.z())
        {
            FatalErrorInFunction
                << "Rank-deficient wall-cell stencil for " << name()
                << " on mesh " << mesh.name() << ": cell " << celli
                << " at " << mesh.C()[celli] << ", normal-matrix "
                << "eigenvalues " << eVals << exit(FatalError);
        }

        if (i == 0 || eVals.x()/eVals.z() < minRatio)
        {
            minRatio = eVals.x()/eVals.z();
            minRatioCentre = mesh.C()[celli];
        }

        const symmTensor ddInv(inv(dd) - dd0);
        for (label j = 1; j < lsv.size(); ++j)
        {
            lsv[j] = ddInv & lsv[j];
            lsv[0] -= lsv[j];
        }
    }

    const label nWallCells = returnReduce(wallCells_.size(), sumOp<label>());
    const scalar globalMinRatio = returnReduce(minRatio, minOp<scalar>());
    const label minRatioProc =
        returnReduce
        (
            wallCells_.size() && minRatio == globalMinRatio
          ? Pstream::myProcNo()
          : Pstream::nProcs(),
            minOp<label>()
        );
    if (Pstream::myProcNo() != minRatioProc)
    {
        minRatioCentre = Zero;
    }
    reduce(minRatioCentre, sumOp<point>());

    Info<< "sealedWallLeastSquares: " << name() << " on mesh " << mesh.name()
        << ", " << nWallCells << " wall cells";
    if (nWallCells)
    {
        Info<< ", smallest eigenvalue ratio " << globalMinRatio
            << " at cell centre " << minRatioCentre;
    }
    Info<< endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::sealedWallLeastSquaresVectors::movePoints()
{
    calcVectors();
    return true;
}


// ************************************************************************* //
