/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2017 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "decayingTaylorGreenVortexMeshMotion.H"
#include "motionInterpolation.H"
#include "motionDiffusivity.H"
#include "fvmLaplacian.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcDiv.H"
#include "fvcGrad.H"
#include "surfaceInterpolate.H"
#include "fvcLaplacian.H"
#include "mapPolyMesh.H"
#include "fvOptions.H"
#include "IFstream.H"

namespace Foam
{
    defineTypeNameAndDebug(decayingTaylorGreenVortexMeshMotion, 0);

    addToRunTimeSelectionTable
    (
        motionSolver,
        decayingTaylorGreenVortexMeshMotion,
        dictionary
    );

    addToRunTimeSelectionTable
    (
        displacementMotionSolver,
        decayingTaylorGreenVortexMeshMotion,
        displacement
    );
}


// * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::decayingTaylorGreenVortexMeshMotion::decayingTaylorGreenVortexMeshMotion
(
    const polyMesh& mesh,
    const IOdictionary& dict
)
:
    displacementMotionSolver(mesh, dict, typeName),
    fvMotionSolver(mesh),
    A_
    (
        readScalar
        (
            dict.subDict("decayingTaylorGreenVortexMeshMotionCoeffs")
                .lookup("scaleFactor")
        )
    ),
    orthogonalMeshMotion_
    (
        dict.subDict("decayingTaylorGreenVortexMeshMotionCoeffs")
            .lookup("orthogonalMeshMotion")
    ),
    xMin_
    (
        readScalar
        (
            dict.subDict("decayingTaylorGreenVortexMeshMotionCoeffs")
                .lookup("xMin")
        )
    ),
    xMax_
    (
        readScalar
        (
            dict.subDict("decayingTaylorGreenVortexMeshMotionCoeffs")
                .lookup("xMax")
        )
    ),
    yMin_
    (
        readScalar
        (
            dict.subDict("decayingTaylorGreenVortexMeshMotionCoeffs")
                .lookup("yMin")
        )
    ),
    yMax_
    (
        readScalar
        (
            dict.subDict("decayingTaylorGreenVortexMeshMotionCoeffs")
                .lookup("yMax")
        )
    ),
    period_
    (
        readScalar
        (
            dict.subDict("decayingTaylorGreenVortexMeshMotionCoeffs")
                .lookup("period")
        )
    )
{}


Foam::decayingTaylorGreenVortexMeshMotion::decayingTaylorGreenVortexMeshMotion
(
    const polyMesh& mesh,
    const IOdictionary& dict,
    const pointVectorField& pointDisplacement,
    const pointIOField& points0
)
:
    displacementMotionSolver(mesh, dict, pointDisplacement, points0, typeName),
    fvMotionSolver(mesh),
    A_
    (
        readScalar
        (
            dict.subDict("decayingTaylorGreenVortexMeshMotionCoeffs")
                .lookup("scaleFactor")
        )
    ),
    orthogonalMeshMotion_
    (
        dict.subDict("decayingTaylorGreenVortexMeshMotionCoeffs")
            .lookup("orthogonalMeshMotion")
    ),
    xMin_
    (
        readScalar
        (
            dict.subDict("decayingTaylorGreenVortexMeshMotionCoeffs")
                .lookup("xMin")
        )
    ),
    xMax_
    (
        readScalar
        (
            dict.subDict("decayingTaylorGreenVortexMeshMotionCoeffs")
                .lookup("xMax")
        )
    ),
    yMin_
    (
        readScalar
        (
            dict.subDict("decayingTaylorGreenVortexMeshMotionCoeffs")
                .lookup("yMin")
        )
    ),
    yMax_
    (
        readScalar
        (
            dict.subDict("decayingTaylorGreenVortexMeshMotionCoeffs")
                .lookup("yMax")
        )
    ),
    period_
    (
        readScalar
        (
            dict.subDict("decayingTaylorGreenVortexMeshMotionCoeffs")
                .lookup("period")
        )
    )
{}


// * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

Foam::decayingTaylorGreenVortexMeshMotion::
~decayingTaylorGreenVortexMeshMotion()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField>
Foam::decayingTaylorGreenVortexMeshMotion::curPoints() const
{
    tmp<pointField> tcurPoints(points0() + pointDisplacement_.primitiveField());
    twoDCorrectPoints(tcurPoints.ref());
    return tcurPoints;
}


void Foam::decayingTaylorGreenVortexMeshMotion::solve()
{
    // Update internal motion-solver state to current mesh position
    movePoints(fvMesh_.points());

    // Initial/reference mesh points
    const vectorField& pts0 = this->points0();

    const scalar pi = constant::mathematical::pi;
    const scalar t  = time().value();

    const scalar Lx = xMax_ - xMin_;
    const scalar Ly = yMax_ - yMin_;

    const scalarField x(pts0.component(vector::X));
    const scalarField y(pts0.component(vector::Y));

    const scalarField xi((x - xMin_)/Lx);
    const scalarField eta((y - yMin_)/Ly);

    const scalar timeFactor = A_*Foam::sin(2.0*pi*t/period_);

    if (orthogonalMeshMotion_)
    {
        pointDisplacement_.primitiveFieldRef() =
              timeFactor*Foam::sin(pi*xi)*vector(1, 0, 0)
            + timeFactor*Foam::sin(pi*eta)*vector(0, 1, 0);
    }
    else
    {
        pointDisplacement_.primitiveFieldRef() =
            timeFactor
           *Foam::sin(pi*xi)
           *Foam::sin(pi*eta)
           *vector(1, 1, 0);
    }

    twoDCorrectPoints(pointDisplacement_.primitiveFieldRef());
}


void Foam::decayingTaylorGreenVortexMeshMotion::updateMesh
(
    const mapPolyMesh& mpm
)
{
    displacementMotionSolver::updateMesh(mpm);
}

// ************************************************************************* //

// ************************************************************************* //
