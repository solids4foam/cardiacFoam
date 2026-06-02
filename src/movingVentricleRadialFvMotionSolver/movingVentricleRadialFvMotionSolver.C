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

#include "movingVentricleRadialFvMotionSolver.H"
#include "motionInterpolation.H"
#include "motionDiffusivity.H"
#include "fvmLaplacian.H"
#include "addToRunTimeSelectionTable.H"
#include "mapPolyMesh.H"
#include "PatchTools.H"
#include "fvOptions.H"
#include "syncTools.H"
#include <cmath>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(movingVentricleRadialFvMotionSolver, 0);

    addToRunTimeSelectionTable
    (
        motionSolver,
        movingVentricleRadialFvMotionSolver,
        dictionary
    );

    addToRunTimeSelectionTable
    (
        displacementMotionSolver,
        movingVentricleRadialFvMotionSolver,
        displacement
    );
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::movingVentricleRadialFvMotionSolver::initialise(const dictionary& dict)
{
    if (patchID_ == -1)
    {
        FatalErrorInFunction
            << "Did not find the patch = " << word(dict.lookup("patchName"))
            << endl;
    }

    Info<< typeName << ": minimum volume = " << volumeMin_ << endl;

    // ---------------- Radial motion settings ----------------
    centrePoint_ =
        motionSolver::coeffDict().getOrDefault<point>("centrePoint", point::zero);

    planeP1_ =
        motionSolver::coeffDict().get<point>("planeP1");

    planeP2_ =
        motionSolver::coeffDict().get<point>("planeP2");

    planeP3_ =
        motionSolver::coeffDict().get<point>("planeP3");

    planeNormal_ = (planeP2_ - planeP1_) ^ (planeP3_ - planeP1_);
    const scalar planeMag = mag(planeNormal_);

    if (planeMag < VSMALL)
    {
        FatalErrorInFunction
            << "planeP1, planeP2 and planeP3 are collinear or too close"
            << exit(FatalError);
    }

    // Unit normal
    planeNormal_ /= planeMag;

    radialExponent_ =
        motionSolver::coeffDict().getOrDefault<scalar>("radialExponent", 1.0);

    radialProfile_ =
        motionSolver::coeffDict().getOrDefault<word>
        (
            "radialProfile",
            "apexWeighted"
        );

    peakLocation_ =
        motionSolver::coeffDict().getOrDefault<scalar>("peakLocation", 0.5);

    peakSharpness_ =
        motionSolver::coeffDict().getOrDefault<scalar>("peakSharpness", 1.0);

    minProfileWeight_ =
        motionSolver::coeffDict().getOrDefault<scalar>("minProfileWeight", 0.0);

    apexProfileWeight_ =
        motionSolver::coeffDict().getOrDefault<scalar>("apexProfileWeight", 0.3);

    dispScale_ =
        motionSolver::coeffDict().getOrDefault<scalar>("dispScale", 1e-2);

    if
    (
        radialProfile_ != "apexWeighted"
     && radialProfile_ != "midPeak"
     && radialProfile_ != "shiftedPeak"
     && radialProfile_ != "baseZeroMidPeak"
    )
    {
        FatalErrorInFunction
            << "Unknown radialProfile " << radialProfile_ << nl
            << "Valid options are apexWeighted, midPeak, shiftedPeak, "
            << "baseZeroMidPeak."
            << exit(FatalError);
    }

    if (peakLocation_ <= 0.0 || peakLocation_ >= 1.0)
    {
        FatalErrorInFunction
            << "peakLocation must lie strictly between 0 and 1. "
            << "Current value: " << peakLocation_
            << exit(FatalError);
    }

    if (peakSharpness_ <= VSMALL)
    {
        FatalErrorInFunction
            << "peakSharpness must be positive. "
            << "Current value: " << peakSharpness_
            << exit(FatalError);
    }

    if (minProfileWeight_ < 0.0 || minProfileWeight_ >= 1.0)
    {
        FatalErrorInFunction
            << "minProfileWeight must lie in the range [0, 1). "
            << "Current value: " << minProfileWeight_
            << exit(FatalError);
    }

    if (apexProfileWeight_ < 0.0 || apexProfileWeight_ > 1.0)
    {
        FatalErrorInFunction
            << "apexProfileWeight must lie in the range [0, 1]. "
            << "Current value: " << apexProfileWeight_
            << exit(FatalError);
    }

    // Read time vs flow rate series
    interpolationTable<scalar> flowTable
    (
        motionSolver::coeffDict().subDict("timeVsFlowSeries")
    );

    // Extract times and flow values
    times_.setSize(flowTable.size());
    flowRates_.setSize(flowTable.size());

    {
        label i = 0;
        forAllConstIter(interpolationTable<scalar>, flowTable, iter)
        {
            const Tuple2<scalar, scalar>& tp = *iter;
            times_[i] = tp.first();
            flowRates_[i] = tp.second();
            ++i;
        }
    }

    // Integrate the flow rate using trapezoidal rule to get the volume versus
    // time
    volumes_.setSize(times_.size(), 0.0);
    for (label i = 1; i < times_.size(); ++i)
    {
        const scalar dt = times_[i] - times_[i-1];
        const scalar avgFlow = 0.5*(flowRates_[i] + flowRates_[i - 1]);
        volumes_[i] = volumes_[i - 1] + avgFlow*dt;
    }

    // Offset by minimum volume
    forAll(volumes_, i)
    {
        volumes_[i] += volumeMin_;
    }

    Info<< typeName << ": using timeVsFlowSeries with "
        << times_.size() << " samples, volumeMin = "
        << volumeMin_ << endl;

    // Set the initial patch points
    initialPatchPoints_ = fvMesh_.boundaryMesh()[patchID_].localPoints();

    // Set the initial point normals
    initialPatchPointNormals_ =
        PatchTools().pointNormals(fvMesh_, fvMesh_.boundaryMesh()[patchID_]);

    // Determine normalization distance directly from the moving patch points
    maxNormalDistance_ = -GREAT;

    forAll(initialPatchPoints_, pI)
    {
        const point& p0 = initialPatchPoints_[pI];
        const scalar s = (p0 - centrePoint_) & planeNormal_;
        maxNormalDistance_ = max(maxNormalDistance_, s);
    }

        reduce(maxNormalDistance_, maxOp<scalar>());

    if (maxNormalDistance_ < VSMALL)
    {
        FatalErrorInFunction
            << "Computed maxNormalDistance_ is too small or negative. "
            << "Check centrePoint and plane point ordering."
            << nl
            << "Normal direction is set by (planeP2-planeP1)^(planeP3-planeP1)."
            << exit(FatalError);
    }

    Info<< typeName
        << ": plane normal = " << planeNormal_
        << ", maxNormalDistance = " << maxNormalDistance_
        << endl;

    Info<< typeName
        << ": radialProfile = " << radialProfile_
        << ", peakLocation = " << peakLocation_
        << ", peakSharpness = " << peakSharpness_
        << ", minProfileWeight = " << minProfileWeight_
        << ", apexProfileWeight = " << apexProfileWeight_
        << endl;
}


Foam::scalar Foam::movingVentricleRadialFvMotionSolver::linearInterp
(
    const scalarField& x, const scalarField& y, const scalar xq
) const
{
    if (xq <= x.first())
    {
        return y.first();
    }

    if (xq >= x.last())
    {
        return y.last();
    }

    label i = 1;
    for (; i < x.size(); ++i)
    {
        if (x[i] >= xq) break;
    }

    const scalar x0 = x[i - 1], x1 = x[i];
    const scalar y0 = y[i - 1], y1 = y[i];
    const scalar w  = (xq - x0)/max(VSMALL, x1 - x0);

    return y0 + w*(y1 - y0);
}


Foam::scalar Foam::movingVentricleRadialFvMotionSolver::profileWeight
(
    const scalar t
) const
{
    const scalar tc = max(0.0, min(1.0, t));
    scalar w = 0.0;

    if (radialProfile_ == "apexWeighted")
    {
        w = pow(tc, radialExponent_);
    }
    else if (radialProfile_ == "midPeak")
    {
        const scalar bell = 4.0*tc*(1.0 - tc);
        w = pow(max(0.0, bell), peakSharpness_);
    }
    else if (radialProfile_ == "shiftedPeak")
    {
        const scalar span =
            max(peakLocation_, 1.0 - peakLocation_);
        const scalar xi = mag(tc - peakLocation_)/max(VSMALL, span);
        const scalar smoothBell = 0.5*(1.0 + cos(constant::mathematical::pi*xi));
        w = pow(max(0.0, smoothBell), peakSharpness_);
    }
    else if (radialProfile_ == "baseZeroMidPeak")
    {
        if (tc <= peakLocation_)
        {
            const scalar u = tc/max(VSMALL, peakLocation_);
            const scalar smoothRise = 3.0*u*u - 2.0*u*u*u;
            w = pow(max(0.0, smoothRise), peakSharpness_);
        }
        else
        {
            const scalar u =
                (tc - peakLocation_)/max(VSMALL, 1.0 - peakLocation_);
            const scalar smoothDecay = 1.0 - (3.0*u*u - 2.0*u*u*u);
            w = apexProfileWeight_ + (1.0 - apexProfileWeight_)*smoothDecay;
        }
    }
    else
    {
        FatalErrorInFunction
            << "Unknown radialProfile " << radialProfile_
            << exit(FatalError);
    }

    return minProfileWeight_ + (1.0 - minProfileWeight_)*w;
}


Foam::scalar Foam::movingVentricleRadialFvMotionSolver::calculateVolume
(
    const fvMesh& mesh,
    const pointField& initialPoints,
    const labelList& ventricleMeshPoints,
    const vectorField& displacement
) const
{
    // We are going to do something a bit hacky: we will use const_cast to move
    // the mesh by the given displacement field, then calculate the new
    // volume, and finally reset the mesh before returning the volume

    // Take a copy of the old time mesh points
    const pointField oldPoints(mesh.points());

    // Full point displacement field over local mesh points
    vectorField pointDisp(mesh.nPoints(), vector::zero);
    scalarField pointCount(mesh.nPoints(), 0.0);

    // Insert displacement on moving patch points
    forAll(ventricleMeshPoints, pI)
    {
        const label pointID = ventricleMeshPoints[pI];
        pointDisp[pointID] += displacement[pI];
        pointCount[pointID] += 1.0;
    }

    // Synchronize point displacements and ownership counts across processors
    syncTools::syncPointList
    (
        const_cast<fvMesh&>(mesh),
        pointDisp,
        plusEqOp<vector>(),
        vector::zero
    );

    syncTools::syncPointList
    (
        const_cast<fvMesh&>(mesh),
        pointCount,
        plusEqOp<scalar>(),
        0.0
    );

    // Average on shared points so all processors get identical point motion
    forAll(pointDisp, pointI)
    {
        if (pointCount[pointI] > SMALL)
        {
            pointDisp[pointI] /= pointCount[pointI];
        }
    }

   // Build moved point coordinates
    pointField movedPoints(initialPoints);

    forAll(movedPoints, pointI)
    {
        movedPoints[pointI] += pointDisp[pointI];
    }


    // Move the mesh
    const_cast<fvMesh&>(mesh).movePoints(movedPoints);

    // Calculate the volume
    const scalar totalVolume = gSum(mesh.cellVolumes());

    // Reset the mesh
    const_cast<fvMesh&>(mesh).movePoints(oldPoints);

    return totalVolume;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::movingVentricleRadialFvMotionSolver::movingVentricleRadialFvMotionSolver
(
    const polyMesh& mesh,
    const IOdictionary& dict
)
:
    displacementMotionSolver(mesh, dict, typeName),
    fvMotionSolver(mesh),
    cellDisplacement_
    (
        IOobject
        (
            "cellDisplacement",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvMesh_,
        dimensionedVector(pointDisplacement_.dimensions(), Zero),
        cellMotionBoundaryTypes<vector>(pointDisplacement_.boundaryField())
    ),
    interpolationPtr_
    (
        coeffDict().found("interpolation")
      ? motionInterpolation::New(fvMesh_, coeffDict().lookup("interpolation"))
      : motionInterpolation::New(fvMesh_)
    ),
    diffusivityPtr_
    (
        motionDiffusivity::New(fvMesh_, coeffDict().lookup("diffusivity"))
    ),
    frozenPointsZone_
    (
        coeffDict().found("frozenPointsZone")
      ? fvMesh_.pointZones().findZoneID
        (
            coeffDict().get<word>("frozenPointsZone")
        )
      : -1
    ),
    patchID_
    (
        mesh.boundaryMesh().findPatchID
        (
            word(motionSolver::coeffDict().lookup("patchName"))
        )
    ),
    volumeMin_(gSum(mesh.cellVolumes())),
    useNewtonRaphson_
    (
        motionSolver::coeffDict().lookupOrDefault<Switch>
        (
            "useNewtonRaphson", true
        )
    ),
    relTol_
    (
        motionSolver::coeffDict().lookupOrDefault<scalar>("relTol", 1e-6)
    ),
    times_(),
    flowRates_(),
    volumes_(),
    initialPatchPoints_(),
    initialPatchPointNormals_(),
    centrePoint_(point::zero),
    planeP1_(point::zero),
    planeP2_(point::zero),
    planeP3_(point::zero),
    planeNormal_(vector::zero),
    maxNormalDistance_(0.0),
    radialExponent_(1.0),
    radialProfile_("apexWeighted"),
    peakLocation_(0.5),
    peakSharpness_(1.0),
    minProfileWeight_(0.0),
    apexProfileWeight_(0.3),
    dispScale_(1e-2)
{
    initialise(dict);
}


Foam::movingVentricleRadialFvMotionSolver::movingVentricleRadialFvMotionSolver
(
    const polyMesh& mesh,
    const IOdictionary& dict,
    const pointVectorField& pointDisplacement,
    const pointIOField& points0
)
:
    displacementMotionSolver(mesh, dict, pointDisplacement, points0, typeName),
    fvMotionSolver(mesh),
    cellDisplacement_
    (
        IOobject
        (
            "cellDisplacement",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvMesh_,
        dimensionedVector(pointDisplacement_.dimensions(), Zero),
        cellMotionBoundaryTypes<vector>(pointDisplacement_.boundaryField())
    ),
    interpolationPtr_
    (
        coeffDict().found("interpolation")
      ? motionInterpolation::New(fvMesh_, coeffDict().lookup("interpolation"))
      : motionInterpolation::New(fvMesh_)
    ),
    diffusivityPtr_
    (
        motionDiffusivity::New(fvMesh_, coeffDict().lookup("diffusivity"))
    ),
    frozenPointsZone_
    (
        coeffDict().found("frozenPointsZone")
      ? fvMesh_.pointZones().findZoneID
        (
            coeffDict().get<word>("frozenPointsZone")
        )
      : -1
    ),
    patchID_
    (
        mesh.boundaryMesh().findPatchID
        (
            word(motionSolver::coeffDict().lookup("patchName"))
        )
    ),
    volumeMin_(gSum(mesh.cellVolumes())),
    useNewtonRaphson_
    (
        motionSolver::coeffDict().lookupOrDefault<Switch>
        (
            "useNewtonRaphson", true
        )
    ),
    relTol_
    (
        motionSolver::coeffDict().lookupOrDefault<scalar>("relTol", 1e-3)
    ),
    times_(),
    flowRates_(),
    volumes_(),
    initialPatchPoints_(),
    initialPatchPointNormals_(),
    centrePoint_(point::zero),
    planeP1_(point::zero),
    planeP2_(point::zero),
    planeP3_(point::zero),
    planeNormal_(vector::zero),
    maxNormalDistance_(0.0),
    radialExponent_(1.0),
    radialProfile_("apexWeighted"),
    peakLocation_(0.5),
    peakSharpness_(1.0),
    minProfileWeight_(0.0),
    apexProfileWeight_(0.3),
    dispScale_(1e-2)
{
    initialise(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::movingVentricleRadialFvMotionSolver::~movingVentricleRadialFvMotionSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::motionDiffusivity&
Foam::movingVentricleRadialFvMotionSolver::diffusivity()
{
    if (!diffusivityPtr_)
    {
        diffusivityPtr_ = motionDiffusivity::New
        (
            fvMesh_,
            coeffDict().lookup("diffusivity")
        );
    }

    return *diffusivityPtr_;
}


Foam::tmp<Foam::pointField>
Foam::movingVentricleRadialFvMotionSolver::curPoints() const
{
    interpolationPtr_->interpolate
    (
        cellDisplacement_,
        pointDisplacement_
    );

    tmp<pointField> tcurPoints
    (
        points0() + pointDisplacement_.primitiveField()
    );
    pointField& curPoints = tcurPoints.ref();

    // Implement frozen points
    if (frozenPointsZone_ != -1)
    {
        const pointZone& pz = fvMesh_.pointZones()[frozenPointsZone_];

        forAll(pz, i)
        {
            curPoints[pz[i]] = points0()[pz[i]];
        }
    }

    twoDCorrectPoints(curPoints);

    return tcurPoints;
}


void Foam::movingVentricleRadialFvMotionSolver::solve()
{
    // The points have moved so before interpolation update
    // the motionSolver accordingly
    movePoints(fvMesh_.points());

    const Time& runTime = fvMesh_.time();

    // Take a reference to the patch and its point indices
    const polyPatch& ventriclePatch = fvMesh_.boundaryMesh()[patchID_];
    const labelList& ventricleMeshPoints = ventriclePatch.meshPoints();

    // Lookup the current volume
    const scalar targetVolume = linearInterp(times_, volumes_, runTime.value());

    // Determine the displacement scaling factor to achieve the target volume
    if (!useNewtonRaphson_)
    {
        FatalErrorInFunction
            << "useNewtonRaphson=false but this solver was requested to run "
            << "without bisection fallback."
            << exit(FatalError);
    }

    // ---------- Helper: build radial patch displacement given scale factor a ----------
    auto makePatchDisp = [&](scalar a) -> vectorField
    {
        vectorField patchPointDisp(ventriclePatch.nPoints(), vector::zero);

        // global displacement amplitude
        const scalar sGlobal = dispScale_*a;

        forAll(patchPointDisp, pI)
        {
            const point& p0 = initialPatchPoints_[pI];

            // signed distance along the user-defined normal direction
            const scalar s = (p0 - centrePoint_) & planeNormal_;

            scalar t = s/max(VSMALL, maxNormalDistance_);
            t = max(0.0, min(1.0, t));

            const scalar w = profileWeight(t);

            // radial direction from centrePoint
            vector r = (p0 - centrePoint_);
            const scalar rMag = mag(r);

            if (rMag > VSMALL)
            {
                patchPointDisp[pI] = (r/rMag) * (sGlobal*w);
            }
        }

        return patchPointDisp;
    };

    // Apply scale 'a', enforce BCs, return volume
    auto evalVolume = [&](scalar a) -> scalar
    {
        vectorField patchPointDisp = makePatchDisp(a);


        // Push to boundary, enforce continuity across procs
        pointDisplacement_.boundaryFieldRef()[patchID_] == patchPointDisp;
        pointDisplacement_.correctBoundaryConditions();
        patchPointDisp =
            pointDisplacement_.boundaryFieldRef()
            [
                patchID_
            ].patchInternalField();

        return
            calculateVolume
            (
                fvMesh_, points0(), ventricleMeshPoints, patchPointDisp
            );
    };

    // Newton-Raphosn using forward finite differencing to approximate the
    // deriviative
    scalar a = 0.0;
    scalar V = evalVolume(a);

    const scalar tol = relTol_*targetVolume;
    const label  maxIter = 20;
    label iter = 0;
    scalar res = mag(V - targetVolume);

    Info<< "Find the displacement scale factor using the Newton-Raphson"
        << " method" << nl
        << "    Target volume = " << targetVolume << endl;

    while (res > tol && ++iter < maxIter)
    {
        // Forward-difference slope: dV/da ≈ (V(a+h) - V(a)) / h
        const scalar h = max(1e-8*max(1.0, mag(a)), 1e-8);
        const scalar Vp = evalVolume(a + h);
        const scalar dV = (Vp - V) / h;

        if (mag(dV) < VSMALL)
        {
            Info<< "  Newton: tiny slope, stopping.\n";
            break;
        }

        scalar aNew = a - (V - targetVolume)/dV;
            
        // Light damping: ensure residual decreases; at most 5 halvings
        scalar Vnew = evalVolume(aNew);
        label ls = 0;
        while
        (
            mag(Vnew - targetVolume) > 0.9*mag(V - targetVolume)
         && ls++ < 5
        )
        {
            aNew = 0.5*(a + aNew);
            Vnew = evalVolume(aNew);
        }

        a = aNew;
        V = Vnew;
        
        // Update the residual
        res = mag(V - targetVolume);

        // Print info
        Info<< "    " << iter
            << ": scale factor = " << a
            << ", volume = " << V
            << ", |res| = " << res << nl;
    }

    if (iter == maxIter)
    {
        FatalErrorInFunction
            << "Max iterations reached in the Newton-Raphson loop!"
            << exit(FatalError);
    }

    // Final: keep the converged displacement in patchPointDisp
    // Final: set converged displacement on the patch (and keep it)
    {
        vectorField patchPointDisp = makePatchDisp(a);

        pointDisplacement_.boundaryFieldRef()[patchID_] == patchPointDisp;
        pointDisplacement_.correctBoundaryConditions();
    }

    // Update the diffusivity field
    diffusivity().correct();

    // Update the pointDisplacement boundary conditions: these will be used
    // by the cellDisplacement boundary conditions
    pointDisplacement_.boundaryFieldRef().updateCoeffs();

    // Solve for cellDisplacement
    Info<< "Solving the mesh motion for cellDisplacement" << endl;
    fvVectorMatrix DEqn
    (
        fvm::laplacian
        (
            dimensionedScalar("viscosity", dimViscosity, 1.0)
           *diffusivity().operator()(),
            cellDisplacement_,
            "laplacian(diffusivity,cellDisplacement)"
        )
    );

    DEqn.solveSegregatedOrCoupled();
}


void Foam::movingVentricleRadialFvMotionSolver::updateMesh
(
    const mapPolyMesh& mpm
)
{
    displacementMotionSolver::updateMesh(mpm);
}


// ************************************************************************* //
