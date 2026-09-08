/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
-------------------------------------------------------------------------------
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

Class
    Foam::movingVentricleRadialFvMotionSolver

Group
    grpMeshMotionSolvers

Description
    Mesh motion solver for an fvMesh for a moving ventricle case, where the mesh
    is moved such that the ventricle volume follows a prescribed time curve.

SourceFiles
    movingVentricleRadialFvMotionSolver.C

Author
    Philip Cardiff, UCD.
    Chanikya Valeti, UCD.
\*---------------------------------------------------------------------------*/

#ifndef movingVentricleRadialFvMotionSolver_H
#define movingVentricleRadialFvMotionSolver_H

#include "displacementMotionSolver.H"
#include "fvMotionSolver.H"
#include "motionInterpolation.H"
#include "motionDiffusivity.H"
#include "interpolationTable.H"
#include "fvMesh.H"
#include "pointFields.H"
#include "volFields.H"

namespace Foam
{
class motionInterpolation;
class motionDiffusivity;
/*---------------------------------------------------------------------------*\
             Class movingVentricleRadialFvMotionSolver Declaration
\*---------------------------------------------------------------------------*/

class movingVentricleRadialFvMotionSolver
:
    public displacementMotionSolver,
    public fvMotionSolver
{
    // Private data

        //- Cell-centre motion field
        mutable volVectorField cellDisplacement_;

        //- Interpolation used to transfer cell displacement to the points
        autoPtr<motionInterpolation> interpolationPtr_;

        //- Diffusivity used to control the motion
        autoPtr<motionDiffusivity> diffusivityPtr_;

        //- Frozen points (that are not on patches). -1 or points that are
        //  fixed to be at points0_ location
        label frozenPointsZone_;

        //- Patch index ventricle moving patch
        const label patchID_;

        //- Minimum volume (assumed to be the starting volume)
        const scalar volumeMin_;

        //- Use Newton-Raphson to determine the displacement scale factor
        const Switch useNewtonRaphson_;

        //- Relative tolerance on the volume calculation loop
        const scalar relTol_;

        //- Times from flow rate series
        scalarField times_;

        //- Flow rate series
        scalarField flowRates_;

        //- Target volumes corresponding to times
        scalarField volumes_;

        //- Store the initial patch points
        vectorField initialPatchPoints_;

        //- Store the initial patch point unit normals
        vectorField initialPatchPointNormals_;

        // Radial motion controls
        point centrePoint_;

        // Explicit base plane from three points
        // Set the normal direction by the right-hand rule.
        point planeP1_;
        point planeP2_;
        point planeP3_;
        vector planeNormal_;

        // Maximum signed distance of moving patch points from centrePoint
        // along planeNormal_, used for normalization
        scalar maxNormalDistance_;

        scalar radialExponent_;

        // Radial weighting profile along the normalized ventricle axis
        word radialProfile_;

        // Peak location for non-monotonic weighting profiles
        scalar peakLocation_;

        // Sharpness control for peaked weighting profiles
        scalar peakSharpness_;

        // Minimum normalized deformation weight applied everywhere
        scalar minProfileWeight_;

        // Apex weight relative to the peak for asymmetric peaked profiles
        scalar apexProfileWeight_;

        // Overall displacement scale multiplier
        scalar dispScale_;


    // Private Member Functions

        //- Initialisation steps given a dictionary
        void initialise(const dictionary& dict);

        //- Linear interpolation helper
        scalar linearInterp
        (
            const scalarField& x,
            const scalarField& y,
            const scalar xq
        ) const;

        //- Evaluate the radial weighting profile for a normalized coordinate
        scalar profileWeight(const scalar t) const;

        //- Calculate the domain volume given the patch displacement
        scalar calculateVolume
        (
            const fvMesh& mesh,
            const pointField& initialPoints,
            const labelList& ventricleMeshPoints,
            const vectorField& displacement
        ) const;

        //- No copy construct
        movingVentricleRadialFvMotionSolver
        (
            const movingVentricleRadialFvMotionSolver&
        ) = delete;

        //- No copy assignment
        void operator=(const movingVentricleRadialFvMotionSolver&) = delete;


public:

    //- Runtime type information
    TypeName("movingVentricleRadialFvMotionSolver");


    // Constructors

        //- Construct from polyMesh and IOdictionary
        movingVentricleRadialFvMotionSolver
        (
            const polyMesh&,
            const IOdictionary&
        );

        //- Construct from components
        movingVentricleRadialFvMotionSolver
        (
            const polyMesh& mesh,
            const IOdictionary& dict,
            const pointVectorField& pointDisplacement,
            const pointIOField& points0
        );


    //- Destructor
    ~movingVentricleRadialFvMotionSolver();


    // Member Functions

        //- Return reference to the diffusivity field
        motionDiffusivity& diffusivity();

        //- Return point location obtained from the current motion field
        virtual tmp<pointField> curPoints() const;

        //- Solve for motion
        virtual void solve();

        //- Update topology
        virtual void updateMesh(const mapPolyMesh&);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
