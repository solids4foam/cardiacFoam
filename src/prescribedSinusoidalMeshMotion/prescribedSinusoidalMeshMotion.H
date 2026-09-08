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

Class
    Foam::prescribedSinusoidalMeshMotion

Group
    grpMeshMotionSolvers

Description
    Mesh motion solver which prescribes a sinusoidal motion, in space and in
    time, over a user-specified box in the x-y plane. The purpose of this
    solver is for testing the performance of a moving mesh fluid solver.

    The point displacement is given in terms of the normalised coordinates

        xi  = (x - xMin)/(xMax - xMin)
        eta = (y - yMin)/(yMax - yMin)

    and the time factor

        f(t) = scaleFactor*sin(2*pi*t/period)

    as

        d = f(t)*sin(pi*xi)*sin(pi*eta)*(1 1 0)

    or, when orthogonalMeshMotion is enabled, as the orthogonal motion

        d = f(t)*sin(pi*xi)*(1 0 0) + f(t)*sin(pi*eta)*(0 1 0)

    All of scaleFactor, orthogonalMeshMotion, xMin, xMax, yMin, yMax and
    period are required entries; none has a default.

SourceFiles
    prescribedSinusoidalMeshMotion.C

Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#ifndef prescribedSinusoidalMeshMotion_H
#define prescribedSinusoidalMeshMotion_H

#include "displacementMotionSolver.H"
#include "fvMotionSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward class declarations
class motionInterpolation;
class motionDiffusivity;

/*---------------------------------------------------------------------------*\
             Class prescribedSinusoidalMeshMotion Declaration
\*---------------------------------------------------------------------------*/

class prescribedSinusoidalMeshMotion
:
    public displacementMotionSolver,
    public fvMotionSolver
{
    // Private data

        //- Motion scale factor
        const scalar A_;

        //- Flag to indicate orthogonal mesh motion or more general motion
        const Switch orthogonalMeshMotion_;

          //- Channel bounds and period
        const scalar xMin_;
        const scalar xMax_;
        const scalar yMin_;
        const scalar yMax_;
        const scalar period_;


    // Private Member Functions

        //- No copy construct
        prescribedSinusoidalMeshMotion
        (
            const prescribedSinusoidalMeshMotion&
        ) = delete;

        //- No copy assignment
        void operator=(const prescribedSinusoidalMeshMotion&) = delete;


public:

    //- Runtime type information
    TypeName("prescribedSinusoidalMeshMotion");


    // Constructors

        //- Construct from polyMesh and IOdictionary
        prescribedSinusoidalMeshMotion
        (
            const polyMesh&,
            const IOdictionary&
        );

        //- Construct from components
        prescribedSinusoidalMeshMotion
        (
            const polyMesh& mesh,
            const IOdictionary& dict,
            const pointVectorField& pointDisplacement,
            const pointIOField& points0
        );


    //- Destructor
    ~prescribedSinusoidalMeshMotion();


    // Member Functions

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
