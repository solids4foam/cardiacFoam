/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2025 AUTHOR,AFFILIATION
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

#include <math.h>
#include "tentusscher_noble_noble_panfilov_2004.H"
#include "ionicModelCellML.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(ionicModelCellML, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ionicModelCellML::ionicModelCellML(const dictionary& dict)
:
    ODESystem(),
    odeSolver_(ODESolver::New(*this, dict)),
    STATES_(17, 0.0),
    CONSTANTS_(46, 0.0),
    ALGEBRAIC_(69, 0.0),
    RATES_(17, 0.0)
{
    Info<< nl << "Calling initConsts" << endl;
    initConsts(CONSTANTS_.data(), RATES_.data(), STATES_.data());

    Info<< nl
        << "CONSTANTS = " << CONSTANTS_ << nl
        << "STATES = " << STATES_ << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ionicModelCellML::~ionicModelCellML()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


void Foam::ionicModelCellML::solve
(
    const scalar tOld, const scalar deltaT
) const
{
    // Define the start time, end time and initial time step
    const scalar tStart = tOld;
    const scalar tEnd = tOld + deltaT;
    scalar step = deltaT;

    // Take a reference to the initial state
    scalarField& yStart = STATES_;

    // Solve the ODE system
    odeSolver_->solve(tStart, tEnd, yStart, step);
}


void Foam::ionicModelCellML::computeVariables(const scalar t) const
{
    ::computeVariables
    (
        t,
        CONSTANTS_.data(),
        RATES_.data(),
        STATES_.data(),
        ALGEBRAIC_.data()
    );
}


void Foam::ionicModelCellML::write(const scalar t, OFstream& output) const
{
    // Write the results
    Info<< "Writing to " << output.name() << endl;
    output
        << t;
    forAll(STATES_, i)
    {
        output
            << " " << STATES_[i];
    }
    forAll(ALGEBRAIC_, i)
    {
        output
            << " " << ALGEBRAIC_[i];
    }
    forAll(RATES_, i)
    {
        output
            << " " << RATES_[i];
    }
    output
        << endl;
}


void Foam::ionicModelCellML::derivatives
(
    const scalar t,
    const scalarField& y,
    scalarField& dydt
) const
{
    // VOI is the curren time
    const scalar VOI = t;

    // Transfer y to STATES
    STATES_ = y;

    // Calculate the rates using the cellML header file
    computeRates
    (
        VOI,
        CONSTANTS_.data(),
        RATES_.data(),
        STATES_.data(),
        ALGEBRAIC_.data()
    );

    // Transfer RATES to dydt
    dydt = RATES_;
}


void Foam::ionicModelCellML::jacobian
(
    const scalar t,
    const scalarField& y,
    scalarField& dfdt,
    scalarSquareMatrix& dfdy
) const
{
    notImplemented("void jacobian(...)");
}

// ************************************************************************* //
