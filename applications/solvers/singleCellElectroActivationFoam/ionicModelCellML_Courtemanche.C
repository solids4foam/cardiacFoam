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
#include "Courtemanche.H"
#include "ionicModelCellML.H"
#include <string>


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
    STATES_(NUM_STATES, 0.0),
    CONSTANTS_(NUM_CONSTANTS, 0.0),
    ALGEBRAIC_(NUM_ALGEBRAIC, 0.0),
    RATES_(NUM_STATES, 0.0),
    step_(readScalar(dict.lookup("initialODEStep")))

    
{ 
    word tissueName;
    dict.lookup("tissue") >> tissueName;

    Info << "Tissue Name " << tissueName << endl;
    tissue_ = (tissueName == "epicardialCells") 
                    ? 1 
            : (tissueName == "mCells") 
                    ? 2            
            : (tissueName == "endocardialCells")    
                    ? 3
            : -1;  //invalid flag 

    if (tissue_ == -1)
    {
        FatalErrorInFunction
            << "Unknown tissue: " << tissueName
            << nl << exit(FatalError);

    }

    Info << "Tissue flag set to: " << tissue_ << endl;
    //see if I need to add flog in function as well.
    Info<< nl << "Calling initConsts" << endl;
    CourtemancheinitConsts(CONSTANTS_.data(), RATES_.data(), STATES_.data(), tissue_);

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

    // Take a reference to the initial state
    scalarField& yStart = STATES_;

    scalar step = min(step_, deltaT);

    // Solve the ODE system
    Info<< "step = " << step_ << endl;
    odeSolver_->solve(tStart, tEnd, yStart, step);

    // Update step
    step_ = step;
}


void Foam::ionicModelCellML::computeVariables(const scalar t) const
{
    ::CourtemanchecomputeVariables
    (
        t,
        CONSTANTS_.data(),
        RATES_.data(),
        STATES_.data(),
        ALGEBRAIC_.data(),
        tissue_
    );
}

void Foam::ionicModelCellML::writeHeader(OFstream& output) const
{
    output << "time";

    for (int i = 0; i < NUM_STATES; ++i)
    {
        output << " " << STATES_NAMES[i];
    }
    for (int i = 0; i < NUM_ALGEBRAIC; ++i)
    {
        output << " " << ALGEBRAIC_NAMES[i];
    }
    for (int i = 0; i < NUM_STATES; ++i)
    {
        output << "RATES_" <<STATES_NAMES[i];
    }

    output 
        << endl;
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
    ::CourtemanchecomputeVariables
    (
        VOI,
        CONSTANTS_.data(),
        RATES_.data(),
        STATES_.data(),
        ALGEBRAIC_.data(),
        tissue_
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
