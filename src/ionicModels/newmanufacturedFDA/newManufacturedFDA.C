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
#include "newManufacturedFDA.H"
#include "newManufacturedFDA_2014.H"
#include "addToRunTimeSelectionTable.H"
#include "newionicModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(newmanufacturedFDA, 0);
    addToRunTimeSelectionTable
    (
        newionicModel, newmanufacturedFDA, dictionary
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::newmanufacturedFDA::newmanufacturedFDA
(
    const dictionary& dict,
    const label num,
    const scalar initialDeltaT,
    const Switch solveVmWithinODESolver
)
:
    newionicModel(dict, num, initialDeltaT, solveVmWithinODESolver),
    STATES_(num),
    STATES_OLD_(num),
    CONSTANTS_(47, 0.0),
    ALGEBRAIC_(num),
    RATES_(num)
{

    //see if I need to add flog in function as well.
    Info<< nl << "Calling TNNP test Constants" << endl;
    forAll(STATES_, i)
    {
        STATES_.set(i, new scalarField(17, 0.0));
        STATES_OLD_.set(i, new scalarField(17, 0.0));
        ALGEBRAIC_.set(i, new scalarField(69, 0.0));
        RATES_.set(i, new scalarField(17, 0.0));



        // Initialise the constants (repeatedly! it's ok...) and the rates and
        // states
        newManufacturedFDAinitConsts
        (
            CONSTANTS_.data(),
            RATES_[i].data(),
            STATES_[i].data(),
            tissue()
        );

        // Initialise old time STATES
        STATES_OLD_[i] = STATES_[i];
    }

    Info<< nl
        << "CONSTANTS = " << CONSTANTS_ << endl
        << "STATES = " << STATES_[0] << endl;
}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::newmanufacturedFDA::~newmanufacturedFDA()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::List<Foam::word> Foam::newmanufacturedFDA::supportedTissues() const
{
    return {"epicardialCells", "mCells", "endocardialCells"};
}


void Foam::newmanufacturedFDA::calculateCurrent
(
    const scalar stepStartTime,
    const scalar deltaT,
    const scalarField& Vm,
    scalarField& Iion,
    Field<Field<double>>& STATESext

    
)
{
    const scalar tStart = stepStartTime*1000;
    const scalar tEnd = (stepStartTime + deltaT)*1000;
    label monitorCell = 0; 
    const label n = STATES_.size();

    if (Vm.size() != n ||  Iion.size() != n || STATESext.size() != n)
    {
        FatalErrorInFunction
            << "Wrong field sizes" << abort(FatalError);
    }

    STATES_ = STATES_OLD_;

    forAll(STATES_, cellI)
    {
        scalarField& STATESI    = STATES_[cellI];
        scalarField& ALGEBRAICI = ALGEBRAIC_[cellI];
        scalarField& RATESI     = RATES_[cellI];
        Field<double>& extS     = STATESext[cellI];

        
        

         // Set Vm
        // Overwrite Vm (always state 0)
        STATESI[0] = Vm[cellI]*1000;

        newManufacturedFDAcomputeVariables
        (
            stepStartTime,
            CONSTANTS_.data(),
            RATESI.data(),
            STATESI.data(),
            ALGEBRAICI.data(),
            tissue(),
            solveVmWithinODESolver()
        );

        Iion[cellI] = ALGEBRAICI[50] + ALGEBRAICI[57] + ALGEBRAICI[51] + ALGEBRAICI[52]
        + ALGEBRAICI[55] + ALGEBRAICI[58] + ALGEBRAICI[53] + ALGEBRAICI[54]
        + ALGEBRAICI[59] + ALGEBRAICI[56] + ALGEBRAICI[61] + ALGEBRAICI[60];

        if (cellI == monitorCell)
            {
            Info<< "integrationPtI = " << cellI
                    << " | t = " << tStart
                    << " → " << tEnd
                    << " | Vm = " << STATESI[0]
                    << " | Iion = " << ALGEBRAICI[50] + ALGEBRAICI[57] + ALGEBRAICI[51] + ALGEBRAICI[52]
                    + ALGEBRAICI[55] + ALGEBRAICI[58] + ALGEBRAICI[53] + ALGEBRAICI[54]
                    + ALGEBRAICI[59] + ALGEBRAICI[56] + ALGEBRAICI[61] + ALGEBRAICI[60]
                    << endl;
            }

        // Copy back to STATESext
        for (label s = 0; s < 17; ++s)
        {
            extS[s] = STATESI[s];
        }
    }
}


void Foam::newmanufacturedFDA::calculateGating
(
    const scalar stepStartTime,
    const scalar deltaT,
    const scalarField& Vm,
    scalarField& Iion,
    Field<Field<double>>& STATESext
)
{
    // stepStartTime in seconds → milliseconds
    const scalar tStart = stepStartTime*1000;

    // Start from last accepted states
    STATES_ = STATES_OLD_;

    const label monitorCell = 0;

    forAll(STATES_, cellI)
    {
        scalarField& STATESI    = STATES_[cellI];
        scalarField& ALGEBRAICI = ALGEBRAIC_[cellI];
        scalarField& RATESI     = RATES_[cellI];
        Field<double>& extS     = STATESext[cellI];

        // Always set Vm in mV
        STATESI[0] = Vm[cellI]*1000;

        // Compute rates at tStart (ms)
        newManufacturedFDAcomputeRates
        (
            tStart,
            CONSTANTS_.data(),
            RATESI.data(),
            STATESI.data(),
            ALGEBRAICI.data(),
            tissue(),
            solveVmWithinODESolver()
        );

        // Explicit Euler in time:
        // deltaT [s] * 1000 [ms/s] * d/dt[per ms] → dimensionless
        for (label s = 1; s < 17; ++s)
        {
            STATESI[s] += deltaT * RATESI[s];
        }

        if (cellI == monitorCell && debug)
        {
            Info<< "[calculateGating] cell=" << cellI
                << " t(ms)=" << tStart
                << " Vm(mV)=" << STATESI[0]
                << " m=" << STATESI[7]
                << " h=" << STATESI[8]
                << " j=" << STATESI[9]
                << " Cai=" << STATESI[3]
                << nl;
        }

        // Optional: recompute Iion with updated states (if needed)
        Iion[cellI] =
            ALGEBRAICI[50] + ALGEBRAICI[57] + ALGEBRAICI[51] + ALGEBRAICI[52]
          + ALGEBRAICI[55] + ALGEBRAICI[58] + ALGEBRAICI[53] + ALGEBRAICI[54]
          + ALGEBRAICI[59] + ALGEBRAICI[56] + ALGEBRAICI[61] + ALGEBRAICI[60];

        // Copy back to external STATES
        for (label s = 0; s < 17; ++s)
        {
            extS[s] = STATESI[s];
        }
    }

    STATES_OLD_ = STATES_;
}

void Foam::newmanufacturedFDA::solveODE
(
    const scalar stepStartTime,
    const scalar deltaT,
    const scalarField& Vm,
    scalarField& Iion,
    Field<Field<double>>& STATESext
)
{
    // Reset from last accepted time level
    STATES_ = STATES_OLD_;

    const scalar tStart = stepStartTime*1000;          // ms
    const scalar tEnd   = (stepStartTime + deltaT)*1000; // ms
    const label monitorCell = 0;

    forAll(STATES_, cellI)
    {
        scalarField& STATESI    = STATES_[cellI];
        scalarField& ALGEBRAICI = ALGEBRAIC_[cellI];
        scalarField& RATESI     = RATES_[cellI];
        Field<double>& extS     = STATESext[cellI];

        scalar& step = newionicModel::step()[cellI];
        step = min(step, deltaT*1000);    // local ODE step in ms

        // Set Vm in mV
        STATESI[0] = Vm[cellI]*1000;

        if (cellI == monitorCell && debug)
        {
            Info<< "[solveODE] BEFORE solve: cell=" << cellI
                << " tStart(ms)=" << tStart
                << " Vm(mV)=" << STATESI[0]
                << " m=" << STATESI[7]
                << " h=" << STATESI[8]
                << nl;
        }

        // Advance ODE system in time
        odeSolver().solve(tStart, tEnd, STATESI, step);

        // Update algebraic variables & Iion at new time level (tEnd)
        newManufacturedFDAcomputeVariables
        (
            tEnd,
            CONSTANTS_.data(),
            RATESI.data(),
            STATESI.data(),
            ALGEBRAICI.data(),
            tissue(),
            solveVmWithinODESolver()
        );

        Iion[cellI] =
            ALGEBRAICI[50] + ALGEBRAICI[57] + ALGEBRAICI[51] + ALGEBRAICI[52]
          + ALGEBRAICI[55] + ALGEBRAICI[58] + ALGEBRAICI[53] + ALGEBRAICI[54]
          + ALGEBRAICI[59] + ALGEBRAICI[56] + ALGEBRAICI[61] + ALGEBRAICI[60];

        if (cellI == monitorCell && debug)
        {
            Info<< "[solveODE] AFTER solve: cell=" << cellI
                << " tEnd(ms)=" << tEnd
                << " Vm(mV)=" << STATESI[0]
                << " m=" << STATESI[7]
                << " h=" << STATESI[8]
                << " Iion=" << Iion[cellI]
                << nl;
        }

        // Copy updated states back out
        for (label s = 0; s < 17; ++s)
        {
            extS[s] = STATESI[s];
        }
    }

    STATES_OLD_ = STATES_;
}


//RIGHT-HAND SIDE: DERIVATIVES

void Foam::newmanufacturedFDA::derivatives
(
    const scalar t,
    const scalarField& y,
    scalarField& dydt
) const
{
    scalarField ALGEBRAIC_TMP(69, 0.0);

    ::newManufacturedFDAcomputeRates
    (
        t,
        CONSTANTS_.data(),
        dydt.data(),
        const_cast<scalarField&>(y).data(),
        ALGEBRAIC_TMP.data(),
        tissue(),
        solveVmWithinODESolver() 
    );
}




// * * * * * * * * * * * * * * Member Functions - postProcessing  * * * * * * * * * * * * * * //

void Foam::newmanufacturedFDA::writeHeader(OFstream& output) const
{

    output << "time Vm";

    for (int i = 0; i < 17; ++i)
        output << " " << newmanufacturedFDASTATES_NAMES[i];

    for (int i = 0; i < 69; ++i)
        output << " " << newmanufacturedFDAALGEBRAIC_NAMES[i];

    for (int i = 0; i < 17; ++i)
        output << " RATES_" << newmanufacturedFDASTATES_NAMES[i];

    output << endl;
}


void Foam::newmanufacturedFDA::write(const scalar t, OFstream& output) const

{

    scalar Vm = STATES_[0][0];

    output
        << t << " " << Vm;

    // States
    forAll(STATES_[0], j)
    {
        output << " " << STATES_[0][j];
    }

    // Algebraic variables
    forAll(ALGEBRAIC_[0], j)
    {
        output << " " << ALGEBRAIC_[0][j];
    }

    // Rates
    forAll(RATES_[0], j)
    {
        output << " " << RATES_[0][j];
    }


output << endl;
}



// ************************************************************************* //