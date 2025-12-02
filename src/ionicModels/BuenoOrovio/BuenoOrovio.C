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
#include "BuenoOrovio.H"
#include "BuenoOrovio_2008.H"
#include "addToRunTimeSelectionTable.H"
#include "ionicModel.H"
#include "ionicModelIO.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(BuenoOrovio, 0);
    addToRunTimeSelectionTable
    (
        ionicModel, BuenoOrovio, dictionary
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::BuenoOrovio::BuenoOrovio
(
    const dictionary& dict,
    const label num,
    const scalar initialDeltaT,
    const Switch solveVmWithinODESolver
)
:
    ionicModel(dict, num, initialDeltaT, solveVmWithinODESolver),
    STATES_(num),
    STATES_OLD_(num),
    CONSTANTS_(NUM_CONSTANTS, 0.0),
    ALGEBRAIC_(num),
    RATES_(num)
{
    Info<< nl << "Initialize Bueno Orovio constants:" << nl;
    forAll(STATES_, i)
    {
        STATES_.set(i,      new scalarField(NUM_STATES,     0.0));
        STATES_OLD_.set(i,  new scalarField(NUM_STATES,     0.0));
        ALGEBRAIC_.set(i,   new scalarField(NUM_ALGEBRAIC,  0.0));
        RATES_.set(i,       new scalarField(NUM_STATES,     0.0));

        // Initialise constants, states and rates from generated code
        BuenoOrovioinitConsts
        (
            CONSTANTS_.data(),
            RATES_[i].data(), 
            STATES_[i].data(),
            tissue(),dict
        );

        ionicModelIO::loadStimulusConstants
        (
            dict, CONSTANTS_,
            stim_start, stim_period_S1,
            stim_duration, stim_amplitude,
            nstim1, stim_period_S2, nstim2
        );
    }
    Info<< CONSTANTS_ << nl;

    label i0 = rand() % STATES_.size();
    Info<< "initial states:" << nl;
    Info<< STATES_[i0] << nl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::BuenoOrovio::~BuenoOrovio()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::List<Foam::word> Foam::BuenoOrovio::supportedTissueTypes() const
{
    // All three tissue variants are supported in the generated code
    return {"endocardialCells", "mCells", "epicardialCells"};
}


// ------------------------------------------------------------------------- //
//  Explicit split: calculateCurrent (Iion only, no state update)
// ------------------------------------------------------------------------- //

void Foam::BuenoOrovio::calculateCurrent
(
    const scalar stepStartTime,
    const scalar deltaT,
    const scalarField& Vm,
    scalarField& Im,
    Field<Field<scalar>>& states
)
{
    // stepStartTime, deltaT are in seconds; BuenoOrovio uses ms
    const scalar tStart = stepStartTime * 1000.0;
    const label monitorCell = 500;

    if (Im.size() != Vm.size())
    {
        FatalErrorInFunction
            << "Im.size() != Vm.size()" << abort(FatalError);
    }

    // We do NOT modify the gating states here – just compute Iion
    forAll(STATES_, integrationPtI)
    {
        scalarField& STATESI    = STATES_[integrationPtI];
        scalarField& ALGEBRAICI = ALGEBRAIC_[integrationPtI];
        scalarField& RATESI     = RATES_[integrationPtI];

        // Vm in the CellML code is in mV
        STATESI[0] = (Vm[integrationPtI] * 1000.0 + 84)/85.7;

        if (integrationPtI == monitorCell)
        {
            Info<< "calculateCurrent: Vm[mV] before computeVariables = "
                << STATESI[0] << nl;
        }

        ::BuenoOroviocomputeVariables
        (
            tStart,
            CONSTANTS_.data(),
            RATESI.data(),
            STATESI.data(),
            ALGEBRAICI.data(),
            tissue(),
            solveVmWithinODESolver()
        );

        if (integrationPtI == monitorCell)
        {
            Info<< "calculateCurrent: iPt = " << integrationPtI
                << " | t = " << tStart
                << " | Vm = " << (STATESI[0] * 85.7 - 84)
                << " | u = " << STATESI[u]
                << " | v = " << STATESI[v]
                << " | w = " << STATESI[w]
                << " | s= " << STATESI[s]
                << " | Jion = " << ALGEBRAICI[Jion]
                << nl;
        }

        // Jion  is the total ionic current density used by the PDE
        Im[integrationPtI] = ALGEBRAICI[Jion] * 85.7;

        // Optionally export internal states to the external buffer
        if (states[integrationPtI].size() >= NUM_STATES)
        {
            for (label s = 0; s < NUM_STATES; ++s)
            {
                states[integrationPtI][s] = STATESI[s];
            }
        }
    }
}

// ------------------------------------------------------------------------- //
//  Implicit path: solve ODE system via OpenFOAM ODESolver
//  (used in the implicit Vm branch of newelectroFoam)
// ------------------------------------------------------------------------- //

void Foam::BuenoOrovio::solveODE
(
    const scalar stepStartTime,
    const scalar deltaT,
    const scalarField& Vm,
    scalarField& Im,
    Field<Field<scalar>>& states
)
{
    const scalar tStart = stepStartTime * 1000.0;
    const scalar tEnd   = (stepStartTime + deltaT) * 1000.0;
    const label monitorCell = 0;

    forAll(STATES_, integrationPtI)
    {
        scalarField& STATESI    = STATES_[integrationPtI];
        scalarField& ALGEBRAICI = ALGEBRAIC_[integrationPtI];
        scalarField& RATESI     = RATES_[integrationPtI];

        
        // Vm fed into the cell model in mV
        if (!solveVmWithinODESolver())
        {
            STATESI[0] = (Vm[integrationPtI] * 1000.0 + 84)/85.7;
        }
        // Per-cell adaptive time step (in ms) for the ODE solver
        scalar& h = ionicModel::step()[integrationPtI];

        // Clamp ODE step
        h = min(h, deltaT * 1000.0);

        if (integrationPtI == monitorCell)
        {
            Info<< "solveODE: i=" << integrationPtI
                << " | t = " << tStart << " → " << tEnd
                << " | step = " << h
                << " | Vm = " << (STATESI[0] * 85.7 - 84)
                << nl;
        }

        // Advance the ODE system
        odeSolver().solve(tStart, tEnd, STATESI, h);

        // Update ALGEBRAIC (incl. Jion) and RATES at tEnd
        ::BuenoOroviocomputeVariables
        (
            tEnd,
            CONSTANTS_.data(),
            RATESI.data(),
            STATESI.data(),
            ALGEBRAICI.data(),
            tissue(),
            solveVmWithinODESolver()
        );

        if (integrationPtI == monitorCell)
        {
            Info<< "solveODE: i=" << integrationPtI
                << " | Vm = " << STATESI[0] * 85.7 - 84
                << " | Jion = " << ALGEBRAICI[Jion]
                << " | dVdt = " << RATESI[0]
                << nl;
        }

        // Total ionic current density used by PDE
        Im[integrationPtI] = ALGEBRAICI[Jion] * 85.7;

        // Export states if requested
        if (states[integrationPtI].size() >= NUM_STATES)
        {
            for (label s = 0; s < NUM_STATES; ++s)
            {
                states[integrationPtI][s] = STATESI[s];
            }
        }
    }
}

// ------------------------------------------------------------------------- //
//  ODE RHS (used by OpenFOAM ODE solvers)
// ------------------------------------------------------------------------- //

void Foam::BuenoOrovio::derivatives
(
    const scalar t,
    const scalarField& y,
    scalarField& dydt
) const
{
    // Must match NUM_ALGEBRAIC from the generated BuenoOrovio code
    scalarField ALGEBRAIC_TMP(NUM_ALGEBRAIC, 0.0);

    ::BuenoOroviocomputeVariables
    (
        t,
        CONSTANTS_.data(),
        dydt.data(),                              // RATES (output)
        const_cast<scalarField&>(y).data(),       // STATES (input)
        ALGEBRAIC_TMP.data(),                     // ALGEBRAIC (scratch)
        tissue(),
        solveVmWithinODESolver()
    );
}

// ------------------------------------------------------------------------- //
//  Writing logic in singleCell and 3D simulations
// ------------------------------------------------------------------------- //
Foam::wordList Foam::BuenoOrovio::exportedFieldNames() const
    {
        return ionicModelIO::exportedFieldNames
        (
            variableExport_,
            BuenoOrovioSTATES_NAMES, NUM_STATES,
            BuenoOrovioALGEBRAIC_NAMES, NUM_ALGEBRAIC
        );
    } 

void Foam::BuenoOrovio::exportStates
(
    const Field<Field<scalar>>&,
    List<volScalarField*>& outFields
)
{
    ionicModelIO::exportStateFields
    (
        STATES_,ALGEBRAIC_,
        exportedFieldNames(),
        BuenoOrovioSTATES_NAMES,NUM_STATES,
        BuenoOrovioALGEBRAIC_NAMES,NUM_ALGEBRAIC,
        outFields
    );
}


void Foam::BuenoOrovio::writeHeader(OFstream& os) const
{
    ionicModelIO::writeHeader(
        os,
        BuenoOrovioSTATES_NAMES,NUM_STATES,
        BuenoOrovioALGEBRAIC_NAMES,NUM_ALGEBRAIC
    );
}

static Foam::scalar BO_vm(const Foam::scalarField& S)
{
    return S[0] * 85.7 - 84.0;
}
void Foam::BuenoOrovio::write(const scalar t, OFstream& os) const
{
    ionicModelIO::write(
        t,os,
        STATES_,ALGEBRAIC_,
        RATES_,
        BO_vm
    );
}


// ------------------------------------------------------------------------- //
//  State sync (used by manufactured / restart logic)
// ------------------------------------------------------------------------- //

void Foam::BuenoOrovio::updateStatesOld
(
    const Field<Field<scalar>>& states
) const
{
    forAll(states, i)
    {
        STATES_OLD_[i] = STATES_[i];
    }
}
void Foam::BuenoOrovio::resetStatesToStatesOld
(
    Field<Field<scalar>>& states
) const
{
    forAll(states, i)
    {
        STATES_[i] = STATES_OLD_[i];
    }
}






