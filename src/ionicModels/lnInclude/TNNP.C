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
#include "TNNP.H"
#include "TNNP_2004.H"
#include "addToRunTimeSelectionTable.H"
#include "ionicModel.H"
#include"ionicModelIO.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(TNNP, 0);
    addToRunTimeSelectionTable
    (
        ionicModel, TNNP, dictionary
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::TNNP::TNNP
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
    Info<< nl << "Initialize TNNP constants:" << nl;
    forAll(STATES_, i)
    {
        STATES_.set(i,      new scalarField(NUM_STATES,     0.0));
        STATES_OLD_.set(i,  new scalarField(NUM_STATES,     0.0));
        ALGEBRAIC_.set(i,   new scalarField(NUM_ALGEBRAIC,  0.0));
        RATES_.set(i,       new scalarField(NUM_STATES,     0.0));

        // Initialise constants, states and rates from generated code
        TNNPinitConsts
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

Foam::TNNP::~TNNP()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::List<Foam::word> Foam::TNNP::supportedTissueTypes() const
{
    // All three tissue variants are supported in the generated code
    return {"endocardialCells", "mCells", "epicardialCells"};
}


// ------------------------------------------------------------------------- //
//  Explicit split: calculateCurrent (Iion only, no state update)
// ------------------------------------------------------------------------- //

void Foam::TNNP::calculateCurrent
(
    const scalar stepStartTime,
    const scalar deltaT,
    const scalarField& Vm,
    scalarField& Im,
    Field<Field<scalar>>& states
)
{
    // stepStartTime, deltaT are in seconds; TNNP uses ms
    const scalar tStart = stepStartTime * 1000.0;
    const label monitorCell = min(500, STATES_.size()-1);


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
        STATESI[0] = Vm[integrationPtI] * 1000.0;

        if (integrationPtI == monitorCell)
        {
            Info<< "calculateCurrent: Vm[mV] before computeVariables = "
                << STATESI[0] << nl;
        }

        ::TNNPcomputeVariables
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
                << " | Vm = " << STATESI[0]
                << " | Iion_cm = " << ALGEBRAICI[Iion_cm]
                << nl;
        }

        // Iion_cm (index 69) is the total ionic current density used by the PDE
        Im[integrationPtI] = ALGEBRAICI[Iion_cm];

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

void Foam::TNNP::solveODE
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
    const label monitorCell = min(500, STATES_.size()-1);

    forAll(STATES_, integrationPtI)
    {
        scalarField& STATESI    = STATES_[integrationPtI];
        scalarField& ALGEBRAICI = ALGEBRAIC_[integrationPtI];
        scalarField& RATESI     = RATES_[integrationPtI];

        // Per-cell adaptive time step (in ms) for the ODE solver
        scalar& h = ionicModel::step()[integrationPtI];

        // Vm fed into the cell model in mV
        if (!solveVmWithinODESolver())
        {
            STATESI[0] = Vm[integrationPtI] * 1000.0;
        }
        // Clamp ODE step
        h = min(h, deltaT * 1000.0);

        if (integrationPtI == monitorCell)
        {
            Info<< "solveODE: i=" << integrationPtI
                << " | t = " << tStart << " → " << tEnd
                << " | step = " << h
                << " | Vm = " << STATESI[0]
                << nl;
        }

        // Advance the ODE system
        odeSolver().solve(tStart, tEnd, STATESI, h);

        // Update ALGEBRAIC (incl. Iion_cm) and RATES at tEnd
        ::TNNPcomputeVariables
        (
            tEnd,
            CONSTANTS_.data(),
            RATESI.data(),
            STATESI.data(),
            ALGEBRAICI.data(),
            tissue(),
            solveVmWithinODESolver()
        );

        ::TNNPcomputeRates
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
                << " | Vm = " << STATESI[0]
                << " | Iion_cm = " << ALGEBRAICI[Iion_cm]
                << " | dVdt = " << RATESI[0]
                << nl;
        }

        // Total ionic current density used by PDE
        Im[integrationPtI] = ALGEBRAICI[Iion_cm];

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

//--------------------------------------------------//
void Foam::TNNP::derivatives
(
    const scalar t,
    const scalarField& y,
    scalarField& dydt
) const
{
    // Must match NUM_ALGEBRAIC from the generated TNNP code
    scalarField ALGEBRAIC_TMP(NUM_ALGEBRAIC, 0.0);

    ::TNNPcomputeRates
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
//  IO
// ------------------------------------------------------------------------- //
Foam::wordList Foam::TNNP::exportedFieldNames() const
    {
        return ionicModelIO::exportedFieldNames
        (
            variableExport_,
            TNNP_STATES_NAMES, NUM_STATES,
            TNNP_ALGEBRAIC_NAMES, NUM_ALGEBRAIC
        );
    } 

void Foam::TNNP::exportStates
(
    const Field<Field<scalar>>&,
    List<volScalarField*>& outFields
)
{
    ionicModelIO::exportStateFields
    (
        STATES_,ALGEBRAIC_,
        exportedFieldNames(),
        TNNP_STATES_NAMES,NUM_STATES,
        TNNP_ALGEBRAIC_NAMES,NUM_ALGEBRAIC,
        outFields
    );
}


void Foam::TNNP::writeHeader(OFstream& os) const
{
    ionicModelIO::writeHeader(
        os,
        TNNP_STATES_NAMES, NUM_STATES,
        TNNP_ALGEBRAIC_NAMES, NUM_ALGEBRAIC
    );
}


static Foam::scalar TNNP_Vm(const Foam::scalarField& S)
{
    // Convert normalized u → physical Vm(mV)
    return S[0];
}

void Foam::TNNP::write(const scalar t, OFstream& os) const
{
    ionicModelIO::write(
        t,
        os,
        STATES_,
        ALGEBRAIC_,
        RATES_,
        TNNP_Vm
    );
}


// ------------------------------------------------------------------------- //
//  State sync (used by manufactured / restart logic)
// ------------------------------------------------------------------------- //

void Foam::TNNP::updateStatesOld
(
    const Field<Field<scalar>>& states
) const
{
    forAll(states, i)
    {
        STATES_OLD_[i] = STATES_[i];
    }
}

void Foam::TNNP::resetStatesToStatesOld
(
    Field<Field<scalar>>& states
) const
{
    forAll(states, i)
    {
        STATES_[i] = STATES_OLD_[i];
    }
}

// ************************************************************************* //

