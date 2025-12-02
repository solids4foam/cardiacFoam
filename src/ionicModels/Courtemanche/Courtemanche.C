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
#include "Courtemanche_1998.H"
#include "Courtemanche.H"
#include "addToRunTimeSelectionTable.H"
#include "ionicModel.H"
#include "ionicModelIO.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(Courtemanche, 0);
    addToRunTimeSelectionTable
    (
        ionicModel, Courtemanche, dictionary
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Courtemanche::Courtemanche
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
    // 🔑 First, set tissue using base logic + overrides
    ionicModel::setTissueFromDict();
    Info<< nl << "Calling Courtemanche initConsts" << endl;
    forAll(STATES_, i)
    {
        STATES_.set(i,      new scalarField(NUM_STATES,     0.0));
        STATES_OLD_.set(i,  new scalarField(NUM_STATES,     0.0));
        ALGEBRAIC_.set(i,   new scalarField(NUM_ALGEBRAIC,  0.0));
        RATES_.set(i,       new scalarField(NUM_STATES,     0.0));

        // Initialise constants, states and rates from generated code
        CourtemancheinitConsts
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

Foam::Courtemanche::~Courtemanche()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::List<Foam::word> Foam::Courtemanche::supportedTissueTypes() const

{
    return {"myocyte"};
}

void Foam::Courtemanche::calculateCurrent
(
    const scalar stepStartTime,
    const scalar deltaT,
    const scalarField& Vm,
    scalarField& Im,
    Field<Field<scalar>>& states
)
{
    // stepStartTime, deltaT are in seconds; Courtemanche uses ms
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
        STATESI[0] = Vm[integrationPtI] * 1000.0;

        if (integrationPtI == monitorCell)
        {
            Info<< "calculateCurrent: Vm[mV] before computeVariables = "
                << STATESI[0] << nl;
        }

        ::CourtemanchecomputeVariables
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
//  Solve ODE with mixed singleCell implementation and 1D-3D condition
// ------------------------------------------------------------------------- //
void Foam::Courtemanche::solveODE
(
    const scalar stepStartTime,
    const scalar deltaT,
    const scalarField& Vm,
    scalarField& Im,
    Field<Field<scalar>>& states
)
{
    const scalar tStart = stepStartTime * 1000;
    const scalar tEnd   = (stepStartTime + deltaT) * 1000;
    const label monitorCell = 0;

    forAll(STATES_, integrationPtI)
    {
        scalarField& STATESI    = STATES_[integrationPtI];
        scalarField& ALGEBRAICI = ALGEBRAIC_[integrationPtI];
        scalarField& RATESI     = RATES_[integrationPtI];

        scalar& step = ionicModel::step()[integrationPtI];

        if (!solveVmWithinODESolver())
        {
            //Vm comming from the PDE
            STATESI[membrane_V] = Vm[integrationPtI]*1000;
        }
        
        // Clamp time step
        step = min(step, deltaT * 1000);

        
        // Advance ODE system
        odeSolver().solve(tStart, tEnd, STATESI, step);

        // Update algebraics after solve
        ::CourtemanchecomputeVariables
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
            Info<< "solveODE: integrationPtI=" << integrationPtI
                << " | t = " << tStart << " → " << tEnd
                << " | step = " << step
                << " | Vm = " << STATESI[membrane_V]
                << " | Iion = " << ALGEBRAICI[Iion_cm]
                << " | Iext = " << ALGEBRAICI[AV_I_stim]
                << " | dVdt = " << RATESI[membrane_V]
                << endl;
        }


        Im[integrationPtI] = ALGEBRAICI[Iion_cm];

        // Export some states back if requested
        if (states[integrationPtI].size() >= NUM_STATES)
        {
            for (label s = 0; s < NUM_STATES; ++s)
            {
                states[integrationPtI][s] = STATESI[s];
            }
        }
    }
}

void Foam::Courtemanche::derivatives
(
    const scalar t,
    const scalarField& y,
    scalarField& dydt
) const
{
    // Must match NUM_ALGEBRAIC from the generated code
    scalarField ALGEBRAIC_TMP(NUM_ALGEBRAIC, 0.0);

    ::CourtemanchecomputeVariables
    (
        t,
        CONSTANTS_.data(),
        dydt.data(),                              // RATES
        const_cast<scalarField&>(y).data(),       // STATES
        ALGEBRAIC_TMP.data(),                     // ALGEBRAIC
        tissue(),
        solveVmWithinODESolver()
    );

}


// ------------------------------------------------------------------------- //
//  Writing logic in singleCell and 3D simulations
// ------------------------------------------------------------------------- //
Foam::wordList Foam::Courtemanche::exportedFieldNames() const
    {
        return ionicModelIO::exportedFieldNames
        (
            variableExport_,
            CourtemancheSTATES_NAMES, NUM_STATES,
            CourtemancheALGEBRAIC_NAMES, NUM_ALGEBRAIC
        );
    } 

void Foam::Courtemanche::exportStates
(
    const Field<Field<scalar>>&,
    List<volScalarField*>& outFields
)
{
    ionicModelIO::exportStateFields
    (
        STATES_,ALGEBRAIC_,
        exportedFieldNames(),
        CourtemancheSTATES_NAMES,NUM_STATES,
        CourtemancheALGEBRAIC_NAMES,NUM_ALGEBRAIC,
        outFields
    );
}

void Foam::Courtemanche::writeHeader(OFstream& os) const
{
    ionicModelIO::writeHeader(
        os,
        CourtemancheSTATES_NAMES, NUM_STATES,
        CourtemancheALGEBRAIC_NAMES, NUM_ALGEBRAIC
    );
}


static Foam::scalar CO_VM(const Foam::scalarField& S)
{
    return S[0];
}

void Foam::Courtemanche::write(const scalar t, OFstream& os) const
{
    ionicModelIO::write(
        t,
        os,
        STATES_,
        ALGEBRAIC_,
        RATES_,
        CO_VM
    );
}


// ------------------------------------------------------------------------- //
//  State sync
// ------------------------------------------------------------------------- //

void Foam::Courtemanche::updateStatesOld
(
    const Field<Field<scalar>>& states
) const
{
    forAll(states, i)
        STATES_OLD_[i] = STATES_[i];
}

void Foam::Courtemanche::resetStatesToStatesOld
(
    Field<Field<scalar>>& states
) const
{
    forAll(states, i)
        STATES_[i] = STATES_OLD_[i];
}



