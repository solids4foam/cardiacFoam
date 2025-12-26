/*---------------------------------------------------------------------------*\
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

\*---------------------------------------------------------------------------*/

#include <math.h>
#include "ORd.H"
#include "ORd_2011.H"
#include "addToRunTimeSelectionTable.H"
#include "ionicModel.H"
#include "ionicModelIO.H"
#include "stimulusIO.H"
#include "volFields.H"
#include "HashTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(ORd, 0);
    addToRunTimeSelectionTable
    (
        ionicModel, ORd, dictionary
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ORd::ORd
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
    Info<< nl << "Calling ORd initConsts" << endl;
    forAll(STATES_, i)
    {
        STATES_.set(i,      new scalarField(NUM_STATES,     0.0));
        STATES_OLD_.set(i,  new scalarField(NUM_STATES,     0.0));
        ALGEBRAIC_.set(i,   new scalarField(NUM_ALGEBRAIC,  0.0));
        RATES_.set(i,       new scalarField(NUM_STATES,     0.0));

        // Initialise constants, states and rates from generated code
        ORdinitConsts
        (
            CONSTANTS_.data(),
            RATES_[i].data(),
            STATES_[i].data(),
            tissue(),dict
        );
        if (!utilitiesMode())
        {
            stimulusIO::loadStimulusProtocol
            (
                dict, CONSTANTS_, stim_start, stim_period_S1,stim_duration,
                stim_amplitude, nstim1, stim_period_S2, nstim2
            );
        }
    }
    Info<< CONSTANTS_ << nl;

    label i0 = rand() % STATES_.size();
    Info<< "initial states:" << nl;
    Info<< STATES_[i0] << nl;

}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ORd::~ORd()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::List<Foam::word> Foam::ORd::supportedTissueTypes() const
{
    return {"myocyte"};
}


void Foam::ORd::calculateCurrent
(
    const scalar stepStartTime,
    const scalar deltaT,
    const scalarField& Vm,
    scalarField& Im,
    Field<Field<scalar>>& states
)
{
    const scalar tStart = stepStartTime * 1000.0;
    if (Im.size() != Vm.size())
    {
        FatalErrorInFunction
            << "Im.size() != Vm.size()" << abort(FatalError);
    }
    forAll(STATES_, integrationPtI)
    {
        scalarField& STATESI    = STATES_[integrationPtI];
        scalarField& ALGEBRAICI = ALGEBRAIC_[integrationPtI];
        scalarField& RATESI     = RATES_[integrationPtI];

        // Update voltage for this integration point
        STATESI[membrane_V] = Vm[integrationPtI] * 1000;

        ::ORdcomputeVariables
        (
            tStart,
            CONSTANTS_.data(),
            RATESI.data(),
            STATESI.data(),
            ALGEBRAICI.data(),
            tissue(),
            solveVmWithinODESolver()
        );
        // Jion  is the total ionic current density used by the PDE
        Im[integrationPtI] = ALGEBRAICI[Iion_cm];

        //copy internal STATES to memory external state buffer.
        //----can easily be expanded for all variables------//
        //copyInternalToExternal(STATES_, states, NUM_STATES);
    }
}


// ------------------------------------------------------------------------- //
//  Solve ODE with mixed singleCell implementation and 1D-3D condition
// ------------------------------------------------------------------------- //
void Foam::ORd::solveODE
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

        // If Vm is solved by the PDE, feed that Vm (in mV) into the cell model
        if (!solveVmWithinODESolver())
        {
            STATESI[membrane_V] = Vm[integrationPtI]*1000.0;
        }

        // Clamp time step (ms)
        step = min(step, deltaT * 1000.0);
        // Advance ODE system for all states
        odeSolver().solve(tStart, tEnd, STATESI, step);

        // Update algebraics and rates at tEnd (includes Iion and I_stim)
        ::ORdcomputeVariables
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
        {debugPrintFields(integrationPtI, tStart, tEnd, step);}

        // Total ionic current density used by PDE
        Im[integrationPtI] = ALGEBRAICI[Iion_cm] ;

        //----can easily be expanded for all variables------//
        // copyInternalToExternal(STATES_, states, NUM_STATES);
    }
}

void Foam::ORd::derivatives
(
    const scalar t,
    const scalarField& y,
    scalarField& dydt
) const
{
    // Must match NUM_ALGEBRAIC from the generated ORd code
    scalarField ALGEBRAIC_TMP(NUM_ALGEBRAIC, 0.0);

    ::ORdcomputeVariables
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

void Foam::ORd::updateStatesOld(const Field<Field<scalar>>&) const
{
    saveStateSnapshot(STATES_, STATES_OLD_);
}

void Foam::ORd::resetStatesToStatesOld(Field<Field<scalar>>&) const
{
    restoreStateSnapshot(STATES_, STATES_OLD_);
}

// ------------------------------------------------------------------------- //
//  Writing logic in singleCell and 3D simulations

//Writing functions for singleCell implementation
Foam::wordList Foam::ORd::exportedFieldNames() const
    {
        return ionicModelIO::exportedFieldNames
        (
            variableExport_,
            ORdSTATES_NAMES, NUM_STATES,
            ORdALGEBRAIC_NAMES, NUM_ALGEBRAIC
        );
    }

    Foam::wordList Foam::ORd::debugPrintedNames() const
    {
        return ionicModelIO::exportedFieldNames
        (
            debugVarNames_,
            ORdSTATES_NAMES, NUM_STATES,
            ORdALGEBRAIC_NAMES, NUM_ALGEBRAIC
        );
    }

void Foam::ORd::exportStates
(
    const Field<Field<scalar>>&,
    PtrList<volScalarField>& outFields
)
{
    ionicModelIO::exportStateFields
    (
        STATES_,ALGEBRAIC_,
        exportedFieldNames(),
        ORdSTATES_NAMES,NUM_STATES,
        ORdALGEBRAIC_NAMES,NUM_ALGEBRAIC,
        outFields
    );
}

void Foam::ORd::debugPrintFields
(
    label cellI,
    scalar t1,
    scalar t2,
    scalar step
) const
{
    ionicModelIO::debugPrintFields
    (
        STATES_, ALGEBRAIC_,
        debugPrintedNames(),
        ORdSTATES_NAMES, NUM_STATES,
        ORdALGEBRAIC_NAMES, NUM_ALGEBRAIC,
        cellI,t1,t2,step
    );
}



void Foam::ORd::writeHeader(OFstream& os) const
{
    const wordList names = exportedFieldNames();

    if (!names.empty())
    {
        ionicModelIO::writeSelectedHeader(os, names);
    }
    else
    {
        ionicModelIO::writeHeader
        (
            os,
            ORdSTATES_NAMES, NUM_STATES,
            ORdALGEBRAIC_NAMES, NUM_ALGEBRAIC
        );
    }
}

static Foam::scalar ORd_vm(const Foam::scalarField& S)
{
    return S[0];
}

void Foam::ORd::write(const scalar t, OFstream& os) const
{
    const wordList names = exportedFieldNames();

    if (!names.empty())
    {
        ionicModelIO::writeSelected
        (
            t, os,
            STATES_, ALGEBRAIC_,
            names,
            ORdSTATES_NAMES, NUM_STATES,
            ORdALGEBRAIC_NAMES, NUM_ALGEBRAIC
        );
    }
    else
    {
        ionicModelIO::write
        (
            t,
            os,
            STATES_, ALGEBRAIC_, RATES_,
            ORd_vm
        );
    }
}
void Foam::ORd::sweepCurrent
(
    const word& currentName,
    scalar Vmin,
    scalar Vmax,
    label nPts,
    const fileName& outputFile
) const
{
    // Retrieve dependency variables
    const auto& depMap = ORdDependencyMap();

    if (!depMap.found(currentName))
    {
        FatalErrorInFunction
            << "Unknown current: " << currentName << nl
            << "Available currents: " << depMap.toc() << nl
            << exit(FatalError);
    }

    const wordList& deps = depMap[currentName];
    OFstream os(outputFile);
    // Write sweep header: V,<deps...>
    ionicModelIO::writeSweepHeader(os, deps);

    // Working arrays from integration point 0
    scalarField STATESI = STATES_[0];
    scalarField RATESI(NUM_STATES, 0.0);
    scalarField ALGI(NUM_ALGEBRAIC, 0.0);

    // Voltage sweep
    for (label i = 0; i < nPts; ++i)
    {
        scalar V = Vmin + (Vmax - Vmin) * scalar(i) / (nPts - 1);

        // Reset all states to baseline
        STATESI = STATES_[0];

        // Overwrite membrane voltage (dimensionless in BO2008)
        STATESI[0] = V;

        ::ORdcomputeVariables
        (
            0.0,                       // VOI
            CONSTANTS_.data(),
            RATESI.data(),
            STATESI.data(),
            ALGI.data(),
            tissue(),
            solveVmWithinODESolver()
        );

        ionicModelIO::writeOneSweepRow
        (
            os, V, deps,STATESI,ALGI,
            ORdSTATES_NAMES, NUM_STATES,
            ORdALGEBRAIC_NAMES, NUM_ALGEBRAIC
        );
    }
    Info<< "Sweep for " << currentName
        << " written to " << outputFile << nl;
}
