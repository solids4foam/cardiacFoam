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
#include "Gaur.H"
#include "Gaur_2021.H"
#include "addToRunTimeSelectionTable.H"
#include "ionicModel.H"
#include "ionicModelIO.H"
//Only needs strings for the header writing
//#include <string>
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(Gaur, 0);
    addToRunTimeSelectionTable
    (
        ionicModel, Gaur, dictionary
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Gaur::Gaur
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
    Info<< nl << "Calling Gaur initConsts" << endl;
    forAll(STATES_, i)
    {
        STATES_.set(i,      new scalarField(NUM_STATES,     0.0));
        STATES_OLD_.set(i,  new scalarField(NUM_STATES,     0.0));
        ALGEBRAIC_.set(i,   new scalarField(NUM_ALGEBRAIC,  0.0));
        RATES_.set(i,       new scalarField(NUM_STATES,     0.0));

        // Initialise constants, states and rates from generated code
        GaurinitConsts
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

Foam::Gaur::~Gaur()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::List<Foam::word> Foam::Gaur::supportedTissueTypes() const
{
    return {"myocyte"};
}


void Foam::Gaur::calculateCurrent
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
        STATESI[cell_v] = Vm[integrationPtI] * 1000;

        ::GaurcomputeVariables
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
void Foam::Gaur::solveODE
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
            STATESI[cell_v] = Vm[integrationPtI]*1000.0;
        }

        // Clamp time step (ms)
        step = min(step, deltaT * 1000.0);
        // Advance ODE system for all states
        odeSolver().solve(tStart, tEnd, STATESI, step);

        // Update algebraics and rates at tEnd (includes Iion and I_stim)
        ::GaurcomputeVariables
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

void Foam::Gaur::derivatives
(
    const scalar t,
    const scalarField& y,
    scalarField& dydt
) const
{
    // Must match NUM_ALGEBRAIC from the generated Gaur code
    scalarField ALGEBRAIC_TMP(NUM_ALGEBRAIC, 0.0);

    ::GaurcomputeVariables
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

void Foam::Gaur::updateStatesOld(const Field<Field<scalar>>&) const
{
    saveStateSnapshot(STATES_, STATES_OLD_);
}

void Foam::Gaur::resetStatesToStatesOld(Field<Field<scalar>>&) const
{
    restoreStateSnapshot(STATES_, STATES_OLD_);
}

// ------------------------------------------------------------------------- //
//  Writing logic in singleCell and 3D simulations

//Writing functions for singleCell implementation
Foam::wordList Foam::Gaur::exportedFieldNames() const
    {
        return ionicModelIO::exportedFieldNames
        (
            variableExport_,
            GaurSTATES_NAMES, NUM_STATES,
            GaurALGEBRAIC_NAMES, NUM_ALGEBRAIC
        );
    }

    Foam::wordList Foam::Gaur::debugPrintedNames() const
    {
        return ionicModelIO::exportedFieldNames
        (
            debugVarNames_,
            GaurSTATES_NAMES, NUM_STATES,
            GaurALGEBRAIC_NAMES, NUM_ALGEBRAIC
        );
    }

void Foam::Gaur::exportStates
(
    const Field<Field<scalar>>&,
    PtrList<volScalarField>& outFields
)
{
    ionicModelIO::exportStateFields
    (
        STATES_,ALGEBRAIC_,
        exportedFieldNames(),
        GaurSTATES_NAMES,NUM_STATES,
        GaurALGEBRAIC_NAMES,NUM_ALGEBRAIC,
        outFields
    );
}

void Foam::Gaur::debugPrintFields
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
        GaurSTATES_NAMES, NUM_STATES,
        GaurALGEBRAIC_NAMES, NUM_ALGEBRAIC,
        cellI,t1,t2,step
    );
}



//Writing functions for singleCell implementation
void Foam::Gaur::writeHeader(OFstream& os) const
{
    ionicModelIO::writeHeader(
        os,
        GaurSTATES_NAMES,NUM_STATES,
        GaurALGEBRAIC_NAMES,NUM_ALGEBRAIC
    );
}

static Foam::scalar Gaur_vm(const Foam::scalarField& S)
{
    return S[0];
}
void Foam::Gaur::write(const scalar t, OFstream& os) const
{
    ionicModelIO::write(
        t,os,
        STATES_,ALGEBRAIC_,
        RATES_,
        Gaur_vm
    );
}



