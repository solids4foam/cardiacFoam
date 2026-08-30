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

#include "TWorld.H"
#include "TWorld_2025.H"
#include "addToRunTimeSelectionTable.H"
#include "ionicModel.H"
#include "ionicModelIO.H"
#include "stimulusIO.H"
#include "volFields.H"
#include "HashTable.H"

#include <math.h>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(TWorld, 0);
    addToRunTimeSelectionTable
    (
        ionicModel, TWorld, dictionary
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::TWorld::TWorld
(
    const dictionary& dict,
    const label num,
    const scalar initialDeltaT,
    const Switch solveVmWithinODESolver
)
:
    configuredIonicModel(dict, num, initialDeltaT, solveVmWithinODESolver),
    STATES_(num),
    CONSTANTS_(NUM_CONSTANTS, 0.0),
    ALGEBRAIC_(num),
    RATES_(num)
{

    // First, set tissue using base logic and model overrides
    ionicModel::setTissueFromDict();
    ionicModel::setSexFromDict();
    forAll(STATES_, i)
    {
        STATES_.set(i,      new scalarField(NUM_STATES,     0.0));
        ALGEBRAIC_.set(i,   new scalarField(NUM_ALGEBRAIC,  0.0));
        RATES_.set(i,       new scalarField(NUM_STATES,     0.0));

        TWorldinitConsts
        (
            CONSTANTS_.data(),
            RATES_[i].data(),
            STATES_[i].data(),
            tissue(),
            sex(),
            dict
        );

        if (!utilitiesMode())
        {
            setStimulusProtocolFromDict(dict);
        }
    }

    applyIonicConstantOverrides();

}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::TWorld::~TWorld()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::List<Foam::word> Foam::TWorld::supportedTissueTypes() const
{
    return {"epicardialCells", "mCells", "endocardialCells"};
}

Foam::List<Foam::word> Foam::TWorld::supportedSexTypes() const
{
    return {"neutral", "male", "female"};
}


// ------------------------------------------------------------------------- //
Foam::scalarField& Foam::TWorld::constants(const label integrationPtI) const
{
    if (!HETEROGENEOUS_CONSTANTS_.empty())
    {
        return HETEROGENEOUS_CONSTANTS_[integrationPtI];
    }
    return CONSTANTS_;
}


Foam::scalarField Foam::TWorld::constantsForTissue
(
    const label tissueFlag
) const
{
    scalarField constants(NUM_CONSTANTS, 0.0);
    scalarField rates(NUM_STATES, 0.0);
    scalarField states(NUM_STATES, 0.0);

    TWorldinitConsts
    (
        constants.data(), rates.data(), states.data(), tissueFlag, sex(), dict()
    );

    ionicModelIO::applyConstantOverrides
    (
        constants,
        TWorldCONSTANTS_NAMES,
        NUM_CONSTANTS,
        dict(),
        type(),
        tissueFlag
    );

    return constants;
}


//  Solve the cell ODE over [tStart, tEnd], converting time bounds to ms for the model
// ------------------------------------------------------------------------- //
void Foam::TWorld::solveODE
(
    const scalar stepStartTime,
    const scalar deltaT,
    const scalarField& Vm,
    scalarField& Im
)
{
    const scalar tStart = stepStartTime * 1000;
    const scalar tEnd   = (stepStartTime + deltaT) * 1000;
    const label sampleCell = sampleIntegrationPoint(STATES_.size());

    forAll(STATES_, integrationPtI)
    {
        scalarField& STATESI    = STATES_[integrationPtI];
        scalarField& ALGEBRAICI = ALGEBRAIC_[integrationPtI];
        scalarField& RATESI     = RATES_[integrationPtI];

        scalar& step = ionicModel::step()[integrationPtI];

        if (!solveVmWithinODESolver())
        {
            STATESI[v] = Vm[integrationPtI]*1000.0;
        }

        step = min(step, deltaT * 1000.0);
        activeIntegrationPoint_ = integrationPtI;
        setActiveVmRate(integrationPtI);
        odeSolver().solve(tStart, tEnd, STATESI, step);

        ::TWorldcomputeVariables
        (
            tEnd,
            constants(integrationPtI).data(),
            RATESI.data(),
            STATESI.data(),
            ALGEBRAICI.data(),
            solveVmWithinODESolver(),
            stimulusProtocol()
        );
        if (integrationPtI == sampleCell)
        {debugPrintFields(integrationPtI, tStart, tEnd, step);}

        Im[integrationPtI] = ALGEBRAICI[Iion_cm] ;
    }

    clearVmRate();
}


void Foam::TWorld::evaluateIonicCurrent
(
    const scalar t,
    const scalarField& Vm,
    scalarField& Im
)
{
    scalarField S(NUM_STATES, 0.0);
    scalarField A(NUM_ALGEBRAIC, 0.0);
    scalarField R(NUM_STATES, 0.0);

    forAll(STATES_, integrationPtI)
    {
        S = STATES_[integrationPtI];
        S[v] = Vm[integrationPtI]*1000.0;
        A = 0.0;
        R = 0.0;

        ::TWorldcomputeVariables
        (
            t,
            constants(integrationPtI).data(),
            R.data(),
            S.data(),
            A.data(),
            solveVmWithinODESolver(),
            stimulusProtocol()
        );

        Im[integrationPtI] = A[Iion_cm];
    }
}

void Foam::TWorld::derivatives
(
    const scalar t,
    const scalarField& y,
    scalarField& dydt
) const
{
    scalarField ALGEBRAIC_TMP(NUM_ALGEBRAIC, 0.0);

    ::TWorldcomputeVariables
    (
        t,
        constants(activeIntegrationPoint_).data(),
        dydt.data(),                              // RATES (output)
        const_cast<scalarField&>(y).data(),       // STATES (input)
        ALGEBRAIC_TMP.data(),                     // ALGEBRAIC (scratch)
        solveVmWithinODESolver(),
        stimulusProtocol()
    );

    if (!solveVmWithinODESolver())
    {
        dydt[v] = activeVmRate();
    }
}

// ------------------------------------------------------------------------- //
//  Writing logic in singleCell and 3D simulations

const char* const* Foam::TWorld::ioStateNames() const
{
    return TWorldSTATES_NAMES;
}

const char* const* Foam::TWorld::ioConstantNames() const
{
    return TWorldCONSTANTS_NAMES;
}

const char* const* Foam::TWorld::ioAlgebraicNames() const
{
    return TWorldALGEBRAIC_NAMES;
}

void Foam::TWorld::sweepCurrent
(
    const word& currentName,
    scalar Vmin,
    scalar Vmax,
    label nPts,
    const fileName& outputFile
) const
{
    const auto& depMap = TWorldDependencyMap();

    if (!depMap.found(currentName))
    {
        FatalErrorInFunction
            << "Unknown current: " << currentName << nl
            << "Available currents: " << depMap.toc() << nl
            << exit(FatalError);
    }

    const wordList& deps = depMap[currentName];
    OFstream os(outputFile);
    ionicModelIO::writeSweepHeader(os, deps);

    scalarField STATESI = STATES_[0];
    scalarField RATESI(NUM_STATES, 0.0);
    scalarField ALGI(NUM_ALGEBRAIC, 0.0);
    ionicModelIO::SelectedMapCache sweepPlanCache;

    for (label i = 0; i < nPts; ++i)
    {
        scalar V = Vmin + (Vmax - Vmin) * scalar(i) / (nPts - 1);

        STATESI = STATES_[0];

        STATESI[v] = V;

        ::TWorldcomputeVariables
        (
            0.0,                       // VOI
            CONSTANTS_.data(),
            RATESI.data(),
            STATESI.data(),
            ALGI.data(),
            solveVmWithinODESolver(),
            stimulusProtocol()
        );

        ionicModelIO::writeOneSweepRow
        (
            os, V, deps, STATESI, ALGI,
            TWorldSTATES_NAMES, NUM_STATES,
            TWorldALGEBRAIC_NAMES, NUM_ALGEBRAIC,
            RATESI,
            sweepPlanCache
        );
    }
}

Foam::wordList Foam::TWorld::availableSweepCurrents() const
{
    return TWorldDependencyMap().toc();
}
