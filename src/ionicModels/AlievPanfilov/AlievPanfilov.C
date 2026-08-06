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

#include "AlievPanfilov.H"
#include "AlievPanfilov_1996.H"
#include "HashTable.H"
#include "addToRunTimeSelectionTable.H"
#include "ionicModel.H"
#include "ionicModelFamilyInfo.H"
#include "ionicModelIO.H"
#include "stimulusIO.H"
#include "volFields.H"

#include <math.h>

namespace
{
    Foam::scalar alievPanfilovTransformedVm(const Foam::scalarField& S)
    {
        return S[u]*100.0 - 80.0;
    }
}

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(AlievPanfilov, 0);
    addToRunTimeSelectionTable
    (
        ionicModel, AlievPanfilov, dictionary
    );

    const ionicModelFamilyInfo& AlievPanfilovFamilyInfo()
    {
        static const ionicModelFamilyInfo info
        {
            NUM_CONSTANTS,
            NUM_STATES,
            NUM_ALGEBRAIC,
            AlievPanfilovCONSTANTS_NAMES,
            AlievPanfilovSTATES_NAMES,
            AlievPanfilovALGEBRAIC_NAMES,
            u,
            1000.0/12.9,
            10.0,
            0.8,
            &alievPanfilovTransformedVm
        };
        return info;
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::AlievPanfilov::AlievPanfilov
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
    ionicModel::setTissueFromDict();
    forAll(STATES_, i)
    {
        STATES_.set(i,      new scalarField(NUM_STATES,     0.0));
        ALGEBRAIC_.set(i,   new scalarField(NUM_ALGEBRAIC,  0.0));
        RATES_.set(i,       new scalarField(NUM_STATES,     0.0));

        AlievPanfilovinitConsts
        (
            CONSTANTS_.data(),
            RATES_[i].data(),
            STATES_[i].data(),
            tissue(),dict
        );

        if (!utilitiesMode())
        {
            setStimulusProtocolFromDict(dict);
        }
    }

    applyIonicConstantOverrides();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::AlievPanfilov::~AlievPanfilov()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::List<Foam::word> Foam::AlievPanfilov::supportedTissueTypes() const
{
    return {"epicardialCells", "mCells", "endocardialCells", "myocyte"};
}


Foam::scalarField& Foam::AlievPanfilov::constants(const label integrationPtI) const
{
    if (!HETEROGENEOUS_CONSTANTS_.empty())
    {
        return HETEROGENEOUS_CONSTANTS_[integrationPtI];
    }
    return CONSTANTS_;
}


Foam::scalarField Foam::AlievPanfilov::constantsForTissue
(
    const label tissueFlag
) const
{
    scalarField constants(NUM_CONSTANTS, 0.0);
    scalarField rates(NUM_STATES, 0.0);
    scalarField states(NUM_STATES, 0.0);

    AlievPanfilovinitConsts
    (
        constants.data(), rates.data(), states.data(), tissueFlag, dict()
    );

    ionicModelIO::applyConstantOverrides
    (
        constants, AlievPanfilovCONSTANTS_NAMES, NUM_CONSTANTS, dict(), type(),
        tissueFlag
    );

    return constants;
}


Foam::scalarField Foam::AlievPanfilov::initialStatesForTissue
(
    const label tissueFlag
) const
{
    scalarField constants(NUM_CONSTANTS, 0.0);
    scalarField rates(NUM_STATES, 0.0);
    scalarField states(NUM_STATES, 0.0);

    AlievPanfilovinitConsts
    (
        constants.data(), rates.data(), states.data(), tissueFlag, dict()
    );

    return states;
}


void Foam::AlievPanfilov::solveODE
(
    const scalar stepStartTime,
    const scalar deltaT,
    const scalarField& Vm,
    scalarField& Im
)
{
    const scalar tStart = stepStartTime * 1000.0 / 12.9;
    const scalar tEnd   = (stepStartTime + deltaT) * 1000.0/12.9;
    const label sampleCell = sampleIntegrationPoint(STATES_.size());

    forAll(STATES_, integrationPtI)
    {
        scalarField& STATESI    = STATES_[integrationPtI];
        scalarField& ALGEBRAICI = ALGEBRAIC_[integrationPtI];
        scalarField& RATESI     = RATES_[integrationPtI];


        if (!solveVmWithinODESolver())
        {
            STATESI[0] = (Vm[integrationPtI] * 1000.0 + 80)/100;
        }
        scalar& step = ionicModel::step()[integrationPtI];

        step = min(step, deltaT * 1000.0/12.9);
        activeIntegrationPoint_ = integrationPtI;
        setActiveVmRate(integrationPtI, 1000.0/100.0, 1000.0);
        if (integrationPtI == sampleCell)
            {debugPrintFields(integrationPtI, tStart, tEnd, step);}

        odeSolver().solve(tStart, tEnd, STATESI, step);

        ::AlievPanfilovcomputeVariables
        (
            tEnd,
            constants(integrationPtI).data(),
            RATESI.data(),
            STATESI.data(),
            ALGEBRAICI.data(),
            tissue(),
            solveVmWithinODESolver()
        ,
            stimulusProtocol()
        );

        if (integrationPtI == sampleCell)
            {debugPrintFields(integrationPtI, tStart, tEnd, step);}

        Im[integrationPtI] = ALGEBRAICI[Iion_cm] * 100;
    }

    clearVmRate();
}


void Foam::AlievPanfilov::evaluateIonicCurrent
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
        S[0] = (Vm[integrationPtI] * 1000.0 + 80)/100;
        A = 0.0;
        R = 0.0;

        ::AlievPanfilovcomputeVariables
        (
            t,
            constants(integrationPtI).data(),
            R.data(),
            S.data(),
            A.data(),
            tissue(),
            solveVmWithinODESolver()
        ,
            stimulusProtocol()
        );

        Im[integrationPtI] = A[Iion_cm] * 100;
    }
}


void Foam::AlievPanfilov::derivatives
(
    const scalar t,
    const scalarField& y,
    scalarField& dydt
) const
{
    scalarField ALGEBRAIC_TMP(NUM_ALGEBRAIC, 0.0);

    ::AlievPanfilovcomputeVariables
    (
        t,
        constants(activeIntegrationPoint_).data(),
        dydt.data(),                              // RATES (output)
        const_cast<scalarField&>(y).data(),       // STATES (input)
        ALGEBRAIC_TMP.data(),                     // ALGEBRAIC (scratch)
        tissue(),
        solveVmWithinODESolver()
    ,
            stimulusProtocol()
        );

    if (!solveVmWithinODESolver())
    {
        dydt[0] = activeVmRate();
    }
}

// ------------------------------------------------------------------------- //
//  Writing logic in singleCell and 3D simulations

const char* const* Foam::AlievPanfilov::ioStateNames() const
{
    return AlievPanfilovSTATES_NAMES;
}

const char* const* Foam::AlievPanfilov::ioConstantNames() const
{
    return AlievPanfilovCONSTANTS_NAMES;
}

const char* const* Foam::AlievPanfilov::ioAlgebraicNames() const
{
    return AlievPanfilovALGEBRAIC_NAMES;
}

void Foam::AlievPanfilov::sweepCurrent
(
    const word& currentName,
    scalar Vmin,
    scalar Vmax,
    label nPts,
    const fileName& outputFile
) const
{
    const auto& depMap = AlievPanfilovDependencyMap();

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

        STATESI[u] = V;

        ::AlievPanfilovcomputeVariables
        (
            0.0,                       // VOI
            CONSTANTS_.data(),
            RATESI.data(),
            STATESI.data(),
            ALGI.data(),
            tissue(),
            solveVmWithinODESolver()
        ,
            stimulusProtocol()
        );

        ionicModelIO::writeOneSweepRow
        (
            os, V, deps,STATESI,ALGI,
            AlievPanfilovSTATES_NAMES, NUM_STATES,
            AlievPanfilovALGEBRAIC_NAMES, NUM_ALGEBRAIC,
            RATESI,
            sweepPlanCache
        );
    }
}

Foam::wordList Foam::AlievPanfilov::availableSweepCurrents() const
{
    return AlievPanfilovDependencyMap().toc();
}


//----------------------Coupling Signals-------------------------------- //

bool Foam::AlievPanfilov::hasSignal(const CouplingSignal s) const
{
    return ionicModel::hasSignal(s);
}

Foam::scalar Foam::AlievPanfilov::signal
(
    const label i,
    const CouplingSignal s
) const
{
    return ionicModel::signal(i, s);
}
