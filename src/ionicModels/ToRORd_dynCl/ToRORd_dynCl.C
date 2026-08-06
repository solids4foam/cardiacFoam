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

#include "ToRORd_dynCl.H"
#include "ToRORd_dynCl_2023.H"
#include "HashTable.H"
#include "addToRunTimeSelectionTable.H"
#include "ionicModel.H"
#include "ionicHeterogeneityOrchestrator.H"
#include "ionicModelIO.H"
#include "stimulusIO.H"
#include "volFields.H"

#include <math.h>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(ToRORd_dynCl, 0);
    addToRunTimeSelectionTable
    (
        ionicModel, ToRORd_dynCl, dictionary
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ToRORd_dynCl::ToRORd_dynCl
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

    // First, set tissue using base logic and overrides
    ionicModel::setTissueFromDict();
    forAll(STATES_, i)
    {
        STATES_.set(i,      new scalarField(NUM_STATES,     0.0));
        ALGEBRAIC_.set(i,   new scalarField(NUM_ALGEBRAIC,  0.0));
        RATES_.set(i,       new scalarField(NUM_STATES,     0.0));

        ToRORd_dynClinitConsts
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

Foam::ToRORd_dynCl::~ToRORd_dynCl()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::List<Foam::word> Foam::ToRORd_dynCl::supportedTissueTypes() const
{
    return {"endocardialCells", "mCells", "epicardialCells"};
}


// ------------------------------------------------------------------------- //
//  Solve the cell ODE over [tStart, tEnd], converting time bounds to ms for the model
// ------------------------------------------------------------------------- //
Foam::scalarField& Foam::ToRORd_dynCl::constants
(
    const label integrationPtI
) const
{
    if (!HETEROGENEOUS_CONSTANTS_.empty())
    {
        return HETEROGENEOUS_CONSTANTS_[integrationPtI];
    }
    return CONSTANTS_;
}

void Foam::ToRORd_dynCl::configureIonicHeterogeneity
(
    const scalarField& transmuralDistance,
    const dictionary& heterogeneityDict
)
{
    ionicHeterogeneityOrchestrator::configureTransmuralBandHeterogeneity
    (
        *this, transmuralDistance, heterogeneityDict, HETEROGENEOUS_CONSTANTS_,
        &HETEROGENEOUS_INITIAL_STATES_
    );

    if (!HETEROGENEOUS_INITIAL_STATES_.empty())
    {
        for (label cellI = 0; cellI < STATES_.size(); ++cellI)
        {
            const scalarField& st = HETEROGENEOUS_INITIAL_STATES_[cellI];
            for (label stateI = 0; stateI < NUM_STATES; ++stateI)
            {
                STATES_[cellI][stateI] = st[stateI];
            }
        }
    }
}


Foam::scalarField Foam::ToRORd_dynCl::constantsForTissue
(
    const label tissueFlag
) const
{
    scalarField constants(NUM_CONSTANTS, 0.0);
    scalarField rates(NUM_STATES, 0.0);
    scalarField states(NUM_STATES, 0.0);

    ToRORd_dynClinitConsts
    (
        constants.data(), rates.data(), states.data(), tissueFlag, dict()
    );

    ionicModelIO::applyConstantOverrides
    (
        constants,
        ToRORd_dynClCONSTANTS_NAMES,
        NUM_CONSTANTS,
        dict(),
        type(),
        tissueFlag
    );

    return constants;
}


Foam::scalarField Foam::ToRORd_dynCl::initialStatesForTissue
(
    const label tissueFlag
) const
{
    scalarField constants(NUM_CONSTANTS, 0.0);
    scalarField rates(NUM_STATES, 0.0);
    scalarField states(NUM_STATES, 0.0);

    ToRORd_dynClinitConsts
    (
        constants.data(),
        rates.data(),
        states.data(),
        tissueFlag,
        dict()
    );

    return states;
}




void Foam::ToRORd_dynCl::solveODE
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
            STATESI[V] = Vm[integrationPtI]*1000.0;
        }

        step = min(step, deltaT * 1000.0);
        activeIntegrationPoint_ = integrationPtI;
        setActiveVmRate(integrationPtI);
        odeSolver().solve(tStart, tEnd, STATESI, step);

        {
            const scalar sum =
                STATESI[IKr_C1] + STATESI[IKr_C2] + STATESI[IKr_C3]
              + STATESI[IKr_I]  + STATESI[IKr_O];
            const scalar inv = 1.0/sum;
            STATESI[IKr_C1] *= inv;
            STATESI[IKr_C2] *= inv;
            STATESI[IKr_C3] *= inv;
            STATESI[IKr_I]  *= inv;
            STATESI[IKr_O]  *= inv;
        }

        ::ToRORd_dynClcomputeVariables
        (
            tEnd,
            constants(integrationPtI).data(),
            RATESI.data(),
            STATESI.data(),
            ALGEBRAICI.data(),
            solveVmWithinODESolver()
        ,
            stimulusProtocol()
        );
        if (integrationPtI == sampleCell)
        {debugPrintFields(integrationPtI, tStart, tEnd, step);}

        Im[integrationPtI] = ALGEBRAICI[Iion_cm] ;
    }

    clearVmRate();
}


void Foam::ToRORd_dynCl::evaluateIonicCurrent
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
        S[V] = Vm[integrationPtI]*1000.0;
        A = 0.0;
        R = 0.0;

        ::ToRORd_dynClcomputeVariables
        (
            t,
            constants(integrationPtI).data(),
            R.data(),
            S.data(),
            A.data(),
            solveVmWithinODESolver()
        ,
            stimulusProtocol()
        );

        Im[integrationPtI] = A[Iion_cm];
    }
}

void Foam::ToRORd_dynCl::derivatives
(
    const scalar t,
    const scalarField& y,
    scalarField& dydt
) const
{
    scalarField ALGEBRAIC_TMP(NUM_ALGEBRAIC, 0.0);

    ::ToRORd_dynClcomputeVariables
    (
        t,
        constants(activeIntegrationPoint_).data(),
        dydt.data(),                              // RATES (output)
        const_cast<scalarField&>(y).data(),       // STATES (input)
        ALGEBRAIC_TMP.data(),                     // ALGEBRAIC (scratch)
        solveVmWithinODESolver()
    ,
            stimulusProtocol()
        );

    if (!solveVmWithinODESolver())
    {
        dydt[V] = activeVmRate();
    }
}

// ------------------------------------------------------------------------- //
//  Writing logic in singleCell and 3D simulations

const char* const* Foam::ToRORd_dynCl::ioStateNames() const
{
    return ToRORd_dynClSTATES_NAMES;
}

const char* const* Foam::ToRORd_dynCl::ioConstantNames() const
{
    return ToRORd_dynClCONSTANTS_NAMES;
}

const char* const* Foam::ToRORd_dynCl::ioAlgebraicNames() const
{
    return ToRORd_dynClALGEBRAIC_NAMES;
}

void Foam::ToRORd_dynCl::sweepCurrent
(
    const word& currentName,
    scalar Vmin,
    scalar Vmax,
    label nPts,
    const fileName& outputFile
) const
{
    const auto& depMap = ToRORd_dynClDependencyMap();

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

        STATESI[0] = V;

        ::ToRORd_dynClcomputeVariables
        (
            0.0,                       // VOI
            CONSTANTS_.data(),
            RATESI.data(),
            STATESI.data(),
            ALGI.data(),
            solveVmWithinODESolver()
        ,
            stimulusProtocol()
        );

        ionicModelIO::writeOneSweepRow
        (
            os, V, deps,STATESI,ALGI,
            ToRORd_dynClSTATES_NAMES, NUM_STATES,
            ToRORd_dynClALGEBRAIC_NAMES, NUM_ALGEBRAIC,
            RATESI,
            sweepPlanCache
        );
    }
}

Foam::wordList Foam::ToRORd_dynCl::availableSweepCurrents() const
{
    return ToRORd_dynClDependencyMap().toc();
}
