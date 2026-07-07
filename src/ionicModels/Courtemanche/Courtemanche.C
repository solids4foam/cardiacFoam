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

#include "Courtemanche.H"
#include "Courtemanche_1998.H"
#include "HashTable.H"
#include "addToRunTimeSelectionTable.H"
#include "ionicModel.H"
#include "ionicModelFamilyInfo.H"
#include "ionicModelIO.H"
#include "stimulusIO.H"
#include "volFields.H"

#include <math.h>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(Courtemanche, 0);
    addToRunTimeSelectionTable
    (
        ionicModel, Courtemanche, dictionary
    );

    const ionicModelFamilyInfo& CourtemancheFamilyInfo()
    {
        static const ionicModelFamilyInfo info
        {
            NUM_CONSTANTS,
            NUM_STATES,
            NUM_ALGEBRAIC,
            CourtemancheCONSTANTS_NAMES,
            CourtemancheSTATES_NAMES,
            CourtemancheALGEBRAIC_NAMES,
            membrane_V,
            1000.0,
            1000.0,
            0.0,
            nullptr
        };
        return info;
    }
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

        CourtemancheinitConsts
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

Foam::Courtemanche::~Courtemanche()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::List<Foam::word> Foam::Courtemanche::supportedTissueTypes() const
{
    return {"epicardialCells", "mCells", "endocardialCells", "myocyte"};
}


Foam::scalarField& Foam::Courtemanche::constants(const label integrationPtI) const
{
    if (!HETEROGENEOUS_CONSTANTS_.empty())
    {
        return HETEROGENEOUS_CONSTANTS_[integrationPtI];
    }
    return CONSTANTS_;
}


Foam::scalarField Foam::Courtemanche::constantsForTissue
(
    const label tissueFlag
) const
{
    scalarField constants(NUM_CONSTANTS, 0.0);
    scalarField rates(NUM_STATES, 0.0);
    scalarField states(NUM_STATES, 0.0);

    CourtemancheinitConsts
    (
        constants.data(), rates.data(), states.data(), tissueFlag, dict()
    );

    ionicModelIO::applyConstantOverrides
    (
        constants, CourtemancheCONSTANTS_NAMES, NUM_CONSTANTS, dict(), type(),
        tissueFlag
    );

    return constants;
}


Foam::scalarField Foam::Courtemanche::initialStatesForTissue
(
    const label tissueFlag
) const
{
    scalarField constants(NUM_CONSTANTS, 0.0);
    scalarField rates(NUM_STATES, 0.0);
    scalarField states(NUM_STATES, 0.0);

    CourtemancheinitConsts
    (
        constants.data(), rates.data(), states.data(), tissueFlag, dict()
    );

    return states;
}


void Foam::Courtemanche::configureIonicHeterogeneity
(
    const scalarField& transmuralDistance,
    const dictionary& heterogeneityDict
)
{
    configureTransmuralBandHeterogeneity
    (
        transmuralDistance, heterogeneityDict, HETEROGENEOUS_CONSTANTS_
    );
}


void Foam::Courtemanche::configureApexBaseBandsHeterogeneity
(
    const scalarField& apexDist,
    const dictionary& dict
)
{
    configureApexBaseBandsHeterogeneityImpl
    (
        apexDist, dict, HETEROGENEOUS_CONSTANTS_
    );
}


// ------------------------------------------------------------------------- //
//  Solve the cell ODE over [tStart, tEnd], converting time bounds to ms for the model
// ------------------------------------------------------------------------- //
void Foam::Courtemanche::solveODE
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
            STATESI[0] = Vm[integrationPtI]*1000.0;
        }

        step = min(step, deltaT * 1000.0);
        activeIntegrationPoint_ = integrationPtI;
        odeSolver().solve(tStart, tEnd, STATESI, step);

        ::CourtemanchecomputeVariables
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

        Im[integrationPtI] = ALGEBRAICI[Iion_cm] ;

    }
}

void Foam::Courtemanche::derivatives
(
    const scalar t,
    const scalarField& y,
    scalarField& dydt
) const
{
    scalarField ALGEBRAIC_TMP(NUM_ALGEBRAIC, 0.0);

    ::CourtemanchecomputeVariables
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
}

// ------------------------------------------------------------------------- //
//  Writing logic in singleCell and 3D simulations

const char* const* Foam::Courtemanche::ioStateNames() const
{
    return CourtemancheSTATES_NAMES;
}

const char* const* Foam::Courtemanche::ioConstantNames() const
{
    return CourtemancheCONSTANTS_NAMES;
}

const char* const* Foam::Courtemanche::ioAlgebraicNames() const
{
    return CourtemancheALGEBRAIC_NAMES;
}

void Foam::Courtemanche::sweepCurrent
(
    const word& currentName,
    scalar Vmin,
    scalar Vmax,
    label nPts,
    const fileName& outputFile
) const
{
    const auto& depMap = CourtemancheDependencyMap();

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

        ::CourtemanchecomputeVariables
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
            CourtemancheSTATES_NAMES, NUM_STATES,
            CourtemancheALGEBRAIC_NAMES, NUM_ALGEBRAIC,
            RATESI,
            sweepPlanCache
        );
    }
}

Foam::wordList Foam::Courtemanche::availableSweepCurrents() const
{
    return CourtemancheDependencyMap().toc();
}
