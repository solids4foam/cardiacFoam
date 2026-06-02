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

#include "BuenoOrovio.H"
#include "BuenoOrovio_2008.H"
#include "HashTable.H"
#include "addToRunTimeSelectionTable.H"
#include "ionicHeterogeneity.H"
#include "ionicModel.H"
#include "ionicModelIO.H"
#include "stimulusIO.H"
#include "volFields.H"

#include <math.h>


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(BuenoOrovio, 0);
    addToRunTimeSelectionTable
    (
        ionicModel, BuenoOrovio, dictionary
    );
}


namespace
{

Foam::scalarField constantsForTissue
(
    const Foam::label tissueFlag,
    const Foam::dictionary& dict
)
{
    Foam::scalarField constants(NUM_CONSTANTS, 0.0);
    Foam::scalarField rates(NUM_STATES, 0.0);
    Foam::scalarField states(NUM_STATES, 0.0);

    BuenoOrovioinitConsts
    (
        constants.data(),
        rates.data(),
        states.data(),
        tissueFlag,
        dict
    );

    Foam::ionicModelIO::applyConstantOverrides
    (
        constants,
        BuenoOrovioCONSTANTS_NAMES,
        NUM_CONSTANTS,
        dict,
        "BuenoOrovio"
    );

    return constants;
}

} // End anonymous namespace

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
    CONSTANTS_(NUM_CONSTANTS, 0.0),
    HETEROGENEOUS_CONSTANTS_(),
    activeIntegrationPoint_(0),
    ALGEBRAIC_(num),
    RATES_(num)

{
    ionicModel::setTissueFromDict();
    forAll(STATES_, i)
    {
        STATES_.set(i,      new scalarField(NUM_STATES,     0.0));
        ALGEBRAIC_.set(i,   new scalarField(NUM_ALGEBRAIC,  0.0));
        RATES_.set(i,       new scalarField(NUM_STATES,     0.0));

        BuenoOrovioinitConsts
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

Foam::BuenoOrovio::~BuenoOrovio()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalarField& Foam::BuenoOrovio::constants
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


Foam::List<Foam::word> Foam::BuenoOrovio::supportedTissueTypes() const
{
    return {"endocardialCells", "mCells", "epicardialCells"};
}


void Foam::BuenoOrovio::configureIonicHeterogeneity
(
    const scalarField& transmuralDistance,
    const dictionary& heterogeneityDict
)
{
    const word mode =
        heterogeneityDict.lookupOrDefault<word>("mode", "transmuralBands");

    if (mode != "transmuralBands")
    {
        FatalErrorInFunction
            << "Unsupported BuenoOrovio ionicHeterogeneity mode '" << mode
            << "'. Supported mode: transmuralBands."
            << exit(FatalError);
    }

    const word smoothing =
        heterogeneityDict.lookupOrDefault<word>("smoothing", "smoothstep");
    const word transitionMode =
        heterogeneityDict.lookupOrDefault<word>("transitionMode", "blend");

    if (transmuralDistance.size() != STATES_.size())
    {
        FatalErrorInFunction
            << "Transmural distance field has " << transmuralDistance.size()
            << " values, but BuenoOrovio was configured with "
            << STATES_.size() << " integration points."
            << exit(FatalError);
    }

    const scalar endoMInterface =
        heterogeneityDict.lookupOrDefault<scalar>("endoMInterface", 0.3);
    const scalar mEpiInterface =
        heterogeneityDict.lookupOrDefault<scalar>("mEpiInterface", 0.7);
    const scalar transitionWidth =
        heterogeneityDict.lookupOrDefault<scalar>("transitionWidth", 0.1);

    ionicHeterogeneity::validateTransmuralBandConfig
    (
        endoMInterface,
        mEpiInterface,
        transitionWidth,
        smoothing,
        transitionMode
    );

    const scalarField endoConstants =
        ::constantsForTissue(3, dict());
    const scalarField mCellConstants =
        ::constantsForTissue(2, dict());
    const scalarField epiConstants =
        ::constantsForTissue(1, dict());

    HETEROGENEOUS_CONSTANTS_.clear();
    HETEROGENEOUS_CONSTANTS_.setSize(STATES_.size());

    forAll(transmuralDistance, integrationPtI)
    {
        const scalar rawT = transmuralDistance[integrationPtI];

        if (rawT < -SMALL || rawT > 1.0 + SMALL)
        {
            FatalErrorInFunction
                << "Transmural distance value t=" << rawT
                << " at integration point " << integrationPtI
                << " is outside the expected [0, 1] range."
                << exit(FatalError);
        }

        const scalar t = min(max(rawT, scalar(0.0)), scalar(1.0));
        scalarField mappedConstants(NUM_CONSTANTS, 0.0);
        const ionicHeterogeneity::TransmuralBandWeights weights =
            ionicHeterogeneity::transmuralBandWeights
            (
                t,
                endoMInterface,
                mEpiInterface,
                transitionWidth,
                smoothing,
                transitionMode
            );

        forAll(mappedConstants, constantI)
        {
            mappedConstants[constantI] =
                weights.endo*endoConstants[constantI]
              + weights.mCell*mCellConstants[constantI]
              + weights.epi*epiConstants[constantI];
        }

        HETEROGENEOUS_CONSTANTS_.set
        (
            integrationPtI,
            new scalarField(mappedConstants)
        );
    }

    Info<< "Configured BuenoOrovio transmural ionic heterogeneity using "
        << "endo/M interface " << endoMInterface
        << ", M/epi interface " << mEpiInterface
        << ", transitionWidth " << transitionWidth
        << ", transitionMode " << transitionMode
        << ", smoothing " << smoothing << "." << nl << endl;
}


void Foam::BuenoOrovio::solveODE
(
    const scalar stepStartTime,
    const scalar deltaT,
    const scalarField& Vm,
    scalarField& Im
)
{
    const scalar tStart = stepStartTime * 1000.0;
    const scalar tEnd   = (stepStartTime + deltaT) * 1000.0;
    const label sampleCell = sampleIntegrationPoint(STATES_.size());

    forAll(STATES_, integrationPtI)
    {
        scalarField& STATESI    = STATES_[integrationPtI];
        scalarField& ALGEBRAICI = ALGEBRAIC_[integrationPtI];
        scalarField& RATESI     = RATES_[integrationPtI];


        if (!solveVmWithinODESolver())
        {
            STATESI[0] = (Vm[integrationPtI] * 1000.0 + 84)/85.7;
        }
        scalar& step = ionicModel::step()[integrationPtI];

        step = min(step, deltaT * 1000.0);
        if (integrationPtI == sampleCell)
            {debugPrintFields(integrationPtI, tStart, tEnd, step);}

        activeIntegrationPoint_ = integrationPtI;
        odeSolver().solve(tStart, tEnd, STATESI, step);

        ::BuenoOroviocomputeVariables
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

        Im[integrationPtI] = ALGEBRAICI[Jion] * 85.7;

        //---------Currently with no use. -------------//

    }
}


void Foam::BuenoOrovio::derivatives
(
    const scalar t,
    const scalarField& y,
    scalarField& dydt
) const
{
    scalarField ALGEBRAIC_TMP(NUM_ALGEBRAIC, 0.0);

    ::BuenoOroviocomputeVariables
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

const char* const* Foam::BuenoOrovio::ioStateNames() const
{
    return BuenoOrovioSTATES_NAMES;
}

const char* const* Foam::BuenoOrovio::ioConstantNames() const
{
    return BuenoOrovioCONSTANTS_NAMES;
}

const char* const* Foam::BuenoOrovio::ioAlgebraicNames() const
{
    return BuenoOrovioALGEBRAIC_NAMES;
}

void Foam::BuenoOrovio::sweepCurrent
(
    const word& currentName,
    scalar Vmin,
    scalar Vmax,
    label nPts,
    const fileName& outputFile
) const
{
    const auto& depMap = BuenoOrovioDependencyMap();

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

        ::BuenoOroviocomputeVariables
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
            BuenoOrovioSTATES_NAMES, NUM_STATES,
            BuenoOrovioALGEBRAIC_NAMES, NUM_ALGEBRAIC,
            RATESI,
            sweepPlanCache
        );
    }

}

Foam::wordList Foam::BuenoOrovio::availableSweepCurrents() const
{
    return BuenoOrovioDependencyMap().toc();
}



// Coupling signals

bool Foam::BuenoOrovio::hasSignal(const CouplingSignal s) const
{
    return ionicModel::hasSignal(s);
}

Foam::scalar Foam::BuenoOrovio::signal
(
    const label i,
    const CouplingSignal s
) const
{
    return ionicModel::signal(i, s);
}
