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

#include "BuenoOrovioGPU.H"
#include "BuenoOrovio_2008.H"
#include "../BuenoOrovio/BuenoOrovioModelInfo.H"
#include "gpuMath.H"
#include "addToRunTimeSelectionTable.H"
#include "ionicModelIO.H"

namespace Foam
{
    defineTypeNameAndDebug(BuenoOrovioGPU, 0);
    addToRunTimeSelectionTable
    (
        ionicModel,
        BuenoOrovioGPU,
        dictionary
    );
}

Foam::BuenoOrovioGPU::BuenoOrovioGPU
(
    const dictionary& dict,
    const label num,
    const scalar initialDeltaT,
    const Switch solveVmWithinODESolver
)
:
    configuredBatchedIonicModel
    (
        dict,
        num,
        initialDeltaT,
        solveVmWithinODESolver,
        NUM_STATES,
        NUM_ALGEBRAIC,
        BuenoOrovioFamilyInfo()
    )
{
    ionicModel::setTissueFromDict();

    double initialRates[NUM_STATES] = {0.0};
    double initialStates[NUM_STATES] = {0.0};

    BuenoOrovioinitConsts
    (
        CONSTANTS_.data(),
        initialRates,
        initialStates,
        tissue(),
        dict
    );

    for (label cellI = 0; cellI < nCells(); ++cellI)
    {
        for (label stateI = 0; stateI < NUM_STATES; ++stateI)
        {
            state(cellI, stateI) = initialStates[stateI];
            rate(cellI, stateI) = initialRates[stateI];
        }
    }

    configurePersistentAlgebraics();
    syncAllToIO();

    if (!utilitiesMode())
    {
        setStimulusProtocolFromDict(dict);
    }
}

Foam::BuenoOrovioGPU::~BuenoOrovioGPU()
{}

Foam::List<Foam::word> Foam::BuenoOrovioGPU::supportedTissueTypes() const
{
    return {"endocardialCells", "mCells", "epicardialCells"};
}

void Foam::BuenoOrovioGPU::evaluateState
(
    const scalar modelTime,
    const scalarUList& stateValues,
    scalarUList& rateValues,
    scalarUList& algebraicValues
) const
{
    BuenoOroviocomputeVariables
    (
        modelTime,
        const_cast<double*>(CONSTANTS_.cdata()),
        rateValues.data(),
        const_cast<scalarUList&>(stateValues).data(),
        algebraicValues.data(),
        tissue(),
        solveVmWithinODESolver(),
        stimulusProtocol()
    );
}

bool Foam::BuenoOrovioGPU::rushLarsenParameters
(
    const label stateI,
    const scalarUList& stateValues,
    const scalarUList& rateValues,
    const scalarUList& algebraicValues,
    scalar& steadyState,
    scalar& tau
) const
{
    using Foam::smoothHeaviside;

    const scalar cellV = stateValues[u];

    switch (stateI)
    {
        case v:
        {
            const scalar hV = smoothHeaviside(cellV - CONSTANTS_[thetaV]);
            const scalar invTau =
                (1.0 - hV)/algebraicValues[tauVMinus]
              + hV/CONSTANTS_[tauVPlus];

            if (invTau <= VSMALL)
            {
                return false;
            }

            tau = 1.0/invTau;
            steadyState = stateValues[v] + rateValues[v]*tau;
            return true;
        }

        case w:
        {
            const scalar hW = smoothHeaviside(cellV - CONSTANTS_[thetaW]);
            const scalar invTau =
                (1.0 - hW)/algebraicValues[tauWMinus]
              + hW/CONSTANTS_[tauWPlus];

            if (invTau <= VSMALL)
            {
                return false;
            }

            tau = 1.0/invTau;
            steadyState = stateValues[w] + rateValues[w]*tau;
            return true;
        }

        case s:
        {
            tau = algebraicValues[tauS];

            if (tau <= VSMALL)
            {
                return false;
            }

            steadyState = stateValues[s] + rateValues[s]*tau;
            return true;
        }

        default:
            return false;
    }
}

void Foam::BuenoOrovioGPU::derivatives
(
    const scalar t,
    const scalarField& y,
    scalarField& dydt
) const
{
    scalarField algebraics(NUM_ALGEBRAIC, 0.0);

    evaluateState(t, y, dydt, algebraics);
}

void Foam::BuenoOrovioGPU::sweepCurrent
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

    scalarField statesLocal(NUM_STATES, 0.0);
    scalarField ratesLocal(NUM_STATES, 0.0);
    scalarField algebraicsLocal(NUM_ALGEBRAIC, 0.0);
    ionicModelIO::SelectedMapCache sweepPlanCache;

    for (label stateI = 0; stateI < NUM_STATES; ++stateI)
    {
        statesLocal[stateI] = state(0, stateI);
    }

    for (label i = 0; i < nPts; ++i)
    {
        const scalar V = Vmin + (Vmax - Vmin) * scalar(i) / (nPts - 1);
        statesLocal[u] = V;
        ratesLocal = 0.0;
        algebraicsLocal = 0.0;

        evaluateState(0.0, statesLocal, ratesLocal, algebraicsLocal);

        ionicModelIO::writeOneSweepRow
        (
            os,
            V,
            deps,
            statesLocal,
            algebraicsLocal,
            BuenoOrovioSTATES_NAMES,
            NUM_STATES,
            BuenoOrovioALGEBRAIC_NAMES,
            NUM_ALGEBRAIC,
            ratesLocal,
            sweepPlanCache
        );
    }
}

Foam::wordList Foam::BuenoOrovioGPU::availableSweepCurrents() const
{
    return BuenoOrovioDependencyMap().toc();
}

// ************************************************************************* //
