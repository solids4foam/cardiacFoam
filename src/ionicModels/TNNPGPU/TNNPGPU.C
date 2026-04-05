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

#include "TNNPGPU.H"
#include "TNNP_2004.H"
#include "../TNNP/TNNPModelInfo.H"
#include "addToRunTimeSelectionTable.H"
#include "ionicModelIO.H"
#include <array>

namespace Foam
{
    defineTypeNameAndDebug(TNNPGPU, 0);
    addToRunTimeSelectionTable
    (
        ionicModel,
        TNNPGPU,
        dictionary
    );
}

namespace
{
    template<class Access>
    void stabilizeTNNPStateValues(Access access)
    {
        auto clampRange =
            [&](const Foam::label stateI, const Foam::scalar upper)
        {
            if (access(stateI) < 0.0)
            {
                access(stateI) = 0.0;
            }
            else if (access(stateI) > upper)
            {
                access(stateI) = upper;
            }
        };

        if (access(K_i) < Foam::SMALL)
        {
            access(K_i) = Foam::SMALL;
        }

        if (access(Na_i) < Foam::SMALL)
        {
            access(Na_i) = Foam::SMALL;
        }

        if (access(Ca_i) < Foam::SMALL)
        {
            access(Ca_i) = Foam::SMALL;
        }

        if (access(Ca_SR) < Foam::SMALL)
        {
            access(Ca_SR) = Foam::SMALL;
        }

        clampRange(Xr1, 1.0);
        clampRange(Xr2, 1.0);
        clampRange(Xs, 1.0);
        clampRange(m, 1.0);
        clampRange(h, 1.0);
        clampRange(j, 1.0);
        clampRange(d, 1.0);
        clampRange(f, 1.0);
        clampRange(fCa, 1.0);
        clampRange(s, 1.1);
        clampRange(r, 1.0);
        clampRange(g, 1.0);
    }

    enum class TNNPRLTauSource : unsigned char
    {
        none,
        lookup,
        constant
    };

    struct TNNPRLDispatchEntry
    {
        TNNPRLTauSource tauSource;
        Foam::label tauIndex;
        Foam::label supportIndex;
    };

    inline TNNPRLDispatchEntry rlNone()
    {
        return {TNNPRLTauSource::none, -1, -1};
    }

    inline TNNPRLDispatchEntry rlLookup
    (
        const Foam::label tauIndex,
        const Foam::label supportIndex
    )
    {
        return {TNNPRLTauSource::lookup, tauIndex, supportIndex};
    }

    inline TNNPRLDispatchEntry rlConstant(const Foam::label constantIndex)
    {
        return {TNNPRLTauSource::constant, constantIndex, -1};
    }

    const std::array<TNNPRLDispatchEntry, NUM_STATES> TNNPRushLarsenDispatch = []()
    {
        std::array<TNNPRLDispatchEntry, NUM_STATES> entries{};
        entries.fill(rlNone());

        entries[Xr1] = rlLookup(tau_xr1, TNNPCS_tau_xr1);
        entries[Xr2] = rlLookup(tau_xr2, TNNPCS_tau_xr2);
        entries[Xs] = rlLookup(tau_xs, TNNPCS_tau_xs);
        entries[m] = rlLookup(tau_m, TNNPCS_tau_m);
        entries[h] = rlLookup(tau_h, TNNPCS_tau_h);
        entries[j] = rlLookup(tau_j, TNNPCS_tau_j);
        entries[d] = rlLookup(tau_d, TNNPCS_tau_d);
        entries[f] = rlLookup(tau_f, TNNPCS_tau_f);
        entries[fCa] = rlConstant(tau_fCa);
        entries[s] = rlLookup(tau_s, TNNPCS_tau_s);
        entries[r] = rlLookup(tau_r, TNNPCS_tau_r);
        entries[g] = rlConstant(tau_g);

        return entries;
    }();
}

Foam::TNNPGPU::TNNPGPU
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
        TNNPFamilyInfo()
    )
{
    ionicModel::setTissueFromDict();
    setHotPathSupportSize(NUM_TNNP_COMPACT_SUPPORT);

    double initialRates[NUM_STATES] = {0.0};
    double initialStates[NUM_STATES] = {0.0};

    TNNPinitConsts
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

Foam::TNNPGPU::~TNNPGPU()
{}

Foam::List<Foam::word> Foam::TNNPGPU::supportedTissueTypes() const
{
    return {"endocardialCells", "mCells", "epicardialCells"};
}

void Foam::TNNPGPU::solveODE
(
    const scalar stepStartTime,
    const scalar deltaT,
    const scalarField& Vm,
    scalarField& Im
)
{
    solveODEImpl(*this, stepStartTime, deltaT, Vm, Im);
}

Foam::scalar Foam::TNNPGPU::ionicCurrentFromHotPathSupport
(
    const scalarUList& supportValues
) const
{
    return supportValues[TNNPCS_Iion_cm];
}

void Foam::TNNPGPU::stabilizeCellState(const label cellI) const
{
    stabilizeTNNPStateValues
    (
        [&](const label stateI) -> scalar&
        {
            return state(cellI, stateI);
        }
    );
}

void Foam::TNNPGPU::stabilizeStateValues(scalarUList& stateValues) const
{
    stabilizeTNNPStateValues
    (
        [&](const label stateI) -> scalar&
        {
            return stateValues[stateI];
        }
    );
}

void Foam::TNNPGPU::evaluateState
(
    const scalar modelTime,
    const scalarUList& stateValues,
    scalarUList& rateValues,
    scalarUList& algebraicValues
) const
{
    TNNPcomputeRates
    (
        modelTime,
        CONSTANTS_.data(),
        rateValues.data(),
        const_cast<scalarUList&>(stateValues).data(),
        algebraicValues.data(),
        tissue(),
        solveVmWithinODESolver(),
        stimulusProtocol()
    );
}

void Foam::TNNPGPU::evaluateHotPathState
(
    const scalar modelTime,
    const scalarUList& stateValues,
    scalarUList& rateValues,
    scalarUList& supportValues
) const
{
    TNNPcomputeRatesCompact
    (
        modelTime,
        CONSTANTS_.data(),
        rateValues.data(),
        const_cast<scalarUList&>(stateValues).data(),
        supportValues.data(),
        tissue(),
        solveVmWithinODESolver(),
        stimulusProtocol()
    );
}

bool Foam::TNNPGPU::rushLarsenParameters
(
    const label stateI,
    const scalarUList& stateValues,
    const scalarUList& rateValues,
    const scalarUList& algebraicValues,
    scalar& steadyState,
    scalar& tau
) const
{
    if (stateI < 0 || stateI >= NUM_STATES)
    {
        return false;
    }

    const TNNPRLDispatchEntry& entry = TNNPRushLarsenDispatch[stateI];

    switch (entry.tauSource)
    {
        case TNNPRLTauSource::lookup:
            tau = algebraicValues[entry.tauIndex];
            break;

        case TNNPRLTauSource::constant:
            tau = CONSTANTS_[entry.tauIndex];
            break;

        case TNNPRLTauSource::none:
        default:
            return false;
    }

    if (tau <= VSMALL)
    {
        return false;
    }

    steadyState = stateValues[stateI] + rateValues[stateI]*tau;
    return std::isfinite(steadyState) && std::isfinite(tau);
}

bool Foam::TNNPGPU::rushLarsenParametersFromHotPathSupport
(
    const label stateI,
    const scalarUList& stateValues,
    const scalarUList& rateValues,
    const scalarUList& supportValues,
    scalar& steadyState,
    scalar& tau
) const
{
    if (stateI < 0 || stateI >= NUM_STATES)
    {
        return false;
    }

    const TNNPRLDispatchEntry& entry = TNNPRushLarsenDispatch[stateI];

    switch (entry.tauSource)
    {
        case TNNPRLTauSource::lookup:
            tau = supportValues[entry.supportIndex];
            break;

        case TNNPRLTauSource::constant:
            tau = CONSTANTS_[entry.tauIndex];
            break;

        case TNNPRLTauSource::none:
        default:
            return false;
    }

    if (tau <= VSMALL)
    {
        return false;
    }

    steadyState = stateValues[stateI] + rateValues[stateI]*tau;
    return std::isfinite(steadyState) && std::isfinite(tau);
}

void Foam::TNNPGPU::derivatives
(
    const scalar t,
    const scalarField& y,
    scalarField& dydt
) const
{
    scalarField algebraics(NUM_ALGEBRAIC, 0.0);

    evaluateState(t, y, dydt, algebraics);
}

void Foam::TNNPGPU::sweepCurrent
(
    const word& currentName,
    scalar Vmin,
    scalar Vmax,
    label nPts,
    const fileName& outputFile
) const
{
    const auto& depMap = TNNPDependencyMap();

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
        const scalar voltage = Vmin + (Vmax - Vmin) * scalar(i) / (nPts - 1);
        statesLocal[V] = voltage;
        ratesLocal = 0.0;
        algebraicsLocal = 0.0;

        TNNPcomputeVariables
        (
            0.0,
            CONSTANTS_.data(),
            ratesLocal.data(),
            statesLocal.data(),
            algebraicsLocal.data(),
            tissue(),
            solveVmWithinODESolver(),
            stimulusProtocol()
        );

        ionicModelIO::writeOneSweepRow
        (
            os,
            voltage,
            deps,
            statesLocal,
            algebraicsLocal,
            TNNP_STATES_NAMES,
            NUM_STATES,
            TNNP_ALGEBRAIC_NAMES,
            NUM_ALGEBRAIC,
            ratesLocal,
            sweepPlanCache
        );
    }
}

Foam::wordList Foam::TNNPGPU::availableSweepCurrents() const
{
    return TNNPDependencyMap().toc();
}

// ************************************************************************* //
