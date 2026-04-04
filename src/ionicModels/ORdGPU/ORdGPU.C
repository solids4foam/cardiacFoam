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

#include "ORdGPU.H"
#include "ORd_2011Compact.H"
#include "ORd_2011Optimized.H"
#include "../ORd/ORdModelInfo.H"
#include "addToRunTimeSelectionTable.H"
#include "ionicModelIO.H"
#include <array>

namespace Foam
{
    defineTypeNameAndDebug(ORdGPU, 0);
    addToRunTimeSelectionTable
    (
        ionicModel,
        ORdGPU,
        dictionary
    );
}

namespace
{
    template<class Access>
    void stabilizeORdStateValues(Access access)
    {
        auto clamp01 = [&](const Foam::label stateI)
        {
            if (access(stateI) < 0.0)
            {
                access(stateI) = 0.0;
            }
            else if (access(stateI) > 1.0)
            {
                access(stateI) = 1.0;
            }
        };

        auto clampPositive = [&](const Foam::label stateI)
        {
            if (access(stateI) < Foam::SMALL)
            {
                access(stateI) = Foam::SMALL;
            }
        };

        clampPositive(CaMK_CaMKt);
        clampPositive(Na_i);
        clampPositive(Na_ss);
        clampPositive(K_i);
        clampPositive(K_ss);
        clampPositive(Ca_ss);
        clampPositive(Ca_nsr);
        clampPositive(Ca_jsr);
        clampPositive(Ca_i);

        clamp01(INa_m);
        clamp01(INa_hf);
        clamp01(INa_hs);
        clamp01(INa_j);
        clamp01(INa_hsp);
        clamp01(INa_jp);
        clamp01(INaL_mL);
        clamp01(INaL_hL);
        clamp01(INaL_hLp);
        clamp01(Ito_a);
        clamp01(Ito_iF);
        clamp01(Ito_iS);
        clamp01(Ito_ap);
        clamp01(Ito_iFp);
        clamp01(Ito_iSp);
        clamp01(ICaL_d);
        clamp01(ICaL_ff);
        clamp01(ICaL_fs);
        clamp01(ICaL_fcaf);
        clamp01(ICaL_fcas);
        clamp01(ICaL_jca);
        clamp01(ICaL_ffp);
        clamp01(ICaL_fcafp);
        clamp01(ICaL_nca);
        clamp01(IKr_xrf);
        clamp01(IKr_xrs);
        clamp01(IKs_xs1);
        clamp01(IKs_xs2);
        clamp01(IK1_xk1);

        if (access(ryr_Jrelnp) < 0.0)
        {
            access(ryr_Jrelnp) = 0.0;
        }

        if (access(ryr_Jrelp) < 0.0)
        {
            access(ryr_Jrelp) = 0.0;
        }
    }

    enum class ORdRLTauSource : unsigned char
    {
        none,
        lookup,
        constant,
        inverseLookup
    };

    struct ORdRLDispatchEntry
    {
        ORdRLTauSource tauSource;
        Foam::label tauIndex;
        Foam::label steadyStateIndex;
    };

    inline ORdRLDispatchEntry rlNone()
    {
        return {ORdRLTauSource::none, -1, -1};
    }

    inline ORdRLDispatchEntry rlLookup(const Foam::label tauIndex, const Foam::label steadyStateIndex)
    {
        return {ORdRLTauSource::lookup, tauIndex, steadyStateIndex};
    }

    inline ORdRLDispatchEntry rlConstant(const Foam::label constantIndex, const Foam::label steadyStateIndex)
    {
        return {ORdRLTauSource::constant, constantIndex, steadyStateIndex};
    }

    inline ORdRLDispatchEntry rlInverseLookup(const Foam::label tauIndex)
    {
        return {ORdRLTauSource::inverseLookup, tauIndex, -1};
    }

    const std::array<ORdRLDispatchEntry, NUM_STATES> ORdRushLarsenDispatch = []()
    {
        std::array<ORdRLDispatchEntry, NUM_STATES> entries{};
        entries.fill(rlNone());

        entries[IK1_xk1] = rlLookup(AV_txk1, AV_xk1ss);
        entries[IKr_xrf] = rlLookup(AV_txrf, AV_xrss);
        entries[IKr_xrs] = rlLookup(AV_txrs, AV_xrss);
        entries[IKs_xs1] = rlLookup(AV_txs1, AV_xs1ss);
        entries[IKs_xs2] = rlLookup(AV_txs2, AV_xs2ss);
        entries[INa_hf] = rlLookup(AV_thf, AV_hss);
        entries[INa_hs] = rlLookup(AV_ths, AV_hss);
        entries[INa_hsp] = rlLookup(AV_thsp, AV_hssp);
        entries[INa_m] = rlLookup(AV_tm, AV_mss);
        entries[INa_j] = rlLookup(AV_tj, AV_jss);
        entries[INa_jp] = rlLookup(AV_tjp, AV_jss);
        entries[Ito_a] = rlLookup(AV_ta, AV_ass);
        entries[Ito_ap] = rlLookup(AV_ta, AV_assp);
        entries[Ito_iF] = rlLookup(AV_tiF, AV_iss);
        entries[Ito_iS] = rlLookup(AV_tiS, AV_iss);
        entries[Ito_iFp] = rlLookup(AV_tiFp, AV_iss);
        entries[Ito_iSp] = rlLookup(AV_tiSp, AV_iss);
        entries[INaL_hL] = rlConstant(AC_thL, AV_hLss);
        entries[INaL_mL] = rlLookup(AV_tmL, AV_mLss);
        entries[INaL_hLp] = rlConstant(AC_thLp, AV_hLssp);
        entries[ICaL_d] = rlLookup(AV_td, AV_dss);
        entries[ICaL_ff] = rlLookup(AV_tff, AV_fss);
        entries[ICaL_fs] = rlLookup(AV_tfs, AV_fss);
        entries[ICaL_fcaf] = rlLookup(AV_tfcaf, AV_fcass);
        entries[ICaL_fcafp] = rlLookup(AV_tfcafp, AV_fcass);
        entries[ICaL_fcas] = rlLookup(AV_tfcas, AV_fcass);
        entries[ICaL_jca] = rlConstant(AC_tjca, AV_fcass);
        entries[ICaL_ffp] = rlLookup(AV_tffp, AV_fss);
        entries[ICaL_nca] = rlInverseLookup(AV_km2n);
        entries[ryr_Jrelnp] = rlLookup(AV_tau_rel, AV_Jrel_inf);
        entries[ryr_Jrelp] = rlLookup(AV_tau_relp, AV_Jrel_infp);

        return entries;
    }();

    enum class ORdCompactRLTauSource : unsigned char
    {
        none,
        lookup,
        constant,
        inverseLookup
    };

    struct ORdCompactRLDispatchEntry
    {
        ORdCompactRLTauSource tauSource;
        Foam::label tauIndex;
        Foam::label steadyStateIndex;
    };

    inline ORdCompactRLDispatchEntry compactRlNone()
    {
        return {ORdCompactRLTauSource::none, -1, -1};
    }

    inline ORdCompactRLDispatchEntry compactRlLookup
    (
        const Foam::label tauIndex,
        const Foam::label steadyStateIndex
    )
    {
        return {ORdCompactRLTauSource::lookup, tauIndex, steadyStateIndex};
    }

    inline ORdCompactRLDispatchEntry compactRlConstant
    (
        const Foam::label constantIndex,
        const Foam::label steadyStateIndex
    )
    {
        return
        {
            ORdCompactRLTauSource::constant,
            constantIndex,
            steadyStateIndex
        };
    }

    inline ORdCompactRLDispatchEntry compactRlInverseLookup
    (
        const Foam::label tauIndex
    )
    {
        return {ORdCompactRLTauSource::inverseLookup, tauIndex, -1};
    }

    const std::array<ORdCompactRLDispatchEntry, NUM_STATES>
        ORdCompactRushLarsenDispatch = []()
    {
        std::array<ORdCompactRLDispatchEntry, NUM_STATES> entries{};
        entries.fill(compactRlNone());

        entries[IK1_xk1] = compactRlLookup(ORdCS_txk1, ORdCS_xk1ss);
        entries[IKr_xrf] = compactRlLookup(ORdCS_txrf, ORdCS_xrss);
        entries[IKr_xrs] = compactRlLookup(ORdCS_txrs, ORdCS_xrss);
        entries[IKs_xs1] = compactRlLookup(ORdCS_txs1, ORdCS_xs1ss);
        entries[IKs_xs2] = compactRlLookup(ORdCS_txs2, ORdCS_xs2ss);
        entries[INa_hf] = compactRlLookup(ORdCS_thf, ORdCS_hss);
        entries[INa_hs] = compactRlLookup(ORdCS_ths, ORdCS_hss);
        entries[INa_hsp] = compactRlLookup(ORdCS_thsp, ORdCS_hssp);
        entries[INa_m] = compactRlLookup(ORdCS_tm, ORdCS_mss);
        entries[INa_j] = compactRlLookup(ORdCS_tj, ORdCS_jss);
        entries[INa_jp] = compactRlLookup(ORdCS_tjp, ORdCS_jss);
        entries[Ito_a] = compactRlLookup(ORdCS_ta, ORdCS_ass);
        entries[Ito_ap] = compactRlLookup(ORdCS_ta, ORdCS_assp);
        entries[Ito_iF] = compactRlLookup(ORdCS_tiF, ORdCS_iss);
        entries[Ito_iS] = compactRlLookup(ORdCS_tiS, ORdCS_iss);
        entries[Ito_iFp] = compactRlLookup(ORdCS_tiFp, ORdCS_iss);
        entries[Ito_iSp] = compactRlLookup(ORdCS_tiSp, ORdCS_iss);
        entries[INaL_hL] = compactRlConstant(AC_thL, ORdCS_hLss);
        entries[INaL_mL] = compactRlLookup(ORdCS_tmL, ORdCS_mLss);
        entries[INaL_hLp] = compactRlConstant(AC_thLp, ORdCS_hLssp);
        entries[ICaL_d] = compactRlLookup(ORdCS_td, ORdCS_dss);
        entries[ICaL_ff] = compactRlLookup(ORdCS_tff, ORdCS_fss);
        entries[ICaL_fs] = compactRlLookup(ORdCS_tfs, ORdCS_fss);
        entries[ICaL_fcaf] = compactRlLookup(ORdCS_tfcaf, ORdCS_fcass);
        entries[ICaL_fcafp] = compactRlLookup(ORdCS_tfcafp, ORdCS_fcass);
        entries[ICaL_fcas] = compactRlLookup(ORdCS_tfcas, ORdCS_fcass);
        entries[ICaL_jca] = compactRlConstant(AC_tjca, ORdCS_fcass);
        entries[ICaL_ffp] = compactRlLookup(ORdCS_tffp, ORdCS_fss);
        entries[ICaL_nca] = compactRlInverseLookup(ORdCS_km2n);
        entries[ryr_Jrelnp] = compactRlLookup(ORdCS_tau_rel, ORdCS_Jrel_inf);
        entries[ryr_Jrelp] = compactRlLookup(ORdCS_tau_relp, ORdCS_Jrel_infp);

        return entries;
    }();
}

Foam::ORdGPU::ORdGPU
(
    const dictionary& dict,
    const label num,
    const scalar initialDeltaT,
    const Switch solveVmWithinODESolver,
    const FullEvaluationBackend fullEvaluationBackend
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
        ORdFamilyInfo()
    ),
    fullEvaluationBackend_(fullEvaluationBackend)
{
    ionicModel::setTissueFromDict();
    setHotPathSupportSize(NUM_ORD_COMPACT_SUPPORT);

    double initialRates[NUM_STATES] = {0.0};
    double initialStates[NUM_STATES] = {0.0};

    ORdinitConsts
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

Foam::ORdGPU::ORdGPU
(
    const dictionary& dict,
    const label num,
    const scalar initialDeltaT,
    const Switch solveVmWithinODESolver
)
:
    ORdGPU
    (
        dict,
        num,
        initialDeltaT,
        solveVmWithinODESolver,
        FullEvaluationBackend::generated
    )
{}

Foam::ORdGPU::~ORdGPU()
{}

Foam::List<Foam::word> Foam::ORdGPU::supportedTissueTypes() const
{
    return {"myocyte"};
}

void Foam::ORdGPU::solveODE
(
    const scalar stepStartTime,
    const scalar deltaT,
    const scalarField& Vm,
    scalarField& Im
)
{
    solveODEImpl(*this, stepStartTime, deltaT, Vm, Im);
}

Foam::scalar Foam::ORdGPU::ionicCurrentFromHotPathSupport
(
    const scalarUList& supportValues
) const
{
    return supportValues.size() == 1
        ? supportValues[0]
        : supportValues[ORdCS_Iion_cm];
}

void Foam::ORdGPU::stabilizeCellState(const label cellI) const
{
    stabilizeORdStateValues
    (
        [&](const label stateI) -> scalar&
        {
            return state(cellI, stateI);
        }
    );
}

void Foam::ORdGPU::stabilizeStateValues(scalarUList& stateValues) const
{
    stabilizeORdStateValues
    (
        [&](const label stateI) -> scalar&
        {
            return stateValues[stateI];
        }
    );
}

void Foam::ORdGPU::evaluateState
(
    const scalar modelTime,
    const scalarUList& stateValues,
    scalarUList& rateValues,
    scalarUList& algebraicValues
) const
{
    switch (fullEvaluationBackend_)
    {
        case FullEvaluationBackend::optimized:
            ORdcomputeVariablesOptimized
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
            break;

        case FullEvaluationBackend::generated:
        default:
            ORdcomputeVariables
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
            break;
    }
}

void Foam::ORdGPU::evaluateHotPathState
(
    const scalar modelTime,
    const scalarUList& stateValues,
    scalarUList& rateValues,
    scalarUList& supportValues
) const
{
    ORdcomputeVariablesCompact
    (
        modelTime,
        CONSTANTS_.data(),
        rateValues.data(),
        const_cast<scalarUList&>(stateValues).data(),
        supportValues.data(),
        supportValues.size(),
        tissue(),
        solveVmWithinODESolver(),
        stimulusProtocol()
    );
}

bool Foam::ORdGPU::rushLarsenParameters
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

    const ORdRLDispatchEntry& entry = ORdRushLarsenDispatch[stateI];

    switch (entry.tauSource)
    {
        case ORdRLTauSource::lookup:
            tau = algebraicValues[entry.tauIndex];
            steadyState = algebraicValues[entry.steadyStateIndex];
            return tau > VSMALL;

        case ORdRLTauSource::constant:
            tau = CONSTANTS_[entry.tauIndex];
            steadyState = algebraicValues[entry.steadyStateIndex];
            return tau > VSMALL;

        case ORdRLTauSource::inverseLookup:
            tau = 1.0/max(algebraicValues[entry.tauIndex], scalar(VSMALL));
            steadyState = stateValues[ICaL_nca] + rateValues[ICaL_nca]*tau;
            return std::isfinite(steadyState) && tau > VSMALL;

        case ORdRLTauSource::none:
        default:
            return false;
    }
}

bool Foam::ORdGPU::rushLarsenParametersFromHotPathSupport
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

    const ORdCompactRLDispatchEntry& entry = ORdCompactRushLarsenDispatch[stateI];

    switch (entry.tauSource)
    {
        case ORdCompactRLTauSource::lookup:
            tau = supportValues[entry.tauIndex];
            steadyState = supportValues[entry.steadyStateIndex];
            return tau > VSMALL;

        case ORdCompactRLTauSource::constant:
            tau = CONSTANTS_[entry.tauIndex];
            steadyState = supportValues[entry.steadyStateIndex];
            return tau > VSMALL;

        case ORdCompactRLTauSource::inverseLookup:
            tau = 1.0/max(supportValues[entry.tauIndex], scalar(VSMALL));
            steadyState = stateValues[ICaL_nca] + rateValues[ICaL_nca]*tau;
            return std::isfinite(steadyState) && tau > VSMALL;

        case ORdCompactRLTauSource::none:
        default:
            return false;
    }
}

void Foam::ORdGPU::derivatives
(
    const scalar t,
    const scalarField& y,
    scalarField& dydt
) const
{
    scalarField algebraics(NUM_ALGEBRAIC, 0.0);
    evaluateState(t, y, dydt, algebraics);
}

void Foam::ORdGPU::sweepCurrent
(
    const word& currentName,
    scalar Vmin,
    scalar Vmax,
    label nPts,
    const fileName& outputFile
) const
{
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
        statesLocal[membrane_V] = voltage;
        ratesLocal = 0.0;
        algebraicsLocal = 0.0;

        ORdcomputeVariables
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
            ORdSTATES_NAMES,
            NUM_STATES,
            ORdALGEBRAIC_NAMES,
            NUM_ALGEBRAIC,
            ratesLocal,
            sweepPlanCache
        );
    }
}

Foam::wordList Foam::ORdGPU::availableSweepCurrents() const
{
    return ORdDependencyMap().toc();
}

// ************************************************************************* //
