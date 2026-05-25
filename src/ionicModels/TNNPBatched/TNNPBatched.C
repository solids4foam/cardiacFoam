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

#include "TNNPBatched.H"
#include "TNNP_2004Batch.H"
#include "TNNP_2004Names.H"

namespace Foam
{
    const ionicModelFamilyInfo& TNNPFamilyInfo()
    {
        static const ionicModelFamilyInfo info
        {
            NUM_CONSTANTS,
            NUM_STATES,
            NUM_ALGEBRAIC,
            TNNP_CONSTANTS_NAMES,
            TNNP_STATES_NAMES,
            TNNP_ALGEBRAIC_NAMES,
            V,
            1000.0,
            1000.0,
            0.0,
            nullptr
        };
        return info;
    }
}
#include "addToRunTimeSelectionTable.H"
#include "ionicModelIO.H"
#include "stimulusIO.H"
#include "Pstream.H"
#include <array>

#ifdef HAS_CUDA
#include <cuda_runtime.h>

namespace Foam
{
    void launchTnnpBatchKernel
    (
        double t,
        const double* d_CONSTANTS,
        int N,
        const double* d_STATES,
        double* d_RATES,
        double* d_SUPPORT,
        int tissueFlag,
        bool solveVm,
        StimulusProtocolPOD stimulus
    );
    void launchTnnpEulerStepKernel
    (
        double* d_STATES,
        const double* d_RATES,
        double dt,
        int N,
        int nStates,
        bool solveVm,
        int vmStateI
    );
    void launchTnnpScaleIonKernel
    (
        const double* d_SUPPORT,
        double* d_Im,
        double scale,
        int N,
        int IionBase
    );
}

#endif // HAS_CUDA

namespace Foam
{
    defineTypeNameAndDebug(TNNPBatched, 0);
    addToRunTimeSelectionTable
    (
        ionicModel,
        TNNPBatched,
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

        entries[Xr1] = rlLookup(tau_xr1, Foam::TNNP_BATCH_SUPPORT_tau_xr1);
        entries[Xr2] = rlLookup(tau_xr2, Foam::TNNP_BATCH_SUPPORT_tau_xr2);
        entries[Xs] = rlLookup(tau_xs, Foam::TNNP_BATCH_SUPPORT_tau_xs);
        entries[m] = rlLookup(tau_m, Foam::TNNP_BATCH_SUPPORT_tau_m);
        entries[h] = rlLookup(tau_h, Foam::TNNP_BATCH_SUPPORT_tau_h);
        entries[j] = rlLookup(tau_j, Foam::TNNP_BATCH_SUPPORT_tau_j);
        entries[d] = rlLookup(tau_d, Foam::TNNP_BATCH_SUPPORT_tau_d);
        entries[f] = rlLookup(tau_f, Foam::TNNP_BATCH_SUPPORT_tau_f);
        entries[fCa] = rlConstant(tau_fCa);
        entries[s] = rlLookup(tau_s, Foam::TNNP_BATCH_SUPPORT_tau_s);
        entries[r] = rlLookup(tau_r, Foam::TNNP_BATCH_SUPPORT_tau_r);
        entries[g] = rlConstant(tau_g);

        return entries;
    }();
}

Foam::TNNPBatched::TNNPBatched
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
    ),
    useSoAEvaluator_
    (
        dict.lookupOrDefault<Switch>("useSoAEvaluator", false)
    ),
    stimulusPOD_()
#ifdef HAS_CUDA
  , useDevice_(false)
#endif
{
    ionicModel::setTissueFromDict();
    setHotPathSupportSize(NUM_TNNP_BATCH_SUPPORT);

#ifdef HAS_CUDA
    if (useSoAEvaluator_)
    {
        int nDevices = 0;
        cudaError_t err = cudaGetDeviceCount(&nDevices);
        if (err == cudaSuccess && nDevices > 0)
        {
            const int rank = Pstream::myProcNo();
            CARDIAC_CUDA_CHECK(cudaSetDevice(rank % nDevices));
            useDevice_ = true;
            Info<< "TNNPBatched: rank " << rank << " using CUDA device "
                << (rank % nDevices) << " of " << nDevices << nl;
        }
        else
        {
            WarningInFunction
                << "TNNPBatched: useSoAEvaluator is on but no CUDA "
                << "device is visible (" << cudaGetErrorString(err)
                << "); falling back to host SIMD path." << nl;
        }
    }
#endif

    if (useSoAEvaluator_)
    {
        const word integrator =
            dict.lookupOrDefault<word>("batchedIntegrator", "rlGatesHeunRest");
        if (integrator != "euler")
        {
            WarningInFunction
                << "useSoAEvaluator is enabled for " << type()
                << " but `batchedIntegrator " << integrator
                << "` is not honored on the batched path; reverting to "
                << "explicit Euler. Set `batchedIntegrator euler;` to "
                << "silence this warning." << nl;
        }
    }

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

Foam::TNNPBatched::~TNNPBatched()
{
#ifdef HAS_CUDA
    if (useDevice_)
    {
        cuda_.free();
    }
#endif
}




void Foam::TNNPBatched::prepareIOAccess
(
    const wordList& requestedNames,
    const bool needsAlgebraics
) const
{
#ifdef HAS_CUDA
    if (useDevice_ && cuda_.deviceDirty)
    {
        cuda_.syncStatesDeviceToHost
        (
            statesSoAData(),
            static_cast<std::size_t>(NUM_STATES),
            static_cast<std::size_t>(nCells())
        );
        cuda_.syncSupportDeviceToHost
        (
            supportSoAData(),
            static_cast<std::size_t>(NUM_TNNP_BATCH_SUPPORT),
            static_cast<std::size_t>(nCells())
        );
    }
#endif
    configuredBatchedIonicModel::prepareIOAccess(requestedNames, needsAlgebraics);
}


void Foam::TNNPBatched::importFields
(
    const volScalarField& Vm,
    const wordList& fieldNames,
    const PtrList<volScalarField>& inFields
)
{
    configuredBatchedIonicModel::importFields(Vm, fieldNames, inFields);
#ifdef HAS_CUDA
    if (useDevice_)
    {
        cuda_.hostDirty = true;
    }
#endif
}


Foam::List<Foam::word> Foam::TNNPBatched::supportedTissueTypes() const
{
    return {"endocardialCells", "mCells", "epicardialCells"};
}

void Foam::TNNPBatched::solveODE
(
    const scalar stepStartTime,
    const scalar deltaT,
    const scalarField& Vm,
    scalarField& Im
)
{
    if (useSoAEvaluator_ && !utilitiesMode())
    {
        solveBatched(stepStartTime, deltaT, Vm, Im);
        return;
    }

    solveODEImpl(*this, stepStartTime, deltaT, Vm, Im);
}


void Foam::TNNPBatched::solveBatched
(
    const scalar stepStartTime,
    const scalar deltaT,
    const scalarField& Vm,
    scalarField& Im
)
{

    const label N = nCells();
    const label nSub = nSubsteps();
    const scalar dtModel = deltaT*timeScaleFactor();
    const scalar dtSubstep = dtModel/scalar(nSub);
    const scalar tStart = stepStartTime*timeScaleFactor();
    const bool   solveVm = solveVmWithinODESolver();

    stimulusPOD_ = stimulusIO::toPOD(stimulusProtocol());

    if (!solveVm)
    {
        for (label cellI = 0; cellI < N; ++cellI)
        {
            state(cellI, V) = vmToState(Vm[cellI]);
        }
    }

    scalar* const STATES_SoA  = statesSoAData();
    scalar* const RATES_SoA   = ratesSoAData();
    scalar* const SUPPORT_SoA = supportSoAData();
    const scalar* const CONSTS = CONSTANTS_.cdata();

    markIODirty();
    setIOEvaluationModelTime(tStart);

#ifdef HAS_CUDA
    if (useDevice_)
    {
        const int tFlag = static_cast<int>(tissue());

        cuda_.allocate
        (
            static_cast<std::size_t>(N),
            static_cast<std::size_t>(NUM_STATES),
            static_cast<std::size_t>(NUM_TNNP_BATCH_SUPPORT),
            static_cast<std::size_t>(CONSTANTS_.size())
        );
        cuda_.uploadConstants(CONSTANTS_.cdata(), CONSTANTS_.size());

        if (cuda_.hostDirty)
        {
            cuda_.syncStatesHostToDevice
            (
                statesSoAData(),
                static_cast<std::size_t>(NUM_STATES),
                static_cast<std::size_t>(N)
            );
        }

        for (label sub = 0; sub < nSub; ++sub)
        {
            const scalar tSub = tStart + scalar(sub)*dtSubstep;
            launchTnnpBatchKernel
            (
                tSub, cuda_.d_constants, static_cast<int>(N),
                cuda_.d_states, cuda_.d_rates, cuda_.d_support,
                tFlag, solveVm, stimulusPOD_
            );
            launchTnnpEulerStepKernel
            (
                cuda_.d_states, cuda_.d_rates,
                static_cast<double>(dtSubstep),
                static_cast<int>(N),
                static_cast<int>(NUM_STATES),
                solveVm,
                static_cast<int>(V)
            );
        }

        launchTnnpBatchKernel
        (
            tStart + dtModel, cuda_.d_constants, static_cast<int>(N),
            cuda_.d_states, cuda_.d_rates, cuda_.d_support,
            tFlag, solveVm, stimulusPOD_
        );
        launchTnnpScaleIonKernel
        (
            cuda_.d_support, cuda_.d_Im, 1.0,        // TNNP Iion is direct
            static_cast<int>(N),
            static_cast<int>(TNNP_BATCH_SUPPORT_Iion_cm)
        );
        cuda_.downloadIm(Im.data(), static_cast<std::size_t>(N));
        cuda_.deviceDirty = true;
        setIOEvaluationModelTime(tStart + dtModel);
        return;
    }
#endif // HAS_CUDA

    for (label sub = 0; sub < nSub; ++sub)
    {
        const scalar tSub = tStart + scalar(sub)*dtSubstep;
        TNNPComputeRatesBatch
        (
            tSub, CONSTS, static_cast<int>(N), 0, static_cast<int>(N),
            STATES_SoA, RATES_SoA, SUPPORT_SoA,
            solveVm, stimulusPOD_
        );

        for (label stateI = 0; stateI < NUM_STATES; ++stateI)
        {
            if (!solveVm && stateI == V) continue;
            const label base = stateI*N;
#ifdef _OPENMP
            #pragma omp simd
#endif
            for (label cellI = 0; cellI < N; ++cellI)
            {
                STATES_SoA[base + cellI] +=
                    dtSubstep*RATES_SoA[base + cellI];
            }
        }
    }

    TNNPComputeRatesBatch
    (
        tStart + dtModel, CONSTS,
        static_cast<int>(N), 0, static_cast<int>(N),
        STATES_SoA, RATES_SoA, SUPPORT_SoA,
        solveVm, stimulusPOD_
    );

    if (SUPPORT_SoA != nullptr)
    {
        const label IionBase = TNNP_BATCH_SUPPORT_Iion_cm*N;
        for (label cellI = 0; cellI < N; ++cellI)
        {
            Im[cellI] = SUPPORT_SoA[IionBase + cellI];
        }
    }

    setIOEvaluationModelTime(tStart + dtModel);
}

Foam::scalar Foam::TNNPBatched::ionicCurrentFromHotPathSupport
(
    const scalarUList& supportValues
) const
{
    return supportValues[TNNP_BATCH_SUPPORT_Iion_cm];
}

void Foam::TNNPBatched::stabilizeCellState(const label cellI) const
{
    stabilizeTNNPStateValues
    (
        [&](const label stateI) -> scalar&
        {
            return state(cellI, stateI);
        }
    );
}

void Foam::TNNPBatched::stabilizeStateValues(scalarUList& stateValues) const
{
    stabilizeTNNPStateValues
    (
        [&](const label stateI) -> scalar&
        {
            return stateValues[stateI];
        }
    );
}

void Foam::TNNPBatched::evaluateState
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
        solveVmWithinODESolver(),
        stimulusProtocol()
    );
}

void Foam::TNNPBatched::evaluateHotPathState
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
                solveVmWithinODESolver(),
        stimulusProtocol()
    );
}

bool Foam::TNNPBatched::rushLarsenParameters
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

bool Foam::TNNPBatched::rushLarsenParametersFromHotPathSupport
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

void Foam::TNNPBatched::derivatives
(
    const scalar t,
    const scalarField& y,
    scalarField& dydt
) const
{
    scalarField algebraics(NUM_ALGEBRAIC, 0.0);

    evaluateState(t, y, dydt, algebraics);
}



// ************************************************************************* //
