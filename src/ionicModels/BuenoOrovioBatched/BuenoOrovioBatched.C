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

#include "BuenoOrovioBatched.H"
#include "BuenoOrovio_2008.H"
#include "gpuMath.H"
#include "batchedRushLarsenEntry.H"
#include <array>

namespace
{
    Foam::scalar buenoOrovioTransformedVm(const Foam::scalarField& S)
    {
        return S[0] * 85.7 - 84.0;
    }

    const std::array<Foam::batchedRushLarsenEntry, NUM_STATES>
    BuenoOrovioRushLarsenDispatch = []()
    {
        std::array<Foam::batchedRushLarsenEntry, NUM_STATES> t{};
        t.fill(Foam::rlNone());

        t[v] = Foam::rlSupport
        (
            Foam::BO_BATCH_SUPPORT_tau_v,
            Foam::BO_BATCH_SUPPORT_gInf_v
        );
        t[w] = Foam::rlSupport
        (
            Foam::BO_BATCH_SUPPORT_tau_w,
            Foam::BO_BATCH_SUPPORT_gInf_w
        );
        t[s] = Foam::rlSupport
        (
            Foam::BO_BATCH_SUPPORT_tau_s,
            Foam::BO_BATCH_SUPPORT_gInf_s
        );

        return t;
    }();
}

namespace Foam
{
    const ionicModelFamilyInfo& BuenoOrovioFamilyInfo()
    {
        static const ionicModelFamilyInfo info
        {
            NUM_CONSTANTS,
            NUM_STATES,
            NUM_ALGEBRAIC,
            BuenoOrovioCONSTANTS_NAMES,
            BuenoOrovioSTATES_NAMES,
            BuenoOrovioALGEBRAIC_NAMES,
            u,
            1000.0,
            1000.0/85.7,
            84.0/85.7,
            &buenoOrovioTransformedVm
        };
        return info;
    }
}

#include "addToRunTimeSelectionTable.H"
#include "ionicModelIO.H"
#include "stimulusIO.H"
#include "Pstream.H"

#ifdef HAS_CUDA
#include <cuda_runtime.h>

namespace Foam
{
    void launchBuenoBatchKernel
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

    void launchEulerStepKernel
    (
        double* d_STATES,
        const double* d_RATES,
        double dt,
        int N,
        int nStates,
        bool solveVm,
        int vmStateI
    );

    void launchScaleIonKernel
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
    bool useBuenoOrovioCompactSupport(const dictionary& dict)
    {
        const word modelName =
            dict.lookupOrDefault<word>("ionicModel", word::null);

        return modelName == "BuenoOroviocompactBatched"
            || dict.lookupOrDefault<Switch>("useCompactSupport", false);
    }

    defineTypeNameAndDebug(BuenoOrovioBatched, 0);
    addToRunTimeSelectionTable
    (
        ionicModel,
        BuenoOrovioBatched,
        dictionary
    );
    defineTypeNameAndDebug(BuenoOroviocompactBatched, 0);
    addToRunTimeSelectionTable
    (
        ionicModel,
        BuenoOroviocompactBatched,
        dictionary
    );
}

Foam::BuenoOrovioBatched::BuenoOrovioBatched
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
    ),
    useSoAEvaluator_
    (
        dict.lookupOrDefault<Switch>("useSoAEvaluator", false)
    ),
    useCompactSupport_(useBuenoOrovioCompactSupport(dict)),
#ifdef HAS_CUDA
    useDevice_(false),
#endif
    stimulusPOD_()
{
    ionicModel::setTissueFromDict();

#ifdef HAS_CUDA
    if (useSoAEvaluator_)
    {
        int nDevices = 0;
        cudaError_t err = cudaGetDeviceCount(&nDevices);
        if (err == cudaSuccess && nDevices > 0)
        {
            const int rank = Pstream::myProcNo();
            const int devId = rank % nDevices;
            CARDIAC_CUDA_CHECK(cudaSetDevice(devId));
            useDevice_ = true;
            Info<< "BuenoOrovioBatched: rank " << rank << " using CUDA device "
                << devId << " of " << nDevices << nl;
        }
        else
        {
            WarningInFunction
                << "BuenoOrovioBatched: useSoAEvaluator is on but "
                << "no CUDA device is visible (" << cudaGetErrorString(err)
                << "); falling back to the host SIMD path." << nl;
        }
    }
#endif

    if (useSoAEvaluator_ || useCompactSupport_)
    {
        setHotPathSupportSize(NUM_BO_BATCH_SUPPORT);
    }

    if (useSoAEvaluator_)
    {

        const word integrator =
            dict.lookupOrDefault<word>("batchedIntegrator", "euler");
        if (integrator != "euler")
        {
            WarningInFunction
                << "useSoAEvaluator is enabled for "
                << type() << " but `batchedIntegrator " << integrator
                << "` is not yet honored on the batched path; reverting "
                << "to explicit Euler for the SoA solver. Set "
                << "`batchedIntegrator euler;` to silence this warning, "
                << "or unset `useSoAEvaluator` to use the cell-major "
                << "integrator path." << nl;
        }
    }

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

Foam::BuenoOroviocompactBatched::BuenoOroviocompactBatched
(
    const dictionary& dict,
    const label num,
    const scalar initialDeltaT,
    const Switch solveVmWithinODESolver
)
:
    BuenoOrovioBatched(dict, num, initialDeltaT, solveVmWithinODESolver)
{}

Foam::BuenoOrovioBatched::~BuenoOrovioBatched()
{
#ifdef HAS_CUDA
    if (useDevice_)
    {
        cuda_.free();
    }
#endif
}




void Foam::BuenoOrovioBatched::prepareIOAccess
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
            static_cast<std::size_t>(NUM_BO_BATCH_SUPPORT),
            static_cast<std::size_t>(nCells())
        );
    }
#endif
    configuredBatchedIonicModel::prepareIOAccess(requestedNames, needsAlgebraics);
}


void Foam::BuenoOrovioBatched::importFields
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


void Foam::BuenoOrovioBatched::solveODE
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


void Foam::BuenoOrovioBatched::solveBatched
(
    const scalar stepStartTime,
    const scalar deltaT,
    const scalarField& Vm,
    scalarField& Im
)
{

    const label N = nCells();
    const label nSub = nSubsteps();
    const scalar dtModel = deltaT * timeScaleFactor();
    const scalar dtSubstep = dtModel/scalar(nSub);
    const scalar tStart = stepStartTime * timeScaleFactor();
    const bool   solveVm = solveVmWithinODESolver();

    stimulusPOD_ = stimulusIO::toPOD(stimulusProtocol());

    if (!solveVm)
    {
        for (label cellI = 0; cellI < N; ++cellI)
        {
            state(cellI, u) = vmToState(Vm[cellI]);
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
            static_cast<std::size_t>(NUM_BO_BATCH_SUPPORT),
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
            launchBuenoBatchKernel
            (
                tSub, cuda_.d_constants,
                static_cast<int>(N),
                cuda_.d_states, cuda_.d_rates, cuda_.d_support,
                tFlag, solveVm, stimulusPOD_
            );
            launchEulerStepKernel
            (
                cuda_.d_states, cuda_.d_rates,
                static_cast<double>(dtSubstep),
                static_cast<int>(N),
                static_cast<int>(NUM_STATES),
                solveVm,
                static_cast<int>(u)
            );
        }

        launchBuenoBatchKernel
        (
            tStart + dtModel, cuda_.d_constants,
            static_cast<int>(N),
            cuda_.d_states, cuda_.d_rates, cuda_.d_support,
            tFlag, solveVm, stimulusPOD_
        );

        launchScaleIonKernel
        (
            cuda_.d_support, cuda_.d_Im, 85.7,
            static_cast<int>(N),
            static_cast<int>(BO_BATCH_SUPPORT_Iion)
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

        BuenoOrovioComputeVariablesBatch
        (
            tSub,
            CONSTS,
            static_cast<int>(N),
            0, static_cast<int>(N),
            STATES_SoA,
            RATES_SoA,
            SUPPORT_SoA,
            solveVm,
            stimulusPOD_
        );

        for (label stateI = 0; stateI < NUM_STATES; ++stateI)
        {
            if (!solveVm && stateI == u)
            {
                continue;
            }

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

    const scalar tEnd = tStart + dtModel;
    BuenoOrovioComputeVariablesBatch
    (
        tEnd,
        CONSTS,
        static_cast<int>(N),
        0, static_cast<int>(N),
        STATES_SoA,
        RATES_SoA,
        SUPPORT_SoA,
        solveVm,
        stimulusPOD_
    );

    if (SUPPORT_SoA != nullptr)
    {
        const label IionBase = BO_BATCH_SUPPORT_Iion*N;
        for (label cellI = 0; cellI < N; ++cellI)
        {
            Im[cellI] = SUPPORT_SoA[IionBase + cellI]*85.7;
        }
    }
    else
    {
        for (label cellI = 0; cellI < N; ++cellI)
        {
            Im[cellI] = 0.0;
        }
    }

    setIOEvaluationModelTime(tEnd);
}

Foam::List<Foam::word> Foam::BuenoOrovioBatched::supportedTissueTypes() const
{
    return {"endocardialCells", "mCells", "epicardialCells"};
}

void Foam::BuenoOrovioBatched::evaluateState
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
        CONSTANTS_.data(),
        rateValues.data(),
        const_cast<scalarUList&>(stateValues).data(),
        algebraicValues.data(),
        tissue(),
        solveVmWithinODESolver(),
        stimulusProtocol()
    );
}

Foam::scalar Foam::BuenoOrovioBatched::ionicCurrentFromHotPathSupport
(
    const scalarUList& supportValues
) const
{
    return useCompactSupport_
      ? supportValues[BO_BATCH_SUPPORT_Iion]*85.7
      : ionicCurrentFromEvaluation(supportValues);
}

void Foam::BuenoOrovioBatched::evaluateHotPathState
(
    const scalar modelTime,
    const scalarUList& stateValues,
    scalarUList& rateValues,
    scalarUList& supportValues
) const
{
    if (!useCompactSupport_)
    {
        evaluateState(modelTime, stateValues, rateValues, supportValues);
        return;
    }

    using Foam::smoothHeaviside;

    scalarField algebraics(NUM_ALGEBRAIC, 0.0);
    evaluateState(modelTime, stateValues, rateValues, algebraics);

    const scalar cellV = stateValues[u];

    const scalar hV = smoothHeaviside(cellV - CONSTANTS_[thetaV]);
    const scalar invTauV =
        (1.0 - hV)/algebraics[tauVMinus]
      + hV/CONSTANTS_[tauVPlus];
    supportValues[BO_BATCH_SUPPORT_tau_v] = 1.0/invTauV;
    supportValues[BO_BATCH_SUPPORT_gInf_v] =
        (1.0 - hV)*algebraics[vInfty]/(algebraics[tauVMinus]*invTauV);

    const scalar hW = smoothHeaviside(cellV - CONSTANTS_[thetaW]);
    const scalar invTauW =
        (1.0 - hW)/algebraics[tauWMinus]
      + hW/CONSTANTS_[tauWPlus];
    supportValues[BO_BATCH_SUPPORT_tau_w] = 1.0/invTauW;
    supportValues[BO_BATCH_SUPPORT_gInf_w] =
        (1.0 - hW)*algebraics[wInfty]/(algebraics[tauWMinus]*invTauW);

    supportValues[BO_BATCH_SUPPORT_tau_s] = algebraics[tauS];
    supportValues[BO_BATCH_SUPPORT_gInf_s] =
        0.5*(1.0 + std::tanh(CONSTANTS_[kS]*(cellV - CONSTANTS_[uS])));
    supportValues[BO_BATCH_SUPPORT_Iion] = algebraics[Jion];
}

bool Foam::BuenoOrovioBatched::rushLarsenParameters
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
            if (invTau <= VSMALL) return false;
            tau = 1.0/invTau;

            steadyState =
                (1.0 - hV) * algebraicValues[vInfty]
              / (algebraicValues[tauVMinus] * invTau);

            return std::isfinite(steadyState);
        }

        case w:
        {
            const scalar hW = smoothHeaviside(cellV - CONSTANTS_[thetaW]);
            const scalar invTau =
                (1.0 - hW)/algebraicValues[tauWMinus]
              + hW/CONSTANTS_[tauWPlus];
            if (invTau <= VSMALL) return false;
            tau = 1.0/invTau;

            steadyState =
                (1.0 - hW) * algebraicValues[wInfty]
              / (algebraicValues[tauWMinus] * invTau);

            return std::isfinite(steadyState);
        }

        case s:
        {
            tau = algebraicValues[tauS];
            if (tau <= VSMALL) return false;

            steadyState = 0.5*(1.0 + std::tanh(
                CONSTANTS_[kS]*(cellV - CONSTANTS_[uS])
            ));
            return std::isfinite(steadyState);
        }

        default:
            return false;
    }
}

bool Foam::BuenoOrovioBatched::rushLarsenParametersFromHotPathSupport
(
    const label stateI,
    const scalarUList& stateValues,
    const scalarUList& rateValues,
    const scalarUList& supportValues,
    scalar& steadyState,
    scalar& tau
) const
{
    if (!useCompactSupport_)
    {
        return rushLarsenParameters
        (
            stateI,
            stateValues,
            rateValues,
            supportValues,
            steadyState,
            tau
        );
    }

    if (stateI < 0 || stateI >= NUM_STATES) return false;

    return resolveSupportRushLarsenEntry
    (
        BuenoOrovioRushLarsenDispatch[stateI],
        CONSTANTS_,
        supportValues,
        VSMALL,
        steadyState,
        tau
    );
}

void Foam::BuenoOrovioBatched::derivatives
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
