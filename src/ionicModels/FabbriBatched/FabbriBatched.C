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

#include "FabbriBatched.H"
#include "Fabbri_2017Batch.H"
#include "batchedRushLarsenEntry.H"
#include <array>
#include <cmath>

namespace Foam
{
    const ionicModelFamilyInfo& FabbriFamilyInfo();
}

namespace
{
    const std::array<Foam::batchedRushLarsenEntry, NUM_STATES>
    FabbriRushLarsenDispatch = []()
    {
        std::array<Foam::batchedRushLarsenEntry, NUM_STATES> t{};
        t.fill(Foam::rlNone());

        t[If_y_gate_y]           = Foam::rlScalarAlgAndSupport(AV_tau_y,     AV_y_infinity,     Foam::FABBRI_BATCH_SUPPORT_tau_y,     Foam::FABBRI_BATCH_SUPPORT_gInf_y);
        t[INa_m_gate_m]          = Foam::rlScalarAlgAndSupport(AV_tau_m,     AV_m_infinity,     Foam::FABBRI_BATCH_SUPPORT_tau_m,     Foam::FABBRI_BATCH_SUPPORT_gInf_m);
        t[INa_h_gate_h]          = Foam::rlScalarAlgAndSupport(AV_tau_h,     AV_h_infinity,     Foam::FABBRI_BATCH_SUPPORT_tau_h,     Foam::FABBRI_BATCH_SUPPORT_gInf_h);
        t[ICaL_dL_gate_dL]       = Foam::rlScalarAlgAndSupport(AV_tau_dL,    AV_dL_infinity,    Foam::FABBRI_BATCH_SUPPORT_tau_dL,    Foam::FABBRI_BATCH_SUPPORT_gInf_dL);
        t[ICaL_fL_gate_fL]       = Foam::rlScalarAlgAndSupport(AV_tau_fL,    AV_fL_infinity,    Foam::FABBRI_BATCH_SUPPORT_tau_fL,    Foam::FABBRI_BATCH_SUPPORT_gInf_fL);
        t[ICaL_fCa_gate_fCa]     = Foam::rlScalarAlgAndSupport(AV_tau_fCa,   AV_fCa_infinity,   Foam::FABBRI_BATCH_SUPPORT_tau_fCa,   Foam::FABBRI_BATCH_SUPPORT_gInf_fCa);
        t[ICaT_dT_gate_dT]       = Foam::rlScalarAlgAndSupport(AV_tau_dT,    AV_dT_infinity,    Foam::FABBRI_BATCH_SUPPORT_tau_dT,    Foam::FABBRI_BATCH_SUPPORT_gInf_dT);
        t[ICaT_fT_gate_fT]       = Foam::rlScalarAlgAndSupport(AV_tau_fT,    AV_fT_infinity,    Foam::FABBRI_BATCH_SUPPORT_tau_fT,    Foam::FABBRI_BATCH_SUPPORT_gInf_fT);
        t[IKur_rKur_gate_r_Kur]  = Foam::rlScalarAlgAndSupport(AV_tau_r_Kur, AV_r_Kur_infinity, Foam::FABBRI_BATCH_SUPPORT_tau_r_Kur, Foam::FABBRI_BATCH_SUPPORT_gInf_r_Kur);
        t[IKur_sKur_gate_s_Kur]  = Foam::rlScalarAlgAndSupport(AV_tau_s_Kur, AV_s_Kur_infinity, Foam::FABBRI_BATCH_SUPPORT_tau_s_Kur, Foam::FABBRI_BATCH_SUPPORT_gInf_s_Kur);
        t[Ito_q_gate_q]          = Foam::rlScalarAlgAndSupport(AV_tau_q,     AV_q_infinity,     Foam::FABBRI_BATCH_SUPPORT_tau_q,     Foam::FABBRI_BATCH_SUPPORT_gInf_q);
        t[Ito_r_gate_r]          = Foam::rlScalarAlgAndSupport(AV_tau_r,     AV_r_infinity,     Foam::FABBRI_BATCH_SUPPORT_tau_r,     Foam::FABBRI_BATCH_SUPPORT_gInf_r);
        t[IKr_pa_gate_paS]       = Foam::rlScalarAlgAndSupport(AV_tau_paS,   AV_pa_infinity,    Foam::FABBRI_BATCH_SUPPORT_tau_paS,   Foam::FABBRI_BATCH_SUPPORT_gInf_paS);
        t[IKr_pa_gate_paF]       = Foam::rlScalarAlgAndSupport(AV_tau_paF,   AV_pa_infinity,    Foam::FABBRI_BATCH_SUPPORT_tau_paF,   Foam::FABBRI_BATCH_SUPPORT_gInf_paF);
        t[IKr_pi_gate_piy]       = Foam::rlScalarAlgAndSupport(AV_tau_pi,    AV_pi_infinity,    Foam::FABBRI_BATCH_SUPPORT_tau_pi,    Foam::FABBRI_BATCH_SUPPORT_gInf_pi);
        t[IKs_n_gate_n]          = Foam::rlScalarAlgAndSupport(AV_tau_n,     AV_n_infinity,     Foam::FABBRI_BATCH_SUPPORT_tau_n,     Foam::FABBRI_BATCH_SUPPORT_gInf_n);
        t[IKACh_a_gate_a]        = Foam::rlScalarAlgAndSupport(AV_tau_a,     AV_a_infinity,     Foam::FABBRI_BATCH_SUPPORT_tau_a,     Foam::FABBRI_BATCH_SUPPORT_gInf_a);

        return t;
    }();
}

#include "addToRunTimeSelectionTable.H"
#include "Pstream.H"
#include "stimulusIO.H"

#ifdef HAS_CUDA
#include <cuda_runtime.h>

namespace Foam
{
    void launchFabbriBatchKernel
    (
        double t,
        const double* d_CONSTANTS,
        int N,
        const double* d_STATES,
        double* d_RATES,
        double* d_SUPPORT,
        bool solveVm,
        StimulusProtocolPOD stimulus
    );

    void launchFabbriEulerStepKernel
    (
        double* d_STATES,
        const double* d_RATES,
        double dt,
        int N,
        int nStates,
        bool solveVm,
        int vmStateI
    );

    void launchFabbriRushLarsenStepKernel
    (
        double* d_STATES,
        const double* d_RATES,
        const double* d_SUPPORT,
        double dt,
        int N,
        int nStates,
        bool solveVm,
        int vmStateI
    );

    void launchFabbriScaleIonKernel
    (
        const double* d_SUPPORT,
        double* d_Im,
        double scale,
        int N,
        int IionSlot
    );
}

#endif // HAS_CUDA

namespace Foam
{
    defineTypeNameAndDebug(FabbriBatched, 0);
    defineTypeNameAndDebug(FabbricompactBatched, 0);
    addToRunTimeSelectionTable
    (
        ionicModel,
        FabbricompactBatched,
        dictionary
    );
}

Foam::FabbriBatched::FabbriBatched
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
        FabbriFamilyInfo()
    ),
#ifdef HAS_CUDA
    useDevice_(false),
#endif
    stimulusPOD_()
{
    ionicModel::setTissueFromDict();

    setHotPathSupportSize(NUM_FABBRI_BATCH_SUPPORT);

#ifdef HAS_CUDA
    {
        int nDevices = 0;
        cudaError_t err = cudaGetDeviceCount(&nDevices);
        if (err == cudaSuccess && nDevices > 0)
        {
            const int rank = Pstream::myProcNo();
            const int devId = rank % nDevices;
            CARDIAC_CUDA_CHECK(cudaSetDevice(devId));
            useDevice_ = true;
            Info<< "FabbriBatched: rank " << rank
                << " using CUDA device " << devId
                << " of " << nDevices << nl;
        }
    }
#endif

    double initialRates[NUM_STATES] = {0.0};
    double initialStates[NUM_STATES] = {0.0};

    FabbriinitConsts
    (
        CONSTANTS_.data(),
        initialRates,
        initialStates,
        tissue(),
        dict
    );

    applyIonicConstantOverrides();

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

Foam::FabbricompactBatched::FabbricompactBatched
(
    const dictionary& dict,
    const label num,
    const scalar initialDeltaT,
    const Switch solveVmWithinODESolver
)
:
    FabbriBatched(dict, num, initialDeltaT, solveVmWithinODESolver)
{}

Foam::FabbriBatched::~FabbriBatched()
{
#ifdef HAS_CUDA
    if (useDevice_)
    {
        cuda_.free();
    }
#endif
}


Foam::List<Foam::word> Foam::FabbriBatched::supportedTissueTypes() const
{
    return {"myocyte"};
}

void Foam::FabbriBatched::evaluateState
(
    const scalar modelTime,
    const scalarUList& stateValues,
    scalarUList& rateValues,
    scalarUList& algebraicValues
) const
{
    FabbricomputeVariables
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

Foam::scalar Foam::FabbriBatched::ionicCurrentFromHotPathSupport
(
    const scalarUList& supportValues
) const
{
    return supportValues[FABBRI_BATCH_SUPPORT_Iion_cm];
}

void Foam::FabbriBatched::evaluateHotPathState
(
    const scalar modelTime,
    const scalarUList& stateValues,
    scalarUList& rateValues,
    scalarUList& supportValues
) const
{
    scalarField algebraics(NUM_ALGEBRAIC, 0.0);
    evaluateState(modelTime, stateValues, rateValues, algebraics);

    for (const auto& entry : FabbriRushLarsenDispatch)
    {
        projectScalarRushLarsenEntryToSupport
        (
            entry,
            CONSTANTS_,
            algebraics,
            supportValues
        );
    }

    supportValues[FABBRI_BATCH_SUPPORT_Iion_cm] = algebraics[Iion_cm];
}

bool Foam::FabbriBatched::rushLarsenParametersFromHotPathSupport
(
    const label stateI,
    const scalarUList& stateValues,
    const scalarUList& rateValues,
    const scalarUList& supportValues,
    scalar& steadyState,
    scalar& tau
) const
{
    if (stateI < 0 || stateI >= NUM_STATES) return false;

    const auto& entry = FabbriRushLarsenDispatch[stateI];
    return resolveSupportRushLarsenEntry
    (
        entry,
        CONSTANTS_,
        supportValues,
        VSMALL,
        steadyState,
        tau
    );
}

void Foam::FabbriBatched::prepareIOAccess
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
            static_cast<std::size_t>(NUM_FABBRI_BATCH_SUPPORT),
            static_cast<std::size_t>(nCells())
        );
    }
#endif

    configuredBatchedIonicModel::prepareIOAccess
    (
        requestedNames,
        needsAlgebraics
    );
}


void Foam::FabbriBatched::importFields
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


void Foam::FabbriBatched::solveODE
(
    const scalar stepStartTime,
    const scalar deltaT,
    const scalarField& Vm,
    scalarField& Im
)
{
    solveODEImpl(*this, stepStartTime, deltaT, Vm, Im);
}

void Foam::FabbricompactBatched::solveODE
(
    const scalar stepStartTime,
    const scalar deltaT,
    const scalarField& Vm,
    scalarField& Im
)
{
#ifdef HAS_CUDA
    if (useDevice_ && !utilitiesMode())
    {
        solveOnDevice(stepStartTime, deltaT, Vm, Im);
        return;
    }
#endif
    solveODEImpl(*this, stepStartTime, deltaT, Vm, Im);
}

#ifdef HAS_CUDA
void Foam::FabbricompactBatched::solveOnDevice
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
    const bool solveVm = solveVmWithinODESolver();
    scalarField vmStateSlice;

    stimulusPOD_ = stimulusIO::toPOD(stimulusProtocol());

    if (!solveVm)
    {
        vmStateSlice.setSize(N);
        for (label cellI = 0; cellI < N; ++cellI)
        {
            const scalar vmState = vmToState(Vm[cellI]);
            state(cellI, membrane_V) = vmState;
            vmStateSlice[cellI] = vmState;
        }
    }

    markIODirty();
    setIOEvaluationModelTime(tStart);

    cuda_.allocate
    (
        static_cast<std::size_t>(N),
        static_cast<std::size_t>(NUM_STATES),
        static_cast<std::size_t>(nHotPathSupport()),
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
    else if (!solveVm)
    {
        cuda_.syncStateSliceHostToDevice
        (
            vmStateSlice.cdata(),
            static_cast<std::size_t>(membrane_V),
            static_cast<std::size_t>(N)
        );
    }

    for (label sub = 0; sub < nSub; ++sub)
    {
        const scalar tSub = tStart + scalar(sub)*dtSubstep;
        launchFabbriBatchKernel
        (
            tSub, cuda_.d_constants, static_cast<int>(N),
            cuda_.d_states, cuda_.d_rates, cuda_.d_support,
            solveVm, stimulusPOD_
        );
        launchFabbriRushLarsenStepKernel
        (
            cuda_.d_states, cuda_.d_rates, cuda_.d_support,
            static_cast<double>(dtSubstep),
            static_cast<int>(N),
            static_cast<int>(NUM_STATES),
            solveVm,
            static_cast<int>(membrane_V)
        );
    }

    launchFabbriBatchKernel
    (
        tStart + dtModel, cuda_.d_constants, static_cast<int>(N),
        cuda_.d_states, cuda_.d_rates, cuda_.d_support,
        solveVm, stimulusPOD_
    );
    launchFabbriScaleIonKernel
    (
        cuda_.d_support, cuda_.d_Im, 1.0,
        static_cast<int>(N),
        static_cast<int>(FABBRI_BATCH_SUPPORT_Iion_cm)
    );
    cuda_.downloadIm(Im.data(), static_cast<std::size_t>(N));
    cuda_.deviceDirty = true;
    setIOEvaluationModelTime(tStart + dtModel);
}
#endif // HAS_CUDA


void Foam::FabbriBatched::derivatives
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
