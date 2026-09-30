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
#include "BuenoOrovio_2008Batch.H"
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
        const double* d_CELL_CONSTANTS,
        bool useCellConstants,
        int N,
        const double* d_STATES,
        double* d_RATES,
        double* d_SUPPORT,
        int tissueFlag,
        bool solveVm,
        StimulusProtocolPOD stimulus
    );

    void launchBuenoEulerStepKernel
    (
        double* d_STATES,
        const double* d_RATES,
        double dt,
        int N,
        int nStates,
        bool solveVm,
        int vmStateI
    );

    void launchBuenoRushLarsenStepKernel
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

    void launchBuenoScaleIonKernel
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
    defineTypeNameAndDebug(BuenoOrovioBatched, 0);
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
#ifdef HAS_CUDA
    useDevice_(false),
#endif
    stimulusPOD_()
{
    ionicModel::setTissueFromDict();

    setHotPathSupportSize(NUM_BO_BATCH_SUPPORT);

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
            Info<< "BuenoOrovioBatched: rank " << rank << " using CUDA device "
                << devId << " of " << nDevices << nl;
        }
    }
#endif

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

void Foam::BuenoOroviocompactBatched::solveODE
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
    stimulusPOD_ = stimulusIO::toPOD(stimulusProtocol());
    solveODEImpl(*this, stepStartTime, deltaT, Vm, Im);
}

#ifdef HAS_CUDA
void Foam::BuenoOroviocompactBatched::solveOnDevice
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
    const bool useVmExtrapolant = hasVmRate() && !solveVm;
    scalarField vmStateStart;
    scalarField vmStateRate;

    if (useVmExtrapolant)
    {
        vmStateStart.setSize(N);
        vmStateRate.setSize(N);
        for (label cellI = 0; cellI < N; ++cellI)
        {
            const scalar vmStart = vmToState(Vm[cellI]);
            const scalar vmEnd =
                vmToState(Vm[cellI] + VmRateSI(cellI)*deltaT);
            vmStateStart[cellI] = vmStart;
            vmStateRate[cellI] =
                (mag(dtModel) > VSMALL)
              ? (vmEnd - vmStart)/dtModel
              : 0.0;
        }
    }
    const int tFlag = static_cast<int>(tissue());

    stimulusPOD_ = stimulusIO::toPOD(stimulusProtocol());

    if (!solveVm)
    {
        for (label cellI = 0; cellI < N; ++cellI)
        {
            state(cellI, u) = vmToState(Vm[cellI]);
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
    if (useVmExtrapolant)
    {
        cuda_.uploadVmExtrapolant
        (
            vmStateStart.cdata(),
            vmStateRate.cdata(),
            static_cast<std::size_t>(N)
        );
    }
    scalarField flattenedCellConstants;
    if (hasHeterogeneousConstants())
    {
        flattenHeterogeneousConstants(flattenedCellConstants);
        cuda_.uploadCellConstants
        (
            flattenedCellConstants.cdata(),
            static_cast<std::size_t>(N),
            static_cast<std::size_t>(CONSTANTS_.size())
        );
    }

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
        const std::size_t vmStateI =
            static_cast<std::size_t>(voltageStateIndex());
        cuda_.syncStateSliceHostToDevice
        (
            statesSoAData() + vmStateI*static_cast<std::size_t>(N),
            vmStateI,
            static_cast<std::size_t>(N)
        );
    }

    for (label sub = 0; sub < nSub; ++sub)
    {
        const scalar tSub = tStart + scalar(sub)*dtSubstep;
        if (useVmExtrapolant)
        {
            cuda_.applyVmExtrapolant
            (
                static_cast<double>(tSub),
                static_cast<double>(tStart),
                static_cast<std::size_t>(voltageStateIndex()),
                static_cast<std::size_t>(N)
            );
        }
        launchBuenoBatchKernel
        (
            tSub,
            cuda_.d_constants,
            cuda_.d_cellConstants,
            hasHeterogeneousConstants(),
            static_cast<int>(N),
            cuda_.d_states, cuda_.d_rates, cuda_.d_support,
            tFlag, solveVm, stimulusPOD_
        );
        if (useEulerIntegrator())
        {
            launchBuenoEulerStepKernel
            (
                cuda_.d_states, cuda_.d_rates,
                static_cast<double>(dtSubstep),
                static_cast<int>(N),
                static_cast<int>(NUM_STATES),
                solveVm,
                static_cast<int>(u)
            );
        }
        else
        {
            launchBuenoRushLarsenStepKernel
            (
                cuda_.d_states, cuda_.d_rates, cuda_.d_support,
                static_cast<double>(dtSubstep),
                static_cast<int>(N),
                static_cast<int>(NUM_STATES),
                solveVm,
                static_cast<int>(u)
            );
        }
    }

    if (useVmExtrapolant)
    {
        cuda_.applyVmExtrapolant
        (
            static_cast<double>(tStart + dtModel),
            static_cast<double>(tStart),
            static_cast<std::size_t>(voltageStateIndex()),
            static_cast<std::size_t>(N)
        );
    }
    launchBuenoBatchKernel
    (
        tStart + dtModel,
        cuda_.d_constants,
        cuda_.d_cellConstants,
        hasHeterogeneousConstants(),
        static_cast<int>(N),
        cuda_.d_states, cuda_.d_rates, cuda_.d_support,
        tFlag, solveVm, stimulusPOD_
    );
    launchBuenoScaleIonKernel
    (
        cuda_.d_support, cuda_.d_Im, 85.7,
        static_cast<int>(N),
        static_cast<int>(BO_BATCH_SUPPORT_Iion)
    );
    cuda_.downloadIm(Im.data(), static_cast<std::size_t>(N));
    cuda_.deviceDirty = true;
    setIOEvaluationModelTime(tStart + dtModel);
}
#endif // HAS_CUDA

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
#ifdef HAS_CUDA
    if (useDevice_ && cuda_.allocated && gpuSelectionNeedsRates(requestedNames))
    {
        cuda_.syncRatesDeviceToHost
        (
            ratesSoAData(),
            static_cast<std::size_t>(NUM_STATES),
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
    stimulusPOD_ = stimulusIO::toPOD(stimulusProtocol());
    solveODEImpl(*this, stepStartTime, deltaT, Vm, Im);
}

Foam::List<Foam::word> Foam::BuenoOrovioBatched::supportedTissueTypes() const
{
    return {"endocardialCells", "mCells", "epicardialCells"};
}


Foam::scalarField Foam::BuenoOrovioBatched::constantsForTissue
(
    const label tissueFlag
) const
{
    scalarField constants(NUM_CONSTANTS, 0.0);
    scalarField rates(NUM_STATES, 0.0);
    scalarField states(NUM_STATES, 0.0);

    BuenoOrovioinitConsts
    (
        constants.data(),
        rates.data(),
        states.data(),
        tissueFlag,
        dict()
    );

    ionicModelIO::applyConstantOverrides
    (
        constants,
        BuenoOrovioCONSTANTS_NAMES,
        NUM_CONSTANTS,
        dict(),
        type(),
        tissueFlag
    );

    return constants;
}

Foam::scalarField Foam::BuenoOrovioBatched::initialStatesForTissue
(
    const label tissueFlag
) const
{
    scalarField constants(NUM_CONSTANTS, 0.0);
    scalarField rates(NUM_STATES, 0.0);
    scalarField states(NUM_STATES, 0.0);

    BuenoOrovioinitConsts
    (
        constants.data(),
        rates.data(),
        states.data(),
        tissueFlag,
        dict()
    );

    return states;
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


void Foam::BuenoOrovioBatched::evaluateState
(
    const label cellI,
    const scalar modelTime,
    const scalarUList& stateValues,
    scalarUList& rateValues,
    scalarUList& algebraicValues
) const
{
    scalarField& cellConstants = constants(cellI);

    BuenoOroviocomputeVariables
    (
        modelTime,
        cellConstants.data(),
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
    // The scalar model converts reduced Jion to tissue current as 85.7*Jion.
    return 85.7*supportValues[BO_BATCH_SUPPORT_Iion];
}

void Foam::BuenoOrovioBatched::evaluateHotPathState
(
    const scalar modelTime,
    const scalarUList& stateValues,
    scalarUList& rateValues,
    scalarUList& supportValues
) const
{
    BuenoOrovioComputeVariablesBatch
    (
        modelTime,
        CONSTANTS_.data(),
        1,
        0,
        1,
        const_cast<scalarUList&>(stateValues).data(),
        rateValues.data(),
        supportValues.data(),
        solveVmWithinODESolver(),
        stimulusPOD_
    );
}


void Foam::BuenoOrovioBatched::evaluateHotPathState
(
    const label cellI,
    const scalar modelTime,
    const scalarUList& stateValues,
    scalarUList& rateValues,
    scalarUList& supportValues
) const
{
    scalarField& cellConstants = constants(cellI);
    BuenoOrovioComputeVariablesBatch
    (
        modelTime,
        cellConstants.data(),
        1,
        0,
        1,
        const_cast<scalarUList&>(stateValues).data(),
        rateValues.data(),
        supportValues.data(),
        solveVmWithinODESolver(),
        stimulusPOD_
    );
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


bool Foam::BuenoOrovioBatched::rushLarsenParametersFromHotPathSupport
(
    const label cellI,
    const label stateI,
    const scalarUList& stateValues,
    const scalarUList& rateValues,
    const scalarUList& supportValues,
    scalar& steadyState,
    scalar& tau
) const
{
    if (stateI < 0 || stateI >= NUM_STATES) return false;

    return resolveSupportRushLarsenEntry
    (
        BuenoOrovioRushLarsenDispatch[stateI],
        constants(cellI),
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
