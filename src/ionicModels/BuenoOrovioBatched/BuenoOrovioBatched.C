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

    for (label sub = 0; sub < nSub; ++sub)
    {
        const scalar tSub = tStart + scalar(sub)*dtSubstep;
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
    return supportValues[BO_BATCH_SUPPORT_Iion];
}

void Foam::BuenoOrovioBatched::evaluateHotPathState
(
    const scalar modelTime,
    const scalarUList& stateValues,
    scalarUList& rateValues,
    scalarUList& supportValues
) const
{
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


void Foam::BuenoOrovioBatched::evaluateHotPathState
(
    const label cellI,
    const scalar modelTime,
    const scalarUList& stateValues,
    scalarUList& rateValues,
    scalarUList& supportValues
) const
{
    using Foam::smoothHeaviside;

    scalarField& cellConstants = constants(cellI);
    scalarField algebraics(NUM_ALGEBRAIC, 0.0);
    evaluateState(cellI, modelTime, stateValues, rateValues, algebraics);

    const scalar cellV = stateValues[u];

    const scalar hV = smoothHeaviside(cellV - cellConstants[thetaV]);
    const scalar invTauV =
        (1.0 - hV)/algebraics[tauVMinus]
      + hV/cellConstants[tauVPlus];
    supportValues[BO_BATCH_SUPPORT_tau_v] = 1.0/invTauV;
    supportValues[BO_BATCH_SUPPORT_gInf_v] =
        (1.0 - hV)*algebraics[vInfty]/(algebraics[tauVMinus]*invTauV);

    const scalar hW = smoothHeaviside(cellV - cellConstants[thetaW]);
    const scalar invTauW =
        (1.0 - hW)/algebraics[tauWMinus]
      + hW/cellConstants[tauWPlus];
    supportValues[BO_BATCH_SUPPORT_tau_w] = 1.0/invTauW;
    supportValues[BO_BATCH_SUPPORT_gInf_w] =
        (1.0 - hW)*algebraics[wInfty]/(algebraics[tauWMinus]*invTauW);

    supportValues[BO_BATCH_SUPPORT_tau_s] = algebraics[tauS];
    supportValues[BO_BATCH_SUPPORT_gInf_s] =
        0.5*(1.0 + std::tanh(cellConstants[kS]*(cellV - cellConstants[uS])));
    supportValues[BO_BATCH_SUPPORT_Iion] = algebraics[Jion];
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
