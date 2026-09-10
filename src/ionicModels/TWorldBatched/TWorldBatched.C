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

#include "TWorldBatched.H"
#include "TWorld_2025Batch.H"
#include "TWorld_2025.H"
#include "batchedRushLarsenEntry.H"
#include <array>

namespace Foam
{
    const ionicModelFamilyInfo& TWorldFamilyInfo()
    {
        static const ionicModelFamilyInfo info
        {
            NUM_CONSTANTS,
            NUM_STATES,
            NUM_ALGEBRAIC,
            TWorldCONSTANTS_NAMES,
            TWorldSTATES_NAMES,
            TWorldALGEBRAIC_NAMES,
            v,
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

namespace
{
    const std::array<Foam::batchedRushLarsenEntry, NUM_STATES>
    TWorldRushLarsenDispatch = []()
    {
        std::array<Foam::batchedRushLarsenEntry, NUM_STATES> t{};
        t.fill(Foam::rlNone());

        // CaMK fractions (constant tau)
        t[camk_f_ICaL] = Foam::rlScalarConstTauAlgInfAndSupport(AC_tau_cal, AV_CaMK_Phos_ss_ICaL, Foam::TWORLD_BATCH_SUPPORT_tau_camk_ICaL, Foam::TWORLD_BATCH_SUPPORT_gInf_camk_ICaL);
        t[camk_f_PLB]  = Foam::rlScalarConstTauAlgInfAndSupport(AC_tau_plb, AV_CaMK_Phos_ss_PLB,   Foam::TWORLD_BATCH_SUPPORT_tau_camk_PLB,  Foam::TWORLD_BATCH_SUPPORT_gInf_camk_PLB);
        t[camk_f_RyR]  = Foam::rlScalarConstTauAlgInfAndSupport(AC_tau_ryr, AV_CaMK_Phos_ss_RyR,   Foam::TWORLD_BATCH_SUPPORT_tau_camk_RyR,  Foam::TWORLD_BATCH_SUPPORT_gInf_camk_RyR);

        // INa
        t[m]    = Foam::rlScalarAlgAndSupport(AV_tm,  AV_mss,    Foam::TWORLD_BATCH_SUPPORT_tau_m,    Foam::TWORLD_BATCH_SUPPORT_gInf_m);
        t[m_P]  = Foam::rlScalarAlgAndSupport(AV_tm,  AV_mss_P,  Foam::TWORLD_BATCH_SUPPORT_tau_m_P,  Foam::TWORLD_BATCH_SUPPORT_gInf_m_P);
        t[h]    = Foam::rlScalarAlgAndSupport(AV_th,  AV_hss,    Foam::TWORLD_BATCH_SUPPORT_tau_h,    Foam::TWORLD_BATCH_SUPPORT_gInf_h);
        t[h_P]  = Foam::rlScalarAlgAndSupport(AV_th,  AV_hss_P,  Foam::TWORLD_BATCH_SUPPORT_tau_h_P,  Foam::TWORLD_BATCH_SUPPORT_gInf_h_P);
        t[hp]   = Foam::rlScalarAlgAndSupport(AV_th,  AV_hssp,   Foam::TWORLD_BATCH_SUPPORT_tau_hp,   Foam::TWORLD_BATCH_SUPPORT_gInf_hp);
        t[hp_P] = Foam::rlScalarAlgAndSupport(AV_th,  AV_hssp_P, Foam::TWORLD_BATCH_SUPPORT_tau_hp_P, Foam::TWORLD_BATCH_SUPPORT_gInf_hp_P);
        t[j]    = Foam::rlScalarAlgAndSupport(AV_tj,  AV_jss,    Foam::TWORLD_BATCH_SUPPORT_tau_j,    Foam::TWORLD_BATCH_SUPPORT_gInf_j);
        t[j_P]  = Foam::rlScalarAlgAndSupport(AV_tj,  AV_jss_P,  Foam::TWORLD_BATCH_SUPPORT_tau_j_P,  Foam::TWORLD_BATCH_SUPPORT_gInf_j_P);
        t[jp]   = Foam::rlScalarAlgAndSupport(AV_tjp, AV_jss,    Foam::TWORLD_BATCH_SUPPORT_tau_jp,   Foam::TWORLD_BATCH_SUPPORT_gInf_jp);
        t[jp_P] = Foam::rlScalarAlgAndSupport(AV_tjp, AV_jssp_P, Foam::TWORLD_BATCH_SUPPORT_tau_jp_P, Foam::TWORLD_BATCH_SUPPORT_gInf_jp_P);

        // INaL
        t[mL]  = Foam::rlScalarAlgAndSupport(AV_tmL, AV_mLss,  Foam::TWORLD_BATCH_SUPPORT_tau_mL,  Foam::TWORLD_BATCH_SUPPORT_gInf_mL);
        t[hL]  = Foam::rlScalarConstTauAlgInfAndSupport(AC_thL,  AV_hLss,  Foam::TWORLD_BATCH_SUPPORT_tau_hL,  Foam::TWORLD_BATCH_SUPPORT_gInf_hL);
        t[hLp] = Foam::rlScalarConstTauAlgInfAndSupport(AC_thLp, AV_hLssp, Foam::TWORLD_BATCH_SUPPORT_tau_hLp, Foam::TWORLD_BATCH_SUPPORT_gInf_hLp);

        // ICaL standard
        t[d]     = Foam::rlScalarAlgAndSupport(AV_td,     AV_dss,   Foam::TWORLD_BATCH_SUPPORT_tau_d,     Foam::TWORLD_BATCH_SUPPORT_gInf_d);
        t[ff]    = Foam::rlScalarAlgAndSupport(AV_tff,    AV_fss,   Foam::TWORLD_BATCH_SUPPORT_tau_ff,    Foam::TWORLD_BATCH_SUPPORT_gInf_ff);
        t[fs]    = Foam::rlScalarAlgAndSupport(AV_tfs,    AV_fss,   Foam::TWORLD_BATCH_SUPPORT_tau_fs,    Foam::TWORLD_BATCH_SUPPORT_gInf_fs);
        t[fcaf]  = Foam::rlScalarAlgAndSupport(AV_tfcaf,  AV_fcass, Foam::TWORLD_BATCH_SUPPORT_tau_fcaf,  Foam::TWORLD_BATCH_SUPPORT_gInf_fcaf);
        t[fcas]  = Foam::rlScalarAlgAndSupport(AV_tfcas,  AV_fcass, Foam::TWORLD_BATCH_SUPPORT_tau_fcas,  Foam::TWORLD_BATCH_SUPPORT_gInf_fcas);
        t[jca]   = Foam::rlScalarConstTauAlgInfAndSupport(AC_tjca,  AV_jcass, Foam::TWORLD_BATCH_SUPPORT_tau_jca,   Foam::TWORLD_BATCH_SUPPORT_gInf_jca);
        t[ffp]   = Foam::rlScalarAlgAndSupport(AV_tffp,   AV_fss,   Foam::TWORLD_BATCH_SUPPORT_tau_ffp,   Foam::TWORLD_BATCH_SUPPORT_gInf_ffp);
        t[fcafp] = Foam::rlScalarAlgAndSupport(AV_tfcafp, AV_fcass, Foam::TWORLD_BATCH_SUPPORT_tau_fcafp, Foam::TWORLD_BATCH_SUPPORT_gInf_fcafp);

        // ICaL PKA/BP variants
        t[d_P]    = Foam::rlScalarAlgAndSupport(AV_td,    AV_dPss,    Foam::TWORLD_BATCH_SUPPORT_tau_d_P,    Foam::TWORLD_BATCH_SUPPORT_gInf_d_P);
        t[ff_P]   = Foam::rlScalarAlgAndSupport(AV_tff,   AV_fss_P,   Foam::TWORLD_BATCH_SUPPORT_tau_ff_P,   Foam::TWORLD_BATCH_SUPPORT_gInf_ff_P);
        t[fs_P]   = Foam::rlScalarAlgAndSupport(AV_tfs,   AV_fss_P,   Foam::TWORLD_BATCH_SUPPORT_tau_fs_P,   Foam::TWORLD_BATCH_SUPPORT_gInf_fs_P);
        t[fcaf_P] = Foam::rlScalarAlgAndSupport(AV_tfcaf, AV_fcass_P, Foam::TWORLD_BATCH_SUPPORT_tau_fcaf_P, Foam::TWORLD_BATCH_SUPPORT_gInf_fcaf_P);
        t[fcas_P] = Foam::rlScalarAlgAndSupport(AV_tfcas, AV_fcass_P, Foam::TWORLD_BATCH_SUPPORT_tau_fcas_P, Foam::TWORLD_BATCH_SUPPORT_gInf_fcas_P);
        t[fBPf]   = Foam::rlScalarAlgAndSupport(AV_tffp,  AV_fBPss,   Foam::TWORLD_BATCH_SUPPORT_tau_fBPf,   Foam::TWORLD_BATCH_SUPPORT_gInf_fBPf);
        t[fcaBPf] = Foam::rlScalarAlgAndSupport(AV_tfcafp,AV_fcaBPss, Foam::TWORLD_BATCH_SUPPORT_tau_fcaBPf, Foam::TWORLD_BATCH_SUPPORT_gInf_fcaBPf);

        // Ito
        t[xtos]   = Foam::rlScalarAlgAndSupport(AV_tauxtos,   AV_xtoss,   Foam::TWORLD_BATCH_SUPPORT_tau_xtos,   Foam::TWORLD_BATCH_SUPPORT_gInf_xtos);
        t[xtos_p] = Foam::rlScalarAlgAndSupport(AV_tauxtos,   AV_xtoss_p, Foam::TWORLD_BATCH_SUPPORT_tau_xtos_p, Foam::TWORLD_BATCH_SUPPORT_gInf_xtos_p);
        t[xtof]   = Foam::rlScalarAlgAndSupport(AV_tauxtof,   AV_xtoss,   Foam::TWORLD_BATCH_SUPPORT_tau_xtof,   Foam::TWORLD_BATCH_SUPPORT_gInf_xtof);
        t[xtof_p] = Foam::rlScalarAlgAndSupport(AV_tauxtof,   AV_xtoss_p, Foam::TWORLD_BATCH_SUPPORT_tau_xtof_p, Foam::TWORLD_BATCH_SUPPORT_gInf_xtof_p);
        t[ytos]   = Foam::rlScalarAlgAndSupport(AV_tauytos,   AV_ytoss,   Foam::TWORLD_BATCH_SUPPORT_tau_ytos,   Foam::TWORLD_BATCH_SUPPORT_gInf_ytos);
        t[ytos_p] = Foam::rlScalarAlgAndSupport(AV_tauytos_p, AV_ytoss,   Foam::TWORLD_BATCH_SUPPORT_tau_ytos_p, Foam::TWORLD_BATCH_SUPPORT_gInf_ytos_p);
        t[ytof]   = Foam::rlScalarAlgAndSupport(AV_tauytof,   AV_ytoss,   Foam::TWORLD_BATCH_SUPPORT_tau_ytof,   Foam::TWORLD_BATCH_SUPPORT_gInf_ytof);
        t[ytof_p] = Foam::rlScalarAlgAndSupport(AV_tauytof_p, AV_ytoss,   Foam::TWORLD_BATCH_SUPPORT_tau_ytof_p, Foam::TWORLD_BATCH_SUPPORT_gInf_ytof_p);

        // IKs
        t[xs_junc] = Foam::rlScalarAlgAndSupport(AV_tauxs_junc, AV_xsss_junc, Foam::TWORLD_BATCH_SUPPORT_tau_xs_junc, Foam::TWORLD_BATCH_SUPPORT_gInf_xs_junc);
        t[xs_sl]   = Foam::rlScalarAlgAndSupport(AV_tauxs_sl,   AV_xsss_sl,   Foam::TWORLD_BATCH_SUPPORT_tau_xs_sl,   Foam::TWORLD_BATCH_SUPPORT_gInf_xs_sl);

        // Jrel ICaL-dependent
        t[jrel_icaldep_act] = Foam::rlScalarAlgAndSupport(AV_tau_rel, AV_Jrel_inf,        Foam::TWORLD_BATCH_SUPPORT_tau_jrel_act, Foam::TWORLD_BATCH_SUPPORT_gInf_jrel_act);
        t[jrel_icaldep_f1]  = Foam::rlScalarConstTauAlgInfAndSupport(AC_tauInact,  AV_Jrel_inact_inf,  Foam::TWORLD_BATCH_SUPPORT_tau_jrel_f1,  Foam::TWORLD_BATCH_SUPPORT_gInf_jrel_f1);
        t[jrel_icaldep_f2]  = Foam::rlScalarConstTauAlgInfAndSupport(AC_tauInact2, AV_Jrel_inact_inf2, Foam::TWORLD_BATCH_SUPPORT_tau_jrel_f2,  Foam::TWORLD_BATCH_SUPPORT_gInf_jrel_f2);

        return t;
    }();
}

#ifdef HAS_CUDA
#include <cuda_runtime.h>

namespace Foam
{
    void launchTWorldBatchKernel
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

    void launchTWorldEulerStepKernel
    (
        double* d_STATES,
        const double* d_RATES,
        double dt,
        int N,
        int nStates,
        bool solveVm,
        int vmStateI
    );

    void launchTWorldRushLarsenStepKernel
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

    void launchTWorldScaleIonKernel
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
    defineTypeNameAndDebug(TWorldBatched, 0);
    defineTypeNameAndDebug(TWorldcompactBatched, 0);
    addToRunTimeSelectionTable
    (
        ionicModel,
        TWorldcompactBatched,
        dictionary
    );
}


Foam::TWorldBatched::TWorldBatched
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
        TWorldFamilyInfo()
    ),
    stimulusPOD_()
#ifdef HAS_CUDA
  , useDevice_(false)
#endif
{
    ionicModel::setTissueFromDict();
    ionicModel::setSexFromDict();

    setHotPathSupportSize(NUM_TWORLD_BATCH_SUPPORT);

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
            Info<< "TWorldBatched: rank " << rank << " using CUDA device "
                << devId << " of " << nDevices << nl;
        }
    }
#endif

    double initialRates[NUM_STATES] = {0.0};
    double initialStates[NUM_STATES] = {0.0};

    TWorldinitConsts
    (
        CONSTANTS_.data(),
        initialRates,
        initialStates,
        tissue(),
        sex(),
        dict
    );

    applyIonicConstantOverrides();

    for (label cellI = 0; cellI < nCells(); ++cellI)
    {
        for (label stateI = 0; stateI < NUM_STATES; ++stateI)
        {
            state(cellI, stateI) = initialStates[stateI];
            rate(cellI, stateI)  = initialRates[stateI];
        }
    }

    configurePersistentAlgebraics();
    syncAllToIO();

    if (!utilitiesMode())
    {
        setStimulusProtocolFromDict(dict);
    }
}


Foam::TWorldcompactBatched::TWorldcompactBatched
(
    const dictionary& dict,
    const label num,
    const scalar initialDeltaT,
    const Switch solveVmWithinODESolver
)
:
    TWorldBatched(dict, num, initialDeltaT, solveVmWithinODESolver)
{}


Foam::TWorldBatched::~TWorldBatched()
{
#ifdef HAS_CUDA
    if (useDevice_)
    {
        cuda_.free();
    }
#endif
}


void Foam::TWorldBatched::prepareIOAccess
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
            static_cast<std::size_t>(NUM_TWORLD_BATCH_SUPPORT),
            static_cast<std::size_t>(nCells())
        );
    }
#endif
    configuredBatchedIonicModel::prepareIOAccess(requestedNames, needsAlgebraics);
}


void Foam::TWorldBatched::importFields
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


void Foam::TWorldBatched::solveODE
(
    const scalar stepStartTime,
    const scalar deltaT,
    const scalarField& Vm,
    scalarField& Im
)
{
    solveODEImpl(*this, stepStartTime, deltaT, Vm, Im);
}


void Foam::TWorldcompactBatched::solveODE
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
void Foam::TWorldcompactBatched::solveOnDevice
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
            state(cellI, v) = vmToState(Vm[cellI]);
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
        launchTWorldBatchKernel
        (
            tSub,
            cuda_.d_constants,
            cuda_.d_cellConstants,
            hasHeterogeneousConstants(),
            static_cast<int>(N),
            cuda_.d_states, cuda_.d_rates, cuda_.d_support,
            tFlag, solveVm, stimulusPOD_
        );
        launchTWorldRushLarsenStepKernel
        (
            cuda_.d_states, cuda_.d_rates, cuda_.d_support,
            static_cast<double>(dtSubstep),
            static_cast<int>(N),
            static_cast<int>(NUM_STATES),
            solveVm,
            static_cast<int>(v)
        );
    }

    launchTWorldBatchKernel
    (
        tStart + dtModel,
        cuda_.d_constants,
        cuda_.d_cellConstants,
        hasHeterogeneousConstants(),
        static_cast<int>(N),
        cuda_.d_states, cuda_.d_rates, cuda_.d_support,
        tFlag, solveVm, stimulusPOD_
    );
    launchTWorldScaleIonKernel
    (
        cuda_.d_support, cuda_.d_Im, 1.0,
        static_cast<int>(N),
        static_cast<int>(TWORLD_BATCH_SUPPORT_Iion_cm)
    );
    cuda_.downloadIm(Im.data(), static_cast<std::size_t>(N));
    cuda_.deviceDirty = true;
    setIOEvaluationModelTime(tStart + dtModel);
}
#endif // HAS_CUDA


Foam::List<Foam::word> Foam::TWorldBatched::supportedTissueTypes() const
{
    return {"epicardialCells", "mCells", "endocardialCells"};
}

Foam::List<Foam::word> Foam::TWorldBatched::supportedSexTypes() const
{
    return {"neutral", "male", "female"};
}


Foam::scalarField Foam::TWorldBatched::constantsForTissue
(
    const label tissueFlag
) const
{
    scalarField constants(NUM_CONSTANTS, 0.0);
    scalarField rates(NUM_STATES, 0.0);
    scalarField states(NUM_STATES, 0.0);

    TWorldinitConsts
    (
        constants.data(),
        rates.data(),
        states.data(),
        tissueFlag,
        sex(),
        dict()
    );

    ionicModelIO::applyConstantOverrides
    (
        constants,
        TWorldCONSTANTS_NAMES,
        NUM_CONSTANTS,
        dict(),
        type(),
        tissueFlag
    );

    return constants;
}

Foam::scalarField Foam::TWorldBatched::initialStatesForTissue
(
    const label tissueFlag
) const
{
    scalarField constants(NUM_CONSTANTS, 0.0);
    scalarField rates(NUM_STATES, 0.0);
    scalarField states(NUM_STATES, 0.0);

    TWorldinitConsts
    (
        constants.data(),
        rates.data(),
        states.data(),
        tissueFlag,
        sex(),
        dict()
    );

    return states;
}


void Foam::TWorldBatched::evaluateState
(
    const scalar modelTime,
    const scalarUList& stateValues,
    scalarUList& rateValues,
    scalarUList& algebraicValues
) const
{
    TWorldcomputeVariables
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


void Foam::TWorldBatched::evaluateState
(
    const label cellI,
    const scalar modelTime,
    const scalarUList& stateValues,
    scalarUList& rateValues,
    scalarUList& algebraicValues
) const
{
    scalarField& cellConstants = constants(cellI);

    TWorldcomputeVariables
    (
        modelTime,
        cellConstants.data(),
        rateValues.data(),
        const_cast<scalarUList&>(stateValues).data(),
        algebraicValues.data(),
        solveVmWithinODESolver(),
        stimulusProtocol()
    );
}


Foam::scalar Foam::TWorldBatched::ionicCurrentFromHotPathSupport
(
    const scalarUList& supportValues
) const
{
    return supportValues[Foam::TWORLD_BATCH_SUPPORT_Iion_cm];
}


void Foam::TWorldBatched::evaluateHotPathState
(
    const scalar modelTime,
    const scalarUList& stateValues,
    scalarUList& rateValues,
    scalarUList& supportValues
) const
{
    scalarField algebraics(NUM_ALGEBRAIC, 0.0);
    evaluateState(modelTime, stateValues, rateValues, algebraics);

    for (label stateI = 0; stateI < NUM_STATES; ++stateI)
    {
        projectScalarRushLarsenEntryToSupport
        (
            TWorldRushLarsenDispatch[stateI],
            CONSTANTS_,
            algebraics,
            supportValues
        );
    }
    supportValues[Foam::TWORLD_BATCH_SUPPORT_Iion_cm] = algebraics[Iion_cm];
}


void Foam::TWorldBatched::evaluateHotPathState
(
    const label cellI,
    const scalar modelTime,
    const scalarUList& stateValues,
    scalarUList& rateValues,
    scalarUList& supportValues
) const
{
    scalarField& cellConstants = constants(cellI);
    scalarField algebraics(NUM_ALGEBRAIC, 0.0);
    evaluateState(cellI, modelTime, stateValues, rateValues, algebraics);

    for (label stateI = 0; stateI < NUM_STATES; ++stateI)
    {
        projectScalarRushLarsenEntryToSupport
        (
            TWorldRushLarsenDispatch[stateI],
            cellConstants,
            algebraics,
            supportValues
        );
    }
    supportValues[Foam::TWORLD_BATCH_SUPPORT_Iion_cm] = algebraics[Iion_cm];
}


void Foam::TWorldBatched::derivatives
(
    const scalar t,
    const scalarField& y,
    scalarField& dydt
) const
{
    scalarField algebraics(NUM_ALGEBRAIC, 0.0);
    evaluateState(t, y, dydt, algebraics);
}


bool Foam::TWorldBatched::rushLarsenParametersFromHotPathSupport
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
        TWorldRushLarsenDispatch[stateI],
        CONSTANTS_,
        supportValues,
        VSMALL,
        steadyState,
        tau
    );
}


bool Foam::TWorldBatched::rushLarsenParametersFromHotPathSupport
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
        TWorldRushLarsenDispatch[stateI],
        constants(cellI),
        supportValues,
        VSMALL,
        steadyState,
        tau
    );
}


// ************************************************************************* //
