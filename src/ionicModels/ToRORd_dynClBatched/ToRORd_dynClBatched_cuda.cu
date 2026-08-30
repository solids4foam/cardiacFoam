/*---------------------------------------------------------------------------*\
License
    This file is part of cardiacFoam.
\*---------------------------------------------------------------------------*/

#include <cuda_runtime.h>
#include <cstdio>

#include "ToRORd_dynCl_2020Batch.H"

// ---- file-local kernel error check ----------------------------------------
#define CUDA_LAUNCH_CHECK()                                                    \
    do {                                                                       \
        cudaError_t _err = cudaGetLastError();                                 \
        if (_err != cudaSuccess)                                               \
        {                                                                      \
            fprintf                                                            \
            (                                                                  \
                stderr,                                                        \
                "[cardiacFoam CUDA] kernel error at %s:%d — %s\n",            \
                __FILE__, __LINE__, cudaGetErrorString(_err)                   \
            );                                                                 \
            abort();                                                           \
        }                                                                      \
    } while (0)
// ---------------------------------------------------------------------------

namespace Foam
{
namespace
{
    constexpr int blockSize = 256;

    __global__ void toRORd_dynClBatchKernel
    (
        const double t,
        const double* __restrict__ CONSTANTS,
        const double* __restrict__ CELL_CONSTANTS,
        const bool useCellConstants,
        const int N,
        const double* __restrict__ STATES,
        double* __restrict__ RATES,
        double* __restrict__ SUPPORT,
        const bool solveVm,
        const StimulusProtocolPOD stimulus
    )
    {
        const int cellI = blockIdx.x*blockDim.x + threadIdx.x;
        if (cellI >= N)
        {
            return;
        }

        const double* cellConstants =
            useCellConstants ? CELL_CONSTANTS + cellI*NUM_CONSTANTS : CONSTANTS;

        ToRORd_dynClComputeVariablesBatch
        (
            t, cellConstants, N, cellI, cellI + 1,
            STATES, RATES, SUPPORT,
            solveVm, stimulus
        );
    }


    __global__ void toRORd_dynClEulerStepKernel
    (
        double* __restrict__ STATES,
        const double* __restrict__ RATES,
        const double dt,
        const int N,
        const int nStates,
        const bool solveVm,
        const int vmStateI
    )
    {
        const int cellI = blockIdx.x*blockDim.x + threadIdx.x;
        if (cellI >= N)
        {
            return;
        }

        for (int stateI = 0; stateI < nStates; ++stateI)
        {
            if (!solveVm && stateI == vmStateI)
            {
                continue;
            }

            const int idx = stateI*N + cellI;
            STATES[idx] += dt*RATES[idx];
        }
    }


    __global__ void toRORd_dynClScaleIonKernel
    (
        const double* __restrict__ SUPPORT,
        double* __restrict__ Im,
        const double scale,
        const int N,
        const int IionSlot
    )
    {
        const int cellI = blockIdx.x*blockDim.x + threadIdx.x;
        if (cellI >= N)
        {
            return;
        }

        Im[cellI] = scale*SUPPORT[IionSlot*N + cellI];
    }


    inline int nBlocks(const int N)
    {
        return (N + blockSize - 1)/blockSize;
    }
}


void launchToRORd_dynClBatchKernel
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
)
{
    (void)tissueFlag;
    toRORd_dynClBatchKernel<<<nBlocks(N), blockSize>>>
    (
        t, d_CONSTANTS, d_CELL_CONSTANTS, useCellConstants, N,
        d_STATES, d_RATES, d_SUPPORT,
        solveVm, stimulus
    );
    CUDA_LAUNCH_CHECK();
}


void launchToRORd_dynClEulerStepKernel
(
    double* d_STATES,
    const double* d_RATES,
    double dt,
    int N,
    int nStates,
    bool solveVm,
    int vmStateI
)
{
    toRORd_dynClEulerStepKernel<<<nBlocks(N), blockSize>>>
    (
        d_STATES, d_RATES, dt,
        N, nStates, solveVm, vmStateI
    );
    CUDA_LAUNCH_CHECK();
}


void launchToRORd_dynClScaleIonKernel
(
    const double* d_SUPPORT,
    double* d_Im,
    double scale,
    int N,
    int IionSlot
)
{
    toRORd_dynClScaleIonKernel<<<nBlocks(N), blockSize>>>
    (
        d_SUPPORT, d_Im, scale, N, IionSlot
    );
    CUDA_LAUNCH_CHECK();
}


namespace
{
    __global__ void toRORd_dynClRushLarsenStepKernel
    (
        double* __restrict__ STATES,
        const double* __restrict__ RATES,
        const double* __restrict__ SUPPORT,
        const double dt,
        const int N,
        const int nStates,
        const bool solveVm,
        const int vmStateI
    )
    {
        const int cellI = blockIdx.x*blockDim.x + threadIdx.x;
        if (cellI >= N)
        {
            return;
        }

        #define TORORD_RL(si, tSlot, iSlot)                                  \
        {                                                                     \
            const double _x   = STATES[(si)*N + cellI];                      \
            const double _inf = SUPPORT[(iSlot)*N + cellI];                  \
            const double _tau = SUPPORT[(tSlot)*N + cellI];                  \
            STATES[(si)*N + cellI] = _inf + (_x - _inf)*exp(-dt/_tau);       \
        }

        TORORD_RL(INa_m,      TORORD_DYNCL_BATCH_SUPPORT_tau_INa_m,      TORORD_DYNCL_BATCH_SUPPORT_gInf_INa_m)
        TORORD_RL(INa_h,      TORORD_DYNCL_BATCH_SUPPORT_tau_INa_h,      TORORD_DYNCL_BATCH_SUPPORT_gInf_INa_h)
        TORORD_RL(INa_j,      TORORD_DYNCL_BATCH_SUPPORT_tau_INa_j,      TORORD_DYNCL_BATCH_SUPPORT_gInf_INa_j)
        TORORD_RL(INa_hp,     TORORD_DYNCL_BATCH_SUPPORT_tau_INa_hp,     TORORD_DYNCL_BATCH_SUPPORT_gInf_INa_hp)
        TORORD_RL(INa_jp,     TORORD_DYNCL_BATCH_SUPPORT_tau_INa_jp,     TORORD_DYNCL_BATCH_SUPPORT_gInf_INa_jp)
        TORORD_RL(INaL_mL,    TORORD_DYNCL_BATCH_SUPPORT_tau_INaL_mL,    TORORD_DYNCL_BATCH_SUPPORT_gInf_INaL_mL)
        TORORD_RL(INaL_hL,    TORORD_DYNCL_BATCH_SUPPORT_tau_INaL_hL,    TORORD_DYNCL_BATCH_SUPPORT_gInf_INaL_hL)
        TORORD_RL(INaL_hLp,   TORORD_DYNCL_BATCH_SUPPORT_tau_INaL_hLp,   TORORD_DYNCL_BATCH_SUPPORT_gInf_INaL_hLp)
        TORORD_RL(Ito_a,      TORORD_DYNCL_BATCH_SUPPORT_tau_Ito_a,      TORORD_DYNCL_BATCH_SUPPORT_gInf_Ito_a)
        TORORD_RL(Ito_iF,     TORORD_DYNCL_BATCH_SUPPORT_tau_Ito_iF,     TORORD_DYNCL_BATCH_SUPPORT_gInf_Ito_iF)
        TORORD_RL(Ito_iS,     TORORD_DYNCL_BATCH_SUPPORT_tau_Ito_iS,     TORORD_DYNCL_BATCH_SUPPORT_gInf_Ito_iS)
        TORORD_RL(Ito_ap,     TORORD_DYNCL_BATCH_SUPPORT_tau_Ito_ap,     TORORD_DYNCL_BATCH_SUPPORT_gInf_Ito_ap)
        TORORD_RL(Ito_iFp,    TORORD_DYNCL_BATCH_SUPPORT_tau_Ito_iFp,    TORORD_DYNCL_BATCH_SUPPORT_gInf_Ito_iFp)
        TORORD_RL(Ito_iSp,    TORORD_DYNCL_BATCH_SUPPORT_tau_Ito_iSp,    TORORD_DYNCL_BATCH_SUPPORT_gInf_Ito_iSp)
        TORORD_RL(ICaL_d,     TORORD_DYNCL_BATCH_SUPPORT_tau_ICaL_d,     TORORD_DYNCL_BATCH_SUPPORT_gInf_ICaL_d)
        TORORD_RL(ICaL_ff,    TORORD_DYNCL_BATCH_SUPPORT_tau_ICaL_ff,    TORORD_DYNCL_BATCH_SUPPORT_gInf_ICaL_ff)
        TORORD_RL(ICaL_fs,    TORORD_DYNCL_BATCH_SUPPORT_tau_ICaL_fs,    TORORD_DYNCL_BATCH_SUPPORT_gInf_ICaL_fs)
        TORORD_RL(ICaL_fcaf,  TORORD_DYNCL_BATCH_SUPPORT_tau_ICaL_fcaf,  TORORD_DYNCL_BATCH_SUPPORT_gInf_ICaL_fcaf)
        TORORD_RL(ICaL_fcas,  TORORD_DYNCL_BATCH_SUPPORT_tau_ICaL_fcas,  TORORD_DYNCL_BATCH_SUPPORT_gInf_ICaL_fcas)
        TORORD_RL(ICaL_jca,   TORORD_DYNCL_BATCH_SUPPORT_tau_ICaL_jca,   TORORD_DYNCL_BATCH_SUPPORT_gInf_ICaL_jca)
        TORORD_RL(ICaL_ffp,   TORORD_DYNCL_BATCH_SUPPORT_tau_ICaL_ffp,   TORORD_DYNCL_BATCH_SUPPORT_gInf_ICaL_ffp)
        TORORD_RL(ICaL_fcafp, TORORD_DYNCL_BATCH_SUPPORT_tau_ICaL_fcafp, TORORD_DYNCL_BATCH_SUPPORT_gInf_ICaL_fcafp)
        TORORD_RL(IKs_xs1,    TORORD_DYNCL_BATCH_SUPPORT_tau_IKs_xs1,    TORORD_DYNCL_BATCH_SUPPORT_gInf_IKs_xs1)
        TORORD_RL(IKs_xs2,    TORORD_DYNCL_BATCH_SUPPORT_tau_IKs_xs2,    TORORD_DYNCL_BATCH_SUPPORT_gInf_IKs_xs2)
        TORORD_RL(Jrel_np,    TORORD_DYNCL_BATCH_SUPPORT_tau_Jrel_np,    TORORD_DYNCL_BATCH_SUPPORT_gInf_Jrel_np)
        TORORD_RL(Jrel_p,     TORORD_DYNCL_BATCH_SUPPORT_tau_Jrel_p,     TORORD_DYNCL_BATCH_SUPPORT_gInf_Jrel_p)

        #undef TORORD_RL

        for (int si = 0; si < nStates; ++si)
        {
            if (si == INa_m     || si == INa_h     || si == INa_j     ||
                si == INa_hp    || si == INa_jp    || si == INaL_mL   ||
                si == INaL_hL   || si == INaL_hLp  || si == Ito_a    ||
                si == Ito_iF    || si == Ito_iS    || si == Ito_ap   ||
                si == Ito_iFp   || si == Ito_iSp   || si == ICaL_d   ||
                si == ICaL_ff   || si == ICaL_fs   || si == ICaL_fcaf ||
                si == ICaL_fcas || si == ICaL_jca  || si == ICaL_ffp  ||
                si == ICaL_fcafp|| si == IKs_xs1   || si == IKs_xs2  ||
                si == Jrel_np   || si == Jrel_p) continue;
            if (!solveVm && si == vmStateI) continue;
            const int idx = si*N + cellI;
            STATES[idx] += dt*RATES[idx];
        }
    }
}


void launchToRORd_dynClRushLarsenStepKernel
(
    double* d_STATES,
    const double* d_RATES,
    const double* d_SUPPORT,
    double dt,
    int N,
    int nStates,
    bool solveVm,
    int vmStateI
)
{
    toRORd_dynClRushLarsenStepKernel<<<nBlocks(N), blockSize>>>
    (
        d_STATES, d_RATES, d_SUPPORT,
        dt, N, nStates, solveVm, vmStateI
    );
    CUDA_LAUNCH_CHECK();
}

} // End namespace Foam

// ************************************************************************* //
