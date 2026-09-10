/*---------------------------------------------------------------------------*\
License
    This file is part of cardiacFoam.
\*---------------------------------------------------------------------------*/

#include <cuda_runtime.h>
#include <cstdio>

#include "Trovato_2020Batch.H"

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

    __global__ void trovatoBatchKernel
    (
        const double t,
        const double* __restrict__ CONSTANTS,
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

        TrovatoComputeVariablesBatch
        (
            t, CONSTANTS, N, cellI, cellI + 1,
            STATES, RATES, SUPPORT,
            solveVm, stimulus
        );
    }


    __global__ void trovatoEulerStepKernel
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


    __global__ void trovatoScaleIonKernel
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


void launchTrovatoBatchKernel
(
    double t,
    const double* d_CONSTANTS,
    int N,
    const double* d_STATES,
    double* d_RATES,
    double* d_SUPPORT,
    bool solveVm,
    StimulusProtocolPOD stimulus
)
{
    trovatoBatchKernel<<<nBlocks(N), blockSize>>>
    (
        t, d_CONSTANTS, N,
        d_STATES, d_RATES, d_SUPPORT,
        solveVm, stimulus
    );
    CUDA_LAUNCH_CHECK();
}


void launchTrovatoEulerStepKernel
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
    trovatoEulerStepKernel<<<nBlocks(N), blockSize>>>
    (
        d_STATES, d_RATES, dt,
        N, nStates, solveVm, vmStateI
    );
    CUDA_LAUNCH_CHECK();
}


void launchTrovatoScaleIonKernel
(
    const double* d_SUPPORT,
    double* d_Im,
    double scale,
    int N,
    int IionSlot
)
{
    trovatoScaleIonKernel<<<nBlocks(N), blockSize>>>
    (
        d_SUPPORT, d_Im, scale, N, IionSlot
    );
    CUDA_LAUNCH_CHECK();
}


namespace
{
    __global__ void trovatoRushLarsenStepKernel
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
        if (cellI >= N) { return; }

        #define TROVATO_RL(si, tSlot, iSlot)                                 \
        {                                                                     \
            const double _x   = STATES[(si)*N + cellI];                      \
            const double _inf = SUPPORT[(iSlot)*N + cellI];                  \
            const double _tau = SUPPORT[(tSlot)*N + cellI];                  \
            STATES[(si)*N + cellI] = _inf + (_x - _inf)*exp(-dt/_tau);       \
        }

        TROVATO_RL(ICaT_b,      TROVATO_BATCH_SUPPORT_tau_ICaT_b,     TROVATO_BATCH_SUPPORT_gInf_ICaT_b)
        TROVATO_RL(ICaT_g,      TROVATO_BATCH_SUPPORT_tau_ICaT_g,     TROVATO_BATCH_SUPPORT_gInf_ICaT_g)
        TROVATO_RL(IK1_xk1,     TROVATO_BATCH_SUPPORT_tau_IK1_xk1,    TROVATO_BATCH_SUPPORT_gInf_IK1_xk1)
        TROVATO_RL(IKr_xrf,     TROVATO_BATCH_SUPPORT_tau_IKr_xrf,    TROVATO_BATCH_SUPPORT_gInf_IKr_xrf)
        TROVATO_RL(IKr_xrs,     TROVATO_BATCH_SUPPORT_tau_IKr_xrs,    TROVATO_BATCH_SUPPORT_gInf_IKr_xrs)
        TROVATO_RL(IKs_xs1,     TROVATO_BATCH_SUPPORT_tau_IKs_xs1,    TROVATO_BATCH_SUPPORT_gInf_IKs_xs1)
        TROVATO_RL(IKs_xs2,     TROVATO_BATCH_SUPPORT_tau_IKs_xs2,    TROVATO_BATCH_SUPPORT_gInf_IKs_xs2)
        TROVATO_RL(INa_hf,      TROVATO_BATCH_SUPPORT_tau_INa_hf,     TROVATO_BATCH_SUPPORT_gInf_INa_hf)
        TROVATO_RL(INa_hs,      TROVATO_BATCH_SUPPORT_tau_INa_hs,     TROVATO_BATCH_SUPPORT_gInf_INa_hs)
        TROVATO_RL(INa_hsp,     TROVATO_BATCH_SUPPORT_tau_INa_hsp,    TROVATO_BATCH_SUPPORT_gInf_INa_hsp)
        TROVATO_RL(INa_m,       TROVATO_BATCH_SUPPORT_tau_INa_m,      TROVATO_BATCH_SUPPORT_gInf_INa_m)
        TROVATO_RL(INa_j,       TROVATO_BATCH_SUPPORT_tau_INa_j,      TROVATO_BATCH_SUPPORT_gInf_INa_j)
        TROVATO_RL(INa_jp,      TROVATO_BATCH_SUPPORT_tau_INa_jp,     TROVATO_BATCH_SUPPORT_gInf_INa_jp)
        TROVATO_RL(If_y,        TROVATO_BATCH_SUPPORT_tau_If_y,       TROVATO_BATCH_SUPPORT_gInf_If_y)
        TROVATO_RL(Ito_a,       TROVATO_BATCH_SUPPORT_tau_Ito_a,      TROVATO_BATCH_SUPPORT_gInf_Ito_a)
        TROVATO_RL(Ito_i1,      TROVATO_BATCH_SUPPORT_tau_Ito_i1,     TROVATO_BATCH_SUPPORT_gInf_Ito_i1)
        TROVATO_RL(Ito_i2,      TROVATO_BATCH_SUPPORT_tau_Ito_i2,     TROVATO_BATCH_SUPPORT_gInf_Ito_i2)
        TROVATO_RL(INaL_hL,     TROVATO_BATCH_SUPPORT_tau_INaL_hL,    TROVATO_BATCH_SUPPORT_gInf_INaL_hL)
        TROVATO_RL(INaL_mL,     TROVATO_BATCH_SUPPORT_tau_INaL_mL,    TROVATO_BATCH_SUPPORT_gInf_INaL_mL)
        TROVATO_RL(INaL_hLp,    TROVATO_BATCH_SUPPORT_tau_INaL_hLp,   TROVATO_BATCH_SUPPORT_gInf_INaL_hLp)
        TROVATO_RL(ICaL_d,      TROVATO_BATCH_SUPPORT_tau_ICaL_d,     TROVATO_BATCH_SUPPORT_gInf_ICaL_d)
        TROVATO_RL(ICaL_ff,     TROVATO_BATCH_SUPPORT_tau_ICaL_ff,    TROVATO_BATCH_SUPPORT_gInf_ICaL_ff)
        TROVATO_RL(ICaL_fs,     TROVATO_BATCH_SUPPORT_tau_ICaL_fs,    TROVATO_BATCH_SUPPORT_gInf_ICaL_fs)
        TROVATO_RL(ICaL_fcaf,   TROVATO_BATCH_SUPPORT_tau_ICaL_fcaf,  TROVATO_BATCH_SUPPORT_gInf_ICaL_fcaf)
        TROVATO_RL(ICaL_fcafp,  TROVATO_BATCH_SUPPORT_tau_ICaL_fcafp, TROVATO_BATCH_SUPPORT_gInf_ICaL_fcafp)
        TROVATO_RL(ICaL_fcas,   TROVATO_BATCH_SUPPORT_tau_ICaL_fcas,  TROVATO_BATCH_SUPPORT_gInf_ICaL_fcas)
        TROVATO_RL(ICaL_ffp,    TROVATO_BATCH_SUPPORT_tau_ICaL_ffp,   TROVATO_BATCH_SUPPORT_gInf_ICaL_ffp)
        TROVATO_RL(ICaL_jca,    TROVATO_BATCH_SUPPORT_tau_ICaL_jca,   TROVATO_BATCH_SUPPORT_gInf_ICaL_jca)
        TROVATO_RL(ryr_Jrel1,   TROVATO_BATCH_SUPPORT_tau_ryr_Jrel1,  TROVATO_BATCH_SUPPORT_gInf_ryr_Jrel1)
        TROVATO_RL(ryr_Jrel2,   TROVATO_BATCH_SUPPORT_tau_ryr_Jrel2,  TROVATO_BATCH_SUPPORT_gInf_ryr_Jrel2)

        #undef TROVATO_RL

        for (int si = 0; si < nStates; ++si)
        {
            if (si == ICaT_b    || si == ICaT_g    || si == IK1_xk1   ||
                si == IKr_xrf   || si == IKr_xrs   || si == IKs_xs1   ||
                si == IKs_xs2   || si == INa_hf    || si == INa_hs    ||
                si == INa_hsp   || si == INa_m     || si == INa_j     ||
                si == INa_jp    || si == If_y      || si == Ito_a     ||
                si == Ito_i1    || si == Ito_i2    || si == INaL_hL   ||
                si == INaL_mL   || si == INaL_hLp  || si == ICaL_d    ||
                si == ICaL_ff   || si == ICaL_fs   || si == ICaL_fcaf  ||
                si == ICaL_fcafp|| si == ICaL_fcas  || si == ICaL_ffp  ||
                si == ICaL_jca  || si == ryr_Jrel1  || si == ryr_Jrel2) continue;
            if (!solveVm && si == vmStateI) continue;
            const int idx = si*N + cellI;
            STATES[idx] += dt*RATES[idx];
        }
    }
}


void launchTrovatoRushLarsenStepKernel
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
    trovatoRushLarsenStepKernel<<<nBlocks(N), blockSize>>>
    (
        d_STATES, d_RATES, d_SUPPORT,
        dt, N, nStates, solveVm, vmStateI
    );
    CUDA_LAUNCH_CHECK();
}

} // End namespace Foam

// ************************************************************************* //
