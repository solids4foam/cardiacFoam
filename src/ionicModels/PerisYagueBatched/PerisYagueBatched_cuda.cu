/*---------------------------------------------------------------------------*\
License
    This file is part of cardiacFoam.
\*---------------------------------------------------------------------------*/

#include <cuda_runtime.h>
#include <cstdio>
#include <cstdlib>

#include "PerisYague_2022Batch.H"

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
            std::abort();                                                           \
        }                                                                      \
    } while (0)
// ---------------------------------------------------------------------------

namespace Foam
{
namespace
{
    constexpr int blockSize = 256;

    __global__ void perisYagueBatchKernel
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

        PerisYagueComputeVariablesBatch
        (
            t, CONSTANTS, N, cellI, cellI + 1,
            STATES, RATES, SUPPORT,
            solveVm, stimulus
        );
    }


    __global__ void perisYagueEulerStepKernel
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


    __global__ void perisYagueScaleIonKernel
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


void launchPerisYagueBatchKernel
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
    perisYagueBatchKernel<<<nBlocks(N), blockSize>>>
    (
        t, d_CONSTANTS, N,
        d_STATES, d_RATES, d_SUPPORT,
        solveVm, stimulus
    );
    CUDA_LAUNCH_CHECK();
}


void launchPerisYagueEulerStepKernel
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
    perisYagueEulerStepKernel<<<nBlocks(N), blockSize>>>
    (
        d_STATES, d_RATES, dt,
        N, nStates, solveVm, vmStateI
    );
    CUDA_LAUNCH_CHECK();
}


void launchPerisYagueScaleIonKernel
(
    const double* d_SUPPORT,
    double* d_Im,
    double scale,
    int N,
    int IionSlot
)
{
    perisYagueScaleIonKernel<<<nBlocks(N), blockSize>>>
    (
        d_SUPPORT, d_Im, scale, N, IionSlot
    );
    CUDA_LAUNCH_CHECK();
}


namespace
{
    __global__ void perisYagueRushLarsenStepKernel
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

        #define PERISYAGUE_RL(si, tSlot, iSlot)                              \
        {                                                                     \
            const double _x   = STATES[(si)*N + cellI];                      \
            const double _inf = SUPPORT[(iSlot)*N + cellI];                  \
            const double _tau = SUPPORT[(tSlot)*N + cellI];                  \
            STATES[(si)*N + cellI] = _inf + (_x - _inf)*::exp(-dt/_tau);       \
        }

        PERISYAGUE_RL(ina_m,     PERISYAGUE_BATCH_SUPPORT_tau_m,      PERISYAGUE_BATCH_SUPPORT_gInf_m)
        PERISYAGUE_RL(ina_h,     PERISYAGUE_BATCH_SUPPORT_tau_h,      PERISYAGUE_BATCH_SUPPORT_gInf_h)
        PERISYAGUE_RL(ina_j,     PERISYAGUE_BATCH_SUPPORT_tau_j,      PERISYAGUE_BATCH_SUPPORT_gInf_j)
        PERISYAGUE_RL(ikr_xr,    PERISYAGUE_BATCH_SUPPORT_tau_xr,     PERISYAGUE_BATCH_SUPPORT_gInf_xr)
        PERISYAGUE_RL(iks_xs,    PERISYAGUE_BATCH_SUPPORT_tau_xs,     PERISYAGUE_BATCH_SUPPORT_gInf_xs)
        PERISYAGUE_RL(ikur_ua,   PERISYAGUE_BATCH_SUPPORT_tau_ua,     PERISYAGUE_BATCH_SUPPORT_gInf_ua)
        PERISYAGUE_RL(ikur_uif,  PERISYAGUE_BATCH_SUPPORT_tau_uif,    PERISYAGUE_BATCH_SUPPORT_gInf_uif)
        PERISYAGUE_RL(ikur_uis,  PERISYAGUE_BATCH_SUPPORT_tau_uis,    PERISYAGUE_BATCH_SUPPORT_gInf_uis)
        PERISYAGUE_RL(ical_d,    PERISYAGUE_BATCH_SUPPORT_tau_d,      PERISYAGUE_BATCH_SUPPORT_gInf_d)
        PERISYAGUE_RL(ical_f,    PERISYAGUE_BATCH_SUPPORT_tau_f,      PERISYAGUE_BATCH_SUPPORT_gInf_f)
        PERISYAGUE_RL(ical_fCa,  PERISYAGUE_BATCH_SUPPORT_tau_fCa,    PERISYAGUE_BATCH_SUPPORT_gInf_fCa)
        PERISYAGUE_RL(iclca_qCa, PERISYAGUE_BATCH_SUPPORT_tau_qCa,    PERISYAGUE_BATCH_SUPPORT_gInf_qCa)
        PERISYAGUE_RL(ryr_u,     PERISYAGUE_BATCH_SUPPORT_tau_ryr_u,  PERISYAGUE_BATCH_SUPPORT_gInf_ryr_u)
        PERISYAGUE_RL(ryr_w,     PERISYAGUE_BATCH_SUPPORT_tau_ryr_w,  PERISYAGUE_BATCH_SUPPORT_gInf_ryr_w)

        #undef PERISYAGUE_RL

        for (int si = 0; si < nStates; ++si)
        {
            if (si == ina_m    || si == ina_h    || si == ina_j    ||
                si == ikr_xr   || si == iks_xs   || si == ikur_ua  ||
                si == ikur_uif || si == ikur_uis  || si == ical_d  ||
                si == ical_f   || si == ical_fCa  || si == iclca_qCa||
                si == ryr_u    || si == ryr_w) continue;
            if (!solveVm && si == vmStateI) continue;
            const int idx = si*N + cellI;
            STATES[idx] += dt*RATES[idx];
        }
    }
}


void launchPerisYagueRushLarsenStepKernel
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
    perisYagueRushLarsenStepKernel<<<nBlocks(N), blockSize>>>
    (
        d_STATES, d_RATES, d_SUPPORT,
        dt, N, nStates, solveVm, vmStateI
    );
    CUDA_LAUNCH_CHECK();
}

} // End namespace Foam

// ************************************************************************* //
