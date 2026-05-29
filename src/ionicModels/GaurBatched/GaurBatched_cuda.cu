/*---------------------------------------------------------------------------*\
License
    This file is part of cardiacFoam.
\*---------------------------------------------------------------------------*/

#include <cuda_runtime.h>
#include <cstdio>

#include "Gaur_2021Batch.H"

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

    __global__ void gaurBatchKernel
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

        GaurComputeVariablesBatch
        (
            t, CONSTANTS, N, cellI, cellI + 1,
            STATES, RATES, SUPPORT,
            solveVm, stimulus
        );
    }


    __global__ void gaurEulerStepKernel
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


    __global__ void gaurScaleIonKernel
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


void launchGaurBatchKernel
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
)
{
    // GPU path requires homogeneous single-tissue mesh.
    // Per-cell tissue heterogeneity is not yet implemented on GPU.
    (void)tissueFlag;
    gaurBatchKernel<<<nBlocks(N), blockSize>>>
    (
        t, d_CONSTANTS, N,
        d_STATES, d_RATES, d_SUPPORT,
        solveVm, stimulus
    );
    CUDA_LAUNCH_CHECK();
}


void launchGaurEulerStepKernel
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
    gaurEulerStepKernel<<<nBlocks(N), blockSize>>>
    (
        d_STATES, d_RATES, dt,
        N, nStates, solveVm, vmStateI
    );
    CUDA_LAUNCH_CHECK();
}


void launchGaurScaleIonKernel
(
    const double* d_SUPPORT,
    double* d_Im,
    double scale,
    int N,
    int IionSlot
)
{
    gaurScaleIonKernel<<<nBlocks(N), blockSize>>>
    (
        d_SUPPORT, d_Im, scale, N, IionSlot
    );
    CUDA_LAUNCH_CHECK();
}


namespace
{
    __global__ void gaurRushLarsenStepKernel
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

        #define GAUR_RL(si, tSlot, iSlot)                                    \
        {                                                                     \
            const double _x   = STATES[(si)*N + cellI];                      \
            const double _inf = SUPPORT[(iSlot)*N + cellI];                  \
            const double _tau = SUPPORT[(tSlot)*N + cellI];                  \
            STATES[(si)*N + cellI] = _inf + (_x - _inf)*exp(-dt/_tau);       \
        }

        GAUR_RL(I_Na_m,     GAUR_BATCH_SUPPORT_tau_m,      GAUR_BATCH_SUPPORT_gInf_m)
        GAUR_RL(I_Na_h,     GAUR_BATCH_SUPPORT_tau_h,      GAUR_BATCH_SUPPORT_gInf_h)
        GAUR_RL(I_Na_j,     GAUR_BATCH_SUPPORT_tau_j,      GAUR_BATCH_SUPPORT_gInf_j)
        GAUR_RL(INaL_ml,    GAUR_BATCH_SUPPORT_tau_ml,     GAUR_BATCH_SUPPORT_gInf_ml)
        GAUR_RL(INaL_hl,    GAUR_BATCH_SUPPORT_tau_hl,     GAUR_BATCH_SUPPORT_gInf_hl)
        GAUR_RL(ICaL_d,     GAUR_BATCH_SUPPORT_tau_d,      GAUR_BATCH_SUPPORT_gInf_d)
        GAUR_RL(ICaL_fca,   GAUR_BATCH_SUPPORT_tau_fca,    GAUR_BATCH_SUPPORT_gInf_fca)
        GAUR_RL(ICaL_ff,    GAUR_BATCH_SUPPORT_tau_ff,     GAUR_BATCH_SUPPORT_gInf_ff)
        GAUR_RL(ICaL_fs,    GAUR_BATCH_SUPPORT_tau_fs,     GAUR_BATCH_SUPPORT_gInf_fs)
        GAUR_RL(IKr_xr,     GAUR_BATCH_SUPPORT_tau_xr,     GAUR_BATCH_SUPPORT_gInf_xr)
        GAUR_RL(IKs_xs1,    GAUR_BATCH_SUPPORT_tau_xs1,    GAUR_BATCH_SUPPORT_gInf_xs1)
        GAUR_RL(IKs_xs2,    GAUR_BATCH_SUPPORT_tau_xs2,    GAUR_BATCH_SUPPORT_gInf_xs2)
        GAUR_RL(CICR_Jrel1, GAUR_BATCH_SUPPORT_tau_Jrel1,  GAUR_BATCH_SUPPORT_gInf_Jrel1)
        GAUR_RL(CICR_Jrel2, GAUR_BATCH_SUPPORT_tau_Jrel2,  GAUR_BATCH_SUPPORT_gInf_Jrel2)

        #undef GAUR_RL

        for (int si = 0; si < nStates; ++si)
        {
            if (si == I_Na_m   || si == I_Na_h   || si == I_Na_j   ||
                si == INaL_ml  || si == INaL_hl  || si == ICaL_d   ||
                si == ICaL_fca || si == ICaL_ff  || si == ICaL_fs  ||
                si == IKr_xr   || si == IKs_xs1  || si == IKs_xs2  ||
                si == CICR_Jrel1 || si == CICR_Jrel2) continue;
            if (!solveVm && si == vmStateI) continue;
            const int idx = si*N + cellI;
            STATES[idx] += dt*RATES[idx];
        }
    }
}


void launchGaurRushLarsenStepKernel
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
    gaurRushLarsenStepKernel<<<nBlocks(N), blockSize>>>
    (
        d_STATES, d_RATES, d_SUPPORT,
        dt, N, nStates, solveVm, vmStateI
    );
    CUDA_LAUNCH_CHECK();
}

} // End namespace Foam

// ************************************************************************* //
