/*---------------------------------------------------------------------------*\
License
    This file is part of cardiacFoam.
\*---------------------------------------------------------------------------*/

#include <cuda_runtime.h>
#include <cstdio>

#include "TNNP_2004Batch.H"

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

    __global__ void tnnpBatchKernel
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

        TNNPComputeVariablesBatch
        (
            t, CONSTANTS, N, cellI, cellI + 1,
            STATES, RATES, SUPPORT,
            solveVm, stimulus
        );
    }


    __global__ void tnnpEulerStepKernel
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


    __global__ void tnnpScaleIonKernel
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
)
{
    // GPU path requires homogeneous single-tissue mesh.
    // Per-cell tissue heterogeneity is not yet implemented on GPU.
    (void)tissueFlag;
    tnnpBatchKernel<<<nBlocks(N), blockSize>>>
    (
        t, d_CONSTANTS, N,
        d_STATES, d_RATES, d_SUPPORT,
        solveVm, stimulus
    );
    CUDA_LAUNCH_CHECK();
}


void launchTnnpEulerStepKernel
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
    tnnpEulerStepKernel<<<nBlocks(N), blockSize>>>
    (
        d_STATES, d_RATES, dt,
        N, nStates, solveVm, vmStateI
    );
    CUDA_LAUNCH_CHECK();
}


void launchTnnpScaleIonKernel
(
    const double* d_SUPPORT,
    double* d_Im,
    double scale,
    int N,
    int IionSlot
)
{
    tnnpScaleIonKernel<<<nBlocks(N), blockSize>>>
    (
        d_SUPPORT, d_Im, scale, N, IionSlot
    );
    CUDA_LAUNCH_CHECK();
}


namespace
{
    __global__ void tnnpRushLarsenStepKernel
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

        #define TNNP_RL(si, tSlot, iSlot)                                    \
        {                                                                     \
            const double _x   = STATES[(si)*N + cellI];                      \
            const double _inf = SUPPORT[(iSlot)*N + cellI];                  \
            const double _tau = SUPPORT[(tSlot)*N + cellI];                  \
            STATES[(si)*N + cellI] = _inf + (_x - _inf)*exp(-dt/_tau);       \
        }

        TNNP_RL(Xr1, TNNP_BATCH_SUPPORT_tau_xr1, TNNP_BATCH_SUPPORT_gInf_xr1)
        TNNP_RL(Xr2, TNNP_BATCH_SUPPORT_tau_xr2, TNNP_BATCH_SUPPORT_gInf_xr2)
        TNNP_RL(Xs,  TNNP_BATCH_SUPPORT_tau_xs,  TNNP_BATCH_SUPPORT_gInf_xs)
        TNNP_RL(m,   TNNP_BATCH_SUPPORT_tau_m,   TNNP_BATCH_SUPPORT_gInf_m)
        TNNP_RL(h,   TNNP_BATCH_SUPPORT_tau_h,   TNNP_BATCH_SUPPORT_gInf_h)
        TNNP_RL(j,   TNNP_BATCH_SUPPORT_tau_j,   TNNP_BATCH_SUPPORT_gInf_j)
        TNNP_RL(d,   TNNP_BATCH_SUPPORT_tau_d,   TNNP_BATCH_SUPPORT_gInf_d)
        TNNP_RL(f,   TNNP_BATCH_SUPPORT_tau_f,   TNNP_BATCH_SUPPORT_gInf_f)
        TNNP_RL(s,   TNNP_BATCH_SUPPORT_tau_s,   TNNP_BATCH_SUPPORT_gInf_s)
        TNNP_RL(r,   TNNP_BATCH_SUPPORT_tau_r,   TNNP_BATCH_SUPPORT_gInf_r)

        #undef TNNP_RL

        for (int si = 0; si < nStates; ++si)
        {
            if (si == Xr1 || si == Xr2 || si == Xs || si == m ||
                si == h  || si == j   || si == d  || si == f  ||
                si == s  || si == r) continue;
            if (!solveVm && si == vmStateI) continue;
            const int idx = si*N + cellI;
            STATES[idx] += dt*RATES[idx];
        }
    }
}


void launchTnnpRushLarsenStepKernel
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
    tnnpRushLarsenStepKernel<<<nBlocks(N), blockSize>>>
    (
        d_STATES, d_RATES, d_SUPPORT,
        dt, N, nStates, solveVm, vmStateI
    );
    CUDA_LAUNCH_CHECK();
}

} // End namespace Foam

// ************************************************************************* //
