/*---------------------------------------------------------------------------*\
License
    This file is part of cardiacFoam.
\*---------------------------------------------------------------------------*/

#include <cuda_runtime.h>
#include <cstdio>

#include "Grandi_2011Batch.H"

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

    __global__ void grandiBatchKernel
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

        GrandiComputeVariablesBatch
        (
            t, CONSTANTS, N, cellI, cellI + 1,
            STATES, RATES, SUPPORT,
            solveVm, stimulus
        );
    }


    __global__ void grandiEulerStepKernel
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


    __global__ void grandiScaleIonKernel
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


void launchGrandiBatchKernel
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
    grandiBatchKernel<<<nBlocks(N), blockSize>>>
    (
        t, d_CONSTANTS, N,
        d_STATES, d_RATES, d_SUPPORT,
        solveVm, stimulus
    );
    CUDA_LAUNCH_CHECK();
}


void launchGrandiEulerStepKernel
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
    grandiEulerStepKernel<<<nBlocks(N), blockSize>>>
    (
        d_STATES, d_RATES, dt,
        N, nStates, solveVm, vmStateI
    );
    CUDA_LAUNCH_CHECK();
}


void launchGrandiScaleIonKernel
(
    const double* d_SUPPORT,
    double* d_Im,
    double scale,
    int N,
    int IionSlot
)
{
    grandiScaleIonKernel<<<nBlocks(N), blockSize>>>
    (
        d_SUPPORT, d_Im, scale, N, IionSlot
    );
    CUDA_LAUNCH_CHECK();
}


namespace
{
    __global__ void grandiRushLarsenStepKernel
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

        #define GRANDI_RL(si, tSlot, iSlot)                                  \
        {                                                                     \
            const double _x   = STATES[(si)*N + cellI];                      \
            const double _inf = SUPPORT[(iSlot)*N + cellI];                  \
            const double _tau = SUPPORT[(tSlot)*N + cellI];                  \
            STATES[(si)*N + cellI] = _inf + (_x - _inf)*exp(-dt/_tau);       \
        }

        GRANDI_RL(Ikr_xr,        GRANDI_BATCH_SUPPORT_tau_Ikr_xr,        GRANDI_BATCH_SUPPORT_gInf_Ikr_xr)
        GRANDI_RL(Iks_xs,        GRANDI_BATCH_SUPPORT_tau_Iks_xs,        GRANDI_BATCH_SUPPORT_gInf_Iks_xs)
        GRANDI_RL(Ikur_ikur_r,   GRANDI_BATCH_SUPPORT_tau_Ikur_ikur_r,   GRANDI_BATCH_SUPPORT_gInf_Ikur_ikur_r)
        GRANDI_RL(Ikur_s,        GRANDI_BATCH_SUPPORT_tau_Ikur_s,        GRANDI_BATCH_SUPPORT_gInf_Ikur_s)
        GRANDI_RL(Ina_h,         GRANDI_BATCH_SUPPORT_tau_Ina_h,         GRANDI_BATCH_SUPPORT_gInf_Ina_h)
        GRANDI_RL(Ina_j,         GRANDI_BATCH_SUPPORT_tau_Ina_j,         GRANDI_BATCH_SUPPORT_gInf_Ina_j)
        GRANDI_RL(Ina_m,         GRANDI_BATCH_SUPPORT_tau_Ina_m,         GRANDI_BATCH_SUPPORT_gInf_Ina_m)
        GRANDI_RL(Inal_hl,       GRANDI_BATCH_SUPPORT_tau_Inal_hl,       GRANDI_BATCH_SUPPORT_gInf_Inal_hl)
        GRANDI_RL(Inal_ml,       GRANDI_BATCH_SUPPORT_tau_Inal_ml,       GRANDI_BATCH_SUPPORT_gInf_Inal_ml)
        GRANDI_RL(Ical_d,        GRANDI_BATCH_SUPPORT_tau_Ical_d,        GRANDI_BATCH_SUPPORT_gInf_Ical_d)
        GRANDI_RL(Ical_f,        GRANDI_BATCH_SUPPORT_tau_Ical_f,        GRANDI_BATCH_SUPPORT_gInf_Ical_f)
        GRANDI_RL(Ical_fCaB_jn,  GRANDI_BATCH_SUPPORT_tau_Ical_fCaB_jn, GRANDI_BATCH_SUPPORT_gInf_Ical_fCaB_jn)
        GRANDI_RL(Ical_fCaB_sl,  GRANDI_BATCH_SUPPORT_tau_Ical_fCaB_sl, GRANDI_BATCH_SUPPORT_gInf_Ical_fCaB_sl)
        GRANDI_RL(Ito_x,         GRANDI_BATCH_SUPPORT_tau_Ito_x,         GRANDI_BATCH_SUPPORT_gInf_Ito_x)
        GRANDI_RL(Ito_y,         GRANDI_BATCH_SUPPORT_tau_Ito_y,         GRANDI_BATCH_SUPPORT_gInf_Ito_y)

        #undef GRANDI_RL

        for (int si = 0; si < nStates; ++si)
        {
            if (si == Ikr_xr      || si == Iks_xs      || si == Ikur_ikur_r ||
                si == Ikur_s      || si == Ina_h       || si == Ina_j       ||
                si == Ina_m       || si == Inal_hl     || si == Inal_ml     ||
                si == Ical_d      || si == Ical_f      || si == Ical_fCaB_jn||
                si == Ical_fCaB_sl|| si == Ito_x       || si == Ito_y) continue;
            if (!solveVm && si == vmStateI) continue;
            const int idx = si*N + cellI;
            STATES[idx] += dt*RATES[idx];
        }
    }
}


void launchGrandiRushLarsenStepKernel
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
    grandiRushLarsenStepKernel<<<nBlocks(N), blockSize>>>
    (
        d_STATES, d_RATES, d_SUPPORT,
        dt, N, nStates, solveVm, vmStateI
    );
    CUDA_LAUNCH_CHECK();
}

} // End namespace Foam

// ************************************************************************* //
