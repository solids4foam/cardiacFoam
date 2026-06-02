/*---------------------------------------------------------------------------*\
License
    This file is part of cardiacFoam.
\*---------------------------------------------------------------------------*/

#include <cuda_runtime.h>
#include <cstdio>

#include "Courtemanche_1998Batch.H"

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

    __global__ void courtemancheBatchKernel
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

        CourtemancheComputeVariablesBatch
        (
            t, CONSTANTS, N, cellI, cellI + 1,
            STATES, RATES, SUPPORT,
            solveVm, stimulus
        );
    }


    __global__ void courtemancheEulerStepKernel
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


    __global__ void courtemancheScaleIonKernel
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


void launchCourtemancheBatchKernel
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
    courtemancheBatchKernel<<<nBlocks(N), blockSize>>>
    (
        t, d_CONSTANTS, N,
        d_STATES, d_RATES, d_SUPPORT,
        solveVm, stimulus
    );
    CUDA_LAUNCH_CHECK();
}


void launchCourtemancheEulerStepKernel
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
    courtemancheEulerStepKernel<<<nBlocks(N), blockSize>>>
    (
        d_STATES, d_RATES, dt,
        N, nStates, solveVm, vmStateI
    );
    CUDA_LAUNCH_CHECK();
}


void launchCourtemancheScaleIonKernel
(
    const double* d_SUPPORT,
    double* d_Im,
    double scale,
    int N,
    int IionSlot
)
{
    courtemancheScaleIonKernel<<<nBlocks(N), blockSize>>>
    (
        d_SUPPORT, d_Im, scale, N, IionSlot
    );
    CUDA_LAUNCH_CHECK();
}


namespace
{
    __global__ void courtemancheRushLarsenStepKernel
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

        #define COURTEMANCHE_RL(si, tSlot, iSlot)                            \
        {                                                                     \
            const double _x   = STATES[(si)*N + cellI];                      \
            const double _inf = SUPPORT[(iSlot)*N + cellI];                  \
            const double _tau = SUPPORT[(tSlot)*N + cellI];                  \
            STATES[(si)*N + cellI] = _inf + (_x - _inf)*exp(-dt/_tau);       \
        }

        COURTEMANCHE_RL(ina_m,    COURTEMANCHE_BATCH_SUPPORT_tau_m,       COURTEMANCHE_BATCH_SUPPORT_gInf_m)
        COURTEMANCHE_RL(ina_h,    COURTEMANCHE_BATCH_SUPPORT_tau_h,       COURTEMANCHE_BATCH_SUPPORT_gInf_h)
        COURTEMANCHE_RL(ina_j,    COURTEMANCHE_BATCH_SUPPORT_tau_j,       COURTEMANCHE_BATCH_SUPPORT_gInf_j)
        COURTEMANCHE_RL(ical_d,   COURTEMANCHE_BATCH_SUPPORT_tau_d,       COURTEMANCHE_BATCH_SUPPORT_gInf_d)
        COURTEMANCHE_RL(ical_f,   COURTEMANCHE_BATCH_SUPPORT_tau_f,       COURTEMANCHE_BATCH_SUPPORT_gInf_f)
        COURTEMANCHE_RL(ical_fCa, COURTEMANCHE_BATCH_SUPPORT_tau_fCa,     COURTEMANCHE_BATCH_SUPPORT_gInf_fCa)
        COURTEMANCHE_RL(ito_oa,   COURTEMANCHE_BATCH_SUPPORT_tau_oa,      COURTEMANCHE_BATCH_SUPPORT_gInf_oa)
        COURTEMANCHE_RL(ito_oi,   COURTEMANCHE_BATCH_SUPPORT_tau_oi,      COURTEMANCHE_BATCH_SUPPORT_gInf_oi)
        COURTEMANCHE_RL(ikur_ua,  COURTEMANCHE_BATCH_SUPPORT_tau_ua,      COURTEMANCHE_BATCH_SUPPORT_gInf_ua)
        COURTEMANCHE_RL(ikur_ui,  COURTEMANCHE_BATCH_SUPPORT_tau_ui,      COURTEMANCHE_BATCH_SUPPORT_gInf_ui)
        COURTEMANCHE_RL(ikr_xr,   COURTEMANCHE_BATCH_SUPPORT_tau_xr,      COURTEMANCHE_BATCH_SUPPORT_gInf_xr)
        COURTEMANCHE_RL(iks_xs,   COURTEMANCHE_BATCH_SUPPORT_tau_xs,      COURTEMANCHE_BATCH_SUPPORT_gInf_xs)
        COURTEMANCHE_RL(cajsr_u,  COURTEMANCHE_BATCH_SUPPORT_tau_cajsr_u, COURTEMANCHE_BATCH_SUPPORT_gInf_cajsr_u)
        COURTEMANCHE_RL(cajsr_w,  COURTEMANCHE_BATCH_SUPPORT_tau_cajsr_w, COURTEMANCHE_BATCH_SUPPORT_gInf_cajsr_w)

        #undef COURTEMANCHE_RL

        for (int si = 0; si < nStates; ++si)
        {
            if (si == ina_m  || si == ina_h   || si == ina_j   ||
                si == ical_d || si == ical_f   || si == ical_fCa ||
                si == ito_oa || si == ito_oi   || si == ikur_ua  ||
                si == ikur_ui|| si == ikr_xr   || si == iks_xs   ||
                si == cajsr_u|| si == cajsr_w) continue;
            if (!solveVm && si == vmStateI) continue;
            const int idx = si*N + cellI;
            STATES[idx] += dt*RATES[idx];
        }
    }
}


void launchCourtemancheRushLarsenStepKernel
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
    courtemancheRushLarsenStepKernel<<<nBlocks(N), blockSize>>>
    (
        d_STATES, d_RATES, d_SUPPORT,
        dt, N, nStates, solveVm, vmStateI
    );
    CUDA_LAUNCH_CHECK();
}

} // End namespace Foam

// ************************************************************************* //
