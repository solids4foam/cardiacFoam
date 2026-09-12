/*---------------------------------------------------------------------------*\
License
    This file is part of cardiacFoam.
\*---------------------------------------------------------------------------*/

#include <cuda_runtime.h>
#include <cstdio>

#include "Stewart_2009Batch.H"

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

    __global__ void stewartBatchKernel
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

        StewartComputeVariablesBatch
        (
            t, CONSTANTS, N, cellI, cellI + 1,
            STATES, RATES, SUPPORT,
            solveVm, stimulus
        );
    }


    __global__ void stewartEulerStepKernel
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


    __global__ void stewartScaleIonKernel
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


void launchStewartBatchKernel
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
    stewartBatchKernel<<<nBlocks(N), blockSize>>>
    (
        t, d_CONSTANTS, N,
        d_STATES, d_RATES, d_SUPPORT,
        solveVm, stimulus
    );
    CUDA_LAUNCH_CHECK();
}


void launchStewartEulerStepKernel
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
    stewartEulerStepKernel<<<nBlocks(N), blockSize>>>
    (
        d_STATES, d_RATES, dt,
        N, nStates, solveVm, vmStateI
    );
    CUDA_LAUNCH_CHECK();
}


void launchStewartScaleIonKernel
(
    const double* d_SUPPORT,
    double* d_Im,
    double scale,
    int N,
    int IionSlot
)
{
    stewartScaleIonKernel<<<nBlocks(N), blockSize>>>
    (
        d_SUPPORT, d_Im, scale, N, IionSlot
    );
    CUDA_LAUNCH_CHECK();
}


namespace
{
    __global__ void stewartRushLarsenStepKernel
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

        #define STEWART_RL(si, tSlot, iSlot)                                 \
        {                                                                     \
            const double _x   = STATES[(si)*N + cellI];                      \
            const double _inf = SUPPORT[(iSlot)*N + cellI];                  \
            const double _tau = SUPPORT[(tSlot)*N + cellI];                  \
            STATES[(si)*N + cellI] = _inf + (_x - _inf)*exp(-dt/_tau);       \
        }

        STEWART_RL(Ihyperpolarization_activated_current_y_gate_y,              STEWART_BATCH_SUPPORT_tau_y,    STEWART_BATCH_SUPPORT_gInf_y)
        STEWART_RL(Irapid_time_dependent_potassium_current_Xr1_gate_Xr1,       STEWART_BATCH_SUPPORT_tau_xr1,  STEWART_BATCH_SUPPORT_gInf_xr1)
        STEWART_RL(Irapid_time_dependent_potassium_current_Xr2_gate_Xr2,       STEWART_BATCH_SUPPORT_tau_xr2,  STEWART_BATCH_SUPPORT_gInf_xr2)
        STEWART_RL(Islow_time_dependent_potassium_current_Xs_gate_Xs,          STEWART_BATCH_SUPPORT_tau_xs,   STEWART_BATCH_SUPPORT_gInf_xs)
        STEWART_RL(Ifast_sodium_current_m_gate_m,                              STEWART_BATCH_SUPPORT_tau_m,    STEWART_BATCH_SUPPORT_gInf_m)
        STEWART_RL(Ifast_sodium_current_h_gate_h,                              STEWART_BATCH_SUPPORT_tau_h,    STEWART_BATCH_SUPPORT_gInf_h)
        STEWART_RL(Ifast_sodium_current_j_gate_j,                              STEWART_BATCH_SUPPORT_tau_j,    STEWART_BATCH_SUPPORT_gInf_j)
        STEWART_RL(IL_type_Ca_current_d_gate_d,                                STEWART_BATCH_SUPPORT_tau_d,    STEWART_BATCH_SUPPORT_gInf_d)
        STEWART_RL(IL_type_Ca_current_f_gate_f,                                STEWART_BATCH_SUPPORT_tau_f,    STEWART_BATCH_SUPPORT_gInf_f)
        STEWART_RL(IL_type_Ca_current_f2_gate_f2,                              STEWART_BATCH_SUPPORT_tau_f2,   STEWART_BATCH_SUPPORT_gInf_f2)
        STEWART_RL(IL_type_Ca_current_fCass_gate_fCass,                        STEWART_BATCH_SUPPORT_tau_fCass,STEWART_BATCH_SUPPORT_gInf_fCass)
        STEWART_RL(Itransient_outward_current_s_gate_s,                        STEWART_BATCH_SUPPORT_tau_s,    STEWART_BATCH_SUPPORT_gInf_s)
        STEWART_RL(Itransient_outward_current_r_gate_r,                        STEWART_BATCH_SUPPORT_tau_r,    STEWART_BATCH_SUPPORT_gInf_r)

        #undef STEWART_RL

        for (int si = 0; si < nStates; ++si)
        {
            if (si == Ihyperpolarization_activated_current_y_gate_y         ||
                si == Irapid_time_dependent_potassium_current_Xr1_gate_Xr1  ||
                si == Irapid_time_dependent_potassium_current_Xr2_gate_Xr2  ||
                si == Islow_time_dependent_potassium_current_Xs_gate_Xs     ||
                si == Ifast_sodium_current_m_gate_m                         ||
                si == Ifast_sodium_current_h_gate_h                         ||
                si == Ifast_sodium_current_j_gate_j                         ||
                si == IL_type_Ca_current_d_gate_d                           ||
                si == IL_type_Ca_current_f_gate_f                           ||
                si == IL_type_Ca_current_f2_gate_f2                         ||
                si == IL_type_Ca_current_fCass_gate_fCass                   ||
                si == Itransient_outward_current_s_gate_s                   ||
                si == Itransient_outward_current_r_gate_r) continue;
            if (!solveVm && si == vmStateI) continue;
            const int idx = si*N + cellI;
            STATES[idx] += dt*RATES[idx];
        }
    }
}


void launchStewartRushLarsenStepKernel
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
    stewartRushLarsenStepKernel<<<nBlocks(N), blockSize>>>
    (
        d_STATES, d_RATES, d_SUPPORT,
        dt, N, nStates, solveVm, vmStateI
    );
    CUDA_LAUNCH_CHECK();
}

} // End namespace Foam

// ************************************************************************* //
