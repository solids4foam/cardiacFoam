/*---------------------------------------------------------------------------*\
License
    This file is part of cardiacFoam.
\*---------------------------------------------------------------------------*/

#include <cuda_runtime.h>
#include <cstdio>
#include <cstdlib>

#include "Fabbri_2017Batch.H"

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

    __global__ void fabbriBatchKernel
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

        FabbriComputeVariablesBatch
        (
            t, CONSTANTS, N, cellI, cellI + 1,
            STATES, RATES, SUPPORT,
            solveVm, stimulus
        );
    }


    __global__ void fabbriEulerStepKernel
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


    __global__ void fabbriScaleIonKernel
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


void launchFabbriBatchKernel
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
    fabbriBatchKernel<<<nBlocks(N), blockSize>>>
    (
        t, d_CONSTANTS, N,
        d_STATES, d_RATES, d_SUPPORT,
        solveVm, stimulus
    );
    CUDA_LAUNCH_CHECK();
}


void launchFabbriEulerStepKernel
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
    fabbriEulerStepKernel<<<nBlocks(N), blockSize>>>
    (
        d_STATES, d_RATES, dt,
        N, nStates, solveVm, vmStateI
    );
    CUDA_LAUNCH_CHECK();
}


void launchFabbriScaleIonKernel
(
    const double* d_SUPPORT,
    double* d_Im,
    double scale,
    int N,
    int IionSlot
)
{
    fabbriScaleIonKernel<<<nBlocks(N), blockSize>>>
    (
        d_SUPPORT, d_Im, scale, N, IionSlot
    );
    CUDA_LAUNCH_CHECK();
}


namespace
{
    __global__ void fabbriRushLarsenStepKernel
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

        #define FABBRI_RL(si, tSlot, iSlot)                                  \
        {                                                                     \
            const double _x   = STATES[(si)*N + cellI];                      \
            const double _inf = SUPPORT[(iSlot)*N + cellI];                  \
            const double _tau = SUPPORT[(tSlot)*N + cellI];                  \
            STATES[(si)*N + cellI] = _inf + (_x - _inf)*::exp(-dt/_tau);       \
        }

        FABBRI_RL(If_y_gate_y,          FABBRI_BATCH_SUPPORT_tau_y,       FABBRI_BATCH_SUPPORT_gInf_y)
        FABBRI_RL(INa_m_gate_m,         FABBRI_BATCH_SUPPORT_tau_m,       FABBRI_BATCH_SUPPORT_gInf_m)
        FABBRI_RL(INa_h_gate_h,         FABBRI_BATCH_SUPPORT_tau_h,       FABBRI_BATCH_SUPPORT_gInf_h)
        FABBRI_RL(ICaL_dL_gate_dL,      FABBRI_BATCH_SUPPORT_tau_dL,      FABBRI_BATCH_SUPPORT_gInf_dL)
        FABBRI_RL(ICaL_fL_gate_fL,      FABBRI_BATCH_SUPPORT_tau_fL,      FABBRI_BATCH_SUPPORT_gInf_fL)
        FABBRI_RL(ICaL_fCa_gate_fCa,    FABBRI_BATCH_SUPPORT_tau_fCa,     FABBRI_BATCH_SUPPORT_gInf_fCa)
        FABBRI_RL(ICaT_dT_gate_dT,      FABBRI_BATCH_SUPPORT_tau_dT,      FABBRI_BATCH_SUPPORT_gInf_dT)
        FABBRI_RL(ICaT_fT_gate_fT,      FABBRI_BATCH_SUPPORT_tau_fT,      FABBRI_BATCH_SUPPORT_gInf_fT)
        FABBRI_RL(IKur_rKur_gate_r_Kur, FABBRI_BATCH_SUPPORT_tau_r_Kur,   FABBRI_BATCH_SUPPORT_gInf_r_Kur)
        FABBRI_RL(IKur_sKur_gate_s_Kur, FABBRI_BATCH_SUPPORT_tau_s_Kur,   FABBRI_BATCH_SUPPORT_gInf_s_Kur)
        FABBRI_RL(Ito_q_gate_q,         FABBRI_BATCH_SUPPORT_tau_q,       FABBRI_BATCH_SUPPORT_gInf_q)
        FABBRI_RL(Ito_r_gate_r,         FABBRI_BATCH_SUPPORT_tau_r,       FABBRI_BATCH_SUPPORT_gInf_r)
        FABBRI_RL(IKr_pa_gate_paS,      FABBRI_BATCH_SUPPORT_tau_paS,     FABBRI_BATCH_SUPPORT_gInf_paS)
        FABBRI_RL(IKr_pa_gate_paF,      FABBRI_BATCH_SUPPORT_tau_paF,     FABBRI_BATCH_SUPPORT_gInf_paF)
        FABBRI_RL(IKr_pi_gate_piy,      FABBRI_BATCH_SUPPORT_tau_pi,      FABBRI_BATCH_SUPPORT_gInf_pi)
        FABBRI_RL(IKs_n_gate_n,         FABBRI_BATCH_SUPPORT_tau_n,       FABBRI_BATCH_SUPPORT_gInf_n)
        FABBRI_RL(IKACh_a_gate_a,       FABBRI_BATCH_SUPPORT_tau_a,       FABBRI_BATCH_SUPPORT_gInf_a)

        #undef FABBRI_RL

        for (int si = 0; si < nStates; ++si)
        {
            if (si == If_y_gate_y          || si == INa_m_gate_m         ||
                si == INa_h_gate_h         || si == ICaL_dL_gate_dL      ||
                si == ICaL_fL_gate_fL      || si == ICaL_fCa_gate_fCa    ||
                si == ICaT_dT_gate_dT      || si == ICaT_fT_gate_fT      ||
                si == IKur_rKur_gate_r_Kur || si == IKur_sKur_gate_s_Kur ||
                si == Ito_q_gate_q         || si == Ito_r_gate_r         ||
                si == IKr_pa_gate_paS      || si == IKr_pa_gate_paF      ||
                si == IKr_pi_gate_piy      || si == IKs_n_gate_n         ||
                si == IKACh_a_gate_a) continue;
            if (!solveVm && si == vmStateI) continue;
            const int idx = si*N + cellI;
            STATES[idx] += dt*RATES[idx];
        }
    }
}


void launchFabbriRushLarsenStepKernel
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
    fabbriRushLarsenStepKernel<<<nBlocks(N), blockSize>>>
    (
        d_STATES, d_RATES, d_SUPPORT,
        dt, N, nStates, solveVm, vmStateI
    );
    CUDA_LAUNCH_CHECK();
}

} // End namespace Foam

// ************************************************************************* //
