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

#include <cuda_runtime.h>
#include <cstdio>

#include "TWorld_2025Batch.H"

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

    __global__ void tWorldBatchKernel
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

        TWorldComputeVariablesBatch
        (
            t, cellConstants, N, cellI, cellI + 1,
            STATES, RATES, SUPPORT,
            solveVm, stimulus
        );
    }


    __global__ void tWorldEulerStepKernel
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


    __global__ void tWorldScaleIonKernel
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
)
{
    (void)tissueFlag;
    tWorldBatchKernel<<<nBlocks(N), blockSize>>>
    (
        t, d_CONSTANTS, d_CELL_CONSTANTS, useCellConstants, N,
        d_STATES, d_RATES, d_SUPPORT,
        solveVm, stimulus
    );
    CUDA_LAUNCH_CHECK();
}


void launchTWorldEulerStepKernel
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
    tWorldEulerStepKernel<<<nBlocks(N), blockSize>>>
    (
        d_STATES, d_RATES, dt,
        N, nStates, solveVm, vmStateI
    );
    CUDA_LAUNCH_CHECK();
}


void launchTWorldScaleIonKernel
(
    const double* d_SUPPORT,
    double* d_Im,
    double scale,
    int N,
    int IionSlot
)
{
    tWorldScaleIonKernel<<<nBlocks(N), blockSize>>>
    (
        d_SUPPORT, d_Im, scale, N, IionSlot
    );
    CUDA_LAUNCH_CHECK();
}


namespace
{
    __global__ void tWorldRushLarsenStepKernel
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

        #define TW_RL(si, tSlot, iSlot)                                       \
        {                                                                      \
            const double _x   = STATES[(si)*N + cellI];                       \
            const double _inf = SUPPORT[(iSlot)*N + cellI];                   \
            const double _tau = SUPPORT[(tSlot)*N + cellI];                   \
            STATES[(si)*N + cellI] = _inf + (_x - _inf)*exp(-dt/_tau);        \
        }

        TW_RL(camk_f_ICaL, TWORLD_BATCH_SUPPORT_tau_camk_ICaL, TWORLD_BATCH_SUPPORT_gInf_camk_ICaL)
        TW_RL(camk_f_PLB,  TWORLD_BATCH_SUPPORT_tau_camk_PLB,  TWORLD_BATCH_SUPPORT_gInf_camk_PLB)
        TW_RL(camk_f_RyR,  TWORLD_BATCH_SUPPORT_tau_camk_RyR,  TWORLD_BATCH_SUPPORT_gInf_camk_RyR)
        TW_RL(m,    TWORLD_BATCH_SUPPORT_tau_m,    TWORLD_BATCH_SUPPORT_gInf_m)
        TW_RL(m_P,  TWORLD_BATCH_SUPPORT_tau_m_P,  TWORLD_BATCH_SUPPORT_gInf_m_P)
        TW_RL(h,    TWORLD_BATCH_SUPPORT_tau_h,    TWORLD_BATCH_SUPPORT_gInf_h)
        TW_RL(h_P,  TWORLD_BATCH_SUPPORT_tau_h_P,  TWORLD_BATCH_SUPPORT_gInf_h_P)
        TW_RL(hp,   TWORLD_BATCH_SUPPORT_tau_hp,   TWORLD_BATCH_SUPPORT_gInf_hp)
        TW_RL(hp_P, TWORLD_BATCH_SUPPORT_tau_hp_P, TWORLD_BATCH_SUPPORT_gInf_hp_P)
        TW_RL(j,    TWORLD_BATCH_SUPPORT_tau_j,    TWORLD_BATCH_SUPPORT_gInf_j)
        TW_RL(j_P,  TWORLD_BATCH_SUPPORT_tau_j_P,  TWORLD_BATCH_SUPPORT_gInf_j_P)
        TW_RL(jp,   TWORLD_BATCH_SUPPORT_tau_jp,   TWORLD_BATCH_SUPPORT_gInf_jp)
        TW_RL(jp_P, TWORLD_BATCH_SUPPORT_tau_jp_P, TWORLD_BATCH_SUPPORT_gInf_jp_P)
        TW_RL(mL,   TWORLD_BATCH_SUPPORT_tau_mL,   TWORLD_BATCH_SUPPORT_gInf_mL)
        TW_RL(hL,   TWORLD_BATCH_SUPPORT_tau_hL,   TWORLD_BATCH_SUPPORT_gInf_hL)
        TW_RL(hLp,  TWORLD_BATCH_SUPPORT_tau_hLp,  TWORLD_BATCH_SUPPORT_gInf_hLp)
        TW_RL(d,      TWORLD_BATCH_SUPPORT_tau_d,      TWORLD_BATCH_SUPPORT_gInf_d)
        TW_RL(ff,     TWORLD_BATCH_SUPPORT_tau_ff,     TWORLD_BATCH_SUPPORT_gInf_ff)
        TW_RL(fs,     TWORLD_BATCH_SUPPORT_tau_fs,     TWORLD_BATCH_SUPPORT_gInf_fs)
        TW_RL(fcaf,   TWORLD_BATCH_SUPPORT_tau_fcaf,   TWORLD_BATCH_SUPPORT_gInf_fcaf)
        TW_RL(fcas,   TWORLD_BATCH_SUPPORT_tau_fcas,   TWORLD_BATCH_SUPPORT_gInf_fcas)
        TW_RL(jca,    TWORLD_BATCH_SUPPORT_tau_jca,    TWORLD_BATCH_SUPPORT_gInf_jca)
        TW_RL(ffp,    TWORLD_BATCH_SUPPORT_tau_ffp,    TWORLD_BATCH_SUPPORT_gInf_ffp)
        TW_RL(fcafp,  TWORLD_BATCH_SUPPORT_tau_fcafp,  TWORLD_BATCH_SUPPORT_gInf_fcafp)
        TW_RL(d_P,    TWORLD_BATCH_SUPPORT_tau_d_P,    TWORLD_BATCH_SUPPORT_gInf_d_P)
        TW_RL(ff_P,   TWORLD_BATCH_SUPPORT_tau_ff_P,   TWORLD_BATCH_SUPPORT_gInf_ff_P)
        TW_RL(fs_P,   TWORLD_BATCH_SUPPORT_tau_fs_P,   TWORLD_BATCH_SUPPORT_gInf_fs_P)
        TW_RL(fcaf_P, TWORLD_BATCH_SUPPORT_tau_fcaf_P, TWORLD_BATCH_SUPPORT_gInf_fcaf_P)
        TW_RL(fcas_P, TWORLD_BATCH_SUPPORT_tau_fcas_P, TWORLD_BATCH_SUPPORT_gInf_fcas_P)
        TW_RL(fBPf,   TWORLD_BATCH_SUPPORT_tau_fBPf,   TWORLD_BATCH_SUPPORT_gInf_fBPf)
        TW_RL(fcaBPf, TWORLD_BATCH_SUPPORT_tau_fcaBPf, TWORLD_BATCH_SUPPORT_gInf_fcaBPf)
        TW_RL(xtos,   TWORLD_BATCH_SUPPORT_tau_xtos,   TWORLD_BATCH_SUPPORT_gInf_xtos)
        TW_RL(xtos_p, TWORLD_BATCH_SUPPORT_tau_xtos_p, TWORLD_BATCH_SUPPORT_gInf_xtos_p)
        TW_RL(xtof,   TWORLD_BATCH_SUPPORT_tau_xtof,   TWORLD_BATCH_SUPPORT_gInf_xtof)
        TW_RL(xtof_p, TWORLD_BATCH_SUPPORT_tau_xtof_p, TWORLD_BATCH_SUPPORT_gInf_xtof_p)
        TW_RL(ytos,   TWORLD_BATCH_SUPPORT_tau_ytos,   TWORLD_BATCH_SUPPORT_gInf_ytos)
        TW_RL(ytos_p, TWORLD_BATCH_SUPPORT_tau_ytos_p, TWORLD_BATCH_SUPPORT_gInf_ytos_p)
        TW_RL(ytof,   TWORLD_BATCH_SUPPORT_tau_ytof,   TWORLD_BATCH_SUPPORT_gInf_ytof)
        TW_RL(ytof_p, TWORLD_BATCH_SUPPORT_tau_ytof_p, TWORLD_BATCH_SUPPORT_gInf_ytof_p)
        TW_RL(xs_junc, TWORLD_BATCH_SUPPORT_tau_xs_junc, TWORLD_BATCH_SUPPORT_gInf_xs_junc)
        TW_RL(xs_sl,   TWORLD_BATCH_SUPPORT_tau_xs_sl,   TWORLD_BATCH_SUPPORT_gInf_xs_sl)
        TW_RL(jrel_icaldep_act, TWORLD_BATCH_SUPPORT_tau_jrel_act, TWORLD_BATCH_SUPPORT_gInf_jrel_act)
        TW_RL(jrel_icaldep_f1,  TWORLD_BATCH_SUPPORT_tau_jrel_f1,  TWORLD_BATCH_SUPPORT_gInf_jrel_f1)
        TW_RL(jrel_icaldep_f2,  TWORLD_BATCH_SUPPORT_tau_jrel_f2,  TWORLD_BATCH_SUPPORT_gInf_jrel_f2)

        #undef TW_RL

        for (int si = 0; si < nStates; ++si)
        {
            if (si == camk_f_ICaL || si == camk_f_PLB  || si == camk_f_RyR  ||
                si == m            || si == m_P          || si == h           ||
                si == h_P          || si == hp           || si == hp_P        ||
                si == j            || si == j_P          || si == jp          ||
                si == jp_P         || si == mL           || si == hL          ||
                si == hLp          || si == d            || si == ff          ||
                si == fs           || si == fcaf         || si == fcas        ||
                si == jca          || si == ffp          || si == fcafp       ||
                si == d_P          || si == ff_P         || si == fs_P        ||
                si == fcaf_P       || si == fcas_P       || si == fBPf        ||
                si == fcaBPf       || si == xtos         || si == xtos_p      ||
                si == xtof         || si == xtof_p       || si == ytos        ||
                si == ytos_p       || si == ytof         || si == ytof_p      ||
                si == xs_junc      || si == xs_sl        ||
                si == jrel_icaldep_act || si == jrel_icaldep_f1 || si == jrel_icaldep_f2)
            {
                continue;
            }
            if (!solveVm && si == vmStateI) continue;
            const int idx = si*N + cellI;
            STATES[idx] += dt*RATES[idx];
        }
    }
}


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
)
{
    tWorldRushLarsenStepKernel<<<nBlocks(N), blockSize>>>
    (
        d_STATES, d_RATES, d_SUPPORT,
        dt, N, nStates, solveVm, vmStateI
    );
    CUDA_LAUNCH_CHECK();
}

} // End namespace Foam

// ************************************************************************* //
