/*---------------------------------------------------------------------------*\
License
    This file is part of cardiacFoam.
\*---------------------------------------------------------------------------*/

#include <cuda_runtime.h>

#include "TNNP_2004Batch.H"

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

        TNNPComputeRatesBatch
        (
            t, CONSTANTS, N, cellI, cellI + 1,
            STATES, RATES, SUPPORT,
            solveVm, stimulus
        );
    }


    __global__ void tnnpFullBatchKernel
    (
        const double t,
        const double* __restrict__ CONSTANTS,
        const int N,
        const double* __restrict__ STATES,
        double* __restrict__ RATES,
        double* __restrict__ ALGEBRAIC,
        const bool solveVm,
        const StimulusProtocolPOD stimulus
    )
    {
        const int cellI = blockIdx.x*blockDim.x + threadIdx.x;
        if (cellI >= N)
        {
            return;
        }

        TNNPComputeRatesFullBatch
        (
            t, CONSTANTS, N, cellI, cellI + 1,
            STATES, RATES, ALGEBRAIC,
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


    __global__ void tnnpStabilizeKernel
    (
        double* __restrict__ STATES,
        const int N
    )
    {
        const int cellI = blockIdx.x*blockDim.x + threadIdx.x;
        if (cellI >= N)
        {
            return;
        }

        constexpr double small = 1.0e-15;

        #define TNNP_STATE(stateI) STATES[(stateI)*N + cellI]
        #define TNNP_CLAMP_RANGE(stateI, upper)                                \
            do                                                                 \
            {                                                                  \
                if (TNNP_STATE(stateI) < 0.0)                                  \
                {                                                              \
                    TNNP_STATE(stateI) = 0.0;                                  \
                }                                                              \
                else if (TNNP_STATE(stateI) > (upper))                        \
                {                                                              \
                    TNNP_STATE(stateI) = (upper);                              \
                }                                                              \
            } while (0)

        if (TNNP_STATE(K_i) < small)
        {
            TNNP_STATE(K_i) = small;
        }

        if (TNNP_STATE(Na_i) < small)
        {
            TNNP_STATE(Na_i) = small;
        }

        if (TNNP_STATE(Ca_i) < small)
        {
            TNNP_STATE(Ca_i) = small;
        }

        if (TNNP_STATE(Ca_SR) < small)
        {
            TNNP_STATE(Ca_SR) = small;
        }

        TNNP_CLAMP_RANGE(Xr1, 1.0);
        TNNP_CLAMP_RANGE(Xr2, 1.0);
        TNNP_CLAMP_RANGE(Xs, 1.0);
        TNNP_CLAMP_RANGE(m, 1.0);
        TNNP_CLAMP_RANGE(h, 1.0);
        TNNP_CLAMP_RANGE(j, 1.0);
        TNNP_CLAMP_RANGE(d, 1.0);
        TNNP_CLAMP_RANGE(f, 1.0);
        TNNP_CLAMP_RANGE(fCa, 1.0);
        TNNP_CLAMP_RANGE(s, 1.1);
        TNNP_CLAMP_RANGE(r, 1.0);
        TNNP_CLAMP_RANGE(g, 1.0);

        #undef TNNP_CLAMP_RANGE
        #undef TNNP_STATE
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
    (void)tissueFlag;
    tnnpBatchKernel<<<nBlocks(N), blockSize>>>
    (
        t, d_CONSTANTS, N,
        d_STATES, d_RATES, d_SUPPORT,
        solveVm, stimulus
    );
}


void launchTnnpFullBatchKernel
(
    double t,
    const double* d_CONSTANTS,
    int N,
    const double* d_STATES,
    double* d_RATES,
    double* d_ALGEBRAIC,
    int tissueFlag,
    bool solveVm,
    StimulusProtocolPOD stimulus
)
{
    (void)tissueFlag;
    tnnpFullBatchKernel<<<nBlocks(N), blockSize>>>
    (
        t, d_CONSTANTS, N,
        d_STATES, d_RATES, d_ALGEBRAIC,
        solveVm, stimulus
    );
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
}


void launchTnnpStabilizeKernel
(
    double* d_STATES,
    int N
)
{
    tnnpStabilizeKernel<<<nBlocks(N), blockSize>>>(d_STATES, N);
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
}

} // End namespace Foam

// ************************************************************************* //
