/*---------------------------------------------------------------------------*\
License
    This file is part of cardiacFoam.
\*---------------------------------------------------------------------------*/

#include <cuda_runtime.h>

#include "BuenoOrovio_2008Batch.H"

namespace Foam
{
namespace
{
    constexpr int blockSize = 256;

    __global__ void buenoBatchKernel
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

        BuenoOrovioComputeVariablesBatch
        (
            t, CONSTANTS, N, cellI, cellI + 1,
            STATES, RATES, SUPPORT,
            solveVm, stimulus
        );
    }


    __global__ void buenoEulerStepKernel
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


    __global__ void buenoScaleIonKernel
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


void launchBuenoBatchKernel
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
    buenoBatchKernel<<<nBlocks(N), blockSize>>>
    (
        t, d_CONSTANTS, N,
        d_STATES, d_RATES, d_SUPPORT,
        solveVm, stimulus
    );
}


void launchEulerStepKernel
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
    buenoEulerStepKernel<<<nBlocks(N), blockSize>>>
    (
        d_STATES, d_RATES, dt,
        N, nStates, solveVm, vmStateI
    );
}


void launchScaleIonKernel
(
    const double* d_SUPPORT,
    double* d_Im,
    double scale,
    int N,
    int IionSlot
)
{
    buenoScaleIonKernel<<<nBlocks(N), blockSize>>>
    (
        d_SUPPORT, d_Im, scale, N, IionSlot
    );
}

} // End namespace Foam

// ************************************************************************* //
