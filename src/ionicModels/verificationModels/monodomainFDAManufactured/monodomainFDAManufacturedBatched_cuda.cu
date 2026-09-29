/*---------------------------------------------------------------------------*\
License
    This file is part of cardiacFoam.
\*---------------------------------------------------------------------------*/

#include <cuda_runtime.h>
#include <cstdio>
#include <cstdlib>

#include "monodomainFDAManufacturedBatch.H"

namespace Foam
{
namespace
{
    constexpr int blockSize = 256;

    inline int nBlocks(const int N)
    {
        return (N + blockSize - 1)/blockSize;
    }

    __global__ void monodomainFDAManufacturedBatchKernel
    (
        const double* __restrict__ CONSTANTS,
        const int N,
        const double* __restrict__ STATES,
        double* __restrict__ RATES,
        double* __restrict__ SUPPORT
    )
    {
        const int cellI = blockIdx.x*blockDim.x + threadIdx.x;
        if (cellI >= N)
        {
            return;
        }

        monodomainFDAManufacturedComputeBatch
        (
            N, cellI, cellI + 1, CONSTANTS,
            STATES, RATES, SUPPORT
        );
    }


    __global__ void monodomainFDAManufacturedEulerStepKernel
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


    __global__ void monodomainFDAManufacturedScaleIonKernel
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


    __global__ void monodomainFDAManufacturedSetVmKernel
    (
        double* __restrict__ STATES,
        const double* __restrict__ VM_START,
        const double* __restrict__ VM_RATE,
        const double modelTimeOffset,
        const int N,
        const int vmStateI
    )
    {
        const int cellI = blockIdx.x*blockDim.x + threadIdx.x;
        if (cellI >= N)
        {
            return;
        }

        STATES[vmStateI*N + cellI] =
            VM_START[cellI] + modelTimeOffset*VM_RATE[cellI];
    }
}


void launchMonodomainFDAManufacturedBatchKernel
(
    const double* d_CONSTANTS,
    const int N,
    const double* d_STATES,
    double* d_RATES,
    double* d_SUPPORT
)
{
    monodomainFDAManufacturedBatchKernel<<<nBlocks(N), blockSize>>>
    (
        d_CONSTANTS, N, d_STATES, d_RATES, d_SUPPORT
    );
    const cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess)
    {
        fprintf
        (
            stderr,
            "[cardiacFoam CUDA] manufactured batch kernel failed: %s\n",
            cudaGetErrorString(err)
        );
        std::abort();
    }
}


void launchMonodomainFDAManufacturedEulerStepKernel
(
    double* d_STATES,
    const double* d_RATES,
    const double dt,
    const int N,
    const int nStates,
    const bool solveVm,
    const int vmStateI
)
{
    monodomainFDAManufacturedEulerStepKernel<<<nBlocks(N), blockSize>>>
    (
        d_STATES, d_RATES, dt, N, nStates, solveVm, vmStateI
    );
    const cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess)
    {
        fprintf
        (
            stderr,
            "[cardiacFoam CUDA] manufactured Euler kernel failed: %s\n",
            cudaGetErrorString(err)
        );
        std::abort();
    }
}


void launchMonodomainFDAManufacturedScaleIonKernel
(
    const double* d_SUPPORT,
    double* d_Im,
    const double scale,
    const int N,
    const int IionSlot
)
{
    monodomainFDAManufacturedScaleIonKernel<<<nBlocks(N), blockSize>>>
    (
        d_SUPPORT, d_Im, scale, N, IionSlot
    );
    const cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess)
    {
        fprintf
        (
            stderr,
            "[cardiacFoam CUDA] manufactured current kernel failed: %s\n",
            cudaGetErrorString(err)
        );
        std::abort();
    }
}


void launchMonodomainFDAManufacturedSetVmKernel
(
    double* d_STATES,
    const double* d_VM_START,
    const double* d_VM_RATE,
    const double modelTimeOffset,
    const int N,
    const int vmStateI
)
{
    monodomainFDAManufacturedSetVmKernel<<<nBlocks(N), blockSize>>>
    (
        d_STATES, d_VM_START, d_VM_RATE, modelTimeOffset, N, vmStateI
    );
    const cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess)
    {
        fprintf
        (
            stderr,
            "[cardiacFoam CUDA] manufactured Vm extrapolation kernel failed: %s\n",
            cudaGetErrorString(err)
        );
        std::abort();
    }
}

} // End namespace Foam

// ************************************************************************* //
