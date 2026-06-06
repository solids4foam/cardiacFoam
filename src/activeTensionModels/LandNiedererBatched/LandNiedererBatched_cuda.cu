/*---------------------------------------------------------------------------*\
License
    This file is part of cardiacFoam.

    cardiacFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

Description
    CUDA launch wrappers for LandNiederer batched tension model.

Author
    Simao Nieto de Castro, UCD.
\*---------------------------------------------------------------------------*/

#include <cuda_runtime.h>
#include <cstdio>
#include "../LandNiederer/LandNiederer_2017.H"

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

namespace Foam
{
namespace
{
    constexpr int blockSize = 256;

    __global__ void landNiedererBatchKernel
    (
        const double* __restrict__ CONSTANTS,
        const int N,
        const int nStates,
        const int nAlgebraics,
        const double* __restrict__ STATES,
        double* __restrict__ RATES,
        double* __restrict__ ALGEBRAIC
    )
    {
        const int cellI = blockIdx.x*blockDim.x + threadIdx.x;
        if (cellI >= N) return;

        // AoS Thread-Local Registers to interface with standard ODE function
        double localStates[6];
        double localRates[6];
        double localAlgebraics[2];

        // Gather from Global SoA
        for (int i = 0; i < 6; ++i)
        {
            localStates[i] = STATES[i * N + cellI];
        }

        // Compute using the standard biophysics core
        LandNiederer2017computeVariables
        (
            0.0,
            const_cast<double*>(CONSTANTS),
            localRates,
            localStates,
            localAlgebraics
        );

        // Scatter to Global SoA
        for (int i = 0; i < 6; ++i)
        {
            RATES[i * N + cellI] = localRates[i];
        }
        for (int i = 0; i < 2; ++i)
        {
            ALGEBRAIC[i * N + cellI] = localAlgebraics[i];
        }
    }

    inline int nBlocks(const int N)
    {
        return (N + blockSize - 1) / blockSize;
    }
}

void launchLandNiedererBatchKernel
(
    const double* d_CONSTANTS,
    int N,
    const double* d_STATES,
    double* d_RATES,
    double* d_ALGEBRAIC
)
{
    // LandNiederer has 6 states and 2 algebraics.
    landNiedererBatchKernel<<<nBlocks(N), blockSize>>>
    (
        d_CONSTANTS,
        N, 6, 2,
        d_STATES,
        d_RATES,
        d_ALGEBRAIC
    );
    CUDA_LAUNCH_CHECK();
}

} // End namespace Foam
