/*---------------------------------------------------------------------------*\
License
    This file is part of cardiacFoam.

    cardiacFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

Description
    CUDA launch wrappers for the LandNiedererTWorld batched tension model.

Author
    Simao Nieto de Castro, UCD.
\*---------------------------------------------------------------------------*/

#include <cuda_runtime.h>
#include <cstdio>
#include "../LandNiedererTWorld/LandNiedererTWorld_2025.H"

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
        double localStates[NUM_STATES];
        double localRates[NUM_STATES];
        double localAlgebraics[NUM_ALGEBRAIC];

        // Gather from Global SoA
        for (int i = 0; i < NUM_STATES; ++i)
        {
            localStates[i] = STATES[i * N + cellI];
        }
        for (int i = 0; i < NUM_ALGEBRAIC; ++i)
        {
            localAlgebraics[i] = ALGEBRAIC[i * N + cellI];
        }

        // Convert lambda_rate from s^-1 to ms^-1 for the biophysics core
        localAlgebraics[AV_lambda_rate] *= 1e-3;

        // Compute using the standard biophysics core
        LandNiedererTWorld2025computeVariables
        (
            0.0,
            const_cast<double*>(CONSTANTS),
            localRates,
            localStates,
            localAlgebraics
        );

        // Convert computed rates from ms^-1 back to s^-1
        for (int i = 0; i < NUM_STATES; ++i)
        {
            localRates[i] *= 1000.0;
        }

        // Scatter to Global SoA
        for (int i = 0; i < NUM_STATES; ++i)
        {
            RATES[i * N + cellI] = localRates[i];
        }
        for (int i = 0; i < NUM_ALGEBRAIC; ++i)
        {
            ALGEBRAIC[i * N + cellI] = localAlgebraics[i];
        }
    }

    inline int nBlocks(const int N)
    {
        return (N + blockSize - 1) / blockSize;
    }
}

void launchLandNiedererTWorldBatchKernel
(
    const double* d_CONSTANTS,
    int N,
    const double* d_STATES,
    double* d_RATES,
    double* d_ALGEBRAIC
)
{
    // LandNiedererTWorld has NUM_STATES states and NUM_ALGEBRAIC algebraics.
    landNiedererBatchKernel<<<nBlocks(N), blockSize>>>
    (
        d_CONSTANTS,
        N, NUM_STATES, NUM_ALGEBRAIC,
        d_STATES,
        d_RATES,
        d_ALGEBRAIC
    );
    CUDA_LAUNCH_CHECK();
}

} // End namespace Foam
