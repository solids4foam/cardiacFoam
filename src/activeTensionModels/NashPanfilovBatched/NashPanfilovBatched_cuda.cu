/*---------------------------------------------------------------------------*\
License
    This file is part of cardiacFoam.

    cardiacFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

Description
    CUDA launch wrappers for NashPanfilov batched tension model.

Author
    Simao Nieto de Castro, UCD.
\*---------------------------------------------------------------------------*/

#include <cuda_runtime.h>
#include <cstdio>
#include "NashPanfilovBatch.H"

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

    __global__ void nashPanfilovBatchKernel
    (
        const double* __restrict__ driveSignals,
        const double* __restrict__ CONSTANTS,
        const int N,
        const double* __restrict__ STATES,
        double* __restrict__ RATES,
        double* __restrict__ ALGEBRAIC
    )
    {
        const int cellI = blockIdx.x*blockDim.x + threadIdx.x;
        if (cellI >= N) return;

        NashPanfilovComputeVariablesBatch
        (
            driveSignals[cellI],
            CONSTANTS,
            N, cellI,
            STATES, RATES, ALGEBRAIC
        );
    }

    inline int nBlocks(const int N)
    {
        return (N + blockSize - 1) / blockSize;
    }
}

void launchNashPanfilovBatchKernel
(
    const double* d_driveSignals,
    const double* d_CONSTANTS,
    int N,
    const double* d_STATES,
    double* d_RATES,
    double* d_ALGEBRAIC
)
{
    nashPanfilovBatchKernel<<<nBlocks(N), blockSize>>>
    (
        d_driveSignals,
        d_CONSTANTS,
        N,
        d_STATES,
        d_RATES,
        d_ALGEBRAIC
    );
    CUDA_LAUNCH_CHECK();
}

} // End namespace Foam
