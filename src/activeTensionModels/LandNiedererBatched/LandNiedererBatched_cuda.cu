/*---------------------------------------------------------------------------*\
License
    This file is part of cardiacFoam.

    cardiacFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.
\*---------------------------------------------------------------------------*/

#include <cuda_runtime.h>
#include <cstdio>

#include "../LandNiederer/LandNiederer_2017.H"

#define CUDA_LAND_NIEDERER_CHECK()                                           \
    do                                                                        \
    {                                                                         \
        const cudaError_t error = cudaGetLastError();                         \
        if (error != cudaSuccess)                                             \
        {                                                                     \
            fprintf                                                           \
            (                                                                 \
                stderr,                                                       \
                "[cardiacFoam CUDA] kernel error at %s:%d: %s\\n",         \
                __FILE__, __LINE__, cudaGetErrorString(error)                \
            );                                                                \
            abort();                                                          \
        }                                                                     \
    } while (0)

namespace Foam
{
namespace
{

constexpr int landNiedererBlockSize = 256;

__global__ void landNiedererBatchedKernel
(
    const double* __restrict__ constants,
    const int nCells,
    const double* __restrict__ states,
    double* __restrict__ rates,
    double* __restrict__ algebraics
)
{
    const int cellI = blockIdx.x * blockDim.x + threadIdx.x;
    if (cellI >= nCells)
    {
        return;
    }

    double localStates[NUM_STATES];
    double localRates[NUM_STATES];
    double localAlgebraics[NUM_ALGEBRAIC];

    for (int stateI = 0; stateI < NUM_STATES; ++stateI)
    {
        localStates[stateI] = states[stateI * nCells + cellI];
    }
    for (int algebraicI = 0; algebraicI < NUM_ALGEBRAIC; ++algebraicI)
    {
        localAlgebraics[algebraicI] = algebraics[algebraicI * nCells + cellI];
    }

    // Convert lambda rate from s^-1 to ms^-1.
    localAlgebraics[AV_lambda_rate] *= 1.0e-3;

    LandNiederer2017computeVariables
    (
        0.0,
        const_cast<double*>(constants),
        localRates,
        localStates,
        localAlgebraics
    );

    for (int stateI = 0; stateI < NUM_STATES; ++stateI)
    {
        rates[stateI * nCells + cellI] = localRates[stateI] * 1000.0;
    }
    for (int algebraicI = 0; algebraicI < NUM_ALGEBRAIC; ++algebraicI)
    {
        algebraics[algebraicI * nCells + cellI] = localAlgebraics[algebraicI];
    }
}

} // End unnamed namespace


void launchLandNiedererBatchedKernel
(
    const double* deviceConstants,
    const int nCells,
    const double* deviceStates,
    double* deviceRates,
    double* deviceAlgebraics
)
{
    const int nBlocks =
        (nCells + landNiedererBlockSize - 1) / landNiedererBlockSize;
    landNiedererBatchedKernel<<<nBlocks, landNiedererBlockSize>>>
    (
        deviceConstants,
        nCells,
        deviceStates,
        deviceRates,
        deviceAlgebraics
    );
    CUDA_LAND_NIEDERER_CHECK();
}

} // End namespace Foam

// ************************************************************************* //
