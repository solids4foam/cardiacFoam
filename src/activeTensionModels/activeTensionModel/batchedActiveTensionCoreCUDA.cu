/*---------------------------------------------------------------------------*\\
License
    This file is part of cardiacFoam.
\\*---------------------------------------------------------------------------*/

#include <cuda_runtime.h>
#include "batchedActiveTensionCoreCUDA.H"

namespace Foam
{
namespace
{
constexpr int blockSize = 256;

__global__ void batchedActiveTensionEulerKernel
(
    double* states,
    const double* rates,
    const int nCells,
    const int nStates,
    const double dtModel
)
{
    const int cellI = blockIdx.x*blockDim.x + threadIdx.x;
    if (cellI >= nCells) return;

    for (int stateI = 0; stateI < nStates; ++stateI)
    {
        const std::size_t offset = static_cast<std::size_t>(stateI)*nCells + cellI;
        states[offset] += dtModel*rates[offset];
    }
}
}

void launchBatchedActiveTensionEulerKernel
(
    double* d_states,
    const double* d_rates,
    const std::size_t nCells,
    const std::size_t nStates,
    const double dtModel
)
{
    batchedActiveTensionEulerKernel<<<
        (static_cast<int>(nCells) + blockSize - 1)/blockSize, blockSize
    >>>
    (
        d_states, d_rates, static_cast<int>(nCells), static_cast<int>(nStates),
        dtModel
    );
    CARDIAC_ACTIVE_TENSION_CUDA_CHECK(cudaGetLastError());
}
} // End namespace Foam

