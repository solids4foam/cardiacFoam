/*---------------------------------------------------------------------------*\
License
    This file is part of cardiacFoam.
\*---------------------------------------------------------------------------*/

#include <cuda_runtime.h>
#include "batchedIonicCoreCUDA.H"

namespace Foam
{
namespace
{
    constexpr int blockSize = 256;

    __global__ void batchedVmExtrapolantKernel
    (
        double* states,
        const double* vmStateStart,
        const double* vmStateRate,
        const int nCells,
        const int vmStateI,
        const double modelTime,
        const double modelStartTime
    )
    {
        const int cellI = blockIdx.x*blockDim.x + threadIdx.x;
        if (cellI >= nCells) return;

        states[static_cast<std::size_t>(vmStateI)*nCells + cellI] =
            vmStateStart[cellI]
          + vmStateRate[cellI]*(modelTime - modelStartTime);
    }
}

void launchBatchedVmExtrapolantKernel
(
    double* d_states,
    const double* d_vmStateStart,
    const double* d_vmStateRate,
    const std::size_t nCells,
    const std::size_t vmStateI,
    const double modelTime,
    const double modelStartTime
)
{
    batchedVmExtrapolantKernel<<<
        (static_cast<int>(nCells) + blockSize - 1)/blockSize,
        blockSize
    >>>
    (
        d_states,
        d_vmStateStart,
        d_vmStateRate,
        static_cast<int>(nCells),
        static_cast<int>(vmStateI),
        modelTime,
        modelStartTime
    );
    CARDIAC_CUDA_CHECK(cudaGetLastError());
}
}

// ************************************************************************* //
