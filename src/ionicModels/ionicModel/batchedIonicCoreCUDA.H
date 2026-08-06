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

Class
    Foam::batchedIonicCoreCUDA

Description
    Shared RAII helper for the device-side mirror of a batchedIonicCore. Owns
    five double* device buffers (states, rates, support, constants, Im) and
    the dirty flags that drive host<->device synchronization.

    Header-only and guarded by HAS_CUDA. When CUDA is not enabled at compile
    time, including this header is a no-op.

Author
    Simao Nieto de Castro, UCD.
\*---------------------------------------------------------------------------*/

#ifndef batchedIonicCoreCUDA_H
#define batchedIonicCoreCUDA_H

#ifdef HAS_CUDA

#include <cuda_runtime.h>
#include <cstddef>
#include "error.H"

namespace Foam
{

#define CARDIAC_CUDA_CHECK(expr)                                               \
    do {                                                                       \
        cudaError_t _err = (expr);                                             \
        if (_err != cudaSuccess) {                                             \
            FatalErrorInFunction                                               \
                << "CUDA error: " << cudaGetErrorString(_err)                  \
                << " at " << #expr                                             \
                << exit(FatalError);                                           \
        }                                                                      \
    } while (0)

class batchedIonicCoreCUDA
{
public:
    double* d_states    = nullptr;
    double* d_rates     = nullptr;
    double* d_support   = nullptr;
    double* d_constants = nullptr;
    double* d_cellConstants = nullptr;
    double* d_Im        = nullptr;

    bool hostDirty   = true;
    bool deviceDirty = false;
    bool allocated   = false;

    batchedIonicCoreCUDA() = default;

    batchedIonicCoreCUDA(const batchedIonicCoreCUDA&) = delete;
    batchedIonicCoreCUDA& operator=(const batchedIonicCoreCUDA&) = delete;

    ~batchedIonicCoreCUDA()
    {
        free();
    }

    void allocate
    (
        const std::size_t nCells,
        const std::size_t nStates,
        const std::size_t nSupport,
        const std::size_t nConstants
    )
    {
        if (allocated) return;
        const std::size_t Nbytes       = sizeof(double)*nStates*nCells;
        const std::size_t supportBytes = sizeof(double)*nSupport*nCells;
        const std::size_t constsBytes  = sizeof(double)*nConstants;
        const std::size_t cellConstsBytes = sizeof(double)*nConstants*nCells;
        const std::size_t imBytes      = sizeof(double)*nCells;

        CARDIAC_CUDA_CHECK(cudaMalloc(&d_states,    Nbytes));
        CARDIAC_CUDA_CHECK(cudaMalloc(&d_rates,     Nbytes));
        CARDIAC_CUDA_CHECK(cudaMalloc(&d_support,   supportBytes));
        CARDIAC_CUDA_CHECK(cudaMalloc(&d_constants, constsBytes));
        CARDIAC_CUDA_CHECK(cudaMalloc(&d_cellConstants, cellConstsBytes));
        CARDIAC_CUDA_CHECK(cudaMalloc(&d_Im,        imBytes));
        allocated = true;
    }

    void uploadConstants(const double* h_constants, const std::size_t nConstants)
    {
        CARDIAC_CUDA_CHECK(cudaMemcpy
        (
            d_constants, h_constants,
            sizeof(double)*nConstants, cudaMemcpyHostToDevice
        ));
    }

    void uploadCellConstants
    (
        const double* h_cellConstants,
        const std::size_t nCells,
        const std::size_t nConstants
    )
    {
        CARDIAC_CUDA_CHECK(cudaMemcpy
        (
            d_cellConstants, h_cellConstants,
            sizeof(double)*nCells*nConstants, cudaMemcpyHostToDevice
        ));
    }

    void syncStatesHostToDevice(const double* h_states, const std::size_t nStates, const std::size_t nCells)
    {
        CARDIAC_CUDA_CHECK(cudaMemcpy
        (
            d_states, h_states,
            sizeof(double)*nStates*nCells, cudaMemcpyHostToDevice
        ));
        hostDirty = false;
    }

    void syncStatesDeviceToHost(double* h_states, const std::size_t nStates, const std::size_t nCells)
    {
        CARDIAC_CUDA_CHECK(cudaMemcpy
        (
            h_states, d_states,
            sizeof(double)*nStates*nCells, cudaMemcpyDeviceToHost
        ));
        deviceDirty = false;
    }

    void syncSupportDeviceToHost(double* h_support, const std::size_t nSupport, const std::size_t nCells)
    {
        if (!h_support || nSupport == 0) return;
        CARDIAC_CUDA_CHECK(cudaMemcpy
        (
            h_support, d_support,
            sizeof(double)*nSupport*nCells, cudaMemcpyDeviceToHost
        ));
    }

    void downloadIm(double* h_Im, const std::size_t nCells)
    {
        CARDIAC_CUDA_CHECK(cudaMemcpy
        (
            h_Im, d_Im,
            sizeof(double)*nCells, cudaMemcpyDeviceToHost
        ));
    }

    void free()
    {
        if (!allocated) return;
        auto safeFree = [](double*& p)
        {
            if (p) { cudaFree(p); p = nullptr; }
        };
        safeFree(d_states);
        safeFree(d_rates);
        safeFree(d_support);
        safeFree(d_constants);
        safeFree(d_cellConstants);
        safeFree(d_Im);
        allocated = false;
    }
};

} // End namespace Foam

#endif // HAS_CUDA

#endif

// ************************************************************************* //
