/*---------------------------------------------------------------------------*\\
License
    This file is part of cardiacFoam.
\\*---------------------------------------------------------------------------*/

#ifndef batchedActiveTensionCoreCUDA_H
#define batchedActiveTensionCoreCUDA_H

#ifdef HAS_CUDA

#include <cuda_runtime.h>
#include <cstddef>
#include "error.H"

namespace Foam
{

#define CARDIAC_ACTIVE_TENSION_CUDA_CHECK(expr)                               \
    do {                                                                       \
        const cudaError_t _error = (expr);                                     \
        if (_error != cudaSuccess)                                             \
        {                                                                      \
            FatalErrorInFunction                                               \
                << "CUDA error: " << cudaGetErrorString(_error)               \
                << " at " << #expr << exit(FatalError);                       \
        }                                                                      \
    } while (0)

// Explicit Euler update for the state-major (SoA) batched active-tension
// buffers.  The model-specific wrappers only evaluate the ODE right-hand side.
void launchBatchedActiveTensionEulerKernel
(
    double* d_states,
    const double* d_rates,
    std::size_t nCells,
    std::size_t nStates,
    double dtModel
);

class batchedActiveTensionCoreCUDA
{
public:
    double* d_states = nullptr;
    double* d_rates = nullptr;
    double* d_algebraics = nullptr;
    double* d_constants = nullptr;
    double* d_driveSignals = nullptr;
    bool allocated = false;
    bool constantsUploaded = false;
    bool stateResident = false;

    batchedActiveTensionCoreCUDA() = default;
    batchedActiveTensionCoreCUDA(const batchedActiveTensionCoreCUDA&) = delete;
    batchedActiveTensionCoreCUDA& operator=(const batchedActiveTensionCoreCUDA&) = delete;

    ~batchedActiveTensionCoreCUDA()
    {
        free();
    }

    void allocate
    (
        const std::size_t nCells,
        const std::size_t nStates,
        const std::size_t nAlgebraics,
        const std::size_t nConstants
    )
    {
        if (allocated) return;

        CARDIAC_ACTIVE_TENSION_CUDA_CHECK
        (
            cudaMalloc(&d_states, sizeof(double)*nCells*nStates)
        );
        CARDIAC_ACTIVE_TENSION_CUDA_CHECK
        (
            cudaMalloc(&d_rates, sizeof(double)*nCells*nStates)
        );
        CARDIAC_ACTIVE_TENSION_CUDA_CHECK
        (
            cudaMalloc(&d_algebraics, sizeof(double)*nCells*nAlgebraics)
        );
        CARDIAC_ACTIVE_TENSION_CUDA_CHECK
        (
            cudaMalloc(&d_constants, sizeof(double)*nConstants)
        );
        CARDIAC_ACTIVE_TENSION_CUDA_CHECK
        (
            cudaMalloc(&d_driveSignals, sizeof(double)*nCells)
        );
        allocated = true;
    }

    void uploadConstants(const double* host, const std::size_t nConstants)
    {
        CARDIAC_ACTIVE_TENSION_CUDA_CHECK
        (
            cudaMemcpy
            (
                d_constants, host, sizeof(double)*nConstants,
                cudaMemcpyHostToDevice
            )
        );
        constantsUploaded = true;
    }

    void uploadInitialFields
    (
        const double* states,
        const double* algebraics,
        const double* driveSignals,
        const std::size_t nCells,
        const std::size_t nStates,
        const std::size_t nAlgebraics
    )
    {
        CARDIAC_ACTIVE_TENSION_CUDA_CHECK
        (
            cudaMemcpy(d_states, states, sizeof(double)*nCells*nStates,
                       cudaMemcpyHostToDevice)
        );
        CARDIAC_ACTIVE_TENSION_CUDA_CHECK
        (
            cudaMemcpy(d_algebraics, algebraics,
                       sizeof(double)*nCells*nAlgebraics,
                       cudaMemcpyHostToDevice)
        );
        CARDIAC_ACTIVE_TENSION_CUDA_CHECK
        (
            cudaMemcpy(d_driveSignals, driveSignals, sizeof(double)*nCells,
                       cudaMemcpyHostToDevice)
        );
    }

    void uploadDriveSignals(const double* driveSignals, const std::size_t nCells)
    {
        CARDIAC_ACTIVE_TENSION_CUDA_CHECK
        (
            cudaMemcpy(d_driveSignals, driveSignals, sizeof(double)*nCells,
                       cudaMemcpyHostToDevice)
        );
    }

    void uploadAlgebraicRow
    (
        const double* host,
        const std::size_t algebraicI,
        const std::size_t nCells
    )
    {
        CARDIAC_ACTIVE_TENSION_CUDA_CHECK
        (
            cudaMemcpy
            (
                d_algebraics + algebraicI*nCells,
                host,
                sizeof(double)*nCells,
                cudaMemcpyHostToDevice
            )
        );
    }

    void downloadOutputs
    (
        double* states,
        double* rates,
        double* algebraics,
        const std::size_t nCells,
        const std::size_t nStates,
        const std::size_t nAlgebraics
    ) const
    {
        CARDIAC_ACTIVE_TENSION_CUDA_CHECK(cudaDeviceSynchronize());
        CARDIAC_ACTIVE_TENSION_CUDA_CHECK
        (
            cudaMemcpy(states, d_states, sizeof(double)*nCells*nStates,
                       cudaMemcpyDeviceToHost)
        );
        CARDIAC_ACTIVE_TENSION_CUDA_CHECK
        (
            cudaMemcpy(rates, d_rates, sizeof(double)*nCells*nStates,
                       cudaMemcpyDeviceToHost)
        );
        CARDIAC_ACTIVE_TENSION_CUDA_CHECK
        (
            cudaMemcpy(algebraics, d_algebraics,
                       sizeof(double)*nCells*nAlgebraics,
                       cudaMemcpyDeviceToHost)
        );
    }

    void downloadStateRow
    (
        double* host,
        const std::size_t stateI,
        const std::size_t nCells
    ) const
    {
        CARDIAC_ACTIVE_TENSION_CUDA_CHECK
        (
            cudaMemcpy
            (
                host,
                d_states + stateI*nCells,
                sizeof(double)*nCells,
                cudaMemcpyDeviceToHost
            )
        );
    }

    void downloadAlgebraicRow
    (
        double* host,
        const std::size_t algebraicI,
        const std::size_t nCells
    ) const
    {
        CARDIAC_ACTIVE_TENSION_CUDA_CHECK
        (
            cudaMemcpy
            (
                host,
                d_algebraics + algebraicI*nCells,
                sizeof(double)*nCells,
                cudaMemcpyDeviceToHost
            )
        );
    }

    void free()
    {
        if (d_states) cudaFree(d_states);
        if (d_rates) cudaFree(d_rates);
        if (d_algebraics) cudaFree(d_algebraics);
        if (d_constants) cudaFree(d_constants);
        if (d_driveSignals) cudaFree(d_driveSignals);
        d_states = d_rates = d_algebraics = d_constants = d_driveSignals = nullptr;
        allocated = false;
        constantsUploaded = false;
        stateResident = false;
    }
};

} // End namespace Foam

#endif
#endif
