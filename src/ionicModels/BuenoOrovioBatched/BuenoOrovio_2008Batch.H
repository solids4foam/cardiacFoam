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

Description
    State-major (SoA) batched evaluator for the reduced Bueno-Orovio 2008
    ionic model.

    Why this exists
    ---------------
    The existing `BuenoOroviocomputeVariables` in BuenoOrovio_2008.H is a
    *cell-major* evaluator: it takes one cell's state vector at a time
    via per-cell scratch buffers (`scalarField`-shaped) and writes the
    result into a per-cell rate / algebraic vector. That layout cannot be
    vectorised across cells by the compiler and cannot be launched as a
    one-thread-per-cell device kernel.

    `BuenoOrovioComputeVariablesBatch` is the same equation set rewritten
    so the *outer* index is the cell and every state/rate access is a
    simple SoA stride: `STATES[stateI*N + cell]`. This shape:

      * Lets the host build vectorise across cells with `#pragma omp simd`
        (AVX2/AVX-512 lanes process the math in parallel).
      * Is the exact shape a CUDA / HIP / OpenMP-target kernel needs:
        wrap the body in `for (cell = begin + tid; cell < end; cell += stride)`
        or launch one thread per cell, and the same source compiles for
        the device.
      * Takes a POD stimulus (`StimulusProtocolPOD`) - no Foam types
        cross the device boundary.

    For now this is host-only. The `CARDIAC_HOST_DEVICE` qualifier is
    set up so adding `__device__` later is a build-flag change, not a
    code rewrite. See `gpu_ionic_model_audit.md` Section 6.1.

    Numerical equivalence
    ---------------------
    The body uses the same hard-step `heaviside` that
    `BuenoOroviocomputeVariables` uses, so a per-cell run of this batched
    function should match the legacy evaluator to floating-point exactness
    (modulo associativity reordering the compiler may pick).
\*---------------------------------------------------------------------------*/

#ifndef BuenoOrovio_2008Batch_H
#define BuenoOrovio_2008Batch_H

#include "BuenoOrovio_2008Names.H"
#include "StimulusProtocolPOD.H"
#include <cmath>

// Local hard-step Heaviside that is host/device callable. Kept identical
// in semantics to `Foam::heaviside` (see BuenoOrovio/heaviside.H) so the
// batched and legacy paths match bit-for-bit on the same input.
CARDIAC_HOST_DEVICE inline double bo_heaviside(const double x)
{
    return (x > 0.0) ? 1.0 : 0.0;
}


namespace Foam
{

//- Compact-support index map for the BuenoOrovio batched evaluator.
//  Compact-support layout: keeps only the slots needed
//  on the hot path. `BO_BATCH_SUPPORT_Iion` is the only required slot
//  for a Vm-only solve; the gate-tau values are exposed for an
//  eventual Rush-Larsen-from-support path on this model.
enum BO_BATCH_SUPPORT_INDEX
{
    BO_BATCH_SUPPORT_tauVMinus,
    BO_BATCH_SUPPORT_tauWMinus,
    BO_BATCH_SUPPORT_tauS,
    BO_BATCH_SUPPORT_Iion,

    NUM_BO_BATCH_SUPPORT
};


//- State-major batched evaluator. One call processes cells
//  `[beginCell, endCell)`. State, rate and (optional) support arrays
//  are SoA: index = `stateI * N + cellI`. Pass `SUPPORT = nullptr` if
//  no compact-support output is required (e.g. Vm-only Euler step).
//
//  Parameters
//    t            : current model time (ms).
//    CONSTANTS    : the model's constants array, identical to the one
//                   built by `BuenoOrovioinitConsts(...)`. Read-only on
//                   device; pass via constant memory or by value through
//                   a kernel-launch wrapper if desired.
//    N            : total number of cells in the SoA buffers (i.e. the
//                   stride between consecutive states for the same cell).
//    beginCell,
//    endCell      : half-open cell range to process.
//    STATES_SoA   : input states. Layout STATES_SoA[stateI*N + cellI].
//    RATES_SoA    : output rates, same SoA layout. Always written.
//    SUPPORT_SoA  : optional output, layout SUPPORT_SoA[supportI*N + cellI].
//                   Pass `nullptr` if not needed.
//    solveVm      : if true, write `RATES[u]` from `-Jion - Istim`. If
//                   false (Vm coupled externally), `RATES[u]` is left
//                   untouched.
//    stimulus     : POD stimulus protocol; passed by value.
//
//  Tissue type is not a per-evaluator argument: it is baked into CONSTANTS
//  by `BuenoOrovioinitConsts` at construction time. This matches the
//  GaurComputeVariablesBatch / TNNPComputeRatesBatch / ToRORd_dynClComputeVariablesBatch
//  signatures.
CARDIAC_HOST_DEVICE inline void BuenoOrovioComputeVariablesBatch
(
    const double t,
    const double* __restrict__ CONSTANTS,
    const int N,
    const int beginCell,
    const int endCell,
    const double* __restrict__ STATES_SoA,
    double*       __restrict__ RATES_SoA,
    double*       __restrict__ SUPPORT_SoA,
    const bool solveVm,
    const StimulusProtocolPOD& stimulus
)
{
    // Stimulus is loop-invariant in cell - compute once.
    const double Istim = computeStimulusPOD(t, stimulus);

    // Hoist constant reads outside the cell loop. The compiler can hoist
    // these on its own under `-O2`/`-O3`, but doing it explicitly makes
    // the SIMD step robust under O1 / debug builds and signals intent.
    const double c_uO         = CONSTANTS[uO];
    const double c_uU         = CONSTANTS[uU];
    const double c_thetaV     = CONSTANTS[thetaV];
    const double c_thetaW     = CONSTANTS[thetaW];
    const double c_thetaVm    = CONSTANTS[thetaVMinus];
    const double c_thetaO     = CONSTANTS[thetaO];
    const double c_tauV1m     = CONSTANTS[tauV1Minus];
    const double c_tauV2m     = CONSTANTS[tauV2Minus];
    const double c_tauVPlus   = CONSTANTS[tauVPlus];
    const double c_tauW1m     = CONSTANTS[tauW1Minus];
    const double c_tauW2m     = CONSTANTS[tauW2Minus];
    const double c_kWm        = CONSTANTS[kWMinus];
    const double c_uWm        = CONSTANTS[uWMinus];
    const double c_tauWPlus   = CONSTANTS[tauWPlus];
    const double c_tauFi      = CONSTANTS[tauFi];
    const double c_tauO1      = CONSTANTS[tauO1];
    const double c_tauO2      = CONSTANTS[tauO2];
    const double c_tauSo1     = CONSTANTS[tauSo1];
    const double c_tauSo2     = CONSTANTS[tauSo2];
    const double c_kSo        = CONSTANTS[kSo];
    const double c_uSo        = CONSTANTS[uSo];
    const double c_tauS1      = CONSTANTS[tauS1];
    const double c_tauS2      = CONSTANTS[tauS2];
    const double c_kS         = CONSTANTS[kS];
    const double c_uS         = CONSTANTS[uS];
    const double c_tauSi      = CONSTANTS[tauSi];
    const double c_tauWInfty  = CONSTANTS[tauWInfty];
    const double c_wInftyStar = CONSTANTS[wInftyStar];

    // SoA strides for each state. Resolved at compile time on the device
    // when `u`, `v`, etc. are integer enum values.
    const int U_BASE = u * N;
    const int V_BASE = v * N;
    const int W_BASE = w * N;
    const int S_BASE = s * N;

#ifdef _OPENMP
    #pragma omp simd
#endif
    for (int cellI = beginCell; cellI < endCell; ++cellI)
    {
        // SoA loads - one scalar per state per cell.
        const double Su = STATES_SoA[U_BASE + cellI];
        const double Sv = STATES_SoA[V_BASE + cellI];
        const double Sw = STATES_SoA[W_BASE + cellI];
        const double Ss = STATES_SoA[S_BASE + cellI];

        // Heaviside switches - these are the divergence sources on a
        // SIMT device. On CPU SIMD they predicate cleanly.
        const double hV  = bo_heaviside(Su - c_thetaV);
        const double hW  = bo_heaviside(Su - c_thetaW);
        const double hVm = bo_heaviside(Su - c_thetaVm);
        const double hO  = bo_heaviside(Su - c_thetaO);

        // Time constants
        const double tauSoLocal =
            c_tauSo1 + (c_tauSo2 - c_tauSo1)
            * (1.0 + std::tanh(c_kSo * (Su - c_uSo))) * 0.5;

        const double tauOLocal = (1.0 - hO) * c_tauO1 + hO * c_tauO2;

        const double tauVMinusLocal = (1.0 - hVm) * c_tauV1m + hVm * c_tauV2m;

        const double tauWMinusLocal =
            c_tauW1m + (c_tauW2m - c_tauW1m)
            * (1.0 + std::tanh(c_kWm * (Su - c_uWm))) * 0.5;

        const double tauSLocal = (1.0 - hW) * c_tauS1 + hW * c_tauS2;

        // Currents
        const double Jfi =
            -Sv * hV * (Su - c_thetaV) * (c_uU - Su) / c_tauFi;

        const double Jso =
            (Su - c_uO) * (1.0 - hW) / tauOLocal
          + hW / tauSoLocal;

        const double Jsi = -hW * Sw * Ss / c_tauSi;

        const double Jion = Jfi + Jso + Jsi;

        // Steady-state targets
        const double vInftyLocal = bo_heaviside(c_thetaVm - Su);

        const double wInftyLocal =
            (1.0 - hO) * (1.0 - Su / c_tauWInfty)
          + hO * c_wInftyStar;

        // Rates (writes go straight to SoA; rate of u handled below).
        RATES_SoA[V_BASE + cellI] =
            (1.0 - hV) * (vInftyLocal - Sv) / tauVMinusLocal
          - hV * Sv / c_tauVPlus;

        RATES_SoA[W_BASE + cellI] =
            (1.0 - hW) * (wInftyLocal - Sw) / tauWMinusLocal
          - hW * Sw / c_tauWPlus;

        RATES_SoA[S_BASE + cellI] =
            ((1.0 + std::tanh(c_kS * (Su - c_uS))) * 0.5 - Ss) / tauSLocal;

        if (solveVm)
        {
            RATES_SoA[U_BASE + cellI] = -Jion - Istim;
        }

        // Optional compact-support output. Layout: SUPPORT[supportI*N + cellI].
        if (SUPPORT_SoA != nullptr)
        {
            SUPPORT_SoA[BO_BATCH_SUPPORT_tauVMinus * N + cellI] = tauVMinusLocal;
            SUPPORT_SoA[BO_BATCH_SUPPORT_tauWMinus * N + cellI] = tauWMinusLocal;
            SUPPORT_SoA[BO_BATCH_SUPPORT_tauS      * N + cellI] = tauSLocal;
            // BuenoOrovio's `Iion` is unitless (Jion); the host wrapper
            // applies the `*85.7` scale (see BuenoOrovioBatched.H::
            // ionicCurrentFromEvaluation). We store the raw Jion here
            // and let the wrapper scale at extraction time, matching the
            // existing host-path contract.
            SUPPORT_SoA[BO_BATCH_SUPPORT_Iion      * N + cellI] = Jion;
        }
    }
}

} // End namespace Foam

#endif

// ************************************************************************* //
