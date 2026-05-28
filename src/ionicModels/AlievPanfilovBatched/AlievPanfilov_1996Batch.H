/*---------------------------------------------------------------------------*\
License
    This file is part of cardiacFoam.

Description
    State-major SoA batched evaluator for the Aliev-Panfilov 1996 ionic model.
    This model uses the reduced voltage state `u`; Vm/Im scaling is handled
    by the wrapper metadata and extraction path.
\*---------------------------------------------------------------------------*/

#pragma once

#include "AlievPanfilov_1996Names.H"
#include "StimulusProtocolPOD.H"
#include <cmath>

namespace Foam
{

enum ALIEVPANFILOV_BATCH_SUPPORT_INDEX
{
    ALIEVPANFILOV_BATCH_SUPPORT_tau_recovery_r,
    ALIEVPANFILOV_BATCH_SUPPORT_gInf_recovery_r,
    ALIEVPANFILOV_BATCH_SUPPORT_Iion_cm,

    NUM_ALIEVPANFILOV_BATCH_SUPPORT
};

static_assert(NUM_ALIEVPANFILOV_BATCH_SUPPORT == 3, "AlievPanfilov compact support layout drift");

CARDIAC_HOST_DEVICE inline void AlievPanfilovComputeVariablesBatch
(
    const double VOI,
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
#ifdef _OPENMP
    #pragma omp simd
#endif
    for (int cellI = beginCell; cellI < endCell; ++cellI)
    {
        const double state_u = STATES_SoA[u * N + cellI];
        const double state_recovery_r =
            STATES_SoA[recovery_r * N + cellI];

        const double eps =
            CONSTANTS[AC_epsilon]
          + CONSTANTS[AC_mu1]*state_recovery_r
           /(CONSTANTS[AC_mu2] + state_u);

        const double recoverySteady =
            -CONSTANTS[AC_k]*state_u*(state_u - (CONSTANTS[AC_a] + 1.0));

        const double rate_recovery_r =
            eps*(recoverySteady - state_recovery_r);

        const double Iion_cm =
            CONSTANTS[AC_k]*state_u*(state_u - CONSTANTS[AC_a])
           *(state_u - 1.0)
          + state_recovery_r*state_u;

        const double Istim = computeStimulusPOD(VOI, stimulus);
        const double rate_u = solveVm ? (-Iion_cm - Istim) : 0.0;

        RATES_SoA[u * N + cellI] = rate_u;
        RATES_SoA[recovery_r * N + cellI] = rate_recovery_r;

        if (SUPPORT_SoA != nullptr)
        {
            SUPPORT_SoA
            [
                ALIEVPANFILOV_BATCH_SUPPORT_tau_recovery_r * N + cellI
            ] = 1.0/eps;

            SUPPORT_SoA
            [
                ALIEVPANFILOV_BATCH_SUPPORT_gInf_recovery_r * N + cellI
            ] = recoverySteady;

            SUPPORT_SoA
            [
                ALIEVPANFILOV_BATCH_SUPPORT_Iion_cm * N + cellI
            ] = Iion_cm;
        }
    }
}

} // End namespace Foam

// ************************************************************************* //
