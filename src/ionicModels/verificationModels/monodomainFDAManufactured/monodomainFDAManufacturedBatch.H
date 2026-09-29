/*---------------------------------------------------------------------------*\
License
    This file is part of cardiacFoam.
\*---------------------------------------------------------------------------*/

#pragma once

#include "monodomainFDAManufactured_2014Names.H"
#include "StimulusProtocolPOD.H"
#include <cmath>

namespace Foam
{

enum MONODOMAIN_FDA_BATCH_SUPPORT_INDEX
{
    MONODOMAIN_FDA_BATCH_SUPPORT_Iion,
    NUM_MONODOMAIN_FDA_BATCH_SUPPORT
};

static_assert
(
    NUM_MONODOMAIN_FDA_BATCH_SUPPORT == 1,
    "monodomainFDAManufactured batch support layout drift"
);

CARDIAC_HOST_DEVICE inline void monodomainFDAManufacturedComputeBatch
(
    const int N,
    const int beginCell,
    const int endCell,
    const double* __restrict__ CONSTANTS,
    const double* __restrict__ STATES_SoA,
    double* __restrict__ RATES_SoA,
    double* __restrict__ SUPPORT_SoA
)
{
    const int V_BASE = V*N;
    const int u1_BASE = u1*N;
    const int u2_BASE = u2*N;
    const int u3_BASE = u3*N;

    const double cm = CONSTANTS[Cm];
    const double beta = CONSTANTS[Beta];
    const double chi = CONSTANTS[Chi];

    for (int cellI = beginCell; cellI < endCell; ++cellI)
    {
        const double v = STATES_SoA[V_BASE + cellI];
        const double s1 = STATES_SoA[u1_BASE + cellI];
        const double s2 = STATES_SoA[u2_BASE + cellI];
        const double s3 = STATES_SoA[u3_BASE + cellI];
        const double q = s1 + s3 - v;
        const double s2sq = s2*s2;
        const double vMinusU3 = v - s3;

        RATES_SoA[V_BASE + cellI] = 0.0;
        RATES_SoA[u1_BASE + cellI] = q*q*s2sq + 0.5*q*s2sq*vMinusU3;
        RATES_SoA[u2_BASE + cellI] = -q*s2sq*s2;
        RATES_SoA[u3_BASE + cellI] = 0.0;

        SUPPORT_SoA[MONODOMAIN_FDA_BATCH_SUPPORT_Iion*N + cellI] =
            -0.5*cm*q*s2sq*vMinusU3 + (beta*vMinusU3)/chi;
    }
}

} // End namespace Foam
