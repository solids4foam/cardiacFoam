/*---------------------------------------------------------------------------*\
License
    This file is part of cardiacFoam.

    cardiacFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

Description
    Batched, inline implementation of the GoktepeKuhl model ODEs.

Author
    Simao Nieto de Castro, UCD.
\*---------------------------------------------------------------------------*/

#ifndef GoktepeKuhlBatch_H
#define GoktepeKuhlBatch_H

#include "../GoktepeKuhl/GoktepeKuhl_2004Names.H"

#ifndef CARDIAC_HOST_DEVICE
#define CARDIAC_HOST_DEVICE
#endif

#include <cmath>

CARDIAC_HOST_DEVICE inline void GoktepeKuhlComputeVariablesBatch
(
    const double driveSignal,
    const double* CONSTANTS,
    const int N,
    const int cellI,
    const double* STATES,
    double* RATES,
    double* ALGEBRAIC
)
{
    double uSignal = (driveSignal - CONSTANTS[AC_Vr]) / 100.0;
    if (uSignal < 0.0) uSignal = 0.0;
    if (uSignal > 1.0) uSignal = 1.0;

    ALGEBRAIC[AV_Vm * N + cellI] = driveSignal;
    ALGEBRAIC[AV_u * N + cellI]  = uSignal;

    ALGEBRAIC[AV_e * N + cellI] =
        CONSTANTS[AC_e0] + (CONSTANTS[AC_eInfty] - CONSTANTS[AC_e0])
            * exp(-exp(-CONSTANTS[AC_eXi] * (uSignal - CONSTANTS[AC_Vshift])));

    RATES[Ta * N + cellI] = ALGEBRAIC[AV_e * N + cellI] * (CONSTANTS[AC_kTa] * uSignal - STATES[Ta * N + cellI]);
}

#endif
