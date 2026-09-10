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
    Batched, inline implementation of the NashPanfilov model ODEs.

Author
    Simao Nieto de Castro, UCD.
\*---------------------------------------------------------------------------*/

#ifndef NashPanfilovBatch_H
#define NashPanfilovBatch_H

#include "../NashPanfilov/NashPanfilov_2004Names.H"

#ifndef CARDIAC_HOST_DEVICE
#define CARDIAC_HOST_DEVICE
#endif

CARDIAC_HOST_DEVICE inline void NashPanfilovComputeVariablesBatch
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
    double uSignal =
        (driveSignal - CONSTANTS[AC_Vr])
      / (CONSTANTS[AC_Vp] - CONSTANTS[AC_Vr]);

    if (uSignal < 0.0) uSignal = 0.0;
    if (uSignal > 1.0) uSignal = 1.0;

    ALGEBRAIC[AV_u * N + cellI] = uSignal;

    ALGEBRAIC[AV_e * N + cellI] = ((uSignal < 0.05) ? 1.0 * CONSTANTS[AC_e0] : 10.0 * CONSTANTS[AC_e0]);
    RATES[Ta * N + cellI] = ALGEBRAIC[AV_e * N + cellI] * (CONSTANTS[AC_kTa] * uSignal - STATES[Ta * N + cellI]);
}

#endif
