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
    Index enums and function prototypes for the Goktepe-Kuhl 2004
    active tension model.

Author
    Simao Nieto de Castro, UCD.
\*---------------------------------------------------------------------------*/

#pragma once

enum STATES_INDEX {
    Ta,
    NUM_STATES
};

enum ALGEBRAIC_INDEX {
    AV_e,
    AV_Vm,
    AV_u,
    NUM_ALGEBRAIC
};

enum CONSTANTS_INDEX {
    AC_Vr,
    AC_eInfty,
    AC_e0,
    AC_eXi,
    AC_Vshift,
    AC_kTa,
    NUM_CONSTANTS
};


inline void
GoktepeKuhlinitConsts(double* CONSTANTS, double* RATES, double* STATES);

#ifndef CARDIAC_HOST_DEVICE
#define CARDIAC_HOST_DEVICE
#endif

CARDIAC_HOST_DEVICE inline void
GoktepeKuhlcomputeVariables
(
    double VOI,
    double* CONSTANTS,
    double* RATES,
    double* STATES,
    double* ALGEBRAIC
);
