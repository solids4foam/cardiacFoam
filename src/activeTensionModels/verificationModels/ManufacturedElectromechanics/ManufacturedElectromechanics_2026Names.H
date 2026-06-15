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
    Index enums and function prototypes for the manufactured
    electromechanics active tension model.

Author
    Simao Nieto de Castro, UCD.
\*---------------------------------------------------------------------------*/

#pragma once

enum MANUFACTURED_ELECTROMECHANICS_STATES_INDEX
{
    MMS_STATE_Ta,
    NUM_STATES
};

enum MANUFACTURED_ELECTROMECHANICS_ALGEBRAIC_INDEX
{
    MMS_AV_Vm,
    MMS_AV_lambda,
    MMS_AV_voltageActivation,
    MMS_AV_lengthFactor,
    NUM_ALGEBRAIC
};

enum MANUFACTURED_ELECTROMECHANICS_CONSTANTS_INDEX
{
    MMS_AC_Tmax,
    MMS_AC_V0,
    MMS_AC_gamma,
    NUM_CONSTANTS
};

inline void ManufacturedElectromechanics2026initConsts
(
    double* CONSTANTS,
    double* RATES,
    double* STATES
);

#ifndef CARDIAC_HOST_DEVICE
#define CARDIAC_HOST_DEVICE
#endif

CARDIAC_HOST_DEVICE inline void ManufacturedElectromechanics2026computeVariables
(
    double VOI,
    double* CONSTANTS,
    double* RATES,
    double* STATES,
    double* ALGEBRAIC
);

// ************************************************************************* //
