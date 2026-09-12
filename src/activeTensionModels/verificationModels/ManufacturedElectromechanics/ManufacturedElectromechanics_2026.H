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
    Self-contained algebraic equations for the manufactured electromechanics
    active tension model.

Author
    Simao Nieto de Castro, UCD.
\*---------------------------------------------------------------------------*/

#include "ManufacturedElectromechanics_2026Names.H"

static const char* ManufacturedElectromechanicsSTATES_NAMES[NUM_STATES] =
{
    "Ta"
};

static const char* ManufacturedElectromechanicsCONSTANTS_NAMES[NUM_CONSTANTS] =
{
    "Tmax",
    "V0",
    "gamma"
};

static const char* ManufacturedElectromechanicsALGEBRAIC_NAMES[NUM_ALGEBRAIC] =
{
    "Vm",
    "lambda",
    "voltageActivation",
    "lengthFactor"
};

inline void ManufacturedElectromechanics2026initConsts
(
    double* CONSTANTS,
    double* RATES,
    double* STATES
)
{
    CONSTANTS[MMS_AC_Tmax] = 1.0;
    CONSTANTS[MMS_AC_V0] = 1.0;
    CONSTANTS[MMS_AC_gamma] = 1.0;

    STATES[MMS_STATE_Ta] = 0.0;
    RATES[MMS_STATE_Ta] = 0.0;
}


#ifndef CARDIAC_HOST_DEVICE
#define CARDIAC_HOST_DEVICE
#endif

CARDIAC_HOST_DEVICE inline void ManufacturedElectromechanics2026computeVariables
(
    double,
    double* CONSTANTS,
    double* RATES,
    double* STATES,
    double* ALGEBRAIC
)
{
    const double Vm = ALGEBRAIC[MMS_AV_Vm];
    const double lambda = ALGEBRAIC[MMS_AV_lambda];
    const double V02 = CONSTANTS[MMS_AC_V0]*CONSTANTS[MMS_AC_V0];
    const double Vm2 = Vm*Vm;

    ALGEBRAIC[MMS_AV_voltageActivation] = Vm2/(V02 + Vm2);
    ALGEBRAIC[MMS_AV_lengthFactor] =
        1.0 + CONSTANTS[MMS_AC_gamma]*(lambda - 1.0);

    STATES[MMS_STATE_Ta] =
        CONSTANTS[MMS_AC_Tmax]
       *ALGEBRAIC[MMS_AV_voltageActivation]
       *ALGEBRAIC[MMS_AV_lengthFactor];

    RATES[MMS_STATE_Ta] = 0.0;
}

// ************************************************************************* //
