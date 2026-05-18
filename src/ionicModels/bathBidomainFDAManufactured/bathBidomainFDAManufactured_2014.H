/*---------------------------------------------------------------------------*\
License
    This file is part of cardiacFoam.

    cardiacFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    cardiacFoam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
    Public License for more details.

    You should have received a copy of the GNU General Public License
    along with cardiacFoam.  If not, see <http://www.gnu.org/licenses/>.

Description
    Constants and algebra for the bath-bidomain FDA manufactured solution.

    The bath variant differs from the heart-only FDA manufactured ionic
    model in exactly one place: the Beta scaling. The heart-only solution
    is built to make the (single-domain) bidomain balance vanish in a
    closed myocardium with zero-flux walls; the bath variant adds the
    surface-current "u3 = -alpha*x/se" linear drift, which shifts the
    dimensional balance and requires betaScale = -1.1 in every dimension.

    The dV/dt = 0 invariant and the (u1, u2) dynamics are unchanged. STATES,
    ALGEBRAIC and CONSTANTS layouts are owned by this header; see
    bathBidomainFDAManufactured_2014Names.H for the enum.
\*---------------------------------------------------------------------------*/

#include "bathBidomainFDAManufactured_2014Names.H"

const char* bathBidomainFDAManufacturedSTATES_NAMES[BATH_BIDOMAIN_NUM_STATES] =
{
    "V",
    "u1",
    "u2",
    "u3"
};

const char* bathBidomainFDAManufacturedALGEBRAIC_NAMES[BATH_BIDOMAIN_NUM_ALGEBRAIC] =
{
    "Iion"
};

const char* bathBidomainFDAManufacturedCONSTANTS_NAMES[BATH_BIDOMAIN_NUM_CONSTANTS] =
{
    "Cm",
    "Beta",
    "Chi"
};

void
bathBidomainFDAManufacturedInitConsts
(
    double* CONSTANTS,
    double* RATES,
    double* STATES,
    int dimension,
    double k
)
{
    (void)RATES;
    (void)STATES;
    (void)dimension;

    CONSTANTS[BathBidomainCm] = 2.0;
    CONSTANTS[BathBidomainChi] = 3.0;

    // The bath FDA manufactured solution carries u3 = -alpha*x/se in every
    // mesh dimension. Beta * (1 - k) absorbs the geometric mismatch from
    // this linear drift; the scaling is independent of dimensionality.
    const double betaScale = -1.1;

    CONSTANTS[BathBidomainBeta] = betaScale*(1.0 - k);
}


void
bathBidomainFDAManufacturedComputeVariables
(
    double VOI,
    double* CONSTANTS,
    double* RATES,
    double* STATES,
    double* ALGEBRAIC,
    int dimension,
    bool solveVmWithinODESolver
)
{
    (void)VOI;
    (void)dimension;
    (void)solveVmWithinODESolver;

    RATES[BathBidomainV] = 0.0;
    RATES[BathBidomainU1] =
        pow(STATES[BathBidomainU1] + STATES[BathBidomainU3] - STATES[BathBidomainV], 2.0)
      * pow(STATES[BathBidomainU2], 2.0)
      + 0.5
      * (STATES[BathBidomainU1] + STATES[BathBidomainU3] - STATES[BathBidomainV])
      * pow(STATES[BathBidomainU2], 2.0)
      * (STATES[BathBidomainV] - STATES[BathBidomainU3]);
    RATES[BathBidomainU2] =
        - (STATES[BathBidomainU1] + STATES[BathBidomainU3] - STATES[BathBidomainV])
        * pow(STATES[BathBidomainU2], 3);
    RATES[BathBidomainU3] = 0.0;

    ALGEBRAIC[BathBidomainIion] =
        -0.5*CONSTANTS[BathBidomainCm]
       * (STATES[BathBidomainU1] + STATES[BathBidomainU3] - STATES[BathBidomainV])
       * pow(STATES[BathBidomainU2], 2)
       * (STATES[BathBidomainV] - STATES[BathBidomainU3])
       + (CONSTANTS[BathBidomainBeta]
       * (STATES[BathBidomainV] - STATES[BathBidomainU3]))
       / CONSTANTS[BathBidomainChi];
}

// ************************************************************************* //
