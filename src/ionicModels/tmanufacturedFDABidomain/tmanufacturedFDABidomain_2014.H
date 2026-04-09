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

\*---------------------------------------------------------------------------*/

#include "tmanufacturedFDABidomain_2014Names.H"

const char* tmanufacturedFDABidomainSTATES_NAMES[BIDOMAIN_NUM_STATES] = {
    "V",
    "u1",
    "u2",
    "u3"
};
const char* tmanufacturedFDABidomainALGEBRAIC_NAMES[BIDOMAIN_NUM_ALGEBRAIC] = {
    "Iion"
};

const char* tmanufacturedFDABidomainCONSTANTS_NAMES[BIDOMAIN_NUM_CONSTANTS] = {
    "Cm",
    "Beta",
    "Chi",
};

void
tmanufacturedFDABidomainInitConsts
(
    double* CONSTANTS,
    double* RATES,
    double* STATES,
    int dimension,
    double k
)
{
    CONSTANTS[BidomainCm] = 2.0;
    CONSTANTS[BidomainChi] = 3.0;

    double betaBase = 0.0;

    if (dimension == 1)
    {
        betaBase = -1.1;
    }
    else if (dimension == 2)
    {
        betaBase = -5.9;
    }
    else if (dimension == 3)
    {
        betaBase = -8.6;
    }
    else
    {
        FatalErrorInFunction
            << "Invalid mesh dimension detected!" << exit(Foam::FatalError);
    }

    CONSTANTS[BidomainBeta] = betaBase*(1.0 - k);
}


void
tmanufacturedFDABidomainComputeVariables
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

    RATES[BidomainV] = 0.0;
    RATES[BidomainU1] =
        pow(STATES[BidomainU1] + STATES[BidomainU3] - STATES[BidomainV], 2.0)
      * pow(STATES[BidomainU2], 2.0)
      + 0.5
      * (STATES[BidomainU1] + STATES[BidomainU3] - STATES[BidomainV])
      * pow(STATES[BidomainU2], 2.0)
      * (STATES[BidomainV] - STATES[BidomainU3]);
    RATES[BidomainU2] =
        - (STATES[BidomainU1] + STATES[BidomainU3] - STATES[BidomainV])
        * pow(STATES[BidomainU2], 3);
    RATES[BidomainU3] = 0.0;

    ALGEBRAIC[BidomainIion] =
        -0.5*CONSTANTS[BidomainCm]
       * (STATES[BidomainU1] + STATES[BidomainU3] - STATES[BidomainV])
       * pow(STATES[BidomainU2], 2)
       * (STATES[BidomainV] - STATES[BidomainU3])
       + (CONSTANTS[BidomainBeta]
       * (STATES[BidomainV] - STATES[BidomainU3]))
       / CONSTANTS[BidomainChi];
}
