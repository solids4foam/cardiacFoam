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

#pragma once

enum BIDOMAIN_CONSTANTS_INDEX{
    BidomainCm,
    BidomainBeta,
    BidomainChi,
    BIDOMAIN_NUM_CONSTANTS
};
enum BIDOMAIN_STATES_INDEX{
    BidomainV,
    BidomainU1,
    BidomainU2,
    BidomainU3,
    BIDOMAIN_NUM_STATES
};
enum BIDOMAIN_ALGEBRAIC_INDEX{
    BidomainIion,
    BIDOMAIN_NUM_ALGEBRAIC
};

void
tmanufacturedFDABidomainInitConsts
(
    double* CONSTANTS,
    double* RATES,
    double* STATES,
    int dimension,
    double k
);

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
);
