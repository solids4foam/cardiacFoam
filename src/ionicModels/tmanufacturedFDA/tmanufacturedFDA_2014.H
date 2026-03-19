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

\*---------------------------------------------------------------------------*/

#include "tmanufacturedFDA_2014Names.H"

const char* tmanufacturedFDASTATES_NAMES[NUM_STATES] = {
    "V", //
    "u1", //
    "u2", //
    "u3"  //

};
const char* tmanufacturedFDAALGEBRAIC_NAMES[NUM_ALGEBRAIC] = {
    "Iion"
};

const char* tmanufacturedFDACONSTANTS_NAMES[NUM_CONSTANTS] = {
    "Cm",
    "Beta",
    "Chi",
};

// Dim 1
// Dim 2
// Dim 3
void
tmanufacturedFDAinitConsts(double* CONSTANTS, double* RATES, double *STATES, int dimension)
{

    CONSTANTS[Cm]  = 2.0;
    CONSTANTS[Chi] = 3.0;

    if (dimension == 1)
    {
        CONSTANTS[Beta] = -1.1;
    }
    else if (dimension == 2)
    {
        CONSTANTS[Beta] = -5.9;
    }
    else if (dimension == 3)
    {
        CONSTANTS[Beta] = -8.6;
    }
    else
    {
        FatalErrorInFunction
            << "Invalid mesh dimension detected!" << exit(Foam::FatalError);
    }

}



void
tmanufacturedFDAcomputeVariables(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC,int dimension, bool solveVmWithinODESolver)
{
    // Local variables
    RATES[V] = 0.0;
    RATES[u1] = pow(STATES[u1] + STATES[u3] - STATES[V], 2.0) * pow(STATES[u2], 2.0) + 0.5 * (STATES[u1] + STATES[u3] - STATES[V]) * pow(STATES[u2], 2.0) * (STATES[V] - STATES[u3]);
    RATES[u2] = - (STATES[u1] + STATES[u3] - STATES[V]) *  pow(STATES[u2], 3);
    RATES[u3] = 0;
    ALGEBRAIC[Iion] =  - (STATES[u1] + STATES[u3] - STATES[V])  * pow(STATES[u2], 2) * (STATES[V] - STATES[u3]) + (CONSTANTS[Beta] * (STATES[V] - STATES[u3])) / CONSTANTS[Chi];

}
