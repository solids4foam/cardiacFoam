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

#include "GoktepeKuhl_2004Names.H"

static const char* GoktepeKuhlSTATES_NAMES[NUM_STATES] = {
    "Ta"
};

static const char* GoktepeKuhlCONSTANTS_NAMES[NUM_CONSTANTS] = {
    "AC_Vr",
    "AC_eInfty",
    "AC_e0",
    "AC_eXi",
    "AC_Vshift",
    "AC_kTa"
};

static const char* GoktepeKuhlALGEBRAIC_NAMES[NUM_ALGEBRAIC] = {
    "AV_e",
    "AV_Vm",
    "AV_u"
};


static inline const Foam::HashTable<Foam::wordList>&
GoktepeKuhlDependencyMap()
{
   static const Foam::HashTable<Foam::wordList> dep
    (
        {
        }
    );

    return dep;
}
void
GoktepeKuhlinitConsts(double* CONSTANTS, double* RATES, double* STATES)
{

    CONSTANTS[AC_Vr] = -80.0;

    // Smooth damping tension parameters
    CONSTANTS[AC_eInfty] = 10.0;     // mV-1
    CONSTANTS[AC_e0] = 1.0;         // mV-1
    CONSTANTS[AC_eXi] = 1.0;        // mV-1
    CONSTANTS[AC_Vshift] = 0.0;

    CONSTANTS[AC_kTa] = 47.9;      // Kpa mV−1

    /* initial conditions */
    STATES[Ta] = 0.0;
}


void
GoktepeKuhlcomputeVariables(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC)
{
    /* membrane_potential */

    //Tension Calculation
    ALGEBRAIC[AV_e] =
        CONSTANTS[AC_e0] + (CONSTANTS[AC_eInfty] - CONSTANTS[AC_e0])
            * exp(-exp(-CONSTANTS[AC_eXi] * (ALGEBRAIC[AV_u] - CONSTANTS[AC_Vshift])));

    RATES[Ta] = ALGEBRAIC[AV_e] * (CONSTANTS[AC_kTa] * ALGEBRAIC[AV_u] - STATES[Ta]);



}
