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

#include "NashPanfilov_2004Names.H"

static const char* NashPanfilovSTATES_NAMES[NUM_STATES] = {
    "Ta"
};

static const char* NashPanfilovCONSTANTS_NAMES[NUM_CONSTANTS] = {
    "AC_Vp",
    "AC_Vr",
    "AC_Vth",
    "AC_e0",
    "AC_kTa",
    "AC_a"
};

static const char* NashPanfilovALGEBRAIC_NAMES[NUM_ALGEBRAIC] = {
    "AV_u",
    "AV_e",

};


static inline const Foam::HashTable<Foam::wordList>&
NashPanfilovDependencyMap()
{
   static const Foam::HashTable<Foam::wordList> dep
    (
        {
        }
    );

    return dep;
}


void
NashPanfilovinitConsts(double* CONSTANTS, double* RATES, double* STATES)
{
    CONSTANTS[AC_Vp] = 20.0;
    CONSTANTS[AC_Vr] = -80.0;
    CONSTANTS[AC_Vth] = -70.0;
    CONSTANTS[AC_e0] = 1.0;
    CONSTANTS[AC_kTa] = 47.9;

    STATES[Ta] = 0.0;
}


void
NashPanfilovcomputeVariables(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC)
{
    ALGEBRAIC[AV_e] = ((ALGEBRAIC[AV_u] < 0.05) ? 1.0 * CONSTANTS[AC_e0] : 10.0 * CONSTANTS[AC_e0]);
    RATES[Ta] = ALGEBRAIC[AV_e] * (CONSTANTS[AC_kTa] * ALGEBRAIC[AV_u] - STATES[Ta]);
}
