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

#pragma once
#include "stimulusIO.H"

enum STATES_INDEX {
    u,
    recovery_r,
    NUM_STATES
};

enum ALGEBRAIC_INDEX {
    AV_eps,
    Istim,
    Iion_cm,
    NUM_ALGEBRAIC
};

enum CONSTANTS_INDEX {
    AC_Vp,
    AC_Vr,
    AC_Vth,
    AC_epsilon,
    AC_k,
    AC_mu1,
    AC_mu2,
    AC_a,

    NUM_CONSTANTS};


void
AlievPanfilovinitConsts(double* CONSTANTS,double* RATES,double* STATES,int tissueFlag,const Foam::dictionary& stimulus);

void
AlievPanfilovcomputeVariables(double VOI,double* CONSTANTS,double* RATES,double* STATES,double* ALGEBRAIC,int tissueFlag,bool solveVmWithinODESolver, const Foam::StimulusProtocol& stimulus);









