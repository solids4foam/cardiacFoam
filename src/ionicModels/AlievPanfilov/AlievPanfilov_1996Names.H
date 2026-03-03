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









