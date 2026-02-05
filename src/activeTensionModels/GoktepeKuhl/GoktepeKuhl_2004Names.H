#pragma once

enum STATES_INDEX {
    Ta,
    AV_u,
    NUM_STATES
};

enum ALGEBRAIC_INDEX { 
    AV_e,
    AV_Vm,
    NUM_ALGEBRAIC
};

enum CONSTANTS_INDEX {
    AC_Vr,
    AC_eInfty,
    AC_e0,
    AC_eXi,
    AC_Vshift,
    AC_kTa,
    NUM_CONSTANTS
};


void
GoktepeKuhlinitConsts(double* CONSTANTS,double* RATES,double* STATES);

void
GoktepeKuhlcomputeVariables(double VOI,double* CONSTANTS,double* RATES,double* STATES,double* ALGEBRAIC);
