#pragma once

enum STATES_INDEX {

    Ta,
    NUM_STATES
};

enum ALGEBRAIC_INDEX {
    AV_u,
    AV_e,

    NUM_ALGEBRAIC
};

enum CONSTANTS_INDEX {

    AC_Vp,
    AC_Vr,
    AC_Vth,
    AC_e0,
    AC_kTa,
    AC_a,
    NUM_CONSTANTS
};


void
NashPanfilovinitConsts(double* CONSTANTS,double* RATES,double* STATES);

void
NashPanfilovcomputeVariables(double VOI,double* CONSTANTS,double* RATES,double* STATES,double* ALGEBRAIC);
