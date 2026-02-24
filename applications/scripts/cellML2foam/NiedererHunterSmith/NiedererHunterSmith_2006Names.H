#pragma once

enum STATES_INDEX {
    Ca_i_TRPN,
    tropomyosin_z,
    Cross_Bridges_Q1,
    Cross_Bridges_Q2,
    Cross_Bridges_Q3,
    NUM_STATES
};

enum ALGEBRAIC_INDEX {
    AV_ExtensionRatio,
    AV_lambda,
    AV_lambda_prev,
    AV_overlap,
    AV_Q,
    AV_Ca_i,
    AV_Ca_50,
    AV_beta_Tm,
    AV_Ca_b,
    AV_Ca_TRPN_50,
    AV_alpha_Tm,
    AV_z_max,
    AV_T_Base,
    AV_T_0,
    AV_Tension,
    AV_k_off,
    AV_J_TRPN,
    NUM_ALGEBRAIC
};

enum CONSTANTS_INDEX {
    AC_dExtensionRatiodt,
    AC_beta_0,
    AC_A_1,
    AC_A_2,
    AC_A_3,
    AC_a,
    AC_alpha_1,
    AC_alpha_2,
    AC_alpha_3,
    AC_T_ref,
    AC_Ca_50ref,
    AC_K_z,
    AC_alpha_0,
    AC_alpha_r1,
    AC_alpha_r2,
    AC_beta_1,
    AC_n_Hill,
    AC_n_Rel,
    AC_z_p,
    AC_K_1,
    AC_K_2,
    AC_Ca_TRPN_Max,
    AC_gamma_trpn,
    AC_k_Ref_off,
    AC_k_on,
    NUM_CONSTANTS
};


void
NiedererHunterSmithinitConsts(double* CONSTANTS,double* RATES,double* STATES);

void
NiedererHunterSmithcomputeVariables(double VOI,double* CONSTANTS,double* RATES,double* STATES,double* ALGEBRAIC,bool solveVmWithinODESolver);
