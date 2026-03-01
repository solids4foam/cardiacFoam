#pragma once
#include "dictionary.H"



enum CONSTANT_INDEX{    AC_CMDN_max,    AC_CSQN_max,    AC_Ca_up_max,    AC_Cao,    AC_Cm,
    AC_ECaL,    AC_F,    AC_FRT,    AC_INaCa_max,    AC_INaK_max,    AC_I_diff,    AC_I_up_max,
    AC_IpCa_max,    AC_KQ10,    AC_K_rel,    AC_K_up,    AC_KmCa,    AC_KmKo,    AC_KmNa,
    AC_KmNai,    AC_Km_CMDN,    AC_Km_CSQN,    AC_Km_TRPN,    AC_Ko,    AC_Nao,    AC_R,
    AC_RTF,    AC_T,    AC_TRPN_max,    AC_V_cell,    AC_V_i,    AC_V_rel,    AC_V_up,
    AC_c1,    AC_c2,   AC_cajsr_u_tau,    AC_g,    AC_gCaL,    AC_gK1,    AC_gKr,    AC_gKs,
    AC_gKur_base,   AC_gNa,    AC_gbCa,    AC_gbNa,    AC_gto,    AC_ical_fCa_tau,    AC_ksat,    AC_sigma,    AC_tau_tr,

    stim_start, stim_duration, stim_amplitude,
    stim_period_S1, nstim1,
    stim_period_S2, nstim2,

    NUM_CONSTANTS};

enum STATE_INDEX {
    membrane_V,
    sodium_Nai,
    potassium_Ki,
    calcium_Cai,
    calcium_CaUp,
    calcium_CaRel,
    ina_m,
    ina_h,
    ina_j,
    ito_oa,
    ito_oi,
    ikur_ua,
    ikur_ui,
    ikr_xr,
    iks_xs,
    ical_d,
    ical_f,
    ical_fCa,
    cajsr_u,
    cajsr_v,
    cajsr_w,

    NUM_STATES
};

enum ALGEBRAIC_INDEX
{
    // --- Currents ---
    AV_ICaL,
    AV_IK1,
    AV_IKr,
    AV_IKs,
    AV_IKur,
    AV_INa,
    AV_INaCa,
    AV_INaK,
    AV_IbCa,
    AV_IbNa,
    AV_IpCa,
    AV_Ito,

    // --- INa gating ---
    AV_ina_m_inf,
    AV_ina_m_tau,
    AV_ina_h_inf,
    AV_ina_h_tau,
    AV_ina_j_inf,
    AV_ina_j_tau,

    // --- ICaL gating ---
    AV_ical_d_inf,
    AV_ical_d_tau,
    AV_ical_f_inf,
    AV_ical_f_tau,

    // --- Ito gating ---
    AV_ito_oa_inf,
    AV_ito_oa_tau,
    AV_ito_oi_inf,
    AV_ito_oi_tau,

    // --- IKur gating ---
    AV_ikur_ua_inf,
    AV_ikur_ua_tau,
    AV_ikur_ui_inf,
    AV_ikur_ui_tau,
    AV_gKur,

    // --- IKr gating ---
    AV_ikr_xr_inf,
    AV_ikr_xr_tau,

    // --- IKs gating ---
    AV_iks_xs_inf,
    AV_iks_xs_tau,

    // --- Other V-dependent ---
    AV_fNaK,

    // --- CaJSR gating ---
    AV_cajsr_w_inf,
    AV_cajsr_w_tau,

    // --- Totals / stimulus ---
    Iion_cm,
    Istim,

    NUM_ALGEBRAIC
};

void
CourtemancheinitConsts(double* CONSTANTS, double* RATES, double *STATES, int tissueFlag, const Foam::dictionary& stimulus);
void
CourtemanchecomputeVariables(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC,int tissueFlag, bool solveVmWithinODESolver);


