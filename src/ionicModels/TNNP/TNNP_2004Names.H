#pragma once
#include "stimulusIO.H"


enum CONSTANTS_INDEX {
    R,
    T,
    F,
    Cm,
    V_c,
    P_kna,
    K_o,
    Na_o,
    Ca_o,
    g_K1,
    g_Kr,
    g_Ks,
    g_Na,
    g_bNa,
    g_CaL,
    g_bCa,
    g_to,
    P_NaK,
    K_mk,
    K_mNa,
    K_NaCa,
    K_sat,
    alpha_NCX,
    gamma_NCX,
    Km_Ca,
    Km_Nai,
    g_pCa,
    K_pCa,
    g_pK,
    tau_g,
    a_rel,
    b_rel,
    c_rel,
    K_up,
    V_leak,
    Vmax_up,
    Buf_c,
    K_buf_c,
    Buf_sr,
    K_buf_sr,
    V_sr,
    tau_fCa,

    NUM_CONSTANTS};


enum STATES_INDEX {
    V,      // membrane voltage
    K_i,    // intracellular potassium
    Na_i,   // intracellular sodium
    Ca_i,   // intracellular calcium
    Xr1,
    Xr2,
    Xs,
    m,
    h,
    j,
    d,
    f,
    fCa,
    s,
    r,
    Ca_SR,
    g,
    NUM_STATES
};


enum ALGEBRAIC_INDEX {
    Istim,          
    xr1_inf,
    xr2_inf,
    xs_inf,
    m_inf,
    h_inf,
    j_inf,
    d_inf,
    f_inf,
    alpha_fCa,
    s_inf,
    r_inf,
    g_inf,
    E_Na,
    alpha_xr1,
    alpha_xr2,
    alpha_xs,
    alpha_m,
    alpha_h,
    alpha_j,
    alpha_d,
    tau_f,
    beta_fCa,
    tau_s,
    tau_r,
    d_g,
    E_K,
    beta_xr1,
    beta_xr2,
    beta_xs,
    beta_m,
    beta_h,
    beta_j,
    gama_fCa,
    beta_d,
    E_Ks,
    tau_xr1,
    tau_xr2,
    tau_xs,
    tau_m,
    tau_h,
    tau_j,
    gamma_d,
    fCa_inf,
    E_Ca,
    tau_d,
    d_fCa,
    alpha_K1,
    beta_K1,
    xK1_inf,
    i_K1,
    i_Kr,
    i_Ks,
    i_Na,
    i_b_Na,
    i_CaL,
    i_b_Ca,
    i_to,
    i_NaK,
    i_NaCa,
    i_p_Ca,
    i_p_K,
    i_rel,
    i_up,
    i_leak,
    ddt_Ca_i_total,
    ddt_Ca_sr_total,
    f_JCa_i_free,
    f_JCa_sr_free,
    Iion_cm,

    NUM_ALGEBRAIC
};



void
TNNPinitConsts(double* CONSTANTS, double* RATES, double *STATES, int tissueFlag, const Foam::dictionary& stimulus);
void
TNNPcomputeVariables(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC,int tissueFlag, bool solveVmWithinODESolver, const Foam::StimulusProtocol& stimulus);
void
TNNPcomputeRates(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC,int tissueFlag, bool solveVmWithinODESolver, const Foam::StimulusProtocol& stimulus);

