#pragma once


enum CONSTANTS_INDEX{
    R,              // 0  : gas constant
    T,              // 1  : temperature
    F,              // 2  : Faraday constant
    Cm,             // 3  : membrane capacitance (pF)
    V_c,            // 4  : cell volume (μm^3)

    // Stimulus (S1)
    stim_start,         // 5
    stim_period_S1,     // 6
    stim_duration,      // 7
    stim_amplitude,     // 8

    // Reversal potentials & ion concentrations
    P_kna,          // 9
    K_o,            // 10
    Na_o,           // 11
    Ca_o,           // 12

    // Conductances
    g_K1,           // 13
    g_Kr,           // 14
    g_Ks,           // 15
    g_Na,           // 16
    g_bNa,          // 17
    g_CaL,          // 18
    g_bCa,          // 19
    g_to,           // 20

    // Na/K Pump
    P_NaK,          // 21
    K_mk,           // 22
    K_mNa,          // 23

    // NCX exchanger
    K_NaCa,         // 24
    K_sat,          // 25
    alpha_NCX,      // 26
    gamma_NCX,      // 27
    Km_Ca,          // 28
    Km_Nai,         // 29

    // Calcium pump + potassium pump
    g_pCa,          // 30
    K_pCa,          // 31
    g_pK,           // 32

    // Calcium dynamics (SR release / uptake)
    tau_g,          // 33
    a_rel,          // 34
    b_rel,          // 35
    c_rel,          // 36
    K_up,           // 37
    V_leak,         // 38
    Vmax_up,        // 39

    // Buffers
    Buf_c,          // 40
    K_buf_c,        // 41
    Buf_sr,         // 42
    K_buf_sr,       // 43
    V_sr,           // 44
    tau_fCa,        // 45

    // Added multi-stimulus parameters (S1/S2 protocol)
    nstim1,         // 46
    stim_period_S2, // 47
    nstim2,         // 48

    NUM_CONSTANTS
};


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
TNNPcomputeVariables(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC,int tissueFlag, bool solveVmWithinODESolver);
void
TNNPcomputeRates(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC,int tissueFlag, bool solveVmWithinODESolver);

