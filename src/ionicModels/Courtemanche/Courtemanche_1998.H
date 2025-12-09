#include "Courtemanche_1998Names.H"
#include "dictionary.H"


const char* CourtemancheSTATES_NAMES[NUM_STATES] = {
    "membrane_V",       // Membrane voltage
    "sodium_Nai",       // Intracellular sodium concentration
    "potassium_Ki",     // Intracellular potassium concentration
    "calcium_Cai",      // Intracellular calcium concentration
    "calcium_CaUp",     // Calcium in SR uptake compartment
    "calcium_CaRel",    // Calcium in SR release compartment
    "ina_m",            // INa activation gate
    "ina_h",            // INa fast inactivation gate
    "ina_j",            // INa slow inactivation gate
    "ito_oa",           // Ito activation gate
    "ito_oi",           // Ito inactivation gate
    "ikur_ua",          // IKur activation gate
    "ikur_ui",          // IKur inactivation gate
    "ikr_xr",           // IKr activation gate
    "iks_xs",           // IKs activation gate
    "ical_d",           // ICaL activation gate
    "ical_f",           // ICaL voltage-dependent inactivation
    "ical_fCa",         // ICaL calcium-dependent inactivation
    "cajsr_u",          // RyR u-state (SR calcium release)
    "cajsr_v",          // RyR v-state
    "cajsr_w"           // RyR w-state
};


const char* CourtemancheALGEBRAIC_NAMES[NUM_ALGEBRAIC] = {
    "AV_ICaL",     // L-type calcium current
    "AV_IK1",      // Inward rectifier potassium current
    "AV_IKr",      // Rapid delayed rectifier potassium current
    "AV_IKs",      // Slow delayed rectifier potassium current
    "AV_IKur",     // Ultrarapid delayed rectifier potassium current
    "AV_INa",      // Fast sodium current
    "AV_INaCa",    // Sodium-calcium exchanger current
    "AV_INaK",     // Sodium-potassium pump current
    "AV_IbCa",     // Background calcium current
    "AV_IbNa",     // Background sodium current
    "AV_IpCa",     // Calcium pump current
    "AV_Ito",      // Transient outward potassium current
    "Iion_cm",    // Total ionic current
    "AV_I_stim"    // External stimulus current
};

void
CourtemancheinitConsts(double* CONSTANTS, double* RATES, double *STATES, int tissueFlag,  const Foam::dictionary& stimulus)
{

    /* ca_buffers */
    CONSTANTS[AC_CMDN_max] = 0.05;
    CONSTANTS[AC_CSQN_max] = 10.0;
    CONSTANTS[AC_Km_CMDN] = 0.00238;
    CONSTANTS[AC_Km_CSQN] = 0.8;
    CONSTANTS[AC_Km_TRPN] = 0.0005;
    CONSTANTS[AC_TRPN_max] = 0.07;
    
    /* cansr */
    CONSTANTS[AC_Ca_up_max] = 15.0;
    CONSTANTS[AC_I_up_max] = 0.005;
    CONSTANTS[AC_K_up] = 0.00092;
    
    /* extra */
    CONSTANTS[AC_Cao] = 1.8;
    CONSTANTS[AC_Ko] = 5.4;
    CONSTANTS[AC_Nao] = 140.0;
    
    /* geom */
    CONSTANTS[AC_Cm] = 100.0;
    CONSTANTS[AC_V_cell] = 20100.0;
    CONSTANTS[AC_V_i] = CONSTANTS[AC_V_cell] * 0.68;
    CONSTANTS[AC_V_rel] = 0.0048 * CONSTANTS[AC_V_cell];
    CONSTANTS[AC_V_up] = 0.0552 * CONSTANTS[AC_V_cell];
    
    /* ical */
    CONSTANTS[AC_ECaL] = 65.0;
    CONSTANTS[AC_gCaL] = 0.12375;
    CONSTANTS[AC_ical_fCa_tau] = 2.0;
    
    /* ipca */
    CONSTANTS[AC_IpCa_max] = 0.275;
    
    /* itr */
    CONSTANTS[AC_tau_tr] = 180.0;
    
    /* phys */
    CONSTANTS[AC_F] = 96.4867;
    CONSTANTS[AC_R] = 8.3143;
    CONSTANTS[AC_T] = 310.0;
    CONSTANTS[AC_RTF] = CONSTANTS[AC_R] * CONSTANTS[AC_T] / CONSTANTS[AC_F];
    CONSTANTS[AC_FRT] = 1.0 / CONSTANTS[AC_RTF];
    
    /* temp */
    CONSTANTS[AC_KQ10] = 3.0;
    
    /* inaca */
    CONSTANTS[AC_INaCa_max] = 1600.0;
    CONSTANTS[AC_KmCa] = 1.38;
    CONSTANTS[AC_KmNa] = 87.5;
    CONSTANTS[AC_g] = 0.35;
    CONSTANTS[AC_ksat] = 0.1;
    
    /* inak */
    CONSTANTS[AC_INaK_max] = 5.99338739999999981e-01;
    CONSTANTS[AC_KmKo] = 1.5;
    CONSTANTS[AC_KmNai] = 10.0;
    CONSTANTS[AC_sigma] = (exp(CONSTANTS[AC_Nao] / 67.3) - 1.0) / 7.0;
    
    // /* stimulus */
    // CONSTANTS[AC_amplitude] = 2.0 * -4618.0;
    // CONSTANTS[AC_stimulus_start]     = 50.0;     // Start time of pulse in ms
    // CONSTANTS[AC_stimulus_duration]  = 0.5;      // Pulse duration in ms
    // CONSTANTS[AC_stimulus_period]    = 1000;    // Time between pulses in ms

    
    /* cajsr */
    CONSTANTS[AC_K_rel] = 30.0;
    CONSTANTS[AC_c1] = 3.41749999999999983e-13;
    CONSTANTS[AC_c2] = 1.367e-15;
    CONSTANTS[AC_cajsr_u_tau] = 8.0;
    
    /* ib */
    CONSTANTS[AC_gbCa] = 0.001131;
    CONSTANTS[AC_gbNa] = 6.74437500000000015e-04;
    
    /* ik1 */
    CONSTANTS[AC_gK1] = 0.09;
    
    /* ikr */
    CONSTANTS[AC_gKr] = 2.94117649999999994e-02;
    
    /* iks */
    CONSTANTS[AC_gKs] = 1.29411759999999987e-01;
    
    /* ikur */
    CONSTANTS[AC_gKur_base] = 0.005;
    
    /* ina */
    CONSTANTS[AC_gNa] = 7.8;
    
    /* ito */
    CONSTANTS[AC_gto] = 0.1652;
    
    /* membrane */
    CONSTANTS[AC_I_diff] = 0.0;

    STATES[membrane_V] = -8.19463303822041098e+01;
    STATES[sodium_Nai] =  1.38169746305367962e+01;
    STATES[potassium_Ki] =  1.36355229902154434e+02;
    STATES[calcium_Cai] =  1.23092247890489894e-04;
    STATES[calcium_CaUp] =  1.54668119199095355e+00;
    STATES[calcium_CaRel] =  1.07650740580354909e+00;
    STATES[ina_m] =  2.56385228666526068e-03;
    STATES[ina_h] =  9.70298907063270155e-01;
    STATES[ina_j] =  9.81123905023234988e-01;
    STATES[ito_oa] =  2.91755626557170314e-02;
    STATES[ito_oi] =  9.99342865333055497e-01;
    STATES[ikur_ua] =  4.58838038240151104e-03;
    STATES[ikur_ui] =  9.91468962753066063e-01;
    STATES[ikr_xr] =  8.33819909884048389e-04;
    STATES[iks_xs] =  1.86683180787284714e-02;
    STATES[ical_d] =  1.24231529593716656e-04;
    STATES[ical_f] =  9.51907788168154578e-01;
    STATES[ical_fCa] =  7.39682838459564729e-01;
    STATES[cajsr_u] = -1.97647749727073971e-40;
    STATES[cajsr_v] = 1.0;
    STATES[cajsr_w] =  9.99233799248152699e-01;

}


void
CourtemanchecomputeVariables(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC,int tissueFlag, bool solveVmWithinODESolver)
{
    double AV_B1;
    double AV_B2;
    double AV_Ca_CMDN;
    double AV_Ca_CSQN;
    double AV_Ca_TRPN;
    double AV_ECa;
    double AV_EK;
    double AV_ENa;
    double AV_Fn;
    double AV_I_rel;
    double AV_I_tr;
    double AV_I_up;
    double AV_I_up_leak;
    double AV_cajsr_u_inf;
    double AV_cajsr_v_inf;
    double AV_cajsr_v_tau;
    double AV_cajsr_w_inf;
    double AV_cajsr_w_tau;
    double AV_fNaK;
    double AV_gKur;
    double AV_ical_d_inf;
    double AV_ical_d_tau;
    double AV_ical_fCa_inf;
    double AV_ical_f_inf;
    double AV_ical_f_tau;
    double AV_ikr_xr_inf;
    double AV_ikr_xr_tau;
    double AV_ikr_xr_tau_alpha;
    double AV_ikr_xr_tau_beta;
    double AV_iks_xs_inf;
    double AV_iks_xs_tau;
    double AV_iks_xs_tau_alpha;
    double AV_iks_xs_tau_beta;
    double AV_ikur_ua_inf;
    double AV_ikur_ua_tau;
    double AV_ikur_ua_tau_alpha;
    double AV_ikur_ua_tau_beta;
    double AV_ikur_ui_inf;
    double AV_ikur_ui_tau;
    double AV_ikur_ui_tau_alpha;
    double AV_ikur_ui_tau_beta;
    double AV_ina_h_alpha;
    double AV_ina_h_beta;
    double AV_ina_h_inf;
    double AV_ina_h_tau;
    double AV_ina_j_alpha;
    double AV_ina_j_beta;
    double AV_ina_j_inf;
    double AV_ina_j_tau;
    double AV_ina_m_alpha;
    double AV_ina_m_beta;
    double AV_ina_m_inf;
    double AV_ina_m_tau;
    double AV_ito_oa_alpha;
    double AV_ito_oa_beta;
    double AV_ito_oa_inf;
    double AV_ito_oa_tau;
    double AV_ito_oi_alpha;
    double AV_ito_oi_beta;
    double AV_ito_oi_inf;
    double AV_ito_oi_tau;
    /* ca_buffers */
    AV_Ca_CMDN = CONSTANTS[AC_CMDN_max] * STATES[calcium_Cai] / (STATES[calcium_Cai] + CONSTANTS[AC_Km_CMDN]);
    AV_Ca_CSQN = CONSTANTS[AC_CSQN_max] * STATES[calcium_CaRel] / (STATES[calcium_CaRel] + CONSTANTS[AC_Km_CSQN]);
    AV_Ca_TRPN = CONSTANTS[AC_TRPN_max] * STATES[calcium_Cai] / (STATES[calcium_Cai] + CONSTANTS[AC_Km_TRPN]);
    
    /* cansr */
    AV_I_up = CONSTANTS[AC_I_up_max] / (1.0 + CONSTANTS[AC_K_up] / STATES[calcium_Cai]);
    AV_I_up_leak = CONSTANTS[AC_I_up_max] * STATES[calcium_CaUp] / CONSTANTS[AC_Ca_up_max];
    
    //std::cout << "CURRENTS " << __LINE__ << std::endl;
    
    /* ical */
    AV_ical_d_inf = 1.0 / (1.0 + exp((STATES[membrane_V] + 10.0) / -8.0));
    AV_ical_d_tau = ((fabs(STATES[membrane_V] + 10.0) < 1e-06) ? 1.0 / (6.24 * 2.0 * 0.035) : (1.0 - exp((STATES[membrane_V] + 10.0) / -6.24)) / (0.035 * (STATES[membrane_V] + 10.0) * (1.0 + exp((STATES[membrane_V] + 10.0) / -6.24))));
    RATES[ical_d] = (AV_ical_d_inf - STATES[ical_d]) / AV_ical_d_tau;
    AV_ical_f_inf = 1.0 / (1.0 + exp((STATES[membrane_V] + 28.0) / 6.9));
    AV_ical_f_tau = 9.0 / (0.0197 * exp(-pow(0.0337, 2.0) * pow(STATES[membrane_V] + 10.0, 2.0)) + 0.02);
    RATES[ical_f] = (AV_ical_f_inf - STATES[ical_f]) / AV_ical_f_tau;
    AV_ical_fCa_inf = 1.0 / (1.0 + STATES[calcium_Cai] / 0.00035);
    RATES[ical_fCa] = (AV_ical_fCa_inf - STATES[ical_fCa]) / CONSTANTS[AC_ical_fCa_tau];
    ALGEBRAIC[AV_ICaL] = CONSTANTS[AC_gCaL] * STATES[ical_d] * STATES[ical_f] * STATES[ical_fCa] * (STATES[membrane_V] - CONSTANTS[AC_ECaL]);
    
    /* ipca */
    ALGEBRAIC[AV_IpCa] = CONSTANTS[AC_IpCa_max] * STATES[calcium_Cai] / (0.0005 + STATES[calcium_Cai]);
    
    /* itr */
    AV_I_tr = (STATES[calcium_CaUp] - STATES[calcium_CaRel]) / CONSTANTS[AC_tau_tr];
    
    /* inaca */
    ALGEBRAIC[AV_INaCa] = CONSTANTS[AC_INaCa_max] * (exp(CONSTANTS[AC_g] * STATES[membrane_V] * CONSTANTS[AC_FRT]) * pow(STATES[sodium_Nai], 3.0) * CONSTANTS[AC_Cao] - exp((CONSTANTS[AC_g] - 1.0) * STATES[membrane_V] * CONSTANTS[AC_FRT]) * pow(CONSTANTS[AC_Nao], 3.0) * STATES[calcium_Cai]) / ((pow(CONSTANTS[AC_KmNa], 3.0) + pow(CONSTANTS[AC_Nao], 3.0)) * (CONSTANTS[AC_KmCa] + CONSTANTS[AC_Cao]) * (1.0 + CONSTANTS[AC_ksat] * exp((CONSTANTS[AC_g] - 1.0) * STATES[membrane_V] * CONSTANTS[AC_FRT])));
    
    /* inak */
    AV_fNaK = 1.0 / (1.0 + 0.1245 * exp(-0.1 * STATES[membrane_V] * CONSTANTS[AC_FRT]) + 0.0365 * CONSTANTS[AC_sigma] * exp(-STATES[membrane_V] * CONSTANTS[AC_FRT]));
    ALGEBRAIC[AV_INaK] = CONSTANTS[AC_INaK_max] * AV_fNaK * CONSTANTS[AC_Ko] / (CONSTANTS[AC_Ko] + CONSTANTS[AC_KmKo]) / (1.0 + pow(CONSTANTS[AC_KmNai] / STATES[sodium_Nai], 1.5));
    
    /* nernst */
    AV_ECa = 0.5 * CONSTANTS[AC_RTF] * log(CONSTANTS[AC_Cao] / STATES[calcium_Cai]);
    AV_EK = CONSTANTS[AC_RTF] * log(CONSTANTS[AC_Ko] / STATES[potassium_Ki]);
    AV_ENa = CONSTANTS[AC_RTF] * log(CONSTANTS[AC_Nao] / STATES[sodium_Nai]);
    
    
    /* cajsr */
    AV_I_rel = CONSTANTS[AC_K_rel] * pow(STATES[cajsr_u], 2.0) * STATES[cajsr_v] * STATES[cajsr_w] * (STATES[calcium_CaRel] - STATES[calcium_Cai]);
    AV_cajsr_w_inf = 1.0 - 1.0 / (1.0 + exp(-(STATES[membrane_V] - 40.0) / 17.0));
    AV_cajsr_w_tau = 6.0 * ((fabs(STATES[membrane_V] - 7.9) < 1e-06) ? 2.0 / 13.0 : (1.0 - exp(-(STATES[membrane_V] - 7.9) / 5.0)) / ((1.0 + 0.3 * exp(-(STATES[membrane_V] - 7.9) / 5.0)) * (STATES[membrane_V] - 7.9)));
    RATES[cajsr_w] = (AV_cajsr_w_inf - STATES[cajsr_w]) / AV_cajsr_w_tau;
    AV_Fn = 1e-12 * CONSTANTS[AC_V_rel] * AV_I_rel - 5e-13 / CONSTANTS[AC_F] * (0.5 * ALGEBRAIC[AV_ICaL] - 0.2 * ALGEBRAIC[AV_INaCa]) * CONSTANTS[AC_Cm];
    AV_cajsr_u_inf = 1.0 / (1.0 + exp(-(AV_Fn - CONSTANTS[AC_c1]) / CONSTANTS[AC_c2]));
    RATES[cajsr_u] = (AV_cajsr_u_inf - STATES[cajsr_u]) / CONSTANTS[AC_cajsr_u_tau];
    AV_cajsr_v_inf = 1.0 - pow(1.0 + exp(-(AV_Fn - 0.2 * CONSTANTS[AC_c1]) / CONSTANTS[AC_c2]), -1.0);
    AV_cajsr_v_tau = 1.91 + 2.09 / (1.0 + exp(-(AV_Fn - CONSTANTS[AC_c1]) / CONSTANTS[AC_c2]));
    RATES[cajsr_v] = (AV_cajsr_v_inf - STATES[cajsr_v]) / AV_cajsr_v_tau;
    
    /* ib */
    ALGEBRAIC[AV_IbCa] = CONSTANTS[AC_gbCa] * (STATES[membrane_V] - AV_ECa);
    ALGEBRAIC[AV_IbNa] = CONSTANTS[AC_gbNa] * (STATES[membrane_V] - AV_ENa);
    
    /* ik1 */
    ALGEBRAIC[AV_IK1] = CONSTANTS[AC_gK1] * (STATES[membrane_V] - AV_EK) / (1.0 + exp(0.07 * (STATES[membrane_V] + 80.0)));
    
    /* ikr */
    ALGEBRAIC[AV_IKr] = CONSTANTS[AC_gKr] * STATES[ikr_xr] * (STATES[membrane_V] - AV_EK) / (1.0 + exp((STATES[membrane_V] + 15.0) / 22.4));
    AV_ikr_xr_inf = 1.0 / (1.0 + exp((STATES[membrane_V] + 14.1) / -6.5));
    AV_ikr_xr_tau_alpha = 0.0003 * ((fabs(STATES[membrane_V] + 14.1) < 1e-06) ? 5.0 : (STATES[membrane_V] + 14.1) / (1.0 - exp((STATES[membrane_V] + 14.1) / -5.0)));
    AV_ikr_xr_tau_beta = 7.38980000000000030e-05 * ((fabs(STATES[membrane_V] - 3.3328) < 1e-07) ? 5.1237 : (STATES[membrane_V] - 3.3328) / (exp((STATES[membrane_V] - 3.3328) / 5.1237) - 1.0));
    AV_ikr_xr_tau = 1.0 / (AV_ikr_xr_tau_alpha + AV_ikr_xr_tau_beta);
    RATES[ikr_xr] = (AV_ikr_xr_inf - STATES[ikr_xr]) / AV_ikr_xr_tau;
    
    /* iks */
    ALGEBRAIC[AV_IKs] = CONSTANTS[AC_gKs] * pow(STATES[iks_xs], 2.0) * (STATES[membrane_V] - AV_EK);
    AV_iks_xs_inf = 1.0 / sqrt(1.0 + exp((STATES[membrane_V] - 19.9) / -12.7));
    AV_iks_xs_tau_alpha = 4e-05 * ((fabs(STATES[membrane_V] - 19.9) < 1e-06) ? 17.0 : (STATES[membrane_V] - 19.9) / (1.0 - exp((STATES[membrane_V] - 19.9) / -17.0)));
    AV_iks_xs_tau_beta = 3.5e-05 * ((fabs(STATES[membrane_V] - 19.9) < 1e-06) ? 9.0 : (STATES[membrane_V] - 19.9) / (exp((STATES[membrane_V] - 19.9) / 9.0) - 1.0));
    AV_iks_xs_tau = 0.5 / (AV_iks_xs_tau_alpha + AV_iks_xs_tau_beta);
    RATES[iks_xs] = (AV_iks_xs_inf - STATES[iks_xs]) / AV_iks_xs_tau;
    
    /* ikur */
    AV_gKur = CONSTANTS[AC_gKur_base] * (1.0 + 10.0 / (1.0 + exp((STATES[membrane_V] - 15.0) / -13.0)));
    AV_ikur_ua_inf = 1.0 / (1.0 + exp((STATES[membrane_V] + 30.3) / -9.6));
    AV_ikur_ua_tau_alpha = 0.65 / (exp((STATES[membrane_V] + 10.0) / -8.5) + exp((STATES[membrane_V] - 30.0) / -59.0));
    AV_ikur_ua_tau_beta = 0.65 / (2.5 + exp((STATES[membrane_V] + 82.0) / 17.0));
    AV_ikur_ua_tau = 1.0 / (AV_ikur_ua_tau_alpha + AV_ikur_ua_tau_beta) / CONSTANTS[AC_KQ10];
    RATES[ikur_ua] = (AV_ikur_ua_inf - STATES[ikur_ua]) / AV_ikur_ua_tau;
    AV_ikur_ui_inf = 1.0 / (1.0 + exp((STATES[membrane_V] - 99.45) / 27.48));
    AV_ikur_ui_tau_alpha = 1.0 / (21.0 + exp((STATES[membrane_V] - 185.0) / -28.0));
    AV_ikur_ui_tau_beta = 1.0 / exp((STATES[membrane_V] - 158.0) / -16.0);
    AV_ikur_ui_tau = 1.0 / (AV_ikur_ui_tau_alpha + AV_ikur_ui_tau_beta) / CONSTANTS[AC_KQ10];
    RATES[ikur_ui] = (AV_ikur_ui_inf - STATES[ikur_ui]) / AV_ikur_ui_tau;
    ALGEBRAIC[AV_IKur] = AV_gKur * pow(STATES[ikur_ua], 3.0) * STATES[ikur_ui] * (STATES[membrane_V] - AV_EK);
    
    /* ina */
    ALGEBRAIC[AV_INa] = CONSTANTS[AC_gNa] * pow(STATES[ina_m], 3.0) * STATES[ina_h] * STATES[ina_j] * (STATES[membrane_V] - AV_ENa);
    AV_ina_h_alpha = ((STATES[membrane_V] < -40.0) ? 0.135 * exp((STATES[membrane_V] + 80.0) / -6.8) : 0.0);
    AV_ina_h_beta = ((STATES[membrane_V] < -40.0) ? 3.56 * exp(0.079 * STATES[membrane_V]) + 310000.0 * exp(0.35 * STATES[membrane_V]) : 1.0 / (0.13 * (1.0 + exp((STATES[membrane_V] + 10.66) / -11.1))));
    AV_ina_h_inf = AV_ina_h_alpha / (AV_ina_h_alpha + AV_ina_h_beta);
    AV_ina_h_tau = 1.0 / (AV_ina_h_alpha + AV_ina_h_beta);
    RATES[ina_h] = (AV_ina_h_inf - STATES[ina_h]) / AV_ina_h_tau;
    AV_ina_j_alpha = ((STATES[membrane_V] < -40.0) ? (-127140.0 * exp(0.2444 * STATES[membrane_V]) - 3.474e-05 * exp(-0.04391 * STATES[membrane_V])) * (STATES[membrane_V] + 37.78) / (1.0 + exp(0.311 * (STATES[membrane_V] + 79.23))) : 0.0);
    AV_ina_j_beta = ((STATES[membrane_V] < -40.0) ? 0.1212 * exp(-0.01052 * STATES[membrane_V]) / (1.0 + exp(-0.1378 * (STATES[membrane_V] + 40.14))) : 0.3 * exp(-2.535e-07 * STATES[membrane_V]) / (1.0 + exp(-0.1 * (STATES[membrane_V] + 32.0))));
    AV_ina_j_inf = AV_ina_j_alpha / (AV_ina_j_alpha + AV_ina_j_beta);
    AV_ina_j_tau = 1.0 / (AV_ina_j_alpha + AV_ina_j_beta);
    RATES[ina_j] = (AV_ina_j_inf - STATES[ina_j]) / AV_ina_j_tau;
    AV_ina_m_alpha = ((STATES[membrane_V] == -47.13) ? 3.2 : 0.32 * (STATES[membrane_V] + 47.13) / (1.0 - exp(-0.1 * (STATES[membrane_V] + 47.13))));
    AV_ina_m_beta = 0.08 * exp(-STATES[membrane_V] / 11.0);
    AV_ina_m_inf = AV_ina_m_alpha / (AV_ina_m_alpha + AV_ina_m_beta);
    AV_ina_m_tau = 1.0 / (AV_ina_m_alpha + AV_ina_m_beta);
    RATES[ina_m] = (AV_ina_m_inf - STATES[ina_m]) / AV_ina_m_tau;
    
    /* ito */
    ALGEBRAIC[AV_Ito] = CONSTANTS[AC_gto] * pow(STATES[ito_oa], 3.0) * STATES[ito_oi] * (STATES[membrane_V] - AV_EK);
    AV_ito_oa_alpha = 0.65 / (exp((STATES[membrane_V] + 10.0) / -8.5) + exp((STATES[membrane_V] - 30.0) / -59.0));
    AV_ito_oa_beta = 0.65 / (2.5 + exp((STATES[membrane_V] + 82.0) / 17.0));
    AV_ito_oa_inf = 1.0 / (1.0 + exp((STATES[membrane_V] + 20.47) / -17.54));
    AV_ito_oa_tau = 1.0 / (AV_ito_oa_alpha + AV_ito_oa_beta) / CONSTANTS[AC_KQ10];
    RATES[ito_oa] = (AV_ito_oa_inf - STATES[ito_oa]) / AV_ito_oa_tau;
    AV_ito_oi_alpha = 1.0 / (18.53 + exp((STATES[membrane_V] + 113.7) / 10.95));
    AV_ito_oi_beta = 1.0 / (35.56 + exp((STATES[membrane_V] + 1.26) / -7.44));
    AV_ito_oi_inf = 1.0 / (1.0 + exp((STATES[membrane_V] + 43.1) / 5.3));
    AV_ito_oi_tau = 1.0 / (AV_ito_oi_alpha + AV_ito_oi_beta) / CONSTANTS[AC_KQ10];
    RATES[ito_oi] = (AV_ito_oi_inf - STATES[ito_oi]) / AV_ito_oi_tau;
    
    /* calcium */
    RATES[calcium_CaRel] = (AV_I_tr - AV_I_rel) * pow(1.0 + CONSTANTS[AC_CSQN_max] * CONSTANTS[AC_Km_CSQN] / pow(STATES[calcium_CaRel] + CONSTANTS[AC_Km_CSQN], 2.0), -1.0);
    RATES[calcium_CaUp] = AV_I_up - (AV_I_up_leak + AV_I_tr * CONSTANTS[AC_V_rel] / CONSTANTS[AC_V_up]);
    AV_B1 = (2.0 * ALGEBRAIC[AV_INaCa] - (ALGEBRAIC[AV_IpCa] + ALGEBRAIC[AV_ICaL] + ALGEBRAIC[AV_IbCa])) * CONSTANTS[AC_Cm] / (2.0 * CONSTANTS[AC_V_i] * CONSTANTS[AC_F]) + (CONSTANTS[AC_V_up] * (AV_I_up_leak - AV_I_up) + AV_I_rel * CONSTANTS[AC_V_rel]) / CONSTANTS[AC_V_i];
    AV_B2 = 1.0 + CONSTANTS[AC_TRPN_max] * CONSTANTS[AC_Km_TRPN] / pow(STATES[calcium_Cai] + CONSTANTS[AC_Km_TRPN], 2.0) + CONSTANTS[AC_CMDN_max] * CONSTANTS[AC_Km_CMDN] / pow(STATES[calcium_Cai] + CONSTANTS[AC_Km_CMDN], 2.0);
    RATES[calcium_Cai] = AV_B1 / AV_B2;
    
    /* membrane */
    ALGEBRAIC[Iion_cm] = ALGEBRAIC[AV_INa] + ALGEBRAIC[AV_IK1] + ALGEBRAIC[AV_Ito] + ALGEBRAIC[AV_IKur] + ALGEBRAIC[AV_IKr] + ALGEBRAIC[AV_IKs] + ALGEBRAIC[AV_ICaL] + ALGEBRAIC[AV_IpCa] + ALGEBRAIC[AV_INaK] + ALGEBRAIC[AV_INaCa] + ALGEBRAIC[AV_IbNa] + ALGEBRAIC[AV_IbCa];
    //RATES[membrane_V] = -(ALGEBRAIC[AV_I_ion] + CONSTANTS[AC_I_diff] + ALGEBRAIC[AV_I_stim]);
    
    /* potassium */
    RATES[potassium_Ki] = (2.0 * ALGEBRAIC[AV_INaK] - (ALGEBRAIC[AV_IK1] + ALGEBRAIC[AV_Ito] + ALGEBRAIC[AV_IKur] + ALGEBRAIC[AV_IKr] + ALGEBRAIC[AV_IKs] + ALGEBRAIC[AV_I_stim])) * CONSTANTS[AC_Cm] / (CONSTANTS[AC_V_i] * CONSTANTS[AC_F]);
    
    /* sodium */
    RATES[sodium_Nai] = (-3.0 * ALGEBRAIC[AV_INaK] - (3.0 * ALGEBRAIC[AV_INaCa] + ALGEBRAIC[AV_IbNa] + ALGEBRAIC[AV_INa])) * CONSTANTS[AC_Cm] / (CONSTANTS[AC_V_i] * CONSTANTS[AC_F]);
    

    // S1–S2 STIMULUS LOGIC
    //std::cout << "CURRENTS " << __LINE__ << std::endl;
// ------------------------------------------------------------
// S1–S2 stimulus logic (Courtemanche version)
// ------------------------------------------------------------
if (solveVmWithinODESolver)
{
    double I_stim = 0.0;

    const double tStart  = CONSTANTS[stim_start];
    const double T1      = CONSTANTS[stim_period_S1];
    const double tDur    = CONSTANTS[stim_duration];
    const double Istim   = CONSTANTS[stim_amplitude];
    const int    n1      = static_cast<int>(CONSTANTS[nstim1]);

    const double T2      = CONSTANTS[stim_period_S2];
    const int    n2      = static_cast<int>(CONSTANTS[nstim2]);

    // ---------------------- S1 ----------------------
    if (T1 > 0.0 && n1 > 0)
    {
        // avoid VOI / 0
        double phase_S1 = VOI - floor(VOI/T1)*T1;

        bool in_S1 =
        (
            phase_S1 >= tStart &&
            phase_S1 <= tStart + tDur &&
            VOI <= tStart + T1*n1
        );

        if (in_S1)
        {
            I_stim = Istim;
        }
    }

    // ---------------------- S2 ----------------------
    if (I_stim == 0.0 && T2 > 0.0 && n2 > 0 && n1 > 0)
    {
        double tS1End = tStart + T1*n1;

        bool in_S2_total_window =
        (
            VOI >= tS1End &&
            VOI <= tS1End + n2*T2
        );

        if (in_S2_total_window)
        {
            double phase = fmod(VOI - tS1End, T2);

            if (phase <= tDur)
            {
                I_stim = Istim;
            }
        }
    }

    // Store algebraic external current
    ALGEBRAIC[AV_I_stim] = I_stim;

    // Vm equation
    RATES[membrane_V] = -(ALGEBRAIC[Iion_cm] + ALGEBRAIC[AV_I_stim]);

}
else
{
    RATES[membrane_V] = 0.0;
}

}

