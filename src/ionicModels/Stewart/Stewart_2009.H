#include "Stewart_2009Names.H"
#include "stimulusIO.H"

static const char* StewartSTATES_NAMES[NUM_STATES] = {
    "membrane_V",
    "Ihyperpolarization_activated_current_y_gate_y",
    "Irapid_time_dependent_potassium_current_Xr1_gate_Xr1",
    "Irapid_time_dependent_potassium_current_Xr2_gate_Xr2",
    "Islow_time_dependent_potassium_current_Xs_gate_Xs",
    "Ifast_sodium_current_m_gate_m",
    "Ifast_sodium_current_h_gate_h",
    "Ifast_sodium_current_j_gate_j",
    "IL_type_Ca_current_d_gate_d",
    "IL_type_Ca_current_f_gate_f",
    "IL_type_Ca_current_f2_gate_f2",
    "IL_type_Ca_current_fCass_gate_fCass",
    "Itransient_outward_current_s_gate_s",
    "Itransient_outward_current_r_gate_r",
    "Ca_i",
    "Ca_SR",
    "Ca_ss",
    "R_prime",
    "Na_i",
    "K_i",
};

static const char* StewartALGEBRAIC_NAMES[NUM_ALGEBRAIC] = {
    "AV_alpha_d",
    "AV_beta_d",
    "AV_d_inf",
    "AV_gamma_d",
    "AV_tau_d",
    "AV_f2_inf",
    "AV_tau_f2",
    "AV_fCass_inf",
    "AV_tau_fCass",
    "AV_f_inf",
    "AV_tau_f",
    "AV_i_p_Ca",
    "AV_alpha_h",
    "AV_beta_h",
    "AV_h_inf",
    "AV_tau_h",
    "AV_alpha_j",
    "AV_beta_j",
    "AV_j_inf",
    "AV_tau_j",
    "AV_alpha_m",
    "AV_beta_m",
    "AV_m_inf",
    "AV_tau_m",
    "AV_alpha_y",
    "AV_beta_y",
    "AV_y_inf",
    "AV_tau_y",
    "AV_alpha_xr1",
    "AV_beta_xr1",
    "AV_xr1_inf",
    "AV_tau_xr1",
    "AV_alpha_xr2",
    "AV_beta_xr2",
    "AV_xr2_inf",
    "AV_tau_xr2",
    "AV_alpha_xs",
    "AV_beta_xs",
    "AV_xs_inf",
    "AV_tau_xs",
    "AV_r_inf",
    "AV_tau_r",
    "AV_s_inf",
    "AV_tau_s",
    "AV_Ca_i_bufc",
    "AV_Ca_sr_bufsr",
    "AV_Ca_ss_bufss",
    "AV_i_leak",
    "AV_i_up",
    "AV_i_xfer",
    "AV_kcasr",
    "AV_k1",
    "AV_k2",
    "AV_O",
    "AV_i_rel",
    "AV_xK1_inf",
    "AV_E_Ca",
    "AV_E_K",
    "AV_i_NaK",
    "AV_a",
    "AV_i_sus",
    "AV_i_to",
    "AV_i_CaL",
    "AV_i_b_Ca",
    "AV_i_f_K",
    "AV_i_K1",
    "AV_i_p_K",
    "AV_i_Kr",
    "AV_E_Ks",
    "AV_E_Na",
    "AV_i_NaCa",
    "AV_i_Na",
    "AV_i_f_Na",
    "AV_i_Ks",
    "AV_i_b_Na",
    "AV_i_f",
    "Istim",
    "Iion_cm",
};


static inline const Foam::HashTable<Foam::wordList>&
StewartDependencyMap()
{
   static const Foam::HashTable<Foam::wordList> dep
    (
        {
        }
    );

    return dep;
}

inline Foam::scalar computeIstim(Foam::scalar t, const double* C)
{
    return Foam::stimulusIO::computeStimulus
    (
        t,
        C[stim_start],
        C[stim_period_S1],
        C[stim_duration],
        C[stim_amplitude],
        Foam::label(C[nstim1]),
        C[stim_period_S2],
        Foam::label(C[nstim2])
    );
}
void
StewartinitConsts(double* CONSTANTS, double* RATES, double* STATES, int tissueFlag, const Foam::dictionary& stimulus)
{
/* calcium_pump_current */
    CONSTANTS[AC_K_pCa] = 0.0005;
    CONSTANTS[AC_g_pCa] = 0.1238;

    /* L_type_Ca_current */
    CONSTANTS[AC_g_CaL] = 3.98e-05;

    /* calcium_background_current */
    CONSTANTS[AC_g_bca] = 0.000592;

    /* calcium_dynamics */
    CONSTANTS[AC_Buf_c] = 0.2;
    CONSTANTS[AC_Buf_sr] = 10.0;
    CONSTANTS[AC_Buf_ss] = 0.4;
    CONSTANTS[AC_Ca_o] = 2.0;
    CONSTANTS[AC_EC] = 1.5;
    CONSTANTS[AC_K_buf_c] = 0.001;
    CONSTANTS[AC_K_buf_sr] = 0.3;
    CONSTANTS[AC_K_buf_ss] = 0.00025;
    CONSTANTS[AC_K_up] = 0.00025;
    CONSTANTS[AC_V_leak] = 0.00036;
    CONSTANTS[AC_V_rel] = 0.102;
    CONSTANTS[AC_V_sr] = 0.001094;
    CONSTANTS[AC_V_ss] = 5.468e-05;
    CONSTANTS[AC_V_xfer] = 0.0038;
    CONSTANTS[AC_Vmax_up] = 0.006375;
    CONSTANTS[AC_k1_prime] = 0.15;
    CONSTANTS[AC_k2_prime] = 0.045;
    CONSTANTS[AC_k3] = 0.06;
    CONSTANTS[AC_k4] = 0.005;
    CONSTANTS[AC_max_sr] = 2.5;
    CONSTANTS[AC_min_sr] = 1.0;

    /* fast_sodium_current */
    CONSTANTS[AC_g_Na] = 130.5744;

    /* hyperpolarization_activated_current */
    CONSTANTS[AC_g_f_K] = 0.0234346;
    CONSTANTS[AC_g_f_Na] = 0.0145654;

    /* inward_rectifier_potassium_current */
    CONSTANTS[AC_g_K1] = 0.065;

    /* membrane */
    CONSTANTS[AC_Cm] = 0.185;
    CONSTANTS[AC_F] = 9.64853414999999950e+04;
    CONSTANTS[AC_R] = 8314.472;
    CONSTANTS[AC_T] = 310.0;
    CONSTANTS[AC_V_c] = 0.016404;

    /* potassium_dynamics */
    CONSTANTS[AC_K_o] = 5.4;

    /* potassium_pump_current */
    CONSTANTS[AC_g_pK] = 0.0146;

    /* rapid_time_dependent_potassium_current */
    CONSTANTS[AC_g_Kr] = 0.0918;

    /* reversal_potentials */
    CONSTANTS[AC_P_kna] = 0.03;

    /* slow_time_dependent_potassium_current */
    CONSTANTS[AC_g_Ks] = 0.2352;

    /* sodium_background_current */
    CONSTANTS[AC_g_bna] = 0.00029;

    /* sodium_calcium_exchanger_current */
    CONSTANTS[AC_K_NaCa] = 1000.0;
    CONSTANTS[AC_K_sat] = 0.1;
    CONSTANTS[AC_Km_Ca] = 1.38;
    CONSTANTS[AC_Km_Nai] = 87.5;
    CONSTANTS[AC_alpha] = 2.5;
    CONSTANTS[AC_sodium_calcium_exchanger_current_gamma] = 0.35;

    /* sodium_dynamics */
    CONSTANTS[AC_Na_o] = 140.0;

    /* sodium_potassium_pump_current */
    CONSTANTS[AC_K_mNa] = 40.0;
    CONSTANTS[AC_K_mk] = 1.0;
    CONSTANTS[AC_P_NaK] = 2.724;

    /* sustained_outward_current */
    CONSTANTS[AC_g_sus] = 0.0227;

    /* transient_outward_current */
    CONSTANTS[AC_g_to] = 0.08184;

    /* --- Initial values --- */

    STATES[membrane_V] = -6.91370441635924067e+01;
    STATES[Ihyperpolarization_activated_current_y_gate_y] =  4.57562667986601973e-02;
    STATES[Irapid_time_dependent_potassium_current_Xr1_gate_Xr1] =  5.50281999719088016e-03;
    STATES[Irapid_time_dependent_potassium_current_Xr2_gate_Xr2] =  3.13213286437995009e-01;
    STATES[Islow_time_dependent_potassium_current_Xs_gate_Xs] =  9.53708522974788982e-03;
    STATES[Ifast_sodium_current_m_gate_m] =  4.17391656294996971e-02;
    STATES[Ifast_sodium_current_h_gate_h] =  1.90678733735144990e-01;
    STATES[Ifast_sodium_current_j_gate_j] =  2.38219836154029002e-01;
    STATES[IL_type_Ca_current_d_gate_d] =  2.87906256206415015e-04;
    STATES[IL_type_Ca_current_f_gate_f] =  9.89328560287986991e-01;
    STATES[IL_type_Ca_current_f2_gate_f2] =  9.95474890442184956e-01;
    STATES[IL_type_Ca_current_fCass_gate_fCass] =  9.99955429598213041e-01;
    STATES[Itransient_outward_current_s_gate_s] =  9.63861017995010005e-01;
    STATES[Itransient_outward_current_r_gate_r] =  1.03618091196912006e-03;
    STATES[Ca_i] =  1.01878186157052002e-04;
    STATES[Ca_SR] =  3.10836886659417022e+00;
    STATES[Ca_ss] =  4.46818714055411022e-04;
    STATES[R_prime] =  9.91580051907845039e-01;
    STATES[Na_i] =  8.80420286531673035e+00;
    STATES[K_i] =  1.36781894160227012e+02;
}


void
StewartcomputeVariables(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC, int tissueFlag, bool solveVmWithinODESolver)
{
/* L_type_Ca_current_d_gate */
    ALGEBRAIC[AV_alpha_d] = 1.4 / (1.0 + exp((-35.0 - STATES[membrane_V]) / 13.0)) + 0.25;
    ALGEBRAIC[AV_beta_d] = 1.4 / (1.0 + exp((STATES[membrane_V] + 5.0) / 5.0));
    ALGEBRAIC[AV_d_inf] = 1.0 / (1.0 + exp((-8.0 - STATES[membrane_V]) / 7.5));
    ALGEBRAIC[AV_gamma_d] = 1.0 / (1.0 + exp((50.0 - STATES[membrane_V]) / 20.0));
    ALGEBRAIC[AV_tau_d] = 1.0 * ALGEBRAIC[AV_alpha_d] * ALGEBRAIC[AV_beta_d] + ALGEBRAIC[AV_gamma_d];
    RATES[IL_type_Ca_current_d_gate_d] = (ALGEBRAIC[AV_d_inf] - STATES[IL_type_Ca_current_d_gate_d]) / ALGEBRAIC[AV_tau_d];

    /* L_type_Ca_current_f2_gate */
    ALGEBRAIC[AV_f2_inf] = 0.67 / (1.0 + exp((STATES[membrane_V] + 35.0) / 7.0)) + 0.33;
    ALGEBRAIC[AV_tau_f2] = 562.0 * exp(-pow(STATES[membrane_V] + 27.0, 2.0) / 240.0) + 31.0 / (1.0 + exp((25.0 - STATES[membrane_V]) / 10.0)) + 80.0 / (1.0 + exp((STATES[membrane_V] + 30.0) / 10.0));
    RATES[IL_type_Ca_current_f2_gate_f2] = (ALGEBRAIC[AV_f2_inf] - STATES[IL_type_Ca_current_f2_gate_f2]) / ALGEBRAIC[AV_tau_f2];

    /* L_type_Ca_current_fCass_gate */
    ALGEBRAIC[AV_fCass_inf] = 0.6 / (1.0 + pow(STATES[Ca_ss] / 0.05, 2.0)) + 0.4;
    ALGEBRAIC[AV_tau_fCass] = 80.0 / (1.0 + pow(STATES[Ca_ss] / 0.05, 2.0)) + 2.0;
    RATES[IL_type_Ca_current_fCass_gate_fCass] = (ALGEBRAIC[AV_fCass_inf] - STATES[IL_type_Ca_current_fCass_gate_fCass]) / ALGEBRAIC[AV_tau_fCass];

    /* L_type_Ca_current_f_gate */
    ALGEBRAIC[AV_f_inf] = 1.0 / (1.0 + exp((STATES[membrane_V] + 20.0) / 7.0));
    ALGEBRAIC[AV_tau_f] = 1102.5 * exp(-pow(STATES[membrane_V] + 27.0, 2.0) / 225.0) + 200.0 / (1.0 + exp((13.0 - STATES[membrane_V]) / 10.0)) + 180.0 / (1.0 + exp((STATES[membrane_V] + 30.0) / 10.0)) + 20.0;
    RATES[IL_type_Ca_current_f_gate_f] = (ALGEBRAIC[AV_f_inf] - STATES[IL_type_Ca_current_f_gate_f]) / ALGEBRAIC[AV_tau_f];

    /* calcium_pump_current */
    ALGEBRAIC[AV_i_p_Ca] = CONSTANTS[AC_g_pCa] * STATES[Ca_i] / (STATES[Ca_i] + CONSTANTS[AC_K_pCa]);


    /* fast_sodium_current_h_gate */
    ALGEBRAIC[AV_alpha_h] = ((STATES[membrane_V] < -40.0) ? 0.057 * exp(-(STATES[membrane_V] + 80.0) / 6.8) : 0.0);
    ALGEBRAIC[AV_beta_h] = ((STATES[membrane_V] < -40.0) ? 2.7 * exp(0.079 * STATES[membrane_V]) + 310000.0 * exp(0.3485 * STATES[membrane_V]) : 0.77 / (0.13 * (1.0 + exp((STATES[membrane_V] + 10.66) / -11.1))));
    ALGEBRAIC[AV_h_inf] = 1.0 / pow(1.0 + exp((STATES[membrane_V] + 71.55) / 7.43), 2.0);
    ALGEBRAIC[AV_tau_h] = 1.0 / (ALGEBRAIC[AV_alpha_h] + ALGEBRAIC[AV_beta_h]);
    RATES[Ifast_sodium_current_h_gate_h] = (ALGEBRAIC[AV_h_inf] - STATES[Ifast_sodium_current_h_gate_h]) / ALGEBRAIC[AV_tau_h];

    /* fast_sodium_current_j_gate */
    ALGEBRAIC[AV_alpha_j] = ((STATES[membrane_V] < -40.0) ? (-25428.0 * exp(0.2444 * STATES[membrane_V]) - 6.948e-06 * exp(-0.04391 * STATES[membrane_V])) * (STATES[membrane_V] + 37.78) / 1.0 / (1.0 + exp(0.311 * (STATES[membrane_V] + 79.23))) : 0.0);
    ALGEBRAIC[AV_beta_j] = ((STATES[membrane_V] < -40.0) ? 0.02424 * exp(-0.01052 * STATES[membrane_V]) / (1.0 + exp(-0.1378 * (STATES[membrane_V] + 40.14))) : 0.6 * exp(0.057 * STATES[membrane_V]) / (1.0 + exp(-0.1 * (STATES[membrane_V] + 32.0))));
    ALGEBRAIC[AV_j_inf] = 1.0 / pow(1.0 + exp((STATES[membrane_V] + 71.55) / 7.43), 2.0);
    ALGEBRAIC[AV_tau_j] = 1.0 / (ALGEBRAIC[AV_alpha_j] + ALGEBRAIC[AV_beta_j]);
    RATES[Ifast_sodium_current_j_gate_j] = (ALGEBRAIC[AV_j_inf] - STATES[Ifast_sodium_current_j_gate_j]) / ALGEBRAIC[AV_tau_j];

    /* fast_sodium_current_m_gate */
    ALGEBRAIC[AV_alpha_m] = 1.0 / (1.0 + exp((-60.0 - STATES[membrane_V]) / 5.0));
    ALGEBRAIC[AV_beta_m] = 0.1 / (1.0 + exp((STATES[membrane_V] + 35.0) / 5.0)) + 0.1 / (1.0 + exp((STATES[membrane_V] - 50.0) / 200.0));
    ALGEBRAIC[AV_m_inf] = 1.0 / pow(1.0 + exp((-56.86 - STATES[membrane_V]) / 9.03), 2.0);
    ALGEBRAIC[AV_tau_m] = 1.0 * ALGEBRAIC[AV_alpha_m] * ALGEBRAIC[AV_beta_m];
    RATES[Ifast_sodium_current_m_gate_m] = (ALGEBRAIC[AV_m_inf] - STATES[Ifast_sodium_current_m_gate_m]) / ALGEBRAIC[AV_tau_m];

    /* hyperpolarization_activated_current_y_gate */
    ALGEBRAIC[AV_alpha_y] = 1.0 * exp(-2.9 - 0.04 * STATES[membrane_V]);
    ALGEBRAIC[AV_beta_y] = 1.0 * exp(3.6 + 0.11 * STATES[membrane_V]);
    ALGEBRAIC[AV_y_inf] = 1.0 / (1.0 + exp((STATES[membrane_V] + 80.6) / 6.8));
    ALGEBRAIC[AV_tau_y] = 4000.0 / (ALGEBRAIC[AV_alpha_y] + ALGEBRAIC[AV_beta_y]);
    RATES[Ihyperpolarization_activated_current_y_gate_y] = (ALGEBRAIC[AV_y_inf] - STATES[Ihyperpolarization_activated_current_y_gate_y]) / ALGEBRAIC[AV_tau_y];

    /* rapid_time_dependent_potassium_current_Xr1_gate */
    ALGEBRAIC[AV_alpha_xr1] = 450.0 / (1.0 + exp((-45.0 - STATES[membrane_V]) / 10.0));
    ALGEBRAIC[AV_beta_xr1] = 6.0 / (1.0 + exp((STATES[membrane_V] + 30.0) / 11.5));
    ALGEBRAIC[AV_xr1_inf] = 1.0 / (1.0 + exp((-26.0 - STATES[membrane_V]) / 7.0));
    ALGEBRAIC[AV_tau_xr1] = 1.0 * ALGEBRAIC[AV_alpha_xr1] * ALGEBRAIC[AV_beta_xr1];
    RATES[Irapid_time_dependent_potassium_current_Xr1_gate_Xr1] = (ALGEBRAIC[AV_xr1_inf] - STATES[Irapid_time_dependent_potassium_current_Xr1_gate_Xr1]) / ALGEBRAIC[AV_tau_xr1];

    /* rapid_time_dependent_potassium_current_Xr2_gate */
    ALGEBRAIC[AV_alpha_xr2] = 3.0 / (1.0 + exp((-60.0 - STATES[membrane_V]) / 20.0));
    ALGEBRAIC[AV_beta_xr2] = 1.12 / (1.0 + exp((STATES[membrane_V] - 60.0) / 20.0));
    ALGEBRAIC[AV_xr2_inf] = 1.0 / (1.0 + exp((STATES[membrane_V] + 88.0) / 24.0));
    ALGEBRAIC[AV_tau_xr2] = 1.0 * ALGEBRAIC[AV_alpha_xr2] * ALGEBRAIC[AV_beta_xr2];
    RATES[Irapid_time_dependent_potassium_current_Xr2_gate_Xr2] = (ALGEBRAIC[AV_xr2_inf] - STATES[Irapid_time_dependent_potassium_current_Xr2_gate_Xr2]) / ALGEBRAIC[AV_tau_xr2];

    /* slow_time_dependent_potassium_current_Xs_gate */
    ALGEBRAIC[AV_alpha_xs] = 1400.0 / sqrt(1.0 + exp((5.0 - STATES[membrane_V]) / 6.0));
    ALGEBRAIC[AV_beta_xs] = 1.0 / (1.0 + exp((STATES[membrane_V] - 35.0) / 15.0));
    ALGEBRAIC[AV_xs_inf] = 1.0 / (1.0 + exp((-5.0 - STATES[membrane_V]) / 14.0));
    ALGEBRAIC[AV_tau_xs] = 1.0 * ALGEBRAIC[AV_alpha_xs] * ALGEBRAIC[AV_beta_xs] + 80.0;
    RATES[Islow_time_dependent_potassium_current_Xs_gate_Xs] = (ALGEBRAIC[AV_xs_inf] - STATES[Islow_time_dependent_potassium_current_Xs_gate_Xs]) / ALGEBRAIC[AV_tau_xs];

    /* transient_outward_current_r_gate */
    ALGEBRAIC[AV_r_inf] = 1.0 / (1.0 + exp((20.0 - STATES[membrane_V]) / 13.0));
    ALGEBRAIC[AV_tau_r] = 10.45 * exp(-pow(STATES[membrane_V] + 40.0, 2.0) / 1800.0) + 7.3;
    RATES[Itransient_outward_current_r_gate_r] = (ALGEBRAIC[AV_r_inf] - STATES[Itransient_outward_current_r_gate_r]) / ALGEBRAIC[AV_tau_r];

    /* transient_outward_current_s_gate */
    ALGEBRAIC[AV_s_inf] = 1.0 / (1.0 + exp((STATES[membrane_V] + 27.0) / 13.0));
    ALGEBRAIC[AV_tau_s] = 85.0 * exp(-pow(STATES[membrane_V] + 25.0, 2.0) / 320.0) + 5.0 / (1.0 + exp((STATES[membrane_V] - 40.0) / 5.0)) + 42.0;
    RATES[Itransient_outward_current_s_gate_s] = (ALGEBRAIC[AV_s_inf] - STATES[Itransient_outward_current_s_gate_s]) / ALGEBRAIC[AV_tau_s];

    /* calcium_dynamics */
    ALGEBRAIC[AV_Ca_i_bufc] = 1.0 / (1.0 + CONSTANTS[AC_Buf_c] * CONSTANTS[AC_K_buf_c] / pow(STATES[Ca_i] + CONSTANTS[AC_K_buf_c], 2.0));
    ALGEBRAIC[AV_Ca_sr_bufsr] = 1.0 / (1.0 + CONSTANTS[AC_Buf_sr] * CONSTANTS[AC_K_buf_sr] / pow(STATES[Ca_SR] + CONSTANTS[AC_K_buf_sr], 2.0));
    ALGEBRAIC[AV_Ca_ss_bufss] = 1.0 / (1.0 + CONSTANTS[AC_Buf_ss] * CONSTANTS[AC_K_buf_ss] / pow(STATES[Ca_ss] + CONSTANTS[AC_K_buf_ss], 2.0));
    ALGEBRAIC[AV_i_leak] = CONSTANTS[AC_V_leak] * (STATES[Ca_SR] - STATES[Ca_i]);
    ALGEBRAIC[AV_i_up] = CONSTANTS[AC_Vmax_up] / (1.0 + pow(CONSTANTS[AC_K_up], 2.0) / pow(STATES[Ca_i], 2.0));
    ALGEBRAIC[AV_i_xfer] = CONSTANTS[AC_V_xfer] * (STATES[Ca_ss] - STATES[Ca_i]);
    ALGEBRAIC[AV_kcasr] = CONSTANTS[AC_max_sr] - (CONSTANTS[AC_max_sr] - CONSTANTS[AC_min_sr]) / (1.0 + pow(CONSTANTS[AC_EC] / STATES[Ca_SR], 2.0));
    ALGEBRAIC[AV_k1] = CONSTANTS[AC_k1_prime] / ALGEBRAIC[AV_kcasr];
    ALGEBRAIC[AV_k2] = CONSTANTS[AC_k2_prime] * ALGEBRAIC[AV_kcasr];
    ALGEBRAIC[AV_O] = ALGEBRAIC[AV_k1] * pow(STATES[Ca_ss], 2.0) * STATES[R_prime] / (CONSTANTS[AC_k3] + ALGEBRAIC[AV_k1] * pow(STATES[Ca_ss], 2.0));
    RATES[R_prime] = -ALGEBRAIC[AV_k2] * STATES[Ca_ss] * STATES[R_prime] + CONSTANTS[AC_k4] * (1.0 - STATES[R_prime]);
    ALGEBRAIC[AV_i_rel] = CONSTANTS[AC_V_rel] * ALGEBRAIC[AV_O] * (STATES[Ca_SR] - STATES[Ca_ss]);
    RATES[Ca_SR] = ALGEBRAIC[AV_Ca_sr_bufsr] * (ALGEBRAIC[AV_i_up] - (ALGEBRAIC[AV_i_rel] + ALGEBRAIC[AV_i_leak]));

    /* inward_rectifier_potassium_current */
    ALGEBRAIC[AV_xK1_inf] = 1.0 / (1.0 + exp(0.1 * (STATES[membrane_V] + 75.44)));

    /* reversal_potentials */
    ALGEBRAIC[AV_E_Ca] = 0.5 * CONSTANTS[AC_R] * CONSTANTS[AC_T] / CONSTANTS[AC_F] * log(CONSTANTS[AC_Ca_o] / STATES[Ca_i]);
    ALGEBRAIC[AV_E_K] = CONSTANTS[AC_R] * CONSTANTS[AC_T] / CONSTANTS[AC_F] * log(CONSTANTS[AC_K_o] / STATES[K_i]);

    /* sodium_potassium_pump_current */
    ALGEBRAIC[AV_i_NaK] = CONSTANTS[AC_P_NaK] * CONSTANTS[AC_K_o] / (CONSTANTS[AC_K_o] + CONSTANTS[AC_K_mk]) * STATES[Na_i] / (STATES[Na_i] + CONSTANTS[AC_K_mNa]) / (1.0 + 0.1245 * exp(-0.1 * STATES[membrane_V] * CONSTANTS[AC_F] / (CONSTANTS[AC_R] * CONSTANTS[AC_T])) + 0.0353 * exp(-STATES[membrane_V] * CONSTANTS[AC_F] / (CONSTANTS[AC_R] * CONSTANTS[AC_T])));

    /* sustained_outward_current */
    ALGEBRAIC[AV_a] = 1.0 / (1.0 + exp((5.0 - STATES[membrane_V]) / 17.0));
    ALGEBRAIC[AV_i_sus] = CONSTANTS[AC_g_sus] * ALGEBRAIC[AV_a] * (STATES[membrane_V] - ALGEBRAIC[AV_E_K]);

    /* transient_outward_current */
    ALGEBRAIC[AV_i_to] = CONSTANTS[AC_g_to] * STATES[Itransient_outward_current_r_gate_r] * STATES[Itransient_outward_current_s_gate_s] * (STATES[membrane_V] - ALGEBRAIC[AV_E_K]);

    /* *remaining* */
    ALGEBRAIC[AV_i_CaL] = CONSTANTS[AC_g_CaL] * STATES[IL_type_Ca_current_d_gate_d] * STATES[IL_type_Ca_current_f_gate_f] * STATES[IL_type_Ca_current_f2_gate_f2] * STATES[IL_type_Ca_current_fCass_gate_fCass] * 4.0 * (STATES[membrane_V] - 15.0) * pow(CONSTANTS[AC_F],
                                                                                                                      2.0) / (CONSTANTS[AC_R] * CONSTANTS[AC_T]) * (0.25 * STATES[Ca_ss] * exp(2.0 * (STATES[membrane_V] - 15.0) * CONSTANTS[AC_F] / (CONSTANTS[AC_R] * CONSTANTS[AC_T])) - CONSTANTS[AC_Ca_o]) / (exp(2.0 * (STATES[membrane_V] - 15.0) * CONSTANTS[AC_F] / (CONSTANTS[AC_R] * CONSTANTS[AC_T])) - 1.0);
    ALGEBRAIC[AV_i_b_Ca] = CONSTANTS[AC_g_bca] * (STATES[membrane_V] - ALGEBRAIC[AV_E_Ca]);
    ALGEBRAIC[AV_i_f_K] = STATES[Ihyperpolarization_activated_current_y_gate_y] * CONSTANTS[AC_g_f_K] * (STATES[membrane_V] - ALGEBRAIC[AV_E_K]);
    ALGEBRAIC[AV_i_K1] = CONSTANTS[AC_g_K1] * ALGEBRAIC[AV_xK1_inf] * (STATES[membrane_V] - 8.0 - ALGEBRAIC[AV_E_K]);
    ALGEBRAIC[AV_i_p_K] = CONSTANTS[AC_g_pK] * (STATES[membrane_V] - ALGEBRAIC[AV_E_K]) / (1.0 + exp((25.0 - STATES[membrane_V]) / 5.98));
    ALGEBRAIC[AV_i_Kr] = CONSTANTS[AC_g_Kr] * sqrt(CONSTANTS[AC_K_o] / 5.4) * STATES[Irapid_time_dependent_potassium_current_Xr1_gate_Xr1] * STATES[Irapid_time_dependent_potassium_current_Xr2_gate_Xr2] * (STATES[membrane_V] - ALGEBRAIC[AV_E_K]);
    ALGEBRAIC[AV_E_Ks] = CONSTANTS[AC_R] * CONSTANTS[AC_T] / CONSTANTS[AC_F] * log((CONSTANTS[AC_K_o] + CONSTANTS[AC_P_kna] * CONSTANTS[AC_Na_o]) / (STATES[K_i] + CONSTANTS[AC_P_kna] * STATES[Na_i]));
    ALGEBRAIC[AV_E_Na] = CONSTANTS[AC_R] * CONSTANTS[AC_T] / CONSTANTS[AC_F] * log(CONSTANTS[AC_Na_o] / STATES[Na_i]);
    ALGEBRAIC[AV_i_NaCa] = CONSTANTS[AC_K_NaCa] * (exp(CONSTANTS[AC_sodium_calcium_exchanger_current_gamma] * STATES[membrane_V] * CONSTANTS[AC_F] / (CONSTANTS[AC_R] * CONSTANTS[AC_T])) * pow(STATES[Na_i], 3.0) * CONSTANTS[AC_Ca_o] - exp((CONSTANTS[AC_sodium_calcium_exchanger_current_gamma] - 1.0) * STATES[membrane_V] * CONSTANTS[AC_F] / (CONSTANTS[AC_R] * CONSTANTS[AC_T])) * pow(CONSTANTS[AC_Na_o], 3.0) * STATES[Ca_i] * CONSTANTS[AC_alpha]) / ((pow(CONSTANTS[AC_Km_Nai], 3.0) + pow(CONSTANTS[AC_Na_o], 3.0)) * (CONSTANTS[AC_Km_Ca] + CONSTANTS[AC_Ca_o]) * (1.0 + CONSTANTS[AC_K_sat] * exp((CONSTANTS[AC_sodium_calcium_exchanger_current_gamma] - 1.0) * STATES[membrane_V] * CONSTANTS[AC_F] / (CONSTANTS[AC_R] * CONSTANTS[AC_T]))));
    RATES[Ca_i] = ALGEBRAIC[AV_Ca_i_bufc] * ((ALGEBRAIC[AV_i_leak] - ALGEBRAIC[AV_i_up]) * CONSTANTS[AC_V_sr] / CONSTANTS[AC_V_c] + ALGEBRAIC[AV_i_xfer] - 1.0 * (ALGEBRAIC[AV_i_b_Ca] + ALGEBRAIC[AV_i_p_Ca] - 2.0 * ALGEBRAIC[AV_i_NaCa]) * CONSTANTS[AC_Cm] / (2.0 * 1.0 * CONSTANTS[AC_V_c] * CONSTANTS[AC_F]));
    RATES[Ca_ss] = ALGEBRAIC[AV_Ca_ss_bufss] * (-1.0 * ALGEBRAIC[AV_i_CaL] * CONSTANTS[AC_Cm] / (2.0 * 1.0 * CONSTANTS[AC_V_ss] * CONSTANTS[AC_F]) + ALGEBRAIC[AV_i_rel] * CONSTANTS[AC_V_sr] / CONSTANTS[AC_V_ss] - ALGEBRAIC[AV_i_xfer] * CONSTANTS[AC_V_c] / CONSTANTS[AC_V_ss]);
    ALGEBRAIC[AV_i_Na] = CONSTANTS[AC_g_Na] * pow(STATES[Ifast_sodium_current_m_gate_m], 3.0) * STATES[Ifast_sodium_current_h_gate_h] * STATES[Ifast_sodium_current_j_gate_j] * (STATES[membrane_V] - ALGEBRAIC[AV_E_Na]);
    ALGEBRAIC[AV_i_f_Na] = STATES[Ihyperpolarization_activated_current_y_gate_y] * CONSTANTS[AC_g_f_Na] * (STATES[membrane_V] - ALGEBRAIC[AV_E_Na]);
    ALGEBRAIC[AV_i_Ks] = CONSTANTS[AC_g_Ks] * pow(STATES[Islow_time_dependent_potassium_current_Xs_gate_Xs], 2.0) * (STATES[membrane_V] - ALGEBRAIC[AV_E_Ks]);
    ALGEBRAIC[AV_i_b_Na] = CONSTANTS[AC_g_bna] * (STATES[membrane_V] - ALGEBRAIC[AV_E_Na]);
    ALGEBRAIC[AV_i_f] = ALGEBRAIC[AV_i_f_Na] + ALGEBRAIC[AV_i_f_K];
    RATES[K_i] = -1.0 * (ALGEBRAIC[AV_i_K1] + ALGEBRAIC[AV_i_to] + ALGEBRAIC[AV_i_f_K] + ALGEBRAIC[AV_i_sus] + ALGEBRAIC[AV_i_Kr] + ALGEBRAIC[AV_i_Ks] + ALGEBRAIC[AV_i_p_K] - 2.0 * ALGEBRAIC[AV_i_NaK]) / (1.0 * CONSTANTS[AC_V_c] * CONSTANTS[AC_F]) * CONSTANTS[AC_Cm];
    RATES[Na_i] = -1.0 * (ALGEBRAIC[AV_i_Na] + ALGEBRAIC[AV_i_b_Na] + ALGEBRAIC[AV_i_f_Na] + 3.0 * ALGEBRAIC[AV_i_NaK] + 3.0 * ALGEBRAIC[AV_i_NaCa]) / (1.0 * CONSTANTS[AC_V_c] * CONSTANTS[AC_F]) * CONSTANTS[AC_Cm];
    ALGEBRAIC[Iion_cm]=(ALGEBRAIC[AV_i_K1] + ALGEBRAIC[AV_i_to] + ALGEBRAIC[AV_i_sus] + ALGEBRAIC[AV_i_Kr] + ALGEBRAIC[AV_i_Ks] + ALGEBRAIC[AV_i_CaL] + ALGEBRAIC[AV_i_NaK] + ALGEBRAIC[AV_i_Na] + ALGEBRAIC[AV_i_b_Na] + ALGEBRAIC[AV_i_NaCa] + ALGEBRAIC[AV_i_b_Ca] + ALGEBRAIC[AV_i_p_K] + ALGEBRAIC[AV_i_p_Ca] + ALGEBRAIC[AV_i_f]);


    ALGEBRAIC[Istim] = computeIstim(VOI, CONSTANTS);

    if (solveVmWithinODESolver)
    {
        RATES[membrane_V] = -ALGEBRAIC[Iion_cm] - ALGEBRAIC[Istim];
    }

}
