#include "Fabbri_2017Names.H"
#include "stimulusIO.H"

static const char* FabbriSTATES_NAMES[NUM_STATES] = {
    "membrane_V",
    "Na_i",
    "If_y_gate_y",
    "INa_m_gate_m",
    "INa_h_gate_h",
    "ICaL_dL_gate_dL",
    "ICaL_fL_gate_fL",
    "ICaL_fCa_gate_fCa",
    "ICaT_dT_gate_dT",
    "ICaT_fT_gate_fT",
    "SR_R",
    "SR_O",
    "SR_I",
    "SR_RI",
    "Buffer_fTMM",
    "Buffer_fCMi",
    "Buffer_fCMs",
    "Buffer_fTC",
    "Buffer_fTMC",
    "Buffer_fCQ",
    "Ca_i",
    "Ca_nsr",
    "Ca_jsr",
    "Ca_sub",
    "IKur_rKur_gate_r_Kur",
    "IKur_sKur_gate_s_Kur",
    "Ito_q_gate_q",
    "Ito_r_gate_r",
    "IKr_pa_gate_paS",
    "IKr_pa_gate_paF",
    "IKr_pi_gate_piy",
    "IKs_n_gate_n",
    "IKACh_a_gate_a",
};

static const char* FabbriALGEBRAIC_NAMES[NUM_ALGEBRAIC] = {
    "AV_P_tot",
    "AV_diff",
    "AV_j_SRCarel",
    "AV_kCaSR",
    "AV_kiSRCa",
    "AV_koSRCa",
    "AV_delta_fCMi",
    "AV_delta_fCMs",
    "AV_delta_fCQ",
    "AV_delta_fTC",
    "AV_delta_fTMC",
    "AV_delta_fTMM",
    "AV_fCa_infinity",
    "AV_tau_fCa",
    "AV_j_Ca_dif",
    "AV_j_tr",
    "AV_j_up",
    "AV_V_clamp",
    "AV_V",
    "AV_Nai",
    "AV_i_siCa",
    "AV_i_siK",
    "AV_i_siNa",
    "AV_i_CaL_i_CaL",
    "AV_adVm",
    "AV_bdVm",
    "AV_alpha_dL",
    "AV_beta_dL",
    "AV_dL_infinity",
    "AV_tau_dL",
    "AV_tau_fL",
    "AV_fL_infinity",
    "AV_i_CaT_i_CaT",
    "AV_dT_infinity",
    "AV_tau_dT",
    "AV_fT_infinity",
    "AV_tau_fT",
    "AV_beta_a",
    "AV_a_infinity",
    "AV_tau_a",
    "AV_alfapaF",
    "AV_betapaF",
    "AV_pa_infinity",
    "AV_tau_paF",
    "AV_tau_paS",
    "AV_pi_infinity",
    "AV_tau_pi",
    "AV_E_Ks",
    "AV_i_Ks_i_Ks",
    "AV_alpha_n",
    "AV_beta_n",
    "AV_n_infinity",
    "AV_tau_n",
    "AV_r_Kur_infinity",
    "AV_tau_r_Kur",
    "AV_s_Kur_infinity",
    "AV_tau_s_Kur",
    "AV_E_mh",
    "AV_i_Na_",
    "AV_i_Na_L",
    "AV_i_Na_i_Na",
    "AV_di",
    "AV_i_NaCa_do",
    "AV_k32",
    "AV_k41",
    "AV_k43",
    "AV_k12",
    "AV_k14",
    "AV_k21",
    "AV_k23",
    "AV_x1",
    "AV_x2",
    "AV_x3",
    "AV_x4",
    "AV_i_NaCa_i_NaCa",
    "AV_alpha_h",
    "AV_beta_h",
    "AV_h_infinity",
    "AV_tau_h",
    "AV_E0_m",
    "AV_beta_m",
    "AV_m_infinity",
    "AV_alpha_m",
    "AV_tau_m",
    "AV_tau_y",
    "AV_y_infinity",
    "AV_q_infinity",
    "AV_tau_q",
    "AV_r_infinity",
    "AV_tau_r",
    "AV_E_Ca",
    "AV_E_Na",
    "AV_i_KACh_i_KACh",
    "AV_i_Kr_i_Kr",
    "AV_i_Kur_i_Kur",
    "AV_i_NaK_i_NaK",
    "AV_i_fK",
    "AV_i_fNa",
    "AV_i_to_i_to",
    "AV_i_f_i_f",
    "Istim",
    "Iion_cm",
};


static inline const Foam::HashTable<Foam::wordList>&
FabbriDependencyMap()
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
FabbriinitConsts(double* CONSTANTS, double* RATES, double* STATES, int tissueFlag, const Foam::dictionary& stimulus)
{
/* Ca_SR_release */
    CONSTANTS[AC_EC50_SR] = 0.45;
    CONSTANTS[AC_HSR] = 2.5;
    CONSTANTS[AC_MaxSR] = 15.0;
    CONSTANTS[AC_MinSR] = 1.0;
    CONSTANTS[AC_kiCa] = 500.0;
    CONSTANTS[AC_kim] = 5.0;
    CONSTANTS[AC_koCa] = 10000.0;
    CONSTANTS[AC_kom] = 660.0;
    CONSTANTS[AC_ks] = 1.48041085099999994e+08;
    
    /* Ca_buffering */
    CONSTANTS[AC_CM_tot] = 0.045;
    CONSTANTS[AC_CQ_tot] = 10.0;
    CONSTANTS[AC_Mgi] = 2.5;
    CONSTANTS[AC_TC_tot] = 0.031;
    CONSTANTS[AC_TMC_tot] = 0.062;
    CONSTANTS[AC_kb_CM] = 542.0;
    CONSTANTS[AC_kb_CQ] = 445.0;
    CONSTANTS[AC_kb_TC] = 446.0;
    CONSTANTS[AC_kb_TMC] = 7.51;
    CONSTANTS[AC_kb_TMM] = 751.0;
    CONSTANTS[AC_kf_CM] = 1642000.0;
    CONSTANTS[AC_kf_CQ] = 175.4;
    CONSTANTS[AC_kf_TC] = 88800.0;
    CONSTANTS[AC_kf_TMC] = 227700.0;
    CONSTANTS[AC_kf_TMM] = 2277.0;
    
    /* Cell_parameters */
    CONSTANTS[AC_L_cell] = 67.0;
    CONSTANTS[AC_L_sub] = 0.02;
    CONSTANTS[AC_R_cell] = 3.9;
    CONSTANTS[AC_V_i_part] = 0.46;
    CONSTANTS[AC_V_jsr_part] = 0.0012;
    CONSTANTS[AC_V_nsr_part] = 0.0116;
    CONSTANTS[AC_V_cell] = 1e-09 * 3.14159265358979312e+00 * pow(CONSTANTS[AC_R_cell],
                                                                 2.0) * CONSTANTS[AC_L_cell];
    CONSTANTS[AC_V_sub] = 1e-09 * 2.0 * 3.14159265358979312e+00 * CONSTANTS[AC_L_sub] * (CONSTANTS[AC_R_cell] - CONSTANTS[AC_L_sub] / 2.0) * CONSTANTS[AC_L_cell];
    CONSTANTS[AC_V_i] = CONSTANTS[AC_V_i_part] * CONSTANTS[AC_V_cell] - CONSTANTS[AC_V_sub];
    CONSTANTS[AC_V_jsr] = CONSTANTS[AC_V_jsr_part] * CONSTANTS[AC_V_cell];
    CONSTANTS[AC_V_nsr] = CONSTANTS[AC_V_nsr_part] * CONSTANTS[AC_V_cell];
    
    /* Rate_modulation_experiments */
    CONSTANTS[AC_ACh] = 0.0;
    CONSTANTS[AC_Iso_1_uM] = 0.0;
    
    /* i_CaL_fCa_gate */
    CONSTANTS[AC_Km_fCa] = 0.000338;
    CONSTANTS[AC_alpha_fCa] = 0.0075;
    
    /* Ca_intracellular_fluxes */
    CONSTANTS[AC_K_up] = 2.86113000000000003e-04;
    CONSTANTS[AC_P_up_basal] = 5.0;
    CONSTANTS[AC_b_up] = ((CONSTANTS[AC_Iso_1_uM] > 0.0) ? -0.25 : ((CONSTANTS[AC_ACh] > 0.0) ? 0.7 * CONSTANTS[AC_ACh] / (9e-05 + CONSTANTS[AC_ACh]) : 0.0));
    CONSTANTS[AC_slope_up] = 5e-05;
    CONSTANTS[AC_tau_dif_Ca] = 5.469e-05;
    CONSTANTS[AC_tau_tr] = 0.04;
    CONSTANTS[AC_P_up] = CONSTANTS[AC_P_up_basal] * (1.0 - CONSTANTS[AC_b_up]);
    
    /* Voltage_clamp */
    CONSTANTS[AC_V_holding] = -45.0;
    CONSTANTS[AC_V_test] = -35.0;
    CONSTANTS[AC_t_holding] = 0.5;
    CONSTANTS[AC_t_test] = 0.5;
    
    /* Ionic_values */
    CONSTANTS[AC_Cao] = 1.8;
    CONSTANTS[AC_Ki] = 140.0;
    CONSTANTS[AC_Ko] = 5.4;
    CONSTANTS[AC_Nao] = 140.0;
    
    /* Membrane */
    //CONSTANTS[AC_C] = 5.7e-05;
    CONSTANTS[AC_F] = 9.64853414999999950e+04;
    CONSTANTS[AC_Membrane_R] = 8314.472;
    CONSTANTS[AC_T] = 310.0;
    CONSTANTS[AC_clamp_mode] = 0.0;
    CONSTANTS[AC_RTONF] = CONSTANTS[AC_Membrane_R] * CONSTANTS[AC_T] / CONSTANTS[AC_F];
    
    /* Nai_concentration */
    CONSTANTS[AC_Nai_clamp] = 1.0;
    
    /* i_CaL */
    CONSTANTS[AC_ACh_block] = 0.31 * CONSTANTS[AC_ACh] / (CONSTANTS[AC_ACh] + 9e-05);
    CONSTANTS[AC_i_CaL_Iso_increase] = ((CONSTANTS[AC_Iso_1_uM] > 0.0) ? 1.23 : 1.0);
    CONSTANTS[AC_P_CaL] = 0.4578;
    
    /* i_CaL_dL_gate */
    CONSTANTS[AC_Iso_shift_dL] = ((CONSTANTS[AC_Iso_1_uM] > 0.0) ? -8.0 : 0.0);
    CONSTANTS[AC_Iso_slope_dL] = ((CONSTANTS[AC_Iso_1_uM] > 0.0) ? -27.0 : 0.0);
    CONSTANTS[AC_V_dL] = -16.4508;
    CONSTANTS[AC_k_dL] = 4.3371;
    
    /* i_CaL_fL_gate */
    CONSTANTS[AC_k_fL] = 0.0;
    CONSTANTS[AC_shift_fL] = 0.0;
    
    /* i_CaT */
    CONSTANTS[AC_P_CaT] = 0.04132;
    
    /* i_CaT_fT_gate */
    CONSTANTS[AC_offset_fT] = 0.0;
    
    /* i_KACh */
    CONSTANTS[AC_ACh_on] = 1.0;
    CONSTANTS[AC_g_KACh] = 0.00345;
    
    /* i_KACh_a_gate */
    CONSTANTS[AC_alpha_a] = (3.5988 - 0.025641) / (1.0 + 1.21550000000000005e-06 / pow(1.0 * CONSTANTS[AC_ACh], 1.6951)) + 0.025641;
    
    /* i_Kr */
    CONSTANTS[AC_g_Kr] = 0.00424;
    
    /* i_Ks */
    CONSTANTS[AC_g_Ks_] = 0.00065;
    CONSTANTS[AC_g_Ks] = ((CONSTANTS[AC_Iso_1_uM] > 0.0) ? 1.2 * CONSTANTS[AC_g_Ks_] : CONSTANTS[AC_g_Ks_]);
    
    /* i_Ks_n_gate */
    CONSTANTS[AC_i_Ks_n_gate_Iso_shift] = ((CONSTANTS[AC_Iso_1_uM] > 0.0) ? -14.0 : 0.0);
    
    /* i_Kur */
    CONSTANTS[AC_g_Kur] = 0.0001539;
    
    /* i_Na */
    CONSTANTS[AC_g_Na] = 0.0223;
    CONSTANTS[AC_g_Na_L] = 0.0;
    
    /* i_NaCa */
    CONSTANTS[AC_K1ni] = 395.3;
    CONSTANTS[AC_K1no] = 1628.0;
    CONSTANTS[AC_K2ni] = 2.289;
    CONSTANTS[AC_K2no] = 561.4;
    CONSTANTS[AC_K3ni] = 26.44;
    CONSTANTS[AC_K3no] = 4.663;
    CONSTANTS[AC_K_NaCa] = 3.343;
    CONSTANTS[AC_Kci] = 0.0207;
    CONSTANTS[AC_Kcni] = 26.44;
    CONSTANTS[AC_Kco] = 3.663;
    CONSTANTS[AC_Qci] = 0.1369;
    CONSTANTS[AC_Qco] = 0.0;
    CONSTANTS[AC_Qn] = 0.4315;
    CONSTANTS[AC_blockade_NaCa] = 0.0;
    CONSTANTS[AC_k34] = CONSTANTS[AC_Nao] / (CONSTANTS[AC_K3no] + CONSTANTS[AC_Nao]);
    
    /* i_NaK */
    CONSTANTS[AC_i_NaK_Iso_increase] = ((CONSTANTS[AC_Iso_1_uM] > 0.0) ? 1.2 : 1.0);
    CONSTANTS[AC_Km_Kp] = 1.4;
    CONSTANTS[AC_Km_Nap] = 14.0;
    CONSTANTS[AC_i_NaK_max] = 0.08105;
    
    /* i_Na_m_gate */
    CONSTANTS[AC_delta_m] = 1e-05;
    
    /* i_f */
    CONSTANTS[AC_Km_f] = 45.0;
    CONSTANTS[AC_alpha] = 0.5927;
    CONSTANTS[AC_blockade] = 0.0;
    CONSTANTS[AC_g_f] = 0.00427;
    CONSTANTS[AC_G_f] = CONSTANTS[AC_g_f] / (CONSTANTS[AC_Ko] / (CONSTANTS[AC_Ko] + CONSTANTS[AC_Km_f]));
    CONSTANTS[AC_G_f_K] = CONSTANTS[AC_G_f] / (CONSTANTS[AC_alpha] + 1.0);
    CONSTANTS[AC_G_f_Na] = CONSTANTS[AC_alpha] * CONSTANTS[AC_G_f_K];
    CONSTANTS[AC_g_f_K] = CONSTANTS[AC_G_f_K] * CONSTANTS[AC_Ko] / (CONSTANTS[AC_Ko] + CONSTANTS[AC_Km_f]);
    CONSTANTS[AC_g_f_Na] = CONSTANTS[AC_G_f_Na] * CONSTANTS[AC_Ko] / (CONSTANTS[AC_Ko] + CONSTANTS[AC_Km_f]);
    
    /* i_f_y_gate */
    CONSTANTS[AC_ACh_shift] = ((CONSTANTS[AC_ACh] > 0.0) ? -1.0 - 9.898 * pow(1.0 * CONSTANTS[AC_ACh], 0.618) / (pow(1.0 * CONSTANTS[AC_ACh], 0.618) + 1.22422999999999998e-03) : 0.0);
    CONSTANTS[AC_i_f_y_gate_Iso_shift] = ((CONSTANTS[AC_Iso_1_uM] > 0.0) ? 7.5 : 0.0);
    CONSTANTS[AC_y_shift] = 0.0;
    
    /* i_to */
    CONSTANTS[AC_g_to] = 0.0035;
    
    /* *remaining* */
    CONSTANTS[AC_E_K] = CONSTANTS[AC_RTONF] * log(CONSTANTS[AC_Ko] / CONSTANTS[AC_Ki]);

    /* --- Initial values --- */

    STATES[membrane_V] = -4.77871680000000012e+01;
    STATES[Na_i] = 5.0;
    STATES[If_y_gate_y] = 0.009508;
    STATES[INa_m_gate_m] = 0.447724;
    STATES[INa_h_gate_h] = 0.003058;
    STATES[ICaL_dL_gate_dL] = 0.001921;
    STATES[ICaL_fL_gate_fL] = 0.846702;
    STATES[ICaL_fCa_gate_fCa] = 0.844449;
    STATES[ICaT_dT_gate_dT] = 0.268909;
    STATES[ICaT_fT_gate_fT] = 0.020484;
    STATES[SR_R] = 0.9308;
    STATES[SR_O] =  6.18151200000000032e-09;
    STATES[SR_I] =  4.59562199999999991e-10;
    STATES[SR_RI] = 0.069199;
    STATES[Buffer_fTMM] = 0.653777;
    STATES[Buffer_fCMi] = 0.217311;
    STATES[Buffer_fCMs] = 0.158521;
    STATES[Buffer_fTC] = 0.017929;
    STATES[Buffer_fTMC] = 0.259947;
    STATES[Buffer_fCQ] = 0.138975;
    STATES[Ca_i] =  9.15640999999999997e-06;
    STATES[Ca_nsr] = 0.435148;
    STATES[Ca_jsr] = 0.409551;
    STATES[Ca_sub] =  6.22610399999999950e-05;
    STATES[IKur_rKur_gate_r_Kur] = 0.011845;
    STATES[IKur_sKur_gate_s_Kur] = 0.845304;
    STATES[Ito_q_gate_q] = 0.430836;
    STATES[Ito_r_gate_r] = 0.014523;
    STATES[IKr_pa_gate_paS] = 0.283185;
    STATES[IKr_pa_gate_paF] = 0.011068;
    STATES[IKr_pi_gate_piy] = 0.709051;
    STATES[IKs_n_gate_n] = 0.1162;
    STATES[IKACh_a_gate_a] = 0.00277;
}


void
FabbricomputeVariables(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC, int tissueFlag, bool solveVmWithinODESolver)
{
/* Ca_SR_release */
    ALGEBRAIC[AV_P_tot] = STATES[SR_R] + STATES[SR_O] + STATES[SR_I] + STATES[SR_RI];
    ALGEBRAIC[AV_diff] = STATES[Ca_jsr] - STATES[Ca_sub];
    ALGEBRAIC[AV_j_SRCarel] = CONSTANTS[AC_ks] * STATES[SR_O] * (STATES[Ca_jsr] - STATES[Ca_sub]);
    ALGEBRAIC[AV_kCaSR] = CONSTANTS[AC_MaxSR] - (CONSTANTS[AC_MaxSR] - CONSTANTS[AC_MinSR]) / (1.0 + pow(CONSTANTS[AC_EC50_SR] / STATES[Ca_jsr], CONSTANTS[AC_HSR]));
    ALGEBRAIC[AV_kiSRCa] = CONSTANTS[AC_kiCa] * ALGEBRAIC[AV_kCaSR];
    ALGEBRAIC[AV_koSRCa] = CONSTANTS[AC_koCa] / ALGEBRAIC[AV_kCaSR];
    RATES[SR_I] = ALGEBRAIC[AV_kiSRCa] * STATES[Ca_sub] * STATES[SR_O] - CONSTANTS[AC_kim] * STATES[SR_I] - (CONSTANTS[AC_kom] * STATES[SR_I] - ALGEBRAIC[AV_koSRCa] * pow(STATES[Ca_sub], 2.0) * STATES[SR_RI]);
    RATES[SR_O] = ALGEBRAIC[AV_koSRCa] * pow(STATES[Ca_sub], 2.0) * STATES[SR_R] - CONSTANTS[AC_kom] * STATES[SR_O] - (ALGEBRAIC[AV_kiSRCa] * STATES[Ca_sub] * STATES[SR_O] - CONSTANTS[AC_kim] * STATES[SR_I]);
    RATES[SR_R] = CONSTANTS[AC_kim] * STATES[SR_RI] - ALGEBRAIC[AV_kiSRCa] * STATES[Ca_sub] * STATES[SR_R] - (ALGEBRAIC[AV_koSRCa] * pow(STATES[Ca_sub], 2.0) * STATES[SR_R] - CONSTANTS[AC_kom] * STATES[SR_O]);
    RATES[SR_RI] = CONSTANTS[AC_kom] * STATES[SR_I] - ALGEBRAIC[AV_koSRCa] * pow(STATES[Ca_sub],
                                                      2.0) * STATES[SR_RI] - (CONSTANTS[AC_kim] * STATES[SR_RI] - ALGEBRAIC[AV_kiSRCa] * STATES[Ca_sub] * STATES[SR_R]);
    
    /* Ca_buffering */
    ALGEBRAIC[AV_delta_fCMi] = CONSTANTS[AC_kf_CM] * STATES[Ca_i] * (1.0 - STATES[Buffer_fCMi]) - CONSTANTS[AC_kb_CM] * STATES[Buffer_fCMi];
    ALGEBRAIC[AV_delta_fCMs] = CONSTANTS[AC_kf_CM] * STATES[Ca_sub] * (1.0 - STATES[Buffer_fCMs]) - CONSTANTS[AC_kb_CM] * STATES[Buffer_fCMs];
    ALGEBRAIC[AV_delta_fCQ] = CONSTANTS[AC_kf_CQ] * STATES[Ca_jsr] * (1.0 - STATES[Buffer_fCQ]) - CONSTANTS[AC_kb_CQ] * STATES[Buffer_fCQ];
    ALGEBRAIC[AV_delta_fTC] = CONSTANTS[AC_kf_TC] * STATES[Ca_i] * (1.0 - STATES[Buffer_fTC]) - CONSTANTS[AC_kb_TC] * STATES[Buffer_fTC];
    ALGEBRAIC[AV_delta_fTMC] = CONSTANTS[AC_kf_TMC] * STATES[Ca_i] * (1.0 - (STATES[Buffer_fTMC] + STATES[Buffer_fTMM])) - CONSTANTS[AC_kb_TMC] * STATES[Buffer_fTMC];
    ALGEBRAIC[AV_delta_fTMM] = CONSTANTS[AC_kf_TMM] * CONSTANTS[AC_Mgi] * (1.0 - (STATES[Buffer_fTMC] + STATES[Buffer_fTMM])) - CONSTANTS[AC_kb_TMM] * STATES[Buffer_fTMM];
    RATES[Buffer_fCMi] = ALGEBRAIC[AV_delta_fCMi];
    RATES[Buffer_fCMs] = ALGEBRAIC[AV_delta_fCMs];
    RATES[Buffer_fCQ] = ALGEBRAIC[AV_delta_fCQ];
    RATES[Buffer_fTC] = ALGEBRAIC[AV_delta_fTC];
    RATES[Buffer_fTMC] = ALGEBRAIC[AV_delta_fTMC];
    RATES[Buffer_fTMM] = ALGEBRAIC[AV_delta_fTMM];
    

    
    /* i_CaL_fCa_gate */
    ALGEBRAIC[AV_fCa_infinity] = CONSTANTS[AC_Km_fCa] / (CONSTANTS[AC_Km_fCa] + STATES[Ca_sub]);
    ALGEBRAIC[AV_tau_fCa] = 0.001 * ALGEBRAIC[AV_fCa_infinity] / CONSTANTS[AC_alpha_fCa];
    RATES[ICaL_fCa_gate_fCa] = (ALGEBRAIC[AV_fCa_infinity] - STATES[ICaL_fCa_gate_fCa]) / ALGEBRAIC[AV_tau_fCa];
    
    /* Ca_intracellular_fluxes */
    ALGEBRAIC[AV_j_Ca_dif] = (STATES[Ca_sub] - STATES[Ca_i]) / CONSTANTS[AC_tau_dif_Ca];
    ALGEBRAIC[AV_j_tr] = (STATES[Ca_nsr] - STATES[Ca_jsr]) / CONSTANTS[AC_tau_tr];
    ALGEBRAIC[AV_j_up] = CONSTANTS[AC_P_up] / (1.0 + exp((-STATES[Ca_i] + CONSTANTS[AC_K_up]) / CONSTANTS[AC_slope_up]));
    
    /* Voltage_clamp */
    ALGEBRAIC[AV_V_clamp] = (((VOI > CONSTANTS[AC_t_holding]) && (VOI < CONSTANTS[AC_t_holding] + CONSTANTS[AC_t_test])) ? CONSTANTS[AC_V_test] : CONSTANTS[AC_V_holding]);
    
    /* Ca_dynamics */
    RATES[Ca_jsr] = ALGEBRAIC[AV_j_tr] - (ALGEBRAIC[AV_j_SRCarel] + CONSTANTS[AC_CQ_tot] * ALGEBRAIC[AV_delta_fCQ]);
    RATES[Ca_nsr] = ALGEBRAIC[AV_j_up] - ALGEBRAIC[AV_j_tr] * CONSTANTS[AC_V_jsr] / CONSTANTS[AC_V_nsr];
    RATES[Ca_i] = 1.0 * (ALGEBRAIC[AV_j_Ca_dif] * CONSTANTS[AC_V_sub] - ALGEBRAIC[AV_j_up] * CONSTANTS[AC_V_nsr]) / CONSTANTS[AC_V_i] - (CONSTANTS[AC_CM_tot] * ALGEBRAIC[AV_delta_fCMi] + CONSTANTS[AC_TC_tot] * ALGEBRAIC[AV_delta_fTC] + CONSTANTS[AC_TMC_tot] * ALGEBRAIC[AV_delta_fTMC]);
    
    /* Membrane */
    ALGEBRAIC[AV_V] = ((CONSTANTS[AC_clamp_mode] >= 1.0) ? ALGEBRAIC[AV_V_clamp] : STATES[membrane_V]);
    
    /* Nai_concentration */
    ALGEBRAIC[AV_Nai] = STATES[Na_i];
    
    /* i_CaL */
    ALGEBRAIC[AV_i_siCa] = 2.0 * CONSTANTS[AC_P_CaL] * (ALGEBRAIC[AV_V] - 0.0) / (CONSTANTS[AC_RTONF] * (1.0 - exp(-1.0 * (ALGEBRAIC[AV_V] - 0.0) * 2.0 / CONSTANTS[AC_RTONF]))) * (STATES[Ca_sub] - CONSTANTS[AC_Cao] * exp(-2.0 * (ALGEBRAIC[AV_V] - 0.0) / CONSTANTS[AC_RTONF])) * STATES[ICaL_dL_gate_dL] * STATES[ICaL_fL_gate_fL] * STATES[ICaL_fCa_gate_fCa];
    ALGEBRAIC[AV_i_siK] = 0.000365 * CONSTANTS[AC_P_CaL] * (ALGEBRAIC[AV_V] - 0.0) / (CONSTANTS[AC_RTONF] * (1.0 - exp(-1.0 * (ALGEBRAIC[AV_V] - 0.0) / CONSTANTS[AC_RTONF]))) * (CONSTANTS[AC_Ki] - CONSTANTS[AC_Ko] * exp(-1.0 * (ALGEBRAIC[AV_V] - 0.0) / CONSTANTS[AC_RTONF])) * STATES[ICaL_dL_gate_dL] * STATES[ICaL_fL_gate_fL] * STATES[ICaL_fCa_gate_fCa];
    ALGEBRAIC[AV_i_siNa] = 1.85e-05 * CONSTANTS[AC_P_CaL] * (ALGEBRAIC[AV_V] - 0.0) / (CONSTANTS[AC_RTONF] * (1.0 - exp(-1.0 * (ALGEBRAIC[AV_V] - 0.0) / CONSTANTS[AC_RTONF]))) * (ALGEBRAIC[AV_Nai] - CONSTANTS[AC_Nao] * exp(-1.0 * (ALGEBRAIC[AV_V] - 0.0) / CONSTANTS[AC_RTONF])) * STATES[ICaL_dL_gate_dL] * STATES[ICaL_fL_gate_fL] * STATES[ICaL_fCa_gate_fCa];
    ALGEBRAIC[AV_i_CaL_i_CaL] = (ALGEBRAIC[AV_i_siCa] + ALGEBRAIC[AV_i_siK] + ALGEBRAIC[AV_i_siNa]) * (1.0 - CONSTANTS[AC_ACh_block]) * 1.0 * CONSTANTS[AC_i_CaL_Iso_increase];
    
    /* i_CaL_dL_gate */
    ALGEBRAIC[AV_adVm] = ((ALGEBRAIC[AV_V] == -41.8) ? -41.80001 : ((ALGEBRAIC[AV_V] == 0.0) ? 0.0 : ((ALGEBRAIC[AV_V] == -6.8) ? -6.80001 : ALGEBRAIC[AV_V])));
    ALGEBRAIC[AV_bdVm] = ((ALGEBRAIC[AV_V] == -1.8) ? -1.80001 : ALGEBRAIC[AV_V]);
    ALGEBRAIC[AV_alpha_dL] = -0.02839 * (ALGEBRAIC[AV_adVm] + 41.8) / (exp(-(ALGEBRAIC[AV_adVm] + 41.8) / 2.5) - 1.0) - 0.0849 * (ALGEBRAIC[AV_adVm] + 6.8) / (exp(-(ALGEBRAIC[AV_adVm] + 6.8) / 4.8) - 1.0);
    ALGEBRAIC[AV_beta_dL] = 0.01143 * (ALGEBRAIC[AV_bdVm] + 1.8) / (exp((ALGEBRAIC[AV_bdVm] + 1.8) / 2.5) - 1.0);
    ALGEBRAIC[AV_dL_infinity] = 1.0 / (1.0 + exp(-(ALGEBRAIC[AV_V] - CONSTANTS[AC_V_dL] - CONSTANTS[AC_Iso_shift_dL]) / (CONSTANTS[AC_k_dL] * (1.0 + CONSTANTS[AC_Iso_slope_dL] / 100.0))));
    ALGEBRAIC[AV_tau_dL] = 0.001 / (ALGEBRAIC[AV_alpha_dL] + ALGEBRAIC[AV_beta_dL]);
    RATES[ICaL_dL_gate_dL] = (ALGEBRAIC[AV_dL_infinity] - STATES[ICaL_dL_gate_dL]) / ALGEBRAIC[AV_tau_dL];
    
    /* i_CaL_fL_gate */
    ALGEBRAIC[AV_tau_fL] = 0.001 * (44.3 + 230.0 * exp(-pow((ALGEBRAIC[AV_V] + 36.0) / 10.0, 2.0)));
    ALGEBRAIC[AV_fL_infinity] = 1.0 / (1.0 + exp((ALGEBRAIC[AV_V] + 37.4 + CONSTANTS[AC_shift_fL]) / (5.3 + CONSTANTS[AC_k_fL])));
    RATES[ICaL_fL_gate_fL] = (ALGEBRAIC[AV_fL_infinity] - STATES[ICaL_fL_gate_fL]) / ALGEBRAIC[AV_tau_fL];
    
    /* i_CaT */
    ALGEBRAIC[AV_i_CaT_i_CaT] = 2.0 * CONSTANTS[AC_P_CaT] * ALGEBRAIC[AV_V] / (CONSTANTS[AC_RTONF] * (1.0 - exp(-1.0 * ALGEBRAIC[AV_V] * 2.0 / CONSTANTS[AC_RTONF]))) * (STATES[Ca_sub] - CONSTANTS[AC_Cao] * exp(-2.0 * ALGEBRAIC[AV_V] / CONSTANTS[AC_RTONF])) * STATES[ICaT_dT_gate_dT] * STATES[ICaT_fT_gate_fT];
    
    /* i_CaT_dT_gate */
    ALGEBRAIC[AV_dT_infinity] = 1.0 / (1.0 + exp(-(ALGEBRAIC[AV_V] + 38.3) / 5.5));
    ALGEBRAIC[AV_tau_dT] = 0.001 / (1.068 * exp((ALGEBRAIC[AV_V] + 38.3) / 30.0) + 1.068 * exp(-(ALGEBRAIC[AV_V] + 38.3) / 30.0));
    RATES[ICaT_dT_gate_dT] = (ALGEBRAIC[AV_dT_infinity] - STATES[ICaT_dT_gate_dT]) / ALGEBRAIC[AV_tau_dT];
    
    /* i_CaT_fT_gate */
    ALGEBRAIC[AV_fT_infinity] = 1.0 / (1.0 + exp((ALGEBRAIC[AV_V] + 58.7) / 3.8));
    ALGEBRAIC[AV_tau_fT] = 1.0 / (16.67 * exp(-(ALGEBRAIC[AV_V] + 75.0) / 83.3) + 16.67 * exp((ALGEBRAIC[AV_V] + 75.0) / 15.38)) + CONSTANTS[AC_offset_fT];
    RATES[ICaT_fT_gate_fT] = (ALGEBRAIC[AV_fT_infinity] - STATES[ICaT_fT_gate_fT]) / ALGEBRAIC[AV_tau_fT];
    
    /* i_KACh_a_gate */
    ALGEBRAIC[AV_beta_a] = 10.0 * exp(0.0133 * (ALGEBRAIC[AV_V] + 40.0));
    ALGEBRAIC[AV_a_infinity] = CONSTANTS[AC_alpha_a] / (CONSTANTS[AC_alpha_a] + ALGEBRAIC[AV_beta_a]);
    ALGEBRAIC[AV_tau_a] = 1.0 / (CONSTANTS[AC_alpha_a] + ALGEBRAIC[AV_beta_a]);
    RATES[IKACh_a_gate_a] = (ALGEBRAIC[AV_a_infinity] - STATES[IKACh_a_gate_a]) / ALGEBRAIC[AV_tau_a];
    
    /* i_Kr_pa_gate */
    ALGEBRAIC[AV_alfapaF] = 1.0 / (1.0 + exp(-(ALGEBRAIC[AV_V] + 23.2) / 6.6)) / (8.46553540000000049e-01 / (37.2 * exp(ALGEBRAIC[AV_V] / 11.9) + 0.96 * exp(-ALGEBRAIC[AV_V] / 18.5)));
    ALGEBRAIC[AV_betapaF] = 4.0 * ((37.2 * exp(ALGEBRAIC[AV_V] / 15.9) + 0.96 * exp(-ALGEBRAIC[AV_V] / 22.5)) / 8.46553540000000049e-01 - 1.0 / (1.0 + exp(-(ALGEBRAIC[AV_V] + 23.2) / 10.6)) / (8.46553540000000049e-01 / (37.2 * exp(ALGEBRAIC[AV_V] / 15.9) + 0.96 * exp(-ALGEBRAIC[AV_V] / 22.5))));
    ALGEBRAIC[AV_pa_infinity] = 1.0 / (1.0 + exp(-(ALGEBRAIC[AV_V] + 10.0144) / 7.6607));
    ALGEBRAIC[AV_tau_paF] = 1.0 / (30.0 * exp(ALGEBRAIC[AV_V] / 10.0) + exp(-ALGEBRAIC[AV_V] / 12.0));
    ALGEBRAIC[AV_tau_paS] = 8.46553540000000049e-01 / (4.2 * exp(ALGEBRAIC[AV_V] / 17.0) + 0.15 * exp(-ALGEBRAIC[AV_V] / 21.6));
    RATES[IKr_pa_gate_paF] = (ALGEBRAIC[AV_pa_infinity] - STATES[IKr_pa_gate_paF]) / ALGEBRAIC[AV_tau_paF];
    RATES[IKr_pa_gate_paS] = (ALGEBRAIC[AV_pa_infinity] - STATES[IKr_pa_gate_paS]) / ALGEBRAIC[AV_tau_paS];
    
    /* i_Kr_pi_gate */
    ALGEBRAIC[AV_pi_infinity] = 1.0 / (1.0 + exp((ALGEBRAIC[AV_V] + 28.6) / 17.1));
    ALGEBRAIC[AV_tau_pi] = 1.0 / (100.0 * exp(-ALGEBRAIC[AV_V] / 54.645) + 656.0 * exp(ALGEBRAIC[AV_V] / 106.157));
    RATES[IKr_pi_gate_piy] = (ALGEBRAIC[AV_pi_infinity] - STATES[IKr_pi_gate_piy]) / ALGEBRAIC[AV_tau_pi];
    
    /* i_Ks */
    ALGEBRAIC[AV_E_Ks] = CONSTANTS[AC_RTONF] * log((CONSTANTS[AC_Ko] + 0.12 * CONSTANTS[AC_Nao]) / (CONSTANTS[AC_Ki] + 0.12 * ALGEBRAIC[AV_Nai]));
    ALGEBRAIC[AV_i_Ks_i_Ks] = CONSTANTS[AC_g_Ks] * (ALGEBRAIC[AV_V] - ALGEBRAIC[AV_E_Ks]) * pow(STATES[IKs_n_gate_n], 2.0);
    
    /* i_Ks_n_gate */
    ALGEBRAIC[AV_alpha_n] = 28.0 / (1.0 + exp(-(ALGEBRAIC[AV_V] - 40.0 - CONSTANTS[AC_i_Ks_n_gate_Iso_shift]) / 3.0));
    ALGEBRAIC[AV_beta_n] = 1.0 * exp(-(ALGEBRAIC[AV_V] - CONSTANTS[AC_i_Ks_n_gate_Iso_shift] - 5.0) / 25.0);
    ALGEBRAIC[AV_n_infinity] = sqrt(1.0 / (1.0 + exp(-(ALGEBRAIC[AV_V] + 0.6383 - CONSTANTS[AC_i_Ks_n_gate_Iso_shift]) / 10.7071)));
    ALGEBRAIC[AV_tau_n] = 1.0 / (ALGEBRAIC[AV_alpha_n] + ALGEBRAIC[AV_beta_n]);
    RATES[IKs_n_gate_n] = (ALGEBRAIC[AV_n_infinity] - STATES[IKs_n_gate_n]) / ALGEBRAIC[AV_tau_n];
    
    /* i_Kur_rKur_gate */
    ALGEBRAIC[AV_r_Kur_infinity] = 1.0 / (1.0 + exp((ALGEBRAIC[AV_V] + 6.0) / -8.6));
    ALGEBRAIC[AV_tau_r_Kur] = 0.009 / (1.0 + exp((ALGEBRAIC[AV_V] + 5.0) / 12.0)) + 0.0005;
    RATES[IKur_rKur_gate_r_Kur] = (ALGEBRAIC[AV_r_Kur_infinity] - STATES[IKur_rKur_gate_r_Kur]) / ALGEBRAIC[AV_tau_r_Kur];
    
    /* i_Kur_sKur_gate */
    ALGEBRAIC[AV_s_Kur_infinity] = 1.0 / (1.0 + exp((ALGEBRAIC[AV_V] + 7.5) / 10.0));
    ALGEBRAIC[AV_tau_s_Kur] = 0.59 / (1.0 + exp((ALGEBRAIC[AV_V] + 60.0) / 10.0)) + 3.05;
    RATES[IKur_sKur_gate_s_Kur] = (ALGEBRAIC[AV_s_Kur_infinity] - STATES[IKur_sKur_gate_s_Kur]) / ALGEBRAIC[AV_tau_s_Kur];
    
    /* i_Na */
    ALGEBRAIC[AV_E_mh] = CONSTANTS[AC_RTONF] * log((CONSTANTS[AC_Nao] + 0.12 * CONSTANTS[AC_Ko]) / (ALGEBRAIC[AV_Nai] + 0.12 * CONSTANTS[AC_Ki]));
    ALGEBRAIC[AV_i_Na_] = CONSTANTS[AC_g_Na] * pow(STATES[INa_m_gate_m], 3.0) * STATES[INa_h_gate_h] * (ALGEBRAIC[AV_V] - ALGEBRAIC[AV_E_mh]);
    ALGEBRAIC[AV_i_Na_L] = CONSTANTS[AC_g_Na_L] * pow(STATES[INa_m_gate_m], 3.0) * (ALGEBRAIC[AV_V] - ALGEBRAIC[AV_E_mh]);
    ALGEBRAIC[AV_i_Na_i_Na] = ALGEBRAIC[AV_i_Na_] + ALGEBRAIC[AV_i_Na_L];
    
    /* i_NaCa */
    ALGEBRAIC[AV_di] = 1.0 + STATES[Ca_sub] / CONSTANTS[AC_Kci] * (1.0 + exp(-CONSTANTS[AC_Qci] * ALGEBRAIC[AV_V] / CONSTANTS[AC_RTONF]) + ALGEBRAIC[AV_Nai] / CONSTANTS[AC_Kcni]) + ALGEBRAIC[AV_Nai] / CONSTANTS[AC_K1ni] * (1.0 + ALGEBRAIC[AV_Nai] / CONSTANTS[AC_K2ni] * (1.0 + ALGEBRAIC[AV_Nai] / CONSTANTS[AC_K3ni]));
    ALGEBRAIC[AV_i_NaCa_do] = 1.0 + CONSTANTS[AC_Cao] / CONSTANTS[AC_Kco] * (1.0 + exp(CONSTANTS[AC_Qco] * ALGEBRAIC[AV_V] / CONSTANTS[AC_RTONF])) + CONSTANTS[AC_Nao] / CONSTANTS[AC_K1no] * (1.0 + CONSTANTS[AC_Nao] / CONSTANTS[AC_K2no] * (1.0 + CONSTANTS[AC_Nao] / CONSTANTS[AC_K3no]));
    ALGEBRAIC[AV_k32] = exp(CONSTANTS[AC_Qn] * ALGEBRAIC[AV_V] / (2.0 * CONSTANTS[AC_RTONF]));
    ALGEBRAIC[AV_k41] = exp(-CONSTANTS[AC_Qn] * ALGEBRAIC[AV_V] / (2.0 * CONSTANTS[AC_RTONF]));
    ALGEBRAIC[AV_k43] = ALGEBRAIC[AV_Nai] / (CONSTANTS[AC_K3ni] + ALGEBRAIC[AV_Nai]);
    ALGEBRAIC[AV_k12] = STATES[Ca_sub] / CONSTANTS[AC_Kci] * exp(-CONSTANTS[AC_Qci] * ALGEBRAIC[AV_V] / CONSTANTS[AC_RTONF]) / ALGEBRAIC[AV_di];
    ALGEBRAIC[AV_k14] = ALGEBRAIC[AV_Nai] / CONSTANTS[AC_K1ni] * ALGEBRAIC[AV_Nai] / CONSTANTS[AC_K2ni] * (1.0 + ALGEBRAIC[AV_Nai] / CONSTANTS[AC_K3ni]) * exp(CONSTANTS[AC_Qn] * ALGEBRAIC[AV_V] / (2.0 * CONSTANTS[AC_RTONF])) / ALGEBRAIC[AV_di];
    ALGEBRAIC[AV_k21] = CONSTANTS[AC_Cao] / CONSTANTS[AC_Kco] * exp(CONSTANTS[AC_Qco] * ALGEBRAIC[AV_V] / CONSTANTS[AC_RTONF]) / ALGEBRAIC[AV_i_NaCa_do];
    ALGEBRAIC[AV_k23] = CONSTANTS[AC_Nao] / CONSTANTS[AC_K1no] * CONSTANTS[AC_Nao] / CONSTANTS[AC_K2no] * (1.0 + CONSTANTS[AC_Nao] / CONSTANTS[AC_K3no]) * exp(-CONSTANTS[AC_Qn] * ALGEBRAIC[AV_V] / (2.0 * CONSTANTS[AC_RTONF])) / ALGEBRAIC[AV_i_NaCa_do];
    ALGEBRAIC[AV_x1] = ALGEBRAIC[AV_k41] * CONSTANTS[AC_k34] * (ALGEBRAIC[AV_k23] + ALGEBRAIC[AV_k21]) + ALGEBRAIC[AV_k21] * ALGEBRAIC[AV_k32] * (ALGEBRAIC[AV_k43] + ALGEBRAIC[AV_k41]);
    ALGEBRAIC[AV_x2] = ALGEBRAIC[AV_k32] * ALGEBRAIC[AV_k43] * (ALGEBRAIC[AV_k14] + ALGEBRAIC[AV_k12]) + ALGEBRAIC[AV_k41] * ALGEBRAIC[AV_k12] * (CONSTANTS[AC_k34] + ALGEBRAIC[AV_k32]);
    ALGEBRAIC[AV_x3] = ALGEBRAIC[AV_k14] * ALGEBRAIC[AV_k43] * (ALGEBRAIC[AV_k23] + ALGEBRAIC[AV_k21]) + ALGEBRAIC[AV_k12] * ALGEBRAIC[AV_k23] * (ALGEBRAIC[AV_k43] + ALGEBRAIC[AV_k41]);
    ALGEBRAIC[AV_x4] = ALGEBRAIC[AV_k23] * CONSTANTS[AC_k34] * (ALGEBRAIC[AV_k14] + ALGEBRAIC[AV_k12]) + ALGEBRAIC[AV_k14] * ALGEBRAIC[AV_k21] * (CONSTANTS[AC_k34] + ALGEBRAIC[AV_k32]);
    ALGEBRAIC[AV_i_NaCa_i_NaCa] = (1.0 - CONSTANTS[AC_blockade_NaCa]) * CONSTANTS[AC_K_NaCa] * (ALGEBRAIC[AV_x2] * ALGEBRAIC[AV_k21] - ALGEBRAIC[AV_x1] * ALGEBRAIC[AV_k12]) / (ALGEBRAIC[AV_x1] + ALGEBRAIC[AV_x2] + ALGEBRAIC[AV_x3] + ALGEBRAIC[AV_x4]);
    
    /* i_Na_h_gate */
    ALGEBRAIC[AV_alpha_h] = 20.0 * exp(-0.125 * (ALGEBRAIC[AV_V] + 75.0));
    ALGEBRAIC[AV_beta_h] = 2000.0 / (320.0 * exp(-0.1 * (ALGEBRAIC[AV_V] + 75.0)) + 1.0);
    ALGEBRAIC[AV_h_infinity] = 1.0 / (1.0 + exp((ALGEBRAIC[AV_V] + 69.804) / 4.4565));
    ALGEBRAIC[AV_tau_h] = 1.0 / (ALGEBRAIC[AV_alpha_h] + ALGEBRAIC[AV_beta_h]);
    RATES[INa_h_gate_h] = (ALGEBRAIC[AV_h_infinity] - STATES[INa_h_gate_h]) / ALGEBRAIC[AV_tau_h];
    
    /* i_Na_m_gate */
    ALGEBRAIC[AV_E0_m] = ALGEBRAIC[AV_V] + 41.0;
    ALGEBRAIC[AV_beta_m] = 8000.0 * exp(-0.056 * (ALGEBRAIC[AV_V] + 66.0));
    ALGEBRAIC[AV_m_infinity] = 1.0 / (1.0 + exp(-(ALGEBRAIC[AV_V] + 42.0504) / 8.3106));
    ALGEBRAIC[AV_alpha_m] = ((fabs(ALGEBRAIC[AV_E0_m]) < CONSTANTS[AC_delta_m]) ? 2000.0 : 200.0 * ALGEBRAIC[AV_E0_m] / (1.0 - exp(-0.1 * ALGEBRAIC[AV_E0_m])));
    ALGEBRAIC[AV_tau_m] = 1.0 / (ALGEBRAIC[AV_alpha_m] + ALGEBRAIC[AV_beta_m]);
    RATES[INa_m_gate_m] = (ALGEBRAIC[AV_m_infinity] - STATES[INa_m_gate_m]) / ALGEBRAIC[AV_tau_m];
    
    /* i_f_y_gate */
    ALGEBRAIC[AV_tau_y] = 1.0 / (0.36 * (ALGEBRAIC[AV_V] + 148.8 - CONSTANTS[AC_ACh_shift] - CONSTANTS[AC_i_f_y_gate_Iso_shift]) / (exp(0.066 * (ALGEBRAIC[AV_V] + 148.8 - CONSTANTS[AC_ACh_shift] - CONSTANTS[AC_i_f_y_gate_Iso_shift])) - 1.0) + 0.1 * (ALGEBRAIC[AV_V] + 87.3 - CONSTANTS[AC_ACh_shift] - CONSTANTS[AC_i_f_y_gate_Iso_shift]) / (1.0 - exp(-0.2 * (ALGEBRAIC[AV_V] + 87.3 - CONSTANTS[AC_ACh_shift] - CONSTANTS[AC_i_f_y_gate_Iso_shift])))) - 0.054;
    ALGEBRAIC[AV_y_infinity] = ((ALGEBRAIC[AV_V] < -(80.0 - CONSTANTS[AC_ACh_shift] - CONSTANTS[AC_i_f_y_gate_Iso_shift] - CONSTANTS[AC_y_shift])) ? 0.01329 + 0.99921 / (1.0 + exp((ALGEBRAIC[AV_V] + 97.134 - CONSTANTS[AC_ACh_shift] - CONSTANTS[AC_i_f_y_gate_Iso_shift] - CONSTANTS[AC_y_shift]) / 8.1752)) : 0.0002501 * exp(-(ALGEBRAIC[AV_V] - CONSTANTS[AC_ACh_shift] - CONSTANTS[AC_i_f_y_gate_Iso_shift] - CONSTANTS[AC_y_shift]) / 12.861));
    RATES[If_y_gate_y] = (ALGEBRAIC[AV_y_infinity] - STATES[If_y_gate_y]) / ALGEBRAIC[AV_tau_y];
    
    /* i_to_q_gate */
    ALGEBRAIC[AV_q_infinity] = 1.0 / (1.0 + exp((ALGEBRAIC[AV_V] + 49.0) / 13.0));
    ALGEBRAIC[AV_tau_q] = 0.001 * 0.6 * (65.17 / (0.57 * exp(-0.08 * (ALGEBRAIC[AV_V] + 44.0)) + 0.065 * exp(0.1 * (ALGEBRAIC[AV_V] + 45.93))) + 10.1);
    RATES[Ito_q_gate_q] = (ALGEBRAIC[AV_q_infinity] - STATES[Ito_q_gate_q]) / ALGEBRAIC[AV_tau_q];
    
    /* i_to_r_gate */
    ALGEBRAIC[AV_r_infinity] = 1.0 / (1.0 + exp(-(ALGEBRAIC[AV_V] - 19.3) / 15.0));
    ALGEBRAIC[AV_tau_r] = 0.001 * 0.66 * 1.4 * (15.59 / (1.037 * exp(0.09 * (ALGEBRAIC[AV_V] + 30.61)) + 0.369 * exp(-0.12 * (ALGEBRAIC[AV_V] + 23.84))) + 2.98);
    RATES[Ito_r_gate_r] = (ALGEBRAIC[AV_r_infinity] - STATES[Ito_r_gate_r]) / ALGEBRAIC[AV_tau_r];
    
    /* *remaining* */
    RATES[Ca_sub] = ALGEBRAIC[AV_j_SRCarel] * CONSTANTS[AC_V_jsr] / CONSTANTS[AC_V_sub] - ((ALGEBRAIC[AV_i_siCa] + ALGEBRAIC[AV_i_CaT_i_CaT] - 2.0 * ALGEBRAIC[AV_i_NaCa_i_NaCa]) / (2.0 * CONSTANTS[AC_F] * CONSTANTS[AC_V_sub]) + ALGEBRAIC[AV_j_Ca_dif] + CONSTANTS[AC_CM_tot] * ALGEBRAIC[AV_delta_fCMs]);
    ALGEBRAIC[AV_E_Ca] = 0.5 * CONSTANTS[AC_RTONF] * log(CONSTANTS[AC_Cao] / STATES[Ca_sub]);
    ALGEBRAIC[AV_E_Na] = CONSTANTS[AC_RTONF] * log(CONSTANTS[AC_Nao] / ALGEBRAIC[AV_Nai]);
    ALGEBRAIC[AV_i_KACh_i_KACh] = ((CONSTANTS[AC_ACh] > 0.0) ? CONSTANTS[AC_ACh_on] * CONSTANTS[AC_g_KACh] * (ALGEBRAIC[AV_V] - CONSTANTS[AC_E_K]) * (1.0 + exp((ALGEBRAIC[AV_V] + 20.0) / 20.0)) * STATES[IKACh_a_gate_a] : 0.0);
    ALGEBRAIC[AV_i_Kr_i_Kr] = CONSTANTS[AC_g_Kr] * (ALGEBRAIC[AV_V] - CONSTANTS[AC_E_K]) * (0.9 * STATES[IKr_pa_gate_paF] + 0.1 * STATES[IKr_pa_gate_paS]) * STATES[IKr_pi_gate_piy];
    ALGEBRAIC[AV_i_Kur_i_Kur] = CONSTANTS[AC_g_Kur] * STATES[IKur_rKur_gate_r_Kur] * STATES[IKur_sKur_gate_s_Kur] * (ALGEBRAIC[AV_V] - CONSTANTS[AC_E_K]);
    ALGEBRAIC[AV_i_NaK_i_NaK] = CONSTANTS[AC_i_NaK_Iso_increase] * CONSTANTS[AC_i_NaK_max] * pow(1.0 + pow(CONSTANTS[AC_Km_Kp] / CONSTANTS[AC_Ko], 1.2),
                                                                                      -1.0) * pow(1.0 + pow(CONSTANTS[AC_Km_Nap] / ALGEBRAIC[AV_Nai], 1.3),
                                                                                                  -1.0) * pow(1.0 + exp(-(ALGEBRAIC[AV_V] - ALGEBRAIC[AV_E_Na] + 110.0) / 20.0),
                                                                                                              -1.0);
    ALGEBRAIC[AV_i_fK] = STATES[If_y_gate_y] * CONSTANTS[AC_g_f_K] * (ALGEBRAIC[AV_V] - CONSTANTS[AC_E_K]) * (1.0 - CONSTANTS[AC_blockade]);
    ALGEBRAIC[AV_i_fNa] = STATES[If_y_gate_y] * CONSTANTS[AC_g_f_Na] * (ALGEBRAIC[AV_V] - ALGEBRAIC[AV_E_Na]) * (1.0 - CONSTANTS[AC_blockade]);
    ALGEBRAIC[AV_i_to_i_to] = CONSTANTS[AC_g_to] * (ALGEBRAIC[AV_V] - CONSTANTS[AC_E_K]) * STATES[Ito_q_gate_q] * STATES[Ito_r_gate_r];
    RATES[Na_i] = (1.0 - CONSTANTS[AC_Nai_clamp]) * -1.0 * (ALGEBRAIC[AV_i_Na_i_Na] + ALGEBRAIC[AV_i_fNa] + ALGEBRAIC[AV_i_siNa] + 3.0 * ALGEBRAIC[AV_i_NaK_i_NaK] + 3.0 * ALGEBRAIC[AV_i_NaCa_i_NaCa]) / (1.0 * (CONSTANTS[AC_V_i] + CONSTANTS[AC_V_sub]) * CONSTANTS[AC_F]);
    ALGEBRAIC[AV_i_f_i_f] = ALGEBRAIC[AV_i_fNa] + ALGEBRAIC[AV_i_fK];
    ALGEBRAIC[Iion_cm] = ALGEBRAIC[AV_i_f_i_f] + ALGEBRAIC[AV_i_Kr_i_Kr] + ALGEBRAIC[AV_i_Ks_i_Ks] + ALGEBRAIC[AV_i_to_i_to] + ALGEBRAIC[AV_i_NaK_i_NaK] + ALGEBRAIC[AV_i_NaCa_i_NaCa] + ALGEBRAIC[AV_i_Na_i_Na] + ALGEBRAIC[AV_i_CaL_i_CaL] + ALGEBRAIC[AV_i_CaT_i_CaT] + ALGEBRAIC[AV_i_KACh_i_KACh] + ALGEBRAIC[AV_i_Kur_i_Kur];
    
    //RATES[membrane_V] = -ALGEBRAIC[AV_i_tot] / CONSTANTS[AC_C];
    

    

    ALGEBRAIC[Istim] = computeIstim(VOI, CONSTANTS);

    if (solveVmWithinODESolver)
    {
        RATES[membrane_V] = -ALGEBRAIC[Iion_cm] - ALGEBRAIC[Istim];
    }

}
