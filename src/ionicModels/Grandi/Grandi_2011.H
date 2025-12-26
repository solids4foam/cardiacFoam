#include "Grandi_2011Names.H"
#include "stimulusIO.H"

static const char* GrandiSTATES_NAMES[NUM_STATES] = {
    "membrane_V",
    "Ca_i",
    "Ca_jn",
    "Ca_sl",
    "Ca_sr",
    "Csqn",
    "Ical_d",
    "Ical_f",
    "Ical_fCaB_jn",
    "Ical_fCaB_sl",
    "Ikr_xr",
    "Iks_xs",
    "Ikur_ikur_r",
    "Ikur_s",
    "Ina_h",
    "Ina_j",
    "Ina_m",
    "Inal_hl",
    "Inal_ml",
    "Ito_x",
    "Ito_y",
    
    "K_i",
    "ryr_i",
    "ryr_o",
    "ryr_ryr_r",
    "Na_i",
    "Na_jn",
    "Na_sl",

    "buffCa_CaM",
    "buffCa_Myoc",
    "buffCa_Myom",
    "buffCa_SLH_jn",
    "buffCa_SLH_sl",
    "buffCa_SLL_jn",
    "buffCa_SLL_sl",
    "buffCa_SRB",
    "buffCa_TnCHc",
    "buffCa_TnCHm",
    "buffCa_TnCL",
    "buffNa_NaB_jn",
    "buffNa_NaB_sl",
};

static const char* GrandiALGEBRAIC_NAMES[NUM_ALGEBRAIC] = {
    "AV_ikr_xr_inf",
    "AV_ikr_xr_tau",
    "AV_ikur_r_inf",
    "AV_ikur_r_tau",
    "AV_ikur_s_inf",
    "AV_ikur_s_tau",
    "AV_ina_h_a",
    "AV_ina_h_b",
    "AV_ina_h_inf",
    "AV_ina_h_tau",
    "AV_ina_j_a",
    "AV_ina_j_b",
    "AV_ina_j_inf",
    "AV_ina_j_tau",
    "AV_ina_m_inf",
    "AV_ina_m_tau",
    "AV_inal_hl_inf",
    "AV_inal_ml_a",
    "AV_inal_ml_b",
    "AV_ito_x_inf",
    "AV_ito_x_tau",
    "AV_ito_y_inf",
    "AV_ito_y_tau",
    "AV_J_CaB_cytosol",
    "AV_J_CaB_jntion",
    "AV_J_CaB_sl",
    "AV_ical_d_inf",
    "AV_ical_d_tau",
    "AV_ical_f_inf",
    "AV_ical_f_tau",
    "AV_iks_xs_inf",
    "AV_iks_xs_tau",
    "AV_fnak",
    "AV_INaK_jn",
    "AV_INaK_sl",
    "AV_INaK",
    "AV_ECa_jn",
    "AV_ECa_sl",
    "AV_EK",
    "AV_ENa_jn",
    "AV_ENa_sl",
    "AV_RI",
    "AV_J_SRCarel",
    "AV_J_SR_leak",
    "AV_J_serca",
    "AV_kCaSR",
    "AV_kiSRCa",
    "AV_koSRCa",
    "AV_ICaB_jn",
    "AV_ICaB_sl",
    "AV_ICaB",
    "AV_ibarca_jn",
    "AV_ibarca_sl",
    "AV_ibark",
    "AV_ibarna_jn",
    "AV_ibarna_sl",
    "AV_ICaL_Ca_jn",
    "AV_ICaL_Ca_sl",
    "AV_ICaL_K",
    "AV_ICaL_Na_jn",
    "AV_ICaL_Na_sl",
    "AV_ICaL_Ca",
    "AV_ICaL_Na",
    "AV_ICaL",
    "AV_IClB",
    "AV_IClCa_jn",
    "AV_IClCa_sl",
    "AV_IClCa",
    "AV_ik1_inf_a",
    "AV_ik1_inf_b",
    "AV_kp",
    "AV_IKp_jn",
    "AV_IKp_sl",
    "AV_IKp",
    "AV_rr",
    "AV_IKr",
    "AV_EKs",
    "AV_IKs_jn",
    "AV_IKs_sl",
    "AV_IKs",
    "AV_IKur",
    "AV_INa_jn",
    "AV_INa_sl",
    "AV_INa",
    "AV_INaB_jn",
    "AV_INaB_sl",
    "AV_INaB",
    "AV_INaL_jn",
    "AV_INaL_sl",
    "AV_INaL",
    "AV_Ito",
    "AV_ik1_inf",
    "AV_IK1",
    "AV_IK_tot",
    "AV_Ka_jn",
    "AV_Ka_sl",
    "AV_s1_jn",
    "AV_s2_jn",
    "AV_s3_jn",
    "AV_s1_sl",
    "AV_s2_sl",
    "AV_s3_sl",
    "AV_ipca_IpCa_jn_a",
    "AV_ipca_IpCa_sl_a",
    "AV_ICl_tot",
    "AV_INaCa_jn",
    "AV_INaCa_sl",
    "AV_IpCa_jn",
    "AV_IpCa_sl",
    "AV_ICa_tot_jn",
    "AV_ICa_tot_sl",
    "AV_INaCa",
    "AV_IpCa",
    "AV_INa_tot_jn",
    "AV_INa_tot_sl",
    "AV_ICa_tot",
    "AV_INa_tot",
    "Istim",
    "Iion_cm",
};


static inline const Foam::HashTable<Foam::wordList>&
GrandiDependencyMap()
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
GrandiinitConsts(double* CONSTANTS, double* RATES, double* STATES, int tissueFlag, const Foam::dictionary& stimulus)
{
/* buffna */
    CONSTANTS[AC_Bmax_Na_jn] = 7.561;
    CONSTANTS[AC_Bmax_Na_sl] = 1.65;
    CONSTANTS[AC_koff_na] = 0.001;
    CONSTANTS[AC_kon_na] = 0.0001;
    
    /* cell */
    CONSTANTS[AC_AF] = 0.0;
    CONSTANTS[AC_C] = 1.1e-10;
    CONSTANTS[AC_ISO] = 0.0;
    CONSTANTS[AC_RA] = 0.0;
    
    /* geom */
    CONSTANTS[AC_JCa_jnsl] = 8.24130542277896849e-13;
    CONSTANTS[AC_JCa_slmyo] = 3.72425607984805052e-12;
    CONSTANTS[AC_JNa_jnsl] = 1.83127823220607955e-14;
    CONSTANTS[AC_JNa_slmyo] = 1.63862792221979433e-12;
    CONSTANTS[AC_cell_length] = 100.0;
    CONSTANTS[AC_cell_radius] = 10.25;
    CONSTANTS[AC_pi] = 3.14159265358979312e+00;
    CONSTANTS[AC_Vcell] = CONSTANTS[AC_pi] * pow(CONSTANTS[AC_cell_radius],
                                                 2.0) * CONSTANTS[AC_cell_length] * 1e-15;
    CONSTANTS[AC_Vjn] = 0.0539 * 0.01 * CONSTANTS[AC_Vcell];
    CONSTANTS[AC_Vmyo] = 0.65 * CONSTANTS[AC_Vcell];
    CONSTANTS[AC_Vsl] = 0.02 * CONSTANTS[AC_Vcell];
    CONSTANTS[AC_Vsr] = 0.035 * CONSTANTS[AC_Vcell];
    
    /* inal_hl */
    CONSTANTS[AC_inal_hl_tau] = 600.0;
    
    /* ion */
    CONSTANTS[AC_Ca_o] = 1.8;
    CONSTANTS[AC_Cl_i] = 15.0;
    CONSTANTS[AC_Cl_o] = 150.0;
    CONSTANTS[AC_K_o] = 5.4;
    CONSTANTS[AC_Mg_i] = 1.0;
    CONSTANTS[AC_Na_o] = 140.0;
    
    /* junc */
    CONSTANTS[AC_Fjn] = 0.11;
    CONSTANTS[AC_Fjn_CaL] = 0.9;
    CONSTANTS[AC_Fsl] = 1.0 - CONSTANTS[AC_Fjn];
    CONSTANTS[AC_Fsl_CaL] = 1.0 - CONSTANTS[AC_Fjn_CaL];
    
    /* phys */
    CONSTANTS[AC_F] = 96485.0;
    CONSTANTS[AC_R] = 8314.0;
    CONSTANTS[AC_T] = 310.0;
    CONSTANTS[AC_FRT] = CONSTANTS[AC_F] / CONSTANTS[AC_R] / CONSTANTS[AC_T];
    CONSTANTS[AC_Q] = (CONSTANTS[AC_T] - 310.0) / 10.0;
    
    /* buffca */
    CONSTANTS[AC_Bmax_CaM] = 0.024;
    CONSTANTS[AC_Bmax_SLhighj] = 0.00165 * CONSTANTS[AC_Vmyo] / CONSTANTS[AC_Vjn] * 0.1;
    CONSTANTS[AC_Bmax_SLhighsl] = CONSTANTS[AC_Vmyo] / CONSTANTS[AC_Vsl] * 0.0134;
    CONSTANTS[AC_Bmax_SLlowj] = 0.0046 * CONSTANTS[AC_Vmyo] / CONSTANTS[AC_Vjn] * 0.1;
    CONSTANTS[AC_Bmax_SLlowsl] = CONSTANTS[AC_Vmyo] / CONSTANTS[AC_Vsl] * 0.0374;
    CONSTANTS[AC_Bmax_SR] = 0.0171;
    CONSTANTS[AC_Bmax_TnChigh] = 0.14;
    CONSTANTS[AC_Bmax_TnClow] = 0.07;
    CONSTANTS[AC_Bmax_myosin] = 0.14;
    CONSTANTS[AC_koff_cam] = 0.238;
    CONSTANTS[AC_koff_myoca] = 0.00046;
    CONSTANTS[AC_koff_myomg] = 5.7e-05;
    CONSTANTS[AC_koff_slh] = 0.03;
    CONSTANTS[AC_koff_sll] = 1.3;
    CONSTANTS[AC_koff_sr] = 0.06;
    CONSTANTS[AC_koff_tnchca] = 3.2e-05;
    CONSTANTS[AC_koff_tnchmg] = 0.00333;
    CONSTANTS[AC_koff_tncl] = (1.0 + 0.5 * CONSTANTS[AC_ISO]) * 0.0196;
    CONSTANTS[AC_kon_cam] = 34.0;
    CONSTANTS[AC_kon_myoca] = 13.8;
    CONSTANTS[AC_kon_myomg] = 0.0157;
    CONSTANTS[AC_kon_slh] = 100.0;
    CONSTANTS[AC_kon_sll] = 100.0;
    CONSTANTS[AC_kon_sr] = 100.0;
    CONSTANTS[AC_kon_tnchca] = 2.37;
    CONSTANTS[AC_kon_tnchmg] = 0.003;
    CONSTANTS[AC_kon_tncl] = 32.7;
    
    /* calcium_Csqn */
    CONSTANTS[AC_Bmax_Csqn] = CONSTANTS[AC_Vmyo] / CONSTANTS[AC_Vsr] * 0.14;
    CONSTANTS[AC_koff_csqn] = 65.0;
    CONSTANTS[AC_kon_csqn] = 100.0;
    
    /* inak */
    CONSTANTS[AC_IbarNaK] = 1.26;
    CONSTANTS[AC_KmKo] = 1.5;
    CONSTANTS[AC_KmNaip] = 11.0 * (1.0 - 0.25 * CONSTANTS[AC_ISO]);
    CONSTANTS[AC_sigma] = (exp(CONSTANTS[AC_Na_o] / 67.3) - 1.0) / 7.0;
    
    /* nernst */
    CONSTANTS[AC_ECl] = 1.0 / CONSTANTS[AC_FRT] * log(CONSTANTS[AC_Cl_i] / CONSTANTS[AC_Cl_o]);
    
    /* ryr */
    CONSTANTS[AC_J_SR_leak_max] = 5.348e-06;
    CONSTANTS[AC_Kmf] = (2.5 - 1.25 * CONSTANTS[AC_ISO]) * 0.000246;
    CONSTANTS[AC_Kmr] = 1.7;
    CONSTANTS[AC_MaxSR] = 15.0;
    CONSTANTS[AC_MinSR] = 1.0;
    CONSTANTS[AC_Q10SRCaP] = 2.6;
    CONSTANTS[AC_Vmax_SRCaP] = 0.0053114;
    CONSTANTS[AC_ec50SR] = 0.45;
    CONSTANTS[AC_hillSRCaP] = 1.787;
    CONSTANTS[AC_kiCa] = 0.5;
    CONSTANTS[AC_kim] = 0.005;
    CONSTANTS[AC_koCa] = (10.0 + 20.0 * CONSTANTS[AC_AF] + 10.0 * CONSTANTS[AC_ISO] * (1.0 - CONSTANTS[AC_AF])) * 1.0;
    CONSTANTS[AC_kom] = 0.06;
    CONSTANTS[AC_ks] = 25.0;
    
    /* stimulus */
    CONSTANTS[AC_amplitude] = -12.5;
    CONSTANTS[AC_stim_duration] = 5.0;
    CONSTANTS[AC_stim_offset] = 50.0;
    CONSTANTS[AC_stim_period] = 1000.0;
    
    /* icab */
    CONSTANTS[AC_gCaB] = 6.06430000000000033e-04;
    
    /* ical */
    CONSTANTS[AC_Q10CaL] = 1.8;
    CONSTANTS[AC_f_conducting] = 1.0;
    CONSTANTS[AC_fcaCaMSL] = 0.0;
    CONSTANTS[AC_fcaCaj] = 0.0;
    CONSTANTS[AC_pCa_max] = 0.00027;
    CONSTANTS[AC_pK_max] = 1.35e-07;
    CONSTANTS[AC_pNa_max] = 7.5e-09;
    CONSTANTS[AC_pCa] = (1.0 + 0.5 * CONSTANTS[AC_ISO]) * (1.0 - 0.5 * CONSTANTS[AC_AF]) * CONSTANTS[AC_pCa_max];
    CONSTANTS[AC_pK] = (1.0 + 0.5 * CONSTANTS[AC_ISO]) * (1.0 - 0.5 * CONSTANTS[AC_AF]) * CONSTANTS[AC_pK_max];
    CONSTANTS[AC_pNa] = (1.0 + 0.5 * CONSTANTS[AC_ISO]) * (1.0 - 0.5 * CONSTANTS[AC_AF]) * CONSTANTS[AC_pNa_max];
    
    /* iclb */
    CONSTANTS[AC_gClB] = 0.009;
    
    /* iclca */
    CONSTANTS[AC_GClCa] = 0.0548;
    CONSTANTS[AC_KdClCa] = 0.1;
    
    /* ikp */
    CONSTANTS[AC_gKp] = 0.002;
    
    /* ikr */
    CONSTANTS[AC_gKr_max] = 0.035;
    CONSTANTS[AC_gKr] = CONSTANTS[AC_gKr_max] * sqrt(CONSTANTS[AC_K_o] / 5.4);
    
    /* iks */
    CONSTANTS[AC_gKs_max] = 0.0035;
    CONSTANTS[AC_pNaK] = 0.01833;
    CONSTANTS[AC_gKs_jn] = (1.0 + CONSTANTS[AC_AF] + 2.0 * CONSTANTS[AC_ISO]) * CONSTANTS[AC_gKs_max];
    CONSTANTS[AC_gKs_sl] = (1.0 + CONSTANTS[AC_AF] + 2.0 * CONSTANTS[AC_ISO]) * CONSTANTS[AC_gKs_max];
    
    /* ikur */
    CONSTANTS[AC_gKur_max] = 0.045;
    CONSTANTS[AC_gKur] = (1.0 - 0.5 * CONSTANTS[AC_AF]) * (1.0 + 2.0 * CONSTANTS[AC_ISO]) * (1.0 + 0.2 * CONSTANTS[AC_RA]) * CONSTANTS[AC_gKur_max];
    
    /* ina */
    CONSTANTS[AC_gNa_max] = 23.0;
    CONSTANTS[AC_gNa] = CONSTANTS[AC_gNa_max] * (1.0 - 0.1 * CONSTANTS[AC_AF]);
    
    /* inab */
    CONSTANTS[AC_gNaB] = 0.000597;
    
    /* inal */
    CONSTANTS[AC_gNaL_max] = 0.0025;
    CONSTANTS[AC_gNaL] = CONSTANTS[AC_gNaL_max] * CONSTANTS[AC_AF];
    
    /* ito */
    CONSTANTS[AC_gto_max] = 0.165;
    CONSTANTS[AC_gto] = (1.0 - 0.7 * CONSTANTS[AC_AF]) * CONSTANTS[AC_gto_max];
    
    /* ik1 */
    CONSTANTS[AC_gK1_max] = 0.0525;
    CONSTANTS[AC_gK1] = (1.0 + CONSTANTS[AC_AF]) * sqrt(CONSTANTS[AC_K_o] / 5.4) * CONSTANTS[AC_gK1_max];
    
    /* inaca */
    CONSTANTS[AC_IbarNCX_max] = 3.15;
    CONSTANTS[AC_Kdact] = 0.000384;
    CONSTANTS[AC_KmCai] = 0.00359;
    CONSTANTS[AC_KmCao] = 1.3;
    CONSTANTS[AC_KmNai] = 12.29;
    CONSTANTS[AC_KmNao] = 87.5;
    CONSTANTS[AC_Q10NCX] = 1.57;
    CONSTANTS[AC_ksat] = 0.27;
    CONSTANTS[AC_nu] = 0.35;
    CONSTANTS[AC_IbarNCX] = (1.0 + 0.4 * CONSTANTS[AC_AF]) * CONSTANTS[AC_IbarNCX_max];
    
    /* ipca */
    CONSTANTS[AC_IbarSLCaP] = 0.0471;
    CONSTANTS[AC_KmPCa] = 0.0005;
    CONSTANTS[AC_Q10SLCaP] = 2.35;
    
    /* ipca_IpCa_jn */
    CONSTANTS[AC_ipca_IpCa_jn_b] = pow(CONSTANTS[AC_KmPCa] * 1.0, 1.6);
    
    /* ipca_IpCa_sl */
    CONSTANTS[AC_ipca_IpCa_sl_b] = pow(CONSTANTS[AC_KmPCa] * 1.0, 1.6);

    /* --- Initial values --- */

    STATES[buffCa_CaM] =  7.02128101897185673e-04;
    STATES[buffCa_Myoc] =  3.94923428392655786e-03;
    STATES[buffCa_Myom] =  1.35538532457244482e-01;
    STATES[buffCa_SLH_jn] =  1.03674364292988680e-01;
    STATES[buffCa_SLH_sl] =  1.90759804527589089e-01;
    STATES[buffCa_SLL_jn] =  1.35640688636079511e-02;
    STATES[buffCa_SLL_sl] =  2.14063418881809235e-02;
    STATES[buffCa_SRB] =  4.45327242854324807e-03;
    STATES[buffCa_TnCHc] =  1.27856586024588575e-01;
    STATES[buffCa_TnCHm] =  5.69999505293381902e-03;
    STATES[buffCa_TnCL] =  1.83143535034222225e-02;
    STATES[buffNa_NaB_jn] =  3.61396062660070427e+00;
    STATES[buffNa_NaB_sl] =  7.88607791910409195e-01;
    STATES[Ca_i] =  2.10808768153058460e-04;
    STATES[Ca_jn] =  3.25814677291117296e-04;
    STATES[Ca_sl] =  2.33018340557575125e-04;
    STATES[Ca_sr] =  5.02305826642838293e-01;
    STATES[Csqn] =  1.13337536953687845e+00;
    STATES[Ical_d] =  2.16850216379767157e-05;
    STATES[Ical_f] =  9.98384427312367095e-01;
    STATES[Ical_fCaB_jn] =  4.49572164109603364e-02;
    STATES[Ical_fCaB_sl] =  3.28512098597005947e-02;
    STATES[Ikr_xr] =  1.31290096227093382e-03;
    STATES[Iks_xs] =  7.49436760722081534e-03;
    STATES[Ikur_ikur_r] =  3.93548562883350357e-04;
    STATES[Ikur_s] =  9.58234428284286399e-01;
    STATES[Ina_h] =  3.15482710277587786e-01;
    STATES[Ina_j] =  2.48034071360795916e-01;
    STATES[Ina_m] =  1.89326933812916480e-02;
    STATES[Inal_hl] =  3.79829335413739144e-02;
    STATES[Inal_ml] =  1.01974216400706526e-02;
    STATES[Ito_x] =  1.37939236359928058e-03;
    STATES[Ito_y] =  9.45874848392074696e-01;
    STATES[membrane_V] = -7.34336366728778671e+01;
    STATES[K_i] = 120.0;
    STATES[ryr_i] =  5.01323282772066123e-07;
    STATES[ryr_o] =  2.01567245823636694e-06;
    STATES[ryr_ryr_r] =  8.00819151705148946e-01;
    STATES[Na_i] =  9.15199678386256998e+00;
    STATES[Na_jn] =  9.15153381546177336e+00;
    STATES[Na_sl] =  9.15182798281732346e+00;
}


void
GrandicomputeVariables(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC, int tissueFlag, bool solveVmWithinODESolver)
{
/* buffna */
    RATES[buffNa_NaB_jn] = CONSTANTS[AC_kon_na] * STATES[Na_jn] * (CONSTANTS[AC_Bmax_Na_jn] - STATES[buffNa_NaB_jn]) - CONSTANTS[AC_koff_na] * STATES[buffNa_NaB_jn];
    RATES[buffNa_NaB_sl] = CONSTANTS[AC_kon_na] * STATES[Na_sl] * (CONSTANTS[AC_Bmax_Na_sl] - STATES[buffNa_NaB_sl]) - CONSTANTS[AC_koff_na] * STATES[buffNa_NaB_sl];
    
    /* ikr_xr */
    ALGEBRAIC[AV_ikr_xr_inf] = 1.0 / (1.0 + exp(-(STATES[membrane_V] + 10.0) / 5.0));
    ALGEBRAIC[AV_ikr_xr_tau] = +550.0 / (1.0 + exp((-22.0 - STATES[membrane_V]) / 9.0)) * 6.0 / (1.0 + exp((STATES[membrane_V] + 11.0) / 9.0)) + 230.0 / (1.0 + exp((STATES[membrane_V] + 40.0) / 20.0));
    
    /* ikur_ikur_r */
    ALGEBRAIC[AV_ikur_r_inf] = 1.0 / (1.0 + exp((STATES[membrane_V] + 6.0) / -8.6));
    ALGEBRAIC[AV_ikur_r_tau] = 9.0 / (1.0 + exp((STATES[membrane_V] + 5.0) / 12.0)) + 0.5;
    
    /* ikur_s */
    ALGEBRAIC[AV_ikur_s_inf] = 1.0 / (1.0 + exp((STATES[membrane_V] + 7.5) / 10.0));
    ALGEBRAIC[AV_ikur_s_tau] = 590.0 / (1.0 + exp((STATES[membrane_V] + 60.0) / 10.0)) + 3050.0;
    
    /* ina_h */
    ALGEBRAIC[AV_ina_h_a] = ((STATES[membrane_V] >= -40.0) ? 0.0 : 0.057 * exp(-(STATES[membrane_V] + 80.0) / 6.8)) * 1.0;
    ALGEBRAIC[AV_ina_h_b] = ((STATES[membrane_V] >= -40.0) ? 0.77 / (0.13 * (1.0 + exp(-(STATES[membrane_V] + 10.66) / 11.1))) : 2.7 * exp(0.079 * STATES[membrane_V]) + 3.1 * pow(10.0, 5.0) * exp(0.3485 * STATES[membrane_V])) * 1.0;
    ALGEBRAIC[AV_ina_h_inf] = 1.0 / pow(1.0 + exp((STATES[membrane_V] + 71.55) / 7.43), 2.0);
    ALGEBRAIC[AV_ina_h_tau] = 1.0 / (ALGEBRAIC[AV_ina_h_a] + ALGEBRAIC[AV_ina_h_b]);
    
    /* ina_j */
    ALGEBRAIC[AV_ina_j_a] = ((STATES[membrane_V] >= -40.0) ? 0.0 : (-2.5428 * pow(10.0, 4.0) * exp(0.2444 * STATES[membrane_V]) - 6.948 * pow(10.0, -6.0) * exp(-0.04391 * STATES[membrane_V])) * (STATES[membrane_V] + 37.78) / (1.0 + exp(0.311 * (STATES[membrane_V] + 79.23))) * 1.0) * 1.0;
    ALGEBRAIC[AV_ina_j_b] = ((STATES[membrane_V] >= -40.0) ? 0.6 * exp(0.057 * STATES[membrane_V]) / (1.0 + exp(-0.1 * (STATES[membrane_V] + 32.0))) : 0.02424 * exp(-0.01052 * STATES[membrane_V]) / (1.0 + exp(-0.1378 * (STATES[membrane_V] + 40.14)))) * 1.0;
    ALGEBRAIC[AV_ina_j_inf] = 1.0 / pow(1.0 + exp((STATES[membrane_V] + 71.55) / 7.43), 2.0);
    ALGEBRAIC[AV_ina_j_tau] = 1.0 / (ALGEBRAIC[AV_ina_j_a] + ALGEBRAIC[AV_ina_j_b]);
    
    /* ina_m */
    ALGEBRAIC[AV_ina_m_inf] = 1.0 / pow(1.0 + exp(-(56.86 + STATES[membrane_V]) / 9.03), 2.0);
    ALGEBRAIC[AV_ina_m_tau] = 0.1292 * exp(-pow((STATES[membrane_V] + 45.79) / 15.54, 2.0)) + 0.06487 * exp(-pow((STATES[membrane_V] - 4.823) / 51.12, 2.0));
    
    /* inal_hl */
    ALGEBRAIC[AV_inal_hl_inf] = 1.0 / (1.0 + exp((STATES[membrane_V] + 91.0) / 6.1));
    
    /* inal_ml */
    ALGEBRAIC[AV_inal_ml_a] = 0.32 * (STATES[membrane_V] + 47.13) / (1.0 - exp(-0.1 * (STATES[membrane_V] + 47.13)));
    ALGEBRAIC[AV_inal_ml_b] = 0.08 * exp(-STATES[membrane_V] / 11.0);
    
    /* ito_x */
    ALGEBRAIC[AV_ito_x_inf] = 1.0 / (1.0 + exp(-(STATES[membrane_V] + 1.0) / 11.0));
    ALGEBRAIC[AV_ito_x_tau] = 3.5 * exp(-pow(STATES[membrane_V] / 30.0, 2.0)) + 1.5;
    
    /* ito_y */
    ALGEBRAIC[AV_ito_y_inf] = 1.0 / (1.0 + exp((STATES[membrane_V] + 40.5) / 11.5));
    ALGEBRAIC[AV_ito_y_tau] = 25.635 * exp(-pow((STATES[membrane_V] + 52.45) / 15.8827, 2.0)) + 24.14;
    
    /* buffca */
    RATES[buffCa_CaM] = CONSTANTS[AC_kon_cam] * STATES[Ca_i] * (CONSTANTS[AC_Bmax_CaM] - STATES[buffCa_CaM]) - CONSTANTS[AC_koff_cam] * STATES[buffCa_CaM];
    RATES[buffCa_Myoc] = CONSTANTS[AC_kon_myoca] * STATES[Ca_i] * (CONSTANTS[AC_Bmax_myosin] - STATES[buffCa_Myoc] - STATES[buffCa_Myom]) - CONSTANTS[AC_koff_myoca] * STATES[buffCa_Myoc];
    RATES[buffCa_Myom] = CONSTANTS[AC_kon_myomg] * CONSTANTS[AC_Mg_i] * (CONSTANTS[AC_Bmax_myosin] - STATES[buffCa_Myoc] - STATES[buffCa_Myom]) - CONSTANTS[AC_koff_myomg] * STATES[buffCa_Myom];
    RATES[buffCa_SLH_jn] = CONSTANTS[AC_kon_slh] * STATES[Ca_jn] * (CONSTANTS[AC_Bmax_SLhighj] - STATES[buffCa_SLH_jn]) - CONSTANTS[AC_koff_slh] * STATES[buffCa_SLH_jn];
    RATES[buffCa_SLH_sl] = CONSTANTS[AC_kon_slh] * STATES[Ca_sl] * (CONSTANTS[AC_Bmax_SLhighsl] - STATES[buffCa_SLH_sl]) - CONSTANTS[AC_koff_slh] * STATES[buffCa_SLH_sl];
    RATES[buffCa_SLL_jn] = CONSTANTS[AC_kon_sll] * STATES[Ca_jn] * (CONSTANTS[AC_Bmax_SLlowj] - STATES[buffCa_SLL_jn]) - CONSTANTS[AC_koff_sll] * STATES[buffCa_SLL_jn];
    RATES[buffCa_SLL_sl] = CONSTANTS[AC_kon_sll] * STATES[Ca_sl] * (CONSTANTS[AC_Bmax_SLlowsl] - STATES[buffCa_SLL_sl]) - CONSTANTS[AC_koff_sll] * STATES[buffCa_SLL_sl];
    RATES[buffCa_SRB] = CONSTANTS[AC_kon_sr] * STATES[Ca_i] * (CONSTANTS[AC_Bmax_SR] - STATES[buffCa_SRB]) - CONSTANTS[AC_koff_sr] * STATES[buffCa_SRB];
    RATES[buffCa_TnCHc] = CONSTANTS[AC_kon_tnchca] * STATES[Ca_i] * (CONSTANTS[AC_Bmax_TnChigh] - STATES[buffCa_TnCHc] - STATES[buffCa_TnCHm]) - CONSTANTS[AC_koff_tnchca] * STATES[buffCa_TnCHc];
    RATES[buffCa_TnCHm] = CONSTANTS[AC_kon_tnchmg] * CONSTANTS[AC_Mg_i] * (CONSTANTS[AC_Bmax_TnChigh] - STATES[buffCa_TnCHc] - STATES[buffCa_TnCHm]) - CONSTANTS[AC_koff_tnchmg] * STATES[buffCa_TnCHm];
    RATES[buffCa_TnCL] = CONSTANTS[AC_kon_tncl] * STATES[Ca_i] * (CONSTANTS[AC_Bmax_TnClow] - STATES[buffCa_TnCL]) - CONSTANTS[AC_koff_tncl] * STATES[buffCa_TnCL];
    ALGEBRAIC[AV_J_CaB_cytosol] = RATES[buffCa_TnCL] + RATES[buffCa_TnCHc] + RATES[buffCa_TnCHm] + RATES[buffCa_CaM] + RATES[buffCa_Myoc] + RATES[buffCa_Myom] + RATES[buffCa_SRB];
    ALGEBRAIC[AV_J_CaB_jntion] = RATES[buffCa_SLL_jn] + RATES[buffCa_SLH_jn];
    ALGEBRAIC[AV_J_CaB_sl] = RATES[buffCa_SLL_sl] + RATES[buffCa_SLH_sl];
    
    /* ical_d */
    ALGEBRAIC[AV_ical_d_inf] = 1.0 / (1.0 + exp(-(STATES[membrane_V] + 3.0 * CONSTANTS[AC_ISO] + 9.0) / 6.0));
    ALGEBRAIC[AV_ical_d_tau] = ALGEBRAIC[AV_ical_d_inf] * (1.0 - exp(-(STATES[membrane_V] + 3.0 * CONSTANTS[AC_ISO] + 9.0) / 6.0)) / (0.035 * (STATES[membrane_V] + 3.0 * CONSTANTS[AC_ISO] + 9.0));
    
    /* ical_f */
    ALGEBRAIC[AV_ical_f_inf] = 1.0 / (1.0 + exp((STATES[membrane_V] + 3.0 * CONSTANTS[AC_ISO] + 30.0) / 7.0)) + 0.2 / (1.0 + exp((50.0 - STATES[membrane_V] - 3.0 * CONSTANTS[AC_ISO]) / 20.0));
    ALGEBRAIC[AV_ical_f_tau] = 1.0 / (0.0197 * exp(-pow(0.0337 * (STATES[membrane_V] + 3.0 * CONSTANTS[AC_ISO] + 25.0), 2.0)) + 0.02);
    
    /* iks_xs */
    ALGEBRAIC[AV_iks_xs_inf] = 1.0 / (1.0 + exp(-(STATES[membrane_V] + 40.0 * CONSTANTS[AC_ISO] + 3.8) / 14.25));
    ALGEBRAIC[AV_iks_xs_tau] = 990.1 / (1.0 + exp(-(STATES[membrane_V] + 40.0 * CONSTANTS[AC_ISO] + 2.436) / 14.12));
    
    /* inak */
    ALGEBRAIC[AV_fnak] = 1.0 / (1.0 + 0.1245 * exp(-0.1 * STATES[membrane_V] * CONSTANTS[AC_FRT]) + 0.0365 * CONSTANTS[AC_sigma] * exp(-STATES[membrane_V] * CONSTANTS[AC_FRT]));
    ALGEBRAIC[AV_INaK_jn] = CONSTANTS[AC_Fjn] * CONSTANTS[AC_IbarNaK] * ALGEBRAIC[AV_fnak] * CONSTANTS[AC_K_o] / (1.0 + pow(CONSTANTS[AC_KmNaip] / STATES[Na_jn], 4.0)) / (CONSTANTS[AC_K_o] + CONSTANTS[AC_KmKo]);
    ALGEBRAIC[AV_INaK_sl] = CONSTANTS[AC_Fsl] * CONSTANTS[AC_IbarNaK] * ALGEBRAIC[AV_fnak] * CONSTANTS[AC_K_o] / (1.0 + pow(CONSTANTS[AC_KmNaip] / STATES[Na_sl], 4.0)) / (CONSTANTS[AC_K_o] + CONSTANTS[AC_KmKo]);
    ALGEBRAIC[AV_INaK] = ALGEBRAIC[AV_INaK_jn] + ALGEBRAIC[AV_INaK_sl];
    
    /* nernst */
    ALGEBRAIC[AV_ECa_jn] = 1.0 / CONSTANTS[AC_FRT] / 2.0 * log(CONSTANTS[AC_Ca_o] / STATES[Ca_jn]);
    ALGEBRAIC[AV_ECa_sl] = 1.0 / CONSTANTS[AC_FRT] / 2.0 * log(CONSTANTS[AC_Ca_o] / STATES[Ca_sl]);
    ALGEBRAIC[AV_EK] = 1.0 / CONSTANTS[AC_FRT] * log(CONSTANTS[AC_K_o] / STATES[K_i]);
    ALGEBRAIC[AV_ENa_jn] = 1.0 / CONSTANTS[AC_FRT] * log(CONSTANTS[AC_Na_o] / STATES[Na_jn]);
    ALGEBRAIC[AV_ENa_sl] = 1.0 / CONSTANTS[AC_FRT] * log(CONSTANTS[AC_Na_o] / STATES[Na_sl]);
    
    /* ryr */
    ALGEBRAIC[AV_RI] = 1.0 - STATES[ryr_ryr_r] - STATES[ryr_o] - STATES[ryr_i];
    ALGEBRAIC[AV_J_SRCarel] = CONSTANTS[AC_ks] * STATES[ryr_o] * (STATES[Ca_sr] - STATES[Ca_jn]);
    ALGEBRAIC[AV_J_SR_leak] = (1.0 + 0.25 * CONSTANTS[AC_AF]) * (STATES[Ca_sr] - STATES[Ca_jn]) * CONSTANTS[AC_J_SR_leak_max];
    ALGEBRAIC[AV_J_serca] = pow(CONSTANTS[AC_Q10SRCaP], CONSTANTS[AC_Q]) * CONSTANTS[AC_Vmax_SRCaP] * (pow(STATES[Ca_i] / CONSTANTS[AC_Kmf], CONSTANTS[AC_hillSRCaP]) - pow(STATES[Ca_sr] / CONSTANTS[AC_Kmr], CONSTANTS[AC_hillSRCaP])) / (1.0 + pow(STATES[Ca_i] / CONSTANTS[AC_Kmf], CONSTANTS[AC_hillSRCaP]) + pow(STATES[Ca_sr] / CONSTANTS[AC_Kmr], CONSTANTS[AC_hillSRCaP]));
    ALGEBRAIC[AV_kCaSR] = CONSTANTS[AC_MaxSR] - (CONSTANTS[AC_MaxSR] - CONSTANTS[AC_MinSR]) / (1.0 + pow(CONSTANTS[AC_ec50SR] / STATES[Ca_sr], 2.5));
    ALGEBRAIC[AV_kiSRCa] = CONSTANTS[AC_kiCa] * ALGEBRAIC[AV_kCaSR];
    ALGEBRAIC[AV_koSRCa] = CONSTANTS[AC_koCa] / ALGEBRAIC[AV_kCaSR];
    RATES[ryr_i] = ALGEBRAIC[AV_kiSRCa] * STATES[Ca_jn] * STATES[ryr_o] - CONSTANTS[AC_kim] * STATES[ryr_i] - (CONSTANTS[AC_kom] * STATES[ryr_i] - ALGEBRAIC[AV_koSRCa] * pow(STATES[Ca_jn], 2.0) * ALGEBRAIC[AV_RI]);
    RATES[ryr_o] = ALGEBRAIC[AV_koSRCa] * pow(STATES[Ca_jn], 2.0) * STATES[ryr_ryr_r] - CONSTANTS[AC_kom] * STATES[ryr_o] - (ALGEBRAIC[AV_kiSRCa] * STATES[Ca_jn] * STATES[ryr_o] - CONSTANTS[AC_kim] * STATES[ryr_i]);
    RATES[ryr_ryr_r] = CONSTANTS[AC_kim] * ALGEBRAIC[AV_RI] - ALGEBRAIC[AV_kiSRCa] * STATES[Ca_jn] * STATES[ryr_ryr_r] - (ALGEBRAIC[AV_koSRCa] * pow(STATES[Ca_jn], 2.0) * STATES[ryr_ryr_r] - CONSTANTS[AC_kom] * STATES[ryr_o]);
    

    /* icab */
    ALGEBRAIC[AV_ICaB_jn] = CONSTANTS[AC_Fjn] * CONSTANTS[AC_gCaB] * (STATES[membrane_V] - ALGEBRAIC[AV_ECa_jn]);
    ALGEBRAIC[AV_ICaB_sl] = CONSTANTS[AC_Fsl] * CONSTANTS[AC_gCaB] * (STATES[membrane_V] - ALGEBRAIC[AV_ECa_sl]);
    ALGEBRAIC[AV_ICaB] = ALGEBRAIC[AV_ICaB_jn] + ALGEBRAIC[AV_ICaB_sl];
    
    /* ical */
    RATES[Ical_d] = (ALGEBRAIC[AV_ical_d_inf] - STATES[Ical_d]) / ALGEBRAIC[AV_ical_d_tau];
    RATES[Ical_f] = (ALGEBRAIC[AV_ical_f_inf] - STATES[Ical_f]) / ALGEBRAIC[AV_ical_f_tau];
    RATES[Ical_fCaB_jn] = (1.7 * STATES[Ca_jn] * (1.0 - STATES[Ical_fCaB_jn]) - 0.0119 * STATES[Ical_fCaB_jn]) * 1.0;
    RATES[Ical_fCaB_sl] = (1.7 * STATES[Ca_sl] * (1.0 - STATES[Ical_fCaB_sl]) - 0.0119 * STATES[Ical_fCaB_sl]) * 1.0;
    ALGEBRAIC[AV_ibarca_jn] = CONSTANTS[AC_pCa] * 4.0 * (STATES[membrane_V] * CONSTANTS[AC_F] * CONSTANTS[AC_FRT]) * (0.341 * STATES[Ca_jn] * exp(2.0 * STATES[membrane_V] * CONSTANTS[AC_FRT]) - 0.341 * CONSTANTS[AC_Ca_o]) / (exp(2.0 * STATES[membrane_V] * CONSTANTS[AC_FRT]) - 1.0) * CONSTANTS[AC_f_conducting];
    ALGEBRAIC[AV_ibarca_sl] = CONSTANTS[AC_pCa] * 4.0 * (STATES[membrane_V] * CONSTANTS[AC_F] * CONSTANTS[AC_FRT]) * (0.341 * STATES[Ca_sl] * exp(2.0 * STATES[membrane_V] * CONSTANTS[AC_FRT]) - 0.341 * CONSTANTS[AC_Ca_o]) / (exp(2.0 * STATES[membrane_V] * CONSTANTS[AC_FRT]) - 1.0) * CONSTANTS[AC_f_conducting];
    ALGEBRAIC[AV_ibark] = CONSTANTS[AC_pK] * STATES[membrane_V] * CONSTANTS[AC_F] * CONSTANTS[AC_FRT] * (0.75 * STATES[K_i] * exp(STATES[membrane_V] * CONSTANTS[AC_FRT]) - 0.75 * CONSTANTS[AC_K_o]) / (exp(STATES[membrane_V] * CONSTANTS[AC_FRT]) - 1.0) * CONSTANTS[AC_f_conducting];
    ALGEBRAIC[AV_ibarna_jn] = CONSTANTS[AC_pNa] * STATES[membrane_V] * CONSTANTS[AC_F] * CONSTANTS[AC_FRT] * (0.75 * STATES[Na_jn] * exp(STATES[membrane_V] * CONSTANTS[AC_FRT]) - 0.75 * CONSTANTS[AC_Na_o]) / (exp(STATES[membrane_V] * CONSTANTS[AC_FRT]) - 1.0) * CONSTANTS[AC_f_conducting];
    ALGEBRAIC[AV_ibarna_sl] = CONSTANTS[AC_pNa] * STATES[membrane_V] * CONSTANTS[AC_F] * CONSTANTS[AC_FRT] * (0.75 * STATES[Na_sl] * exp(STATES[membrane_V] * CONSTANTS[AC_FRT]) - 0.75 * CONSTANTS[AC_Na_o]) / (exp(STATES[membrane_V] * CONSTANTS[AC_FRT]) - 1.0) * CONSTANTS[AC_f_conducting];
    ALGEBRAIC[AV_ICaL_Ca_jn] = CONSTANTS[AC_Fjn_CaL] * ALGEBRAIC[AV_ibarca_jn] * STATES[Ical_d] * STATES[Ical_f] * (1.0 - STATES[Ical_fCaB_jn] + CONSTANTS[AC_fcaCaj]) * pow(CONSTANTS[AC_Q10CaL],
                                                                                                                                     CONSTANTS[AC_Q]) * 0.45;
    ALGEBRAIC[AV_ICaL_Ca_sl] = CONSTANTS[AC_Fsl_CaL] * ALGEBRAIC[AV_ibarca_sl] * STATES[Ical_d] * STATES[Ical_f] * (1.0 - STATES[Ical_fCaB_sl] + CONSTANTS[AC_fcaCaMSL]) * pow(CONSTANTS[AC_Q10CaL],
                                                                                                                                       CONSTANTS[AC_Q]) * 0.45;
    ALGEBRAIC[AV_ICaL_K] = ALGEBRAIC[AV_ibark] * STATES[Ical_d] * STATES[Ical_f] * (CONSTANTS[AC_Fjn_CaL] * (CONSTANTS[AC_fcaCaj] + (1.0 - STATES[Ical_fCaB_jn])) + CONSTANTS[AC_Fsl_CaL] * (CONSTANTS[AC_fcaCaMSL] + (1.0 - STATES[Ical_fCaB_sl]))) * pow(CONSTANTS[AC_Q10CaL],
                                                                                                                                                                                                         CONSTANTS[AC_Q]) * 0.45;
    ALGEBRAIC[AV_ICaL_Na_jn] = CONSTANTS[AC_Fjn_CaL] * ALGEBRAIC[AV_ibarna_jn] * STATES[Ical_d] * STATES[Ical_f] * (1.0 - STATES[Ical_fCaB_jn] + CONSTANTS[AC_fcaCaj]) * pow(CONSTANTS[AC_Q10CaL],
                                                                                                                                     CONSTANTS[AC_Q]) * 0.45;
    ALGEBRAIC[AV_ICaL_Na_sl] = CONSTANTS[AC_Fsl_CaL] * ALGEBRAIC[AV_ibarna_sl] * STATES[Ical_d] * STATES[Ical_f] * (1.0 - STATES[Ical_fCaB_sl] + CONSTANTS[AC_fcaCaMSL]) * pow(CONSTANTS[AC_Q10CaL],
                                                                                                                                       CONSTANTS[AC_Q]) * 0.45;
    ALGEBRAIC[AV_ICaL_Ca] = ALGEBRAIC[AV_ICaL_Ca_jn] + ALGEBRAIC[AV_ICaL_Ca_sl];
    ALGEBRAIC[AV_ICaL_Na] = ALGEBRAIC[AV_ICaL_Na_jn] + ALGEBRAIC[AV_ICaL_Na_sl];
    ALGEBRAIC[AV_ICaL] = ALGEBRAIC[AV_ICaL_Ca] + ALGEBRAIC[AV_ICaL_K] + ALGEBRAIC[AV_ICaL_Na];
    
    /* iclb */
    ALGEBRAIC[AV_IClB] = CONSTANTS[AC_gClB] * (STATES[membrane_V] - CONSTANTS[AC_ECl]);
    
    /* iclca */
    ALGEBRAIC[AV_IClCa_jn] = CONSTANTS[AC_Fjn] * CONSTANTS[AC_GClCa] / (1.0 + CONSTANTS[AC_KdClCa] / STATES[Ca_jn]) * (STATES[membrane_V] - CONSTANTS[AC_ECl]);
    ALGEBRAIC[AV_IClCa_sl] = CONSTANTS[AC_Fsl] * CONSTANTS[AC_GClCa] / (1.0 + CONSTANTS[AC_KdClCa] / STATES[Ca_sl]) * (STATES[membrane_V] - CONSTANTS[AC_ECl]);
    ALGEBRAIC[AV_IClCa] = ALGEBRAIC[AV_IClCa_jn] + ALGEBRAIC[AV_IClCa_sl];
    
    /* ik1_ik1_inf */
    ALGEBRAIC[AV_ik1_inf_a] = 1.0 * (1.02 / (1.0 + exp(0.2385 * (STATES[membrane_V] - ALGEBRAIC[AV_EK] - 59.215))));
    ALGEBRAIC[AV_ik1_inf_b] = 1.0 * ((0.49124 * exp(0.08032 * (STATES[membrane_V] + 5.476 - ALGEBRAIC[AV_EK])) + exp(0.06175 * (STATES[membrane_V] - ALGEBRAIC[AV_EK] - 594.31))) / (1.0 + exp(-0.5143 * (STATES[membrane_V] - ALGEBRAIC[AV_EK] + 4.753))));
    
    /* ikp */
    ALGEBRAIC[AV_kp] = 1.0 / (1.0 + exp(7.488 - STATES[membrane_V] / 5.98));
    ALGEBRAIC[AV_IKp_jn] = CONSTANTS[AC_Fjn] * CONSTANTS[AC_gKp] * ALGEBRAIC[AV_kp] * (STATES[membrane_V] - ALGEBRAIC[AV_EK]);
    ALGEBRAIC[AV_IKp_sl] = CONSTANTS[AC_Fsl] * CONSTANTS[AC_gKp] * ALGEBRAIC[AV_kp] * (STATES[membrane_V] - ALGEBRAIC[AV_EK]);
    ALGEBRAIC[AV_IKp] = ALGEBRAIC[AV_IKp_jn] + ALGEBRAIC[AV_IKp_sl];
    
    /* ikr */
    RATES[Ikr_xr] = (ALGEBRAIC[AV_ikr_xr_inf] - STATES[Ikr_xr]) / ALGEBRAIC[AV_ikr_xr_tau];
    ALGEBRAIC[AV_rr] = 1.0 / (1.0 + exp((STATES[membrane_V] + 74.0) / 24.0));
    ALGEBRAIC[AV_IKr] = CONSTANTS[AC_gKr] * STATES[Ikr_xr] * ALGEBRAIC[AV_rr] * (STATES[membrane_V] - ALGEBRAIC[AV_EK]);
    
    /* iks */
    RATES[Iks_xs] = (ALGEBRAIC[AV_iks_xs_inf] - STATES[Iks_xs]) / ALGEBRAIC[AV_iks_xs_tau];
    ALGEBRAIC[AV_EKs] = 1.0 / CONSTANTS[AC_FRT] * log((CONSTANTS[AC_K_o] + CONSTANTS[AC_pNaK] * CONSTANTS[AC_Na_o]) / (STATES[K_i] + CONSTANTS[AC_pNaK] * STATES[Na_i]));
    ALGEBRAIC[AV_IKs_jn] = CONSTANTS[AC_Fjn] * CONSTANTS[AC_gKs_jn] * pow(STATES[Iks_xs], 2.0) * (STATES[membrane_V] - ALGEBRAIC[AV_EKs]);
    ALGEBRAIC[AV_IKs_sl] = CONSTANTS[AC_Fsl] * CONSTANTS[AC_gKs_sl] * pow(STATES[Iks_xs], 2.0) * (STATES[membrane_V] - ALGEBRAIC[AV_EKs]);
    ALGEBRAIC[AV_IKs] = ALGEBRAIC[AV_IKs_jn] + ALGEBRAIC[AV_IKs_sl];
    
    /* ikur */
    RATES[Ikur_ikur_r] = (ALGEBRAIC[AV_ikur_r_inf] - STATES[Ikur_ikur_r]) / ALGEBRAIC[AV_ikur_r_tau];
    RATES[Ikur_s] = (ALGEBRAIC[AV_ikur_s_inf] - STATES[Ikur_s]) / ALGEBRAIC[AV_ikur_s_tau];
    ALGEBRAIC[AV_IKur] = CONSTANTS[AC_gKur] * STATES[Ikur_ikur_r] * STATES[Ikur_s] * (STATES[membrane_V] - ALGEBRAIC[AV_EK]);
    
    /* ina */
    RATES[Ina_h] = (ALGEBRAIC[AV_ina_h_inf] - STATES[Ina_h]) / ALGEBRAIC[AV_ina_h_tau];
    RATES[Ina_j] = (ALGEBRAIC[AV_ina_j_inf] - STATES[Ina_j]) / ALGEBRAIC[AV_ina_j_tau];
    RATES[Ina_m] = (ALGEBRAIC[AV_ina_m_inf] - STATES[Ina_m]) / ALGEBRAIC[AV_ina_m_tau];
    ALGEBRAIC[AV_INa_jn] = CONSTANTS[AC_Fjn] * CONSTANTS[AC_gNa] * pow(STATES[Ina_m], 3.0) * STATES[Ina_h] * STATES[Ina_j] * (STATES[membrane_V] - ALGEBRAIC[AV_ENa_jn]);
    ALGEBRAIC[AV_INa_sl] = CONSTANTS[AC_Fsl] * CONSTANTS[AC_gNa] * pow(STATES[Ina_m], 3.0) * STATES[Ina_h] * STATES[Ina_j] * (STATES[membrane_V] - ALGEBRAIC[AV_ENa_sl]);
    ALGEBRAIC[AV_INa] = ALGEBRAIC[AV_INa_jn] + ALGEBRAIC[AV_INa_sl];
    
    /* inab */
    ALGEBRAIC[AV_INaB_jn] = CONSTANTS[AC_Fjn] * CONSTANTS[AC_gNaB] * (STATES[membrane_V] - ALGEBRAIC[AV_ENa_jn]);
    ALGEBRAIC[AV_INaB_sl] = CONSTANTS[AC_Fsl] * CONSTANTS[AC_gNaB] * (STATES[membrane_V] - ALGEBRAIC[AV_ENa_sl]);
    ALGEBRAIC[AV_INaB] = ALGEBRAIC[AV_INaB_jn] + ALGEBRAIC[AV_INaB_sl];
    
    /* inal */
    RATES[Inal_hl] = (ALGEBRAIC[AV_inal_hl_inf] - STATES[Inal_hl]) / CONSTANTS[AC_inal_hl_tau];
    RATES[Inal_ml] = ALGEBRAIC[AV_inal_ml_a] * (1.0 - STATES[Inal_ml]) - ALGEBRAIC[AV_inal_ml_b] * STATES[Inal_ml];
    ALGEBRAIC[AV_INaL_jn] = CONSTANTS[AC_Fjn] * CONSTANTS[AC_gNaL] * pow(STATES[Inal_ml], 3.0) * STATES[Inal_hl] * (STATES[membrane_V] - ALGEBRAIC[AV_ENa_jn]);
    ALGEBRAIC[AV_INaL_sl] = CONSTANTS[AC_Fsl] * CONSTANTS[AC_gNaL] * pow(STATES[Inal_ml], 3.0) * STATES[Inal_hl] * (STATES[membrane_V] - ALGEBRAIC[AV_ENa_sl]);
    ALGEBRAIC[AV_INaL] = ALGEBRAIC[AV_INaL_jn] + ALGEBRAIC[AV_INaL_sl];
    
    /* ito */
    RATES[Ito_x] = (ALGEBRAIC[AV_ito_x_inf] - STATES[Ito_x]) / ALGEBRAIC[AV_ito_x_tau];
    RATES[Ito_y] = (ALGEBRAIC[AV_ito_y_inf] - STATES[Ito_y]) / ALGEBRAIC[AV_ito_y_tau];
    ALGEBRAIC[AV_Ito] = CONSTANTS[AC_gto] * STATES[Ito_x] * STATES[Ito_y] * (STATES[membrane_V] - ALGEBRAIC[AV_EK]);
    
    /* ik1 */
    ALGEBRAIC[AV_ik1_inf] = ALGEBRAIC[AV_ik1_inf_a] / (ALGEBRAIC[AV_ik1_inf_a] + ALGEBRAIC[AV_ik1_inf_b]);
    ALGEBRAIC[AV_IK1] = CONSTANTS[AC_gK1] * ALGEBRAIC[AV_ik1_inf] * (STATES[membrane_V] - ALGEBRAIC[AV_EK]);
    
    
    /* calcium */
    RATES[Ca_i] = -ALGEBRAIC[AV_J_serca] * CONSTANTS[AC_Vsr] / CONSTANTS[AC_Vmyo] - ALGEBRAIC[AV_J_CaB_cytosol] + CONSTANTS[AC_JCa_slmyo] / CONSTANTS[AC_Vmyo] * (STATES[Ca_sl] - STATES[Ca_i]);
    RATES[Csqn] = CONSTANTS[AC_kon_csqn] * STATES[Ca_sr] * (CONSTANTS[AC_Bmax_Csqn] - STATES[Csqn]) - CONSTANTS[AC_koff_csqn] * STATES[Csqn];
    RATES[Ca_sr] = ALGEBRAIC[AV_J_serca] - (ALGEBRAIC[AV_J_SR_leak] * CONSTANTS[AC_Vmyo] / CONSTANTS[AC_Vsr] + ALGEBRAIC[AV_J_SRCarel]) - RATES[Csqn];
    
    /* inaca */
    ALGEBRAIC[AV_Ka_jn] = 1.0 / (1.0 + pow(CONSTANTS[AC_Kdact] / STATES[Ca_jn], 2.0));
    ALGEBRAIC[AV_Ka_sl] = 1.0 / (1.0 + pow(CONSTANTS[AC_Kdact] / STATES[Ca_sl], 2.0));
    
    /* inaca_INaCa_jn */
    ALGEBRAIC[AV_s1_jn] = exp(CONSTANTS[AC_nu] * STATES[membrane_V] * CONSTANTS[AC_FRT]) * pow(STATES[Na_jn], 3.0) * CONSTANTS[AC_Ca_o];
    ALGEBRAIC[AV_s2_jn] = exp((CONSTANTS[AC_nu] - 1.0) * STATES[membrane_V] * CONSTANTS[AC_FRT]) * pow(CONSTANTS[AC_Na_o],
                                                                                    3.0) * STATES[Ca_jn];
    ALGEBRAIC[AV_s3_jn] = CONSTANTS[AC_KmCai] * pow(CONSTANTS[AC_Na_o], 3.0) * (1.0 + pow(STATES[Na_jn] / CONSTANTS[AC_KmNai], 3.0)) + pow(CONSTANTS[AC_KmNao],
                                                                                                                             3.0) * STATES[Ca_jn] * (1.0 + STATES[Ca_jn] / CONSTANTS[AC_KmCai]) + CONSTANTS[AC_KmCao] * pow(STATES[Na_jn],
                                                                                                                                                                         3.0) + pow(STATES[Na_jn],
                                                                                                                                                                                    3.0) * CONSTANTS[AC_Ca_o] + pow(CONSTANTS[AC_Na_o],
                                                                                                                                                                                                                    3.0) * STATES[Ca_jn];
    
    /* inaca_INaCa_sl */
    ALGEBRAIC[AV_s1_sl] = exp(CONSTANTS[AC_nu] * STATES[membrane_V] * CONSTANTS[AC_FRT]) * pow(STATES[Na_sl], 3.0) * CONSTANTS[AC_Ca_o];
    ALGEBRAIC[AV_s2_sl] = exp((CONSTANTS[AC_nu] - 1.0) * STATES[membrane_V] * CONSTANTS[AC_FRT]) * pow(CONSTANTS[AC_Na_o],
                                                                                    3.0) * STATES[Ca_sl];
    ALGEBRAIC[AV_s3_sl] = CONSTANTS[AC_KmCai] * pow(CONSTANTS[AC_Na_o], 3.0) * (1.0 + pow(STATES[Na_sl] / CONSTANTS[AC_KmNai], 3.0)) + pow(CONSTANTS[AC_KmNao],
                                                                                                                             3.0) * STATES[Ca_sl] * (1.0 + STATES[Ca_sl] / CONSTANTS[AC_KmCai]) + CONSTANTS[AC_KmCao] * pow(STATES[Na_sl],
                                                                                                                                                                         3.0) + pow(STATES[Na_sl],
                                                                                                                                                                                    3.0) * CONSTANTS[AC_Ca_o] + pow(CONSTANTS[AC_Na_o],
                                                                                                                                                                                                                    3.0) * STATES[Ca_sl];
    /* ipca_IpCa_jn */
    ALGEBRAIC[AV_ipca_IpCa_jn_a] = pow(STATES[Ca_jn] * 1.0, 1.6);
    
    /* ipca_IpCa_sl */
    ALGEBRAIC[AV_ipca_IpCa_sl_a] = pow(STATES[Ca_sl] * 1.0, 1.6);
    
    /* sodium */
    RATES[Na_i] = CONSTANTS[AC_JNa_slmyo] / CONSTANTS[AC_Vmyo] * (STATES[Na_sl] - STATES[Na_i]);
    
    /* *remaining* */
    ALGEBRAIC[AV_INaCa_jn] = CONSTANTS[AC_Fjn] * CONSTANTS[AC_IbarNCX] * pow(CONSTANTS[AC_Q10NCX],
                                                                  CONSTANTS[AC_Q]) * ALGEBRAIC[AV_Ka_jn] * (ALGEBRAIC[AV_s1_jn] - ALGEBRAIC[AV_s2_jn]) / ALGEBRAIC[AV_s3_jn] / (1.0 + CONSTANTS[AC_ksat] * exp((CONSTANTS[AC_nu] - 1.0) * STATES[membrane_V] * CONSTANTS[AC_FRT]));
    ALGEBRAIC[AV_INaCa_sl] = CONSTANTS[AC_Fsl] * CONSTANTS[AC_IbarNCX] * pow(CONSTANTS[AC_Q10NCX],
                                                                  CONSTANTS[AC_Q]) * ALGEBRAIC[AV_Ka_sl] * (ALGEBRAIC[AV_s1_sl] - ALGEBRAIC[AV_s2_sl]) / ALGEBRAIC[AV_s3_sl] / (1.0 + CONSTANTS[AC_ksat] * exp((CONSTANTS[AC_nu] - 1.0) * STATES[membrane_V] * CONSTANTS[AC_FRT]));
    ALGEBRAIC[AV_IpCa_jn] = CONSTANTS[AC_Fjn] * pow(CONSTANTS[AC_Q10SLCaP],
                                         CONSTANTS[AC_Q]) * CONSTANTS[AC_IbarSLCaP] * ALGEBRAIC[AV_ipca_IpCa_jn_a] / (ALGEBRAIC[AV_ipca_IpCa_jn_a] + CONSTANTS[AC_ipca_IpCa_jn_b]);
    ALGEBRAIC[AV_IpCa_sl] = CONSTANTS[AC_Fsl] * pow(CONSTANTS[AC_Q10SLCaP],
                                         CONSTANTS[AC_Q]) * CONSTANTS[AC_IbarSLCaP] * ALGEBRAIC[AV_ipca_IpCa_sl_a] / (ALGEBRAIC[AV_ipca_IpCa_sl_a] + CONSTANTS[AC_ipca_IpCa_sl_b]);
    ALGEBRAIC[AV_ICa_tot_jn] = ALGEBRAIC[AV_ICaL_Ca_jn] + ALGEBRAIC[AV_ICaB_jn] + ALGEBRAIC[AV_IpCa_jn] - 2.0 * ALGEBRAIC[AV_INaCa_jn];
    ALGEBRAIC[AV_ICa_tot_sl] = ALGEBRAIC[AV_ICaL_Ca_sl] + ALGEBRAIC[AV_ICaB_sl] + ALGEBRAIC[AV_IpCa_sl] - 2.0 * ALGEBRAIC[AV_INaCa_sl];
    ALGEBRAIC[AV_INaCa] = ALGEBRAIC[AV_INaCa_jn] + ALGEBRAIC[AV_INaCa_sl];
    ALGEBRAIC[AV_IpCa] = ALGEBRAIC[AV_IpCa_jn] + ALGEBRAIC[AV_IpCa_sl];
    ALGEBRAIC[AV_INa_tot_jn] = ALGEBRAIC[AV_INa_jn] + ALGEBRAIC[AV_INaB_jn] + 3.0 * ALGEBRAIC[AV_INaCa_jn] + 3.0 * ALGEBRAIC[AV_INaK_jn] + ALGEBRAIC[AV_ICaL_Na_jn] + ALGEBRAIC[AV_INaL_jn];
    ALGEBRAIC[AV_INa_tot_sl] = ALGEBRAIC[AV_INa_sl] + ALGEBRAIC[AV_INaB_sl] + 3.0 * ALGEBRAIC[AV_INaCa_sl] + 3.0 * ALGEBRAIC[AV_INaK_sl] + ALGEBRAIC[AV_ICaL_Na_sl] + ALGEBRAIC[AV_INaL_sl];
    RATES[Ca_jn] = -ALGEBRAIC[AV_ICa_tot_jn] * CONSTANTS[AC_C] / (CONSTANTS[AC_Vjn] * 2.0 * CONSTANTS[AC_F]) + CONSTANTS[AC_JCa_jnsl] / CONSTANTS[AC_Vjn] * (STATES[Ca_sl] - STATES[Ca_jn]) - ALGEBRAIC[AV_J_CaB_jntion] + ALGEBRAIC[AV_J_SRCarel] * CONSTANTS[AC_Vsr] / CONSTANTS[AC_Vjn] + ALGEBRAIC[AV_J_SR_leak] * CONSTANTS[AC_Vmyo] / CONSTANTS[AC_Vjn];
    RATES[Ca_sl] = -ALGEBRAIC[AV_ICa_tot_sl] * CONSTANTS[AC_C] / (CONSTANTS[AC_Vsl] * 2.0 * CONSTANTS[AC_F]) + CONSTANTS[AC_JCa_jnsl] / CONSTANTS[AC_Vsl] * (STATES[Ca_jn] - STATES[Ca_sl]) + CONSTANTS[AC_JCa_slmyo] / CONSTANTS[AC_Vsl] * (STATES[Ca_i] - STATES[Ca_sl]) - ALGEBRAIC[AV_J_CaB_sl];
    
    RATES[Na_jn] = -ALGEBRAIC[AV_INa_tot_jn] * CONSTANTS[AC_C] / (CONSTANTS[AC_Vjn] * CONSTANTS[AC_F]) + CONSTANTS[AC_JNa_jnsl] / CONSTANTS[AC_Vjn] * (STATES[Na_sl] - STATES[Na_jn]) - RATES[buffNa_NaB_jn];
    RATES[Na_sl] = -ALGEBRAIC[AV_INa_tot_sl] * CONSTANTS[AC_C] / (CONSTANTS[AC_Vsl] * CONSTANTS[AC_F]) + CONSTANTS[AC_JNa_jnsl] / CONSTANTS[AC_Vsl] * (STATES[Na_jn] - STATES[Na_sl]) + CONSTANTS[AC_JNa_slmyo] / CONSTANTS[AC_Vsl] * (STATES[Na_i] - STATES[Na_sl]) - RATES[buffNa_NaB_sl];

    /* potassium */
    RATES[K_i] = 0.0;
    ALGEBRAIC[AV_IK_tot] = ALGEBRAIC[AV_Ito] + ALGEBRAIC[AV_IKr] + ALGEBRAIC[AV_IKs] + ALGEBRAIC[AV_IK1] - 2.0 * ALGEBRAIC[AV_INaK] + ALGEBRAIC[AV_ICaL_K] + ALGEBRAIC[AV_IKp] + ALGEBRAIC[AV_IKur];
    ALGEBRAIC[AV_ICl_tot] = ALGEBRAIC[AV_IClCa] + ALGEBRAIC[AV_IClB];
    ALGEBRAIC[AV_ICa_tot] = ALGEBRAIC[AV_ICa_tot_jn] + ALGEBRAIC[AV_ICa_tot_sl];
    ALGEBRAIC[AV_INa_tot] = ALGEBRAIC[AV_INa_tot_jn] + ALGEBRAIC[AV_INa_tot_sl];
    ALGEBRAIC[Iion_cm] = ALGEBRAIC[AV_INa_tot] + ALGEBRAIC[AV_ICl_tot] + ALGEBRAIC[AV_ICa_tot] + ALGEBRAIC[AV_IK_tot];

    
    ALGEBRAIC[Istim] = computeIstim(VOI, CONSTANTS);

    if (solveVmWithinODESolver)
    {
        RATES[membrane_V] = -ALGEBRAIC[Iion_cm] - ALGEBRAIC[Istim];
    }

}
