
enum CONSTANTS_INDEX{ 
    AC_Acap, AC_Ageo, AC_BSLmax, AC_BSRmax, AC_CaMKo, AC_ECl, AC_F, AC_GK1,
    AC_GKb, AC_GKr, AC_GKs, AC_GNa, AC_GNaL, AC_Gncx, AC_GpCa, AC_Gto,
    AC_H, AC_Khp, AC_Kki, AC_Kko, AC_KmBSL, AC_KmBSR, AC_KmCaAct, AC_KmCaM,
    AC_KmCaMK, AC_Kmgatp, AC_Knai0, AC_Knao0, AC_Knap, AC_Kxkur, AC_L, AC_MgADP,
    AC_MgATP, AC_PCa, AC_PCaK, AC_PCaNa, AC_PCab, AC_PKNa, AC_PNab, AC_Pnak,
    AC_R, AC_SOICR, AC_T, AC_a2, AC_a4, AC_aCaMK, AC_b1, AC_bCaMK,
    AC_cao, AC_cli, AC_clo, AC_cmdnmax, AC_csqnmax, AC_delta, AC_eP, AC_grelbarjsrol,
    AC_h10, AC_h101, AC_h11, AC_h1111, AC_h12, AC_h121, AC_k1, AC_k11,
    AC_k1m, AC_k1p, AC_k2, AC_k21, AC_k2m, AC_k2p, AC_k3m, AC_k3p2,
    AC_k4m, AC_k4p2, AC_k5, AC_k51, AC_kasymm, AC_kcaoff, AC_kcaon, AC_kmcmdn,
    AC_kmcsqn, AC_kmtrpn, AC_kna1, AC_kna2, AC_kna3, AC_ko, AC_nao, AC_pi,
    AC_qca, AC_qna, AC_rad, AC_tau_gap, AC_tau_hl, AC_tauoff, AC_tauon,
    AC_trpnmax, AC_vcell, AC_vcsr, AC_vhalf_d, AC_vhalff, AC_vjsr, AC_vmyo, AC_vmyo1,
    AC_vmyo1frac, AC_vmyo2, AC_vnsr, AC_vnsr1, AC_vnsr2, AC_vss, AC_wca, AC_wna,
    AC_wnaca, AC_zca, AC_zk, AC_zna, 
    AC_stimulus_start,AC_stimulus_period , AC_stimulus_duration, AC_stimulus_amplitude, NUM_CONSTANTS
}; 

enum STATES_INDEX{
    cell_v,
    ionic_concentrations_nai,
    ionic_concentrations_nass,
    ionic_concentrations_ki,
    ionic_concentrations_kss,
    ionic_concentrations_cai,
    ionic_concentrations_cai2,
    ionic_concentrations_cass,
    ionic_concentrations_cansr,
    ionic_concentrations_cajsr,
    ionic_concentrations_cacsr,
    I_Na_m,
    I_Na_h,
    I_Na_j,
    INaL_ml,
    INaL_hl,
    ICaL_d,
    ICaL_fca,
    IKr_xr,
    IKs_xs1,
    IKs_xs2,
    ITo_aa,
    CICR_Jrel2,
    CICR_Jrel1,
    CaMK_CaMKt,
    CICR_tjsrol,
    CICR_A,
    ICaL_fs,
    ICaL_ff,
    NUM_STATES
};

const char* STATES_NAMES[NUM_STATES] = {
    "cell_v",
    "ionic_concentrations_nai",
    "ionic_concentrations_nass",
    "ionic_concentrations_ki",
    "ionic_concentrations_kss",
    "ionic_concentrations_cai",
    "ionic_concentrations_cai2",
    "ionic_concentrations_cass",
    "ionic_concentrations_cansr",
    "ionic_concentrations_cajsr",
    "ionic_concentrations_cacsr",
    "I_Na_m",
    "I_Na_h",
    "I_Na_j",
    "INaL_ml",
    "INaL_hl",
    "ICaL_d",
    "ICaL_fca",
    "IKr_xr",
    "IKs_xs1",
    "IKs_xs2",
    "ITo_aa",
    "CICR_Jrel2",
    "CICR_Jrel1",
    "CaMK_CaMKt",
    "CICR_tjsrol",
    "CICR_A",
    "ICaL_fs",
    "ICaL_ff",
};


enum ALGEBRAIC_INDEX{  
    AV_Iion,
    AV_I_stim,  
    AV_ICaK,    
    AV_ICaL_ICaL,    
    AV_ICaNa,    
    AV_ICab_ICab,    
    AV_IK1_IK1,    
    AV_IKb_IKb,    
    AV_IKr_IKr,    
    AV_IKs_IKs,    
    AV_INa,    
    AV_INaCa_i_INaCa_i,    
    AV_INaCa_ss_INaCa_ss, 
    AV_INaCa,   
    AV_INaK_INaK,    
    AV_INaL_INaL,    
    AV_INab_INab,   
     AV_ITo_ITo,    
     AV_IpCa_IpCa,
     NUM_ALGEBRAIC
};

static const char* ALGEBRAIC_NAMES[NUM_ALGEBRAIC] = {
    "AV_Iion",
    "AV_I_stim",
    "AV_ICaK",
    "AV_ICaL_ICaL",
    "AV_ICaNa",
    "AV_ICab_ICab",
    "AV_IK1_IK1",
    "AV_IKb_IKb",
    "AV_IKr_IKr",
    "AV_IKs_IKs",
    "AV_INa",
    "AV_INaCa_i_INaCa_i",
    "AV_INaCa_ss_INaCa_ss",
    "AV_INaCa",
    "AV_INaK_INaK",
    "AV_INaL_INaL",
    "AV_INab_INab",
    "AV_ITo_ITo",
    "AV_IpCa_IpCa"
};



void
GaurinitConsts(double* CONSTANTS, double* RATES, double *STATES, int tissueFlag)
{ 
    CONSTANTS[AC_stimulus_start]     = 20.0;     // Start time of pulse in ms
    CONSTANTS[AC_stimulus_duration]  = 0.5;      // Pulse duration in ms
    CONSTANTS[AC_stimulus_period]    = 500.0;    // Time between pulses in ms
    CONSTANTS[AC_stimulus_amplitude] = -80.0;    // Amplitude of stimulus
        
    CONSTANTS[AC_F] = 96485.0;
    CONSTANTS[AC_R] = 8314.0;
    CONSTANTS[AC_T] = 310.0;

    /* cell */

    CONSTANTS[AC_cao] = 5.4;
    CONSTANTS[AC_cli] = 19.53;
    CONSTANTS[AC_clo] = 100.0;
    CONSTANTS[AC_ko] = 5.4;
    CONSTANTS[AC_nao] = 140.0;

    CONSTANTS[AC_L] = 0.017;
    CONSTANTS[AC_pi] = 3.14;
    CONSTANTS[AC_rad] = 0.0011;
    CONSTANTS[AC_vmyo1frac] = 0.4;
    CONSTANTS[AC_Ageo] = 2.0 * CONSTANTS[AC_pi] * CONSTANTS[AC_rad] * CONSTANTS[AC_rad] + 2.0 * CONSTANTS[AC_pi] * CONSTANTS[AC_rad] * CONSTANTS[AC_L];
    CONSTANTS[AC_vcell] = 1000.0 * CONSTANTS[AC_pi] * CONSTANTS[AC_rad] * CONSTANTS[AC_rad] * CONSTANTS[AC_L];
    CONSTANTS[AC_Acap] = 2.0 * CONSTANTS[AC_Ageo];
    CONSTANTS[AC_vjsr] = 0.0048 * CONSTANTS[AC_vcell] / 2.0;
    CONSTANTS[AC_vmyo] = 0.68 * CONSTANTS[AC_vcell];

    CONSTANTS[AC_vnsr] = 0.0552 * CONSTANTS[AC_vcell];
    CONSTANTS[AC_vss] = 0.02 * CONSTANTS[AC_vcell] / 2.0;
    CONSTANTS[AC_vcsr] = CONSTANTS[AC_vjsr];
    CONSTANTS[AC_vmyo1] = CONSTANTS[AC_vmyo] * CONSTANTS[AC_vmyo1frac];
    CONSTANTS[AC_vnsr1] = CONSTANTS[AC_vnsr] / 2.0;
    CONSTANTS[AC_vnsr2] = CONSTANTS[AC_vnsr] / 2.0;
    CONSTANTS[AC_vmyo2] = CONSTANTS[AC_vmyo] - CONSTANTS[AC_vmyo1];
    
    /* CaMK */
    CONSTANTS[AC_CaMKo] = 0.05;
    CONSTANTS[AC_ECl] = CONSTANTS[AC_R] * CONSTANTS[AC_T] / CONSTANTS[AC_F] * log(CONSTANTS[AC_cli] / CONSTANTS[AC_clo]);
    CONSTANTS[AC_KmCaM] = 0.0015;
    CONSTANTS[AC_KmCaMK] = 0.065;
    CONSTANTS[AC_PKNa] = 0.01833;
    CONSTANTS[AC_aCaMK] = 0.05;
    CONSTANTS[AC_bCaMK] = 0.00068;
    
    /* I_Na */
    CONSTANTS[AC_GNa] = 14.838;
    
    /* INaL */
    CONSTANTS[AC_GNaL] = 0.0075;
    CONSTANTS[AC_tau_hl] = 600.0;
    
    /* ICaL */
    CONSTANTS[AC_PCa] = 0.0002;
    CONSTANTS[AC_vhalf_d] = 3.4;
    CONSTANTS[AC_vhalff] = -22.9;
    CONSTANTS[AC_zca] = 2.0;
    CONSTANTS[AC_PCaK] = 0.0003574 * CONSTANTS[AC_PCa];
    CONSTANTS[AC_PCaNa] = 0.00125 * CONSTANTS[AC_PCa];
    
    /* IKr */
    CONSTANTS[AC_GKr] = 0.5 * 0.015;
    
    /* IKs */
    CONSTANTS[AC_GKs] = 0.021;
    
    /* IK1 */
    CONSTANTS[AC_GK1] = 0.75 * 0.2;
    
    /* ITo */
    CONSTANTS[AC_Gto] = 0.4 * 0.5;
    
    /* CICR */
    CONSTANTS[AC_SOICR] = 10.0;
    CONSTANTS[AC_grelbarjsrol] = 2.0;
    CONSTANTS[AC_tau_gap] = 3.0;
    CONSTANTS[AC_tauoff] = 5.0;
    CONSTANTS[AC_tauon] = 0.5;
    
    /* INaCa_i */
    CONSTANTS[AC_Gncx] = 1.0 * 0.0008;
    CONSTANTS[AC_KmCaAct] = 0.00015;
    CONSTANTS[AC_kasymm] = 12.5;
    CONSTANTS[AC_kcaoff] = 5000.0;
    CONSTANTS[AC_kcaon] = 1500000.0;
    CONSTANTS[AC_kna1] = 15.0;
    CONSTANTS[AC_kna2] = 5.0;
    CONSTANTS[AC_kna3] = 88.12;
    CONSTANTS[AC_qca] = 0.167;
    CONSTANTS[AC_qna] = 0.5224;
    CONSTANTS[AC_wca] = 60000.0;
    CONSTANTS[AC_wna] = 60000.0;
    CONSTANTS[AC_wnaca] = 5000.0;
    CONSTANTS[AC_zna] = 1.0;
    CONSTANTS[AC_h10] = CONSTANTS[AC_kasymm] + 1.0 + CONSTANTS[AC_nao] / CONSTANTS[AC_kna1] * (1.0 + CONSTANTS[AC_nao] / CONSTANTS[AC_kna2]);
    CONSTANTS[AC_k2] = CONSTANTS[AC_kcaoff];
    CONSTANTS[AC_k5] = CONSTANTS[AC_kcaoff];
    CONSTANTS[AC_h11] = CONSTANTS[AC_nao] * CONSTANTS[AC_nao] / (CONSTANTS[AC_h10] * CONSTANTS[AC_kna1] * CONSTANTS[AC_kna2]);
    CONSTANTS[AC_h12] = 1.0 / CONSTANTS[AC_h10];
    CONSTANTS[AC_k1] = CONSTANTS[AC_h12] * CONSTANTS[AC_cao] * CONSTANTS[AC_kcaon];
    
    /* INaCa_ss */
    CONSTANTS[AC_h101] = CONSTANTS[AC_kasymm] + 1.0 + CONSTANTS[AC_nao] / CONSTANTS[AC_kna1] * (1.0 + CONSTANTS[AC_nao] / CONSTANTS[AC_kna2]);
    CONSTANTS[AC_k21] = CONSTANTS[AC_kcaoff];
    CONSTANTS[AC_k51] = CONSTANTS[AC_kcaoff];
    CONSTANTS[AC_h1111] = CONSTANTS[AC_nao] * CONSTANTS[AC_nao] / (CONSTANTS[AC_h101] * CONSTANTS[AC_kna1] * CONSTANTS[AC_kna2]);
    CONSTANTS[AC_h121] = 1.0 / CONSTANTS[AC_h101];
    CONSTANTS[AC_k11] = CONSTANTS[AC_h121] * CONSTANTS[AC_cao] * CONSTANTS[AC_kcaon];
    
    /* INaK */
    CONSTANTS[AC_H] = 1e-07;
    CONSTANTS[AC_Khp] = 1.698e-07;
    CONSTANTS[AC_Kki] = 0.5;
    CONSTANTS[AC_Kko] = 0.3582;
    CONSTANTS[AC_Kmgatp] = 1.698e-07;
    CONSTANTS[AC_Knai0] = 9.073;
    CONSTANTS[AC_Knao0] = 27.78;
    CONSTANTS[AC_Knap] = 224.0;
    CONSTANTS[AC_Kxkur] = 292.0;
    CONSTANTS[AC_MgADP] = 0.05;
    CONSTANTS[AC_MgATP] = 9.8;
    CONSTANTS[AC_Pnak] = 1.0 * 30.0;
    CONSTANTS[AC_delta] = -0.155;
    CONSTANTS[AC_eP] = 4.2;
    CONSTANTS[AC_k1m] = 182.4;
    CONSTANTS[AC_k1p] = 949.5;
    CONSTANTS[AC_k2m] = 39.4;
    CONSTANTS[AC_k2p] = 687.2;
    CONSTANTS[AC_k3m] = 79300.0;
    CONSTANTS[AC_k3p2] = 1899.0;
    CONSTANTS[AC_k4m] = 40.0;
    CONSTANTS[AC_k4p2] = 639.0;
    CONSTANTS[AC_zk] = 1.0;
    CONSTANTS[AC_a2] = CONSTANTS[AC_k2p];
    CONSTANTS[AC_a4] = CONSTANTS[AC_k4p2] * CONSTANTS[AC_MgATP] / CONSTANTS[AC_Kmgatp] / (1.0 + CONSTANTS[AC_MgATP] / CONSTANTS[AC_Kmgatp]);
    CONSTANTS[AC_b1] = CONSTANTS[AC_k1m] * CONSTANTS[AC_MgADP];
    
    /* IKb */
    CONSTANTS[AC_GKb] = 0.003;
    
    /* INab */
    CONSTANTS[AC_PNab] = 3.75e-10;
    
    /* ICab */
    CONSTANTS[AC_PCab] = 1.0 * 2.5e-08;
    
    /* IpCa */
    CONSTANTS[AC_GpCa] = 0.0575;
    
    /* SR_uptake */
    CONSTANTS[AC_BSLmax] = 1.124;
    CONSTANTS[AC_BSRmax] = 0.047;
    CONSTANTS[AC_KmBSL] = 0.0087;
    CONSTANTS[AC_KmBSR] = 0.00087;
    CONSTANTS[AC_cmdnmax] = 0.05;
    CONSTANTS[AC_csqnmax] = 10.0;
    CONSTANTS[AC_kmcmdn] = 0.00238;
    CONSTANTS[AC_kmcsqn] = 0.8;
    CONSTANTS[AC_kmtrpn] = 0.0005;
    CONSTANTS[AC_trpnmax] = 0.07;


   
    
    STATES[cell_v] = -87.3;
    STATES[ionic_concentrations_nai] = 6.43;
    STATES[ionic_concentrations_nass] = 6.43;
    STATES[ionic_concentrations_ki] = 140.76;
    STATES[ionic_concentrations_kss] = 140.76;
    STATES[ionic_concentrations_cai] = 6.75e-05;
    STATES[ionic_concentrations_cai2] = 6.81e-05;
    STATES[ionic_concentrations_cass] = 6.75e-05;
    STATES[ionic_concentrations_cansr] = 1.37;
    STATES[ionic_concentrations_cajsr] = 1.31;
    STATES[ionic_concentrations_cacsr] = 1.3;
    STATES[CaMK_CaMKt] =  8.23910999999999914e-03;
    STATES[CICR_Jrel2] = 8.59e-10;
    STATES[CICR_Jrel1] = 9.71e-22;
    STATES[ITo_aa] = 0.995;
    STATES[ICaL_d] = 4.42e-07;
    STATES[ICaL_fca] = 0.988;
    STATES[ICaL_ff] = 1.0;
    STATES[ICaL_fs] = 0.0;
    STATES[I_Na_m] = 0.0022;
    STATES[I_Na_h] = 0.631;
    STATES[I_Na_j] = 0.631;
    STATES[INaL_ml] = 0.001;
    STATES[INaL_hl] = 0.301;
    STATES[IKr_xr] = 0.153;
    STATES[IKs_xs1] = 0.054;
    STATES[IKs_xs2] = 0.046;

    STATES[ CICR_tjsrol] = 100.0; 
    STATES[CICR_A] = 100.0; 
}


void
GaurcomputeVariables(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC,int tissueFlag)
{

    double AV_Bcacsr;
    double AV_Bcai;
    double AV_Bcai2;
    double AV_Bcajsr;
    double AV_Bcass;
    double AV_CaMK_f;
    double AV_CaMKa;
    double AV_CaMKb;
    double AV_E1;
    double AV_E11;
    double AV_E12;
    double AV_E2;
    double AV_E21;
    double AV_E22;
    double AV_E3;
    double AV_E31;
    double AV_E32;
    double AV_E4;
    double AV_E41;
    double AV_E42;
    double AV_EK;
    double AV_EKs;
    double AV_ENa;
    double AV_Jdiff;
    double AV_JdiffK;
    double AV_JdiffNa;
    double AV_Jgap;
    double AV_Jleak;
    double AV_JnakK;
    double AV_JnakNa;
    double AV_JncxCa;
    double AV_JncxCa1;
    double AV_JncxNa;
    double AV_JncxNa1;
    double AV_Jrel;
    double AV_Jrel1_inf;
    double AV_Jrel2_inf;
    double AV_Jrelol;
    double AV_Jtr;
    double AV_Jtr2;
    double AV_Jup;
    double AV_Jup2;
    double AV_Jupnp;
    double AV_Jupnp2;
    double AV_Jupp;
    double AV_Jupp2;
    double AV_Knai;
    double AV_Knao;
    double AV_KsCa;
    double AV_P;
    double AV_PhiCaK;
    double AV_PhiCaL;
    double AV_PhiCaNa;
    double AV_Rel1;
    double AV_Rel2;
    double AV_a1;
    double AV_a3;
    double AV_aa_h;
    double AV_aa_j;
    double AV_aa_m;
    double AV_allo;
    double AV_allo1;
    double AV_alpha_aa;
    double AV_aml;
    double AV_b2;
    double AV_b3;
    double AV_b4;
    double AV_bb_h;
    double AV_bb_j;
    double AV_bb_m;
    double AV_beta_aa;
    double AV_bml;
    double AV_cai_avg;
    double AV_d_inf;
    double AV_diff_A;
    double AV_diff_CaMKt;
    double AV_diff_cacsr;
    double AV_diff_cai;
    double AV_diff_cai2;
    double AV_diff_cajsr;
    double AV_diff_cansr;
    double AV_diff_cass;
    double AV_diff_ki;
    double AV_diff_kss;
    double AV_diff_nai;
    double AV_diff_nass;
    double AV_diff_tjsrol;
    double AV_f;
    double AV_fJupp;
    double AV_f_inf;
    double AV_fca_inf;
    double AV_ff_inf;
    double AV_fs_inf;
    double AV_greljsrol;
    double AV_h1;
    double AV_h111;
    double AV_h2;
    double AV_h21;
    double AV_h3;
    double AV_h31;
    double AV_h4;
    double AV_h41;
    double AV_h5;
    double AV_h51;
    double AV_h6;
    double AV_h61;
    double AV_h7;
    double AV_h71;
    double AV_h8;
    double AV_h81;
    double AV_h9;
    double AV_h91;
    double AV_h_inf;
    double AV_hca;
    double AV_hl_inf;
    double AV_hna;
    double AV_j_inf;
    double AV_k3;
    double AV_k31;
    double AV_k3p;
    double AV_k3p1;
    double AV_k3pp;
    double AV_k3pp1;
    double AV_k4;
    double AV_k41;
    double AV_k4p;
    double AV_k4p1;
    double AV_k4pp;
    double AV_k4pp1;
    double AV_k6;
    double AV_k61;
    double AV_k7;
    double AV_k71;
    double AV_k8;
    double AV_k81;
    double AV_kito2;
    double AV_m_inf;
    double AV_ml_inf;
    double AV_rito2;
    double AV_rk1;
    double AV_rkr;
    double AV_tau_Jrel1;
    double AV_tau_Jrel2;
    double AV_tau_d;
    double AV_tau_fca;
    double AV_tau_ff;
    double AV_tau_fs;
    double AV_tau_h;
    double AV_tau_j;
    double AV_tau_m;
    double AV_tau_ml;
    double AV_tau_xr;
    double AV_tau_xs1;
    double AV_tau_xs2;
    double AV_trel1factor;
    double AV_trel2factor;
    double AV_vfrt;
    double AV_x1;
    double AV_x11;
    double AV_x12;
    double AV_x2;
    double AV_x21;
    double AV_x22;
    double AV_x3;
    double AV_x31;
    double AV_x32;
    double AV_x4;
    double AV_x41;
    double AV_x42;
    double AV_xkb;
    double AV_xr_inf;
    double AV_xs1_inf;
    double AV_xs2_inf;
  
 
    /* diffusion */
    AV_Jdiff = (STATES[ionic_concentrations_cass] - STATES[ionic_concentrations_cai]) / 0.2;
    AV_JdiffK = (STATES[ionic_concentrations_kss] - STATES[ionic_concentrations_ki]) / 2.0;
    AV_JdiffNa = (STATES[ionic_concentrations_nass] - STATES[ionic_concentrations_nai]) / 2.0;
    
    /* stimulus */
    ALGEBRAIC[AV_I_stim] = (
        (VOI - floor(VOI / CONSTANTS[AC_stimulus_period]) * CONSTANTS[AC_stimulus_period] >= CONSTANTS[AC_stimulus_start]) &&
        (VOI - floor(VOI / CONSTANTS[AC_stimulus_period]) * CONSTANTS[AC_stimulus_period] <= CONSTANTS[AC_stimulus_start] + CONSTANTS[AC_stimulus_duration])
        ? CONSTANTS[AC_stimulus_amplitude]
        : 0.0
    );
    
    
    /* CaMK */
    AV_EK = CONSTANTS[AC_R] * CONSTANTS[AC_T] / CONSTANTS[AC_F] * log(CONSTANTS[AC_ko] / STATES[ionic_concentrations_ki]);
    AV_ENa = CONSTANTS[AC_R] * CONSTANTS[AC_T] / CONSTANTS[AC_F] * log(CONSTANTS[AC_nao] / STATES[ionic_concentrations_nai]);
    AV_vfrt = STATES[cell_v] * CONSTANTS[AC_F] / (CONSTANTS[AC_R] * CONSTANTS[AC_T]);
    AV_CaMKb = CONSTANTS[AC_CaMKo] * (1.0 - STATES[CaMK_CaMKt]) / (1.0 + CONSTANTS[AC_KmCaM] / STATES[ionic_concentrations_cass]);
    AV_EKs = CONSTANTS[AC_R] * CONSTANTS[AC_T] / CONSTANTS[AC_F] * log((CONSTANTS[AC_ko] + CONSTANTS[AC_PKNa] * CONSTANTS[AC_nao]) / (STATES[ionic_concentrations_ki] + CONSTANTS[AC_PKNa] * STATES[ionic_concentrations_nai]));
    AV_CaMKa = AV_CaMKb + STATES[CaMK_CaMKt];
    AV_diff_CaMKt = CONSTANTS[AC_aCaMK] * AV_CaMKb * (AV_CaMKb + STATES[CaMK_CaMKt]) - CONSTANTS[AC_bCaMK] * STATES[CaMK_CaMKt];
    AV_CaMK_f = 1.0 / (1.0 + CONSTANTS[AC_KmCaMK] / AV_CaMKa);
    RATES[CaMK_CaMKt] = AV_diff_CaMKt;
    
    /* I_Na */
    AV_aa_h = ((STATES[cell_v] >= -40.0) ? 0.0 : 0.057 * exp(-(STATES[cell_v] + 80.0) / 6.8));
    AV_aa_j = ((STATES[cell_v] >= -40.0) ? 0.0 : (-25428.0 * exp(0.2444 * STATES[cell_v]) - 6.948e-06 * exp(-0.04391 * STATES[cell_v])) * (STATES[cell_v] + 37.78) / (1.0 + exp(0.311 * (STATES[cell_v] + 79.23))));
    AV_aa_m = 1.0 / (1.0 + exp((-60.0 - STATES[cell_v]) / 5.0));
    AV_bb_h = ((STATES[cell_v] >= -40.0) ? 0.77 / (0.13 * (1.0 + exp(-(STATES[cell_v] + 10.66) / 11.1))) : 2.7 * exp(0.079 * STATES[cell_v]) + 310000.0 * exp(0.3485 * STATES[cell_v]));
    AV_bb_j = ((STATES[cell_v] >= -40.0) ? 0.6 * exp(0.057 * STATES[cell_v]) / (1.0 + exp(-0.1 * (STATES[cell_v] + 32.0))) : 0.02424 * exp(-0.01052 * STATES[cell_v]) / (1.0 + exp(-0.1378 * (STATES[cell_v] + 40.14))));
    AV_bb_m = 0.1 / (1.0 + exp((STATES[cell_v] + 35.0) / 5.0)) + 0.1 / (1.0 + exp((STATES[cell_v] - 50.0) / 200.0));
    AV_h_inf = 1.0 / ((1.0 + exp((STATES[cell_v] + 71.55) / 7.43)) * (1.0 + exp((STATES[cell_v] + 71.55) / 7.43)));
    AV_m_inf = 1.0 / ((1.0 + exp((-56.86 - STATES[cell_v]) / 9.03)) * (1.0 + exp((-56.86 - STATES[cell_v]) / 9.03)));
    ALGEBRAIC[AV_INa] = CONSTANTS[AC_GNa] * STATES[I_Na_m] * STATES[I_Na_m] * STATES[I_Na_m] * STATES[I_Na_h] * STATES[I_Na_j] * (STATES[cell_v] - AV_ENa);
    AV_j_inf = AV_h_inf;
    AV_tau_h = 1.0 / (AV_aa_h + AV_bb_h);
    AV_tau_j = 1.0 / (AV_aa_j + AV_bb_j);
    AV_tau_m = AV_aa_m * AV_bb_m;
    RATES[I_Na_h] = (AV_h_inf - STATES[I_Na_h]) / AV_tau_h;
    RATES[I_Na_j] = (AV_j_inf - STATES[I_Na_j]) / AV_tau_j;
    RATES[I_Na_m] = (AV_m_inf - STATES[I_Na_m]) / AV_tau_m;
    
    /* INaL */
    AV_aml = 0.32 * ((STATES[cell_v] == -47.13) ? 10.0 : -(STATES[cell_v] + 47.13) / (exp(-0.1 * (STATES[cell_v] + 47.13)) - 1.0));
    AV_bml = 0.08 * exp(-STATES[cell_v] / 11.0);
    AV_hl_inf = 1.0 / (1.0 + exp((STATES[cell_v] + 91.0) / 6.1));
    ALGEBRAIC[AV_INaL_INaL] = CONSTANTS[AC_GNaL] * STATES[INaL_ml] * STATES[INaL_ml] * STATES[INaL_ml] * STATES[INaL_hl] * (STATES[cell_v] - AV_ENa);
    AV_ml_inf = AV_aml / (AV_aml + AV_bml);
    AV_tau_ml = 1.0 / (AV_aml + AV_bml);
    RATES[INaL_hl] = (AV_hl_inf - STATES[INaL_hl]) / CONSTANTS[AC_tau_hl];
    RATES[INaL_ml] = (AV_ml_inf - STATES[INaL_ml]) / AV_tau_ml;
    
    /* ICaL */
    AV_PhiCaK = CONSTANTS[AC_F] * (0.75 * STATES[ionic_concentrations_kss] * exp(AV_vfrt) - 0.75 * CONSTANTS[AC_ko]) * ((STATES[cell_v] == 0.0) ? 1.0 : AV_vfrt / (exp(AV_vfrt) - 1.0));
    AV_PhiCaL = 2.0 * CONSTANTS[AC_F] * (STATES[ionic_concentrations_cass] * exp(2.0 * (STATES[cell_v] * CONSTANTS[AC_F] / (CONSTANTS[AC_R] * CONSTANTS[AC_T]))) - 0.341 * CONSTANTS[AC_cao]) * ((STATES[cell_v] == 0.0) ? 1.0 : 2.0 * AV_vfrt / (exp(2.0 * AV_vfrt) - 1.0));
    AV_PhiCaNa = CONSTANTS[AC_F] * (0.75 * STATES[ionic_concentrations_nass] * exp(AV_vfrt) - 0.75 * CONSTANTS[AC_nao]) * ((STATES[cell_v] == 0.0) ? 1.0 : AV_vfrt / (exp(AV_vfrt) - 1.0));
    AV_f = STATES[ICaL_ff] * STATES[ICaL_fs];
    AV_tau_d = 0.6 + 1.0 / (exp(-0.05 * (STATES[cell_v] + 6.0)) + exp(0.09 * (STATES[cell_v] + 14.0)));
    AV_tau_fca = 10.0 * AV_CaMK_f + 0.5 + 1.0 / (1.0 + STATES[ionic_concentrations_cass] / 0.003);
    AV_tau_ff = 7.0 + 1.0 / (0.0045 * exp(-(STATES[cell_v] + 20.0) / 10.0) + 0.0045 * exp((STATES[cell_v] + 20.0) / 10.0));
    AV_tau_fs = 70.0 + 1.0 / (3.5e-05 * exp(-(STATES[cell_v] + 5.0) / 4.0) + 3.5e-05 * exp((STATES[cell_v] + 5.0) / 6.0));
    ALGEBRAIC[AV_ICaL_ICaL] = CONSTANTS[AC_PCa] * AV_PhiCaL * STATES[ICaL_d] * AV_f * STATES[ICaL_fca];
    AV_d_inf = 1.0 / (1.0 + exp(-(STATES[cell_v] - CONSTANTS[AC_vhalf_d]) / 6.2));
    AV_f_inf = 1.0 / (1.0 + exp((STATES[cell_v] - CONSTANTS[AC_vhalff]) / 4.9)) + 0.35 / (1.0 + exp((45.0 - STATES[cell_v]) / 20.0));
    ALGEBRAIC[AV_ICaK] = CONSTANTS[AC_PCaK] * AV_PhiCaK * STATES[ICaL_d] * AV_f * STATES[ICaL_fca];
    ALGEBRAIC[AV_ICaNa] = CONSTANTS[AC_PCaNa] * AV_PhiCaNa * STATES[ICaL_d] * AV_f * STATES[ICaL_fca];
    AV_fca_inf = 0.3 / (1.0 - ((ALGEBRAIC[AV_ICaL_ICaL] > 0.0) ? 0.0 : ALGEBRAIC[AV_ICaL_ICaL] / 0.05)) + 0.55 / (1.0 + STATES[ionic_concentrations_cass] / 0.003) + 0.15;
    AV_ff_inf = AV_f_inf;
    AV_fs_inf = AV_f_inf;
    RATES[ICaL_d] = (AV_d_inf - STATES[ICaL_d]) / AV_tau_d;
    RATES[ICaL_fca] = (AV_fca_inf - STATES[ICaL_fca]) / AV_tau_fca;
    RATES[ICaL_ff] = (AV_ff_inf - STATES[ICaL_ff]) / AV_tau_ff;
    RATES[ICaL_fs] = (AV_fs_inf - STATES[ICaL_fs]) / AV_tau_fs;
    
    /* IKr */
    AV_rkr = 1.0 / (1.0 + exp((STATES[cell_v] + 22.0) / 15.0));
    AV_tau_xr = 12.98 + 1.0 / (0.3652 * exp((STATES[cell_v] - 31.66) / 3.869) + 4.123e-05 * exp(-(STATES[cell_v] - 47.78) / 20.38));
    AV_xr_inf = 1.0 / (1.0 + exp(-(STATES[cell_v] + 56.8) / 17.8));
    ALGEBRAIC[AV_IKr_IKr] = CONSTANTS[AC_GKr] * sqrt(CONSTANTS[AC_ko] / 5.4) * STATES[IKr_xr] * AV_rkr * (STATES[cell_v] - AV_EK);
    RATES[IKr_xr] = (AV_xr_inf - STATES[IKr_xr]) / AV_tau_xr;
    
    /* IKs */
    AV_KsCa = 1.0 + 0.6 / (1.0 + pow(3.8e-05 / STATES[ionic_concentrations_cai], 1.4));
    AV_tau_xs1 = 300.0 + 1.0 / (1e-06 * exp((STATES[cell_v] + 50.0) / 20.0) + 0.04 * exp(-(STATES[cell_v] + 50.0) / 20.0));
    AV_tau_xs2 = 1.0 / (0.01 * exp((STATES[cell_v] - 50.0) / 20.0) + 0.0193 * exp(-(STATES[cell_v] + 66.54) / 31.0));
    AV_xs1_inf = 1.0 / (1.0 + exp(-(STATES[cell_v] - 25.1) / 37.1));
    ALGEBRAIC[AV_IKs_IKs] = CONSTANTS[AC_GKs] * AV_KsCa * STATES[IKs_xs1] * STATES[IKs_xs2] * (STATES[cell_v] - AV_EKs);
    AV_xs2_inf = AV_xs1_inf;
    RATES[IKs_xs1] = (AV_xs1_inf - STATES[IKs_xs1]) / AV_tau_xs1;
    RATES[IKs_xs2] = (AV_xs2_inf - STATES[IKs_xs2]) / AV_tau_xs2;
    
    /* IK1 */
    AV_rk1 = 1.0 / (1.0 + exp((STATES[cell_v] + 79.3 - 2.6 * CONSTANTS[AC_ko]) / 19.6));
    ALGEBRAIC[AV_IK1_IK1] = CONSTANTS[AC_GK1] * sqrt(CONSTANTS[AC_ko] / 5.4) * AV_rk1 * (STATES[cell_v] - AV_EK);
    
    /* ITo */
    AV_alpha_aa = 0.025 / (1.0 + exp((STATES[cell_v] + 58.0) / 5.0));
    AV_beta_aa = 1.0 / (5.0 * (1.0 + exp((STATES[cell_v] + 19.0) / -9.0)));
    AV_rito2 = 1.0 / (1.0 + exp(-(STATES[cell_v] + 10.0) / 5.0));
    RATES[ITo_aa] = AV_alpha_aa * (1.0 - STATES[ITo_aa]) - AV_beta_aa * STATES[ITo_aa];
    
    /* CICR */
    AV_diff_A = ((STATES[CICR_tjsrol] < 5.0) ? -STATES[CICR_A] / 0.1 : 1.0);
    AV_trel1factor = 1.0 * 20.0 * (1.0 + 1.0 / (1.0 + pow(CONSTANTS[AC_KmCaMK] / AV_CaMKa, 8.0))) / (1.0 + pow(0.5 / STATES[ionic_concentrations_cajsr], 8.0));
    AV_trel2factor = 50.0 * (1.0 + 1.0 / (1.0 + pow(CONSTANTS[AC_KmCaMK] / AV_CaMKa, 8.0))) / (1.0 + pow(0.5 / STATES[ionic_concentrations_cacsr], 8.0));
    AV_Jgap = (STATES[ionic_concentrations_cai] - STATES[ionic_concentrations_cai2]) / CONSTANTS[AC_tau_gap];
    AV_diff_tjsrol = (((STATES[ionic_concentrations_cajsr] > CONSTANTS[AC_SOICR]) && (STATES[CICR_A] > 45.0)) ? -STATES[CICR_tjsrol] / 0.001 : 1.0);
    AV_greljsrol = CONSTANTS[AC_grelbarjsrol] * (1.0 - exp(-STATES[CICR_tjsrol] / CONSTANTS[AC_tauon])) * exp(-STATES[CICR_tjsrol] / CONSTANTS[AC_tauoff]);
    AV_tau_Jrel1 = ((0.001 > AV_trel1factor) ? 0.001 : AV_trel1factor);
    AV_tau_Jrel2 = ((0.001 > AV_trel2factor) ? 0.001 : AV_trel2factor);
    RATES[CICR_A] = AV_diff_A;
    AV_Jrelol = AV_greljsrol * (STATES[ionic_concentrations_cajsr] - STATES[ionic_concentrations_cass]);
    RATES[CICR_tjsrol] = AV_diff_tjsrol;
    AV_Jrel = STATES[CICR_Jrel1] + AV_Jrelol;
    
    /* INaCa_i */
    AV_allo = 1.0 / (1.0 + CONSTANTS[AC_KmCaAct] / STATES[ionic_concentrations_cai] * CONSTANTS[AC_KmCaAct] / STATES[ionic_concentrations_cai]);
    AV_h4 = 1.0 + STATES[ionic_concentrations_nai] / CONSTANTS[AC_kna1] * (1.0 + STATES[ionic_concentrations_nai] / CONSTANTS[AC_kna2]);
    AV_hca = exp(CONSTANTS[AC_qca] * STATES[cell_v] * CONSTANTS[AC_F] / (CONSTANTS[AC_R] * CONSTANTS[AC_T]));
    AV_hna = exp(CONSTANTS[AC_qna] * STATES[cell_v] * CONSTANTS[AC_F] / (CONSTANTS[AC_R] * CONSTANTS[AC_T]));
    AV_h1 = 1.0 + STATES[ionic_concentrations_nai] / CONSTANTS[AC_kna3] * (1.0 + AV_hna);
    AV_h5 = STATES[ionic_concentrations_nai] * STATES[ionic_concentrations_nai] / (AV_h4 * CONSTANTS[AC_kna1] * CONSTANTS[AC_kna2]);
    AV_h6 = 1.0 / AV_h4;
    AV_h7 = 1.0 + CONSTANTS[AC_nao] / CONSTANTS[AC_kna3] * (1.0 + 1.0 / AV_hna);
    AV_h2 = STATES[ionic_concentrations_nai] * AV_hna / (CONSTANTS[AC_kna3] * AV_h1);
    AV_h3 = 1.0 / AV_h1;
    AV_h8 = CONSTANTS[AC_nao] / (CONSTANTS[AC_kna3] * AV_hna * AV_h7);
    AV_h9 = 1.0 / AV_h7;
    AV_k6 = AV_h6 * STATES[ionic_concentrations_cai] * CONSTANTS[AC_kcaon];
    AV_k3p = AV_h9 * CONSTANTS[AC_wca];
    AV_k3pp = AV_h8 * CONSTANTS[AC_wnaca];
    AV_k4p = AV_h3 * CONSTANTS[AC_wca] / AV_hca;
    AV_k4pp = AV_h2 * CONSTANTS[AC_wnaca];
    AV_k7 = AV_h5 * AV_h2 * CONSTANTS[AC_wna];
    AV_k8 = AV_h8 * CONSTANTS[AC_h11] * CONSTANTS[AC_wna];
    AV_k3 = AV_k3p + AV_k3pp;
    AV_k4 = AV_k4p + AV_k4pp;
    AV_x1 = CONSTANTS[AC_k2] * AV_k4 * (AV_k7 + AV_k6) + CONSTANTS[AC_k5] * AV_k7 * (CONSTANTS[AC_k2] + AV_k3);
    AV_x2 = CONSTANTS[AC_k1] * AV_k7 * (AV_k4 + CONSTANTS[AC_k5]) + AV_k4 * AV_k6 * (CONSTANTS[AC_k1] + AV_k8);
    AV_x3 = CONSTANTS[AC_k1] * AV_k3 * (AV_k7 + AV_k6) + AV_k8 * AV_k6 * (CONSTANTS[AC_k2] + AV_k3);
    AV_x4 = CONSTANTS[AC_k2] * AV_k8 * (AV_k4 + CONSTANTS[AC_k5]) + AV_k3 * CONSTANTS[AC_k5] * (CONSTANTS[AC_k1] + AV_k8);
    AV_E1 = AV_x1 / (AV_x1 + AV_x2 + AV_x3 + AV_x4);
    AV_E2 = AV_x2 / (AV_x1 + AV_x2 + AV_x3 + AV_x4);
    AV_E3 = AV_x3 / (AV_x1 + AV_x2 + AV_x3 + AV_x4);
    AV_E4 = AV_x4 / (AV_x1 + AV_x2 + AV_x3 + AV_x4);
    AV_JncxCa = AV_E2 * CONSTANTS[AC_k2] - AV_E1 * CONSTANTS[AC_k1];
    AV_JncxNa = 3.0 * (AV_E4 * AV_k7 - AV_E1 * AV_k8) + AV_E3 * AV_k4pp - AV_E2 * AV_k3pp;
    ALGEBRAIC[AV_INaCa_i_INaCa_i] = 0.8 * CONSTANTS[AC_Gncx] * AV_allo * (CONSTANTS[AC_zna] * AV_JncxNa + CONSTANTS[AC_zca] * AV_JncxCa);
    
    /* INaCa_ss */
    AV_allo1 = 1.0 / (1.0 + CONSTANTS[AC_KmCaAct] / STATES[ionic_concentrations_cass] * CONSTANTS[AC_KmCaAct] / STATES[ionic_concentrations_cass]);
    AV_h111 = 1.0 + STATES[ionic_concentrations_nass] / CONSTANTS[AC_kna3] * (1.0 + AV_hna);
    AV_h41 = 1.0 + STATES[ionic_concentrations_nass] / CONSTANTS[AC_kna1] * (1.0 + STATES[ionic_concentrations_nass] / CONSTANTS[AC_kna2]);
    AV_h71 = 1.0 + CONSTANTS[AC_nao] / CONSTANTS[AC_kna3] * (1.0 + 1.0 / AV_hna);
    AV_h21 = STATES[ionic_concentrations_nass] * AV_hna / (CONSTANTS[AC_kna3] * AV_h111);
    AV_h31 = 1.0 / AV_h111;
    AV_h51 = STATES[ionic_concentrations_nass] * STATES[ionic_concentrations_nass] / (AV_h41 * CONSTANTS[AC_kna1] * CONSTANTS[AC_kna2]);
    AV_h61 = 1.0 / AV_h41;
    AV_h81 = CONSTANTS[AC_nao] / (CONSTANTS[AC_kna3] * AV_hna * AV_h71);
    AV_h91 = 1.0 / AV_h71;
    AV_k3p1 = AV_h91 * CONSTANTS[AC_wca];
    AV_k3pp1 = AV_h81 * CONSTANTS[AC_wnaca];
    AV_k4p1 = AV_h31 * CONSTANTS[AC_wca] / AV_hca;
    AV_k4pp1 = AV_h21 * CONSTANTS[AC_wnaca];
    AV_k61 = AV_h61 * STATES[ionic_concentrations_cass] * CONSTANTS[AC_kcaon];
    AV_k71 = AV_h51 * AV_h21 * CONSTANTS[AC_wna];
    AV_k81 = AV_h81 * CONSTANTS[AC_h1111] * CONSTANTS[AC_wna];
    AV_k31 = AV_k3p1 + AV_k3pp1;
    AV_k41 = AV_k4p1 + AV_k4pp1;
    AV_x11 = CONSTANTS[AC_k21] * AV_k41 * (AV_k71 + AV_k61) + CONSTANTS[AC_k51] * AV_k71 * (CONSTANTS[AC_k21] + AV_k31);
    AV_x21 = CONSTANTS[AC_k11] * AV_k71 * (AV_k41 + CONSTANTS[AC_k51]) + AV_k41 * AV_k61 * (CONSTANTS[AC_k11] + AV_k81);
    AV_x31 = CONSTANTS[AC_k11] * AV_k31 * (AV_k71 + AV_k61) + AV_k81 * AV_k61 * (CONSTANTS[AC_k21] + AV_k31);
    AV_x41 = CONSTANTS[AC_k21] * AV_k81 * (AV_k41 + CONSTANTS[AC_k51]) + AV_k31 * CONSTANTS[AC_k51] * (CONSTANTS[AC_k11] + AV_k81);
    AV_E11 = AV_x11 / (AV_x11 + AV_x21 + AV_x31 + AV_x41);
    AV_E21 = AV_x21 / (AV_x11 + AV_x21 + AV_x31 + AV_x41);
    AV_E31 = AV_x31 / (AV_x11 + AV_x21 + AV_x31 + AV_x41);
    AV_E41 = AV_x41 / (AV_x11 + AV_x21 + AV_x31 + AV_x41);
    AV_JncxCa1 = AV_E21 * CONSTANTS[AC_k21] - AV_E11 * CONSTANTS[AC_k11];
    AV_JncxNa1 = 3.0 * (AV_E41 * AV_k71 - AV_E11 * AV_k81) + AV_E31 * AV_k4pp1 - AV_E21 * AV_k3pp1;
    ALGEBRAIC[AV_INaCa_ss_INaCa_ss] = 0.2 * CONSTANTS[AC_Gncx] * AV_allo1 * (CONSTANTS[AC_zna] * AV_JncxNa1 + CONSTANTS[AC_zca] * AV_JncxCa1);
    ALGEBRAIC[AV_INaCa] = ALGEBRAIC[AV_INaCa_i_INaCa_i] + ALGEBRAIC[AV_INaCa_ss_INaCa_ss];
    
    /* INaK */
    AV_Knai = CONSTANTS[AC_Knai0] * exp(CONSTANTS[AC_delta] * STATES[cell_v] * CONSTANTS[AC_F] / (3.0 * CONSTANTS[AC_R] * CONSTANTS[AC_T]));
    AV_Knao = CONSTANTS[AC_Knao0] * exp((1.0 - CONSTANTS[AC_delta]) * STATES[cell_v] * CONSTANTS[AC_F] / (3.0 * CONSTANTS[AC_R] * CONSTANTS[AC_T]));
    AV_P = CONSTANTS[AC_eP] / (1.0 + CONSTANTS[AC_H] / CONSTANTS[AC_Khp] + STATES[ionic_concentrations_nai] / CONSTANTS[AC_Knap] + STATES[ionic_concentrations_ki] / CONSTANTS[AC_Kxkur]);
    AV_a1 = CONSTANTS[AC_k1p] * (STATES[ionic_concentrations_nai] / AV_Knai * STATES[ionic_concentrations_nai] / AV_Knai * (STATES[ionic_concentrations_nai] / AV_Knai)) / ((1.0 + STATES[ionic_concentrations_nai] / AV_Knai) * (1.0 + STATES[ionic_concentrations_nai] / AV_Knai) * (1.0 + STATES[ionic_concentrations_nai] / AV_Knai) + (1.0 + STATES[ionic_concentrations_ki] / CONSTANTS[AC_Kki]) * (1.0 + STATES[ionic_concentrations_ki] / CONSTANTS[AC_Kki]) - 1.0);
    AV_a3 = CONSTANTS[AC_k3p2] * (CONSTANTS[AC_ko] / CONSTANTS[AC_Kko] * CONSTANTS[AC_ko] / CONSTANTS[AC_Kko]) / ((1.0 + CONSTANTS[AC_nao] / AV_Knao) * (1.0 + CONSTANTS[AC_nao] / AV_Knao) * (1.0 + CONSTANTS[AC_nao] / AV_Knao) + (1.0 + CONSTANTS[AC_ko] / CONSTANTS[AC_Kko]) * (1.0 + CONSTANTS[AC_ko] / CONSTANTS[AC_Kko]) - 1.0);
    AV_b2 = CONSTANTS[AC_k2m] * (CONSTANTS[AC_nao] / AV_Knao * CONSTANTS[AC_nao] / AV_Knao * (CONSTANTS[AC_nao] / AV_Knao)) / ((1.0 + CONSTANTS[AC_nao] / AV_Knao) * (1.0 + CONSTANTS[AC_nao] / AV_Knao) * (1.0 + CONSTANTS[AC_nao] / AV_Knao) + (1.0 + CONSTANTS[AC_ko] / CONSTANTS[AC_Kko]) * (1.0 + CONSTANTS[AC_ko] / CONSTANTS[AC_Kko]) - 1.0);
    AV_b3 = CONSTANTS[AC_k3m] * AV_P * CONSTANTS[AC_H] / (1.0 + CONSTANTS[AC_MgATP] / CONSTANTS[AC_Kmgatp]);
    AV_b4 = CONSTANTS[AC_k4m] * (STATES[ionic_concentrations_ki] / CONSTANTS[AC_Kki] * STATES[ionic_concentrations_ki] / CONSTANTS[AC_Kki]) / ((1.0 + STATES[ionic_concentrations_nai] / AV_Knai) * (1.0 + STATES[ionic_concentrations_nai] / AV_Knai) * (1.0 + STATES[ionic_concentrations_nai] / AV_Knai) + (1.0 + STATES[ionic_concentrations_ki] / CONSTANTS[AC_Kki]) * (1.0 + STATES[ionic_concentrations_ki] / CONSTANTS[AC_Kki]) - 1.0);
    AV_x12 = CONSTANTS[AC_a4] * AV_a1 * CONSTANTS[AC_a2] + AV_b2 * AV_b4 * AV_b3 + CONSTANTS[AC_a2] * AV_b4 * AV_b3 + AV_b3 * AV_a1 * CONSTANTS[AC_a2];
    AV_x22 = AV_b2 * CONSTANTS[AC_b1] * AV_b4 + AV_a1 * CONSTANTS[AC_a2] * AV_a3 + AV_a3 * CONSTANTS[AC_b1] * AV_b4 + CONSTANTS[AC_a2] * AV_a3 * AV_b4;
    AV_x32 = CONSTANTS[AC_a2] * AV_a3 * CONSTANTS[AC_a4] + AV_b3 * AV_b2 * CONSTANTS[AC_b1] + AV_b2 * CONSTANTS[AC_b1] * CONSTANTS[AC_a4] + AV_a3 * CONSTANTS[AC_a4] * CONSTANTS[AC_b1];
    AV_x42 = AV_b4 * AV_b3 * AV_b2 + AV_a3 * CONSTANTS[AC_a4] * AV_a1 + AV_b2 * CONSTANTS[AC_a4] * AV_a1 + AV_b3 * AV_b2 * AV_a1;
    AV_E12 = AV_x12 / (AV_x12 + AV_x22 + AV_x32 + AV_x42);
    AV_E22 = AV_x22 / (AV_x12 + AV_x22 + AV_x32 + AV_x42);
    AV_E32 = AV_x32 / (AV_x12 + AV_x22 + AV_x32 + AV_x42);
    AV_E42 = AV_x42 / (AV_x12 + AV_x22 + AV_x32 + AV_x42);
    AV_JnakK = 2.0 * (AV_E42 * CONSTANTS[AC_b1] - AV_E32 * AV_a1);
    AV_JnakNa = 3.0 * (AV_E12 * AV_a3 - AV_E22 * AV_b3);
    ALGEBRAIC[AV_INaK_INaK] = CONSTANTS[AC_Pnak] * (CONSTANTS[AC_zna] * AV_JnakNa + CONSTANTS[AC_zk] * AV_JnakK);
    
    /* IKb */
    AV_xkb = 1.0 / (1.0 + exp(-(STATES[cell_v] - 14.48) / 18.34));
    ALGEBRAIC[AV_IKb_IKb] = CONSTANTS[AC_GKb] * AV_xkb * (STATES[cell_v] - AV_EK);
    
    /* INab */
    ALGEBRAIC[AV_INab_INab] = CONSTANTS[AC_PNab] * CONSTANTS[AC_F] * (STATES[ionic_concentrations_nai] * exp(AV_vfrt) - CONSTANTS[AC_nao]) * ((STATES[cell_v] == 0.0) ? 1.0 : AV_vfrt / (exp(AV_vfrt) - 1.0));
    
    /* ICab */
    ALGEBRAIC[AV_ICab_ICab] = CONSTANTS[AC_PCab] * 2.0 * CONSTANTS[AC_F] * (STATES[ionic_concentrations_cai] * exp(2.0 * AV_vfrt) - 0.341 * CONSTANTS[AC_cao]) * ((STATES[cell_v] == 0.0) ? 1.0 : 2.0 * AV_vfrt / (exp(2.0 * AV_vfrt) - 1.0));
    
    /* IpCa */
    ALGEBRAIC[AV_IpCa_IpCa] = CONSTANTS[AC_GpCa] * STATES[ionic_concentrations_cai] / (0.0005 + STATES[ionic_concentrations_cai]);
    
    /* SR_uptake */
    AV_Jleak = 0.005 * STATES[ionic_concentrations_cansr] / 15.0;
    AV_Jtr = (STATES[ionic_concentrations_cansr] - STATES[ionic_concentrations_cajsr]) / 100.0;
    AV_Jtr2 = (STATES[ionic_concentrations_cansr] - STATES[ionic_concentrations_cacsr]) / 100.0;
    AV_Jupnp = 1.0 * 0.005 * STATES[ionic_concentrations_cai] / (STATES[ionic_concentrations_cai] + 0.001);
    AV_Jupnp2 = 1.0 * 0.005 * STATES[ionic_concentrations_cai2] / (STATES[ionic_concentrations_cai2] + 0.001);
    AV_Jupp = 1.0 * 3.0 * 0.005 * STATES[ionic_concentrations_cai] / (STATES[ionic_concentrations_cai] + 0.001 - 0.0002);
    AV_Jupp2 = 1.0 * 3.0 * 0.005 * STATES[ionic_concentrations_cai2] / (STATES[ionic_concentrations_cai2] + 0.001 - 0.0002);
    AV_fJupp = 1.0 / (1.0 + CONSTANTS[AC_KmCaMK] / AV_CaMKa);
    AV_Jup = (1.0 - AV_fJupp) * AV_Jupnp + AV_fJupp * AV_Jupp - AV_Jleak;
    AV_Jup2 = (1.0 - AV_fJupp) * AV_Jupnp2 + AV_fJupp * AV_Jupp2 - AV_Jleak;
    
    /* ionic_concentrations */
    AV_Bcacsr = 1.0 / (1.0 + CONSTANTS[AC_csqnmax] * CONSTANTS[AC_kmcsqn] / ((CONSTANTS[AC_kmcsqn] + STATES[ionic_concentrations_cacsr]) * (CONSTANTS[AC_kmcsqn] + STATES[ionic_concentrations_cacsr])));
    AV_Bcai = 1.0 / (1.0 + CONSTANTS[AC_cmdnmax] * CONSTANTS[AC_kmcmdn] / ((CONSTANTS[AC_kmcmdn] + STATES[ionic_concentrations_cai]) * (CONSTANTS[AC_kmcmdn] + STATES[ionic_concentrations_cai])) + CONSTANTS[AC_trpnmax] * CONSTANTS[AC_kmtrpn] / ((CONSTANTS[AC_kmtrpn] + STATES[ionic_concentrations_cai]) * (CONSTANTS[AC_kmtrpn] + STATES[ionic_concentrations_cai])));
    AV_Bcai2 = 1.0 / (1.0 + CONSTANTS[AC_cmdnmax] * CONSTANTS[AC_kmcmdn] / ((CONSTANTS[AC_kmcmdn] + STATES[ionic_concentrations_cai2]) * (CONSTANTS[AC_kmcmdn] + STATES[ionic_concentrations_cai2])) + CONSTANTS[AC_trpnmax] * CONSTANTS[AC_kmtrpn] / ((CONSTANTS[AC_kmtrpn] + STATES[ionic_concentrations_cai2]) * (CONSTANTS[AC_kmtrpn] + STATES[ionic_concentrations_cai2])));
    AV_Bcajsr = 1.0 / (1.0 + CONSTANTS[AC_csqnmax] * CONSTANTS[AC_kmcsqn] / ((CONSTANTS[AC_kmcsqn] + STATES[ionic_concentrations_cajsr]) * (CONSTANTS[AC_kmcsqn] + STATES[ionic_concentrations_cajsr])));
    AV_Bcass = 1.0 / (1.0 + CONSTANTS[AC_BSRmax] * CONSTANTS[AC_KmBSR] / ((CONSTANTS[AC_KmBSR] + STATES[ionic_concentrations_cass]) * (CONSTANTS[AC_KmBSR] + STATES[ionic_concentrations_cass])) + CONSTANTS[AC_BSLmax] * CONSTANTS[AC_KmBSL] / ((CONSTANTS[AC_KmBSL] + STATES[ionic_concentrations_cass]) * (CONSTANTS[AC_KmBSL] + STATES[ionic_concentrations_cass])));
    AV_cai_avg = CONSTANTS[AC_vmyo1frac] * STATES[ionic_concentrations_cai] + (1.0 - CONSTANTS[AC_vmyo1frac]) * STATES[ionic_concentrations_cai2];
    AV_diff_cansr = AV_Jup * (CONSTANTS[AC_vnsr1] / CONSTANTS[AC_vnsr]) + AV_Jup2 * (CONSTANTS[AC_vnsr2] / CONSTANTS[AC_vnsr]) - AV_Jtr * CONSTANTS[AC_vjsr] / CONSTANTS[AC_vnsr] - AV_Jtr2 * CONSTANTS[AC_vcsr] / CONSTANTS[AC_vnsr];
    AV_diff_ki = -(ALGEBRAIC[AV_IKr_IKr] + ALGEBRAIC[AV_IKs_IKs] + ALGEBRAIC[AV_IK1_IK1] + ALGEBRAIC[AV_IKb_IKb] - 2.0 * ALGEBRAIC[AV_INaK_INaK]) * CONSTANTS[AC_Acap] / (2.0 * CONSTANTS[AC_F] * CONSTANTS[AC_vmyo]) + AV_JdiffK * CONSTANTS[AC_vss] / CONSTANTS[AC_vmyo];
    AV_diff_kss = -ALGEBRAIC[AV_ICaK] * CONSTANTS[AC_Acap] / (2.0 * CONSTANTS[AC_F] * CONSTANTS[AC_vss]) - AV_JdiffK;
    AV_diff_nai = -(ALGEBRAIC[AV_INa] + ALGEBRAIC[AV_INaL_INaL] + 3.0 * ALGEBRAIC[AV_INaCa_i_INaCa_i] + 3.0 * ALGEBRAIC[AV_INaK_INaK] + ALGEBRAIC[AV_INab_INab]) * CONSTANTS[AC_Acap] / (2.0 * CONSTANTS[AC_F] * CONSTANTS[AC_vmyo]) + AV_JdiffNa * CONSTANTS[AC_vss] / CONSTANTS[AC_vmyo];
    AV_diff_nass = -(ALGEBRAIC[AV_ICaNa] + 3.0 * ALGEBRAIC[AV_INaCa_ss_INaCa_ss]) * CONSTANTS[AC_Acap] / (2.0 * CONSTANTS[AC_F] * CONSTANTS[AC_vss]) - AV_JdiffNa;
    RATES[ionic_concentrations_cansr] = AV_diff_cansr;
    RATES[ionic_concentrations_ki] = AV_diff_ki;
    RATES[ionic_concentrations_kss] = AV_diff_kss;
    RATES[ionic_concentrations_nai] = AV_diff_nai;
    RATES[ionic_concentrations_nass] = AV_diff_nass;
    AV_diff_cacsr = AV_Bcacsr * (AV_Jtr2 - STATES[CICR_Jrel2]);
    AV_diff_cai = AV_Bcai * (-(ALGEBRAIC[AV_IpCa_IpCa] + ALGEBRAIC[AV_ICab_ICab] - 2.0 * ALGEBRAIC[AV_INaCa_i_INaCa_i]) * CONSTANTS[AC_Acap] / (2.0 * 2.0 * CONSTANTS[AC_F] * CONSTANTS[AC_vmyo1]) - AV_Jup * CONSTANTS[AC_vnsr1] / CONSTANTS[AC_vmyo1] + AV_Jdiff * CONSTANTS[AC_vss] / CONSTANTS[AC_vmyo1] - AV_Jgap);
    AV_diff_cai2 = AV_Bcai2 * (STATES[CICR_Jrel2] * (CONSTANTS[AC_vcsr] / CONSTANTS[AC_vmyo2]) + AV_Jgap * CONSTANTS[AC_vmyo2] / CONSTANTS[AC_vmyo1] - AV_Jup2 * CONSTANTS[AC_vnsr2] / CONSTANTS[AC_vmyo2]);
    AV_diff_cajsr = AV_Bcajsr * (AV_Jtr - AV_Jrel);
    AV_diff_cass = AV_Bcass * (-(ALGEBRAIC[AV_ICaL_ICaL] - 2.0 * ALGEBRAIC[AV_INaCa_ss_INaCa_ss]) * CONSTANTS[AC_Acap] / (2.0 * 2.0 * CONSTANTS[AC_F] * CONSTANTS[AC_vss]) + AV_Jrel * CONSTANTS[AC_vjsr] / CONSTANTS[AC_vss] - AV_Jdiff);
    RATES[ionic_concentrations_cacsr] = AV_diff_cacsr;
    RATES[ionic_concentrations_cai] = AV_diff_cai;
    RATES[ionic_concentrations_cai2] = AV_diff_cai2;
    RATES[ionic_concentrations_cajsr] = AV_diff_cajsr;
    RATES[ionic_concentrations_cass] = AV_diff_cass;
    
    /* *remaining* */
    AV_kito2 = 1.0 - 1.0 / (1.0 + AV_Jrel / 0.4 * AV_Jrel / 0.4);
    AV_Rel1 = (-ALGEBRAIC[AV_ICaL_ICaL] + 2.0 * ALGEBRAIC[AV_INaCa_ss_INaCa_ss]) * (CONSTANTS[AC_Acap] / (2.0 * CONSTANTS[AC_vss] * CONSTANTS[AC_F])) + AV_Jrel * (CONSTANTS[AC_vjsr] / CONSTANTS[AC_vss]) - AV_Jdiff;
    AV_Rel2 = AV_Jgap + STATES[CICR_Jrel2] * (CONSTANTS[AC_vcsr] / CONSTANTS[AC_vmyo2]) - AV_Jup2 * (CONSTANTS[AC_vnsr2] / CONSTANTS[AC_vmyo2]);
    ALGEBRAIC[AV_ITo_ITo] = CONSTANTS[AC_Gto] * STATES[ITo_aa] * AV_rito2 * AV_kito2 * (STATES[cell_v] - CONSTANTS[AC_ECl]);
    AV_Jrel1_inf = ((AV_Rel1 > 0.0) ? 1.0 * 60.0 * AV_Rel1 * (1.0 + 1.0 / (1.0 + pow(CONSTANTS[AC_KmCaMK] / AV_CaMKa, 8.0))) / (1.0 + pow(0.75 / STATES[ionic_concentrations_cajsr], 8.0)) : 0.0);
    AV_Jrel2_inf = ((AV_Rel2 > 0.0) ? 1.0 * (250.0 * AV_Rel2 * (1.0 + 1.0 / (1.0 + pow(CONSTANTS[AC_KmCaMK] / AV_CaMKa, 8.0)))) / (1.0 + pow(0.75 / STATES[ionic_concentrations_cacsr], 8.0)) : 0.0);
    RATES[CICR_Jrel1] = (AV_Jrel1_inf - STATES[CICR_Jrel1]) / AV_tau_Jrel1;
    RATES[CICR_Jrel2] = (AV_Jrel2_inf - STATES[CICR_Jrel2]) / AV_tau_Jrel2;
    ALGEBRAIC[AV_Iion] = ALGEBRAIC[AV_INa] + ALGEBRAIC[AV_INaL_INaL] + ALGEBRAIC[AV_ICaL_ICaL] + ALGEBRAIC[AV_ICaNa] + ALGEBRAIC[AV_ICaK] + ALGEBRAIC[AV_IKr_IKr] + ALGEBRAIC[AV_IKs_IKs] + ALGEBRAIC[AV_IK1_IK1] + ALGEBRAIC[AV_ITo_ITo] + ALGEBRAIC[AV_INaCa] + ALGEBRAIC[AV_INaK_INaK] + ALGEBRAIC[AV_INab_INab] + ALGEBRAIC[AV_IKb_IKb] + ALGEBRAIC[AV_IpCa_IpCa] + ALGEBRAIC[AV_ICab_ICab];
    RATES[cell_v] = -ALGEBRAIC[AV_Iion] - ALGEBRAIC[AV_I_stim]; 
}
