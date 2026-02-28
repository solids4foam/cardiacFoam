#include "ToRORd_dynCl_2023Names.H"
#include "stimulusIO.H"

static const char* ToRORd_dynClSTATES_NAMES[NUM_STATES] = {
    "V",
    "CaMKt",
    "Nai",
    "Nass",
    "Ki",
    "Kss",
    "Cass",
    "Cansr",
    "Cajsr",
    "Cai",
    "Cli",
    "Clss",
    "INa_m",
    "INa_h",
    "INa_j",
    "INa_hp",
    "INa_jp",
    "INaL_mL",
    "INaL_hL",
    "INaL_hLp",
    "Ito_a",
    "Ito_iF",
    "Ito_iS",
    "Ito_ap",
    "Ito_iFp",
    "Ito_iSp",
    "ICaL_d",
    "ICaL_ff",
    "ICaL_fs",
    "ICaL_fcaf",
    "ICaL_fcas",
    "ICaL_jca",
    "ICaL_ffp",
    "ICaL_fcafp",
    "ICaL_nca_ss",
    "ICaL_nca_i",
    "IKr_C1",
    "IKr_C2",
    "IKr_C3",
    "IKr_I",
    "IKr_O",
    "IKs_xs1",
    "IKs_xs2",
    "Jrel_np",
    "Jrel_p",
};

static const char* ToRORd_dynClALGEBRAIC_NAMES[NUM_ALGEBRAIC] = {
    "AV_CaMKb",
    "AV_CaMKa",
    "AV_IpCa_IpCa",
    "AV_Jdiff",
    "AV_JdiffCl",
    "AV_JdiffK",
    "AV_JdiffNa",
    "AV_time",
    "AV_Jtr",
    "AV_Jleak",
    "AV_fJupp",
    "AV_Jupnp",
    "AV_Jupp",
    "AV_Jup",
    "AV_ECl",
    "AV_EClss",
    "AV_EK",
    "AV_ENa",
    "AV_EKs",
    "AV_IClCa_junc",
    "AV_IClCa_sl",
    "AV_IClb",
    "AV_IClCa",
    "AV_aK1",
    "AV_bK1",
    "AV_K1ss",
    "AV_IK1_IK1",
    "AV_xkb",
    "AV_IKb_IKb",
    "AV_KsCa",
    "AV_txs1",
    "AV_txs2",
    "AV_xs1ss",
    "AV_xs2ss",
    "AV_IKs_IKs",
    "AV_ah",
    "AV_aj",
    "AV_bh",
    "AV_bj",
    "AV_fINap",
    "AV_hss",
    "AV_hssp",
    "AV_mss",
    "AV_tm",
    "AV_INa_INa",
    "AV_jss",
    "AV_th",
    "AV_tj",
    "AV_tjp",
    "AV_fINaLp",
    "AV_hLss",
    "AV_hLssp",
    "AV_mLss",
    "AV_tmL",
    "AV_INaL_INaL",
    "AV_I_katp_I_katp",
    "AV_fItop",
    "AV_AiF",
    "AV_ass",
    "AV_assp",
    "AV_delta_epi",
    "AV_dti_develop",
    "AV_dti_recover",
    "AV_iss",
    "AV_ta",
    "AV_tiF_b",
    "AV_tiS_b",
    "AV_AiS",
    "AV_tiF",
    "AV_tiS",
    "AV_i",
    "AV_ip",
    "AV_tiFp",
    "AV_tiSp",
    "AV_Ito_Ito",
    "AV_Afcaf",
    "AV_Ii",
    "AV_Iss",
    "AV_dss",
    "AV_fICaLp",
    "AV_fss",
    "AV_jcass",
    "AV_km2n",
    "AV_tfcaf",
    "AV_tfcas",
    "AV_tff",
    "AV_tfs",
    "AV_Afcas",
    "AV_anca_i",
    "AV_anca_ss",
    "AV_fcass",
    "AV_td",
    "AV_tfcafp",
    "AV_tffp",
    "AV_f",
    "AV_fca",
    "AV_fcap",
    "AV_fp",
    "AV_gamma_cai",
    "AV_gamma_cass",
    "AV_gamma_ki",
    "AV_gamma_kss",
    "AV_gamma_nai",
    "AV_gamma_nass",
    "AV_IKr_IKr",
    "AV_allo_i",
    "AV_allo_ss",
    "AV_h4_i",
    "AV_h4_ss",
    "AV_h5_i",
    "AV_h5_ss",
    "AV_h6_i",
    "AV_h6_ss",
    "AV_k6_i",
    "AV_k6_ss",
    "AV_P",
    "AV_b3",
    "AV_Bcajsr",
    "AV_Bcass",
    "AV_Bcai",
    "AV_vffrt",
    "AV_vfrt",
    "AV_fJrelp",
    "AV_Jrel",
    "AV_tau_rel_b",
    "AV_tau_rel",
    "AV_tau_relp_b",
    "AV_tau_relp",
    "AV_PhiCaK_i",
    "AV_PhiCaK_ss",
    "AV_PhiCaL_i",
    "AV_PhiCaL_ss",
    "AV_PhiCaNa_i",
    "AV_PhiCaNa_ss",
    "AV_ICab_ICab",
    "AV_alpha",
    "AV_alpha_2",
    "AV_alpha_C2ToI",
    "AV_alpha_i",
    "AV_beta",
    "AV_beta_2",
    "AV_beta_i",
    "AV_hca",
    "AV_hna",
    "AV_Knai",
    "AV_Knao",
    "AV_INab_INab",
    "AV_ICaK_i",
    "AV_ICaK_ss",
    "AV_ICaL_i",
    "AV_ICaL_ss",
    "AV_ICaNa_i",
    "AV_ICaNa_ss",
    "AV_beta_ItoC2",
    "AV_h1_i",
    "AV_h1_ss",
    "AV_h7_i",
    "AV_h7_ss",
    "AV_a1",
    "AV_a3",
    "AV_b2",
    "AV_b4",
    "AV_ICaK",
    "AV_ICaL_ICaL",
    "AV_ICaNa",
    "AV_h2_i",
    "AV_h2_ss",
    "AV_h3_i",
    "AV_h3_ss",
    "AV_h8_i",
    "AV_h8_ss",
    "AV_h9_i",
    "AV_h9_ss",
    "AV_x1",
    "AV_x2",
    "AV_x3",
    "AV_x4",
    "AV_Jrel_inf_b",
    "AV_Jrel_infp_b",
    "AV_k3p_i",
    "AV_k3p_ss",
    "AV_k3pp_i",
    "AV_k3pp_ss",
    "AV_k4p_i",
    "AV_k4p_ss",
    "AV_k4pp_i",
    "AV_k4pp_ss",
    "AV_k7_i",
    "AV_k7_ss",
    "AV_k8_i",
    "AV_k8_ss",
    "AV_E1",
    "AV_E2",
    "AV_E3",
    "AV_E4",
    "AV_Jrel_inf",
    "AV_Jrel_infp",
    "AV_k3_i",
    "AV_k3_ss",
    "AV_k4_i",
    "AV_k4_ss",
    "AV_JnakK",
    "AV_JnakNa",
    "AV_x1_i",
    "AV_x1_ss",
    "AV_x2_i",
    "AV_x2_ss",
    "AV_x3_i",
    "AV_x3_ss",
    "AV_x4_i",
    "AV_x4_ss",
    "AV_INaK_INaK",
    "AV_E1_i",
    "AV_E1_ss",
    "AV_E2_i",
    "AV_E2_ss",
    "AV_E3_i",
    "AV_E3_ss",
    "AV_E4_i",
    "AV_E4_ss",
    "AV_JncxCa_i",
    "AV_JncxCa_ss",
    "AV_JncxNa_i",
    "AV_JncxNa_ss",
    "AV_INaCa_i",
    "AV_INaCa_ss",
    "Istim",
    "Iion_cm"
};


static inline const Foam::HashTable<Foam::wordList>&
ToRORd_dynClDependencyMap()
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

inline Foam::scalar tissueSelect
(
    const Foam::label tissueFlag,
    const Foam::scalar v1,
    const Foam::scalar v2,
    const Foam::scalar v3
)
{
    return (tissueFlag == 3) ? v1 :
           (tissueFlag == 2) ? v2 : v3;
}

void
ToRORd_dynClinitConsts(double* CONSTANTS, double* RATES, double* STATES, int tissueFlag, const Foam::dictionary& stimulus)
{
/* CaMK */
    CONSTANTS[AC_CaMKo] = 0.05;
    CONSTANTS[AC_KmCaM] = 0.0015;
    CONSTANTS[AC_KmCaMK] = 0.15;
    CONSTANTS[AC_aCaMK] = 0.05;
    CONSTANTS[AC_bCaMK] = 0.00068;
    
    /* IpCa */
    CONSTANTS[AC_GpCa] = 0.0005;
    CONSTANTS[AC_KmCap] = 0.0005;
    
    /* cell_geometry */
    CONSTANTS[AC_L] = 0.01;
    CONSTANTS[AC_rad] = 0.0011;
    CONSTANTS[AC_Ageo] = 2.0 * 3.14 * CONSTANTS[AC_rad] * CONSTANTS[AC_rad] + 2.0 * 3.14 * CONSTANTS[AC_rad] * CONSTANTS[AC_L];
    CONSTANTS[AC_vcell] = 1000.0 * 3.14 * CONSTANTS[AC_rad] * CONSTANTS[AC_rad] * CONSTANTS[AC_L];
    CONSTANTS[AC_Acap] = 2.0 * CONSTANTS[AC_Ageo];
    CONSTANTS[AC_vjsr] = 0.0048 * CONSTANTS[AC_vcell];
    CONSTANTS[AC_vmyo] = 0.68 * CONSTANTS[AC_vcell];
    CONSTANTS[AC_vnsr] = 0.0552 * CONSTANTS[AC_vcell];
    CONSTANTS[AC_vss] = 0.02 * CONSTANTS[AC_vcell];
    
    /* diff */
    CONSTANTS[AC_tauCa] = 0.2;
    CONSTANTS[AC_tauCl] = 2.0;
    CONSTANTS[AC_tauK] = 2.0;
    CONSTANTS[AC_tauNa] = 2.0;

    
    /* extracellular */
    CONSTANTS[AC_cao] = 1.8;
    CONSTANTS[AC_clo] = 150.0;
    CONSTANTS[AC_ko] = 5.0;
    CONSTANTS[AC_nao] = 140.0;
    
    /* physical_constants */
    CONSTANTS[AC_F] = 96485.0;
    CONSTANTS[AC_R] = 8314.0;
    CONSTANTS[AC_T] = 310.0;
    CONSTANTS[AC_zca] = 2.0;
    CONSTANTS[AC_zcl] = -1.0;
    CONSTANTS[AC_zk] = 1.0;
    CONSTANTS[AC_zna] = 1.0;
    
    /* SERCA */
    CONSTANTS[AC_Jup_b] = 1.0;
    CONSTANTS[AC_upScale] = ((tissueFlag == 1) ? 1.3 : 1.0);
    
    /* reversal_potentials */
    CONSTANTS[AC_PKNa] = 0.01833;
    
    /* ICl */
    CONSTANTS[AC_Fjunc] = 1.0;
    CONSTANTS[AC_GClCa] = 0.2843;
    CONSTANTS[AC_GClb] = 0.00198;
    CONSTANTS[AC_KdClCa] = 0.1;
    
    /* IK1 */
    CONSTANTS[AC_GK1_b] = 0.6992;
    CONSTANTS[AC_GK1] = ((tissueFlag == 1) ? CONSTANTS[AC_GK1_b] * 1.2 : ((tissueFlag == 2.0) ? CONSTANTS[AC_GK1_b] * 1.3 : CONSTANTS[AC_GK1_b]));
    
    /* IKb */
    CONSTANTS[AC_GKb_b] = 0.0189;
    CONSTANTS[AC_GKb] = ((tissueFlag == 1) ? CONSTANTS[AC_GKb_b] * 0.6 : CONSTANTS[AC_GKb_b]);
    
    /* IKs */
    CONSTANTS[AC_GKs_b] = 0.0011;
    CONSTANTS[AC_GKs] = ((tissueFlag == 1) ? CONSTANTS[AC_GKs_b] * 1.4 : CONSTANTS[AC_GKs_b]);
    
    /* INa */
    CONSTANTS[AC_GNa] = 11.7802;
    
    /* INaL */
    CONSTANTS[AC_GNaL_b] = 0.0279;
    CONSTANTS[AC_thL] = 200.0;
    CONSTANTS[AC_GNaL] = ((tissueFlag == 1) ? CONSTANTS[AC_GNaL_b] * 0.6 : CONSTANTS[AC_GNaL_b]);
    CONSTANTS[AC_thLp] = 3.0 * CONSTANTS[AC_thL];
    
    /* I_katp */
    CONSTANTS[AC_A_atp] = 2.0;
    CONSTANTS[AC_K_atp] = 0.25;
    CONSTANTS[AC_K_o_n] = 5.0;
    CONSTANTS[AC_fkatp] = 0.0;
    CONSTANTS[AC_gkatp] = 4.3195;
    CONSTANTS[AC_akik] = pow(CONSTANTS[AC_ko] / CONSTANTS[AC_K_o_n], 0.24);
    CONSTANTS[AC_bkik] = 1.0 / (1.0 + pow(CONSTANTS[AC_A_atp] / CONSTANTS[AC_K_atp], 2.0));
    
    /* Ito */
    CONSTANTS[AC_EKshift] = 0.0;
    CONSTANTS[AC_Gto_b] = 0.16;
    CONSTANTS[AC_Gto] = ((tissueFlag == 1) ? CONSTANTS[AC_Gto_b] * 2.0 : ((tissueFlag == 2.0) ? CONSTANTS[AC_Gto_b] * 2.0 : CONSTANTS[AC_Gto_b]));
    
    /* ICaL */
    CONSTANTS[AC_Aff] = 0.6;
    CONSTANTS[AC_ICaL_fractionSS] = 0.8;
    CONSTANTS[AC_Io] = 0.5 * (CONSTANTS[AC_nao] + CONSTANTS[AC_ko] + CONSTANTS[AC_clo] + 4.0 * CONSTANTS[AC_cao]) / 1000.0;
    CONSTANTS[AC_Kmn] = 0.002;
    CONSTANTS[AC_PCa_b] = 8.37570000000000046e-05;
    CONSTANTS[AC_dielConstant] = 74.0;
    CONSTANTS[AC_k2n] = 500.0;
    CONSTANTS[AC_offset] = 0.0;
    CONSTANTS[AC_tjca] = 72.5;
    CONSTANTS[AC_vShift] = 0.0;
    CONSTANTS[AC_Afs] = 1.0 - CONSTANTS[AC_Aff];
    CONSTANTS[AC_PCa] = ((tissueFlag == 1) ? CONSTANTS[AC_PCa_b] * 1.2 : ((tissueFlag == 2.0) ? CONSTANTS[AC_PCa_b] * 2.0 : CONSTANTS[AC_PCa_b]));
    CONSTANTS[AC_constA] = 1820000.0 * pow(CONSTANTS[AC_dielConstant] * CONSTANTS[AC_T],
                                           -1.5);
    CONSTANTS[AC_PCaK] = 0.0003574 * CONSTANTS[AC_PCa];
    CONSTANTS[AC_PCaNa] = 0.00125 * CONSTANTS[AC_PCa];
    CONSTANTS[AC_PCap] = 1.1 * CONSTANTS[AC_PCa];
    CONSTANTS[AC_gamma_cao] = exp(-CONSTANTS[AC_constA] * 4.0 * (sqrt(CONSTANTS[AC_Io]) / (1.0 + sqrt(CONSTANTS[AC_Io])) - 0.3 * CONSTANTS[AC_Io]));
    CONSTANTS[AC_gamma_ko] = exp(-CONSTANTS[AC_constA] * 1.0 * (sqrt(CONSTANTS[AC_Io]) / (1.0 + sqrt(CONSTANTS[AC_Io])) - 0.3 * CONSTANTS[AC_Io]));
    CONSTANTS[AC_gamma_nao] = exp(-CONSTANTS[AC_constA] * 1.0 * (sqrt(CONSTANTS[AC_Io]) / (1.0 + sqrt(CONSTANTS[AC_Io])) - 0.3 * CONSTANTS[AC_Io]));
    CONSTANTS[AC_PCaKp] = 0.0003574 * CONSTANTS[AC_PCap];
    CONSTANTS[AC_PCaNap] = 0.00125 * CONSTANTS[AC_PCap];
    
    /* ICab */
    CONSTANTS[AC_PCab] = 5.91940000000000005e-08;
    
    /* IKr */
    CONSTANTS[AC_GKr_b] = 0.0321;
    CONSTANTS[AC_alpha_1] = 0.154375;
    CONSTANTS[AC_beta_1] = 0.1911;
    CONSTANTS[AC_GKr] = ((tissueFlag == 1) ? CONSTANTS[AC_GKr_b] * 1.3 : ((tissueFlag == 2.0) ? CONSTANTS[AC_GKr_b] * 0.8 : CONSTANTS[AC_GKr_b]));
    
    /* INaCa */
    CONSTANTS[AC_Gncx_b] = 0.0034;
    CONSTANTS[AC_INaCa_fractionSS] = 0.35;
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
    CONSTANTS[AC_Gncx] = ((tissueFlag == 1) ? CONSTANTS[AC_Gncx_b] * 1.1 : ((tissueFlag == 2.0) ? CONSTANTS[AC_Gncx_b] * 1.4 : CONSTANTS[AC_Gncx_b]));
    CONSTANTS[AC_h10_i] = CONSTANTS[AC_kasymm] + 1.0 + CONSTANTS[AC_nao] / CONSTANTS[AC_kna1] * (1.0 + CONSTANTS[AC_nao] / CONSTANTS[AC_kna2]);
    CONSTANTS[AC_h10_ss] = CONSTANTS[AC_kasymm] + 1.0 + CONSTANTS[AC_nao] / CONSTANTS[AC_kna1] * (1.0 + CONSTANTS[AC_nao] / CONSTANTS[AC_kna2]);
    CONSTANTS[AC_k2_i] = CONSTANTS[AC_kcaoff];
    CONSTANTS[AC_k2_ss] = CONSTANTS[AC_kcaoff];
    CONSTANTS[AC_k5_i] = CONSTANTS[AC_kcaoff];
    CONSTANTS[AC_k5_ss] = CONSTANTS[AC_kcaoff];
    CONSTANTS[AC_h11_i] = CONSTANTS[AC_nao] * CONSTANTS[AC_nao] / (CONSTANTS[AC_h10_i] * CONSTANTS[AC_kna1] * CONSTANTS[AC_kna2]);
    CONSTANTS[AC_h11_ss] = CONSTANTS[AC_nao] * CONSTANTS[AC_nao] / (CONSTANTS[AC_h10_ss] * CONSTANTS[AC_kna1] * CONSTANTS[AC_kna2]);
    CONSTANTS[AC_h12_i] = 1.0 / CONSTANTS[AC_h10_i];
    CONSTANTS[AC_h12_ss] = 1.0 / CONSTANTS[AC_h10_ss];
    CONSTANTS[AC_k1_i] = CONSTANTS[AC_h12_i] * CONSTANTS[AC_cao] * CONSTANTS[AC_kcaon];
    CONSTANTS[AC_k1_ss] = CONSTANTS[AC_h12_ss] * CONSTANTS[AC_cao] * CONSTANTS[AC_kcaon];
    
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
    CONSTANTS[AC_Pnak_b] = 15.4509;
    CONSTANTS[AC_delta] = -0.155;
    CONSTANTS[AC_eP] = 4.2;
    CONSTANTS[AC_k1m] = 182.4;
    CONSTANTS[AC_k1p] = 949.5;
    CONSTANTS[AC_k2m] = 39.4;
    CONSTANTS[AC_k2p] = 687.2;
    CONSTANTS[AC_k3m] = 79300.0;
    CONSTANTS[AC_k3p] = 1899.0;
    CONSTANTS[AC_k4m] = 40.0;
    CONSTANTS[AC_k4p] = 639.0;
    CONSTANTS[AC_Pnak] = ((tissueFlag == 1) ? CONSTANTS[AC_Pnak_b] * 0.9 : ((tissueFlag == 2.0) ? CONSTANTS[AC_Pnak_b] * 0.7 : CONSTANTS[AC_Pnak_b]));
    CONSTANTS[AC_a2] = CONSTANTS[AC_k2p];
    CONSTANTS[AC_a4] = CONSTANTS[AC_k4p] * CONSTANTS[AC_MgATP] / CONSTANTS[AC_Kmgatp] / (1.0 + CONSTANTS[AC_MgATP] / CONSTANTS[AC_Kmgatp]);
    CONSTANTS[AC_b1] = CONSTANTS[AC_k1m] * CONSTANTS[AC_MgADP];
    
    /* INab */
    CONSTANTS[AC_PNab] = 1.9239e-09;
    
    /* intracellular_ions */
    CONSTANTS[AC_BSLmax] = 1.124;
    CONSTANTS[AC_BSRmax] = 0.047;
    CONSTANTS[AC_KmBSL] = 0.0087;
    CONSTANTS[AC_KmBSR] = 0.00087;
    CONSTANTS[AC_cmdnmax_b] = 0.05;
    CONSTANTS[AC_csqnmax] = 10.0;
    CONSTANTS[AC_kmcmdn] = 0.00238;
    CONSTANTS[AC_kmcsqn] = 0.8;
    CONSTANTS[AC_kmtrpn] = 0.0005;
    CONSTANTS[AC_trpnmax] = 0.07;
    CONSTANTS[AC_cmdnmax] = ((tissueFlag == 1) ? CONSTANTS[AC_cmdnmax_b] * 1.3 : CONSTANTS[AC_cmdnmax_b]);
    
    
    /* ryr */
    CONSTANTS[AC_Jrel_b] = 1.5378;
    CONSTANTS[AC_bt] = 4.75;
    CONSTANTS[AC_cajsr_half] = 1.7;
    CONSTANTS[AC_a_rel] = 0.5 * CONSTANTS[AC_bt] / 1.0;
    CONSTANTS[AC_btp] = 1.25 * CONSTANTS[AC_bt];
    CONSTANTS[AC_a_relp] = 0.5 * CONSTANTS[AC_btp] / 1.0;

    /* --- Initial values --- */

    STATES[V]          = tissueSelect(tissueFlag, -89.74808, -91.33918, -90.74563);
    STATES[CaMKt]      = tissueSelect(tissueFlag, 1.10e-02, 2.01882e-02, 1.27e-02);
    
    STATES[Nai]        = tissueSelect(tissueFlag, 12.39736, 15.94867, 13.40062);
    STATES[Nass]       = tissueSelect(tissueFlag, 12.3977,  15.94922, 13.40094);
    STATES[Ki]         = tissueSelect(tissueFlag, 147.7115, 156.7131, 152.3639);
    STATES[Kss]        = tissueSelect(tissueFlag, 147.7114, 156.7130, 152.3638);
    
    STATES[Cass]       = tissueSelect(tissueFlag, 6.50e-05, 6.64e-05, 5.75e-05);
    STATES[Cansr]      = tissueSelect(tissueFlag, 1.528001, 2.012225, 1.806794);
    STATES[Cajsr]      = tissueSelect(tissueFlag, 1.525693, 2.016415, 1.805047);
    STATES[Cai]        = tissueSelect(tissueFlag, 7.45e-05, 8.30e-05, 6.62e-05);
    
    STATES[Cli]        = tissueSelect(tissueFlag, 29.20698, 48.91277, 34.31721);
    STATES[Clss]       = tissueSelect(tissueFlag, 29.20696, 48.91274, 34.31719);
    
    /* INa */
    STATES[INa_m]      = tissueSelect(tissueFlag, 6.52e-04, 4.62e-04, 5.25e-04);
    STATES[INa_h]      = tissueSelect(tissueFlag, 0.8473267, 0.8739077, 0.8645148);
    STATES[INa_j]      = tissueSelect(tissueFlag, 0.8471657, 0.8737841, 0.8644571);
    STATES[INa_hp]     = tissueSelect(tissueFlag, 0.7018454, 0.7478972, 0.7313656);
    STATES[INa_jp]     = tissueSelect(tissueFlag, 0.8469014, 0.8735375, 0.8643527);
    
    /* INaL */
    STATES[INaL_mL]    = tissueSelect(tissueFlag, 1.35e-04, 9.99e-05, 1.12e-04);
    STATES[INaL_hL]    = tissueSelect(tissueFlag, 0.5566017, 0.5986118, 0.5916536);
    STATES[INaL_hLp]   = tissueSelect(tissueFlag, 0.3115491, 0.3339899, 0.3476812);
    
    /* Ito */
    STATES[Ito_a]      = tissueSelect(tissueFlag, 8.90e-04, 7.99e-04, 8.32e-04);
    STATES[Ito_iF]     = tissueSelect(tissueFlag, 0.9996716, 0.9997514, 0.9997242);
    STATES[Ito_iS]     = tissueSelect(tissueFlag, 0.5988908, 0.5702538, 0.9997235);
    STATES[Ito_ap]     = tissueSelect(tissueFlag, 4.53e-04, 4.07e-04, 4.24e-04);
    STATES[Ito_iFp]    = tissueSelect(tissueFlag, 0.9996716, 0.9997514, 0.9997242);
    STATES[Ito_iSp]    = tissueSelect(tissueFlag, 0.6620692, 0.6351927, 0.9997241);
    
    /* ICaL */
    STATES[ICaL_d]     = tissueSelect(tissueFlag, 1.59e-31, -8.33e-30, -2.49e-36);
    STATES[ICaL_ff]    = tissueSelect(tissueFlag, 1.0, 1.0, 1.0);
    STATES[ICaL_fs]    = tissueSelect(tissueFlag, 0.9401791, 0.9183587, 0.9510602);
    STATES[ICaL_fcaf]  = tissueSelect(tissueFlag, 1.0, 1.0, 1.0);
    STATES[ICaL_fcas]  = tissueSelect(tissueFlag, 0.9999014, 0.9997540, 0.9999377);
    STATES[ICaL_jca]   = tissueSelect(tissueFlag, 0.9999846, 0.9999743, 0.9999886);
    STATES[ICaL_ffp]   = tissueSelect(tissueFlag, 1.0, 1.0, 1.0);
    STATES[ICaL_fcafp] = tissueSelect(tissueFlag, 1.0, 1.0, 1.0);
    STATES[ICaL_nca_ss]= tissueSelect(tissueFlag, 4.90e-04, 5.34e-04, 3.05e-04);
    STATES[ICaL_nca_i] = tissueSelect(tissueFlag, 8.33e-04, 1.26e-03, 5.27e-04);
    
    /* IKr */
    STATES[IKr_C1]     = tissueSelect(tissueFlag, 0.9982511, 0.9983451, 0.9984733);
    STATES[IKr_C2]     = tissueSelect(tissueFlag, 7.94e-04, 7.09e-04, 7.39e-04);
    STATES[IKr_C3]     = tissueSelect(tissueFlag, 6.53e-04, 5.91e-04, 6.03e-04);
    STATES[IKr_I]      = tissueSelect(tissueFlag, 9.80e-06, 1.06e-05, 5.68e-06);
    STATES[IKr_O]      = tissueSelect(tissueFlag, 2.92e-04, 3.45e-04, 1.79e-04);
    
    /* IKs */
    STATES[IKs_xs1]    = tissueSelect(tissueFlag, 0.243959, 0.2642293, 0.2233584);
    STATES[IKs_xs2]    = tissueSelect(tissueFlag, 1.59e-04, 1.33e-04, 1.42e-04);
    
    /* RyR */
    STATES[Jrel_np]    = tissueSelect(tissueFlag, 1.81e-22, -1.30e-21, 6.78e-25);
    STATES[Jrel_p]     = tissueSelect(tissueFlag, 4.36e-21, -7.61e-20, -1.58e-23);
    

}
void
ToRORd_dynClcomputeVariables(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC, int tissueFlag, bool solveVmWithinODESolver)
{
/* CaMK */
    ALGEBRAIC[AV_CaMKb] = CONSTANTS[AC_CaMKo] * (1.0 - STATES[CaMKt]) / (1.0 + CONSTANTS[AC_KmCaM] / STATES[Cass]);
    ALGEBRAIC[AV_CaMKa] = ALGEBRAIC[AV_CaMKb] + STATES[CaMKt];
    RATES[CaMKt] = CONSTANTS[AC_aCaMK] * ALGEBRAIC[AV_CaMKb] * (ALGEBRAIC[AV_CaMKb] + STATES[CaMKt]) - CONSTANTS[AC_bCaMK] * STATES[CaMKt];
    
    /* IpCa */
    ALGEBRAIC[AV_IpCa_IpCa] = CONSTANTS[AC_GpCa] * STATES[Cai] / (CONSTANTS[AC_KmCap] + STATES[Cai]);
    
    /* diff */
    ALGEBRAIC[AV_Jdiff] = (STATES[Cass] - STATES[Cai]) / CONSTANTS[AC_tauCa];
    ALGEBRAIC[AV_JdiffCl] = (STATES[Clss] - STATES[Cli]) / CONSTANTS[AC_tauNa];
    ALGEBRAIC[AV_JdiffK] = (STATES[Kss] - STATES[Ki]) / CONSTANTS[AC_tauK];
    ALGEBRAIC[AV_JdiffNa] = (STATES[Nass] - STATES[Nai]) / CONSTANTS[AC_tauNa];
    
    
    /* trans_flux */
    ALGEBRAIC[AV_Jtr] = (STATES[Cansr] - STATES[Cajsr]) / 60.0;
    
    /* SERCA */
    ALGEBRAIC[AV_Jleak] = 0.0048825 * STATES[Cansr] / 15.0;
    ALGEBRAIC[AV_fJupp] = 1.0 / (1.0 + CONSTANTS[AC_KmCaMK] / ALGEBRAIC[AV_CaMKa]);
    ALGEBRAIC[AV_Jupnp] = CONSTANTS[AC_upScale] * 0.005425 * STATES[Cai] / (STATES[Cai] + 0.00092);
    ALGEBRAIC[AV_Jupp] = CONSTANTS[AC_upScale] * 2.75 * 0.005425 * STATES[Cai] / (STATES[Cai] + 0.00092 - 0.00017);
    ALGEBRAIC[AV_Jup] = CONSTANTS[AC_Jup_b] * ((1.0 - ALGEBRAIC[AV_fJupp]) * ALGEBRAIC[AV_Jupnp] + ALGEBRAIC[AV_fJupp] * ALGEBRAIC[AV_Jupp] - ALGEBRAIC[AV_Jleak]);
    
    /* reversal_potentials */
    ALGEBRAIC[AV_ECl] = CONSTANTS[AC_R] * CONSTANTS[AC_T] / (CONSTANTS[AC_zcl] * CONSTANTS[AC_F]) * log(CONSTANTS[AC_clo] / STATES[Cli]);
    ALGEBRAIC[AV_EClss] = CONSTANTS[AC_R] * CONSTANTS[AC_T] / (CONSTANTS[AC_zcl] * CONSTANTS[AC_F]) * log(CONSTANTS[AC_clo] / STATES[Clss]);
    ALGEBRAIC[AV_EK] = CONSTANTS[AC_R] * CONSTANTS[AC_T] / (CONSTANTS[AC_zk] * CONSTANTS[AC_F]) * log(CONSTANTS[AC_ko] / STATES[Ki]);
    ALGEBRAIC[AV_ENa] = CONSTANTS[AC_R] * CONSTANTS[AC_T] / (CONSTANTS[AC_zna] * CONSTANTS[AC_F]) * log(CONSTANTS[AC_nao] / STATES[Nai]);
    ALGEBRAIC[AV_EKs] = CONSTANTS[AC_R] * CONSTANTS[AC_T] / (CONSTANTS[AC_zk] * CONSTANTS[AC_F]) * log((CONSTANTS[AC_ko] + CONSTANTS[AC_PKNa] * CONSTANTS[AC_nao]) / (STATES[Ki] + CONSTANTS[AC_PKNa] * STATES[Nai]));
    
    /* ICl */
    ALGEBRAIC[AV_IClCa_junc] = CONSTANTS[AC_Fjunc] * CONSTANTS[AC_GClCa] / (1.0 + CONSTANTS[AC_KdClCa] / STATES[Cass]) * (STATES[V] - ALGEBRAIC[AV_EClss]);
    ALGEBRAIC[AV_IClCa_sl] = (1.0 - CONSTANTS[AC_Fjunc]) * CONSTANTS[AC_GClCa] / (1.0 + CONSTANTS[AC_KdClCa] / STATES[Cai]) * (STATES[V] - ALGEBRAIC[AV_ECl]);
    ALGEBRAIC[AV_IClb] = CONSTANTS[AC_GClb] * (STATES[V] - ALGEBRAIC[AV_ECl]);
    ALGEBRAIC[AV_IClCa] = ALGEBRAIC[AV_IClCa_junc] + ALGEBRAIC[AV_IClCa_sl];
    
    /* IK1 */
    ALGEBRAIC[AV_aK1] = 4.094 / (1.0 + exp(0.1217 * (STATES[V] - ALGEBRAIC[AV_EK] - 49.934)));
    ALGEBRAIC[AV_bK1] = (15.72 * exp(0.0674 * (STATES[V] - ALGEBRAIC[AV_EK] - 3.257)) + exp(0.0618 * (STATES[V] - ALGEBRAIC[AV_EK] - 594.31))) / (1.0 + exp(-0.1629 * (STATES[V] - ALGEBRAIC[AV_EK] + 14.207)));
    ALGEBRAIC[AV_K1ss] = ALGEBRAIC[AV_aK1] / (ALGEBRAIC[AV_aK1] + ALGEBRAIC[AV_bK1]);
    ALGEBRAIC[AV_IK1_IK1] = CONSTANTS[AC_GK1] * sqrt(CONSTANTS[AC_ko] / 5.0) * ALGEBRAIC[AV_K1ss] * (STATES[V] - ALGEBRAIC[AV_EK]);
    
    /* IKb */
    ALGEBRAIC[AV_xkb] = 1.0 / (1.0 + exp(-(STATES[V] - 10.8968) / 23.9871));
    ALGEBRAIC[AV_IKb_IKb] = CONSTANTS[AC_GKb] * ALGEBRAIC[AV_xkb] * (STATES[V] - ALGEBRAIC[AV_EK]);
    
    /* IKs */
    ALGEBRAIC[AV_KsCa] = 1.0 + 0.6 / (1.0 + pow(3.8e-05 / STATES[Cai], 1.4));
    ALGEBRAIC[AV_txs1] = 817.3 + 1.0 / (0.0002326 * exp((STATES[V] + 48.28) / 17.8) + 0.001292 * exp(-(STATES[V] + 210.0) / 230.0));
    ALGEBRAIC[AV_txs2] = 1.0 / (0.01 * exp((STATES[V] - 50.0) / 20.0) + 0.0193 * exp(-(STATES[V] + 66.54) / 31.0));
    ALGEBRAIC[AV_xs1ss] = 1.0 / (1.0 + exp(-(STATES[V] + 11.6) / 8.932));
    ALGEBRAIC[AV_xs2ss] = ALGEBRAIC[AV_xs1ss];
    RATES[IKs_xs1] = (ALGEBRAIC[AV_xs1ss] - STATES[IKs_xs1]) / ALGEBRAIC[AV_txs1];
    ALGEBRAIC[AV_IKs_IKs] = CONSTANTS[AC_GKs] * ALGEBRAIC[AV_KsCa] * STATES[IKs_xs1] * STATES[IKs_xs2] * (STATES[V] - ALGEBRAIC[AV_EKs]);
    RATES[IKs_xs2] = (ALGEBRAIC[AV_xs2ss] - STATES[IKs_xs2]) / ALGEBRAIC[AV_txs2];
    
    /* INa */
    ALGEBRAIC[AV_ah] = ((STATES[V] >= -40.0) ? 0.0 : 0.057 * exp(-(STATES[V] + 80.0) / 6.8));
    ALGEBRAIC[AV_aj] = ((STATES[V] >= -40.0) ? 0.0 : (-25428.0 * exp(0.2444 * STATES[V]) - 6.948e-06 * exp(-0.04391 * STATES[V])) * (STATES[V] + 37.78) / (1.0 + exp(0.311 * (STATES[V] + 79.23))));
    ALGEBRAIC[AV_bh] = ((STATES[V] >= -40.0) ? 0.77 / (0.13 * (1.0 + exp(-(STATES[V] + 10.66) / 11.1))) : 2.7 * exp(0.079 * STATES[V]) + 310000.0 * exp(0.3485 * STATES[V]));
    ALGEBRAIC[AV_bj] = ((STATES[V] >= -40.0) ? 0.6 * exp(0.057 * STATES[V]) / (1.0 + exp(-0.1 * (STATES[V] + 32.0))) : 0.02424 * exp(-0.01052 * STATES[V]) / (1.0 + exp(-0.1378 * (STATES[V] + 40.14))));
    ALGEBRAIC[AV_fINap] = 1.0 / (1.0 + CONSTANTS[AC_KmCaMK] / ALGEBRAIC[AV_CaMKa]);
    ALGEBRAIC[AV_hss] = 1.0 / pow(1.0 + exp((STATES[V] + 71.55) / 7.43), 2.0);
    ALGEBRAIC[AV_hssp] = 1.0 / pow(1.0 + exp((STATES[V] + 77.55) / 7.43), 2.0);
    ALGEBRAIC[AV_mss] = 1.0 / pow(1.0 + exp(-(STATES[V] + 56.86) / 9.03), 2.0);
    ALGEBRAIC[AV_tm] = 0.1292 * exp(-pow((STATES[V] + 45.79) / 15.54, 2.0)) + 0.06487 * exp(-pow((STATES[V] - 4.823) / 51.12, 2.0));
    ALGEBRAIC[AV_INa_INa] = CONSTANTS[AC_GNa] * (STATES[V] - ALGEBRAIC[AV_ENa]) * pow(STATES[INa_m], 3.0) * ((1.0 - ALGEBRAIC[AV_fINap]) * STATES[INa_h] * STATES[INa_j] + ALGEBRAIC[AV_fINap] * STATES[INa_hp] * STATES[INa_jp]);
    ALGEBRAIC[AV_jss] = ALGEBRAIC[AV_hss];
    ALGEBRAIC[AV_th] = 1.0 / (ALGEBRAIC[AV_ah] + ALGEBRAIC[AV_bh]);
    ALGEBRAIC[AV_tj] = 1.0 / (ALGEBRAIC[AV_aj] + ALGEBRAIC[AV_bj]);
    RATES[INa_m] = (ALGEBRAIC[AV_mss] - STATES[INa_m]) / ALGEBRAIC[AV_tm];
    ALGEBRAIC[AV_tjp] = 1.46 * ALGEBRAIC[AV_tj];
    RATES[INa_h] = (ALGEBRAIC[AV_hss] - STATES[INa_h]) / ALGEBRAIC[AV_th];
    RATES[INa_hp] = (ALGEBRAIC[AV_hssp] - STATES[INa_hp]) / ALGEBRAIC[AV_th];
    RATES[INa_j] = (ALGEBRAIC[AV_jss] - STATES[INa_j]) / ALGEBRAIC[AV_tj];
    RATES[INa_jp] = (ALGEBRAIC[AV_jss] - STATES[INa_jp]) / ALGEBRAIC[AV_tjp];
    
    /* INaL */
    ALGEBRAIC[AV_fINaLp] = 1.0 / (1.0 + CONSTANTS[AC_KmCaMK] / ALGEBRAIC[AV_CaMKa]);
    ALGEBRAIC[AV_hLss] = 1.0 / (1.0 + exp((STATES[V] + 87.61) / 7.488));
    ALGEBRAIC[AV_hLssp] = 1.0 / (1.0 + exp((STATES[V] + 93.81) / 7.488));
    ALGEBRAIC[AV_mLss] = 1.0 / (1.0 + exp(-(STATES[V] + 42.85) / 5.264));
    ALGEBRAIC[AV_tmL] = 0.1292 * exp(-pow((STATES[V] + 45.79) / 15.54, 2.0)) + 0.06487 * exp(-pow((STATES[V] - 4.823) / 51.12, 2.0));
    RATES[INaL_hL] = (ALGEBRAIC[AV_hLss] - STATES[INaL_hL]) / CONSTANTS[AC_thL];
    RATES[INaL_mL] = (ALGEBRAIC[AV_mLss] - STATES[INaL_mL]) / ALGEBRAIC[AV_tmL];
    ALGEBRAIC[AV_INaL_INaL] = CONSTANTS[AC_GNaL] * (STATES[V] - ALGEBRAIC[AV_ENa]) * STATES[INaL_mL] * ((1.0 - ALGEBRAIC[AV_fINaLp]) * STATES[INaL_hL] + ALGEBRAIC[AV_fINaLp] * STATES[INaL_hLp]);
    RATES[INaL_hLp] = (ALGEBRAIC[AV_hLssp] - STATES[INaL_hLp]) / CONSTANTS[AC_thLp];
    
    /* I_katp */
    ALGEBRAIC[AV_I_katp_I_katp] = CONSTANTS[AC_fkatp] * CONSTANTS[AC_gkatp] * CONSTANTS[AC_akik] * CONSTANTS[AC_bkik] * (STATES[V] - ALGEBRAIC[AV_EK]);
    
    /* Ito */
    ALGEBRAIC[AV_fItop] = 1.0 / (1.0 + CONSTANTS[AC_KmCaMK] / ALGEBRAIC[AV_CaMKa]);
    ALGEBRAIC[AV_AiF] = 1.0 / (1.0 + exp((STATES[V] + CONSTANTS[AC_EKshift] - 213.6) / 151.2));
    ALGEBRAIC[AV_ass] = 1.0 / (1.0 + exp(-(STATES[V] + CONSTANTS[AC_EKshift] - 14.34) / 14.82));
    ALGEBRAIC[AV_assp] = 1.0 / (1.0 + exp(-(STATES[V] + CONSTANTS[AC_EKshift] - 24.34) / 14.82));
    ALGEBRAIC[AV_delta_epi] = ((tissueFlag == 1) ? 1.0 - 0.95 / (1.0 + exp((STATES[V] + CONSTANTS[AC_EKshift] + 70.0) / 5.0)) : 1.0);
    ALGEBRAIC[AV_dti_develop] = 1.354 + 0.0001 / (exp((STATES[V] + CONSTANTS[AC_EKshift] - 167.4) / 15.89) + exp(-(STATES[V] + CONSTANTS[AC_EKshift] - 12.23) / 0.2154));
    ALGEBRAIC[AV_dti_recover] = 1.0 - 0.5 / (1.0 + exp((STATES[V] + CONSTANTS[AC_EKshift] + 70.0) / 20.0));
    ALGEBRAIC[AV_iss] = 1.0 / (1.0 + exp((STATES[V] + CONSTANTS[AC_EKshift] + 43.94) / 5.711));
    ALGEBRAIC[AV_ta] = 1.0515 / (1.0 / (1.2089 * (1.0 + exp(-(STATES[V] + CONSTANTS[AC_EKshift] - 18.4099) / 29.3814))) + 3.5 / (1.0 + exp((STATES[V] + CONSTANTS[AC_EKshift] + 100.0) / 29.3814)));
    ALGEBRAIC[AV_tiF_b] = 4.562 + 1.0 / (0.3933 * exp(-(STATES[V] + CONSTANTS[AC_EKshift] + 100.0) / 100.0) + 0.08004 * exp((STATES[V] + CONSTANTS[AC_EKshift] + 50.0) / 16.59));
    ALGEBRAIC[AV_tiS_b] = 23.62 + 1.0 / (0.001416 * exp(-(STATES[V] + CONSTANTS[AC_EKshift] + 96.52) / 59.05) + 1.78e-08 * exp((STATES[V] + CONSTANTS[AC_EKshift] + 114.1) / 8.079));
    ALGEBRAIC[AV_AiS] = 1.0 - ALGEBRAIC[AV_AiF];
    ALGEBRAIC[AV_tiF] = ALGEBRAIC[AV_tiF_b] * ALGEBRAIC[AV_delta_epi];
    ALGEBRAIC[AV_tiS] = ALGEBRAIC[AV_tiS_b] * ALGEBRAIC[AV_delta_epi];
    RATES[Ito_a] = (ALGEBRAIC[AV_ass] - STATES[Ito_a]) / ALGEBRAIC[AV_ta];
    RATES[Ito_ap] = (ALGEBRAIC[AV_assp] - STATES[Ito_ap]) / ALGEBRAIC[AV_ta];
    ALGEBRAIC[AV_i] = ALGEBRAIC[AV_AiF] * STATES[Ito_iF] + ALGEBRAIC[AV_AiS] * STATES[Ito_iS];
    ALGEBRAIC[AV_ip] = ALGEBRAIC[AV_AiF] * STATES[Ito_iFp] + ALGEBRAIC[AV_AiS] * STATES[Ito_iSp];
    ALGEBRAIC[AV_tiFp] = ALGEBRAIC[AV_dti_develop] * ALGEBRAIC[AV_dti_recover] * ALGEBRAIC[AV_tiF];
    ALGEBRAIC[AV_tiSp] = ALGEBRAIC[AV_dti_develop] * ALGEBRAIC[AV_dti_recover] * ALGEBRAIC[AV_tiS];
    RATES[Ito_iF] = (ALGEBRAIC[AV_iss] - STATES[Ito_iF]) / ALGEBRAIC[AV_tiF];
    RATES[Ito_iS] = (ALGEBRAIC[AV_iss] - STATES[Ito_iS]) / ALGEBRAIC[AV_tiS];
    ALGEBRAIC[AV_Ito_Ito] = CONSTANTS[AC_Gto] * (STATES[V] - ALGEBRAIC[AV_EK]) * ((1.0 - ALGEBRAIC[AV_fItop]) * STATES[Ito_a] * ALGEBRAIC[AV_i] + ALGEBRAIC[AV_fItop] * STATES[Ito_ap] * ALGEBRAIC[AV_ip]);
    RATES[Ito_iFp] = (ALGEBRAIC[AV_iss] - STATES[Ito_iFp]) / ALGEBRAIC[AV_tiFp];
    RATES[Ito_iSp] = (ALGEBRAIC[AV_iss] - STATES[Ito_iSp]) / ALGEBRAIC[AV_tiSp];
    
    /* ICaL */
    ALGEBRAIC[AV_Afcaf] = 0.3 + 0.6 / (1.0 + exp((STATES[V] - 10.0) / 10.0));
    ALGEBRAIC[AV_Ii] = 0.5 * (STATES[Nai] + STATES[Ki] + STATES[Cli] + 4.0 * STATES[Cai]) / 1000.0;
    ALGEBRAIC[AV_Iss] = 0.5 * (STATES[Nass] + STATES[Kss] + STATES[Clss] + 4.0 * STATES[Cass]) / 1000.0;
    ALGEBRAIC[AV_dss] = ((STATES[V] >= 31.4978) ? 1.0 : 1.0763 * exp(-1.007 * exp(-0.0829 * STATES[V])));
    ALGEBRAIC[AV_fICaLp] = 1.0 / (1.0 + CONSTANTS[AC_KmCaMK] / ALGEBRAIC[AV_CaMKa]);
    ALGEBRAIC[AV_fss] = 1.0 / (1.0 + exp((STATES[V] + 19.58) / 3.696));
    ALGEBRAIC[AV_jcass] = 1.0 / (1.0 + exp((STATES[V] + 18.08) / 2.7916));
    ALGEBRAIC[AV_km2n] = STATES[ICaL_jca] * 1.0;
    ALGEBRAIC[AV_tfcaf] = 7.0 + 1.0 / (0.04 * exp(-(STATES[V] - 4.0) / 7.0) + 0.04 * exp((STATES[V] - 4.0) / 7.0));
    ALGEBRAIC[AV_tfcas] = 100.0 + 1.0 / (0.00012 * exp(-STATES[V] / 3.0) + 0.00012 * exp(STATES[V] / 7.0));
    ALGEBRAIC[AV_tff] = 7.0 + 1.0 / (0.0045 * exp(-(STATES[V] + 20.0) / 10.0) + 0.0045 * exp((STATES[V] + 20.0) / 10.0));
    ALGEBRAIC[AV_tfs] = 1000.0 + 1.0 / (3.5e-05 * exp(-(STATES[V] + 5.0) / 4.0) + 3.5e-05 * exp((STATES[V] + 5.0) / 6.0));
    ALGEBRAIC[AV_Afcas] = 1.0 - ALGEBRAIC[AV_Afcaf];
    ALGEBRAIC[AV_anca_i] = 1.0 / (CONSTANTS[AC_k2n] / ALGEBRAIC[AV_km2n] + pow(1.0 + CONSTANTS[AC_Kmn] / STATES[Cai], 4.0));
    ALGEBRAIC[AV_anca_ss] = 1.0 / (CONSTANTS[AC_k2n] / ALGEBRAIC[AV_km2n] + pow(1.0 + CONSTANTS[AC_Kmn] / STATES[Cass], 4.0));
    ALGEBRAIC[AV_fcass] = ALGEBRAIC[AV_fss];
    ALGEBRAIC[AV_td] = CONSTANTS[AC_offset] + 0.6 + 1.0 / (exp(-0.05 * (STATES[V] + CONSTANTS[AC_vShift] + 6.0)) + exp(0.09 * (STATES[V] + CONSTANTS[AC_vShift] + 14.0)));
    ALGEBRAIC[AV_tfcafp] = 2.5 * ALGEBRAIC[AV_tfcaf];
    ALGEBRAIC[AV_tffp] = 2.5 * ALGEBRAIC[AV_tff];
    RATES[ICaL_ff] = (ALGEBRAIC[AV_fss] - STATES[ICaL_ff]) / ALGEBRAIC[AV_tff];
    RATES[ICaL_fs] = (ALGEBRAIC[AV_fss] - STATES[ICaL_fs]) / ALGEBRAIC[AV_tfs];
    RATES[ICaL_jca] = (ALGEBRAIC[AV_jcass] - STATES[ICaL_jca]) / CONSTANTS[AC_tjca];
    ALGEBRAIC[AV_f] = CONSTANTS[AC_Aff] * STATES[ICaL_ff] + CONSTANTS[AC_Afs] * STATES[ICaL_fs];
    ALGEBRAIC[AV_fca] = ALGEBRAIC[AV_Afcaf] * STATES[ICaL_fcaf] + ALGEBRAIC[AV_Afcas] * STATES[ICaL_fcas];
    ALGEBRAIC[AV_fcap] = ALGEBRAIC[AV_Afcaf] * STATES[ICaL_fcafp] + ALGEBRAIC[AV_Afcas] * STATES[ICaL_fcas];
    ALGEBRAIC[AV_fp] = CONSTANTS[AC_Aff] * STATES[ICaL_ffp] + CONSTANTS[AC_Afs] * STATES[ICaL_fs];
    ALGEBRAIC[AV_gamma_cai] = exp(-CONSTANTS[AC_constA] * 4.0 * (sqrt(ALGEBRAIC[AV_Ii]) / (1.0 + sqrt(ALGEBRAIC[AV_Ii])) - 0.3 * ALGEBRAIC[AV_Ii]));
    ALGEBRAIC[AV_gamma_cass] = exp(-CONSTANTS[AC_constA] * 4.0 * (sqrt(ALGEBRAIC[AV_Iss]) / (1.0 + sqrt(ALGEBRAIC[AV_Iss])) - 0.3 * ALGEBRAIC[AV_Iss]));
    ALGEBRAIC[AV_gamma_ki] = exp(-CONSTANTS[AC_constA] * 1.0 * (sqrt(ALGEBRAIC[AV_Ii]) / (1.0 + sqrt(ALGEBRAIC[AV_Ii])) - 0.3 * ALGEBRAIC[AV_Ii]));
    ALGEBRAIC[AV_gamma_kss] = exp(-CONSTANTS[AC_constA] * 1.0 * (sqrt(ALGEBRAIC[AV_Iss]) / (1.0 + sqrt(ALGEBRAIC[AV_Iss])) - 0.3 * ALGEBRAIC[AV_Iss]));
    ALGEBRAIC[AV_gamma_nai] = exp(-CONSTANTS[AC_constA] * 1.0 * (sqrt(ALGEBRAIC[AV_Ii]) / (1.0 + sqrt(ALGEBRAIC[AV_Ii])) - 0.3 * ALGEBRAIC[AV_Ii]));
    ALGEBRAIC[AV_gamma_nass] = exp(-CONSTANTS[AC_constA] * 1.0 * (sqrt(ALGEBRAIC[AV_Iss]) / (1.0 + sqrt(ALGEBRAIC[AV_Iss])) - 0.3 * ALGEBRAIC[AV_Iss]));
    RATES[ICaL_d] = (ALGEBRAIC[AV_dss] - STATES[ICaL_d]) / ALGEBRAIC[AV_td];
    RATES[ICaL_fcaf] = (ALGEBRAIC[AV_fcass] - STATES[ICaL_fcaf]) / ALGEBRAIC[AV_tfcaf];
    RATES[ICaL_fcafp] = (ALGEBRAIC[AV_fcass] - STATES[ICaL_fcafp]) / ALGEBRAIC[AV_tfcafp];
    RATES[ICaL_fcas] = (ALGEBRAIC[AV_fcass] - STATES[ICaL_fcas]) / ALGEBRAIC[AV_tfcas];
    RATES[ICaL_ffp] = (ALGEBRAIC[AV_fss] - STATES[ICaL_ffp]) / ALGEBRAIC[AV_tffp];
    RATES[ICaL_nca_i] = ALGEBRAIC[AV_anca_i] * CONSTANTS[AC_k2n] - STATES[ICaL_nca_i] * ALGEBRAIC[AV_km2n];
    RATES[ICaL_nca_ss] = ALGEBRAIC[AV_anca_ss] * CONSTANTS[AC_k2n] - STATES[ICaL_nca_ss] * ALGEBRAIC[AV_km2n];
    
    /* IKr */
    ALGEBRAIC[AV_IKr_IKr] = CONSTANTS[AC_GKr] * sqrt(CONSTANTS[AC_ko] / 5.0) * STATES[IKr_O] * (STATES[V] - ALGEBRAIC[AV_EK]);
    
    /* INaCa */
    ALGEBRAIC[AV_allo_i] = 1.0 / (1.0 + pow(CONSTANTS[AC_KmCaAct] / STATES[Cai], 2.0));
    ALGEBRAIC[AV_allo_ss] = 1.0 / (1.0 + pow(CONSTANTS[AC_KmCaAct] / STATES[Cass], 2.0));
    ALGEBRAIC[AV_h4_i] = 1.0 + STATES[Nai] / CONSTANTS[AC_kna1] * (1.0 + STATES[Nai] / CONSTANTS[AC_kna2]);
    ALGEBRAIC[AV_h4_ss] = 1.0 + STATES[Nass] / CONSTANTS[AC_kna1] * (1.0 + STATES[Nass] / CONSTANTS[AC_kna2]);
    ALGEBRAIC[AV_h5_i] = STATES[Nai] * STATES[Nai] / (ALGEBRAIC[AV_h4_i] * CONSTANTS[AC_kna1] * CONSTANTS[AC_kna2]);
    ALGEBRAIC[AV_h5_ss] = STATES[Nass] * STATES[Nass] / (ALGEBRAIC[AV_h4_ss] * CONSTANTS[AC_kna1] * CONSTANTS[AC_kna2]);
    ALGEBRAIC[AV_h6_i] = 1.0 / ALGEBRAIC[AV_h4_i];
    ALGEBRAIC[AV_h6_ss] = 1.0 / ALGEBRAIC[AV_h4_ss];
    ALGEBRAIC[AV_k6_i] = ALGEBRAIC[AV_h6_i] * STATES[Cai] * CONSTANTS[AC_kcaon];
    ALGEBRAIC[AV_k6_ss] = ALGEBRAIC[AV_h6_ss] * STATES[Cass] * CONSTANTS[AC_kcaon];
    
    /* INaK */
    ALGEBRAIC[AV_P] = CONSTANTS[AC_eP] / (1.0 + CONSTANTS[AC_H] / CONSTANTS[AC_Khp] + STATES[Nai] / CONSTANTS[AC_Knap] + STATES[Ki] / CONSTANTS[AC_Kxkur]);
    ALGEBRAIC[AV_b3] = CONSTANTS[AC_k3m] * ALGEBRAIC[AV_P] * CONSTANTS[AC_H] / (1.0 + CONSTANTS[AC_MgATP] / CONSTANTS[AC_Kmgatp]);
    
    /* intracellular_ions */
    RATES[Cansr] = ALGEBRAIC[AV_Jup] - ALGEBRAIC[AV_Jtr] * CONSTANTS[AC_vjsr] / CONSTANTS[AC_vnsr];
    RATES[Cli] = (ALGEBRAIC[AV_IClb] + ALGEBRAIC[AV_IClCa_sl]) * CONSTANTS[AC_Acap] / (CONSTANTS[AC_F] * CONSTANTS[AC_vmyo]) + ALGEBRAIC[AV_JdiffCl] * CONSTANTS[AC_vss] / CONSTANTS[AC_vmyo];
    RATES[Clss] = ALGEBRAIC[AV_IClCa_junc] * CONSTANTS[AC_Acap] / (CONSTANTS[AC_F] * CONSTANTS[AC_vss]) - ALGEBRAIC[AV_JdiffCl];
    ALGEBRAIC[AV_Bcajsr] = 1.0 / (1.0 + CONSTANTS[AC_csqnmax] * CONSTANTS[AC_kmcsqn] / pow(CONSTANTS[AC_kmcsqn] + STATES[Cajsr], 2.0));
    ALGEBRAIC[AV_Bcass] = 1.0 / (1.0 + CONSTANTS[AC_BSRmax] * CONSTANTS[AC_KmBSR] / pow(CONSTANTS[AC_KmBSR] + STATES[Cass], 2.0) + CONSTANTS[AC_BSLmax] * CONSTANTS[AC_KmBSL] / pow(CONSTANTS[AC_KmBSL] + STATES[Cass], 2.0));
    ALGEBRAIC[AV_Bcai] = 1.0 / (1.0 + CONSTANTS[AC_cmdnmax] * CONSTANTS[AC_kmcmdn] / pow(CONSTANTS[AC_kmcmdn] + STATES[Cai], 2.0) + CONSTANTS[AC_trpnmax] * CONSTANTS[AC_kmtrpn] / pow(CONSTANTS[AC_kmtrpn] + STATES[Cai], 2.0));
    
    /* membrane */
    ALGEBRAIC[AV_vffrt] = STATES[V] * CONSTANTS[AC_F] * CONSTANTS[AC_F] / (CONSTANTS[AC_R] * CONSTANTS[AC_T]);
    ALGEBRAIC[AV_vfrt] = STATES[V] * CONSTANTS[AC_F] / (CONSTANTS[AC_R] * CONSTANTS[AC_T]);

    
    /* ryr */
    ALGEBRAIC[AV_fJrelp] = 1.0 / (1.0 + CONSTANTS[AC_KmCaMK] / ALGEBRAIC[AV_CaMKa]);
    ALGEBRAIC[AV_Jrel] = CONSTANTS[AC_Jrel_b] * ((1.0 - ALGEBRAIC[AV_fJrelp]) * STATES[Jrel_np] + ALGEBRAIC[AV_fJrelp] * STATES[Jrel_p]);
    ALGEBRAIC[AV_tau_rel_b] = CONSTANTS[AC_bt] / (1.0 + 0.0123 / STATES[Cajsr]);
    ALGEBRAIC[AV_tau_rel] = ((ALGEBRAIC[AV_tau_rel_b] < 0.001) ? 0.001 : ALGEBRAIC[AV_tau_rel_b]);
    ALGEBRAIC[AV_tau_relp_b] = CONSTANTS[AC_btp] / (1.0 + 0.0123 / STATES[Cajsr]);
    ALGEBRAIC[AV_tau_relp] = ((ALGEBRAIC[AV_tau_relp_b] < 0.001) ? 0.001 : ALGEBRAIC[AV_tau_relp_b]);
    
    /* *remaining* */
    ALGEBRAIC[AV_PhiCaK_i] = 1.0 * ALGEBRAIC[AV_vffrt] * (ALGEBRAIC[AV_gamma_ki] * STATES[Ki] * exp(1.0 * ALGEBRAIC[AV_vfrt]) - CONSTANTS[AC_gamma_ko] * CONSTANTS[AC_ko]) / (exp(1.0 * ALGEBRAIC[AV_vfrt]) - 1.0);
    ALGEBRAIC[AV_PhiCaK_ss] = 1.0 * ALGEBRAIC[AV_vffrt] * (ALGEBRAIC[AV_gamma_kss] * STATES[Kss] * exp(1.0 * ALGEBRAIC[AV_vfrt]) - CONSTANTS[AC_gamma_ko] * CONSTANTS[AC_ko]) / (exp(1.0 * ALGEBRAIC[AV_vfrt]) - 1.0);
    ALGEBRAIC[AV_PhiCaL_i] = 4.0 * ALGEBRAIC[AV_vffrt] * (ALGEBRAIC[AV_gamma_cai] * STATES[Cai] * exp(2.0 * ALGEBRAIC[AV_vfrt]) - CONSTANTS[AC_gamma_cao] * CONSTANTS[AC_cao]) / (exp(2.0 * ALGEBRAIC[AV_vfrt]) - 1.0);
    ALGEBRAIC[AV_PhiCaL_ss] = 4.0 * ALGEBRAIC[AV_vffrt] * (ALGEBRAIC[AV_gamma_cass] * STATES[Cass] * exp(2.0 * ALGEBRAIC[AV_vfrt]) - CONSTANTS[AC_gamma_cao] * CONSTANTS[AC_cao]) / (exp(2.0 * ALGEBRAIC[AV_vfrt]) - 1.0);
    ALGEBRAIC[AV_PhiCaNa_i] = 1.0 * ALGEBRAIC[AV_vffrt] * (ALGEBRAIC[AV_gamma_nai] * STATES[Nai] * exp(1.0 * ALGEBRAIC[AV_vfrt]) - CONSTANTS[AC_gamma_nao] * CONSTANTS[AC_nao]) / (exp(1.0 * ALGEBRAIC[AV_vfrt]) - 1.0);
    ALGEBRAIC[AV_PhiCaNa_ss] = 1.0 * ALGEBRAIC[AV_vffrt] * (ALGEBRAIC[AV_gamma_nass] * STATES[Nass] * exp(1.0 * ALGEBRAIC[AV_vfrt]) - CONSTANTS[AC_gamma_nao] * CONSTANTS[AC_nao]) / (exp(1.0 * ALGEBRAIC[AV_vfrt]) - 1.0);
    ALGEBRAIC[AV_ICab_ICab] = CONSTANTS[AC_PCab] * 4.0 * ALGEBRAIC[AV_vffrt] * (ALGEBRAIC[AV_gamma_cai] * STATES[Cai] * exp(2.0 * ALGEBRAIC[AV_vfrt]) - CONSTANTS[AC_gamma_cao] * CONSTANTS[AC_cao]) / (exp(2.0 * ALGEBRAIC[AV_vfrt]) - 1.0);
    ALGEBRAIC[AV_alpha] = 0.1161 * exp(0.299 * ALGEBRAIC[AV_vfrt]);
    ALGEBRAIC[AV_alpha_2] = 0.0578 * exp(0.971 * ALGEBRAIC[AV_vfrt]);
    ALGEBRAIC[AV_alpha_C2ToI] = 5.2e-05 * exp(1.525 * ALGEBRAIC[AV_vfrt]);
    ALGEBRAIC[AV_alpha_i] = 0.2533 * exp(0.5953 * ALGEBRAIC[AV_vfrt]);
    ALGEBRAIC[AV_beta] = 0.2442 * exp(-1.604 * ALGEBRAIC[AV_vfrt]);
    ALGEBRAIC[AV_beta_2] = 0.000349 * exp(-1.062 * ALGEBRAIC[AV_vfrt]);
    ALGEBRAIC[AV_beta_i] = 0.06525 * exp(-0.8209 * ALGEBRAIC[AV_vfrt]);
    ALGEBRAIC[AV_hca] = exp(CONSTANTS[AC_qca] * ALGEBRAIC[AV_vfrt]);
    ALGEBRAIC[AV_hna] = exp(CONSTANTS[AC_qna] * ALGEBRAIC[AV_vfrt]);
    ALGEBRAIC[AV_Knai] = CONSTANTS[AC_Knai0] * exp(CONSTANTS[AC_delta] * ALGEBRAIC[AV_vfrt] / 3.0);
    ALGEBRAIC[AV_Knao] = CONSTANTS[AC_Knao0] * exp((1.0 - CONSTANTS[AC_delta]) * ALGEBRAIC[AV_vfrt] / 3.0);
    ALGEBRAIC[AV_INab_INab] = CONSTANTS[AC_PNab] * ALGEBRAIC[AV_vffrt] * (STATES[Nai] * exp(ALGEBRAIC[AV_vfrt]) - CONSTANTS[AC_nao]) / (exp(ALGEBRAIC[AV_vfrt]) - 1.0);
    RATES[Cajsr] = ALGEBRAIC[AV_Bcajsr] * (ALGEBRAIC[AV_Jtr] - ALGEBRAIC[AV_Jrel]);
    ALGEBRAIC[AV_ICaK_i] = (1.0 - CONSTANTS[AC_ICaL_fractionSS]) * ((1.0 - ALGEBRAIC[AV_fICaLp]) * CONSTANTS[AC_PCaK] * ALGEBRAIC[AV_PhiCaK_i] * STATES[ICaL_d] * (ALGEBRAIC[AV_f] * (1.0 - STATES[ICaL_nca_i]) + STATES[ICaL_jca] * ALGEBRAIC[AV_fca] * STATES[ICaL_nca_i]) + ALGEBRAIC[AV_fICaLp] * CONSTANTS[AC_PCaKp] * ALGEBRAIC[AV_PhiCaK_i] * STATES[ICaL_d] * (ALGEBRAIC[AV_fp] * (1.0 - STATES[ICaL_nca_i]) + STATES[ICaL_jca] * ALGEBRAIC[AV_fcap] * STATES[ICaL_nca_i]));
    ALGEBRAIC[AV_ICaK_ss] = CONSTANTS[AC_ICaL_fractionSS] * ((1.0 - ALGEBRAIC[AV_fICaLp]) * CONSTANTS[AC_PCaK] * ALGEBRAIC[AV_PhiCaK_ss] * STATES[ICaL_d] * (ALGEBRAIC[AV_f] * (1.0 - STATES[ICaL_nca_ss]) + STATES[ICaL_jca] * ALGEBRAIC[AV_fca] * STATES[ICaL_nca_ss]) + ALGEBRAIC[AV_fICaLp] * CONSTANTS[AC_PCaKp] * ALGEBRAIC[AV_PhiCaK_ss] * STATES[ICaL_d] * (ALGEBRAIC[AV_fp] * (1.0 - STATES[ICaL_nca_ss]) + STATES[ICaL_jca] * ALGEBRAIC[AV_fcap] * STATES[ICaL_nca_ss]));
    ALGEBRAIC[AV_ICaL_i] = (1.0 - CONSTANTS[AC_ICaL_fractionSS]) * ((1.0 - ALGEBRAIC[AV_fICaLp]) * CONSTANTS[AC_PCa] * ALGEBRAIC[AV_PhiCaL_i] * STATES[ICaL_d] * (ALGEBRAIC[AV_f] * (1.0 - STATES[ICaL_nca_i]) + STATES[ICaL_jca] * ALGEBRAIC[AV_fca] * STATES[ICaL_nca_i]) + ALGEBRAIC[AV_fICaLp] * CONSTANTS[AC_PCap] * ALGEBRAIC[AV_PhiCaL_i] * STATES[ICaL_d] * (ALGEBRAIC[AV_fp] * (1.0 - STATES[ICaL_nca_i]) + STATES[ICaL_jca] * ALGEBRAIC[AV_fcap] * STATES[ICaL_nca_i]));
    ALGEBRAIC[AV_ICaL_ss] = CONSTANTS[AC_ICaL_fractionSS] * ((1.0 - ALGEBRAIC[AV_fICaLp]) * CONSTANTS[AC_PCa] * ALGEBRAIC[AV_PhiCaL_ss] * STATES[ICaL_d] * (ALGEBRAIC[AV_f] * (1.0 - STATES[ICaL_nca_ss]) + STATES[ICaL_jca] * ALGEBRAIC[AV_fca] * STATES[ICaL_nca_ss]) + ALGEBRAIC[AV_fICaLp] * CONSTANTS[AC_PCap] * ALGEBRAIC[AV_PhiCaL_ss] * STATES[ICaL_d] * (ALGEBRAIC[AV_fp] * (1.0 - STATES[ICaL_nca_ss]) + STATES[ICaL_jca] * ALGEBRAIC[AV_fcap] * STATES[ICaL_nca_ss]));
    ALGEBRAIC[AV_ICaNa_i] = (1.0 - CONSTANTS[AC_ICaL_fractionSS]) * ((1.0 - ALGEBRAIC[AV_fICaLp]) * CONSTANTS[AC_PCaNa] * ALGEBRAIC[AV_PhiCaNa_i] * STATES[ICaL_d] * (ALGEBRAIC[AV_f] * (1.0 - STATES[ICaL_nca_i]) + STATES[ICaL_jca] * ALGEBRAIC[AV_fca] * STATES[ICaL_nca_i]) + ALGEBRAIC[AV_fICaLp] * CONSTANTS[AC_PCaNap] * ALGEBRAIC[AV_PhiCaNa_i] * STATES[ICaL_d] * (ALGEBRAIC[AV_fp] * (1.0 - STATES[ICaL_nca_i]) + STATES[ICaL_jca] * ALGEBRAIC[AV_fcap] * STATES[ICaL_nca_i]));
    ALGEBRAIC[AV_ICaNa_ss] = CONSTANTS[AC_ICaL_fractionSS] * ((1.0 - ALGEBRAIC[AV_fICaLp]) * CONSTANTS[AC_PCaNa] * ALGEBRAIC[AV_PhiCaNa_ss] * STATES[ICaL_d] * (ALGEBRAIC[AV_f] * (1.0 - STATES[ICaL_nca_ss]) + STATES[ICaL_jca] * ALGEBRAIC[AV_fca] * STATES[ICaL_nca_ss]) + ALGEBRAIC[AV_fICaLp] * CONSTANTS[AC_PCaNap] * ALGEBRAIC[AV_PhiCaNa_ss] * STATES[ICaL_d] * (ALGEBRAIC[AV_fp] * (1.0 - STATES[ICaL_nca_ss]) + STATES[ICaL_jca] * ALGEBRAIC[AV_fcap] * STATES[ICaL_nca_ss]));
    ALGEBRAIC[AV_beta_ItoC2] = ALGEBRAIC[AV_beta_2] * ALGEBRAIC[AV_beta_i] * ALGEBRAIC[AV_alpha_C2ToI] / (ALGEBRAIC[AV_alpha_2] * ALGEBRAIC[AV_alpha_i]);
    RATES[IKr_C2] = ALGEBRAIC[AV_alpha] * STATES[IKr_C3] + CONSTANTS[AC_beta_1] * STATES[IKr_C1] - (ALGEBRAIC[AV_beta] + CONSTANTS[AC_alpha_1]) * STATES[IKr_C2];
    RATES[IKr_C3] = ALGEBRAIC[AV_beta] * STATES[IKr_C2] - ALGEBRAIC[AV_alpha] * STATES[IKr_C3];
    RATES[IKr_O] = ALGEBRAIC[AV_alpha_2] * STATES[IKr_C1] + ALGEBRAIC[AV_beta_i] * STATES[IKr_I] - (ALGEBRAIC[AV_beta_2] + ALGEBRAIC[AV_alpha_i]) * STATES[IKr_O];
    ALGEBRAIC[AV_h1_i] = 1.0 + STATES[Nai] / CONSTANTS[AC_kna3] * (1.0 + ALGEBRAIC[AV_hna]);
    ALGEBRAIC[AV_h1_ss] = 1.0 + STATES[Nass] / CONSTANTS[AC_kna3] * (1.0 + ALGEBRAIC[AV_hna]);
    ALGEBRAIC[AV_h7_i] = 1.0 + CONSTANTS[AC_nao] / CONSTANTS[AC_kna3] * (1.0 + 1.0 / ALGEBRAIC[AV_hna]);
    ALGEBRAIC[AV_h7_ss] = 1.0 + CONSTANTS[AC_nao] / CONSTANTS[AC_kna3] * (1.0 + 1.0 / ALGEBRAIC[AV_hna]);
    ALGEBRAIC[AV_a1] = CONSTANTS[AC_k1p] * pow(STATES[Nai] / ALGEBRAIC[AV_Knai],
                                               3.0) / (pow(1.0 + STATES[Nai] / ALGEBRAIC[AV_Knai], 3.0) + pow(1.0 + STATES[Ki] / CONSTANTS[AC_Kki], 2.0) - 1.0);
    ALGEBRAIC[AV_a3] = CONSTANTS[AC_k3p] * pow(CONSTANTS[AC_ko] / CONSTANTS[AC_Kko], 2.0) / (pow(1.0 + CONSTANTS[AC_nao] / ALGEBRAIC[AV_Knao], 3.0) + pow(1.0 + CONSTANTS[AC_ko] / CONSTANTS[AC_Kko], 2.0) - 1.0);
    ALGEBRAIC[AV_b2] = CONSTANTS[AC_k2m] * pow(CONSTANTS[AC_nao] / ALGEBRAIC[AV_Knao],
                                               3.0) / (pow(1.0 + CONSTANTS[AC_nao] / ALGEBRAIC[AV_Knao], 3.0) + pow(1.0 + CONSTANTS[AC_ko] / CONSTANTS[AC_Kko], 2.0) - 1.0);
    ALGEBRAIC[AV_b4] = CONSTANTS[AC_k4m] * pow(STATES[Ki] / CONSTANTS[AC_Kki], 2.0) / (pow(1.0 + STATES[Nai] / ALGEBRAIC[AV_Knai], 3.0) + pow(1.0 + STATES[Ki] / CONSTANTS[AC_Kki], 2.0) - 1.0);
    ALGEBRAIC[AV_ICaK] = ALGEBRAIC[AV_ICaK_ss] + ALGEBRAIC[AV_ICaK_i];
    ALGEBRAIC[AV_ICaL_ICaL] = ALGEBRAIC[AV_ICaL_ss] + ALGEBRAIC[AV_ICaL_i];
    ALGEBRAIC[AV_ICaNa] = ALGEBRAIC[AV_ICaNa_ss] + ALGEBRAIC[AV_ICaNa_i];
    RATES[IKr_C1] = CONSTANTS[AC_alpha_1] * STATES[IKr_C2] + ALGEBRAIC[AV_beta_2] * STATES[IKr_O] + ALGEBRAIC[AV_beta_ItoC2] * STATES[IKr_I] - (CONSTANTS[AC_beta_1] + ALGEBRAIC[AV_alpha_2] + ALGEBRAIC[AV_alpha_C2ToI]) * STATES[IKr_C1];
    RATES[IKr_I] = ALGEBRAIC[AV_alpha_C2ToI] * STATES[IKr_C1] + ALGEBRAIC[AV_alpha_i] * STATES[IKr_O] - (ALGEBRAIC[AV_beta_ItoC2] + ALGEBRAIC[AV_beta_i]) * STATES[IKr_I];
    ALGEBRAIC[AV_h2_i] = STATES[Nai] * ALGEBRAIC[AV_hna] / (CONSTANTS[AC_kna3] * ALGEBRAIC[AV_h1_i]);
    ALGEBRAIC[AV_h2_ss] = STATES[Nass] * ALGEBRAIC[AV_hna] / (CONSTANTS[AC_kna3] * ALGEBRAIC[AV_h1_ss]);
    ALGEBRAIC[AV_h3_i] = 1.0 / ALGEBRAIC[AV_h1_i];
    ALGEBRAIC[AV_h3_ss] = 1.0 / ALGEBRAIC[AV_h1_ss];
    ALGEBRAIC[AV_h8_i] = CONSTANTS[AC_nao] / (CONSTANTS[AC_kna3] * ALGEBRAIC[AV_hna] * ALGEBRAIC[AV_h7_i]);
    ALGEBRAIC[AV_h8_ss] = CONSTANTS[AC_nao] / (CONSTANTS[AC_kna3] * ALGEBRAIC[AV_hna] * ALGEBRAIC[AV_h7_ss]);
    ALGEBRAIC[AV_h9_i] = 1.0 / ALGEBRAIC[AV_h7_i];
    ALGEBRAIC[AV_h9_ss] = 1.0 / ALGEBRAIC[AV_h7_ss];
    ALGEBRAIC[AV_x1] = CONSTANTS[AC_a4] * ALGEBRAIC[AV_a1] * CONSTANTS[AC_a2] + ALGEBRAIC[AV_b2] * ALGEBRAIC[AV_b4] * ALGEBRAIC[AV_b3] + CONSTANTS[AC_a2] * ALGEBRAIC[AV_b4] * ALGEBRAIC[AV_b3] + ALGEBRAIC[AV_b3] * ALGEBRAIC[AV_a1] * CONSTANTS[AC_a2];
    ALGEBRAIC[AV_x2] = ALGEBRAIC[AV_b2] * CONSTANTS[AC_b1] * ALGEBRAIC[AV_b4] + ALGEBRAIC[AV_a1] * CONSTANTS[AC_a2] * ALGEBRAIC[AV_a3] + ALGEBRAIC[AV_a3] * CONSTANTS[AC_b1] * ALGEBRAIC[AV_b4] + CONSTANTS[AC_a2] * ALGEBRAIC[AV_a3] * ALGEBRAIC[AV_b4];
    ALGEBRAIC[AV_x3] = CONSTANTS[AC_a2] * ALGEBRAIC[AV_a3] * CONSTANTS[AC_a4] + ALGEBRAIC[AV_b3] * ALGEBRAIC[AV_b2] * CONSTANTS[AC_b1] + ALGEBRAIC[AV_b2] * CONSTANTS[AC_b1] * CONSTANTS[AC_a4] + ALGEBRAIC[AV_a3] * CONSTANTS[AC_a4] * CONSTANTS[AC_b1];
    ALGEBRAIC[AV_x4] = ALGEBRAIC[AV_b4] * ALGEBRAIC[AV_b3] * ALGEBRAIC[AV_b2] + ALGEBRAIC[AV_a3] * CONSTANTS[AC_a4] * ALGEBRAIC[AV_a1] + ALGEBRAIC[AV_b2] * CONSTANTS[AC_a4] * ALGEBRAIC[AV_a1] + ALGEBRAIC[AV_b3] * ALGEBRAIC[AV_b2] * ALGEBRAIC[AV_a1];
    RATES[Kss] = -ALGEBRAIC[AV_ICaK_ss] * CONSTANTS[AC_Acap] / (CONSTANTS[AC_F] * CONSTANTS[AC_vss]) - ALGEBRAIC[AV_JdiffK];
    ALGEBRAIC[AV_Jrel_inf_b] = -CONSTANTS[AC_a_rel] * ALGEBRAIC[AV_ICaL_ss] / 1.0 / (1.0 + pow(CONSTANTS[AC_cajsr_half] / STATES[Cajsr], 8.0));
    ALGEBRAIC[AV_Jrel_infp_b] = -CONSTANTS[AC_a_relp] * ALGEBRAIC[AV_ICaL_ss] / 1.0 / (1.0 + pow(CONSTANTS[AC_cajsr_half] / STATES[Cajsr], 8.0));
    ALGEBRAIC[AV_k3p_i] = ALGEBRAIC[AV_h9_i] * CONSTANTS[AC_wca];
    ALGEBRAIC[AV_k3p_ss] = ALGEBRAIC[AV_h9_ss] * CONSTANTS[AC_wca];
    ALGEBRAIC[AV_k3pp_i] = ALGEBRAIC[AV_h8_i] * CONSTANTS[AC_wnaca];
    ALGEBRAIC[AV_k3pp_ss] = ALGEBRAIC[AV_h8_ss] * CONSTANTS[AC_wnaca];
    ALGEBRAIC[AV_k4p_i] = ALGEBRAIC[AV_h3_i] * CONSTANTS[AC_wca] / ALGEBRAIC[AV_hca];
    ALGEBRAIC[AV_k4p_ss] = ALGEBRAIC[AV_h3_ss] * CONSTANTS[AC_wca] / ALGEBRAIC[AV_hca];
    ALGEBRAIC[AV_k4pp_i] = ALGEBRAIC[AV_h2_i] * CONSTANTS[AC_wnaca];
    ALGEBRAIC[AV_k4pp_ss] = ALGEBRAIC[AV_h2_ss] * CONSTANTS[AC_wnaca];
    ALGEBRAIC[AV_k7_i] = ALGEBRAIC[AV_h5_i] * ALGEBRAIC[AV_h2_i] * CONSTANTS[AC_wna];
    ALGEBRAIC[AV_k7_ss] = ALGEBRAIC[AV_h5_ss] * ALGEBRAIC[AV_h2_ss] * CONSTANTS[AC_wna];
    ALGEBRAIC[AV_k8_i] = ALGEBRAIC[AV_h8_i] * CONSTANTS[AC_h11_i] * CONSTANTS[AC_wna];
    ALGEBRAIC[AV_k8_ss] = ALGEBRAIC[AV_h8_ss] * CONSTANTS[AC_h11_ss] * CONSTANTS[AC_wna];
    ALGEBRAIC[AV_E1] = ALGEBRAIC[AV_x1] / (ALGEBRAIC[AV_x1] + ALGEBRAIC[AV_x2] + ALGEBRAIC[AV_x3] + ALGEBRAIC[AV_x4]);
    ALGEBRAIC[AV_E2] = ALGEBRAIC[AV_x2] / (ALGEBRAIC[AV_x1] + ALGEBRAIC[AV_x2] + ALGEBRAIC[AV_x3] + ALGEBRAIC[AV_x4]);
    ALGEBRAIC[AV_E3] = ALGEBRAIC[AV_x3] / (ALGEBRAIC[AV_x1] + ALGEBRAIC[AV_x2] + ALGEBRAIC[AV_x3] + ALGEBRAIC[AV_x4]);
    ALGEBRAIC[AV_E4] = ALGEBRAIC[AV_x4] / (ALGEBRAIC[AV_x1] + ALGEBRAIC[AV_x2] + ALGEBRAIC[AV_x3] + ALGEBRAIC[AV_x4]);
    ALGEBRAIC[AV_Jrel_inf] = ((tissueFlag == 2.0) ? ALGEBRAIC[AV_Jrel_inf_b] * 1.7 : ALGEBRAIC[AV_Jrel_inf_b]);
    ALGEBRAIC[AV_Jrel_infp] = ((tissueFlag == 2.0) ? ALGEBRAIC[AV_Jrel_infp_b] * 1.7 : ALGEBRAIC[AV_Jrel_infp_b]);
    ALGEBRAIC[AV_k3_i] = ALGEBRAIC[AV_k3p_i] + ALGEBRAIC[AV_k3pp_i];
    ALGEBRAIC[AV_k3_ss] = ALGEBRAIC[AV_k3p_ss] + ALGEBRAIC[AV_k3pp_ss];
    ALGEBRAIC[AV_k4_i] = ALGEBRAIC[AV_k4p_i] + ALGEBRAIC[AV_k4pp_i];
    ALGEBRAIC[AV_k4_ss] = ALGEBRAIC[AV_k4p_ss] + ALGEBRAIC[AV_k4pp_ss];
    ALGEBRAIC[AV_JnakK] = 2.0 * (ALGEBRAIC[AV_E4] * CONSTANTS[AC_b1] - ALGEBRAIC[AV_E3] * ALGEBRAIC[AV_a1]);
    ALGEBRAIC[AV_JnakNa] = 3.0 * (ALGEBRAIC[AV_E1] * ALGEBRAIC[AV_a3] - ALGEBRAIC[AV_E2] * ALGEBRAIC[AV_b3]);
    RATES[Jrel_np] = (ALGEBRAIC[AV_Jrel_inf] - STATES[Jrel_np]) / ALGEBRAIC[AV_tau_rel];
    RATES[Jrel_p] = (ALGEBRAIC[AV_Jrel_infp] - STATES[Jrel_p]) / ALGEBRAIC[AV_tau_relp];
    ALGEBRAIC[AV_x1_i] = CONSTANTS[AC_k2_i] * ALGEBRAIC[AV_k4_i] * (ALGEBRAIC[AV_k7_i] + ALGEBRAIC[AV_k6_i]) + CONSTANTS[AC_k5_i] * ALGEBRAIC[AV_k7_i] * (CONSTANTS[AC_k2_i] + ALGEBRAIC[AV_k3_i]);
    ALGEBRAIC[AV_x1_ss] = CONSTANTS[AC_k2_ss] * ALGEBRAIC[AV_k4_ss] * (ALGEBRAIC[AV_k7_ss] + ALGEBRAIC[AV_k6_ss]) + CONSTANTS[AC_k5_ss] * ALGEBRAIC[AV_k7_ss] * (CONSTANTS[AC_k2_ss] + ALGEBRAIC[AV_k3_ss]);
    ALGEBRAIC[AV_x2_i] = CONSTANTS[AC_k1_i] * ALGEBRAIC[AV_k7_i] * (ALGEBRAIC[AV_k4_i] + CONSTANTS[AC_k5_i]) + ALGEBRAIC[AV_k4_i] * ALGEBRAIC[AV_k6_i] * (CONSTANTS[AC_k1_i] + ALGEBRAIC[AV_k8_i]);
    ALGEBRAIC[AV_x2_ss] = CONSTANTS[AC_k1_ss] * ALGEBRAIC[AV_k7_ss] * (ALGEBRAIC[AV_k4_ss] + CONSTANTS[AC_k5_ss]) + ALGEBRAIC[AV_k4_ss] * ALGEBRAIC[AV_k6_ss] * (CONSTANTS[AC_k1_ss] + ALGEBRAIC[AV_k8_ss]);
    ALGEBRAIC[AV_x3_i] = CONSTANTS[AC_k1_i] * ALGEBRAIC[AV_k3_i] * (ALGEBRAIC[AV_k7_i] + ALGEBRAIC[AV_k6_i]) + ALGEBRAIC[AV_k8_i] * ALGEBRAIC[AV_k6_i] * (CONSTANTS[AC_k2_i] + ALGEBRAIC[AV_k3_i]);
    ALGEBRAIC[AV_x3_ss] = CONSTANTS[AC_k1_ss] * ALGEBRAIC[AV_k3_ss] * (ALGEBRAIC[AV_k7_ss] + ALGEBRAIC[AV_k6_ss]) + ALGEBRAIC[AV_k8_ss] * ALGEBRAIC[AV_k6_ss] * (CONSTANTS[AC_k2_ss] + ALGEBRAIC[AV_k3_ss]);
    ALGEBRAIC[AV_x4_i] = CONSTANTS[AC_k2_i] * ALGEBRAIC[AV_k8_i] * (ALGEBRAIC[AV_k4_i] + CONSTANTS[AC_k5_i]) + ALGEBRAIC[AV_k3_i] * CONSTANTS[AC_k5_i] * (CONSTANTS[AC_k1_i] + ALGEBRAIC[AV_k8_i]);
    ALGEBRAIC[AV_x4_ss] = CONSTANTS[AC_k2_ss] * ALGEBRAIC[AV_k8_ss] * (ALGEBRAIC[AV_k4_ss] + CONSTANTS[AC_k5_ss]) + ALGEBRAIC[AV_k3_ss] * CONSTANTS[AC_k5_ss] * (CONSTANTS[AC_k1_ss] + ALGEBRAIC[AV_k8_ss]);
    ALGEBRAIC[AV_INaK_INaK] = CONSTANTS[AC_Pnak] * (CONSTANTS[AC_zna] * ALGEBRAIC[AV_JnakNa] + CONSTANTS[AC_zk] * ALGEBRAIC[AV_JnakK]);
    ALGEBRAIC[AV_E1_i] = ALGEBRAIC[AV_x1_i] / (ALGEBRAIC[AV_x1_i] + ALGEBRAIC[AV_x2_i] + ALGEBRAIC[AV_x3_i] + ALGEBRAIC[AV_x4_i]);
    ALGEBRAIC[AV_E1_ss] = ALGEBRAIC[AV_x1_ss] / (ALGEBRAIC[AV_x1_ss] + ALGEBRAIC[AV_x2_ss] + ALGEBRAIC[AV_x3_ss] + ALGEBRAIC[AV_x4_ss]);
    ALGEBRAIC[AV_E2_i] = ALGEBRAIC[AV_x2_i] / (ALGEBRAIC[AV_x1_i] + ALGEBRAIC[AV_x2_i] + ALGEBRAIC[AV_x3_i] + ALGEBRAIC[AV_x4_i]);
    ALGEBRAIC[AV_E2_ss] = ALGEBRAIC[AV_x2_ss] / (ALGEBRAIC[AV_x1_ss] + ALGEBRAIC[AV_x2_ss] + ALGEBRAIC[AV_x3_ss] + ALGEBRAIC[AV_x4_ss]);
    ALGEBRAIC[AV_E3_i] = ALGEBRAIC[AV_x3_i] / (ALGEBRAIC[AV_x1_i] + ALGEBRAIC[AV_x2_i] + ALGEBRAIC[AV_x3_i] + ALGEBRAIC[AV_x4_i]);
    ALGEBRAIC[AV_E3_ss] = ALGEBRAIC[AV_x3_ss] / (ALGEBRAIC[AV_x1_ss] + ALGEBRAIC[AV_x2_ss] + ALGEBRAIC[AV_x3_ss] + ALGEBRAIC[AV_x4_ss]);
    ALGEBRAIC[AV_E4_i] = ALGEBRAIC[AV_x4_i] / (ALGEBRAIC[AV_x1_i] + ALGEBRAIC[AV_x2_i] + ALGEBRAIC[AV_x3_i] + ALGEBRAIC[AV_x4_i]);
    ALGEBRAIC[AV_E4_ss] = ALGEBRAIC[AV_x4_ss] / (ALGEBRAIC[AV_x1_ss] + ALGEBRAIC[AV_x2_ss] + ALGEBRAIC[AV_x3_ss] + ALGEBRAIC[AV_x4_ss]);
    RATES[Ki] = -(ALGEBRAIC[AV_Ito_Ito] + ALGEBRAIC[AV_IKr_IKr] + ALGEBRAIC[AV_IKs_IKs] + ALGEBRAIC[AV_IK1_IK1] + ALGEBRAIC[AV_IKb_IKb] + ALGEBRAIC[AV_I_katp_I_katp] + ALGEBRAIC[Istim] - 2.0 * ALGEBRAIC[AV_INaK_INaK] + ALGEBRAIC[AV_ICaK_i]) * CONSTANTS[AC_Acap] / (CONSTANTS[AC_F] * CONSTANTS[AC_vmyo]) + ALGEBRAIC[AV_JdiffK] * CONSTANTS[AC_vss] / CONSTANTS[AC_vmyo];
    ALGEBRAIC[AV_JncxCa_i] = ALGEBRAIC[AV_E2_i] * CONSTANTS[AC_k2_i] - ALGEBRAIC[AV_E1_i] * CONSTANTS[AC_k1_i];
    ALGEBRAIC[AV_JncxCa_ss] = ALGEBRAIC[AV_E2_ss] * CONSTANTS[AC_k2_ss] - ALGEBRAIC[AV_E1_ss] * CONSTANTS[AC_k1_ss];
    ALGEBRAIC[AV_JncxNa_i] = 3.0 * (ALGEBRAIC[AV_E4_i] * ALGEBRAIC[AV_k7_i] - ALGEBRAIC[AV_E1_i] * ALGEBRAIC[AV_k8_i]) + ALGEBRAIC[AV_E3_i] * ALGEBRAIC[AV_k4pp_i] - ALGEBRAIC[AV_E2_i] * ALGEBRAIC[AV_k3pp_i];
    ALGEBRAIC[AV_JncxNa_ss] = 3.0 * (ALGEBRAIC[AV_E4_ss] * ALGEBRAIC[AV_k7_ss] - ALGEBRAIC[AV_E1_ss] * ALGEBRAIC[AV_k8_ss]) + ALGEBRAIC[AV_E3_ss] * ALGEBRAIC[AV_k4pp_ss] - ALGEBRAIC[AV_E2_ss] * ALGEBRAIC[AV_k3pp_ss];
    ALGEBRAIC[AV_INaCa_i] = (1.0 - CONSTANTS[AC_INaCa_fractionSS]) * CONSTANTS[AC_Gncx] * ALGEBRAIC[AV_allo_i] * (CONSTANTS[AC_zna] * ALGEBRAIC[AV_JncxNa_i] + CONSTANTS[AC_zca] * ALGEBRAIC[AV_JncxCa_i]);
    ALGEBRAIC[AV_INaCa_ss] = CONSTANTS[AC_INaCa_fractionSS] * CONSTANTS[AC_Gncx] * ALGEBRAIC[AV_allo_ss] * (CONSTANTS[AC_zna] * ALGEBRAIC[AV_JncxNa_ss] + CONSTANTS[AC_zca] * ALGEBRAIC[AV_JncxCa_ss]);
    RATES[Cai] = ALGEBRAIC[AV_Bcai] * (-(ALGEBRAIC[AV_ICaL_i] + ALGEBRAIC[AV_IpCa_IpCa] + ALGEBRAIC[AV_ICab_ICab] - 2.0 * ALGEBRAIC[AV_INaCa_i]) * CONSTANTS[AC_Acap] / (2.0 * CONSTANTS[AC_F] * CONSTANTS[AC_vmyo]) - ALGEBRAIC[AV_Jup] * CONSTANTS[AC_vnsr] / CONSTANTS[AC_vmyo] + ALGEBRAIC[AV_Jdiff] * CONSTANTS[AC_vss] / CONSTANTS[AC_vmyo]);
    RATES[Cass] = ALGEBRAIC[AV_Bcass] * (-(ALGEBRAIC[AV_ICaL_ss] - 2.0 * ALGEBRAIC[AV_INaCa_ss]) * CONSTANTS[AC_Acap] / (2.0 * CONSTANTS[AC_F] * CONSTANTS[AC_vss]) + ALGEBRAIC[AV_Jrel] * CONSTANTS[AC_vjsr] / CONSTANTS[AC_vss] - ALGEBRAIC[AV_Jdiff]);
    RATES[Nai] = -(ALGEBRAIC[AV_INa_INa] + ALGEBRAIC[AV_INaL_INaL] + 3.0 * ALGEBRAIC[AV_INaCa_i] + ALGEBRAIC[AV_ICaNa_i] + 3.0 * ALGEBRAIC[AV_INaK_INaK] + ALGEBRAIC[AV_INab_INab]) * CONSTANTS[AC_Acap] / (CONSTANTS[AC_F] * CONSTANTS[AC_vmyo]) + ALGEBRAIC[AV_JdiffNa] * CONSTANTS[AC_vss] / CONSTANTS[AC_vmyo];
    RATES[Nass] = -(ALGEBRAIC[AV_ICaNa_ss] + 3.0 * ALGEBRAIC[AV_INaCa_ss]) * CONSTANTS[AC_Acap] / (CONSTANTS[AC_F] * CONSTANTS[AC_vss]) - ALGEBRAIC[AV_JdiffNa];
    
    
    ALGEBRAIC[Iion_cm]= (ALGEBRAIC[AV_INa_INa] + ALGEBRAIC[AV_INaL_INaL] + ALGEBRAIC[AV_Ito_Ito] + ALGEBRAIC[AV_ICaL_ICaL] + ALGEBRAIC[AV_ICaNa] + ALGEBRAIC[AV_ICaK] + ALGEBRAIC[AV_IKr_IKr] + ALGEBRAIC[AV_IKs_IKs] + ALGEBRAIC[AV_IK1_IK1] + ALGEBRAIC[AV_INaCa_i] + ALGEBRAIC[AV_INaCa_ss] + ALGEBRAIC[AV_INaK_INaK] + ALGEBRAIC[AV_INab_INab] + ALGEBRAIC[AV_IKb_IKb] + ALGEBRAIC[AV_IpCa_IpCa] + ALGEBRAIC[AV_ICab_ICab] + ALGEBRAIC[AV_IClCa] + ALGEBRAIC[AV_IClb] + ALGEBRAIC[AV_I_katp_I_katp]);
    
    ALGEBRAIC[Istim] = computeIstim(VOI, CONSTANTS);

    if (solveVmWithinODESolver)
    {
        RATES[V] = -ALGEBRAIC[Iion_cm] - ALGEBRAIC[Istim];
    }

}
