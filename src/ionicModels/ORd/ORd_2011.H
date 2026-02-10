#include "ORd_2011Names.H"
#include "stimulusIO.H"

static const char* ORdSTATES_NAMES[NUM_STATES] = {
    "membrane_V",
    "CaMK_CaMKt",
    "Na_i",
    "Na_ss",
    "K_i",
    "K_ss",
    "Ca_ss",
    "Ca_nsr",
    "Ca_jsr",
    "Ca_i",
    "INa_m",
    "INa_hf",
    "INa_hs",
    "INa_j",
    "INa_hsp",
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
    "ICaL_nca",
    "IKr_xrf",
    "IKr_xrs",
    "IKs_xs1",
    "IKs_xs2",
    "IK1_xk1",
    "ryr_Jrelnp",
    "ryr_Jrelp",
};

static const char* ORdALGEBRAIC_NAMES[NUM_ALGEBRAIC] = {
    "AV_CaMKb",
    "AV_CaMKa",
    "AV_IpCa_IpCa",
    "AV_Jdiff",
    "AV_JdiffK",
    "AV_JdiffNa",
    "AV_Jtr",
    "AV_allo_i",
    "AV_allo_ss",
    "AV_h4_i",
    "AV_h4_ss",
    "AV_hca",
    "AV_hna",
    "AV_h1_i",
    "AV_h1_ss",
    "AV_h5_i",
    "AV_h5_ss",
    "AV_h6_i",
    "AV_h6_ss",
    "AV_h7_i",
    "AV_h7_ss",
    "AV_h2_i",
    "AV_h2_ss",
    "AV_h3_i",
    "AV_h3_ss",
    "AV_h8_i",
    "AV_h8_ss",
    "AV_h9_i",
    "AV_h9_ss",
    "AV_k6_i",
    "AV_k6_ss",
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
    "AV_k3_i",
    "AV_k3_ss",
    "AV_k4_i",
    "AV_k4_ss",
    "AV_x1_i",
    "AV_x1_ss",
    "AV_x2_i",
    "AV_x2_ss",
    "AV_x3_i",
    "AV_x3_ss",
    "AV_x4_i",
    "AV_x4_ss",
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
    "AV_INaCa_i_INaCa_i",
    "AV_INaCa_ss",
    "AV_Knai",
    "AV_Knao",
    "AV_P",
    "AV_a1",
    "AV_a3",
    "AV_b2",
    "AV_b3",
    "AV_b4",
    "AV_x1",
    "AV_x2",
    "AV_x3",
    "AV_x4",
    "AV_E1",
    "AV_E2",
    "AV_E3",
    "AV_E4",
    "AV_JnakK",
    "AV_JnakNa",
    "AV_INaK_INaK",
    "AV_Jleak",
    "AV_fJupp",
    "AV_Jupnp",
    "AV_Jupp",
    "AV_Jup",
    "AV_EK",
    "AV_ENa",
    "AV_EKs",
    "AV_rk1",
    "AV_txk1",
    "AV_xk1ss",
    "AV_IK1_IK1",
    "AV_xkb",
    "AV_IKb_IKb",
    "AV_Axrf",
    "AV_rkr",
    "AV_txrf",
    "AV_txrs",
    "AV_xrss",
    "AV_Axrs",
    "AV_xr",
    "AV_IKr_IKr",
    "AV_KsCa",
    "AV_txs1",
    "AV_txs2",
    "AV_xs1ss",
    "AV_xs2ss",
    "AV_IKs_IKs",
    "AV_fINap",
    "AV_hssp",
    "AV_thf",
    "AV_ths",
    "AV_tj",
    "AV_hss",
    "AV_mss",
    "AV_thsp",
    "AV_tjp",
    "AV_tm",
    "AV_h",
    "AV_hp",
    "AV_jss",
    "AV_INa_INa",
    "AV_AiF",
    "AV_ass",
    "AV_assp",
    "AV_delta_epi",
    "AV_dti_develop",
    "AV_dti_recover",
    "AV_fItop",
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
    "AV_fINaLp",
    "AV_hLss",
    "AV_hLssp",
    "AV_mLss",
    "AV_tmL",
    "AV_INaL_INaL",
    "AV_Afcaf",
    "AV_dss",
    "AV_fICaLp",
    "AV_fss",
    "AV_km2n",
    "AV_td",
    "AV_tfcaf",
    "AV_tfcas",
    "AV_tff",
    "AV_tfs",
    "AV_Afcas",
    "AV_anca",
    "AV_fcass",
    "AV_tfcafp",
    "AV_tffp",
    "AV_f",
    "AV_fca",
    "AV_fcap",
    "AV_fp",
    "AV_Bcajsr",
    "AV_Bcass",
    "AV_Bcai",
    "AV_vffrt",
    "AV_vfrt",
    "AV_fJrelp",
    "AV_Jrel",
    "AV_tau_rel_temp",
    "AV_tau_rel",
    "AV_tau_relp_temp",
    "AV_tau_relp",
    "AV_PhiCaK",
    "AV_PhiCaL",
    "AV_PhiCaNa",
    "AV_ICab_ICab",
    "AV_INab_INab",
    "AV_ICaK",
    "AV_ICaL_ICaL",
    "AV_ICaNa",
    "AV_Jrel_inf_temp",
    "AV_Jrel_temp",
    "AV_Jrel_inf",
    "AV_Jrel_infp",
    "Istim",
    "Iion_cm",
};


static inline const Foam::HashTable<Foam::wordList>&
ORdDependencyMap()
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
ORdinitConsts(double* CONSTANTS, double* RATES, double* STATES, int tissueFlag, const Foam::dictionary& stimulus)
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

    /* environment */
    CONSTANTS[AC_celltype] = 0.0;

    /* extracellular */
    CONSTANTS[AC_cao] = 1.8;
    CONSTANTS[AC_ko] = 5.4;
    CONSTANTS[AC_nao] = 140.0;

    /* physical_constants */
    CONSTANTS[AC_F] = 96485.0;
    CONSTANTS[AC_R] = 8314.0;
    CONSTANTS[AC_T] = 310.0;
    CONSTANTS[AC_zca] = 2.0;
    CONSTANTS[AC_zk] = 1.0;
    CONSTANTS[AC_zna] = 1.0;

    /* INaCa_i */
    CONSTANTS[AC_Gncx_b] = 0.0008;
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
    CONSTANTS[AC_Gncx] = ((tissueFlag == 1) ? CONSTANTS[AC_Gncx_b] * 1.1 : ((tissueFlag == 2) ? CONSTANTS[AC_Gncx_b] * 1.4 : CONSTANTS[AC_Gncx_b]));
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


    /* reversal_potentials */
    CONSTANTS[AC_PKNa] = 0.01833;

    /* IK1 */
    CONSTANTS[AC_GK1_b] = 0.1908;
    CONSTANTS[AC_GK1] = ((tissueFlag == 1) ? CONSTANTS[AC_GK1_b] * 1.2 : ((tissueFlag == 2) ? CONSTANTS[AC_GK1_b] * 1.3 : CONSTANTS[AC_GK1_b]));

    /* IKb */
    CONSTANTS[AC_GKb_b] = 0.003;
    CONSTANTS[AC_GKb] = ((tissueFlag == 1) ? CONSTANTS[AC_GKb_b] * 0.6 : CONSTANTS[AC_GKb_b]);

    /* IKr */
    CONSTANTS[AC_GKr_b] = 0.046;
    CONSTANTS[AC_GKr] = ((tissueFlag == 1) ? CONSTANTS[AC_GKr_b] * 1.3 : ((tissueFlag == 2) ? CONSTANTS[AC_GKr_b] * 0.8 : CONSTANTS[AC_GKr_b]));

    /* IKs */
    CONSTANTS[AC_GKs_b] = 0.0034;
    CONSTANTS[AC_GKs] = ((tissueFlag == 1) ? CONSTANTS[AC_GKs_b] * 1.4 : CONSTANTS[AC_GKs_b]);

    /* INa */
    CONSTANTS[AC_Ahf] = 0.99;
    CONSTANTS[AC_GNa] = 75.0;
    CONSTANTS[AC_hssV1] = 82.9;
    CONSTANTS[AC_hssV2] = 6.086;
    CONSTANTS[AC_mssV1] = 39.57;
    CONSTANTS[AC_mssV2] = 9.871;
    CONSTANTS[AC_mtD1] = 6.765;
    CONSTANTS[AC_mtD2] = 8.552;
    CONSTANTS[AC_mtV1] = 11.64;
    CONSTANTS[AC_mtV2] = 34.77;
    CONSTANTS[AC_mtV3] = 77.42;
    CONSTANTS[AC_mtV4] = 5.955;
    CONSTANTS[AC_Ahs] = 1.0 - CONSTANTS[AC_Ahf];

    /* Ito */
    CONSTANTS[AC_Gto_b] = 0.02;
    CONSTANTS[AC_Gto] = ((tissueFlag == 1) ? CONSTANTS[AC_Gto_b] * 4.0 : ((tissueFlag == 2) ? CONSTANTS[AC_Gto_b] * 4.0 : CONSTANTS[AC_Gto_b]));

    /* INaL */
    CONSTANTS[AC_GNaL_b] = 0.0075;
    CONSTANTS[AC_thL] = 200.0;
    CONSTANTS[AC_GNaL] = ((tissueFlag == 1) ? CONSTANTS[AC_GNaL_b] * 0.6 : CONSTANTS[AC_GNaL_b]);
    CONSTANTS[AC_thLp] = 3.0 * CONSTANTS[AC_thL];

    /* ICaL */
    CONSTANTS[AC_Aff] = 0.6;
    CONSTANTS[AC_Kmn] = 0.002;
    CONSTANTS[AC_PCa_b] = 0.0001;
    CONSTANTS[AC_k2n] = 1000.0;
    CONSTANTS[AC_tjca] = 75.0;
    CONSTANTS[AC_Afs] = 1.0 - CONSTANTS[AC_Aff];
    CONSTANTS[AC_PCa] = ((tissueFlag == 1) ? CONSTANTS[AC_PCa_b] * 1.2 : ((tissueFlag == 2) ? CONSTANTS[AC_PCa_b] * 2.5 : CONSTANTS[AC_PCa_b]));
    CONSTANTS[AC_PCaK] = 0.0003574 * CONSTANTS[AC_PCa];
    CONSTANTS[AC_PCaNa] = 0.00125 * CONSTANTS[AC_PCa];
    CONSTANTS[AC_PCap] = 1.1 * CONSTANTS[AC_PCa];
    CONSTANTS[AC_PCaKp] = 0.0003574 * CONSTANTS[AC_PCap];
    CONSTANTS[AC_PCaNap] = 0.00125 * CONSTANTS[AC_PCap];

    /* ICab */
    CONSTANTS[AC_PCab] = 2.5e-08;

    /* INab */
    CONSTANTS[AC_PNab] = 3.75e-10;

    /* intracellular_ions */
    CONSTANTS[AC_BSLmax] = 1.124;
    CONSTANTS[AC_BSRmax] = 0.047;
    CONSTANTS[AC_KmBSL] = 0.0087;
    CONSTANTS[AC_KmBSR] = 0.00087;
    CONSTANTS[AC_cm] = 1.0;
    CONSTANTS[AC_cmdnmax_b] = 0.05;
    CONSTANTS[AC_csqnmax] = 10.0;
    CONSTANTS[AC_kmcmdn] = 0.00238;
    CONSTANTS[AC_kmcsqn] = 0.8;
    CONSTANTS[AC_kmtrpn] = 0.0005;
    CONSTANTS[AC_trpnmax] = 0.07;
    CONSTANTS[AC_cmdnmax] = ((tissueFlag == 1) ? CONSTANTS[AC_cmdnmax_b] * 1.3 : CONSTANTS[AC_cmdnmax_b]);



    /* ryr */
    CONSTANTS[AC_bt] = 4.75;
    CONSTANTS[AC_a_rel] = 0.5 * CONSTANTS[AC_bt];
    CONSTANTS[AC_btp] = 1.25 * CONSTANTS[AC_bt];
    CONSTANTS[AC_a_relp] = 0.5 * CONSTANTS[AC_btp];

    /* --- Initial values --- */

    STATES[membrane_V] = -87.0;
    STATES[CaMK_CaMKt] = 0.0;
    STATES[Na_i] = 7.0;
    STATES[Na_ss] = 7.0;
    STATES[K_i] = 145.0;
    STATES[K_ss] = 145.0;
    STATES[Ca_ss] = 0.0001;
    STATES[Ca_nsr] = 1.2;
    STATES[Ca_jsr] = 1.2;
    STATES[Ca_i] = 0.0001;
    STATES[INa_m] = 0.0;
    STATES[INa_hf] = 1.0;
    STATES[INa_hs] = 1.0;
    STATES[INa_j] = 1.0;
    STATES[INa_hsp] = 1.0;
    STATES[INa_jp] = 1.0;
    STATES[INaL_mL] = 0.0;
    STATES[INaL_hL] = 1.0;
    STATES[INaL_hLp] = 1.0;
    STATES[Ito_a] = 0.0;
    STATES[Ito_iF] = 1.0;
    STATES[Ito_iS] = 1.0;
    STATES[Ito_ap] = 0.0;
    STATES[Ito_iFp] = 1.0;
    STATES[Ito_iSp] = 1.0;
    STATES[ICaL_d] = 0.0;
    STATES[ICaL_ff] = 1.0;
    STATES[ICaL_fs] = 1.0;
    STATES[ICaL_fcaf] = 1.0;
    STATES[ICaL_fcas] = 1.0;
    STATES[ICaL_jca] = 1.0;
    STATES[ICaL_ffp] = 1.0;
    STATES[ICaL_fcafp] = 1.0;
    STATES[ICaL_nca] = 0.0;
    STATES[IKr_xrf] = 0.0;
    STATES[IKr_xrs] = 0.0;
    STATES[IKs_xs1] = 0.0;
    STATES[IKs_xs2] = 0.0;
    STATES[IK1_xk1] = 1.0;
    STATES[ryr_Jrelnp] = 0.0;
    STATES[ryr_Jrelp] = 0.0;
}


void
ORdcomputeVariables(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC, int tissueFlag, bool solveVmWithinODESolver)
{
/* CaMK */
    ALGEBRAIC[AV_CaMKb] = CONSTANTS[AC_CaMKo] * (1.0 - STATES[CaMK_CaMKt]) / (1.0 + CONSTANTS[AC_KmCaM] / STATES[Ca_ss]);
    ALGEBRAIC[AV_CaMKa] = ALGEBRAIC[AV_CaMKb] + STATES[CaMK_CaMKt];
    RATES[CaMK_CaMKt] = CONSTANTS[AC_aCaMK] * ALGEBRAIC[AV_CaMKb] * (ALGEBRAIC[AV_CaMKb] + STATES[CaMK_CaMKt]) - CONSTANTS[AC_bCaMK] * STATES[CaMK_CaMKt];

    /* IpCa */
    ALGEBRAIC[AV_IpCa_IpCa] = CONSTANTS[AC_GpCa] * STATES[Ca_i] / (CONSTANTS[AC_KmCap] + STATES[Ca_i]);

    /* diff */
    ALGEBRAIC[AV_Jdiff] = (STATES[Ca_ss] - STATES[Ca_i]) / 0.2;
    ALGEBRAIC[AV_JdiffK] = (STATES[K_ss] - STATES[K_i]) / 2.0;
    ALGEBRAIC[AV_JdiffNa] = (STATES[Na_ss] - STATES[Na_i]) / 2.0;


    /* trans_flux */
    ALGEBRAIC[AV_Jtr] = (STATES[Ca_nsr] - STATES[Ca_jsr]) / 100.0;

    /* INaCa_i */
    ALGEBRAIC[AV_allo_i] = 1.0 / (1.0 + pow(CONSTANTS[AC_KmCaAct] / STATES[Ca_i], 2.0));
    ALGEBRAIC[AV_allo_ss] = 1.0 / (1.0 + pow(CONSTANTS[AC_KmCaAct] / STATES[Ca_ss], 2.0));
    ALGEBRAIC[AV_h4_i] = 1.0 + STATES[Na_i] / CONSTANTS[AC_kna1] * (1.0 + STATES[Na_i] / CONSTANTS[AC_kna2]);
    ALGEBRAIC[AV_h4_ss] = 1.0 + STATES[Na_ss] / CONSTANTS[AC_kna1] * (1.0 + STATES[Na_ss] / CONSTANTS[AC_kna2]);
    ALGEBRAIC[AV_hca] = exp(CONSTANTS[AC_qca] * STATES[membrane_V] * CONSTANTS[AC_F] / (CONSTANTS[AC_R] * CONSTANTS[AC_T]));
    ALGEBRAIC[AV_hna] = exp(CONSTANTS[AC_qna] * STATES[membrane_V] * CONSTANTS[AC_F] / (CONSTANTS[AC_R] * CONSTANTS[AC_T]));
    ALGEBRAIC[AV_h1_i] = 1.0 + STATES[Na_i] / CONSTANTS[AC_kna3] * (1.0 + ALGEBRAIC[AV_hna]);
    ALGEBRAIC[AV_h1_ss] = 1.0 + STATES[Na_ss] / CONSTANTS[AC_kna3] * (1.0 + ALGEBRAIC[AV_hna]);
    ALGEBRAIC[AV_h5_i] = STATES[Na_i] * STATES[Na_i] / (ALGEBRAIC[AV_h4_i] * CONSTANTS[AC_kna1] * CONSTANTS[AC_kna2]);
    ALGEBRAIC[AV_h5_ss] = STATES[Na_ss] * STATES[Na_ss] / (ALGEBRAIC[AV_h4_ss] * CONSTANTS[AC_kna1] * CONSTANTS[AC_kna2]);
    ALGEBRAIC[AV_h6_i] = 1.0 / ALGEBRAIC[AV_h4_i];
    ALGEBRAIC[AV_h6_ss] = 1.0 / ALGEBRAIC[AV_h4_ss];
    ALGEBRAIC[AV_h7_i] = 1.0 + CONSTANTS[AC_nao] / CONSTANTS[AC_kna3] * (1.0 + 1.0 / ALGEBRAIC[AV_hna]);
    ALGEBRAIC[AV_h7_ss] = 1.0 + CONSTANTS[AC_nao] / CONSTANTS[AC_kna3] * (1.0 + 1.0 / ALGEBRAIC[AV_hna]);
    ALGEBRAIC[AV_h2_i] = STATES[Na_i] * ALGEBRAIC[AV_hna] / (CONSTANTS[AC_kna3] * ALGEBRAIC[AV_h1_i]);
    ALGEBRAIC[AV_h2_ss] = STATES[Na_ss] * ALGEBRAIC[AV_hna] / (CONSTANTS[AC_kna3] * ALGEBRAIC[AV_h1_ss]);
    ALGEBRAIC[AV_h3_i] = 1.0 / ALGEBRAIC[AV_h1_i];
    ALGEBRAIC[AV_h3_ss] = 1.0 / ALGEBRAIC[AV_h1_ss];
    ALGEBRAIC[AV_h8_i] = CONSTANTS[AC_nao] / (CONSTANTS[AC_kna3] * ALGEBRAIC[AV_hna] * ALGEBRAIC[AV_h7_i]);
    ALGEBRAIC[AV_h8_ss] = CONSTANTS[AC_nao] / (CONSTANTS[AC_kna3] * ALGEBRAIC[AV_hna] * ALGEBRAIC[AV_h7_ss]);
    ALGEBRAIC[AV_h9_i] = 1.0 / ALGEBRAIC[AV_h7_i];
    ALGEBRAIC[AV_h9_ss] = 1.0 / ALGEBRAIC[AV_h7_ss];
    ALGEBRAIC[AV_k6_i] = ALGEBRAIC[AV_h6_i] * STATES[Ca_i] * CONSTANTS[AC_kcaon];
    ALGEBRAIC[AV_k6_ss] = ALGEBRAIC[AV_h6_ss] * STATES[Ca_ss] * CONSTANTS[AC_kcaon];
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
    ALGEBRAIC[AV_k3_i] = ALGEBRAIC[AV_k3p_i] + ALGEBRAIC[AV_k3pp_i];
    ALGEBRAIC[AV_k3_ss] = ALGEBRAIC[AV_k3p_ss] + ALGEBRAIC[AV_k3pp_ss];
    ALGEBRAIC[AV_k4_i] = ALGEBRAIC[AV_k4p_i] + ALGEBRAIC[AV_k4pp_i];
    ALGEBRAIC[AV_k4_ss] = ALGEBRAIC[AV_k4p_ss] + ALGEBRAIC[AV_k4pp_ss];
    ALGEBRAIC[AV_x1_i] = CONSTANTS[AC_k2_i] * ALGEBRAIC[AV_k4_i] * (ALGEBRAIC[AV_k7_i] + ALGEBRAIC[AV_k6_i]) + CONSTANTS[AC_k5_i] * ALGEBRAIC[AV_k7_i] * (CONSTANTS[AC_k2_i] + ALGEBRAIC[AV_k3_i]);
    ALGEBRAIC[AV_x1_ss] = CONSTANTS[AC_k2_ss] * ALGEBRAIC[AV_k4_ss] * (ALGEBRAIC[AV_k7_ss] + ALGEBRAIC[AV_k6_ss]) + CONSTANTS[AC_k5_ss] * ALGEBRAIC[AV_k7_ss] * (CONSTANTS[AC_k2_ss] + ALGEBRAIC[AV_k3_ss]);
    ALGEBRAIC[AV_x2_i] = CONSTANTS[AC_k1_i] * ALGEBRAIC[AV_k7_i] * (ALGEBRAIC[AV_k4_i] + CONSTANTS[AC_k5_i]) + ALGEBRAIC[AV_k4_i] * ALGEBRAIC[AV_k6_i] * (CONSTANTS[AC_k1_i] + ALGEBRAIC[AV_k8_i]);
    ALGEBRAIC[AV_x2_ss] = CONSTANTS[AC_k1_ss] * ALGEBRAIC[AV_k7_ss] * (ALGEBRAIC[AV_k4_ss] + CONSTANTS[AC_k5_ss]) + ALGEBRAIC[AV_k4_ss] * ALGEBRAIC[AV_k6_ss] * (CONSTANTS[AC_k1_ss] + ALGEBRAIC[AV_k8_ss]);
    ALGEBRAIC[AV_x3_i] = CONSTANTS[AC_k1_i] * ALGEBRAIC[AV_k3_i] * (ALGEBRAIC[AV_k7_i] + ALGEBRAIC[AV_k6_i]) + ALGEBRAIC[AV_k8_i] * ALGEBRAIC[AV_k6_i] * (CONSTANTS[AC_k2_i] + ALGEBRAIC[AV_k3_i]);
    ALGEBRAIC[AV_x3_ss] = CONSTANTS[AC_k1_ss] * ALGEBRAIC[AV_k3_ss] * (ALGEBRAIC[AV_k7_ss] + ALGEBRAIC[AV_k6_ss]) + ALGEBRAIC[AV_k8_ss] * ALGEBRAIC[AV_k6_ss] * (CONSTANTS[AC_k2_ss] + ALGEBRAIC[AV_k3_ss]);
    ALGEBRAIC[AV_x4_i] = CONSTANTS[AC_k2_i] * ALGEBRAIC[AV_k8_i] * (ALGEBRAIC[AV_k4_i] + CONSTANTS[AC_k5_i]) + ALGEBRAIC[AV_k3_i] * CONSTANTS[AC_k5_i] * (CONSTANTS[AC_k1_i] + ALGEBRAIC[AV_k8_i]);
    ALGEBRAIC[AV_x4_ss] = CONSTANTS[AC_k2_ss] * ALGEBRAIC[AV_k8_ss] * (ALGEBRAIC[AV_k4_ss] + CONSTANTS[AC_k5_ss]) + ALGEBRAIC[AV_k3_ss] * CONSTANTS[AC_k5_ss] * (CONSTANTS[AC_k1_ss] + ALGEBRAIC[AV_k8_ss]);
    ALGEBRAIC[AV_E1_i] = ALGEBRAIC[AV_x1_i] / (ALGEBRAIC[AV_x1_i] + ALGEBRAIC[AV_x2_i] + ALGEBRAIC[AV_x3_i] + ALGEBRAIC[AV_x4_i]);
    ALGEBRAIC[AV_E1_ss] = ALGEBRAIC[AV_x1_ss] / (ALGEBRAIC[AV_x1_ss] + ALGEBRAIC[AV_x2_ss] + ALGEBRAIC[AV_x3_ss] + ALGEBRAIC[AV_x4_ss]);
    ALGEBRAIC[AV_E2_i] = ALGEBRAIC[AV_x2_i] / (ALGEBRAIC[AV_x1_i] + ALGEBRAIC[AV_x2_i] + ALGEBRAIC[AV_x3_i] + ALGEBRAIC[AV_x4_i]);
    ALGEBRAIC[AV_E2_ss] = ALGEBRAIC[AV_x2_ss] / (ALGEBRAIC[AV_x1_ss] + ALGEBRAIC[AV_x2_ss] + ALGEBRAIC[AV_x3_ss] + ALGEBRAIC[AV_x4_ss]);
    ALGEBRAIC[AV_E3_i] = ALGEBRAIC[AV_x3_i] / (ALGEBRAIC[AV_x1_i] + ALGEBRAIC[AV_x2_i] + ALGEBRAIC[AV_x3_i] + ALGEBRAIC[AV_x4_i]);
    ALGEBRAIC[AV_E3_ss] = ALGEBRAIC[AV_x3_ss] / (ALGEBRAIC[AV_x1_ss] + ALGEBRAIC[AV_x2_ss] + ALGEBRAIC[AV_x3_ss] + ALGEBRAIC[AV_x4_ss]);
    ALGEBRAIC[AV_E4_i] = ALGEBRAIC[AV_x4_i] / (ALGEBRAIC[AV_x1_i] + ALGEBRAIC[AV_x2_i] + ALGEBRAIC[AV_x3_i] + ALGEBRAIC[AV_x4_i]);
    ALGEBRAIC[AV_E4_ss] = ALGEBRAIC[AV_x4_ss] / (ALGEBRAIC[AV_x1_ss] + ALGEBRAIC[AV_x2_ss] + ALGEBRAIC[AV_x3_ss] + ALGEBRAIC[AV_x4_ss]);
    ALGEBRAIC[AV_JncxCa_i] = ALGEBRAIC[AV_E2_i] * CONSTANTS[AC_k2_i] - ALGEBRAIC[AV_E1_i] * CONSTANTS[AC_k1_i];
    ALGEBRAIC[AV_JncxCa_ss] = ALGEBRAIC[AV_E2_ss] * CONSTANTS[AC_k2_ss] - ALGEBRAIC[AV_E1_ss] * CONSTANTS[AC_k1_ss];
    ALGEBRAIC[AV_JncxNa_i] = 3.0 * (ALGEBRAIC[AV_E4_i] * ALGEBRAIC[AV_k7_i] - ALGEBRAIC[AV_E1_i] * ALGEBRAIC[AV_k8_i]) + ALGEBRAIC[AV_E3_i] * ALGEBRAIC[AV_k4pp_i] - ALGEBRAIC[AV_E2_i] * ALGEBRAIC[AV_k3pp_i];
    ALGEBRAIC[AV_JncxNa_ss] = 3.0 * (ALGEBRAIC[AV_E4_ss] * ALGEBRAIC[AV_k7_ss] - ALGEBRAIC[AV_E1_ss] * ALGEBRAIC[AV_k8_ss]) + ALGEBRAIC[AV_E3_ss] * ALGEBRAIC[AV_k4pp_ss] - ALGEBRAIC[AV_E2_ss] * ALGEBRAIC[AV_k3pp_ss];
    ALGEBRAIC[AV_INaCa_i_INaCa_i] = 0.8 * CONSTANTS[AC_Gncx] * ALGEBRAIC[AV_allo_i] * (CONSTANTS[AC_zna] * ALGEBRAIC[AV_JncxNa_i] + CONSTANTS[AC_zca] * ALGEBRAIC[AV_JncxCa_i]);
    ALGEBRAIC[AV_INaCa_ss] = 0.2 * CONSTANTS[AC_Gncx] * ALGEBRAIC[AV_allo_ss] * (CONSTANTS[AC_zna] * ALGEBRAIC[AV_JncxNa_ss] + CONSTANTS[AC_zca] * ALGEBRAIC[AV_JncxCa_ss]);

    /* INaK */
    ALGEBRAIC[AV_Knai] = CONSTANTS[AC_Knai0] * exp(CONSTANTS[AC_delta] * STATES[membrane_V] * CONSTANTS[AC_F] / (3.0 * CONSTANTS[AC_R] * CONSTANTS[AC_T]));
    ALGEBRAIC[AV_Knao] = CONSTANTS[AC_Knao0] * exp((1.0 - CONSTANTS[AC_delta]) * STATES[membrane_V] * CONSTANTS[AC_F] / (3.0 * CONSTANTS[AC_R] * CONSTANTS[AC_T]));
    ALGEBRAIC[AV_P] = CONSTANTS[AC_eP] / (1.0 + CONSTANTS[AC_H] / CONSTANTS[AC_Khp] + STATES[Na_i] / CONSTANTS[AC_Knap] + STATES[K_i] / CONSTANTS[AC_Kxkur]);
    ALGEBRAIC[AV_a1] = CONSTANTS[AC_k1p] * pow(STATES[Na_i] / ALGEBRAIC[AV_Knai],
                                               3.0) / (pow(1.0 + STATES[Na_i] / ALGEBRAIC[AV_Knai], 3.0) + pow(1.0 + STATES[K_i] / CONSTANTS[AC_Kki], 2.0) - 1.0);
    ALGEBRAIC[AV_a3] = CONSTANTS[AC_k3p] * pow(CONSTANTS[AC_ko] / CONSTANTS[AC_Kko], 2.0) / (pow(1.0 + CONSTANTS[AC_nao] / ALGEBRAIC[AV_Knao], 3.0) + pow(1.0 + CONSTANTS[AC_ko] / CONSTANTS[AC_Kko], 2.0) - 1.0);
    ALGEBRAIC[AV_b2] = CONSTANTS[AC_k2m] * pow(CONSTANTS[AC_nao] / ALGEBRAIC[AV_Knao],
                                               3.0) / (pow(1.0 + CONSTANTS[AC_nao] / ALGEBRAIC[AV_Knao], 3.0) + pow(1.0 + CONSTANTS[AC_ko] / CONSTANTS[AC_Kko], 2.0) - 1.0);
    ALGEBRAIC[AV_b3] = CONSTANTS[AC_k3m] * ALGEBRAIC[AV_P] * CONSTANTS[AC_H] / (1.0 + CONSTANTS[AC_MgATP] / CONSTANTS[AC_Kmgatp]);
    ALGEBRAIC[AV_b4] = CONSTANTS[AC_k4m] * pow(STATES[K_i] / CONSTANTS[AC_Kki], 2.0) / (pow(1.0 + STATES[Na_i] / ALGEBRAIC[AV_Knai], 3.0) + pow(1.0 + STATES[K_i] / CONSTANTS[AC_Kki], 2.0) - 1.0);
    ALGEBRAIC[AV_x1] = CONSTANTS[AC_a4] * ALGEBRAIC[AV_a1] * CONSTANTS[AC_a2] + ALGEBRAIC[AV_b2] * ALGEBRAIC[AV_b4] * ALGEBRAIC[AV_b3] + CONSTANTS[AC_a2] * ALGEBRAIC[AV_b4] * ALGEBRAIC[AV_b3] + ALGEBRAIC[AV_b3] * ALGEBRAIC[AV_a1] * CONSTANTS[AC_a2];
    ALGEBRAIC[AV_x2] = ALGEBRAIC[AV_b2] * CONSTANTS[AC_b1] * ALGEBRAIC[AV_b4] + ALGEBRAIC[AV_a1] * CONSTANTS[AC_a2] * ALGEBRAIC[AV_a3] + ALGEBRAIC[AV_a3] * CONSTANTS[AC_b1] * ALGEBRAIC[AV_b4] + CONSTANTS[AC_a2] * ALGEBRAIC[AV_a3] * ALGEBRAIC[AV_b4];
    ALGEBRAIC[AV_x3] = CONSTANTS[AC_a2] * ALGEBRAIC[AV_a3] * CONSTANTS[AC_a4] + ALGEBRAIC[AV_b3] * ALGEBRAIC[AV_b2] * CONSTANTS[AC_b1] + ALGEBRAIC[AV_b2] * CONSTANTS[AC_b1] * CONSTANTS[AC_a4] + ALGEBRAIC[AV_a3] * CONSTANTS[AC_a4] * CONSTANTS[AC_b1];
    ALGEBRAIC[AV_x4] = ALGEBRAIC[AV_b4] * ALGEBRAIC[AV_b3] * ALGEBRAIC[AV_b2] + ALGEBRAIC[AV_a3] * CONSTANTS[AC_a4] * ALGEBRAIC[AV_a1] + ALGEBRAIC[AV_b2] * CONSTANTS[AC_a4] * ALGEBRAIC[AV_a1] + ALGEBRAIC[AV_b3] * ALGEBRAIC[AV_b2] * ALGEBRAIC[AV_a1];
    ALGEBRAIC[AV_E1] = ALGEBRAIC[AV_x1] / (ALGEBRAIC[AV_x1] + ALGEBRAIC[AV_x2] + ALGEBRAIC[AV_x3] + ALGEBRAIC[AV_x4]);
    ALGEBRAIC[AV_E2] = ALGEBRAIC[AV_x2] / (ALGEBRAIC[AV_x1] + ALGEBRAIC[AV_x2] + ALGEBRAIC[AV_x3] + ALGEBRAIC[AV_x4]);
    ALGEBRAIC[AV_E3] = ALGEBRAIC[AV_x3] / (ALGEBRAIC[AV_x1] + ALGEBRAIC[AV_x2] + ALGEBRAIC[AV_x3] + ALGEBRAIC[AV_x4]);
    ALGEBRAIC[AV_E4] = ALGEBRAIC[AV_x4] / (ALGEBRAIC[AV_x1] + ALGEBRAIC[AV_x2] + ALGEBRAIC[AV_x3] + ALGEBRAIC[AV_x4]);
    ALGEBRAIC[AV_JnakK] = 2.0 * (ALGEBRAIC[AV_E4] * CONSTANTS[AC_b1] - ALGEBRAIC[AV_E3] * ALGEBRAIC[AV_a1]);
    ALGEBRAIC[AV_JnakNa] = 3.0 * (ALGEBRAIC[AV_E1] * ALGEBRAIC[AV_a3] - ALGEBRAIC[AV_E2] * ALGEBRAIC[AV_b3]);
    ALGEBRAIC[AV_INaK_INaK] = CONSTANTS[AC_Pnak] * (CONSTANTS[AC_zna] * ALGEBRAIC[AV_JnakNa] + CONSTANTS[AC_zk] * ALGEBRAIC[AV_JnakK]);

    /* SERCA */
    ALGEBRAIC[AV_Jleak] = 0.0039375 * STATES[Ca_nsr] / 15.0;
    ALGEBRAIC[AV_fJupp] = 1.0 / (1.0 + CONSTANTS[AC_KmCaMK] / ALGEBRAIC[AV_CaMKa]);
    ALGEBRAIC[AV_Jupnp] = CONSTANTS[AC_upScale] * 0.004375 * STATES[Ca_i] / (STATES[Ca_i] + 0.00092);
    ALGEBRAIC[AV_Jupp] = CONSTANTS[AC_upScale] * 2.75 * 0.004375 * STATES[Ca_i] / (STATES[Ca_i] + 0.00092 - 0.00017);
    ALGEBRAIC[AV_Jup] = (1.0 - ALGEBRAIC[AV_fJupp]) * ALGEBRAIC[AV_Jupnp] + ALGEBRAIC[AV_fJupp] * ALGEBRAIC[AV_Jupp] - ALGEBRAIC[AV_Jleak];

    /* reversal_potentials */
    ALGEBRAIC[AV_EK] = CONSTANTS[AC_R] * CONSTANTS[AC_T] / CONSTANTS[AC_F] * log(CONSTANTS[AC_ko] / STATES[K_i]);
    ALGEBRAIC[AV_ENa] = CONSTANTS[AC_R] * CONSTANTS[AC_T] / CONSTANTS[AC_F] * log(CONSTANTS[AC_nao] / STATES[Na_i]);
    ALGEBRAIC[AV_EKs] = CONSTANTS[AC_R] * CONSTANTS[AC_T] / CONSTANTS[AC_F] * log((CONSTANTS[AC_ko] + CONSTANTS[AC_PKNa] * CONSTANTS[AC_nao]) / (STATES[K_i] + CONSTANTS[AC_PKNa] * STATES[Na_i]));

    /* IK1 */
    ALGEBRAIC[AV_rk1] = 1.0 / (1.0 + exp((STATES[membrane_V] + 105.8 - 2.6 * CONSTANTS[AC_ko]) / 9.493));
    ALGEBRAIC[AV_txk1] = 122.2 / (exp(-(STATES[membrane_V] + 127.2) / 20.36) + exp((STATES[membrane_V] + 236.8) / 69.33));
    ALGEBRAIC[AV_xk1ss] = 1.0 / (1.0 + exp(-(STATES[membrane_V] + 2.5538 * CONSTANTS[AC_ko] + 144.59) / (1.5692 * CONSTANTS[AC_ko] + 3.8115)));
    RATES[IK1_xk1] = (ALGEBRAIC[AV_xk1ss] - STATES[IK1_xk1]) / ALGEBRAIC[AV_txk1];
    ALGEBRAIC[AV_IK1_IK1] = CONSTANTS[AC_GK1] * sqrt(CONSTANTS[AC_ko]) * ALGEBRAIC[AV_rk1] * STATES[IK1_xk1] * (STATES[membrane_V] - ALGEBRAIC[AV_EK]);

    /* IKb */
    ALGEBRAIC[AV_xkb] = 1.0 / (1.0 + exp(-(STATES[membrane_V] - 14.48) / 18.34));
    ALGEBRAIC[AV_IKb_IKb] = CONSTANTS[AC_GKb] * ALGEBRAIC[AV_xkb] * (STATES[membrane_V] - ALGEBRAIC[AV_EK]);

    /* IKr */
    ALGEBRAIC[AV_Axrf] = 1.0 / (1.0 + exp((STATES[membrane_V] + 54.81) / 38.21));
    ALGEBRAIC[AV_rkr] = 1.0 / (1.0 + exp((STATES[membrane_V] + 55.0) / 75.0)) * 1.0 / (1.0 + exp((STATES[membrane_V] - 10.0) / 30.0));
    ALGEBRAIC[AV_txrf] = 12.98 + 1.0 / (0.3652 * exp((STATES[membrane_V] - 31.66) / 3.869) + 4.123e-05 * exp(-(STATES[membrane_V] - 47.78) / 20.38));
    ALGEBRAIC[AV_txrs] = 1.865 + 1.0 / (0.06629 * exp((STATES[membrane_V] - 34.7) / 7.355) + 1.128e-05 * exp(-(STATES[membrane_V] - 29.74) / 25.94));
    ALGEBRAIC[AV_xrss] = 1.0 / (1.0 + exp(-(STATES[membrane_V] + 8.337) / 6.789));
    ALGEBRAIC[AV_Axrs] = 1.0 - ALGEBRAIC[AV_Axrf];
    RATES[IKr_xrf] = (ALGEBRAIC[AV_xrss] - STATES[IKr_xrf]) / ALGEBRAIC[AV_txrf];
    RATES[IKr_xrs] = (ALGEBRAIC[AV_xrss] - STATES[IKr_xrs]) / ALGEBRAIC[AV_txrs];
    ALGEBRAIC[AV_xr] = ALGEBRAIC[AV_Axrf] * STATES[IKr_xrf] + ALGEBRAIC[AV_Axrs] * STATES[IKr_xrs];
    ALGEBRAIC[AV_IKr_IKr] = CONSTANTS[AC_GKr] * sqrt(CONSTANTS[AC_ko] / 5.4) * ALGEBRAIC[AV_xr] * ALGEBRAIC[AV_rkr] * (STATES[membrane_V] - ALGEBRAIC[AV_EK]);

    /* IKs */
    ALGEBRAIC[AV_KsCa] = 1.0 + 0.6 / (1.0 + pow(3.8e-05 / STATES[Ca_i], 1.4));
    ALGEBRAIC[AV_txs1] = 817.3 + 1.0 / (0.0002326 * exp((STATES[membrane_V] + 48.28) / 17.8) + 0.001292 * exp(-(STATES[membrane_V] + 210.0) / 230.0));
    ALGEBRAIC[AV_txs2] = 1.0 / (0.01 * exp((STATES[membrane_V] - 50.0) / 20.0) + 0.0193 * exp(-(STATES[membrane_V] + 66.54) / 31.0));
    ALGEBRAIC[AV_xs1ss] = 1.0 / (1.0 + exp(-(STATES[membrane_V] + 11.6) / 8.932));
    ALGEBRAIC[AV_xs2ss] = ALGEBRAIC[AV_xs1ss];
    RATES[IKs_xs1] = (ALGEBRAIC[AV_xs1ss] - STATES[IKs_xs1]) / ALGEBRAIC[AV_txs1];
    ALGEBRAIC[AV_IKs_IKs] = CONSTANTS[AC_GKs] * ALGEBRAIC[AV_KsCa] * STATES[IKs_xs1] * STATES[IKs_xs2] * (STATES[membrane_V] - ALGEBRAIC[AV_EKs]);
    RATES[IKs_xs2] = (ALGEBRAIC[AV_xs2ss] - STATES[IKs_xs2]) / ALGEBRAIC[AV_txs2];

    /* INa */
    ALGEBRAIC[AV_fINap] = 1.0 / (1.0 + CONSTANTS[AC_KmCaMK] / ALGEBRAIC[AV_CaMKa]);
    ALGEBRAIC[AV_hssp] = 1.0 / (1.0 + exp((STATES[membrane_V] + 89.1) / 6.086));
    ALGEBRAIC[AV_thf] = 1.0 / (1.432e-05 * exp(-(STATES[membrane_V] + 1.196) / 6.285) + 6.149 * exp((STATES[membrane_V] + 0.5096) / 20.27));
    ALGEBRAIC[AV_ths] = 1.0 / (0.009794 * exp(-(STATES[membrane_V] + 17.95) / 28.05) + 0.3343 * exp((STATES[membrane_V] + 5.73) / 56.66));
    ALGEBRAIC[AV_tj] = 2.038 + 1.0 / (0.02136 * exp(-(STATES[membrane_V] + 100.6) / 8.281) + 0.3052 * exp((STATES[membrane_V] + 0.9941) / 38.45));
    ALGEBRAIC[AV_hss] = 1.0 / (1.0 + exp((STATES[membrane_V] + CONSTANTS[AC_hssV1]) / CONSTANTS[AC_hssV2]));
    ALGEBRAIC[AV_mss] = 1.0 / (1.0 + exp(-(STATES[membrane_V] + CONSTANTS[AC_mssV1]) / CONSTANTS[AC_mssV2]));
    ALGEBRAIC[AV_thsp] = 3.0 * ALGEBRAIC[AV_ths];
    ALGEBRAIC[AV_tjp] = 1.46 * ALGEBRAIC[AV_tj];
    ALGEBRAIC[AV_tm] = 1.0 / (CONSTANTS[AC_mtD1] * exp((STATES[membrane_V] + CONSTANTS[AC_mtV1]) / CONSTANTS[AC_mtV2]) + CONSTANTS[AC_mtD2] * exp(-(STATES[membrane_V] + CONSTANTS[AC_mtV3]) / CONSTANTS[AC_mtV4]));
    ALGEBRAIC[AV_h] = CONSTANTS[AC_Ahf] * STATES[INa_hf] + CONSTANTS[AC_Ahs] * STATES[INa_hs];
    ALGEBRAIC[AV_hp] = CONSTANTS[AC_Ahf] * STATES[INa_hf] + CONSTANTS[AC_Ahs] * STATES[INa_hsp];
    ALGEBRAIC[AV_jss] = ALGEBRAIC[AV_hss];
    RATES[INa_hf] = (ALGEBRAIC[AV_hss] - STATES[INa_hf]) / ALGEBRAIC[AV_thf];
    RATES[INa_hs] = (ALGEBRAIC[AV_hss] - STATES[INa_hs]) / ALGEBRAIC[AV_ths];
    RATES[INa_hsp] = (ALGEBRAIC[AV_hssp] - STATES[INa_hsp]) / ALGEBRAIC[AV_thsp];
    RATES[INa_m] = (ALGEBRAIC[AV_mss] - STATES[INa_m]) / ALGEBRAIC[AV_tm];
    ALGEBRAIC[AV_INa_INa] = CONSTANTS[AC_GNa] * (STATES[membrane_V] - ALGEBRAIC[AV_ENa]) * pow(STATES[INa_m], 3.0) * ((1.0 - ALGEBRAIC[AV_fINap]) * ALGEBRAIC[AV_h] * STATES[INa_j] + ALGEBRAIC[AV_fINap] * ALGEBRAIC[AV_hp] * STATES[INa_jp]);
    RATES[INa_j] = (ALGEBRAIC[AV_jss] - STATES[INa_j]) / ALGEBRAIC[AV_tj];
    RATES[INa_jp] = (ALGEBRAIC[AV_jss] - STATES[INa_jp]) / ALGEBRAIC[AV_tjp];

    /* Ito */
    ALGEBRAIC[AV_AiF] = 1.0 / (1.0 + exp((STATES[membrane_V] - 213.6) / 151.2));
    ALGEBRAIC[AV_ass] = 1.0 / (1.0 + exp(-(STATES[membrane_V] - 14.34) / 14.82));
    ALGEBRAIC[AV_assp] = 1.0 / (1.0 + exp(-(STATES[membrane_V] - 24.34) / 14.82));
    ALGEBRAIC[AV_delta_epi] = ((tissueFlag == 1) ? 1.0 - 0.95 / (1.0 + exp((STATES[membrane_V] + 70.0) / 5.0)) : 1.0);
    ALGEBRAIC[AV_dti_develop] = 1.354 + 0.0001 / (exp((STATES[membrane_V] - 167.4) / 15.89) + exp(-(STATES[membrane_V] - 12.23) / 0.2154));
    ALGEBRAIC[AV_dti_recover] = 1.0 - 0.5 / (1.0 + exp((STATES[membrane_V] + 70.0) / 20.0));
    ALGEBRAIC[AV_fItop] = 1.0 / (1.0 + CONSTANTS[AC_KmCaMK] / ALGEBRAIC[AV_CaMKa]);
    ALGEBRAIC[AV_iss] = 1.0 / (1.0 + exp((STATES[membrane_V] + 43.94) / 5.711));
    ALGEBRAIC[AV_ta] = 1.0515 / (1.0 / (1.2089 * (1.0 + exp(-(STATES[membrane_V] - 18.4099) / 29.3814))) + 3.5 / (1.0 + exp((STATES[membrane_V] + 100.0) / 29.3814)));
    ALGEBRAIC[AV_tiF_b] = 4.562 + 1.0 / (0.3933 * exp(-(STATES[membrane_V] + 100.0) / 100.0) + 0.08004 * exp((STATES[membrane_V] + 50.0) / 16.59));
    ALGEBRAIC[AV_tiS_b] = 23.62 + 1.0 / (0.001416 * exp(-(STATES[membrane_V] + 96.52) / 59.05) + 1.78e-08 * exp((STATES[membrane_V] + 114.1) / 8.079));
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
    ALGEBRAIC[AV_Ito_Ito] = CONSTANTS[AC_Gto] * (STATES[membrane_V] - ALGEBRAIC[AV_EK]) * ((1.0 - ALGEBRAIC[AV_fItop]) * STATES[Ito_a] * ALGEBRAIC[AV_i] + ALGEBRAIC[AV_fItop] * STATES[Ito_ap] * ALGEBRAIC[AV_ip]);
    RATES[Ito_iFp] = (ALGEBRAIC[AV_iss] - STATES[Ito_iFp]) / ALGEBRAIC[AV_tiFp];
    RATES[Ito_iSp] = (ALGEBRAIC[AV_iss] - STATES[Ito_iSp]) / ALGEBRAIC[AV_tiSp];

    /* INaL */
    ALGEBRAIC[AV_fINaLp] = 1.0 / (1.0 + CONSTANTS[AC_KmCaMK] / ALGEBRAIC[AV_CaMKa]);
    ALGEBRAIC[AV_hLss] = 1.0 / (1.0 + exp((STATES[membrane_V] + 87.61) / 7.488));
    ALGEBRAIC[AV_hLssp] = 1.0 / (1.0 + exp((STATES[membrane_V] + 93.81) / 7.488));
    ALGEBRAIC[AV_mLss] = 1.0 / (1.0 + exp(-(STATES[membrane_V] + 42.85) / 5.264));
    ALGEBRAIC[AV_tmL] = ALGEBRAIC[AV_tm];
    RATES[INaL_hL] = (ALGEBRAIC[AV_hLss] - STATES[INaL_hL]) / CONSTANTS[AC_thL];
    RATES[INaL_mL] = (ALGEBRAIC[AV_mLss] - STATES[INaL_mL]) / ALGEBRAIC[AV_tmL];
    ALGEBRAIC[AV_INaL_INaL] = CONSTANTS[AC_GNaL] * (STATES[membrane_V] - ALGEBRAIC[AV_ENa]) * STATES[INaL_mL] * ((1.0 - ALGEBRAIC[AV_fINaLp]) * STATES[INaL_hL] + ALGEBRAIC[AV_fINaLp] * STATES[INaL_hLp]);
    RATES[INaL_hLp] = (ALGEBRAIC[AV_hLssp] - STATES[INaL_hLp]) / CONSTANTS[AC_thLp];

    /* ICaL */
    ALGEBRAIC[AV_Afcaf] = 0.3 + 0.6 / (1.0 + exp((STATES[membrane_V] - 10.0) / 10.0));
    ALGEBRAIC[AV_dss] = 1.0 / (1.0 + exp(-(STATES[membrane_V] + 3.94) / 4.23));
    ALGEBRAIC[AV_fICaLp] = 1.0 / (1.0 + CONSTANTS[AC_KmCaMK] / ALGEBRAIC[AV_CaMKa]);
    ALGEBRAIC[AV_fss] = 1.0 / (1.0 + exp((STATES[membrane_V] + 19.58) / 3.696));
    ALGEBRAIC[AV_km2n] = STATES[ICaL_jca] * 1.0;
    ALGEBRAIC[AV_td] = 0.6 + 1.0 / (exp(-0.05 * (STATES[membrane_V] + 6.0)) + exp(0.09 * (STATES[membrane_V] + 14.0)));
    ALGEBRAIC[AV_tfcaf] = 7.0 + 1.0 / (0.04 * exp(-(STATES[membrane_V] - 4.0) / 7.0) + 0.04 * exp((STATES[membrane_V] - 4.0) / 7.0));
    ALGEBRAIC[AV_tfcas] = 100.0 + 1.0 / (0.00012 * exp(-STATES[membrane_V] / 3.0) + 0.00012 * exp(STATES[membrane_V] / 7.0));
    ALGEBRAIC[AV_tff] = 7.0 + 1.0 / (0.0045 * exp(-(STATES[membrane_V] + 20.0) / 10.0) + 0.0045 * exp((STATES[membrane_V] + 20.0) / 10.0));
    ALGEBRAIC[AV_tfs] = 1000.0 + 1.0 / (3.5e-05 * exp(-(STATES[membrane_V] + 5.0) / 4.0) + 3.5e-05 * exp((STATES[membrane_V] + 5.0) / 6.0));
    ALGEBRAIC[AV_Afcas] = 1.0 - ALGEBRAIC[AV_Afcaf];
    ALGEBRAIC[AV_anca] = 1.0 / (CONSTANTS[AC_k2n] / ALGEBRAIC[AV_km2n] + pow(1.0 + CONSTANTS[AC_Kmn] / STATES[Ca_ss], 4.0));
    ALGEBRAIC[AV_fcass] = ALGEBRAIC[AV_fss];
    ALGEBRAIC[AV_tfcafp] = 2.5 * ALGEBRAIC[AV_tfcaf];
    ALGEBRAIC[AV_tffp] = 2.5 * ALGEBRAIC[AV_tff];
    RATES[ICaL_d] = (ALGEBRAIC[AV_dss] - STATES[ICaL_d]) / ALGEBRAIC[AV_td];
    RATES[ICaL_ff] = (ALGEBRAIC[AV_fss] - STATES[ICaL_ff]) / ALGEBRAIC[AV_tff];
    RATES[ICaL_fs] = (ALGEBRAIC[AV_fss] - STATES[ICaL_fs]) / ALGEBRAIC[AV_tfs];
    ALGEBRAIC[AV_f] = CONSTANTS[AC_Aff] * STATES[ICaL_ff] + CONSTANTS[AC_Afs] * STATES[ICaL_fs];
    ALGEBRAIC[AV_fca] = ALGEBRAIC[AV_Afcaf] * STATES[ICaL_fcaf] + ALGEBRAIC[AV_Afcas] * STATES[ICaL_fcas];
    ALGEBRAIC[AV_fcap] = ALGEBRAIC[AV_Afcaf] * STATES[ICaL_fcafp] + ALGEBRAIC[AV_Afcas] * STATES[ICaL_fcas];
    ALGEBRAIC[AV_fp] = CONSTANTS[AC_Aff] * STATES[ICaL_ffp] + CONSTANTS[AC_Afs] * STATES[ICaL_fs];
    RATES[ICaL_fcaf] = (ALGEBRAIC[AV_fcass] - STATES[ICaL_fcaf]) / ALGEBRAIC[AV_tfcaf];
    RATES[ICaL_fcafp] = (ALGEBRAIC[AV_fcass] - STATES[ICaL_fcafp]) / ALGEBRAIC[AV_tfcafp];
    RATES[ICaL_fcas] = (ALGEBRAIC[AV_fcass] - STATES[ICaL_fcas]) / ALGEBRAIC[AV_tfcas];
    RATES[ICaL_ffp] = (ALGEBRAIC[AV_fss] - STATES[ICaL_ffp]) / ALGEBRAIC[AV_tffp];
    RATES[ICaL_jca] = (ALGEBRAIC[AV_fcass] - STATES[ICaL_jca]) / CONSTANTS[AC_tjca];
    RATES[ICaL_nca] = ALGEBRAIC[AV_anca] * CONSTANTS[AC_k2n] - STATES[ICaL_nca] * ALGEBRAIC[AV_km2n];

    /* intracellular_ions */
    RATES[Ca_nsr] = ALGEBRAIC[AV_Jup] - ALGEBRAIC[AV_Jtr] * CONSTANTS[AC_vjsr] / CONSTANTS[AC_vnsr];
    ALGEBRAIC[AV_Bcajsr] = 1.0 / (1.0 + CONSTANTS[AC_csqnmax] * CONSTANTS[AC_kmcsqn] / pow(CONSTANTS[AC_kmcsqn] + STATES[Ca_jsr], 2.0));
    ALGEBRAIC[AV_Bcass] = 1.0 / (1.0 + CONSTANTS[AC_BSRmax] * CONSTANTS[AC_KmBSR] / pow(CONSTANTS[AC_KmBSR] + STATES[Ca_ss], 2.0) + CONSTANTS[AC_BSLmax] * CONSTANTS[AC_KmBSL] / pow(CONSTANTS[AC_KmBSL] + STATES[Ca_ss], 2.0));
    ALGEBRAIC[AV_Bcai] = 1.0 / (1.0 + CONSTANTS[AC_cmdnmax] * CONSTANTS[AC_kmcmdn] / pow(CONSTANTS[AC_kmcmdn] + STATES[Ca_i], 2.0) + CONSTANTS[AC_trpnmax] * CONSTANTS[AC_kmtrpn] / pow(CONSTANTS[AC_kmtrpn] + STATES[Ca_i], 2.0));

    /* membrane */
    ALGEBRAIC[AV_vffrt] = STATES[membrane_V] * CONSTANTS[AC_F] * CONSTANTS[AC_F] / (CONSTANTS[AC_R] * CONSTANTS[AC_T]);
    ALGEBRAIC[AV_vfrt] = STATES[membrane_V] * CONSTANTS[AC_F] / (CONSTANTS[AC_R] * CONSTANTS[AC_T]);

    /* ryr */
    ALGEBRAIC[AV_fJrelp] = 1.0 / (1.0 + CONSTANTS[AC_KmCaMK] / ALGEBRAIC[AV_CaMKa]);
    ALGEBRAIC[AV_Jrel] = (1.0 - ALGEBRAIC[AV_fJrelp]) * STATES[ryr_Jrelnp] + ALGEBRAIC[AV_fJrelp] * STATES[ryr_Jrelp];
    ALGEBRAIC[AV_tau_rel_temp] = CONSTANTS[AC_bt] / (1.0 + 0.0123 / STATES[Ca_jsr]);
    ALGEBRAIC[AV_tau_rel] = ((ALGEBRAIC[AV_tau_rel_temp] < 0.001) ? 0.001 : ALGEBRAIC[AV_tau_rel_temp]);
    ALGEBRAIC[AV_tau_relp_temp] = CONSTANTS[AC_btp] / (1.0 + 0.0123 / STATES[Ca_jsr]);
    ALGEBRAIC[AV_tau_relp] = ((ALGEBRAIC[AV_tau_relp_temp] < 0.001) ? 0.001 : ALGEBRAIC[AV_tau_relp_temp]);

    /* *remaining* */
    ALGEBRAIC[AV_PhiCaK] = 1.0 * ALGEBRAIC[AV_vffrt] * (0.75 * STATES[K_ss] * exp(1.0 * ALGEBRAIC[AV_vfrt]) - 0.75 * CONSTANTS[AC_ko]) / (exp(1.0 * ALGEBRAIC[AV_vfrt]) - 1.0);
    ALGEBRAIC[AV_PhiCaL] = 4.0 * ALGEBRAIC[AV_vffrt] * (STATES[Ca_ss] * exp(2.0 * ALGEBRAIC[AV_vfrt]) - 0.341 * CONSTANTS[AC_cao]) / (exp(2.0 * ALGEBRAIC[AV_vfrt]) - 1.0);
    ALGEBRAIC[AV_PhiCaNa] = 1.0 * ALGEBRAIC[AV_vffrt] * (0.75 * STATES[Na_ss] * exp(1.0 * ALGEBRAIC[AV_vfrt]) - 0.75 * CONSTANTS[AC_nao]) / (exp(1.0 * ALGEBRAIC[AV_vfrt]) - 1.0);
    ALGEBRAIC[AV_ICab_ICab] = CONSTANTS[AC_PCab] * 4.0 * ALGEBRAIC[AV_vffrt] * (STATES[Ca_i] * exp(2.0 * ALGEBRAIC[AV_vfrt]) - 0.341 * CONSTANTS[AC_cao]) / (exp(2.0 * ALGEBRAIC[AV_vfrt]) - 1.0);
    ALGEBRAIC[AV_INab_INab] = CONSTANTS[AC_PNab] * ALGEBRAIC[AV_vffrt] * (STATES[Na_i] * exp(ALGEBRAIC[AV_vfrt]) - CONSTANTS[AC_nao]) / (exp(ALGEBRAIC[AV_vfrt]) - 1.0);
    RATES[Ca_jsr] = ALGEBRAIC[AV_Bcajsr] * (ALGEBRAIC[AV_Jtr] - ALGEBRAIC[AV_Jrel]);
    RATES[K_i] = -(ALGEBRAIC[AV_Ito_Ito] + ALGEBRAIC[AV_IKr_IKr] + ALGEBRAIC[AV_IKs_IKs] + ALGEBRAIC[AV_IK1_IK1] + ALGEBRAIC[AV_IKb_IKb] + ALGEBRAIC[Istim] - 2.0 * ALGEBRAIC[AV_INaK_INaK]) * CONSTANTS[AC_cm] * CONSTANTS[AC_Acap] / (CONSTANTS[AC_F] * CONSTANTS[AC_vmyo]) + ALGEBRAIC[AV_JdiffK] * CONSTANTS[AC_vss] / CONSTANTS[AC_vmyo];
    ALGEBRAIC[AV_ICaK] = (1.0 - ALGEBRAIC[AV_fICaLp]) * CONSTANTS[AC_PCaK] * ALGEBRAIC[AV_PhiCaK] * STATES[ICaL_d] * (ALGEBRAIC[AV_f] * (1.0 - STATES[ICaL_nca]) + STATES[ICaL_jca] * ALGEBRAIC[AV_fca] * STATES[ICaL_nca]) + ALGEBRAIC[AV_fICaLp] * CONSTANTS[AC_PCaKp] * ALGEBRAIC[AV_PhiCaK] * STATES[ICaL_d] * (ALGEBRAIC[AV_fp] * (1.0 - STATES[ICaL_nca]) + STATES[ICaL_jca] * ALGEBRAIC[AV_fcap] * STATES[ICaL_nca]);
    ALGEBRAIC[AV_ICaL_ICaL] = (1.0 - ALGEBRAIC[AV_fICaLp]) * CONSTANTS[AC_PCa] * ALGEBRAIC[AV_PhiCaL] * STATES[ICaL_d] * (ALGEBRAIC[AV_f] * (1.0 - STATES[ICaL_nca]) + STATES[ICaL_jca] * ALGEBRAIC[AV_fca] * STATES[ICaL_nca]) + ALGEBRAIC[AV_fICaLp] * CONSTANTS[AC_PCap] * ALGEBRAIC[AV_PhiCaL] * STATES[ICaL_d] * (ALGEBRAIC[AV_fp] * (1.0 - STATES[ICaL_nca]) + STATES[ICaL_jca] * ALGEBRAIC[AV_fcap] * STATES[ICaL_nca]);
    ALGEBRAIC[AV_ICaNa] = (1.0 - ALGEBRAIC[AV_fICaLp]) * CONSTANTS[AC_PCaNa] * ALGEBRAIC[AV_PhiCaNa] * STATES[ICaL_d] * (ALGEBRAIC[AV_f] * (1.0 - STATES[ICaL_nca]) + STATES[ICaL_jca] * ALGEBRAIC[AV_fca] * STATES[ICaL_nca]) + ALGEBRAIC[AV_fICaLp] * CONSTANTS[AC_PCaNap] * ALGEBRAIC[AV_PhiCaNa] * STATES[ICaL_d] * (ALGEBRAIC[AV_fp] * (1.0 - STATES[ICaL_nca]) + STATES[ICaL_jca] * ALGEBRAIC[AV_fcap] * STATES[ICaL_nca]);
    RATES[Ca_i] = ALGEBRAIC[AV_Bcai] * (-(ALGEBRAIC[AV_IpCa_IpCa] + ALGEBRAIC[AV_ICab_ICab] - 2.0 * ALGEBRAIC[AV_INaCa_i_INaCa_i]) * CONSTANTS[AC_cm] * CONSTANTS[AC_Acap] / (2.0 * CONSTANTS[AC_F] * CONSTANTS[AC_vmyo]) - ALGEBRAIC[AV_Jup] * CONSTANTS[AC_vnsr] / CONSTANTS[AC_vmyo] + ALGEBRAIC[AV_Jdiff] * CONSTANTS[AC_vss] / CONSTANTS[AC_vmyo]);
    RATES[Na_i] = -(ALGEBRAIC[AV_INa_INa] + ALGEBRAIC[AV_INaL_INaL] + 3.0 * ALGEBRAIC[AV_INaCa_i_INaCa_i] + 3.0 * ALGEBRAIC[AV_INaK_INaK] + ALGEBRAIC[AV_INab_INab]) * CONSTANTS[AC_Acap] * CONSTANTS[AC_cm] / (CONSTANTS[AC_F] * CONSTANTS[AC_vmyo]) + ALGEBRAIC[AV_JdiffNa] * CONSTANTS[AC_vss] / CONSTANTS[AC_vmyo];
    RATES[Ca_ss] = ALGEBRAIC[AV_Bcass] * (-(ALGEBRAIC[AV_ICaL_ICaL] - 2.0 * ALGEBRAIC[AV_INaCa_ss]) * CONSTANTS[AC_cm] * CONSTANTS[AC_Acap] / (2.0 * CONSTANTS[AC_F] * CONSTANTS[AC_vss]) + ALGEBRAIC[AV_Jrel] * CONSTANTS[AC_vjsr] / CONSTANTS[AC_vss] - ALGEBRAIC[AV_Jdiff]);
    RATES[K_ss] = -ALGEBRAIC[AV_ICaK] * CONSTANTS[AC_cm] * CONSTANTS[AC_Acap] / (CONSTANTS[AC_F] * CONSTANTS[AC_vss]) - ALGEBRAIC[AV_JdiffK];
    RATES[Na_ss] = -(ALGEBRAIC[AV_ICaNa] + 3.0 * ALGEBRAIC[AV_INaCa_ss]) * CONSTANTS[AC_cm] * CONSTANTS[AC_Acap] / (CONSTANTS[AC_F] * CONSTANTS[AC_vss]) - ALGEBRAIC[AV_JdiffNa];
    ALGEBRAIC[AV_Jrel_inf_temp] = CONSTANTS[AC_a_rel] * -ALGEBRAIC[AV_ICaL_ICaL] / (1.0 + 1.0 * pow(1.5 / STATES[Ca_jsr], 8.0));
    ALGEBRAIC[AV_Jrel_temp] = CONSTANTS[AC_a_relp] * -ALGEBRAIC[AV_ICaL_ICaL] / (1.0 + pow(1.5 / STATES[Ca_jsr], 8.0));
    ALGEBRAIC[AV_Jrel_inf] = ((tissueFlag == 2) ? ALGEBRAIC[AV_Jrel_inf_temp] * 1.7 : ALGEBRAIC[AV_Jrel_inf_temp]);
    ALGEBRAIC[AV_Jrel_infp] = ((tissueFlag == 2) ? ALGEBRAIC[AV_Jrel_temp] * 1.7 : ALGEBRAIC[AV_Jrel_temp]);
    RATES[ryr_Jrelnp] = (ALGEBRAIC[AV_Jrel_inf] - STATES[ryr_Jrelnp]) / ALGEBRAIC[AV_tau_rel];
    RATES[ryr_Jrelp] = (ALGEBRAIC[AV_Jrel_infp] - STATES[ryr_Jrelp]) / ALGEBRAIC[AV_tau_relp];


    ALGEBRAIC[Iion_cm] = ALGEBRAIC[AV_INa_INa] + ALGEBRAIC[AV_INaL_INaL] + ALGEBRAIC[AV_Ito_Ito] + ALGEBRAIC[AV_ICaL_ICaL] + ALGEBRAIC[AV_ICaNa] + ALGEBRAIC[AV_ICaK] + ALGEBRAIC[AV_IKr_IKr] + ALGEBRAIC[AV_IKs_IKs] + ALGEBRAIC[AV_IK1_IK1] + ALGEBRAIC[AV_INaCa_i_INaCa_i] + ALGEBRAIC[AV_INaCa_ss] + ALGEBRAIC[AV_INaK_INaK] + ALGEBRAIC[AV_INab_INab] + ALGEBRAIC[AV_IKb_IKb] + ALGEBRAIC[AV_IpCa_IpCa] + ALGEBRAIC[AV_ICab_ICab];


    ALGEBRAIC[Istim] = computeIstim(VOI, CONSTANTS);

    if (solveVmWithinODESolver)
    {
        RATES[membrane_V] = -ALGEBRAIC[Iion_cm] - ALGEBRAIC[Istim];
    }

}
