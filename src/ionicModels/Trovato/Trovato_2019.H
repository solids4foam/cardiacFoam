#include "Trovato_2019Names.H"
#include "stimulusIO.H"

static const char* TrovatoSTATES_NAMES[NUM_STATES] = {
    "membrane_v",
    "CaMKt",
    "Nai",
    "Nasl",
    "Nass",
    "Ki",
    "Kss",
    "Ksl",
    "Cai",
    "Cass",
    "Casl",
    "Cansr",
    "Cajsr",
    "Cacsr",
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
    "Ito_i1",
    "Ito_i2",
    "ICaL_d",
    "ICaL_ff",
    "ICaL_fs",
    "ICaL_fcaf",
    "ICaL_fcas",
    "ICaL_jca",
    "ICaL_ffp",
    "ICaL_fcafp",
    "ICaL_nca",
    "ICaT_b",
    "ICaT_g",
    "IKr_xrf",
    "IKr_xrs",
    "IKs_xs1",
    "IKs_xs2",
    "If_y",
    "IK1_xk1",
    "ryr_Jrel1",
    "ryr_Jrel2",
    "IP3_u",
};

static const char* TrovatoALGEBRAIC_NAMES[NUM_ALGEBRAIC] = {
    "AV_CaMKb",
    "AV_CaMKa",
    "AV_POip3",
    "AV_Jip3",
    "AV_IpCa_IpCa",
    "AV_Jdiff",
    "AV_JdiffK",
    "AV_JdiffNa",
    "AV_Jgap",
    "AV_JgapK",
    "AV_JgapNa",
    "AV_Jtr1",
    "AV_Jtr2",
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
    "AV_dkmplb",
    "AV_dqupcamk",
    "AV_Jup1",
    "AV_Jup2",
    "AV_ECa",
    "AV_EK",
    "AV_ENa",
    "AV_EKs",
    "AV_bss",
    "AV_gss",
    "AV_taub",
    "AV_taug",
    "AV_ICaT_ICaT",
    "AV_rk1",
    "AV_txk1",
    "AV_xk1ss",
    "AV_IK1_IK1",
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
    "AV_IKs_IKs",
    "AV_xs2ss",
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
    "AV_tauy",
    "AV_yss",
    "AV_IfK",
    "AV_IfNa",
    "AV_If_If",
    "AV_asus",
    "AV_Isus_Isus",
    "AV_ass",
    "AV_iss",
    "AV_taua",
    "AV_tauif",
    "AV_tauis",
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
    "AV_Bcacsr",
    "AV_Bcai",
    "AV_Bcajsr",
    "AV_Bcasl",
    "AV_Bcass",
    "AV_vffrt",
    "AV_vfrt",
    "AV_REL2",
    "AV_ireltau",
    "AV_ireltau2",
    "AV_irelss2",
    "AV_PhiCaK",
    "AV_PhiCaL",
    "AV_PhiCaNa",
    "AV_ICab_ICab",
    "AV_INab_INab",
    "AV_ICaK",
    "AV_ICaL_ICaL",
    "AV_ICaNa",
    "AV_REL",
    "AV_irelss",
    "Istim",
    "Iion_cm",
};


static inline const Foam::HashTable<Foam::wordList>&
TrovatoDependencyMap()
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
TrovatoinitConsts(double* CONSTANTS, double* RATES, double* STATES, int tissueFlag, const Foam::dictionary& stimulus)
{
/* CaMK */
    CONSTANTS[AC_CaMKo] = 0.05;
    CONSTANTS[AC_KmCaM] = 0.0015;
    CONSTANTS[AC_KmCaMK] = 0.15;
    CONSTANTS[AC_aCaMK] = 0.05;
    CONSTANTS[AC_bCaMK] = 0.00068;

    /* IP3 */
    CONSTANTS[AC_IP3_IP3] = 0.0001;
    CONSTANTS[AC_k0] = 96000.0;
    CONSTANTS[AC_k0a] = 9.6;
    CONSTANTS[AC_k1] = 150000.0;
    CONSTANTS[AC_k1a] = 16.5;
    CONSTANTS[AC_k2] = 1800.0;
    CONSTANTS[AC_k2a] = 0.21;
    CONSTANTS[AC_tauip3] = 3.7;

    /* IpCa */
    CONSTANTS[AC_GpCa] = 0.0005;
    CONSTANTS[AC_KmCap] = 0.0005;

    /* cell_geometry */
    CONSTANTS[AC_L] = 0.0164;
    CONSTANTS[AC_greekpi] = 3.14159265000000021e+00;
    CONSTANTS[AC_rad] = 0.00175;
    CONSTANTS[AC_Ageo] = 2.0 * CONSTANTS[AC_greekpi] * CONSTANTS[AC_rad] * CONSTANTS[AC_rad] + 2.0 * CONSTANTS[AC_greekpi] * CONSTANTS[AC_rad] * CONSTANTS[AC_L];
    CONSTANTS[AC_vcell] = 1000.0 * 3.14159265000000021e+00 * CONSTANTS[AC_rad] * CONSTANTS[AC_rad] * CONSTANTS[AC_L];
    CONSTANTS[AC_Acap] = 2.0 * CONSTANTS[AC_Ageo];
    CONSTANTS[AC_vcsr] = 0.008 * CONSTANTS[AC_vcell];
    CONSTANTS[AC_vjsr] = 0.002 * CONSTANTS[AC_vcell];
    CONSTANTS[AC_vmyo] = 0.6 * CONSTANTS[AC_vcell];
    CONSTANTS[AC_vnsr] = 0.04 * CONSTANTS[AC_vcell];
    CONSTANTS[AC_vsl] = 0.15 * CONSTANTS[AC_vcell];
    CONSTANTS[AC_vss] = 0.02 * CONSTANTS[AC_vcell];

    /* diff */
    CONSTANTS[AC_gaptau] = 12.0;
    CONSTANTS[AC_sstau] = 0.2;

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
    CONSTANTS[AC_Gncx] = 9.57089999999999959e-04;
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
    CONSTANTS[AC_Pnak] = 32.4872;
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
    CONSTANTS[AC_a2] = CONSTANTS[AC_k2p];
    CONSTANTS[AC_a4] = CONSTANTS[AC_k4p] * CONSTANTS[AC_MgATP] / CONSTANTS[AC_Kmgatp] / (1.0 + CONSTANTS[AC_MgATP] / CONSTANTS[AC_Kmgatp]);
    CONSTANTS[AC_b1] = CONSTANTS[AC_k1m] * CONSTANTS[AC_MgADP];

    /* SERCA */
    CONSTANTS[AC_dkmplbbar] = 0.00017;
    CONSTANTS[AC_dqupcamkbar] = 0.75;
    CONSTANTS[AC_kmup] = 0.00028;
    CONSTANTS[AC_nsrbar] = 15.0;

    /* reversal_potentials */
    CONSTANTS[AC_PKNa] = 0.01833;

    /* ICaT */
    CONSTANTS[AC_GCaT] = 0.0754;

    /* IK1 */
    CONSTANTS[AC_GK1] = 0.0455;

    /* IKr */
    CONSTANTS[AC_GKr] = 0.0342;

    /* IKs */
    CONSTANTS[AC_GKs] = 0.0029;

    /* INa */
    CONSTANTS[AC_Ahf] = 0.99;
    CONSTANTS[AC_GNa] = 39.4572;
    CONSTANTS[AC_hssV1] = 78.5;
    CONSTANTS[AC_hssV2] = 6.22;
    CONSTANTS[AC_mssV1] = 48.4264;
    CONSTANTS[AC_mssV2] = 7.5653;
    CONSTANTS[AC_mtD1] = 6.765;
    CONSTANTS[AC_mtD2] = 8.552;
    CONSTANTS[AC_mtV1] = 11.64;
    CONSTANTS[AC_mtV2] = 34.77;
    CONSTANTS[AC_mtV3] = 77.42;
    CONSTANTS[AC_mtV4] = 5.955;
    CONSTANTS[AC_Ahs] = 1.0 - CONSTANTS[AC_Ahf];

    /* If */
    CONSTANTS[AC_GfK] = 0.0232;
    CONSTANTS[AC_GfNa] = 0.0116;

    /* Isus */
    CONSTANTS[AC_Gsus] = 0.0301;

    /* Ito */
    CONSTANTS[AC_Gto] = 0.192;

    /* INaL */
    CONSTANTS[AC_GNaL] = 0.0189;
    CONSTANTS[AC_thL] = 200.0;
    CONSTANTS[AC_thLp] = 3.0 * CONSTANTS[AC_thL];

    /* ICaL */
    CONSTANTS[AC_Aff] = 0.6;
    CONSTANTS[AC_Kmn] = 0.002;
    CONSTANTS[AC_PCa] = 7.76769999999999951e-05;
    CONSTANTS[AC_k2n] = 1000.0;
    CONSTANTS[AC_tjca] = 75.0;
    CONSTANTS[AC_Afs] = 1.0 - CONSTANTS[AC_Aff];
    CONSTANTS[AC_PCaK] = 0.0003574 * CONSTANTS[AC_PCa];
    CONSTANTS[AC_PCaNa] = 0.00125 * CONSTANTS[AC_PCa];
    CONSTANTS[AC_PCap] = 1.1 * CONSTANTS[AC_PCa];
    CONSTANTS[AC_PCaKp] = 0.0003574 * CONSTANTS[AC_PCap];
    CONSTANTS[AC_PCaNap] = 0.00125 * CONSTANTS[AC_PCap];

    /* ICab */
    CONSTANTS[AC_PCab] = 2.5e-08;

    /* INab */
    CONSTANTS[AC_PNab] = 9.375e-10;

    /* intracellular_ions */
    CONSTANTS[AC_BSLmax] = 0.4777;
    CONSTANTS[AC_BSRmax] = 0.019975;
    CONSTANTS[AC_KmBSL] = 0.0087;
    CONSTANTS[AC_KmBSR] = 0.00087;
    CONSTANTS[AC_cm] = 1.0;
    CONSTANTS[AC_cmdnmax] = 0.1125;
    CONSTANTS[AC_cmdnmaxsl] = 0.0125;
    CONSTANTS[AC_csqnmax] = 2.88;
    CONSTANTS[AC_csqnmaxsl] = 1.2;
    CONSTANTS[AC_kmcmdn] = 0.00238;
    CONSTANTS[AC_kmcsqn] = 0.8;
    CONSTANTS[AC_kmtrpn] = 0.0005;
    CONSTANTS[AC_trpnmax] = 0.0315;
    CONSTANTS[AC_trpnmaxsl] = 0.0035;



    /* --- Initial values --- */

    STATES[membrane_v] = -8.66814002878592049e+01;
    STATES[CaMKt] =  5.05983330678751002e-03;
    STATES[Nai] =  8.23183964616931974e+00;
    STATES[Nasl] =  8.23153516580562084e+00;
    STATES[Nass] =  8.23154325237267948e+00;
    STATES[Ki] =  1.43767359809131989e+02;
    STATES[Kss] =  1.43767768218104010e+02;
    STATES[Ksl] =  1.43767769906215989e+02;
    STATES[Cai] =  4.36004404734282013e-05;
    STATES[Cass] =  1.01777993438817996e-04;
    STATES[Casl] =  1.02004317781146999e-04;
    STATES[Cansr] =  1.26350902016857991e+00;
    STATES[Cajsr] =  1.24811940209535011e+00;
    STATES[Cacsr] =  1.26516959198517998e+00;
    STATES[INa_m] =  6.32661703915808008e-03;
    STATES[INa_hf] =  7.88611739889677033e-01;
    STATES[INa_hs] =  7.88545979951331022e-01;
    STATES[INa_j] =  7.90474358603665994e-01;
    STATES[INa_hsp] =  5.79693514309867042e-01;
    STATES[INa_jp] =  7.90947058236416978e-01;
    STATES[INaL_mL] =  2.41925773627233007e-04;
    STATES[INaL_hL] =  4.63574582508217981e-01;
    STATES[INaL_hLp] =  2.40216198686475008e-01;
    STATES[Ito_a] =  2.72851144435704018e-04;
    STATES[Ito_i1] =  6.49604795721571038e-01;
    STATES[Ito_i2] =  9.89965695822495051e-01;
    STATES[ICaL_d] =  6.97735089296891969e-09;
    STATES[ICaL_ff] =  9.99999968230738001e-01;
    STATES[ICaL_fs] =  9.26692153319136014e-01;
    STATES[ICaL_fcaf] =  9.99999968195730005e-01;
    STATES[ICaL_fcas] =  9.99999905741935979e-01;
    STATES[ICaL_jca] =  9.99978907334661993e-01;
    STATES[ICaL_ffp] =  9.99999968365902991e-01;
    STATES[ICaL_fcafp] =  9.99999968278239004e-01;
    STATES[ICaL_nca] =  5.47252500964926034e-03;
    STATES[ICaT_b] =  3.04250912559619012e-04;
    STATES[ICaT_g] =  9.94214357917907021e-01;
    STATES[IKr_xrf] =  3.31691184084271975e-04;
    STATES[IKr_xrs] =  5.68716473334160977e-01;
    STATES[IKs_xs1] =  1.91165248085394007e-01;
    STATES[IKs_xs2] =  2.22677365291218995e-04;
    STATES[If_y] =  2.33119011214907995e-01;
    STATES[IK1_xk1] =  9.97084813729908981e-01;
    STATES[ryr_Jrel1] =  1.08240945806961999e-04;
    STATES[ryr_Jrel2] =  1.25045800437316998e-69;
    STATES[IP3_u] =  4.66236137183557997e-01;
}


void
TrovatocomputeVariables(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC, int tissueFlag, bool solveVmWithinODESolver)
{
/* CaMK */
    ALGEBRAIC[AV_CaMKb] = CONSTANTS[AC_CaMKo] * (1.0 - STATES[CaMKt]) / (1.0 + CONSTANTS[AC_KmCaM] / STATES[Cass]);
    ALGEBRAIC[AV_CaMKa] = ALGEBRAIC[AV_CaMKb] + STATES[CaMKt];
    RATES[CaMKt] = CONSTANTS[AC_aCaMK] * ALGEBRAIC[AV_CaMKb] * (ALGEBRAIC[AV_CaMKb] + STATES[CaMKt]) - CONSTANTS[AC_bCaMK] * STATES[CaMKt];

    /* IP3 */
    ALGEBRAIC[AV_POip3] = CONSTANTS[AC_tauip3] * CONSTANTS[AC_IP3_IP3] * STATES[Cass] * (1.0 - STATES[IP3_u]) / ((1.0 + CONSTANTS[AC_IP3_IP3] * CONSTANTS[AC_k0] / CONSTANTS[AC_k0a]) * (1.0 + STATES[Cass] * CONSTANTS[AC_k1] / CONSTANTS[AC_k1a]));
    RATES[IP3_u] = STATES[Cass] * CONSTANTS[AC_k2] * (1.0 - STATES[IP3_u]) - CONSTANTS[AC_k2a] * STATES[IP3_u];
    ALGEBRAIC[AV_Jip3] = 10.92 * ALGEBRAIC[AV_POip3] * (STATES[Cajsr] - STATES[Cass]);

    /* IpCa */
    ALGEBRAIC[AV_IpCa_IpCa] = CONSTANTS[AC_GpCa] * STATES[Casl] / (CONSTANTS[AC_KmCap] + STATES[Casl]);

    /* diff */
    ALGEBRAIC[AV_Jdiff] = (STATES[Cass] - STATES[Casl]) / CONSTANTS[AC_sstau];
    ALGEBRAIC[AV_JdiffK] = (STATES[Kss] - STATES[Ksl]) / CONSTANTS[AC_sstau];
    ALGEBRAIC[AV_JdiffNa] = (STATES[Nass] - STATES[Nasl]) / CONSTANTS[AC_sstau];
    ALGEBRAIC[AV_Jgap] = (STATES[Casl] - STATES[Cai]) / CONSTANTS[AC_gaptau];
    ALGEBRAIC[AV_JgapK] = (STATES[Ksl] - STATES[Ki]) / CONSTANTS[AC_gaptau];
    ALGEBRAIC[AV_JgapNa] = (STATES[Nasl] - STATES[Nai]) / CONSTANTS[AC_gaptau];



    /* trans_flux */
    ALGEBRAIC[AV_Jtr1] = (STATES[Cansr] - STATES[Cajsr]) / 120.0;
    ALGEBRAIC[AV_Jtr2] = (STATES[Cansr] - STATES[Cacsr]) / 120.0;

    /* INaCa_i */
    ALGEBRAIC[AV_allo_i] = 1.0 / (1.0 + pow(CONSTANTS[AC_KmCaAct] / STATES[Casl], 2.0));
    ALGEBRAIC[AV_allo_ss] = 1.0 / (1.0 + pow(CONSTANTS[AC_KmCaAct] / STATES[Cass], 2.0));
    ALGEBRAIC[AV_h4_i] = 1.0 + STATES[Nasl] / CONSTANTS[AC_kna1] * (1.0 + STATES[Nasl] / CONSTANTS[AC_kna2]);
    ALGEBRAIC[AV_h4_ss] = 1.0 + STATES[Nass] / CONSTANTS[AC_kna1] * (1.0 + STATES[Nass] / CONSTANTS[AC_kna2]);
    ALGEBRAIC[AV_hca] = exp(CONSTANTS[AC_qca] * STATES[membrane_v] * CONSTANTS[AC_F] / (CONSTANTS[AC_R] * CONSTANTS[AC_T]));
    ALGEBRAIC[AV_hna] = exp(CONSTANTS[AC_qna] * STATES[membrane_v] * CONSTANTS[AC_F] / (CONSTANTS[AC_R] * CONSTANTS[AC_T]));
    ALGEBRAIC[AV_h1_i] = 1.0 + STATES[Nasl] / CONSTANTS[AC_kna3] * (1.0 + ALGEBRAIC[AV_hna]);
    ALGEBRAIC[AV_h1_ss] = 1.0 + STATES[Nass] / CONSTANTS[AC_kna3] * (1.0 + ALGEBRAIC[AV_hna]);
    ALGEBRAIC[AV_h5_i] = STATES[Nasl] * STATES[Nasl] / (ALGEBRAIC[AV_h4_i] * CONSTANTS[AC_kna1] * CONSTANTS[AC_kna2]);
    ALGEBRAIC[AV_h5_ss] = STATES[Nass] * STATES[Nass] / (ALGEBRAIC[AV_h4_ss] * CONSTANTS[AC_kna1] * CONSTANTS[AC_kna2]);
    ALGEBRAIC[AV_h6_i] = 1.0 / ALGEBRAIC[AV_h4_i];
    ALGEBRAIC[AV_h6_ss] = 1.0 / ALGEBRAIC[AV_h4_ss];
    ALGEBRAIC[AV_h7_i] = 1.0 + CONSTANTS[AC_nao] / CONSTANTS[AC_kna3] * (1.0 + 1.0 / ALGEBRAIC[AV_hna]);
    ALGEBRAIC[AV_h7_ss] = 1.0 + CONSTANTS[AC_nao] / CONSTANTS[AC_kna3] * (1.0 + 1.0 / ALGEBRAIC[AV_hna]);
    ALGEBRAIC[AV_h2_i] = STATES[Nasl] * ALGEBRAIC[AV_hna] / (CONSTANTS[AC_kna3] * ALGEBRAIC[AV_h1_i]);
    ALGEBRAIC[AV_h2_ss] = STATES[Nass] * ALGEBRAIC[AV_hna] / (CONSTANTS[AC_kna3] * ALGEBRAIC[AV_h1_ss]);
    ALGEBRAIC[AV_h3_i] = 1.0 / ALGEBRAIC[AV_h1_i];
    ALGEBRAIC[AV_h3_ss] = 1.0 / ALGEBRAIC[AV_h1_ss];
    ALGEBRAIC[AV_h8_i] = CONSTANTS[AC_nao] / (CONSTANTS[AC_kna3] * ALGEBRAIC[AV_hna] * ALGEBRAIC[AV_h7_i]);
    ALGEBRAIC[AV_h8_ss] = CONSTANTS[AC_nao] / (CONSTANTS[AC_kna3] * ALGEBRAIC[AV_hna] * ALGEBRAIC[AV_h7_ss]);
    ALGEBRAIC[AV_h9_i] = 1.0 / ALGEBRAIC[AV_h7_i];
    ALGEBRAIC[AV_h9_ss] = 1.0 / ALGEBRAIC[AV_h7_ss];
    ALGEBRAIC[AV_k6_i] = ALGEBRAIC[AV_h6_i] * STATES[Casl] * CONSTANTS[AC_kcaon];
    ALGEBRAIC[AV_k6_ss] = ALGEBRAIC[AV_h6_ss] * STATES[Cass] * CONSTANTS[AC_kcaon];
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
    ALGEBRAIC[AV_Knai] = CONSTANTS[AC_Knai0] * exp(CONSTANTS[AC_delta] * STATES[membrane_v] * CONSTANTS[AC_F] / (3.0 * CONSTANTS[AC_R] * CONSTANTS[AC_T]));
    ALGEBRAIC[AV_Knao] = CONSTANTS[AC_Knao0] * exp((1.0 - CONSTANTS[AC_delta]) * STATES[membrane_v] * CONSTANTS[AC_F] / (3.0 * CONSTANTS[AC_R] * CONSTANTS[AC_T]));
    ALGEBRAIC[AV_P] = CONSTANTS[AC_eP] / (1.0 + CONSTANTS[AC_H] / CONSTANTS[AC_Khp] + STATES[Nasl] / CONSTANTS[AC_Knap] + STATES[Ksl] / CONSTANTS[AC_Kxkur]);
    ALGEBRAIC[AV_a1] = CONSTANTS[AC_k1p] * pow(STATES[Nasl] / ALGEBRAIC[AV_Knai],
                                               3.0) / (pow(1.0 + STATES[Nasl] / ALGEBRAIC[AV_Knai], 3.0) + pow(1.0 + STATES[Ksl] / CONSTANTS[AC_Kki], 2.0) - 1.0);
    ALGEBRAIC[AV_a3] = CONSTANTS[AC_k3p] * pow(CONSTANTS[AC_ko] / CONSTANTS[AC_Kko], 2.0) / (pow(1.0 + CONSTANTS[AC_nao] / ALGEBRAIC[AV_Knao], 3.0) + pow(1.0 + CONSTANTS[AC_ko] / CONSTANTS[AC_Kko], 2.0) - 1.0);
    ALGEBRAIC[AV_b2] = CONSTANTS[AC_k2m] * pow(CONSTANTS[AC_nao] / ALGEBRAIC[AV_Knao],
                                               3.0) / (pow(1.0 + CONSTANTS[AC_nao] / ALGEBRAIC[AV_Knao], 3.0) + pow(1.0 + CONSTANTS[AC_ko] / CONSTANTS[AC_Kko], 2.0) - 1.0);
    ALGEBRAIC[AV_b3] = CONSTANTS[AC_k3m] * ALGEBRAIC[AV_P] * CONSTANTS[AC_H] / (1.0 + CONSTANTS[AC_MgATP] / CONSTANTS[AC_Kmgatp]);
    ALGEBRAIC[AV_b4] = CONSTANTS[AC_k4m] * pow(STATES[Ksl] / CONSTANTS[AC_Kki], 2.0) / (pow(1.0 + STATES[Nasl] / ALGEBRAIC[AV_Knai], 3.0) + pow(1.0 + STATES[Ksl] / CONSTANTS[AC_Kki], 2.0) - 1.0);
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
    ALGEBRAIC[AV_dkmplb] = CONSTANTS[AC_dkmplbbar] * ALGEBRAIC[AV_CaMKa] / (ALGEBRAIC[AV_CaMKa] + CONSTANTS[AC_KmCaMK]);
    ALGEBRAIC[AV_dqupcamk] = CONSTANTS[AC_dqupcamkbar] * ALGEBRAIC[AV_CaMKa] / (ALGEBRAIC[AV_CaMKa] + CONSTANTS[AC_KmCaMK]);
    ALGEBRAIC[AV_Jup1] = 0.0002 * (ALGEBRAIC[AV_dqupcamk] + 1.0) / (1.0 + (CONSTANTS[AC_kmup] - ALGEBRAIC[AV_dkmplb]) / STATES[Casl]) - 0.00105 * STATES[Cansr] / CONSTANTS[AC_nsrbar];
    ALGEBRAIC[AV_Jup2] = 0.0026 * (ALGEBRAIC[AV_dqupcamk] + 1.0) / (1.0 + (CONSTANTS[AC_kmup] - ALGEBRAIC[AV_dkmplb]) / STATES[Cai]) - 0.0042 * STATES[Cansr] / CONSTANTS[AC_nsrbar];

    /* reversal_potentials */
    ALGEBRAIC[AV_ECa] = CONSTANTS[AC_R] * CONSTANTS[AC_T] / (2.0 * CONSTANTS[AC_F]) * log(CONSTANTS[AC_cao] / STATES[Casl]);
    ALGEBRAIC[AV_EK] = CONSTANTS[AC_R] * CONSTANTS[AC_T] / CONSTANTS[AC_F] * log(CONSTANTS[AC_ko] / STATES[Ksl]);
    ALGEBRAIC[AV_ENa] = CONSTANTS[AC_R] * CONSTANTS[AC_T] / CONSTANTS[AC_F] * log(CONSTANTS[AC_nao] / STATES[Nasl]);
    ALGEBRAIC[AV_EKs] = CONSTANTS[AC_R] * CONSTANTS[AC_T] / CONSTANTS[AC_F] * log((CONSTANTS[AC_ko] + CONSTANTS[AC_PKNa] * CONSTANTS[AC_nao]) / (STATES[Ksl] + CONSTANTS[AC_PKNa] * STATES[Nasl]));

    /* ICaT */
    ALGEBRAIC[AV_bss] = 1.0 / (1.0 + exp(-(STATES[membrane_v] + 30.0) / 7.0));
    ALGEBRAIC[AV_gss] = 1.0 / (1.0 + exp((STATES[membrane_v] + 61.0) / 5.0));
    ALGEBRAIC[AV_taub] = 1.0 / (1.068 * exp((STATES[membrane_v] + 16.3) / 30.0) + 1.068 * exp(-(STATES[membrane_v] + 16.3) / 30.0));
    ALGEBRAIC[AV_taug] = 1.0 / (0.015 * exp((STATES[membrane_v] + 71.7) / 15.4) + 0.015 * exp(-(STATES[membrane_v] + 71.7) / 83.3));
    ALGEBRAIC[AV_ICaT_ICaT] = CONSTANTS[AC_GCaT] * STATES[ICaT_b] * STATES[ICaT_g] * (STATES[membrane_v] - ALGEBRAIC[AV_ECa]);
    RATES[ICaT_b] = (ALGEBRAIC[AV_bss] - STATES[ICaT_b]) / ALGEBRAIC[AV_taub];
    RATES[ICaT_g] = (ALGEBRAIC[AV_gss] - STATES[ICaT_g]) / ALGEBRAIC[AV_taug];

    /* IK1 */
    ALGEBRAIC[AV_rk1] = 1.0 / (1.0 + exp((STATES[membrane_v] + 116.0 - 5.5 * CONSTANTS[AC_ko]) / 11.0));
    ALGEBRAIC[AV_txk1] = 122.2 / (exp(-(STATES[membrane_v] + 127.2) / 20.36) + exp((STATES[membrane_v] + 236.8) / 69.33));
    ALGEBRAIC[AV_xk1ss] = 1.0 / (1.0 + exp(-(STATES[membrane_v] + 2.5538 * CONSTANTS[AC_ko] + 144.59) / (1.5692 * CONSTANTS[AC_ko] + 3.8115)));
    ALGEBRAIC[AV_IK1_IK1] = CONSTANTS[AC_GK1] * 2.3238 * sqrt(CONSTANTS[AC_ko] / 5.4) * ALGEBRAIC[AV_rk1] * STATES[IK1_xk1] * (STATES[membrane_v] - ALGEBRAIC[AV_EK]);
    RATES[IK1_xk1] = (ALGEBRAIC[AV_xk1ss] - STATES[IK1_xk1]) / ALGEBRAIC[AV_txk1];

    /* IKr */
    ALGEBRAIC[AV_Axrf] = 1.0 / (1.0 + exp((STATES[membrane_v] + 54.81) / 38.21));
    ALGEBRAIC[AV_rkr] = 1.0 / (1.0 + exp((STATES[membrane_v] + 55.0) / (0.32 * 75.0))) * 1.0 / (1.0 + exp((STATES[membrane_v] - 10.0) / (0.32 * 30.0)));
    ALGEBRAIC[AV_txrf] = 12.98 + 1.0 / (0.3652 * exp((STATES[membrane_v] + 17.6 - 31.66) / 3.869) + 4.123e-05 * exp(-(STATES[membrane_v] + 17.6 - 47.78) / 20.38));
    ALGEBRAIC[AV_txrs] = 1.865 + 1.0 / (0.06629 * exp((STATES[membrane_v] + 17.2 - 34.7) / 7.355) + 1.128e-05 * exp(-(STATES[membrane_v] + 17.2 - 29.74) / 25.94));
    ALGEBRAIC[AV_xrss] = 1.0 / (1.0 + exp(-(STATES[membrane_v] + 8.337) / 6.789));
    ALGEBRAIC[AV_Axrs] = 1.0 - ALGEBRAIC[AV_Axrf];
    RATES[IKr_xrf] = (ALGEBRAIC[AV_xrss] - STATES[IKr_xrf]) / ALGEBRAIC[AV_txrf];
    RATES[IKr_xrs] = (ALGEBRAIC[AV_xrss] - STATES[IKr_xrs]) / ALGEBRAIC[AV_txrs];
    ALGEBRAIC[AV_xr] = ALGEBRAIC[AV_Axrf] * STATES[IKr_xrf] + ALGEBRAIC[AV_Axrs] * STATES[IKr_xrs];
    ALGEBRAIC[AV_IKr_IKr] = CONSTANTS[AC_GKr] * sqrt(CONSTANTS[AC_ko] / 5.4) * ALGEBRAIC[AV_xr] * ALGEBRAIC[AV_rkr] * (STATES[membrane_v] - ALGEBRAIC[AV_EK]);

    /* IKs */
    ALGEBRAIC[AV_KsCa] = 1.0 + 0.6 / (1.0 + pow(3.8e-05 / STATES[Casl], 1.4));
    ALGEBRAIC[AV_txs1] = 817.3 + 1.0 / (0.0002326 * exp((STATES[membrane_v] + 48.28) / 17.8) + 0.001292 * exp(-(STATES[membrane_v] + 210.0) / 230.0));
    ALGEBRAIC[AV_txs2] = 1.0 / (0.01 * exp((STATES[membrane_v] - 50.0) / 20.0) + 0.0193 * exp(-(STATES[membrane_v] + 66.54) / 31.0));
    ALGEBRAIC[AV_xs1ss] = 1.0 / (1.0 + exp(-(STATES[membrane_v] + 11.6) / 8.932));
    ALGEBRAIC[AV_IKs_IKs] = CONSTANTS[AC_GKs] * ALGEBRAIC[AV_KsCa] * STATES[IKs_xs1] * STATES[IKs_xs2] * (STATES[membrane_v] - ALGEBRAIC[AV_EKs]);
    ALGEBRAIC[AV_xs2ss] = ALGEBRAIC[AV_xs1ss];
    RATES[IKs_xs1] = (ALGEBRAIC[AV_xs1ss] - STATES[IKs_xs1]) / ALGEBRAIC[AV_txs1];
    RATES[IKs_xs2] = (ALGEBRAIC[AV_xs2ss] - STATES[IKs_xs2]) / ALGEBRAIC[AV_txs2];

    /* INa */
    ALGEBRAIC[AV_fINap] = 1.0 / (1.0 + CONSTANTS[AC_KmCaMK] / ALGEBRAIC[AV_CaMKa]);
    ALGEBRAIC[AV_hssp] = 1.0 / (1.0 + exp((STATES[membrane_v] + 84.7) / 6.22));
    ALGEBRAIC[AV_thf] = 1.0 / (3.686e-06 * exp(-(STATES[membrane_v] + 3.8875) / 7.8579) + 16.0 * exp((STATES[membrane_v] - 0.4963) / 9.1843));
    ALGEBRAIC[AV_ths] = 1.0 / (0.009794 * exp(-(STATES[membrane_v] + 17.95) / 28.05) + 0.3343 * exp((STATES[membrane_v] + 5.73) / 56.66));
    ALGEBRAIC[AV_tj] = 4.859 + 1.0 / (0.8628 * exp(-(STATES[membrane_v] + 116.7258) / 7.6005) + 1.1096 * exp((STATES[membrane_v] + 6.2719) / 9.0358));
    ALGEBRAIC[AV_hss] = 1.0 / (1.0 + exp((STATES[membrane_v] + CONSTANTS[AC_hssV1]) / CONSTANTS[AC_hssV2]));
    ALGEBRAIC[AV_mss] = 1.0 / (1.0 + exp(-(STATES[membrane_v] + CONSTANTS[AC_mssV1]) / CONSTANTS[AC_mssV2]));
    ALGEBRAIC[AV_thsp] = 3.0 * ALGEBRAIC[AV_ths];
    ALGEBRAIC[AV_tjp] = 1.46 * ALGEBRAIC[AV_tj];
    ALGEBRAIC[AV_tm] = 1.0 / (CONSTANTS[AC_mtD1] * exp((STATES[membrane_v] + CONSTANTS[AC_mtV1]) / CONSTANTS[AC_mtV2]) + CONSTANTS[AC_mtD2] * exp(-(STATES[membrane_v] + CONSTANTS[AC_mtV3]) / CONSTANTS[AC_mtV4]));
    ALGEBRAIC[AV_h] = CONSTANTS[AC_Ahf] * STATES[INa_hf] + CONSTANTS[AC_Ahs] * STATES[INa_hs];
    ALGEBRAIC[AV_hp] = CONSTANTS[AC_Ahf] * STATES[INa_hf] + CONSTANTS[AC_Ahs] * STATES[INa_hsp];
    ALGEBRAIC[AV_jss] = ALGEBRAIC[AV_hss];
    RATES[INa_hf] = (ALGEBRAIC[AV_hss] - STATES[INa_hf]) / ALGEBRAIC[AV_thf];
    RATES[INa_hs] = (ALGEBRAIC[AV_hss] - STATES[INa_hs]) / ALGEBRAIC[AV_ths];
    RATES[INa_hsp] = (ALGEBRAIC[AV_hssp] - STATES[INa_hsp]) / ALGEBRAIC[AV_thsp];
    RATES[INa_m] = (ALGEBRAIC[AV_mss] - STATES[INa_m]) / ALGEBRAIC[AV_tm];
    ALGEBRAIC[AV_INa_INa] = CONSTANTS[AC_GNa] * (STATES[membrane_v] - ALGEBRAIC[AV_ENa]) * pow(STATES[INa_m], 3.0) * ((1.0 - ALGEBRAIC[AV_fINap]) * ALGEBRAIC[AV_h] * STATES[INa_j] + ALGEBRAIC[AV_fINap] * ALGEBRAIC[AV_hp] * STATES[INa_jp]);
    RATES[INa_j] = (ALGEBRAIC[AV_jss] - STATES[INa_j]) / ALGEBRAIC[AV_tj];
    RATES[INa_jp] = (ALGEBRAIC[AV_jss] - STATES[INa_jp]) / ALGEBRAIC[AV_tjp];

    /* If */
    ALGEBRAIC[AV_tauy] = 2000.0 / (exp((STATES[membrane_v] + 57.0) / 60.0) + exp(-(STATES[membrane_v] + 132.0) / 10.0));
    ALGEBRAIC[AV_yss] = 1.0 / (1.0 + exp((STATES[membrane_v] + 87.0) / 9.5));
    ALGEBRAIC[AV_IfK] = CONSTANTS[AC_GfK] * STATES[If_y] * STATES[If_y] * (STATES[membrane_v] - ALGEBRAIC[AV_EK]);
    ALGEBRAIC[AV_IfNa] = CONSTANTS[AC_GfNa] * STATES[If_y] * STATES[If_y] * (STATES[membrane_v] - ALGEBRAIC[AV_ENa]);
    RATES[If_y] = (ALGEBRAIC[AV_yss] - STATES[If_y]) / ALGEBRAIC[AV_tauy];
    ALGEBRAIC[AV_If_If] = ALGEBRAIC[AV_IfNa] + ALGEBRAIC[AV_IfK];

    /* Isus */
    ALGEBRAIC[AV_asus] = 1.0 / (1.0 + exp(-(STATES[membrane_v] - 12.0) / 16.0));
    ALGEBRAIC[AV_Isus_Isus] = CONSTANTS[AC_Gsus] * ALGEBRAIC[AV_asus] * (STATES[membrane_v] - ALGEBRAIC[AV_EK]);

    /* Ito */
    ALGEBRAIC[AV_ass] = 1.0 / (1.0 + exp((20.0 - STATES[membrane_v]) / 13.0));
    ALGEBRAIC[AV_iss] = 1.0 / (1.0 + exp((27.0 + STATES[membrane_v]) / 13.0));
    ALGEBRAIC[AV_taua] = 1.0515 / (1.0 / (1.2089 * (1.0 + exp(-(STATES[membrane_v] - 18.4099) / 29.3814))) + 3.5 / (1.0 + exp((STATES[membrane_v] + 100.0) / 29.3814)));
    ALGEBRAIC[AV_tauif] = 6.162 + 1.0 / (0.3933 * exp(-(STATES[membrane_v] + 100.0) / 100.0) + 0.08004 * exp((STATES[membrane_v] - 8.0) / 8.59));
    ALGEBRAIC[AV_tauis] = 43.0 + 1.0 / (0.001416 * exp(-(STATES[membrane_v] + 96.52) / 59.05) + 1.78e-08 * exp((STATES[membrane_v] + 114.1) / 8.079));
    ALGEBRAIC[AV_Ito_Ito] = CONSTANTS[AC_Gto] * STATES[Ito_a] * STATES[Ito_i1] * STATES[Ito_i2] * (STATES[membrane_v] - ALGEBRAIC[AV_EK]);
    RATES[Ito_a] = (ALGEBRAIC[AV_ass] - STATES[Ito_a]) / ALGEBRAIC[AV_taua];
    RATES[Ito_i1] = (ALGEBRAIC[AV_iss] - STATES[Ito_i1]) / ALGEBRAIC[AV_tauis];
    RATES[Ito_i2] = (ALGEBRAIC[AV_iss] - STATES[Ito_i2]) / ALGEBRAIC[AV_tauif];

    /* INaL */
    ALGEBRAIC[AV_fINaLp] = 1.0 / (1.0 + CONSTANTS[AC_KmCaMK] / ALGEBRAIC[AV_CaMKa]);
    ALGEBRAIC[AV_hLss] = 1.0 / (1.0 + exp((STATES[membrane_v] + 87.61) / 7.488));
    ALGEBRAIC[AV_hLssp] = 1.0 / (1.0 + exp((STATES[membrane_v] + 93.81) / 7.488));
    ALGEBRAIC[AV_mLss] = 1.0 / (1.0 + exp(-(STATES[membrane_v] + 42.85) / 5.264));
    ALGEBRAIC[AV_tmL] = ALGEBRAIC[AV_tm];
    ALGEBRAIC[AV_INaL_INaL] = CONSTANTS[AC_GNaL] * (STATES[membrane_v] - ALGEBRAIC[AV_ENa]) * STATES[INaL_mL] * ((1.0 - ALGEBRAIC[AV_fINaLp]) * STATES[INaL_hL] + ALGEBRAIC[AV_fINaLp] * STATES[INaL_hLp]);
    RATES[INaL_hL] = (ALGEBRAIC[AV_hLss] - STATES[INaL_hL]) / CONSTANTS[AC_thL];
    RATES[INaL_mL] = (ALGEBRAIC[AV_mLss] - STATES[INaL_mL]) / ALGEBRAIC[AV_tmL];
    RATES[INaL_hLp] = (ALGEBRAIC[AV_hLssp] - STATES[INaL_hLp]) / CONSTANTS[AC_thLp];

    /* ICaL */
    ALGEBRAIC[AV_Afcaf] = 0.3 + 0.6 / (1.0 + exp((STATES[membrane_v] - 10.0) / 10.0));
    ALGEBRAIC[AV_dss] = 1.0 / (1.0 + exp(-(STATES[membrane_v] + 3.94 + 3.3) / 4.23));
    ALGEBRAIC[AV_fICaLp] = 1.0 / (1.0 + CONSTANTS[AC_KmCaMK] / ALGEBRAIC[AV_CaMKa]);
    ALGEBRAIC[AV_fss] = 1.0 / (1.0 + exp((STATES[membrane_v] + 19.58 + 3.3) / 3.696));
    ALGEBRAIC[AV_km2n] = STATES[ICaL_jca] * 1.0;
    ALGEBRAIC[AV_td] = 0.6 + 1.0 / (exp(-0.05 * (STATES[membrane_v] + 6.0)) + exp(0.09 * (STATES[membrane_v] + 14.0)));
    ALGEBRAIC[AV_tfcaf] = 0.72 * (7.0 + 1.0 / (0.04 * exp(-(STATES[membrane_v] + 15.19 - 4.0) / 7.0) + 0.04 * exp((STATES[membrane_v] + 15.19 - 4.0) / 7.0)));
    ALGEBRAIC[AV_tfcas] = 0.49 * (100.0 + 1.0 / (0.00012 * exp(-(STATES[membrane_v] + 15.19) / 3.0) + 0.00012 * exp((STATES[membrane_v] + 15.19) / 7.0)));
    ALGEBRAIC[AV_tff] = 7.0 + 1.0 / (0.0045 * exp(-(STATES[membrane_v] + 20.0 + 15.19) / 10.0) + 0.0045 * exp((STATES[membrane_v] + 20.0 + 15.19) / 10.0));
    ALGEBRAIC[AV_tfs] = 1000.0 + 1.0 / (3.5e-05 * exp(-(STATES[membrane_v] + 5.0 + 15.19) / 4.0) + 3.5e-05 * exp((STATES[membrane_v] + 5.0 + 15.19) / 6.0));
    ALGEBRAIC[AV_Afcas] = 1.0 - ALGEBRAIC[AV_Afcaf];
    ALGEBRAIC[AV_anca] = 1.0 / (CONSTANTS[AC_k2n] / ALGEBRAIC[AV_km2n] + pow(1.0 + CONSTANTS[AC_Kmn] / STATES[Cass], 4.0));
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
    RATES[Cansr] = ALGEBRAIC[AV_Jup1] + ALGEBRAIC[AV_Jup2] - (ALGEBRAIC[AV_Jtr1] * CONSTANTS[AC_vjsr] / CONSTANTS[AC_vnsr] + ALGEBRAIC[AV_Jtr2] * CONSTANTS[AC_vcsr] / CONSTANTS[AC_vnsr]);
    RATES[Ki] = ALGEBRAIC[AV_JgapK] * CONSTANTS[AC_vsl] / CONSTANTS[AC_vmyo];
    RATES[Nai] = ALGEBRAIC[AV_JgapNa] * CONSTANTS[AC_vsl] / CONSTANTS[AC_vmyo];
    ALGEBRAIC[AV_Bcacsr] = 1.0 / (1.0 + CONSTANTS[AC_csqnmax] * CONSTANTS[AC_kmcsqn] / pow(CONSTANTS[AC_kmcsqn] + STATES[Cacsr], 2.0));
    ALGEBRAIC[AV_Bcai] = 1.0 / (1.0 + CONSTANTS[AC_cmdnmax] * CONSTANTS[AC_kmcmdn] / pow(CONSTANTS[AC_kmcmdn] + STATES[Cai], 2.0) + CONSTANTS[AC_trpnmax] * CONSTANTS[AC_kmtrpn] / pow(CONSTANTS[AC_kmtrpn] + STATES[Cai], 2.0));
    ALGEBRAIC[AV_Bcajsr] = 1.0 / (1.0 + CONSTANTS[AC_csqnmaxsl] * CONSTANTS[AC_kmcsqn] / pow(CONSTANTS[AC_kmcsqn] + STATES[Cajsr], 2.0));
    ALGEBRAIC[AV_Bcasl] = 1.0 / (1.0 + CONSTANTS[AC_cmdnmaxsl] * CONSTANTS[AC_kmcmdn] / pow(CONSTANTS[AC_kmcmdn] + STATES[Casl], 2.0) + CONSTANTS[AC_trpnmaxsl] * CONSTANTS[AC_kmtrpn] / pow(CONSTANTS[AC_kmtrpn] + STATES[Casl], 2.0));
    ALGEBRAIC[AV_Bcass] = 1.0 / (1.0 + CONSTANTS[AC_BSRmax] * CONSTANTS[AC_KmBSR] / pow(CONSTANTS[AC_KmBSR] + STATES[Cass], 2.0) + CONSTANTS[AC_BSLmax] * CONSTANTS[AC_KmBSL] / pow(CONSTANTS[AC_KmBSL] + STATES[Cass], 2.0));
    RATES[Cacsr] = ALGEBRAIC[AV_Bcacsr] * (ALGEBRAIC[AV_Jtr2] - STATES[ryr_Jrel2]);
    RATES[Cai] = ALGEBRAIC[AV_Bcai] * (ALGEBRAIC[AV_Jgap] * CONSTANTS[AC_vsl] / CONSTANTS[AC_vmyo] + STATES[ryr_Jrel2] * CONSTANTS[AC_vcsr] / CONSTANTS[AC_vmyo] - ALGEBRAIC[AV_Jup2] * CONSTANTS[AC_vnsr] / CONSTANTS[AC_vmyo]);
    RATES[Cajsr] = ALGEBRAIC[AV_Bcajsr] * (ALGEBRAIC[AV_Jtr1] - (STATES[ryr_Jrel1] + ALGEBRAIC[AV_Jip3]));

    /* membrane */
    ALGEBRAIC[AV_vffrt] = STATES[membrane_v] * CONSTANTS[AC_F] * CONSTANTS[AC_F] / (CONSTANTS[AC_R] * CONSTANTS[AC_T]);
    ALGEBRAIC[AV_vfrt] = STATES[membrane_v] * CONSTANTS[AC_F] / (CONSTANTS[AC_R] * CONSTANTS[AC_T]);


    /* ryr */
    ALGEBRAIC[AV_REL2] = ALGEBRAIC[AV_Jgap] * CONSTANTS[AC_vsl] / CONSTANTS[AC_vmyo] + -ALGEBRAIC[AV_Jup2] * CONSTANTS[AC_vnsr] / CONSTANTS[AC_vmyo] + STATES[ryr_Jrel2] * CONSTANTS[AC_vcsr] / CONSTANTS[AC_vmyo];
    ALGEBRAIC[AV_ireltau] = 2.0 * (1.0 + 1.0 * 1.0 / (1.0 + pow(0.28 / ALGEBRAIC[AV_CaMKa], 8.0))) / (1.0 + 0.0123 / STATES[Cajsr]);
    ALGEBRAIC[AV_ireltau2] = 6.0 * (1.0 + 1.0 * 1.0 / (1.0 + pow(0.28 / ALGEBRAIC[AV_CaMKa], 8.0))) / (1.0 + 0.0123 / STATES[Cacsr]);
    ALGEBRAIC[AV_irelss2] = ((ALGEBRAIC[AV_REL2] > 0.0) ? 91.0 * (1.0 + 1.0 * 1.0 / (1.0 + pow(0.28 / ALGEBRAIC[AV_CaMKa], 8.0))) * ALGEBRAIC[AV_REL2] / (1.0 + pow(1.0 / STATES[Cacsr], 8.0)) : 0.0);
    RATES[ryr_Jrel2] = (ALGEBRAIC[AV_irelss2] - STATES[ryr_Jrel2]) / ALGEBRAIC[AV_ireltau2];

    /* *remaining* */
    ALGEBRAIC[AV_PhiCaK] = 1.0 * ALGEBRAIC[AV_vffrt] * (0.75 * STATES[Kss] * exp(1.0 * ALGEBRAIC[AV_vfrt]) - 0.75 * CONSTANTS[AC_ko]) / (exp(1.0 * ALGEBRAIC[AV_vfrt]) - 1.0);
    ALGEBRAIC[AV_PhiCaL] = 4.0 * ALGEBRAIC[AV_vffrt] * (STATES[Cass] * exp(2.0 * ALGEBRAIC[AV_vfrt]) - 0.341 * CONSTANTS[AC_cao]) / (exp(2.0 * ALGEBRAIC[AV_vfrt]) - 1.0);
    ALGEBRAIC[AV_PhiCaNa] = 1.0 * ALGEBRAIC[AV_vffrt] * (0.75 * STATES[Nass] * exp(1.0 * ALGEBRAIC[AV_vfrt]) - 0.75 * CONSTANTS[AC_nao]) / (exp(1.0 * ALGEBRAIC[AV_vfrt]) - 1.0);
    ALGEBRAIC[AV_ICab_ICab] = CONSTANTS[AC_PCab] * 4.0 * ALGEBRAIC[AV_vffrt] * (STATES[Casl] * exp(2.0 * ALGEBRAIC[AV_vfrt]) - 0.341 * CONSTANTS[AC_cao]) / (exp(2.0 * ALGEBRAIC[AV_vfrt]) - 1.0);
    ALGEBRAIC[AV_INab_INab] = CONSTANTS[AC_PNab] * ALGEBRAIC[AV_vffrt] * (STATES[Nasl] * exp(ALGEBRAIC[AV_vfrt]) - CONSTANTS[AC_nao]) / (exp(ALGEBRAIC[AV_vfrt]) - 1.0);
    RATES[Ksl] = -(ALGEBRAIC[AV_Ito_Ito] + ALGEBRAIC[AV_Isus_Isus] + ALGEBRAIC[AV_IKr_IKr] + ALGEBRAIC[AV_IKs_IKs] + ALGEBRAIC[AV_IfK] + ALGEBRAIC[AV_IK1_IK1] + ALGEBRAIC[Istim] - 2.0 * ALGEBRAIC[AV_INaK_INaK]) * CONSTANTS[AC_cm] * CONSTANTS[AC_Acap] / (CONSTANTS[AC_F] * CONSTANTS[AC_vsl]) + ALGEBRAIC[AV_JdiffK] * CONSTANTS[AC_vss] / CONSTANTS[AC_vsl] + -ALGEBRAIC[AV_JgapK];
    ALGEBRAIC[AV_ICaK] = (1.0 - ALGEBRAIC[AV_fICaLp]) * CONSTANTS[AC_PCaK] * ALGEBRAIC[AV_PhiCaK] * STATES[ICaL_d] * (ALGEBRAIC[AV_f] * (1.0 - STATES[ICaL_nca]) + STATES[ICaL_jca] * ALGEBRAIC[AV_fca] * STATES[ICaL_nca]) + ALGEBRAIC[AV_fICaLp] * CONSTANTS[AC_PCaKp] * ALGEBRAIC[AV_PhiCaK] * STATES[ICaL_d] * (ALGEBRAIC[AV_fp] * (1.0 - STATES[ICaL_nca]) + STATES[ICaL_jca] * ALGEBRAIC[AV_fcap] * STATES[ICaL_nca]);
    ALGEBRAIC[AV_ICaL_ICaL] = (1.0 - ALGEBRAIC[AV_fICaLp]) * CONSTANTS[AC_PCa] * ALGEBRAIC[AV_PhiCaL] * STATES[ICaL_d] * (ALGEBRAIC[AV_f] * (1.0 - STATES[ICaL_nca]) + STATES[ICaL_jca] * ALGEBRAIC[AV_fca] * STATES[ICaL_nca]) + ALGEBRAIC[AV_fICaLp] * CONSTANTS[AC_PCap] * ALGEBRAIC[AV_PhiCaL] * STATES[ICaL_d] * (ALGEBRAIC[AV_fp] * (1.0 - STATES[ICaL_nca]) + STATES[ICaL_jca] * ALGEBRAIC[AV_fcap] * STATES[ICaL_nca]);
    ALGEBRAIC[AV_ICaNa] = (1.0 - ALGEBRAIC[AV_fICaLp]) * CONSTANTS[AC_PCaNa] * ALGEBRAIC[AV_PhiCaNa] * STATES[ICaL_d] * (ALGEBRAIC[AV_f] * (1.0 - STATES[ICaL_nca]) + STATES[ICaL_jca] * ALGEBRAIC[AV_fca] * STATES[ICaL_nca]) + ALGEBRAIC[AV_fICaLp] * CONSTANTS[AC_PCaNap] * ALGEBRAIC[AV_PhiCaNa] * STATES[ICaL_d] * (ALGEBRAIC[AV_fp] * (1.0 - STATES[ICaL_nca]) + STATES[ICaL_jca] * ALGEBRAIC[AV_fcap] * STATES[ICaL_nca]);
    RATES[Casl] = ALGEBRAIC[AV_Bcasl] * (-(ALGEBRAIC[AV_IpCa_IpCa] + ALGEBRAIC[AV_ICab_ICab] + ALGEBRAIC[AV_ICaT_ICaT] - 2.0 * ALGEBRAIC[AV_INaCa_i_INaCa_i]) * CONSTANTS[AC_cm] * CONSTANTS[AC_Acap] / (2.0 * CONSTANTS[AC_F] * CONSTANTS[AC_vsl]) - ALGEBRAIC[AV_Jup1] * CONSTANTS[AC_vnsr] / CONSTANTS[AC_vsl] + -ALGEBRAIC[AV_Jgap] + ALGEBRAIC[AV_Jdiff] * CONSTANTS[AC_vss] / CONSTANTS[AC_vsl]);
    RATES[Nasl] = -(ALGEBRAIC[AV_INa_INa] + ALGEBRAIC[AV_INaL_INaL] + 3.0 * ALGEBRAIC[AV_INaCa_i_INaCa_i] + 3.0 * ALGEBRAIC[AV_INaK_INaK] + ALGEBRAIC[AV_IfNa] + ALGEBRAIC[AV_INab_INab]) * CONSTANTS[AC_Acap] * CONSTANTS[AC_cm] / (CONSTANTS[AC_F] * CONSTANTS[AC_vsl]) + ALGEBRAIC[AV_JdiffNa] * CONSTANTS[AC_vss] / CONSTANTS[AC_vsl] + -ALGEBRAIC[AV_JgapNa];
    RATES[Cass] = ALGEBRAIC[AV_Bcass] * (-(ALGEBRAIC[AV_ICaL_ICaL] - 2.0 * ALGEBRAIC[AV_INaCa_ss]) * CONSTANTS[AC_cm] * CONSTANTS[AC_Acap] / (2.0 * CONSTANTS[AC_F] * CONSTANTS[AC_vss]) + (STATES[ryr_Jrel1] + ALGEBRAIC[AV_Jip3]) * CONSTANTS[AC_vjsr] / CONSTANTS[AC_vss] - ALGEBRAIC[AV_Jdiff]);
    RATES[Kss] = -ALGEBRAIC[AV_ICaK] * CONSTANTS[AC_cm] * CONSTANTS[AC_Acap] / (CONSTANTS[AC_F] * CONSTANTS[AC_vss]) - ALGEBRAIC[AV_JdiffK];
    RATES[Nass] = -(ALGEBRAIC[AV_ICaNa] + 3.0 * ALGEBRAIC[AV_INaCa_ss]) * CONSTANTS[AC_cm] * CONSTANTS[AC_Acap] / (CONSTANTS[AC_F] * CONSTANTS[AC_vss]) - ALGEBRAIC[AV_JdiffNa];
    ALGEBRAIC[AV_REL] = -(ALGEBRAIC[AV_ICaL_ICaL] * CONSTANTS[AC_Acap] / (CONSTANTS[AC_vss] * CONSTANTS[AC_F] * 2.0) + -(STATES[ryr_Jrel1] + ALGEBRAIC[AV_Jip3]) * CONSTANTS[AC_vjsr] / CONSTANTS[AC_vss] + ALGEBRAIC[AV_Jdiff]);
    ALGEBRAIC[AV_irelss] = ((ALGEBRAIC[AV_REL] > 0.0) ? 15.0 * (1.0 + 1.0 * 1.0 / (1.0 + pow(0.28 / ALGEBRAIC[AV_CaMKa], 8.0))) * ALGEBRAIC[AV_REL] / (1.0 + pow(1.0 / STATES[Cajsr], 8.0)) : 0.0);
    RATES[ryr_Jrel1] = (ALGEBRAIC[AV_irelss] - STATES[ryr_Jrel1]) / ALGEBRAIC[AV_ireltau];

    ALGEBRAIC[Iion_cm] = ALGEBRAIC[AV_INa_INa] + ALGEBRAIC[AV_INaL_INaL] + ALGEBRAIC[AV_Ito_Ito] + ALGEBRAIC[AV_Isus_Isus] + ALGEBRAIC[AV_ICaL_ICaL] + ALGEBRAIC[AV_ICaT_ICaT] + ALGEBRAIC[AV_ICaNa] + ALGEBRAIC[AV_ICaK] + ALGEBRAIC[AV_IKr_IKr] + ALGEBRAIC[AV_IKs_IKs] + ALGEBRAIC[AV_If_If] + ALGEBRAIC[AV_IK1_IK1] + ALGEBRAIC[AV_INaCa_i_INaCa_i] + ALGEBRAIC[AV_INaCa_ss] + ALGEBRAIC[AV_INaK_INaK] + ALGEBRAIC[AV_INab_INab] + ALGEBRAIC[AV_IpCa_IpCa] + ALGEBRAIC[AV_ICab_ICab];





    ALGEBRAIC[Istim] = computeIstim(VOI, CONSTANTS);

    if (solveVmWithinODESolver)
    {
        RATES[membrane_v] = -ALGEBRAIC[Iion_cm] - ALGEBRAIC[Istim];
    }

}
