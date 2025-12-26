#include "Gaur_2021Names.H"
#include "stimulusIO.H"

static const char* GaurSTATES_NAMES[NUM_STATES] = {
    "cell_v",
    "nai",
    "nass",
    "ki",
    "kss",
    "cai",
    "cai2",
    "cass",
    "cansr",
    "cajsr",
    "cacsr",
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

static const char* GaurALGEBRAIC_NAMES[NUM_ALGEBRAIC] = {
    "AV_Jdiff",
    "AV_JdiffK",
    "AV_JdiffNa",
    "AV_EK",
    "AV_ENa",
    "AV_vffrt",
    "AV_vfrt",
    "AV_CaMKb",
    "AV_EKs",
    "AV_CaMKa",
    "AV_diff_CaMKt",
    "AV_CaMK_f",
    "AV_aa_h",
    "AV_aa_j",
    "AV_aa_m",
    "AV_bb_h",
    "AV_bb_j",
    "AV_bb_m",
    "AV_h_inf",
    "AV_m_inf",
    "AV_INa",
    "AV_j_inf",
    "AV_tau_h",
    "AV_tau_j",
    "AV_tau_m",
    "AV_aml",
    "AV_bml",
    "AV_hl_inf",
    "AV_INaL_INaL",
    "AV_ml_inf",
    "AV_tau_ml",
    "AV_PhiCaK",
    "AV_PhiCaL",
    "AV_PhiCaNa",
    "AV_f",
    "AV_tau_d",
    "AV_tau_fca",
    "AV_tau_ff",
    "AV_tau_fs",
    "AV_ICaL_ICaL",
    "AV_d_inf",
    "AV_f_inf",
    "AV_ICaK",
    "AV_ICaNa",
    "AV_fca_inf",
    "AV_ff_inf",
    "AV_fs_inf",
    "AV_rkr",
    "AV_tau_xr",
    "AV_xr_inf",
    "AV_IKr_IKr",
    "AV_KsCa",
    "AV_tau_xs1",
    "AV_tau_xs2",
    "AV_xs1_inf",
    "AV_IKs_IKs",
    "AV_xs2_inf",
    "AV_rk1",
    "AV_IK1_IK1",
    "AV_alpha_aa",
    "AV_beta_aa",
    "AV_rito2",
    "AV_diff_A",
    "AV_trel1factor",
    "AV_trel2factor",
    "AV_Jgap",
    "AV_diff_tjsrol",
    "AV_greljsrol",
    "AV_tau_Jrel1",
    "AV_tau_Jrel2",
    "AV_Jrelol",
    "AV_Jrel",
    "AV_allo",
    "AV_h4",
    "AV_hca",
    "AV_hna",
    "AV_h1",
    "AV_h5",
    "AV_h6",
    "AV_h7",
    "AV_h2",
    "AV_h3",
    "AV_h8",
    "AV_h9",
    "AV_k6",
    "AV_k3p",
    "AV_k3pp",
    "AV_k4p",
    "AV_k4pp",
    "AV_k7",
    "AV_k8",
    "AV_k3",
    "AV_k4",
    "AV_x1",
    "AV_x2",
    "AV_x3",
    "AV_x4",
    "AV_E1",
    "AV_E2",
    "AV_E3",
    "AV_E4",
    "AV_JncxCa",
    "AV_JncxNa",
    "AV_INaCa_i_INaCa_i",
    "AV_allo1",
    "AV_h111",
    "AV_h41",
    "AV_h71",
    "AV_h21",
    "AV_h31",
    "AV_h51",
    "AV_h61",
    "AV_h81",
    "AV_h91",
    "AV_k3p1",
    "AV_k3pp1",
    "AV_k4p1",
    "AV_k4pp1",
    "AV_k61",
    "AV_k71",
    "AV_k81",
    "AV_k31",
    "AV_k41",
    "AV_x11",
    "AV_x21",
    "AV_x31",
    "AV_x41",
    "AV_E11",
    "AV_E21",
    "AV_E31",
    "AV_E41",
    "AV_JncxCa1",
    "AV_JncxNa1",
    "AV_INaCa_ss_INaCa_ss",
    "AV_INaCa",
    "AV_Knai",
    "AV_Knao",
    "AV_P",
    "AV_a1",
    "AV_a3",
    "AV_b2",
    "AV_b3",
    "AV_b4",
    "AV_x12",
    "AV_x22",
    "AV_x32",
    "AV_x42",
    "AV_E12",
    "AV_E22",
    "AV_E32",
    "AV_E42",
    "AV_JnakK",
    "AV_JnakNa",
    "AV_INaK_INaK",
    "AV_xkb",
    "AV_IKb_IKb",
    "AV_INab_INab",
    "AV_ICab_ICab",
    "AV_IpCa_IpCa",
    "AV_Jleak",
    "AV_Jtr",
    "AV_Jtr2",
    "AV_Jupnp",
    "AV_Jupnp2",
    "AV_Jupp",
    "AV_Jupp2",
    "AV_fJupp",
    "AV_Jup",
    "AV_Jup2",
    "AV_Bcacsr",
    "AV_Bcai",
    "AV_Bcai2",
    "AV_Bcajsr",
    "AV_Bcass",
    "AV_cai_avg",
    "AV_diff_cansr",
    "AV_diff_ki",
    "AV_diff_kss",
    "AV_diff_nai",
    "AV_diff_nass",
    "AV_diff_cacsr",
    "AV_diff_cai",
    "AV_diff_cai2",
    "AV_diff_cajsr",
    "AV_diff_cass",
    "AV_kito2",
    "AV_Rel1",
    "AV_Rel2",
    "AV_ITo_ITo",
    "AV_Jrel1_inf",
    "AV_Jrel2_inf",
    "Istim",
    "Iion_cm",
};


static inline const Foam::HashTable<Foam::wordList>&
GaurDependencyMap()
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
GaurinitConsts(double* CONSTANTS, double* RATES, double* STATES, int tissueFlag, const Foam::dictionary& stimulus)
{

    /* cell */
    CONSTANTS[AC_F] = 96485.0;
    CONSTANTS[AC_L] = 0.017;
    CONSTANTS[AC_R] = 8314.0;
    CONSTANTS[AC_T] = 310.0;
    CONSTANTS[AC_cao] = 1.0 * 1.8;
    CONSTANTS[AC_cli] = 19.53;
    CONSTANTS[AC_clo] = 100.0;
    CONSTANTS[AC_ko] = 5.4;
    CONSTANTS[AC_nao] = 140.0;
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

    /* --- Initial values --- */

STATES[cell_v] = -87.3;
    STATES[nai] = 6.43;
    STATES[nass] = 6.43;
    STATES[ki] = 140.76;
    STATES[kss] = 140.76;
    STATES[cai] = 6.75e-05;
    STATES[cai2] = 6.81e-05;
    STATES[cass] = 6.75e-05;
    STATES[cansr] = 1.37;
    STATES[cajsr] = 1.31;
    STATES[cacsr] = 1.3;
    STATES[I_Na_m] = 0.0022;
    STATES[I_Na_h] = 0.631;
    STATES[I_Na_j] = 0.631;
    STATES[INaL_ml] = 0.001;
    STATES[INaL_hl] = 0.301;
    STATES[ICaL_d] = 4.42e-07;
    STATES[ICaL_fca] = 0.988;
    STATES[IKr_xr] = 0.153;
    STATES[IKs_xs1] = 0.054;
    STATES[IKs_xs2] = 0.046;
    STATES[ITo_aa] = 0.995;
    STATES[CICR_Jrel2] = 8.59e-10;
    STATES[CICR_Jrel1] = 9.71e-22;
    STATES[CaMK_CaMKt] =  8.23910999999999914e-03;
    STATES[CICR_tjsrol] = 100.0;
    STATES[CICR_A] = 100.0;
    STATES[ICaL_fs] = 0.0;
    STATES[ICaL_ff] = 0.0;
}


void
GaurcomputeVariables(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC, int tissueFlag, bool solveVmWithinODESolver)
{
/* diffusion */
    ALGEBRAIC[AV_Jdiff] = (STATES[cass] - STATES[cai]) / 0.2;
    ALGEBRAIC[AV_JdiffK] = (STATES[kss] - STATES[ki]) / 2.0;
    ALGEBRAIC[AV_JdiffNa] = (STATES[nass] - STATES[nai]) / 2.0;
    

    
    /* CaMK */
    ALGEBRAIC[AV_EK] = CONSTANTS[AC_R] * CONSTANTS[AC_T] / CONSTANTS[AC_F] * log(CONSTANTS[AC_ko] / STATES[ki]);
    ALGEBRAIC[AV_ENa] = CONSTANTS[AC_R] * CONSTANTS[AC_T] / CONSTANTS[AC_F] * log(CONSTANTS[AC_nao] / STATES[nai]);
    ALGEBRAIC[AV_vffrt] = STATES[cell_v] * CONSTANTS[AC_F] * CONSTANTS[AC_F] / (CONSTANTS[AC_R] * CONSTANTS[AC_T]);
    ALGEBRAIC[AV_vfrt] = STATES[cell_v] * CONSTANTS[AC_F] / (CONSTANTS[AC_R] * CONSTANTS[AC_T]);
    ALGEBRAIC[AV_CaMKb] = CONSTANTS[AC_CaMKo] * (1.0 - STATES[CaMK_CaMKt]) / (1.0 + CONSTANTS[AC_KmCaM] / STATES[cass]);
    ALGEBRAIC[AV_EKs] = CONSTANTS[AC_R] * CONSTANTS[AC_T] / CONSTANTS[AC_F] * log((CONSTANTS[AC_ko] + CONSTANTS[AC_PKNa] * CONSTANTS[AC_nao]) / (STATES[ki] + CONSTANTS[AC_PKNa] * STATES[nai]));
    ALGEBRAIC[AV_CaMKa] = ALGEBRAIC[AV_CaMKb] + STATES[CaMK_CaMKt];
    ALGEBRAIC[AV_diff_CaMKt] = CONSTANTS[AC_aCaMK] * ALGEBRAIC[AV_CaMKb] * (ALGEBRAIC[AV_CaMKb] + STATES[CaMK_CaMKt]) - CONSTANTS[AC_bCaMK] * STATES[CaMK_CaMKt];
    ALGEBRAIC[AV_CaMK_f] = 1.0 / (1.0 + CONSTANTS[AC_KmCaMK] / ALGEBRAIC[AV_CaMKa]);
    RATES[CaMK_CaMKt] = ALGEBRAIC[AV_diff_CaMKt];
    
    /* I_Na */
    ALGEBRAIC[AV_aa_h] = ((STATES[cell_v] >= -40.0) ? 0.0 : 0.057 * exp(-(STATES[cell_v] + 80.0) / 6.8));
    ALGEBRAIC[AV_aa_j] = ((STATES[cell_v] >= -40.0) ? 0.0 : (-25428.0 * exp(0.2444 * STATES[cell_v]) - 6.948e-06 * exp(-0.04391 * STATES[cell_v])) * (STATES[cell_v] + 37.78) / (1.0 + exp(0.311 * (STATES[cell_v] + 79.23))));
    ALGEBRAIC[AV_aa_m] = 1.0 / (1.0 + exp((-60.0 - STATES[cell_v]) / 5.0));
    ALGEBRAIC[AV_bb_h] = ((STATES[cell_v] >= -40.0) ? 0.77 / (0.13 * (1.0 + exp(-(STATES[cell_v] + 10.66) / 11.1))) : 2.7 * exp(0.079 * STATES[cell_v]) + 310000.0 * exp(0.3485 * STATES[cell_v]));
    ALGEBRAIC[AV_bb_j] = ((STATES[cell_v] >= -40.0) ? 0.6 * exp(0.057 * STATES[cell_v]) / (1.0 + exp(-0.1 * (STATES[cell_v] + 32.0))) : 0.02424 * exp(-0.01052 * STATES[cell_v]) / (1.0 + exp(-0.1378 * (STATES[cell_v] + 40.14))));
    ALGEBRAIC[AV_bb_m] = 0.1 / (1.0 + exp((STATES[cell_v] + 35.0) / 5.0)) + 0.1 / (1.0 + exp((STATES[cell_v] - 50.0) / 200.0));
    ALGEBRAIC[AV_h_inf] = 1.0 / ((1.0 + exp((STATES[cell_v] + 71.55) / 7.43)) * (1.0 + exp((STATES[cell_v] + 71.55) / 7.43)));
    ALGEBRAIC[AV_m_inf] = 1.0 / ((1.0 + exp((-56.86 - STATES[cell_v]) / 9.03)) * (1.0 + exp((-56.86 - STATES[cell_v]) / 9.03)));
    ALGEBRAIC[AV_INa] = CONSTANTS[AC_GNa] * STATES[I_Na_m] * STATES[I_Na_m] * STATES[I_Na_m] * STATES[I_Na_h] * STATES[I_Na_j] * (STATES[cell_v] - ALGEBRAIC[AV_ENa]);
    ALGEBRAIC[AV_j_inf] = ALGEBRAIC[AV_h_inf];
    ALGEBRAIC[AV_tau_h] = 1.0 / (ALGEBRAIC[AV_aa_h] + ALGEBRAIC[AV_bb_h]);
    ALGEBRAIC[AV_tau_j] = 1.0 / (ALGEBRAIC[AV_aa_j] + ALGEBRAIC[AV_bb_j]);
    ALGEBRAIC[AV_tau_m] = ALGEBRAIC[AV_aa_m] * ALGEBRAIC[AV_bb_m];
    RATES[I_Na_h] = (ALGEBRAIC[AV_h_inf] - STATES[I_Na_h]) / ALGEBRAIC[AV_tau_h];
    RATES[I_Na_j] = (ALGEBRAIC[AV_j_inf] - STATES[I_Na_j]) / ALGEBRAIC[AV_tau_j];
    RATES[I_Na_m] = (ALGEBRAIC[AV_m_inf] - STATES[I_Na_m]) / ALGEBRAIC[AV_tau_m];
    
    /* INaL */
    ALGEBRAIC[AV_aml] = 0.32 * ((STATES[cell_v] == -47.13) ? 10.0 : -(STATES[cell_v] + 47.13) / (exp(-0.1 * (STATES[cell_v] + 47.13)) - 1.0));
    ALGEBRAIC[AV_bml] = 0.08 * exp(-STATES[cell_v] / 11.0);
    ALGEBRAIC[AV_hl_inf] = 1.0 / (1.0 + exp((STATES[cell_v] + 91.0) / 6.1));
    ALGEBRAIC[AV_INaL_INaL] = CONSTANTS[AC_GNaL] * STATES[INaL_ml] * STATES[INaL_ml] * STATES[INaL_ml] * STATES[INaL_hl] * (STATES[cell_v] - ALGEBRAIC[AV_ENa]);
    ALGEBRAIC[AV_ml_inf] = ALGEBRAIC[AV_aml] / (ALGEBRAIC[AV_aml] + ALGEBRAIC[AV_bml]);
    ALGEBRAIC[AV_tau_ml] = 1.0 / (ALGEBRAIC[AV_aml] + ALGEBRAIC[AV_bml]);
    RATES[INaL_hl] = (ALGEBRAIC[AV_hl_inf] - STATES[INaL_hl]) / CONSTANTS[AC_tau_hl];
    RATES[INaL_ml] = (ALGEBRAIC[AV_ml_inf] - STATES[INaL_ml]) / ALGEBRAIC[AV_tau_ml];
    
    /* ICaL */
    ALGEBRAIC[AV_PhiCaK] = CONSTANTS[AC_F] * (0.75 * STATES[kss] * exp(ALGEBRAIC[AV_vfrt]) - 0.75 * CONSTANTS[AC_ko]) * ((STATES[cell_v] == 0.0) ? 1.0 : ALGEBRAIC[AV_vfrt] / (exp(ALGEBRAIC[AV_vfrt]) - 1.0));
    ALGEBRAIC[AV_PhiCaL] = 2.0 * CONSTANTS[AC_F] * (STATES[cass] * exp(2.0 * (STATES[cell_v] * CONSTANTS[AC_F] / (CONSTANTS[AC_R] * CONSTANTS[AC_T]))) - 0.341 * CONSTANTS[AC_cao]) * ((STATES[cell_v] == 0.0) ? 1.0 : 2.0 * ALGEBRAIC[AV_vfrt] / (exp(2.0 * ALGEBRAIC[AV_vfrt]) - 1.0));
    ALGEBRAIC[AV_PhiCaNa] = CONSTANTS[AC_F] * (0.75 * STATES[nass] * exp(ALGEBRAIC[AV_vfrt]) - 0.75 * CONSTANTS[AC_nao]) * ((STATES[cell_v] == 0.0) ? 1.0 : ALGEBRAIC[AV_vfrt] / (exp(ALGEBRAIC[AV_vfrt]) - 1.0));
    ALGEBRAIC[AV_f] = STATES[ICaL_ff] * STATES[ICaL_fs];
    ALGEBRAIC[AV_tau_d] = 0.6 + 1.0 / (exp(-0.05 * (STATES[cell_v] + 6.0)) + exp(0.09 * (STATES[cell_v] + 14.0)));
    ALGEBRAIC[AV_tau_fca] = 10.0 * ALGEBRAIC[AV_CaMK_f] + 0.5 + 1.0 / (1.0 + STATES[cass] / 0.003);
    ALGEBRAIC[AV_tau_ff] = 7.0 + 1.0 / (0.0045 * exp(-(STATES[cell_v] + 20.0) / 10.0) + 0.0045 * exp((STATES[cell_v] + 20.0) / 10.0));
    ALGEBRAIC[AV_tau_fs] = 70.0 + 1.0 / (3.5e-05 * exp(-(STATES[cell_v] + 5.0) / 4.0) + 3.5e-05 * exp((STATES[cell_v] + 5.0) / 6.0));
    ALGEBRAIC[AV_ICaL_ICaL] = CONSTANTS[AC_PCa] * ALGEBRAIC[AV_PhiCaL] * STATES[ICaL_d] * ALGEBRAIC[AV_f] * STATES[ICaL_fca];
    ALGEBRAIC[AV_d_inf] = 1.0 / (1.0 + exp(-(STATES[cell_v] - CONSTANTS[AC_vhalf_d]) / 6.2));
    ALGEBRAIC[AV_f_inf] = 1.0 / (1.0 + exp((STATES[cell_v] - CONSTANTS[AC_vhalff]) / 4.9)) + 0.35 / (1.0 + exp((45.0 - STATES[cell_v]) / 20.0));
    ALGEBRAIC[AV_ICaK] = CONSTANTS[AC_PCaK] * ALGEBRAIC[AV_PhiCaK] * STATES[ICaL_d] * ALGEBRAIC[AV_f] * STATES[ICaL_fca];
    ALGEBRAIC[AV_ICaNa] = CONSTANTS[AC_PCaNa] * ALGEBRAIC[AV_PhiCaNa] * STATES[ICaL_d] * ALGEBRAIC[AV_f] * STATES[ICaL_fca];
    ALGEBRAIC[AV_fca_inf] = 0.3 / (1.0 - ((ALGEBRAIC[AV_ICaL_ICaL] > 0.0) ? 0.0 : ALGEBRAIC[AV_ICaL_ICaL] / 0.05)) + 0.55 / (1.0 + STATES[cass] / 0.003) + 0.15;
    ALGEBRAIC[AV_ff_inf] = ALGEBRAIC[AV_f_inf];
    ALGEBRAIC[AV_fs_inf] = ALGEBRAIC[AV_f_inf];
    RATES[ICaL_d] = (ALGEBRAIC[AV_d_inf] - STATES[ICaL_d]) / ALGEBRAIC[AV_tau_d];
    RATES[ICaL_fca] = (ALGEBRAIC[AV_fca_inf] - STATES[ICaL_fca]) / ALGEBRAIC[AV_tau_fca];
    RATES[ICaL_ff] = (ALGEBRAIC[AV_ff_inf] - STATES[ICaL_ff]) / ALGEBRAIC[AV_tau_ff];
    RATES[ICaL_fs] = (ALGEBRAIC[AV_fs_inf] - STATES[ICaL_fs]) / ALGEBRAIC[AV_tau_fs];
    
    /* IKr */
    ALGEBRAIC[AV_rkr] = 1.0 / (1.0 + exp((STATES[cell_v] + 22.0) / 15.0));
    ALGEBRAIC[AV_tau_xr] = 12.98 + 1.0 / (0.3652 * exp((STATES[cell_v] - 31.66) / 3.869) + 4.123e-05 * exp(-(STATES[cell_v] - 47.78) / 20.38));
    ALGEBRAIC[AV_xr_inf] = 1.0 / (1.0 + exp(-(STATES[cell_v] + 56.8) / 17.8));
    ALGEBRAIC[AV_IKr_IKr] = CONSTANTS[AC_GKr] * sqrt(CONSTANTS[AC_ko] / 5.4) * STATES[IKr_xr] * ALGEBRAIC[AV_rkr] * (STATES[cell_v] - ALGEBRAIC[AV_EK]);
    RATES[IKr_xr] = (ALGEBRAIC[AV_xr_inf] - STATES[IKr_xr]) / ALGEBRAIC[AV_tau_xr];
    
    /* IKs */
    ALGEBRAIC[AV_KsCa] = 1.0 + 0.6 / (1.0 + pow(3.8e-05 / STATES[cai], 1.4));
    ALGEBRAIC[AV_tau_xs1] = 300.0 + 1.0 / (1e-06 * exp((STATES[cell_v] + 50.0) / 20.0) + 0.04 * exp(-(STATES[cell_v] + 50.0) / 20.0));
    ALGEBRAIC[AV_tau_xs2] = 1.0 / (0.01 * exp((STATES[cell_v] - 50.0) / 20.0) + 0.0193 * exp(-(STATES[cell_v] + 66.54) / 31.0));
    ALGEBRAIC[AV_xs1_inf] = 1.0 / (1.0 + exp(-(STATES[cell_v] - 25.1) / 37.1));
    ALGEBRAIC[AV_IKs_IKs] = CONSTANTS[AC_GKs] * ALGEBRAIC[AV_KsCa] * STATES[IKs_xs1] * STATES[IKs_xs2] * (STATES[cell_v] - ALGEBRAIC[AV_EKs]);
    ALGEBRAIC[AV_xs2_inf] = ALGEBRAIC[AV_xs1_inf];
    RATES[IKs_xs1] = (ALGEBRAIC[AV_xs1_inf] - STATES[IKs_xs1]) / ALGEBRAIC[AV_tau_xs1];
    RATES[IKs_xs2] = (ALGEBRAIC[AV_xs2_inf] - STATES[IKs_xs2]) / ALGEBRAIC[AV_tau_xs2];
    
    /* IK1 */
    ALGEBRAIC[AV_rk1] = 1.0 / (1.0 + exp((STATES[cell_v] + 79.3 - 2.6 * CONSTANTS[AC_ko]) / 19.6));
    ALGEBRAIC[AV_IK1_IK1] = CONSTANTS[AC_GK1] * sqrt(CONSTANTS[AC_ko] / 5.4) * ALGEBRAIC[AV_rk1] * (STATES[cell_v] - ALGEBRAIC[AV_EK]);
    
    /* ITo */
    ALGEBRAIC[AV_alpha_aa] = 0.025 / (1.0 + exp((STATES[cell_v] + 58.0) / 5.0));
    ALGEBRAIC[AV_beta_aa] = 1.0 / (5.0 * (1.0 + exp((STATES[cell_v] + 19.0) / -9.0)));
    ALGEBRAIC[AV_rito2] = 1.0 / (1.0 + exp(-(STATES[cell_v] + 10.0) / 5.0));
    RATES[ITo_aa] = ALGEBRAIC[AV_alpha_aa] * (1.0 - STATES[ITo_aa]) - ALGEBRAIC[AV_beta_aa] * STATES[ITo_aa];
    
    /* CICR */
    ALGEBRAIC[AV_diff_A] = ((STATES[CICR_tjsrol] < 5.0) ? -STATES[CICR_A] / 0.1 : 1.0);
    ALGEBRAIC[AV_trel1factor] = 1.0 * 20.0 * (1.0 + 1.0 / (1.0 + pow(CONSTANTS[AC_KmCaMK] / ALGEBRAIC[AV_CaMKa], 8.0))) / (1.0 + pow(0.5 / STATES[cajsr], 8.0));
    ALGEBRAIC[AV_trel2factor] = 50.0 * (1.0 + 1.0 / (1.0 + pow(CONSTANTS[AC_KmCaMK] / ALGEBRAIC[AV_CaMKa], 8.0))) / (1.0 + pow(0.5 / STATES[cacsr], 8.0));
    ALGEBRAIC[AV_Jgap] = (STATES[cai] - STATES[cai2]) / CONSTANTS[AC_tau_gap];
    ALGEBRAIC[AV_diff_tjsrol] = (((STATES[cajsr] > CONSTANTS[AC_SOICR]) && (STATES[CICR_A] > 45.0)) ? -STATES[CICR_tjsrol] / 0.001 : 1.0);
    ALGEBRAIC[AV_greljsrol] = CONSTANTS[AC_grelbarjsrol] * (1.0 - exp(-STATES[CICR_tjsrol] / CONSTANTS[AC_tauon])) * exp(-STATES[CICR_tjsrol] / CONSTANTS[AC_tauoff]);
    ALGEBRAIC[AV_tau_Jrel1] = ((0.001 > ALGEBRAIC[AV_trel1factor]) ? 0.001 : ALGEBRAIC[AV_trel1factor]);
    ALGEBRAIC[AV_tau_Jrel2] = ((0.001 > ALGEBRAIC[AV_trel2factor]) ? 0.001 : ALGEBRAIC[AV_trel2factor]);
    RATES[CICR_A] = ALGEBRAIC[AV_diff_A];
    ALGEBRAIC[AV_Jrelol] = ALGEBRAIC[AV_greljsrol] * (STATES[cajsr] - STATES[cass]);
    RATES[CICR_tjsrol] = ALGEBRAIC[AV_diff_tjsrol];
    ALGEBRAIC[AV_Jrel] = STATES[CICR_Jrel1] + ALGEBRAIC[AV_Jrelol];
    
    /* INaCa_i */
    ALGEBRAIC[AV_allo] = 1.0 / (1.0 + CONSTANTS[AC_KmCaAct] / STATES[cai] * CONSTANTS[AC_KmCaAct] / STATES[cai]);
    ALGEBRAIC[AV_h4] = 1.0 + STATES[nai] / CONSTANTS[AC_kna1] * (1.0 + STATES[nai] / CONSTANTS[AC_kna2]);
    ALGEBRAIC[AV_hca] = exp(CONSTANTS[AC_qca] * STATES[cell_v] * CONSTANTS[AC_F] / (CONSTANTS[AC_R] * CONSTANTS[AC_T]));
    ALGEBRAIC[AV_hna] = exp(CONSTANTS[AC_qna] * STATES[cell_v] * CONSTANTS[AC_F] / (CONSTANTS[AC_R] * CONSTANTS[AC_T]));
    ALGEBRAIC[AV_h1] = 1.0 + STATES[nai] / CONSTANTS[AC_kna3] * (1.0 + ALGEBRAIC[AV_hna]);
    ALGEBRAIC[AV_h5] = STATES[nai] * STATES[nai] / (ALGEBRAIC[AV_h4] * CONSTANTS[AC_kna1] * CONSTANTS[AC_kna2]);
    ALGEBRAIC[AV_h6] = 1.0 / ALGEBRAIC[AV_h4];
    ALGEBRAIC[AV_h7] = 1.0 + CONSTANTS[AC_nao] / CONSTANTS[AC_kna3] * (1.0 + 1.0 / ALGEBRAIC[AV_hna]);
    ALGEBRAIC[AV_h2] = STATES[nai] * ALGEBRAIC[AV_hna] / (CONSTANTS[AC_kna3] * ALGEBRAIC[AV_h1]);
    ALGEBRAIC[AV_h3] = 1.0 / ALGEBRAIC[AV_h1];
    ALGEBRAIC[AV_h8] = CONSTANTS[AC_nao] / (CONSTANTS[AC_kna3] * ALGEBRAIC[AV_hna] * ALGEBRAIC[AV_h7]);
    ALGEBRAIC[AV_h9] = 1.0 / ALGEBRAIC[AV_h7];
    ALGEBRAIC[AV_k6] = ALGEBRAIC[AV_h6] * STATES[cai] * CONSTANTS[AC_kcaon];
    ALGEBRAIC[AV_k3p] = ALGEBRAIC[AV_h9] * CONSTANTS[AC_wca];
    ALGEBRAIC[AV_k3pp] = ALGEBRAIC[AV_h8] * CONSTANTS[AC_wnaca];
    ALGEBRAIC[AV_k4p] = ALGEBRAIC[AV_h3] * CONSTANTS[AC_wca] / ALGEBRAIC[AV_hca];
    ALGEBRAIC[AV_k4pp] = ALGEBRAIC[AV_h2] * CONSTANTS[AC_wnaca];
    ALGEBRAIC[AV_k7] = ALGEBRAIC[AV_h5] * ALGEBRAIC[AV_h2] * CONSTANTS[AC_wna];
    ALGEBRAIC[AV_k8] = ALGEBRAIC[AV_h8] * CONSTANTS[AC_h11] * CONSTANTS[AC_wna];
    ALGEBRAIC[AV_k3] = ALGEBRAIC[AV_k3p] + ALGEBRAIC[AV_k3pp];
    ALGEBRAIC[AV_k4] = ALGEBRAIC[AV_k4p] + ALGEBRAIC[AV_k4pp];
    ALGEBRAIC[AV_x1] = CONSTANTS[AC_k2] * ALGEBRAIC[AV_k4] * (ALGEBRAIC[AV_k7] + ALGEBRAIC[AV_k6]) + CONSTANTS[AC_k5] * ALGEBRAIC[AV_k7] * (CONSTANTS[AC_k2] + ALGEBRAIC[AV_k3]);
    ALGEBRAIC[AV_x2] = CONSTANTS[AC_k1] * ALGEBRAIC[AV_k7] * (ALGEBRAIC[AV_k4] + CONSTANTS[AC_k5]) + ALGEBRAIC[AV_k4] * ALGEBRAIC[AV_k6] * (CONSTANTS[AC_k1] + ALGEBRAIC[AV_k8]);
    ALGEBRAIC[AV_x3] = CONSTANTS[AC_k1] * ALGEBRAIC[AV_k3] * (ALGEBRAIC[AV_k7] + ALGEBRAIC[AV_k6]) + ALGEBRAIC[AV_k8] * ALGEBRAIC[AV_k6] * (CONSTANTS[AC_k2] + ALGEBRAIC[AV_k3]);
    ALGEBRAIC[AV_x4] = CONSTANTS[AC_k2] * ALGEBRAIC[AV_k8] * (ALGEBRAIC[AV_k4] + CONSTANTS[AC_k5]) + ALGEBRAIC[AV_k3] * CONSTANTS[AC_k5] * (CONSTANTS[AC_k1] + ALGEBRAIC[AV_k8]);
    ALGEBRAIC[AV_E1] = ALGEBRAIC[AV_x1] / (ALGEBRAIC[AV_x1] + ALGEBRAIC[AV_x2] + ALGEBRAIC[AV_x3] + ALGEBRAIC[AV_x4]);
    ALGEBRAIC[AV_E2] = ALGEBRAIC[AV_x2] / (ALGEBRAIC[AV_x1] + ALGEBRAIC[AV_x2] + ALGEBRAIC[AV_x3] + ALGEBRAIC[AV_x4]);
    ALGEBRAIC[AV_E3] = ALGEBRAIC[AV_x3] / (ALGEBRAIC[AV_x1] + ALGEBRAIC[AV_x2] + ALGEBRAIC[AV_x3] + ALGEBRAIC[AV_x4]);
    ALGEBRAIC[AV_E4] = ALGEBRAIC[AV_x4] / (ALGEBRAIC[AV_x1] + ALGEBRAIC[AV_x2] + ALGEBRAIC[AV_x3] + ALGEBRAIC[AV_x4]);
    ALGEBRAIC[AV_JncxCa] = ALGEBRAIC[AV_E2] * CONSTANTS[AC_k2] - ALGEBRAIC[AV_E1] * CONSTANTS[AC_k1];
    ALGEBRAIC[AV_JncxNa] = 3.0 * (ALGEBRAIC[AV_E4] * ALGEBRAIC[AV_k7] - ALGEBRAIC[AV_E1] * ALGEBRAIC[AV_k8]) + ALGEBRAIC[AV_E3] * ALGEBRAIC[AV_k4pp] - ALGEBRAIC[AV_E2] * ALGEBRAIC[AV_k3pp];
    ALGEBRAIC[AV_INaCa_i_INaCa_i] = 0.8 * CONSTANTS[AC_Gncx] * ALGEBRAIC[AV_allo] * (CONSTANTS[AC_zna] * ALGEBRAIC[AV_JncxNa] + CONSTANTS[AC_zca] * ALGEBRAIC[AV_JncxCa]);
    
    /* INaCa_ss */
    ALGEBRAIC[AV_allo1] = 1.0 / (1.0 + CONSTANTS[AC_KmCaAct] / STATES[cass] * CONSTANTS[AC_KmCaAct] / STATES[cass]);
    ALGEBRAIC[AV_h111] = 1.0 + STATES[nass] / CONSTANTS[AC_kna3] * (1.0 + ALGEBRAIC[AV_hna]);
    ALGEBRAIC[AV_h41] = 1.0 + STATES[nass] / CONSTANTS[AC_kna1] * (1.0 + STATES[nass] / CONSTANTS[AC_kna2]);
    ALGEBRAIC[AV_h71] = 1.0 + CONSTANTS[AC_nao] / CONSTANTS[AC_kna3] * (1.0 + 1.0 / ALGEBRAIC[AV_hna]);
    ALGEBRAIC[AV_h21] = STATES[nass] * ALGEBRAIC[AV_hna] / (CONSTANTS[AC_kna3] * ALGEBRAIC[AV_h111]);
    ALGEBRAIC[AV_h31] = 1.0 / ALGEBRAIC[AV_h111];
    ALGEBRAIC[AV_h51] = STATES[nass] * STATES[nass] / (ALGEBRAIC[AV_h41] * CONSTANTS[AC_kna1] * CONSTANTS[AC_kna2]);
    ALGEBRAIC[AV_h61] = 1.0 / ALGEBRAIC[AV_h41];
    ALGEBRAIC[AV_h81] = CONSTANTS[AC_nao] / (CONSTANTS[AC_kna3] * ALGEBRAIC[AV_hna] * ALGEBRAIC[AV_h71]);
    ALGEBRAIC[AV_h91] = 1.0 / ALGEBRAIC[AV_h71];
    ALGEBRAIC[AV_k3p1] = ALGEBRAIC[AV_h91] * CONSTANTS[AC_wca];
    ALGEBRAIC[AV_k3pp1] = ALGEBRAIC[AV_h81] * CONSTANTS[AC_wnaca];
    ALGEBRAIC[AV_k4p1] = ALGEBRAIC[AV_h31] * CONSTANTS[AC_wca] / ALGEBRAIC[AV_hca];
    ALGEBRAIC[AV_k4pp1] = ALGEBRAIC[AV_h21] * CONSTANTS[AC_wnaca];
    ALGEBRAIC[AV_k61] = ALGEBRAIC[AV_h61] * STATES[cass] * CONSTANTS[AC_kcaon];
    ALGEBRAIC[AV_k71] = ALGEBRAIC[AV_h51] * ALGEBRAIC[AV_h21] * CONSTANTS[AC_wna];
    ALGEBRAIC[AV_k81] = ALGEBRAIC[AV_h81] * CONSTANTS[AC_h1111] * CONSTANTS[AC_wna];
    ALGEBRAIC[AV_k31] = ALGEBRAIC[AV_k3p1] + ALGEBRAIC[AV_k3pp1];
    ALGEBRAIC[AV_k41] = ALGEBRAIC[AV_k4p1] + ALGEBRAIC[AV_k4pp1];
    ALGEBRAIC[AV_x11] = CONSTANTS[AC_k21] * ALGEBRAIC[AV_k41] * (ALGEBRAIC[AV_k71] + ALGEBRAIC[AV_k61]) + CONSTANTS[AC_k51] * ALGEBRAIC[AV_k71] * (CONSTANTS[AC_k21] + ALGEBRAIC[AV_k31]);
    ALGEBRAIC[AV_x21] = CONSTANTS[AC_k11] * ALGEBRAIC[AV_k71] * (ALGEBRAIC[AV_k41] + CONSTANTS[AC_k51]) + ALGEBRAIC[AV_k41] * ALGEBRAIC[AV_k61] * (CONSTANTS[AC_k11] + ALGEBRAIC[AV_k81]);
    ALGEBRAIC[AV_x31] = CONSTANTS[AC_k11] * ALGEBRAIC[AV_k31] * (ALGEBRAIC[AV_k71] + ALGEBRAIC[AV_k61]) + ALGEBRAIC[AV_k81] * ALGEBRAIC[AV_k61] * (CONSTANTS[AC_k21] + ALGEBRAIC[AV_k31]);
    ALGEBRAIC[AV_x41] = CONSTANTS[AC_k21] * ALGEBRAIC[AV_k81] * (ALGEBRAIC[AV_k41] + CONSTANTS[AC_k51]) + ALGEBRAIC[AV_k31] * CONSTANTS[AC_k51] * (CONSTANTS[AC_k11] + ALGEBRAIC[AV_k81]);
    ALGEBRAIC[AV_E11] = ALGEBRAIC[AV_x11] / (ALGEBRAIC[AV_x11] + ALGEBRAIC[AV_x21] + ALGEBRAIC[AV_x31] + ALGEBRAIC[AV_x41]);
    ALGEBRAIC[AV_E21] = ALGEBRAIC[AV_x21] / (ALGEBRAIC[AV_x11] + ALGEBRAIC[AV_x21] + ALGEBRAIC[AV_x31] + ALGEBRAIC[AV_x41]);
    ALGEBRAIC[AV_E31] = ALGEBRAIC[AV_x31] / (ALGEBRAIC[AV_x11] + ALGEBRAIC[AV_x21] + ALGEBRAIC[AV_x31] + ALGEBRAIC[AV_x41]);
    ALGEBRAIC[AV_E41] = ALGEBRAIC[AV_x41] / (ALGEBRAIC[AV_x11] + ALGEBRAIC[AV_x21] + ALGEBRAIC[AV_x31] + ALGEBRAIC[AV_x41]);
    ALGEBRAIC[AV_JncxCa1] = ALGEBRAIC[AV_E21] * CONSTANTS[AC_k21] - ALGEBRAIC[AV_E11] * CONSTANTS[AC_k11];
    ALGEBRAIC[AV_JncxNa1] = 3.0 * (ALGEBRAIC[AV_E41] * ALGEBRAIC[AV_k71] - ALGEBRAIC[AV_E11] * ALGEBRAIC[AV_k81]) + ALGEBRAIC[AV_E31] * ALGEBRAIC[AV_k4pp1] - ALGEBRAIC[AV_E21] * ALGEBRAIC[AV_k3pp1];
    ALGEBRAIC[AV_INaCa_ss_INaCa_ss] = 0.2 * CONSTANTS[AC_Gncx] * ALGEBRAIC[AV_allo1] * (CONSTANTS[AC_zna] * ALGEBRAIC[AV_JncxNa1] + CONSTANTS[AC_zca] * ALGEBRAIC[AV_JncxCa1]);
    ALGEBRAIC[AV_INaCa] = ALGEBRAIC[AV_INaCa_i_INaCa_i] + ALGEBRAIC[AV_INaCa_ss_INaCa_ss];
    
    /* INaK */
    ALGEBRAIC[AV_Knai] = CONSTANTS[AC_Knai0] * exp(CONSTANTS[AC_delta] * STATES[cell_v] * CONSTANTS[AC_F] / (3.0 * CONSTANTS[AC_R] * CONSTANTS[AC_T]));
    ALGEBRAIC[AV_Knao] = CONSTANTS[AC_Knao0] * exp((1.0 - CONSTANTS[AC_delta]) * STATES[cell_v] * CONSTANTS[AC_F] / (3.0 * CONSTANTS[AC_R] * CONSTANTS[AC_T]));
    ALGEBRAIC[AV_P] = CONSTANTS[AC_eP] / (1.0 + CONSTANTS[AC_H] / CONSTANTS[AC_Khp] + STATES[nai] / CONSTANTS[AC_Knap] + STATES[ki] / CONSTANTS[AC_Kxkur]);
    ALGEBRAIC[AV_a1] = CONSTANTS[AC_k1p] * (STATES[nai] / ALGEBRAIC[AV_Knai] * STATES[nai] / ALGEBRAIC[AV_Knai] * (STATES[nai] / ALGEBRAIC[AV_Knai])) / ((1.0 + STATES[nai] / ALGEBRAIC[AV_Knai]) * (1.0 + STATES[nai] / ALGEBRAIC[AV_Knai]) * (1.0 + STATES[nai] / ALGEBRAIC[AV_Knai]) + (1.0 + STATES[ki] / CONSTANTS[AC_Kki]) * (1.0 + STATES[ki] / CONSTANTS[AC_Kki]) - 1.0);
    ALGEBRAIC[AV_a3] = CONSTANTS[AC_k3p2] * (CONSTANTS[AC_ko] / CONSTANTS[AC_Kko] * CONSTANTS[AC_ko] / CONSTANTS[AC_Kko]) / ((1.0 + CONSTANTS[AC_nao] / ALGEBRAIC[AV_Knao]) * (1.0 + CONSTANTS[AC_nao] / ALGEBRAIC[AV_Knao]) * (1.0 + CONSTANTS[AC_nao] / ALGEBRAIC[AV_Knao]) + (1.0 + CONSTANTS[AC_ko] / CONSTANTS[AC_Kko]) * (1.0 + CONSTANTS[AC_ko] / CONSTANTS[AC_Kko]) - 1.0);
    ALGEBRAIC[AV_b2] = CONSTANTS[AC_k2m] * (CONSTANTS[AC_nao] / ALGEBRAIC[AV_Knao] * CONSTANTS[AC_nao] / ALGEBRAIC[AV_Knao] * (CONSTANTS[AC_nao] / ALGEBRAIC[AV_Knao])) / ((1.0 + CONSTANTS[AC_nao] / ALGEBRAIC[AV_Knao]) * (1.0 + CONSTANTS[AC_nao] / ALGEBRAIC[AV_Knao]) * (1.0 + CONSTANTS[AC_nao] / ALGEBRAIC[AV_Knao]) + (1.0 + CONSTANTS[AC_ko] / CONSTANTS[AC_Kko]) * (1.0 + CONSTANTS[AC_ko] / CONSTANTS[AC_Kko]) - 1.0);
    ALGEBRAIC[AV_b3] = CONSTANTS[AC_k3m] * ALGEBRAIC[AV_P] * CONSTANTS[AC_H] / (1.0 + CONSTANTS[AC_MgATP] / CONSTANTS[AC_Kmgatp]);
    ALGEBRAIC[AV_b4] = CONSTANTS[AC_k4m] * (STATES[ki] / CONSTANTS[AC_Kki] * STATES[ki] / CONSTANTS[AC_Kki]) / ((1.0 + STATES[nai] / ALGEBRAIC[AV_Knai]) * (1.0 + STATES[nai] / ALGEBRAIC[AV_Knai]) * (1.0 + STATES[nai] / ALGEBRAIC[AV_Knai]) + (1.0 + STATES[ki] / CONSTANTS[AC_Kki]) * (1.0 + STATES[ki] / CONSTANTS[AC_Kki]) - 1.0);
    ALGEBRAIC[AV_x12] = CONSTANTS[AC_a4] * ALGEBRAIC[AV_a1] * CONSTANTS[AC_a2] + ALGEBRAIC[AV_b2] * ALGEBRAIC[AV_b4] * ALGEBRAIC[AV_b3] + CONSTANTS[AC_a2] * ALGEBRAIC[AV_b4] * ALGEBRAIC[AV_b3] + ALGEBRAIC[AV_b3] * ALGEBRAIC[AV_a1] * CONSTANTS[AC_a2];
    ALGEBRAIC[AV_x22] = ALGEBRAIC[AV_b2] * CONSTANTS[AC_b1] * ALGEBRAIC[AV_b4] + ALGEBRAIC[AV_a1] * CONSTANTS[AC_a2] * ALGEBRAIC[AV_a3] + ALGEBRAIC[AV_a3] * CONSTANTS[AC_b1] * ALGEBRAIC[AV_b4] + CONSTANTS[AC_a2] * ALGEBRAIC[AV_a3] * ALGEBRAIC[AV_b4];
    ALGEBRAIC[AV_x32] = CONSTANTS[AC_a2] * ALGEBRAIC[AV_a3] * CONSTANTS[AC_a4] + ALGEBRAIC[AV_b3] * ALGEBRAIC[AV_b2] * CONSTANTS[AC_b1] + ALGEBRAIC[AV_b2] * CONSTANTS[AC_b1] * CONSTANTS[AC_a4] + ALGEBRAIC[AV_a3] * CONSTANTS[AC_a4] * CONSTANTS[AC_b1];
    ALGEBRAIC[AV_x42] = ALGEBRAIC[AV_b4] * ALGEBRAIC[AV_b3] * ALGEBRAIC[AV_b2] + ALGEBRAIC[AV_a3] * CONSTANTS[AC_a4] * ALGEBRAIC[AV_a1] + ALGEBRAIC[AV_b2] * CONSTANTS[AC_a4] * ALGEBRAIC[AV_a1] + ALGEBRAIC[AV_b3] * ALGEBRAIC[AV_b2] * ALGEBRAIC[AV_a1];
    ALGEBRAIC[AV_E12] = ALGEBRAIC[AV_x12] / (ALGEBRAIC[AV_x12] + ALGEBRAIC[AV_x22] + ALGEBRAIC[AV_x32] + ALGEBRAIC[AV_x42]);
    ALGEBRAIC[AV_E22] = ALGEBRAIC[AV_x22] / (ALGEBRAIC[AV_x12] + ALGEBRAIC[AV_x22] + ALGEBRAIC[AV_x32] + ALGEBRAIC[AV_x42]);
    ALGEBRAIC[AV_E32] = ALGEBRAIC[AV_x32] / (ALGEBRAIC[AV_x12] + ALGEBRAIC[AV_x22] + ALGEBRAIC[AV_x32] + ALGEBRAIC[AV_x42]);
    ALGEBRAIC[AV_E42] = ALGEBRAIC[AV_x42] / (ALGEBRAIC[AV_x12] + ALGEBRAIC[AV_x22] + ALGEBRAIC[AV_x32] + ALGEBRAIC[AV_x42]);
    ALGEBRAIC[AV_JnakK] = 2.0 * (ALGEBRAIC[AV_E42] * CONSTANTS[AC_b1] - ALGEBRAIC[AV_E32] * ALGEBRAIC[AV_a1]);
    ALGEBRAIC[AV_JnakNa] = 3.0 * (ALGEBRAIC[AV_E12] * ALGEBRAIC[AV_a3] - ALGEBRAIC[AV_E22] * ALGEBRAIC[AV_b3]);
    ALGEBRAIC[AV_INaK_INaK] = CONSTANTS[AC_Pnak] * (CONSTANTS[AC_zna] * ALGEBRAIC[AV_JnakNa] + CONSTANTS[AC_zk] * ALGEBRAIC[AV_JnakK]);
    
    /* IKb */
    ALGEBRAIC[AV_xkb] = 1.0 / (1.0 + exp(-(STATES[cell_v] - 14.48) / 18.34));
    ALGEBRAIC[AV_IKb_IKb] = CONSTANTS[AC_GKb] * ALGEBRAIC[AV_xkb] * (STATES[cell_v] - ALGEBRAIC[AV_EK]);
    
    /* INab */
    ALGEBRAIC[AV_INab_INab] = CONSTANTS[AC_PNab] * CONSTANTS[AC_F] * (STATES[nai] * exp(ALGEBRAIC[AV_vfrt]) - CONSTANTS[AC_nao]) * ((STATES[cell_v] == 0.0) ? 1.0 : ALGEBRAIC[AV_vfrt] / (exp(ALGEBRAIC[AV_vfrt]) - 1.0));
    
    /* ICab */
    ALGEBRAIC[AV_ICab_ICab] = CONSTANTS[AC_PCab] * 2.0 * CONSTANTS[AC_F] * (STATES[cai] * exp(2.0 * ALGEBRAIC[AV_vfrt]) - 0.341 * CONSTANTS[AC_cao]) * ((STATES[cell_v] == 0.0) ? 1.0 : 2.0 * ALGEBRAIC[AV_vfrt] / (exp(2.0 * ALGEBRAIC[AV_vfrt]) - 1.0));
    
    /* IpCa */
    ALGEBRAIC[AV_IpCa_IpCa] = CONSTANTS[AC_GpCa] * STATES[cai] / (0.0005 + STATES[cai]);
    
    /* SR_uptake */
    ALGEBRAIC[AV_Jleak] = 0.005 * STATES[cansr] / 15.0;
    ALGEBRAIC[AV_Jtr] = (STATES[cansr] - STATES[cajsr]) / 100.0;
    ALGEBRAIC[AV_Jtr2] = (STATES[cansr] - STATES[cacsr]) / 100.0;
    ALGEBRAIC[AV_Jupnp] = 1.0 * 0.005 * STATES[cai] / (STATES[cai] + 0.001);
    ALGEBRAIC[AV_Jupnp2] = 1.0 * 0.005 * STATES[cai2] / (STATES[cai2] + 0.001);
    ALGEBRAIC[AV_Jupp] = 1.0 * 3.0 * 0.005 * STATES[cai] / (STATES[cai] + 0.001 - 0.0002);
    ALGEBRAIC[AV_Jupp2] = 1.0 * 3.0 * 0.005 * STATES[cai2] / (STATES[cai2] + 0.001 - 0.0002);
    ALGEBRAIC[AV_fJupp] = 1.0 / (1.0 + CONSTANTS[AC_KmCaMK] / ALGEBRAIC[AV_CaMKa]);
    ALGEBRAIC[AV_Jup] = (1.0 - ALGEBRAIC[AV_fJupp]) * ALGEBRAIC[AV_Jupnp] + ALGEBRAIC[AV_fJupp] * ALGEBRAIC[AV_Jupp] - ALGEBRAIC[AV_Jleak];
    ALGEBRAIC[AV_Jup2] = (1.0 - ALGEBRAIC[AV_fJupp]) * ALGEBRAIC[AV_Jupnp2] + ALGEBRAIC[AV_fJupp] * ALGEBRAIC[AV_Jupp2] - ALGEBRAIC[AV_Jleak];
    
    /* ionic_concentrations */
    ALGEBRAIC[AV_Bcacsr] = 1.0 / (1.0 + CONSTANTS[AC_csqnmax] * CONSTANTS[AC_kmcsqn] / ((CONSTANTS[AC_kmcsqn] + STATES[cacsr]) * (CONSTANTS[AC_kmcsqn] + STATES[cacsr])));
    ALGEBRAIC[AV_Bcai] = 1.0 / (1.0 + CONSTANTS[AC_cmdnmax] * CONSTANTS[AC_kmcmdn] / ((CONSTANTS[AC_kmcmdn] + STATES[cai]) * (CONSTANTS[AC_kmcmdn] + STATES[cai])) + CONSTANTS[AC_trpnmax] * CONSTANTS[AC_kmtrpn] / ((CONSTANTS[AC_kmtrpn] + STATES[cai]) * (CONSTANTS[AC_kmtrpn] + STATES[cai])));
    ALGEBRAIC[AV_Bcai2] = 1.0 / (1.0 + CONSTANTS[AC_cmdnmax] * CONSTANTS[AC_kmcmdn] / ((CONSTANTS[AC_kmcmdn] + STATES[cai2]) * (CONSTANTS[AC_kmcmdn] + STATES[cai2])) + CONSTANTS[AC_trpnmax] * CONSTANTS[AC_kmtrpn] / ((CONSTANTS[AC_kmtrpn] + STATES[cai2]) * (CONSTANTS[AC_kmtrpn] + STATES[cai2])));
    ALGEBRAIC[AV_Bcajsr] = 1.0 / (1.0 + CONSTANTS[AC_csqnmax] * CONSTANTS[AC_kmcsqn] / ((CONSTANTS[AC_kmcsqn] + STATES[cajsr]) * (CONSTANTS[AC_kmcsqn] + STATES[cajsr])));
    ALGEBRAIC[AV_Bcass] = 1.0 / (1.0 + CONSTANTS[AC_BSRmax] * CONSTANTS[AC_KmBSR] / ((CONSTANTS[AC_KmBSR] + STATES[cass]) * (CONSTANTS[AC_KmBSR] + STATES[cass])) + CONSTANTS[AC_BSLmax] * CONSTANTS[AC_KmBSL] / ((CONSTANTS[AC_KmBSL] + STATES[cass]) * (CONSTANTS[AC_KmBSL] + STATES[cass])));
    ALGEBRAIC[AV_cai_avg] = CONSTANTS[AC_vmyo1frac] * STATES[cai] + (1.0 - CONSTANTS[AC_vmyo1frac]) * STATES[cai2];
    ALGEBRAIC[AV_diff_cansr] = ALGEBRAIC[AV_Jup] * (CONSTANTS[AC_vnsr1] / CONSTANTS[AC_vnsr]) + ALGEBRAIC[AV_Jup2] * (CONSTANTS[AC_vnsr2] / CONSTANTS[AC_vnsr]) - ALGEBRAIC[AV_Jtr] * CONSTANTS[AC_vjsr] / CONSTANTS[AC_vnsr] - ALGEBRAIC[AV_Jtr2] * CONSTANTS[AC_vcsr] / CONSTANTS[AC_vnsr];
    ALGEBRAIC[AV_diff_ki] = -(ALGEBRAIC[AV_IKr_IKr] + ALGEBRAIC[AV_IKs_IKs] + ALGEBRAIC[AV_IK1_IK1] + ALGEBRAIC[AV_IKb_IKb] - 2.0 * ALGEBRAIC[AV_INaK_INaK]) * CONSTANTS[AC_Acap] / (2.0 * CONSTANTS[AC_F] * CONSTANTS[AC_vmyo]) + ALGEBRAIC[AV_JdiffK] * CONSTANTS[AC_vss] / CONSTANTS[AC_vmyo];
    ALGEBRAIC[AV_diff_kss] = -ALGEBRAIC[AV_ICaK] * CONSTANTS[AC_Acap] / (2.0 * CONSTANTS[AC_F] * CONSTANTS[AC_vss]) - ALGEBRAIC[AV_JdiffK];
    ALGEBRAIC[AV_diff_nai] = -(ALGEBRAIC[AV_INa] + ALGEBRAIC[AV_INaL_INaL] + 3.0 * ALGEBRAIC[AV_INaCa_i_INaCa_i] + 3.0 * ALGEBRAIC[AV_INaK_INaK] + ALGEBRAIC[AV_INab_INab]) * CONSTANTS[AC_Acap] / (2.0 * CONSTANTS[AC_F] * CONSTANTS[AC_vmyo]) + ALGEBRAIC[AV_JdiffNa] * CONSTANTS[AC_vss] / CONSTANTS[AC_vmyo];
    ALGEBRAIC[AV_diff_nass] = -(ALGEBRAIC[AV_ICaNa] + 3.0 * ALGEBRAIC[AV_INaCa_ss_INaCa_ss]) * CONSTANTS[AC_Acap] / (2.0 * CONSTANTS[AC_F] * CONSTANTS[AC_vss]) - ALGEBRAIC[AV_JdiffNa];
    RATES[cansr] = ALGEBRAIC[AV_diff_cansr];
    RATES[ki] = ALGEBRAIC[AV_diff_ki];
    RATES[kss] = ALGEBRAIC[AV_diff_kss];
    RATES[nai] = ALGEBRAIC[AV_diff_nai];
    RATES[nass] = ALGEBRAIC[AV_diff_nass];
    ALGEBRAIC[AV_diff_cacsr] = ALGEBRAIC[AV_Bcacsr] * (ALGEBRAIC[AV_Jtr2] - STATES[CICR_Jrel2]);
    ALGEBRAIC[AV_diff_cai] = ALGEBRAIC[AV_Bcai] * (-(ALGEBRAIC[AV_IpCa_IpCa] + ALGEBRAIC[AV_ICab_ICab] - 2.0 * ALGEBRAIC[AV_INaCa_i_INaCa_i]) * CONSTANTS[AC_Acap] / (2.0 * 2.0 * CONSTANTS[AC_F] * CONSTANTS[AC_vmyo1]) - ALGEBRAIC[AV_Jup] * CONSTANTS[AC_vnsr1] / CONSTANTS[AC_vmyo1] + ALGEBRAIC[AV_Jdiff] * CONSTANTS[AC_vss] / CONSTANTS[AC_vmyo1] - ALGEBRAIC[AV_Jgap]);
    ALGEBRAIC[AV_diff_cai2] = ALGEBRAIC[AV_Bcai2] * (STATES[CICR_Jrel2] * (CONSTANTS[AC_vcsr] / CONSTANTS[AC_vmyo2]) + ALGEBRAIC[AV_Jgap] * CONSTANTS[AC_vmyo2] / CONSTANTS[AC_vmyo1] - ALGEBRAIC[AV_Jup2] * CONSTANTS[AC_vnsr2] / CONSTANTS[AC_vmyo2]);
    ALGEBRAIC[AV_diff_cajsr] = ALGEBRAIC[AV_Bcajsr] * (ALGEBRAIC[AV_Jtr] - ALGEBRAIC[AV_Jrel]);
    ALGEBRAIC[AV_diff_cass] = ALGEBRAIC[AV_Bcass] * (-(ALGEBRAIC[AV_ICaL_ICaL] - 2.0 * ALGEBRAIC[AV_INaCa_ss_INaCa_ss]) * CONSTANTS[AC_Acap] / (2.0 * 2.0 * CONSTANTS[AC_F] * CONSTANTS[AC_vss]) + ALGEBRAIC[AV_Jrel] * CONSTANTS[AC_vjsr] / CONSTANTS[AC_vss] - ALGEBRAIC[AV_Jdiff]);
    RATES[cacsr] = ALGEBRAIC[AV_diff_cacsr];
    RATES[cai] = ALGEBRAIC[AV_diff_cai];
    RATES[cai2] = ALGEBRAIC[AV_diff_cai2];
    RATES[cajsr] = ALGEBRAIC[AV_diff_cajsr];
    RATES[cass] = ALGEBRAIC[AV_diff_cass];
    
    /* *remaining* */
    ALGEBRAIC[AV_kito2] = 1.0 - 1.0 / (1.0 + ALGEBRAIC[AV_Jrel] / 0.4 * ALGEBRAIC[AV_Jrel] / 0.4);
    ALGEBRAIC[AV_Rel1] = (-ALGEBRAIC[AV_ICaL_ICaL] + 2.0 * ALGEBRAIC[AV_INaCa_ss_INaCa_ss]) * (CONSTANTS[AC_Acap] / (2.0 * CONSTANTS[AC_vss] * CONSTANTS[AC_F])) + ALGEBRAIC[AV_Jrel] * (CONSTANTS[AC_vjsr] / CONSTANTS[AC_vss]) - ALGEBRAIC[AV_Jdiff];
    ALGEBRAIC[AV_Rel2] = ALGEBRAIC[AV_Jgap] + STATES[CICR_Jrel2] * (CONSTANTS[AC_vcsr] / CONSTANTS[AC_vmyo2]) - ALGEBRAIC[AV_Jup2] * (CONSTANTS[AC_vnsr2] / CONSTANTS[AC_vmyo2]);
    ALGEBRAIC[AV_ITo_ITo] = CONSTANTS[AC_Gto] * STATES[ITo_aa] * ALGEBRAIC[AV_rito2] * ALGEBRAIC[AV_kito2] * (STATES[cell_v] - CONSTANTS[AC_ECl]);
    ALGEBRAIC[AV_Jrel1_inf] = ((ALGEBRAIC[AV_Rel1] > 0.0) ? 1.0 * 60.0 * ALGEBRAIC[AV_Rel1] * (1.0 + 1.0 / (1.0 + pow(CONSTANTS[AC_KmCaMK] / ALGEBRAIC[AV_CaMKa], 8.0))) / (1.0 + pow(0.75 / STATES[cajsr], 8.0)) : 0.0);
    ALGEBRAIC[AV_Jrel2_inf] = ((ALGEBRAIC[AV_Rel2] > 0.0) ? 1.0 * (250.0 * ALGEBRAIC[AV_Rel2] * (1.0 + 1.0 / (1.0 + pow(CONSTANTS[AC_KmCaMK] / ALGEBRAIC[AV_CaMKa], 8.0)))) / (1.0 + pow(0.75 / STATES[cacsr], 8.0)) : 0.0);
    RATES[CICR_Jrel1] = (ALGEBRAIC[AV_Jrel1_inf] - STATES[CICR_Jrel1]) / ALGEBRAIC[AV_tau_Jrel1];
    RATES[CICR_Jrel2] = (ALGEBRAIC[AV_Jrel2_inf] - STATES[CICR_Jrel2]) / ALGEBRAIC[AV_tau_Jrel2];
    
    ALGEBRAIC[Iion_cm] = ALGEBRAIC[AV_INa] + ALGEBRAIC[AV_INaL_INaL] + ALGEBRAIC[AV_ICaL_ICaL] + ALGEBRAIC[AV_ICaNa] + ALGEBRAIC[AV_ICaK] + ALGEBRAIC[AV_IKr_IKr] + ALGEBRAIC[AV_IKs_IKs] + ALGEBRAIC[AV_IK1_IK1] + ALGEBRAIC[AV_ITo_ITo] + ALGEBRAIC[AV_INaCa_i_INaCa_i] + ALGEBRAIC[AV_INaCa_ss_INaCa_ss] + ALGEBRAIC[AV_INaK_INaK] + ALGEBRAIC[AV_INab_INab] + ALGEBRAIC[AV_IKb_IKb] + ALGEBRAIC[AV_IpCa_IpCa] + ALGEBRAIC[AV_ICab_ICab];
    
    ALGEBRAIC[Istim] = computeIstim(VOI, CONSTANTS);

    if (solveVmWithinODESolver)
    {
        RATES[0] = -ALGEBRAIC[Iion_cm] - ALGEBRAIC[Istim];
    }

}
