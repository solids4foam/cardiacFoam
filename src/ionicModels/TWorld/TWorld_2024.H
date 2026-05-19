/*---------------------------------------------------------------------------*\
License
    This file is part of cardiacFoam.

    cardiacFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    cardiacFoam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with cardiacFoam.  If not, see <http://www.gnu.org/licenses/>.

Description
    Self-contained ODE equations for the TWorld 2024 cardiac
    cell model. Converted from CellML.

Author
    Simao Nieto de Castro, UCD.
\*---------------------------------------------------------------------------*/

#include "TWorld_2024Names.H"
#include "stimulusIO.H"

static const char* TWorldSTATES_NAMES[NUM_STATES] = {
    "v",
    "camk_trap",
    "camk_f_ICaL",
    "camk_f_RyR",
    "camk_f_PLB",
    "casig_serca_trap",
    "buffers_NaBj",
    "buffers_NaBsl",
    "buffers_TnClow",
    "buffers_TnCHc",
    "buffers_TnCHm",
    "buffers_CaM",
    "buffers_Myosin_ca",
    "buffers_Myosin_mg",
    "buffers_SRB",
    "buffers_SLLj",
    "buffers_SLLsl",
    "buffers_SLHj",
    "buffers_SLHsl",
    "buffers_Csqn",
    "naj",
    "nasl",
    "nai",
    "ki",
    "cli",
    "casr",
    "caj",
    "casl",
    "cai",
    "m",
    "h",
    "j",
    "hp",
    "jp",
    "m_P",
    "h_P",
    "j_P",
    "hp_P",
    "jp_P",
    "mL",
    "hL",
    "hLp",
    "d",
    "ff",
    "fs",
    "fcaf",
    "fcas",
    "jca",
    "ffp",
    "fcafp",
    "nca",
    "nca_i",
    "d_P",
    "ff_P",
    "fs_P",
    "fcaf_P",
    "fcas_P",
    "fBPf",
    "fcaBPf",
    "ical_pureCDI_junc",
    "ical_pureCDI_sl",
    "xtos",
    "ytos",
    "xtof",
    "ytof",
    "xtos_p",
    "xtof_p",
    "ytos_p",
    "ytof_p",
    "C0",
    "C1",
    "C2",
    "I",
    "O",
    "xs_junc",
    "xs_sl",
    "jrel_icaldep_act",
    "jrel_icaldep_f1",
    "jrel_icaldep_f2",
    "ryr_R",
    "ryr_O",
    "ryr_I",
    "ryr_CaRI",
    "ryr_R_p",
    "ryr_O_p",
    "ryr_I_p",
    "ryr_CaRI_p",
    "contraction_TmBlocked",
    "contraction_XW",
    "contraction_XS",
    "contraction_ZETAS",
    "contraction_ZETAW",
    "contraction_Ca_TRPN",
};

static const char* TWorldALGEBRAIC_NAMES[NUM_ALGEBRAIC] = {
    "AV_CaMK_Phos_ss_ICaL",
    "AV_CaMK_Phos_ss_PLB",
    "AV_CaMK_Phos_ss_RyR",
    "AV_CaMK_active",
    "AV_bound",
    "AV_bound_serca",
    "AV_casig_SERCA_act",
    "AV_Afcaf",
    "AV_Afcas",
    "AV_ICaK",
    "AV_ICaK_i_BP",
    "AV_ICaK_i_CaMK",
    "AV_ICaK_i_NP",
    "AV_ICaK_i_PKA",
    "AV_ICaK_junc",
    "AV_ICaK_sl",
    "AV_ICaK_ss_BP",
    "AV_ICaK_ss_CaMK",
    "AV_ICaK_ss_NP",
    "AV_ICaK_ss_PKA",
    "AV_ICaL_ICaL",
    "AV_ICaL_i_BP",
    "AV_ICaL_i_CaMK",
    "AV_ICaL_i_NP",
    "AV_ICaL_i_PKA",
    "AV_ICaL_junc",
    "AV_ICaL_sl",
    "AV_ICaL_ss_BP",
    "AV_ICaL_ss_CaMK",
    "AV_ICaL_ss_NP",
    "AV_ICaL_ss_PKA",
    "AV_ICaNa",
    "AV_ICaNa_i_BP",
    "AV_ICaNa_i_CaMK",
    "AV_ICaNa_i_NP",
    "AV_ICaNa_i_PKA",
    "AV_ICaNa_junc",
    "AV_ICaNa_sl",
    "AV_ICaNa_ss_BP",
    "AV_ICaNa_ss_CaMK",
    "AV_ICaNa_ss_NP",
    "AV_ICaNa_ss_PKA",
    "AV_ICa_tot",
    "AV_Ii",
    "AV_Iss",
    "AV_PhiCaK_i",
    "AV_PhiCaK_ss",
    "AV_PhiCaL_i",
    "AV_PhiCaL_ss",
    "AV_PhiCaNa_i",
    "AV_PhiCaNa_ss",
    "AV_anca",
    "AV_anca_i",
    "AV_dPss",
    "AV_dss",
    "AV_f",
    "AV_fBP",
    "AV_fBPss",
    "AV_fICaL_BP",
    "AV_fICaL_CaMKonly",
    "AV_fICaL_PKAonly",
    "AV_fICaLp",
    "AV_f_P",
    "AV_fca",
    "AV_fcaBP",
    "AV_fcaBPss",
    "AV_fcap",
    "AV_fcap_P",
    "AV_fcass",
    "AV_fcass_P",
    "AV_fp",
    "AV_fss",
    "AV_fss_P",
    "AV_gamma_cai",
    "AV_gamma_cass",
    "AV_gamma_ki",
    "AV_gamma_kss",
    "AV_gamma_nai",
    "AV_gamma_nass",
    "AV_jcass",
    "AV_km2n",
    "AV_sigmoidTransition",
    "AV_sigmoidTransition2",
    "AV_tauTransition",
    "AV_tauTransition2",
    "AV_td",
    "AV_tfcaf",
    "AV_tfcafp",
    "AV_tfcas",
    "AV_tff",
    "AV_tffp",
    "AV_tfs",
    "AV_ICab_ICab",
    "AV_ICab_junc",
    "AV_ICab_sl",
    "AV_IClCa",
    "AV_IClCa_junc",
    "AV_IClCa_sl",
    "AV_IClb",
    "AV_IK1_IK1",
    "AV_K1ss",
    "AV_aK1",
    "AV_bK1",
    "AV_IKb_IKb",
    "AV_xkb",
    "AV_IKr_IKr",
    "AV_IKr_alpha",
    "AV_alpha_2",
    "AV_alpha_C2ToI",
    "AV_alpha_i",
    "AV_IKr_beta",
    "AV_beta_2",
    "AV_beta_ItoC2",
    "AV_beta_i",
    "AV_IKs_IKs",
    "AV_IKs_junc",
    "AV_IKs_sl",
    "AV_VsTs_Ca_junc",
    "AV_VsTs_Ca_sl",
    "AV_VsXs_Ca_junc",
    "AV_VsXs_Ca_sl",
    "AV_gks_junc",
    "AV_gks_sl",
    "AV_tauxs_junc",
    "AV_tauxs_sl",
    "AV_xsss_junc",
    "AV_xsss_sl",
    "AV_INa_INa",
    "AV_INaBase",
    "AV_INaBase_BP",
    "AV_INaBase_CaMK",
    "AV_INaBase_NP",
    "AV_INaBase_PKA",
    "AV_INaj",
    "AV_INasl",
    "AV_ah",
    "AV_aj",
    "AV_bh",
    "AV_bj",
    "AV_fINa_BP",
    "AV_fINa_CaMKonly",
    "AV_fINa_PKAonly",
    "AV_fINap",
    "AV_hss",
    "AV_hss_P",
    "AV_hssp",
    "AV_hssp_P",
    "AV_jss",
    "AV_jss_P",
    "AV_jssp_P",
    "AV_mss",
    "AV_mss_P",
    "AV_th",
    "AV_tj",
    "AV_tjp",
    "AV_tm",
    "AV_E1_i",
    "AV_E1_ss",
    "AV_E2_i",
    "AV_E2_ss",
    "AV_E3_i",
    "AV_E3_ss",
    "AV_E4_i",
    "AV_E4_ss",
    "AV_INaCa_INaCa",
    "AV_INaCa_i",
    "AV_INaCa_ss",
    "AV_JncxCa_i",
    "AV_JncxCa_ss",
    "AV_JncxNa_i",
    "AV_JncxNa_ss",
    "AV_allo_i",
    "AV_allo_ss",
    "AV_h1_i",
    "AV_h1_ss",
    "AV_h2_i",
    "AV_h2_ss",
    "AV_h3_i",
    "AV_h3_ss",
    "AV_h4_i",
    "AV_h4_ss",
    "AV_h5_i",
    "AV_h5_ss",
    "AV_h6_i",
    "AV_h6_ss",
    "AV_h7_i",
    "AV_h7_ss",
    "AV_h8_i",
    "AV_h8_ss",
    "AV_h9_i",
    "AV_h9_ss",
    "AV_hca",
    "AV_hna",
    "AV_k3_i",
    "AV_k3_ss",
    "AV_k3p_i",
    "AV_k3p_ss",
    "AV_k3pp_i",
    "AV_k3pp_ss",
    "AV_k4_i",
    "AV_k4_ss",
    "AV_k4p_i",
    "AV_k4p_ss",
    "AV_k4pp_i",
    "AV_k4pp_ss",
    "AV_k6_i",
    "AV_k6_ss",
    "AV_k7_i",
    "AV_k7_ss",
    "AV_k8_i",
    "AV_k8_ss",
    "AV_x1_i",
    "AV_x1_ss",
    "AV_x2_i",
    "AV_x2_ss",
    "AV_x3_i",
    "AV_x3_ss",
    "AV_x4_i",
    "AV_x4_ss",
    "AV_INaK_INaK",
    "AV_INaK_junc_PKA",
    "AV_INaK_junc_noPKA",
    "AV_INaK_sl_PKA",
    "AV_INaK_sl_noPKA",
    "AV_INaKj",
    "AV_INaKsl",
    "AV_fnak",
    "AV_GNaL",
    "AV_INaL_INaL",
    "AV_INaLj",
    "AV_INaLsl",
    "AV_fINaLp",
    "AV_hLss",
    "AV_hLssp",
    "AV_mLss",
    "AV_tmL",
    "AV_INab_INab",
    "AV_INabj",
    "AV_INabsl",
    "AV_IpCa_IpCa",
    "AV_IpCa_junc",
    "AV_IpCa_sl",
    "AV_Ito_Ito",
    "AV_Itof",
    "AV_Itos",
    "AV_dti_develop",
    "AV_dti_recover",
    "AV_fItop",
    "AV_tauxtof",
    "AV_tauxtos",
    "AV_tauytof",
    "AV_tauytof_p",
    "AV_tauytos",
    "AV_tauytos_p",
    "AV_xtoss",
    "AV_xtoss_p",
    "AV_ytoss",
    "AV_Jserca",
    "AV_Jserca_np",
    "AV_Jserca_p",
    "AV_Vmax_mult",
    "AV_phosphorylationTotal",
    "AV_Ta",
    "AV_XU",
    "AV_d_contraction_Ca_TRPN",
    "AV_gamma_rate",
    "AV_gamma_rate_w",
    "AV_xb_su",
    "AV_xb_su_gamma",
    "AV_xb_uw",
    "AV_xb_ws",
    "AV_xb_wu",
    "AV_xb_wu_gamma",
    "AV_time",
    "AV_ICa_tot_junc",
    "AV_ICa_tot_sl",
    "AV_ICl_tot",
    "AV_IK_tot",
    "AV_INa_tot_junc",
    "AV_INa_tot_sl",
    "AV_J_CaB_cytosol",
    "AV_J_CaB_junction",
    "AV_J_CaB_sl",
    "AV_d_buffers_CaM",
    "AV_d_buffers_Csqn",
    "AV_d_buffers_Myosin_ca",
    "AV_d_buffers_NaBj",
    "AV_d_buffers_NaBsl",
    "AV_d_buffers_SLHj",
    "AV_d_buffers_SLHsl",
    "AV_d_buffers_SLLj",
    "AV_d_buffers_SLLsl",
    "AV_d_buffers_SRB",
    "AV_d_buffers_TnCHc",
    "AV_d_buffers_TnClow",
    "AV_vffrt",
    "AV_vfrt",
    "AV_ECaj",
    "AV_ECasl",
    "AV_ECl",
    "AV_EK",
    "AV_EKs",
    "AV_ENaj",
    "AV_ENasl",
    "AV_CaMKIILeakMultiplier",
    "AV_ICaL_junc_positive",
    "AV_ICaL_junc_sigmoided",
    "AV_J_SRCarel",
    "AV_J_SRCarel_np",
    "AV_J_SRCarel_p",
    "AV_J_SRleak",
    "AV_Jrel_ICaLdep",
    "AV_Jrel_inact_inf",
    "AV_Jrel_inact_inf2",
    "AV_Jrel_inf",
    "AV_RI_to_CI",
    "AV_RIcleft",
    "AV_RIcleft_p",
    "AV_kCaSR",
    "AV_kiSRCa",
    "AV_koSRCa",
    "AV_nonlinearModifier",
    "AV_sigmoidBaseCaI",
    "AV_tau_rel",
    "Istim",
    "Iion_cm",
};


static const char* TWorldCONSTANTS_NAMES[NUM_CONSTANTS] = {
    "AC_CaMK0",
    "AC_K_Phos_CaMK",
    "AC_Km_CaMK_Ca",
    "AC_PP1_tot",
    "AC_Whole_cell_PP1",
    "AC_CaMK_alpha",
    "AC_alpha_serca",
    "AC_CaMK_beta",
    "AC_tau_cal",
    "AC_tau_plb",
    "AC_tau_ryr",
    "AC_Aff",
    "AC_Afs",
    "AC_ICaL_fractionSS",
    "AC_Io",
    "AC_Kmn",
    "AC_PCa",
    "AC_PCaK",
    "AC_PCaK_P",
    "AC_PCaKp",
    "AC_PCaNa",
    "AC_PCaNa_P",
    "AC_PCaNap",
    "AC_PCa_P",
    "AC_PCa_Pb",
    "AC_PCa_b",
    "AC_PCap",
    "AC_constA",
    "AC_dielConstant",
    "AC_fICaLP",
    "AC_fICaL_P",
    "AC_gamma_cao",
    "AC_gamma_ko",
    "AC_gamma_nao",
    "AC_k2n",
    "AC_rateRecovery",
    "AC_tjca",
    "AC_GCab",
    "AC_GCab_b",
    "AC_GClCa",
    "AC_GClCa_b",
    "AC_GClb",
    "AC_GClb_b",
    "AC_KdClCa",
    "AC_GK1",
    "AC_GK1_b",
    "AC_IK1_celltype_factor",
    "AC_IK1_sex_factor",
    "AC_GKb",
    "AC_GKb_b",
    "AC_GKr",
    "AC_GKr_b",
    "AC_alpha_1",
    "AC_IKr_beta_1",
    "AC_IKr_celltype_factor",
    "AC_IKr_sex_factor",
    "AC_GKs",
    "AC_GKs_b",
    "AC_P_g_0",
    "AC_P_g_max",
    "AC_P_tau_0",
    "AC_P_tau_max",
    "AC_P_vh_0",
    "AC_P_vh_max",
    "AC_IKs_celltype_factor",
    "AC_gKs_factor",
    "AC_kPKA_IKs",
    "AC_IKs_sex_factor",
    "AC_GNa",
    "AC_GNa_P",
    "AC_GNa_b",
    "AC_fINa_P",
    "AC_Gncx",
    "AC_Gncx_b",
    "AC_INaCa_fractionSS",
    "AC_KmCaAct",
    "AC_INaCa_celltype_factor",
    "AC_h10_i",
    "AC_h10_ss",
    "AC_h11_i",
    "AC_h11_ss",
    "AC_h12_i",
    "AC_h12_ss",
    "AC_k1_i",
    "AC_k1_ss",
    "AC_k2_i",
    "AC_k2_ss",
    "AC_k5_i",
    "AC_k5_ss",
    "AC_kasymm",
    "AC_kcaoff",
    "AC_kcaon",
    "AC_kna1",
    "AC_kna2",
    "AC_kna3",
    "AC_qca",
    "AC_qna",
    "AC_INaCa_sex_factor",
    "AC_wca",
    "AC_wna",
    "AC_wnaca",
    "AC_IbarNaK",
    "AC_IbarNaK_b",
    "AC_KmKo",
    "AC_KmNaip",
    "AC_KmNaip_PKA",
    "AC_GNaL_b",
    "AC_thL",
    "AC_thLp",
    "AC_GNab",
    "AC_GNab_b",
    "AC_IbarSLCaP",
    "AC_IbarSLCaP_b",
    "AC_KmPCa",
    "AC_Q10SLCaP",
    "AC_Gto_fast",
    "AC_Gto_slow",
    "AC_fICaL_PKA",
    "AC_fIKs_PKA",
    "AC_fINaK_PKA",
    "AC_fINa_PKA",
    "AC_fMyBPC_PKA",
    "AC_fPLB_PKA",
    "AC_fTnI_PKA",
    "AC_Km_SERCA_Ca",
    "AC_Kmf",
    "AC_Kmf_p",
    "AC_Kmr",
    "AC_Max_Vmax_SERCA_Ca",
    "AC_Q10SRCaP",
    "AC_Vmax_SRCaP",
    "AC_Vmax_SRCaP_b",
    "AC_hillSRCaP",
    "AC_cellLength",
    "AC_cellRadius",
    "AC_vcell",
    "AC_vjunc",
    "AC_vmyo",
    "AC_vsl",
    "AC_vsr",
    "AC_A",
    "AC_Lfac",
    "AC_PKAForceMultiplier",
    "AC_TOT_A",
    "AC_TRPN_n",
    "AC_Tref",
    "AC_XSSS",
    "AC_XWSS",
    "AC_beta_0",
    "AC_contraction_beta_1",
    "AC_ca50",
    "AC_cds",
    "AC_cdw",
    "AC_dr",
    "AC_fPKA_TnI",
    "AC_fracTnIpo",
    "AC_contraction_gamma",
    "AC_gamma_wu",
    "AC_k_su",
    "AC_k_uw",
    "AC_k_ws",
    "AC_k_wu",
    "AC_koff",
    "AC_ktm_block",
    "AC_ktm_unblock",
    "AC_lambda",
    "AC_lambda0",
    "AC_lambda_max",
    "AC_lambda_min",
    "AC_lambda_rate",
    "AC_mu",
    "AC_nperm",
    "AC_nu",
    "AC_perm50",
    "AC_phi",
    "AC_wfrac",
    "AC_sex",
    "AC_cao",
    "AC_clo",
    "AC_ko",
    "AC_nao",
    "AC_Bmax_CaM",
    "AC_Bmax_Csqn",
    "AC_Bmax_Naj",
    "AC_Bmax_Nasl",
    "AC_Bmax_SLhighj",
    "AC_Bmax_SLhighsl",
    "AC_Bmax_SLlowj",
    "AC_Bmax_SLlowsl",
    "AC_Bmax_SR",
    "AC_Bmax_TnChigh",
    "AC_Bmax_TnClow",
    "AC_Bmax_myosin",
    "AC_J_ca_juncsl",
    "AC_J_ca_slmyo",
    "AC_J_na_juncsl",
    "AC_J_na_slmyo",
    "AC_koff_cam",
    "AC_koff_csqn",
    "AC_koff_myoca",
    "AC_koff_myomg",
    "AC_koff_na",
    "AC_koff_slh",
    "AC_koff_sll",
    "AC_koff_sr",
    "AC_koff_tnchca",
    "AC_koff_tnchmg",
    "AC_kon_cam",
    "AC_kon_csqn",
    "AC_kon_myoca",
    "AC_kon_myomg",
    "AC_kon_na",
    "AC_kon_slh",
    "AC_kon_sll",
    "AC_kon_sr",
    "AC_kon_tnchca",
    "AC_kon_tnchmg",
    "AC_mgi",
    "AC_Cmem",
    "AC_Fjunc",
    "AC_Fsl",
    "AC_ICaLPCa_multiplier",
    "AC_ICab_multiplier",
    "AC_IClCa_multiplier",
    "AC_IClb_multiplier",
    "AC_IK1_multiplier",
    "AC_IKb_multiplier",
    "AC_IKr_multiplier",
    "AC_IKs_multiplier",
    "AC_INaCa_multiplier",
    "AC_INaK_multiplier",
    "AC_INaL_multiplier",
    "AC_INa_multiplier",
    "AC_INab_multiplier",
    "AC_IpCa_multiplier",
    "AC_Itof_multiplier",
    "AC_Itos_multiplier",
    "AC_Jrel_multiplier",
    "AC_Jup_multiplier",
    "AC_F",
    "AC_Qpow",
    "AC_R",
    "AC_T",
    "AC_zca",
    "AC_zcl",
    "AC_zk",
    "AC_zna",
    "AC_PNaK",
    "AC_CI_to_RI",
    "AC_MaxSR",
    "AC_MinSR",
    "AC_a_rel",
    "AC_baseRateCaI",
    "AC_bt",
    "AC_caExpFactor",
    "AC_caExpFactor2",
    "AC_caTransFactor",
    "AC_caTransFactor2",
    "AC_caTransFactor2p",
    "AC_directRelMidpoint",
    "AC_ec50SR",
    "AC_ecCaI",
    "AC_kiCa",
    "AC_kim",
    "AC_koCa",
    "AC_kom",
    "AC_ks",
    "AC_maxCaI",
    "AC_minCaI",
    "AC_steepnessCaI",
    "AC_steepnessCaSR",
    "AC_tauInact",
    "AC_tauInact2",
    "gnalTissueScale",
};

static inline const Foam::HashTable<Foam::wordList>&
TWorldDependencyMap()
{
   static const Foam::HashTable<Foam::wordList> dep
    (
        {
        {"INa", Foam::wordList({"AV_INa_INa","AV_INaj","AV_INasl","m","h","j","hp","jp","m_P","h_P","j_P","AV_mss","AV_hss","AV_jss","AV_tm","AV_th","AV_tj"})},
        {"I_Na", Foam::wordList({"AV_INa_INa","AV_INaj","AV_INasl","m","h","j","hp","jp","m_P","h_P","j_P","AV_mss","AV_hss","AV_jss","AV_tm","AV_th","AV_tj"})},
        {"AV_INa_INa", Foam::wordList({"AV_INa_INa","AV_INaj","AV_INasl","m","h","j","hp","jp","m_P","h_P","j_P","AV_mss","AV_hss","AV_jss","AV_tm","AV_th","AV_tj"})},

        {"INaL", Foam::wordList({"AV_INaL_INaL","AV_INaLj","AV_INaLsl","mL","hL","hLp","AV_GNaL","AV_mLss","AV_hLss","AV_hLssp","AV_tmL"})},
        {"I_NaL", Foam::wordList({"AV_INaL_INaL","AV_INaLj","AV_INaLsl","mL","hL","hLp","AV_GNaL","AV_mLss","AV_hLss","AV_hLssp","AV_tmL"})},
        {"AV_INaL_INaL", Foam::wordList({"AV_INaL_INaL","AV_INaLj","AV_INaLsl","mL","hL","hLp","AV_GNaL","AV_mLss","AV_hLss","AV_hLssp","AV_tmL"})},

        {"Ito", Foam::wordList({"AV_Ito_Ito","AV_Itof","AV_Itos","xtof","ytof","xtos","ytos","xtof_p","ytof_p","xtos_p","ytos_p","AV_xtoss","AV_xtoss_p","AV_ytoss","AV_tauxtof","AV_tauytof","AV_tauxtos","AV_tauytos"})},
        {"I_to", Foam::wordList({"AV_Ito_Ito","AV_Itof","AV_Itos","xtof","ytof","xtos","ytos","xtof_p","ytof_p","xtos_p","ytos_p","AV_xtoss","AV_xtoss_p","AV_ytoss","AV_tauxtof","AV_tauytof","AV_tauxtos","AV_tauytos"})},
        {"AV_Ito_Ito", Foam::wordList({"AV_Ito_Ito","AV_Itof","AV_Itos","xtof","ytof","xtos","ytos","xtof_p","ytof_p","xtos_p","ytos_p","AV_xtoss","AV_xtoss_p","AV_ytoss","AV_tauxtof","AV_tauytof","AV_tauxtos","AV_tauytos"})},

        {"ICaL", Foam::wordList({"AV_ICaL_ICaL","AV_ICaL_junc","AV_ICaL_sl","d","d_P","ff","fs","fcaf","fcas","jca","nca","nca_i","AV_dss","AV_dPss","AV_fss","AV_fss_P","AV_td","AV_tff","AV_tfs","AV_tfcaf","AV_tfcas"})},
        {"I_CaL", Foam::wordList({"AV_ICaL_ICaL","AV_ICaL_junc","AV_ICaL_sl","d","d_P","ff","fs","fcaf","fcas","jca","nca","nca_i","AV_dss","AV_dPss","AV_fss","AV_fss_P","AV_td","AV_tff","AV_tfs","AV_tfcaf","AV_tfcas"})},
        {"AV_ICaL_ICaL", Foam::wordList({"AV_ICaL_ICaL","AV_ICaL_junc","AV_ICaL_sl","d","d_P","ff","fs","fcaf","fcas","jca","nca","nca_i","AV_dss","AV_dPss","AV_fss","AV_fss_P","AV_td","AV_tff","AV_tfs","AV_tfcaf","AV_tfcas"})},

        {"ICaNa", Foam::wordList({"AV_ICaNa","AV_ICaNa_junc","AV_ICaNa_sl","AV_PhiCaNa_i","AV_PhiCaNa_ss"})},
        {"I_CaNa", Foam::wordList({"AV_ICaNa","AV_ICaNa_junc","AV_ICaNa_sl","AV_PhiCaNa_i","AV_PhiCaNa_ss"})},
        {"AV_ICaNa", Foam::wordList({"AV_ICaNa","AV_ICaNa_junc","AV_ICaNa_sl","AV_PhiCaNa_i","AV_PhiCaNa_ss"})},

        {"ICaK", Foam::wordList({"AV_ICaK","AV_ICaK_junc","AV_ICaK_sl","AV_PhiCaK_i","AV_PhiCaK_ss"})},
        {"I_CaK", Foam::wordList({"AV_ICaK","AV_ICaK_junc","AV_ICaK_sl","AV_PhiCaK_i","AV_PhiCaK_ss"})},
        {"AV_ICaK", Foam::wordList({"AV_ICaK","AV_ICaK_junc","AV_ICaK_sl","AV_PhiCaK_i","AV_PhiCaK_ss"})},

        {"IKr", Foam::wordList({"AV_IKr_IKr","C0","C1","C2","O","I","AV_IKr_alpha","AV_IKr_beta"})},
        {"I_Kr", Foam::wordList({"AV_IKr_IKr","C0","C1","C2","O","I","AV_IKr_alpha","AV_IKr_beta"})},
        {"AV_IKr_IKr", Foam::wordList({"AV_IKr_IKr","C0","C1","C2","O","I","AV_IKr_alpha","AV_IKr_beta"})},

        {"IKs", Foam::wordList({"AV_IKs_IKs","AV_IKs_junc","AV_IKs_sl","xs_junc","xs_sl","AV_xsss_junc","AV_xsss_sl","AV_tauxs_junc","AV_tauxs_sl"})},
        {"I_Ks", Foam::wordList({"AV_IKs_IKs","AV_IKs_junc","AV_IKs_sl","xs_junc","xs_sl","AV_xsss_junc","AV_xsss_sl","AV_tauxs_junc","AV_tauxs_sl"})},
        {"AV_IKs_IKs", Foam::wordList({"AV_IKs_IKs","AV_IKs_junc","AV_IKs_sl","xs_junc","xs_sl","AV_xsss_junc","AV_xsss_sl","AV_tauxs_junc","AV_tauxs_sl"})},

        {"IK1", Foam::wordList({"AV_IK1_IK1","AV_K1ss","AV_aK1","AV_bK1"})},
        {"I_K1", Foam::wordList({"AV_IK1_IK1","AV_K1ss","AV_aK1","AV_bK1"})},
        {"AV_IK1_IK1", Foam::wordList({"AV_IK1_IK1","AV_K1ss","AV_aK1","AV_bK1"})},

        {"INaCa", Foam::wordList({"AV_INaCa_INaCa","AV_INaCa_i","AV_INaCa_ss","AV_JncxNa_i","AV_JncxCa_i","AV_JncxNa_ss","AV_JncxCa_ss"})},
        {"I_NaCa", Foam::wordList({"AV_INaCa_INaCa","AV_INaCa_i","AV_INaCa_ss","AV_JncxNa_i","AV_JncxCa_i","AV_JncxNa_ss","AV_JncxCa_ss"})},
        {"AV_INaCa_INaCa", Foam::wordList({"AV_INaCa_INaCa","AV_INaCa_i","AV_INaCa_ss","AV_JncxNa_i","AV_JncxCa_i","AV_JncxNa_ss","AV_JncxCa_ss"})},

        {"INaK", Foam::wordList({"AV_INaK_INaK","AV_INaKj","AV_INaKsl","AV_fnak"})},
        {"I_NaK", Foam::wordList({"AV_INaK_INaK","AV_INaKj","AV_INaKsl","AV_fnak"})},
        {"AV_INaK_INaK", Foam::wordList({"AV_INaK_INaK","AV_INaKj","AV_INaKsl","AV_fnak"})},

        {"IKb", Foam::wordList({"AV_IKb_IKb","AV_xkb"})},
        {"I_Kb", Foam::wordList({"AV_IKb_IKb","AV_xkb"})},
        {"AV_IKb_IKb", Foam::wordList({"AV_IKb_IKb","AV_xkb"})},

        {"INab", Foam::wordList({"AV_INab_INab","AV_INabj","AV_INabsl"})},
        {"I_Nab", Foam::wordList({"AV_INab_INab","AV_INabj","AV_INabsl"})},
        {"AV_INab_INab", Foam::wordList({"AV_INab_INab","AV_INabj","AV_INabsl"})},

        {"ICab", Foam::wordList({"AV_ICab_ICab","AV_ICab_junc","AV_ICab_sl"})},
        {"I_Cab", Foam::wordList({"AV_ICab_ICab","AV_ICab_junc","AV_ICab_sl"})},
        {"AV_ICab_ICab", Foam::wordList({"AV_ICab_ICab","AV_ICab_junc","AV_ICab_sl"})},

        {"IpCa", Foam::wordList({"AV_IpCa_IpCa","AV_IpCa_junc","AV_IpCa_sl"})},
        {"I_pCa", Foam::wordList({"AV_IpCa_IpCa","AV_IpCa_junc","AV_IpCa_sl"})},
        {"AV_IpCa_IpCa", Foam::wordList({"AV_IpCa_IpCa","AV_IpCa_junc","AV_IpCa_sl"})},

        {"IClCa", Foam::wordList({"AV_IClCa","AV_IClCa_junc","AV_IClCa_sl"})},
        {"I_ClCa", Foam::wordList({"AV_IClCa","AV_IClCa_junc","AV_IClCa_sl"})},
        {"AV_IClCa", Foam::wordList({"AV_IClCa","AV_IClCa_junc","AV_IClCa_sl"})},

        {"IClb", Foam::wordList({"AV_IClb"})},
        {"I_Clb", Foam::wordList({"AV_IClb"})},
        {"AV_IClb", Foam::wordList({"AV_IClb"})},

        {"Iion", Foam::wordList({"AV_INa_tot_junc","AV_INa_tot_sl","AV_ICa_tot_junc","AV_ICa_tot_sl","AV_IK_tot","AV_ICl_tot","Iion_cm"})},
        {"I_ion", Foam::wordList({"AV_INa_tot_junc","AV_INa_tot_sl","AV_ICa_tot_junc","AV_ICa_tot_sl","AV_IK_tot","AV_ICl_tot","Iion_cm"})},
        {"Iion_cm", Foam::wordList({"AV_INa_tot_junc","AV_INa_tot_sl","AV_ICa_tot_junc","AV_ICa_tot_sl","AV_IK_tot","AV_ICl_tot","Iion_cm"})}
        }
    );

    return dep;
}

inline Foam::scalar computeIstim(Foam::scalar t, const Foam::StimulusProtocol& stimulus)
{
    return Foam::stimulusIO::computeStimulus(t, stimulus);
}

void
TWorldinitConsts(double* CONSTANTS, double* RATES, double* STATES, int tissueFlag, const Foam::dictionary& stimulus)
{
/* CaMK */
    CONSTANTS[AC_CaMK0] = 0.1;
    CONSTANTS[AC_K_Phos_CaMK] = 0.35;
    CONSTANTS[AC_Km_CaMK_Ca] = 0.0075;
    CONSTANTS[AC_Whole_cell_PP1] = 0.13698;
    CONSTANTS[AC_CaMK_alpha] = 0.05;
    CONSTANTS[AC_alpha_serca] = 0.05;
    CONSTANTS[AC_CaMK_beta] = 0.00068;
    CONSTANTS[AC_tau_plb] = 100000.0;
    CONSTANTS[AC_tau_ryr] = 10000.0;
    CONSTANTS[AC_PP1_tot] = CONSTANTS[AC_Whole_cell_PP1];
    CONSTANTS[AC_tau_cal] = CONSTANTS[AC_tau_ryr];
    
    /* PKA */
    CONSTANTS[AC_fICaL_PKA] = 0.0;
    CONSTANTS[AC_fIKs_PKA] = 0.0;
    CONSTANTS[AC_fINaK_PKA] = 0.0;
    CONSTANTS[AC_fINa_PKA] = 0.0;
    CONSTANTS[AC_fMyBPC_PKA] = 0.0;
    CONSTANTS[AC_fPLB_PKA] = 0.0;
    CONSTANTS[AC_fTnI_PKA] = 0.0;
    
    /* cell_geometry */
    CONSTANTS[AC_cellLength] = 100.0;
    CONSTANTS[AC_cellRadius] = 10.25;
    CONSTANTS[AC_vcell] = 1e-15 * 3.14159265358979312e+00 * pow(CONSTANTS[AC_cellRadius], 2.0) * CONSTANTS[AC_cellLength];
    CONSTANTS[AC_vjunc] = 0.000539 * CONSTANTS[AC_vcell];
    CONSTANTS[AC_vmyo] = 0.65 * CONSTANTS[AC_vcell];
    CONSTANTS[AC_vsl] = 0.02 * CONSTANTS[AC_vcell];
    CONSTANTS[AC_vsr] = 0.035 * CONSTANTS[AC_vcell];
    
    /* environment */
    CONSTANTS[AC_sex] = 0.0;
    
    /* extracellular */
    CONSTANTS[AC_cao] = 1.8;
    CONSTANTS[AC_clo] = 150.0;
    CONSTANTS[AC_ko] = 5.0;
    CONSTANTS[AC_nao] = 140.0;
    
    /* multipliers */
    CONSTANTS[AC_ICaLPCa_multiplier] = 1.0;
    CONSTANTS[AC_ICab_multiplier] = 1.0;
    CONSTANTS[AC_IClCa_multiplier] = 1.0;
    CONSTANTS[AC_IClb_multiplier] = 1.0;
    CONSTANTS[AC_IK1_multiplier] = 1.0;
    CONSTANTS[AC_IKb_multiplier] = 1.0;
    CONSTANTS[AC_IKr_multiplier] = 1.0;
    CONSTANTS[AC_IKs_multiplier] = 1.0;
    CONSTANTS[AC_INaCa_multiplier] = 1.0;
    CONSTANTS[AC_INaK_multiplier] = 1.0;
    CONSTANTS[AC_INaL_multiplier] = 1.0;
    CONSTANTS[AC_INa_multiplier] = 1.0;
    CONSTANTS[AC_INab_multiplier] = 1.0;
    CONSTANTS[AC_IpCa_multiplier] = 1.0;
    CONSTANTS[AC_Itof_multiplier] = 1.0;
    CONSTANTS[AC_Itos_multiplier] = 1.0;
    CONSTANTS[AC_Jrel_multiplier] = 1.0;
    CONSTANTS[AC_Jup_multiplier] = 1.0;
    
    /* physical_constants */
    CONSTANTS[AC_F] = 96485.0;
    CONSTANTS[AC_R] = 8314.0;
    CONSTANTS[AC_T] = 310.0;
    CONSTANTS[AC_zca] = 2.0;
    CONSTANTS[AC_zcl] = -1.0;
    CONSTANTS[AC_zk] = 1.0;
    CONSTANTS[AC_zna] = 1.0;
    CONSTANTS[AC_Qpow] = (CONSTANTS[AC_T] - 310.0) / 10.0;
    
    /* INaCa */
    CONSTANTS[AC_Gncx_b] = 0.00179;
    CONSTANTS[AC_INaCa_fractionSS] = 0.31;
    CONSTANTS[AC_KmCaAct] = 0.00015;
    CONSTANTS[AC_INaCa_celltype_factor] = ((tissueFlag == 1) ? 1.1 : ((tissueFlag == 2) ? 1.4 : 1.0));
    CONSTANTS[AC_kasymm] = 19.4258;
    CONSTANTS[AC_kcaoff] = 3853.2;
    CONSTANTS[AC_kcaon] = 3416400.0;
    CONSTANTS[AC_kna1] = 11.9712;
    CONSTANTS[AC_kna2] = 2.76;
    CONSTANTS[AC_kna3] = 88.767;
    CONSTANTS[AC_qca] = 0.0955;
    CONSTANTS[AC_qna] = 0.6718;
    CONSTANTS[AC_INaCa_sex_factor] = ((CONSTANTS[AC_sex] == 1.0) ? 1.0 / 1.0724 : ((CONSTANTS[AC_sex] == 2.0) ? 1.0724 : 1.0));
    CONSTANTS[AC_wca] = 51756.0;
    CONSTANTS[AC_wna] = 32978.0;
    CONSTANTS[AC_wnaca] = 2776.3;
    CONSTANTS[AC_Gncx] = CONSTANTS[AC_Gncx_b] * CONSTANTS[AC_INaCa_celltype_factor] * CONSTANTS[AC_INaCa_sex_factor] * CONSTANTS[AC_INaCa_multiplier];
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
    
    /* SERCA */
    CONSTANTS[AC_Km_SERCA_Ca] = 0.4;
    CONSTANTS[AC_Kmf] = 3.06720000000000018e-04;
    CONSTANTS[AC_Kmr] = 2.31442;
    CONSTANTS[AC_Max_Vmax_SERCA_Ca] = 1.11142;
    CONSTANTS[AC_Q10SRCaP] = 2.6;
    CONSTANTS[AC_Vmax_SRCaP_b] = 0.00543;
    CONSTANTS[AC_hillSRCaP] = 1.02809;
    CONSTANTS[AC_Kmf_p] = CONSTANTS[AC_Kmf] * 0.5;
    CONSTANTS[AC_Vmax_SRCaP] = ((tissueFlag == 1) ? 1.2 * CONSTANTS[AC_Vmax_SRCaP_b] * CONSTANTS[AC_Jup_multiplier] : CONSTANTS[AC_Vmax_SRCaP_b] * CONSTANTS[AC_Jup_multiplier]);
    
    /* contraction */
    CONSTANTS[AC_PKAForceMultiplier] = 1.0 + 0.26 * CONSTANTS[AC_fMyBPC_PKA];
    CONSTANTS[AC_TOT_A] = 25.0;
    CONSTANTS[AC_TRPN_n] = 1.65;
    CONSTANTS[AC_Tref] = 80.0;
    CONSTANTS[AC_beta_0] = 2.3;
    CONSTANTS[AC_contraction_beta_1] = -2.4;
    CONSTANTS[AC_dr] = 0.25;
    CONSTANTS[AC_fracTnIpo] = 0.0031;
    CONSTANTS[AC_contraction_gamma] = 0.0085;
    CONSTANTS[AC_gamma_wu] = 0.615;
    CONSTANTS[AC_koff] = 0.07854;
    CONSTANTS[AC_ktm_unblock] = 0.02626;
    CONSTANTS[AC_lambda] = 1.0;
    CONSTANTS[AC_lambda_max] = 1.2;
    CONSTANTS[AC_lambda_min] = 0.87;
    CONSTANTS[AC_lambda_rate] = 0.0;
    CONSTANTS[AC_mu] = 3.94046;
    CONSTANTS[AC_nperm] = 2.036;
    CONSTANTS[AC_nu] = 10.15996;
    CONSTANTS[AC_perm50] = 0.35;
    CONSTANTS[AC_phi] = 2.23;
    CONSTANTS[AC_wfrac] = 0.5;
    CONSTANTS[AC_A] = CONSTANTS[AC_TOT_A] * CONSTANTS[AC_dr] / ((1.0 - CONSTANTS[AC_dr]) * CONSTANTS[AC_wfrac] + CONSTANTS[AC_dr]);
    CONSTANTS[AC_XSSS] = CONSTANTS[AC_dr] * 0.5;
    CONSTANTS[AC_XWSS] = (1.0 - CONSTANTS[AC_dr]) * CONSTANTS[AC_wfrac] * 0.5;
    CONSTANTS[AC_fPKA_TnI] = 1.45 - 0.45 * (1.0 - CONSTANTS[AC_fTnI_PKA]) / (1.0 - CONSTANTS[AC_fracTnIpo]);
    CONSTANTS[AC_k_uw] = 0.026 * CONSTANTS[AC_nu];
    CONSTANTS[AC_k_ws] = 0.004 * (1.0 + CONSTANTS[AC_fMyBPC_PKA] / 2.0) * CONSTANTS[AC_mu];
    CONSTANTS[AC_lambda0] = ((CONSTANTS[AC_lambda_max] <= CONSTANTS[AC_lambda]) ? CONSTANTS[AC_lambda_max] : CONSTANTS[AC_lambda]);
    CONSTANTS[AC_Lfac] = ((0.0 >= 1.0 + CONSTANTS[AC_beta_0] * (CONSTANTS[AC_lambda0] + ((CONSTANTS[AC_lambda_min] <= CONSTANTS[AC_lambda0]) ? CONSTANTS[AC_lambda_min] : CONSTANTS[AC_lambda0]) - (1.0 + CONSTANTS[AC_lambda_min]))) ? 0.0 : 1.0 + CONSTANTS[AC_beta_0] * (CONSTANTS[AC_lambda0] + ((CONSTANTS[AC_lambda_min] <= CONSTANTS[AC_lambda0]) ? CONSTANTS[AC_lambda_min] : CONSTANTS[AC_lambda0]) - (1.0 + CONSTANTS[AC_lambda_min])));
    CONSTANTS[AC_ca50] = CONSTANTS[AC_fPKA_TnI] * 0.7645 + CONSTANTS[AC_contraction_beta_1] * ((0.2 <= CONSTANTS[AC_lambda] - 1.0) ? 0.2 : CONSTANTS[AC_lambda] - 1.0);
    CONSTANTS[AC_cds] = CONSTANTS[AC_phi] * CONSTANTS[AC_k_ws] * CONSTANTS[AC_wfrac] * (1.0 - CONSTANTS[AC_dr]) / CONSTANTS[AC_dr];
    CONSTANTS[AC_cdw] = CONSTANTS[AC_phi] * CONSTANTS[AC_k_uw] * (1.0 - CONSTANTS[AC_wfrac]) / CONSTANTS[AC_wfrac];
    CONSTANTS[AC_k_su] = CONSTANTS[AC_k_ws] * (1.0 / CONSTANTS[AC_dr] - 1.0) * CONSTANTS[AC_wfrac];
    CONSTANTS[AC_k_wu] = CONSTANTS[AC_k_uw] * (1.0 / CONSTANTS[AC_wfrac] - 1.0) - CONSTANTS[AC_k_ws];
    CONSTANTS[AC_ktm_block] = CONSTANTS[AC_ktm_unblock] * pow(CONSTANTS[AC_perm50], CONSTANTS[AC_nperm]) * 0.5 / (0.5 - CONSTANTS[AC_XSSS] - CONSTANTS[AC_XWSS]);
    
    /* reversal_potentials */
    CONSTANTS[AC_PNaK] = 0.01833;
    
    /* ICl */
    CONSTANTS[AC_GClCa_b] = 0.01615;
    CONSTANTS[AC_GClb_b] = 0.00241;
    CONSTANTS[AC_KdClCa] = 0.1;
    CONSTANTS[AC_GClCa] = CONSTANTS[AC_IClCa_multiplier] * CONSTANTS[AC_GClCa_b];
    CONSTANTS[AC_GClb] = CONSTANTS[AC_IClb_multiplier] * CONSTANTS[AC_GClb_b];
    
    /* IK1 */
    CONSTANTS[AC_GK1_b] = 0.6992;
    CONSTANTS[AC_IK1_celltype_factor] = ((tissueFlag == 1) ? 1.1 : ((tissueFlag == 2) ? 1.3 : 1.0));
    CONSTANTS[AC_IK1_sex_factor] = ((CONSTANTS[AC_sex] == 1.0) ? 1.07 : ((CONSTANTS[AC_sex] == 2.0) ? 0.92 : 1.0));
    CONSTANTS[AC_GK1] = CONSTANTS[AC_GK1_b] * CONSTANTS[AC_IK1_celltype_factor] * CONSTANTS[AC_IK1_sex_factor] * CONSTANTS[AC_IK1_multiplier];
    
    /* IKb */
    CONSTANTS[AC_GKb_b] = 0.010879;
    CONSTANTS[AC_GKb] = CONSTANTS[AC_GKb_b] * CONSTANTS[AC_IKb_multiplier];
    
    /* Ito */
    CONSTANTS[AC_Gto_fast] = ((tissueFlag == 1) ? 0.29856 * CONSTANTS[AC_Itof_multiplier] : ((tissueFlag == 2) ? 0.14928 * CONSTANTS[AC_Itof_multiplier] : 0.01276 * CONSTANTS[AC_Itof_multiplier]));
    CONSTANTS[AC_Gto_slow] = ((tissueFlag == 1) ? 0.02036 * CONSTANTS[AC_Itos_multiplier] : ((tissueFlag == 2) ? 0.04632 * CONSTANTS[AC_Itos_multiplier] : 0.0721 * CONSTANTS[AC_Itos_multiplier]));
    
    /* ICaL */
    CONSTANTS[AC_Aff] = 0.52477;
    CONSTANTS[AC_ICaL_fractionSS] = 0.8;
    CONSTANTS[AC_Io] = 0.5 * (CONSTANTS[AC_nao] + CONSTANTS[AC_ko] + CONSTANTS[AC_clo] + 4.0 * CONSTANTS[AC_cao]) / 1000.0;
    CONSTANTS[AC_Kmn] = 0.00222;
    CONSTANTS[AC_PCa_b] = 1.57679999999999995e-04;
    CONSTANTS[AC_dielConstant] = 74.0;
    CONSTANTS[AC_fICaLP] = CONSTANTS[AC_fICaL_PKA];
    CONSTANTS[AC_k2n] = 957.85903;
    CONSTANTS[AC_rateRecovery] = 0.02313;
    CONSTANTS[AC_tjca] = 66.0;
    CONSTANTS[AC_Afs] = 1.0 - CONSTANTS[AC_Aff];
    CONSTANTS[AC_PCa] = ((tissueFlag == 1) ? CONSTANTS[AC_PCa_b] * CONSTANTS[AC_ICaLPCa_multiplier] * 1.025 : ((tissueFlag == 2) ? CONSTANTS[AC_PCa_b] * CONSTANTS[AC_ICaLPCa_multiplier] * 1.1 : CONSTANTS[AC_PCa_b] * CONSTANTS[AC_ICaLPCa_multiplier]));
    CONSTANTS[AC_constA] = 1820000.0 * pow(CONSTANTS[AC_dielConstant] * CONSTANTS[AC_T], -1.5);
    CONSTANTS[AC_fICaL_P] = CONSTANTS[AC_fICaLP];
    CONSTANTS[AC_PCaK] = 2.21139954662870998e-04 * CONSTANTS[AC_PCa];
    CONSTANTS[AC_PCaNa] = 7.73432969581949454e-04 * CONSTANTS[AC_PCa];
    CONSTANTS[AC_PCa_Pb] = 1.9 * CONSTANTS[AC_PCa];
    CONSTANTS[AC_PCap] = 1.1 * CONSTANTS[AC_PCa];
    CONSTANTS[AC_gamma_cao] = pow(10.0, -CONSTANTS[AC_constA] * 4.0 * (sqrt(CONSTANTS[AC_Io]) / (1.0 + sqrt(CONSTANTS[AC_Io])) - 0.3 * CONSTANTS[AC_Io]));
    CONSTANTS[AC_gamma_ko] = pow(10.0, -CONSTANTS[AC_constA] * 1.0 * (sqrt(CONSTANTS[AC_Io]) / (1.0 + sqrt(CONSTANTS[AC_Io])) - 0.3 * CONSTANTS[AC_Io]));
    CONSTANTS[AC_gamma_nao] = pow(10.0, -CONSTANTS[AC_constA] * 1.0 * (sqrt(CONSTANTS[AC_Io]) / (1.0 + sqrt(CONSTANTS[AC_Io])) - 0.3 * CONSTANTS[AC_Io]));
    CONSTANTS[AC_PCaKp] = 2.21139954662870998e-04 * CONSTANTS[AC_PCap];
    CONSTANTS[AC_PCaNap] = 7.73432969581949454e-04 * CONSTANTS[AC_PCap];
    CONSTANTS[AC_PCa_P] = ((tissueFlag == 1) ? CONSTANTS[AC_PCa_Pb] * 1.025 : ((tissueFlag == 2) ? CONSTANTS[AC_PCa_Pb] * 1.1 : CONSTANTS[AC_PCa_Pb]));
    CONSTANTS[AC_PCaK_P] = 2.21139954662870998e-04 * CONSTANTS[AC_PCa_P];
    CONSTANTS[AC_PCaNa_P] = 7.73432969581949454e-04 * CONSTANTS[AC_PCa_P];
    
    /* ICab */
    CONSTANTS[AC_GCab_b] = 5.15574999999999996e-04;
    CONSTANTS[AC_GCab] = CONSTANTS[AC_ICab_multiplier] * CONSTANTS[AC_GCab_b];
    
    /* IKr */
    CONSTANTS[AC_GKr_b] = 0.043;
    CONSTANTS[AC_alpha_1] = 0.154375;
    CONSTANTS[AC_IKr_beta_1] = 0.1911;
    CONSTANTS[AC_IKr_celltype_factor] = ((tissueFlag == 1) ? 1.25 : ((tissueFlag == 2) ? 0.7 : 1.0));
    CONSTANTS[AC_IKr_sex_factor] = ((CONSTANTS[AC_sex] == 1.0) ? 1.1 : ((CONSTANTS[AC_sex] == 2.0) ? 0.87 : 1.0));
    CONSTANTS[AC_GKr] = CONSTANTS[AC_GKr_b] * CONSTANTS[AC_IKr_celltype_factor] * CONSTANTS[AC_IKr_sex_factor] * CONSTANTS[AC_IKr_multiplier];
    
    /* IKs */
    CONSTANTS[AC_GKs_b] = 2.97002;
    CONSTANTS[AC_IKs_celltype_factor] = ((tissueFlag == 1) ? 1.4 : ((tissueFlag == 2) ? 0.5 : 1.0));
    CONSTANTS[gnalTissueScale] = (tissueFlag == 1) ? 0.7 : 1.0;
    CONSTANTS[AC_gKs_factor] = 0.01;
    CONSTANTS[AC_kPKA_IKs] = CONSTANTS[AC_fIKs_PKA];
    CONSTANTS[AC_IKs_sex_factor] = ((CONSTANTS[AC_sex] == 1.0) ? 1.09 : ((CONSTANTS[AC_sex] == 2.0) ? 0.905 : 1.0));
    CONSTANTS[AC_GKs] = CONSTANTS[AC_GKs_b] * CONSTANTS[AC_IKs_celltype_factor] * CONSTANTS[AC_IKs_sex_factor] * CONSTANTS[AC_IKs_multiplier];
    CONSTANTS[AC_P_g_0] = CONSTANTS[AC_gKs_factor] * (0.2 + 0.2 * CONSTANTS[AC_kPKA_IKs]);
    CONSTANTS[AC_P_g_max] = CONSTANTS[AC_gKs_factor] * (0.8 + 7.0 * CONSTANTS[AC_kPKA_IKs]);
    CONSTANTS[AC_P_tau_0] = 26.0 + 9.0 * CONSTANTS[AC_kPKA_IKs];
    CONSTANTS[AC_P_tau_max] = 40.0 + 4.0 * CONSTANTS[AC_kPKA_IKs];
    CONSTANTS[AC_P_vh_0] = -1.0 - 10.0 * CONSTANTS[AC_kPKA_IKs];
    CONSTANTS[AC_P_vh_max] = -12.0 - 9.0 * CONSTANTS[AC_kPKA_IKs];
    
    /* INa */
    CONSTANTS[AC_GNa_b] = 22.08788;
    CONSTANTS[AC_fINa_P] = CONSTANTS[AC_fINa_PKA];
    CONSTANTS[AC_GNa] = CONSTANTS[AC_GNa_b] * CONSTANTS[AC_INa_multiplier];
    CONSTANTS[AC_GNa_P] = 1.25 * CONSTANTS[AC_GNa];
    
    /* INaK */
    CONSTANTS[AC_IbarNaK_b] = 2.10774;
    CONSTANTS[AC_KmKo] = 1.5;
    CONSTANTS[AC_KmNaip] = 11.0;
    CONSTANTS[AC_KmNaip_PKA] = 8.4615;
    CONSTANTS[AC_IbarNaK] = CONSTANTS[AC_INaK_multiplier] * CONSTANTS[AC_IbarNaK_b];
    
    /* INaL */
    CONSTANTS[AC_GNaL_b] = 0.04229;
    CONSTANTS[AC_thL] = 145.0;
    CONSTANTS[AC_thLp] = 3.0 * CONSTANTS[AC_thL];
    
    /* INab */
    CONSTANTS[AC_GNab_b] = 0.000594;
    CONSTANTS[AC_GNab] = CONSTANTS[AC_INab_multiplier] * CONSTANTS[AC_GNab_b];
    
    /* IpCa */
    CONSTANTS[AC_IbarSLCaP_b] = 0.02064;
    CONSTANTS[AC_KmPCa] = 0.0005;
    CONSTANTS[AC_Q10SLCaP] = 2.35;
    CONSTANTS[AC_IbarSLCaP] = ((CONSTANTS[AC_sex] == 1.0) ? CONSTANTS[AC_IbarSLCaP_b] * 0.8 * CONSTANTS[AC_IpCa_multiplier] : ((CONSTANTS[AC_sex] == 2.0) ? CONSTANTS[AC_IbarSLCaP_b] * 1.28 * CONSTANTS[AC_IpCa_multiplier] : CONSTANTS[AC_IbarSLCaP_b] * CONSTANTS[AC_IpCa_multiplier]));
    
    /* intracellular_ions */
    CONSTANTS[AC_Bmax_CaM] = ((CONSTANTS[AC_sex] == 1.0) ? 0.024 / 1.1 : ((CONSTANTS[AC_sex] == 2.0) ? 0.024 * 1.1 : 0.024));
    CONSTANTS[AC_Bmax_Csqn] = 1.36552139999999989e-01 * CONSTANTS[AC_vmyo] / CONSTANTS[AC_vsr];
    CONSTANTS[AC_Bmax_Naj] = 7.561;
    CONSTANTS[AC_Bmax_Nasl] = 1.65;
    CONSTANTS[AC_Bmax_SLhighj] = 1.75754999999999987e-04 * CONSTANTS[AC_vmyo] / CONSTANTS[AC_vjunc];
    CONSTANTS[AC_Bmax_SLhighsl] = 1.21542300000000002e-02 * CONSTANTS[AC_vmyo] / CONSTANTS[AC_vsl];
    CONSTANTS[AC_Bmax_SLlowj] = 4.89983000000000042e-04 * CONSTANTS[AC_vmyo] / CONSTANTS[AC_vjunc];
    CONSTANTS[AC_Bmax_SLlowsl] = 0.033923 * CONSTANTS[AC_vmyo] / CONSTANTS[AC_vsl];
    CONSTANTS[AC_Bmax_SR] = 1.78585399999999991e-02;
    CONSTANTS[AC_Bmax_TnChigh] = 0.14;
    CONSTANTS[AC_Bmax_TnClow] = 0.07;
    CONSTANTS[AC_Bmax_myosin] = 0.14;
    CONSTANTS[AC_J_ca_juncsl] = 3.26067463358168778e-13;
    CONSTANTS[AC_J_ca_slmyo] = 1.34127367348033708e-11;
    CONSTANTS[AC_J_na_juncsl] = 1.83127823220607987e-14;
    CONSTANTS[AC_J_na_slmyo] = 1.63862792221979393e-12;
    CONSTANTS[AC_koff_cam] = 0.238;
    CONSTANTS[AC_koff_csqn] = 65.0;
    CONSTANTS[AC_koff_myoca] = 0.00046;
    CONSTANTS[AC_koff_myomg] = 5.7e-05;
    CONSTANTS[AC_koff_na] = 0.001;
    CONSTANTS[AC_koff_slh] = 0.03;
    CONSTANTS[AC_koff_sll] = 1.3;
    CONSTANTS[AC_koff_sr] = 0.06;
    CONSTANTS[AC_koff_tnchca] = 3.2e-05;
    CONSTANTS[AC_koff_tnchmg] = 0.00333;
    CONSTANTS[AC_kon_cam] = 34.0;
    CONSTANTS[AC_kon_csqn] = 100.0;
    CONSTANTS[AC_kon_myoca] = 13.8;
    CONSTANTS[AC_kon_myomg] = 0.0157;
    CONSTANTS[AC_kon_na] = 0.0001;
    CONSTANTS[AC_kon_slh] = 100.0;
    CONSTANTS[AC_kon_sll] = 100.0;
    CONSTANTS[AC_kon_sr] = 100.0;
    CONSTANTS[AC_kon_tnchca] = 2.37;
    CONSTANTS[AC_kon_tnchmg] = 0.003;
    CONSTANTS[AC_mgi] = 0.5;
    
    /* membrane */
    CONSTANTS[AC_Cmem] = 1.381e-10;
    CONSTANTS[AC_Fjunc] = 0.11;
    CONSTANTS[AC_Fsl] = 1.0 - CONSTANTS[AC_Fjunc];
    
    /* ryr */
    CONSTANTS[AC_CI_to_RI] = 0.00248;
    CONSTANTS[AC_MaxSR] = 15.0;
    CONSTANTS[AC_MinSR] = 1.0;
    CONSTANTS[AC_baseRateCaI] = 3.02320000000000019e-04;
    CONSTANTS[AC_bt] = 12.4767;
    CONSTANTS[AC_caExpFactor] = 0.68655;
    CONSTANTS[AC_caExpFactor2] = 2.06273;
    CONSTANTS[AC_caTransFactor] = 0.94428;
    CONSTANTS[AC_caTransFactor2] = 0.52967;
    CONSTANTS[AC_directRelMidpoint] = 0.95271;
    CONSTANTS[AC_ec50SR] = 0.75385;
    CONSTANTS[AC_ecCaI] = 0.001;
    CONSTANTS[AC_kiCa] = 0.39871;
    CONSTANTS[AC_kim] = 0.04311;
    CONSTANTS[AC_koCa] = 23.87221;
    CONSTANTS[AC_kom] = 0.16219;
    CONSTANTS[AC_ks] = CONSTANTS[AC_Jrel_multiplier] * 26.6;
    CONSTANTS[AC_maxCaI] = 30.13294;
    CONSTANTS[AC_minCaI] = 0.93249;
    CONSTANTS[AC_steepnessCaI] = 5.93447;
    CONSTANTS[AC_steepnessCaSR] = 5.09473;
    CONSTANTS[AC_tauInact] = 64.11202;
    CONSTANTS[AC_tauInact2] = 119.48978;
    CONSTANTS[AC_a_rel] = 1.25 * CONSTANTS[AC_bt] / 1.0;
    CONSTANTS[AC_caTransFactor2p] = CONSTANTS[AC_caTransFactor2] * 1.5;

    /* --- Initial values --- */

STATES[v] = -8.70655895407835061e+01;
    STATES[camk_trap] =  4.14775695414328974e-02;
    STATES[camk_f_ICaL] =  1.49310462149601003e-01;
    STATES[camk_f_RyR] =  5.83425856770706033e-02;
    STATES[camk_f_PLB] =  1.02519516890396001e-02;
    STATES[casig_serca_trap] =  4.14775695414004025e-02;
    STATES[buffers_NaBj] =  3.55888399344769990e+00;
    STATES[buffers_NaBsl] =  7.76467094025442006e-01;
    STATES[buffers_TnClow] =  1.04491410305814006e-02;
    STATES[buffers_TnCHc] =  1.24415836833344007e-01;
    STATES[buffers_TnCHm] =  4.75139089328611031e-03;
    STATES[buffers_CaM] =  3.26761159683871009e-04;
    STATES[buffers_Myosin_ca] =  4.35725142843640997e-03;
    STATES[buffers_Myosin_mg] =  1.34587917639654997e-01;
    STATES[buffers_SRB] =  2.47759889601233997e-03;
    STATES[buffers_SLLj] =  4.15486833433121022e-03;
    STATES[buffers_SLLsl] =  8.10233315864154063e-03;
    STATES[buffers_SLHj] =  4.97605868295565024e-02;
    STATES[buffers_SLHsl] =  9.60185900338078047e-02;
    STATES[buffers_Csqn] =  1.32358809148371992e+00;
    STATES[naj] =  8.88424477947198987e+00;
    STATES[nasl] =  8.88341217642103054e+00;
    STATES[nai] =  8.88365804124162040e+00;
    STATES[ki] =  1.34213097948671987e+02;
    STATES[cli] =  1.55745509675962008e+01;
    STATES[casr] =  7.09622552333154966e-01;
    STATES[caj] =  9.20580586542241966e-05;
    STATES[casl] =  9.62428218317655005e-05;
    STATES[cai] =  9.66105402693763000e-05;
    STATES[m] =  1.15997496446638006e-03;
    STATES[h] =  7.91660456106285992e-01;
    STATES[j] =  7.91587530878067991e-01;
    STATES[hp] =  6.12409807882460022e-01;
    STATES[jp] =  7.91052149435448992e-01;
    STATES[m_P] =  3.34028683421829994e-03;
    STATES[h_P] =  6.47375310489835987e-01;
    STATES[j_P] =  6.47279109941896968e-01;
    STATES[hp_P] =  4.19158183843564014e-01;
    STATES[jp_P] =  4.18669666153645981e-01;
    STATES[mL] =  2.24902633151065996e-04;
    STATES[hL] =  4.78572465113524026e-01;
    STATES[hLp] =  2.58793079229037992e-01;
    STATES[d] =  6.00992839633468041e-26;
    STATES[ff] =  9.99999988243592042e-01;
    STATES[fs] =  9.11094950010905014e-01;
    STATES[fcaf] =  9.99999988251085048e-01;
    STATES[fcas] =  9.91477881027688968e-01;
    STATES[jca] =  9.99993204834444982e-01;
    STATES[ffp] =  9.99999987868131046e-01;
    STATES[fcafp] =  9.99999987980700999e-01;
    STATES[nca] =  5.28788208273222975e-03;
    STATES[nca_i] =  6.21547539490541028e-03;
    STATES[d_P] = 0.0;
    STATES[ff_P] =  9.99999987992091000e-01;
    STATES[fs_P] =  9.30411822666438959e-01;
    STATES[fcaf_P] =  9.99999987985881966e-01;
    STATES[fcas_P] =  9.99145398038995047e-01;
    STATES[fBPf] =  9.99999987992091000e-01;
    STATES[fcaBPf] =  9.99999987985881966e-01;
    STATES[ical_pureCDI_junc] =  9.98342275864489981e-01;
    STATES[ical_pureCDI_sl] =  9.98263185911480000e-01;
    STATES[xtos] =  2.86101402392893993e-04;
    STATES[ytos] =  7.68184915918178945e-01;
    STATES[xtof] =  2.86099374494749989e-04;
    STATES[ytof] =  9.99998646875180031e-01;
    STATES[xtos_p] =  1.32590983159417000e-04;
    STATES[xtof_p] =  1.32590043205026988e-04;
    STATES[ytos_p] =  7.77117031645023992e-01;
    STATES[ytof_p] =  9.99998646879750042e-01;
    STATES[C0] =  9.97558791218494001e-01;
    STATES[C1] =  9.60122333286174043e-04;
    STATES[C2] =  8.04204051120587025e-04;
    STATES[I] =  3.55901681542090999e-05;
    STATES[O] =  6.41292278097939979e-04;
    STATES[xs_junc] =  3.11875010014238986e-02;
    STATES[xs_sl] =  3.12048105743730000e-02;
    STATES[jrel_icaldep_act] = -9.23949112019060968e-15;
    STATES[jrel_icaldep_f1] =  9.99994847409610976e-01;
    STATES[jrel_icaldep_f2] =  9.98725931920675047e-01;
    STATES[ryr_R] =  8.21525956842405036e-01;
    STATES[ryr_O] =  3.33929041065478997e-08;
    STATES[ryr_I] =  4.57362306340243018e-09;
    STATES[ryr_CaRI] =  6.53313479138347974e-02;
    STATES[ryr_R_p] =  8.21663945287199993e-01;
    STATES[ryr_O_p] =  5.00978004918836978e-08;
    STATES[ryr_I_p] =  6.86118926495390029e-09;
    STATES[ryr_CaRI_p] =  6.51828099424918983e-02;
    STATES[contraction_TmBlocked] =  9.92378492940848012e-01;
    STATES[contraction_XW] =  2.85370717430549021e-03;
    STATES[contraction_XS] =  1.91507884532543006e-03;
    STATES[contraction_ZETAS] = 0.0;
    STATES[contraction_ZETAW] = 0.0;
    STATES[contraction_Ca_TRPN] =  3.19778784167467980e-02;
}


void
TWorldcomputeVariables(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC, bool solveVmWithinODESolver, const Foam::StimulusProtocol& stimulus)
{
/* CaMK */
    ALGEBRAIC[AV_bound] = CONSTANTS[AC_CaMK0] * (1.0 - STATES[camk_trap]) / (1.0 + CONSTANTS[AC_Km_CaMK_Ca] / STATES[caj]);
    ALGEBRAIC[AV_bound_serca] = CONSTANTS[AC_CaMK0] * (1.0 - STATES[casig_serca_trap]) / (1.0 + CONSTANTS[AC_Km_CaMK_Ca] / STATES[caj]);
    ALGEBRAIC[AV_CaMK_active] = ALGEBRAIC[AV_bound] + STATES[camk_trap];
    ALGEBRAIC[AV_casig_SERCA_act] = ALGEBRAIC[AV_bound_serca] + STATES[casig_serca_trap];
    ALGEBRAIC[AV_CaMK_Phos_ss_ICaL] = ALGEBRAIC[AV_CaMK_active] / (ALGEBRAIC[AV_CaMK_active] + CONSTANTS[AC_K_Phos_CaMK]);
    ALGEBRAIC[AV_CaMK_Phos_ss_PLB] = ALGEBRAIC[AV_CaMK_active] / (ALGEBRAIC[AV_CaMK_active] + 10.0);
    ALGEBRAIC[AV_CaMK_Phos_ss_RyR] = ALGEBRAIC[AV_CaMK_active] / (ALGEBRAIC[AV_CaMK_active] + 1.0);
    RATES[camk_trap] = CONSTANTS[AC_CaMK_alpha] * ALGEBRAIC[AV_bound] * ALGEBRAIC[AV_CaMK_active] - CONSTANTS[AC_CaMK_beta] * STATES[camk_trap] * (0.1 + 0.9 * CONSTANTS[AC_PP1_tot] / 0.1371);
    RATES[casig_serca_trap] = CONSTANTS[AC_alpha_serca] * ALGEBRAIC[AV_bound_serca] * ALGEBRAIC[AV_casig_SERCA_act] - CONSTANTS[AC_CaMK_beta] * STATES[casig_serca_trap] * (0.1 + 0.9 * CONSTANTS[AC_PP1_tot] / 0.1371);
    RATES[camk_f_ICaL] = (ALGEBRAIC[AV_CaMK_Phos_ss_ICaL] - STATES[camk_f_ICaL]) / CONSTANTS[AC_tau_cal];
    RATES[camk_f_PLB] = (ALGEBRAIC[AV_CaMK_Phos_ss_PLB] - STATES[camk_f_PLB]) / CONSTANTS[AC_tau_plb];
    RATES[camk_f_RyR] = (ALGEBRAIC[AV_CaMK_Phos_ss_RyR] - STATES[camk_f_RyR]) / CONSTANTS[AC_tau_ryr];
    
    /* environment */
    ALGEBRAIC[AV_time] = VOI;
    
    /* INaCa */
    ALGEBRAIC[AV_allo_i] = 1.0 / (1.0 + pow(CONSTANTS[AC_KmCaAct] / STATES[casl], 2.0));
    ALGEBRAIC[AV_allo_ss] = 1.0 / (1.0 + pow(CONSTANTS[AC_KmCaAct] / STATES[caj], 2.0));
    ALGEBRAIC[AV_h4_i] = 1.0 + STATES[nai] / CONSTANTS[AC_kna1] * (1.0 + STATES[nai] / CONSTANTS[AC_kna2]);
    ALGEBRAIC[AV_h4_ss] = 1.0 + STATES[nasl] / CONSTANTS[AC_kna1] * (1.0 + STATES[nasl] / CONSTANTS[AC_kna2]);
    ALGEBRAIC[AV_hca] = exp(CONSTANTS[AC_qca] * (STATES[v] - 8.3117) * CONSTANTS[AC_F] / (CONSTANTS[AC_R] * CONSTANTS[AC_T]));
    ALGEBRAIC[AV_hna] = exp(CONSTANTS[AC_qna] * (STATES[v] - 8.3117) * CONSTANTS[AC_F] / (CONSTANTS[AC_R] * CONSTANTS[AC_T]));
    ALGEBRAIC[AV_h1_i] = 1.0 + STATES[nai] / CONSTANTS[AC_kna3] * (1.0 + ALGEBRAIC[AV_hna]);
    ALGEBRAIC[AV_h1_ss] = 1.0 + STATES[nasl] / CONSTANTS[AC_kna3] * (1.0 + ALGEBRAIC[AV_hna]);
    ALGEBRAIC[AV_h5_i] = STATES[nai] * STATES[nai] / (ALGEBRAIC[AV_h4_i] * CONSTANTS[AC_kna1] * CONSTANTS[AC_kna2]);
    ALGEBRAIC[AV_h5_ss] = STATES[nasl] * STATES[nasl] / (ALGEBRAIC[AV_h4_ss] * CONSTANTS[AC_kna1] * CONSTANTS[AC_kna2]);
    ALGEBRAIC[AV_h6_i] = 1.0 / ALGEBRAIC[AV_h4_i];
    ALGEBRAIC[AV_h6_ss] = 1.0 / ALGEBRAIC[AV_h4_ss];
    ALGEBRAIC[AV_h7_i] = 1.0 + CONSTANTS[AC_nao] / CONSTANTS[AC_kna3] * (1.0 + 1.0 / ALGEBRAIC[AV_hna]);
    ALGEBRAIC[AV_h7_ss] = 1.0 + CONSTANTS[AC_nao] / CONSTANTS[AC_kna3] * (1.0 + 1.0 / ALGEBRAIC[AV_hna]);
    ALGEBRAIC[AV_h2_i] = STATES[nai] * ALGEBRAIC[AV_hna] / (CONSTANTS[AC_kna3] * ALGEBRAIC[AV_h1_i]);
    ALGEBRAIC[AV_h2_ss] = STATES[nasl] * ALGEBRAIC[AV_hna] / (CONSTANTS[AC_kna3] * ALGEBRAIC[AV_h1_ss]);
    ALGEBRAIC[AV_h3_i] = 1.0 / ALGEBRAIC[AV_h1_i];
    ALGEBRAIC[AV_h3_ss] = 1.0 / ALGEBRAIC[AV_h1_ss];
    ALGEBRAIC[AV_h8_i] = CONSTANTS[AC_nao] / (CONSTANTS[AC_kna3] * ALGEBRAIC[AV_hna] * ALGEBRAIC[AV_h7_i]);
    ALGEBRAIC[AV_h8_ss] = CONSTANTS[AC_nao] / (CONSTANTS[AC_kna3] * ALGEBRAIC[AV_hna] * ALGEBRAIC[AV_h7_ss]);
    ALGEBRAIC[AV_h9_i] = 1.0 / ALGEBRAIC[AV_h7_i];
    ALGEBRAIC[AV_h9_ss] = 1.0 / ALGEBRAIC[AV_h7_ss];
    ALGEBRAIC[AV_k6_i] = ALGEBRAIC[AV_h6_i] * STATES[cai] * CONSTANTS[AC_kcaon];
    ALGEBRAIC[AV_k6_ss] = ALGEBRAIC[AV_h6_ss] * STATES[casl] * CONSTANTS[AC_kcaon];
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
    ALGEBRAIC[AV_INaCa_i] = (1.0 - CONSTANTS[AC_INaCa_fractionSS]) * CONSTANTS[AC_Gncx] * ALGEBRAIC[AV_allo_i] * (CONSTANTS[AC_zna] * ALGEBRAIC[AV_JncxNa_i] + CONSTANTS[AC_zca] * ALGEBRAIC[AV_JncxCa_i]);
    ALGEBRAIC[AV_INaCa_ss] = CONSTANTS[AC_INaCa_fractionSS] * CONSTANTS[AC_Gncx] * ALGEBRAIC[AV_allo_ss] * (CONSTANTS[AC_zna] * ALGEBRAIC[AV_JncxNa_ss] + CONSTANTS[AC_zca] * ALGEBRAIC[AV_JncxCa_ss]);
    ALGEBRAIC[AV_INaCa_INaCa] = ALGEBRAIC[AV_INaCa_i] + ALGEBRAIC[AV_INaCa_ss];
    
    /* SERCA */
    ALGEBRAIC[AV_phosphorylationTotal] = STATES[camk_f_PLB] + CONSTANTS[AC_fPLB_PKA] - STATES[camk_f_PLB] * CONSTANTS[AC_fPLB_PKA];
    ALGEBRAIC[AV_Vmax_mult] = 1.0 + CONSTANTS[AC_Max_Vmax_SERCA_Ca] / (1.0 + pow(CONSTANTS[AC_Km_SERCA_Ca] / ALGEBRAIC[AV_casig_SERCA_act], 2.0));
    ALGEBRAIC[AV_Jserca_np] = pow(CONSTANTS[AC_Q10SRCaP], CONSTANTS[AC_Qpow]) * CONSTANTS[AC_Vmax_SRCaP] * ALGEBRAIC[AV_Vmax_mult] * (pow(STATES[cai] / CONSTANTS[AC_Kmf], CONSTANTS[AC_hillSRCaP]) - pow(STATES[casr] / CONSTANTS[AC_Kmr], CONSTANTS[AC_hillSRCaP])) / (1.0 + pow(STATES[cai] / CONSTANTS[AC_Kmf], CONSTANTS[AC_hillSRCaP]) + pow(STATES[casr] / CONSTANTS[AC_Kmr], CONSTANTS[AC_hillSRCaP]));
    ALGEBRAIC[AV_Jserca_p] = pow(CONSTANTS[AC_Q10SRCaP], CONSTANTS[AC_Qpow]) * CONSTANTS[AC_Vmax_SRCaP] * ALGEBRAIC[AV_Vmax_mult] * (pow(STATES[cai] / CONSTANTS[AC_Kmf_p], CONSTANTS[AC_hillSRCaP]) - pow(STATES[casr] / CONSTANTS[AC_Kmr], CONSTANTS[AC_hillSRCaP])) / (1.0 + pow(STATES[cai] / CONSTANTS[AC_Kmf_p], CONSTANTS[AC_hillSRCaP]) + pow(STATES[casr] / CONSTANTS[AC_Kmr], CONSTANTS[AC_hillSRCaP]));
    ALGEBRAIC[AV_Jserca] = ALGEBRAIC[AV_Jserca_np] * (1.0 - ALGEBRAIC[AV_phosphorylationTotal]) + ALGEBRAIC[AV_Jserca_p] * ALGEBRAIC[AV_phosphorylationTotal];
    
    /* contraction */
    ALGEBRAIC[AV_XU] = 1.0 - STATES[contraction_TmBlocked] - STATES[contraction_XW] - STATES[contraction_XS];
    ALGEBRAIC[AV_gamma_rate] = ((STATES[contraction_ZETAS] > 0.0) ? CONSTANTS[AC_contraction_gamma] * STATES[contraction_ZETAS] : ((STATES[contraction_ZETAS] < -1.0) ? CONSTANTS[AC_contraction_gamma] * (-STATES[contraction_ZETAS] - 1.0) : 0.0));
    ALGEBRAIC[AV_gamma_rate_w] = CONSTANTS[AC_gamma_wu] * fabs(STATES[contraction_ZETAW]);
    ALGEBRAIC[AV_xb_su_gamma] = ALGEBRAIC[AV_gamma_rate] * STATES[contraction_XS];
    ALGEBRAIC[AV_xb_uw] = CONSTANTS[AC_k_uw] * ALGEBRAIC[AV_XU] * (1.0 + CONSTANTS[AC_fMyBPC_PKA] / 2.0);
    ALGEBRAIC[AV_xb_ws] = CONSTANTS[AC_k_ws] * STATES[contraction_XW];
    ALGEBRAIC[AV_xb_wu_gamma] = ALGEBRAIC[AV_gamma_rate_w] * STATES[contraction_XW];
    ALGEBRAIC[AV_Ta] = CONSTANTS[AC_Lfac] * CONSTANTS[AC_PKAForceMultiplier] * CONSTANTS[AC_Tref] / CONSTANTS[AC_dr] * ((STATES[contraction_ZETAS] + 1.0) * STATES[contraction_XS] + STATES[contraction_ZETAW] * STATES[contraction_XW]);
    ALGEBRAIC[AV_d_contraction_Ca_TRPN] = CONSTANTS[AC_koff] * (pow(STATES[cai] * 1000.0 / CONSTANTS[AC_ca50], CONSTANTS[AC_TRPN_n]) * (1.0 - STATES[contraction_Ca_TRPN]) - STATES[contraction_Ca_TRPN]);
    ALGEBRAIC[AV_xb_su] = CONSTANTS[AC_k_su] * STATES[contraction_XS];
    ALGEBRAIC[AV_xb_wu] = CONSTANTS[AC_k_wu] * STATES[contraction_XW];
    RATES[contraction_TmBlocked] = CONSTANTS[AC_ktm_block] * ((100.0 <= pow(STATES[contraction_Ca_TRPN], -CONSTANTS[AC_nperm] / 2.0)) ? 100.0 : pow(STATES[contraction_Ca_TRPN], -CONSTANTS[AC_nperm] / 2.0)) * ALGEBRAIC[AV_XU] - CONSTANTS[AC_ktm_unblock] * pow(STATES[contraction_Ca_TRPN], CONSTANTS[AC_nperm] / 2.0) * STATES[contraction_TmBlocked];
    RATES[contraction_ZETAS] = CONSTANTS[AC_A] * CONSTANTS[AC_lambda_rate] - CONSTANTS[AC_cds] * STATES[contraction_ZETAS];
    RATES[contraction_ZETAW] = CONSTANTS[AC_A] * CONSTANTS[AC_lambda_rate] - CONSTANTS[AC_cdw] * STATES[contraction_ZETAW];
    RATES[contraction_Ca_TRPN] = ALGEBRAIC[AV_d_contraction_Ca_TRPN];
    RATES[contraction_XS] = ALGEBRAIC[AV_xb_ws] - ALGEBRAIC[AV_xb_su] - ALGEBRAIC[AV_xb_su_gamma];
    RATES[contraction_XW] = ALGEBRAIC[AV_xb_uw] - ALGEBRAIC[AV_xb_wu] - ALGEBRAIC[AV_xb_ws] - ALGEBRAIC[AV_xb_wu_gamma];
    
    /* reversal_potentials */
    ALGEBRAIC[AV_ECaj] = CONSTANTS[AC_R] * CONSTANTS[AC_T] / (CONSTANTS[AC_zca] * CONSTANTS[AC_F]) * log(CONSTANTS[AC_cao] / STATES[caj]);
    ALGEBRAIC[AV_ECasl] = CONSTANTS[AC_R] * CONSTANTS[AC_T] / (CONSTANTS[AC_zca] * CONSTANTS[AC_F]) * log(CONSTANTS[AC_cao] / STATES[casl]);
    ALGEBRAIC[AV_ECl] = CONSTANTS[AC_R] * CONSTANTS[AC_T] / (CONSTANTS[AC_zcl] * CONSTANTS[AC_F]) * log(CONSTANTS[AC_clo] / STATES[cli]);
    ALGEBRAIC[AV_EK] = CONSTANTS[AC_R] * CONSTANTS[AC_T] / (CONSTANTS[AC_zk] * CONSTANTS[AC_F]) * log(CONSTANTS[AC_ko] / STATES[ki]);
    ALGEBRAIC[AV_ENaj] = CONSTANTS[AC_R] * CONSTANTS[AC_T] / (CONSTANTS[AC_zna] * CONSTANTS[AC_F]) * log(CONSTANTS[AC_nao] / STATES[naj]);
    ALGEBRAIC[AV_ENasl] = CONSTANTS[AC_R] * CONSTANTS[AC_T] / (CONSTANTS[AC_zna] * CONSTANTS[AC_F]) * log(CONSTANTS[AC_nao] / STATES[nasl]);
    ALGEBRAIC[AV_EKs] = CONSTANTS[AC_R] * CONSTANTS[AC_T] / (CONSTANTS[AC_zk] * CONSTANTS[AC_F]) * log((CONSTANTS[AC_ko] + CONSTANTS[AC_PNaK] * CONSTANTS[AC_nao]) / (STATES[ki] + CONSTANTS[AC_PNaK] * STATES[nai]));
    
    /* ICl */
    ALGEBRAIC[AV_IClCa_junc] = 0.5 * CONSTANTS[AC_GClCa] / (1.0 + CONSTANTS[AC_KdClCa] / STATES[caj]) * (STATES[v] - ALGEBRAIC[AV_ECl]);
    ALGEBRAIC[AV_IClCa_sl] = 0.5 * CONSTANTS[AC_GClCa] / (1.0 + CONSTANTS[AC_KdClCa] / STATES[casl]) * (STATES[v] - ALGEBRAIC[AV_ECl]);
    ALGEBRAIC[AV_IClb] = CONSTANTS[AC_GClb] * (STATES[v] - ALGEBRAIC[AV_ECl]);
    ALGEBRAIC[AV_IClCa] = ALGEBRAIC[AV_IClCa_junc] + ALGEBRAIC[AV_IClCa_sl];
    
    /* IK1 */
    ALGEBRAIC[AV_aK1] = 4.094 / (1.0 + exp(0.1217 * (STATES[v] - ALGEBRAIC[AV_EK] - 49.934)));
    ALGEBRAIC[AV_bK1] = (15.72 * exp(0.0674 * (STATES[v] - ALGEBRAIC[AV_EK] - 3.257)) + exp(0.0618 * (STATES[v] - ALGEBRAIC[AV_EK] - 594.31))) / (1.0 + exp(-0.1629 * (STATES[v] - ALGEBRAIC[AV_EK] + 14.207)));
    ALGEBRAIC[AV_K1ss] = ALGEBRAIC[AV_aK1] / (ALGEBRAIC[AV_aK1] + ALGEBRAIC[AV_bK1]);
    ALGEBRAIC[AV_IK1_IK1] = CONSTANTS[AC_GK1] * sqrt(CONSTANTS[AC_ko] / 5.0) * ALGEBRAIC[AV_K1ss] * (STATES[v] - ALGEBRAIC[AV_EK]);
    
    /* IKb */
    ALGEBRAIC[AV_xkb] = 1.0 / (1.0 + exp(-(STATES[v] - 10.8968) / 23.9871));
    ALGEBRAIC[AV_IKb_IKb] = CONSTANTS[AC_GKb] * ALGEBRAIC[AV_xkb] * (STATES[v] - ALGEBRAIC[AV_EK]);
    
    /* Ito */
    ALGEBRAIC[AV_dti_develop] = 1.354 + 0.0001 / (exp((STATES[v] - 167.4) / 15.89) + exp(-(STATES[v] - 12.23) / 0.2154));
    ALGEBRAIC[AV_dti_recover] = 1.0 - 0.5 / (1.0 + exp((STATES[v] + 70.0) / 20.0));
    ALGEBRAIC[AV_fItop] = STATES[camk_f_RyR];
    ALGEBRAIC[AV_tauxtof] = 8.5 * exp(-pow((STATES[v] + 45.0) / 50.0, 2.0)) + 0.5;
    ALGEBRAIC[AV_tauxtos] = 9.0 / (1.0 + exp((STATES[v] + 3.0) / 15.0)) + 0.5;
    ALGEBRAIC[AV_tauytof] = 85.0 * exp(-pow(STATES[v] + 40.0, 2.0) / 220.0) + 7.0;
    ALGEBRAIC[AV_tauytos] = 800.0 / (1.0 + exp((STATES[v] + 60.0) / 10.0)) + 30.0;
    ALGEBRAIC[AV_xtoss] = 1.0 / (1.0 + exp(-(STATES[v] - 19.0) / 13.0));
    ALGEBRAIC[AV_xtoss_p] = 1.0 / (1.0 + exp(-(STATES[v] - 29.0) / 13.0));
    ALGEBRAIC[AV_ytoss] = 1.0 / (1.0 + exp((STATES[v] + 19.5) / 5.0));
    ALGEBRAIC[AV_Itof] = CONSTANTS[AC_Gto_fast] * (STATES[v] - ALGEBRAIC[AV_EK]) * ((1.0 - ALGEBRAIC[AV_fItop]) * STATES[xtof] * STATES[ytof] + ALGEBRAIC[AV_fItop] * STATES[xtof_p] * STATES[ytof_p]);
    ALGEBRAIC[AV_Itos] = CONSTANTS[AC_Gto_slow] * (STATES[v] - ALGEBRAIC[AV_EK]) * ((1.0 - ALGEBRAIC[AV_fItop]) * STATES[xtos] * STATES[ytos] + ALGEBRAIC[AV_fItop] * STATES[xtos_p] * STATES[ytos_p]);
    ALGEBRAIC[AV_tauytof_p] = ALGEBRAIC[AV_tauytof] * ALGEBRAIC[AV_dti_develop] * ALGEBRAIC[AV_dti_recover];
    ALGEBRAIC[AV_tauytos_p] = ALGEBRAIC[AV_tauytos] * ALGEBRAIC[AV_dti_develop] * ALGEBRAIC[AV_dti_recover];
    RATES[xtof] = (ALGEBRAIC[AV_xtoss] - STATES[xtof]) / ALGEBRAIC[AV_tauxtof];
    RATES[xtof_p] = (ALGEBRAIC[AV_xtoss_p] - STATES[xtof_p]) / ALGEBRAIC[AV_tauxtof];
    RATES[xtos] = (ALGEBRAIC[AV_xtoss] - STATES[xtos]) / ALGEBRAIC[AV_tauxtos];
    RATES[xtos_p] = (ALGEBRAIC[AV_xtoss_p] - STATES[xtos_p]) / ALGEBRAIC[AV_tauxtos];
    RATES[ytof] = (ALGEBRAIC[AV_ytoss] - STATES[ytof]) / ALGEBRAIC[AV_tauytof];
    RATES[ytos] = (ALGEBRAIC[AV_ytoss] - STATES[ytos]) / ALGEBRAIC[AV_tauytos];
    ALGEBRAIC[AV_Ito_Ito] = ALGEBRAIC[AV_Itos] + ALGEBRAIC[AV_Itof];
    RATES[ytof_p] = (ALGEBRAIC[AV_ytoss] - STATES[ytof_p]) / ALGEBRAIC[AV_tauytof_p];
    RATES[ytos_p] = (ALGEBRAIC[AV_ytoss] - STATES[ytos_p]) / ALGEBRAIC[AV_tauytos_p];
    
    /* ICaL */
    ALGEBRAIC[AV_Afcaf] = 0.3 + 0.6 / (1.0 + exp((STATES[v] - 9.24247) / 27.96201));
    ALGEBRAIC[AV_Ii] = 0.5 * (STATES[nai] + STATES[ki] + STATES[cli] + 4.0 * STATES[cai]) / 1000.0;
    ALGEBRAIC[AV_Iss] = 0.5 * (STATES[nasl] + STATES[ki] + STATES[cli] + 4.0 * STATES[casl]) / 1000.0;
    ALGEBRAIC[AV_dPss] = ((1.0323 * exp(-1.0553 * exp(-0.081 * (STATES[v] + 12.62483))) <= 1.0) ? 1.0323 * exp(-1.0553 * exp(-0.081 * (STATES[v] + 12.62483))) : 1.0);
    ALGEBRAIC[AV_dss] = ((1.0763 * exp(-1.007 * exp(-0.0829 * (STATES[v] + 3.62483))) <= 1.0) ? 1.0763 * exp(-1.007 * exp(-0.0829 * (STATES[v] + 3.62483))) : 1.0);
    ALGEBRAIC[AV_fICaLp] = STATES[camk_f_ICaL];
    ALGEBRAIC[AV_fss] = 1.0 / (1.0 + exp((STATES[v] + 19.58) / 3.696));
    ALGEBRAIC[AV_fss_P] = 1.0 / (1.0 + exp((STATES[v] + 25.58) / 3.696));
    ALGEBRAIC[AV_jcass] = 1.0 / (1.0 + exp((STATES[v] + 17.66945) / 3.21501));
    ALGEBRAIC[AV_km2n] = STATES[jca] * 0.84191;
    ALGEBRAIC[AV_sigmoidTransition] = 1.0 - 1.0 / (1.0 + 1.86532 * STATES[caj] / 0.032);
    ALGEBRAIC[AV_sigmoidTransition2] = 1.0 - 1.0 / (1.0 + 1.86532 * STATES[casl] / 0.032);
    ALGEBRAIC[AV_td] = 1.5 + 1.0 / (exp(-0.05 * (STATES[v] + 6.0)) + exp(0.09 * (STATES[v] + 14.0)));
    ALGEBRAIC[AV_tfcaf] = 13.50673 + 1.0 / (0.1542 * exp(-(STATES[v] - 1.31611) / 11.3396) + 0.1542 * exp((STATES[v] - 1.31611) / 11.3396));
    ALGEBRAIC[AV_tfcas] = 177.95813 + 1.0 / (4.73955000000000025e-04 * exp(-(STATES[v] - 0.79049) / 0.81777) + 4.73955000000000025e-04 * exp((STATES[v] + 2.40474) / 1.90812));
    ALGEBRAIC[AV_tff] = 6.17111 + 1.0 / (0.00126 * exp(-(STATES[v] + 26.63596) / 9.69961) + 0.00126 * exp((STATES[v] + 26.63596) / 9.69961));
    ALGEBRAIC[AV_tfs] = 2.71922488999999996e+03 + 1.0 / (7.19411000000000014e-05 * exp(-(STATES[v] + 5.74631) / 10.8769) + 7.19411000000000014e-05 * exp((STATES[v] + 5.74631) / 16.31535));
    ALGEBRAIC[AV_Afcas] = 1.0 - ALGEBRAIC[AV_Afcaf];
    ALGEBRAIC[AV_anca] = 1.0 / (CONSTANTS[AC_k2n] / ALGEBRAIC[AV_km2n] + pow(1.0 + CONSTANTS[AC_Kmn] / STATES[caj], 3.80763));
    ALGEBRAIC[AV_anca_i] = 1.0 / (CONSTANTS[AC_k2n] / ALGEBRAIC[AV_km2n] + pow(1.0 + CONSTANTS[AC_Kmn] / STATES[casl], 3.80763));
    ALGEBRAIC[AV_fBPss] = ALGEBRAIC[AV_fss_P];
    ALGEBRAIC[AV_fcass] = ALGEBRAIC[AV_fss];
    ALGEBRAIC[AV_fcass_P] = ALGEBRAIC[AV_fss_P];
    ALGEBRAIC[AV_tauTransition] = 1.0967 + (1.0 - ALGEBRAIC[AV_sigmoidTransition]) * 141.4299;
    ALGEBRAIC[AV_tauTransition2] = 1.0967 + (1.0 - ALGEBRAIC[AV_sigmoidTransition2]) * 141.4299;
    ALGEBRAIC[AV_tfcafp] = 2.5 * ALGEBRAIC[AV_tfcaf];
    ALGEBRAIC[AV_tffp] = 2.5 * ALGEBRAIC[AV_tff];
    RATES[d] = (ALGEBRAIC[AV_dss] - STATES[d]) / ALGEBRAIC[AV_td];
    RATES[d_P] = (ALGEBRAIC[AV_dPss] - STATES[d_P]) / ALGEBRAIC[AV_td];
    RATES[ff] = (ALGEBRAIC[AV_fss] - STATES[ff]) / ALGEBRAIC[AV_tff];
    RATES[ff_P] = (ALGEBRAIC[AV_fss_P] - STATES[ff_P]) / ALGEBRAIC[AV_tff];
    RATES[fs] = (ALGEBRAIC[AV_fss] - STATES[fs]) / ALGEBRAIC[AV_tfs];
    RATES[fs_P] = (ALGEBRAIC[AV_fss_P] - STATES[fs_P]) / ALGEBRAIC[AV_tfs];
    RATES[jca] = (ALGEBRAIC[AV_jcass] - STATES[jca]) / CONSTANTS[AC_tjca];
    ALGEBRAIC[AV_f] = CONSTANTS[AC_Aff] * STATES[ff] + CONSTANTS[AC_Afs] * STATES[fs];
    ALGEBRAIC[AV_fBP] = CONSTANTS[AC_Aff] * STATES[fBPf] + CONSTANTS[AC_Afs] * STATES[fs_P];
    ALGEBRAIC[AV_fICaL_BP] = ALGEBRAIC[AV_fICaLp] * CONSTANTS[AC_fICaL_P];
    ALGEBRAIC[AV_f_P] = CONSTANTS[AC_Aff] * STATES[ff_P] + CONSTANTS[AC_Afs] * STATES[fs_P];
    ALGEBRAIC[AV_fca] = ALGEBRAIC[AV_Afcaf] * STATES[fcaf] + ALGEBRAIC[AV_Afcas] * STATES[fcas];
    ALGEBRAIC[AV_fcaBP] = ALGEBRAIC[AV_Afcaf] * STATES[fcaBPf] + ALGEBRAIC[AV_Afcas] * STATES[fcas_P];
    ALGEBRAIC[AV_fcaBPss] = ALGEBRAIC[AV_fcass_P];
    ALGEBRAIC[AV_fcap] = ALGEBRAIC[AV_Afcaf] * STATES[fcafp] + ALGEBRAIC[AV_Afcas] * STATES[fcas];
    ALGEBRAIC[AV_fcap_P] = ALGEBRAIC[AV_Afcaf] * STATES[fcaf_P] + ALGEBRAIC[AV_Afcas] * STATES[fcas_P];
    ALGEBRAIC[AV_fp] = CONSTANTS[AC_Aff] * STATES[ffp] + CONSTANTS[AC_Afs] * STATES[fs];
    ALGEBRAIC[AV_gamma_cai] = pow(10.0, -CONSTANTS[AC_constA] * 4.0 * (sqrt(ALGEBRAIC[AV_Ii]) / (1.0 + sqrt(ALGEBRAIC[AV_Ii])) - 0.3 * ALGEBRAIC[AV_Ii]));
    ALGEBRAIC[AV_gamma_cass] = pow(10.0, -CONSTANTS[AC_constA] * 4.0 * (sqrt(ALGEBRAIC[AV_Iss]) / (1.0 + sqrt(ALGEBRAIC[AV_Iss])) - 0.3 * ALGEBRAIC[AV_Iss]));
    ALGEBRAIC[AV_gamma_ki] = pow(10.0, -CONSTANTS[AC_constA] * 1.0 * (sqrt(ALGEBRAIC[AV_Ii]) / (1.0 + sqrt(ALGEBRAIC[AV_Ii])) - 0.3 * ALGEBRAIC[AV_Ii]));
    ALGEBRAIC[AV_gamma_kss] = pow(10.0, -CONSTANTS[AC_constA] * 1.0 * (sqrt(ALGEBRAIC[AV_Iss]) / (1.0 + sqrt(ALGEBRAIC[AV_Iss])) - 0.3 * ALGEBRAIC[AV_Iss]));
    ALGEBRAIC[AV_gamma_nai] = pow(10.0, -CONSTANTS[AC_constA] * 1.0 * (sqrt(ALGEBRAIC[AV_Ii]) / (1.0 + sqrt(ALGEBRAIC[AV_Ii])) - 0.3 * ALGEBRAIC[AV_Ii]));
    ALGEBRAIC[AV_gamma_nass] = pow(10.0, -CONSTANTS[AC_constA] * 1.0 * (sqrt(ALGEBRAIC[AV_Iss]) / (1.0 + sqrt(ALGEBRAIC[AV_Iss])) - 0.3 * ALGEBRAIC[AV_Iss]));
    RATES[fBPf] = (ALGEBRAIC[AV_fBPss] - STATES[fBPf]) / ALGEBRAIC[AV_tffp];
    RATES[fcaf] = (ALGEBRAIC[AV_fcass] - STATES[fcaf]) / ALGEBRAIC[AV_tfcaf];
    RATES[fcaf_P] = (ALGEBRAIC[AV_fcass_P] - STATES[fcaf_P]) / ALGEBRAIC[AV_tfcaf];
    RATES[fcafp] = (ALGEBRAIC[AV_fcass] - STATES[fcafp]) / ALGEBRAIC[AV_tfcafp];
    RATES[fcas] = (ALGEBRAIC[AV_fcass] - STATES[fcas]) / ALGEBRAIC[AV_tfcas];
    RATES[fcas_P] = (ALGEBRAIC[AV_fcass_P] - STATES[fcas_P]) / ALGEBRAIC[AV_tfcas];
    RATES[ffp] = (ALGEBRAIC[AV_fss] - STATES[ffp]) / ALGEBRAIC[AV_tffp];
    RATES[ical_pureCDI_junc] = -STATES[ical_pureCDI_junc] * ALGEBRAIC[AV_sigmoidTransition] / ALGEBRAIC[AV_tauTransition] + (1.0 - STATES[ical_pureCDI_junc]) * CONSTANTS[AC_rateRecovery];
    RATES[ical_pureCDI_sl] = -STATES[ical_pureCDI_sl] * ALGEBRAIC[AV_sigmoidTransition2] / ALGEBRAIC[AV_tauTransition2] + (1.0 - STATES[ical_pureCDI_sl]) * CONSTANTS[AC_rateRecovery];
    RATES[nca] = ALGEBRAIC[AV_anca] * CONSTANTS[AC_k2n] - STATES[nca] * ALGEBRAIC[AV_km2n];
    RATES[nca_i] = ALGEBRAIC[AV_anca_i] * CONSTANTS[AC_k2n] - STATES[nca_i] * ALGEBRAIC[AV_km2n];
    ALGEBRAIC[AV_fICaL_CaMKonly] = ALGEBRAIC[AV_fICaLp] - ALGEBRAIC[AV_fICaL_BP];
    ALGEBRAIC[AV_fICaL_PKAonly] = CONSTANTS[AC_fICaL_P] - ALGEBRAIC[AV_fICaL_BP];
    RATES[fcaBPf] = (ALGEBRAIC[AV_fcaBPss] - STATES[fcaBPf]) / ALGEBRAIC[AV_tfcafp];
    
    /* IKr */
    ALGEBRAIC[AV_IKr_IKr] = CONSTANTS[AC_GKr] * sqrt(CONSTANTS[AC_ko] / 5.0) * STATES[O] * (STATES[v] - ALGEBRAIC[AV_EK]);
    
    /* IKs */
    ALGEBRAIC[AV_VsTs_Ca_junc] = CONSTANTS[AC_P_tau_0] + (CONSTANTS[AC_P_tau_max] - CONSTANTS[AC_P_tau_0]) / (1.0 + pow(0.00015 / STATES[caj], 3.0));
    ALGEBRAIC[AV_VsTs_Ca_sl] = CONSTANTS[AC_P_tau_0] + (CONSTANTS[AC_P_tau_max] - CONSTANTS[AC_P_tau_0]) / (1.0 + pow(0.00015 / STATES[casl], 3.0));
    ALGEBRAIC[AV_VsXs_Ca_junc] = CONSTANTS[AC_P_vh_0] + (CONSTANTS[AC_P_vh_max] - CONSTANTS[AC_P_vh_0]) / (1.0 + pow(0.00035 / STATES[caj], 4.0));
    ALGEBRAIC[AV_VsXs_Ca_sl] = CONSTANTS[AC_P_vh_0] + (CONSTANTS[AC_P_vh_max] - CONSTANTS[AC_P_vh_0]) / (1.0 + pow(0.00035 / STATES[casl], 4.0));
    ALGEBRAIC[AV_gks_junc] = CONSTANTS[AC_P_g_0] + (CONSTANTS[AC_P_g_max] - CONSTANTS[AC_P_g_0]) / (1.0 + pow(0.00015 / STATES[caj], 1.3));
    ALGEBRAIC[AV_gks_sl] = CONSTANTS[AC_P_g_0] + (CONSTANTS[AC_P_g_max] - CONSTANTS[AC_P_g_0]) / (1.0 + pow(0.00015 / STATES[casl], 1.3));
    ALGEBRAIC[AV_tauxs_junc] = 2.0 * (50.0 + (50.0 + 350.0 * exp(-pow(STATES[v] + 30.0, 2.0) / 4000.0)) / (1.0 + exp(-(STATES[v] + ALGEBRAIC[AV_VsTs_Ca_junc]) / 10.0)));
    ALGEBRAIC[AV_tauxs_sl] = 2.0 * (50.0 + (50.0 + 350.0 * exp(-pow(STATES[v] + 30.0, 2.0) / 4000.0)) / (1.0 + exp(-(STATES[v] + ALGEBRAIC[AV_VsTs_Ca_sl]) / 10.0)));
    ALGEBRAIC[AV_xsss_junc] = 1.0 / (1.0 + exp(-(STATES[v] - ALGEBRAIC[AV_VsXs_Ca_junc]) / 25.0));
    ALGEBRAIC[AV_xsss_sl] = 1.0 / (1.0 + exp(-(STATES[v] - ALGEBRAIC[AV_VsXs_Ca_sl]) / 25.0));
    RATES[xs_junc] = (ALGEBRAIC[AV_xsss_junc] - STATES[xs_junc]) / ALGEBRAIC[AV_tauxs_junc];
    RATES[xs_sl] = (ALGEBRAIC[AV_xsss_sl] - STATES[xs_sl]) / ALGEBRAIC[AV_tauxs_sl];
    
    /* INa */
    ALGEBRAIC[AV_ah] = ((STATES[v] >= -40.0) ? 0.0 : 0.057 * exp(-(STATES[v] + 80.0) / 6.8));
    ALGEBRAIC[AV_aj] = ((STATES[v] >= -40.0) ? 0.0 : (-25428.0 * exp(0.2444 * STATES[v]) - 6.948e-06 * exp(-0.04391 * STATES[v])) * (STATES[v] + 37.78) / (1.0 + exp(0.311 * (STATES[v] + 79.23))));
    ALGEBRAIC[AV_bh] = ((STATES[v] >= -40.0) ? 0.77 / (0.13 * (1.0 + exp(-(STATES[v] + 10.66) / 11.1))) : 2.7 * exp(0.079 * STATES[v]) + 310000.0 * exp(0.3485 * STATES[v]));
    ALGEBRAIC[AV_bj] = ((STATES[v] >= -40.0) ? 0.6 * exp(0.057 * STATES[v]) / (1.0 + exp(-0.1 * (STATES[v] + 32.0))) : 0.02424 * exp(-0.01052 * STATES[v]) / (1.0 + exp(-0.1378 * (STATES[v] + 40.14))));
    ALGEBRAIC[AV_fINap] = STATES[camk_f_RyR];
    ALGEBRAIC[AV_hss] = 1.0 / pow(1.0 + exp((STATES[v] + 71.55) / 7.43), 2.0);
    ALGEBRAIC[AV_hss_P] = 1.0 / pow(1.0 + exp((STATES[v] + 76.55) / 7.43), 2.0);
    ALGEBRAIC[AV_hssp] = 1.0 / pow(1.0 + exp((STATES[v] + 77.55) / 7.43), 2.0);
    ALGEBRAIC[AV_hssp_P] = 1.0 / pow(1.0 + exp((STATES[v] + 82.55) / 7.43), 2.0);
    ALGEBRAIC[AV_mss] = 1.0 / pow(1.0 + exp(-(STATES[v] + 56.86) / 9.03), 2.0);
    ALGEBRAIC[AV_mss_P] = 1.0 / pow(1.0 + exp(-(STATES[v] + 61.86) / 9.03), 2.0);
    ALGEBRAIC[AV_tm] = 0.1292 * exp(-pow((STATES[v] + 45.79) / 15.54, 2.0)) + 0.06487 * exp(-pow((STATES[v] - 4.823) / 51.12, 2.0));
    ALGEBRAIC[AV_fINa_BP] = ALGEBRAIC[AV_fINap] * CONSTANTS[AC_fINa_P];
    ALGEBRAIC[AV_jss] = ALGEBRAIC[AV_hss];
    ALGEBRAIC[AV_jss_P] = ALGEBRAIC[AV_hss_P];
    ALGEBRAIC[AV_jssp_P] = ALGEBRAIC[AV_hssp_P];
    ALGEBRAIC[AV_th] = 1.0 / (ALGEBRAIC[AV_ah] + ALGEBRAIC[AV_bh]);
    ALGEBRAIC[AV_tj] = 1.0 / (ALGEBRAIC[AV_aj] + ALGEBRAIC[AV_bj]);
    RATES[m] = (ALGEBRAIC[AV_mss] - STATES[m]) / ALGEBRAIC[AV_tm];
    RATES[m_P] = (ALGEBRAIC[AV_mss_P] - STATES[m_P]) / ALGEBRAIC[AV_tm];
    ALGEBRAIC[AV_INaBase_CaMK] = CONSTANTS[AC_GNa] * pow(STATES[m], 3.0) * STATES[hp] * STATES[jp];
    ALGEBRAIC[AV_INaBase_NP] = CONSTANTS[AC_GNa] * pow(STATES[m], 3.0) * STATES[h] * STATES[j];
    ALGEBRAIC[AV_fINa_CaMKonly] = ALGEBRAIC[AV_fINap] - ALGEBRAIC[AV_fINa_BP];
    ALGEBRAIC[AV_fINa_PKAonly] = CONSTANTS[AC_fINa_P] - ALGEBRAIC[AV_fINa_BP];
    ALGEBRAIC[AV_tjp] = 1.46 * ALGEBRAIC[AV_tj];
    RATES[h] = (ALGEBRAIC[AV_hss] - STATES[h]) / ALGEBRAIC[AV_th];
    RATES[h_P] = (ALGEBRAIC[AV_hss_P] - STATES[h_P]) / ALGEBRAIC[AV_th];
    RATES[hp] = (ALGEBRAIC[AV_hssp] - STATES[hp]) / ALGEBRAIC[AV_th];
    RATES[hp_P] = (ALGEBRAIC[AV_hssp_P] - STATES[hp_P]) / ALGEBRAIC[AV_th];
    RATES[j] = (ALGEBRAIC[AV_jss] - STATES[j]) / ALGEBRAIC[AV_tj];
    RATES[j_P] = (ALGEBRAIC[AV_jss_P] - STATES[j_P]) / ALGEBRAIC[AV_tj];
    ALGEBRAIC[AV_INaBase_BP] = CONSTANTS[AC_GNa_P] * pow(STATES[m_P], 3.0) * STATES[hp_P] * STATES[jp_P];
    ALGEBRAIC[AV_INaBase_PKA] = CONSTANTS[AC_GNa_P] * pow(STATES[m_P], 3.0) * STATES[h_P] * STATES[j_P];
    RATES[jp] = (ALGEBRAIC[AV_jss] - STATES[jp]) / ALGEBRAIC[AV_tjp];
    RATES[jp_P] = (ALGEBRAIC[AV_jssp_P] - STATES[jp_P]) / ALGEBRAIC[AV_tjp];
    ALGEBRAIC[AV_INaBase] = (1.0 - ALGEBRAIC[AV_fINa_CaMKonly] - ALGEBRAIC[AV_fINa_PKAonly] - ALGEBRAIC[AV_fINa_BP]) * ALGEBRAIC[AV_INaBase_NP] + ALGEBRAIC[AV_fINa_CaMKonly] * ALGEBRAIC[AV_INaBase_CaMK] + ALGEBRAIC[AV_fINa_PKAonly] * ALGEBRAIC[AV_INaBase_PKA] + ALGEBRAIC[AV_fINa_BP] * ALGEBRAIC[AV_INaBase_BP];
    
    /* INaK */
    ALGEBRAIC[AV_fnak] = 0.75 + (0.00375 - (140.0 - CONSTANTS[AC_nao]) / 50.0 * 0.001) * STATES[v];
    
    /* INaL */
    ALGEBRAIC[AV_fINaLp] = STATES[camk_f_RyR];
    ALGEBRAIC[AV_hLss] = 1.0 / (1.0 + exp((STATES[v] + 87.61) / 7.488));
    ALGEBRAIC[AV_hLssp] = 1.0 / (1.0 + exp((STATES[v] + 93.81) / 7.488));
    ALGEBRAIC[AV_mLss] = 1.0 / (1.0 + exp(-(STATES[v] + 42.85) / 5.264));
    ALGEBRAIC[AV_tmL] = 0.1292 * exp(-pow((STATES[v] + 45.79) / 15.54, 2.0)) + 0.06487 * exp(-pow((STATES[v] - 4.823) / 51.12, 2.0));
    // AV_GNaL: epi-only 0.7 factor lives in CONSTANTS[gnalTissueScale]
    // (1.0 in non-epi). Branch-free hot path.
    ALGEBRAIC[AV_GNaL] = CONSTANTS[AC_GNaL_b] * CONSTANTS[AC_INaL_multiplier]
                       * (1.0 + ALGEBRAIC[AV_fINaLp]) * CONSTANTS[gnalTissueScale];
    RATES[hL] = (ALGEBRAIC[AV_hLss] - STATES[hL]) / CONSTANTS[AC_thL];
    RATES[mL] = (ALGEBRAIC[AV_mLss] - STATES[mL]) / ALGEBRAIC[AV_tmL];
    RATES[hLp] = (ALGEBRAIC[AV_hLssp] - STATES[hLp]) / CONSTANTS[AC_thLp];
    
    /* intracellular_ions */
    ALGEBRAIC[AV_ICl_tot] = ALGEBRAIC[AV_IClCa] + ALGEBRAIC[AV_IClb];
    RATES[buffers_Myosin_mg] = CONSTANTS[AC_kon_myomg] * CONSTANTS[AC_mgi] * (CONSTANTS[AC_Bmax_myosin] - STATES[buffers_Myosin_ca] - STATES[buffers_Myosin_mg]) - CONSTANTS[AC_koff_myomg] * STATES[buffers_Myosin_mg];
    RATES[buffers_TnCHm] = CONSTANTS[AC_kon_tnchmg] * CONSTANTS[AC_mgi] * (CONSTANTS[AC_Bmax_TnChigh] - STATES[buffers_TnCHc] - STATES[buffers_TnCHm]) - CONSTANTS[AC_koff_tnchmg] * STATES[buffers_TnCHm];
    RATES[nai] = CONSTANTS[AC_J_na_slmyo] / CONSTANTS[AC_vmyo] * (STATES[nasl] - STATES[nai]);
    ALGEBRAIC[AV_d_buffers_CaM] = CONSTANTS[AC_kon_cam] * STATES[cai] * (CONSTANTS[AC_Bmax_CaM] - STATES[buffers_CaM]) - CONSTANTS[AC_koff_cam] * STATES[buffers_CaM];
    ALGEBRAIC[AV_d_buffers_Csqn] = CONSTANTS[AC_kon_csqn] * STATES[casr] * (CONSTANTS[AC_Bmax_Csqn] - STATES[buffers_Csqn]) - CONSTANTS[AC_koff_csqn] * STATES[buffers_Csqn];
    ALGEBRAIC[AV_d_buffers_Myosin_ca] = CONSTANTS[AC_kon_myoca] * STATES[cai] * (CONSTANTS[AC_Bmax_myosin] - STATES[buffers_Myosin_ca] - STATES[buffers_Myosin_mg]) - CONSTANTS[AC_koff_myoca] * STATES[buffers_Myosin_ca];
    ALGEBRAIC[AV_d_buffers_NaBj] = CONSTANTS[AC_kon_na] * STATES[naj] * (CONSTANTS[AC_Bmax_Naj] - STATES[buffers_NaBj]) - CONSTANTS[AC_koff_na] * STATES[buffers_NaBj];
    ALGEBRAIC[AV_d_buffers_NaBsl] = CONSTANTS[AC_kon_na] * STATES[nasl] * (CONSTANTS[AC_Bmax_Nasl] - STATES[buffers_NaBsl]) - CONSTANTS[AC_koff_na] * STATES[buffers_NaBsl];
    ALGEBRAIC[AV_d_buffers_SLHj] = CONSTANTS[AC_kon_slh] * STATES[caj] * (CONSTANTS[AC_Bmax_SLhighj] - STATES[buffers_SLHj]) - CONSTANTS[AC_koff_slh] * STATES[buffers_SLHj];
    ALGEBRAIC[AV_d_buffers_SLHsl] = CONSTANTS[AC_kon_slh] * STATES[casl] * (CONSTANTS[AC_Bmax_SLhighsl] - STATES[buffers_SLHsl]) - CONSTANTS[AC_koff_slh] * STATES[buffers_SLHsl];
    ALGEBRAIC[AV_d_buffers_SLLj] = CONSTANTS[AC_kon_sll] * STATES[caj] * (CONSTANTS[AC_Bmax_SLlowj] - STATES[buffers_SLLj]) - CONSTANTS[AC_koff_sll] * STATES[buffers_SLLj];
    ALGEBRAIC[AV_d_buffers_SLLsl] = CONSTANTS[AC_kon_sll] * STATES[casl] * (CONSTANTS[AC_Bmax_SLlowsl] - STATES[buffers_SLLsl]) - CONSTANTS[AC_koff_sll] * STATES[buffers_SLLsl];
    ALGEBRAIC[AV_d_buffers_SRB] = CONSTANTS[AC_kon_sr] * STATES[cai] * (CONSTANTS[AC_Bmax_SR] - STATES[buffers_SRB]) - CONSTANTS[AC_koff_sr] * STATES[buffers_SRB];
    ALGEBRAIC[AV_d_buffers_TnCHc] = CONSTANTS[AC_kon_tnchca] * STATES[cai] * (CONSTANTS[AC_Bmax_TnChigh] - STATES[buffers_TnCHc] - STATES[buffers_TnCHm]) - CONSTANTS[AC_koff_tnchca] * STATES[buffers_TnCHc];
    ALGEBRAIC[AV_d_buffers_TnClow] = CONSTANTS[AC_Bmax_TnClow] * ALGEBRAIC[AV_d_contraction_Ca_TRPN];
    RATES[buffers_CaM] = ALGEBRAIC[AV_d_buffers_CaM];
    RATES[buffers_Csqn] = ALGEBRAIC[AV_d_buffers_Csqn];
    RATES[buffers_Myosin_ca] = ALGEBRAIC[AV_d_buffers_Myosin_ca];
    RATES[buffers_NaBj] = ALGEBRAIC[AV_d_buffers_NaBj];
    RATES[buffers_NaBsl] = ALGEBRAIC[AV_d_buffers_NaBsl];
    RATES[buffers_SLHj] = ALGEBRAIC[AV_d_buffers_SLHj];
    RATES[buffers_SLHsl] = ALGEBRAIC[AV_d_buffers_SLHsl];
    RATES[buffers_SLLj] = ALGEBRAIC[AV_d_buffers_SLLj];
    RATES[buffers_SLLsl] = ALGEBRAIC[AV_d_buffers_SLLsl];
    RATES[buffers_SRB] = ALGEBRAIC[AV_d_buffers_SRB];
    RATES[buffers_TnCHc] = ALGEBRAIC[AV_d_buffers_TnCHc];
    RATES[buffers_TnClow] = ALGEBRAIC[AV_d_buffers_TnClow];
    ALGEBRAIC[AV_J_CaB_cytosol] = ALGEBRAIC[AV_d_buffers_TnClow] + ALGEBRAIC[AV_d_buffers_TnCHc] + ALGEBRAIC[AV_d_buffers_CaM] + ALGEBRAIC[AV_d_buffers_Myosin_ca] + ALGEBRAIC[AV_d_buffers_SRB];
    ALGEBRAIC[AV_J_CaB_junction] = ALGEBRAIC[AV_d_buffers_SLLj] + ALGEBRAIC[AV_d_buffers_SLHj];
    ALGEBRAIC[AV_J_CaB_sl] = ALGEBRAIC[AV_d_buffers_SLLsl] + ALGEBRAIC[AV_d_buffers_SLHsl];
    RATES[cai] = -ALGEBRAIC[AV_Jserca] * CONSTANTS[AC_vsr] / CONSTANTS[AC_vmyo] - ALGEBRAIC[AV_J_CaB_cytosol] + CONSTANTS[AC_J_ca_slmyo] / CONSTANTS[AC_vmyo] * (STATES[casl] - STATES[cai]);
    
    /* membrane */
    ALGEBRAIC[AV_vffrt] = STATES[v] * CONSTANTS[AC_F] * CONSTANTS[AC_F] / (CONSTANTS[AC_R] * CONSTANTS[AC_T]);
    ALGEBRAIC[AV_vfrt] = STATES[v] * CONSTANTS[AC_F] / (CONSTANTS[AC_R] * CONSTANTS[AC_T]);
    
    /* ryr */
    ALGEBRAIC[AV_CaMKIILeakMultiplier] = 1.0 + 2.0 * STATES[camk_f_RyR];
    ALGEBRAIC[AV_Jrel_ICaLdep] = 0.00174 * STATES[jrel_icaldep_act] * STATES[jrel_icaldep_f1] * STATES[jrel_icaldep_f2];
    ALGEBRAIC[AV_RIcleft] = 1.0 - STATES[ryr_R] - STATES[ryr_O] - STATES[ryr_I] - STATES[ryr_CaRI];
    ALGEBRAIC[AV_RIcleft_p] = 1.0 - STATES[ryr_R_p] - STATES[ryr_O_p] - STATES[ryr_I_p] - STATES[ryr_CaRI_p];
    ALGEBRAIC[AV_nonlinearModifier] = 0.2144 * exp(1.83 * STATES[casr]);
    ALGEBRAIC[AV_J_SRCarel_np] = CONSTANTS[AC_ks] * STATES[ryr_O] * (STATES[casr] - STATES[caj]) + ALGEBRAIC[AV_Jrel_ICaLdep];
    ALGEBRAIC[AV_J_SRCarel_p] = CONSTANTS[AC_ks] * STATES[ryr_O_p] * (STATES[casr] - STATES[caj]) + ALGEBRAIC[AV_Jrel_ICaLdep];
    ALGEBRAIC[AV_J_SRleak] = 1.59306000000000000e-06 * (STATES[casr] - STATES[caj]) * ALGEBRAIC[AV_nonlinearModifier] * ALGEBRAIC[AV_CaMKIILeakMultiplier];
    ALGEBRAIC[AV_kCaSR] = CONSTANTS[AC_MaxSR] - (CONSTANTS[AC_MaxSR] - CONSTANTS[AC_MinSR]) / (1.0 + pow(CONSTANTS[AC_ec50SR] / STATES[casr], CONSTANTS[AC_steepnessCaSR]));
    ALGEBRAIC[AV_sigmoidBaseCaI] = CONSTANTS[AC_minCaI] + (CONSTANTS[AC_maxCaI] - CONSTANTS[AC_minCaI]) / (1.0 + pow(CONSTANTS[AC_ecCaI] / STATES[caj], CONSTANTS[AC_steepnessCaI]));
    ALGEBRAIC[AV_tau_rel] = ((CONSTANTS[AC_bt] / (1.0 + 0.0123 / STATES[casr]) >= 0.001) ? CONSTANTS[AC_bt] / (1.0 + 0.0123 / STATES[casr]) : 0.001);
    ALGEBRAIC[AV_J_SRCarel] = ALGEBRAIC[AV_J_SRCarel_p] * STATES[camk_f_RyR] + ALGEBRAIC[AV_J_SRCarel_np] * (1.0 - STATES[camk_f_RyR]);
    ALGEBRAIC[AV_RI_to_CI] = CONSTANTS[AC_baseRateCaI] * ALGEBRAIC[AV_sigmoidBaseCaI];
    ALGEBRAIC[AV_kiSRCa] = CONSTANTS[AC_kiCa] * ALGEBRAIC[AV_kCaSR];
    ALGEBRAIC[AV_koSRCa] = CONSTANTS[AC_koCa] / ALGEBRAIC[AV_kCaSR];
    RATES[ryr_CaRI] = ALGEBRAIC[AV_RI_to_CI] * ALGEBRAIC[AV_RIcleft] - CONSTANTS[AC_CI_to_RI] * STATES[ryr_CaRI];
    RATES[ryr_CaRI_p] = ALGEBRAIC[AV_RI_to_CI] * ALGEBRAIC[AV_RIcleft_p] - CONSTANTS[AC_CI_to_RI] * STATES[ryr_CaRI_p];
    RATES[ryr_I] = ALGEBRAIC[AV_kiSRCa] * CONSTANTS[AC_caTransFactor] * pow(STATES[caj] / 1.0, CONSTANTS[AC_caExpFactor]) * STATES[ryr_O] - CONSTANTS[AC_kim] * STATES[ryr_I] - (CONSTANTS[AC_kom] * STATES[ryr_I] - CONSTANTS[AC_caTransFactor2] * ALGEBRAIC[AV_koSRCa] * pow(STATES[caj] / 1.0, CONSTANTS[AC_caExpFactor2]) * ALGEBRAIC[AV_RIcleft]);
    RATES[ryr_I_p] = ALGEBRAIC[AV_kiSRCa] * CONSTANTS[AC_caTransFactor] * pow(STATES[caj] / 1.0, CONSTANTS[AC_caExpFactor]) * STATES[ryr_O_p] - CONSTANTS[AC_kim] * STATES[ryr_I_p] - (CONSTANTS[AC_kom] * STATES[ryr_I_p] - CONSTANTS[AC_caTransFactor2p] * ALGEBRAIC[AV_koSRCa] * pow(STATES[caj] / 1.0, CONSTANTS[AC_caExpFactor2]) * ALGEBRAIC[AV_RIcleft_p]);
    RATES[ryr_O] = CONSTANTS[AC_caTransFactor2] * ALGEBRAIC[AV_koSRCa] * pow(STATES[caj] / 1.0, CONSTANTS[AC_caExpFactor2]) * STATES[ryr_R] - CONSTANTS[AC_kom] * STATES[ryr_O] - (ALGEBRAIC[AV_kiSRCa] * CONSTANTS[AC_caTransFactor] * pow(STATES[caj] / 1.0, CONSTANTS[AC_caExpFactor]) * STATES[ryr_O] - CONSTANTS[AC_kim] * STATES[ryr_I]);
    RATES[ryr_O_p] = CONSTANTS[AC_caTransFactor2p] * ALGEBRAIC[AV_koSRCa] * pow(STATES[caj] / 1.0, CONSTANTS[AC_caExpFactor2]) * STATES[ryr_R_p] - CONSTANTS[AC_kom] * STATES[ryr_O_p] - (ALGEBRAIC[AV_kiSRCa] * CONSTANTS[AC_caTransFactor] * pow(STATES[caj] / 1.0, CONSTANTS[AC_caExpFactor]) * STATES[ryr_O_p] - CONSTANTS[AC_kim] * STATES[ryr_I_p]);
    RATES[ryr_R] = CONSTANTS[AC_kim] * ALGEBRAIC[AV_RIcleft] - ALGEBRAIC[AV_kiSRCa] * CONSTANTS[AC_caTransFactor] * pow(STATES[caj] / 1.0, CONSTANTS[AC_caExpFactor]) * STATES[ryr_R] - (CONSTANTS[AC_caTransFactor2] * ALGEBRAIC[AV_koSRCa] * pow(STATES[caj] / 1.0, CONSTANTS[AC_caExpFactor2]) * STATES[ryr_R] - CONSTANTS[AC_kom] * STATES[ryr_O]);
    RATES[ryr_R_p] = CONSTANTS[AC_kim] * ALGEBRAIC[AV_RIcleft_p] - ALGEBRAIC[AV_kiSRCa] * CONSTANTS[AC_caTransFactor] * pow(STATES[caj] / 1.0, CONSTANTS[AC_caExpFactor]) * STATES[ryr_R_p] - (CONSTANTS[AC_caTransFactor2p] * ALGEBRAIC[AV_koSRCa] * pow(STATES[caj] / 1.0, CONSTANTS[AC_caExpFactor2]) * STATES[ryr_R_p] - CONSTANTS[AC_kom] * STATES[ryr_O_p]);
    
    /* *remaining* */
    ALGEBRAIC[AV_PhiCaK_i] = 1.0 * ALGEBRAIC[AV_vffrt] * (ALGEBRAIC[AV_gamma_ki] * STATES[ki] * exp(1.0 * ALGEBRAIC[AV_vfrt]) - CONSTANTS[AC_gamma_ko] * CONSTANTS[AC_ko]) / (exp(1.0 * ALGEBRAIC[AV_vfrt]) - 1.0);
    ALGEBRAIC[AV_PhiCaK_ss] = 1.0 * ALGEBRAIC[AV_vffrt] * (ALGEBRAIC[AV_gamma_kss] * STATES[ki] * exp(1.0 * ALGEBRAIC[AV_vfrt]) - CONSTANTS[AC_gamma_ko] * CONSTANTS[AC_ko]) / (exp(1.0 * ALGEBRAIC[AV_vfrt]) - 1.0);
    ALGEBRAIC[AV_PhiCaL_i] = 4.0 * ALGEBRAIC[AV_vffrt] * (ALGEBRAIC[AV_gamma_cai] * STATES[cai] * exp(2.0 * ALGEBRAIC[AV_vfrt]) - CONSTANTS[AC_gamma_cao] * CONSTANTS[AC_cao]) / (exp(2.0 * ALGEBRAIC[AV_vfrt]) - 1.0);
    ALGEBRAIC[AV_PhiCaL_ss] = 4.0 * ALGEBRAIC[AV_vffrt] * (ALGEBRAIC[AV_gamma_cass] * STATES[casl] * exp(2.0 * ALGEBRAIC[AV_vfrt]) - CONSTANTS[AC_gamma_cao] * CONSTANTS[AC_cao]) / (exp(2.0 * ALGEBRAIC[AV_vfrt]) - 1.0);
    ALGEBRAIC[AV_PhiCaNa_i] = 1.0 * ALGEBRAIC[AV_vffrt] * (ALGEBRAIC[AV_gamma_nai] * STATES[nai] * exp(1.0 * ALGEBRAIC[AV_vfrt]) - CONSTANTS[AC_gamma_nao] * CONSTANTS[AC_nao]) / (exp(1.0 * ALGEBRAIC[AV_vfrt]) - 1.0);
    ALGEBRAIC[AV_PhiCaNa_ss] = 1.0 * ALGEBRAIC[AV_vffrt] * (ALGEBRAIC[AV_gamma_nass] * STATES[nasl] * exp(1.0 * ALGEBRAIC[AV_vfrt]) - CONSTANTS[AC_gamma_nao] * CONSTANTS[AC_nao]) / (exp(1.0 * ALGEBRAIC[AV_vfrt]) - 1.0);
    ALGEBRAIC[AV_ICab_junc] = CONSTANTS[AC_Fjunc] * CONSTANTS[AC_GCab] * (STATES[v] - ALGEBRAIC[AV_ECaj]);
    ALGEBRAIC[AV_ICab_sl] = CONSTANTS[AC_Fsl] * CONSTANTS[AC_GCab] * (STATES[v] - ALGEBRAIC[AV_ECasl]);
    ALGEBRAIC[AV_IKr_alpha] = 0.1161 * exp(0.299 * ALGEBRAIC[AV_vfrt]);
    ALGEBRAIC[AV_alpha_2] = 0.0578 * exp(0.971 * ALGEBRAIC[AV_vfrt]);
    ALGEBRAIC[AV_alpha_C2ToI] = 5.2e-05 * exp(1.525 * ALGEBRAIC[AV_vfrt]);
    ALGEBRAIC[AV_alpha_i] = 0.2533 * exp(0.5953 * ALGEBRAIC[AV_vfrt]);
    ALGEBRAIC[AV_IKr_beta] = 0.2442 * exp(-1.604 * ALGEBRAIC[AV_vfrt]);
    ALGEBRAIC[AV_beta_2] = 0.000349 * exp(-1.062 * ALGEBRAIC[AV_vfrt]);
    ALGEBRAIC[AV_beta_i] = 0.04568 * exp(-0.8209 * ALGEBRAIC[AV_vfrt]);
    ALGEBRAIC[AV_IKs_junc] = CONSTANTS[AC_Fjunc] * CONSTANTS[AC_GKs] * ALGEBRAIC[AV_gks_junc] * pow(STATES[xs_junc], 2.0) * (STATES[v] - ALGEBRAIC[AV_EKs]);
    ALGEBRAIC[AV_IKs_sl] = CONSTANTS[AC_Fsl] * CONSTANTS[AC_GKs] * ALGEBRAIC[AV_gks_sl] * pow(STATES[xs_sl], 2.0) * (STATES[v] - ALGEBRAIC[AV_EKs]);
    ALGEBRAIC[AV_INaj] = CONSTANTS[AC_Fjunc] * ALGEBRAIC[AV_INaBase] * (STATES[v] - ALGEBRAIC[AV_ENaj]);
    ALGEBRAIC[AV_INasl] = CONSTANTS[AC_Fsl] * ALGEBRAIC[AV_INaBase] * (STATES[v] - ALGEBRAIC[AV_ENasl]);
    ALGEBRAIC[AV_INaK_junc_PKA] = CONSTANTS[AC_Fjunc] * CONSTANTS[AC_IbarNaK] * ALGEBRAIC[AV_fnak] / (1.0 + pow(CONSTANTS[AC_KmNaip_PKA] / STATES[naj], 4.0)) * CONSTANTS[AC_ko] / (CONSTANTS[AC_ko] + CONSTANTS[AC_KmKo]);
    ALGEBRAIC[AV_INaK_junc_noPKA] = CONSTANTS[AC_Fjunc] * CONSTANTS[AC_IbarNaK] * ALGEBRAIC[AV_fnak] / (1.0 + pow(CONSTANTS[AC_KmNaip] / STATES[naj], 4.0)) * CONSTANTS[AC_ko] / (CONSTANTS[AC_ko] + CONSTANTS[AC_KmKo]);
    ALGEBRAIC[AV_INaK_sl_PKA] = CONSTANTS[AC_Fsl] * CONSTANTS[AC_IbarNaK] * ALGEBRAIC[AV_fnak] / (1.0 + pow(CONSTANTS[AC_KmNaip_PKA] / STATES[nasl], 4.0)) * CONSTANTS[AC_ko] / (CONSTANTS[AC_ko] + CONSTANTS[AC_KmKo]);
    ALGEBRAIC[AV_INaK_sl_noPKA] = CONSTANTS[AC_Fsl] * CONSTANTS[AC_IbarNaK] * ALGEBRAIC[AV_fnak] / (1.0 + pow(CONSTANTS[AC_KmNaip] / STATES[nasl], 4.0)) * CONSTANTS[AC_ko] / (CONSTANTS[AC_ko] + CONSTANTS[AC_KmKo]);
    ALGEBRAIC[AV_INaLj] = CONSTANTS[AC_Fjunc] * ALGEBRAIC[AV_GNaL] * (STATES[v] - ALGEBRAIC[AV_ENaj]) * STATES[mL] * ((1.0 - ALGEBRAIC[AV_fINaLp]) * STATES[hL] + ALGEBRAIC[AV_fINaLp] * STATES[hLp]);
    ALGEBRAIC[AV_INaLsl] = CONSTANTS[AC_Fsl] * ALGEBRAIC[AV_GNaL] * (STATES[v] - ALGEBRAIC[AV_ENasl]) * STATES[mL] * ((1.0 - ALGEBRAIC[AV_fINaLp]) * STATES[hL] + ALGEBRAIC[AV_fINaLp] * STATES[hLp]);
    ALGEBRAIC[AV_INabj] = CONSTANTS[AC_Fjunc] * CONSTANTS[AC_GNab] * (STATES[v] - ALGEBRAIC[AV_ENaj]);
    ALGEBRAIC[AV_INabsl] = CONSTANTS[AC_Fsl] * CONSTANTS[AC_GNab] * (STATES[v] - ALGEBRAIC[AV_ENasl]);
    ALGEBRAIC[AV_IpCa_junc] = CONSTANTS[AC_Fjunc] * pow(CONSTANTS[AC_Q10SLCaP], CONSTANTS[AC_Qpow]) * CONSTANTS[AC_IbarSLCaP] * pow(STATES[caj], 1.6) / (pow(CONSTANTS[AC_KmPCa], 1.6) + pow(STATES[caj], 1.6));
    ALGEBRAIC[AV_IpCa_sl] = CONSTANTS[AC_Fsl] * pow(CONSTANTS[AC_Q10SLCaP], CONSTANTS[AC_Qpow]) * CONSTANTS[AC_IbarSLCaP] * pow(STATES[casl], 1.6) / (pow(CONSTANTS[AC_KmPCa], 1.6) + pow(STATES[casl], 1.6));
    RATES[casr] = ALGEBRAIC[AV_Jserca] - (ALGEBRAIC[AV_J_SRleak] * CONSTANTS[AC_vmyo] / CONSTANTS[AC_vsr] + ALGEBRAIC[AV_J_SRCarel]) - ALGEBRAIC[AV_d_buffers_Csqn];
    RATES[cli] = ALGEBRAIC[AV_ICl_tot] * CONSTANTS[AC_Cmem] / (CONSTANTS[AC_vmyo] * CONSTANTS[AC_F]);
    ALGEBRAIC[AV_ICaK_i_BP] = CONSTANTS[AC_PCaK_P] * ALGEBRAIC[AV_PhiCaK_i] * STATES[d_P] * (ALGEBRAIC[AV_fBP] * (1.0 - STATES[nca_i]) + STATES[jca] * ALGEBRAIC[AV_fcaBP] * STATES[nca_i]);
    ALGEBRAIC[AV_ICaK_i_CaMK] = CONSTANTS[AC_PCaKp] * ALGEBRAIC[AV_PhiCaK_i] * STATES[d] * (ALGEBRAIC[AV_fp] * (1.0 - STATES[nca_i]) + STATES[jca] * ALGEBRAIC[AV_fcap] * STATES[nca_i]);
    ALGEBRAIC[AV_ICaK_i_NP] = CONSTANTS[AC_PCaK] * ALGEBRAIC[AV_PhiCaK_i] * STATES[d] * (ALGEBRAIC[AV_f] * (1.0 - STATES[nca_i]) + STATES[jca] * ALGEBRAIC[AV_fca] * STATES[nca_i]);
    ALGEBRAIC[AV_ICaK_i_PKA] = CONSTANTS[AC_PCaK_P] * ALGEBRAIC[AV_PhiCaK_i] * STATES[d_P] * (ALGEBRAIC[AV_f_P] * (1.0 - STATES[nca_i]) + STATES[jca] * ALGEBRAIC[AV_fcap_P] * STATES[nca_i]);
    ALGEBRAIC[AV_ICaK_ss_BP] = CONSTANTS[AC_PCaK_P] * ALGEBRAIC[AV_PhiCaK_ss] * STATES[d_P] * (ALGEBRAIC[AV_fBP] * (1.0 - STATES[nca]) + STATES[jca] * ALGEBRAIC[AV_fcaBP] * STATES[nca]);
    ALGEBRAIC[AV_ICaK_ss_CaMK] = CONSTANTS[AC_PCaKp] * ALGEBRAIC[AV_PhiCaK_ss] * STATES[d] * (ALGEBRAIC[AV_fp] * (1.0 - STATES[nca]) + STATES[jca] * ALGEBRAIC[AV_fcap] * STATES[nca]);
    ALGEBRAIC[AV_ICaK_ss_NP] = CONSTANTS[AC_PCaK] * ALGEBRAIC[AV_PhiCaK_ss] * STATES[d] * (ALGEBRAIC[AV_f] * (1.0 - STATES[nca]) + STATES[jca] * ALGEBRAIC[AV_fca] * STATES[nca]);
    ALGEBRAIC[AV_ICaK_ss_PKA] = CONSTANTS[AC_PCaK_P] * ALGEBRAIC[AV_PhiCaK_ss] * STATES[d_P] * (ALGEBRAIC[AV_f_P] * (1.0 - STATES[nca]) + STATES[jca] * ALGEBRAIC[AV_fcap_P] * STATES[nca]);
    ALGEBRAIC[AV_ICaL_i_BP] = CONSTANTS[AC_PCa_P] * ALGEBRAIC[AV_PhiCaL_i] * STATES[d_P] * (ALGEBRAIC[AV_fBP] * (1.0 - STATES[nca_i]) + STATES[jca] * ALGEBRAIC[AV_fcaBP] * STATES[nca_i]);
    ALGEBRAIC[AV_ICaL_i_CaMK] = CONSTANTS[AC_PCap] * ALGEBRAIC[AV_PhiCaL_i] * STATES[d] * (ALGEBRAIC[AV_fp] * (1.0 - STATES[nca_i]) + STATES[jca] * ALGEBRAIC[AV_fcap] * STATES[nca_i]);
    ALGEBRAIC[AV_ICaL_i_NP] = CONSTANTS[AC_PCa] * ALGEBRAIC[AV_PhiCaL_i] * STATES[d] * (ALGEBRAIC[AV_f] * (1.0 - STATES[nca_i]) + STATES[jca] * ALGEBRAIC[AV_fca] * STATES[nca_i]);
    ALGEBRAIC[AV_ICaL_i_PKA] = CONSTANTS[AC_PCa_P] * ALGEBRAIC[AV_PhiCaL_i] * STATES[d_P] * (ALGEBRAIC[AV_f_P] * (1.0 - STATES[nca_i]) + STATES[jca] * ALGEBRAIC[AV_fcap_P] * STATES[nca_i]);
    ALGEBRAIC[AV_ICaL_ss_BP] = CONSTANTS[AC_PCa_P] * ALGEBRAIC[AV_PhiCaL_ss] * STATES[d_P] * (ALGEBRAIC[AV_fBP] * (1.0 - STATES[nca]) + STATES[jca] * ALGEBRAIC[AV_fcaBP] * STATES[nca]);
    ALGEBRAIC[AV_ICaL_ss_CaMK] = CONSTANTS[AC_PCap] * ALGEBRAIC[AV_PhiCaL_ss] * STATES[d] * (ALGEBRAIC[AV_fp] * (1.0 - STATES[nca]) + STATES[jca] * ALGEBRAIC[AV_fcap] * STATES[nca]);
    ALGEBRAIC[AV_ICaL_ss_NP] = CONSTANTS[AC_PCa] * ALGEBRAIC[AV_PhiCaL_ss] * STATES[d] * (ALGEBRAIC[AV_f] * (1.0 - STATES[nca]) + STATES[jca] * ALGEBRAIC[AV_fca] * STATES[nca]);
    ALGEBRAIC[AV_ICaL_ss_PKA] = CONSTANTS[AC_PCa_P] * ALGEBRAIC[AV_PhiCaL_ss] * STATES[d_P] * (ALGEBRAIC[AV_f_P] * (1.0 - STATES[nca]) + STATES[jca] * ALGEBRAIC[AV_fcap_P] * STATES[nca]);
    ALGEBRAIC[AV_ICaNa_i_BP] = CONSTANTS[AC_PCaNa_P] * ALGEBRAIC[AV_PhiCaNa_i] * STATES[d_P] * (ALGEBRAIC[AV_fBP] * (1.0 - STATES[nca_i]) + STATES[jca] * ALGEBRAIC[AV_fcaBP] * STATES[nca_i]);
    ALGEBRAIC[AV_ICaNa_i_CaMK] = CONSTANTS[AC_PCaNap] * ALGEBRAIC[AV_PhiCaNa_i] * STATES[d] * (ALGEBRAIC[AV_fp] * (1.0 - STATES[nca_i]) + STATES[jca] * ALGEBRAIC[AV_fcap] * STATES[nca_i]);
    ALGEBRAIC[AV_ICaNa_i_NP] = CONSTANTS[AC_PCaNa] * ALGEBRAIC[AV_PhiCaNa_i] * STATES[d] * (ALGEBRAIC[AV_f] * (1.0 - STATES[nca_i]) + STATES[jca] * ALGEBRAIC[AV_fca] * STATES[nca_i]);
    ALGEBRAIC[AV_ICaNa_i_PKA] = CONSTANTS[AC_PCaNa_P] * ALGEBRAIC[AV_PhiCaNa_i] * STATES[d_P] * (ALGEBRAIC[AV_f_P] * (1.0 - STATES[nca_i]) + STATES[jca] * ALGEBRAIC[AV_fcap_P] * STATES[nca_i]);
    ALGEBRAIC[AV_ICaNa_ss_BP] = CONSTANTS[AC_PCaNa_P] * ALGEBRAIC[AV_PhiCaNa_ss] * STATES[d_P] * (ALGEBRAIC[AV_fBP] * (1.0 - STATES[nca]) + STATES[jca] * ALGEBRAIC[AV_fcaBP] * STATES[nca]);
    ALGEBRAIC[AV_ICaNa_ss_CaMK] = CONSTANTS[AC_PCaNap] * ALGEBRAIC[AV_PhiCaNa_ss] * STATES[d] * (ALGEBRAIC[AV_fp] * (1.0 - STATES[nca]) + STATES[jca] * ALGEBRAIC[AV_fcap] * STATES[nca]);
    ALGEBRAIC[AV_ICaNa_ss_NP] = CONSTANTS[AC_PCaNa] * ALGEBRAIC[AV_PhiCaNa_ss] * STATES[d] * (ALGEBRAIC[AV_f] * (1.0 - STATES[nca]) + STATES[jca] * ALGEBRAIC[AV_fca] * STATES[nca]);
    ALGEBRAIC[AV_ICaNa_ss_PKA] = CONSTANTS[AC_PCaNa_P] * ALGEBRAIC[AV_PhiCaNa_ss] * STATES[d_P] * (ALGEBRAIC[AV_f_P] * (1.0 - STATES[nca]) + STATES[jca] * ALGEBRAIC[AV_fcap_P] * STATES[nca]);
    ALGEBRAIC[AV_ICab_ICab] = ALGEBRAIC[AV_ICab_junc] + ALGEBRAIC[AV_ICab_sl];
    ALGEBRAIC[AV_beta_ItoC2] = ALGEBRAIC[AV_beta_2] * ALGEBRAIC[AV_beta_i] * ALGEBRAIC[AV_alpha_C2ToI] / (ALGEBRAIC[AV_alpha_2] * ALGEBRAIC[AV_alpha_i]);
    RATES[C0] = ALGEBRAIC[AV_IKr_beta] * STATES[C1] - ALGEBRAIC[AV_IKr_alpha] * STATES[C0];
    RATES[C1] = ALGEBRAIC[AV_IKr_alpha] * STATES[C0] + CONSTANTS[AC_IKr_beta_1] * STATES[C2] - (ALGEBRAIC[AV_IKr_beta] + CONSTANTS[AC_alpha_1]) * STATES[C1];
    RATES[O] = ALGEBRAIC[AV_alpha_2] * STATES[C2] + ALGEBRAIC[AV_beta_i] * STATES[I] - (ALGEBRAIC[AV_beta_2] + ALGEBRAIC[AV_alpha_i]) * STATES[O];
    ALGEBRAIC[AV_IKs_IKs] = ALGEBRAIC[AV_IKs_junc] + ALGEBRAIC[AV_IKs_sl];
    ALGEBRAIC[AV_INa_INa] = ALGEBRAIC[AV_INaj] + ALGEBRAIC[AV_INasl];
    ALGEBRAIC[AV_INaKj] = (1.0 - CONSTANTS[AC_fINaK_PKA]) * ALGEBRAIC[AV_INaK_junc_noPKA] + CONSTANTS[AC_fINaK_PKA] * ALGEBRAIC[AV_INaK_junc_PKA];
    ALGEBRAIC[AV_INaKsl] = (1.0 - CONSTANTS[AC_fINaK_PKA]) * ALGEBRAIC[AV_INaK_sl_noPKA] + CONSTANTS[AC_fINaK_PKA] * ALGEBRAIC[AV_INaK_sl_PKA];
    ALGEBRAIC[AV_INaL_INaL] = ALGEBRAIC[AV_INaLj] + ALGEBRAIC[AV_INaLsl];
    ALGEBRAIC[AV_INab_INab] = ALGEBRAIC[AV_INabj] + ALGEBRAIC[AV_INabsl];
    ALGEBRAIC[AV_IpCa_IpCa] = ALGEBRAIC[AV_IpCa_junc] + ALGEBRAIC[AV_IpCa_sl];
    ALGEBRAIC[AV_ICaK_junc] = ((1.0 - ALGEBRAIC[AV_fICaL_CaMKonly] - ALGEBRAIC[AV_fICaL_PKAonly] - ALGEBRAIC[AV_fICaL_BP]) * ALGEBRAIC[AV_ICaK_ss_NP] + ALGEBRAIC[AV_fICaL_CaMKonly] * ALGEBRAIC[AV_ICaK_ss_CaMK] + ALGEBRAIC[AV_fICaL_PKAonly] * ALGEBRAIC[AV_ICaK_ss_PKA] + ALGEBRAIC[AV_fICaL_BP] * ALGEBRAIC[AV_ICaK_ss_BP]) * STATES[ical_pureCDI_junc] * CONSTANTS[AC_ICaL_fractionSS];
    ALGEBRAIC[AV_ICaK_sl] = ((1.0 - ALGEBRAIC[AV_fICaL_CaMKonly] - ALGEBRAIC[AV_fICaL_PKAonly] - ALGEBRAIC[AV_fICaL_BP]) * ALGEBRAIC[AV_ICaK_i_NP] + ALGEBRAIC[AV_fICaL_CaMKonly] * ALGEBRAIC[AV_ICaK_i_CaMK] + ALGEBRAIC[AV_fICaL_PKAonly] * ALGEBRAIC[AV_ICaK_i_PKA] + ALGEBRAIC[AV_fICaL_BP] * ALGEBRAIC[AV_ICaK_i_BP]) * STATES[ical_pureCDI_sl] * (1.0 - CONSTANTS[AC_ICaL_fractionSS]);
    ALGEBRAIC[AV_ICaL_junc] = ((1.0 - ALGEBRAIC[AV_fICaL_CaMKonly] - ALGEBRAIC[AV_fICaL_PKAonly] - ALGEBRAIC[AV_fICaL_BP]) * ALGEBRAIC[AV_ICaL_ss_NP] + ALGEBRAIC[AV_fICaL_CaMKonly] * ALGEBRAIC[AV_ICaL_ss_CaMK] + ALGEBRAIC[AV_fICaL_PKAonly] * ALGEBRAIC[AV_ICaL_ss_PKA] + ALGEBRAIC[AV_fICaL_BP] * ALGEBRAIC[AV_ICaL_ss_BP]) * STATES[ical_pureCDI_junc] * CONSTANTS[AC_ICaL_fractionSS];
    ALGEBRAIC[AV_ICaL_sl] = ((1.0 - ALGEBRAIC[AV_fICaL_CaMKonly] - ALGEBRAIC[AV_fICaL_PKAonly] - ALGEBRAIC[AV_fICaL_BP]) * ALGEBRAIC[AV_ICaL_i_NP] + ALGEBRAIC[AV_fICaL_CaMKonly] * ALGEBRAIC[AV_ICaL_i_CaMK] + ALGEBRAIC[AV_fICaL_PKAonly] * ALGEBRAIC[AV_ICaL_i_PKA] + ALGEBRAIC[AV_fICaL_BP] * ALGEBRAIC[AV_ICaL_i_BP]) * STATES[ical_pureCDI_sl] * (1.0 - CONSTANTS[AC_ICaL_fractionSS]);
    ALGEBRAIC[AV_ICaNa_junc] = ((1.0 - ALGEBRAIC[AV_fICaL_CaMKonly] - ALGEBRAIC[AV_fICaL_PKAonly] - ALGEBRAIC[AV_fICaL_BP]) * ALGEBRAIC[AV_ICaNa_ss_NP] + ALGEBRAIC[AV_fICaL_CaMKonly] * ALGEBRAIC[AV_ICaNa_ss_CaMK] + ALGEBRAIC[AV_fICaL_PKAonly] * ALGEBRAIC[AV_ICaNa_ss_PKA] + ALGEBRAIC[AV_fICaL_BP] * ALGEBRAIC[AV_ICaNa_ss_BP]) * STATES[ical_pureCDI_junc] * CONSTANTS[AC_ICaL_fractionSS];
    ALGEBRAIC[AV_ICaNa_sl] = ((1.0 - ALGEBRAIC[AV_fICaL_CaMKonly] - ALGEBRAIC[AV_fICaL_PKAonly] - ALGEBRAIC[AV_fICaL_BP]) * ALGEBRAIC[AV_ICaNa_i_NP] + ALGEBRAIC[AV_fICaL_CaMKonly] * ALGEBRAIC[AV_ICaNa_i_CaMK] + ALGEBRAIC[AV_fICaL_PKAonly] * ALGEBRAIC[AV_ICaNa_i_PKA] + ALGEBRAIC[AV_fICaL_BP] * ALGEBRAIC[AV_ICaNa_i_BP]) * STATES[ical_pureCDI_sl] * (1.0 - CONSTANTS[AC_ICaL_fractionSS]);
    RATES[C2] = CONSTANTS[AC_alpha_1] * STATES[C1] + ALGEBRAIC[AV_beta_2] * STATES[O] + ALGEBRAIC[AV_beta_ItoC2] * STATES[I] - (CONSTANTS[AC_IKr_beta_1] + ALGEBRAIC[AV_alpha_2] + ALGEBRAIC[AV_alpha_C2ToI]) * STATES[C2];
    RATES[I] = ALGEBRAIC[AV_alpha_C2ToI] * STATES[C2] + ALGEBRAIC[AV_alpha_i] * STATES[O] - (ALGEBRAIC[AV_beta_ItoC2] + ALGEBRAIC[AV_beta_i]) * STATES[I];
    ALGEBRAIC[AV_INaK_INaK] = ALGEBRAIC[AV_INaKj] + ALGEBRAIC[AV_INaKsl];
    ALGEBRAIC[AV_ICaK] = ALGEBRAIC[AV_ICaK_junc] + ALGEBRAIC[AV_ICaK_sl];
    ALGEBRAIC[AV_ICaL_ICaL] = ALGEBRAIC[AV_ICaL_junc] + ALGEBRAIC[AV_ICaL_sl];
    ALGEBRAIC[AV_ICaNa] = ALGEBRAIC[AV_ICaNa_junc] + ALGEBRAIC[AV_ICaNa_sl];
    ALGEBRAIC[AV_ICa_tot_junc] = ALGEBRAIC[AV_ICaL_junc] + ALGEBRAIC[AV_ICab_junc] + ALGEBRAIC[AV_IpCa_junc] - 2.0 * ALGEBRAIC[AV_INaCa_ss];
    ALGEBRAIC[AV_ICa_tot_sl] = ALGEBRAIC[AV_ICaL_sl] + ALGEBRAIC[AV_ICab_sl] + ALGEBRAIC[AV_IpCa_sl] - 2.0 * ALGEBRAIC[AV_INaCa_i];
    ALGEBRAIC[AV_INa_tot_junc] = ALGEBRAIC[AV_INaj] + ALGEBRAIC[AV_INaLj] + ALGEBRAIC[AV_INabj] + 3.0 * ALGEBRAIC[AV_INaCa_ss] + 3.0 * ALGEBRAIC[AV_INaKj] + ALGEBRAIC[AV_ICaNa_junc];
    ALGEBRAIC[AV_INa_tot_sl] = ALGEBRAIC[AV_INasl] + ALGEBRAIC[AV_INaLsl] + ALGEBRAIC[AV_INabsl] + 3.0 * ALGEBRAIC[AV_INaCa_i] + 3.0 * ALGEBRAIC[AV_INaKsl] + ALGEBRAIC[AV_ICaNa_sl];
    ALGEBRAIC[AV_ICaL_junc_positive] = fabs(ALGEBRAIC[AV_ICaL_junc]);
    ALGEBRAIC[AV_ICa_tot] = ALGEBRAIC[AV_ICaL_ICaL] + ALGEBRAIC[AV_ICaK] + ALGEBRAIC[AV_ICaNa];
    RATES[caj] = -ALGEBRAIC[AV_ICa_tot_junc] * CONSTANTS[AC_Cmem] / (CONSTANTS[AC_vjunc] * 2.0 * CONSTANTS[AC_F]) + CONSTANTS[AC_J_ca_juncsl] / CONSTANTS[AC_vjunc] * (STATES[casl] - STATES[caj]) - ALGEBRAIC[AV_J_CaB_junction] + ALGEBRAIC[AV_J_SRCarel] * CONSTANTS[AC_vsr] / CONSTANTS[AC_vjunc] + ALGEBRAIC[AV_J_SRleak] * CONSTANTS[AC_vmyo] / CONSTANTS[AC_vjunc];
    RATES[casl] = -ALGEBRAIC[AV_ICa_tot_sl] * CONSTANTS[AC_Cmem] / (CONSTANTS[AC_vsl] * 2.0 * CONSTANTS[AC_F]) + CONSTANTS[AC_J_ca_juncsl] / CONSTANTS[AC_vsl] * (STATES[caj] - STATES[casl]) + CONSTANTS[AC_J_ca_slmyo] / CONSTANTS[AC_vsl] * (STATES[cai] - STATES[casl]) - ALGEBRAIC[AV_J_CaB_sl];
    RATES[naj] = -ALGEBRAIC[AV_INa_tot_junc] * CONSTANTS[AC_Cmem] / (CONSTANTS[AC_vjunc] * CONSTANTS[AC_F]) + CONSTANTS[AC_J_na_juncsl] / CONSTANTS[AC_vjunc] * (STATES[nasl] - STATES[naj]) - ALGEBRAIC[AV_d_buffers_NaBj];
    RATES[nasl] = -ALGEBRAIC[AV_INa_tot_sl] * CONSTANTS[AC_Cmem] / (CONSTANTS[AC_vsl] * CONSTANTS[AC_F]) + CONSTANTS[AC_J_na_juncsl] / CONSTANTS[AC_vsl] * (STATES[naj] - STATES[nasl]) + CONSTANTS[AC_J_na_slmyo] / CONSTANTS[AC_vsl] * (STATES[nai] - STATES[nasl]) - ALGEBRAIC[AV_d_buffers_NaBsl];
    ALGEBRAIC[AV_IK_tot] = ALGEBRAIC[AV_Ito_Ito] + ALGEBRAIC[AV_IKr_IKr] + ALGEBRAIC[AV_IKs_IKs] + ALGEBRAIC[AV_IK1_IK1] - 2.0 * ALGEBRAIC[AV_INaK_INaK] + ALGEBRAIC[AV_ICaK] + ALGEBRAIC[AV_IKb_IKb];
    ALGEBRAIC[AV_ICaL_junc_sigmoided] = 1.0 - 1.0 / (1.0 + pow(ALGEBRAIC[AV_ICaL_junc_positive] / 0.45, 4.5));
    RATES[ki] = -ALGEBRAIC[AV_IK_tot] * CONSTANTS[AC_Cmem] / (CONSTANTS[AC_vmyo] * CONSTANTS[AC_F]);
    ALGEBRAIC[Iion_cm] = ALGEBRAIC[AV_INa_tot_junc] + ALGEBRAIC[AV_INa_tot_sl] + ALGEBRAIC[AV_ICa_tot_junc] + ALGEBRAIC[AV_ICa_tot_sl] + ALGEBRAIC[AV_IK_tot] + ALGEBRAIC[AV_ICl_tot];
    ALGEBRAIC[AV_Jrel_inact_inf] = 1.0 / (1.0 + ALGEBRAIC[AV_ICaL_junc_sigmoided] / 0.001);
    ALGEBRAIC[AV_Jrel_inact_inf2] = 1.0 / (1.0 + ALGEBRAIC[AV_ICaL_junc_sigmoided] / 0.0006);
    ALGEBRAIC[AV_Jrel_inf] = CONSTANTS[AC_a_rel] * ALGEBRAIC[AV_ICaL_junc_sigmoided] / (1.0 + pow(CONSTANTS[AC_directRelMidpoint] / STATES[casr], 7.72672));
    RATES[v] = 0.0;
    RATES[jrel_icaldep_act] = (ALGEBRAIC[AV_Jrel_inf] - STATES[jrel_icaldep_act]) / ALGEBRAIC[AV_tau_rel];
    RATES[jrel_icaldep_f1] = (ALGEBRAIC[AV_Jrel_inact_inf] - STATES[jrel_icaldep_f1]) / CONSTANTS[AC_tauInact];
    RATES[jrel_icaldep_f2] = (ALGEBRAIC[AV_Jrel_inact_inf2] - STATES[jrel_icaldep_f2]) / CONSTANTS[AC_tauInact2];
    

    

    ALGEBRAIC[Istim] = computeIstim(VOI, stimulus);

    if (solveVmWithinODESolver)
    {
        RATES[v] = -ALGEBRAIC[Iion_cm] - ALGEBRAIC[Istim];
    }

}
