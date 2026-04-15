"""
Ionic Model Catalog

A static, build-time-generated catalog of all ionic and active tension models.
This module exposes the exact variables each model supports so an autonomous
agent can plan outputVariables without running the solver.

All variable names are extracted from C++ source files and are guaranteed to be
exact.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Final


@dataclass(frozen=True)
class IonicModelEntry:
    """Metadata for a single ionic model."""

    states: tuple[str, ...]
    """State variables from stateVariableNames()."""

    algebraic: tuple[str, ...]
    """Algebraic variables from algebraicVariableNames()."""

    constants: tuple[str, ...]
    """Constant parameters from constantVariableNames()."""

    recommended_exports: tuple[str, ...]
    """Minimum useful export set (subset of states)."""

    compatible_tissues: tuple[str, ...]
    """Which tissue types this model supports (e.g. epicardialCells, myocyte)."""

    compatible_solvers: tuple[str, ...]
    """Which myocardium solvers can host this model."""

    species: tuple[str, ...]
    """Species the model is designed for (e.g. ('human',), ('pig',), ('generic',))."""

    cardiac_region: tuple[str, ...]
    """Cardiac region the model is designed for (e.g. ('ventricle',), ('atrium',), ('sinoatrial_node',), ('purkinje',), ('manufactured',))."""

    model_type: str
    """Model classification: 'phenomenological', 'ionic', or 'manufactured'."""

    description: str
    """Human-readable description."""

    notes: str = ""
    """Additional notes or warnings."""


@dataclass(frozen=True)
class ActiveTensionModelEntry:
    """Metadata for a single active tension model."""

    states: tuple[str, ...]
    """State variables."""

    algebraic: tuple[str, ...]
    """Algebraic variables."""

    constants: tuple[str, ...]
    """Constant parameters."""

    rates: tuple[str, ...]
    """Rate variables (derivatives of states)."""

    recommended_exports: tuple[str, ...]
    """Minimum useful export set."""

    description: str
    """Human-readable description."""

    notes: str = ""
    """Additional notes or warnings."""


SOLVER_COMPATIBILITY_RULES: Final[tuple[dict, ...]] = (
    {
        "myocardium_solver": "monodomainSolver",
        "purkinje_solver": "monodomain1DSolver",
        "required_coupler": "reactionDiffusionPvjCoupler",
        "valid": True,
    },
    {
        "myocardium_solver": "eikonalSolver",
        "purkinje_solver": "eikonalSolver",
        "required_coupler": "eikonalPvjCoupler",
        "valid": True,
    },
    {
        "myocardium_solver": "monodomainSolver",
        "purkinje_solver": "eikonalSolver",
        "required_coupler": None,
        "valid": False,
        "reason": "Incompatible physics: reaction-diffusion myocardium cannot couple to eikonal Purkinje",
    },
    {
        "myocardium_solver": "eikonalSolver",
        "purkinje_solver": "monodomain1DSolver",
        "required_coupler": None,
        "valid": False,
        "reason": "Incompatible physics: eikonal myocardium cannot couple to reaction-diffusion Purkinje",
    },
    {
        "myocardium_solver": "bidomainSolver",
        "purkinje_solver": "*",
        "required_coupler": None,
        "valid": False,
        "reason": "bidomainSolver does not support Purkinje network coupling",
    },
    {
        "myocardium_solver": "singleCellSolver",
        "purkinje_solver": "*",
        "required_coupler": None,
        "valid": False,
        "reason": "singleCellSolver has no PDE domain; Purkinje coupling not applicable",
    },
)


IONIC_MODEL_CATALOG: Final[dict[str, IonicModelEntry]] = {
    "AlievPanfilov": IonicModelEntry(
        states=("u", "recovery_r"),
        algebraic=("AV_eps", "Istim", "Iion_cm"),
        constants=("AC_Vp", "AC_Vr", "AC_Vth", "AC_epsilon", "AC_k", "AC_mu1", "AC_mu2", "AC_a"),
        recommended_exports=("u", "recovery_r"),
        compatible_tissues=("myocyte",),
        compatible_solvers=("monodomainSolver", "bidomainSolver", "singleCellSolver"),
        species=("generic",),
        cardiac_region=("ventricle",),
        model_type="phenomenological",
        description="Simplified 2-variable phenomenological model; not species-specific (Aliev & Panfilov 1996).",
        notes="Phenomenological; works in any PDE-based or ODE-only solver, same as BuenoOrovio.",
    ),
    "BuenoOrovio": IonicModelEntry(
        states=("u", "v", "w", "s"),
        algebraic=("Jfi", "Jso", "Jsi", "Jion", "Istim", "tauSo", "tauO", "tauVMinus", "tauWMinus", "tauS", "vInfty", "wInfty"),
        constants=("uO", "uU", "thetaV", "thetaW", "thetaVMinus", "thetaO", "tauV1Minus", "tauV2Minus", "tauVPlus", "tauW1Minus", "tauW2Minus", "kWMinus", "uWMinus", "tauWPlus", "tauFi", "tauO1", "tauO2", "tauSo1", "tauSo2", "kSo", "uSo", "tauS1", "tauS2", "kS", "uS", "tauSi", "tauWInfty", "wInftyStar"),
        recommended_exports=("u", "v", "w", "s"),
        compatible_tissues=("epicardialCells", "mCells", "endocardialCells"),
        compatible_solvers=("monodomainSolver", "bidomainSolver", "singleCellSolver"),
        species=("generic",),
        cardiac_region=("ventricle",),
        model_type="phenomenological",
        description="Minimal phenomenological model designed to replicate human/mammalian ventricular action potentials (Bueno-Orovio et al. 2008).",
    ),
    "Courtemanche": IonicModelEntry(
        states=("membrane_V", "sodium_Nai", "potassium_Ki", "calcium_Cai", "calcium_CaUp", "calcium_CaRel", "ina_m", "ina_h", "ina_j", "ito_oa", "ito_oi", "ikur_ua", "ikur_ui", "ikr_xr", "iks_xs", "ical_d", "ical_f", "ical_fCa", "cajsr_u", "cajsr_v", "cajsr_w"),
        algebraic=("AV_ICaL", "AV_IK1", "AV_IKr", "AV_IKs", "AV_IKur", "AV_INa", "AV_INaCa", "AV_INaK", "AV_IbCa", "AV_IbNa", "AV_IpCa", "AV_Ito", "AV_ina_m_inf", "AV_ina_m_tau", "AV_ina_h_inf", "AV_ina_h_tau", "AV_ina_j_inf", "AV_ina_j_tau", "AV_ical_d_inf", "AV_ical_d_tau", "AV_ical_f_inf", "AV_ical_f_tau", "AV_ito_oa_inf", "AV_ito_oa_tau", "AV_ito_oi_inf", "AV_ito_oi_tau", "AV_ikur_ua_inf", "AV_ikur_ua_tau", "AV_ikur_ui_inf", "AV_ikur_ui_tau", "AV_gKur", "AV_ikr_xr_inf", "AV_ikr_xr_tau", "AV_iks_xs_inf", "AV_iks_xs_tau", "AV_fNaK", "AV_cajsr_w_inf", "AV_cajsr_w_tau", "Iion_cm", "Istim"),
        constants=("AC_CMDN_max", "AC_CSQN_max", "AC_Ca_up_max", "AC_Cao", "AC_Cm", "AC_ECaL", "AC_F", "AC_FRT", "AC_INaCa_max", "AC_INaK_max", "AC_I_diff", "AC_I_up_max", "AC_IpCa_max", "AC_KQ10", "AC_K_rel", "AC_K_up", "AC_KmCa", "AC_KmKo", "AC_KmNa", "AC_KmNai", "AC_Km_CMDN", "AC_Km_CSQN", "AC_Km_TRPN", "AC_Ko", "AC_Nao", "AC_R", "AC_RTF", "AC_T", "AC_TRPN_max", "AC_V_cell", "AC_V_i", "AC_V_rel", "AC_V_up", "AC_c1", "AC_c2", "AC_cajsr_u_tau", "AC_g", "AC_gCaL", "AC_gK1", "AC_gKr", "AC_gKs", "AC_gKur_base", "AC_gNa", "AC_gbCa", "AC_gbNa", "AC_gto", "AC_ical_fCa_tau", "AC_ksat", "AC_sigma", "AC_tau_tr"),
        recommended_exports=("membrane_V", "calcium_Cai"),
        compatible_tissues=("myocyte",),
        compatible_solvers=("monodomainSolver", "bidomainSolver", "singleCellSolver"),
        species=("human",),
        cardiac_region=("atrium",),
        model_type="ionic",
        description="Gold-standard ionic model for human atrial simulations (Courtemanche et al. 1998).",
    ),
    "Fabbri": IonicModelEntry(
        states=("membrane_V", "Na_i", "If_y_gate_y", "INa_m_gate_m", "INa_h_gate_h", "ICaL_dL_gate_dL", "ICaL_fL_gate_fL", "ICaL_fCa_gate_fCa", "ICaT_dT_gate_dT", "ICaT_fT_gate_fT", "SR_R", "SR_O", "SR_I", "SR_RI", "Buffer_fTMM", "Buffer_fCMi", "Buffer_fCMs", "Buffer_fTC", "Buffer_fTMC", "Buffer_fCQ", "Ca_i", "Ca_nsr", "Ca_jsr", "Ca_sub", "IKur_rKur_gate_r_Kur", "IKur_sKur_gate_s_Kur", "Ito_q_gate_q", "Ito_r_gate_r", "IKr_pa_gate_paS", "IKr_pa_gate_paF", "IKr_pi_gate_piy", "IKs_n_gate_n", "IKACh_a_gate_a"),
        algebraic=("AV_P_tot", "AV_diff", "AV_j_SRCarel", "AV_kCaSR", "AV_kiSRCa", "AV_koSRCa", "AV_delta_fCMi", "AV_delta_fCMs", "AV_delta_fCQ", "AV_delta_fTC", "AV_delta_fTMC", "AV_delta_fTMM", "AV_fCa_infinity", "AV_tau_fCa", "AV_j_Ca_dif", "AV_j_tr", "AV_j_up", "AV_V_clamp", "AV_V", "AV_Nai", "AV_i_siCa", "AV_i_siK", "AV_i_siNa", "AV_i_CaL_i_CaL", "AV_adVm", "AV_bdVm", "AV_alpha_dL", "AV_beta_dL", "AV_dL_infinity", "AV_tau_dL", "AV_tau_fL", "AV_fL_infinity", "AV_i_CaT_i_CaT", "AV_dT_infinity", "AV_tau_dT", "AV_fT_infinity", "AV_tau_fT", "AV_beta_a", "AV_a_infinity", "AV_tau_a", "AV_alfapaF", "AV_betapaF", "AV_pa_infinity", "AV_tau_paF", "AV_tau_paS", "AV_pi_infinity", "AV_tau_pi", "AV_E_Ks", "AV_i_Ks_i_Ks", "AV_alpha_n", "AV_beta_n", "AV_n_infinity", "AV_tau_n", "AV_r_Kur_infinity", "AV_tau_r_Kur", "AV_s_Kur_infinity", "AV_tau_s_Kur", "AV_E_mh", "AV_i_Na_", "AV_i_Na_L", "AV_i_Na_i_Na", "AV_di", "AV_i_NaCa_do", "AV_k32", "AV_k41", "AV_k43", "AV_k12", "AV_k14", "AV_k21", "AV_k23", "AV_x1", "AV_x2", "AV_x3", "AV_x4", "AV_i_NaCa_i_NaCa", "AV_alpha_h", "AV_beta_h", "AV_h_infinity", "AV_tau_h", "AV_E0_m", "AV_beta_m", "AV_m_infinity", "AV_alpha_m", "AV_tau_m", "AV_tau_y", "AV_y_infinity", "AV_q_infinity", "AV_tau_q", "AV_r_infinity", "AV_tau_r", "AV_E_Ca", "AV_E_Na", "AV_i_KACh_i_KACh", "AV_i_Kr_i_Kr", "AV_i_Kur_i_Kur", "AV_i_NaK_i_NaK", "AV_i_fK", "AV_i_fNa", "AV_i_to_i_to", "AV_i_f_i_f", "Istim", "Iion_cm"),
        constants=("AC_EC50_SR", "AC_HSR", "AC_MaxSR", "AC_MinSR", "AC_kiCa", "AC_kim", "AC_koCa", "AC_kom", "AC_ks", "AC_CM_tot", "AC_CQ_tot", "AC_Mgi", "AC_TC_tot", "AC_TMC_tot", "AC_kb_CM", "AC_kb_CQ", "AC_kb_TC", "AC_kb_TMC", "AC_kb_TMM", "AC_kf_CM", "AC_kf_CQ", "AC_kf_TC", "AC_kf_TMC", "AC_kf_TMM", "AC_L_cell", "AC_L_sub", "AC_R_cell", "AC_V_i_part", "AC_V_jsr_part", "AC_V_nsr_part", "AC_V_cell", "AC_V_sub", "AC_V_i", "AC_V_jsr", "AC_V_nsr", "AC_ACh", "AC_Iso_1_uM", "AC_Km_fCa", "AC_alpha_fCa", "AC_K_up", "AC_P_up_basal", "AC_b_up", "AC_slope_up", "AC_tau_dif_Ca", "AC_tau_tr", "AC_P_up", "AC_V_holding", "AC_V_test", "AC_t_holding", "AC_t_test", "AC_Cao", "AC_Ki", "AC_Ko", "AC_Nao", "AC_C", "AC_F", "AC_Membrane_R", "AC_T", "AC_clamp_mode", "AC_RTONF", "AC_Nai_clamp", "AC_ACh_block", "AC_i_CaL_Iso_increase", "AC_P_CaL", "AC_Iso_shift_dL", "AC_Iso_slope_dL", "AC_V_dL", "AC_k_dL", "AC_k_fL", "AC_shift_fL", "AC_P_CaT", "AC_offset_fT", "AC_ACh_on", "AC_g_KACh", "AC_alpha_a", "AC_g_Kr", "AC_g_Ks_", "AC_g_Ks", "AC_i_Ks_n_gate_Iso_shift", "AC_g_Kur", "AC_g_Na", "AC_g_Na_L", "AC_K1ni", "AC_K1no", "AC_K2ni", "AC_K2no", "AC_K3ni", "AC_K3no", "AC_K_NaCa", "AC_Kci", "AC_Kcni", "AC_Kco", "AC_Qci", "AC_Qco", "AC_Qn", "AC_blockade_NaCa", "AC_k34", "AC_i_NaK_Iso_increase", "AC_Km_Kp", "AC_Km_Nap", "AC_i_NaK_max", "AC_delta_m", "AC_Km_f", "AC_alpha", "AC_blockade", "AC_g_f", "AC_G_f", "AC_G_f_K", "AC_G_f_Na", "AC_g_f_K", "AC_g_f_Na", "AC_ACh_shift", "AC_i_f_y_gate_Iso_shift", "AC_y_shift", "AC_g_to", "AC_E_K"),
        recommended_exports=("membrane_V", "Ca_i"),
        compatible_tissues=("myocyte",),
        compatible_solvers=("monodomainSolver", "bidomainSolver", "singleCellSolver"),
        species=("human",),
        cardiac_region=("sinoatrial_node",),
        model_type="ionic",
        description="Ionic model specifically for human sinoatrial node (pacemaker) cells (Fabbri et al. 2017).",
    ),
    "Gaur": IonicModelEntry(
        states=("cell_v", "nai", "nass", "ki", "kss", "cai", "cai2", "cass", "cansr", "cajsr", "cacsr", "I_Na_m", "I_Na_h", "I_Na_j", "INaL_ml", "INaL_hl", "ICaL_d", "ICaL_fca", "IKr_xr", "IKs_xs1", "IKs_xs2", "ITo_aa", "CICR_Jrel2", "CICR_Jrel1", "CaMK_CaMKt", "CICR_tjsrol", "CICR_A", "ICaL_fs", "ICaL_ff"),
        algebraic=("AV_Jdiff", "AV_JdiffK", "AV_JdiffNa", "AV_EK", "AV_ENa", "AV_vffrt", "AV_vfrt", "AV_CaMKb", "AV_EKs", "AV_CaMKa", "AV_diff_CaMKt", "AV_CaMK_f", "AV_aa_h", "AV_aa_j", "AV_aa_m", "AV_alpha_h", "AV_alpha_j", "AV_alpha_m", "AV_beta_h", "AV_beta_j", "AV_beta_m", "AV_h_inf", "AV_j_inf", "AV_m_inf", "AV_tau_h", "AV_tau_j", "AV_tau_m", "AV_aa", "AV_i_Na", "AV_i_NaL", "AV_d_inf", "AV_d_tau", "AV_f_Ca_inf", "AV_f_Ca_tau", "AV_fp_Ca_inf", "AV_fp_Ca_tau", "AV_f_inf", "AV_f_tau", "AV_fss_inf", "AV_fss_tau", "AV_fs_inf", "AV_fs_tau", "AV_i_CaL", "AV_i_Ca_Lp_Ca", "AV_i_Ca_L_Ca", "AV_i_CaNa_L", "AV_i_CaK_L", "AV_i_pCa", "AV_k1", "AV_k2", "AV_k2prime", "AV_k3", "AV_k4", "AV_k5", "AV_k6", "AV_k7", "AV_x1to2", "AV_x2to3", "AV_x3to4", "AV_x4to5", "AV_x5to6", "AV_x6to1", "AV_x6to7", "AV_x7to6", "AV_E_Ca_CaMK", "AV_fracLCaCaMK", "AV_i_Ca_L_CaMK", "AV_i_Ca_LCaMK_Ca", "AV_i_CaNCX", "AV_i_NaCa", "AV_i_NaCa_i", "AV_i_NaCa_ss", "AV_i_pNaK", "AV_i_Ks", "AV_xs1_inf", "AV_xs1_tau", "AV_xs2_inf", "AV_xs2_tau", "AV_xr_inf", "AV_xr_tau", "AV_i_Kr", "AV_a_inf", "AV_a_tau", "AV_i_to", "AV_iF_Ca", "AV_i_f", "AV_i_fNa", "AV_i_fK", "AV_i_K1", "AV_i_Kb", "AV_E_K", "AV_Istim", "Iion_cm"),
        constants=("AC_nao", "AC_cao", "AC_ko", "AC_Cm", "AC_R", "AC_T", "AC_F", "AC_zNa", "AC_zK", "AC_zCa", "AC_g_Na", "AC_g_NaL", "AC_P_Ca_L", "AC_k_NaCa", "AC_K_NaCa_3Na", "AC_K_NaCa_2Ca", "AC_K_NaCa_Km_Nai", "AC_K_NaCa_Km_Ca", "AC_h_NaCa", "AC_P_NaK", "AC_K_NaK_K", "AC_K_NaK_Na", "AC_g_Ks", "AC_g_Kr", "AC_g_to", "AC_g_f", "AC_V_nsr", "AC_V_ss", "AC_V_myoplasm", "AC_V_jsr", "AC_L_cell", "AC_L_sub", "AC_R_cell", "AC_P_rel_max", "AC_k_rel_inf", "AC_k_rel_tau", "AC_Max_SR", "AC_Min_SR", "AC_ec_50_SR", "AC_HSR", "AC_k_jup", "AC_P_up_max", "AC_K_up", "AC_n_CaMK", "AC_K_m_CaMK", "AC_CaMK_0", "AC_b_CaMK", "AC_B_Nai", "AC_B_Cai", "AC_B_Cass", "AC_K_B_Nai", "AC_K_B_Cai", "AC_K_B_Cass", "AC_alpha_Nass", "AC_beta_Nass", "AC_alpha_K", "AC_beta_K", "AC_alpha_Ca", "AC_beta_Ca", "AC_g_K1", "AC_p_f_i", "AC_r_i", "AC_s_i", "AC_V_half_If", "AC_k_If", "AC_g_b_Na", "AC_g_b_Ca"),
        recommended_exports=("cell_v", "cai"),
        compatible_tissues=("myocyte",),
        compatible_solvers=("monodomainSolver", "bidomainSolver", "singleCellSolver"),
        species=("pig",),
        cardiac_region=("ventricle",),
        model_type="ionic",
        description="Pig ventricular ionic model based on the Gaur-Rudy-Luo formulation (1996).",
    ),
    "Grandi": IonicModelEntry(
        states=("V", "m", "hf", "hs", "j", "xrf", "xrs", "d", "ff", "fs", "fcaf", "fcas", "jca", "nca", "ffp", "fcafp", "fcasp", "xrsp", "xs1", "xs2", "y", "oa", "oi", "r", "u", "Nai", "Cass", "Cajsr", "Cansr"),
        algebraic=("AV_INa", "AV_INaL", "AV_ICaL", "AV_ICaNa", "AV_ICaK", "AV_IpCa", "AV_INaCa", "AV_INaK", "AV_IKr", "AV_IKs", "AV_IK1", "AV_IKb", "AV_INab", "AV_ICab", "AV_Ist", "AV_Iion", "AV_jup", "AV_jtr", "AV_jrel", "AV_Istim", "Iion_cm"),
        constants=("AC_Vc", "AC_Vsr", "AC_Vjsr", "AC_Vss", "AC_ko", "AC_nao", "AC_cao", "AC_F", "AC_R", "AC_T", "AC_zNa", "AC_zCa", "AC_zK", "AC_g_Na", "AC_g_NaL", "AC_pCa_L", "AC_g_K1", "AC_g_Kr", "AC_g_Ks", "AC_g_Kp", "AC_g_Kur", "AC_g_to", "AC_g_f", "AC_g_b_Na", "AC_g_b_Ca", "AC_pNaK", "AC_KmKo", "AC_KmNai", "AC_p_NaCa", "AC_KmCai", "AC_KmCao", "AC_KmNai_ncx", "AC_KmNao", "AC_ksat", "AC_nu", "AC_p_Ca_L", "AC_P_up", "AC_K_up", "AC_n_CICR", "AC_K_rel", "AC_CaMK_0", "AC_b_CaMK", "AC_K_m_CaMK", "AC_n", "AC_p_Ca_p", "AC_K_pCa", "AC_tau_m", "AC_m_ss", "AC_tau_hf", "AC_tau_hs", "AC_h_ss", "AC_tau_j", "AC_j_ss", "AC_tau_xrf", "AC_tau_xrs", "AC_x_ss", "AC_tau_d", "AC_d_ss", "AC_tau_ff", "AC_tau_fs", "AC_f_ss", "AC_tau_fcaf", "AC_tau_fcas", "AC_fc_ss", "AC_tau_jca", "AC_nca_ss", "AC_tau_nca", "AC_K_Nai", "AC_B_Nai", "AC_B_Cass", "AC_K_Cass", "AC_tau_cajsr", "AC_cajsr_ss", "AC_tau_cansr", "AC_cansr_ss"),
        recommended_exports=("V", "Cass"),
        compatible_tissues=("myocyte",),
        compatible_solvers=("monodomainSolver", "bidomainSolver", "singleCellSolver"),
        species=("human",),
        cardiac_region=("ventricle",),
        model_type="ionic",
        description="Human ventricular ionic model with detailed calcium signalling and electrolytes (Grandi et al. 2010).",
    ),
    "ORd": IonicModelEntry(
        states=("V", "m", "hf", "hs", "j", "hsp", "jp", "fLCa", "mL", "hL", "hLp", "a", "iF", "iS", "ap", "iFp", "iSp", "d", "ff", "fs", "fcaf", "fcas", "jca", "nca", "ffp", "fcafp", "fcasp", "xrf", "xrs", "xs1", "xs2", "xk1", "y", "oa", "oi", "r", "u", "Nai", "Cass", "Cajsr", "Cansr", "CaMKt"),
        algebraic=("AV_INa", "AV_INaL", "AV_ICaL", "AV_ICaNa", "AV_ICaK", "AV_IpCa", "AV_INaCa", "AV_INaK", "AV_IKr", "AV_IKs", "AV_IK1", "AV_IKb", "AV_INab", "AV_ICab", "AV_Ist", "AV_Iion", "AV_jup", "AV_jtr", "AV_jrel", "AV_Istim", "Iion_cm"),
        constants=("AC_nao", "AC_cao", "AC_ko", "AC_zNa", "AC_zCa", "AC_zK", "AC_R", "AC_T", "AC_F", "AC_g_Na", "AC_g_NaL", "AC_pCa_L", "AC_g_K1", "AC_g_Kr", "AC_g_Ks", "AC_g_Kp", "AC_g_Kur", "AC_g_to", "AC_g_f", "AC_g_b_Na", "AC_g_b_Ca", "AC_pNaK", "AC_KmKo", "AC_KmNai", "AC_p_NaCa", "AC_KmCai", "AC_KmCao", "AC_KmNai_ncx", "AC_KmNao", "AC_ksat", "AC_nu", "AC_p_Ca_L", "AC_P_up", "AC_K_up", "AC_n_CICR", "AC_K_rel", "AC_CaMK_0", "AC_b_CaMK", "AC_K_m_CaMK", "AC_n", "AC_p_Ca_p", "AC_K_pCa", "AC_Vc", "AC_Vsr", "AC_Vjsr", "AC_Vss"),
        recommended_exports=("V", "Cass"),
        compatible_tissues=("myocyte",),
        compatible_solvers=("monodomainSolver", "bidomainSolver", "singleCellSolver"),
        species=("human",),
        cardiac_region=("ventricle",),
        model_type="ionic",
        description="O'Hara-Rudy; industry-standard human ventricular ionic model (O'Hara et al. 2011).",
    ),
    "Stewart": IonicModelEntry(
        states=("V", "m", "h1", "h2", "j", "d", "f", "f_ca", "r", "s", "xK", "Nai", "Ki", "Cai", "CaRel", "CaUp"),
        algebraic=("AV_INa", "AV_ICa", "AV_IbCa", "AV_IbNa", "AV_IK1", "AV_IKp", "AV_IKr", "AV_IKs", "AV_Ito", "AV_IpCa", "AV_INaK", "AV_INaCa", "AV_Iion", "AV_J_up", "AV_J_tr", "AV_J_rel", "AV_J_xfer", "AV_J_leak", "AV_Istim", "Iion_cm"),
        constants=("AC_g_Na", "AC_g_K1", "AC_g_Kr", "AC_g_Ks", "AC_g_b_Na", "AC_g_b_Ca", "AC_g_to", "AC_g_Kp", "AC_g_Ca", "AC_pNaK", "AC_pCa", "AC_pNaCa", "AC_kNaCa", "AC_KmCai", "AC_KmNai", "AC_KmCao", "AC_KmNao", "AC_ksat", "AC_zNa", "AC_zCa", "AC_nH", "AC_Vc", "AC_VJSR", "AC_Vup", "AC_Vrel", "AC_Vss", "AC_F", "AC_R", "AC_T", "AC_zK", "AC_kao", "AC_kio", "AC_nao", "AC_cao", "AC_Ko", "AC_Nao", "AC_Cao", "AC_CaMK0", "AC_bCaMK", "AC_KmCaMK"),
        recommended_exports=("V", "Cai"),
        compatible_tissues=("myocyte",),
        compatible_solvers=("monodomainSolver", "bidomainSolver", "singleCellSolver"),
        species=("human",),
        cardiac_region=("purkinje",),
        model_type="ionic",
        description="Human Purkinje fibre ionic model for specialised cardiac conduction (Stewart et al. 2009).",
    ),
    "TNNP": IonicModelEntry(
        states=("membrane_V", "sodium_Nai", "potassium_Ki", "calcium_Cai", "calcium_CaUp", "calcium_CaRel", "ina_m", "ina_h", "ina_j", "ito_oa", "ito_oi", "ikur_ua", "ikur_ui", "ikr_xr", "iks_xs", "ical_d", "ical_f", "ical_fCa", "cajsr_u", "cajsr_v", "cajsr_w"),
        algebraic=("AV_ICaL", "AV_IK1", "AV_IKr", "AV_IKs", "AV_IKur", "AV_INa", "AV_INaCa", "AV_INaK", "AV_IbCa", "AV_IbNa", "AV_IpCa", "AV_Ito", "AV_ina_m_inf", "AV_ina_m_tau", "AV_ina_h_inf", "AV_ina_h_tau", "AV_ina_j_inf", "AV_ina_j_tau", "AV_ical_d_inf", "AV_ical_d_tau", "AV_ical_f_inf", "AV_ical_f_tau", "AV_ito_oa_inf", "AV_ito_oa_tau", "AV_ito_oi_inf", "AV_ito_oi_tau", "AV_ikur_ua_inf", "AV_ikur_ua_tau", "AV_ikur_ui_inf", "AV_ikur_ui_tau", "AV_gKur", "AV_ikr_xr_inf", "AV_ikr_xr_tau", "AV_iks_xs_inf", "AV_iks_xs_tau", "AV_fNaK", "AV_cajsr_w_inf", "AV_cajsr_w_tau", "Iion_cm", "Istim"),
        constants=("AC_CMDN_max", "AC_CSQN_max", "AC_Ca_up_max", "AC_Cao", "AC_Cm", "AC_ECaL", "AC_F", "AC_FRT", "AC_INaCa_max", "AC_INaK_max", "AC_I_diff", "AC_I_up_max", "AC_IpCa_max", "AC_KQ10", "AC_K_rel", "AC_K_up", "AC_KmCa", "AC_KmKo", "AC_KmNa", "AC_KmNai", "AC_Km_CMDN", "AC_Km_CSQN", "AC_Km_TRPN", "AC_Ko", "AC_Nao", "AC_R", "AC_RTF", "AC_T", "AC_TRPN_max", "AC_V_cell", "AC_V_i", "AC_V_rel", "AC_V_up", "AC_c1", "AC_c2", "AC_cajsr_u_tau", "AC_g", "AC_gCaL", "AC_gK1", "AC_gKr", "AC_gKs", "AC_gKur_base", "AC_gNa", "AC_gbCa", "AC_gbNa", "AC_gto", "AC_ical_fCa_tau", "AC_ksat", "AC_sigma", "AC_tau_tr"),
        recommended_exports=("membrane_V", "calcium_Cai"),
        compatible_tissues=("epicardialCells", "mCells", "endocardialCells"),
        compatible_solvers=("monodomainSolver", "bidomainSolver", "singleCellSolver"),
        species=("human",),
        cardiac_region=("ventricle",),
        model_type="ionic",
        description="Ten Tusscher-Noble-Noble-Panfilov; widely used human ventricular ionic model for 3-D simulations (ten Tusscher & Panfilov 2006).",
    ),
    "ToRORd_dynCl": IonicModelEntry(
        states=("V", "m", "hf", "hs", "j", "hsp", "jp", "fLCa", "mL", "hL", "hLp", "a", "iF", "iS", "ap", "iFp", "iSp", "d", "ff", "fs", "fcaf", "fcas", "jca", "nca", "ffp", "fcafp", "fcasp", "xrf", "xrs", "xs1", "xs2", "xk1", "y", "oa", "oi", "r", "u", "Nai", "Cass", "Cajsr", "Cansr", "CaMKt", "Cli", "Clss"),
        algebraic=("AV_INa", "AV_INaL", "AV_ICaL", "AV_ICaNa", "AV_ICaK", "AV_IpCa", "AV_INaCa", "AV_INaK", "AV_IKr", "AV_IKs", "AV_IK1", "AV_IKb", "AV_IKCl", "AV_INab", "AV_ICab", "AV_IClb", "AV_Ist", "AV_Iion", "AV_jup", "AV_jtr", "AV_jrel", "AV_Istim", "Iion_cm"),
        constants=("AC_nao", "AC_cao", "AC_ko", "AC_clo", "AC_clss", "AC_zNa", "AC_zCa", "AC_zK", "AC_R", "AC_T", "AC_F", "AC_g_Na", "AC_g_NaL", "AC_pCa_L", "AC_g_K1", "AC_g_Kr", "AC_g_Ks", "AC_g_Kp", "AC_g_Kur", "AC_g_to", "AC_g_f", "AC_g_b_Na", "AC_g_b_Ca", "AC_g_Cl_b", "AC_g_Cl_K", "AC_pNaK", "AC_KmKo", "AC_KmNai", "AC_p_NaCa", "AC_KmCai", "AC_KmCao", "AC_KmNai_ncx", "AC_KmNao", "AC_ksat", "AC_nu", "AC_p_Ca_L", "AC_P_up", "AC_K_up", "AC_n_CICR", "AC_K_rel", "AC_CaMK_0", "AC_b_CaMK", "AC_K_m_CaMK", "AC_n", "AC_p_Ca_p", "AC_K_pCa", "AC_Vc", "AC_Vsr", "AC_Vjsr", "AC_Vss"),
        recommended_exports=("V", "Cass"),
        compatible_tissues=("epicardialCells", "mCells", "endocardialCells"),
        compatible_solvers=("monodomainSolver", "bidomainSolver", "singleCellSolver"),
        species=("human",),
        cardiac_region=("ventricle",),
        model_type="ionic",
        description="Evolution of ORd with dynamic chloride handling (Tomek et al. 2019).",
    ),
    "Trovato": IonicModelEntry(
        states=("Vm", "m", "hf", "hs", "j", "hsp", "jp", "fLCa", "mL", "hL", "hLp", "a", "iF", "iS", "ap", "iFp", "iSp", "d", "ff", "fs", "fcaf", "fcas", "jca", "nca", "ffp", "fcafp", "fcasp", "xrf", "xrs", "xs1", "xs2", "xk1", "y", "oa", "oi", "r", "u", "Nai", "Cass", "Cajsr", "Cansr", "CaMKt"),
        algebraic=("AV_INa", "AV_INaL", "AV_ICaL", "AV_ICaNa", "AV_ICaK", "AV_IpCa", "AV_INaCa", "AV_INaK", "AV_IKr", "AV_IKs", "AV_IK1", "AV_IKb", "AV_INab", "AV_ICab", "AV_Ist", "AV_Iion", "AV_jup", "AV_jtr", "AV_jrel", "AV_Istim", "Iion_cm"),
        constants=("AC_nao", "AC_cao", "AC_ko", "AC_zNa", "AC_zCa", "AC_zK", "AC_R", "AC_T", "AC_F", "AC_g_Na", "AC_g_NaL", "AC_pCa_L", "AC_g_K1", "AC_g_Kr", "AC_g_Ks", "AC_g_Kp", "AC_g_to", "AC_g_f", "AC_g_b_Na", "AC_g_b_Ca", "AC_pNaK", "AC_KmKo", "AC_KmNai", "AC_p_NaCa", "AC_KmCai", "AC_KmCao", "AC_KmNai_ncx", "AC_KmNao", "AC_ksat", "AC_nu", "AC_p_Ca_L", "AC_P_up", "AC_K_up", "AC_n_CICR", "AC_K_rel", "AC_CaMK_0", "AC_b_CaMK", "AC_K_m_CaMK", "AC_n", "AC_p_Ca_p", "AC_K_pCa", "AC_Vc", "AC_Vsr", "AC_Vjsr", "AC_Vss"),
        recommended_exports=("Vm", "Cass"),
        compatible_tissues=("myocyte",),
        compatible_solvers=("monodomainSolver", "bidomainSolver", "singleCellSolver"),
        species=("human",),
        cardiac_region=("ventricle",),
        model_type="ionic",
        description="Recent update to ORd focusing on early after-depolarisations (EADs) (Trovato et al. 2020).",
    ),
    "monodomainFDAManufactured": IonicModelEntry(
        states=("phi",),
        algebraic=("Iion_cm",),
        constants=(),
        recommended_exports=("phi",),
        compatible_tissues=("manufactured",),
        compatible_solvers=("monodomainSolver", "singleCellSolver"),
        species=("generic",),
        cardiac_region=("manufactured",),
        model_type="manufactured",
        description="Manufactured (verification) monodomain ionic model; not physiological.",
    ),
    "bidomainFDAManufactured": IonicModelEntry(
        states=("phi",),
        algebraic=("Iion_cm",),
        constants=("kappa",),
        recommended_exports=("phi",),
        compatible_tissues=("manufactured",),
        compatible_solvers=("bidomainSolver", "singleCellSolver"),
        species=("generic",),
        cardiac_region=("manufactured",),
        model_type="manufactured",
        description="Manufactured (verification) bidomain ionic model; not physiological.",
    ),
}

ACTIVE_TENSION_MODEL_CATALOG: Final[dict[str, ActiveTensionModelEntry]] = {
    "GoktepeKuhl": ActiveTensionModelEntry(
        states=("Ta",),
        algebraic=("AV_e", "AV_Vm", "AV_u"),
        constants=("AC_Vr", "AC_eInfty", "AC_e0", "AC_eXi", "AC_Vshift", "AC_kTa"),
        rates=("Ta",),
        recommended_exports=("Ta",),
        description="Goktepe-Kuhl active tension model (2004).",
    ),
    "NashPanfilov": ActiveTensionModelEntry(
        states=("Ta",),
        algebraic=("AV_u", "AV_e"),
        constants=("AC_Vp", "AC_Vr", "AC_Vth", "AC_e0", "AC_kTa", "AC_a"),
        rates=("Ta",),
        recommended_exports=("Ta",),
        description="Nash-Panfilov active tension model (2004).",
    ),
}


def get_ionic_model_entry(name: str) -> IonicModelEntry:
    """
    Return the catalog entry for the named ionic model.

    Args:
        name: The ionic model name (e.g. 'TNNP', 'BuenoOrovio').

    Returns:
        The IonicModelEntry for that model.

    Raises:
        KeyError: If the model is not in the catalog.
    """
    if name not in IONIC_MODEL_CATALOG:
        raise KeyError(
            f"Unknown ionic model '{name}'. "
            f"Available models: {', '.join(IONIC_MODEL_CATALOG.keys())}"
        )
    return IONIC_MODEL_CATALOG[name]


def list_compatible_ionic_models(
    solver: str, tissue: str | None = None
) -> list[str]:
    """
    Return ionic model names compatible with the given myocardium solver and optional tissue.

    Args:
        solver: The myocardium solver name (e.g. 'monodomainSolver', 'bidomainSolver', 'singleCellSolver').
        tissue: Optional tissue type filter (e.g. 'epicardialCells', 'myocyte').

    Returns:
        List of compatible model names.
    """
    compatible = []
    for model_name, entry in IONIC_MODEL_CATALOG.items():
        if solver not in entry.compatible_solvers:
            continue
        if tissue is not None and tissue not in entry.compatible_tissues:
            continue
        compatible.append(model_name)
    return compatible


def get_active_tension_entry(name: str) -> ActiveTensionModelEntry:
    """
    Return the catalog entry for the named active tension model.

    Args:
        name: The active tension model name (e.g. 'GoktepeKuhl').

    Returns:
        The ActiveTensionModelEntry for that model.

    Raises:
        KeyError: If the model is not in the catalog.
    """
    if name not in ACTIVE_TENSION_MODEL_CATALOG:
        raise KeyError(
            f"Unknown active tension model '{name}'. "
            f"Available models: {', '.join(ACTIVE_TENSION_MODEL_CATALOG.keys())}"
        )
    return ACTIVE_TENSION_MODEL_CATALOG[name]


def list_models_by_region(cardiac_region: str) -> list[str]:
    """Return ionic model names whose cardiac_region tuple contains the given region."""
    return [
        name
        for name, entry in IONIC_MODEL_CATALOG.items()
        if cardiac_region in entry.cardiac_region
    ]


def list_models_by_species(species: str) -> list[str]:
    """Return ionic model names whose species tuple contains the given species."""
    return [
        name
        for name, entry in IONIC_MODEL_CATALOG.items()
        if species in entry.species
    ]
