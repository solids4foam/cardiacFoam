"""
Ionic Model Catalog

A static, build-time-generated catalog of ionic models.
This module exposes the exact variables each model supports so an autonomous
agent can plan ionic outputVariables without running the solver.

All variable names are extracted from C++ source files and are guaranteed to be
exact.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Final

from .active_tension_catalog import (
    ACTIVE_TENSION_MODEL_CATALOG,
    ActiveTensionModelEntry,
    get_active_tension_entry,
)


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

    aliases: tuple[str, ...] = ()
    """Alternative names for this model."""

    recommended_ode_step: float = 1e-5
    """Suggested ODE timestep for stable integration."""

    recommended_stimulus_duration: float | None = 0.002
    """Suggested stimulus pulse duration in milliseconds."""

    recommended_stimulus_intensity: float | None = 80000.0
    """Suggested stimulus current intensity in pA."""

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
        description="Two-variable phenomenological model; not species-specific (Aliev & Panfilov 1996). Computationally very cheap — suited for large-scale qualitative propagation studies and parameter sweeps. Does not reproduce ion concentrations or realistic AP morphology. Dimensionless units: stimulus intensity ~0.5, duration ~1.0 time unit.",
        aliases=("AP model", "Aliev-Panfilov", "AP96"),
        recommended_ode_step=1e-3,
        recommended_stimulus_duration=1.0,
        recommended_stimulus_intensity=0.5,
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
        description="Four-variable phenomenological model reproducing better AP morphology than AlievPanfilov; not species-specific (Bueno-Orovio et al. 2008). Computationally cheap. No ion dynamics. Compatible with epicardial, M-cell, and endocardial tissue variants. Dimensionless units.",
        aliases=("BO model", "Bueno-Orovio", "minimal ventricular model", "BO08"),
        recommended_ode_step=1e-3,
        recommended_stimulus_duration=1.0,
        recommended_stimulus_intensity=0.5,
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
        description="Human atrial ionic model (Courtemanche et al. 1998). 21 states. Suited for atrial fibrillation and atrial remodelling studies. Computationally moderate.",
        aliases=("CRN", "Courtemanche-Ramirez-Nattel", "human atrial model", "CRN98"),
        recommended_ode_step=1e-5,
        recommended_stimulus_duration=0.002,
        recommended_stimulus_intensity=80000.0,
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
        description="Human sinoatrial node model with spontaneous pacing (Fabbri et al. 2017). 33 states. Computationally demanding. Use for SAN automaticity and pacemaker current studies. Not a working myocyte model — does not produce AP morphology typical of atrial or ventricular cells.",
        aliases=("Fabbri-Fantini", "sinoatrial node model", "SAN model", "pacemaker model"),
        recommended_ode_step=1e-5,
        recommended_stimulus_duration=0.002,
        recommended_stimulus_intensity=80000.0,
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
        description="Pig ventricular ionic model (Gaur & Rudy 2011). Use when pig (guinea pig) animal experiments are the reference. Computationally moderate.",
        aliases=("Gaur-Rudy", "guinea pig ventricular model"),
        recommended_ode_step=1e-5,
        recommended_stimulus_duration=0.002,
        recommended_stimulus_intensity=80000.0,
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
        description="Human ventricular ionic model with detailed calcium signalling and electrolytes (Grandi et al. 2010). 29 states. Suited for Ca²⁺ handling and atrial fibrillation drug studies. Computationally demanding.",
        aliases=("Grandi-Pasqualini-Bers", "GPB model"),
        recommended_ode_step=1e-5,
        recommended_stimulus_duration=0.002,
        recommended_stimulus_intensity=80000.0,
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
        description="O'Hara-Rudy; industry-standard 41-state human ventricular ionic model (O'Hara et al. 2011). Preferred for drug-induced arrhythmia and late INa studies. Computationally demanding.",
        aliases=("O'Hara-Rudy", "OR model", "O'Hara 2011", "ORd2011"),
        recommended_ode_step=1e-5,
        recommended_stimulus_duration=0.002,
        recommended_stimulus_intensity=80000.0,
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
        description="Human Purkinje fibre ionic model for specialised cardiac conduction (Stewart et al. 2009). 16 states. Computationally moderate. Use with monodomain1DSolver for Purkinje network simulations. Not a ventricular myocyte model — do not use for ventricular tissue simulations.",
        aliases=("Stewart-Aslanidi-Noble", "human Purkinje model"),
        recommended_ode_step=1e-5,
        recommended_stimulus_duration=0.002,
        recommended_stimulus_intensity=80000.0,
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
        description="Ten Tusscher-Noble-Noble-Panfilov; widely used 21-state human ventricular ionic model (ten Tusscher & Panfilov 2006). Well-validated for 3-D reentry and drug screening. Tissue variants: epicardial, M-cell, endocardial. Computationally moderate.",
        aliases=("ten Tusscher", "TT04", "TT06", "ten Tusscher-Panfilov", "TT2006"),
        recommended_ode_step=1e-5,
        recommended_stimulus_duration=0.002,
        recommended_stimulus_intensity=50000.0,
    ),
    "ToRORd_dynCl": IonicModelEntry(
        states=("V", "CaMKt", "Nai", "Nass", "Ki", "Kss", "Cass", "Cansr", "Cajsr", "Cai", "Cli", "Clss", "INa_m", "INa_h", "INa_j", "INa_hp", "INa_jp", "INaL_mL", "INaL_hL", "INaL_hLp", "Ito_a", "Ito_iF", "Ito_iS", "Ito_ap", "Ito_iFp", "Ito_iSp", "ICaL_d", "ICaL_ff", "ICaL_fs", "ICaL_fcaf", "ICaL_fcas", "ICaL_jca", "ICaL_ffp", "ICaL_fcafp", "ICaL_nca_ss", "ICaL_nca_i", "IKr_C1", "IKr_C2", "IKr_C3", "IKr_I", "IKr_O", "IKs_xs1", "IKs_xs2", "Jrel_np", "Jrel_p"),
        algebraic=("AV_CaMKb", "AV_CaMKa", "AV_IpCa_IpCa", "AV_Jdiff", "AV_JdiffCl", "AV_JdiffK", "AV_JdiffNa", "AV_time", "AV_Jtr", "AV_Jleak", "AV_fJupp", "AV_Jupnp", "AV_Jupp", "AV_Jup", "AV_ECl", "AV_EClss", "AV_EK", "AV_ENa", "AV_EKs", "AV_IClCa_junc", "AV_IClCa_sl", "AV_IClb", "AV_IClCa", "AV_aK1", "AV_bK1", "AV_K1ss", "AV_IK1_IK1", "AV_xkb", "AV_IKb_IKb", "AV_KsCa", "AV_txs1", "AV_txs2", "AV_xs1ss", "AV_xs2ss", "AV_IKs_IKs", "AV_ah", "AV_aj", "AV_bh", "AV_bj", "AV_fINap", "AV_hss", "AV_hssp", "AV_mss", "AV_tm", "AV_INa_INa", "AV_jss", "AV_th", "AV_tj", "AV_tjp", "AV_fINaLp", "AV_hLss", "AV_hLssp", "AV_mLss", "AV_tmL", "AV_INaL_INaL", "AV_I_katp_I_katp", "AV_fItop", "AV_AiF", "AV_ass", "AV_assp", "AV_delta_epi", "AV_dti_develop", "AV_dti_recover", "AV_iss", "AV_ta", "AV_tiF_b", "AV_tiS_b", "AV_AiS", "AV_tiF", "AV_tiS", "AV_i", "AV_ip", "AV_tiFp", "AV_tiSp", "AV_Ito_Ito", "AV_Afcaf", "AV_Ii", "AV_Iss", "AV_dss", "AV_fICaLp", "AV_fss", "AV_jcass", "AV_km2n", "AV_tfcaf", "AV_tfcas", "AV_tff", "AV_tfs", "AV_Afcas", "AV_anca_i", "AV_anca_ss", "AV_fcass", "AV_td", "AV_tfcafp", "AV_tffp", "AV_f", "AV_fca", "AV_fcap", "AV_fp", "AV_gamma_cai", "AV_gamma_cass", "AV_gamma_ki", "AV_gamma_kss", "AV_gamma_nai", "AV_gamma_nass", "AV_IKr_IKr", "AV_allo_i", "AV_allo_ss", "AV_h4_i", "AV_h4_ss", "AV_h5_i", "AV_h5_ss", "AV_h6_i", "AV_h6_ss", "AV_k6_i", "AV_k6_ss", "AV_P", "AV_b3", "AV_Bcajsr", "AV_Bcass", "AV_Bcai", "AV_vffrt", "AV_vfrt", "AV_fJrelp", "AV_Jrel", "AV_tau_rel_b", "AV_tau_rel", "AV_tau_relp_b", "AV_tau_relp", "AV_PhiCaK_i", "AV_PhiCaK_ss", "AV_PhiCaL_i", "AV_PhiCaL_ss", "AV_PhiCaNa_i", "AV_PhiCaNa_ss", "AV_ICab_ICab", "AV_alpha", "AV_alpha_2", "AV_alpha_C2ToI", "AV_alpha_i", "AV_beta", "AV_beta_2", "AV_beta_i", "AV_hca", "AV_hna", "AV_Knai", "AV_Knao", "AV_INab_INab", "AV_ICaK_i", "AV_ICaK_ss", "AV_ICaL_i", "AV_ICaL_ss", "AV_ICaNa_i", "AV_ICaNa_ss", "AV_beta_ItoC2", "AV_h1_i", "AV_h1_ss", "AV_h7_i", "AV_h7_ss", "AV_a1", "AV_a3", "AV_b2", "AV_b4", "AV_ICaK", "AV_ICaL_ICaL", "AV_ICaNa", "AV_h2_i", "AV_h2_ss", "AV_h3_i", "AV_h3_ss", "AV_h8_i", "AV_h8_ss", "AV_h9_i", "AV_h9_ss", "AV_x1", "AV_x2", "AV_x3", "AV_x4", "AV_Jrel_inf_b", "AV_Jrel_infp_b", "AV_k3p_i", "AV_k3p_ss", "AV_k3pp_i", "AV_k3pp_ss", "AV_k4p_i", "AV_k4p_ss", "AV_k4pp_i", "AV_k4pp_ss", "AV_k7_i", "AV_k7_ss", "AV_k8_i", "AV_k8_ss", "AV_E1", "AV_E2", "AV_E3", "AV_E4", "AV_Jrel_inf", "AV_Jrel_infp", "AV_k3_i", "AV_k3_ss", "AV_k4_i", "AV_k4_ss", "AV_JnakK", "AV_JnakNa", "AV_x1_i", "AV_x1_ss", "AV_x2_i", "AV_x2_ss", "AV_x3_i", "AV_x3_ss", "AV_x4_i", "AV_x4_ss", "AV_INaK_INaK", "AV_E1_i", "AV_E1_ss", "AV_E2_i", "AV_E2_ss", "AV_E3_i", "AV_E3_ss", "AV_E4_i", "AV_E4_ss", "AV_JncxCa_i", "AV_JncxCa_ss", "AV_JncxNa_i", "AV_JncxNa_ss", "AV_INaCa_i", "AV_INaCa_ss", "Istim", "Iion_cm"),
        constants=("AC_CaMKo", "AC_KmCaM", "AC_KmCaMK", "AC_aCaMK", "AC_bCaMK", "AC_GpCa", "AC_KmCap", "AC_L", "AC_rad", "AC_Ageo", "AC_vcell", "AC_Acap", "AC_vjsr", "AC_vmyo", "AC_vnsr", "AC_vss", "AC_tauCa", "AC_tauCl", "AC_tauK", "AC_tauNa", "AC_cao", "AC_clo", "AC_ko", "AC_nao", "AC_F", "AC_R", "AC_T", "AC_zca", "AC_zcl", "AC_zk", "AC_zna", "AC_Jup_b", "AC_upScale", "AC_PKNa", "AC_Fjunc", "AC_GClCa", "AC_GClb", "AC_KdClCa", "AC_GK1_b", "AC_GK1", "AC_GKb_b", "AC_GKb", "AC_GKs_b", "AC_GKs", "AC_GNa", "AC_GNaL_b", "AC_thL", "AC_GNaL", "AC_thLp", "AC_A_atp", "AC_K_atp", "AC_K_o_n", "AC_fkatp", "AC_gkatp", "AC_akik", "AC_bkik", "AC_EKshift", "AC_Gto_b", "AC_Gto", "AC_Aff", "AC_ICaL_fractionSS", "AC_Io", "AC_Kmn", "AC_PCa_b", "AC_dielConstant", "AC_k2n", "AC_offset", "AC_tjca", "AC_vShift", "AC_Afs", "AC_PCa", "AC_constA", "AC_PCaK", "AC_PCaNa", "AC_PCap", "AC_gamma_cao", "AC_gamma_ko", "AC_gamma_nao", "AC_PCaKp", "AC_PCaNap", "AC_PCab", "AC_GKr_b", "AC_alpha_1", "AC_beta_1", "AC_GKr", "AC_Gncx_b", "AC_INaCa_fractionSS", "AC_KmCaAct", "AC_kasymm", "AC_kcaoff", "AC_kcaon", "AC_kna1", "AC_kna2", "AC_kna3", "AC_qca", "AC_qna", "AC_wca", "AC_wna", "AC_wnaca", "AC_Gncx", "AC_h10_i", "AC_h10_ss", "AC_k2_i", "AC_k2_ss", "AC_k5_i", "AC_k5_ss", "AC_h11_i", "AC_h11_ss", "AC_h12_i", "AC_h12_ss", "AC_k1_i", "AC_k1_ss", "AC_H", "AC_Khp", "AC_Kki", "AC_Kko", "AC_Kmgatp", "AC_Knai0", "AC_Knao0", "AC_Knap", "AC_Kxkur", "AC_MgADP", "AC_MgATP", "AC_Pnak_b", "AC_delta", "AC_eP", "AC_k1m", "AC_k1p", "AC_k2m", "AC_k2p", "AC_k3m", "AC_k3p", "AC_k4m", "AC_k4p", "AC_Pnak", "AC_a2", "AC_a4", "AC_b1", "AC_PNab", "AC_BSLmax", "AC_BSRmax", "AC_KmBSL", "AC_KmBSR", "AC_cmdnmax_b", "AC_csqnmax", "AC_kmcmdn", "AC_kmcsqn", "AC_kmtrpn", "AC_trpnmax", "AC_cmdnmax", "AC_Jrel_b", "AC_bt", "AC_cajsr_half", "AC_a_rel", "AC_btp", "AC_a_relp", "deltaEpiAmp", "deltaEpiShift", "deltaEpiScale", "jrelTissueScale"),
        recommended_exports=("V", "Cass"),
        compatible_tissues=("epicardialCells", "mCells", "endocardialCells"),
        compatible_solvers=("monodomainSolver", "bidomainSolver", "singleCellSolver"),
        species=("human",),
        cardiac_region=("ventricle",),
        model_type="ionic",
        description="Evolution of ORd with dynamic chloride handling; 45 states (Tomek et al. 2019). Suited for studies involving chloride currents or arrhythmia. Computationally demanding.",
        aliases=("ToR-ORd", "Tomek", "ToRORd", "Tomek-Rodriguez-ORd", "dynamic chloride"),
        recommended_ode_step=1e-5,
        recommended_stimulus_duration=0.002,
        recommended_stimulus_intensity=80000.0,
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
        description="ORd update focusing on early after-depolarisations (EADs); 41 states (Trovato et al. 2020). Computationally demanding — same cost as ORd. Use for EAD and triggered activity studies. Not a general-purpose ventricular replacement for ORd; purpose-built for EAD-prone regimes.",
        aliases=("Trovato 2020", "updated ORd", "EAD model"),
        recommended_ode_step=1e-5,
        recommended_stimulus_duration=0.002,
        recommended_stimulus_intensity=80000.0,
    ),
    "TWorld": IonicModelEntry(
        states=("v", "camk_trap", "camk_f_ICaL", "camk_f_RyR", "camk_f_PLB", "casig_serca_trap", "buffers_NaBj", "buffers_NaBsl", "buffers_TnClow", "buffers_TnCHc", "buffers_TnCHm", "buffers_CaM", "buffers_Myosin_ca", "buffers_Myosin_mg", "buffers_SRB", "buffers_SLLj", "buffers_SLLsl", "buffers_SLHj", "buffers_SLHsl", "buffers_Csqn", "naj", "nasl", "nai", "ki", "cli", "casr", "caj", "casl", "cai", "m", "h", "j", "hp", "jp", "m_P", "h_P", "j_P", "hp_P", "jp_P", "mL", "hL", "hLp", "d", "ff", "fs", "fcaf", "fcas", "jca", "ffp", "fcafp", "nca", "nca_i", "d_P", "ff_P", "fs_P", "fcaf_P", "fcas_P", "fBPf", "fcaBPf", "ical_pureCDI_junc", "ical_pureCDI_sl", "xtos", "ytos", "xtof", "ytof", "xtos_p", "xtof_p", "ytos_p", "ytof_p", "C0", "C1", "C2", "I", "O", "xs_junc", "xs_sl", "jrel_icaldep_act", "jrel_icaldep_f1", "jrel_icaldep_f2", "ryr_R", "ryr_O", "ryr_I", "ryr_CaRI", "ryr_R_p", "ryr_O_p", "ryr_I_p", "ryr_CaRI_p", "contraction_TmBlocked", "contraction_XW", "contraction_XS", "contraction_ZETAS", "contraction_ZETAW", "contraction_Ca_TRPN"),
        algebraic=("AV_CaMK_Phos_ss_ICaL", "AV_CaMK_Phos_ss_PLB", "AV_CaMK_Phos_ss_RyR", "AV_CaMK_active", "AV_bound", "AV_bound_serca", "AV_casig_SERCA_act", "AV_Afcaf", "AV_Afcas", "AV_ICaK", "AV_ICaK_i_BP", "AV_ICaK_i_CaMK", "AV_ICaK_i_NP", "AV_ICaK_i_PKA", "AV_ICaK_junc", "AV_ICaK_sl", "AV_ICaK_ss_BP", "AV_ICaK_ss_CaMK", "AV_ICaK_ss_NP", "AV_ICaK_ss_PKA", "AV_ICaL_ICaL", "AV_ICaL_i_BP", "AV_ICaL_i_CaMK", "AV_ICaL_i_NP", "AV_ICaL_i_PKA", "AV_ICaL_junc", "AV_ICaL_sl", "AV_ICaL_ss_BP", "AV_ICaL_ss_CaMK", "AV_ICaL_ss_NP", "AV_ICaL_ss_PKA", "AV_ICaNa", "AV_ICaNa_i_BP", "AV_ICaNa_i_CaMK", "AV_ICaNa_i_NP", "AV_ICaNa_i_PKA", "AV_ICaNa_junc", "AV_ICaNa_sl", "AV_ICaNa_ss_BP", "AV_ICaNa_ss_CaMK", "AV_ICaNa_ss_NP", "AV_ICaNa_ss_PKA", "AV_ICa_tot", "AV_Ii", "AV_Iss", "AV_PhiCaK_i", "AV_PhiCaK_ss", "AV_PhiCaL_i", "AV_PhiCaL_ss", "AV_PhiCaNa_i", "AV_PhiCaNa_ss", "AV_anca", "AV_anca_i", "AV_dPss", "AV_dss", "AV_f", "AV_fBP", "AV_fBPss", "AV_fICaL_BP", "AV_fICaL_CaMKonly", "AV_fICaL_PKAonly", "AV_fICaLp", "AV_f_P", "AV_fca", "AV_fcaBP", "AV_fcaBPss", "AV_fcap", "AV_fcap_P", "AV_fcass", "AV_fcass_P", "AV_fp", "AV_fss", "AV_fss_P", "AV_gamma_cai", "AV_gamma_cass", "AV_gamma_ki", "AV_gamma_kss", "AV_gamma_nai", "AV_gamma_nass", "AV_jcass", "AV_km2n", "AV_sigmoidTransition", "AV_sigmoidTransition2", "AV_tauTransition", "AV_tauTransition2", "AV_td", "AV_tfcaf", "AV_tfcafp", "AV_tfcas", "AV_tff", "AV_tffp", "AV_tfs", "AV_ICab_ICab", "AV_ICab_junc", "AV_ICab_sl", "AV_IClCa", "AV_IClCa_junc", "AV_IClCa_sl", "AV_IClb", "AV_IK1_IK1", "AV_K1ss", "AV_aK1", "AV_bK1", "AV_IKb_IKb", "AV_xkb", "AV_IKr_IKr", "AV_IKr_alpha", "AV_alpha_2", "AV_alpha_C2ToI", "AV_alpha_i", "AV_IKr_beta", "AV_beta_2", "AV_beta_ItoC2", "AV_beta_i", "AV_IKs_IKs", "AV_IKs_junc", "AV_IKs_sl", "AV_INa_INa", "AV_INaBase", "AV_INaBase_BP", "AV_INaBase_CaMK", "AV_INaBase_NP", "AV_INaBase_PKA", "AV_INaj", "AV_INasl", "AV_ah", "AV_aj", "AV_bh", "AV_bj", "AV_fINa_BP", "AV_fINa_CaMKonly", "AV_fINa_PKAonly", "AV_fINap", "AV_hss", "AV_hss_P", "AV_hssp", "AV_hssp_P", "AV_jss", "AV_jss_P", "AV_jssp_P", "AV_mss", "AV_mss_P", "AV_th", "AV_tj", "AV_tjp", "AV_tm", "AV_INaCa_INaCa", "AV_INaCa_i", "AV_INaCa_ss", "AV_JncxCa_i", "AV_JncxCa_ss", "AV_JncxNa_i", "AV_JncxNa_ss", "AV_allo_i", "AV_allo_ss", "AV_INaK_INaK", "AV_INaKj", "AV_INaKsl", "AV_fnak", "AV_GNaL", "AV_INaL_INaL", "AV_INaLj", "AV_INaLsl", "AV_fINaLp", "AV_hLss", "AV_hLssp", "AV_mLss", "AV_tmL", "AV_INab_INab", "AV_INabj", "AV_INabsl", "AV_IpCa_IpCa", "AV_IpCa_junc", "AV_IpCa_sl", "AV_Ito_Ito", "AV_Itof", "AV_Itos", "AV_Jserca", "AV_Jserca_np", "AV_Jserca_p", "AV_Ta", "AV_XU", "AV_J_SRCarel", "AV_J_SRleak", "AV_Jrel_ICaLdep", "AV_Jrel_inf", "AV_kCaSR", "AV_kiSRCa", "AV_koSRCa", "AV_vffrt", "AV_vfrt", "AV_ECaj", "AV_ECasl", "AV_ECl", "AV_EK", "AV_EKs", "AV_ENaj", "AV_ENasl", "Istim", "Iion_cm"),
        constants=("AC_CaMK0", "AC_K_Phos_CaMK", "AC_Km_CaMK_Ca", "AC_PP1_tot", "AC_Whole_cell_PP1", "AC_CaMK_alpha", "AC_alpha_serca", "AC_CaMK_beta", "AC_tau_cal", "AC_tau_plb", "AC_tau_ryr", "AC_Aff", "AC_Afs", "AC_ICaL_fractionSS", "AC_Io", "AC_Kmn", "AC_PCa", "AC_PCaK", "AC_PCaK_P", "AC_PCaKp", "AC_PCaNa", "AC_PCaNa_P", "AC_PCaNap", "AC_PCa_P", "AC_PCa_Pb", "AC_PCa_b", "AC_PCap", "AC_constA", "AC_dielConstant", "AC_fICaLP", "AC_fICaL_P", "AC_gamma_cao", "AC_gamma_ko", "AC_gamma_nao", "AC_k2n", "AC_rateRecovery", "AC_tjca", "AC_GCab", "AC_GCab_b", "AC_GClCa", "AC_GClCa_b", "AC_GClb", "AC_GClb_b", "AC_KdClCa", "AC_GK1", "AC_GK1_b", "AC_IK1_celltype_factor", "AC_IK1_sex_factor", "AC_GKb", "AC_GKb_b", "AC_GKr", "AC_GKr_b", "AC_alpha_1", "AC_IKr_beta_1", "AC_IKr_celltype_factor", "AC_IKr_sex_factor", "AC_GKs", "AC_GKs_b", "AC_IKs_celltype_factor", "AC_gKs_factor", "AC_kPKA_IKs", "AC_IKs_sex_factor", "AC_GNa", "AC_GNa_P", "AC_GNa_b", "AC_fINa_P", "AC_Gncx", "AC_Gncx_b", "AC_INaCa_fractionSS", "AC_KmCaAct", "AC_INaCa_celltype_factor", "AC_IbarNaK", "AC_IbarNaK_b", "AC_KmKo", "AC_KmNaip", "AC_KmNaip_PKA", "AC_GNaL_b", "AC_thL", "AC_thLp", "AC_GNab", "AC_GNab_b", "AC_IbarSLCaP", "AC_IbarSLCaP_b", "AC_KmPCa", "AC_Q10SLCaP", "AC_Gto_fast", "AC_Gto_slow", "AC_fICaL_PKA", "AC_fIKs_PKA", "AC_fINaK_PKA", "AC_fINa_PKA", "AC_fMyBPC_PKA", "AC_fPLB_PKA", "AC_fTnI_PKA", "AC_Km_SERCA_Ca", "AC_Kmf", "AC_Kmf_p", "AC_Kmr", "AC_Max_Vmax_SERCA_Ca", "AC_Q10SRCaP", "AC_Vmax_SRCaP", "AC_Vmax_SRCaP_b", "AC_hillSRCaP", "AC_cellLength", "AC_cellRadius", "AC_vcell", "AC_vjunc", "AC_vmyo", "AC_vsl", "AC_vsr", "AC_A", "AC_Lfac", "AC_PKAForceMultiplier", "AC_TOT_A", "AC_TRPN_n", "AC_Tref", "AC_XSSS", "AC_XWSS", "AC_beta_0", "AC_contraction_beta_1", "AC_ca50", "AC_cds", "AC_cdw", "AC_dr", "AC_fPKA_TnI", "AC_fracTnIpo", "AC_contraction_gamma", "AC_gamma_wu", "AC_k_su", "AC_k_uw", "AC_k_ws", "AC_k_wu", "AC_koff", "AC_ktm_block", "AC_ktm_unblock", "AC_lambda", "AC_lambda0", "AC_lambda_max", "AC_lambda_min", "AC_lambda_rate", "AC_mu", "AC_nperm", "AC_nu", "AC_perm50", "AC_phi", "AC_wfrac", "AC_sex", "AC_cao", "AC_clo", "AC_ko", "AC_nao", "AC_Bmax_CaM", "AC_Bmax_Csqn", "AC_Bmax_Naj", "AC_Bmax_Nasl", "AC_Bmax_SLhighj", "AC_Bmax_SLhighsl", "AC_Bmax_SLlowj", "AC_Bmax_SLlowsl", "AC_Bmax_SR", "AC_Bmax_TnChigh", "AC_Bmax_TnClow", "AC_Bmax_myosin", "AC_J_ca_juncsl", "AC_J_ca_slmyo", "AC_J_na_juncsl", "AC_J_na_slmyo", "AC_koff_cam", "AC_koff_csqn", "AC_koff_myoca", "AC_koff_myomg", "AC_koff_na", "AC_koff_slh", "AC_koff_sll", "AC_koff_sr", "AC_koff_tnchca", "AC_koff_tnchmg", "AC_kon_cam", "AC_kon_csqn", "AC_kon_myoca", "AC_kon_myomg", "AC_kon_na", "AC_kon_slh", "AC_kon_sll", "AC_kon_sr", "AC_kon_tnchca", "AC_kon_tnchmg", "AC_mgi", "AC_Cmem", "AC_Fjunc", "AC_Fsl", "AC_ICaLPCa_multiplier", "AC_ICab_multiplier", "AC_IClCa_multiplier", "AC_IClb_multiplier", "AC_IK1_multiplier", "AC_IKb_multiplier", "AC_IKr_multiplier", "AC_IKs_multiplier", "AC_INaCa_multiplier", "AC_INaK_multiplier", "AC_INaL_multiplier", "AC_INa_multiplier", "AC_INab_multiplier", "AC_IpCa_multiplier", "AC_Itof_multiplier", "AC_Itos_multiplier", "AC_Jrel_multiplier", "AC_Jup_multiplier", "AC_F", "AC_Qpow", "AC_R", "AC_T", "AC_zca", "AC_zcl", "AC_zk", "AC_zna", "AC_PNaK", "AC_CI_to_RI", "AC_MaxSR", "AC_MinSR", "AC_a_rel", "AC_baseRateCaI", "AC_bt", "AC_caExpFactor", "AC_caExpFactor2", "AC_caTransFactor", "AC_caTransFactor2", "AC_caTransFactor2p", "AC_directRelMidpoint", "AC_ec50SR", "AC_ecCaI", "AC_kiCa", "AC_kim", "AC_koCa", "AC_kom", "AC_ks", "AC_maxCaI", "AC_minCaI", "AC_steepnessCaI", "AC_steepnessCaSR", "AC_tauInact", "AC_tauInact2", "gnalTissueScale"),
        recommended_exports=("v", "cai", "contraction_Ca_TRPN"),
        compatible_tissues=("epicardialCells", "mCells", "endocardialCells"),
        compatible_solvers=("monodomainSolver", "bidomainSolver", "singleCellSolver"),
        species=("human",),
        cardiac_region=("ventricle",),
        model_type="ionic",
        description="T-World 2024 — large human ventricular ionic model with CaMKII signalling, Cl⁻ and Cl(Ca) currents, junction/sl compartments, and integrated active-tension contraction (Land model). 93 states. Computationally very demanding. Use for studies requiring detailed Ca²⁺ handling, CaMKII-dependent remodelling, or coupled electromechanics.",
        aliases=("T-World", "TWorld 2024"),
        recommended_ode_step=1e-5,
        recommended_stimulus_duration=0.002,
        recommended_stimulus_intensity=80000.0,
        notes="Contains active-tension states (contraction_*); the contraction subsystem is part of the ionic ODE system — no separate activeTensionModel needed.",
    ),
    "PerisYague": IonicModelEntry(
        states=("membrane_V", "sodium_Nai", "potassium_Ki", "chloride_Cli", "calcium_Cai", "calcium_CaUp", "calcium_CaRel", "ina_m", "ina_h", "ina_j", "ikur_ua", "ikur_uif", "ikur_uis", "ikr_xr", "iks_xs", "ical_d", "ical_f", "ical_fCa", "iclca_qCa", "ryr_u", "ryr_v", "ryr_w"),
        algebraic=("AV_ICaL", "AV_IK1", "AV_IKr", "AV_IKs", "AV_IKur", "AV_INa", "AV_INaCa", "AV_INaK", "AV_IClCa", "AV_IbCa", "AV_IbNa", "AV_IbK", "AV_IpCa", "AV_ina_m_inf", "AV_ina_m_tau", "AV_ina_h_inf", "AV_ina_h_tau", "AV_ina_j_inf", "AV_ina_j_tau", "AV_ical_d_inf", "AV_ical_d_tau", "AV_ical_f_inf", "AV_ical_f_tau", "AV_ikur_ua_inf", "AV_ikur_ua_tau", "AV_ikur_uif_inf", "AV_ikur_uif_tau", "AV_ikur_uis_inf", "AV_ikur_uis_tau", "AV_gKur", "AV_ikr_xr_inf", "AV_ikr_xr_tau", "AV_iks_xs_inf", "AV_iks_xs_tau", "AV_iclca_qCa_inf", "AV_iclca_qCa_tau", "AV_fNaK", "AV_ryr_w_inf", "AV_ryr_w_tau", "Iion_cm", "Istim"),
        constants=("AC_CMDN_max", "AC_CSQN_max", "AC_Ca_up_max", "AC_Cao", "AC_Cm", "AC_ECaL", "AC_F", "AC_FRT", "AC_INaCa_max", "AC_INaK_max", "AC_I_diff", "AC_I_up_max", "AC_KQ10", "AC_krel", "AC_K_up", "AC_KmCa", "AC_KmKo", "AC_KmNa", "AC_KmNai", "AC_Km_CMDN", "AC_Km_CSQN", "AC_Km_TRPN", "AC_Ko", "AC_Nao", "AC_Clo", "AC_R", "AC_RTF", "AC_T", "AC_TRPN_max", "AC_V_cell", "AC_V_i", "AC_V_rel", "AC_V_up", "AC_c1", "AC_c2", "AC_cajsr_u_tau", "AC_gNa", "AC_gK1", "AC_gKr", "AC_gKs", "AC_gKur_amp", "AC_gClCa", "AC_gbCa", "AC_gbNa", "AC_gbK", "AC_ical_fCa_tau", "AC_gamma", "AC_ksat", "AC_sigma", "AC_tau_tr", "AC_gCaL", "AC_IpCa_max"),
        recommended_exports=("membrane_V", "calcium_Cai"),
        compatible_tissues=("myocyte",),
        compatible_solvers=("monodomainSolver", "bidomainSolver", "singleCellSolver"),
        species=("human",),
        cardiac_region=("atrium",),
        model_type="ionic",
        description="Human atrial ionic model with chloride current IClCa; 22 states (Peris-Yaguë et al. 2022). Extension of TNNP with Cl⁻ handling and two-component IKur inactivation. Suited for atrial Cl⁻-dependent arrhythmia and drug studies. Computationally moderate.",
        aliases=("Peris-Yague", "Peris Yague 2022", "atrial Cl model"),
        recommended_ode_step=1e-5,
        recommended_stimulus_duration=0.002,
        recommended_stimulus_intensity=50000.0,
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
        description="Manufactured solution for monodomain solver convergence testing and verification (FDA benchmark). Computational cost is not applicable — this is a single synthetic variable, not a physiological model. Not for physiological studies.",
        aliases=("FDA manufactured", "monodomain manufactured solution"),
        recommended_ode_step=1e-4,
        recommended_stimulus_duration=None,
        recommended_stimulus_intensity=None,
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
        description="Manufactured solution for bidomain solver convergence testing and verification (FDA benchmark). Computational cost is not applicable — this is a single synthetic variable, not a physiological model. Not for physiological studies.",
        aliases=("FDA manufactured", "bidomain manufactured solution"),
        recommended_ode_step=1e-4,
        recommended_stimulus_duration=None,
        recommended_stimulus_intensity=None,
    ),
    "bathBidomainFDAManufactured": IonicModelEntry(
        states=("V", "u1", "u2", "u3"),
        algebraic=("Iion_cm",),
        constants=("Cm", "Beta", "Chi"),
        recommended_exports=("V", "u1", "u2", "u3"),
        compatible_tissues=("manufactured",),
        compatible_solvers=("bidomainSolver",),
        species=("generic",),
        cardiac_region=("manufactured",),
        model_type="manufactured",
        description="Manufactured solution for bidomain+bath (extracellular potential) solver convergence testing and verification. Requires a potentialDomain block and the manufacturedFDABathBidomainVerifier hook. Not for physiological studies.",
        aliases=("FDA bath manufactured", "bath bidomain manufactured solution"),
        recommended_ode_step=1e-4,
        recommended_stimulus_duration=None,
        recommended_stimulus_intensity=None,
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
