#----------------------------------------------------------------------------#
# License
#     This file is part of cardiacFoam.
#
#     cardiacFoam is free software: you can redistribute it and/or modify it
#     under the terms of the GNU General Public License as published by the
#     Free Software Foundation, either version 3 of the License, or (at your
#     option) any later version.
#
#     cardiacFoam is distributed in the hope that it will be useful, but
#     WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#     General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with cardiacFoam.  If not, see <http://www.gnu.org/licenses/>.
#
# Module
#     ionic_model_catalog
#
# Description
#     Manages metadata and parameter constraints for ionic ODE models.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

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
    """Suggested ODE timestep for stable integration. Model-specific:
    stiff full ionic models (TNNP, ORd) need ~1e-5; phenomenological models
    (AlievPanfilov, BuenoOrovio) tolerate larger steps."""

    recommended_stimulus_duration: float | None = 0.002
    """Default pulse duration for **volumetric (PDE) external stimulus**
    (electroProperties.<solver>Coeffs.externalStimulus.stimulusDuration).
    NOT a single-cell default — single-cell stim_duration is protocol-
    specific and lives in dict_entries.py with its own typical_value."""

    recommended_stimulus_intensity: float | None = 80000.0
    """Default current density for **volumetric (PDE) external stimulus**
    (externalStimulus.stimulusIntensity) in A/m³. NOT a single-cell default
    — single-cell stim_amplitude varies sharply by model class (~60 pA for
    full ionic models, ~0.4 for phenomenological) and lives in
    dict_entries.py with model-class guidance in its notes field."""

    notes: str = ""
    """Additional notes or warnings."""

    supports_heterogeneity: bool = False
    """Whether the model implements transmural ionicHeterogeneity (endo/M/epi
    blend) via configureIonicHeterogeneity. Batched variants inherit this flag
    through dataclasses.replace()."""

    supports_apex_base_heterogeneity: bool = False
    """Whether the model implements apexBaseBands longitudinal heterogeneity
    via configureApexBaseBandsHeterogeneity. Batched variants inherit this flag
    through dataclasses.replace()."""


# SOLVER_COMPATIBILITY_RULES moved to openfoam_driver/solver_coupling.py;
# re-exported here for backward compatibility with consumers that imported
# it from this module. Prefer the new home for new imports.
from .solver_coupling import SOLVER_COMPATIBILITY_RULES  # noqa: E402, F401


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
        supports_heterogeneity=True,
        supports_apex_base_heterogeneity=True,
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
        algebraic=("AV_ICaL", "AV_IK1", "AV_IKr", "AV_IKs", "AV_IKur", "AV_INa", "AV_INaCa", "AV_INaK", "AV_IbCa", "AV_IbNa", "AV_IpCa", "AV_Ito", "AV_ina_m_inf", "AV_ina_m_tau", "AV_ina_h_inf", "AV_ina_h_tau", "AV_ina_j_inf", "AV_ina_j_tau", "AV_ical_d_inf", "AV_ical_d_tau", "AV_ical_f_inf", "AV_ical_f_tau", "AV_ito_oa_inf", "AV_ito_oa_tau", "AV_ito_oi_inf", "AV_ito_oi_tau", "AV_ikur_ua_inf", "AV_ikur_ua_tau", "AV_ikur_ui_inf", "AV_ikur_ui_tau", "AV_gKur", "AV_ikr_xr_inf", "AV_ikr_xr_tau", "AV_iks_xs_inf", "AV_iks_xs_tau", "AV_fNaK", "AV_cajsr_w_inf", "AV_cajsr_w_tau", "AV_ical_fCa_inf", "AV_cajsr_u_inf", "Iion_cm", "Istim"),
        constants=("AC_CMDN_max", "AC_CSQN_max", "AC_Ca_up_max", "AC_Cao", "AC_Cm", "AC_ECaL", "AC_F", "AC_FRT", "AC_INaCa_max", "AC_INaK_max", "AC_I_diff", "AC_I_up_max", "AC_IpCa_max", "AC_KQ10", "AC_K_rel", "AC_K_up", "AC_KmCa", "AC_KmKo", "AC_KmNa", "AC_KmNai", "AC_Km_CMDN", "AC_Km_CSQN", "AC_Km_TRPN", "AC_Ko", "AC_Nao", "AC_R", "AC_RTF", "AC_T", "AC_TRPN_max", "AC_V_cell", "AC_V_i", "AC_V_rel", "AC_V_up", "AC_c1", "AC_c2", "AC_cajsr_u_tau", "AC_g", "AC_gCaL", "AC_gK1", "AC_gKr", "AC_gKs", "AC_gKur_base", "AC_gNa", "AC_gbCa", "AC_gbNa", "AC_gto", "AC_ical_fCa_tau", "AC_ksat", "AC_sigma", "AC_tau_tr"),
        recommended_exports=("membrane_V", "calcium_Cai", "AV_INa", "AV_ICaL", "AV_IKr", "AV_IKs", "AV_Ito", "Iion_cm"),
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
        recommended_exports=("membrane_V", "Ca_i", "AV_i_Na_i_Na", "AV_i_CaL_i_CaL", "AV_i_Kr_i_Kr", "AV_i_f_i_f", "AV_i_to_i_to", "Iion_cm"),
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
        algebraic=("AV_Jdiff", "AV_JdiffK", "AV_JdiffNa", "AV_EK", "AV_ENa", "AV_vffrt", "AV_vfrt", "AV_CaMKb", "AV_EKs", "AV_CaMKa", "AV_diff_CaMKt", "AV_CaMK_f", "AV_aa_h", "AV_aa_j", "AV_aa_m", "AV_bb_h", "AV_bb_j", "AV_bb_m", "AV_h_inf", "AV_m_inf", "AV_INa", "AV_j_inf", "AV_tau_h", "AV_tau_j", "AV_tau_m", "AV_aml", "AV_bml", "AV_hl_inf", "AV_INaL_INaL", "AV_ml_inf", "AV_tau_ml", "AV_PhiCaK", "AV_PhiCaL", "AV_PhiCaNa", "AV_f", "AV_tau_d", "AV_tau_fca", "AV_tau_ff", "AV_tau_fs", "AV_ICaL_ICaL", "AV_d_inf", "AV_f_inf", "AV_ICaK", "AV_ICaNa", "AV_fca_inf", "AV_ff_inf", "AV_fs_inf", "AV_rkr", "AV_tau_xr", "AV_xr_inf", "AV_IKr_IKr", "AV_KsCa", "AV_tau_xs1", "AV_tau_xs2", "AV_xs1_inf", "AV_IKs_IKs", "AV_xs2_inf", "AV_rk1", "AV_IK1_IK1", "AV_alpha_aa", "AV_beta_aa", "AV_rito2", "AV_diff_A", "AV_trel1factor", "AV_trel2factor", "AV_Jgap", "AV_diff_tjsrol", "AV_greljsrol", "AV_tau_Jrel1", "AV_tau_Jrel2", "AV_Jrelol", "AV_Jrel", "AV_allo", "AV_h4", "AV_hca", "AV_hna", "AV_h1", "AV_h5", "AV_h6", "AV_h7", "AV_h2", "AV_h3", "AV_h8", "AV_h9", "AV_k6", "AV_k3p", "AV_k3pp", "AV_k4p", "AV_k4pp", "AV_k7", "AV_k8", "AV_k3", "AV_k4", "AV_x1", "AV_x2", "AV_x3", "AV_x4", "AV_E1", "AV_E2", "AV_E3", "AV_E4", "AV_JncxCa", "AV_JncxNa", "AV_INaCa_i_INaCa_i", "AV_allo1", "AV_h111", "AV_h41", "AV_h71", "AV_h21", "AV_h31", "AV_h51", "AV_h61", "AV_h81", "AV_h91", "AV_k3p1", "AV_k3pp1", "AV_k4p1", "AV_k4pp1", "AV_k61", "AV_k71", "AV_k81", "AV_k31", "AV_k41", "AV_x11", "AV_x21", "AV_x31", "AV_x41", "AV_E11", "AV_E21", "AV_E31", "AV_E41", "AV_JncxCa1", "AV_JncxNa1", "AV_INaCa_ss_INaCa_ss", "AV_INaCa", "AV_Knai", "AV_Knao", "AV_P", "AV_a1", "AV_a3", "AV_b2", "AV_b3", "AV_b4", "AV_x12", "AV_x22", "AV_x32", "AV_x42", "AV_E12", "AV_E22", "AV_E32", "AV_E42", "AV_JnakK", "AV_JnakNa", "AV_INaK_INaK", "AV_xkb", "AV_IKb_IKb", "AV_INab_INab", "AV_ICab_ICab", "AV_IpCa_IpCa", "AV_Jleak", "AV_Jtr", "AV_Jtr2", "AV_Jupnp", "AV_Jupnp2", "AV_Jupp", "AV_Jupp2", "AV_fJupp", "AV_Jup", "AV_Jup2", "AV_Bcacsr", "AV_Bcai", "AV_Bcai2", "AV_Bcajsr", "AV_Bcass", "AV_cai_avg", "AV_diff_cansr", "AV_diff_ki", "AV_diff_kss", "AV_diff_nai", "AV_diff_nass", "AV_diff_cacsr", "AV_diff_cai", "AV_diff_cai2", "AV_diff_cajsr", "AV_diff_cass", "AV_kito2", "AV_Rel1", "AV_Rel2", "AV_ITo_ITo", "AV_Jrel1_inf", "AV_Jrel2_inf", "Istim", "Iion_cm"),
        constants=("AC_F", "AC_L", "AC_R", "AC_T", "AC_cao", "AC_cli", "AC_clo", "AC_ko", "AC_nao", "AC_pi", "AC_rad", "AC_vmyo1frac", "AC_Ageo", "AC_vcell", "AC_Acap", "AC_vjsr", "AC_vmyo", "AC_vnsr", "AC_vss", "AC_vcsr", "AC_vmyo1", "AC_vnsr1", "AC_vnsr2", "AC_vmyo2", "AC_CaMKo", "AC_ECl", "AC_KmCaM", "AC_KmCaMK", "AC_PKNa", "AC_aCaMK", "AC_bCaMK", "AC_GNa", "AC_GNaL", "AC_tau_hl", "AC_PCa", "AC_vhalf_d", "AC_vhalff", "AC_zca", "AC_PCaK", "AC_PCaNa", "AC_GKr", "AC_GKs", "AC_GK1", "AC_Gto", "AC_SOICR", "AC_grelbarjsrol", "AC_tau_gap", "AC_tauoff", "AC_tauon", "AC_Gncx", "AC_KmCaAct", "AC_kasymm", "AC_kcaoff", "AC_kcaon", "AC_kna1", "AC_kna2", "AC_kna3", "AC_qca", "AC_qna", "AC_wca", "AC_wna", "AC_wnaca", "AC_zna", "AC_h10", "AC_k2", "AC_k5", "AC_h11", "AC_h12", "AC_k1", "AC_h101", "AC_k21", "AC_k51", "AC_h1111", "AC_h121", "AC_k11", "AC_H", "AC_Khp", "AC_Kki", "AC_Kko", "AC_Kmgatp", "AC_Knai0", "AC_Knao0", "AC_Knap", "AC_Kxkur", "AC_MgADP", "AC_MgATP", "AC_Pnak", "AC_delta", "AC_eP", "AC_k1m", "AC_k1p", "AC_k2m", "AC_k2p", "AC_k3m", "AC_k3p2", "AC_k4m", "AC_k4p2", "AC_zk", "AC_a2", "AC_a4", "AC_b1", "AC_GKb", "AC_PNab", "AC_PCab", "AC_GpCa", "AC_BSLmax", "AC_BSRmax", "AC_KmBSL", "AC_KmBSR", "AC_cmdnmax", "AC_csqnmax", "AC_kmcmdn", "AC_kmcsqn", "AC_kmtrpn", "AC_trpnmax"),
        recommended_exports=("cell_v", "cai", "AV_INa", "AV_ICaL_ICaL", "AV_IKr_IKr", "AV_IKs_IKs", "AV_ITo_ITo", "Iion_cm"),
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
        states=("membrane_V", "Ca_i", "Ca_jn", "Ca_sl", "Ca_sr", "Csqn", "Ical_d", "Ical_f", "Ical_fCaB_jn", "Ical_fCaB_sl", "Ikr_xr", "Iks_xs", "Ikur_ikur_r", "Ikur_s", "Ina_h", "Ina_j", "Ina_m", "Inal_hl", "Inal_ml", "Ito_x", "Ito_y", "K_i", "ryr_i", "ryr_o", "ryr_ryr_r", "Na_i", "Na_jn", "Na_sl", "buffCa_CaM", "buffCa_Myoc", "buffCa_Myom", "buffCa_SLH_jn", "buffCa_SLH_sl", "buffCa_SLL_jn", "buffCa_SLL_sl", "buffCa_SRB", "buffCa_TnCHc", "buffCa_TnCHm", "buffCa_TnCL", "buffNa_NaB_jn", "buffNa_NaB_sl"),
        algebraic=("AV_ikr_xr_inf", "AV_ikr_xr_tau", "AV_ikur_r_inf", "AV_ikur_r_tau", "AV_ikur_s_inf", "AV_ikur_s_tau", "AV_ina_h_a", "AV_ina_h_b", "AV_ina_h_inf", "AV_ina_h_tau", "AV_ina_j_a", "AV_ina_j_b", "AV_ina_j_inf", "AV_ina_j_tau", "AV_ina_m_inf", "AV_ina_m_tau", "AV_inal_hl_inf", "AV_inal_ml_a", "AV_inal_ml_b", "AV_inal_ml_tau", "AV_inal_ml_inf", "AV_ical_fCaB_jn_tau", "AV_ical_fCaB_jn_inf", "AV_ical_fCaB_sl_tau", "AV_ical_fCaB_sl_inf", "AV_ito_x_inf", "AV_ito_x_tau", "AV_ito_y_inf", "AV_ito_y_tau", "AV_J_CaB_cytosol", "AV_J_CaB_jntion", "AV_J_CaB_sl", "AV_ical_d_inf", "AV_ical_d_tau", "AV_ical_f_inf", "AV_ical_f_tau", "AV_iks_xs_inf", "AV_iks_xs_tau", "AV_fnak", "AV_INaK_jn", "AV_INaK_sl", "AV_INaK", "AV_ECa_jn", "AV_ECa_sl", "AV_EK", "AV_ENa_jn", "AV_ENa_sl", "AV_RI", "AV_J_SRCarel", "AV_J_SR_leak", "AV_J_serca", "AV_kCaSR", "AV_kiSRCa", "AV_koSRCa", "AV_ICaB_jn", "AV_ICaB_sl", "AV_ICaB", "AV_ibarca_jn", "AV_ibarca_sl", "AV_ibark", "AV_ibarna_jn", "AV_ibarna_sl", "AV_ICaL_Ca_jn", "AV_ICaL_Ca_sl", "AV_ICaL_K", "AV_ICaL_Na_jn", "AV_ICaL_Na_sl", "AV_ICaL_Ca", "AV_ICaL_Na", "AV_ICaL", "AV_IClB", "AV_IClCa_jn", "AV_IClCa_sl", "AV_IClCa", "AV_ik1_inf_a", "AV_ik1_inf_b", "AV_kp", "AV_IKp_jn", "AV_IKp_sl", "AV_IKp", "AV_rr", "AV_IKr", "AV_EKs", "AV_IKs_jn", "AV_IKs_sl", "AV_IKs", "AV_IKur", "AV_INa_jn", "AV_INa_sl", "AV_INa", "AV_INaB_jn", "AV_INaB_sl", "AV_INaB", "AV_INaL_jn", "AV_INaL_sl", "AV_INaL", "AV_Ito", "AV_ik1_inf", "AV_IK1", "AV_IK_tot", "AV_Ka_jn", "AV_Ka_sl", "AV_s1_jn", "AV_s2_jn", "AV_s3_jn", "AV_s1_sl", "AV_s2_sl", "AV_s3_sl", "AV_ipca_IpCa_jn_a", "AV_ipca_IpCa_sl_a", "AV_ICl_tot", "AV_INaCa_jn", "AV_INaCa_sl", "AV_IpCa_jn", "AV_IpCa_sl", "AV_ICa_tot_jn", "AV_ICa_tot_sl", "AV_INaCa", "AV_IpCa", "AV_INa_tot_jn", "AV_INa_tot_sl", "AV_ICa_tot", "AV_INa_tot", "Istim", "Iion_cm"),
        constants=("AC_Bmax_Na_jn", "AC_Bmax_Na_sl", "AC_koff_na", "AC_kon_na", "AC_AF", "AC_C", "AC_ISO", "AC_RA", "AC_JCa_jnsl", "AC_JCa_slmyo", "AC_JNa_jnsl", "AC_JNa_slmyo", "AC_cell_length", "AC_cell_radius", "AC_pi", "AC_Vcell", "AC_Vjn", "AC_Vmyo", "AC_Vsl", "AC_Vsr", "AC_inal_hl_tau", "AC_Ca_o", "AC_Cl_i", "AC_Cl_o", "AC_K_o", "AC_Mg_i", "AC_Na_o", "AC_Fjn", "AC_Fjn_CaL", "AC_Fsl", "AC_Fsl_CaL", "AC_F", "AC_R", "AC_T", "AC_FRT", "AC_Q", "AC_Bmax_CaM", "AC_Bmax_SLhighj", "AC_Bmax_SLhighsl", "AC_Bmax_SLlowj", "AC_Bmax_SLlowsl", "AC_Bmax_SR", "AC_Bmax_TnChigh", "AC_Bmax_TnClow", "AC_Bmax_myosin", "AC_koff_cam", "AC_koff_myoca", "AC_koff_myomg", "AC_koff_slh", "AC_koff_sll", "AC_koff_sr", "AC_koff_tnchca", "AC_koff_tnchmg", "AC_koff_tncl", "AC_kon_cam", "AC_kon_myoca", "AC_kon_myomg", "AC_kon_slh", "AC_kon_sll", "AC_kon_sr", "AC_kon_tnchca", "AC_kon_tnchmg", "AC_kon_tncl", "AC_Bmax_Csqn", "AC_koff_csqn", "AC_kon_csqn", "AC_IbarNaK", "AC_KmKo", "AC_KmNaip", "AC_sigma", "AC_ECl", "AC_J_SR_leak_max", "AC_Kmf", "AC_Kmr", "AC_MaxSR", "AC_MinSR", "AC_Q10SRCaP", "AC_Vmax_SRCaP", "AC_ec50SR", "AC_hillSRCaP", "AC_kiCa", "AC_kim", "AC_koCa", "AC_kom", "AC_ks", "AC_amplitude", "AC_stim_duration", "AC_stim_offset", "AC_stim_period", "AC_gCaB", "AC_Q10CaL", "AC_f_conducting", "AC_fcaCaMSL", "AC_fcaCaj", "AC_pCa_max", "AC_pK_max", "AC_pNa_max", "AC_pCa", "AC_pK", "AC_pNa", "AC_gClB", "AC_GClCa", "AC_KdClCa", "AC_gKp", "AC_gKr_max", "AC_gKr", "AC_gKs_max", "AC_pNaK", "AC_gKs_jn", "AC_gKs_sl", "AC_gKur_max", "AC_gKur", "AC_gNa_max", "AC_gNa", "AC_gNaB", "AC_gNaL_max", "AC_gNaL", "AC_gto_max", "AC_gto", "AC_gK1_max", "AC_gK1", "AC_IbarNCX_max", "AC_Kdact", "AC_KmCai", "AC_KmCao", "AC_KmNai", "AC_KmNao", "AC_Q10NCX", "AC_ksat", "AC_nu", "AC_IbarNCX", "AC_IbarSLCaP", "AC_KmPCa", "AC_Q10SLCaP", "AC_ipca_IpCa_jn_b", "AC_ipca_IpCa_sl_b"),
        recommended_exports=("V", "Cass", "AV_INa", "AV_ICaL", "AV_IKr", "AV_IKs", "AV_Ito", "Iion_cm"),
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
    "Stewart": IonicModelEntry(
        states=("membrane_V", "Ihyperpolarization_activated_current_y_gate_y", "Irapid_time_dependent_potassium_current_Xr1_gate_Xr1", "Irapid_time_dependent_potassium_current_Xr2_gate_Xr2", "Islow_time_dependent_potassium_current_Xs_gate_Xs", "Ifast_sodium_current_m_gate_m", "Ifast_sodium_current_h_gate_h", "Ifast_sodium_current_j_gate_j", "IL_type_Ca_current_d_gate_d", "IL_type_Ca_current_f_gate_f", "IL_type_Ca_current_f2_gate_f2", "IL_type_Ca_current_fCass_gate_fCass", "Itransient_outward_current_s_gate_s", "Itransient_outward_current_r_gate_r", "Ca_i", "Ca_SR", "Ca_ss", "R_prime", "Na_i", "K_i"),
        algebraic=("AV_alpha_d", "AV_beta_d", "AV_d_inf", "AV_gamma_d", "AV_tau_d", "AV_f2_inf", "AV_tau_f2", "AV_fCass_inf", "AV_tau_fCass", "AV_f_inf", "AV_tau_f", "AV_i_p_Ca", "AV_alpha_h", "AV_beta_h", "AV_h_inf", "AV_tau_h", "AV_alpha_j", "AV_beta_j", "AV_j_inf", "AV_tau_j", "AV_alpha_m", "AV_beta_m", "AV_m_inf", "AV_tau_m", "AV_alpha_y", "AV_beta_y", "AV_y_inf", "AV_tau_y", "AV_alpha_xr1", "AV_beta_xr1", "AV_xr1_inf", "AV_tau_xr1", "AV_alpha_xr2", "AV_beta_xr2", "AV_xr2_inf", "AV_tau_xr2", "AV_alpha_xs", "AV_beta_xs", "AV_xs_inf", "AV_tau_xs", "AV_r_inf", "AV_tau_r", "AV_s_inf", "AV_tau_s", "AV_Ca_i_bufc", "AV_Ca_sr_bufsr", "AV_Ca_ss_bufss", "AV_i_leak", "AV_i_up", "AV_i_xfer", "AV_kcasr", "AV_k1", "AV_k2", "AV_O", "AV_i_rel", "AV_xK1_inf", "AV_E_Ca", "AV_E_K", "AV_i_NaK", "AV_a", "AV_i_sus", "AV_i_to", "AV_i_CaL", "AV_i_b_Ca", "AV_i_f_K", "AV_i_K1", "AV_i_p_K", "AV_i_Kr", "AV_E_Ks", "AV_E_Na", "AV_i_NaCa", "AV_i_Na", "AV_i_f_Na", "AV_i_Ks", "AV_i_b_Na", "AV_i_f", "Istim", "Iion_cm"),
        constants=("AC_K_pCa", "AC_g_pCa", "AC_g_CaL", "AC_g_bca", "AC_Buf_c", "AC_Buf_sr", "AC_Buf_ss", "AC_Ca_o", "AC_EC", "AC_K_buf_c", "AC_K_buf_sr", "AC_K_buf_ss", "AC_K_up", "AC_V_leak", "AC_V_rel", "AC_V_sr", "AC_V_ss", "AC_V_xfer", "AC_Vmax_up", "AC_k1_prime", "AC_k2_prime", "AC_k3", "AC_k4", "AC_max_sr", "AC_min_sr", "AC_g_Na", "AC_g_f_K", "AC_g_f_Na", "AC_g_K1", "AC_Cm", "AC_F", "AC_R", "AC_T", "AC_V_c", "AC_K_o", "AC_g_pK", "AC_g_Kr", "AC_P_kna", "AC_g_Ks", "AC_g_bna", "AC_K_NaCa", "AC_K_sat", "AC_Km_Ca", "AC_Km_Nai", "AC_alpha", "AC_sodium_calcium_exchanger_current_gamma", "AC_Na_o", "AC_K_mNa", "AC_K_mk", "AC_P_NaK", "AC_g_sus", "AC_g_to"),
        recommended_exports=("V", "Cai", "AV_i_Na", "AV_i_CaL", "AV_i_Kr", "AV_i_f", "AV_i_to", "Iion_cm"),
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
        states=("V", "K_i", "Na_i", "Ca_i", "Xr1", "Xr2", "Xs", "m", "h", "j", "d", "f", "fCa", "s", "r", "Ca_SR", "g"),
        algebraic=("Istim", "xr1_inf", "xr2_inf", "xs_inf", "m_inf", "h_inf", "j_inf", "d_inf", "f_inf", "alpha_fCa", "s_inf", "r_inf", "g_inf", "E_Na", "alpha_xr1", "alpha_xr2", "alpha_xs", "alpha_m", "alpha_h", "alpha_j", "alpha_d", "tau_f", "beta_fCa", "tau_s", "tau_r", "d_g", "E_K", "beta_xr1", "beta_xr2", "beta_xs", "beta_m", "beta_h", "beta_j", "gama_fCa", "beta_d", "E_Ks", "tau_xr1", "tau_xr2", "tau_xs", "tau_m", "tau_h", "tau_j", "gamma_d", "fCa_inf", "E_Ca", "tau_d", "d_fCa", "alpha_K1", "beta_K1", "xK1_inf", "i_K1", "i_Kr", "i_Ks", "i_Na", "i_b_Na", "i_CaL", "i_b_Ca", "i_to", "i_NaK", "i_NaCa", "i_p_Ca", "i_p_K", "i_rel", "i_up", "i_leak", "ddt_Ca_i_total", "ddt_Ca_sr_total", "f_JCa_i_free", "f_JCa_sr_free", "Iion_cm"),
        constants=("R", "T", "F", "Cm", "V_c", "P_kna", "K_o", "Na_o", "Ca_o", "g_K1", "g_Kr", "g_Ks", "g_Na", "g_bNa", "g_CaL", "g_bCa", "g_to", "P_NaK", "K_mk", "K_mNa", "K_NaCa", "K_sat", "alpha_NCX", "gamma_NCX", "Km_Ca", "Km_Nai", "g_pCa", "K_pCa", "g_pK", "tau_g", "a_rel", "b_rel", "c_rel", "K_up", "V_leak", "Vmax_up", "Buf_c", "K_buf_c", "Buf_sr", "K_buf_sr", "V_sr", "tau_fCa", "sInfPrefactor", "sInfShift", "sInfScale", "tauSGaussAmp", "tauSGaussShift", "tauSGaussWidth", "tauSSigmoidAmp", "tauSSigmoidShift", "tauSSigmoidScale", "tauSOffset"),
        recommended_exports=("membrane_V", "calcium_Cai", "i_Na", "i_CaL", "i_Kr", "i_Ks", "i_to", "Iion_cm"),
        compatible_tissues=("epicardialCells", "mCells", "endocardialCells"),
        supports_heterogeneity=True,
        supports_apex_base_heterogeneity=True,
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
        recommended_exports=("V", "Cass", "AV_INa_INa", "AV_ICaL_ICaL", "AV_IKr_IKr", "AV_IKs_IKs", "AV_Ito_Ito", "Iion_cm"),
        compatible_tissues=("epicardialCells", "mCells", "endocardialCells"),
        supports_heterogeneity=True,
        supports_apex_base_heterogeneity=True,
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
        states=("membrane_v", "CaMKt", "Nai", "Nasl", "Nass", "Ki", "Kss", "Ksl", "Cai", "Cass", "Casl", "Cansr", "Cajsr", "Cacsr", "INa_m", "INa_hf", "INa_hs", "INa_j", "INa_hsp", "INa_jp", "INaL_mL", "INaL_hL", "INaL_hLp", "Ito_a", "Ito_i1", "Ito_i2", "ICaL_d", "ICaL_ff", "ICaL_fs", "ICaL_fcaf", "ICaL_fcas", "ICaL_jca", "ICaL_ffp", "ICaL_fcafp", "ICaL_nca", "ICaT_b", "ICaT_g", "IKr_xrf", "IKr_xrs", "IKs_xs1", "IKs_xs2", "If_y", "IK1_xk1", "ryr_Jrel1", "ryr_Jrel2", "IP3_u"),
        algebraic=("AV_CaMKb", "AV_CaMKa", "AV_POip3", "AV_Jip3", "AV_IpCa_IpCa", "AV_Jdiff", "AV_JdiffK", "AV_JdiffNa", "AV_Jgap", "AV_JgapK", "AV_JgapNa", "AV_Jtr1", "AV_Jtr2", "AV_allo_i", "AV_allo_ss", "AV_h4_i", "AV_h4_ss", "AV_hca", "AV_hna", "AV_h1_i", "AV_h1_ss", "AV_h5_i", "AV_h5_ss", "AV_h6_i", "AV_h6_ss", "AV_h7_i", "AV_h7_ss", "AV_h2_i", "AV_h2_ss", "AV_h3_i", "AV_h3_ss", "AV_h8_i", "AV_h8_ss", "AV_h9_i", "AV_h9_ss", "AV_k6_i", "AV_k6_ss", "AV_k3p_i", "AV_k3p_ss", "AV_k3pp_i", "AV_k3pp_ss", "AV_k4p_i", "AV_k4p_ss", "AV_k4pp_i", "AV_k4pp_ss", "AV_k7_i", "AV_k7_ss", "AV_k8_i", "AV_k8_ss", "AV_k3_i", "AV_k3_ss", "AV_k4_i", "AV_k4_ss", "AV_x1_i", "AV_x1_ss", "AV_x2_i", "AV_x2_ss", "AV_x3_i", "AV_x3_ss", "AV_x4_i", "AV_x4_ss", "AV_E1_i", "AV_E1_ss", "AV_E2_i", "AV_E2_ss", "AV_E3_i", "AV_E3_ss", "AV_E4_i", "AV_E4_ss", "AV_JncxCa_i", "AV_JncxCa_ss", "AV_JncxNa_i", "AV_JncxNa_ss", "AV_INaCa_i_INaCa_i", "AV_INaCa_ss", "AV_Knai", "AV_Knao", "AV_P", "AV_a1", "AV_a3", "AV_b2", "AV_b3", "AV_b4", "AV_x1", "AV_x2", "AV_x3", "AV_x4", "AV_E1", "AV_E2", "AV_E3", "AV_E4", "AV_JnakK", "AV_JnakNa", "AV_INaK_INaK", "AV_dkmplb", "AV_dqupcamk", "AV_Jup1", "AV_Jup2", "AV_ECa", "AV_EK", "AV_ENa", "AV_EKs", "AV_bss", "AV_gss", "AV_taub", "AV_taug", "AV_ICaT_ICaT", "AV_rk1", "AV_txk1", "AV_xk1ss", "AV_IK1_IK1", "AV_Axrf", "AV_rkr", "AV_txrf", "AV_txrs", "AV_xrss", "AV_Axrs", "AV_xr", "AV_IKr_IKr", "AV_KsCa", "AV_txs1", "AV_txs2", "AV_xs1ss", "AV_IKs_IKs", "AV_xs2ss", "AV_fINap", "AV_hssp", "AV_thf", "AV_ths", "AV_tj", "AV_hss", "AV_mss", "AV_thsp", "AV_tjp", "AV_tm", "AV_h", "AV_hp", "AV_jss", "AV_INa_INa", "AV_tauy", "AV_yss", "AV_IfK", "AV_IfNa", "AV_If_If", "AV_asus", "AV_Isus_Isus", "AV_ass", "AV_iss", "AV_taua", "AV_tauif", "AV_tauis", "AV_Ito_Ito", "AV_fINaLp", "AV_hLss", "AV_hLssp", "AV_mLss", "AV_tmL", "AV_INaL_INaL", "AV_Afcaf", "AV_dss", "AV_fICaLp", "AV_fss", "AV_km2n", "AV_td", "AV_tfcaf", "AV_tfcas", "AV_tff", "AV_tfs", "AV_Afcas", "AV_anca", "AV_fcass", "AV_tfcafp", "AV_tffp", "AV_f", "AV_fca", "AV_fcap", "AV_fp", "AV_Bcacsr", "AV_Bcai", "AV_Bcajsr", "AV_Bcasl", "AV_Bcass", "AV_vffrt", "AV_vfrt", "AV_REL2", "AV_ireltau", "AV_ireltau2", "AV_irelss2", "AV_PhiCaK", "AV_PhiCaL", "AV_PhiCaNa", "AV_ICab_ICab", "AV_INab_INab", "AV_ICaK", "AV_ICaL_ICaL", "AV_ICaNa", "AV_REL", "AV_irelss", "Istim", "Iion_cm"),
        constants=("AC_CaMKo", "AC_KmCaM", "AC_KmCaMK", "AC_aCaMK", "AC_bCaMK", "AC_IP3_IP3", "AC_k0", "AC_k0a", "AC_k1", "AC_k1a", "AC_k2", "AC_k2a", "AC_tauip3", "AC_GpCa", "AC_KmCap", "AC_L", "AC_greekpi", "AC_rad", "AC_Ageo", "AC_vcell", "AC_Acap", "AC_vcsr", "AC_vjsr", "AC_vmyo", "AC_vnsr", "AC_vsl", "AC_vss", "AC_gaptau", "AC_sstau", "AC_cao", "AC_ko", "AC_nao", "AC_F", "AC_R", "AC_T", "AC_zca", "AC_zk", "AC_zna", "AC_Gncx", "AC_KmCaAct", "AC_kasymm", "AC_kcaoff", "AC_kcaon", "AC_kna1", "AC_kna2", "AC_kna3", "AC_qca", "AC_qna", "AC_wca", "AC_wna", "AC_wnaca", "AC_h10_i", "AC_h10_ss", "AC_k2_i", "AC_k2_ss", "AC_k5_i", "AC_k5_ss", "AC_h11_i", "AC_h11_ss", "AC_h12_i", "AC_h12_ss", "AC_k1_i", "AC_k1_ss", "AC_H", "AC_Khp", "AC_Kki", "AC_Kko", "AC_Kmgatp", "AC_Knai0", "AC_Knao0", "AC_Knap", "AC_Kxkur", "AC_MgADP", "AC_MgATP", "AC_Pnak", "AC_delta", "AC_eP", "AC_k1m", "AC_k1p", "AC_k2m", "AC_k2p", "AC_k3m", "AC_k3p", "AC_k4m", "AC_k4p", "AC_a2", "AC_a4", "AC_b1", "AC_dkmplbbar", "AC_dqupcamkbar", "AC_kmup", "AC_nsrbar", "AC_PKNa", "AC_GCaT", "AC_GK1", "AC_GKr", "AC_GKs", "AC_Ahf", "AC_GNa", "AC_hssV1", "AC_hssV2", "AC_mssV1", "AC_mssV2", "AC_mtD1", "AC_mtD2", "AC_mtV1", "AC_mtV2", "AC_mtV3", "AC_mtV4", "AC_Ahs", "AC_GfK", "AC_GfNa", "AC_Gsus", "AC_Gto", "AC_GNaL", "AC_thL", "AC_thLp", "AC_Aff", "AC_Kmn", "AC_PCa", "AC_k2n", "AC_tjca", "AC_Afs", "AC_PCaK", "AC_PCaNa", "AC_PCap", "AC_PCaKp", "AC_PCaNap", "AC_PCab", "AC_PNab", "AC_BSLmax", "AC_BSRmax", "AC_KmBSL", "AC_KmBSR", "AC_cm", "AC_cmdnmax", "AC_cmdnmaxsl", "AC_csqnmax", "AC_csqnmaxsl", "AC_kmcmdn", "AC_kmcsqn", "AC_kmtrpn", "AC_trpnmax", "AC_trpnmaxsl"),
        recommended_exports=("Vm", "Cass", "AV_INa_INa", "AV_ICaL_ICaL", "AV_IKr_IKr", "AV_IKs_IKs", "AV_Ito_Ito", "Iion_cm"),
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
        algebraic=("AV_CaMK_Phos_ss_ICaL", "AV_CaMK_Phos_ss_PLB", "AV_CaMK_Phos_ss_RyR", "AV_CaMK_active", "AV_bound", "AV_bound_serca", "AV_casig_SERCA_act", "AV_Afcaf", "AV_Afcas", "AV_ICaK", "AV_ICaK_i_BP", "AV_ICaK_i_CaMK", "AV_ICaK_i_NP", "AV_ICaK_i_PKA", "AV_ICaK_junc", "AV_ICaK_sl", "AV_ICaK_ss_BP", "AV_ICaK_ss_CaMK", "AV_ICaK_ss_NP", "AV_ICaK_ss_PKA", "AV_ICaL_ICaL", "AV_ICaL_i_BP", "AV_ICaL_i_CaMK", "AV_ICaL_i_NP", "AV_ICaL_i_PKA", "AV_ICaL_junc", "AV_ICaL_sl", "AV_ICaL_ss_BP", "AV_ICaL_ss_CaMK", "AV_ICaL_ss_NP", "AV_ICaL_ss_PKA", "AV_ICaNa", "AV_ICaNa_i_BP", "AV_ICaNa_i_CaMK", "AV_ICaNa_i_NP", "AV_ICaNa_i_PKA", "AV_ICaNa_junc", "AV_ICaNa_sl", "AV_ICaNa_ss_BP", "AV_ICaNa_ss_CaMK", "AV_ICaNa_ss_NP", "AV_ICaNa_ss_PKA", "AV_ICa_tot", "AV_Ii", "AV_Iss", "AV_PhiCaK_i", "AV_PhiCaK_ss", "AV_PhiCaL_i", "AV_PhiCaL_ss", "AV_PhiCaNa_i", "AV_PhiCaNa_ss", "AV_anca", "AV_anca_i", "AV_dPss", "AV_dss", "AV_f", "AV_fBP", "AV_fBPss", "AV_fICaL_BP", "AV_fICaL_CaMKonly", "AV_fICaL_PKAonly", "AV_fICaLp", "AV_f_P", "AV_fca", "AV_fcaBP", "AV_fcaBPss", "AV_fcap", "AV_fcap_P", "AV_fcass", "AV_fcass_P", "AV_fp", "AV_fss", "AV_fss_P", "AV_gamma_cai", "AV_gamma_cass", "AV_gamma_ki", "AV_gamma_kss", "AV_gamma_nai", "AV_gamma_nass", "AV_jcass", "AV_km2n", "AV_sigmoidTransition", "AV_sigmoidTransition2", "AV_tauTransition", "AV_tauTransition2", "AV_td", "AV_tfcaf", "AV_tfcafp", "AV_tfcas", "AV_tff", "AV_tffp", "AV_tfs", "AV_ICab_ICab", "AV_ICab_junc", "AV_ICab_sl", "AV_IClCa", "AV_IClCa_junc", "AV_IClCa_sl", "AV_IClb", "AV_IK1_IK1", "AV_K1ss", "AV_aK1", "AV_bK1", "AV_IKb_IKb", "AV_xkb", "AV_IKr_IKr", "AV_IKr_alpha", "AV_alpha_2", "AV_alpha_C2ToI", "AV_alpha_i", "AV_IKr_beta", "AV_beta_2", "AV_beta_ItoC2", "AV_beta_i", "AV_IKs_IKs", "AV_IKs_junc", "AV_IKs_sl", "AV_VsTs_Ca_junc", "AV_VsTs_Ca_sl", "AV_VsXs_Ca_junc", "AV_VsXs_Ca_sl", "AV_gks_junc", "AV_gks_sl", "AV_tauxs_junc", "AV_tauxs_sl", "AV_xsss_junc", "AV_xsss_sl", "AV_INa_INa", "AV_INaBase", "AV_INaBase_BP", "AV_INaBase_CaMK", "AV_INaBase_NP", "AV_INaBase_PKA", "AV_INaj", "AV_INasl", "AV_ah", "AV_aj", "AV_bh", "AV_bj", "AV_fINa_BP", "AV_fINa_CaMKonly", "AV_fINa_PKAonly", "AV_fINap", "AV_hss", "AV_hss_P", "AV_hssp", "AV_hssp_P", "AV_jss", "AV_jss_P", "AV_jssp_P", "AV_mss", "AV_mss_P", "AV_th", "AV_tj", "AV_tjp", "AV_tm", "AV_E1_i", "AV_E1_ss", "AV_E2_i", "AV_E2_ss", "AV_E3_i", "AV_E3_ss", "AV_E4_i", "AV_E4_ss", "AV_INaCa_INaCa", "AV_INaCa_i", "AV_INaCa_ss", "AV_JncxCa_i", "AV_JncxCa_ss", "AV_JncxNa_i", "AV_JncxNa_ss", "AV_allo_i", "AV_allo_ss", "AV_h1_i", "AV_h1_ss", "AV_h2_i", "AV_h2_ss", "AV_h3_i", "AV_h3_ss", "AV_h4_i", "AV_h4_ss", "AV_h5_i", "AV_h5_ss", "AV_h6_i", "AV_h6_ss", "AV_h7_i", "AV_h7_ss", "AV_h8_i", "AV_h8_ss", "AV_h9_i", "AV_h9_ss", "AV_hca", "AV_hna", "AV_k3_i", "AV_k3_ss", "AV_k3p_i", "AV_k3p_ss", "AV_k3pp_i", "AV_k3pp_ss", "AV_k4_i", "AV_k4_ss", "AV_k4p_i", "AV_k4p_ss", "AV_k4pp_i", "AV_k4pp_ss", "AV_k6_i", "AV_k6_ss", "AV_k7_i", "AV_k7_ss", "AV_k8_i", "AV_k8_ss", "AV_x1_i", "AV_x1_ss", "AV_x2_i", "AV_x2_ss", "AV_x3_i", "AV_x3_ss", "AV_x4_i", "AV_x4_ss", "AV_INaK_INaK", "AV_INaK_junc_PKA", "AV_INaK_junc_noPKA", "AV_INaK_sl_PKA", "AV_INaK_sl_noPKA", "AV_INaKj", "AV_INaKsl", "AV_fnak", "AV_GNaL", "AV_INaL_INaL", "AV_INaLj", "AV_INaLsl", "AV_fINaLp", "AV_hLss", "AV_hLssp", "AV_mLss", "AV_tmL", "AV_INab_INab", "AV_INabj", "AV_INabsl", "AV_IpCa_IpCa", "AV_IpCa_junc", "AV_IpCa_sl", "AV_Ito_Ito", "AV_Itof", "AV_Itos", "AV_dti_develop", "AV_dti_recover", "AV_fItop", "AV_tauxtof", "AV_tauxtos", "AV_tauytof", "AV_tauytof_p", "AV_tauytos", "AV_tauytos_p", "AV_xtoss", "AV_xtoss_p", "AV_ytoss", "AV_Jserca", "AV_Jserca_np", "AV_Jserca_p", "AV_Vmax_mult", "AV_phosphorylationTotal", "AV_Ta", "AV_XU", "AV_d_contraction_Ca_TRPN", "AV_gamma_rate", "AV_gamma_rate_w", "AV_xb_su", "AV_xb_su_gamma", "AV_xb_uw", "AV_xb_ws", "AV_xb_wu", "AV_xb_wu_gamma", "AV_time", "AV_ICa_tot_junc", "AV_ICa_tot_sl", "AV_ICl_tot", "AV_IK_tot", "AV_INa_tot_junc", "AV_INa_tot_sl", "AV_J_CaB_cytosol", "AV_J_CaB_junction", "AV_J_CaB_sl", "AV_d_buffers_CaM", "AV_d_buffers_Csqn", "AV_d_buffers_Myosin_ca", "AV_d_buffers_NaBj", "AV_d_buffers_NaBsl", "AV_d_buffers_SLHj", "AV_d_buffers_SLHsl", "AV_d_buffers_SLLj", "AV_d_buffers_SLLsl", "AV_d_buffers_SRB", "AV_d_buffers_TnCHc", "AV_d_buffers_TnClow", "AV_vffrt", "AV_vfrt", "AV_ECaj", "AV_ECasl", "AV_ECl", "AV_EK", "AV_EKs", "AV_ENaj", "AV_ENasl", "AV_CaMKIILeakMultiplier", "AV_ICaL_junc_positive", "AV_ICaL_junc_sigmoided", "AV_J_SRCarel", "AV_J_SRCarel_np", "AV_J_SRCarel_p", "AV_J_SRleak", "AV_Jrel_ICaLdep", "AV_Jrel_inact_inf", "AV_Jrel_inact_inf2", "AV_Jrel_inf", "AV_RI_to_CI", "AV_RIcleft", "AV_RIcleft_p", "AV_kCaSR", "AV_kiSRCa", "AV_koSRCa", "AV_nonlinearModifier", "AV_sigmoidBaseCaI", "AV_tau_rel", "Istim", "Iion_cm"),
        constants=("AC_CaMK0", "AC_K_Phos_CaMK", "AC_Km_CaMK_Ca", "AC_PP1_tot", "AC_Whole_cell_PP1", "AC_CaMK_alpha", "AC_alpha_serca", "AC_CaMK_beta", "AC_tau_cal", "AC_tau_plb", "AC_tau_ryr", "AC_Aff", "AC_Afs", "AC_ICaL_fractionSS", "AC_Io", "AC_Kmn", "AC_PCa", "AC_PCaK", "AC_PCaK_P", "AC_PCaKp", "AC_PCaNa", "AC_PCaNa_P", "AC_PCaNap", "AC_PCa_P", "AC_PCa_Pb", "AC_PCa_b", "AC_PCap", "AC_constA", "AC_dielConstant", "AC_fICaLP", "AC_fICaL_P", "AC_gamma_cao", "AC_gamma_ko", "AC_gamma_nao", "AC_k2n", "AC_rateRecovery", "AC_tjca", "AC_GCab", "AC_GCab_b", "AC_GClCa", "AC_GClCa_b", "AC_GClb", "AC_GClb_b", "AC_KdClCa", "AC_GK1", "AC_GK1_b", "AC_IK1_celltype_factor", "AC_IK1_sex_factor", "AC_GKb", "AC_GKb_b", "AC_GKr", "AC_GKr_b", "AC_alpha_1", "AC_IKr_beta_1", "AC_IKr_celltype_factor", "AC_IKr_sex_factor", "AC_GKs", "AC_GKs_b", "AC_P_g_0", "AC_P_g_max", "AC_P_tau_0", "AC_P_tau_max", "AC_P_vh_0", "AC_P_vh_max", "AC_IKs_celltype_factor", "AC_gKs_factor", "AC_kPKA_IKs", "AC_IKs_sex_factor", "AC_GNa", "AC_GNa_P", "AC_GNa_b", "AC_fINa_P", "AC_Gncx", "AC_Gncx_b", "AC_INaCa_fractionSS", "AC_KmCaAct", "AC_INaCa_celltype_factor", "AC_h10_i", "AC_h10_ss", "AC_h11_i", "AC_h11_ss", "AC_h12_i", "AC_h12_ss", "AC_k1_i", "AC_k1_ss", "AC_k2_i", "AC_k2_ss", "AC_k5_i", "AC_k5_ss", "AC_kasymm", "AC_kcaoff", "AC_kcaon", "AC_kna1", "AC_kna2", "AC_kna3", "AC_qca", "AC_qna", "AC_INaCa_sex_factor", "AC_wca", "AC_wna", "AC_wnaca", "AC_IbarNaK", "AC_IbarNaK_b", "AC_KmKo", "AC_KmNaip", "AC_KmNaip_PKA", "AC_GNaL_b", "AC_thL", "AC_thLp", "AC_GNab", "AC_GNab_b", "AC_IbarSLCaP", "AC_IbarSLCaP_b", "AC_KmPCa", "AC_Q10SLCaP", "AC_Gto_fast", "AC_Gto_slow", "AC_fICaL_PKA", "AC_fIKs_PKA", "AC_fINaK_PKA", "AC_fINa_PKA", "AC_fMyBPC_PKA", "AC_fPLB_PKA", "AC_fTnI_PKA", "AC_Km_SERCA_Ca", "AC_Kmf", "AC_Kmf_p", "AC_Kmr", "AC_Max_Vmax_SERCA_Ca", "AC_Q10SRCaP", "AC_Vmax_SRCaP", "AC_Vmax_SRCaP_b", "AC_hillSRCaP", "AC_cellLength", "AC_cellRadius", "AC_vcell", "AC_vjunc", "AC_vmyo", "AC_vsl", "AC_vsr", "AC_A", "AC_Lfac", "AC_PKAForceMultiplier", "AC_TOT_A", "AC_TRPN_n", "AC_Tref", "AC_XSSS", "AC_XWSS", "AC_beta_0", "AC_contraction_beta_1", "AC_ca50", "AC_cds", "AC_cdw", "AC_dr", "AC_fPKA_TnI", "AC_fracTnIpo", "AC_contraction_gamma", "AC_gamma_wu", "AC_k_su", "AC_k_uw", "AC_k_ws", "AC_k_wu", "AC_koff", "AC_ktm_block", "AC_ktm_unblock", "AC_lambda", "AC_lambda0", "AC_lambda_max", "AC_lambda_min", "AC_lambda_rate", "AC_mu", "AC_nperm", "AC_nu", "AC_perm50", "AC_phi", "AC_wfrac", "AC_sex", "AC_cao", "AC_clo", "AC_ko", "AC_nao", "AC_Bmax_CaM", "AC_Bmax_Csqn", "AC_Bmax_Naj", "AC_Bmax_Nasl", "AC_Bmax_SLhighj", "AC_Bmax_SLhighsl", "AC_Bmax_SLlowj", "AC_Bmax_SLlowsl", "AC_Bmax_SR", "AC_Bmax_TnChigh", "AC_Bmax_TnClow", "AC_Bmax_myosin", "AC_J_ca_juncsl", "AC_J_ca_slmyo", "AC_J_na_juncsl", "AC_J_na_slmyo", "AC_koff_cam", "AC_koff_csqn", "AC_koff_myoca", "AC_koff_myomg", "AC_koff_na", "AC_koff_slh", "AC_koff_sll", "AC_koff_sr", "AC_koff_tnchca", "AC_koff_tnchmg", "AC_kon_cam", "AC_kon_csqn", "AC_kon_myoca", "AC_kon_myomg", "AC_kon_na", "AC_kon_slh", "AC_kon_sll", "AC_kon_sr", "AC_kon_tnchca", "AC_kon_tnchmg", "AC_mgi", "AC_Cmem", "AC_Fjunc", "AC_Fsl", "AC_ICaLPCa_multiplier", "AC_ICab_multiplier", "AC_IClCa_multiplier", "AC_IClb_multiplier", "AC_IK1_multiplier", "AC_IKb_multiplier", "AC_IKr_multiplier", "AC_IKs_multiplier", "AC_INaCa_multiplier", "AC_INaK_multiplier", "AC_INaL_multiplier", "AC_INa_multiplier", "AC_INab_multiplier", "AC_IpCa_multiplier", "AC_Itof_multiplier", "AC_Itos_multiplier", "AC_Jrel_multiplier", "AC_Jup_multiplier", "AC_F", "AC_Qpow", "AC_R", "AC_T", "AC_zca", "AC_zcl", "AC_zk", "AC_zna", "AC_PNaK", "AC_CI_to_RI", "AC_MaxSR", "AC_MinSR", "AC_a_rel", "AC_baseRateCaI", "AC_bt", "AC_caExpFactor", "AC_caExpFactor2", "AC_caTransFactor", "AC_caTransFactor2", "AC_caTransFactor2p", "AC_directRelMidpoint", "AC_ec50SR", "AC_ecCaI", "AC_kiCa", "AC_kim", "AC_koCa", "AC_kom", "AC_ks", "AC_maxCaI", "AC_minCaI", "AC_steepnessCaI", "AC_steepnessCaSR", "AC_tauInact", "AC_tauInact2", "gnalTissueScale"),
        recommended_exports=("v", "cai", "contraction_Ca_TRPN", "AV_INa_INa", "AV_ICaL_ICaL", "AV_IKr_IKr", "AV_IKs_IKs", "AV_Ito_Ito", "Iion_cm"),
        compatible_tissues=("epicardialCells", "mCells", "endocardialCells"),
        supports_heterogeneity=True,
        supports_apex_base_heterogeneity=True,
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
        algebraic=("AV_ICaL", "AV_IK1", "AV_IKr", "AV_IKs", "AV_IKur", "AV_INa", "AV_INaCa", "AV_INaK", "AV_IClCa", "AV_IbCa", "AV_IbNa", "AV_IbK", "AV_IpCa", "AV_ina_m_inf", "AV_ina_m_tau", "AV_ina_h_inf", "AV_ina_h_tau", "AV_ina_j_inf", "AV_ina_j_tau", "AV_ical_d_inf", "AV_ical_d_tau", "AV_ical_f_inf", "AV_ical_f_tau", "AV_ikur_ua_inf", "AV_ikur_ua_tau", "AV_ikur_uif_inf", "AV_ikur_uif_tau", "AV_ikur_uis_inf", "AV_ikur_uis_tau", "AV_gKur", "AV_ikr_xr_inf", "AV_ikr_xr_tau", "AV_iks_xs_inf", "AV_iks_xs_tau", "AV_iclca_qCa_inf", "AV_iclca_qCa_tau", "AV_fNaK", "AV_ryr_w_inf", "AV_ryr_w_tau", "AV_ical_fCa_inf", "AV_ryr_u_inf", "Iion_cm", "Istim"),
        constants=("AC_CMDN_max", "AC_CSQN_max", "AC_Ca_up_max", "AC_Cao", "AC_Cm", "AC_ECaL", "AC_F", "AC_FRT", "AC_INaCa_max", "AC_INaK_max", "AC_I_diff", "AC_I_up_max", "AC_KQ10", "AC_krel", "AC_K_up", "AC_KmCa", "AC_KmKo", "AC_KmNa", "AC_KmNai", "AC_Km_CMDN", "AC_Km_CSQN", "AC_Km_TRPN", "AC_Ko", "AC_Nao", "AC_Clo", "AC_R", "AC_RTF", "AC_T", "AC_TRPN_max", "AC_V_cell", "AC_V_i", "AC_V_rel", "AC_V_up", "AC_c1", "AC_c2", "AC_cajsr_u_tau", "AC_gNa", "AC_gK1", "AC_gKr", "AC_gKs", "AC_gKur_amp", "AC_gClCa", "AC_gbCa", "AC_gbNa", "AC_gbK", "AC_ical_fCa_tau", "AC_gamma", "AC_ksat", "AC_sigma", "AC_tau_tr", "AC_gCaL", "AC_IpCa_max"),
        recommended_exports=("membrane_V", "calcium_Cai", "AV_INa", "AV_ICaL", "AV_IKr", "AV_IKs", "AV_IClCa", "Iion_cm"),
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
        states=("V", "u1", "u2", "u3"),
        algebraic=("Iion",),
        constants=("Cm", "Beta", "Chi"),
        recommended_exports=("V", "u1", "u2", "u3"),
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
        states=("V", "u1", "u2", "u3"),
        algebraic=("Iion",),
        constants=("Cm", "Beta", "Chi"),
        recommended_exports=("V", "u1", "u2", "u3"),
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
        algebraic=("Iion",),
        constants=("Cm", "Beta", "Chi"),
        recommended_exports=("V", "u1", "u2", "u3"),
        compatible_tissues=("manufactured",),
        compatible_solvers=("bidomainSolver",),
        species=("generic",),
        cardiac_region=("manufactured",),
        model_type="manufactured",
        description="Manufactured solution for bidomain+bath (extracellular potential) solver convergence testing and verification. Requires bidomainSolverCoeffs.bathPotentialDomain and the manufacturedFDABathBidomainVerifier hook. Not for physiological studies.",
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


BATCHED_MODELS = [
    "AlievPanfilovcompactBatched",
    "BuenoOroviocompactBatched",
    "CourtemanchecompactBatched",
    "FabbricompactBatched",
    "GaurcompactBatched",
    "GrandicompactBatched",
    "PerisYaguecompactBatched",
    "StewartcompactBatched",
    "TNNPcompactBatched",
    "ToRORd_dynClcompactBatched",
    "TrovatocompactBatched",
    "TWorldcompactBatched",
]

for batched_name in BATCHED_MODELS:
    parent_name = batched_name.replace("compactBatched", "")
    if parent_name in IONIC_MODEL_CATALOG:
        parent = IONIC_MODEL_CATALOG[parent_name]
        from dataclasses import replace
        IONIC_MODEL_CATALOG[batched_name] = replace(
            parent,
            description=(
                f"GPU/CUDA variant of {parent_name}. "
                f"Identical ion dynamics — same states, algebraic variables, and "
                f"physiological parameters as {parent_name}. "
                f"Executed on a CUDA-capable GPU using compact Rush-Larsen batched ODE "
                f"integration. Select this instead of {parent_name} when a CUDA GPU is "
                f"available; use the non-batched {parent_name} for CPU-only runs. "
                f"Cannot be mixed with the non-batched variant in the same simulation."
            ),
            aliases=(f"{parent_name} GPU", f"{parent_name} CUDA", f"{parent_name} compact batched"),
            notes=(
                f"Requires a CUDA-capable GPU. "
                f"Physics are identical to {parent_name} — see that entry for "
                f"species, cardiac region, tissue compatibility, and recommended parameters. "
                f"The 'compactBatched' suffix is a hardware-selection signal only; "
                f"all electroProperties keys are the same as the non-batched variant."
            ),
        )
