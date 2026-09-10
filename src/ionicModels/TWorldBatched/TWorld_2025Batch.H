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
    State-major (SoA) batched evaluator for the TWorld 2025 ionic model.
Author
    Simao Nieto de Castro, UCD.
\*---------------------------------------------------------------------------*/

#pragma once

#include "TWorld_2025Names.H"
#include "StimulusProtocolPOD.H"
#include <cmath>

namespace Foam
{

enum TWORLD_BATCH_SUPPORT_INDEX
{
    // CaMK fractions
    TWORLD_BATCH_SUPPORT_tau_camk_ICaL,
    TWORLD_BATCH_SUPPORT_gInf_camk_ICaL,
    TWORLD_BATCH_SUPPORT_tau_camk_PLB,
    TWORLD_BATCH_SUPPORT_gInf_camk_PLB,
    TWORLD_BATCH_SUPPORT_tau_camk_RyR,
    TWORLD_BATCH_SUPPORT_gInf_camk_RyR,
    // INa
    TWORLD_BATCH_SUPPORT_tau_m,
    TWORLD_BATCH_SUPPORT_gInf_m,
    TWORLD_BATCH_SUPPORT_tau_m_P,
    TWORLD_BATCH_SUPPORT_gInf_m_P,
    TWORLD_BATCH_SUPPORT_tau_h,
    TWORLD_BATCH_SUPPORT_gInf_h,
    TWORLD_BATCH_SUPPORT_tau_h_P,
    TWORLD_BATCH_SUPPORT_gInf_h_P,
    TWORLD_BATCH_SUPPORT_tau_hp,
    TWORLD_BATCH_SUPPORT_gInf_hp,
    TWORLD_BATCH_SUPPORT_tau_hp_P,
    TWORLD_BATCH_SUPPORT_gInf_hp_P,
    TWORLD_BATCH_SUPPORT_tau_j,
    TWORLD_BATCH_SUPPORT_gInf_j,
    TWORLD_BATCH_SUPPORT_tau_j_P,
    TWORLD_BATCH_SUPPORT_gInf_j_P,
    TWORLD_BATCH_SUPPORT_tau_jp,
    TWORLD_BATCH_SUPPORT_gInf_jp,
    TWORLD_BATCH_SUPPORT_tau_jp_P,
    TWORLD_BATCH_SUPPORT_gInf_jp_P,
    // INaL
    TWORLD_BATCH_SUPPORT_tau_mL,
    TWORLD_BATCH_SUPPORT_gInf_mL,
    TWORLD_BATCH_SUPPORT_tau_hL,
    TWORLD_BATCH_SUPPORT_gInf_hL,
    TWORLD_BATCH_SUPPORT_tau_hLp,
    TWORLD_BATCH_SUPPORT_gInf_hLp,
    // ICaL standard
    TWORLD_BATCH_SUPPORT_tau_d,
    TWORLD_BATCH_SUPPORT_gInf_d,
    TWORLD_BATCH_SUPPORT_tau_ff,
    TWORLD_BATCH_SUPPORT_gInf_ff,
    TWORLD_BATCH_SUPPORT_tau_fs,
    TWORLD_BATCH_SUPPORT_gInf_fs,
    TWORLD_BATCH_SUPPORT_tau_fcaf,
    TWORLD_BATCH_SUPPORT_gInf_fcaf,
    TWORLD_BATCH_SUPPORT_tau_fcas,
    TWORLD_BATCH_SUPPORT_gInf_fcas,
    TWORLD_BATCH_SUPPORT_tau_jca,
    TWORLD_BATCH_SUPPORT_gInf_jca,
    TWORLD_BATCH_SUPPORT_tau_ffp,
    TWORLD_BATCH_SUPPORT_gInf_ffp,
    TWORLD_BATCH_SUPPORT_tau_fcafp,
    TWORLD_BATCH_SUPPORT_gInf_fcafp,
    // ICaL PKA/BP variants
    TWORLD_BATCH_SUPPORT_tau_d_P,
    TWORLD_BATCH_SUPPORT_gInf_d_P,
    TWORLD_BATCH_SUPPORT_tau_ff_P,
    TWORLD_BATCH_SUPPORT_gInf_ff_P,
    TWORLD_BATCH_SUPPORT_tau_fs_P,
    TWORLD_BATCH_SUPPORT_gInf_fs_P,
    TWORLD_BATCH_SUPPORT_tau_fcaf_P,
    TWORLD_BATCH_SUPPORT_gInf_fcaf_P,
    TWORLD_BATCH_SUPPORT_tau_fcas_P,
    TWORLD_BATCH_SUPPORT_gInf_fcas_P,
    TWORLD_BATCH_SUPPORT_tau_fBPf,
    TWORLD_BATCH_SUPPORT_gInf_fBPf,
    TWORLD_BATCH_SUPPORT_tau_fcaBPf,
    TWORLD_BATCH_SUPPORT_gInf_fcaBPf,
    // Ito
    TWORLD_BATCH_SUPPORT_tau_xtos,
    TWORLD_BATCH_SUPPORT_gInf_xtos,
    TWORLD_BATCH_SUPPORT_tau_xtos_p,
    TWORLD_BATCH_SUPPORT_gInf_xtos_p,
    TWORLD_BATCH_SUPPORT_tau_xtof,
    TWORLD_BATCH_SUPPORT_gInf_xtof,
    TWORLD_BATCH_SUPPORT_tau_xtof_p,
    TWORLD_BATCH_SUPPORT_gInf_xtof_p,
    TWORLD_BATCH_SUPPORT_tau_ytos,
    TWORLD_BATCH_SUPPORT_gInf_ytos,
    TWORLD_BATCH_SUPPORT_tau_ytos_p,
    TWORLD_BATCH_SUPPORT_gInf_ytos_p,
    TWORLD_BATCH_SUPPORT_tau_ytof,
    TWORLD_BATCH_SUPPORT_gInf_ytof,
    TWORLD_BATCH_SUPPORT_tau_ytof_p,
    TWORLD_BATCH_SUPPORT_gInf_ytof_p,
    // IKs
    TWORLD_BATCH_SUPPORT_tau_xs_junc,
    TWORLD_BATCH_SUPPORT_gInf_xs_junc,
    TWORLD_BATCH_SUPPORT_tau_xs_sl,
    TWORLD_BATCH_SUPPORT_gInf_xs_sl,
    // Jrel ICaL-dependent
    TWORLD_BATCH_SUPPORT_tau_jrel_act,
    TWORLD_BATCH_SUPPORT_gInf_jrel_act,
    TWORLD_BATCH_SUPPORT_tau_jrel_f1,
    TWORLD_BATCH_SUPPORT_gInf_jrel_f1,
    TWORLD_BATCH_SUPPORT_tau_jrel_f2,
    TWORLD_BATCH_SUPPORT_gInf_jrel_f2,
    // Iion
    TWORLD_BATCH_SUPPORT_Iion_cm,

    NUM_TWORLD_BATCH_SUPPORT
};

static_assert(NUM_TWORLD_BATCH_SUPPORT == 89, "TWorld compact support layout drift");

CARDIAC_HOST_DEVICE inline void TWorldComputeVariablesBatch
(
    const double t,
    const double* __restrict__ CONSTANTS,
    const int N,
    const int beginCell,
    const int endCell,
    const double* __restrict__ STATES_SoA,
    double*       __restrict__ RATES_SoA,
    double*       __restrict__ SUPPORT_SoA,
    const bool solveVmWithinODESolver,
    const StimulusProtocolPOD& stimulus
)
{
    const double Istim = computeStimulusPOD(t, stimulus);

    #define S(si) STATES_SoA[(si)*N + cellI]
    #define R(si) RATES_SoA[(si)*N + cellI]
    #define SUP(slot) SUPPORT_SoA[(slot)*N + cellI]

    for (int cellI = beginCell; cellI < endCell; ++cellI)
    {
        const double s_v                      = S(v);
        const double s_camk_trap              = S(camk_trap);
        const double s_camk_f_ICaL            = S(camk_f_ICaL);
        const double s_camk_f_RyR             = S(camk_f_RyR);
        const double s_camk_f_PLB             = S(camk_f_PLB);
        const double s_casig_serca_trap        = S(casig_serca_trap);
        const double s_buffers_NaBj           = S(buffers_NaBj);
        const double s_buffers_NaBsl          = S(buffers_NaBsl);
        const double s_buffers_TnCHc          = S(buffers_TnCHc);
        const double s_buffers_TnCHm          = S(buffers_TnCHm);
        const double s_buffers_CaM            = S(buffers_CaM);
        const double s_buffers_Myosin_ca      = S(buffers_Myosin_ca);
        const double s_buffers_Myosin_mg      = S(buffers_Myosin_mg);
        const double s_buffers_SRB            = S(buffers_SRB);
        const double s_buffers_SLLj           = S(buffers_SLLj);
        const double s_buffers_SLLsl          = S(buffers_SLLsl);
        const double s_buffers_SLHj           = S(buffers_SLHj);
        const double s_buffers_SLHsl          = S(buffers_SLHsl);
        const double s_buffers_Csqn           = S(buffers_Csqn);
        const double s_naj                    = S(naj);
        const double s_nasl                   = S(nasl);
        const double s_nai                    = S(nai);
        const double s_ki                     = S(ki);
        const double s_cli                    = S(cli);
        const double s_casr                   = S(casr);
        const double s_caj                    = S(caj);
        const double s_casl                   = S(casl);
        const double s_cai                    = S(cai);
        const double s_m                      = S(m);
        const double s_h                      = S(h);
        const double s_j                      = S(j);
        const double s_hp                     = S(hp);
        const double s_jp                     = S(jp);
        const double s_m_P                    = S(m_P);
        const double s_h_P                    = S(h_P);
        const double s_j_P                    = S(j_P);
        const double s_hp_P                   = S(hp_P);
        const double s_jp_P                   = S(jp_P);
        const double s_mL                     = S(mL);
        const double s_hL                     = S(hL);
        const double s_hLp                    = S(hLp);
        const double s_d                      = S(d);
        const double s_ff                     = S(ff);
        const double s_fs                     = S(fs);
        const double s_fcaf                   = S(fcaf);
        const double s_fcas                   = S(fcas);
        const double s_jca                    = S(jca);
        const double s_ffp                    = S(ffp);
        const double s_fcafp                  = S(fcafp);
        const double s_nca                    = S(nca);
        const double s_nca_i                  = S(nca_i);
        const double s_d_P                    = S(d_P);
        const double s_ff_P                   = S(ff_P);
        const double s_fs_P                   = S(fs_P);
        const double s_fcaf_P                 = S(fcaf_P);
        const double s_fcas_P                 = S(fcas_P);
        const double s_fBPf                   = S(fBPf);
        const double s_fcaBPf                 = S(fcaBPf);
        const double s_ical_pureCDI_junc      = S(ical_pureCDI_junc);
        const double s_ical_pureCDI_sl        = S(ical_pureCDI_sl);
        const double s_xtos                   = S(xtos);
        const double s_ytos                   = S(ytos);
        const double s_xtof                   = S(xtof);
        const double s_ytof                   = S(ytof);
        const double s_xtos_p                 = S(xtos_p);
        const double s_xtof_p                 = S(xtof_p);
        const double s_ytos_p                 = S(ytos_p);
        const double s_ytof_p                 = S(ytof_p);
        const double s_C0                     = S(C0);
        const double s_C1                     = S(C1);
        const double s_C2                     = S(C2);
        const double s_I                      = S(::I);
        const double s_O                      = S(::O);
        const double s_xs_junc                = S(xs_junc);
        const double s_xs_sl                  = S(xs_sl);
        const double s_jrel_icaldep_act       = S(jrel_icaldep_act);
        const double s_jrel_icaldep_f1        = S(jrel_icaldep_f1);
        const double s_jrel_icaldep_f2        = S(jrel_icaldep_f2);
        const double s_ryr_R                  = S(ryr_R);
        const double s_ryr_O                  = S(ryr_O);
        const double s_ryr_I                  = S(ryr_I);
        const double s_ryr_CaRI               = S(ryr_CaRI);
        const double s_ryr_R_p                = S(ryr_R_p);
        const double s_ryr_O_p                = S(ryr_O_p);
        const double s_ryr_I_p                = S(ryr_I_p);
        const double s_ryr_CaRI_p             = S(ryr_CaRI_p);
        const double s_contraction_TmBlocked  = S(contraction_TmBlocked);
        const double s_contraction_XW         = S(contraction_XW);
        const double s_contraction_XS         = S(contraction_XS);
        const double s_contraction_ZETAS      = S(contraction_ZETAS);
        const double s_contraction_ZETAW      = S(contraction_ZETAW);
        const double s_contraction_Ca_TRPN    = S(contraction_Ca_TRPN);

        /* CaMK */
        const double a_bound            = CONSTANTS[AC_CaMK0] * (1.0 - s_camk_trap) / (1.0 + CONSTANTS[AC_Km_CaMK_Ca] / s_caj);
        const double a_bound_serca      = CONSTANTS[AC_CaMK0] * (1.0 - s_casig_serca_trap) / (1.0 + CONSTANTS[AC_Km_CaMK_Ca] / s_caj);
        const double a_CaMK_active      = a_bound + s_camk_trap;
        const double a_casig_SERCA_act  = a_bound_serca + s_casig_serca_trap;
        const double a_CaMK_Phos_ss_ICaL = a_CaMK_active / (a_CaMK_active + CONSTANTS[AC_K_Phos_CaMK]);
        const double a_CaMK_Phos_ss_PLB  = a_CaMK_active / (a_CaMK_active + 10.0);
        const double a_CaMK_Phos_ss_RyR  = a_CaMK_active / (a_CaMK_active + 1.0);
        const double r_camk_trap        = CONSTANTS[AC_CaMK_alpha] * a_bound * a_CaMK_active - CONSTANTS[AC_CaMK_beta] * s_camk_trap * (0.1 + 0.9 * CONSTANTS[AC_PP1_tot] / 0.1371);
        const double r_casig_serca_trap = CONSTANTS[AC_alpha_serca] * a_bound_serca * a_casig_SERCA_act - CONSTANTS[AC_CaMK_beta] * s_casig_serca_trap * (0.1 + 0.9 * CONSTANTS[AC_PP1_tot] / 0.1371);
        const double r_camk_f_ICaL      = (a_CaMK_Phos_ss_ICaL - s_camk_f_ICaL) / CONSTANTS[AC_tau_cal];
        const double r_camk_f_PLB       = (a_CaMK_Phos_ss_PLB  - s_camk_f_PLB)  / CONSTANTS[AC_tau_plb];
        const double r_camk_f_RyR       = (a_CaMK_Phos_ss_RyR  - s_camk_f_RyR)  / CONSTANTS[AC_tau_ryr];

        /* INaCa */
        const double a_allo_i   = 1.0 / (1.0 + pow(CONSTANTS[AC_KmCaAct] / s_casl, 2.0));
        const double a_allo_ss  = 1.0 / (1.0 + pow(CONSTANTS[AC_KmCaAct] / s_caj,  2.0));
        const double a_h4_i     = 1.0 + s_nai  / CONSTANTS[AC_kna1] * (1.0 + s_nai  / CONSTANTS[AC_kna2]);
        const double a_h4_ss    = 1.0 + s_nasl / CONSTANTS[AC_kna1] * (1.0 + s_nasl / CONSTANTS[AC_kna2]);
        const double a_hca      = exp(CONSTANTS[AC_qca] * (s_v - 8.3117) * CONSTANTS[AC_F] / (CONSTANTS[AC_R] * CONSTANTS[AC_T]));
        const double a_hna      = exp(CONSTANTS[AC_qna] * (s_v - 8.3117) * CONSTANTS[AC_F] / (CONSTANTS[AC_R] * CONSTANTS[AC_T]));
        const double a_h1_i     = 1.0 + s_nai  / CONSTANTS[AC_kna3] * (1.0 + a_hna);
        const double a_h1_ss    = 1.0 + s_nasl / CONSTANTS[AC_kna3] * (1.0 + a_hna);
        const double a_h5_i     = s_nai  * s_nai  / (a_h4_i  * CONSTANTS[AC_kna1] * CONSTANTS[AC_kna2]);
        const double a_h5_ss    = s_nasl * s_nasl / (a_h4_ss * CONSTANTS[AC_kna1] * CONSTANTS[AC_kna2]);
        const double a_h6_i     = 1.0 / a_h4_i;
        const double a_h6_ss    = 1.0 / a_h4_ss;
        const double a_h7_i     = 1.0 + CONSTANTS[AC_nao] / CONSTANTS[AC_kna3] * (1.0 + 1.0 / a_hna);
        const double a_h7_ss    = a_h7_i;
        const double a_h2_i     = s_nai  * a_hna / (CONSTANTS[AC_kna3] * a_h1_i);
        const double a_h2_ss    = s_nasl * a_hna / (CONSTANTS[AC_kna3] * a_h1_ss);
        const double a_h3_i     = 1.0 / a_h1_i;
        const double a_h3_ss    = 1.0 / a_h1_ss;
        const double a_h8_i     = CONSTANTS[AC_nao] / (CONSTANTS[AC_kna3] * a_hna * a_h7_i);
        const double a_h8_ss    = CONSTANTS[AC_nao] / (CONSTANTS[AC_kna3] * a_hna * a_h7_ss);
        const double a_h9_i     = 1.0 / a_h7_i;
        const double a_h9_ss    = 1.0 / a_h7_ss;
        const double a_k6_i     = a_h6_i * s_cai  * CONSTANTS[AC_kcaon];
        const double a_k6_ss    = a_h6_ss * s_casl * CONSTANTS[AC_kcaon];
        const double a_k3p_i    = a_h9_i  * CONSTANTS[AC_wca];
        const double a_k3p_ss   = a_h9_ss * CONSTANTS[AC_wca];
        const double a_k3pp_i   = a_h8_i  * CONSTANTS[AC_wnaca];
        const double a_k3pp_ss  = a_h8_ss * CONSTANTS[AC_wnaca];
        const double a_k4p_i    = a_h3_i  * CONSTANTS[AC_wca] / a_hca;
        const double a_k4p_ss   = a_h3_ss * CONSTANTS[AC_wca] / a_hca;
        const double a_k4pp_i   = a_h2_i  * CONSTANTS[AC_wnaca];
        const double a_k4pp_ss  = a_h2_ss * CONSTANTS[AC_wnaca];
        const double a_k7_i     = a_h5_i  * a_h2_i  * CONSTANTS[AC_wna];
        const double a_k7_ss    = a_h5_ss * a_h2_ss * CONSTANTS[AC_wna];
        const double a_k8_i     = a_h8_i  * CONSTANTS[AC_h11_i]  * CONSTANTS[AC_wna];
        const double a_k8_ss    = a_h8_ss * CONSTANTS[AC_h11_ss] * CONSTANTS[AC_wna];
        const double a_k3_i     = a_k3p_i  + a_k3pp_i;
        const double a_k3_ss    = a_k3p_ss + a_k3pp_ss;
        const double a_k4_i     = a_k4p_i  + a_k4pp_i;
        const double a_k4_ss    = a_k4p_ss + a_k4pp_ss;
        const double a_x1_i     = CONSTANTS[AC_k2_i]  * a_k4_i  * (a_k7_i  + a_k6_i)  + CONSTANTS[AC_k5_i]  * a_k7_i  * (CONSTANTS[AC_k2_i]  + a_k3_i);
        const double a_x1_ss    = CONSTANTS[AC_k2_ss] * a_k4_ss * (a_k7_ss + a_k6_ss) + CONSTANTS[AC_k5_ss] * a_k7_ss * (CONSTANTS[AC_k2_ss] + a_k3_ss);
        const double a_x2_i     = CONSTANTS[AC_k1_i]  * a_k7_i  * (a_k4_i  + CONSTANTS[AC_k5_i])  + a_k4_i  * a_k6_i  * (CONSTANTS[AC_k1_i]  + a_k8_i);
        const double a_x2_ss    = CONSTANTS[AC_k1_ss] * a_k7_ss * (a_k4_ss + CONSTANTS[AC_k5_ss]) + a_k4_ss * a_k6_ss * (CONSTANTS[AC_k1_ss] + a_k8_ss);
        const double a_x3_i     = CONSTANTS[AC_k1_i]  * a_k3_i  * (a_k7_i  + a_k6_i)  + a_k8_i  * a_k6_i  * (CONSTANTS[AC_k2_i]  + a_k3_i);
        const double a_x3_ss    = CONSTANTS[AC_k1_ss] * a_k3_ss * (a_k7_ss + a_k6_ss) + a_k8_ss * a_k6_ss * (CONSTANTS[AC_k2_ss] + a_k3_ss);
        const double a_x4_i     = CONSTANTS[AC_k2_i]  * a_k8_i  * (a_k4_i  + CONSTANTS[AC_k5_i])  + a_k3_i  * CONSTANTS[AC_k5_i]  * (CONSTANTS[AC_k1_i]  + a_k8_i);
        const double a_x4_ss    = CONSTANTS[AC_k2_ss] * a_k8_ss * (a_k4_ss + CONSTANTS[AC_k5_ss]) + a_k3_ss * CONSTANTS[AC_k5_ss] * (CONSTANTS[AC_k1_ss] + a_k8_ss);
        const double a_denom_i  = a_x1_i  + a_x2_i  + a_x3_i  + a_x4_i;
        const double a_denom_ss = a_x1_ss + a_x2_ss + a_x3_ss + a_x4_ss;
        const double a_E1_i     = a_x1_i  / a_denom_i;
        const double a_E1_ss    = a_x1_ss / a_denom_ss;
        const double a_E2_i     = a_x2_i  / a_denom_i;
        const double a_E2_ss    = a_x2_ss / a_denom_ss;
        const double a_E3_i     = a_x3_i  / a_denom_i;
        const double a_E3_ss    = a_x3_ss / a_denom_ss;
        const double a_E4_i     = a_x4_i  / a_denom_i;
        const double a_E4_ss    = a_x4_ss / a_denom_ss;
        const double a_JncxCa_i   = a_E2_i  * CONSTANTS[AC_k2_i]  - a_E1_i  * CONSTANTS[AC_k1_i];
        const double a_JncxCa_ss  = a_E2_ss * CONSTANTS[AC_k2_ss] - a_E1_ss * CONSTANTS[AC_k1_ss];
        const double a_JncxNa_i   = 3.0 * (a_E4_i  * a_k7_i  - a_E1_i  * a_k8_i)  + a_E3_i  * a_k4pp_i  - a_E2_i  * a_k3pp_i;
        const double a_JncxNa_ss  = 3.0 * (a_E4_ss * a_k7_ss - a_E1_ss * a_k8_ss) + a_E3_ss * a_k4pp_ss - a_E2_ss * a_k3pp_ss;
        const double a_INaCa_i    = (1.0 - CONSTANTS[AC_INaCa_fractionSS]) * CONSTANTS[AC_Gncx] * a_allo_i  * (CONSTANTS[AC_zna] * a_JncxNa_i  + CONSTANTS[AC_zca] * a_JncxCa_i);
        const double a_INaCa_ss   = CONSTANTS[AC_INaCa_fractionSS]          * CONSTANTS[AC_Gncx] * a_allo_ss * (CONSTANTS[AC_zna] * a_JncxNa_ss + CONSTANTS[AC_zca] * a_JncxCa_ss);

        /* SERCA */
        const double a_phosphorylationTotal = s_camk_f_PLB + CONSTANTS[AC_fPLB_PKA] - s_camk_f_PLB * CONSTANTS[AC_fPLB_PKA];
        const double a_Vmax_mult   = 1.0 + CONSTANTS[AC_Max_Vmax_SERCA_Ca] / (1.0 + pow(CONSTANTS[AC_Km_SERCA_Ca] / a_casig_SERCA_act, 2.0));
        const double a_Jserca_np   = pow(CONSTANTS[AC_Q10SRCaP], CONSTANTS[AC_Qpow]) * CONSTANTS[AC_Vmax_SRCaP] * a_Vmax_mult * (pow(s_cai / CONSTANTS[AC_Kmf],   CONSTANTS[AC_hillSRCaP]) - pow(s_casr / CONSTANTS[AC_Kmr], CONSTANTS[AC_hillSRCaP])) / (1.0 + pow(s_cai / CONSTANTS[AC_Kmf],   CONSTANTS[AC_hillSRCaP]) + pow(s_casr / CONSTANTS[AC_Kmr], CONSTANTS[AC_hillSRCaP]));
        const double a_Jserca_p    = pow(CONSTANTS[AC_Q10SRCaP], CONSTANTS[AC_Qpow]) * CONSTANTS[AC_Vmax_SRCaP] * a_Vmax_mult * (pow(s_cai / CONSTANTS[AC_Kmf_p], CONSTANTS[AC_hillSRCaP]) - pow(s_casr / CONSTANTS[AC_Kmr], CONSTANTS[AC_hillSRCaP])) / (1.0 + pow(s_cai / CONSTANTS[AC_Kmf_p], CONSTANTS[AC_hillSRCaP]) + pow(s_casr / CONSTANTS[AC_Kmr], CONSTANTS[AC_hillSRCaP]));
        const double a_Jserca      = a_Jserca_np * (1.0 - a_phosphorylationTotal) + a_Jserca_p * a_phosphorylationTotal;

        /* contraction */
        const double a_XU              = 1.0 - s_contraction_TmBlocked - s_contraction_XW - s_contraction_XS;
        const double a_gamma_rate      = ((s_contraction_ZETAS > 0.0) ? CONSTANTS[AC_contraction_gamma] * s_contraction_ZETAS : ((s_contraction_ZETAS < -1.0) ? CONSTANTS[AC_contraction_gamma] * (-s_contraction_ZETAS - 1.0) : 0.0));
        const double a_gamma_rate_w    = CONSTANTS[AC_gamma_wu] * fabs(s_contraction_ZETAW);
        const double a_xb_su_gamma     = a_gamma_rate * s_contraction_XS;
        const double a_xb_uw           = CONSTANTS[AC_k_uw] * a_XU * (1.0 + CONSTANTS[AC_fMyBPC_PKA] / 2.0);
        const double a_xb_ws           = CONSTANTS[AC_k_ws] * s_contraction_XW;
        const double a_xb_wu_gamma     = a_gamma_rate_w * s_contraction_XW;
        const double a_d_contraction_Ca_TRPN = CONSTANTS[AC_koff] * (pow(s_cai * 1000.0 / CONSTANTS[AC_ca50], CONSTANTS[AC_TRPN_n]) * (1.0 - s_contraction_Ca_TRPN) - s_contraction_Ca_TRPN);
        const double a_xb_su           = CONSTANTS[AC_k_su] * s_contraction_XS;
        const double a_xb_wu           = CONSTANTS[AC_k_wu] * s_contraction_XW;
        const double r_contraction_TmBlocked = CONSTANTS[AC_ktm_block] * ((100.0 <= pow(s_contraction_Ca_TRPN, -CONSTANTS[AC_nperm] / 2.0)) ? 100.0 : pow(s_contraction_Ca_TRPN, -CONSTANTS[AC_nperm] / 2.0)) * a_XU - CONSTANTS[AC_ktm_unblock] * pow(s_contraction_Ca_TRPN, CONSTANTS[AC_nperm] / 2.0) * s_contraction_TmBlocked;
        const double r_contraction_ZETAS     = CONSTANTS[AC_A] * CONSTANTS[AC_lambda_rate] - CONSTANTS[AC_cds] * s_contraction_ZETAS;
        const double r_contraction_ZETAW     = CONSTANTS[AC_A] * CONSTANTS[AC_lambda_rate] - CONSTANTS[AC_cdw] * s_contraction_ZETAW;
        const double r_contraction_Ca_TRPN   = a_d_contraction_Ca_TRPN;
        const double r_contraction_XS        = a_xb_ws - a_xb_su - a_xb_su_gamma;
        const double r_contraction_XW        = a_xb_uw - a_xb_wu - a_xb_ws - a_xb_wu_gamma;

        /* reversal potentials */
        const double a_ECaj   = CONSTANTS[AC_R] * CONSTANTS[AC_T] / (CONSTANTS[AC_zca] * CONSTANTS[AC_F]) * log(CONSTANTS[AC_cao] / s_caj);
        const double a_ECasl  = CONSTANTS[AC_R] * CONSTANTS[AC_T] / (CONSTANTS[AC_zca] * CONSTANTS[AC_F]) * log(CONSTANTS[AC_cao] / s_casl);
        const double a_ECl    = CONSTANTS[AC_R] * CONSTANTS[AC_T] / (CONSTANTS[AC_zcl] * CONSTANTS[AC_F]) * log(CONSTANTS[AC_clo]  / s_cli);
        const double a_EK     = CONSTANTS[AC_R] * CONSTANTS[AC_T] / (CONSTANTS[AC_zk]  * CONSTANTS[AC_F]) * log(CONSTANTS[AC_ko]   / s_ki);
        const double a_ENaj   = CONSTANTS[AC_R] * CONSTANTS[AC_T] / (CONSTANTS[AC_zna] * CONSTANTS[AC_F]) * log(CONSTANTS[AC_nao]  / s_naj);
        const double a_ENasl  = CONSTANTS[AC_R] * CONSTANTS[AC_T] / (CONSTANTS[AC_zna] * CONSTANTS[AC_F]) * log(CONSTANTS[AC_nao]  / s_nasl);
        const double a_EKs    = CONSTANTS[AC_R] * CONSTANTS[AC_T] / (CONSTANTS[AC_zk]  * CONSTANTS[AC_F]) * log((CONSTANTS[AC_ko] + CONSTANTS[AC_PNaK] * CONSTANTS[AC_nao]) / (s_ki + CONSTANTS[AC_PNaK] * s_nai));

        /* ICl */
        const double a_IClCa_junc = 0.5 * CONSTANTS[AC_GClCa] / (1.0 + CONSTANTS[AC_KdClCa] / s_caj)  * (s_v - a_ECl);
        const double a_IClCa_sl   = 0.5 * CONSTANTS[AC_GClCa] / (1.0 + CONSTANTS[AC_KdClCa] / s_casl) * (s_v - a_ECl);
        const double a_IClb       = CONSTANTS[AC_GClb] * (s_v - a_ECl);
        const double a_IClCa      = a_IClCa_junc + a_IClCa_sl;

        /* IK1 */
        const double a_aK1   = 4.094 / (1.0 + exp(0.1217 * (s_v - a_EK - 49.934)));
        const double a_bK1   = (15.72 * exp(0.0674 * (s_v - a_EK - 3.257)) + exp(0.0618 * (s_v - a_EK - 594.31))) / (1.0 + exp(-0.1629 * (s_v - a_EK + 14.207)));
        const double a_K1ss  = a_aK1 / (a_aK1 + a_bK1);
        const double a_IK1_IK1 = CONSTANTS[AC_GK1] * sqrt(CONSTANTS[AC_ko] / 5.0) * a_K1ss * (s_v - a_EK);

        /* IKb */
        const double a_xkb     = 1.0 / (1.0 + exp(-(s_v - 10.8968) / 23.9871));
        const double a_IKb_IKb = CONSTANTS[AC_GKb] * a_xkb * (s_v - a_EK);

        /* Ito */
        const double a_dti_develop  = 1.354 + 0.0001 / (exp((s_v - 167.4) / 15.89) + exp(-(s_v - 12.23) / 0.2154));
        const double a_dti_recover  = 1.0 - 0.5 / (1.0 + exp((s_v + 70.0) / 20.0));
        const double a_fItop        = s_camk_f_RyR;
        const double a_tauxtof      = 8.5 * exp(-pow((s_v + 45.0) / 50.0, 2.0)) + 0.5;
        const double a_tauxtos      = 9.0 / (1.0 + exp((s_v + 3.0) / 15.0)) + 0.5;
        const double a_tauytof      = 85.0 * exp(-pow(s_v + 40.0, 2.0) / 220.0) + 7.0;
        const double a_tauytos      = 800.0 / (1.0 + exp((s_v + 60.0) / 10.0)) + 30.0;
        const double a_xtoss        = 1.0 / (1.0 + exp(-(s_v - 19.0) / 13.0));
        const double a_xtoss_p      = 1.0 / (1.0 + exp(-(s_v - 29.0) / 13.0));
        const double a_ytoss        = 1.0 / (1.0 + exp((s_v + 19.5) / 5.0));
        const double a_tauytof_p    = a_tauytof  * a_dti_develop * a_dti_recover;
        const double a_tauytos_p    = a_tauytos * a_dti_develop * a_dti_recover;
        const double a_Itof  = CONSTANTS[AC_Gto_fast] * (s_v - a_EK) * ((1.0 - a_fItop) * s_xtof * s_ytof + a_fItop * s_xtof_p * s_ytof_p);
        const double a_Itos  = CONSTANTS[AC_Gto_slow] * (s_v - a_EK) * ((1.0 - a_fItop) * s_xtos * s_ytos + a_fItop * s_xtos_p * s_ytos_p);
        const double a_Ito_Ito = a_Itos + a_Itof;
        const double r_xtof   = (a_xtoss   - s_xtof)   / a_tauxtof;
        const double r_xtof_p = (a_xtoss_p - s_xtof_p) / a_tauxtof;
        const double r_xtos   = (a_xtoss   - s_xtos)   / a_tauxtos;
        const double r_xtos_p = (a_xtoss_p - s_xtos_p) / a_tauxtos;
        const double r_ytof   = (a_ytoss   - s_ytof)   / a_tauytof;
        const double r_ytos   = (a_ytoss   - s_ytos)   / a_tauytos;
        const double r_ytof_p = (a_ytoss   - s_ytof_p) / a_tauytof_p;
        const double r_ytos_p = (a_ytoss   - s_ytos_p) / a_tauytos_p;

        /* ICaL */
        const double a_Afcaf     = 0.3 + 0.6 / (1.0 + exp((s_v - 9.24247) / 27.96201));
        const double a_Ii        = 0.5 * (s_nai  + s_ki + s_cli + 4.0 * s_cai)  / 1000.0;
        const double a_Iss       = 0.5 * (s_nasl + s_ki + s_cli + 4.0 * s_casl) / 1000.0;
        const double a_dPss      = ((1.0323 * exp(-1.0553 * exp(-0.081 * (s_v + 12.62483))) <= 1.0) ? 1.0323 * exp(-1.0553 * exp(-0.081 * (s_v + 12.62483))) : 1.0);
        const double a_dss       = ((1.0763 * exp(-1.007  * exp(-0.0829 * (s_v + 3.62483)))  <= 1.0) ? 1.0763 * exp(-1.007  * exp(-0.0829 * (s_v + 3.62483)))  : 1.0);
        const double a_fICaLp    = s_camk_f_ICaL;
        const double a_fss       = 1.0 / (1.0 + exp((s_v + 19.58) / 3.696));
        const double a_fss_P     = 1.0 / (1.0 + exp((s_v + 25.58) / 3.696));
        const double a_jcass     = 1.0 / (1.0 + exp((s_v + 17.66945) / 3.21501));
        const double a_km2n      = s_jca * 0.84191;
        const double a_sigmoidTransition  = 1.0 - 1.0 / (1.0 + 1.86532 * s_caj  / 0.032);
        const double a_sigmoidTransition2 = 1.0 - 1.0 / (1.0 + 1.86532 * s_casl / 0.032);
        const double a_td        = 1.5 + 1.0 / (exp(-0.05 * (s_v + 6.0)) + exp(0.09 * (s_v + 14.0)));
        const double a_tfcaf     = 13.50673 + 1.0 / (0.1542 * exp(-(s_v - 1.31611) / 11.3396) + 0.1542 * exp((s_v - 1.31611) / 11.3396));
        const double a_tfcas     = 177.95813 + 1.0 / (4.73955000000000025e-04 * exp(-(s_v - 0.79049) / 0.81777) + 4.73955000000000025e-04 * exp((s_v + 2.40474) / 1.90812));
        const double a_tff       = 6.17111 + 1.0 / (0.00126 * exp(-(s_v + 26.63596) / 9.69961) + 0.00126 * exp((s_v + 26.63596) / 9.69961));
        const double a_tfs       = 2.71922488999999996e+03 + 1.0 / (7.19411000000000014e-05 * exp(-(s_v + 5.74631) / 10.8769) + 7.19411000000000014e-05 * exp((s_v + 5.74631) / 16.31535));
        const double a_Afcas     = 1.0 - a_Afcaf;
        const double a_anca      = 1.0 / (CONSTANTS[AC_k2n] / a_km2n + pow(1.0 + CONSTANTS[AC_Kmn] / s_caj,  3.80763));
        const double a_anca_i    = 1.0 / (CONSTANTS[AC_k2n] / a_km2n + pow(1.0 + CONSTANTS[AC_Kmn] / s_casl, 3.80763));
        const double a_fBPss     = a_fss_P;
        const double a_fcass     = a_fss;
        const double a_fcass_P   = a_fss_P;
        const double a_tauTransition  = 1.0967 + (1.0 - a_sigmoidTransition)  * 141.4299;
        const double a_tauTransition2 = 1.0967 + (1.0 - a_sigmoidTransition2) * 141.4299;
        const double a_tfcafp    = 2.5 * a_tfcaf;
        const double a_tffp      = 2.5 * a_tff;
        const double r_d         = (a_dss   - s_d)      / a_td;
        const double r_d_P       = (a_dPss  - s_d_P)    / a_td;
        const double r_ff        = (a_fss   - s_ff)     / a_tff;
        const double r_ff_P      = (a_fss_P - s_ff_P)   / a_tff;
        const double r_fs        = (a_fss   - s_fs)     / a_tfs;
        const double r_fs_P      = (a_fss_P - s_fs_P)   / a_tfs;
        const double r_jca       = (a_jcass - s_jca)    / CONSTANTS[AC_tjca];
        const double a_f         = CONSTANTS[AC_Aff] * s_ff    + CONSTANTS[AC_Afs] * s_fs;
        const double a_fBP       = CONSTANTS[AC_Aff] * s_fBPf  + CONSTANTS[AC_Afs] * s_fs_P;
        const double a_fICaL_BP  = a_fICaLp * CONSTANTS[AC_fICaLP];
        const double a_f_P       = CONSTANTS[AC_Aff] * s_ff_P  + CONSTANTS[AC_Afs] * s_fs_P;
        const double a_fca       = a_Afcaf  * s_fcaf   + a_Afcas * s_fcas;
        const double a_fcaBP     = a_Afcaf  * s_fcaBPf + a_Afcas * s_fcas_P;
        const double a_fcaBPss   = a_fcass_P;
        const double a_fcap      = a_Afcaf  * s_fcafp  + a_Afcas * s_fcas;
        const double a_fcap_P    = a_Afcaf  * s_fcaf_P + a_Afcas * s_fcas_P;
        const double a_fp        = CONSTANTS[AC_Aff] * s_ffp   + CONSTANTS[AC_Afs] * s_fs;
        const double a_gamma_cai  = pow(10.0, -CONSTANTS[AC_constA] * 4.0 * (sqrt(a_Ii)  / (1.0 + sqrt(a_Ii))  - 0.3 * a_Ii));
        const double a_gamma_cass = pow(10.0, -CONSTANTS[AC_constA] * 4.0 * (sqrt(a_Iss) / (1.0 + sqrt(a_Iss)) - 0.3 * a_Iss));
        const double a_gamma_ki   = pow(10.0, -CONSTANTS[AC_constA] * 1.0 * (sqrt(a_Ii)  / (1.0 + sqrt(a_Ii))  - 0.3 * a_Ii));
        const double a_gamma_kss  = pow(10.0, -CONSTANTS[AC_constA] * 1.0 * (sqrt(a_Iss) / (1.0 + sqrt(a_Iss)) - 0.3 * a_Iss));
        const double a_gamma_nai  = a_gamma_ki;
        const double a_gamma_nass = a_gamma_kss;
        const double r_fBPf       = (a_fBPss   - s_fBPf)   / a_tffp;
        const double r_fcaf       = (a_fcass   - s_fcaf)    / a_tfcaf;
        const double r_fcaf_P     = (a_fcass_P - s_fcaf_P)  / a_tfcaf;
        const double r_fcafp      = (a_fcass   - s_fcafp)   / a_tfcafp;
        const double r_fcas       = (a_fcass   - s_fcas)    / a_tfcas;
        const double r_fcas_P     = (a_fcass_P - s_fcas_P)  / a_tfcas;
        const double r_ffp        = (a_fss     - s_ffp)     / a_tffp;
        const double r_ical_pureCDI_junc = -s_ical_pureCDI_junc * a_sigmoidTransition  / a_tauTransition  + (1.0 - s_ical_pureCDI_junc) * CONSTANTS[AC_rateRecovery];
        const double r_ical_pureCDI_sl   = -s_ical_pureCDI_sl  * a_sigmoidTransition2 / a_tauTransition2 + (1.0 - s_ical_pureCDI_sl)  * CONSTANTS[AC_rateRecovery];
        const double r_nca   = a_anca   * CONSTANTS[AC_k2n] - s_nca   * a_km2n;
        const double r_nca_i = a_anca_i * CONSTANTS[AC_k2n] - s_nca_i * a_km2n;
        const double a_fICaL_CaMKonly = a_fICaLp - a_fICaL_BP;
        const double a_fICaL_PKAonly  = CONSTANTS[AC_fICaL_P] - a_fICaL_BP;
        const double r_fcaBPf = (a_fcaBPss - s_fcaBPf) / a_tfcafp;

        /* IKr (Markov — Euler) */
        const double a_IKr_IKr = CONSTANTS[AC_GKr] * sqrt(CONSTANTS[AC_ko] / 5.0) * s_O * (s_v - a_EK);

        /* IKs */
        const double a_VsTs_Ca_junc  = CONSTANTS[AC_P_tau_0] + (CONSTANTS[AC_P_tau_max] - CONSTANTS[AC_P_tau_0]) / (1.0 + pow(0.00015 / s_caj,  3.0));
        const double a_VsTs_Ca_sl    = CONSTANTS[AC_P_tau_0] + (CONSTANTS[AC_P_tau_max] - CONSTANTS[AC_P_tau_0]) / (1.0 + pow(0.00015 / s_casl, 3.0));
        const double a_VsXs_Ca_junc  = CONSTANTS[AC_P_vh_0]  + (CONSTANTS[AC_P_vh_max]  - CONSTANTS[AC_P_vh_0])  / (1.0 + pow(0.00035 / s_caj,  4.0));
        const double a_VsXs_Ca_sl    = CONSTANTS[AC_P_vh_0]  + (CONSTANTS[AC_P_vh_max]  - CONSTANTS[AC_P_vh_0])  / (1.0 + pow(0.00035 / s_casl, 4.0));
        const double a_gks_junc      = CONSTANTS[AC_P_g_0]   + (CONSTANTS[AC_P_g_max]   - CONSTANTS[AC_P_g_0])   / (1.0 + pow(0.00015 / s_caj,  1.3));
        const double a_gks_sl        = CONSTANTS[AC_P_g_0]   + (CONSTANTS[AC_P_g_max]   - CONSTANTS[AC_P_g_0])   / (1.0 + pow(0.00015 / s_casl, 1.3));
        const double a_tauxs_junc    = 2.0 * (50.0 + (50.0 + 350.0 * exp(-pow(s_v + 30.0, 2.0) / 4000.0)) / (1.0 + exp(-(s_v + a_VsTs_Ca_junc) / 10.0)));
        const double a_tauxs_sl      = 2.0 * (50.0 + (50.0 + 350.0 * exp(-pow(s_v + 30.0, 2.0) / 4000.0)) / (1.0 + exp(-(s_v + a_VsTs_Ca_sl)  / 10.0)));
        const double a_xsss_junc     = 1.0 / (1.0 + exp(-(s_v - a_VsXs_Ca_junc) / 25.0));
        const double a_xsss_sl       = 1.0 / (1.0 + exp(-(s_v - a_VsXs_Ca_sl)   / 25.0));
        const double r_xs_junc       = (a_xsss_junc - s_xs_junc) / a_tauxs_junc;
        const double r_xs_sl         = (a_xsss_sl   - s_xs_sl)   / a_tauxs_sl;
        const double a_IKs_junc      = CONSTANTS[AC_Fjunc] * CONSTANTS[AC_GKs] * a_gks_junc * pow(s_xs_junc, 2.0) * (s_v - a_EKs);
        const double a_IKs_sl        = CONSTANTS[AC_Fsl]   * CONSTANTS[AC_GKs] * a_gks_sl   * pow(s_xs_sl,   2.0) * (s_v - a_EKs);
        const double a_IKs_IKs       = a_IKs_junc + a_IKs_sl;

        /* INa */
        const double a_ah    = ((s_v >= -40.0) ? 0.0 : 0.057 * exp(-(s_v + 80.0) / 6.8));
        const double a_aj    = ((s_v >= -40.0) ? 0.0 : (-25428.0 * exp(0.2444 * s_v) - 6.948e-06 * exp(-0.04391 * s_v)) * (s_v + 37.78) / (1.0 + exp(0.311 * (s_v + 79.23))));
        const double a_bh    = ((s_v >= -40.0) ? 0.77 / (0.13 * (1.0 + exp(-(s_v + 10.66) / 11.1))) : 2.7 * exp(0.079 * s_v) + 310000.0 * exp(0.3485 * s_v));
        const double a_bj    = ((s_v >= -40.0) ? 0.6 * exp(0.057 * s_v) / (1.0 + exp(-0.1 * (s_v + 32.0))) : 0.02424 * exp(-0.01052 * s_v) / (1.0 + exp(-0.1378 * (s_v + 40.14))));
        const double a_fINap = s_camk_f_RyR;
        const double a_hss   = 1.0 / pow(1.0 + exp((s_v + 71.55) / 7.43), 2.0);
        const double a_hss_P = 1.0 / pow(1.0 + exp((s_v + 76.55) / 7.43), 2.0);
        const double a_hssp  = 1.0 / pow(1.0 + exp((s_v + 77.55) / 7.43), 2.0);
        const double a_hssp_P = 1.0 / pow(1.0 + exp((s_v + 82.55) / 7.43), 2.0);
        const double a_mss   = 1.0 / pow(1.0 + exp(-(s_v + 56.86) / 9.03), 2.0);
        const double a_mss_P = 1.0 / pow(1.0 + exp(-(s_v + 61.86) / 9.03), 2.0);
        const double a_tm    = 0.1292 * exp(-pow((s_v + 45.79) / 15.54, 2.0)) + 0.06487 * exp(-pow((s_v - 4.823) / 51.12, 2.0));
        const double a_fINa_BP  = a_fINap * CONSTANTS[AC_fINa_P];
        const double a_jss      = a_hss;
        const double a_jss_P    = a_hss_P;
        const double a_jssp_P   = a_hssp_P;
        const double a_th    = 1.0 / (a_ah + a_bh);
        const double a_tj    = 1.0 / (a_aj + a_bj);
        const double r_m     = (a_mss   - s_m)    / a_tm;
        const double r_m_P   = (a_mss_P - s_m_P)  / a_tm;
        const double a_INaBase_CaMK = CONSTANTS[AC_GNa] * pow(s_m,   3.0) * s_hp    * s_jp;
        const double a_INaBase_NP   = CONSTANTS[AC_GNa] * pow(s_m,   3.0) * s_h     * s_j;
        const double a_fINa_CaMKonly = a_fINap - a_fINa_BP;
        const double a_fINa_PKAonly  = CONSTANTS[AC_fINa_P] - a_fINa_BP;
        const double a_tjp   = 1.46 * a_tj;
        const double r_h     = (a_hss   - s_h)    / a_th;
        const double r_h_P   = (a_hss_P - s_h_P)  / a_th;
        const double r_hp    = (a_hssp  - s_hp)   / a_th;
        const double r_hp_P  = (a_hssp_P - s_hp_P) / a_th;
        const double r_j     = (a_jss   - s_j)    / a_tj;
        const double r_j_P   = (a_jss_P - s_j_P)  / a_tj;
        const double a_INaBase_BP  = CONSTANTS[AC_GNa_P] * pow(s_m_P, 3.0) * s_hp_P * s_jp_P;
        const double a_INaBase_PKA = CONSTANTS[AC_GNa_P] * pow(s_m_P, 3.0) * s_h_P  * s_j_P;
        const double r_jp    = (a_jss   - s_jp)   / a_tjp;
        const double r_jp_P  = (a_jssp_P - s_jp_P) / a_tjp;
        const double a_INaBase = (1.0 - a_fINa_CaMKonly - a_fINa_PKAonly - a_fINa_BP) * a_INaBase_NP + a_fINa_CaMKonly * a_INaBase_CaMK + a_fINa_PKAonly * a_INaBase_PKA + a_fINa_BP * a_INaBase_BP;

        /* INaK */
        const double a_fnak = 0.75 + (0.00375 - (140.0 - CONSTANTS[AC_nao]) / 50.0 * 0.001) * s_v;

        /* INaL */
        const double a_fINaLp  = s_camk_f_RyR;
        const double a_hLss    = 1.0 / (1.0 + exp((s_v + 87.61) / 7.488));
        const double a_hLssp   = 1.0 / (1.0 + exp((s_v + 93.81) / 7.488));
        const double a_mLss    = 1.0 / (1.0 + exp(-(s_v + 42.85) / 5.264));
        const double a_tmL     = 0.1292 * exp(-pow((s_v + 45.79) / 15.54, 2.0)) + 0.06487 * exp(-pow((s_v - 4.823) / 51.12, 2.0));
        const double a_GNaL    = CONSTANTS[AC_GNaL_b] * CONSTANTS[AC_INaL_multiplier] * (1.0 + a_fINaLp) * CONSTANTS[gnalTissueScale];
        const double r_hL      = (a_hLss  - s_hL)  / CONSTANTS[AC_thL];
        const double r_mL      = (a_mLss  - s_mL)  / a_tmL;
        const double r_hLp     = (a_hLssp - s_hLp) / CONSTANTS[AC_thLp];

        /* intracellular_ions */
        const double a_ICl_tot = a_IClCa + a_IClb;
        const double r_buffers_Myosin_mg = CONSTANTS[AC_kon_myomg] * CONSTANTS[AC_mgi] * (CONSTANTS[AC_Bmax_myosin] - s_buffers_Myosin_ca - s_buffers_Myosin_mg) - CONSTANTS[AC_koff_myomg] * s_buffers_Myosin_mg;
        const double r_buffers_TnCHm     = CONSTANTS[AC_kon_tnchmg] * CONSTANTS[AC_mgi] * (CONSTANTS[AC_Bmax_TnChigh] - s_buffers_TnCHc - s_buffers_TnCHm) - CONSTANTS[AC_koff_tnchmg] * s_buffers_TnCHm;
        const double r_nai               = CONSTANTS[AC_J_na_slmyo] / CONSTANTS[AC_vmyo] * (s_nasl - s_nai);
        const double a_d_buffers_CaM     = CONSTANTS[AC_kon_cam]    * s_cai  * (CONSTANTS[AC_Bmax_CaM]    - s_buffers_CaM)    - CONSTANTS[AC_koff_cam]    * s_buffers_CaM;
        const double a_d_buffers_Csqn    = CONSTANTS[AC_kon_csqn]   * s_casr * (CONSTANTS[AC_Bmax_Csqn]   - s_buffers_Csqn)   - CONSTANTS[AC_koff_csqn]   * s_buffers_Csqn;
        const double a_d_buffers_Myosin_ca = CONSTANTS[AC_kon_myoca] * s_cai  * (CONSTANTS[AC_Bmax_myosin] - s_buffers_Myosin_ca - s_buffers_Myosin_mg) - CONSTANTS[AC_koff_myoca] * s_buffers_Myosin_ca;
        const double a_d_buffers_NaBj    = CONSTANTS[AC_kon_na]  * s_naj  * (CONSTANTS[AC_Bmax_Naj]    - s_buffers_NaBj)    - CONSTANTS[AC_koff_na]  * s_buffers_NaBj;
        const double a_d_buffers_NaBsl   = CONSTANTS[AC_kon_na]  * s_nasl * (CONSTANTS[AC_Bmax_Nasl]   - s_buffers_NaBsl)   - CONSTANTS[AC_koff_na]  * s_buffers_NaBsl;
        const double a_d_buffers_SLHj    = CONSTANTS[AC_kon_slh] * s_caj  * (CONSTANTS[AC_Bmax_SLhighj] - s_buffers_SLHj)  - CONSTANTS[AC_koff_slh] * s_buffers_SLHj;
        const double a_d_buffers_SLHsl   = CONSTANTS[AC_kon_slh] * s_casl * (CONSTANTS[AC_Bmax_SLhighsl] - s_buffers_SLHsl) - CONSTANTS[AC_koff_slh] * s_buffers_SLHsl;
        const double a_d_buffers_SLLj    = CONSTANTS[AC_kon_sll] * s_caj  * (CONSTANTS[AC_Bmax_SLlowj]  - s_buffers_SLLj)  - CONSTANTS[AC_koff_sll] * s_buffers_SLLj;
        const double a_d_buffers_SLLsl   = CONSTANTS[AC_kon_sll] * s_casl * (CONSTANTS[AC_Bmax_SLlowsl]  - s_buffers_SLLsl) - CONSTANTS[AC_koff_sll] * s_buffers_SLLsl;
        const double a_d_buffers_SRB     = CONSTANTS[AC_kon_sr]  * s_cai  * (CONSTANTS[AC_Bmax_SR]      - s_buffers_SRB)    - CONSTANTS[AC_koff_sr]  * s_buffers_SRB;
        const double a_d_buffers_TnCHc   = CONSTANTS[AC_kon_tnchca] * s_cai * (CONSTANTS[AC_Bmax_TnChigh] - s_buffers_TnCHc - s_buffers_TnCHm) - CONSTANTS[AC_koff_tnchca] * s_buffers_TnCHc;
        const double a_d_buffers_TnClow  = CONSTANTS[AC_Bmax_TnClow] * a_d_contraction_Ca_TRPN;
        const double a_J_CaB_cytosol = a_d_buffers_TnClow + a_d_buffers_TnCHc + a_d_buffers_CaM + a_d_buffers_Myosin_ca + a_d_buffers_SRB;
        const double a_J_CaB_junction = a_d_buffers_SLLj + a_d_buffers_SLHj;
        const double a_J_CaB_sl       = a_d_buffers_SLLsl + a_d_buffers_SLHsl;
        const double r_cai = -a_Jserca * CONSTANTS[AC_vsr] / CONSTANTS[AC_vmyo] - a_J_CaB_cytosol + CONSTANTS[AC_J_ca_slmyo] / CONSTANTS[AC_vmyo] * (s_casl - s_cai);

        /* membrane */
        const double a_vffrt = s_v * CONSTANTS[AC_F] * CONSTANTS[AC_F] / (CONSTANTS[AC_R] * CONSTANTS[AC_T]);
        const double a_vfrt  = s_v * CONSTANTS[AC_F] / (CONSTANTS[AC_R] * CONSTANTS[AC_T]);

        /* ryr */
        const double a_CaMKIILeakMultiplier = 1.0 + 2.0 * s_camk_f_RyR;
        const double a_Jrel_ICaLdep  = 0.00174 * s_jrel_icaldep_act * s_jrel_icaldep_f1 * s_jrel_icaldep_f2;
        const double a_RIcleft       = 1.0 - s_ryr_R   - s_ryr_O   - s_ryr_I   - s_ryr_CaRI;
        const double a_RIcleft_p     = 1.0 - s_ryr_R_p - s_ryr_O_p - s_ryr_I_p - s_ryr_CaRI_p;
        const double a_nonlinearModifier  = 0.2144 * exp(1.83 * s_casr);
        const double a_J_SRCarel_np  = CONSTANTS[AC_ks] * s_ryr_O   * (s_casr - s_caj) + a_Jrel_ICaLdep;
        const double a_J_SRCarel_p   = CONSTANTS[AC_ks] * s_ryr_O_p * (s_casr - s_caj) + a_Jrel_ICaLdep;
        const double a_J_SRleak      = 1.59306000000000000e-06 * (s_casr - s_caj) * a_nonlinearModifier * a_CaMKIILeakMultiplier;
        const double a_kCaSR         = CONSTANTS[AC_MaxSR] - (CONSTANTS[AC_MaxSR] - CONSTANTS[AC_MinSR]) / (1.0 + pow(CONSTANTS[AC_ec50SR] / s_casr, CONSTANTS[AC_steepnessCaSR]));
        const double a_sigmoidBaseCaI = CONSTANTS[AC_minCaI] + (CONSTANTS[AC_maxCaI] - CONSTANTS[AC_minCaI]) / (1.0 + pow(CONSTANTS[AC_ecCaI] / s_caj, CONSTANTS[AC_steepnessCaI]));
        const double a_tau_rel       = ((CONSTANTS[AC_bt] / (1.0 + 0.0123 / s_casr) >= 0.001) ? CONSTANTS[AC_bt] / (1.0 + 0.0123 / s_casr) : 0.001);
        const double a_J_SRCarel     = a_J_SRCarel_p * s_camk_f_RyR + a_J_SRCarel_np * (1.0 - s_camk_f_RyR);
        const double a_RI_to_CI      = CONSTANTS[AC_baseRateCaI] * a_sigmoidBaseCaI;
        const double a_kiSRCa        = CONSTANTS[AC_kiCa] * a_kCaSR;
        const double a_koSRCa        = CONSTANTS[AC_koCa] / a_kCaSR;
        const double r_ryr_CaRI      = a_RI_to_CI * a_RIcleft   - CONSTANTS[AC_CI_to_RI] * s_ryr_CaRI;
        const double r_ryr_CaRI_p    = a_RI_to_CI * a_RIcleft_p - CONSTANTS[AC_CI_to_RI] * s_ryr_CaRI_p;
        const double r_ryr_I         = a_kiSRCa * CONSTANTS[AC_caTransFactor] * pow(s_caj, CONSTANTS[AC_caExpFactor]) * s_ryr_O   - CONSTANTS[AC_kim] * s_ryr_I   - (CONSTANTS[AC_kom] * s_ryr_I   - CONSTANTS[AC_caTransFactor2]  * a_koSRCa * pow(s_caj, CONSTANTS[AC_caExpFactor2]) * a_RIcleft);
        const double r_ryr_I_p       = a_kiSRCa * CONSTANTS[AC_caTransFactor] * pow(s_caj, CONSTANTS[AC_caExpFactor]) * s_ryr_O_p - CONSTANTS[AC_kim] * s_ryr_I_p - (CONSTANTS[AC_kom] * s_ryr_I_p - CONSTANTS[AC_caTransFactor2p] * a_koSRCa * pow(s_caj, CONSTANTS[AC_caExpFactor2]) * a_RIcleft_p);
        const double r_ryr_O         = CONSTANTS[AC_caTransFactor2]  * a_koSRCa * pow(s_caj, CONSTANTS[AC_caExpFactor2]) * s_ryr_R   - CONSTANTS[AC_kom] * s_ryr_O   - (a_kiSRCa * CONSTANTS[AC_caTransFactor] * pow(s_caj, CONSTANTS[AC_caExpFactor]) * s_ryr_O   - CONSTANTS[AC_kim] * s_ryr_I);
        const double r_ryr_O_p       = CONSTANTS[AC_caTransFactor2p] * a_koSRCa * pow(s_caj, CONSTANTS[AC_caExpFactor2]) * s_ryr_R_p - CONSTANTS[AC_kom] * s_ryr_O_p - (a_kiSRCa * CONSTANTS[AC_caTransFactor] * pow(s_caj, CONSTANTS[AC_caExpFactor]) * s_ryr_O_p - CONSTANTS[AC_kim] * s_ryr_I_p);
        const double r_ryr_R         = CONSTANTS[AC_kim] * a_RIcleft   - a_kiSRCa * CONSTANTS[AC_caTransFactor] * pow(s_caj, CONSTANTS[AC_caExpFactor]) * s_ryr_R   - (CONSTANTS[AC_caTransFactor2]  * a_koSRCa * pow(s_caj, CONSTANTS[AC_caExpFactor2]) * s_ryr_R   - CONSTANTS[AC_kom] * s_ryr_O);
        const double r_ryr_R_p       = CONSTANTS[AC_kim] * a_RIcleft_p - a_kiSRCa * CONSTANTS[AC_caTransFactor] * pow(s_caj, CONSTANTS[AC_caExpFactor]) * s_ryr_R_p - (CONSTANTS[AC_caTransFactor2p] * a_koSRCa * pow(s_caj, CONSTANTS[AC_caExpFactor2]) * s_ryr_R_p - CONSTANTS[AC_kom] * s_ryr_O_p);

        /* remaining currents */
        const double a_PhiCaK_i  = 1.0 * a_vffrt * (a_gamma_ki  * s_ki  * exp(1.0 * a_vfrt) - CONSTANTS[AC_gamma_ko] * CONSTANTS[AC_ko]) / (exp(1.0 * a_vfrt) - 1.0);
        const double a_PhiCaK_ss = 1.0 * a_vffrt * (a_gamma_kss * s_ki  * exp(1.0 * a_vfrt) - CONSTANTS[AC_gamma_ko] * CONSTANTS[AC_ko]) / (exp(1.0 * a_vfrt) - 1.0);
        const double a_PhiCaL_i  = 4.0 * a_vffrt * (a_gamma_cai  * s_cai  * exp(2.0 * a_vfrt) - CONSTANTS[AC_gamma_cao] * CONSTANTS[AC_cao]) / (exp(2.0 * a_vfrt) - 1.0);
        const double a_PhiCaL_ss = 4.0 * a_vffrt * (a_gamma_cass * s_casl * exp(2.0 * a_vfrt) - CONSTANTS[AC_gamma_cao] * CONSTANTS[AC_cao]) / (exp(2.0 * a_vfrt) - 1.0);
        const double a_PhiCaNa_i = 1.0 * a_vffrt * (a_gamma_nai  * s_nai  * exp(1.0 * a_vfrt) - CONSTANTS[AC_gamma_nao] * CONSTANTS[AC_nao]) / (exp(1.0 * a_vfrt) - 1.0);
        const double a_PhiCaNa_ss= 1.0 * a_vffrt * (a_gamma_nass * s_nasl * exp(1.0 * a_vfrt) - CONSTANTS[AC_gamma_nao] * CONSTANTS[AC_nao]) / (exp(1.0 * a_vfrt) - 1.0);
        const double a_ICab_junc = CONSTANTS[AC_Fjunc] * CONSTANTS[AC_GCab] * (s_v - a_ECaj);
        const double a_ICab_sl   = CONSTANTS[AC_Fsl]   * CONSTANTS[AC_GCab] * (s_v - a_ECasl);
        const double a_IKr_alpha = 0.1161 * exp(0.299  * a_vfrt);
        const double a_alpha_2   = 0.0578 * exp(0.971  * a_vfrt);
        const double a_alpha_C2ToI = 5.2e-05 * exp(1.525 * a_vfrt);
        const double a_alpha_i   = 0.2533 * exp(0.5953 * a_vfrt);
        const double a_IKr_beta  = 0.2442 * exp(-1.604 * a_vfrt);
        const double a_beta_2    = 0.000349 * exp(-1.062 * a_vfrt);
        const double a_beta_i    = 0.04568  * exp(-0.8209 * a_vfrt);
        const double a_INaj      = CONSTANTS[AC_Fjunc] * a_INaBase * (s_v - a_ENaj);
        const double a_INasl     = CONSTANTS[AC_Fsl]   * a_INaBase * (s_v - a_ENasl);
        const double a_INaK_junc_PKA   = CONSTANTS[AC_Fjunc] * CONSTANTS[AC_IbarNaK] * a_fnak / (1.0 + pow(CONSTANTS[AC_KmNaip_PKA] / s_naj, 4.0)) * CONSTANTS[AC_ko] / (CONSTANTS[AC_ko] + CONSTANTS[AC_KmKo]);
        const double a_INaK_junc_noPKA = CONSTANTS[AC_Fjunc] * CONSTANTS[AC_IbarNaK] * a_fnak / (1.0 + pow(CONSTANTS[AC_KmNaip]     / s_naj, 4.0)) * CONSTANTS[AC_ko] / (CONSTANTS[AC_ko] + CONSTANTS[AC_KmKo]);
        const double a_INaK_sl_PKA     = CONSTANTS[AC_Fsl]   * CONSTANTS[AC_IbarNaK] * a_fnak / (1.0 + pow(CONSTANTS[AC_KmNaip_PKA] / s_nasl, 4.0)) * CONSTANTS[AC_ko] / (CONSTANTS[AC_ko] + CONSTANTS[AC_KmKo]);
        const double a_INaK_sl_noPKA   = CONSTANTS[AC_Fsl]   * CONSTANTS[AC_IbarNaK] * a_fnak / (1.0 + pow(CONSTANTS[AC_KmNaip]     / s_nasl, 4.0)) * CONSTANTS[AC_ko] / (CONSTANTS[AC_ko] + CONSTANTS[AC_KmKo]);
        const double a_INaLj  = CONSTANTS[AC_Fjunc] * a_GNaL * (s_v - a_ENaj)  * s_mL * ((1.0 - a_fINaLp) * s_hL + a_fINaLp * s_hLp);
        const double a_INaLsl = CONSTANTS[AC_Fsl]   * a_GNaL * (s_v - a_ENasl) * s_mL * ((1.0 - a_fINaLp) * s_hL + a_fINaLp * s_hLp);
        const double a_INabj  = CONSTANTS[AC_Fjunc] * CONSTANTS[AC_GNab] * (s_v - a_ENaj);
        const double a_INabsl = CONSTANTS[AC_Fsl]   * CONSTANTS[AC_GNab] * (s_v - a_ENasl);
        const double a_IpCa_junc = CONSTANTS[AC_Fjunc] * pow(CONSTANTS[AC_Q10SLCaP], CONSTANTS[AC_Qpow]) * CONSTANTS[AC_IbarSLCaP] * pow(s_caj,  1.6) / (pow(CONSTANTS[AC_KmPCa], 1.6) + pow(s_caj,  1.6));
        const double a_IpCa_sl   = CONSTANTS[AC_Fsl]   * pow(CONSTANTS[AC_Q10SLCaP], CONSTANTS[AC_Qpow]) * CONSTANTS[AC_IbarSLCaP] * pow(s_casl, 1.6) / (pow(CONSTANTS[AC_KmPCa], 1.6) + pow(s_casl, 1.6));
        const double r_casr = a_Jserca - (a_J_SRleak * CONSTANTS[AC_vmyo] / CONSTANTS[AC_vsr] + a_J_SRCarel) - a_d_buffers_Csqn;
        const double r_cli  = a_ICl_tot * CONSTANTS[AC_Cmem] / (CONSTANTS[AC_vmyo] * CONSTANTS[AC_F]);

        const double a_ICaK_i_BP  = CONSTANTS[AC_PCaK_P]  * a_PhiCaK_i  * s_d_P * (a_fBP  * (1.0 - s_nca_i) + s_jca * a_fcaBP  * s_nca_i);
        const double a_ICaK_i_CaMK= CONSTANTS[AC_PCaKp]   * a_PhiCaK_i  * s_d   * (a_fp   * (1.0 - s_nca_i) + s_jca * a_fcap   * s_nca_i);
        const double a_ICaK_i_NP  = CONSTANTS[AC_PCaK]    * a_PhiCaK_i  * s_d   * (a_f    * (1.0 - s_nca_i) + s_jca * a_fca    * s_nca_i);
        const double a_ICaK_i_PKA = CONSTANTS[AC_PCaK_P]  * a_PhiCaK_i  * s_d_P * (a_f_P  * (1.0 - s_nca_i) + s_jca * a_fcap_P * s_nca_i);
        const double a_ICaK_ss_BP = CONSTANTS[AC_PCaK_P]  * a_PhiCaK_ss * s_d_P * (a_fBP  * (1.0 - s_nca)   + s_jca * a_fcaBP  * s_nca);
        const double a_ICaK_ss_CaMK=CONSTANTS[AC_PCaKp]   * a_PhiCaK_ss * s_d   * (a_fp   * (1.0 - s_nca)   + s_jca * a_fcap   * s_nca);
        const double a_ICaK_ss_NP = CONSTANTS[AC_PCaK]    * a_PhiCaK_ss * s_d   * (a_f    * (1.0 - s_nca)   + s_jca * a_fca    * s_nca);
        const double a_ICaK_ss_PKA= CONSTANTS[AC_PCaK_P]  * a_PhiCaK_ss * s_d_P * (a_f_P  * (1.0 - s_nca)   + s_jca * a_fcap_P * s_nca);
        const double a_ICaL_i_BP  = CONSTANTS[AC_PCa_P]   * a_PhiCaL_i  * s_d_P * (a_fBP  * (1.0 - s_nca_i) + s_jca * a_fcaBP  * s_nca_i);
        const double a_ICaL_i_CaMK= CONSTANTS[AC_PCap]    * a_PhiCaL_i  * s_d   * (a_fp   * (1.0 - s_nca_i) + s_jca * a_fcap   * s_nca_i);
        const double a_ICaL_i_NP  = CONSTANTS[AC_PCa]     * a_PhiCaL_i  * s_d   * (a_f    * (1.0 - s_nca_i) + s_jca * a_fca    * s_nca_i);
        const double a_ICaL_i_PKA = CONSTANTS[AC_PCa_P]   * a_PhiCaL_i  * s_d_P * (a_f_P  * (1.0 - s_nca_i) + s_jca * a_fcap_P * s_nca_i);
        const double a_ICaL_ss_BP = CONSTANTS[AC_PCa_P]   * a_PhiCaL_ss * s_d_P * (a_fBP  * (1.0 - s_nca)   + s_jca * a_fcaBP  * s_nca);
        const double a_ICaL_ss_CaMK=CONSTANTS[AC_PCap]    * a_PhiCaL_ss * s_d   * (a_fp   * (1.0 - s_nca)   + s_jca * a_fcap   * s_nca);
        const double a_ICaL_ss_NP = CONSTANTS[AC_PCa]     * a_PhiCaL_ss * s_d   * (a_f    * (1.0 - s_nca)   + s_jca * a_fca    * s_nca);
        const double a_ICaL_ss_PKA= CONSTANTS[AC_PCa_P]   * a_PhiCaL_ss * s_d_P * (a_f_P  * (1.0 - s_nca)   + s_jca * a_fcap_P * s_nca);
        const double a_ICaNa_i_BP = CONSTANTS[AC_PCaNa_P]  * a_PhiCaNa_i  * s_d_P * (a_fBP  * (1.0 - s_nca_i) + s_jca * a_fcaBP  * s_nca_i);
        const double a_ICaNa_i_CaMK=CONSTANTS[AC_PCaNap]   * a_PhiCaNa_i  * s_d   * (a_fp   * (1.0 - s_nca_i) + s_jca * a_fcap   * s_nca_i);
        const double a_ICaNa_i_NP = CONSTANTS[AC_PCaNa]    * a_PhiCaNa_i  * s_d   * (a_f    * (1.0 - s_nca_i) + s_jca * a_fca    * s_nca_i);
        const double a_ICaNa_i_PKA= CONSTANTS[AC_PCaNa_P]  * a_PhiCaNa_i  * s_d_P * (a_f_P  * (1.0 - s_nca_i) + s_jca * a_fcap_P * s_nca_i);
        const double a_ICaNa_ss_BP= CONSTANTS[AC_PCaNa_P]  * a_PhiCaNa_ss * s_d_P * (a_fBP  * (1.0 - s_nca)   + s_jca * a_fcaBP  * s_nca);
        const double a_ICaNa_ss_CaMK=CONSTANTS[AC_PCaNap]  * a_PhiCaNa_ss * s_d   * (a_fp   * (1.0 - s_nca)   + s_jca * a_fcap   * s_nca);
        const double a_ICaNa_ss_NP= CONSTANTS[AC_PCaNa]    * a_PhiCaNa_ss * s_d   * (a_f    * (1.0 - s_nca)   + s_jca * a_fca    * s_nca);
        const double a_ICaNa_ss_PKA=CONSTANTS[AC_PCaNa_P]  * a_PhiCaNa_ss * s_d_P * (a_f_P  * (1.0 - s_nca)   + s_jca * a_fcap_P * s_nca);
        const double a_beta_ItoC2  = a_beta_2 * a_beta_i * a_alpha_C2ToI / (a_alpha_2 * a_alpha_i);
        const double r_C0  = a_IKr_beta * s_C1 - a_IKr_alpha * s_C0;
        const double r_C1  = a_IKr_alpha * s_C0 + CONSTANTS[AC_IKr_beta_1] * s_C2 - (a_IKr_beta + CONSTANTS[AC_alpha_1]) * s_C1;
        const double r_O   = a_alpha_2 * s_C2 + a_beta_i * s_I - (a_beta_2 + a_alpha_i) * s_O;
        const double a_INaKj  = (1.0 - CONSTANTS[AC_fINaK_PKA]) * a_INaK_junc_noPKA + CONSTANTS[AC_fINaK_PKA] * a_INaK_junc_PKA;
        const double a_INaKsl = (1.0 - CONSTANTS[AC_fINaK_PKA]) * a_INaK_sl_noPKA   + CONSTANTS[AC_fINaK_PKA] * a_INaK_sl_PKA;
        const double a_INaK_INaK = a_INaKj + a_INaKsl;

        const double a_ICaK_junc = ((1.0 - a_fICaL_CaMKonly - a_fICaL_PKAonly - a_fICaL_BP) * a_ICaK_ss_NP   + a_fICaL_CaMKonly * a_ICaK_ss_CaMK + a_fICaL_PKAonly * a_ICaK_ss_PKA + a_fICaL_BP * a_ICaK_ss_BP) * s_ical_pureCDI_junc * CONSTANTS[AC_ICaL_fractionSS];
        const double a_ICaK_sl   = ((1.0 - a_fICaL_CaMKonly - a_fICaL_PKAonly - a_fICaL_BP) * a_ICaK_i_NP    + a_fICaL_CaMKonly * a_ICaK_i_CaMK  + a_fICaL_PKAonly * a_ICaK_i_PKA  + a_fICaL_BP * a_ICaK_i_BP)  * s_ical_pureCDI_sl   * (1.0 - CONSTANTS[AC_ICaL_fractionSS]);
        const double a_ICaL_junc = ((1.0 - a_fICaL_CaMKonly - a_fICaL_PKAonly - a_fICaL_BP) * a_ICaL_ss_NP   + a_fICaL_CaMKonly * a_ICaL_ss_CaMK + a_fICaL_PKAonly * a_ICaL_ss_PKA + a_fICaL_BP * a_ICaL_ss_BP) * s_ical_pureCDI_junc * CONSTANTS[AC_ICaL_fractionSS];
        const double a_ICaL_sl   = ((1.0 - a_fICaL_CaMKonly - a_fICaL_PKAonly - a_fICaL_BP) * a_ICaL_i_NP    + a_fICaL_CaMKonly * a_ICaL_i_CaMK  + a_fICaL_PKAonly * a_ICaL_i_PKA  + a_fICaL_BP * a_ICaL_i_BP)  * s_ical_pureCDI_sl   * (1.0 - CONSTANTS[AC_ICaL_fractionSS]);
        const double a_ICaNa_junc= ((1.0 - a_fICaL_CaMKonly - a_fICaL_PKAonly - a_fICaL_BP) * a_ICaNa_ss_NP  + a_fICaL_CaMKonly * a_ICaNa_ss_CaMK+ a_fICaL_PKAonly * a_ICaNa_ss_PKA+ a_fICaL_BP * a_ICaNa_ss_BP)* s_ical_pureCDI_junc * CONSTANTS[AC_ICaL_fractionSS];
        const double a_ICaNa_sl  = ((1.0 - a_fICaL_CaMKonly - a_fICaL_PKAonly - a_fICaL_BP) * a_ICaNa_i_NP   + a_fICaL_CaMKonly * a_ICaNa_i_CaMK + a_fICaL_PKAonly * a_ICaNa_i_PKA + a_fICaL_BP * a_ICaNa_i_BP)  * s_ical_pureCDI_sl   * (1.0 - CONSTANTS[AC_ICaL_fractionSS]);
        const double r_C2  = CONSTANTS[AC_alpha_1] * s_C1 + a_beta_2 * s_O + a_beta_ItoC2 * s_I - (CONSTANTS[AC_IKr_beta_1] + a_alpha_2 + a_alpha_C2ToI) * s_C2;
        const double r_I   = a_alpha_C2ToI * s_C2 + a_alpha_i * s_O - (a_beta_ItoC2 + a_beta_i) * s_I;
        const double a_ICaK         = a_ICaK_junc + a_ICaK_sl;
        const double a_ICa_tot_junc = a_ICaL_junc + a_ICab_junc + a_IpCa_junc - 2.0 * a_INaCa_ss;
        const double a_ICa_tot_sl   = a_ICaL_sl   + a_ICab_sl   + a_IpCa_sl   - 2.0 * a_INaCa_i;
        const double a_INa_tot_junc = a_INaj + a_INaLj + a_INabj + 3.0 * a_INaCa_ss + 3.0 * a_INaKj  + a_ICaNa_junc;
        const double a_INa_tot_sl   = a_INasl + a_INaLsl + a_INabsl + 3.0 * a_INaCa_i + 3.0 * a_INaKsl + a_ICaNa_sl;
        const double a_ICaL_junc_positive = fabs(a_ICaL_junc);
        const double a_IK_tot       = a_Ito_Ito + a_IKr_IKr + a_IKs_IKs + a_IK1_IK1 - 2.0 * a_INaK_INaK + a_ICaK + a_IKb_IKb;
        const double a_ICaL_junc_sigmoided = 1.0 - 1.0 / (1.0 + pow(a_ICaL_junc_positive / 0.45, 4.5));
        const double r_ki   = -a_IK_tot * CONSTANTS[AC_Cmem] / (CONSTANTS[AC_vmyo] * CONSTANTS[AC_F]);
        const double a_Iion_cm = a_INa_tot_junc + a_INa_tot_sl + a_ICa_tot_junc + a_ICa_tot_sl + a_IK_tot + a_ICl_tot;
        const double a_Jrel_inact_inf  = 1.0 / (1.0 + a_ICaL_junc_sigmoided / 0.001);
        const double a_Jrel_inact_inf2 = 1.0 / (1.0 + a_ICaL_junc_sigmoided / 0.0006);
        const double a_Jrel_inf        = CONSTANTS[AC_a_rel] * a_ICaL_junc_sigmoided / (1.0 + pow(CONSTANTS[AC_directRelMidpoint] / s_casr, 7.72672));
        const double r_v    = solveVmWithinODESolver ? (-a_Iion_cm - Istim) : 0.0;
        const double r_jrel_icaldep_act = (a_Jrel_inf        - s_jrel_icaldep_act) / a_tau_rel;
        const double r_jrel_icaldep_f1  = (a_Jrel_inact_inf  - s_jrel_icaldep_f1)  / CONSTANTS[AC_tauInact];
        const double r_jrel_icaldep_f2  = (a_Jrel_inact_inf2 - s_jrel_icaldep_f2)  / CONSTANTS[AC_tauInact2];
        const double r_caj  = -a_ICa_tot_junc * CONSTANTS[AC_Cmem] / (CONSTANTS[AC_vjunc] * 2.0 * CONSTANTS[AC_F]) + CONSTANTS[AC_J_ca_juncsl] / CONSTANTS[AC_vjunc] * (s_casl - s_caj) - a_J_CaB_junction + a_J_SRCarel * CONSTANTS[AC_vsr] / CONSTANTS[AC_vjunc] + a_J_SRleak * CONSTANTS[AC_vmyo] / CONSTANTS[AC_vjunc];
        const double r_casl = -a_ICa_tot_sl   * CONSTANTS[AC_Cmem] / (CONSTANTS[AC_vsl]   * 2.0 * CONSTANTS[AC_F]) + CONSTANTS[AC_J_ca_juncsl] / CONSTANTS[AC_vsl]   * (s_caj - s_casl) + CONSTANTS[AC_J_ca_slmyo] / CONSTANTS[AC_vsl] * (s_cai - s_casl) - a_J_CaB_sl;
        const double r_naj  = -a_INa_tot_junc * CONSTANTS[AC_Cmem] / (CONSTANTS[AC_vjunc] * CONSTANTS[AC_F]) + CONSTANTS[AC_J_na_juncsl] / CONSTANTS[AC_vjunc] * (s_nasl - s_naj) - a_d_buffers_NaBj;
        const double r_nasl = -a_INa_tot_sl   * CONSTANTS[AC_Cmem] / (CONSTANTS[AC_vsl]   * CONSTANTS[AC_F]) + CONSTANTS[AC_J_na_juncsl] / CONSTANTS[AC_vsl]   * (s_naj - s_nasl) + CONSTANTS[AC_J_na_slmyo] / CONSTANTS[AC_vsl] * (s_nai - s_nasl) - a_d_buffers_NaBsl;

        /* write rates to SoA */
        R(v)                      = r_v;
        R(camk_trap)              = r_camk_trap;
        R(camk_f_ICaL)            = r_camk_f_ICaL;
        R(camk_f_RyR)             = r_camk_f_RyR;
        R(camk_f_PLB)             = r_camk_f_PLB;
        R(casig_serca_trap)       = r_casig_serca_trap;
        R(buffers_NaBj)           = a_d_buffers_NaBj;
        R(buffers_NaBsl)          = a_d_buffers_NaBsl;
        R(buffers_TnClow)         = a_d_buffers_TnClow;
        R(buffers_TnCHc)          = a_d_buffers_TnCHc;
        R(buffers_TnCHm)          = r_buffers_TnCHm;
        R(buffers_CaM)            = a_d_buffers_CaM;
        R(buffers_Myosin_ca)      = a_d_buffers_Myosin_ca;
        R(buffers_Myosin_mg)      = r_buffers_Myosin_mg;
        R(buffers_SRB)            = a_d_buffers_SRB;
        R(buffers_SLLj)           = a_d_buffers_SLLj;
        R(buffers_SLLsl)          = a_d_buffers_SLLsl;
        R(buffers_SLHj)           = a_d_buffers_SLHj;
        R(buffers_SLHsl)          = a_d_buffers_SLHsl;
        R(buffers_Csqn)           = a_d_buffers_Csqn;
        R(naj)                    = r_naj;
        R(nasl)                   = r_nasl;
        R(nai)                    = r_nai;
        R(ki)                     = r_ki;
        R(cli)                    = r_cli;
        R(casr)                   = r_casr;
        R(caj)                    = r_caj;
        R(casl)                   = r_casl;
        R(cai)                    = r_cai;
        R(m)                      = r_m;
        R(h)                      = r_h;
        R(j)                      = r_j;
        R(hp)                     = r_hp;
        R(jp)                     = r_jp;
        R(m_P)                    = r_m_P;
        R(h_P)                    = r_h_P;
        R(j_P)                    = r_j_P;
        R(hp_P)                   = r_hp_P;
        R(jp_P)                   = r_jp_P;
        R(mL)                     = r_mL;
        R(hL)                     = r_hL;
        R(hLp)                    = r_hLp;
        R(d)                      = r_d;
        R(ff)                     = r_ff;
        R(fs)                     = r_fs;
        R(fcaf)                   = r_fcaf;
        R(fcas)                   = r_fcas;
        R(jca)                    = r_jca;
        R(ffp)                    = r_ffp;
        R(fcafp)                  = r_fcafp;
        R(nca)                    = r_nca;
        R(nca_i)                  = r_nca_i;
        R(d_P)                    = r_d_P;
        R(ff_P)                   = r_ff_P;
        R(fs_P)                   = r_fs_P;
        R(fcaf_P)                 = r_fcaf_P;
        R(fcas_P)                 = r_fcas_P;
        R(fBPf)                   = r_fBPf;
        R(fcaBPf)                 = r_fcaBPf;
        R(ical_pureCDI_junc)      = r_ical_pureCDI_junc;
        R(ical_pureCDI_sl)        = r_ical_pureCDI_sl;
        R(xtos)                   = r_xtos;
        R(ytos)                   = r_ytos;
        R(xtof)                   = r_xtof;
        R(ytof)                   = r_ytof;
        R(xtos_p)                 = r_xtos_p;
        R(xtof_p)                 = r_xtof_p;
        R(ytos_p)                 = r_ytos_p;
        R(ytof_p)                 = r_ytof_p;
        R(C0)                     = r_C0;
        R(C1)                     = r_C1;
        R(C2)                     = r_C2;
        R(::I)                    = r_I;
        R(::O)                    = r_O;
        R(xs_junc)                = r_xs_junc;
        R(xs_sl)                  = r_xs_sl;
        R(jrel_icaldep_act)       = r_jrel_icaldep_act;
        R(jrel_icaldep_f1)        = r_jrel_icaldep_f1;
        R(jrel_icaldep_f2)        = r_jrel_icaldep_f2;
        R(ryr_R)                  = r_ryr_R;
        R(ryr_O)                  = r_ryr_O;
        R(ryr_I)                  = r_ryr_I;
        R(ryr_CaRI)               = r_ryr_CaRI;
        R(ryr_R_p)                = r_ryr_R_p;
        R(ryr_O_p)                = r_ryr_O_p;
        R(ryr_I_p)                = r_ryr_I_p;
        R(ryr_CaRI_p)             = r_ryr_CaRI_p;
        R(contraction_TmBlocked)  = r_contraction_TmBlocked;
        R(contraction_XW)         = r_contraction_XW;
        R(contraction_XS)         = r_contraction_XS;
        R(contraction_ZETAS)      = r_contraction_ZETAS;
        R(contraction_ZETAW)      = r_contraction_ZETAW;
        R(contraction_Ca_TRPN)    = r_contraction_Ca_TRPN;

        /* write SUPPORT — tau/gInf for RL gates, then Iion_cm */
        SUP(TWORLD_BATCH_SUPPORT_tau_camk_ICaL)  = CONSTANTS[AC_tau_cal];
        SUP(TWORLD_BATCH_SUPPORT_gInf_camk_ICaL) = a_CaMK_Phos_ss_ICaL;
        SUP(TWORLD_BATCH_SUPPORT_tau_camk_PLB)   = CONSTANTS[AC_tau_plb];
        SUP(TWORLD_BATCH_SUPPORT_gInf_camk_PLB)  = a_CaMK_Phos_ss_PLB;
        SUP(TWORLD_BATCH_SUPPORT_tau_camk_RyR)   = CONSTANTS[AC_tau_ryr];
        SUP(TWORLD_BATCH_SUPPORT_gInf_camk_RyR)  = a_CaMK_Phos_ss_RyR;
        SUP(TWORLD_BATCH_SUPPORT_tau_m)           = a_tm;
        SUP(TWORLD_BATCH_SUPPORT_gInf_m)          = a_mss;
        SUP(TWORLD_BATCH_SUPPORT_tau_m_P)         = a_tm;
        SUP(TWORLD_BATCH_SUPPORT_gInf_m_P)        = a_mss_P;
        SUP(TWORLD_BATCH_SUPPORT_tau_h)           = a_th;
        SUP(TWORLD_BATCH_SUPPORT_gInf_h)          = a_hss;
        SUP(TWORLD_BATCH_SUPPORT_tau_h_P)         = a_th;
        SUP(TWORLD_BATCH_SUPPORT_gInf_h_P)        = a_hss_P;
        SUP(TWORLD_BATCH_SUPPORT_tau_hp)          = a_th;
        SUP(TWORLD_BATCH_SUPPORT_gInf_hp)         = a_hssp;
        SUP(TWORLD_BATCH_SUPPORT_tau_hp_P)        = a_th;
        SUP(TWORLD_BATCH_SUPPORT_gInf_hp_P)       = a_hssp_P;
        SUP(TWORLD_BATCH_SUPPORT_tau_j)           = a_tj;
        SUP(TWORLD_BATCH_SUPPORT_gInf_j)          = a_jss;
        SUP(TWORLD_BATCH_SUPPORT_tau_j_P)         = a_tj;
        SUP(TWORLD_BATCH_SUPPORT_gInf_j_P)        = a_jss_P;
        SUP(TWORLD_BATCH_SUPPORT_tau_jp)          = a_tjp;
        SUP(TWORLD_BATCH_SUPPORT_gInf_jp)         = a_jss;
        SUP(TWORLD_BATCH_SUPPORT_tau_jp_P)        = a_tjp;
        SUP(TWORLD_BATCH_SUPPORT_gInf_jp_P)       = a_jssp_P;
        SUP(TWORLD_BATCH_SUPPORT_tau_mL)          = a_tmL;
        SUP(TWORLD_BATCH_SUPPORT_gInf_mL)         = a_mLss;
        SUP(TWORLD_BATCH_SUPPORT_tau_hL)          = CONSTANTS[AC_thL];
        SUP(TWORLD_BATCH_SUPPORT_gInf_hL)         = a_hLss;
        SUP(TWORLD_BATCH_SUPPORT_tau_hLp)         = CONSTANTS[AC_thLp];
        SUP(TWORLD_BATCH_SUPPORT_gInf_hLp)        = a_hLssp;
        SUP(TWORLD_BATCH_SUPPORT_tau_d)           = a_td;
        SUP(TWORLD_BATCH_SUPPORT_gInf_d)          = a_dss;
        SUP(TWORLD_BATCH_SUPPORT_tau_ff)          = a_tff;
        SUP(TWORLD_BATCH_SUPPORT_gInf_ff)         = a_fss;
        SUP(TWORLD_BATCH_SUPPORT_tau_fs)          = a_tfs;
        SUP(TWORLD_BATCH_SUPPORT_gInf_fs)         = a_fss;
        SUP(TWORLD_BATCH_SUPPORT_tau_fcaf)        = a_tfcaf;
        SUP(TWORLD_BATCH_SUPPORT_gInf_fcaf)       = a_fcass;
        SUP(TWORLD_BATCH_SUPPORT_tau_fcas)        = a_tfcas;
        SUP(TWORLD_BATCH_SUPPORT_gInf_fcas)       = a_fcass;
        SUP(TWORLD_BATCH_SUPPORT_tau_jca)         = CONSTANTS[AC_tjca];
        SUP(TWORLD_BATCH_SUPPORT_gInf_jca)        = a_jcass;
        SUP(TWORLD_BATCH_SUPPORT_tau_ffp)         = a_tffp;
        SUP(TWORLD_BATCH_SUPPORT_gInf_ffp)        = a_fss;
        SUP(TWORLD_BATCH_SUPPORT_tau_fcafp)       = a_tfcafp;
        SUP(TWORLD_BATCH_SUPPORT_gInf_fcafp)      = a_fcass;
        SUP(TWORLD_BATCH_SUPPORT_tau_d_P)         = a_td;
        SUP(TWORLD_BATCH_SUPPORT_gInf_d_P)        = a_dPss;
        SUP(TWORLD_BATCH_SUPPORT_tau_ff_P)        = a_tff;
        SUP(TWORLD_BATCH_SUPPORT_gInf_ff_P)       = a_fss_P;
        SUP(TWORLD_BATCH_SUPPORT_tau_fs_P)        = a_tfs;
        SUP(TWORLD_BATCH_SUPPORT_gInf_fs_P)       = a_fss_P;
        SUP(TWORLD_BATCH_SUPPORT_tau_fcaf_P)      = a_tfcaf;
        SUP(TWORLD_BATCH_SUPPORT_gInf_fcaf_P)     = a_fcass_P;
        SUP(TWORLD_BATCH_SUPPORT_tau_fcas_P)      = a_tfcas;
        SUP(TWORLD_BATCH_SUPPORT_gInf_fcas_P)     = a_fcass_P;
        SUP(TWORLD_BATCH_SUPPORT_tau_fBPf)        = a_tffp;
        SUP(TWORLD_BATCH_SUPPORT_gInf_fBPf)       = a_fBPss;
        SUP(TWORLD_BATCH_SUPPORT_tau_fcaBPf)      = a_tfcafp;
        SUP(TWORLD_BATCH_SUPPORT_gInf_fcaBPf)     = a_fcaBPss;
        SUP(TWORLD_BATCH_SUPPORT_tau_xtos)        = a_tauxtos;
        SUP(TWORLD_BATCH_SUPPORT_gInf_xtos)       = a_xtoss;
        SUP(TWORLD_BATCH_SUPPORT_tau_xtos_p)      = a_tauxtos;
        SUP(TWORLD_BATCH_SUPPORT_gInf_xtos_p)     = a_xtoss_p;
        SUP(TWORLD_BATCH_SUPPORT_tau_xtof)        = a_tauxtof;
        SUP(TWORLD_BATCH_SUPPORT_gInf_xtof)       = a_xtoss;
        SUP(TWORLD_BATCH_SUPPORT_tau_xtof_p)      = a_tauxtof;
        SUP(TWORLD_BATCH_SUPPORT_gInf_xtof_p)     = a_xtoss_p;
        SUP(TWORLD_BATCH_SUPPORT_tau_ytos)        = a_tauytos;
        SUP(TWORLD_BATCH_SUPPORT_gInf_ytos)       = a_ytoss;
        SUP(TWORLD_BATCH_SUPPORT_tau_ytos_p)      = a_tauytos_p;
        SUP(TWORLD_BATCH_SUPPORT_gInf_ytos_p)     = a_ytoss;
        SUP(TWORLD_BATCH_SUPPORT_tau_ytof)        = a_tauytof;
        SUP(TWORLD_BATCH_SUPPORT_gInf_ytof)       = a_ytoss;
        SUP(TWORLD_BATCH_SUPPORT_tau_ytof_p)      = a_tauytof_p;
        SUP(TWORLD_BATCH_SUPPORT_gInf_ytof_p)     = a_ytoss;
        SUP(TWORLD_BATCH_SUPPORT_tau_xs_junc)     = a_tauxs_junc;
        SUP(TWORLD_BATCH_SUPPORT_gInf_xs_junc)    = a_xsss_junc;
        SUP(TWORLD_BATCH_SUPPORT_tau_xs_sl)       = a_tauxs_sl;
        SUP(TWORLD_BATCH_SUPPORT_gInf_xs_sl)      = a_xsss_sl;
        SUP(TWORLD_BATCH_SUPPORT_tau_jrel_act)    = a_tau_rel;
        SUP(TWORLD_BATCH_SUPPORT_gInf_jrel_act)   = a_Jrel_inf;
        SUP(TWORLD_BATCH_SUPPORT_tau_jrel_f1)     = CONSTANTS[AC_tauInact];
        SUP(TWORLD_BATCH_SUPPORT_gInf_jrel_f1)    = a_Jrel_inact_inf;
        SUP(TWORLD_BATCH_SUPPORT_tau_jrel_f2)     = CONSTANTS[AC_tauInact2];
        SUP(TWORLD_BATCH_SUPPORT_gInf_jrel_f2)    = a_Jrel_inact_inf2;
        SUP(TWORLD_BATCH_SUPPORT_Iion_cm)         = a_Iion_cm;
    }

    #undef S
    #undef R
    #undef SUP
}

} // End namespace Foam
