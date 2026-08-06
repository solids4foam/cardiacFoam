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
    State-major (SoA) batched evaluator for the ToRORd_dynCl 2023 ionic model.
Author
    Simao Nieto de Castro, UCD.
\*---------------------------------------------------------------------------*/

#pragma once

#include "ToRORd_dynCl_2023.H"
#include "ToRORd_dynCl_2023Names.H"
#include "StimulusProtocolPOD.H"
#include <cmath>

namespace Foam
{

enum TORORD_DYNCL_BATCH_SUPPORT_INDEX
{
    TORORD_DYNCL_BATCH_SUPPORT_tau_INa_m,
    TORORD_DYNCL_BATCH_SUPPORT_gInf_INa_m,
    TORORD_DYNCL_BATCH_SUPPORT_tau_INa_h,
    TORORD_DYNCL_BATCH_SUPPORT_gInf_INa_h,
    TORORD_DYNCL_BATCH_SUPPORT_tau_INa_j,
    TORORD_DYNCL_BATCH_SUPPORT_gInf_INa_j,
    TORORD_DYNCL_BATCH_SUPPORT_tau_INa_hp,
    TORORD_DYNCL_BATCH_SUPPORT_gInf_INa_hp,
    TORORD_DYNCL_BATCH_SUPPORT_tau_INa_jp,
    TORORD_DYNCL_BATCH_SUPPORT_gInf_INa_jp,
    TORORD_DYNCL_BATCH_SUPPORT_tau_INaL_mL,
    TORORD_DYNCL_BATCH_SUPPORT_gInf_INaL_mL,
    TORORD_DYNCL_BATCH_SUPPORT_tau_INaL_hL,
    TORORD_DYNCL_BATCH_SUPPORT_gInf_INaL_hL,
    TORORD_DYNCL_BATCH_SUPPORT_tau_INaL_hLp,
    TORORD_DYNCL_BATCH_SUPPORT_gInf_INaL_hLp,
    TORORD_DYNCL_BATCH_SUPPORT_tau_Ito_a,
    TORORD_DYNCL_BATCH_SUPPORT_gInf_Ito_a,
    TORORD_DYNCL_BATCH_SUPPORT_tau_Ito_iF,
    TORORD_DYNCL_BATCH_SUPPORT_gInf_Ito_iF,
    TORORD_DYNCL_BATCH_SUPPORT_tau_Ito_iS,
    TORORD_DYNCL_BATCH_SUPPORT_gInf_Ito_iS,
    TORORD_DYNCL_BATCH_SUPPORT_tau_Ito_ap,
    TORORD_DYNCL_BATCH_SUPPORT_gInf_Ito_ap,
    TORORD_DYNCL_BATCH_SUPPORT_tau_Ito_iFp,
    TORORD_DYNCL_BATCH_SUPPORT_gInf_Ito_iFp,
    TORORD_DYNCL_BATCH_SUPPORT_tau_Ito_iSp,
    TORORD_DYNCL_BATCH_SUPPORT_gInf_Ito_iSp,
    TORORD_DYNCL_BATCH_SUPPORT_tau_ICaL_d,
    TORORD_DYNCL_BATCH_SUPPORT_gInf_ICaL_d,
    TORORD_DYNCL_BATCH_SUPPORT_tau_ICaL_ff,
    TORORD_DYNCL_BATCH_SUPPORT_gInf_ICaL_ff,
    TORORD_DYNCL_BATCH_SUPPORT_tau_ICaL_fs,
    TORORD_DYNCL_BATCH_SUPPORT_gInf_ICaL_fs,
    TORORD_DYNCL_BATCH_SUPPORT_tau_ICaL_fcaf,
    TORORD_DYNCL_BATCH_SUPPORT_gInf_ICaL_fcaf,
    TORORD_DYNCL_BATCH_SUPPORT_tau_ICaL_fcas,
    TORORD_DYNCL_BATCH_SUPPORT_gInf_ICaL_fcas,
    TORORD_DYNCL_BATCH_SUPPORT_tau_ICaL_jca,
    TORORD_DYNCL_BATCH_SUPPORT_gInf_ICaL_jca,
    TORORD_DYNCL_BATCH_SUPPORT_tau_ICaL_ffp,
    TORORD_DYNCL_BATCH_SUPPORT_gInf_ICaL_ffp,
    TORORD_DYNCL_BATCH_SUPPORT_tau_ICaL_fcafp,
    TORORD_DYNCL_BATCH_SUPPORT_gInf_ICaL_fcafp,
    TORORD_DYNCL_BATCH_SUPPORT_tau_IKs_xs1,
    TORORD_DYNCL_BATCH_SUPPORT_gInf_IKs_xs1,
    TORORD_DYNCL_BATCH_SUPPORT_tau_IKs_xs2,
    TORORD_DYNCL_BATCH_SUPPORT_gInf_IKs_xs2,
    TORORD_DYNCL_BATCH_SUPPORT_tau_Jrel_np,
    TORORD_DYNCL_BATCH_SUPPORT_gInf_Jrel_np,
    TORORD_DYNCL_BATCH_SUPPORT_tau_Jrel_p,
    TORORD_DYNCL_BATCH_SUPPORT_gInf_Jrel_p,
    TORORD_DYNCL_BATCH_SUPPORT_Iion_cm,

    NUM_TORORD_DYNCL_BATCH_SUPPORT
};

static_assert(NUM_TORORD_DYNCL_BATCH_SUPPORT == 53, "ToRORd_dynCl compact support layout drift");

CARDIAC_HOST_DEVICE inline void ToRORd_dynClComputeVariablesBatch
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

    const double c_CaMKo = CONSTANTS[AC_CaMKo];
    const double c_KmCaM = CONSTANTS[AC_KmCaM];
    const double c_KmCaMK = CONSTANTS[AC_KmCaMK];
    const double c_aCaMK = CONSTANTS[AC_aCaMK];
    const double c_bCaMK = CONSTANTS[AC_bCaMK];
    const double c_GpCa = CONSTANTS[AC_GpCa];
    const double c_KmCap = CONSTANTS[AC_KmCap];
    const double c_Acap = CONSTANTS[AC_Acap];
    const double c_vjsr = CONSTANTS[AC_vjsr];
    const double c_vmyo = CONSTANTS[AC_vmyo];
    const double c_vnsr = CONSTANTS[AC_vnsr];
    const double c_vss = CONSTANTS[AC_vss];
    const double c_tauCa = CONSTANTS[AC_tauCa];
    const double c_tauK = CONSTANTS[AC_tauK];
    const double c_tauNa = CONSTANTS[AC_tauNa];
    const double c_cao = CONSTANTS[AC_cao];
    const double c_clo = CONSTANTS[AC_clo];
    const double c_ko = CONSTANTS[AC_ko];
    const double c_nao = CONSTANTS[AC_nao];
    const double c_F = CONSTANTS[AC_F];
    const double c_R = CONSTANTS[AC_R];
    const double c_T = CONSTANTS[AC_T];
    const double c_zca = CONSTANTS[AC_zca];
    const double c_zcl = CONSTANTS[AC_zcl];
    const double c_zk = CONSTANTS[AC_zk];
    const double c_zna = CONSTANTS[AC_zna];
    const double c_Jup_b = CONSTANTS[AC_Jup_b];
    const double c_upScale = CONSTANTS[AC_upScale];
    const double c_PKNa = CONSTANTS[AC_PKNa];
    const double c_Fjunc = CONSTANTS[AC_Fjunc];
    const double c_GClCa = CONSTANTS[AC_GClCa];
    const double c_GClb = CONSTANTS[AC_GClb];
    const double c_KdClCa = CONSTANTS[AC_KdClCa];
    const double c_GK1 = CONSTANTS[AC_GK1];
    const double c_GKb = CONSTANTS[AC_GKb];
    const double c_GKs = CONSTANTS[AC_GKs];
    const double c_GNa = CONSTANTS[AC_GNa];
    const double c_thL = CONSTANTS[AC_thL];
    const double c_GNaL = CONSTANTS[AC_GNaL];
    const double c_thLp = CONSTANTS[AC_thLp];
    const double c_fkatp = CONSTANTS[AC_fkatp];
    const double c_gkatp = CONSTANTS[AC_gkatp];
    const double c_akik = CONSTANTS[AC_akik];
    const double c_bkik = CONSTANTS[AC_bkik];
    const double c_EKshift = CONSTANTS[AC_EKshift];

    const double c_deltaEpiAmp     = CONSTANTS[deltaEpiAmp];
    const double c_deltaEpiShift   = CONSTANTS[deltaEpiShift];
    const double c_deltaEpiScale   = CONSTANTS[deltaEpiScale];
    const double c_jrelTissueScale = CONSTANTS[jrelTissueScale];

    const double c_Gto = CONSTANTS[AC_Gto];
    const double c_Aff = CONSTANTS[AC_Aff];
    const double c_ICaL_fractionSS = CONSTANTS[AC_ICaL_fractionSS];
    const double c_Kmn = CONSTANTS[AC_Kmn];
    const double c_k2n = CONSTANTS[AC_k2n];
    const double c_offset = CONSTANTS[AC_offset];
    const double c_tjca = CONSTANTS[AC_tjca];
    const double c_vShift = CONSTANTS[AC_vShift];
    const double c_Afs = CONSTANTS[AC_Afs];
    const double c_PCa = CONSTANTS[AC_PCa];
    const double c_constA = CONSTANTS[AC_constA];
    const double c_PCaK = CONSTANTS[AC_PCaK];
    const double c_PCaNa = CONSTANTS[AC_PCaNa];
    const double c_PCap = CONSTANTS[AC_PCap];
    const double c_gamma_cao = CONSTANTS[AC_gamma_cao];
    const double c_gamma_ko = CONSTANTS[AC_gamma_ko];
    const double c_gamma_nao = CONSTANTS[AC_gamma_nao];
    const double c_PCaKp = CONSTANTS[AC_PCaKp];
    const double c_PCaNap = CONSTANTS[AC_PCaNap];
    const double c_PCab = CONSTANTS[AC_PCab];
    const double c_alpha_1 = CONSTANTS[AC_alpha_1];
    const double c_beta_1 = CONSTANTS[AC_beta_1];
    const double c_GKr = CONSTANTS[AC_GKr];
    const double c_INaCa_fractionSS = CONSTANTS[AC_INaCa_fractionSS];
    const double c_KmCaAct = CONSTANTS[AC_KmCaAct];
    const double c_kcaon = CONSTANTS[AC_kcaon];
    const double c_kna1 = CONSTANTS[AC_kna1];
    const double c_kna2 = CONSTANTS[AC_kna2];
    const double c_kna3 = CONSTANTS[AC_kna3];
    const double c_qca = CONSTANTS[AC_qca];
    const double c_qna = CONSTANTS[AC_qna];
    const double c_wca = CONSTANTS[AC_wca];
    const double c_wna = CONSTANTS[AC_wna];
    const double c_wnaca = CONSTANTS[AC_wnaca];
    const double c_Gncx = CONSTANTS[AC_Gncx];
    const double c_k2_i = CONSTANTS[AC_k2_i];
    const double c_k2_ss = CONSTANTS[AC_k2_ss];
    const double c_k5_i = CONSTANTS[AC_k5_i];
    const double c_k5_ss = CONSTANTS[AC_k5_ss];
    const double c_h11_i = CONSTANTS[AC_h11_i];
    const double c_h11_ss = CONSTANTS[AC_h11_ss];
    const double c_k1_i = CONSTANTS[AC_k1_i];
    const double c_k1_ss = CONSTANTS[AC_k1_ss];
    const double c_H = CONSTANTS[AC_H];
    const double c_Khp = CONSTANTS[AC_Khp];
    const double c_Kki = CONSTANTS[AC_Kki];
    const double c_Kko = CONSTANTS[AC_Kko];
    const double c_Kmgatp = CONSTANTS[AC_Kmgatp];
    const double c_Knai0 = CONSTANTS[AC_Knai0];
    const double c_Knao0 = CONSTANTS[AC_Knao0];
    const double c_Knap = CONSTANTS[AC_Knap];
    const double c_Kxkur = CONSTANTS[AC_Kxkur];
    const double c_MgATP = CONSTANTS[AC_MgATP];
    const double c_delta = CONSTANTS[AC_delta];
    const double c_eP = CONSTANTS[AC_eP];
    const double c_k1p = CONSTANTS[AC_k1p];
    const double c_k2m = CONSTANTS[AC_k2m];
    const double c_k3m = CONSTANTS[AC_k3m];
    const double c_k3p = CONSTANTS[AC_k3p];
    const double c_k4m = CONSTANTS[AC_k4m];
    const double c_Pnak = CONSTANTS[AC_Pnak];
    const double c_a2 = CONSTANTS[AC_a2];
    const double c_a4 = CONSTANTS[AC_a4];
    const double c_b1 = CONSTANTS[AC_b1];
    const double c_PNab = CONSTANTS[AC_PNab];
    const double c_BSLmax = CONSTANTS[AC_BSLmax];
    const double c_BSRmax = CONSTANTS[AC_BSRmax];
    const double c_KmBSL = CONSTANTS[AC_KmBSL];
    const double c_KmBSR = CONSTANTS[AC_KmBSR];
    const double c_csqnmax = CONSTANTS[AC_csqnmax];
    const double c_kmcmdn = CONSTANTS[AC_kmcmdn];
    const double c_kmcsqn = CONSTANTS[AC_kmcsqn];
    const double c_kmtrpn = CONSTANTS[AC_kmtrpn];
    const double c_trpnmax = CONSTANTS[AC_trpnmax];
    const double c_cmdnmax = CONSTANTS[AC_cmdnmax];
    const double c_Jrel_b = CONSTANTS[AC_Jrel_b];
    const double c_bt = CONSTANTS[AC_bt];
    const double c_cajsr_half = CONSTANTS[AC_cajsr_half];
    const double c_a_rel = CONSTANTS[AC_a_rel];
    const double c_btp = CONSTANTS[AC_btp];
    const double c_a_relp = CONSTANTS[AC_a_relp];

    const int V_BASE = V * N;
    const int CaMKt_BASE = CaMKt * N;
    const int Nai_BASE = Nai * N;
    const int Nass_BASE = Nass * N;
    const int Ki_BASE = Ki * N;
    const int Kss_BASE = Kss * N;
    const int Cass_BASE = Cass * N;
    const int Cansr_BASE = Cansr * N;
    const int Cajsr_BASE = Cajsr * N;
    const int Cai_BASE = Cai * N;
    const int Cli_BASE = Cli * N;
    const int Clss_BASE = Clss * N;
    const int INa_m_BASE = INa_m * N;
    const int INa_h_BASE = INa_h * N;
    const int INa_j_BASE = INa_j * N;
    const int INa_hp_BASE = INa_hp * N;
    const int INa_jp_BASE = INa_jp * N;
    const int INaL_mL_BASE = INaL_mL * N;
    const int INaL_hL_BASE = INaL_hL * N;
    const int INaL_hLp_BASE = INaL_hLp * N;
    const int Ito_a_BASE = Ito_a * N;
    const int Ito_iF_BASE = Ito_iF * N;
    const int Ito_iS_BASE = Ito_iS * N;
    const int Ito_ap_BASE = Ito_ap * N;
    const int Ito_iFp_BASE = Ito_iFp * N;
    const int Ito_iSp_BASE = Ito_iSp * N;
    const int ICaL_d_BASE = ICaL_d * N;
    const int ICaL_ff_BASE = ICaL_ff * N;
    const int ICaL_fs_BASE = ICaL_fs * N;
    const int ICaL_fcaf_BASE = ICaL_fcaf * N;
    const int ICaL_fcas_BASE = ICaL_fcas * N;
    const int ICaL_jca_BASE = ICaL_jca * N;
    const int ICaL_ffp_BASE = ICaL_ffp * N;
    const int ICaL_fcafp_BASE = ICaL_fcafp * N;
    const int ICaL_nca_ss_BASE = ICaL_nca_ss * N;
    const int ICaL_nca_i_BASE = ICaL_nca_i * N;
    const int IKr_C1_BASE = IKr_C1 * N;
    const int IKr_C2_BASE = IKr_C2 * N;
    const int IKr_C3_BASE = IKr_C3 * N;
    const int IKr_I_BASE = IKr_I * N;
    const int IKr_O_BASE = IKr_O * N;
    const int IKs_xs1_BASE = IKs_xs1 * N;
    const int IKs_xs2_BASE = IKs_xs2 * N;
    const int Jrel_np_BASE = Jrel_np * N;
    const int Jrel_p_BASE = Jrel_p * N;

#ifdef _OPENMP
    #pragma omp simd
#endif
    for (int cellI = beginCell; cellI < endCell; ++cellI)
    {
        const double s_V = STATES_SoA[V_BASE + cellI];
        const double s_CaMKt = STATES_SoA[CaMKt_BASE + cellI];
        const double s_Nai = STATES_SoA[Nai_BASE + cellI];
        const double s_Nass = STATES_SoA[Nass_BASE + cellI];
        const double s_Ki = STATES_SoA[Ki_BASE + cellI];
        const double s_Kss = STATES_SoA[Kss_BASE + cellI];
        const double s_Cass = STATES_SoA[Cass_BASE + cellI];
        const double s_Cansr = STATES_SoA[Cansr_BASE + cellI];
        const double s_Cajsr = STATES_SoA[Cajsr_BASE + cellI];
        const double s_Cai = STATES_SoA[Cai_BASE + cellI];
        const double s_Cli = STATES_SoA[Cli_BASE + cellI];
        const double s_Clss = STATES_SoA[Clss_BASE + cellI];
        const double s_INa_m = STATES_SoA[INa_m_BASE + cellI];
        const double s_INa_h = STATES_SoA[INa_h_BASE + cellI];
        const double s_INa_j = STATES_SoA[INa_j_BASE + cellI];
        const double s_INa_hp = STATES_SoA[INa_hp_BASE + cellI];
        const double s_INa_jp = STATES_SoA[INa_jp_BASE + cellI];
        const double s_INaL_mL = STATES_SoA[INaL_mL_BASE + cellI];
        const double s_INaL_hL = STATES_SoA[INaL_hL_BASE + cellI];
        const double s_INaL_hLp = STATES_SoA[INaL_hLp_BASE + cellI];
        const double s_Ito_a = STATES_SoA[Ito_a_BASE + cellI];
        const double s_Ito_iF = STATES_SoA[Ito_iF_BASE + cellI];
        const double s_Ito_iS = STATES_SoA[Ito_iS_BASE + cellI];
        const double s_Ito_ap = STATES_SoA[Ito_ap_BASE + cellI];
        const double s_Ito_iFp = STATES_SoA[Ito_iFp_BASE + cellI];
        const double s_Ito_iSp = STATES_SoA[Ito_iSp_BASE + cellI];
        const double s_ICaL_d = STATES_SoA[ICaL_d_BASE + cellI];
        const double s_ICaL_ff = STATES_SoA[ICaL_ff_BASE + cellI];
        const double s_ICaL_fs = STATES_SoA[ICaL_fs_BASE + cellI];
        const double s_ICaL_fcaf = STATES_SoA[ICaL_fcaf_BASE + cellI];
        const double s_ICaL_fcas = STATES_SoA[ICaL_fcas_BASE + cellI];
        const double s_ICaL_jca = STATES_SoA[ICaL_jca_BASE + cellI];
        const double s_ICaL_ffp = STATES_SoA[ICaL_ffp_BASE + cellI];
        const double s_ICaL_fcafp = STATES_SoA[ICaL_fcafp_BASE + cellI];
        const double s_ICaL_nca_ss = STATES_SoA[ICaL_nca_ss_BASE + cellI];
        const double s_ICaL_nca_i = STATES_SoA[ICaL_nca_i_BASE + cellI];
        const double s_IKr_C1 = STATES_SoA[IKr_C1_BASE + cellI];
        const double s_IKr_C2 = STATES_SoA[IKr_C2_BASE + cellI];
        const double s_IKr_C3 = STATES_SoA[IKr_C3_BASE + cellI];
        const double s_IKr_I = STATES_SoA[IKr_I_BASE + cellI];
        const double s_IKr_O = STATES_SoA[IKr_O_BASE + cellI];
        const double s_IKs_xs1 = STATES_SoA[IKs_xs1_BASE + cellI];
        const double s_IKs_xs2 = STATES_SoA[IKs_xs2_BASE + cellI];
        const double s_Jrel_np = STATES_SoA[Jrel_np_BASE + cellI];
        const double s_Jrel_p = STATES_SoA[Jrel_p_BASE + cellI];

        /* CaMK */
        const double a_CaMKb = c_CaMKo * (1.0 - s_CaMKt) / (1.0 + c_KmCaM / s_Cass);
        const double a_CaMKa = a_CaMKb + s_CaMKt;
        const double r_CaMKt = c_aCaMK * a_CaMKb * (a_CaMKb + s_CaMKt) - c_bCaMK * s_CaMKt;

        /* IpCa */
        const double a_IpCa_IpCa = c_GpCa * s_Cai / (c_KmCap + s_Cai);

        /* diff */
        const double a_Jdiff = (s_Cass - s_Cai) / c_tauCa;
        const double a_JdiffCl = (s_Clss - s_Cli) / c_tauNa;
        const double a_JdiffK = (s_Kss - s_Ki) / c_tauK;
        const double a_JdiffNa = (s_Nass - s_Nai) / c_tauNa;

        /* trans_flux */
        const double a_Jtr = (s_Cansr - s_Cajsr) / 60.0;

        /* SERCA */
        const double a_Jleak = 0.0048825 * s_Cansr / 15.0;
        const double a_fJupp = 1.0 / (1.0 + c_KmCaMK / a_CaMKa);
        const double a_Jupnp = c_upScale * 0.005425 * s_Cai / (s_Cai + 0.00092);
        const double a_Jupp = c_upScale * 2.75 * 0.005425 * s_Cai / (s_Cai + 0.00092 - 0.00017);
        const double a_Jup = c_Jup_b * ((1.0 - a_fJupp) * a_Jupnp + a_fJupp * a_Jupp - a_Jleak);

        /* reversal_potentials */
        const double a_ECl = c_R * c_T / (c_zcl * c_F) * std::log(c_clo / s_Cli);
        const double a_EClss = c_R * c_T / (c_zcl * c_F) * std::log(c_clo / s_Clss);
        const double a_EK = c_R * c_T / (c_zk * c_F) * std::log(c_ko / s_Ki);
        const double a_ENa = c_R * c_T / (c_zna * c_F) * std::log(c_nao / s_Nai);
        const double a_EKs = c_R * c_T / (c_zk * c_F) * std::log((c_ko + c_PKNa * c_nao) / (s_Ki + c_PKNa * s_Nai));

        /* ICl */
        const double a_IClCa_junc = c_Fjunc * c_GClCa / (1.0 + c_KdClCa / s_Cass) * (s_V - a_EClss);
        const double a_IClCa_sl = (1.0 - c_Fjunc) * c_GClCa / (1.0 + c_KdClCa / s_Cai) * (s_V - a_ECl);
        const double a_IClb = c_GClb * (s_V - a_ECl);
        const double a_IClCa = a_IClCa_junc + a_IClCa_sl;

        /* IK1 */
        const double a_aK1 = 4.094 / (1.0 + std::exp(0.1217 * (s_V - a_EK - 49.934)));
        const double a_bK1 = (15.72 * std::exp(0.0674 * (s_V - a_EK - 3.257)) + std::exp(0.0618 * (s_V - a_EK - 594.31))) / (1.0 + std::exp(-0.1629 * (s_V - a_EK + 14.207)));
        const double a_K1ss = a_aK1 / (a_aK1 + a_bK1);
        const double a_IK1_IK1 = c_GK1 * std::sqrt(c_ko / 5.0) * a_K1ss * (s_V - a_EK);

        /* IKb */
        const double a_xkb = 1.0 / (1.0 + std::exp(-(s_V - 10.8968) / 23.9871));
        const double a_IKb_IKb = c_GKb * a_xkb * (s_V - a_EK);

        /* IKs */
        const double a_KsCa = 1.0 + 0.6 / (1.0 + std::pow(3.8e-05 / s_Cai, 1.4));
        const double a_txs1 = 817.3 + 1.0 / (0.0002326 * std::exp((s_V + 48.28) / 17.8) + 0.001292 * std::exp(-(s_V + 210.0) / 230.0));
        const double a_txs2 = 1.0 / (0.01 * std::exp((s_V - 50.0) / 20.0) + 0.0193 * std::exp(-(s_V + 66.54) / 31.0));
        const double a_xs1ss = 1.0 / (1.0 + std::exp(-(s_V + 11.6) / 8.932));
        const double a_xs2ss = a_xs1ss;
        const double r_IKs_xs1 = (a_xs1ss - s_IKs_xs1) / a_txs1;
        const double a_IKs_IKs = c_GKs * a_KsCa * s_IKs_xs1 * s_IKs_xs2 * (s_V - a_EKs);
        const double r_IKs_xs2 = (a_xs2ss - s_IKs_xs2) / a_txs2;

        /* INa */
        const double a_ah = ((s_V >= -40.0) ? 0.0 : 0.057 * std::exp(-(s_V + 80.0) / 6.8));
        const double a_aj = ((s_V >= -40.0) ? 0.0 : (-25428.0 * std::exp(0.2444 * s_V) - 6.948e-06 * std::exp(-0.04391 * s_V)) * (s_V + 37.78) / (1.0 + std::exp(0.311 * (s_V + 79.23))));
        const double a_bh = ((s_V >= -40.0) ? 0.77 / (0.13 * (1.0 + std::exp(-(s_V + 10.66) / 11.1))) : 2.7 * std::exp(0.079 * s_V) + 310000.0 * std::exp(0.3485 * s_V));
        const double a_bj = ((s_V >= -40.0) ? 0.6 * std::exp(0.057 * s_V) / (1.0 + std::exp(-0.1 * (s_V + 32.0))) : 0.02424 * std::exp(-0.01052 * s_V) / (1.0 + std::exp(-0.1378 * (s_V + 40.14))));
        const double a_fINap = 1.0 / (1.0 + c_KmCaMK / a_CaMKa);
        const double a_hss = 1.0 / std::pow(1.0 + std::exp((s_V + 71.55) / 7.43), 2.0);
        const double a_hssp = 1.0 / std::pow(1.0 + std::exp((s_V + 77.55) / 7.43), 2.0);
        const double a_mss = 1.0 / std::pow(1.0 + std::exp(-(s_V + 56.86) / 9.03), 2.0);
        const double a_tm = 0.1292 * std::exp(-std::pow((s_V + 45.79) / 15.54, 2.0)) + 0.06487 * std::exp(-std::pow((s_V - 4.823) / 51.12, 2.0));
        const double a_INa_INa = c_GNa * (s_V - a_ENa) * std::pow(s_INa_m, 3.0) * ((1.0 - a_fINap) * s_INa_h * s_INa_j + a_fINap * s_INa_hp * s_INa_jp);
        const double a_jss = a_hss;
        const double a_th = 1.0 / (a_ah + a_bh);
        const double a_tj = 1.0 / (a_aj + a_bj);
        const double r_INa_m = (a_mss - s_INa_m) / a_tm;
        const double a_tjp = 1.46 * a_tj;
        const double r_INa_h = (a_hss - s_INa_h) / a_th;
        const double r_INa_hp = (a_hssp - s_INa_hp) / a_th;
        const double r_INa_j = (a_jss - s_INa_j) / a_tj;
        const double r_INa_jp = (a_jss - s_INa_jp) / a_tjp;

        /* INaL */
        const double a_fINaLp = 1.0 / (1.0 + c_KmCaMK / a_CaMKa);
        const double a_hLss = 1.0 / (1.0 + std::exp((s_V + 87.61) / 7.488));
        const double a_hLssp = 1.0 / (1.0 + std::exp((s_V + 93.81) / 7.488));
        const double a_mLss = 1.0 / (1.0 + std::exp(-(s_V + 42.85) / 5.264));
        const double a_tmL = 0.1292 * std::exp(-std::pow((s_V + 45.79) / 15.54, 2.0)) + 0.06487 * std::exp(-std::pow((s_V - 4.823) / 51.12, 2.0));
        const double r_INaL_hL = (a_hLss - s_INaL_hL) / c_thL;
        const double r_INaL_mL = (a_mLss - s_INaL_mL) / a_tmL;
        const double a_INaL_INaL = c_GNaL * (s_V - a_ENa) * s_INaL_mL * ((1.0 - a_fINaLp) * s_INaL_hL + a_fINaLp * s_INaL_hLp);
        const double r_INaL_hLp = (a_hLssp - s_INaL_hLp) / c_thLp;

        /* I_katp */
        const double a_I_katp_I_katp = c_fkatp * c_gkatp * c_akik * c_bkik * (s_V - a_EK);

        /* Ito */
        const double a_fItop = 1.0 / (1.0 + c_KmCaMK / a_CaMKa);
        const double a_AiF = 1.0 / (1.0 + std::exp((s_V + c_EKshift - 213.6) / 151.2));
        const double a_ass = 1.0 / (1.0 + std::exp(-(s_V + c_EKshift - 14.34) / 14.82));
        const double a_assp = 1.0 / (1.0 + std::exp(-(s_V + c_EKshift - 24.34) / 14.82));
        const double a_delta_epi = 1.0 - c_deltaEpiAmp
            / (1.0 + std::exp((s_V + c_EKshift + c_deltaEpiShift) / c_deltaEpiScale));
        const double a_dti_develop = 1.354 + 0.0001 / (std::exp((s_V + c_EKshift - 167.4) / 15.89) + std::exp(-(s_V + c_EKshift - 12.23) / 0.2154));
        const double a_dti_recover = 1.0 - 0.5 / (1.0 + std::exp((s_V + c_EKshift + 70.0) / 20.0));
        const double a_iss = 1.0 / (1.0 + std::exp((s_V + c_EKshift + 43.94) / 5.711));
        const double a_ta = 1.0515 / (1.0 / (1.2089 * (1.0 + std::exp(-(s_V + c_EKshift - 18.4099) / 29.3814))) + 3.5 / (1.0 + std::exp((s_V + c_EKshift + 100.0) / 29.3814)));
        const double a_tiF_b = 4.562 + 1.0 / (0.3933 * std::exp(-(s_V + c_EKshift + 100.0) / 100.0) + 0.08004 * std::exp((s_V + c_EKshift + 50.0) / 16.59));
        const double a_tiS_b = 23.62 + 1.0 / (0.001416 * std::exp(-(s_V + c_EKshift + 96.52) / 59.05) + 1.78e-08 * std::exp((s_V + c_EKshift + 114.1) / 8.079));
        const double a_AiS = 1.0 - a_AiF;
        const double a_tiF = a_tiF_b * a_delta_epi;
        const double a_tiS = a_tiS_b * a_delta_epi;
        const double r_Ito_a = (a_ass - s_Ito_a) / a_ta;
        const double r_Ito_ap = (a_assp - s_Ito_ap) / a_ta;
        const double a_i = a_AiF * s_Ito_iF + a_AiS * s_Ito_iS;
        const double a_ip = a_AiF * s_Ito_iFp + a_AiS * s_Ito_iSp;
        const double a_tiFp = a_dti_develop * a_dti_recover * a_tiF;
        const double a_tiSp = a_dti_develop * a_dti_recover * a_tiS;
        const double r_Ito_iF = (a_iss - s_Ito_iF) / a_tiF;
        const double r_Ito_iS = (a_iss - s_Ito_iS) / a_tiS;
        const double a_Ito_Ito = c_Gto * (s_V - a_EK) * ((1.0 - a_fItop) * s_Ito_a * a_i + a_fItop * s_Ito_ap * a_ip);
        const double r_Ito_iFp = (a_iss - s_Ito_iFp) / a_tiFp;
        const double r_Ito_iSp = (a_iss - s_Ito_iSp) / a_tiSp;

        /* ICaL */
        const double a_Afcaf = 0.3 + 0.6 / (1.0 + std::exp((s_V - 10.0) / 10.0));
        const double a_Ii = 0.5 * (s_Nai + s_Ki + s_Cli + 4.0 * s_Cai) / 1000.0;
        const double a_Iss = 0.5 * (s_Nass + s_Kss + s_Clss + 4.0 * s_Cass) / 1000.0;
        const double a_dss = ((s_V >= 31.4978) ? 1.0 : 1.0763 * std::exp(-1.007 * std::exp(-0.0829 * s_V)));
        const double a_fICaLp = 1.0 / (1.0 + c_KmCaMK / a_CaMKa);
        const double a_fss = 1.0 / (1.0 + std::exp((s_V + 19.58) / 3.696));
        const double a_jcass = 1.0 / (1.0 + std::exp((s_V + 18.08) / 2.7916));
        const double a_km2n = s_ICaL_jca * 1.0;
        const double a_tfcaf = 7.0 + 1.0 / (0.04 * std::exp(-(s_V - 4.0) / 7.0) + 0.04 * std::exp((s_V - 4.0) / 7.0));
        const double a_tfcas = 100.0 + 1.0 / (0.00012 * std::exp(-s_V / 3.0) + 0.00012 * std::exp(s_V / 7.0));
        const double a_tff = 7.0 + 1.0 / (0.0045 * std::exp(-(s_V + 20.0) / 10.0) + 0.0045 * std::exp((s_V + 20.0) / 10.0));
        const double a_tfs = 1000.0 + 1.0 / (3.5e-05 * std::exp(-(s_V + 5.0) / 4.0) + 3.5e-05 * std::exp((s_V + 5.0) / 6.0));
        const double a_Afcas = 1.0 - a_Afcaf;
        const double a_anca_i = 1.0 / (c_k2n / a_km2n + std::pow(1.0 + c_Kmn / s_Cai, 4.0));
        const double a_anca_ss = 1.0 / (c_k2n / a_km2n + std::pow(1.0 + c_Kmn / s_Cass, 4.0));
        const double a_fcass = a_fss;
        const double a_td = c_offset + 0.6 + 1.0 / (std::exp(-0.05 * (s_V + c_vShift + 6.0)) + std::exp(0.09 * (s_V + c_vShift + 14.0)));
        const double a_tfcafp = 2.5 * a_tfcaf;
        const double a_tffp = 2.5 * a_tff;
        const double r_ICaL_ff = (a_fss - s_ICaL_ff) / a_tff;
        const double r_ICaL_fs = (a_fss - s_ICaL_fs) / a_tfs;
        const double r_ICaL_jca = (a_jcass - s_ICaL_jca) / c_tjca;
        const double a_f = c_Aff * s_ICaL_ff + c_Afs * s_ICaL_fs;
        const double a_fca = a_Afcaf * s_ICaL_fcaf + a_Afcas * s_ICaL_fcas;
        const double a_fcap = a_Afcaf * s_ICaL_fcafp + a_Afcas * s_ICaL_fcas;
        const double a_fp = c_Aff * s_ICaL_ffp + c_Afs * s_ICaL_fs;
        const double a_gamma_cai = std::exp(-c_constA * 4.0 * (std::sqrt(a_Ii) / (1.0 + std::sqrt(a_Ii)) - 0.3 * a_Ii));
        const double a_gamma_cass = std::exp(-c_constA * 4.0 * (std::sqrt(a_Iss) / (1.0 + std::sqrt(a_Iss)) - 0.3 * a_Iss));
        const double a_gamma_ki = std::exp(-c_constA * 1.0 * (std::sqrt(a_Ii) / (1.0 + std::sqrt(a_Ii)) - 0.3 * a_Ii));
        const double a_gamma_kss = std::exp(-c_constA * 1.0 * (std::sqrt(a_Iss) / (1.0 + std::sqrt(a_Iss)) - 0.3 * a_Iss));
        const double a_gamma_nai = std::exp(-c_constA * 1.0 * (std::sqrt(a_Ii) / (1.0 + std::sqrt(a_Ii)) - 0.3 * a_Ii));
        const double a_gamma_nass = std::exp(-c_constA * 1.0 * (std::sqrt(a_Iss) / (1.0 + std::sqrt(a_Iss)) - 0.3 * a_Iss));
        const double r_ICaL_d = (a_dss - s_ICaL_d) / a_td;
        const double r_ICaL_fcaf = (a_fcass - s_ICaL_fcaf) / a_tfcaf;
        const double r_ICaL_fcafp = (a_fcass - s_ICaL_fcafp) / a_tfcafp;
        const double r_ICaL_fcas = (a_fcass - s_ICaL_fcas) / a_tfcas;
        const double r_ICaL_ffp = (a_fss - s_ICaL_ffp) / a_tffp;
        const double r_ICaL_nca_i = a_anca_i * c_k2n - s_ICaL_nca_i * a_km2n;
        const double r_ICaL_nca_ss = a_anca_ss * c_k2n - s_ICaL_nca_ss * a_km2n;

        /* IKr */
        const double a_IKr_IKr = c_GKr * std::sqrt(c_ko / 5.0) * s_IKr_O * (s_V - a_EK);

        /* INaCa */
        const double a_allo_i = 1.0 / (1.0 + std::pow(c_KmCaAct / s_Cai, 2.0));
        const double a_allo_ss = 1.0 / (1.0 + std::pow(c_KmCaAct / s_Cass, 2.0));
        const double a_h4_i = 1.0 + s_Nai / c_kna1 * (1.0 + s_Nai / c_kna2);
        const double a_h4_ss = 1.0 + s_Nass / c_kna1 * (1.0 + s_Nass / c_kna2);
        const double a_h5_i = s_Nai * s_Nai / (a_h4_i * c_kna1 * c_kna2);
        const double a_h5_ss = s_Nass * s_Nass / (a_h4_ss * c_kna1 * c_kna2);
        const double a_h6_i = 1.0 / a_h4_i;
        const double a_h6_ss = 1.0 / a_h4_ss;
        const double a_k6_i = a_h6_i * s_Cai * c_kcaon;
        const double a_k6_ss = a_h6_ss * s_Cass * c_kcaon;

        /* INaK */
        const double a_P = c_eP / (1.0 + c_H / c_Khp + s_Nai / c_Knap + s_Ki / c_Kxkur);
        const double a_b3 = c_k3m * a_P * c_H / (1.0 + c_MgATP / c_Kmgatp);

        /* intracellular_ions */
        const double r_Cansr = a_Jup - a_Jtr * c_vjsr / c_vnsr;
        const double r_Cli = (a_IClb + a_IClCa_sl) * c_Acap / (c_F * c_vmyo) + a_JdiffCl * c_vss / c_vmyo;
        const double r_Clss = a_IClCa_junc * c_Acap / (c_F * c_vss) - a_JdiffCl;
        const double a_Bcajsr = 1.0 / (1.0 + c_csqnmax * c_kmcsqn / std::pow(c_kmcsqn + s_Cajsr, 2.0));
        const double a_Bcass = 1.0 / (1.0 + c_BSRmax * c_KmBSR / std::pow(c_KmBSR + s_Cass, 2.0) + c_BSLmax * c_KmBSL / std::pow(c_KmBSL + s_Cass, 2.0));
        const double a_Bcai = 1.0 / (1.0 + c_cmdnmax * c_kmcmdn / std::pow(c_kmcmdn + s_Cai, 2.0) + c_trpnmax * c_kmtrpn / std::pow(c_kmtrpn + s_Cai, 2.0));

        /* membrane */
        const double a_vffrt = s_V * c_F * c_F / (c_R * c_T);
        const double a_vfrt = s_V * c_F / (c_R * c_T);

        /* ryr */
        const double a_fJrelp = 1.0 / (1.0 + c_KmCaMK / a_CaMKa);
        const double a_Jrel = c_Jrel_b * ((1.0 - a_fJrelp) * s_Jrel_np + a_fJrelp * s_Jrel_p);
        const double a_tau_rel_b = c_bt / (1.0 + 0.0123 / s_Cajsr);
        const double a_tau_rel = ((a_tau_rel_b < 0.001) ? 0.001 : a_tau_rel_b);
        const double a_tau_relp_b = c_btp / (1.0 + 0.0123 / s_Cajsr);
        const double a_tau_relp = ((a_tau_relp_b < 0.001) ? 0.001 : a_tau_relp_b);

        /* *remaining* */
        const double a_PhiCaK_i = 1.0 * a_vffrt * (a_gamma_ki * s_Ki * std::exp(1.0 * a_vfrt) - c_gamma_ko * c_ko) / (std::exp(1.0 * a_vfrt) - 1.0);
        const double a_PhiCaK_ss = 1.0 * a_vffrt * (a_gamma_kss * s_Kss * std::exp(1.0 * a_vfrt) - c_gamma_ko * c_ko) / (std::exp(1.0 * a_vfrt) - 1.0);
        const double a_PhiCaL_i = 4.0 * a_vffrt * (a_gamma_cai * s_Cai * std::exp(2.0 * a_vfrt) - c_gamma_cao * c_cao) / (std::exp(2.0 * a_vfrt) - 1.0);
        const double a_PhiCaL_ss = 4.0 * a_vffrt * (a_gamma_cass * s_Cass * std::exp(2.0 * a_vfrt) - c_gamma_cao * c_cao) / (std::exp(2.0 * a_vfrt) - 1.0);
        const double a_PhiCaNa_i = 1.0 * a_vffrt * (a_gamma_nai * s_Nai * std::exp(1.0 * a_vfrt) - c_gamma_nao * c_nao) / (std::exp(1.0 * a_vfrt) - 1.0);
        const double a_PhiCaNa_ss = 1.0 * a_vffrt * (a_gamma_nass * s_Nass * std::exp(1.0 * a_vfrt) - c_gamma_nao * c_nao) / (std::exp(1.0 * a_vfrt) - 1.0);
        const double a_ICab_ICab = c_PCab * 4.0 * a_vffrt * (a_gamma_cai * s_Cai * std::exp(2.0 * a_vfrt) - c_gamma_cao * c_cao) / (std::exp(2.0 * a_vfrt) - 1.0);
        const double a_alpha = 0.1161 * std::exp(0.299 * a_vfrt);
        const double a_alpha_2 = 0.0578 * std::exp(0.971 * a_vfrt);
        const double a_alpha_C2ToI = 5.2e-05 * std::exp(1.525 * a_vfrt);
        const double a_alpha_i = 0.2533 * std::exp(0.5953 * a_vfrt);
        const double a_beta = 0.2442 * std::exp(-1.604 * a_vfrt);
        const double a_beta_2 = 0.000349 * std::exp(-1.062 * a_vfrt);
        const double a_beta_i = 0.06525 * std::exp(-0.8209 * a_vfrt);
        const double a_hca = std::exp(c_qca * a_vfrt);
        const double a_hna = std::exp(c_qna * a_vfrt);
        const double a_Knai = c_Knai0 * std::exp(c_delta * a_vfrt / 3.0);
        const double a_Knao = c_Knao0 * std::exp((1.0 - c_delta) * a_vfrt / 3.0);
        const double a_INab_INab = c_PNab * a_vffrt * (s_Nai * std::exp(a_vfrt) - c_nao) / (std::exp(a_vfrt) - 1.0);
        const double r_Cajsr = a_Bcajsr * (a_Jtr - a_Jrel);
        const double a_ICaK_i = (1.0 - c_ICaL_fractionSS) * ((1.0 - a_fICaLp) * c_PCaK * a_PhiCaK_i * s_ICaL_d * (a_f * (1.0 - s_ICaL_nca_i) + s_ICaL_jca * a_fca * s_ICaL_nca_i) + a_fICaLp * c_PCaKp * a_PhiCaK_i * s_ICaL_d * (a_fp * (1.0 - s_ICaL_nca_i) + s_ICaL_jca * a_fcap * s_ICaL_nca_i));
        const double a_ICaK_ss = c_ICaL_fractionSS * ((1.0 - a_fICaLp) * c_PCaK * a_PhiCaK_ss * s_ICaL_d * (a_f * (1.0 - s_ICaL_nca_ss) + s_ICaL_jca * a_fca * s_ICaL_nca_ss) + a_fICaLp * c_PCaKp * a_PhiCaK_ss * s_ICaL_d * (a_fp * (1.0 - s_ICaL_nca_ss) + s_ICaL_jca * a_fcap * s_ICaL_nca_ss));
        const double a_ICaL_i = (1.0 - c_ICaL_fractionSS) * ((1.0 - a_fICaLp) * c_PCa * a_PhiCaL_i * s_ICaL_d * (a_f * (1.0 - s_ICaL_nca_i) + s_ICaL_jca * a_fca * s_ICaL_nca_i) + a_fICaLp * c_PCap * a_PhiCaL_i * s_ICaL_d * (a_fp * (1.0 - s_ICaL_nca_i) + s_ICaL_jca * a_fcap * s_ICaL_nca_i));
        const double a_ICaL_ss = c_ICaL_fractionSS * ((1.0 - a_fICaLp) * c_PCa * a_PhiCaL_ss * s_ICaL_d * (a_f * (1.0 - s_ICaL_nca_ss) + s_ICaL_jca * a_fca * s_ICaL_nca_ss) + a_fICaLp * c_PCap * a_PhiCaL_ss * s_ICaL_d * (a_fp * (1.0 - s_ICaL_nca_ss) + s_ICaL_jca * a_fcap * s_ICaL_nca_ss));
        const double a_ICaNa_i = (1.0 - c_ICaL_fractionSS) * ((1.0 - a_fICaLp) * c_PCaNa * a_PhiCaNa_i * s_ICaL_d * (a_f * (1.0 - s_ICaL_nca_i) + s_ICaL_jca * a_fca * s_ICaL_nca_i) + a_fICaLp * c_PCaNap * a_PhiCaNa_i * s_ICaL_d * (a_fp * (1.0 - s_ICaL_nca_i) + s_ICaL_jca * a_fcap * s_ICaL_nca_i));
        const double a_ICaNa_ss = c_ICaL_fractionSS * ((1.0 - a_fICaLp) * c_PCaNa * a_PhiCaNa_ss * s_ICaL_d * (a_f * (1.0 - s_ICaL_nca_ss) + s_ICaL_jca * a_fca * s_ICaL_nca_ss) + a_fICaLp * c_PCaNap * a_PhiCaNa_ss * s_ICaL_d * (a_fp * (1.0 - s_ICaL_nca_ss) + s_ICaL_jca * a_fcap * s_ICaL_nca_ss));
        const double a_beta_ItoC2 = a_beta_2 * a_beta_i * a_alpha_C2ToI / (a_alpha_2 * a_alpha_i);
        const double r_IKr_C2 = a_alpha * s_IKr_C3 + c_beta_1 * s_IKr_C1 - (a_beta + c_alpha_1) * s_IKr_C2;
        const double r_IKr_C3 = a_beta * s_IKr_C2 - a_alpha * s_IKr_C3;
        const double r_IKr_O = a_alpha_2 * s_IKr_C1 + a_beta_i * s_IKr_I - (a_beta_2 + a_alpha_i) * s_IKr_O;
        const double a_h1_i = 1.0 + s_Nai / c_kna3 * (1.0 + a_hna);
        const double a_h1_ss = 1.0 + s_Nass / c_kna3 * (1.0 + a_hna);
        const double a_h7_i = 1.0 + c_nao / c_kna3 * (1.0 + 1.0 / a_hna);
        const double a_h7_ss = 1.0 + c_nao / c_kna3 * (1.0 + 1.0 / a_hna);
        const double a_a1 = c_k1p * std::pow(s_Nai / a_Knai, 3.0) / (std::pow(1.0 + s_Nai / a_Knai, 3.0) + std::pow(1.0 + s_Ki / c_Kki, 2.0) - 1.0);
        const double a_a3 = c_k3p * std::pow(c_nao / a_Knao, 2.0) / (std::pow(1.0 + c_nao / a_Knao, 3.0) + std::pow(1.0 + c_ko / c_Kko, 2.0) - 1.0);
        const double a_b2 = c_k2m * std::pow(c_nao / a_Knao, 3.0) / (std::pow(1.0 + c_nao / a_Knao, 3.0) + std::pow(1.0 + c_ko / c_Kko, 2.0) - 1.0);
        const double a_b4 = c_k4m * std::pow(s_Ki / c_Kki, 2.0) / (std::pow(1.0 + s_Nai / a_Knai, 3.0) + std::pow(1.0 + s_Ki / c_Kki, 2.0) - 1.0);
        const double a_ICaK = a_ICaK_ss + a_ICaK_i;
        const double a_ICaL_ICaL = a_ICaL_ss + a_ICaL_i;
        const double a_ICaNa = a_ICaNa_ss + a_ICaNa_i;
        const double r_IKr_C1 = c_alpha_1 * s_IKr_C2 + a_beta_2 * s_IKr_O + a_beta_ItoC2 * s_IKr_I - (c_beta_1 + a_alpha_2 + a_alpha_C2ToI) * s_IKr_C1;
        const double r_IKr_I = a_alpha_C2ToI * s_IKr_C1 + a_alpha_i * s_IKr_O - (a_beta_ItoC2 + a_beta_i) * s_IKr_I;
        const double a_h2_i = s_Nai * a_hna / (c_kna3 * a_h1_i);
        const double a_h2_ss = s_Nass * a_hna / (c_kna3 * a_h1_ss);
        const double a_h3_i = 1.0 / a_h1_i;
        const double a_h3_ss = 1.0 / a_h1_ss;
        const double a_h8_i = c_nao / (c_kna3 * a_hna * a_h7_i);
        const double a_h8_ss = c_nao / (c_kna3 * a_hna * a_h7_ss);
        const double a_h9_i = 1.0 / a_h7_i;
        const double a_h9_ss = 1.0 / a_h7_ss;
        const double a_x1 = c_a4 * a_a1 * c_a2 + a_b2 * a_b4 * a_b3 + c_a2 * a_b4 * a_b3 + a_b3 * a_a1 * c_a2;
        const double a_x2 = a_b2 * c_b1 * a_b4 + a_a1 * c_a2 * a_a3 + a_a3 * c_b1 * a_b4 + c_a2 * a_a3 * a_b4;
        const double a_x3 = c_a2 * a_a3 * c_a4 + a_b3 * a_b2 * c_b1 + a_b2 * c_b1 * c_a4 + a_a3 * c_a4 * c_b1;
        const double a_x4 = a_b4 * a_b3 * a_b2 + a_a3 * c_a4 * a_a1 + a_b2 * c_a4 * a_a1 + a_b3 * a_b2 * a_a1;
        const double r_Kss = -a_ICaK_ss * c_Acap / (c_F * c_vss) - a_JdiffK;
        const double a_Jrel_inf_b = -c_a_rel * a_ICaL_ss / 1.0 / (1.0 + std::pow(c_cajsr_half / s_Cajsr, 8.0));
        const double a_Jrel_infp_b = -c_a_relp * a_ICaL_ss / 1.0 / (1.0 + std::pow(c_cajsr_half / s_Cajsr, 8.0));
        const double a_k3p_i = a_h9_i * c_wca;
        const double a_k3p_ss = a_h9_ss * c_wca;
        const double a_k3pp_i = a_h8_i * c_wnaca;
        const double a_k3pp_ss = a_h8_ss * c_wnaca;
        const double a_k4p_i = a_h3_i * c_wca / a_hca;
        const double a_k4p_ss = a_h3_ss * c_wca / a_hca;
        const double a_k4pp_i = a_h2_i * c_wnaca;
        const double a_k4pp_ss = a_h2_ss * c_wnaca;
        const double a_k7_i = a_h5_i * a_h2_i * c_wna;
        const double a_k7_ss = a_h5_ss * a_h2_ss * c_wna;
        const double a_k8_i = a_h8_i * c_h11_i * c_wna;
        const double a_k8_ss = a_h8_ss * c_h11_ss * c_wna;
        const double a_E1 = a_x1 / (a_x1 + a_x2 + a_x3 + a_x4);
        const double a_E2 = a_x2 / (a_x1 + a_x2 + a_x3 + a_x4);
        const double a_E3 = a_x3 / (a_x1 + a_x2 + a_x3 + a_x4);
        const double a_E4 = a_x4 / (a_x1 + a_x2 + a_x3 + a_x4);
        const double a_Jrel_inf  = a_Jrel_inf_b  * c_jrelTissueScale;
        const double a_Jrel_infp = a_Jrel_infp_b * c_jrelTissueScale;
        const double a_k3_i = a_k3p_i + a_k3pp_i;
        const double a_k3_ss = a_k3p_ss + a_k3pp_ss;
        const double a_k4_i = a_k4p_i + a_k4pp_i;
        const double a_k4_ss = a_k4p_ss + a_k4pp_ss;
        const double a_JnakK = 2.0 * (a_E4 * c_b1 - a_E3 * a_a1);
        const double a_JnakNa = 3.0 * (a_E1 * a_a3 - a_E2 * a_b3);
        const double r_Jrel_np = (a_Jrel_inf - s_Jrel_np) / a_tau_rel;
        const double r_Jrel_p = (a_Jrel_infp - s_Jrel_p) / a_tau_relp;
        const double a_x1_i = c_k2_i * a_k4_i * (a_k7_i + a_k6_i) + c_k5_i * a_k7_i * (c_k2_i + a_k3_i);
        const double a_x1_ss = c_k2_ss * a_k4_ss * (a_k7_ss + a_k6_ss) + c_k5_ss * a_k7_ss * (c_k2_ss + a_k3_ss);
        const double a_x2_i = c_k1_i * a_k7_i * (a_k4_i + c_k5_i) + a_k4_i * a_k6_i * (c_k1_i + a_k8_i);
        const double a_x2_ss = c_k1_ss * a_k7_ss * (a_k4_ss + c_k5_ss) + a_k4_ss * a_k6_ss * (c_k1_ss + a_k8_ss);
        const double a_x3_i = c_k1_i * a_k3_i * (a_k7_i + a_k6_i) + a_k8_i * a_k6_i * (c_k2_i + a_k3_i);
        const double a_x3_ss = c_k1_ss * a_k3_ss * (a_k7_ss + a_k6_ss) + a_k8_ss * a_k6_ss * (c_k2_ss + a_k3_ss);
        const double a_x4_i = c_k2_i * a_k8_i * (a_k4_i + c_k5_i) + a_k3_i * c_k5_i * (c_k1_i + a_k8_i);
        const double a_x4_ss = c_k2_ss * a_k8_ss * (a_k4_ss + c_k5_ss) + a_k3_ss * c_k5_ss * (c_k1_ss + a_k8_ss);
        const double a_INaK_INaK = c_Pnak * (c_zna * a_JnakNa + c_zk * a_JnakK);
        const double a_E1_i = a_x1_i / (a_x1_i + a_x2_i + a_x3_i + a_x4_i);
        const double a_E1_ss = a_x1_ss / (a_x1_ss + a_x2_ss + a_x3_ss + a_x4_ss);
        const double a_E2_i = a_x2_i / (a_x1_i + a_x2_i + a_x3_i + a_x4_i);
        const double a_E2_ss = a_x2_ss / (a_x1_ss + a_x2_ss + a_x3_ss + a_x4_ss);
        const double a_E3_i = a_x3_i / (a_x1_i + a_x2_i + a_x3_i + a_x4_i);
        const double a_E3_ss = a_x3_ss / (a_x1_ss + a_x2_ss + a_x3_ss + a_x4_ss);
        const double a_E4_i = a_x4_i / (a_x1_i + a_x2_i + a_x3_i + a_x4_i);
        const double a_E4_ss = a_x4_ss / (a_x1_ss + a_x2_ss + a_x3_ss + a_x4_ss);
        const double r_Ki = -(a_Ito_Ito + a_IKr_IKr + a_IKs_IKs + a_IK1_IK1 + a_IKb_IKb + a_I_katp_I_katp + Istim - 2.0 * a_INaK_INaK + a_ICaK_i) * c_Acap / (c_F * c_vmyo) + a_JdiffK * c_vss / c_vmyo;
        const double a_JncxCa_i = a_E2_i * c_k2_i - a_E1_i * c_k1_i;
        const double a_JncxCa_ss = a_E2_ss * c_k2_ss - a_E1_ss * c_k1_ss;
        const double a_JncxNa_i = 3.0 * (a_E4_i * a_k7_i - a_E1_i * a_k8_i) + a_E3_i * a_k4pp_i - a_E2_i * a_k3pp_i;
        const double a_JncxNa_ss = 3.0 * (a_E4_ss * a_k7_ss - a_E1_ss * a_k8_ss) + a_E3_ss * a_k4pp_ss - a_E2_ss * a_k3pp_ss;
        const double a_INaCa_i = (1.0 - c_INaCa_fractionSS) * c_Gncx * a_allo_i * (c_zna * a_JncxNa_i + c_zca * a_JncxCa_i);
        const double a_INaCa_ss = c_INaCa_fractionSS * c_Gncx * a_allo_ss * (c_zna * a_JncxNa_ss + c_zca * a_JncxCa_ss);
        const double r_Cai = a_Bcai * (-(a_ICaL_i + a_IpCa_IpCa + a_ICab_ICab - 2.0 * a_INaCa_i) * c_Acap / (2.0 * c_F * c_vmyo) - a_Jup * c_vnsr / c_vmyo + a_Jdiff * c_vss / c_vmyo);
        const double r_Cass = a_Bcass * (-(a_ICaL_ss - 2.0 * a_INaCa_ss) * c_Acap / (2.0 * c_F * c_vss) + a_Jrel * c_vjsr / c_vss - a_Jdiff);
        const double r_Nai = -(a_INa_INa + a_INaL_INaL + 3.0 * a_INaCa_i + a_ICaNa_i + 3.0 * a_INaK_INaK + a_INab_INab) * c_Acap / (c_F * c_vmyo) + a_JdiffNa * c_vss / c_vmyo;
        const double r_Nass = -(a_ICaNa_ss + 3.0 * a_INaCa_ss) * c_Acap / (c_F * c_vss) - a_JdiffNa;

        const double a_Iion_cm = (a_INa_INa + a_INaL_INaL + a_Ito_Ito + a_ICaL_ICaL + a_ICaNa + a_ICaK + a_IKr_IKr + a_IKs_IKs + a_IK1_IK1 + a_INaCa_i + a_INaCa_ss + a_INaK_INaK + a_INab_INab + a_IKb_IKb + a_IpCa_IpCa + a_ICab_ICab + a_IClCa + a_IClb + a_I_katp_I_katp);

        const double r_V = solveVmWithinODESolver ? (-a_Iion_cm - Istim) : 0.0;

        RATES_SoA[V_BASE + cellI] = r_V;
        RATES_SoA[CaMKt_BASE + cellI] = r_CaMKt;
        RATES_SoA[Nai_BASE + cellI] = r_Nai;
        RATES_SoA[Nass_BASE + cellI] = r_Nass;
        RATES_SoA[Ki_BASE + cellI] = r_Ki;
        RATES_SoA[Kss_BASE + cellI] = r_Kss;
        RATES_SoA[Cass_BASE + cellI] = r_Cass;
        RATES_SoA[Cansr_BASE + cellI] = r_Cansr;
        RATES_SoA[Cajsr_BASE + cellI] = r_Cajsr;
        RATES_SoA[Cai_BASE + cellI] = r_Cai;
        RATES_SoA[Cli_BASE + cellI] = r_Cli;
        RATES_SoA[Clss_BASE + cellI] = r_Clss;
        RATES_SoA[INa_m_BASE + cellI] = r_INa_m;
        RATES_SoA[INa_h_BASE + cellI] = r_INa_h;
        RATES_SoA[INa_j_BASE + cellI] = r_INa_j;
        RATES_SoA[INa_hp_BASE + cellI] = r_INa_hp;
        RATES_SoA[INa_jp_BASE + cellI] = r_INa_jp;
        RATES_SoA[INaL_mL_BASE + cellI] = r_INaL_mL;
        RATES_SoA[INaL_hL_BASE + cellI] = r_INaL_hL;
        RATES_SoA[INaL_hLp_BASE + cellI] = r_INaL_hLp;
        RATES_SoA[Ito_a_BASE + cellI] = r_Ito_a;
        RATES_SoA[Ito_iF_BASE + cellI] = r_Ito_iF;
        RATES_SoA[Ito_iS_BASE + cellI] = r_Ito_iS;
        RATES_SoA[Ito_ap_BASE + cellI] = r_Ito_ap;
        RATES_SoA[Ito_iFp_BASE + cellI] = r_Ito_iFp;
        RATES_SoA[Ito_iSp_BASE + cellI] = r_Ito_iSp;
        RATES_SoA[ICaL_d_BASE + cellI] = r_ICaL_d;
        RATES_SoA[ICaL_ff_BASE + cellI] = r_ICaL_ff;
        RATES_SoA[ICaL_fs_BASE + cellI] = r_ICaL_fs;
        RATES_SoA[ICaL_fcaf_BASE + cellI] = r_ICaL_fcaf;
        RATES_SoA[ICaL_fcas_BASE + cellI] = r_ICaL_fcas;
        RATES_SoA[ICaL_jca_BASE + cellI] = r_ICaL_jca;
        RATES_SoA[ICaL_ffp_BASE + cellI] = r_ICaL_ffp;
        RATES_SoA[ICaL_fcafp_BASE + cellI] = r_ICaL_fcafp;
        RATES_SoA[ICaL_nca_ss_BASE + cellI] = r_ICaL_nca_ss;
        RATES_SoA[ICaL_nca_i_BASE + cellI] = r_ICaL_nca_i;
        RATES_SoA[IKr_C1_BASE + cellI] = r_IKr_C1;
        RATES_SoA[IKr_C2_BASE + cellI] = r_IKr_C2;
        RATES_SoA[IKr_C3_BASE + cellI] = r_IKr_C3;
        RATES_SoA[IKr_I_BASE + cellI] = r_IKr_I;
        RATES_SoA[IKr_O_BASE + cellI] = r_IKr_O;
        RATES_SoA[IKs_xs1_BASE + cellI] = r_IKs_xs1;
        RATES_SoA[IKs_xs2_BASE + cellI] = r_IKs_xs2;
        RATES_SoA[Jrel_np_BASE + cellI] = r_Jrel_np;
        RATES_SoA[Jrel_p_BASE + cellI] = r_Jrel_p;

        if (SUPPORT_SoA != nullptr)
        {
            SUPPORT_SoA[TORORD_DYNCL_BATCH_SUPPORT_tau_INa_m * N + cellI] = a_tm;
            SUPPORT_SoA[TORORD_DYNCL_BATCH_SUPPORT_gInf_INa_m * N + cellI] = a_mss;
            SUPPORT_SoA[TORORD_DYNCL_BATCH_SUPPORT_tau_INa_h * N + cellI] = a_th;
            SUPPORT_SoA[TORORD_DYNCL_BATCH_SUPPORT_gInf_INa_h * N + cellI] = a_hss;
            SUPPORT_SoA[TORORD_DYNCL_BATCH_SUPPORT_tau_INa_j * N + cellI] = a_tj;
            SUPPORT_SoA[TORORD_DYNCL_BATCH_SUPPORT_gInf_INa_j * N + cellI] = a_jss;
            SUPPORT_SoA[TORORD_DYNCL_BATCH_SUPPORT_tau_INa_hp * N + cellI] = a_th;
            SUPPORT_SoA[TORORD_DYNCL_BATCH_SUPPORT_gInf_INa_hp * N + cellI] = a_hssp;
            SUPPORT_SoA[TORORD_DYNCL_BATCH_SUPPORT_tau_INa_jp * N + cellI] = a_tjp;
            SUPPORT_SoA[TORORD_DYNCL_BATCH_SUPPORT_gInf_INa_jp * N + cellI] = a_jss;
            SUPPORT_SoA[TORORD_DYNCL_BATCH_SUPPORT_tau_INaL_mL * N + cellI] = a_tmL;
            SUPPORT_SoA[TORORD_DYNCL_BATCH_SUPPORT_gInf_INaL_mL * N + cellI] = a_mLss;
            SUPPORT_SoA[TORORD_DYNCL_BATCH_SUPPORT_tau_INaL_hL * N + cellI] = c_thL;
            SUPPORT_SoA[TORORD_DYNCL_BATCH_SUPPORT_gInf_INaL_hL * N + cellI] = a_hLss;
            SUPPORT_SoA[TORORD_DYNCL_BATCH_SUPPORT_tau_INaL_hLp * N + cellI] = c_thLp;
            SUPPORT_SoA[TORORD_DYNCL_BATCH_SUPPORT_gInf_INaL_hLp * N + cellI] = a_hLssp;
            SUPPORT_SoA[TORORD_DYNCL_BATCH_SUPPORT_tau_Ito_a * N + cellI] = a_ta;
            SUPPORT_SoA[TORORD_DYNCL_BATCH_SUPPORT_gInf_Ito_a * N + cellI] = a_ass;
            SUPPORT_SoA[TORORD_DYNCL_BATCH_SUPPORT_tau_Ito_iF * N + cellI] = a_tiF;
            SUPPORT_SoA[TORORD_DYNCL_BATCH_SUPPORT_gInf_Ito_iF * N + cellI] = a_iss;
            SUPPORT_SoA[TORORD_DYNCL_BATCH_SUPPORT_tau_Ito_iS * N + cellI] = a_tiS;
            SUPPORT_SoA[TORORD_DYNCL_BATCH_SUPPORT_gInf_Ito_iS * N + cellI] = a_iss;
            SUPPORT_SoA[TORORD_DYNCL_BATCH_SUPPORT_tau_Ito_ap * N + cellI] = a_ta;
            SUPPORT_SoA[TORORD_DYNCL_BATCH_SUPPORT_gInf_Ito_ap * N + cellI] = a_assp;
            SUPPORT_SoA[TORORD_DYNCL_BATCH_SUPPORT_tau_Ito_iFp * N + cellI] = a_tiFp;
            SUPPORT_SoA[TORORD_DYNCL_BATCH_SUPPORT_gInf_Ito_iFp * N + cellI] = a_iss;
            SUPPORT_SoA[TORORD_DYNCL_BATCH_SUPPORT_tau_Ito_iSp * N + cellI] = a_tiSp;
            SUPPORT_SoA[TORORD_DYNCL_BATCH_SUPPORT_gInf_Ito_iSp * N + cellI] = a_iss;
            SUPPORT_SoA[TORORD_DYNCL_BATCH_SUPPORT_tau_ICaL_d * N + cellI] = a_td;
            SUPPORT_SoA[TORORD_DYNCL_BATCH_SUPPORT_gInf_ICaL_d * N + cellI] = a_dss;
            SUPPORT_SoA[TORORD_DYNCL_BATCH_SUPPORT_tau_ICaL_ff * N + cellI] = a_tff;
            SUPPORT_SoA[TORORD_DYNCL_BATCH_SUPPORT_gInf_ICaL_ff * N + cellI] = a_fss;
            SUPPORT_SoA[TORORD_DYNCL_BATCH_SUPPORT_tau_ICaL_fs * N + cellI] = a_tfs;
            SUPPORT_SoA[TORORD_DYNCL_BATCH_SUPPORT_gInf_ICaL_fs * N + cellI] = a_fss;
            SUPPORT_SoA[TORORD_DYNCL_BATCH_SUPPORT_tau_ICaL_fcaf * N + cellI] = a_tfcaf;
            SUPPORT_SoA[TORORD_DYNCL_BATCH_SUPPORT_gInf_ICaL_fcaf * N + cellI] = a_fcass;
            SUPPORT_SoA[TORORD_DYNCL_BATCH_SUPPORT_tau_ICaL_fcas * N + cellI] = a_tfcas;
            SUPPORT_SoA[TORORD_DYNCL_BATCH_SUPPORT_gInf_ICaL_fcas * N + cellI] = a_fcass;
            SUPPORT_SoA[TORORD_DYNCL_BATCH_SUPPORT_tau_ICaL_jca * N + cellI] = c_tjca;
            SUPPORT_SoA[TORORD_DYNCL_BATCH_SUPPORT_gInf_ICaL_jca * N + cellI] = a_jcass;
            SUPPORT_SoA[TORORD_DYNCL_BATCH_SUPPORT_tau_ICaL_ffp * N + cellI] = a_tffp;
            SUPPORT_SoA[TORORD_DYNCL_BATCH_SUPPORT_gInf_ICaL_ffp * N + cellI] = a_fss;
            SUPPORT_SoA[TORORD_DYNCL_BATCH_SUPPORT_tau_ICaL_fcafp * N + cellI] = a_tfcafp;
            SUPPORT_SoA[TORORD_DYNCL_BATCH_SUPPORT_gInf_ICaL_fcafp * N + cellI] = a_fcass;
            SUPPORT_SoA[TORORD_DYNCL_BATCH_SUPPORT_tau_IKs_xs1 * N + cellI] = a_txs1;
            SUPPORT_SoA[TORORD_DYNCL_BATCH_SUPPORT_gInf_IKs_xs1 * N + cellI] = a_xs1ss;
            SUPPORT_SoA[TORORD_DYNCL_BATCH_SUPPORT_tau_IKs_xs2 * N + cellI] = a_txs2;
            SUPPORT_SoA[TORORD_DYNCL_BATCH_SUPPORT_gInf_IKs_xs2 * N + cellI] = a_xs2ss;
            SUPPORT_SoA[TORORD_DYNCL_BATCH_SUPPORT_tau_Jrel_np * N + cellI] = a_tau_rel;
            SUPPORT_SoA[TORORD_DYNCL_BATCH_SUPPORT_gInf_Jrel_np * N + cellI] = a_Jrel_inf;
            SUPPORT_SoA[TORORD_DYNCL_BATCH_SUPPORT_tau_Jrel_p * N + cellI] = a_tau_relp;
            SUPPORT_SoA[TORORD_DYNCL_BATCH_SUPPORT_gInf_Jrel_p * N + cellI] = a_Jrel_infp;
            SUPPORT_SoA[TORORD_DYNCL_BATCH_SUPPORT_Iion_cm * N + cellI] = a_Iion_cm;
        }
    }
}

} // End namespace Foam

// ************************************************************************* //
