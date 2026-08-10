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
    State-major (SoA) batched evaluator for the Gaur 2021 ionic model.
\*---------------------------------------------------------------------------*/

#pragma once

#include "Gaur_2021.H"
#include "Gaur_2021Names.H"
#include "StimulusProtocolPOD.H"
#include <cmath>

namespace Foam
{

enum GAUR_BATCH_SUPPORT_INDEX
{
    GAUR_BATCH_SUPPORT_tau_m,
    GAUR_BATCH_SUPPORT_gInf_m,
    GAUR_BATCH_SUPPORT_tau_h,
    GAUR_BATCH_SUPPORT_gInf_h,
    GAUR_BATCH_SUPPORT_tau_j,
    GAUR_BATCH_SUPPORT_gInf_j,
    GAUR_BATCH_SUPPORT_tau_ml,
    GAUR_BATCH_SUPPORT_gInf_ml,
    GAUR_BATCH_SUPPORT_tau_hl,
    GAUR_BATCH_SUPPORT_gInf_hl,
    GAUR_BATCH_SUPPORT_tau_d,
    GAUR_BATCH_SUPPORT_gInf_d,
    GAUR_BATCH_SUPPORT_tau_fca,
    GAUR_BATCH_SUPPORT_gInf_fca,
    GAUR_BATCH_SUPPORT_tau_ff,
    GAUR_BATCH_SUPPORT_gInf_ff,
    GAUR_BATCH_SUPPORT_tau_fs,
    GAUR_BATCH_SUPPORT_gInf_fs,
    GAUR_BATCH_SUPPORT_tau_xr,
    GAUR_BATCH_SUPPORT_gInf_xr,
    GAUR_BATCH_SUPPORT_tau_xs1,
    GAUR_BATCH_SUPPORT_gInf_xs1,
    GAUR_BATCH_SUPPORT_tau_xs2,
    GAUR_BATCH_SUPPORT_gInf_xs2,
    GAUR_BATCH_SUPPORT_tau_aa,
    GAUR_BATCH_SUPPORT_gInf_aa,
    GAUR_BATCH_SUPPORT_tau_Jrel1,
    GAUR_BATCH_SUPPORT_gInf_Jrel1,
    GAUR_BATCH_SUPPORT_tau_Jrel2,
    GAUR_BATCH_SUPPORT_gInf_Jrel2,
    GAUR_BATCH_SUPPORT_Iion_cm,

    NUM_GAUR_BATCH_SUPPORT
};

static_assert(NUM_GAUR_BATCH_SUPPORT == 31, "Gaur compact support layout drift");

CARDIAC_HOST_DEVICE inline void GaurComputeVariablesBatch
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

    const double c_F = CONSTANTS[AC_F];
    const double c_R = CONSTANTS[AC_R];
    const double c_T = CONSTANTS[AC_T];
    const double c_cao = CONSTANTS[AC_cao];
    const double c_ko = CONSTANTS[AC_ko];
    const double c_nao = CONSTANTS[AC_nao];
    const double c_Acap = CONSTANTS[AC_Acap];
    const double c_vjsr = CONSTANTS[AC_vjsr];
    const double c_vmyo = CONSTANTS[AC_vmyo];
    const double c_vnsr = CONSTANTS[AC_vnsr];
    const double c_vss = CONSTANTS[AC_vss];
    const double c_vcsr = CONSTANTS[AC_vcsr];
    const double c_vmyo1 = CONSTANTS[AC_vmyo1];
    const double c_vnsr1 = CONSTANTS[AC_vnsr1];
    const double c_vnsr2 = CONSTANTS[AC_vnsr2];
    const double c_vmyo2 = CONSTANTS[AC_vmyo2];
    const double c_CaMKo = CONSTANTS[AC_CaMKo];
    const double c_ECl = CONSTANTS[AC_ECl];
    const double c_KmCaM = CONSTANTS[AC_KmCaM];
    const double c_KmCaMK = CONSTANTS[AC_KmCaMK];
    const double c_PKNa = CONSTANTS[AC_PKNa];
    const double c_aCaMK = CONSTANTS[AC_aCaMK];
    const double c_bCaMK = CONSTANTS[AC_bCaMK];
    const double c_GNa = CONSTANTS[AC_GNa];
    const double c_GNaL = CONSTANTS[AC_GNaL];
    const double c_tau_hl = CONSTANTS[AC_tau_hl];
    const double c_PCa = CONSTANTS[AC_PCa];
    const double c_vhalf_d = CONSTANTS[AC_vhalf_d];
    const double c_vhalff = CONSTANTS[AC_vhalff];
    const double c_zca = CONSTANTS[AC_zca];
    const double c_PCaK = CONSTANTS[AC_PCaK];
    const double c_PCaNa = CONSTANTS[AC_PCaNa];
    const double c_GKr = CONSTANTS[AC_GKr];
    const double c_GKs = CONSTANTS[AC_GKs];
    const double c_GK1 = CONSTANTS[AC_GK1];
    const double c_Gto = CONSTANTS[AC_Gto];
    const double c_SOICR = CONSTANTS[AC_SOICR];
    const double c_grelbarjsrol = CONSTANTS[AC_grelbarjsrol];
    const double c_tau_gap = CONSTANTS[AC_tau_gap];
    const double c_tauoff = CONSTANTS[AC_tauoff];
    const double c_tauon = CONSTANTS[AC_tauon];
    const double c_Gncx = CONSTANTS[AC_Gncx];
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
    const double c_zna = CONSTANTS[AC_zna];
    const double c_k2 = CONSTANTS[AC_k2];
    const double c_k5 = CONSTANTS[AC_k5];
    const double c_h11 = CONSTANTS[AC_h11];
    const double c_k1 = CONSTANTS[AC_k1];
    const double c_k21 = CONSTANTS[AC_k21];
    const double c_k51 = CONSTANTS[AC_k51];
    const double c_h1111 = CONSTANTS[AC_h1111];
    const double c_k11 = CONSTANTS[AC_k11];
    const double c_H_const = CONSTANTS[AC_H];
    const double c_Khp = CONSTANTS[AC_Khp];
    const double c_Kki = CONSTANTS[AC_Kki];
    const double c_Kko = CONSTANTS[AC_Kko];
    const double c_Kmgatp = CONSTANTS[AC_Kmgatp];
    const double c_Knai0 = CONSTANTS[AC_Knai0];
    const double c_Knao0 = CONSTANTS[AC_Knao0];
    const double c_Knap = CONSTANTS[AC_Knap];
    const double c_Kxkur = CONSTANTS[AC_Kxkur];
    const double c_MgATP = CONSTANTS[AC_MgATP];
    const double c_Pnak = CONSTANTS[AC_Pnak];
    const double c_delta = CONSTANTS[AC_delta];
    const double c_eP = CONSTANTS[AC_eP];
    const double c_k1p = CONSTANTS[AC_k1p];
    const double c_k2m = CONSTANTS[AC_k2m];
    const double c_k3m = CONSTANTS[AC_k3m];
    const double c_k3p2 = CONSTANTS[AC_k3p2];
    const double c_k4m = CONSTANTS[AC_k4m];
    const double c_zk = CONSTANTS[AC_zk];
    const double c_a2 = CONSTANTS[AC_a2];
    const double c_a4 = CONSTANTS[AC_a4];
    const double c_b1 = CONSTANTS[AC_b1];
    const double c_GKb = CONSTANTS[AC_GKb];
    const double c_PNab = CONSTANTS[AC_PNab];
    const double c_PCab = CONSTANTS[AC_PCab];
    const double c_GpCa = CONSTANTS[AC_GpCa];
    const double c_BSLmax = CONSTANTS[AC_BSLmax];
    const double c_BSRmax = CONSTANTS[AC_BSRmax];
    const double c_KmBSL = CONSTANTS[AC_KmBSL];
    const double c_KmBSR = CONSTANTS[AC_KmBSR];
    const double c_cmdnmax = CONSTANTS[AC_cmdnmax];
    const double c_csqnmax = CONSTANTS[AC_csqnmax];
    const double c_kmcmdn = CONSTANTS[AC_kmcmdn];
    const double c_kmcsqn = CONSTANTS[AC_kmcsqn];
    const double c_kmtrpn = CONSTANTS[AC_kmtrpn];
    const double c_trpnmax = CONSTANTS[AC_trpnmax];

    const int cell_v_BASE = cell_v * N;
    const int nai_BASE = nai * N;
    const int nass_BASE = nass * N;
    const int ki_BASE = ki * N;
    const int kss_BASE = kss * N;
    const int cai_BASE = cai * N;
    const int cai2_BASE = cai2 * N;
    const int cass_BASE = cass * N;
    const int cansr_BASE = cansr * N;
    const int cajsr_BASE = cajsr * N;
    const int cacsr_BASE = cacsr * N;
    const int I_Na_m_BASE = I_Na_m * N;
    const int I_Na_h_BASE = I_Na_h * N;
    const int I_Na_j_BASE = I_Na_j * N;
    const int INaL_ml_BASE = INaL_ml * N;
    const int INaL_hl_BASE = INaL_hl * N;
    const int ICaL_d_BASE = ICaL_d * N;
    const int ICaL_fca_BASE = ICaL_fca * N;
    const int IKr_xr_BASE = IKr_xr * N;
    const int IKs_xs1_BASE = IKs_xs1 * N;
    const int IKs_xs2_BASE = IKs_xs2 * N;
    const int ITo_aa_BASE = ITo_aa * N;
    const int CICR_Jrel2_BASE = CICR_Jrel2 * N;
    const int CICR_Jrel1_BASE = CICR_Jrel1 * N;
    const int CaMK_CaMKt_BASE = CaMK_CaMKt * N;
    const int CICR_tjsrol_BASE = CICR_tjsrol * N;
    const int CICR_A_BASE = CICR_A * N;
    const int ICaL_fs_BASE = ICaL_fs * N;
    const int ICaL_ff_BASE = ICaL_ff * N;

#ifdef _OPENMP
    #pragma omp simd
#endif
    for (int cellI = beginCell; cellI < endCell; ++cellI)
    {
        const double s_cell_v = STATES_SoA[cell_v_BASE + cellI];
        const double s_nai = STATES_SoA[nai_BASE + cellI];
        const double s_nass = STATES_SoA[nass_BASE + cellI];
        const double s_ki = STATES_SoA[ki_BASE + cellI];
        const double s_kss = STATES_SoA[kss_BASE + cellI];
        const double s_cai = STATES_SoA[cai_BASE + cellI];
        const double s_cai2 = STATES_SoA[cai2_BASE + cellI];
        const double s_cass = STATES_SoA[cass_BASE + cellI];
        const double s_cansr = STATES_SoA[cansr_BASE + cellI];
        const double s_cajsr = STATES_SoA[cajsr_BASE + cellI];
        const double s_cacsr = STATES_SoA[cacsr_BASE + cellI];
        const double s_I_Na_m = STATES_SoA[I_Na_m_BASE + cellI];
        const double s_I_Na_h = STATES_SoA[I_Na_h_BASE + cellI];
        const double s_I_Na_j = STATES_SoA[I_Na_j_BASE + cellI];
        const double s_INaL_ml = STATES_SoA[INaL_ml_BASE + cellI];
        const double s_INaL_hl = STATES_SoA[INaL_hl_BASE + cellI];
        const double s_ICaL_d = STATES_SoA[ICaL_d_BASE + cellI];
        const double s_ICaL_fca = STATES_SoA[ICaL_fca_BASE + cellI];
        const double s_IKr_xr = STATES_SoA[IKr_xr_BASE + cellI];
        const double s_IKs_xs1 = STATES_SoA[IKs_xs1_BASE + cellI];
        const double s_IKs_xs2 = STATES_SoA[IKs_xs2_BASE + cellI];
        const double s_ITo_aa = STATES_SoA[ITo_aa_BASE + cellI];
        const double s_CICR_Jrel2 = STATES_SoA[CICR_Jrel2_BASE + cellI];
        const double s_CICR_Jrel1 = STATES_SoA[CICR_Jrel1_BASE + cellI];
        const double s_CaMK_CaMKt = STATES_SoA[CaMK_CaMKt_BASE + cellI];
        const double s_CICR_tjsrol = STATES_SoA[CICR_tjsrol_BASE + cellI];
        const double s_CICR_A = STATES_SoA[CICR_A_BASE + cellI];
        const double s_ICaL_fs = STATES_SoA[ICaL_fs_BASE + cellI];
        const double s_ICaL_ff = STATES_SoA[ICaL_ff_BASE + cellI];

        const double Jdiff = (s_cass - s_cai) / 0.2;
        const double JdiffK = (s_kss - s_ki) / 2.0;
        const double JdiffNa = (s_nass - s_nai) / 2.0;

        const double EK = c_R * c_T / c_F * std::log(c_ko / s_ki);
        const double ENa = c_R * c_T / c_F * std::log(c_nao / s_nai);
        const double vfrt = s_cell_v * c_F / (c_R * c_T);
        const double CaMKb = c_CaMKo * (1.0 - s_CaMK_CaMKt) / (1.0 + c_KmCaM / s_cass);
        const double EKs = c_R * c_T / c_F * std::log((c_ko + c_PKNa * c_nao) / (s_ki + c_PKNa * s_nai));
        const double CaMKa = CaMKb + s_CaMK_CaMKt;
        const double diff_CaMKt = c_aCaMK * CaMKb * (CaMKb + s_CaMK_CaMKt) - c_bCaMK * s_CaMK_CaMKt;
        const double CaMK_f = 1.0 / (1.0 + c_KmCaMK / CaMKa);

        const double aa_h = ((s_cell_v >= -40.0) ? 0.0 : 0.057 * std::exp(-(s_cell_v + 80.0) / 6.8));
        const double aa_j = ((s_cell_v >= -40.0) ? 0.0 : (-25428.0 * std::exp(0.2444 * s_cell_v) - 6.948e-06 * std::exp(-0.04391 * s_cell_v)) * (s_cell_v + 37.78) / (1.0 + std::exp(0.311 * (s_cell_v + 79.23))));
        const double aa_m = 1.0 / (1.0 + std::exp((-60.0 - s_cell_v) / 5.0));
        const double bb_h = ((s_cell_v >= -40.0) ? 0.77 / (0.13 * (1.0 + std::exp(-(s_cell_v + 10.66) / 11.1))) : 2.7 * std::exp(0.079 * s_cell_v) + 310000.0 * std::exp(0.3485 * s_cell_v));
        const double bb_j = ((s_cell_v >= -40.0) ? 0.6 * std::exp(0.057 * s_cell_v) / (1.0 + std::exp(-0.1 * (s_cell_v + 32.0))) : 0.02424 * std::exp(-0.01052 * s_cell_v) / (1.0 + std::exp(-0.1378 * (s_cell_v + 40.14))));
        const double bb_m = 0.1 / (1.0 + std::exp((s_cell_v + 35.0) / 5.0)) + 0.1 / (1.0 + std::exp((s_cell_v - 50.0) / 200.0));
        const double h_inf = 1.0 / ((1.0 + std::exp((s_cell_v + 71.55) / 7.43)) * (1.0 + std::exp((s_cell_v + 71.55) / 7.43)));
        const double m_inf = 1.0 / ((1.0 + std::exp((-56.86 - s_cell_v) / 9.03)) * (1.0 + std::exp((-56.86 - s_cell_v) / 9.03)));
        const double INa = c_GNa * s_I_Na_m * s_I_Na_m * s_I_Na_m * s_I_Na_h * s_I_Na_j * (s_cell_v - ENa);
        const double j_inf = h_inf;
        const double tau_h = 1.0 / (aa_h + bb_h);
        const double tau_j = 1.0 / (aa_j + bb_j);
        const double tau_m = aa_m * bb_m;

        const double aml = 0.32 * ((s_cell_v == -47.13) ? 10.0 : -(s_cell_v + 47.13) / (std::exp(-0.1 * (s_cell_v + 47.13)) - 1.0));
        const double bml = 0.08 * std::exp(-s_cell_v / 11.0);
        const double hl_inf = 1.0 / (1.0 + std::exp((s_cell_v + 91.0) / 6.1));
        const double INaL_val = c_GNaL * s_INaL_ml * s_INaL_ml * s_INaL_ml * s_INaL_hl * (s_cell_v - ENa);
        const double ml_inf = aml / (aml + bml);
        const double tau_ml = 1.0 / (aml + bml);

        const double PhiCaK = c_F * (0.75 * s_kss * std::exp(vfrt) - 0.75 * c_ko) * ((s_cell_v == 0.0) ? 1.0 : vfrt / (std::exp(vfrt) - 1.0));
        const double PhiCaL = 2.0 * c_F * (s_cass * std::exp(2.0 * (s_cell_v * c_F / (c_R * c_T))) - 0.341 * c_cao) * ((s_cell_v == 0.0) ? 1.0 : 2.0 * vfrt / (std::exp(2.0 * vfrt) - 1.0));
        const double PhiCaNa = c_F * (0.75 * s_nass * std::exp(vfrt) - 0.75 * c_nao) * ((s_cell_v == 0.0) ? 1.0 : vfrt / (std::exp(vfrt) - 1.0));
        const double f_act = s_ICaL_ff * s_ICaL_fs;
        const double tau_d = 0.6 + 1.0 / (std::exp(-0.05 * (s_cell_v + 6.0)) + std::exp(0.09 * (s_cell_v + 14.0)));
        const double tau_fca = 10.0 * CaMK_f + 0.5 + 1.0 / (1.0 + s_cass / 0.003);
        const double tau_ff = 7.0 + 1.0 / (0.0045 * std::exp(-(s_cell_v + 20.0) / 10.0) + 0.0045 * std::exp((s_cell_v + 20.0) / 10.0));
        const double tau_fs = 70.0 + 1.0 / (3.5e-05 * std::exp(-(s_cell_v + 5.0) / 4.0) + 3.5e-05 * std::exp((s_cell_v + 5.0) / 6.0));
        const double ICaL_val = c_PCa * PhiCaL * s_ICaL_d * f_act * s_ICaL_fca;
        const double d_inf = 1.0 / (1.0 + std::exp(-(s_cell_v - c_vhalf_d) / 6.2));
        const double f_inf = 1.0 / (1.0 + std::exp((s_cell_v - c_vhalff) / 4.9)) + 0.35 / (1.0 + std::exp((45.0 - s_cell_v) / 20.0));
        const double ICaK = c_PCaK * PhiCaK * s_ICaL_d * f_act * s_ICaL_fca;
        const double ICaNa = c_PCaNa * PhiCaNa * s_ICaL_d * f_act * s_ICaL_fca;
        const double fca_inf = 0.3 / (1.0 - ((ICaL_val > 0.0) ? 0.0 : ICaL_val / 0.05)) + 0.55 / (1.0 + s_cass / 0.003) + 0.15;
        const double ff_inf = f_inf;
        const double fs_inf = f_inf;

        const double rkr = 1.0 / (1.0 + std::exp((s_cell_v + 22.0) / 15.0));
        const double tau_xr = 12.98 + 1.0 / (0.3652 * std::exp((s_cell_v - 31.66) / 3.869) + 4.123e-05 * std::exp(-(s_cell_v - 47.78) / 20.38));
        const double xr_inf = 1.0 / (1.0 + std::exp(-(s_cell_v + 56.8) / 17.8));
        const double IKr_val = c_GKr * std::sqrt(c_ko / 5.4) * s_IKr_xr * rkr * (s_cell_v - EK);

        const double KsCa = 1.0 + 0.6 / (1.0 + std::pow(3.8e-05 / s_cai, 1.4));
        const double tau_xs1 = 300.0 + 1.0 / (1e-06 * std::exp((s_cell_v + 50.0) / 20.0) + 0.04 * std::exp(-(s_cell_v + 50.0) / 20.0));
        const double tau_xs2 = 1.0 / (0.01 * std::exp((s_cell_v - 50.0) / 20.0) + 0.0193 * std::exp(-(s_cell_v + 66.54) / 31.0));
        const double xs1_inf = 1.0 / (1.0 + std::exp(-(s_cell_v - 25.1) / 37.1));
        const double IKs_val = c_GKs * KsCa * s_IKs_xs1 * s_IKs_xs2 * (s_cell_v - EKs);
        const double xs2_inf = xs1_inf;

        const double rk1 = 1.0 / (1.0 + std::exp((s_cell_v + 79.3 - 2.6 * c_ko) / 19.6));
        const double IK1_val = c_GK1 * std::sqrt(c_ko / 5.4) * rk1 * (s_cell_v - EK);

        const double alpha_aa = 0.025 / (1.0 + std::exp((s_cell_v + 58.0) / 5.0));
        const double beta_aa = 1.0 / (5.0 * (1.0 + std::exp((s_cell_v + 19.0) / -9.0)));
        const double rito2 = 1.0 / (1.0 + std::exp(-(s_cell_v + 10.0) / 5.0));
        const double diff_A = ((s_CICR_tjsrol < 5.0) ? -s_CICR_A / 0.1 : 1.0);
        const double trel1factor = 1.0 * 20.0 * (1.0 + 1.0 / (1.0 + std::pow(c_KmCaMK / CaMKa, 8.0))) / (1.0 + std::pow(0.5 / s_cajsr, 8.0));
        const double trel2factor = 50.0 * (1.0 + 1.0 / (1.0 + std::pow(c_KmCaMK / CaMKa, 8.0))) / (1.0 + std::pow(0.5 / s_cacsr, 8.0));
        const double Jgap = (s_cai - s_cai2) / c_tau_gap;
        const double diff_tjsrol = (((s_cajsr > c_SOICR) && (s_CICR_A > 45.0)) ? -s_CICR_tjsrol / 0.001 : 1.0);
        const double greljsrol = c_grelbarjsrol * (1.0 - std::exp(-s_CICR_tjsrol / c_tauon)) * std::exp(-s_CICR_tjsrol / c_tauoff);
        const double tau_Jrel1 = ((0.001 > trel1factor) ? 0.001 : trel1factor);
        const double tau_Jrel2 = ((0.001 > trel2factor) ? 0.001 : trel2factor);
        const double Jrelol = greljsrol * (s_cajsr - s_cass);
        const double Jrel = s_CICR_Jrel1 + Jrelol;

        const double allo = 1.0 / (1.0 + c_KmCaAct / s_cai * c_KmCaAct / s_cai);
        const double h4 = 1.0 + s_nai / c_kna1 * (1.0 + s_nai / c_kna2);
        const double hca = std::exp(c_qca * s_cell_v * c_F / (c_R * c_T));
        const double hna = std::exp(c_qna * s_cell_v * c_F / (c_R * c_T));
        const double h1 = 1.0 + s_nai / c_kna3 * (1.0 + hna);
        const double h5 = s_nai * s_nai / (h4 * c_kna1 * c_kna2);
        const double h6 = 1.0 / h4;
        const double h7 = 1.0 + c_nao / c_kna3 * (1.0 + 1.0 / hna);
        const double h2 = s_nai * hna / (c_kna3 * h1);
        const double h3 = 1.0 / h1;
        const double h8 = c_nao / (c_kna3 * hna * h7);
        const double h9 = 1.0 / h7;
        const double k6 = h6 * s_cai * c_kcaon;
        const double k3p = h9 * c_wca;
        const double k3pp = h8 * c_wnaca;
        const double k4p = h3 * c_wca / hca;
        const double k4pp = h2 * c_wnaca;
        const double k7 = h5 * h2 * c_wna;
        const double k8 = h8 * c_h11 * c_wna;
        const double k3 = k3p + k3pp;
        const double k4 = k4p + k4pp;
        const double x1 = c_k2 * k4 * (k7 + k6) + c_k5 * k7 * (c_k2 + k3);
        const double x2 = c_k1 * k7 * (k4 + c_k5) + k4 * k6 * (c_k1 + k8);
        const double x3 = c_k1 * k3 * (k7 + k6) + k8 * k6 * (c_k2 + k3);
        const double x4 = c_k2 * k8 * (k4 + c_k5) + k3 * c_k5 * (c_k1 + k8);
        const double E1 = x1 / (x1 + x2 + x3 + x4);
        const double E2 = x2 / (x1 + x2 + x3 + x4);
        const double E3 = x3 / (x1 + x2 + x3 + x4);
        const double E4 = x4 / (x1 + x2 + x3 + x4);
        const double JncxCa = E2 * c_k2 - E1 * c_k1;
        const double JncxNa = 3.0 * (E4 * k7 - E1 * k8) + E3 * k4pp - E2 * k3pp;
        const double INaCa_i_val = 0.8 * c_Gncx * allo * (c_zna * JncxNa + c_zca * JncxCa);

        const double allo1 = 1.0 / (1.0 + c_KmCaAct / s_cass * c_KmCaAct / s_cass);
        const double h111 = 1.0 + s_nass / c_kna3 * (1.0 + hna);
        const double h41 = 1.0 + s_nass / c_kna1 * (1.0 + s_nass / c_kna2);
        const double h71 = 1.0 + c_nao / c_kna3 * (1.0 + 1.0 / hna);
        const double h21 = s_nass * hna / (c_kna3 * h111);
        const double h31 = 1.0 / h111;
        const double h51 = s_nass * s_nass / (h41 * c_kna1 * c_kna2);
        const double h61 = 1.0 / h41;
        const double h81 = c_nao / (c_kna3 * hna * h71);
        const double h91 = 1.0 / h71;
        const double k3p1 = h91 * c_wca;
        const double k3pp1 = h81 * c_wnaca;
        const double k4p1 = h31 * c_wca / hca;
        const double k4pp1 = h21 * c_wnaca;
        const double k61 = h61 * s_cass * c_kcaon;
        const double k71 = h51 * h21 * c_wna;
        const double k81 = h81 * c_h1111 * c_wna;
        const double k31 = k3p1 + k3pp1;
        const double k41 = k4p1 + k4pp1;
        const double x11 = c_k21 * k41 * (k71 + k61) + c_k51 * k71 * (c_k21 + k31);
        const double x21 = c_k11 * k71 * (k41 + c_k51) + k41 * k61 * (c_k11 + k81);
        const double x31 = c_k11 * k31 * (k71 + k61) + k81 * k61 * (c_k21 + k31);
        const double x41 = c_k21 * k81 * (k41 + c_k51) + k31 * c_k51 * (c_k11 + k81);
        const double E11 = x11 / (x11 + x21 + x31 + x41);
        const double E21 = x21 / (x11 + x21 + x31 + x41);
        const double E31 = x31 / (x11 + x21 + x31 + x41);
        const double E41 = x41 / (x11 + x21 + x31 + x41);
        const double JncxCa1 = E21 * c_k21 - E11 * c_k11;
        const double JncxNa1 = 3.0 * (E41 * k71 - E11 * k81) + E31 * k4pp1 - E21 * k3pp1;
        const double INaCa_ss_val = 0.2 * c_Gncx * allo1 * (c_zna * JncxNa1 + c_zca * JncxCa1);

        const double Knai = c_Knai0 * std::exp(c_delta * s_cell_v * c_F / (3.0 * c_R * c_T));
        const double Knao = c_Knao0 * std::exp((1.0 - c_delta) * s_cell_v * c_F / (3.0 * c_R * c_T));
        const double P_var = c_eP / (1.0 + c_H_const / c_Khp + s_nai / c_Knap + s_ki / c_Kxkur);
        const double a1_var = c_k1p * (s_nai / Knai * s_nai / Knai * (s_nai / Knai)) / ((1.0 + s_nai / Knai) * (1.0 + s_nai / Knai) * (1.0 + s_nai / Knai) + (1.0 + s_ki / c_Kki) * (1.0 + s_ki / c_Kki) - 1.0);
        const double a3_var = c_k3p2 * (c_ko / c_Kko * c_ko / c_Kko) / ((1.0 + c_nao / Knao) * (1.0 + c_nao / Knao) * (1.0 + c_nao / Knao) + (1.0 + c_ko / c_Kko) * (1.0 + c_ko / c_Kko) - 1.0);
        const double b2_var = c_k2m * (c_nao / Knao * c_nao / Knao * (c_nao / Knao)) / ((1.0 + c_nao / Knao) * (1.0 + c_nao / Knao) * (1.0 + c_nao / Knao) + (1.0 + c_ko / c_Kko) * (1.0 + c_ko / c_Kko) - 1.0);
        const double b3_var = c_k3m * P_var * c_H_const / (1.0 + c_MgATP / c_Kmgatp);
        const double b4_var = c_k4m * (s_ki / c_Kki * s_ki / c_Kki) / ((1.0 + s_nai / Knai) * (1.0 + s_nai / Knai) * (1.0 + s_nai / Knai) + (1.0 + s_ki / c_Kki) * (1.0 + s_ki / c_Kki) - 1.0);
        const double x12 = c_a4 * a1_var * c_a2 + b2_var * b4_var * b3_var + c_a2 * b4_var * b3_var + b3_var * a1_var * c_a2;
        const double x22 = b2_var * c_b1 * b4_var + a1_var * c_a2 * a3_var + a3_var * c_b1 * b4_var + c_a2 * a3_var * b4_var;
        const double x32 = c_a2 * a3_var * c_a4 + b3_var * b2_var * c_b1 + b2_var * c_b1 * c_a4 + a3_var * c_a4 * c_b1;
        const double x42 = b4_var * b3_var * b2_var + a3_var * c_a4 * a1_var + b2_var * c_a4 * a1_var + b3_var * b2_var * a1_var;
        const double E12 = x12 / (x12 + x22 + x32 + x42);
        const double E22 = x22 / (x12 + x22 + x32 + x42);
        const double E32 = x32 / (x12 + x22 + x32 + x42);
        const double E42 = x42 / (x12 + x22 + x32 + x42);
        const double JnakK = 2.0 * (E42 * c_b1 - E32 * a1_var);
        const double JnakNa = 3.0 * (E12 * a3_var - E22 * b3_var);
        const double INaK_val = c_Pnak * (c_zna * JnakNa + c_zk * JnakK);

        const double xkb = 1.0 / (1.0 + ::exp(-(s_cell_v - 14.48) / 18.34));
        const double IKb_val = c_GKb * xkb * (s_cell_v - EK);

        const double INab_val = c_PNab * c_F * (s_nai * std::exp(vfrt) - c_nao) * ((s_cell_v == 0.0) ? 1.0 : vfrt / (std::exp(vfrt) - 1.0));
        const double ICab_val = c_PCab * 2.0 * c_F * (s_cai * std::exp(2.0 * vfrt) - 0.341 * c_cao) * ((s_cell_v == 0.0) ? 1.0 : 2.0 * vfrt / (std::exp(2.0 * vfrt) - 1.0));
        const double IpCa_val = c_GpCa * s_cai / (0.0005 + s_cai);

        const double Jleak = 0.005 * s_cansr / 15.0;
        const double Jtr = (s_cansr - s_cajsr) / 100.0;
        const double Jtr2 = (s_cansr - s_cacsr) / 100.0;
        const double Jupnp = 1.0 * 0.005 * s_cai / (s_cai + 0.001);
        const double Jupnp2 = 1.0 * 0.005 * s_cai2 / (s_cai2 + 0.001);
        const double Jupp = 1.0 * 3.0 * 0.005 * s_cai / (s_cai + 0.001 - 0.0002);
        const double Jupp2 = 1.0 * 3.0 * 0.005 * s_cai2 / (s_cai2 + 0.001 - 0.0002);
        const double fJupp = 1.0 / (1.0 + c_KmCaMK / CaMKa);
        const double Jup = (1.0 - fJupp) * Jupnp + fJupp * Jupp - Jleak;
        const double Jup2 = (1.0 - fJupp) * Jupnp2 + fJupp * Jupp2 - Jleak;

        const double Bcacsr = 1.0 / (1.0 + c_csqnmax * c_kmcsqn / ((c_kmcsqn + s_cacsr) * (c_kmcsqn + s_cacsr)));
        const double Bcai = 1.0 / (1.0 + c_cmdnmax * c_kmcmdn / ((c_kmcmdn + s_cai) * (c_kmcmdn + s_cai)) + c_trpnmax * c_kmtrpn / ((c_kmtrpn + s_cai) * (c_kmtrpn + s_cai)));
        const double Bcai2 = 1.0 / (1.0 + c_cmdnmax * c_kmcmdn / ((c_kmcmdn + s_cai2) * (c_kmcmdn + s_cai2)) + c_trpnmax * c_kmtrpn / ((c_kmtrpn + s_cai2) * (c_kmtrpn + s_cai2)));
        const double Bcajsr = 1.0 / (1.0 + c_csqnmax * c_kmcsqn / ((c_kmcsqn + s_cajsr) * (c_kmcsqn + s_cajsr)));
        const double Bcass = 1.0 / (1.0 + c_BSRmax * c_KmBSR / ((c_KmBSR + s_cass) * (c_KmBSR + s_cass)) + c_BSLmax * c_KmBSL / ((c_KmBSL + s_cass) * (c_KmBSL + s_cass)));

        const double diff_cansr = Jup * (c_vnsr1 / c_vnsr) + Jup2 * (c_vnsr2 / c_vnsr) - Jtr * c_vjsr / c_vnsr - Jtr2 * c_vcsr / c_vnsr;
        const double diff_ki = -(IKr_val + IKs_val + IK1_val + IKb_val - 2.0 * INaK_val) * c_Acap / (2.0 * c_F * c_vmyo) + JdiffK * c_vss / c_vmyo;
        const double diff_kss = -ICaK * c_Acap / (2.0 * c_F * c_vss) - JdiffK;
        const double diff_nai = -(INa + INaL_val + 3.0 * INaCa_i_val + 3.0 * INaK_val + INab_val) * c_Acap / (2.0 * c_F * c_vmyo) + JdiffNa * c_vss / c_vmyo;
        const double diff_nass = -(ICaNa + 3.0 * INaCa_ss_val) * c_Acap / (2.0 * c_F * c_vss) - JdiffNa;

        const double diff_cacsr = Bcacsr * (Jtr2 - s_CICR_Jrel2);
        const double diff_cai = Bcai * (-(IpCa_val + ICab_val - 2.0 * INaCa_i_val) * c_Acap / (2.0 * 2.0 * c_F * c_vmyo1) - Jup * c_vnsr1 / c_vmyo1 + Jdiff * c_vss / c_vmyo1 - Jgap);
        const double diff_cai2 = Bcai2 * (s_CICR_Jrel2 * (c_vcsr / c_vmyo2) + Jgap * c_vmyo2 / c_vmyo1 - Jup2 * c_vnsr2 / c_vmyo2);
        const double diff_cajsr = Bcajsr * (Jtr - Jrel);
        const double diff_cass = Bcass * (-(ICaL_val - 2.0 * INaCa_ss_val) * c_Acap / (2.0 * 2.0 * c_F * c_vss) + Jrel * c_vjsr / c_vss - Jdiff);

        const double kito2 = 1.0 - 1.0 / (1.0 + Jrel / 0.4 * Jrel / 0.4);
        const double Rel1 = (-ICaL_val + 2.0 * INaCa_ss_val) * (c_Acap / (2.0 * c_vss * c_F)) + Jrel * (c_vjsr / c_vss) - Jdiff;
        const double Rel2 = Jgap + s_CICR_Jrel2 * (c_vcsr / c_vmyo2) - Jup2 * (c_vnsr2 / c_vmyo2);
        const double ITo_val = c_Gto * s_ITo_aa * rito2 * kito2 * (s_cell_v - c_ECl);
        const double Jrel1_inf = ((Rel1 > 0.0) ? 1.0 * 60.0 * Rel1 * (1.0 + 1.0 / (1.0 + std::pow(c_KmCaMK / CaMKa, 8.0))) / (1.0 + std::pow(0.75 / s_cajsr, 8.0)) : 0.0);
        const double Jrel2_inf = ((Rel2 > 0.0) ? 1.0 * (250.0 * Rel2 * (1.0 + 1.0 / (1.0 + std::pow(c_KmCaMK / CaMKa, 8.0)))) / (1.0 + std::pow(0.75 / s_cacsr, 8.0)) : 0.0);
        const double diff_Jrel1 = (Jrel1_inf - s_CICR_Jrel1) / tau_Jrel1;
        const double diff_Jrel2 = (Jrel2_inf - s_CICR_Jrel2) / tau_Jrel2;

        const double Iion_cm_val = INa + INaL_val + ICaL_val + ICaNa + ICaK + IKr_val + IKs_val + IK1_val + ITo_val + INaCa_i_val + INaCa_ss_val + INaK_val + INab_val + IKb_val + IpCa_val + ICab_val;

        const double r_cell_v = solveVmWithinODESolver ? (-Iion_cm_val - Istim) : 0.0;

        RATES_SoA[cell_v_BASE + cellI] = r_cell_v;
        RATES_SoA[nai_BASE + cellI] = diff_nai;
        RATES_SoA[nass_BASE + cellI] = diff_nass;
        RATES_SoA[ki_BASE + cellI] = diff_ki;
        RATES_SoA[kss_BASE + cellI] = diff_kss;
        RATES_SoA[cai_BASE + cellI] = diff_cai;
        RATES_SoA[cai2_BASE + cellI] = diff_cai2;
        RATES_SoA[cass_BASE + cellI] = diff_cass;
        RATES_SoA[cansr_BASE + cellI] = diff_cansr;
        RATES_SoA[cajsr_BASE + cellI] = diff_cajsr;
        RATES_SoA[cacsr_BASE + cellI] = diff_cacsr;
        RATES_SoA[I_Na_m_BASE + cellI] = tau_m > 0.0 ? (m_inf - s_I_Na_m)/tau_m : 0.0;
        RATES_SoA[I_Na_h_BASE + cellI] = tau_h > 0.0 ? (h_inf - s_I_Na_h)/tau_h : 0.0;
        RATES_SoA[I_Na_j_BASE + cellI] = tau_j > 0.0 ? (j_inf - s_I_Na_j)/tau_j : 0.0;
        RATES_SoA[INaL_ml_BASE + cellI] = tau_ml > 0.0 ? (ml_inf - s_INaL_ml)/tau_ml : 0.0;
        RATES_SoA[INaL_hl_BASE + cellI] = (hl_inf - s_INaL_hl)/c_tau_hl;
        RATES_SoA[ICaL_d_BASE + cellI] = tau_d > 0.0 ? (d_inf - s_ICaL_d)/tau_d : 0.0;
        RATES_SoA[ICaL_fca_BASE + cellI] = tau_fca > 0.0 ? (fca_inf - s_ICaL_fca)/tau_fca : 0.0;
        RATES_SoA[IKr_xr_BASE + cellI] = tau_xr > 0.0 ? (xr_inf - s_IKr_xr)/tau_xr : 0.0;
        RATES_SoA[IKs_xs1_BASE + cellI] = tau_xs1 > 0.0 ? (xs1_inf - s_IKs_xs1)/tau_xs1 : 0.0;
        RATES_SoA[IKs_xs2_BASE + cellI] = tau_xs2 > 0.0 ? (xs2_inf - s_IKs_xs2)/tau_xs2 : 0.0;
        RATES_SoA[ITo_aa_BASE + cellI] = alpha_aa * (1.0 - s_ITo_aa) - beta_aa * s_ITo_aa;
        RATES_SoA[CICR_Jrel2_BASE + cellI] = diff_Jrel2;
        RATES_SoA[CICR_Jrel1_BASE + cellI] = diff_Jrel1;
        RATES_SoA[CaMK_CaMKt_BASE + cellI] = diff_CaMKt;
        RATES_SoA[CICR_tjsrol_BASE + cellI] = diff_tjsrol;
        RATES_SoA[CICR_A_BASE + cellI] = diff_A;
        RATES_SoA[ICaL_fs_BASE + cellI] = tau_fs > 0.0 ? (fs_inf - s_ICaL_fs)/tau_fs : 0.0;
        RATES_SoA[ICaL_ff_BASE + cellI] = tau_ff > 0.0 ? (ff_inf - s_ICaL_ff)/tau_ff : 0.0;

        if (SUPPORT_SoA != nullptr)
        {
            SUPPORT_SoA[GAUR_BATCH_SUPPORT_tau_m    * N + cellI] = tau_m;
            SUPPORT_SoA[GAUR_BATCH_SUPPORT_tau_h    * N + cellI] = tau_h;
            SUPPORT_SoA[GAUR_BATCH_SUPPORT_tau_j    * N + cellI] = tau_j;
            SUPPORT_SoA[GAUR_BATCH_SUPPORT_tau_ml   * N + cellI] = tau_ml;
            SUPPORT_SoA[GAUR_BATCH_SUPPORT_tau_hl   * N + cellI] = c_tau_hl;
            SUPPORT_SoA[GAUR_BATCH_SUPPORT_tau_d    * N + cellI] = tau_d;
            SUPPORT_SoA[GAUR_BATCH_SUPPORT_tau_fca  * N + cellI] = tau_fca;
            SUPPORT_SoA[GAUR_BATCH_SUPPORT_tau_ff   * N + cellI] = tau_ff;
            SUPPORT_SoA[GAUR_BATCH_SUPPORT_tau_fs   * N + cellI] = tau_fs;
            SUPPORT_SoA[GAUR_BATCH_SUPPORT_tau_xr   * N + cellI] = tau_xr;
            SUPPORT_SoA[GAUR_BATCH_SUPPORT_tau_xs1  * N + cellI] = tau_xs1;
            SUPPORT_SoA[GAUR_BATCH_SUPPORT_tau_xs2  * N + cellI] = tau_xs2;
            SUPPORT_SoA[GAUR_BATCH_SUPPORT_tau_aa   * N + cellI] = 1.0/(alpha_aa + beta_aa);
            SUPPORT_SoA[GAUR_BATCH_SUPPORT_tau_Jrel1 * N + cellI] = tau_Jrel1;
            SUPPORT_SoA[GAUR_BATCH_SUPPORT_tau_Jrel2 * N + cellI] = tau_Jrel2;
            SUPPORT_SoA[GAUR_BATCH_SUPPORT_gInf_m    * N + cellI] = m_inf;
            SUPPORT_SoA[GAUR_BATCH_SUPPORT_gInf_h    * N + cellI] = h_inf;
            SUPPORT_SoA[GAUR_BATCH_SUPPORT_gInf_j    * N + cellI] = j_inf;
            SUPPORT_SoA[GAUR_BATCH_SUPPORT_gInf_ml   * N + cellI] = ml_inf;
            SUPPORT_SoA[GAUR_BATCH_SUPPORT_gInf_hl   * N + cellI] = hl_inf;
            SUPPORT_SoA[GAUR_BATCH_SUPPORT_gInf_d    * N + cellI] = d_inf;
            SUPPORT_SoA[GAUR_BATCH_SUPPORT_gInf_fca  * N + cellI] = fca_inf;
            SUPPORT_SoA[GAUR_BATCH_SUPPORT_gInf_ff   * N + cellI] = ff_inf;
            SUPPORT_SoA[GAUR_BATCH_SUPPORT_gInf_fs   * N + cellI] = fs_inf;
            SUPPORT_SoA[GAUR_BATCH_SUPPORT_gInf_xr   * N + cellI] = xr_inf;
            SUPPORT_SoA[GAUR_BATCH_SUPPORT_gInf_xs1  * N + cellI] = xs1_inf;
            SUPPORT_SoA[GAUR_BATCH_SUPPORT_gInf_xs2  * N + cellI] = xs2_inf;
            SUPPORT_SoA[GAUR_BATCH_SUPPORT_gInf_aa   * N + cellI] = alpha_aa/(alpha_aa + beta_aa);
            SUPPORT_SoA[GAUR_BATCH_SUPPORT_gInf_Jrel1 * N + cellI] = Jrel1_inf;
            SUPPORT_SoA[GAUR_BATCH_SUPPORT_gInf_Jrel2 * N + cellI] = Jrel2_inf;
            SUPPORT_SoA[GAUR_BATCH_SUPPORT_Iion_cm  * N + cellI] = Iion_cm_val;
        }
    }
}

} // End namespace Foam

// ************************************************************************* //
