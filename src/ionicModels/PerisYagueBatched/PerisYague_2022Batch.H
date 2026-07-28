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
    State-major (SoA) batched evaluator for the Peris-Yague 2022 ionic model.
Author
    Simao Nieto de Castro, UCD.
\*---------------------------------------------------------------------------*/

#pragma once

#include "PerisYague_2022Names.H"
#include "StimulusProtocolPOD.H"
#include <cmath>

namespace Foam
{

enum PERISYAGUE_BATCH_SUPPORT_INDEX
{
    PERISYAGUE_BATCH_SUPPORT_tau_m,
    PERISYAGUE_BATCH_SUPPORT_gInf_m,
    PERISYAGUE_BATCH_SUPPORT_tau_h,
    PERISYAGUE_BATCH_SUPPORT_gInf_h,
    PERISYAGUE_BATCH_SUPPORT_tau_j,
    PERISYAGUE_BATCH_SUPPORT_gInf_j,
    PERISYAGUE_BATCH_SUPPORT_tau_xr,
    PERISYAGUE_BATCH_SUPPORT_gInf_xr,
    PERISYAGUE_BATCH_SUPPORT_tau_xs,
    PERISYAGUE_BATCH_SUPPORT_gInf_xs,
    PERISYAGUE_BATCH_SUPPORT_tau_ua,
    PERISYAGUE_BATCH_SUPPORT_gInf_ua,
    PERISYAGUE_BATCH_SUPPORT_tau_uif,
    PERISYAGUE_BATCH_SUPPORT_gInf_uif,
    PERISYAGUE_BATCH_SUPPORT_tau_uis,
    PERISYAGUE_BATCH_SUPPORT_gInf_uis,
    PERISYAGUE_BATCH_SUPPORT_tau_d,
    PERISYAGUE_BATCH_SUPPORT_gInf_d,
    PERISYAGUE_BATCH_SUPPORT_tau_f,
    PERISYAGUE_BATCH_SUPPORT_gInf_f,
    PERISYAGUE_BATCH_SUPPORT_tau_fCa,
    PERISYAGUE_BATCH_SUPPORT_gInf_fCa,
    PERISYAGUE_BATCH_SUPPORT_tau_qCa,
    PERISYAGUE_BATCH_SUPPORT_gInf_qCa,
    PERISYAGUE_BATCH_SUPPORT_tau_ryr_u,
    PERISYAGUE_BATCH_SUPPORT_gInf_ryr_u,
    PERISYAGUE_BATCH_SUPPORT_tau_ryr_v,
    PERISYAGUE_BATCH_SUPPORT_gInf_ryr_v,
    PERISYAGUE_BATCH_SUPPORT_tau_ryr_w,
    PERISYAGUE_BATCH_SUPPORT_gInf_ryr_w,
    PERISYAGUE_BATCH_SUPPORT_Iion_cm,

    NUM_PERISYAGUE_BATCH_SUPPORT
};

static_assert(NUM_PERISYAGUE_BATCH_SUPPORT == 31, "PerisYague compact support layout drift");

CARDIAC_HOST_DEVICE inline void PerisYagueComputeVariablesBatch
(
    const double t,
    const double* __restrict__ CONSTANTS,
    const int N,
    const int beginCell,
    const int endCell,
    const double* __restrict__ STATES_SoA,
    double*       __restrict__ RATES_SoA,
    double*       __restrict__ SUPPORT_SoA,
    const bool solveVm,
    const StimulusProtocolPOD& stimulus
)
{
    const double Istim = computeStimulusPOD(t, stimulus);

    const double c_RTF          = CONSTANTS[AC_RTF];
    const double c_FRT          = CONSTANTS[AC_FRT];
    const double c_Nao          = CONSTANTS[AC_Nao];
    const double c_Ko           = CONSTANTS[AC_Ko];
    const double c_Cao          = CONSTANTS[AC_Cao];
    const double c_Clo          = CONSTANTS[AC_Clo];
    const double c_Cm           = CONSTANTS[AC_Cm];
    const double c_F            = CONSTANTS[AC_F];
    const double c_V_i          = CONSTANTS[AC_V_i];
    const double c_V_up         = CONSTANTS[AC_V_up];
    const double c_V_rel        = CONSTANTS[AC_V_rel];

    const double c_gNa          = CONSTANTS[AC_gNa];
    const double c_gK1          = CONSTANTS[AC_gK1];
    const double c_gKr          = CONSTANTS[AC_gKr];
    const double c_gKs          = CONSTANTS[AC_gKs];
    const double c_gKur_amp     = CONSTANTS[AC_gKur_amp];
    const double c_gCaL         = CONSTANTS[AC_gCaL];
    const double c_gbNa         = CONSTANTS[AC_gbNa];
    const double c_gbCa         = CONSTANTS[AC_gbCa];
    const double c_gbK          = CONSTANTS[AC_gbK];
    const double c_gClCa        = CONSTANTS[AC_gClCa];
    const double c_KQ10         = CONSTANTS[AC_KQ10];
    const double c_ical_fCa_tau = CONSTANTS[AC_ical_fCa_tau];

    const double c_INaCa_max    = CONSTANTS[AC_INaCa_max];
    const double c_INaK_max     = CONSTANTS[AC_INaK_max];
    const double c_IpCa_max     = CONSTANTS[AC_IpCa_max];
    const double c_gamma        = CONSTANTS[AC_gamma];
    const double c_KmNa         = CONSTANTS[AC_KmNa];
    const double c_KmKo         = CONSTANTS[AC_KmKo];
    const double c_KmNai        = CONSTANTS[AC_KmNai];
    const double c_KmCa         = CONSTANTS[AC_KmCa];
    const double c_ksat         = CONSTANTS[AC_ksat];
    const double c_sigma        = CONSTANTS[AC_sigma];

    const double c_krel         = CONSTANTS[AC_krel];
    const double c_K_up         = CONSTANTS[AC_K_up];
    const double c_I_up_max     = CONSTANTS[AC_I_up_max];
    const double c_Ca_up_max    = CONSTANTS[AC_Ca_up_max];
    const double c_tau_tr       = CONSTANTS[AC_tau_tr];
    const double c_c1           = CONSTANTS[AC_c1];
    const double c_c2           = CONSTANTS[AC_c2];
    const double c_cajsr_u_tau  = CONSTANTS[AC_cajsr_u_tau];

    const double c_CMDN_max     = CONSTANTS[AC_CMDN_max];
    const double c_Km_CMDN      = CONSTANTS[AC_Km_CMDN];
    const double c_TRPN_max     = CONSTANTS[AC_TRPN_max];
    const double c_Km_TRPN      = CONSTANTS[AC_Km_TRPN];
    const double c_CSQN_max     = CONSTANTS[AC_CSQN_max];
    const double c_Km_CSQN      = CONSTANTS[AC_Km_CSQN];

    const int V_BASE     = membrane_V    * N;
    const int Nai_BASE   = sodium_Nai    * N;
    const int Ki_BASE    = potassium_Ki  * N;
    const int Cli_BASE   = chloride_Cli  * N;
    const int Cai_BASE   = calcium_Cai   * N;
    const int CaUp_BASE  = calcium_CaUp  * N;
    const int CaRel_BASE = calcium_CaRel * N;
    const int m_BASE     = ina_m         * N;
    const int h_BASE     = ina_h         * N;
    const int j_BASE     = ina_j         * N;
    const int ua_BASE    = ikur_ua       * N;
    const int uif_BASE   = ikur_uif      * N;
    const int uis_BASE   = ikur_uis      * N;
    const int xr_BASE    = ikr_xr        * N;
    const int xs_BASE    = iks_xs        * N;
    const int d_BASE     = ical_d        * N;
    const int f_BASE     = ical_f        * N;
    const int fCa_BASE   = ical_fCa      * N;
    const int qCa_BASE   = iclca_qCa     * N;
    const int ryr_u_BASE = ryr_u         * N;
    const int ryr_v_BASE = ryr_v         * N;
    const int ryr_w_BASE = ryr_w         * N;

#ifdef _OPENMP
    #pragma omp simd
#endif
    for (int cellI = beginCell; cellI < endCell; ++cellI)
    {
        const double V     = STATES_SoA[V_BASE     + cellI];
        const double Nai   = STATES_SoA[Nai_BASE   + cellI];
        const double Ki    = STATES_SoA[Ki_BASE    + cellI];
        const double Cli   = STATES_SoA[Cli_BASE   + cellI];
        const double Cai   = STATES_SoA[Cai_BASE   + cellI];
        const double CaUp  = STATES_SoA[CaUp_BASE  + cellI];
        const double CaRel = STATES_SoA[CaRel_BASE + cellI];
        const double s_m   = STATES_SoA[m_BASE     + cellI];
        const double s_h   = STATES_SoA[h_BASE     + cellI];
        const double s_j   = STATES_SoA[j_BASE     + cellI];
        const double s_ua  = STATES_SoA[ua_BASE    + cellI];
        const double s_uif = STATES_SoA[uif_BASE   + cellI];
        const double s_uis = STATES_SoA[uis_BASE   + cellI];
        const double s_xr  = STATES_SoA[xr_BASE    + cellI];
        const double s_xs  = STATES_SoA[xs_BASE    + cellI];
        const double s_d   = STATES_SoA[d_BASE     + cellI];
        const double s_f   = STATES_SoA[f_BASE     + cellI];
        const double s_fCa = STATES_SoA[fCa_BASE   + cellI];
        const double s_qCa = STATES_SoA[qCa_BASE   + cellI];
        const double s_u   = STATES_SoA[ryr_u_BASE + cellI];
        const double s_v   = STATES_SoA[ryr_v_BASE + cellI];
        const double s_w   = STATES_SoA[ryr_w_BASE + cellI];

        const double ENa = c_RTF * std::log(c_Nao / Nai);
        const double EK  = c_RTF * std::log(c_Ko  / Ki);
        const double ECa = 0.5 * c_RTF * std::log(c_Cao / Cai);
        const double ECl = -c_RTF * std::log(c_Clo / Cli);

        const double am = (std::fabs(V + 47.13) < 1e-6)
                          ? 3.2
                          : 0.32 * (V + 47.13) / (1.0 - std::exp(-0.1 * (V + 47.13)));
        const double bm = 0.08 * std::exp(-V / 11.0);
        const double m_inf = am / (am + bm);
        const double m_tau = 1.7 / (am + bm);

        const double ah = (V < -40.0)
                          ? 0.135 * std::exp(-(V + 80.0) / 6.8)
                          : 0.0;
        const double bh = (V < -40.0)
                          ? 3.56 * std::exp(0.079 * V) + 3.1e5 * std::exp(0.35 * V)
                          : 0.13 / (1.0 + std::exp(-(V + 10.66) / 11.1));
        const double h_inf = ah / (ah + bh);
        const double h_tau = 2.0 / (ah + bh);

        const double aj = (V < -40.0)
                          ? (-127140.0 * std::exp(0.2444 * V) - 3.474e-5 * std::exp(-0.04391 * V))
                            * (V + 37.78) / (1.0 + std::exp(0.311 * (V + 79.23)))
                          : 0.0;
        const double bj = (V < -40.0)
                          ? 0.1212 * std::exp(-0.01052 * V) / (1.0 + std::exp(-0.1378 * (V + 40.14)))
                          : 0.3 * std::exp(-2.535e-7 * V) / (1.0 + std::exp(-0.1 * (V + 32.0)));
        const double j_inf = aj / (aj + bj);
        const double j_tau = 2.0 / (aj + bj);

        const double INa = c_gNa * s_m * s_m * s_m * s_h * s_j * (V - ENa);

        const double IK1 = c_gK1 * (V - EK - 5.0) / (1.0 + std::exp(0.063 * (V + 70.0)));

        const double axr = 0.0003 * (V + 14.1) / (1.0 - std::exp(-(V + 14.1) / 5.0));
        const double bxr = 7.3898e-5 * (V - 3.3328) / (std::exp((V - 3.3328) / 5.1237) - 1.0);
        const double xr_inf = axr / (axr + bxr);
        const double xr_tau = 1.0 / (axr + bxr);
        const double IKr = c_gKr * s_xr * (V - EK) / (1.0 + std::exp((V + 15.0) / 22.4));

        const double axs = 4e-5 * (V - 18.80816) / (1.0 - std::exp(-(V - 18.80816) / 17.0));
        const double bxs = 3.5e-5 * (V - 18.80816) / (std::exp((V - 18.80816) / 9.0) - 1.0);
        const double xs_inf = 1.0 / std::sqrt(1.0 + std::exp(-(V - 19.9) / 12.7));
        const double xs_tau = 0.5 / (axs + bxs);
        const double IKs = c_gKs * s_xs * s_xs * (V - EK);

        const double gKur = c_gKur_amp * (0.005 + 0.05 / (1.0 + std::exp(-(V - 15.0) / 13.0)));
        const double ua_inf = 1.0 / (1.0 + std::exp(-(V + 30.3) / 9.6));
        const double aua = 0.65 / (std::exp(-(V + 10.0) / 8.5) + std::exp(-(V - 30.0) / 59.0));
        const double bua = 0.65 / (2.5 + std::exp((V + 82.0) / 17.0));
        const double ua_tau = 1.0 / (c_KQ10 * (aua + bua));

        const double uif_inf = 1.0 / (1.0 + std::exp((V + 17.5) / 5.85));
        const double V_over_50 = V / 50.0;
        const double uif_tau = 400.0 + 1068.0 * std::exp(-V_over_50 * V_over_50);

        const double uis_inf = uif_inf;
        const double V_plus_39_3_over_30 = (V + 39.3) / 30.0;
        const double uis_tau =
            2000.0 + 60000.0 * std::exp(-V_plus_39_3_over_30 * V_plus_39_3_over_30);

        const double IKur =
            gKur * s_ua * s_ua * s_ua * (0.25 * s_uif + 0.75 * s_uis) * (V - EK);

        const double d_inf = 1.0 / (1.0 + std::exp(-(V + 5.0) / 8.0));
        const double V_plus_5 = V + 5.0;
        const double d_tau = (std::fabs(V_plus_5) < 1e-6)
                             ? 1.0
                             : (1.0 - std::exp(-V_plus_5 / 6.24))
                               / (0.035 * V_plus_5 * (1.0 + std::exp(-V_plus_5 / 6.24)));

        const double f_inf = 1.0 / (1.0 + std::exp((V + 28.0) / 6.9));
        const double f_tau = 9.0
                             / (0.0197 * std::exp(-0.03372 * V_plus_5 * V_plus_5) + 0.02);

        const double fCa_inf = 1.0 / (1.0 + Cai / 0.00035);

        const double ICaL = c_gCaL * s_d * s_f * s_fCa * (V - 65.0);

        const double fNaK = 1.0
                            / (1.0 + 0.1245 * std::exp(-0.1 * V * c_FRT)
                               + 0.0365 * c_sigma * std::exp(-V * c_FRT));
        const double INaK =
            c_INaK_max * fNaK * c_Ko
            / (c_Ko + c_KmKo)
            / (1.0 + std::pow(c_KmNai / Nai, 1.5));

        const double Nai3  = Nai * Nai * Nai;
        const double Nao3  = c_Nao * c_Nao * c_Nao;
        const double KmNa3 = c_KmNa * c_KmNa * c_KmNa;
        const double B1_NaCa =
            std::exp(c_gamma * V * c_FRT) * Nai3 * c_Cao
          - std::exp((c_gamma - 1.0) * V * c_FRT) * Nao3 * Cai;
        const double B2_NaCa =
            (KmNa3 + Nao3) * (c_KmCa + c_Cao)
            * (1.0 + c_ksat * std::exp((c_gamma - 1.0) * V * c_FRT));
        const double INaCa = c_INaCa_max * B1_NaCa / B2_NaCa;

        const double IbNa = c_gbNa * (V - ENa);
        const double IbCa = c_gbCa * (V - ECa);
        const double IbK  = c_gbK  * (V - EK);
        const double IpCa = c_IpCa_max * Cai / (0.0005 + Cai);

        const double Irel     = c_krel * s_u * s_u * s_v * s_w * (CaRel - Cai);
        const double Itr      = (CaUp - CaRel) / c_tau_tr;
        const double Iup      = c_I_up_max / (1.0 + c_K_up / Cai);
        const double Iup_leak = c_I_up_max * CaUp / c_Ca_up_max;

        const double Fn =
            1e-12 * c_V_rel * Irel
          - 5e-13 / c_F * (0.5 * ICaL - 0.2 * INaCa) * c_Cm;

        const double inv_qCa_arg = 1.0 / (1.0 + Fn / 1.1e-10);
        const double qCa_inf = 1.0 - inv_qCa_arg * inv_qCa_arg * inv_qCa_arg;
        const double qCa_tau = 2.0;
        const double IClCa   = c_gClCa * s_qCa * (V - ECl);

        const double w_inf = 1.0 - 1.0 / (1.0 + std::exp(-(V - 40.0) / 17.0));
        const double V_minus_7_9 = V - 7.9;
        const double w_tau = 6.0
                             * ((std::fabs(V_minus_7_9) < 1e-6)
                                ? 2.0 / 13.0
                                : (1.0 - std::exp(-V_minus_7_9 / 5.0))
                                  / ((1.0 + 0.3 * std::exp(-V_minus_7_9 / 5.0)) * V_minus_7_9));

        const double u_inf_ryr = 1.0 / (1.0 + std::exp(-(Fn - c_c1) / c_c2));
        const double v_inf_ryr = 1.0 - 1.0 / (1.0 + std::exp(-(Fn - 0.2 * c_c1) / c_c2));
        const double v_tau_ryr = 1.91 + 2.09 / (1.0 + std::exp(-(Fn - c_c1) / c_c2));

        const double Cai_plus_KmTRPN = Cai + c_Km_TRPN;
        const double Cai_plus_KmCMDN = Cai + c_Km_CMDN;
        const double B1 =
            (2.0 * INaCa - (IpCa + ICaL + IbCa)) * c_Cm / (2.0 * c_V_i * c_F)
          + (c_V_up * (Iup_leak - Iup) + Irel * c_V_rel) / c_V_i;
        const double B2 =
            1.0
          + c_TRPN_max * c_Km_TRPN / (Cai_plus_KmTRPN * Cai_plus_KmTRPN)
          + c_CMDN_max * c_Km_CMDN / (Cai_plus_KmCMDN * Cai_plus_KmCMDN);

        const double CaRel_plus_KmCSQN = CaRel + c_Km_CSQN;

        const double Iion_cm_val =
            INa + IK1 + IKur + IKr + IKs + ICaL + IpCa
          + INaK + INaCa + IbNa + IbCa + IClCa;

        RATES_SoA[m_BASE     + cellI] = (m_inf  - s_m)   / m_tau;
        RATES_SoA[h_BASE     + cellI] = (h_inf  - s_h)   / h_tau;
        RATES_SoA[j_BASE     + cellI] = (j_inf  - s_j)   / j_tau;
        RATES_SoA[xr_BASE    + cellI] = (xr_inf - s_xr)  / xr_tau;
        RATES_SoA[xs_BASE    + cellI] = (xs_inf - s_xs)  / xs_tau;
        RATES_SoA[ua_BASE    + cellI] = (ua_inf - s_ua)  / ua_tau;
        RATES_SoA[uif_BASE   + cellI] = (uif_inf - s_uif) / uif_tau;
        RATES_SoA[uis_BASE   + cellI] = (uis_inf - s_uis) / uis_tau;
        RATES_SoA[d_BASE     + cellI] = (d_inf  - s_d)   / d_tau;
        RATES_SoA[f_BASE     + cellI] = (f_inf  - s_f)   / f_tau;
        RATES_SoA[fCa_BASE   + cellI] = (fCa_inf - s_fCa) / c_ical_fCa_tau;
        RATES_SoA[qCa_BASE   + cellI] = (qCa_inf - s_qCa) / qCa_tau;
        RATES_SoA[ryr_u_BASE + cellI] = (u_inf_ryr - s_u) / c_cajsr_u_tau;
        RATES_SoA[ryr_v_BASE + cellI] = (v_inf_ryr - s_v) / v_tau_ryr;
        RATES_SoA[ryr_w_BASE + cellI] = (w_inf  - s_w)   / w_tau;

        RATES_SoA[Cai_BASE   + cellI] = B1 / B2;
        RATES_SoA[CaRel_BASE + cellI] =
            (Itr - Irel)
            / (1.0 + c_CSQN_max * c_Km_CSQN
                     / (CaRel_plus_KmCSQN * CaRel_plus_KmCSQN));
        RATES_SoA[CaUp_BASE  + cellI] = Iup - (Iup_leak + Itr * c_V_rel / c_V_up);

        RATES_SoA[Nai_BASE   + cellI] =
            (-3.0 * INaK - (3.0 * INaCa + IbNa + INa)) * c_Cm / (c_V_i * c_F);
        RATES_SoA[Ki_BASE    + cellI] =
            (2.0 * INaK - (IK1 + IKur + IKr + IKs + IbK)) * c_Cm / (c_V_i * c_F);
        RATES_SoA[Cli_BASE   + cellI] = IClCa * c_Cm / (c_V_i * c_F);

        if (solveVm)
        {
            RATES_SoA[V_BASE + cellI] = -Iion_cm_val - Istim;
        }

        if (SUPPORT_SoA != nullptr)
        {
            SUPPORT_SoA[PERISYAGUE_BATCH_SUPPORT_tau_m     * N + cellI] = m_tau;
            SUPPORT_SoA[PERISYAGUE_BATCH_SUPPORT_tau_h     * N + cellI] = h_tau;
            SUPPORT_SoA[PERISYAGUE_BATCH_SUPPORT_tau_j     * N + cellI] = j_tau;
            SUPPORT_SoA[PERISYAGUE_BATCH_SUPPORT_tau_xr    * N + cellI] = xr_tau;
            SUPPORT_SoA[PERISYAGUE_BATCH_SUPPORT_tau_xs    * N + cellI] = xs_tau;
            SUPPORT_SoA[PERISYAGUE_BATCH_SUPPORT_tau_ua    * N + cellI] = ua_tau;
            SUPPORT_SoA[PERISYAGUE_BATCH_SUPPORT_tau_uif   * N + cellI] = uif_tau;
            SUPPORT_SoA[PERISYAGUE_BATCH_SUPPORT_tau_uis   * N + cellI] = uis_tau;
            SUPPORT_SoA[PERISYAGUE_BATCH_SUPPORT_tau_d     * N + cellI] = d_tau;
            SUPPORT_SoA[PERISYAGUE_BATCH_SUPPORT_tau_f     * N + cellI] = f_tau;
            SUPPORT_SoA[PERISYAGUE_BATCH_SUPPORT_tau_fCa   * N + cellI] = c_ical_fCa_tau;
            SUPPORT_SoA[PERISYAGUE_BATCH_SUPPORT_tau_qCa   * N + cellI] = qCa_tau;
            SUPPORT_SoA[PERISYAGUE_BATCH_SUPPORT_tau_ryr_u * N + cellI] = c_cajsr_u_tau;
            SUPPORT_SoA[PERISYAGUE_BATCH_SUPPORT_tau_ryr_v * N + cellI] = v_tau_ryr;
            SUPPORT_SoA[PERISYAGUE_BATCH_SUPPORT_tau_ryr_w * N + cellI] = w_tau;
            SUPPORT_SoA[PERISYAGUE_BATCH_SUPPORT_gInf_m     * N + cellI] = m_inf;
            SUPPORT_SoA[PERISYAGUE_BATCH_SUPPORT_gInf_h     * N + cellI] = h_inf;
            SUPPORT_SoA[PERISYAGUE_BATCH_SUPPORT_gInf_j     * N + cellI] = j_inf;
            SUPPORT_SoA[PERISYAGUE_BATCH_SUPPORT_gInf_xr    * N + cellI] = xr_inf;
            SUPPORT_SoA[PERISYAGUE_BATCH_SUPPORT_gInf_xs    * N + cellI] = xs_inf;
            SUPPORT_SoA[PERISYAGUE_BATCH_SUPPORT_gInf_ua    * N + cellI] = ua_inf;
            SUPPORT_SoA[PERISYAGUE_BATCH_SUPPORT_gInf_uif   * N + cellI] = uif_inf;
            SUPPORT_SoA[PERISYAGUE_BATCH_SUPPORT_gInf_uis   * N + cellI] = uis_inf;
            SUPPORT_SoA[PERISYAGUE_BATCH_SUPPORT_gInf_d     * N + cellI] = d_inf;
            SUPPORT_SoA[PERISYAGUE_BATCH_SUPPORT_gInf_f     * N + cellI] = f_inf;
            SUPPORT_SoA[PERISYAGUE_BATCH_SUPPORT_gInf_fCa   * N + cellI] = fCa_inf;
            SUPPORT_SoA[PERISYAGUE_BATCH_SUPPORT_gInf_qCa   * N + cellI] = qCa_inf;
            SUPPORT_SoA[PERISYAGUE_BATCH_SUPPORT_gInf_ryr_u * N + cellI] = u_inf_ryr;
            SUPPORT_SoA[PERISYAGUE_BATCH_SUPPORT_gInf_ryr_v * N + cellI] = v_inf_ryr;
            SUPPORT_SoA[PERISYAGUE_BATCH_SUPPORT_gInf_ryr_w * N + cellI] = w_inf;
            SUPPORT_SoA[PERISYAGUE_BATCH_SUPPORT_Iion_cm * N + cellI] = Iion_cm_val;
        }
    }
}

} // End namespace Foam

// ************************************************************************* //
