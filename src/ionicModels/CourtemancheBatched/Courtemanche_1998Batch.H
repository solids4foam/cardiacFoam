/*---------------------------------------------------------------------------*\
License
    This file is part of cardiacFoam.

Description
    State-major SoA batched evaluator for the Courtemanche 1998 ionic model.
\*---------------------------------------------------------------------------*/

#pragma once

#include "Courtemanche_1998Names.H"
#include "StimulusProtocolPOD.H"
#include <cmath>

namespace Foam
{

enum COURTEMANCHE_BATCH_SUPPORT_INDEX
{
    COURTEMANCHE_BATCH_SUPPORT_tau_m,
    COURTEMANCHE_BATCH_SUPPORT_gInf_m,
    COURTEMANCHE_BATCH_SUPPORT_tau_h,
    COURTEMANCHE_BATCH_SUPPORT_gInf_h,
    COURTEMANCHE_BATCH_SUPPORT_tau_j,
    COURTEMANCHE_BATCH_SUPPORT_gInf_j,
    COURTEMANCHE_BATCH_SUPPORT_tau_d,
    COURTEMANCHE_BATCH_SUPPORT_gInf_d,
    COURTEMANCHE_BATCH_SUPPORT_tau_f,
    COURTEMANCHE_BATCH_SUPPORT_gInf_f,
    COURTEMANCHE_BATCH_SUPPORT_tau_fCa,
    COURTEMANCHE_BATCH_SUPPORT_gInf_fCa,
    COURTEMANCHE_BATCH_SUPPORT_tau_oa,
    COURTEMANCHE_BATCH_SUPPORT_gInf_oa,
    COURTEMANCHE_BATCH_SUPPORT_tau_oi,
    COURTEMANCHE_BATCH_SUPPORT_gInf_oi,
    COURTEMANCHE_BATCH_SUPPORT_tau_ua,
    COURTEMANCHE_BATCH_SUPPORT_gInf_ua,
    COURTEMANCHE_BATCH_SUPPORT_tau_ui,
    COURTEMANCHE_BATCH_SUPPORT_gInf_ui,
    COURTEMANCHE_BATCH_SUPPORT_tau_xr,
    COURTEMANCHE_BATCH_SUPPORT_gInf_xr,
    COURTEMANCHE_BATCH_SUPPORT_tau_xs,
    COURTEMANCHE_BATCH_SUPPORT_gInf_xs,
    COURTEMANCHE_BATCH_SUPPORT_tau_cajsr_u,
    COURTEMANCHE_BATCH_SUPPORT_gInf_cajsr_u,
    COURTEMANCHE_BATCH_SUPPORT_tau_cajsr_v,
    COURTEMANCHE_BATCH_SUPPORT_gInf_cajsr_v,
    COURTEMANCHE_BATCH_SUPPORT_tau_cajsr_w,
    COURTEMANCHE_BATCH_SUPPORT_gInf_cajsr_w,
    COURTEMANCHE_BATCH_SUPPORT_Iion_cm,

    NUM_COURTEMANCHE_BATCH_SUPPORT
};

static_assert(NUM_COURTEMANCHE_BATCH_SUPPORT == 31, "Courtemanche compact support layout drift");

CARDIAC_HOST_DEVICE inline void CourtemancheComputeVariablesBatch
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

    const double c_CMDN_max     = CONSTANTS[AC_CMDN_max];
    const double c_CSQN_max     = CONSTANTS[AC_CSQN_max];
    const double c_Ca_up_max    = CONSTANTS[AC_Ca_up_max];
    const double c_Cao          = CONSTANTS[AC_Cao];
    const double c_Cm           = CONSTANTS[AC_Cm];
    const double c_ECaL         = CONSTANTS[AC_ECaL];
    const double c_F            = CONSTANTS[AC_F];
    const double c_FRT          = CONSTANTS[AC_FRT];
    const double c_INaCa_max    = CONSTANTS[AC_INaCa_max];
    const double c_INaK_max     = CONSTANTS[AC_INaK_max];
    const double c_I_up_max     = CONSTANTS[AC_I_up_max];
    const double c_IpCa_max     = CONSTANTS[AC_IpCa_max];
    const double c_KQ10         = CONSTANTS[AC_KQ10];
    const double c_K_rel        = CONSTANTS[AC_K_rel];
    const double c_K_up         = CONSTANTS[AC_K_up];
    const double c_KmCa         = CONSTANTS[AC_KmCa];
    const double c_KmKo         = CONSTANTS[AC_KmKo];
    const double c_KmNa         = CONSTANTS[AC_KmNa];
    const double c_KmNai        = CONSTANTS[AC_KmNai];
    const double c_Km_CMDN      = CONSTANTS[AC_Km_CMDN];
    const double c_Km_CSQN      = CONSTANTS[AC_Km_CSQN];
    const double c_Km_TRPN      = CONSTANTS[AC_Km_TRPN];
    const double c_Ko           = CONSTANTS[AC_Ko];
    const double c_Nao          = CONSTANTS[AC_Nao];
    const double c_RTF          = CONSTANTS[AC_RTF];
    const double c_TRPN_max     = CONSTANTS[AC_TRPN_max];
    const double c_V_i          = CONSTANTS[AC_V_i];
    const double c_V_rel        = CONSTANTS[AC_V_rel];
    const double c_V_up         = CONSTANTS[AC_V_up];
    const double c_c1           = CONSTANTS[AC_c1];
    const double c_c2           = CONSTANTS[AC_c2];
    const double c_cajsr_u_tau  = CONSTANTS[AC_cajsr_u_tau];
    const double c_g            = CONSTANTS[AC_g];
    const double c_gCaL         = CONSTANTS[AC_gCaL];
    const double c_gK1          = CONSTANTS[AC_gK1];
    const double c_gKr          = CONSTANTS[AC_gKr];
    const double c_gKs          = CONSTANTS[AC_gKs];
    const double c_gKur_base    = CONSTANTS[AC_gKur_base];
    const double c_gNa          = CONSTANTS[AC_gNa];
    const double c_gbCa         = CONSTANTS[AC_gbCa];
    const double c_gbNa         = CONSTANTS[AC_gbNa];
    const double c_gto          = CONSTANTS[AC_gto];
    const double c_ical_fCa_tau = CONSTANTS[AC_ical_fCa_tau];
    const double c_ksat         = CONSTANTS[AC_ksat];
    const double c_sigma        = CONSTANTS[AC_sigma];
    const double c_tau_tr       = CONSTANTS[AC_tau_tr];

    const int V_BASE     = membrane_V    * N;
    const int Nai_BASE   = sodium_Nai    * N;
    const int Ki_BASE    = potassium_Ki  * N;
    const int Cai_BASE   = calcium_Cai   * N;
    const int CaUp_BASE  = calcium_CaUp  * N;
    const int CaRel_BASE = calcium_CaRel * N;
    const int m_BASE     = ina_m         * N;
    const int h_BASE     = ina_h         * N;
    const int j_BASE     = ina_j         * N;
    const int oa_BASE    = ito_oa        * N;
    const int oi_BASE    = ito_oi        * N;
    const int ua_BASE    = ikur_ua       * N;
    const int ui_BASE    = ikur_ui       * N;
    const int xr_BASE    = ikr_xr        * N;
    const int xs_BASE    = iks_xs        * N;
    const int d_BASE     = ical_d        * N;
    const int f_BASE     = ical_f        * N;
    const int fCa_BASE   = ical_fCa      * N;
    const int u_BASE     = cajsr_u       * N;
    const int v_BASE     = cajsr_v       * N;
    const int w_BASE     = cajsr_w       * N;

#ifdef _OPENMP
    #pragma omp simd
#endif
    for (int cellI = beginCell; cellI < endCell; ++cellI)
    {
        const double V     = STATES_SoA[V_BASE     + cellI];
        const double Nai   = STATES_SoA[Nai_BASE   + cellI];
        const double Ki    = STATES_SoA[Ki_BASE    + cellI];
        const double Cai   = STATES_SoA[Cai_BASE   + cellI];
        const double CaUp  = STATES_SoA[CaUp_BASE  + cellI];
        const double CaRel = STATES_SoA[CaRel_BASE + cellI];
        const double s_m   = STATES_SoA[m_BASE     + cellI];
        const double s_h   = STATES_SoA[h_BASE     + cellI];
        const double s_j   = STATES_SoA[j_BASE     + cellI];
        const double s_oa  = STATES_SoA[oa_BASE    + cellI];
        const double s_oi  = STATES_SoA[oi_BASE    + cellI];
        const double s_ua  = STATES_SoA[ua_BASE    + cellI];
        const double s_ui  = STATES_SoA[ui_BASE    + cellI];
        const double s_xr  = STATES_SoA[xr_BASE    + cellI];
        const double s_xs  = STATES_SoA[xs_BASE    + cellI];
        const double s_d   = STATES_SoA[d_BASE     + cellI];
        const double s_f   = STATES_SoA[f_BASE     + cellI];
        const double s_fCa = STATES_SoA[fCa_BASE   + cellI];
        const double s_u   = STATES_SoA[u_BASE     + cellI];
        const double s_v   = STATES_SoA[v_BASE     + cellI];
        const double s_w   = STATES_SoA[w_BASE     + cellI];

        const double Iup      = c_I_up_max / (1.0 + c_K_up / Cai);
        const double IupLeak  = c_I_up_max * CaUp / c_Ca_up_max;
        const double Itr      = (CaUp - CaRel) / c_tau_tr;

        const double d_inf = 1.0 / (1.0 + std::exp((V + 10.0) / -8.0));
        const double d_tau =
            (std::fabs(V + 10.0) < 1e-6)
          ? 1.0 / (6.24 * 2.0 * 0.035)
          : (1.0 - std::exp((V + 10.0) / -6.24))
            / (0.035 * (V + 10.0) * (1.0 + std::exp((V + 10.0) / -6.24)));

        const double f_inf = 1.0 / (1.0 + std::exp((V + 28.0) / 6.9));
        const double f_tau =
            9.0 / (0.0197 * std::exp(-std::pow(0.0337, 2.0) * std::pow(V + 10.0, 2.0)) + 0.02);
        const double fCa_inf = 1.0 / (1.0 + Cai / 0.00035);
        const double ICaL = c_gCaL * s_d * s_f * s_fCa * (V - c_ECaL);

        const double IpCa = c_IpCa_max * Cai / (0.0005 + Cai);
        const double INaCa =
            c_INaCa_max
          * (std::exp(c_g * V * c_FRT) * std::pow(Nai, 3.0) * c_Cao
             - std::exp((c_g - 1.0) * V * c_FRT) * std::pow(c_Nao, 3.0) * Cai)
          / ((std::pow(c_KmNa, 3.0) + std::pow(c_Nao, 3.0))
             * (c_KmCa + c_Cao)
             * (1.0 + c_ksat * std::exp((c_g - 1.0) * V * c_FRT)));

        const double fNaK =
            1.0 / (1.0 + 0.1245 * std::exp(-0.1 * V * c_FRT)
                 + 0.0365 * c_sigma * std::exp(-V * c_FRT));
        const double INaK =
            c_INaK_max * fNaK * c_Ko / (c_Ko + c_KmKo)
          / (1.0 + std::pow(c_KmNai / Nai, 1.5));

        const double ECa = 0.5 * c_RTF * std::log(c_Cao / Cai);
        const double EK  = c_RTF * std::log(c_Ko / Ki);
        const double ENa = c_RTF * std::log(c_Nao / Nai);

        const double Irel = c_K_rel * s_u * s_u * s_v * s_w * (CaRel - Cai);
        const double w_inf = 1.0 - 1.0 / (1.0 + std::exp(-(V - 40.0) / 17.0));
        const double w_tau =
            6.0 * ((std::fabs(V - 7.9) < 1e-6)
          ? 2.0 / 13.0
          : (1.0 - std::exp(-(V - 7.9) / 5.0))
            / ((1.0 + 0.3 * std::exp(-(V - 7.9) / 5.0)) * (V - 7.9)));

        const double Fn =
            1e-12 * c_V_rel * Irel
          - 5e-13 / c_F * (0.5 * ICaL - 0.2 * INaCa) * c_Cm;

        const double u_inf = 1.0 / (1.0 + std::exp(-(Fn - c_c1) / c_c2));
        const double v_inf = 1.0 - 1.0 / (1.0 + std::exp(-(Fn - 0.2 * c_c1) / c_c2));
        const double v_tau = 1.91 + 2.09 / (1.0 + std::exp(-(Fn - c_c1) / c_c2));

        const double IbCa = c_gbCa * (V - ECa);
        const double IbNa = c_gbNa * (V - ENa);
        const double IK1 = c_gK1 * (V - EK) / (1.0 + std::exp(0.07 * (V + 80.0)));

        const double IKr = c_gKr * s_xr * (V - EK) / (1.0 + std::exp((V + 15.0) / 22.4));
        const double xr_inf = 1.0 / (1.0 + std::exp((V + 14.1) / -6.5));
        const double xr_alpha =
            0.0003 * ((std::fabs(V + 14.1) < 1e-6)
          ? 5.0 : (V + 14.1) / (1.0 - std::exp((V + 14.1) / -5.0)));
        const double xr_beta =
            7.3898e-05 * ((std::fabs(V - 3.3328) < 1e-7)
          ? 5.1237 : (V - 3.3328) / (std::exp((V - 3.3328) / 5.1237) - 1.0));
        const double xr_tau = 1.0 / (xr_alpha + xr_beta);

        const double IKs = c_gKs * s_xs * s_xs * (V - EK);
        const double xs_inf = 1.0 / std::sqrt(1.0 + std::exp((V - 19.9) / -12.7));
        const double xs_alpha =
            4e-05 * ((std::fabs(V - 19.9) < 1e-6)
          ? 17.0 : (V - 19.9) / (1.0 - std::exp((V - 19.9) / -17.0)));
        const double xs_beta =
            3.5e-05 * ((std::fabs(V - 19.9) < 1e-6)
          ? 9.0 : (V - 19.9) / (std::exp((V - 19.9) / 9.0) - 1.0));
        const double xs_tau = 0.5 / (xs_alpha + xs_beta);

        const double gKur = c_gKur_base * (1.0 + 10.0 / (1.0 + std::exp((V - 15.0) / -13.0)));
        const double ua_inf = 1.0 / (1.0 + std::exp((V + 30.3) / -9.6));
        const double ua_alpha = 0.65 / (std::exp((V + 10.0) / -8.5) + std::exp((V - 30.0) / -59.0));
        const double ua_beta = 0.65 / (2.5 + std::exp((V + 82.0) / 17.0));
        const double ua_tau = 1.0 / (ua_alpha + ua_beta) / c_KQ10;
        const double ui_inf = 1.0 / (1.0 + std::exp((V - 99.45) / 27.48));
        const double ui_alpha = 1.0 / (21.0 + std::exp((V - 185.0) / -28.0));
        const double ui_beta = 1.0 / std::exp((V - 158.0) / -16.0);
        const double ui_tau = 1.0 / (ui_alpha + ui_beta) / c_KQ10;
        const double IKur = gKur * s_ua * s_ua * s_ua * s_ui * (V - EK);

        const double INa = c_gNa * s_m * s_m * s_m * s_h * s_j * (V - ENa);
        const double h_alpha = (V < -40.0) ? 0.135 * std::exp((V + 80.0) / -6.8) : 0.0;
        const double h_beta = (V < -40.0)
          ? 3.56 * std::exp(0.079 * V) + 310000.0 * std::exp(0.35 * V)
          : 1.0 / (0.13 * (1.0 + std::exp((V + 10.66) / -11.1)));
        const double h_inf = h_alpha / (h_alpha + h_beta);
        const double h_tau = 1.0 / (h_alpha + h_beta);
        const double j_alpha = (V < -40.0)
          ? (-127140.0 * std::exp(0.2444 * V) - 3.474e-05 * std::exp(-0.04391 * V))
            * (V + 37.78) / (1.0 + std::exp(0.311 * (V + 79.23)))
          : 0.0;
        const double j_beta = (V < -40.0)
          ? 0.1212 * std::exp(-0.01052 * V) / (1.0 + std::exp(-0.1378 * (V + 40.14)))
          : 0.3 * std::exp(-2.535e-07 * V) / (1.0 + std::exp(-0.1 * (V + 32.0)));
        const double j_inf = j_alpha / (j_alpha + j_beta);
        const double j_tau = 1.0 / (j_alpha + j_beta);
        const double m_alpha = (std::fabs(V + 47.13) < 1e-12)
          ? 3.2 : 0.32 * (V + 47.13) / (1.0 - std::exp(-0.1 * (V + 47.13)));
        const double m_beta = 0.08 * std::exp(-V / 11.0);
        const double m_inf = m_alpha / (m_alpha + m_beta);
        const double m_tau = 1.0 / (m_alpha + m_beta);

        const double Ito = c_gto * s_oa * s_oa * s_oa * s_oi * (V - EK);
        const double oa_alpha = 0.65 / (std::exp((V + 10.0) / -8.5) + std::exp((V - 30.0) / -59.0));
        const double oa_beta = 0.65 / (2.5 + std::exp((V + 82.0) / 17.0));
        const double oa_inf = 1.0 / (1.0 + std::exp((V + 20.47) / -17.54));
        const double oa_tau = 1.0 / (oa_alpha + oa_beta) / c_KQ10;
        const double oi_alpha = 1.0 / (18.53 + std::exp((V + 113.7) / 10.95));
        const double oi_beta = 1.0 / (35.56 + std::exp((V + 1.26) / -7.44));
        const double oi_inf = 1.0 / (1.0 + std::exp((V + 43.1) / 5.3));
        const double oi_tau = 1.0 / (oi_alpha + oi_beta) / c_KQ10;

        const double CaRelKm = CaRel + c_Km_CSQN;
        const double B1 =
            (2.0 * INaCa - (IpCa + ICaL + IbCa)) * c_Cm / (2.0 * c_V_i * c_F)
          + (c_V_up * (IupLeak - Iup) + Irel * c_V_rel) / c_V_i;
        const double B2 =
            1.0
          + c_TRPN_max * c_Km_TRPN / std::pow(Cai + c_Km_TRPN, 2.0)
          + c_CMDN_max * c_Km_CMDN / std::pow(Cai + c_Km_CMDN, 2.0);

        const double Iion = INa + IK1 + Ito + IKur + IKr + IKs + ICaL
          + IpCa + INaK + INaCa + IbNa + IbCa;

        RATES_SoA[d_BASE     + cellI] = (d_inf - s_d) / d_tau;
        RATES_SoA[f_BASE     + cellI] = (f_inf - s_f) / f_tau;
        RATES_SoA[fCa_BASE   + cellI] = (fCa_inf - s_fCa) / c_ical_fCa_tau;
        RATES_SoA[w_BASE     + cellI] = (w_inf - s_w) / w_tau;
        RATES_SoA[u_BASE     + cellI] = (u_inf - s_u) / c_cajsr_u_tau;
        RATES_SoA[v_BASE     + cellI] = (v_inf - s_v) / v_tau;
        RATES_SoA[xr_BASE    + cellI] = (xr_inf - s_xr) / xr_tau;
        RATES_SoA[xs_BASE    + cellI] = (xs_inf - s_xs) / xs_tau;
        RATES_SoA[ua_BASE    + cellI] = (ua_inf - s_ua) / ua_tau;
        RATES_SoA[ui_BASE    + cellI] = (ui_inf - s_ui) / ui_tau;
        RATES_SoA[h_BASE     + cellI] = (h_inf - s_h) / h_tau;
        RATES_SoA[j_BASE     + cellI] = (j_inf - s_j) / j_tau;
        RATES_SoA[m_BASE     + cellI] = (m_inf - s_m) / m_tau;
        RATES_SoA[oa_BASE    + cellI] = (oa_inf - s_oa) / oa_tau;
        RATES_SoA[oi_BASE    + cellI] = (oi_inf - s_oi) / oi_tau;
        RATES_SoA[CaRel_BASE + cellI] =
            (Itr - Irel)
          / (1.0 + c_CSQN_max * c_Km_CSQN / (CaRelKm * CaRelKm));
        RATES_SoA[CaUp_BASE  + cellI] = Iup - (IupLeak + Itr * c_V_rel / c_V_up);
        RATES_SoA[Cai_BASE   + cellI] = B1 / B2;
        RATES_SoA[Ki_BASE    + cellI] =
            (2.0 * INaK - (IK1 + Ito + IKur + IKr + IKs + Istim))
          * c_Cm / (c_V_i * c_F);
        RATES_SoA[Nai_BASE   + cellI] =
            (-3.0 * INaK - (3.0 * INaCa + IbNa + INa)) * c_Cm / (c_V_i * c_F);

        if (solveVm)
        {
            RATES_SoA[V_BASE + cellI] = -Iion - Istim;
        }

        if (SUPPORT_SoA != nullptr)
        {
            SUPPORT_SoA[COURTEMANCHE_BATCH_SUPPORT_tau_m       * N + cellI] = m_tau;
            SUPPORT_SoA[COURTEMANCHE_BATCH_SUPPORT_tau_h       * N + cellI] = h_tau;
            SUPPORT_SoA[COURTEMANCHE_BATCH_SUPPORT_tau_j       * N + cellI] = j_tau;
            SUPPORT_SoA[COURTEMANCHE_BATCH_SUPPORT_tau_d       * N + cellI] = d_tau;
            SUPPORT_SoA[COURTEMANCHE_BATCH_SUPPORT_tau_f       * N + cellI] = f_tau;
            SUPPORT_SoA[COURTEMANCHE_BATCH_SUPPORT_tau_fCa     * N + cellI] = c_ical_fCa_tau;
            SUPPORT_SoA[COURTEMANCHE_BATCH_SUPPORT_tau_oa      * N + cellI] = oa_tau;
            SUPPORT_SoA[COURTEMANCHE_BATCH_SUPPORT_tau_oi      * N + cellI] = oi_tau;
            SUPPORT_SoA[COURTEMANCHE_BATCH_SUPPORT_tau_ua      * N + cellI] = ua_tau;
            SUPPORT_SoA[COURTEMANCHE_BATCH_SUPPORT_tau_ui      * N + cellI] = ui_tau;
            SUPPORT_SoA[COURTEMANCHE_BATCH_SUPPORT_tau_xr      * N + cellI] = xr_tau;
            SUPPORT_SoA[COURTEMANCHE_BATCH_SUPPORT_tau_xs      * N + cellI] = xs_tau;
            SUPPORT_SoA[COURTEMANCHE_BATCH_SUPPORT_tau_cajsr_u * N + cellI] = c_cajsr_u_tau;
            SUPPORT_SoA[COURTEMANCHE_BATCH_SUPPORT_tau_cajsr_v * N + cellI] = v_tau;
            SUPPORT_SoA[COURTEMANCHE_BATCH_SUPPORT_tau_cajsr_w * N + cellI] = w_tau;
            SUPPORT_SoA[COURTEMANCHE_BATCH_SUPPORT_gInf_m       * N + cellI] = m_inf;
            SUPPORT_SoA[COURTEMANCHE_BATCH_SUPPORT_gInf_h       * N + cellI] = h_inf;
            SUPPORT_SoA[COURTEMANCHE_BATCH_SUPPORT_gInf_j       * N + cellI] = j_inf;
            SUPPORT_SoA[COURTEMANCHE_BATCH_SUPPORT_gInf_d       * N + cellI] = d_inf;
            SUPPORT_SoA[COURTEMANCHE_BATCH_SUPPORT_gInf_f       * N + cellI] = f_inf;
            SUPPORT_SoA[COURTEMANCHE_BATCH_SUPPORT_gInf_fCa     * N + cellI] = fCa_inf;
            SUPPORT_SoA[COURTEMANCHE_BATCH_SUPPORT_gInf_oa      * N + cellI] = oa_inf;
            SUPPORT_SoA[COURTEMANCHE_BATCH_SUPPORT_gInf_oi      * N + cellI] = oi_inf;
            SUPPORT_SoA[COURTEMANCHE_BATCH_SUPPORT_gInf_ua      * N + cellI] = ua_inf;
            SUPPORT_SoA[COURTEMANCHE_BATCH_SUPPORT_gInf_ui      * N + cellI] = ui_inf;
            SUPPORT_SoA[COURTEMANCHE_BATCH_SUPPORT_gInf_xr      * N + cellI] = xr_inf;
            SUPPORT_SoA[COURTEMANCHE_BATCH_SUPPORT_gInf_xs      * N + cellI] = xs_inf;
            SUPPORT_SoA[COURTEMANCHE_BATCH_SUPPORT_gInf_cajsr_u * N + cellI] = u_inf;
            SUPPORT_SoA[COURTEMANCHE_BATCH_SUPPORT_gInf_cajsr_v * N + cellI] = v_inf;
            SUPPORT_SoA[COURTEMANCHE_BATCH_SUPPORT_gInf_cajsr_w * N + cellI] = w_inf;
            SUPPORT_SoA[COURTEMANCHE_BATCH_SUPPORT_Iion_cm     * N + cellI] = Iion;
        }
    }
}

} // End namespace Foam

// ************************************************************************* //
