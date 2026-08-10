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
    State-major (SoA) batched evaluator for the TNNP 2004 ionic model.
\*---------------------------------------------------------------------------*/

#pragma once

#include "TNNP_2004.H"
#include "TNNP_2004Names.H"
#include "StimulusProtocolPOD.H"

namespace Foam
{

enum TNNP_BATCH_SUPPORT_INDEX
{
    TNNP_BATCH_SUPPORT_tau_xr1,
    TNNP_BATCH_SUPPORT_gInf_xr1,
    TNNP_BATCH_SUPPORT_tau_xr2,
    TNNP_BATCH_SUPPORT_gInf_xr2,
    TNNP_BATCH_SUPPORT_tau_xs,
    TNNP_BATCH_SUPPORT_gInf_xs,
    TNNP_BATCH_SUPPORT_tau_m,
    TNNP_BATCH_SUPPORT_gInf_m,
    TNNP_BATCH_SUPPORT_tau_h,
    TNNP_BATCH_SUPPORT_gInf_h,
    TNNP_BATCH_SUPPORT_tau_j,
    TNNP_BATCH_SUPPORT_gInf_j,
    TNNP_BATCH_SUPPORT_tau_d,
    TNNP_BATCH_SUPPORT_gInf_d,
    TNNP_BATCH_SUPPORT_tau_f,
    TNNP_BATCH_SUPPORT_gInf_f,
    TNNP_BATCH_SUPPORT_tau_s,
    TNNP_BATCH_SUPPORT_gInf_s,
    TNNP_BATCH_SUPPORT_tau_r,
    TNNP_BATCH_SUPPORT_gInf_r,
    TNNP_BATCH_SUPPORT_Iion_cm,

    NUM_TNNP_BATCH_SUPPORT
};

static_assert(NUM_TNNP_BATCH_SUPPORT == 21, "TNNP compact support layout drift");


CARDIAC_HOST_DEVICE inline void TNNPComputeVariablesBatch
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
    const double c_R              = CONSTANTS[R];
    const double c_T              = CONSTANTS[::T];  // ::T = temperature enum, not Foam::T() transpose
    const double c_F              = CONSTANTS[F];
    const double c_Cm             = CONSTANTS[Cm];
    const double c_V_c            = CONSTANTS[V_c];
    const double c_P_kna          = CONSTANTS[P_kna];
    const double c_K_o            = CONSTANTS[K_o];
    const double c_Na_o           = CONSTANTS[Na_o];
    const double c_Ca_o           = CONSTANTS[Ca_o];
    const double c_g_K1           = CONSTANTS[g_K1];
    const double c_g_Kr           = CONSTANTS[g_Kr];
    const double c_g_Ks           = CONSTANTS[g_Ks];
    const double c_g_Na           = CONSTANTS[g_Na];
    const double c_g_bNa          = CONSTANTS[g_bNa];
    const double c_g_CaL          = CONSTANTS[g_CaL];
    const double c_g_bCa          = CONSTANTS[g_bCa];
    const double c_g_to           = CONSTANTS[g_to];
    const double c_P_NaK          = CONSTANTS[P_NaK];
    const double c_K_mk           = CONSTANTS[K_mk];
    const double c_K_mNa          = CONSTANTS[K_mNa];
    const double c_K_NaCa         = CONSTANTS[K_NaCa];
    const double c_K_sat          = CONSTANTS[K_sat];
    const double c_alpha_NCX      = CONSTANTS[alpha_NCX];
    const double c_gamma_NCX      = CONSTANTS[gamma_NCX];
    const double c_Km_Ca          = CONSTANTS[Km_Ca];
    const double c_Km_Nai         = CONSTANTS[Km_Nai];
    const double c_g_pCa          = CONSTANTS[g_pCa];
    const double c_K_pCa          = CONSTANTS[K_pCa];
    const double c_g_pK           = CONSTANTS[g_pK];
    const double c_tau_g          = CONSTANTS[tau_g];
    const double c_a_rel          = CONSTANTS[a_rel];
    const double c_b_rel          = CONSTANTS[b_rel];
    const double c_c_rel          = CONSTANTS[c_rel];
    const double c_K_up           = CONSTANTS[K_up];
    const double c_V_leak         = CONSTANTS[V_leak];
    const double c_Vmax_up        = CONSTANTS[Vmax_up];
    const double c_Buf_c          = CONSTANTS[Buf_c];
    const double c_K_buf_c        = CONSTANTS[K_buf_c];
    const double c_Buf_sr         = CONSTANTS[Buf_sr];
    const double c_K_buf_sr       = CONSTANTS[K_buf_sr];
    const double c_V_sr           = CONSTANTS[V_sr];
    const double c_tau_fCa        = CONSTANTS[tau_fCa];
    const double c_sInfPrefactor  = CONSTANTS[sInfPrefactor];
    const double c_sInfShift      = CONSTANTS[sInfShift];
    const double c_sInfScale      = CONSTANTS[sInfScale];
    const double c_tauSGaussAmp   = CONSTANTS[tauSGaussAmp];
    const double c_tauSGaussShift = CONSTANTS[tauSGaussShift];
    const double c_tauSGaussWidth = CONSTANTS[tauSGaussWidth];
    const double c_tauSSigmoidAmp   = CONSTANTS[tauSSigmoidAmp];
    const double c_tauSSigmoidShift = CONSTANTS[tauSSigmoidShift];
    const double c_tauSSigmoidScale = CONSTANTS[tauSSigmoidScale];
    const double c_tauSOffset     = CONSTANTS[tauSOffset];

    const double c_RTF    = (c_R * c_T) / c_F;
    const double c_FRT_V  = c_F / (c_R * c_T);  // F/(RT) — used as multiplier with V

    const int V_BASE     = V     * N;
    const int K_i_BASE   = K_i   * N;
    const int Na_i_BASE  = Na_i  * N;
    const int Ca_i_BASE  = Ca_i  * N;
    const int Xr1_BASE   = Xr1   * N;
    const int Xr2_BASE   = Xr2   * N;
    const int Xs_BASE    = Xs    * N;
    const int m_BASE     = m     * N;
    const int h_BASE     = h     * N;
    const int j_BASE     = j     * N;
    const int d_BASE     = d     * N;
    const int f_BASE     = f     * N;
    const int fCa_BASE   = fCa   * N;
    const int s_BASE     = s     * N;
    const int r_BASE     = r     * N;
    const int Ca_SR_BASE = Ca_SR * N;
    const int g_BASE     = g     * N;

#ifdef _OPENMP
    #pragma omp simd
#endif
    for (int cellI = beginCell; cellI < endCell; ++cellI)
    {
        const double s_V     = STATES_SoA[V_BASE     + cellI];
        const double s_K_i   = STATES_SoA[K_i_BASE   + cellI];
        const double s_Na_i  = STATES_SoA[Na_i_BASE  + cellI];
        const double s_Ca_i  = STATES_SoA[Ca_i_BASE  + cellI];
        const double s_Xr1   = STATES_SoA[Xr1_BASE   + cellI];
        const double s_Xr2   = STATES_SoA[Xr2_BASE   + cellI];
        const double s_Xs    = STATES_SoA[Xs_BASE     + cellI];
        const double s_m     = STATES_SoA[m_BASE      + cellI];
        const double s_h     = STATES_SoA[h_BASE      + cellI];
        const double s_j     = STATES_SoA[j_BASE      + cellI];
        const double s_d     = STATES_SoA[d_BASE      + cellI];
        const double s_f     = STATES_SoA[f_BASE      + cellI];
        const double s_fCa   = STATES_SoA[fCa_BASE    + cellI];
        const double s_s     = STATES_SoA[s_BASE      + cellI];
        const double s_r     = STATES_SoA[r_BASE      + cellI];
        const double s_Ca_SR = STATES_SoA[Ca_SR_BASE  + cellI];
        const double s_g     = STATES_SoA[g_BASE      + cellI];

        const double a_m_inf    = 1.00000 / std::pow(1.00000 + std::exp((-56.8600 - s_V) / 9.03000), 2.00000);
        const double a_alpha_m  = 1.00000 / (1.00000 + std::exp((-60.0000 - s_V) / 5.00000));
        const double a_beta_m   = 0.100000 / (1.00000 + std::exp((s_V + 35.0000) / 5.00000))
                                + 0.100000 / (1.00000 + std::exp((s_V - 50.0000) / 200.000));
        const double a_tau_m    = 1.00000 * a_alpha_m * a_beta_m;

        const double a_h_inf    = 1.00000 / std::pow(1.00000 + std::exp((s_V + 71.5500) / 7.43000), 2.00000);
        const double a_alpha_h  = (s_V < -40.0000
            ? 0.0570000 * std::exp(-(s_V + 80.0000) / 6.80000)
            : 0.00000);
        const double a_beta_h   = (s_V < -40.0000
            ? 2.70000 * std::exp(0.0790000 * s_V) + 310000.0 * std::exp(0.348500 * s_V)
            : 0.770000 / (0.130000 * (1.00000 + std::exp((s_V + 10.6600) / -11.1000))));
        const double a_tau_h    = 1.00000 / (a_alpha_h + a_beta_h);

        const double a_j_inf    = 1.00000 / std::pow(1.00000 + std::exp((s_V + 71.5500) / 7.43000), 2.00000);
        const double a_alpha_j  = (s_V < -40.0000
            ? (((-25428.0 * std::exp(0.244400 * s_V)
                 - 6.94800e-06 * std::exp(-0.0439100 * s_V))
                * (s_V + 37.7800)) / 1.00000)
               / (1.00000 + std::exp(0.311000 * (s_V + 79.2300)))
            : 0.00000);
        const double a_beta_j   = (s_V < -40.0000
            ? (0.0242400 * std::exp(-0.0105200 * s_V))
              / (1.00000 + std::exp(-0.137800 * (s_V + 40.1400)))
            : (0.600000 * std::exp(0.0570000 * s_V))
              / (1.00000 + std::exp(-0.100000 * (s_V + 32.0000))));
        const double a_tau_j    = 1.00000 / (a_alpha_j + a_beta_j);

        const double a_d_inf    = 1.00000 / (1.00000 + std::exp((-5.00000 - s_V) / 7.50000));
        const double a_alpha_d  = 1.40000 / (1.00000 + std::exp((-35.0000 - s_V) / 13.0000)) + 0.250000;
        const double a_gama_fCa = 1.40000 / (1.00000 + std::exp((s_V + 5.00000) / 5.00000));
        const double a_gamma_d  = 1.00000 / (1.00000 + std::exp((50.0000 - s_V) / 20.0000));
        const double a_tau_d    = 1.00000 * a_alpha_d * a_gama_fCa + a_gamma_d;

        const double a_f_inf    = 1.00000 / (1.00000 + std::exp((s_V + 20.0000) / 7.00000));
        const double a_tau_f    = 1125.00 * std::exp(-std::pow(s_V + 27.0000, 2.00000) / 240.000)
                                + 80.0000
                                + 165.000 / (1.00000 + std::exp((25.0000 - s_V) / 10.0000));

        const double a_alpha_fCa = 1.00000 / (1.00000 + std::pow(s_Ca_i / 0.000325000, 8.00000));
        const double a_beta_fCa  = 0.100000 / (1.00000 + std::exp((s_Ca_i - 0.000500000) / 0.000100000));
        const double a_beta_d2   = 0.200000 / (1.00000 + std::exp((s_Ca_i - 0.000750000) / 0.000800000));
        const double a_fCa_inf   = (a_alpha_fCa + a_beta_fCa + a_beta_d2 + 0.230000) / 1.46000;
        const double a_d_fCa     = (a_fCa_inf - s_fCa) / c_tau_fCa;
        const double rate_fCa    = (a_fCa_inf > s_fCa && s_V > -60.0000 ? 0.00000 : a_d_fCa);

        const double a_s_inf    = c_sInfPrefactor
                                / (1.00000 + std::exp((s_V + c_sInfShift) / c_sInfScale));
        const double a_tau_s    = c_tauSGaussAmp
                                  * std::exp(-std::pow(s_V + c_tauSGaussShift, 2.00000)
                                             / c_tauSGaussWidth)
                                + c_tauSSigmoidAmp
                                  / (1.00000 + std::exp((s_V + c_tauSSigmoidShift)
                                                        / c_tauSSigmoidScale))
                                + c_tauSOffset;

        const double a_r_inf    = 1.00000 / (1.00000 + std::exp((20.0000 - s_V) / 6.00000));
        const double a_tau_r    = 9.50000 * std::exp(-std::pow(s_V + 40.0000, 2.00000) / 1800.00)
                                + 0.800000;

        const double a_xr1_inf   = 1.00000 / (1.00000 + std::exp((-26.0000 - s_V) / 7.00000));
        const double a_alpha_xr1 = 450.000 / (1.00000 + std::exp((-45.0000 - s_V) / 10.0000));
        const double a_beta_xr1  = 6.00000 / (1.00000 + std::exp((s_V + 30.0000) / 11.5000));
        const double a_tau_xr1   = 1.00000 * a_alpha_xr1 * a_beta_xr1;

        const double a_xr2_inf   = 1.00000 / (1.00000 + std::exp((s_V + 88.0000) / 24.0000));
        const double a_alpha_xr2 = 3.00000 / (1.00000 + std::exp((-60.0000 - s_V) / 20.0000));
        const double a_beta_xr2  = 1.12000 / (1.00000 + std::exp((s_V - 60.0000) / 20.0000));
        const double a_tau_xr2   = 1.00000 * a_alpha_xr2 * a_beta_xr2;

        const double a_xs_inf   = 1.00000 / (1.00000 + std::exp((-5.00000 - s_V) / 14.0000));
        const double a_alpha_xs = 1100.00
                                / std::pow((1.00000 + std::exp((-10.0000 - s_V) / 6.00000)), 0.5);
        const double a_beta_xs  = 1.00000 / (1.00000 + std::exp((s_V - 60.0000) / 20.0000));
        const double a_tau_xs   = 1.00000 * a_alpha_xs * a_beta_xs;

        const double a_g_inf    = (s_Ca_i < 0.000350000
            ? 1.00000 / (1.00000 + std::pow(s_Ca_i / 0.000350000, 6.00000))
            : 1.00000 / (1.00000 + std::pow(s_Ca_i / 0.000350000, 16.0000)));
        const double a_d_g      = (a_g_inf - s_g) / c_tau_g;
        const double rate_g     = (a_g_inf > s_g && s_V > -60.0000 ? 0.00000 : a_d_g);

        const double a_i_NaK = ((c_P_NaK * c_K_o / (c_K_o + c_K_mk)) * s_Na_i
                               / (s_Na_i + c_K_mNa))
                             / (1.00000
                               + 0.124500 * std::exp(-0.100000 * s_V * c_FRT_V)
                               + 0.0353000 * std::exp(-s_V * c_FRT_V));

        const double a_E_Na  = c_RTF * std::log(c_Na_o / s_Na_i);
        const double a_i_Na  = c_g_Na * std::pow(s_m, 3.00000) * s_h * s_j * (s_V - a_E_Na);
        const double a_i_b_Na = c_g_bNa * (s_V - a_E_Na);

        const double a_i_NaCa = (c_K_NaCa
            * (std::exp(c_gamma_NCX * s_V * c_FRT_V) * std::pow(s_Na_i, 3.00000) * c_Ca_o
             - std::exp((c_gamma_NCX - 1.00000) * s_V * c_FRT_V)
               * std::pow(c_Na_o, 3.00000) * s_Ca_i * c_alpha_NCX))
          / ((std::pow(c_Km_Nai, 3.00000) + std::pow(c_Na_o, 3.00000))
             * (c_Km_Ca + c_Ca_o)
             * (1.00000 + c_K_sat * std::exp((c_gamma_NCX - 1.00000) * s_V * c_FRT_V)));

        const double a_E_K   = c_RTF * std::log(c_K_o / s_K_i);

        const double a_alpha_K1 = 0.100000
            / (1.00000 + std::exp(0.0600000 * ((s_V - a_E_K) - 200.000)));
        const double a_beta_K1  = (3.00000 * std::exp(0.000200000 * ((s_V - a_E_K) + 100.000))
                                  + 1.00000 * std::exp(0.100000 * ((s_V - a_E_K) - 10.0000)))
                                / (1.00000 + std::exp(-0.500000 * (s_V - a_E_K)));
        const double a_xK1_inf  = a_alpha_K1 / (a_alpha_K1 + a_beta_K1);
        const double a_i_K1     = c_g_K1 * a_xK1_inf
                                * std::pow(c_K_o / 5.40000, 0.5)
                                * (s_V - a_E_K);

        const double a_i_to  = c_g_to * s_r * s_s * (s_V - a_E_K);
        const double a_i_Kr  = c_g_Kr * std::pow(c_K_o / 5.40000, 0.5)
                             * s_Xr1 * s_Xr2 * (s_V - a_E_K);

        const double a_E_Ks  = c_RTF
                             * std::log((c_K_o + c_P_kna * c_Na_o)
                                      / (s_K_i + c_P_kna * s_Na_i));
        const double a_i_Ks  = c_g_Ks * std::pow(s_Xs, 2.00000) * (s_V - a_E_Ks);

        const double a_i_CaL = (c_g_CaL * s_d * s_f * s_fCa * 4.00000 * s_V
                               * std::pow(c_F, 2.00000) / (c_R * c_T))
                             * (s_Ca_i * std::exp(2.00000 * s_V * c_FRT_V)
                               - 0.341000 * c_Ca_o)
                             / (std::exp(2.00000 * s_V * c_FRT_V) - 1.00000);

        const double a_E_Ca  = 0.500000 * c_RTF * std::log(c_Ca_o / s_Ca_i);
        const double a_i_b_Ca = c_g_bCa * (s_V - a_E_Ca);
        const double a_i_p_K  = c_g_pK * (s_V - a_E_K)
                              / (1.00000 + std::exp((25.0000 - s_V) / 5.98000));
        const double a_i_p_Ca = c_g_pCa * s_Ca_i / (s_Ca_i + c_K_pCa);

        const double a_Iion_cm = a_i_K1 + a_i_to + a_i_Kr
                               + a_i_Ks + a_i_CaL + a_i_NaK
                               + a_i_Na + a_i_b_Na + a_i_NaCa
                               + a_i_b_Ca + a_i_p_K + a_i_p_Ca;

        const double a_Istim = computeStimulusPOD(t, stimulus);

        const double a_i_rel  = ((c_a_rel * std::pow(s_Ca_SR, 2.00000))
                               / (std::pow(c_b_rel, 2.00000) + std::pow(s_Ca_SR, 2.00000))
                               + c_c_rel) * s_d * s_g;
        const double a_i_up   = c_Vmax_up
                              / (1.00000 + std::pow(c_K_up, 2.00000) / std::pow(s_Ca_i, 2.00000));
        const double a_i_leak = c_V_leak * (s_Ca_SR - s_Ca_i);

        const double a_ddt_Ca_i_total =
            ((-((a_i_CaL + a_i_b_Ca + a_i_p_Ca) - 2.00000 * a_i_NaCa)
              / (2.00000 * c_V_c * c_F)) * c_Cm
             + a_i_leak) - a_i_up + a_i_rel;

        const double a_f_JCa_i_free = 1.00000
            / (1.00000 + (c_Buf_c * c_K_buf_c)
               / std::pow(s_Ca_i + c_K_buf_c, 2.00000));

        const double a_ddt_Ca_sr_total = (c_V_c / c_V_sr)
            * (a_i_up - (a_i_rel + a_i_leak));

        const double a_f_JCa_sr_free = 1.00000
            / (1.00000 + (c_Buf_sr * c_K_buf_sr)
               / std::pow(s_Ca_SR + c_K_buf_sr, 2.00000));

        RATES_SoA[m_BASE      + cellI] = (a_m_inf   - s_m)   / a_tau_m;
        RATES_SoA[h_BASE      + cellI] = (a_h_inf   - s_h)   / a_tau_h;
        RATES_SoA[j_BASE      + cellI] = (a_j_inf   - s_j)   / a_tau_j;
        RATES_SoA[d_BASE      + cellI] = (a_d_inf   - s_d)   / a_tau_d;
        RATES_SoA[f_BASE      + cellI] = (a_f_inf   - s_f)   / a_tau_f;
        RATES_SoA[fCa_BASE    + cellI] = rate_fCa;
        RATES_SoA[s_BASE      + cellI] = (a_s_inf   - s_s)   / a_tau_s;
        RATES_SoA[r_BASE      + cellI] = (a_r_inf   - s_r)   / a_tau_r;
        RATES_SoA[Xr1_BASE    + cellI] = (a_xr1_inf - s_Xr1) / a_tau_xr1;
        RATES_SoA[Xr2_BASE    + cellI] = (a_xr2_inf - s_Xr2) / a_tau_xr2;
        RATES_SoA[Xs_BASE     + cellI] = (a_xs_inf  - s_Xs)  / a_tau_xs;
        RATES_SoA[g_BASE      + cellI] = rate_g;

        RATES_SoA[Na_i_BASE   + cellI] = (-(a_i_Na + a_i_b_Na
                                            + 3.00000 * a_i_NaK
                                            + 3.00000 * a_i_NaCa) * c_Cm)
                                       / (c_V_c * c_F);

        RATES_SoA[K_i_BASE    + cellI] = (-((a_i_K1 + a_i_to + a_i_Kr + a_i_Ks
                                             + a_i_p_K + a_Istim)
                                           - 2.00000 * a_i_NaK) * c_Cm)
                                       / (c_V_c * c_F);

        RATES_SoA[Ca_i_BASE   + cellI] = a_ddt_Ca_i_total * a_f_JCa_i_free;
        RATES_SoA[Ca_SR_BASE  + cellI] = a_ddt_Ca_sr_total * a_f_JCa_sr_free;

        RATES_SoA[V_BASE + cellI] = 0.0;
        if (solveVmWithinODESolver)
        {
            RATES_SoA[V_BASE + cellI] = -a_Iion_cm - a_Istim;
        }

        if (SUPPORT_SoA != nullptr)
        {
            SUPPORT_SoA[TNNP_BATCH_SUPPORT_tau_xr1  * N + cellI] = a_tau_xr1;
            SUPPORT_SoA[TNNP_BATCH_SUPPORT_gInf_xr1 * N + cellI] = a_xr1_inf;
            SUPPORT_SoA[TNNP_BATCH_SUPPORT_tau_xr2  * N + cellI] = a_tau_xr2;
            SUPPORT_SoA[TNNP_BATCH_SUPPORT_gInf_xr2 * N + cellI] = a_xr2_inf;
            SUPPORT_SoA[TNNP_BATCH_SUPPORT_tau_xs   * N + cellI] = a_tau_xs;
            SUPPORT_SoA[TNNP_BATCH_SUPPORT_gInf_xs  * N + cellI] = a_xs_inf;
            SUPPORT_SoA[TNNP_BATCH_SUPPORT_tau_m    * N + cellI] = a_tau_m;
            SUPPORT_SoA[TNNP_BATCH_SUPPORT_gInf_m   * N + cellI] = a_m_inf;
            SUPPORT_SoA[TNNP_BATCH_SUPPORT_tau_h    * N + cellI] = a_tau_h;
            SUPPORT_SoA[TNNP_BATCH_SUPPORT_gInf_h   * N + cellI] = a_h_inf;
            SUPPORT_SoA[TNNP_BATCH_SUPPORT_tau_j    * N + cellI] = a_tau_j;
            SUPPORT_SoA[TNNP_BATCH_SUPPORT_gInf_j   * N + cellI] = a_j_inf;
            SUPPORT_SoA[TNNP_BATCH_SUPPORT_tau_d    * N + cellI] = a_tau_d;
            SUPPORT_SoA[TNNP_BATCH_SUPPORT_gInf_d   * N + cellI] = a_d_inf;
            SUPPORT_SoA[TNNP_BATCH_SUPPORT_tau_f    * N + cellI] = a_tau_f;
            SUPPORT_SoA[TNNP_BATCH_SUPPORT_gInf_f   * N + cellI] = a_f_inf;
            SUPPORT_SoA[TNNP_BATCH_SUPPORT_tau_s    * N + cellI] = a_tau_s;
            SUPPORT_SoA[TNNP_BATCH_SUPPORT_gInf_s   * N + cellI] = a_s_inf;
            SUPPORT_SoA[TNNP_BATCH_SUPPORT_tau_r    * N + cellI] = a_tau_r;
            SUPPORT_SoA[TNNP_BATCH_SUPPORT_gInf_r   * N + cellI] = a_r_inf;
            SUPPORT_SoA[TNNP_BATCH_SUPPORT_Iion_cm  * N + cellI] = a_Iion_cm;
        }
    }
}


} // End namespace Foam

// ************************************************************************* //
