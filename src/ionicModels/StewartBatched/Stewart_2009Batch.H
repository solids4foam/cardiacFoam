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
    State-major SoA batched evaluator for the Stewart 2009 ionic model.
    Generated mechanically from the scalar equation body; no Foam runtime
    objects are used inside the hot path.
\*---------------------------------------------------------------------------*/

#pragma once

#include "Stewart_2009Names.H"
#include "StimulusProtocolPOD.H"
#include <cmath>

namespace Foam
{

enum STEWART_BATCH_SUPPORT_INDEX
{
    STEWART_BATCH_SUPPORT_tau_y,
    STEWART_BATCH_SUPPORT_gInf_y,
    STEWART_BATCH_SUPPORT_tau_xr1,
    STEWART_BATCH_SUPPORT_gInf_xr1,
    STEWART_BATCH_SUPPORT_tau_xr2,
    STEWART_BATCH_SUPPORT_gInf_xr2,
    STEWART_BATCH_SUPPORT_tau_xs,
    STEWART_BATCH_SUPPORT_gInf_xs,
    STEWART_BATCH_SUPPORT_tau_m,
    STEWART_BATCH_SUPPORT_gInf_m,
    STEWART_BATCH_SUPPORT_tau_h,
    STEWART_BATCH_SUPPORT_gInf_h,
    STEWART_BATCH_SUPPORT_tau_j,
    STEWART_BATCH_SUPPORT_gInf_j,
    STEWART_BATCH_SUPPORT_tau_d,
    STEWART_BATCH_SUPPORT_gInf_d,
    STEWART_BATCH_SUPPORT_tau_f,
    STEWART_BATCH_SUPPORT_gInf_f,
    STEWART_BATCH_SUPPORT_tau_f2,
    STEWART_BATCH_SUPPORT_gInf_f2,
    STEWART_BATCH_SUPPORT_tau_fCass,
    STEWART_BATCH_SUPPORT_gInf_fCass,
    STEWART_BATCH_SUPPORT_tau_s,
    STEWART_BATCH_SUPPORT_gInf_s,
    STEWART_BATCH_SUPPORT_tau_r,
    STEWART_BATCH_SUPPORT_gInf_r,
    STEWART_BATCH_SUPPORT_Iion_cm,

    NUM_STEWART_BATCH_SUPPORT
};

static_assert(NUM_STEWART_BATCH_SUPPORT == 27, "Stewart compact support layout drift");

CARDIAC_HOST_DEVICE inline void StewartComputeVariablesBatch
(
    const double VOI,
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
    using std::exp;
    using std::fabs;
    using std::log;
    using std::pow;
    using std::sqrt;

#ifdef _OPENMP
    #pragma omp simd
#endif
    for (int cellI = beginCell; cellI < endCell; ++cellI)
    {
        const double state_Ca_SR = STATES_SoA[Ca_SR * N + cellI];
        const double state_Ca_i = STATES_SoA[Ca_i * N + cellI];
        const double state_Ca_ss = STATES_SoA[Ca_ss * N + cellI];
        const double state_IL_type_Ca_current_d_gate_d = STATES_SoA[IL_type_Ca_current_d_gate_d * N + cellI];
        const double state_IL_type_Ca_current_f2_gate_f2 = STATES_SoA[IL_type_Ca_current_f2_gate_f2 * N + cellI];
        const double state_IL_type_Ca_current_fCass_gate_fCass = STATES_SoA[IL_type_Ca_current_fCass_gate_fCass * N + cellI];
        const double state_IL_type_Ca_current_f_gate_f = STATES_SoA[IL_type_Ca_current_f_gate_f * N + cellI];
        const double state_Ifast_sodium_current_h_gate_h = STATES_SoA[Ifast_sodium_current_h_gate_h * N + cellI];
        const double state_Ifast_sodium_current_j_gate_j = STATES_SoA[Ifast_sodium_current_j_gate_j * N + cellI];
        const double state_Ifast_sodium_current_m_gate_m = STATES_SoA[Ifast_sodium_current_m_gate_m * N + cellI];
        const double state_Ihyperpolarization_activated_current_y_gate_y = STATES_SoA[Ihyperpolarization_activated_current_y_gate_y * N + cellI];
        const double state_Irapid_time_dependent_potassium_current_Xr1_gate_Xr1 = STATES_SoA[Irapid_time_dependent_potassium_current_Xr1_gate_Xr1 * N + cellI];
        const double state_Irapid_time_dependent_potassium_current_Xr2_gate_Xr2 = STATES_SoA[Irapid_time_dependent_potassium_current_Xr2_gate_Xr2 * N + cellI];
        const double state_Islow_time_dependent_potassium_current_Xs_gate_Xs = STATES_SoA[Islow_time_dependent_potassium_current_Xs_gate_Xs * N + cellI];
        const double state_Itransient_outward_current_r_gate_r = STATES_SoA[Itransient_outward_current_r_gate_r * N + cellI];
        const double state_Itransient_outward_current_s_gate_s = STATES_SoA[Itransient_outward_current_s_gate_s * N + cellI];
        const double state_K_i = STATES_SoA[K_i * N + cellI];
        const double state_Na_i = STATES_SoA[Na_i * N + cellI];
        const double state_R_prime = STATES_SoA[R_prime * N + cellI];
        const double state_membrane_V = STATES_SoA[membrane_V * N + cellI];

        double alg_AV_Ca_i_bufc = 0.0;
        double alg_AV_Ca_sr_bufsr = 0.0;
        double alg_AV_Ca_ss_bufss = 0.0;
        double alg_AV_E_Ca = 0.0;
        double alg_AV_E_K = 0.0;
        double alg_AV_E_Ks = 0.0;
        double alg_AV_E_Na = 0.0;
        double alg_AV_O = 0.0;
        double alg_AV_a = 0.0;
        double alg_AV_alpha_d = 0.0;
        double alg_AV_alpha_h = 0.0;
        double alg_AV_alpha_j = 0.0;
        double alg_AV_alpha_m = 0.0;
        double alg_AV_alpha_xr1 = 0.0;
        double alg_AV_alpha_xr2 = 0.0;
        double alg_AV_alpha_xs = 0.0;
        double alg_AV_alpha_y = 0.0;
        double alg_AV_beta_d = 0.0;
        double alg_AV_beta_h = 0.0;
        double alg_AV_beta_j = 0.0;
        double alg_AV_beta_m = 0.0;
        double alg_AV_beta_xr1 = 0.0;
        double alg_AV_beta_xr2 = 0.0;
        double alg_AV_beta_xs = 0.0;
        double alg_AV_beta_y = 0.0;
        double alg_AV_d_inf = 0.0;
        double alg_AV_f2_inf = 0.0;
        double alg_AV_fCass_inf = 0.0;
        double alg_AV_f_inf = 0.0;
        double alg_AV_gamma_d = 0.0;
        double alg_AV_h_inf = 0.0;
        double alg_AV_i_CaL = 0.0;
        double alg_AV_i_K1 = 0.0;
        double alg_AV_i_Kr = 0.0;
        double alg_AV_i_Ks = 0.0;
        double alg_AV_i_Na = 0.0;
        double alg_AV_i_NaCa = 0.0;
        double alg_AV_i_NaK = 0.0;
        double alg_AV_i_b_Ca = 0.0;
        double alg_AV_i_b_Na = 0.0;
        double alg_AV_i_f = 0.0;
        double alg_AV_i_f_K = 0.0;
        double alg_AV_i_f_Na = 0.0;
        double alg_AV_i_leak = 0.0;
        double alg_AV_i_p_Ca = 0.0;
        double alg_AV_i_p_K = 0.0;
        double alg_AV_i_rel = 0.0;
        double alg_AV_i_sus = 0.0;
        double alg_AV_i_to = 0.0;
        double alg_AV_i_up = 0.0;
        double alg_AV_i_xfer = 0.0;
        double alg_AV_j_inf = 0.0;
        double alg_AV_k1 = 0.0;
        double alg_AV_k2 = 0.0;
        double alg_AV_kcasr = 0.0;
        double alg_AV_m_inf = 0.0;
        double alg_AV_r_inf = 0.0;
        double alg_AV_s_inf = 0.0;
        double alg_AV_tau_d = 0.0;
        double alg_AV_tau_f = 0.0;
        double alg_AV_tau_f2 = 0.0;
        double alg_AV_tau_fCass = 0.0;
        double alg_AV_tau_h = 0.0;
        double alg_AV_tau_j = 0.0;
        double alg_AV_tau_m = 0.0;
        double alg_AV_tau_r = 0.0;
        double alg_AV_tau_s = 0.0;
        double alg_AV_tau_xr1 = 0.0;
        double alg_AV_tau_xr2 = 0.0;
        double alg_AV_tau_xs = 0.0;
        double alg_AV_tau_y = 0.0;
        double alg_AV_xK1_inf = 0.0;
        double alg_AV_xr1_inf = 0.0;
        double alg_AV_xr2_inf = 0.0;
        double alg_AV_xs_inf = 0.0;
        double alg_AV_y_inf = 0.0;
        double alg_Iion_cm = 0.0;
        double alg_Istim = 0.0;



        /* L_type_Ca_current_d_gate */
            alg_AV_alpha_d = 1.4 / (1.0 + exp((-35.0 - state_membrane_V) / 13.0)) + 0.25;
            alg_AV_beta_d = 1.4 / (1.0 + exp((state_membrane_V + 5.0) / 5.0));
            alg_AV_d_inf = 1.0 / (1.0 + exp((-8.0 - state_membrane_V) / 7.5));
            alg_AV_gamma_d = 1.0 / (1.0 + exp((50.0 - state_membrane_V) / 20.0));
            alg_AV_tau_d = 1.0 * alg_AV_alpha_d * alg_AV_beta_d + alg_AV_gamma_d;
            RATES_SoA[IL_type_Ca_current_d_gate_d * N + cellI] = (alg_AV_d_inf - state_IL_type_Ca_current_d_gate_d) / alg_AV_tau_d;

            /* L_type_Ca_current_f2_gate */
            alg_AV_f2_inf = 0.67 / (1.0 + exp((state_membrane_V + 35.0) / 7.0)) + 0.33;
            alg_AV_tau_f2 = 562.0 * exp(-pow(state_membrane_V + 27.0, 2.0) / 240.0) + 31.0 / (1.0 + exp((25.0 - state_membrane_V) / 10.0)) + 80.0 / (1.0 + exp((state_membrane_V + 30.0) / 10.0));
            RATES_SoA[IL_type_Ca_current_f2_gate_f2 * N + cellI] = (alg_AV_f2_inf - state_IL_type_Ca_current_f2_gate_f2) / alg_AV_tau_f2;

            /* L_type_Ca_current_fCass_gate */
            alg_AV_fCass_inf = 0.6 / (1.0 + pow(state_Ca_ss / 0.05, 2.0)) + 0.4;
            alg_AV_tau_fCass = 80.0 / (1.0 + pow(state_Ca_ss / 0.05, 2.0)) + 2.0;
            RATES_SoA[IL_type_Ca_current_fCass_gate_fCass * N + cellI] = (alg_AV_fCass_inf - state_IL_type_Ca_current_fCass_gate_fCass) / alg_AV_tau_fCass;

            /* L_type_Ca_current_f_gate */
            alg_AV_f_inf = 1.0 / (1.0 + exp((state_membrane_V + 20.0) / 7.0));
            alg_AV_tau_f = 1102.5 * exp(-pow(state_membrane_V + 27.0, 2.0) / 225.0) + 200.0 / (1.0 + exp((13.0 - state_membrane_V) / 10.0)) + 180.0 / (1.0 + exp((state_membrane_V + 30.0) / 10.0)) + 20.0;
            RATES_SoA[IL_type_Ca_current_f_gate_f * N + cellI] = (alg_AV_f_inf - state_IL_type_Ca_current_f_gate_f) / alg_AV_tau_f;

            /* calcium_pump_current */
            alg_AV_i_p_Ca = CONSTANTS[AC_g_pCa] * state_Ca_i / (state_Ca_i + CONSTANTS[AC_K_pCa]);


            /* fast_sodium_current_h_gate */
            alg_AV_alpha_h = ((state_membrane_V < -40.0) ? 0.057 * exp(-(state_membrane_V + 80.0) / 6.8) : 0.0);
            alg_AV_beta_h = ((state_membrane_V < -40.0) ? 2.7 * exp(0.079 * state_membrane_V) + 310000.0 * exp(0.3485 * state_membrane_V) : 0.77 / (0.13 * (1.0 + exp((state_membrane_V + 10.66) / -11.1))));
            alg_AV_h_inf = 1.0 / pow(1.0 + exp((state_membrane_V + 71.55) / 7.43), 2.0);
            alg_AV_tau_h = 1.0 / (alg_AV_alpha_h + alg_AV_beta_h);
            RATES_SoA[Ifast_sodium_current_h_gate_h * N + cellI] = (alg_AV_h_inf - state_Ifast_sodium_current_h_gate_h) / alg_AV_tau_h;

            /* fast_sodium_current_j_gate */
            alg_AV_alpha_j = ((state_membrane_V < -40.0) ? (-25428.0 * exp(0.2444 * state_membrane_V) - 6.948e-06 * exp(-0.04391 * state_membrane_V)) * (state_membrane_V + 37.78) / 1.0 / (1.0 + exp(0.311 * (state_membrane_V + 79.23))) : 0.0);
            alg_AV_beta_j = ((state_membrane_V < -40.0) ? 0.02424 * exp(-0.01052 * state_membrane_V) / (1.0 + exp(-0.1378 * (state_membrane_V + 40.14))) : 0.6 * exp(0.057 * state_membrane_V) / (1.0 + exp(-0.1 * (state_membrane_V + 32.0))));
            alg_AV_j_inf = 1.0 / pow(1.0 + exp((state_membrane_V + 71.55) / 7.43), 2.0);
            alg_AV_tau_j = 1.0 / (alg_AV_alpha_j + alg_AV_beta_j);
            RATES_SoA[Ifast_sodium_current_j_gate_j * N + cellI] = (alg_AV_j_inf - state_Ifast_sodium_current_j_gate_j) / alg_AV_tau_j;

            /* fast_sodium_current_m_gate */
            alg_AV_alpha_m = 1.0 / (1.0 + exp((-60.0 - state_membrane_V) / 5.0));
            alg_AV_beta_m = 0.1 / (1.0 + exp((state_membrane_V + 35.0) / 5.0)) + 0.1 / (1.0 + exp((state_membrane_V - 50.0) / 200.0));
            alg_AV_m_inf = 1.0 / pow(1.0 + exp((-56.86 - state_membrane_V) / 9.03), 2.0);
            alg_AV_tau_m = 1.0 * alg_AV_alpha_m * alg_AV_beta_m;
            RATES_SoA[Ifast_sodium_current_m_gate_m * N + cellI] = (alg_AV_m_inf - state_Ifast_sodium_current_m_gate_m) / alg_AV_tau_m;

            /* hyperpolarization_activated_current_y_gate */
            alg_AV_alpha_y = 1.0 * exp(-2.9 - 0.04 * state_membrane_V);
            alg_AV_beta_y = 1.0 * exp(3.6 + 0.11 * state_membrane_V);
            alg_AV_y_inf = 1.0 / (1.0 + exp((state_membrane_V + 80.6) / 6.8));
            alg_AV_tau_y = 4000.0 / (alg_AV_alpha_y + alg_AV_beta_y);
            RATES_SoA[Ihyperpolarization_activated_current_y_gate_y * N + cellI] = (alg_AV_y_inf - state_Ihyperpolarization_activated_current_y_gate_y) / alg_AV_tau_y;

            /* rapid_time_dependent_potassium_current_Xr1_gate */
            alg_AV_alpha_xr1 = 450.0 / (1.0 + exp((-45.0 - state_membrane_V) / 10.0));
            alg_AV_beta_xr1 = 6.0 / (1.0 + exp((state_membrane_V + 30.0) / 11.5));
            alg_AV_xr1_inf = 1.0 / (1.0 + exp((-26.0 - state_membrane_V) / 7.0));
            alg_AV_tau_xr1 = 1.0 * alg_AV_alpha_xr1 * alg_AV_beta_xr1;
            RATES_SoA[Irapid_time_dependent_potassium_current_Xr1_gate_Xr1 * N + cellI] = (alg_AV_xr1_inf - state_Irapid_time_dependent_potassium_current_Xr1_gate_Xr1) / alg_AV_tau_xr1;

            /* rapid_time_dependent_potassium_current_Xr2_gate */
            alg_AV_alpha_xr2 = 3.0 / (1.0 + exp((-60.0 - state_membrane_V) / 20.0));
            alg_AV_beta_xr2 = 1.12 / (1.0 + exp((state_membrane_V - 60.0) / 20.0));
            alg_AV_xr2_inf = 1.0 / (1.0 + exp((state_membrane_V + 88.0) / 24.0));
            alg_AV_tau_xr2 = 1.0 * alg_AV_alpha_xr2 * alg_AV_beta_xr2;
            RATES_SoA[Irapid_time_dependent_potassium_current_Xr2_gate_Xr2 * N + cellI] = (alg_AV_xr2_inf - state_Irapid_time_dependent_potassium_current_Xr2_gate_Xr2) / alg_AV_tau_xr2;

            /* slow_time_dependent_potassium_current_Xs_gate */
            alg_AV_alpha_xs = 1400.0 / sqrt(1.0 + exp((5.0 - state_membrane_V) / 6.0));
            alg_AV_beta_xs = 1.0 / (1.0 + exp((state_membrane_V - 35.0) / 15.0));
            alg_AV_xs_inf = 1.0 / (1.0 + exp((-5.0 - state_membrane_V) / 14.0));
            alg_AV_tau_xs = 1.0 * alg_AV_alpha_xs * alg_AV_beta_xs + 80.0;
            RATES_SoA[Islow_time_dependent_potassium_current_Xs_gate_Xs * N + cellI] = (alg_AV_xs_inf - state_Islow_time_dependent_potassium_current_Xs_gate_Xs) / alg_AV_tau_xs;

            /* transient_outward_current_r_gate */
            alg_AV_r_inf = 1.0 / (1.0 + exp((20.0 - state_membrane_V) / 13.0));
            alg_AV_tau_r = 10.45 * exp(-pow(state_membrane_V + 40.0, 2.0) / 1800.0) + 7.3;
            RATES_SoA[Itransient_outward_current_r_gate_r * N + cellI] = (alg_AV_r_inf - state_Itransient_outward_current_r_gate_r) / alg_AV_tau_r;

            /* transient_outward_current_s_gate */
            alg_AV_s_inf = 1.0 / (1.0 + exp((state_membrane_V + 27.0) / 13.0));
            alg_AV_tau_s = 85.0 * exp(-pow(state_membrane_V + 25.0, 2.0) / 320.0) + 5.0 / (1.0 + exp((state_membrane_V - 40.0) / 5.0)) + 42.0;
            RATES_SoA[Itransient_outward_current_s_gate_s * N + cellI] = (alg_AV_s_inf - state_Itransient_outward_current_s_gate_s) / alg_AV_tau_s;

            /* calcium_dynamics */
            alg_AV_Ca_i_bufc = 1.0 / (1.0 + CONSTANTS[AC_Buf_c] * CONSTANTS[AC_K_buf_c] / pow(state_Ca_i + CONSTANTS[AC_K_buf_c], 2.0));
            alg_AV_Ca_sr_bufsr = 1.0 / (1.0 + CONSTANTS[AC_Buf_sr] * CONSTANTS[AC_K_buf_sr] / pow(state_Ca_SR + CONSTANTS[AC_K_buf_sr], 2.0));
            alg_AV_Ca_ss_bufss = 1.0 / (1.0 + CONSTANTS[AC_Buf_ss] * CONSTANTS[AC_K_buf_ss] / pow(state_Ca_ss + CONSTANTS[AC_K_buf_ss], 2.0));
            alg_AV_i_leak = CONSTANTS[AC_V_leak] * (state_Ca_SR - state_Ca_i);
            alg_AV_i_up = CONSTANTS[AC_Vmax_up] / (1.0 + pow(CONSTANTS[AC_K_up], 2.0) / pow(state_Ca_i, 2.0));
            alg_AV_i_xfer = CONSTANTS[AC_V_xfer] * (state_Ca_ss - state_Ca_i);
            alg_AV_kcasr = CONSTANTS[AC_max_sr] - (CONSTANTS[AC_max_sr] - CONSTANTS[AC_min_sr]) / (1.0 + pow(CONSTANTS[AC_EC] / state_Ca_SR, 2.0));
            alg_AV_k1 = CONSTANTS[AC_k1_prime] / alg_AV_kcasr;
            alg_AV_k2 = CONSTANTS[AC_k2_prime] * alg_AV_kcasr;
            alg_AV_O = alg_AV_k1 * pow(state_Ca_ss, 2.0) * state_R_prime / (CONSTANTS[AC_k3] + alg_AV_k1 * pow(state_Ca_ss, 2.0));
            RATES_SoA[R_prime * N + cellI] = -alg_AV_k2 * state_Ca_ss * state_R_prime + CONSTANTS[AC_k4] * (1.0 - state_R_prime);
            alg_AV_i_rel = CONSTANTS[AC_V_rel] * alg_AV_O * (state_Ca_SR - state_Ca_ss);
            RATES_SoA[Ca_SR * N + cellI] = alg_AV_Ca_sr_bufsr * (alg_AV_i_up - (alg_AV_i_rel + alg_AV_i_leak));

            /* inward_rectifier_potassium_current */
            alg_AV_xK1_inf = 1.0 / (1.0 + exp(0.1 * (state_membrane_V + 75.44)));

            /* reversal_potentials */
            alg_AV_E_Ca = 0.5 * CONSTANTS[AC_R] * CONSTANTS[AC_T] / CONSTANTS[AC_F] * log(CONSTANTS[AC_Ca_o] / state_Ca_i);
            alg_AV_E_K = CONSTANTS[AC_R] * CONSTANTS[AC_T] / CONSTANTS[AC_F] * log(CONSTANTS[AC_K_o] / state_K_i);

            /* sodium_potassium_pump_current */
            alg_AV_i_NaK = CONSTANTS[AC_P_NaK] * CONSTANTS[AC_K_o] / (CONSTANTS[AC_K_o] + CONSTANTS[AC_K_mk]) * state_Na_i / (state_Na_i + CONSTANTS[AC_K_mNa]) / (1.0 + 0.1245 * exp(-0.1 * state_membrane_V * CONSTANTS[AC_F] / (CONSTANTS[AC_R] * CONSTANTS[AC_T])) + 0.0353 * exp(-state_membrane_V * CONSTANTS[AC_F] / (CONSTANTS[AC_R] * CONSTANTS[AC_T])));

            /* sustained_outward_current */
            alg_AV_a = 1.0 / (1.0 + exp((5.0 - state_membrane_V) / 17.0));
            alg_AV_i_sus = CONSTANTS[AC_g_sus] * alg_AV_a * (state_membrane_V - alg_AV_E_K);

            /* transient_outward_current */
            alg_AV_i_to = CONSTANTS[AC_g_to] * state_Itransient_outward_current_r_gate_r * state_Itransient_outward_current_s_gate_s * (state_membrane_V - alg_AV_E_K);

            /* *remaining* */
            alg_AV_i_CaL = CONSTANTS[AC_g_CaL] * state_IL_type_Ca_current_d_gate_d * state_IL_type_Ca_current_f_gate_f * state_IL_type_Ca_current_f2_gate_f2 * state_IL_type_Ca_current_fCass_gate_fCass * 4.0 * (state_membrane_V - 15.0) * pow(CONSTANTS[AC_F],
                                                                                                                              2.0) / (CONSTANTS[AC_R] * CONSTANTS[AC_T]) * (0.25 * state_Ca_ss * exp(2.0 * (state_membrane_V - 15.0) * CONSTANTS[AC_F] / (CONSTANTS[AC_R] * CONSTANTS[AC_T])) - CONSTANTS[AC_Ca_o]) / (exp(2.0 * (state_membrane_V - 15.0) * CONSTANTS[AC_F] / (CONSTANTS[AC_R] * CONSTANTS[AC_T])) - 1.0);
            alg_AV_i_b_Ca = CONSTANTS[AC_g_bca] * (state_membrane_V - alg_AV_E_Ca);
            alg_AV_i_f_K = state_Ihyperpolarization_activated_current_y_gate_y * CONSTANTS[AC_g_f_K] * (state_membrane_V - alg_AV_E_K);
            alg_AV_i_K1 = CONSTANTS[AC_g_K1] * alg_AV_xK1_inf * (state_membrane_V - 8.0 - alg_AV_E_K);
            alg_AV_i_p_K = CONSTANTS[AC_g_pK] * (state_membrane_V - alg_AV_E_K) / (1.0 + exp((25.0 - state_membrane_V) / 5.98));
            alg_AV_i_Kr = CONSTANTS[AC_g_Kr] * sqrt(CONSTANTS[AC_K_o] / 5.4) * state_Irapid_time_dependent_potassium_current_Xr1_gate_Xr1 * state_Irapid_time_dependent_potassium_current_Xr2_gate_Xr2 * (state_membrane_V - alg_AV_E_K);
            alg_AV_E_Ks = CONSTANTS[AC_R] * CONSTANTS[AC_T] / CONSTANTS[AC_F] * log((CONSTANTS[AC_K_o] + CONSTANTS[AC_P_kna] * CONSTANTS[AC_Na_o]) / (state_K_i + CONSTANTS[AC_P_kna] * state_Na_i));
            alg_AV_E_Na = CONSTANTS[AC_R] * CONSTANTS[AC_T] / CONSTANTS[AC_F] * log(CONSTANTS[AC_Na_o] / state_Na_i);
            alg_AV_i_NaCa = CONSTANTS[AC_K_NaCa] * (exp(CONSTANTS[AC_sodium_calcium_exchanger_current_gamma] * state_membrane_V * CONSTANTS[AC_F] / (CONSTANTS[AC_R] * CONSTANTS[AC_T])) * pow(state_Na_i, 3.0) * CONSTANTS[AC_Ca_o] - exp((CONSTANTS[AC_sodium_calcium_exchanger_current_gamma] - 1.0) * state_membrane_V * CONSTANTS[AC_F] / (CONSTANTS[AC_R] * CONSTANTS[AC_T])) * pow(CONSTANTS[AC_Na_o], 3.0) * state_Ca_i * CONSTANTS[AC_alpha]) / ((pow(CONSTANTS[AC_Km_Nai], 3.0) + pow(CONSTANTS[AC_Na_o], 3.0)) * (CONSTANTS[AC_Km_Ca] + CONSTANTS[AC_Ca_o]) * (1.0 + CONSTANTS[AC_K_sat] * exp((CONSTANTS[AC_sodium_calcium_exchanger_current_gamma] - 1.0) * state_membrane_V * CONSTANTS[AC_F] / (CONSTANTS[AC_R] * CONSTANTS[AC_T]))));
            RATES_SoA[Ca_i * N + cellI] = alg_AV_Ca_i_bufc * ((alg_AV_i_leak - alg_AV_i_up) * CONSTANTS[AC_V_sr] / CONSTANTS[AC_V_c] + alg_AV_i_xfer - 1.0 * (alg_AV_i_b_Ca + alg_AV_i_p_Ca - 2.0 * alg_AV_i_NaCa) * CONSTANTS[AC_Cm] / (2.0 * 1.0 * CONSTANTS[AC_V_c] * CONSTANTS[AC_F]));
            RATES_SoA[Ca_ss * N + cellI] = alg_AV_Ca_ss_bufss * (-1.0 * alg_AV_i_CaL * CONSTANTS[AC_Cm] / (2.0 * 1.0 * CONSTANTS[AC_V_ss] * CONSTANTS[AC_F]) + alg_AV_i_rel * CONSTANTS[AC_V_sr] / CONSTANTS[AC_V_ss] - alg_AV_i_xfer * CONSTANTS[AC_V_c] / CONSTANTS[AC_V_ss]);
            alg_AV_i_Na = CONSTANTS[AC_g_Na] * pow(state_Ifast_sodium_current_m_gate_m, 3.0) * state_Ifast_sodium_current_h_gate_h * state_Ifast_sodium_current_j_gate_j * (state_membrane_V - alg_AV_E_Na);
            alg_AV_i_f_Na = state_Ihyperpolarization_activated_current_y_gate_y * CONSTANTS[AC_g_f_Na] * (state_membrane_V - alg_AV_E_Na);
            alg_AV_i_Ks = CONSTANTS[AC_g_Ks] * pow(state_Islow_time_dependent_potassium_current_Xs_gate_Xs, 2.0) * (state_membrane_V - alg_AV_E_Ks);
            alg_AV_i_b_Na = CONSTANTS[AC_g_bna] * (state_membrane_V - alg_AV_E_Na);
            alg_AV_i_f = alg_AV_i_f_Na + alg_AV_i_f_K;
            RATES_SoA[K_i * N + cellI] = -1.0 * (alg_AV_i_K1 + alg_AV_i_to + alg_AV_i_f_K + alg_AV_i_sus + alg_AV_i_Kr + alg_AV_i_Ks + alg_AV_i_p_K - 2.0 * alg_AV_i_NaK) / (1.0 * CONSTANTS[AC_V_c] * CONSTANTS[AC_F]) * CONSTANTS[AC_Cm];
            RATES_SoA[Na_i * N + cellI] = -1.0 * (alg_AV_i_Na + alg_AV_i_b_Na + alg_AV_i_f_Na + 3.0 * alg_AV_i_NaK + 3.0 * alg_AV_i_NaCa) / (1.0 * CONSTANTS[AC_V_c] * CONSTANTS[AC_F]) * CONSTANTS[AC_Cm];
            alg_Iion_cm=(alg_AV_i_K1 + alg_AV_i_to + alg_AV_i_sus + alg_AV_i_Kr + alg_AV_i_Ks + alg_AV_i_CaL + alg_AV_i_NaK + alg_AV_i_Na + alg_AV_i_b_Na + alg_AV_i_NaCa + alg_AV_i_b_Ca + alg_AV_i_p_K + alg_AV_i_p_Ca + alg_AV_i_f);


            alg_Istim = computeStimulusPOD(VOI, stimulus);

            if (solveVm)
            {
                RATES_SoA[membrane_V * N + cellI] = -alg_Iion_cm - alg_Istim;
            }
            else
            {
                RATES_SoA[membrane_V * N + cellI] = 0.0;
            }


        if (SUPPORT_SoA != nullptr)
        {
            SUPPORT_SoA[STEWART_BATCH_SUPPORT_tau_y * N + cellI] = alg_AV_tau_y;
            SUPPORT_SoA[STEWART_BATCH_SUPPORT_tau_xr1 * N + cellI] = alg_AV_tau_xr1;
            SUPPORT_SoA[STEWART_BATCH_SUPPORT_tau_xr2 * N + cellI] = alg_AV_tau_xr2;
            SUPPORT_SoA[STEWART_BATCH_SUPPORT_tau_xs * N + cellI] = alg_AV_tau_xs;
            SUPPORT_SoA[STEWART_BATCH_SUPPORT_tau_m * N + cellI] = alg_AV_tau_m;
            SUPPORT_SoA[STEWART_BATCH_SUPPORT_tau_h * N + cellI] = alg_AV_tau_h;
            SUPPORT_SoA[STEWART_BATCH_SUPPORT_tau_j * N + cellI] = alg_AV_tau_j;
            SUPPORT_SoA[STEWART_BATCH_SUPPORT_tau_d * N + cellI] = alg_AV_tau_d;
            SUPPORT_SoA[STEWART_BATCH_SUPPORT_tau_f * N + cellI] = alg_AV_tau_f;
            SUPPORT_SoA[STEWART_BATCH_SUPPORT_tau_f2 * N + cellI] = alg_AV_tau_f2;
            SUPPORT_SoA[STEWART_BATCH_SUPPORT_tau_fCass * N + cellI] = alg_AV_tau_fCass;
            SUPPORT_SoA[STEWART_BATCH_SUPPORT_tau_s * N + cellI] = alg_AV_tau_s;
            SUPPORT_SoA[STEWART_BATCH_SUPPORT_tau_r * N + cellI] = alg_AV_tau_r;
            SUPPORT_SoA[STEWART_BATCH_SUPPORT_gInf_y * N + cellI] = alg_AV_y_inf;
            SUPPORT_SoA[STEWART_BATCH_SUPPORT_gInf_xr1 * N + cellI] = alg_AV_xr1_inf;
            SUPPORT_SoA[STEWART_BATCH_SUPPORT_gInf_xr2 * N + cellI] = alg_AV_xr2_inf;
            SUPPORT_SoA[STEWART_BATCH_SUPPORT_gInf_xs * N + cellI] = alg_AV_xs_inf;
            SUPPORT_SoA[STEWART_BATCH_SUPPORT_gInf_m * N + cellI] = alg_AV_m_inf;
            SUPPORT_SoA[STEWART_BATCH_SUPPORT_gInf_h * N + cellI] = alg_AV_h_inf;
            SUPPORT_SoA[STEWART_BATCH_SUPPORT_gInf_j * N + cellI] = alg_AV_j_inf;
            SUPPORT_SoA[STEWART_BATCH_SUPPORT_gInf_d * N + cellI] = alg_AV_d_inf;
            SUPPORT_SoA[STEWART_BATCH_SUPPORT_gInf_f * N + cellI] = alg_AV_f_inf;
            SUPPORT_SoA[STEWART_BATCH_SUPPORT_gInf_f2 * N + cellI] = alg_AV_f2_inf;
            SUPPORT_SoA[STEWART_BATCH_SUPPORT_gInf_fCass * N + cellI] = alg_AV_fCass_inf;
            SUPPORT_SoA[STEWART_BATCH_SUPPORT_gInf_s * N + cellI] = alg_AV_s_inf;
            SUPPORT_SoA[STEWART_BATCH_SUPPORT_gInf_r * N + cellI] = alg_AV_r_inf;
            SUPPORT_SoA[STEWART_BATCH_SUPPORT_Iion_cm * N + cellI] = alg_Iion_cm;
        }
    }
}

} // End namespace Foam

// ************************************************************************* //
