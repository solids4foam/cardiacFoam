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
    State-major (SoA) batched evaluator for the Grandi 2011 ionic model.
\*---------------------------------------------------------------------------*/

#pragma once

#include "Grandi_2011Names.H"
#include "StimulusProtocolPOD.H"
#include <cmath>

namespace Foam
{

enum GRANDI_BATCH_SUPPORT_INDEX
{
    GRANDI_BATCH_SUPPORT_tau_Ikr_xr,
    GRANDI_BATCH_SUPPORT_gInf_Ikr_xr,
    GRANDI_BATCH_SUPPORT_tau_Iks_xs,
    GRANDI_BATCH_SUPPORT_gInf_Iks_xs,
    GRANDI_BATCH_SUPPORT_tau_Ikur_ikur_r,
    GRANDI_BATCH_SUPPORT_gInf_Ikur_ikur_r,
    GRANDI_BATCH_SUPPORT_tau_Ikur_s,
    GRANDI_BATCH_SUPPORT_gInf_Ikur_s,
    GRANDI_BATCH_SUPPORT_tau_Ina_h,
    GRANDI_BATCH_SUPPORT_gInf_Ina_h,
    GRANDI_BATCH_SUPPORT_tau_Ina_j,
    GRANDI_BATCH_SUPPORT_gInf_Ina_j,
    GRANDI_BATCH_SUPPORT_tau_Ina_m,
    GRANDI_BATCH_SUPPORT_gInf_Ina_m,
    GRANDI_BATCH_SUPPORT_tau_Inal_hl,
    GRANDI_BATCH_SUPPORT_gInf_Inal_hl,
    GRANDI_BATCH_SUPPORT_tau_Inal_ml,
    GRANDI_BATCH_SUPPORT_gInf_Inal_ml,
    GRANDI_BATCH_SUPPORT_tau_Ical_d,
    GRANDI_BATCH_SUPPORT_gInf_Ical_d,
    GRANDI_BATCH_SUPPORT_tau_Ical_f,
    GRANDI_BATCH_SUPPORT_gInf_Ical_f,
    GRANDI_BATCH_SUPPORT_tau_Ical_fCaB_jn,
    GRANDI_BATCH_SUPPORT_gInf_Ical_fCaB_jn,
    GRANDI_BATCH_SUPPORT_tau_Ical_fCaB_sl,
    GRANDI_BATCH_SUPPORT_gInf_Ical_fCaB_sl,
    GRANDI_BATCH_SUPPORT_tau_Ito_x,
    GRANDI_BATCH_SUPPORT_gInf_Ito_x,
    GRANDI_BATCH_SUPPORT_tau_Ito_y,
    GRANDI_BATCH_SUPPORT_gInf_Ito_y,
    GRANDI_BATCH_SUPPORT_Iion_cm,

    NUM_GRANDI_BATCH_SUPPORT
};

static_assert(NUM_GRANDI_BATCH_SUPPORT == 31, "Grandi compact support layout drift");

CARDIAC_HOST_DEVICE inline void GrandiComputeVariablesBatch
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
        const double state_Ca_i = STATES_SoA[Ca_i * N + cellI];
        const double state_Ca_jn = STATES_SoA[Ca_jn * N + cellI];
        const double state_Ca_sl = STATES_SoA[Ca_sl * N + cellI];
        const double state_Ca_sr = STATES_SoA[Ca_sr * N + cellI];
        const double state_Csqn = STATES_SoA[Csqn * N + cellI];
        const double state_Ical_d = STATES_SoA[Ical_d * N + cellI];
        const double state_Ical_f = STATES_SoA[Ical_f * N + cellI];
        const double state_Ical_fCaB_jn = STATES_SoA[Ical_fCaB_jn * N + cellI];
        const double state_Ical_fCaB_sl = STATES_SoA[Ical_fCaB_sl * N + cellI];
        const double state_Ikr_xr = STATES_SoA[Ikr_xr * N + cellI];
        const double state_Iks_xs = STATES_SoA[Iks_xs * N + cellI];
        const double state_Ikur_ikur_r = STATES_SoA[Ikur_ikur_r * N + cellI];
        const double state_Ikur_s = STATES_SoA[Ikur_s * N + cellI];
        const double state_Ina_h = STATES_SoA[Ina_h * N + cellI];
        const double state_Ina_j = STATES_SoA[Ina_j * N + cellI];
        const double state_Ina_m = STATES_SoA[Ina_m * N + cellI];
        const double state_Inal_hl = STATES_SoA[Inal_hl * N + cellI];
        const double state_Inal_ml = STATES_SoA[Inal_ml * N + cellI];
        const double state_Ito_x = STATES_SoA[Ito_x * N + cellI];
        const double state_Ito_y = STATES_SoA[Ito_y * N + cellI];
        const double state_K_i = STATES_SoA[K_i * N + cellI];
        const double state_Na_i = STATES_SoA[Na_i * N + cellI];
        const double state_Na_jn = STATES_SoA[Na_jn * N + cellI];
        const double state_Na_sl = STATES_SoA[Na_sl * N + cellI];
        const double state_buffCa_CaM = STATES_SoA[buffCa_CaM * N + cellI];
        const double state_buffCa_Myoc = STATES_SoA[buffCa_Myoc * N + cellI];
        const double state_buffCa_Myom = STATES_SoA[buffCa_Myom * N + cellI];
        const double state_buffCa_SLH_jn = STATES_SoA[buffCa_SLH_jn * N + cellI];
        const double state_buffCa_SLH_sl = STATES_SoA[buffCa_SLH_sl * N + cellI];
        const double state_buffCa_SLL_jn = STATES_SoA[buffCa_SLL_jn * N + cellI];
        const double state_buffCa_SLL_sl = STATES_SoA[buffCa_SLL_sl * N + cellI];
        const double state_buffCa_SRB = STATES_SoA[buffCa_SRB * N + cellI];
        const double state_buffCa_TnCHc = STATES_SoA[buffCa_TnCHc * N + cellI];
        const double state_buffCa_TnCHm = STATES_SoA[buffCa_TnCHm * N + cellI];
        const double state_buffCa_TnCL = STATES_SoA[buffCa_TnCL * N + cellI];
        const double state_buffNa_NaB_jn = STATES_SoA[buffNa_NaB_jn * N + cellI];
        const double state_buffNa_NaB_sl = STATES_SoA[buffNa_NaB_sl * N + cellI];
        const double state_membrane_V = STATES_SoA[membrane_V * N + cellI];
        const double state_ryr_i = STATES_SoA[ryr_i * N + cellI];
        const double state_ryr_o = STATES_SoA[ryr_o * N + cellI];
        const double state_ryr_ryr_r = STATES_SoA[ryr_ryr_r * N + cellI];

        double alg_AV_ECa_jn = 0.0;
        double alg_AV_ECa_sl = 0.0;
        double alg_AV_EK = 0.0;
        double alg_AV_EKs = 0.0;
        double alg_AV_ENa_jn = 0.0;
        double alg_AV_ENa_sl = 0.0;
        double alg_AV_ICaB_jn = 0.0;
        double alg_AV_ICaB_sl = 0.0;
        double alg_AV_ICaL_Ca_jn = 0.0;
        double alg_AV_ICaL_Ca_sl = 0.0;
        double alg_AV_ICaL_K = 0.0;
        double alg_AV_ICaL_Na_jn = 0.0;
        double alg_AV_ICaL_Na_sl = 0.0;
        double alg_AV_ICa_tot = 0.0;
        double alg_AV_ICa_tot_jn = 0.0;
        double alg_AV_ICa_tot_sl = 0.0;
        double alg_AV_IClB = 0.0;
        double alg_AV_IClCa = 0.0;
        double alg_AV_IClCa_jn = 0.0;
        double alg_AV_IClCa_sl = 0.0;
        double alg_AV_ICl_tot = 0.0;
        double alg_AV_IK1 = 0.0;
        double alg_AV_IK_tot = 0.0;
        double alg_AV_IKp = 0.0;
        double alg_AV_IKp_jn = 0.0;
        double alg_AV_IKp_sl = 0.0;
        double alg_AV_IKr = 0.0;
        double alg_AV_IKs = 0.0;
        double alg_AV_IKs_jn = 0.0;
        double alg_AV_IKs_sl = 0.0;
        double alg_AV_IKur = 0.0;
        double alg_AV_INaB_jn = 0.0;
        double alg_AV_INaB_sl = 0.0;
        double alg_AV_INaCa_jn = 0.0;
        double alg_AV_INaCa_sl = 0.0;
        double alg_AV_INaK = 0.0;
        double alg_AV_INaK_jn = 0.0;
        double alg_AV_INaK_sl = 0.0;
        double alg_AV_INaL_jn = 0.0;
        double alg_AV_INaL_sl = 0.0;
        double alg_AV_INa_jn = 0.0;
        double alg_AV_INa_sl = 0.0;
        double alg_AV_INa_tot = 0.0;
        double alg_AV_INa_tot_jn = 0.0;
        double alg_AV_INa_tot_sl = 0.0;
        double alg_AV_IpCa_jn = 0.0;
        double alg_AV_IpCa_sl = 0.0;
        double alg_AV_Ito = 0.0;
        double alg_AV_J_CaB_cytosol = 0.0;
        double alg_AV_J_CaB_jntion = 0.0;
        double alg_AV_J_CaB_sl = 0.0;
        double alg_AV_J_SRCarel = 0.0;
        double alg_AV_J_SR_leak = 0.0;
        double alg_AV_J_serca = 0.0;
        double alg_AV_Ka_jn = 0.0;
        double alg_AV_Ka_sl = 0.0;
        double alg_AV_RI = 0.0;
        double alg_AV_fnak = 0.0;
        double alg_AV_ibarca_jn = 0.0;
        double alg_AV_ibarca_sl = 0.0;
        double alg_AV_ibark = 0.0;
        double alg_AV_ibarna_jn = 0.0;
        double alg_AV_ibarna_sl = 0.0;
        double alg_AV_ical_d_inf = 0.0;
        double alg_AV_ical_d_tau = 0.0;
        double alg_AV_ical_f_inf = 0.0;
        double alg_AV_ical_f_tau = 0.0;
        double alg_AV_ik1_inf = 0.0;
        double alg_AV_ik1_inf_a = 0.0;
        double alg_AV_ik1_inf_b = 0.0;
        double alg_AV_ikr_xr_inf = 0.0;
        double alg_AV_ikr_xr_tau = 0.0;
        double alg_AV_iks_xs_inf = 0.0;
        double alg_AV_iks_xs_tau = 0.0;
        double alg_AV_ikur_r_inf = 0.0;
        double alg_AV_ikur_r_tau = 0.0;
        double alg_AV_ikur_s_inf = 0.0;
        double alg_AV_ikur_s_tau = 0.0;
        double alg_AV_ina_h_a = 0.0;
        double alg_AV_ina_h_b = 0.0;
        double alg_AV_ina_h_inf = 0.0;
        double alg_AV_ina_h_tau = 0.0;
        double alg_AV_ina_j_a = 0.0;
        double alg_AV_ina_j_b = 0.0;
        double alg_AV_ina_j_inf = 0.0;
        double alg_AV_ina_j_tau = 0.0;
        double alg_AV_ina_m_inf = 0.0;
        double alg_AV_ina_m_tau = 0.0;
        double alg_AV_inal_hl_inf = 0.0;
        double alg_AV_inal_ml_a = 0.0;
        double alg_AV_inal_ml_b = 0.0;
        double alg_AV_ipca_IpCa_jn_a = 0.0;
        double alg_AV_ipca_IpCa_sl_a = 0.0;
        double alg_AV_ito_x_inf = 0.0;
        double alg_AV_ito_x_tau = 0.0;
        double alg_AV_ito_y_inf = 0.0;
        double alg_AV_ito_y_tau = 0.0;
        double alg_AV_kCaSR = 0.0;
        double alg_AV_kiSRCa = 0.0;
        double alg_AV_koSRCa = 0.0;
        double alg_AV_kp = 0.0;
        double alg_AV_rr = 0.0;
        double alg_AV_s1_jn = 0.0;
        double alg_AV_s1_sl = 0.0;
        double alg_AV_s2_jn = 0.0;
        double alg_AV_s2_sl = 0.0;
        double alg_AV_s3_jn = 0.0;
        double alg_AV_s3_sl = 0.0;
        double alg_Iion_cm = 0.0;
        double alg_Istim = 0.0;

        double rate_Ca_i = 0.0;
        double rate_Ca_jn = 0.0;
        double rate_Ca_sl = 0.0;
        double rate_Ca_sr = 0.0;
        double rate_Csqn = 0.0;
        double rate_Ical_d = 0.0;
        double rate_Ical_f = 0.0;
        double rate_Ical_fCaB_jn = 0.0;
        double rate_Ical_fCaB_sl = 0.0;
        double rate_Ikr_xr = 0.0;
        double rate_Iks_xs = 0.0;
        double rate_Ikur_ikur_r = 0.0;
        double rate_Ikur_s = 0.0;
        double rate_Ina_h = 0.0;
        double rate_Ina_j = 0.0;
        double rate_Ina_m = 0.0;
        double rate_Inal_hl = 0.0;
        double rate_Inal_ml = 0.0;
        double rate_Ito_x = 0.0;
        double rate_Ito_y = 0.0;
        double rate_K_i = 0.0;
        double rate_Na_i = 0.0;
        double rate_Na_jn = 0.0;
        double rate_Na_sl = 0.0;
        double rate_buffCa_CaM = 0.0;
        double rate_buffCa_Myoc = 0.0;
        double rate_buffCa_Myom = 0.0;
        double rate_buffCa_SLH_jn = 0.0;
        double rate_buffCa_SLH_sl = 0.0;
        double rate_buffCa_SLL_jn = 0.0;
        double rate_buffCa_SLL_sl = 0.0;
        double rate_buffCa_SRB = 0.0;
        double rate_buffCa_TnCHc = 0.0;
        double rate_buffCa_TnCHm = 0.0;
        double rate_buffCa_TnCL = 0.0;
        double rate_buffNa_NaB_jn = 0.0;
        double rate_buffNa_NaB_sl = 0.0;
        double rate_membrane_V = 0.0;
        double rate_ryr_i = 0.0;
        double rate_ryr_o = 0.0;
        double rate_ryr_ryr_r = 0.0;


/* buffna */
    rate_buffNa_NaB_jn = CONSTANTS[AC_kon_na] * state_Na_jn * (CONSTANTS[AC_Bmax_Na_jn] - state_buffNa_NaB_jn) - CONSTANTS[AC_koff_na] * state_buffNa_NaB_jn;
    rate_buffNa_NaB_sl = CONSTANTS[AC_kon_na] * state_Na_sl * (CONSTANTS[AC_Bmax_Na_sl] - state_buffNa_NaB_sl) - CONSTANTS[AC_koff_na] * state_buffNa_NaB_sl;

    /* ikr_xr */
    alg_AV_ikr_xr_inf = 1.0 / (1.0 + exp(-(state_membrane_V + 10.0) / 5.0));
    alg_AV_ikr_xr_tau = +550.0 / (1.0 + exp((-22.0 - state_membrane_V) / 9.0)) * 6.0 / (1.0 + exp((state_membrane_V + 11.0) / 9.0)) + 230.0 / (1.0 + exp((state_membrane_V + 40.0) / 20.0));

    /* ikur_ikur_r */
    alg_AV_ikur_r_inf = 1.0 / (1.0 + exp((state_membrane_V + 6.0) / -8.6));
    alg_AV_ikur_r_tau = 9.0 / (1.0 + exp((state_membrane_V + 5.0) / 12.0)) + 0.5;

    /* ikur_s */
    alg_AV_ikur_s_inf = 1.0 / (1.0 + exp((state_membrane_V + 7.5) / 10.0));
    alg_AV_ikur_s_tau = 590.0 / (1.0 + exp((state_membrane_V + 60.0) / 10.0)) + 3050.0;

    /* ina_h */
    alg_AV_ina_h_a = ((state_membrane_V >= -40.0) ? 0.0 : 0.057 * exp(-(state_membrane_V + 80.0) / 6.8)) * 1.0;
    alg_AV_ina_h_b = ((state_membrane_V >= -40.0) ? 0.77 / (0.13 * (1.0 + exp(-(state_membrane_V + 10.66) / 11.1))) : 2.7 * exp(0.079 * state_membrane_V) + 3.1 * pow(10.0, 5.0) * exp(0.3485 * state_membrane_V)) * 1.0;
    alg_AV_ina_h_inf = 1.0 / pow(1.0 + exp((state_membrane_V + 71.55) / 7.43), 2.0);
    alg_AV_ina_h_tau = 1.0 / (alg_AV_ina_h_a + alg_AV_ina_h_b);

    /* ina_j */
    alg_AV_ina_j_a = ((state_membrane_V >= -40.0) ? 0.0 : (-2.5428 * pow(10.0, 4.0) * exp(0.2444 * state_membrane_V) - 6.948 * pow(10.0, -6.0) * exp(-0.04391 * state_membrane_V)) * (state_membrane_V + 37.78) / (1.0 + exp(0.311 * (state_membrane_V + 79.23))) * 1.0) * 1.0;
    alg_AV_ina_j_b = ((state_membrane_V >= -40.0) ? 0.6 * exp(0.057 * state_membrane_V) / (1.0 + exp(-0.1 * (state_membrane_V + 32.0))) : 0.02424 * exp(-0.01052 * state_membrane_V) / (1.0 + exp(-0.1378 * (state_membrane_V + 40.14)))) * 1.0;
    alg_AV_ina_j_inf = 1.0 / pow(1.0 + exp((state_membrane_V + 71.55) / 7.43), 2.0);
    alg_AV_ina_j_tau = 1.0 / (alg_AV_ina_j_a + alg_AV_ina_j_b);

    /* ina_m */
    alg_AV_ina_m_inf = 1.0 / pow(1.0 + exp(-(56.86 + state_membrane_V) / 9.03), 2.0);
    alg_AV_ina_m_tau = 0.1292 * exp(-pow((state_membrane_V + 45.79) / 15.54, 2.0)) + 0.06487 * exp(-pow((state_membrane_V - 4.823) / 51.12, 2.0));

    /* inal_hl */
    alg_AV_inal_hl_inf = 1.0 / (1.0 + exp((state_membrane_V + 91.0) / 6.1));

    /* inal_ml */
    alg_AV_inal_ml_a = 0.32 * (state_membrane_V + 47.13) / (1.0 - exp(-0.1 * (state_membrane_V + 47.13)));
    alg_AV_inal_ml_b = 0.08 * exp(-state_membrane_V / 11.0);

    /* ito_x */
    alg_AV_ito_x_inf = 1.0 / (1.0 + exp(-(state_membrane_V + 1.0) / 11.0));
    alg_AV_ito_x_tau = 3.5 * exp(-pow(state_membrane_V / 30.0, 2.0)) + 1.5;

    /* ito_y */
    alg_AV_ito_y_inf = 1.0 / (1.0 + exp((state_membrane_V + 40.5) / 11.5));
    alg_AV_ito_y_tau = 25.635 * exp(-pow((state_membrane_V + 52.45) / 15.8827, 2.0)) + 24.14;

    /* buffca */
    rate_buffCa_CaM = CONSTANTS[AC_kon_cam] * state_Ca_i * (CONSTANTS[AC_Bmax_CaM] - state_buffCa_CaM) - CONSTANTS[AC_koff_cam] * state_buffCa_CaM;
    rate_buffCa_Myoc = CONSTANTS[AC_kon_myoca] * state_Ca_i * (CONSTANTS[AC_Bmax_myosin] - state_buffCa_Myoc - state_buffCa_Myom) - CONSTANTS[AC_koff_myoca] * state_buffCa_Myoc;
    rate_buffCa_Myom = CONSTANTS[AC_kon_myomg] * CONSTANTS[AC_Mg_i] * (CONSTANTS[AC_Bmax_myosin] - state_buffCa_Myoc - state_buffCa_Myom) - CONSTANTS[AC_koff_myomg] * state_buffCa_Myom;
    rate_buffCa_SLH_jn = CONSTANTS[AC_kon_slh] * state_Ca_jn * (CONSTANTS[AC_Bmax_SLhighj] - state_buffCa_SLH_jn) - CONSTANTS[AC_koff_slh] * state_buffCa_SLH_jn;
    rate_buffCa_SLH_sl = CONSTANTS[AC_kon_slh] * state_Ca_sl * (CONSTANTS[AC_Bmax_SLhighsl] - state_buffCa_SLH_sl) - CONSTANTS[AC_koff_slh] * state_buffCa_SLH_sl;
    rate_buffCa_SLL_jn = CONSTANTS[AC_kon_sll] * state_Ca_jn * (CONSTANTS[AC_Bmax_SLlowj] - state_buffCa_SLL_jn) - CONSTANTS[AC_koff_sll] * state_buffCa_SLL_jn;
    rate_buffCa_SLL_sl = CONSTANTS[AC_kon_sll] * state_Ca_sl * (CONSTANTS[AC_Bmax_SLlowsl] - state_buffCa_SLL_sl) - CONSTANTS[AC_koff_sll] * state_buffCa_SLL_sl;
    rate_buffCa_SRB = CONSTANTS[AC_kon_sr] * state_Ca_i * (CONSTANTS[AC_Bmax_SR] - state_buffCa_SRB) - CONSTANTS[AC_koff_sr] * state_buffCa_SRB;
    rate_buffCa_TnCHc = CONSTANTS[AC_kon_tnchca] * state_Ca_i * (CONSTANTS[AC_Bmax_TnChigh] - state_buffCa_TnCHc - state_buffCa_TnCHm) - CONSTANTS[AC_koff_tnchca] * state_buffCa_TnCHc;
    rate_buffCa_TnCHm = CONSTANTS[AC_kon_tnchmg] * CONSTANTS[AC_Mg_i] * (CONSTANTS[AC_Bmax_TnChigh] - state_buffCa_TnCHc - state_buffCa_TnCHm) - CONSTANTS[AC_koff_tnchmg] * state_buffCa_TnCHm;
    rate_buffCa_TnCL = CONSTANTS[AC_kon_tncl] * state_Ca_i * (CONSTANTS[AC_Bmax_TnClow] - state_buffCa_TnCL) - CONSTANTS[AC_koff_tncl] * state_buffCa_TnCL;
    alg_AV_J_CaB_cytosol = rate_buffCa_TnCL + rate_buffCa_TnCHc + rate_buffCa_TnCHm + rate_buffCa_CaM + rate_buffCa_Myoc + rate_buffCa_Myom + rate_buffCa_SRB;
    alg_AV_J_CaB_jntion = rate_buffCa_SLL_jn + rate_buffCa_SLH_jn;
    alg_AV_J_CaB_sl = rate_buffCa_SLL_sl + rate_buffCa_SLH_sl;

    /* ical_d */
    alg_AV_ical_d_inf = 1.0 / (1.0 + exp(-(state_membrane_V + 3.0 * CONSTANTS[AC_ISO] + 9.0) / 6.0));
    alg_AV_ical_d_tau = alg_AV_ical_d_inf * (1.0 - exp(-(state_membrane_V + 3.0 * CONSTANTS[AC_ISO] + 9.0) / 6.0)) / (0.035 * (state_membrane_V + 3.0 * CONSTANTS[AC_ISO] + 9.0));

    /* ical_f */
    alg_AV_ical_f_inf = 1.0 / (1.0 + exp((state_membrane_V + 3.0 * CONSTANTS[AC_ISO] + 30.0) / 7.0)) + 0.2 / (1.0 + exp((50.0 - state_membrane_V - 3.0 * CONSTANTS[AC_ISO]) / 20.0));
    alg_AV_ical_f_tau = 1.0 / (0.0197 * exp(-pow(0.0337 * (state_membrane_V + 3.0 * CONSTANTS[AC_ISO] + 25.0), 2.0)) + 0.02);

    /* iks_xs */
    alg_AV_iks_xs_inf = 1.0 / (1.0 + exp(-(state_membrane_V + 40.0 * CONSTANTS[AC_ISO] + 3.8) / 14.25));
    alg_AV_iks_xs_tau = 990.1 / (1.0 + exp(-(state_membrane_V + 40.0 * CONSTANTS[AC_ISO] + 2.436) / 14.12));

    /* inak */
    alg_AV_fnak = 1.0 / (1.0 + 0.1245 * exp(-0.1 * state_membrane_V * CONSTANTS[AC_FRT]) + 0.0365 * CONSTANTS[AC_sigma] * exp(-state_membrane_V * CONSTANTS[AC_FRT]));
    alg_AV_INaK_jn = CONSTANTS[AC_Fjn] * CONSTANTS[AC_IbarNaK] * alg_AV_fnak * CONSTANTS[AC_K_o] / (1.0 + pow(CONSTANTS[AC_KmNaip] / state_Na_jn, 4.0)) / (CONSTANTS[AC_K_o] + CONSTANTS[AC_KmKo]);
    alg_AV_INaK_sl = CONSTANTS[AC_Fsl] * CONSTANTS[AC_IbarNaK] * alg_AV_fnak * CONSTANTS[AC_K_o] / (1.0 + pow(CONSTANTS[AC_KmNaip] / state_Na_sl, 4.0)) / (CONSTANTS[AC_K_o] + CONSTANTS[AC_KmKo]);
    alg_AV_INaK = alg_AV_INaK_jn + alg_AV_INaK_sl;

    /* nernst */
    alg_AV_ECa_jn = 1.0 / CONSTANTS[AC_FRT] / 2.0 * log(CONSTANTS[AC_Ca_o] / state_Ca_jn);
    alg_AV_ECa_sl = 1.0 / CONSTANTS[AC_FRT] / 2.0 * log(CONSTANTS[AC_Ca_o] / state_Ca_sl);
    alg_AV_EK = 1.0 / CONSTANTS[AC_FRT] * log(CONSTANTS[AC_K_o] / state_K_i);
    alg_AV_ENa_jn = 1.0 / CONSTANTS[AC_FRT] * log(CONSTANTS[AC_Na_o] / state_Na_jn);
    alg_AV_ENa_sl = 1.0 / CONSTANTS[AC_FRT] * log(CONSTANTS[AC_Na_o] / state_Na_sl);

    /* ryr */
    alg_AV_RI = 1.0 - state_ryr_ryr_r - state_ryr_o - state_ryr_i;
    alg_AV_J_SRCarel = CONSTANTS[AC_ks] * state_ryr_o * (state_Ca_sr - state_Ca_jn);
    alg_AV_J_SR_leak = (1.0 + 0.25 * CONSTANTS[AC_AF]) * (state_Ca_sr - state_Ca_jn) * CONSTANTS[AC_J_SR_leak_max];
    alg_AV_J_serca = pow(CONSTANTS[AC_Q10SRCaP], CONSTANTS[AC_Q]) * CONSTANTS[AC_Vmax_SRCaP] * (pow(state_Ca_i / CONSTANTS[AC_Kmf], CONSTANTS[AC_hillSRCaP]) - pow(state_Ca_sr / CONSTANTS[AC_Kmr], CONSTANTS[AC_hillSRCaP])) / (1.0 + pow(state_Ca_i / CONSTANTS[AC_Kmf], CONSTANTS[AC_hillSRCaP]) + pow(state_Ca_sr / CONSTANTS[AC_Kmr], CONSTANTS[AC_hillSRCaP]));
    alg_AV_kCaSR = CONSTANTS[AC_MaxSR] - (CONSTANTS[AC_MaxSR] - CONSTANTS[AC_MinSR]) / (1.0 + pow(CONSTANTS[AC_ec50SR] / state_Ca_sr, 2.5));
    alg_AV_kiSRCa = CONSTANTS[AC_kiCa] * alg_AV_kCaSR;
    alg_AV_koSRCa = CONSTANTS[AC_koCa] / alg_AV_kCaSR;
    rate_ryr_i = alg_AV_kiSRCa * state_Ca_jn * state_ryr_o - CONSTANTS[AC_kim] * state_ryr_i - (CONSTANTS[AC_kom] * state_ryr_i - alg_AV_koSRCa * pow(state_Ca_jn, 2.0) * alg_AV_RI);
    rate_ryr_o = alg_AV_koSRCa * pow(state_Ca_jn, 2.0) * state_ryr_ryr_r - CONSTANTS[AC_kom] * state_ryr_o - (alg_AV_kiSRCa * state_Ca_jn * state_ryr_o - CONSTANTS[AC_kim] * state_ryr_i);
    rate_ryr_ryr_r = CONSTANTS[AC_kim] * alg_AV_RI - alg_AV_kiSRCa * state_Ca_jn * state_ryr_ryr_r - (alg_AV_koSRCa * pow(state_Ca_jn, 2.0) * state_ryr_ryr_r - CONSTANTS[AC_kom] * state_ryr_o);


    /* icab */
    alg_AV_ICaB_jn = CONSTANTS[AC_Fjn] * CONSTANTS[AC_gCaB] * (state_membrane_V - alg_AV_ECa_jn);
    alg_AV_ICaB_sl = CONSTANTS[AC_Fsl] * CONSTANTS[AC_gCaB] * (state_membrane_V - alg_AV_ECa_sl);

    /* ical */
    rate_Ical_d = (alg_AV_ical_d_inf - state_Ical_d) / alg_AV_ical_d_tau;
    rate_Ical_f = (alg_AV_ical_f_inf - state_Ical_f) / alg_AV_ical_f_tau;
    rate_Ical_fCaB_jn = (1.7 * state_Ca_jn * (1.0 - state_Ical_fCaB_jn) - 0.0119 * state_Ical_fCaB_jn) * 1.0;
    rate_Ical_fCaB_sl = (1.7 * state_Ca_sl * (1.0 - state_Ical_fCaB_sl) - 0.0119 * state_Ical_fCaB_sl) * 1.0;
    alg_AV_ibarca_jn = CONSTANTS[AC_pCa] * 4.0 * (state_membrane_V * CONSTANTS[AC_F] * CONSTANTS[AC_FRT]) * (0.341 * state_Ca_jn * exp(2.0 * state_membrane_V * CONSTANTS[AC_FRT]) - 0.341 * CONSTANTS[AC_Ca_o]) / (exp(2.0 * state_membrane_V * CONSTANTS[AC_FRT]) - 1.0) * CONSTANTS[AC_f_conducting];
    alg_AV_ibarca_sl = CONSTANTS[AC_pCa] * 4.0 * (state_membrane_V * CONSTANTS[AC_F] * CONSTANTS[AC_FRT]) * (0.341 * state_Ca_sl * exp(2.0 * state_membrane_V * CONSTANTS[AC_FRT]) - 0.341 * CONSTANTS[AC_Ca_o]) / (exp(2.0 * state_membrane_V * CONSTANTS[AC_FRT]) - 1.0) * CONSTANTS[AC_f_conducting];
    alg_AV_ibark = CONSTANTS[AC_pK] * state_membrane_V * CONSTANTS[AC_F] * CONSTANTS[AC_FRT] * (0.75 * state_K_i * exp(state_membrane_V * CONSTANTS[AC_FRT]) - 0.75 * CONSTANTS[AC_K_o]) / (exp(state_membrane_V * CONSTANTS[AC_FRT]) - 1.0) * CONSTANTS[AC_f_conducting];
    alg_AV_ibarna_jn = CONSTANTS[AC_pNa] * state_membrane_V * CONSTANTS[AC_F] * CONSTANTS[AC_FRT] * (0.75 * state_Na_jn * exp(state_membrane_V * CONSTANTS[AC_FRT]) - 0.75 * CONSTANTS[AC_Na_o]) / (exp(state_membrane_V * CONSTANTS[AC_FRT]) - 1.0) * CONSTANTS[AC_f_conducting];
    alg_AV_ibarna_sl = CONSTANTS[AC_pNa] * state_membrane_V * CONSTANTS[AC_F] * CONSTANTS[AC_FRT] * (0.75 * state_Na_sl * exp(state_membrane_V * CONSTANTS[AC_FRT]) - 0.75 * CONSTANTS[AC_Na_o]) / (exp(state_membrane_V * CONSTANTS[AC_FRT]) - 1.0) * CONSTANTS[AC_f_conducting];
    alg_AV_ICaL_Ca_jn = CONSTANTS[AC_Fjn_CaL] * alg_AV_ibarca_jn * state_Ical_d * state_Ical_f * (1.0 - state_Ical_fCaB_jn + CONSTANTS[AC_fcaCaj]) * pow(CONSTANTS[AC_Q10CaL],
                                                                                                                                     CONSTANTS[AC_Q]) * 0.45;
    alg_AV_ICaL_Ca_sl = CONSTANTS[AC_Fsl_CaL] * alg_AV_ibarca_sl * state_Ical_d * state_Ical_f * (1.0 - state_Ical_fCaB_sl + CONSTANTS[AC_fcaCaMSL]) * pow(CONSTANTS[AC_Q10CaL],
                                                                                                                                       CONSTANTS[AC_Q]) * 0.45;
    alg_AV_ICaL_K = alg_AV_ibark * state_Ical_d * state_Ical_f * (CONSTANTS[AC_Fjn_CaL] * (CONSTANTS[AC_fcaCaj] + (1.0 - state_Ical_fCaB_jn)) + CONSTANTS[AC_Fsl_CaL] * (CONSTANTS[AC_fcaCaMSL] + (1.0 - state_Ical_fCaB_sl))) * pow(CONSTANTS[AC_Q10CaL],
                                                                                                                                                                                                         CONSTANTS[AC_Q]) * 0.45;
    alg_AV_ICaL_Na_jn = CONSTANTS[AC_Fjn_CaL] * alg_AV_ibarna_jn * state_Ical_d * state_Ical_f * (1.0 - state_Ical_fCaB_jn + CONSTANTS[AC_fcaCaj]) * pow(CONSTANTS[AC_Q10CaL],
                                                                                                                                     CONSTANTS[AC_Q]) * 0.45;
    alg_AV_ICaL_Na_sl = CONSTANTS[AC_Fsl_CaL] * alg_AV_ibarna_sl * state_Ical_d * state_Ical_f * (1.0 - state_Ical_fCaB_sl + CONSTANTS[AC_fcaCaMSL]) * pow(CONSTANTS[AC_Q10CaL],
                                                                                                                                       CONSTANTS[AC_Q]) * 0.45;

    /* iclb */
    alg_AV_IClB = CONSTANTS[AC_gClB] * (state_membrane_V - CONSTANTS[AC_ECl]);

    /* iclca */
    alg_AV_IClCa_jn = CONSTANTS[AC_Fjn] * CONSTANTS[AC_GClCa] / (1.0 + CONSTANTS[AC_KdClCa] / state_Ca_jn) * (state_membrane_V - CONSTANTS[AC_ECl]);
    alg_AV_IClCa_sl = CONSTANTS[AC_Fsl] * CONSTANTS[AC_GClCa] / (1.0 + CONSTANTS[AC_KdClCa] / state_Ca_sl) * (state_membrane_V - CONSTANTS[AC_ECl]);
    alg_AV_IClCa = alg_AV_IClCa_jn + alg_AV_IClCa_sl;

    /* ik1_ik1_inf */
    alg_AV_ik1_inf_a = 1.0 * (1.02 / (1.0 + exp(0.2385 * (state_membrane_V - alg_AV_EK - 59.215))));
    alg_AV_ik1_inf_b = 1.0 * ((0.49124 * exp(0.08032 * (state_membrane_V + 5.476 - alg_AV_EK)) + exp(0.06175 * (state_membrane_V - alg_AV_EK - 594.31))) / (1.0 + exp(-0.5143 * (state_membrane_V - alg_AV_EK + 4.753))));

    /* ikp */
    alg_AV_kp = 1.0 / (1.0 + exp(7.488 - state_membrane_V / 5.98));
    alg_AV_IKp_jn = CONSTANTS[AC_Fjn] * CONSTANTS[AC_gKp] * alg_AV_kp * (state_membrane_V - alg_AV_EK);
    alg_AV_IKp_sl = CONSTANTS[AC_Fsl] * CONSTANTS[AC_gKp] * alg_AV_kp * (state_membrane_V - alg_AV_EK);
    alg_AV_IKp = alg_AV_IKp_jn + alg_AV_IKp_sl;

    /* ikr */
    rate_Ikr_xr = (alg_AV_ikr_xr_inf - state_Ikr_xr) / alg_AV_ikr_xr_tau;
    alg_AV_rr = 1.0 / (1.0 + exp((state_membrane_V + 74.0) / 24.0));
    alg_AV_IKr = CONSTANTS[AC_gKr] * state_Ikr_xr * alg_AV_rr * (state_membrane_V - alg_AV_EK);

    /* iks */
    rate_Iks_xs = (alg_AV_iks_xs_inf - state_Iks_xs) / alg_AV_iks_xs_tau;
    alg_AV_EKs = 1.0 / CONSTANTS[AC_FRT] * log((CONSTANTS[AC_K_o] + CONSTANTS[AC_pNaK] * CONSTANTS[AC_Na_o]) / (state_K_i + CONSTANTS[AC_pNaK] * state_Na_i));
    alg_AV_IKs_jn = CONSTANTS[AC_Fjn] * CONSTANTS[AC_gKs_jn] * pow(state_Iks_xs, 2.0) * (state_membrane_V - alg_AV_EKs);
    alg_AV_IKs_sl = CONSTANTS[AC_Fsl] * CONSTANTS[AC_gKs_sl] * pow(state_Iks_xs, 2.0) * (state_membrane_V - alg_AV_EKs);
    alg_AV_IKs = alg_AV_IKs_jn + alg_AV_IKs_sl;

    /* ikur */
    rate_Ikur_ikur_r = (alg_AV_ikur_r_inf - state_Ikur_ikur_r) / alg_AV_ikur_r_tau;
    rate_Ikur_s = (alg_AV_ikur_s_inf - state_Ikur_s) / alg_AV_ikur_s_tau;
    alg_AV_IKur = CONSTANTS[AC_gKur] * state_Ikur_ikur_r * state_Ikur_s * (state_membrane_V - alg_AV_EK);

    /* ina */
    rate_Ina_h = (alg_AV_ina_h_inf - state_Ina_h) / alg_AV_ina_h_tau;
    rate_Ina_j = (alg_AV_ina_j_inf - state_Ina_j) / alg_AV_ina_j_tau;
    rate_Ina_m = (alg_AV_ina_m_inf - state_Ina_m) / alg_AV_ina_m_tau;
    alg_AV_INa_jn = CONSTANTS[AC_Fjn] * CONSTANTS[AC_gNa] * pow(state_Ina_m, 3.0) * state_Ina_h * state_Ina_j * (state_membrane_V - alg_AV_ENa_jn);
    alg_AV_INa_sl = CONSTANTS[AC_Fsl] * CONSTANTS[AC_gNa] * pow(state_Ina_m, 3.0) * state_Ina_h * state_Ina_j * (state_membrane_V - alg_AV_ENa_sl);

    /* inab */
    alg_AV_INaB_jn = CONSTANTS[AC_Fjn] * CONSTANTS[AC_gNaB] * (state_membrane_V - alg_AV_ENa_jn);
    alg_AV_INaB_sl = CONSTANTS[AC_Fsl] * CONSTANTS[AC_gNaB] * (state_membrane_V - alg_AV_ENa_sl);

    /* inal */
    rate_Inal_hl = (alg_AV_inal_hl_inf - state_Inal_hl) / CONSTANTS[AC_inal_hl_tau];
    rate_Inal_ml = alg_AV_inal_ml_a * (1.0 - state_Inal_ml) - alg_AV_inal_ml_b * state_Inal_ml;
    alg_AV_INaL_jn = CONSTANTS[AC_Fjn] * CONSTANTS[AC_gNaL] * pow(state_Inal_ml, 3.0) * state_Inal_hl * (state_membrane_V - alg_AV_ENa_jn);
    alg_AV_INaL_sl = CONSTANTS[AC_Fsl] * CONSTANTS[AC_gNaL] * pow(state_Inal_ml, 3.0) * state_Inal_hl * (state_membrane_V - alg_AV_ENa_sl);

    /* ito */
    rate_Ito_x = (alg_AV_ito_x_inf - state_Ito_x) / alg_AV_ito_x_tau;
    rate_Ito_y = (alg_AV_ito_y_inf - state_Ito_y) / alg_AV_ito_y_tau;
    alg_AV_Ito = CONSTANTS[AC_gto] * state_Ito_x * state_Ito_y * (state_membrane_V - alg_AV_EK);

    /* ik1 */
    alg_AV_ik1_inf = alg_AV_ik1_inf_a / (alg_AV_ik1_inf_a + alg_AV_ik1_inf_b);
    alg_AV_IK1 = CONSTANTS[AC_gK1] * alg_AV_ik1_inf * (state_membrane_V - alg_AV_EK);


    /* calcium */
    rate_Ca_i = -alg_AV_J_serca * CONSTANTS[AC_Vsr] / CONSTANTS[AC_Vmyo] - alg_AV_J_CaB_cytosol + CONSTANTS[AC_JCa_slmyo] / CONSTANTS[AC_Vmyo] * (state_Ca_sl - state_Ca_i);
    rate_Csqn = CONSTANTS[AC_kon_csqn] * state_Ca_sr * (CONSTANTS[AC_Bmax_Csqn] - state_Csqn) - CONSTANTS[AC_koff_csqn] * state_Csqn;
    rate_Ca_sr = alg_AV_J_serca - (alg_AV_J_SR_leak * CONSTANTS[AC_Vmyo] / CONSTANTS[AC_Vsr] + alg_AV_J_SRCarel) - rate_Csqn;

    /* inaca */
    alg_AV_Ka_jn = 1.0 / (1.0 + pow(CONSTANTS[AC_Kdact] / state_Ca_jn, 2.0));
    alg_AV_Ka_sl = 1.0 / (1.0 + pow(CONSTANTS[AC_Kdact] / state_Ca_sl, 2.0));

    /* inaca_INaCa_jn */
    alg_AV_s1_jn = exp(CONSTANTS[AC_nu] * state_membrane_V * CONSTANTS[AC_FRT]) * pow(state_Na_jn, 3.0) * CONSTANTS[AC_Ca_o];
    alg_AV_s2_jn = exp((CONSTANTS[AC_nu] - 1.0) * state_membrane_V * CONSTANTS[AC_FRT]) * pow(CONSTANTS[AC_Na_o],
                                                                                    3.0) * state_Ca_jn;
    alg_AV_s3_jn = CONSTANTS[AC_KmCai] * pow(CONSTANTS[AC_Na_o], 3.0) * (1.0 + pow(state_Na_jn / CONSTANTS[AC_KmNai], 3.0)) + pow(CONSTANTS[AC_KmNao],
                                                                                                                             3.0) * state_Ca_jn * (1.0 + state_Ca_jn / CONSTANTS[AC_KmCai]) + CONSTANTS[AC_KmCao] * pow(state_Na_jn,
                                                                                                                                                                         3.0) + pow(state_Na_jn,
                                                                                                                                                                                    3.0) * CONSTANTS[AC_Ca_o] + pow(CONSTANTS[AC_Na_o],
                                                                                                                                                                                                                    3.0) * state_Ca_jn;

    /* inaca_INaCa_sl */
    alg_AV_s1_sl = exp(CONSTANTS[AC_nu] * state_membrane_V * CONSTANTS[AC_FRT]) * pow(state_Na_sl, 3.0) * CONSTANTS[AC_Ca_o];
    alg_AV_s2_sl = exp((CONSTANTS[AC_nu] - 1.0) * state_membrane_V * CONSTANTS[AC_FRT]) * pow(CONSTANTS[AC_Na_o],
                                                                                    3.0) * state_Ca_sl;
    alg_AV_s3_sl = CONSTANTS[AC_KmCai] * pow(CONSTANTS[AC_Na_o], 3.0) * (1.0 + pow(state_Na_sl / CONSTANTS[AC_KmNai], 3.0)) + pow(CONSTANTS[AC_KmNao],
                                                                                                                             3.0) * state_Ca_sl * (1.0 + state_Ca_sl / CONSTANTS[AC_KmCai]) + CONSTANTS[AC_KmCao] * pow(state_Na_sl,
                                                                                                                                                                         3.0) + pow(state_Na_sl,
                                                                                                                                                                                    3.0) * CONSTANTS[AC_Ca_o] + pow(CONSTANTS[AC_Na_o],
                                                                                                                                                                                                                    3.0) * state_Ca_sl;
    /* ipca_IpCa_jn */
    alg_AV_ipca_IpCa_jn_a = pow(state_Ca_jn * 1.0, 1.6);

    /* ipca_IpCa_sl */
    alg_AV_ipca_IpCa_sl_a = pow(state_Ca_sl * 1.0, 1.6);

    /* sodium */
    rate_Na_i = CONSTANTS[AC_JNa_slmyo] / CONSTANTS[AC_Vmyo] * (state_Na_sl - state_Na_i);

    /* *remaining* */
    alg_AV_INaCa_jn = CONSTANTS[AC_Fjn] * CONSTANTS[AC_IbarNCX] * pow(CONSTANTS[AC_Q10NCX],
                                                                  CONSTANTS[AC_Q]) * alg_AV_Ka_jn * (alg_AV_s1_jn - alg_AV_s2_jn) / alg_AV_s3_jn / (1.0 + CONSTANTS[AC_ksat] * exp((CONSTANTS[AC_nu] - 1.0) * state_membrane_V * CONSTANTS[AC_FRT]));
    alg_AV_INaCa_sl = CONSTANTS[AC_Fsl] * CONSTANTS[AC_IbarNCX] * pow(CONSTANTS[AC_Q10NCX],
                                                                  CONSTANTS[AC_Q]) * alg_AV_Ka_sl * (alg_AV_s1_sl - alg_AV_s2_sl) / alg_AV_s3_sl / (1.0 + CONSTANTS[AC_ksat] * exp((CONSTANTS[AC_nu] - 1.0) * state_membrane_V * CONSTANTS[AC_FRT]));
    alg_AV_IpCa_jn = CONSTANTS[AC_Fjn] * pow(CONSTANTS[AC_Q10SLCaP],
                                         CONSTANTS[AC_Q]) * CONSTANTS[AC_IbarSLCaP] * alg_AV_ipca_IpCa_jn_a / (alg_AV_ipca_IpCa_jn_a + CONSTANTS[AC_ipca_IpCa_jn_b]);
    alg_AV_IpCa_sl = CONSTANTS[AC_Fsl] * pow(CONSTANTS[AC_Q10SLCaP],
                                         CONSTANTS[AC_Q]) * CONSTANTS[AC_IbarSLCaP] * alg_AV_ipca_IpCa_sl_a / (alg_AV_ipca_IpCa_sl_a + CONSTANTS[AC_ipca_IpCa_sl_b]);
    alg_AV_ICa_tot_jn = alg_AV_ICaL_Ca_jn + alg_AV_ICaB_jn + alg_AV_IpCa_jn - 2.0 * alg_AV_INaCa_jn;
    alg_AV_ICa_tot_sl = alg_AV_ICaL_Ca_sl + alg_AV_ICaB_sl + alg_AV_IpCa_sl - 2.0 * alg_AV_INaCa_sl;
    alg_AV_INa_tot_jn = alg_AV_INa_jn + alg_AV_INaB_jn + 3.0 * alg_AV_INaCa_jn + 3.0 * alg_AV_INaK_jn + alg_AV_ICaL_Na_jn + alg_AV_INaL_jn;
    alg_AV_INa_tot_sl = alg_AV_INa_sl + alg_AV_INaB_sl + 3.0 * alg_AV_INaCa_sl + 3.0 * alg_AV_INaK_sl + alg_AV_ICaL_Na_sl + alg_AV_INaL_sl;
    rate_Ca_jn = -alg_AV_ICa_tot_jn * CONSTANTS[AC_C] / (CONSTANTS[AC_Vjn] * 2.0 * CONSTANTS[AC_F]) + CONSTANTS[AC_JCa_jnsl] / CONSTANTS[AC_Vjn] * (state_Ca_sl - state_Ca_jn) - alg_AV_J_CaB_jntion + alg_AV_J_SRCarel * CONSTANTS[AC_Vsr] / CONSTANTS[AC_Vjn] + alg_AV_J_SR_leak * CONSTANTS[AC_Vmyo] / CONSTANTS[AC_Vjn];
    rate_Ca_sl = -alg_AV_ICa_tot_sl * CONSTANTS[AC_C] / (CONSTANTS[AC_Vsl] * 2.0 * CONSTANTS[AC_F]) + CONSTANTS[AC_JCa_jnsl] / CONSTANTS[AC_Vsl] * (state_Ca_jn - state_Ca_sl) + CONSTANTS[AC_JCa_slmyo] / CONSTANTS[AC_Vsl] * (state_Ca_i - state_Ca_sl) - alg_AV_J_CaB_sl;

    rate_Na_jn = -alg_AV_INa_tot_jn * CONSTANTS[AC_C] / (CONSTANTS[AC_Vjn] * CONSTANTS[AC_F]) + CONSTANTS[AC_JNa_jnsl] / CONSTANTS[AC_Vjn] * (state_Na_sl - state_Na_jn) - rate_buffNa_NaB_jn;
    rate_Na_sl = -alg_AV_INa_tot_sl * CONSTANTS[AC_C] / (CONSTANTS[AC_Vsl] * CONSTANTS[AC_F]) + CONSTANTS[AC_JNa_jnsl] / CONSTANTS[AC_Vsl] * (state_Na_jn - state_Na_sl) + CONSTANTS[AC_JNa_slmyo] / CONSTANTS[AC_Vsl] * (state_Na_i - state_Na_sl) - rate_buffNa_NaB_sl;

    /* potassium */
    rate_K_i = 0.0;
    alg_AV_IK_tot = alg_AV_Ito + alg_AV_IKr + alg_AV_IKs + alg_AV_IK1 - 2.0 * alg_AV_INaK + alg_AV_ICaL_K + alg_AV_IKp + alg_AV_IKur;
    alg_AV_ICl_tot = alg_AV_IClCa + alg_AV_IClB;
    alg_AV_ICa_tot = alg_AV_ICa_tot_jn + alg_AV_ICa_tot_sl;
    alg_AV_INa_tot = alg_AV_INa_tot_jn + alg_AV_INa_tot_sl;
    alg_Iion_cm = alg_AV_INa_tot + alg_AV_ICl_tot + alg_AV_ICa_tot + alg_AV_IK_tot;


    alg_Istim = computeStimulusPOD(VOI, stimulus);

    if (solveVm)
    {
        rate_membrane_V = -alg_Iion_cm - alg_Istim;
    }

        RATES_SoA[Ca_i * N + cellI] = rate_Ca_i;
        RATES_SoA[Ca_jn * N + cellI] = rate_Ca_jn;
        RATES_SoA[Ca_sl * N + cellI] = rate_Ca_sl;
        RATES_SoA[Ca_sr * N + cellI] = rate_Ca_sr;
        RATES_SoA[Csqn * N + cellI] = rate_Csqn;
        RATES_SoA[Ical_d * N + cellI] = rate_Ical_d;
        RATES_SoA[Ical_f * N + cellI] = rate_Ical_f;
        RATES_SoA[Ical_fCaB_jn * N + cellI] = rate_Ical_fCaB_jn;
        RATES_SoA[Ical_fCaB_sl * N + cellI] = rate_Ical_fCaB_sl;
        RATES_SoA[Ikr_xr * N + cellI] = rate_Ikr_xr;
        RATES_SoA[Iks_xs * N + cellI] = rate_Iks_xs;
        RATES_SoA[Ikur_ikur_r * N + cellI] = rate_Ikur_ikur_r;
        RATES_SoA[Ikur_s * N + cellI] = rate_Ikur_s;
        RATES_SoA[Ina_h * N + cellI] = rate_Ina_h;
        RATES_SoA[Ina_j * N + cellI] = rate_Ina_j;
        RATES_SoA[Ina_m * N + cellI] = rate_Ina_m;
        RATES_SoA[Inal_hl * N + cellI] = rate_Inal_hl;
        RATES_SoA[Inal_ml * N + cellI] = rate_Inal_ml;
        RATES_SoA[Ito_x * N + cellI] = rate_Ito_x;
        RATES_SoA[Ito_y * N + cellI] = rate_Ito_y;
        RATES_SoA[K_i * N + cellI] = rate_K_i;
        RATES_SoA[Na_i * N + cellI] = rate_Na_i;
        RATES_SoA[Na_jn * N + cellI] = rate_Na_jn;
        RATES_SoA[Na_sl * N + cellI] = rate_Na_sl;
        RATES_SoA[buffCa_CaM * N + cellI] = rate_buffCa_CaM;
        RATES_SoA[buffCa_Myoc * N + cellI] = rate_buffCa_Myoc;
        RATES_SoA[buffCa_Myom * N + cellI] = rate_buffCa_Myom;
        RATES_SoA[buffCa_SLH_jn * N + cellI] = rate_buffCa_SLH_jn;
        RATES_SoA[buffCa_SLH_sl * N + cellI] = rate_buffCa_SLH_sl;
        RATES_SoA[buffCa_SLL_jn * N + cellI] = rate_buffCa_SLL_jn;
        RATES_SoA[buffCa_SLL_sl * N + cellI] = rate_buffCa_SLL_sl;
        RATES_SoA[buffCa_SRB * N + cellI] = rate_buffCa_SRB;
        RATES_SoA[buffCa_TnCHc * N + cellI] = rate_buffCa_TnCHc;
        RATES_SoA[buffCa_TnCHm * N + cellI] = rate_buffCa_TnCHm;
        RATES_SoA[buffCa_TnCL * N + cellI] = rate_buffCa_TnCL;
        RATES_SoA[buffNa_NaB_jn * N + cellI] = rate_buffNa_NaB_jn;
        RATES_SoA[buffNa_NaB_sl * N + cellI] = rate_buffNa_NaB_sl;
        RATES_SoA[membrane_V * N + cellI] = rate_membrane_V;
        RATES_SoA[ryr_i * N + cellI] = rate_ryr_i;
        RATES_SoA[ryr_o * N + cellI] = rate_ryr_o;
        RATES_SoA[ryr_ryr_r * N + cellI] = rate_ryr_ryr_r;

        if (SUPPORT_SoA != nullptr)
        {
            SUPPORT_SoA[GRANDI_BATCH_SUPPORT_tau_Ikr_xr * N + cellI] = alg_AV_ikr_xr_tau;
            SUPPORT_SoA[GRANDI_BATCH_SUPPORT_tau_Iks_xs * N + cellI] = alg_AV_iks_xs_tau;
            SUPPORT_SoA[GRANDI_BATCH_SUPPORT_tau_Ikur_ikur_r * N + cellI] = alg_AV_ikur_r_tau;
            SUPPORT_SoA[GRANDI_BATCH_SUPPORT_tau_Ikur_s * N + cellI] = alg_AV_ikur_s_tau;
            SUPPORT_SoA[GRANDI_BATCH_SUPPORT_tau_Ina_h * N + cellI] = alg_AV_ina_h_tau;
            SUPPORT_SoA[GRANDI_BATCH_SUPPORT_tau_Ina_j * N + cellI] = alg_AV_ina_j_tau;
            SUPPORT_SoA[GRANDI_BATCH_SUPPORT_tau_Ina_m * N + cellI] = alg_AV_ina_m_tau;
            SUPPORT_SoA[GRANDI_BATCH_SUPPORT_tau_Inal_hl * N + cellI] = CONSTANTS[AC_inal_hl_tau];
            SUPPORT_SoA[GRANDI_BATCH_SUPPORT_tau_Inal_ml * N + cellI] = 1.0 / (alg_AV_inal_ml_a + alg_AV_inal_ml_b);
            SUPPORT_SoA[GRANDI_BATCH_SUPPORT_tau_Ical_d * N + cellI] = alg_AV_ical_d_tau;
            SUPPORT_SoA[GRANDI_BATCH_SUPPORT_tau_Ical_f * N + cellI] = alg_AV_ical_f_tau;
            SUPPORT_SoA[GRANDI_BATCH_SUPPORT_tau_Ical_fCaB_jn * N + cellI] = 1.0 / (1.7 * state_Ca_jn + 0.0119);
            SUPPORT_SoA[GRANDI_BATCH_SUPPORT_tau_Ical_fCaB_sl * N + cellI] = 1.0 / (1.7 * state_Ca_sl + 0.0119);
            SUPPORT_SoA[GRANDI_BATCH_SUPPORT_tau_Ito_x * N + cellI] = alg_AV_ito_x_tau;
            SUPPORT_SoA[GRANDI_BATCH_SUPPORT_tau_Ito_y * N + cellI] = alg_AV_ito_y_tau;
            SUPPORT_SoA[GRANDI_BATCH_SUPPORT_gInf_Ikr_xr * N + cellI] = alg_AV_ikr_xr_inf;
            SUPPORT_SoA[GRANDI_BATCH_SUPPORT_gInf_Iks_xs * N + cellI] = alg_AV_iks_xs_inf;
            SUPPORT_SoA[GRANDI_BATCH_SUPPORT_gInf_Ikur_ikur_r * N + cellI] = alg_AV_ikur_r_inf;
            SUPPORT_SoA[GRANDI_BATCH_SUPPORT_gInf_Ikur_s * N + cellI] = alg_AV_ikur_s_inf;
            SUPPORT_SoA[GRANDI_BATCH_SUPPORT_gInf_Ina_h * N + cellI] = alg_AV_ina_h_inf;
            SUPPORT_SoA[GRANDI_BATCH_SUPPORT_gInf_Ina_j * N + cellI] = alg_AV_ina_j_inf;
            SUPPORT_SoA[GRANDI_BATCH_SUPPORT_gInf_Ina_m * N + cellI] = alg_AV_ina_m_inf;
            SUPPORT_SoA[GRANDI_BATCH_SUPPORT_gInf_Inal_hl * N + cellI] = alg_AV_inal_hl_inf;
            SUPPORT_SoA[GRANDI_BATCH_SUPPORT_gInf_Inal_ml * N + cellI] = alg_AV_inal_ml_a / (alg_AV_inal_ml_a + alg_AV_inal_ml_b);
            SUPPORT_SoA[GRANDI_BATCH_SUPPORT_gInf_Ical_d * N + cellI] = alg_AV_ical_d_inf;
            SUPPORT_SoA[GRANDI_BATCH_SUPPORT_gInf_Ical_f * N + cellI] = alg_AV_ical_f_inf;
            SUPPORT_SoA[GRANDI_BATCH_SUPPORT_gInf_Ical_fCaB_jn * N + cellI] = 1.7 * state_Ca_jn / (1.7 * state_Ca_jn + 0.0119);
            SUPPORT_SoA[GRANDI_BATCH_SUPPORT_gInf_Ical_fCaB_sl * N + cellI] = 1.7 * state_Ca_sl / (1.7 * state_Ca_sl + 0.0119);
            SUPPORT_SoA[GRANDI_BATCH_SUPPORT_gInf_Ito_x * N + cellI] = alg_AV_ito_x_inf;
            SUPPORT_SoA[GRANDI_BATCH_SUPPORT_gInf_Ito_y * N + cellI] = alg_AV_ito_y_inf;
            SUPPORT_SoA[GRANDI_BATCH_SUPPORT_Iion_cm * N + cellI] = alg_Iion_cm;
        }
    }
}

} // End namespace Foam

// ************************************************************************* //
