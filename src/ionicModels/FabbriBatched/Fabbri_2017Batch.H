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
    State-major (SoA) batched evaluator for the Fabbri 2017 ionic model.
\*---------------------------------------------------------------------------*/

#pragma once

#include "Fabbri_2017Names.H"
#include "StimulusProtocolPOD.H"
#include <cmath>

namespace Foam
{

enum FABBRI_BATCH_SUPPORT_INDEX
{
    FABBRI_BATCH_SUPPORT_tau_y,
    FABBRI_BATCH_SUPPORT_gInf_y,
    FABBRI_BATCH_SUPPORT_tau_m,
    FABBRI_BATCH_SUPPORT_gInf_m,
    FABBRI_BATCH_SUPPORT_tau_h,
    FABBRI_BATCH_SUPPORT_gInf_h,
    FABBRI_BATCH_SUPPORT_tau_dL,
    FABBRI_BATCH_SUPPORT_gInf_dL,
    FABBRI_BATCH_SUPPORT_tau_fL,
    FABBRI_BATCH_SUPPORT_gInf_fL,
    FABBRI_BATCH_SUPPORT_tau_fCa,
    FABBRI_BATCH_SUPPORT_gInf_fCa,
    FABBRI_BATCH_SUPPORT_tau_dT,
    FABBRI_BATCH_SUPPORT_gInf_dT,
    FABBRI_BATCH_SUPPORT_tau_fT,
    FABBRI_BATCH_SUPPORT_gInf_fT,
    FABBRI_BATCH_SUPPORT_tau_r_Kur,
    FABBRI_BATCH_SUPPORT_gInf_r_Kur,
    FABBRI_BATCH_SUPPORT_tau_s_Kur,
    FABBRI_BATCH_SUPPORT_gInf_s_Kur,
    FABBRI_BATCH_SUPPORT_tau_q,
    FABBRI_BATCH_SUPPORT_gInf_q,
    FABBRI_BATCH_SUPPORT_tau_r,
    FABBRI_BATCH_SUPPORT_gInf_r,
    FABBRI_BATCH_SUPPORT_tau_paS,
    FABBRI_BATCH_SUPPORT_gInf_paS,
    FABBRI_BATCH_SUPPORT_tau_paF,
    FABBRI_BATCH_SUPPORT_gInf_paF,
    FABBRI_BATCH_SUPPORT_tau_pi,
    FABBRI_BATCH_SUPPORT_gInf_pi,
    FABBRI_BATCH_SUPPORT_tau_n,
    FABBRI_BATCH_SUPPORT_gInf_n,
    FABBRI_BATCH_SUPPORT_tau_a,
    FABBRI_BATCH_SUPPORT_gInf_a,
    FABBRI_BATCH_SUPPORT_Iion_cm,

    NUM_FABBRI_BATCH_SUPPORT
};

static_assert(NUM_FABBRI_BATCH_SUPPORT == 35, "Fabbri compact support layout drift");

CARDIAC_HOST_DEVICE inline void FabbriComputeVariablesBatch
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
        const double state_Buffer_fCMi = STATES_SoA[Buffer_fCMi * N + cellI];
        const double state_Buffer_fCMs = STATES_SoA[Buffer_fCMs * N + cellI];
        const double state_Buffer_fCQ = STATES_SoA[Buffer_fCQ * N + cellI];
        const double state_Buffer_fTC = STATES_SoA[Buffer_fTC * N + cellI];
        const double state_Buffer_fTMC = STATES_SoA[Buffer_fTMC * N + cellI];
        const double state_Buffer_fTMM = STATES_SoA[Buffer_fTMM * N + cellI];
        const double state_Ca_i = STATES_SoA[Ca_i * N + cellI];
        const double state_Ca_jsr = STATES_SoA[Ca_jsr * N + cellI];
        const double state_Ca_nsr = STATES_SoA[Ca_nsr * N + cellI];
        const double state_Ca_sub = STATES_SoA[Ca_sub * N + cellI];
        const double state_ICaL_dL_gate_dL = STATES_SoA[ICaL_dL_gate_dL * N + cellI];
        const double state_ICaL_fCa_gate_fCa = STATES_SoA[ICaL_fCa_gate_fCa * N + cellI];
        const double state_ICaL_fL_gate_fL = STATES_SoA[ICaL_fL_gate_fL * N + cellI];
        const double state_ICaT_dT_gate_dT = STATES_SoA[ICaT_dT_gate_dT * N + cellI];
        const double state_ICaT_fT_gate_fT = STATES_SoA[ICaT_fT_gate_fT * N + cellI];
        const double state_IKACh_a_gate_a = STATES_SoA[IKACh_a_gate_a * N + cellI];
        const double state_IKr_pa_gate_paF = STATES_SoA[IKr_pa_gate_paF * N + cellI];
        const double state_IKr_pa_gate_paS = STATES_SoA[IKr_pa_gate_paS * N + cellI];
        const double state_IKr_pi_gate_piy = STATES_SoA[IKr_pi_gate_piy * N + cellI];
        const double state_IKs_n_gate_n = STATES_SoA[IKs_n_gate_n * N + cellI];
        const double state_IKur_rKur_gate_r_Kur = STATES_SoA[IKur_rKur_gate_r_Kur * N + cellI];
        const double state_IKur_sKur_gate_s_Kur = STATES_SoA[IKur_sKur_gate_s_Kur * N + cellI];
        const double state_INa_h_gate_h = STATES_SoA[INa_h_gate_h * N + cellI];
        const double state_INa_m_gate_m = STATES_SoA[INa_m_gate_m * N + cellI];
        const double state_If_y_gate_y = STATES_SoA[If_y_gate_y * N + cellI];
        const double state_Ito_q_gate_q = STATES_SoA[Ito_q_gate_q * N + cellI];
        const double state_Ito_r_gate_r = STATES_SoA[Ito_r_gate_r * N + cellI];
        const double state_Na_i = STATES_SoA[Na_i * N + cellI];
        const double state_SR_I = STATES_SoA[SR_I * N + cellI];
        const double state_SR_O = STATES_SoA[SR_O * N + cellI];
        const double state_SR_R = STATES_SoA[SR_R * N + cellI];
        const double state_SR_RI = STATES_SoA[SR_RI * N + cellI];
        const double state_membrane_V = STATES_SoA[membrane_V * N + cellI];

        double alg_AV_E0_m = 0.0;
        double alg_AV_E_Ks = 0.0;
        double alg_AV_E_Na = 0.0;
        double alg_AV_E_mh = 0.0;
        double alg_AV_Nai = 0.0;
        double alg_AV_V = 0.0;
        double alg_AV_V_clamp = 0.0;
        double alg_AV_a_infinity = 0.0;
        double alg_AV_adVm = 0.0;
        double alg_AV_alpha_dL = 0.0;
        double alg_AV_alpha_h = 0.0;
        double alg_AV_alpha_m = 0.0;
        double alg_AV_alpha_n = 0.0;
        double alg_AV_bdVm = 0.0;
        double alg_AV_beta_a = 0.0;
        double alg_AV_beta_dL = 0.0;
        double alg_AV_beta_h = 0.0;
        double alg_AV_beta_m = 0.0;
        double alg_AV_beta_n = 0.0;
        double alg_AV_dL_infinity = 0.0;
        double alg_AV_dT_infinity = 0.0;
        double alg_AV_delta_fCMi = 0.0;
        double alg_AV_delta_fCMs = 0.0;
        double alg_AV_delta_fCQ = 0.0;
        double alg_AV_delta_fTC = 0.0;
        double alg_AV_delta_fTMC = 0.0;
        double alg_AV_delta_fTMM = 0.0;
        double alg_AV_di = 0.0;
        double alg_AV_fCa_infinity = 0.0;
        double alg_AV_fL_infinity = 0.0;
        double alg_AV_fT_infinity = 0.0;
        double alg_AV_h_infinity = 0.0;
        double alg_AV_i_CaL_i_CaL = 0.0;
        double alg_AV_i_CaT_i_CaT = 0.0;
        double alg_AV_i_KACh_i_KACh = 0.0;
        double alg_AV_i_Kr_i_Kr = 0.0;
        double alg_AV_i_Ks_i_Ks = 0.0;
        double alg_AV_i_Kur_i_Kur = 0.0;
        double alg_AV_i_NaCa_do = 0.0;
        double alg_AV_i_NaCa_i_NaCa = 0.0;
        double alg_AV_i_NaK_i_NaK = 0.0;
        double alg_AV_i_Na_ = 0.0;
        double alg_AV_i_Na_L = 0.0;
        double alg_AV_i_Na_i_Na = 0.0;
        double alg_AV_i_fK = 0.0;
        double alg_AV_i_fNa = 0.0;
        double alg_AV_i_f_i_f = 0.0;
        double alg_AV_i_siCa = 0.0;
        double alg_AV_i_siK = 0.0;
        double alg_AV_i_siNa = 0.0;
        double alg_AV_i_to_i_to = 0.0;
        double alg_AV_j_Ca_dif = 0.0;
        double alg_AV_j_SRCarel = 0.0;
        double alg_AV_j_tr = 0.0;
        double alg_AV_j_up = 0.0;
        double alg_AV_k12 = 0.0;
        double alg_AV_k14 = 0.0;
        double alg_AV_k21 = 0.0;
        double alg_AV_k23 = 0.0;
        double alg_AV_k32 = 0.0;
        double alg_AV_k41 = 0.0;
        double alg_AV_k43 = 0.0;
        double alg_AV_kCaSR = 0.0;
        double alg_AV_kiSRCa = 0.0;
        double alg_AV_koSRCa = 0.0;
        double alg_AV_m_infinity = 0.0;
        double alg_AV_n_infinity = 0.0;
        double alg_AV_pa_infinity = 0.0;
        double alg_AV_pi_infinity = 0.0;
        double alg_AV_q_infinity = 0.0;
        double alg_AV_r_Kur_infinity = 0.0;
        double alg_AV_r_infinity = 0.0;
        double alg_AV_s_Kur_infinity = 0.0;
        double alg_AV_tau_a = 0.0;
        double alg_AV_tau_dL = 0.0;
        double alg_AV_tau_dT = 0.0;
        double alg_AV_tau_fCa = 0.0;
        double alg_AV_tau_fL = 0.0;
        double alg_AV_tau_fT = 0.0;
        double alg_AV_tau_h = 0.0;
        double alg_AV_tau_m = 0.0;
        double alg_AV_tau_n = 0.0;
        double alg_AV_tau_paF = 0.0;
        double alg_AV_tau_paS = 0.0;
        double alg_AV_tau_pi = 0.0;
        double alg_AV_tau_q = 0.0;
        double alg_AV_tau_r = 0.0;
        double alg_AV_tau_r_Kur = 0.0;
        double alg_AV_tau_s_Kur = 0.0;
        double alg_AV_tau_y = 0.0;
        double alg_AV_x1 = 0.0;
        double alg_AV_x2 = 0.0;
        double alg_AV_x3 = 0.0;
        double alg_AV_x4 = 0.0;
        double alg_AV_y_infinity = 0.0;
        double alg_Iion_cm = 0.0;
        double alg_Istim = 0.0;



        /* Ca_SR_release */
            alg_AV_j_SRCarel = CONSTANTS[AC_ks] * state_SR_O * (state_Ca_jsr - state_Ca_sub);
            alg_AV_kCaSR = CONSTANTS[AC_MaxSR] - (CONSTANTS[AC_MaxSR] - CONSTANTS[AC_MinSR]) / (1.0 + pow(CONSTANTS[AC_EC50_SR] / (state_Ca_jsr > 1e-15 ? state_Ca_jsr : 1e-15), CONSTANTS[AC_HSR]));
            alg_AV_kiSRCa = CONSTANTS[AC_kiCa] * alg_AV_kCaSR;
            alg_AV_koSRCa = CONSTANTS[AC_koCa] / alg_AV_kCaSR;
            RATES_SoA[SR_I * N + cellI] = alg_AV_kiSRCa * state_Ca_sub * state_SR_O - CONSTANTS[AC_kim] * state_SR_I - (CONSTANTS[AC_kom] * state_SR_I - alg_AV_koSRCa * pow(state_Ca_sub, 2.0) * state_SR_RI);
            RATES_SoA[SR_O * N + cellI] = alg_AV_koSRCa * pow(state_Ca_sub, 2.0) * state_SR_R - CONSTANTS[AC_kom] * state_SR_O - (alg_AV_kiSRCa * state_Ca_sub * state_SR_O - CONSTANTS[AC_kim] * state_SR_I);
            RATES_SoA[SR_R * N + cellI] = CONSTANTS[AC_kim] * state_SR_RI - alg_AV_kiSRCa * state_Ca_sub * state_SR_R - (alg_AV_koSRCa * pow(state_Ca_sub, 2.0) * state_SR_R - CONSTANTS[AC_kom] * state_SR_O);
            RATES_SoA[SR_RI * N + cellI] = CONSTANTS[AC_kom] * state_SR_I - alg_AV_koSRCa * pow(state_Ca_sub,
                                                              2.0) * state_SR_RI - (CONSTANTS[AC_kim] * state_SR_RI - alg_AV_kiSRCa * state_Ca_sub * state_SR_R);

            /* Ca_buffering */
            alg_AV_delta_fCMi = CONSTANTS[AC_kf_CM] * state_Ca_i * (1.0 - state_Buffer_fCMi) - CONSTANTS[AC_kb_CM] * state_Buffer_fCMi;
            alg_AV_delta_fCMs = CONSTANTS[AC_kf_CM] * state_Ca_sub * (1.0 - state_Buffer_fCMs) - CONSTANTS[AC_kb_CM] * state_Buffer_fCMs;
            alg_AV_delta_fCQ = CONSTANTS[AC_kf_CQ] * state_Ca_jsr * (1.0 - state_Buffer_fCQ) - CONSTANTS[AC_kb_CQ] * state_Buffer_fCQ;
            alg_AV_delta_fTC = CONSTANTS[AC_kf_TC] * state_Ca_i * (1.0 - state_Buffer_fTC) - CONSTANTS[AC_kb_TC] * state_Buffer_fTC;
            alg_AV_delta_fTMC = CONSTANTS[AC_kf_TMC] * state_Ca_i * (1.0 - (state_Buffer_fTMC + state_Buffer_fTMM)) - CONSTANTS[AC_kb_TMC] * state_Buffer_fTMC;
            alg_AV_delta_fTMM = CONSTANTS[AC_kf_TMM] * CONSTANTS[AC_Mgi] * (1.0 - (state_Buffer_fTMC + state_Buffer_fTMM)) - CONSTANTS[AC_kb_TMM] * state_Buffer_fTMM;
            RATES_SoA[Buffer_fCMi * N + cellI] = alg_AV_delta_fCMi;
            RATES_SoA[Buffer_fCMs * N + cellI] = alg_AV_delta_fCMs;
            RATES_SoA[Buffer_fCQ * N + cellI] = alg_AV_delta_fCQ;
            RATES_SoA[Buffer_fTC * N + cellI] = alg_AV_delta_fTC;
            RATES_SoA[Buffer_fTMC * N + cellI] = alg_AV_delta_fTMC;
            RATES_SoA[Buffer_fTMM * N + cellI] = alg_AV_delta_fTMM;



            /* i_CaL_fCa_gate */
            alg_AV_fCa_infinity = CONSTANTS[AC_Km_fCa] / (CONSTANTS[AC_Km_fCa] + state_Ca_sub);
            alg_AV_tau_fCa = alg_AV_fCa_infinity / CONSTANTS[AC_alpha_fCa];
            RATES_SoA[ICaL_fCa_gate_fCa * N + cellI] = (alg_AV_fCa_infinity - state_ICaL_fCa_gate_fCa) / alg_AV_tau_fCa;

            /* Ca_intracellular_fluxes */
            alg_AV_j_Ca_dif = (state_Ca_sub - state_Ca_i) / CONSTANTS[AC_tau_dif_Ca];
            alg_AV_j_tr = (state_Ca_nsr - state_Ca_jsr) / CONSTANTS[AC_tau_tr];
            alg_AV_j_up = CONSTANTS[AC_P_up] / (1.0 + exp((-state_Ca_i + CONSTANTS[AC_K_up]) / CONSTANTS[AC_slope_up]));

            /* Voltage_clamp */
            alg_AV_V_clamp = (((VOI > CONSTANTS[AC_t_holding]) && (VOI < CONSTANTS[AC_t_holding] + CONSTANTS[AC_t_test])) ? CONSTANTS[AC_V_test] : CONSTANTS[AC_V_holding]);

            /* Ca_dynamics */
            RATES_SoA[Ca_jsr * N + cellI] = alg_AV_j_tr - (alg_AV_j_SRCarel + CONSTANTS[AC_CQ_tot] * alg_AV_delta_fCQ);
            RATES_SoA[Ca_nsr * N + cellI] = alg_AV_j_up - alg_AV_j_tr * CONSTANTS[AC_V_jsr] / CONSTANTS[AC_V_nsr];
            RATES_SoA[Ca_i * N + cellI] = 1.0 * (alg_AV_j_Ca_dif * CONSTANTS[AC_V_sub] - alg_AV_j_up * CONSTANTS[AC_V_nsr]) / CONSTANTS[AC_V_i] - (CONSTANTS[AC_CM_tot] * alg_AV_delta_fCMi + CONSTANTS[AC_TC_tot] * alg_AV_delta_fTC + CONSTANTS[AC_TMC_tot] * alg_AV_delta_fTMC);

            /* Membrane */
            alg_AV_V = ((CONSTANTS[AC_clamp_mode] >= 1.0) ? alg_AV_V_clamp : state_membrane_V);
            if (alg_AV_V >  500.0) alg_AV_V =  500.0;
            if (alg_AV_V < -500.0) alg_AV_V = -500.0;

            /* Nai_concentration */
            alg_AV_Nai = state_Na_i;

            /* i_CaL */
            {
                const double _V = alg_AV_V;
                const double _R = CONSTANTS[AC_RTONF];
                const double _g2 = (fabs(_V) < 1e-7) ? (_R * 0.5) : (_V / (1.0 - exp(-2.0 * _V / _R)));
                const double _g1 = (fabs(_V) < 1e-7) ? _R         : (_V / (1.0 - exp(-_V / _R)));
                alg_AV_i_siCa = 2.0 * CONSTANTS[AC_P_CaL] * _g2 / _R * (state_Ca_sub - CONSTANTS[AC_Cao] * exp(-2.0 * _V / _R)) * state_ICaL_dL_gate_dL * state_ICaL_fL_gate_fL * state_ICaL_fCa_gate_fCa;
                alg_AV_i_siK  = 0.000365 * CONSTANTS[AC_P_CaL] * _g1 / _R * (CONSTANTS[AC_Ki] - CONSTANTS[AC_Ko] * exp(-_V / _R)) * state_ICaL_dL_gate_dL * state_ICaL_fL_gate_fL * state_ICaL_fCa_gate_fCa;
                alg_AV_i_siNa = 1.85e-05 * CONSTANTS[AC_P_CaL] * _g1 / _R * (alg_AV_Nai - CONSTANTS[AC_Nao] * exp(-_V / _R)) * state_ICaL_dL_gate_dL * state_ICaL_fL_gate_fL * state_ICaL_fCa_gate_fCa;
            }
            alg_AV_i_CaL_i_CaL = (alg_AV_i_siCa + alg_AV_i_siK + alg_AV_i_siNa) * (1.0 - CONSTANTS[AC_ACh_block]) * 1.0 * CONSTANTS[AC_i_CaL_Iso_increase];

            /* i_CaL_dL_gate */
            alg_AV_adVm = ((alg_AV_V == -41.8) ? -41.80001 : ((alg_AV_V == 0.0) ? 0.0 : ((alg_AV_V == -6.8) ? -6.80001 : alg_AV_V)));
            alg_AV_bdVm = ((alg_AV_V == -1.8) ? -1.80001 : alg_AV_V);
            alg_AV_alpha_dL = -0.02839 * (alg_AV_adVm + 41.8) / (exp(-(alg_AV_adVm + 41.8) / 2.5) - 1.0) - 0.0849 * (alg_AV_adVm + 6.8) / (exp(-(alg_AV_adVm + 6.8) / 4.8) - 1.0);
            alg_AV_beta_dL = 0.01143 * (alg_AV_bdVm + 1.8) / (exp((alg_AV_bdVm + 1.8) / 2.5) - 1.0);
            alg_AV_dL_infinity = 1.0 / (1.0 + exp(-(alg_AV_V - CONSTANTS[AC_V_dL] - CONSTANTS[AC_Iso_shift_dL]) / (CONSTANTS[AC_k_dL] * (1.0 + CONSTANTS[AC_Iso_slope_dL] / 100.0))));
            alg_AV_tau_dL = 1.0 / (alg_AV_alpha_dL + alg_AV_beta_dL);
            RATES_SoA[ICaL_dL_gate_dL * N + cellI] = (alg_AV_dL_infinity - state_ICaL_dL_gate_dL) / alg_AV_tau_dL;

            /* i_CaL_fL_gate */
            alg_AV_tau_fL = (44.3 + 230.0 * exp(-pow((alg_AV_V + 36.0) / 10.0, 2.0)));
            alg_AV_fL_infinity = 1.0 / (1.0 + exp((alg_AV_V + 37.4 + CONSTANTS[AC_shift_fL]) / (5.3 + CONSTANTS[AC_k_fL])));
            RATES_SoA[ICaL_fL_gate_fL * N + cellI] = (alg_AV_fL_infinity - state_ICaL_fL_gate_fL) / alg_AV_tau_fL;

            /* i_CaT */
            {
                const double _V = alg_AV_V;
                const double _R = CONSTANTS[AC_RTONF];
                const double _g2 = (fabs(_V) < 1e-7) ? (_R * 0.5) : (_V / (1.0 - exp(-2.0 * _V / _R)));
                alg_AV_i_CaT_i_CaT = 2.0 * CONSTANTS[AC_P_CaT] * _g2 / _R * (state_Ca_sub - CONSTANTS[AC_Cao] * exp(-2.0 * _V / _R)) * state_ICaT_dT_gate_dT * state_ICaT_fT_gate_fT;
            }

            /* i_CaT_dT_gate */
            alg_AV_dT_infinity = 1.0 / (1.0 + exp(-(alg_AV_V + 38.3) / 5.5));
            alg_AV_tau_dT = 1.0 / (1.068 * exp((alg_AV_V + 38.3) / 30.0) + 1.068 * exp(-(alg_AV_V + 38.3) / 30.0));
            RATES_SoA[ICaT_dT_gate_dT * N + cellI] = (alg_AV_dT_infinity - state_ICaT_dT_gate_dT) / alg_AV_tau_dT;

            /* i_CaT_fT_gate */
            alg_AV_fT_infinity = 1.0 / (1.0 + exp((alg_AV_V + 58.7) / 3.8));
            alg_AV_tau_fT = 1000.0 / (16.67 * exp(-(alg_AV_V + 75.0) / 83.3) + 16.67 * exp((alg_AV_V + 75.0) / 15.38)) + CONSTANTS[AC_offset_fT];
            RATES_SoA[ICaT_fT_gate_fT * N + cellI] = (alg_AV_fT_infinity - state_ICaT_fT_gate_fT) / alg_AV_tau_fT;

            /* i_KACh_a_gate */
            alg_AV_beta_a = 10.0 * exp(0.0133 * (alg_AV_V + 40.0));
            alg_AV_a_infinity = CONSTANTS[AC_alpha_a] / (CONSTANTS[AC_alpha_a] + alg_AV_beta_a);
            alg_AV_tau_a = 1000.0 / (CONSTANTS[AC_alpha_a] + alg_AV_beta_a);
            RATES_SoA[IKACh_a_gate_a * N + cellI] = (alg_AV_a_infinity - state_IKACh_a_gate_a) / alg_AV_tau_a;

            /* i_Kr_pa_gate */
            alg_AV_pa_infinity = 1.0 / (1.0 + exp(-(alg_AV_V + 10.0144) / 7.6607));
            alg_AV_tau_paF = 1000.0 / (30.0 * exp(alg_AV_V / 10.0) + exp(-alg_AV_V / 12.0));
            alg_AV_tau_paS = 8.46553540000000049e+02 / (4.2 * exp(alg_AV_V / 17.0) + 0.15 * exp(-alg_AV_V / 21.6));
            RATES_SoA[IKr_pa_gate_paF * N + cellI] = (alg_AV_pa_infinity - state_IKr_pa_gate_paF) / alg_AV_tau_paF;
            RATES_SoA[IKr_pa_gate_paS * N + cellI] = (alg_AV_pa_infinity - state_IKr_pa_gate_paS) / alg_AV_tau_paS;

            /* i_Kr_pi_gate */
            alg_AV_pi_infinity = 1.0 / (1.0 + exp((alg_AV_V + 28.6) / 17.1));
            alg_AV_tau_pi = 1000.0 / (100.0 * exp(-alg_AV_V / 54.645) + 656.0 * exp(alg_AV_V / 106.157));
            RATES_SoA[IKr_pi_gate_piy * N + cellI] = (alg_AV_pi_infinity - state_IKr_pi_gate_piy) / alg_AV_tau_pi;

            /* i_Ks */
            alg_AV_E_Ks = CONSTANTS[AC_RTONF] * log((CONSTANTS[AC_Ko] + 0.12 * CONSTANTS[AC_Nao]) / (CONSTANTS[AC_Ki] + 0.12 * alg_AV_Nai));
            alg_AV_i_Ks_i_Ks = CONSTANTS[AC_g_Ks] * (alg_AV_V - alg_AV_E_Ks) * pow(state_IKs_n_gate_n, 2.0);

            /* i_Ks_n_gate */
            alg_AV_alpha_n = 28.0 / (1.0 + exp(-(alg_AV_V - 40.0 - CONSTANTS[AC_i_Ks_n_gate_Iso_shift]) / 3.0));
            alg_AV_beta_n = 1.0 * exp(-(alg_AV_V - CONSTANTS[AC_i_Ks_n_gate_Iso_shift] - 5.0) / 25.0);
            alg_AV_n_infinity = sqrt(1.0 / (1.0 + exp(-(alg_AV_V + 0.6383 - CONSTANTS[AC_i_Ks_n_gate_Iso_shift]) / 10.7071)));
            alg_AV_tau_n = 1000.0 / (alg_AV_alpha_n + alg_AV_beta_n);
            RATES_SoA[IKs_n_gate_n * N + cellI] = (alg_AV_n_infinity - state_IKs_n_gate_n) / alg_AV_tau_n;

            /* i_Kur_rKur_gate */
            alg_AV_r_Kur_infinity = 1.0 / (1.0 + exp((alg_AV_V + 6.0) / -8.6));
            alg_AV_tau_r_Kur = 9.0 / (1.0 + exp((alg_AV_V + 5.0) / 12.0)) + 0.5;
            RATES_SoA[IKur_rKur_gate_r_Kur * N + cellI] = (alg_AV_r_Kur_infinity - state_IKur_rKur_gate_r_Kur) / alg_AV_tau_r_Kur;

            /* i_Kur_sKur_gate */
            alg_AV_s_Kur_infinity = 1.0 / (1.0 + exp((alg_AV_V + 7.5) / 10.0));
            alg_AV_tau_s_Kur = 590.0 / (1.0 + exp((alg_AV_V + 60.0) / 10.0)) + 3050.0;
            RATES_SoA[IKur_sKur_gate_s_Kur * N + cellI] = (alg_AV_s_Kur_infinity - state_IKur_sKur_gate_s_Kur) / alg_AV_tau_s_Kur;

            /* i_Na */
            alg_AV_E_mh = CONSTANTS[AC_RTONF] * log((CONSTANTS[AC_Nao] + 0.12 * CONSTANTS[AC_Ko]) / (alg_AV_Nai + 0.12 * CONSTANTS[AC_Ki]));
            alg_AV_i_Na_ = CONSTANTS[AC_g_Na] * pow(state_INa_m_gate_m, 3.0) * state_INa_h_gate_h * (alg_AV_V - alg_AV_E_mh);
            alg_AV_i_Na_L = CONSTANTS[AC_g_Na_L] * pow(state_INa_m_gate_m, 3.0) * (alg_AV_V - alg_AV_E_mh);
            alg_AV_i_Na_i_Na = alg_AV_i_Na_ + alg_AV_i_Na_L;

            /* i_NaCa */
            alg_AV_di = 1.0 + state_Ca_sub / CONSTANTS[AC_Kci] * (1.0 + exp(-CONSTANTS[AC_Qci] * alg_AV_V / CONSTANTS[AC_RTONF]) + alg_AV_Nai / CONSTANTS[AC_Kcni]) + alg_AV_Nai / CONSTANTS[AC_K1ni] * (1.0 + alg_AV_Nai / CONSTANTS[AC_K2ni] * (1.0 + alg_AV_Nai / CONSTANTS[AC_K3ni]));
            alg_AV_i_NaCa_do = 1.0 + CONSTANTS[AC_Cao] / CONSTANTS[AC_Kco] * (1.0 + exp(CONSTANTS[AC_Qco] * alg_AV_V / CONSTANTS[AC_RTONF])) + CONSTANTS[AC_Nao] / CONSTANTS[AC_K1no] * (1.0 + CONSTANTS[AC_Nao] / CONSTANTS[AC_K2no] * (1.0 + CONSTANTS[AC_Nao] / CONSTANTS[AC_K3no]));
            alg_AV_k32 = exp(CONSTANTS[AC_Qn] * alg_AV_V / (2.0 * CONSTANTS[AC_RTONF]));
            alg_AV_k41 = exp(-CONSTANTS[AC_Qn] * alg_AV_V / (2.0 * CONSTANTS[AC_RTONF]));
            alg_AV_k43 = alg_AV_Nai / (CONSTANTS[AC_K3ni] + alg_AV_Nai);
            alg_AV_k12 = state_Ca_sub / CONSTANTS[AC_Kci] * exp(-CONSTANTS[AC_Qci] * alg_AV_V / CONSTANTS[AC_RTONF]) / alg_AV_di;
            alg_AV_k14 = alg_AV_Nai / CONSTANTS[AC_K1ni] * alg_AV_Nai / CONSTANTS[AC_K2ni] * (1.0 + alg_AV_Nai / CONSTANTS[AC_K3ni]) * exp(CONSTANTS[AC_Qn] * alg_AV_V / (2.0 * CONSTANTS[AC_RTONF])) / alg_AV_di;
            alg_AV_k21 = CONSTANTS[AC_Cao] / CONSTANTS[AC_Kco] * exp(CONSTANTS[AC_Qco] * alg_AV_V / CONSTANTS[AC_RTONF]) / alg_AV_i_NaCa_do;
            alg_AV_k23 = CONSTANTS[AC_Nao] / CONSTANTS[AC_K1no] * CONSTANTS[AC_Nao] / CONSTANTS[AC_K2no] * (1.0 + CONSTANTS[AC_Nao] / CONSTANTS[AC_K3no]) * exp(-CONSTANTS[AC_Qn] * alg_AV_V / (2.0 * CONSTANTS[AC_RTONF])) / alg_AV_i_NaCa_do;
            alg_AV_x1 = alg_AV_k41 * CONSTANTS[AC_k34] * (alg_AV_k23 + alg_AV_k21) + alg_AV_k21 * alg_AV_k32 * (alg_AV_k43 + alg_AV_k41);
            alg_AV_x2 = alg_AV_k32 * alg_AV_k43 * (alg_AV_k14 + alg_AV_k12) + alg_AV_k41 * alg_AV_k12 * (CONSTANTS[AC_k34] + alg_AV_k32);
            alg_AV_x3 = alg_AV_k14 * alg_AV_k43 * (alg_AV_k23 + alg_AV_k21) + alg_AV_k12 * alg_AV_k23 * (alg_AV_k43 + alg_AV_k41);
            alg_AV_x4 = alg_AV_k23 * CONSTANTS[AC_k34] * (alg_AV_k14 + alg_AV_k12) + alg_AV_k14 * alg_AV_k21 * (CONSTANTS[AC_k34] + alg_AV_k32);
            alg_AV_i_NaCa_i_NaCa = (1.0 - CONSTANTS[AC_blockade_NaCa]) * CONSTANTS[AC_K_NaCa] * (alg_AV_x2 * alg_AV_k21 - alg_AV_x1 * alg_AV_k12) / (alg_AV_x1 + alg_AV_x2 + alg_AV_x3 + alg_AV_x4);

            /* i_Na_h_gate */
            alg_AV_alpha_h = 20.0 * exp(-0.125 * (alg_AV_V + 75.0));
            alg_AV_beta_h = 2000.0 / (320.0 * exp(-0.1 * (alg_AV_V + 75.0)) + 1.0);
            alg_AV_h_infinity = 1.0 / (1.0 + exp((alg_AV_V + 69.804) / 4.4565));
            alg_AV_tau_h = 1000.0 / (alg_AV_alpha_h + alg_AV_beta_h);
            RATES_SoA[INa_h_gate_h * N + cellI] = (alg_AV_h_infinity - state_INa_h_gate_h) / alg_AV_tau_h;

            /* i_Na_m_gate */
            alg_AV_E0_m = alg_AV_V + 41.0;
            alg_AV_beta_m = 8000.0 * exp(-0.056 * (alg_AV_V + 66.0));
            alg_AV_m_infinity = 1.0 / (1.0 + exp(-(alg_AV_V + 42.0504) / 8.3106));
            alg_AV_alpha_m = ((fabs(alg_AV_E0_m) < CONSTANTS[AC_delta_m]) ? 2000.0 : 200.0 * alg_AV_E0_m / (1.0 - exp(-0.1 * alg_AV_E0_m)));
            alg_AV_tau_m = 1000.0 / (alg_AV_alpha_m + alg_AV_beta_m);
            RATES_SoA[INa_m_gate_m * N + cellI] = (alg_AV_m_infinity - state_INa_m_gate_m) / alg_AV_tau_m;

            /* i_f_y_gate */
            {
                const double _x1 = alg_AV_V + 148.8 - CONSTANTS[AC_ACh_shift] - CONSTANTS[AC_i_f_y_gate_Iso_shift];
                const double _x2 = alg_AV_V + 87.3  - CONSTANTS[AC_ACh_shift] - CONSTANTS[AC_i_f_y_gate_Iso_shift];
                const double _t1 = (fabs(_x1) < 1e-7)
                    ? (0.36 / 0.066)
                    : (0.36 * _x1 / (exp(0.066 * _x1) - 1.0));
                const double _t2 = (fabs(_x2) < 1e-7)
                    ? 0.5
                    : ((-0.2 * _x2 > 700.0) ? 0.0
                       : (0.1 * _x2 / (1.0 - exp(-0.2 * _x2))));
                const double _raw_s = 1.0 / (_t1 + _t2) - 0.054;
                const double _raw = _raw_s * 1000.0;
                alg_AV_tau_y = (_raw > 1.0) ? _raw : 1.0;
            }
            alg_AV_y_infinity = ((alg_AV_V < -(80.0 - CONSTANTS[AC_ACh_shift] - CONSTANTS[AC_i_f_y_gate_Iso_shift] - CONSTANTS[AC_y_shift])) ? 0.01329 + 0.99921 / (1.0 + exp((alg_AV_V + 97.134 - CONSTANTS[AC_ACh_shift] - CONSTANTS[AC_i_f_y_gate_Iso_shift] - CONSTANTS[AC_y_shift]) / 8.1752)) : 0.0002501 * exp(-(alg_AV_V - CONSTANTS[AC_ACh_shift] - CONSTANTS[AC_i_f_y_gate_Iso_shift] - CONSTANTS[AC_y_shift]) / 12.861));
            RATES_SoA[If_y_gate_y * N + cellI] = (alg_AV_y_infinity - state_If_y_gate_y) / alg_AV_tau_y;

            /* i_to_q_gate */
            alg_AV_q_infinity = 1.0 / (1.0 + exp((alg_AV_V + 49.0) / 13.0));
            alg_AV_tau_q = 0.6 * (65.17 / (0.57 * exp(-0.08 * (alg_AV_V + 44.0)) + 0.065 * exp(0.1 * (alg_AV_V + 45.93))) + 10.1);
            RATES_SoA[Ito_q_gate_q * N + cellI] = (alg_AV_q_infinity - state_Ito_q_gate_q) / alg_AV_tau_q;

            /* i_to_r_gate */
            alg_AV_r_infinity = 1.0 / (1.0 + exp(-(alg_AV_V - 19.3) / 15.0));
            alg_AV_tau_r = 0.66 * 1.4 * (15.59 / (1.037 * exp(0.09 * (alg_AV_V + 30.61)) + 0.369 * exp(-0.12 * (alg_AV_V + 23.84))) + 2.98);
            RATES_SoA[Ito_r_gate_r * N + cellI] = (alg_AV_r_infinity - state_Ito_r_gate_r) / alg_AV_tau_r;

            /* *remaining* */
            RATES_SoA[Ca_sub * N + cellI] = alg_AV_j_SRCarel * CONSTANTS[AC_V_jsr] / CONSTANTS[AC_V_sub] - ((alg_AV_i_siCa + alg_AV_i_CaT_i_CaT - 2.0 * alg_AV_i_NaCa_i_NaCa) / (2.0 * CONSTANTS[AC_F] * CONSTANTS[AC_V_sub]) + alg_AV_j_Ca_dif + CONSTANTS[AC_CM_tot] * alg_AV_delta_fCMs);
            alg_AV_E_Na = CONSTANTS[AC_RTONF] * log(CONSTANTS[AC_Nao] / alg_AV_Nai);
            alg_AV_i_KACh_i_KACh = ((CONSTANTS[AC_ACh] > 0.0) ? CONSTANTS[AC_ACh_on] * CONSTANTS[AC_g_KACh] * (alg_AV_V - CONSTANTS[AC_E_K]) * (1.0 + exp((alg_AV_V + 20.0) / 20.0)) * state_IKACh_a_gate_a : 0.0);
            alg_AV_i_Kr_i_Kr = CONSTANTS[AC_g_Kr] * (alg_AV_V - CONSTANTS[AC_E_K]) * (0.9 * state_IKr_pa_gate_paF + 0.1 * state_IKr_pa_gate_paS) * state_IKr_pi_gate_piy;
            alg_AV_i_Kur_i_Kur = CONSTANTS[AC_g_Kur] * state_IKur_rKur_gate_r_Kur * state_IKur_sKur_gate_s_Kur * (alg_AV_V - CONSTANTS[AC_E_K]);
            alg_AV_i_NaK_i_NaK = CONSTANTS[AC_i_NaK_Iso_increase] * CONSTANTS[AC_i_NaK_max] * pow(1.0 + pow(CONSTANTS[AC_Km_Kp] / CONSTANTS[AC_Ko], 1.2),
                                                                                              -1.0) * pow(1.0 + pow(CONSTANTS[AC_Km_Nap] / alg_AV_Nai, 1.3),
                                                                                                          -1.0) * pow(1.0 + exp(-(alg_AV_V - alg_AV_E_Na + 110.0) / 20.0),
                                                                                                                      -1.0);
            alg_AV_i_fK = state_If_y_gate_y * CONSTANTS[AC_g_f_K] * (alg_AV_V - CONSTANTS[AC_E_K]) * (1.0 - CONSTANTS[AC_blockade]);
            alg_AV_i_fNa = state_If_y_gate_y * CONSTANTS[AC_g_f_Na] * (alg_AV_V - alg_AV_E_Na) * (1.0 - CONSTANTS[AC_blockade]);
            alg_AV_i_to_i_to = CONSTANTS[AC_g_to] * (alg_AV_V - CONSTANTS[AC_E_K]) * state_Ito_q_gate_q * state_Ito_r_gate_r;
            RATES_SoA[Na_i * N + cellI] = (1.0 - CONSTANTS[AC_Nai_clamp]) * -1.0 * (alg_AV_i_Na_i_Na + alg_AV_i_fNa + alg_AV_i_siNa + 3.0 * alg_AV_i_NaK_i_NaK + 3.0 * alg_AV_i_NaCa_i_NaCa) / (1.0 * (CONSTANTS[AC_V_i] + CONSTANTS[AC_V_sub]) * CONSTANTS[AC_F]);
            alg_AV_i_f_i_f = alg_AV_i_fNa + alg_AV_i_fK;
            alg_Iion_cm = alg_AV_i_f_i_f + alg_AV_i_Kr_i_Kr + alg_AV_i_Ks_i_Ks + alg_AV_i_to_i_to + alg_AV_i_NaK_i_NaK + alg_AV_i_NaCa_i_NaCa + alg_AV_i_Na_i_Na + alg_AV_i_CaL_i_CaL + alg_AV_i_CaT_i_CaT + alg_AV_i_KACh_i_KACh + alg_AV_i_Kur_i_Kur;

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
            SUPPORT_SoA[FABBRI_BATCH_SUPPORT_tau_y * N + cellI] = alg_AV_tau_y;
            SUPPORT_SoA[FABBRI_BATCH_SUPPORT_tau_m * N + cellI] = alg_AV_tau_m;
            SUPPORT_SoA[FABBRI_BATCH_SUPPORT_tau_h * N + cellI] = alg_AV_tau_h;
            SUPPORT_SoA[FABBRI_BATCH_SUPPORT_tau_dL * N + cellI] = alg_AV_tau_dL;
            SUPPORT_SoA[FABBRI_BATCH_SUPPORT_tau_fL * N + cellI] = alg_AV_tau_fL;
            SUPPORT_SoA[FABBRI_BATCH_SUPPORT_tau_fCa * N + cellI] = alg_AV_tau_fCa;
            SUPPORT_SoA[FABBRI_BATCH_SUPPORT_tau_dT * N + cellI] = alg_AV_tau_dT;
            SUPPORT_SoA[FABBRI_BATCH_SUPPORT_tau_fT * N + cellI] = alg_AV_tau_fT;
            SUPPORT_SoA[FABBRI_BATCH_SUPPORT_tau_r_Kur * N + cellI] = alg_AV_tau_r_Kur;
            SUPPORT_SoA[FABBRI_BATCH_SUPPORT_tau_s_Kur * N + cellI] = alg_AV_tau_s_Kur;
            SUPPORT_SoA[FABBRI_BATCH_SUPPORT_tau_q * N + cellI] = alg_AV_tau_q;
            SUPPORT_SoA[FABBRI_BATCH_SUPPORT_tau_r * N + cellI] = alg_AV_tau_r;
            SUPPORT_SoA[FABBRI_BATCH_SUPPORT_tau_paS * N + cellI] = alg_AV_tau_paS;
            SUPPORT_SoA[FABBRI_BATCH_SUPPORT_tau_paF * N + cellI] = alg_AV_tau_paF;
            SUPPORT_SoA[FABBRI_BATCH_SUPPORT_tau_pi * N + cellI] = alg_AV_tau_pi;
            SUPPORT_SoA[FABBRI_BATCH_SUPPORT_tau_n * N + cellI] = alg_AV_tau_n;
            SUPPORT_SoA[FABBRI_BATCH_SUPPORT_tau_a * N + cellI] = alg_AV_tau_a;
            SUPPORT_SoA[FABBRI_BATCH_SUPPORT_gInf_y * N + cellI] = alg_AV_y_infinity;
            SUPPORT_SoA[FABBRI_BATCH_SUPPORT_gInf_m * N + cellI] = alg_AV_m_infinity;
            SUPPORT_SoA[FABBRI_BATCH_SUPPORT_gInf_h * N + cellI] = alg_AV_h_infinity;
            SUPPORT_SoA[FABBRI_BATCH_SUPPORT_gInf_dL * N + cellI] = alg_AV_dL_infinity;
            SUPPORT_SoA[FABBRI_BATCH_SUPPORT_gInf_fL * N + cellI] = alg_AV_fL_infinity;
            SUPPORT_SoA[FABBRI_BATCH_SUPPORT_gInf_fCa * N + cellI] = alg_AV_fCa_infinity;
            SUPPORT_SoA[FABBRI_BATCH_SUPPORT_gInf_dT * N + cellI] = alg_AV_dT_infinity;
            SUPPORT_SoA[FABBRI_BATCH_SUPPORT_gInf_fT * N + cellI] = alg_AV_fT_infinity;
            SUPPORT_SoA[FABBRI_BATCH_SUPPORT_gInf_r_Kur * N + cellI] = alg_AV_r_Kur_infinity;
            SUPPORT_SoA[FABBRI_BATCH_SUPPORT_gInf_s_Kur * N + cellI] = alg_AV_s_Kur_infinity;
            SUPPORT_SoA[FABBRI_BATCH_SUPPORT_gInf_q * N + cellI] = alg_AV_q_infinity;
            SUPPORT_SoA[FABBRI_BATCH_SUPPORT_gInf_r * N + cellI] = alg_AV_r_infinity;
            SUPPORT_SoA[FABBRI_BATCH_SUPPORT_gInf_paS * N + cellI] = alg_AV_pa_infinity;
            SUPPORT_SoA[FABBRI_BATCH_SUPPORT_gInf_paF * N + cellI] = alg_AV_pa_infinity;
            SUPPORT_SoA[FABBRI_BATCH_SUPPORT_gInf_pi * N + cellI] = alg_AV_pi_infinity;
            SUPPORT_SoA[FABBRI_BATCH_SUPPORT_gInf_n * N + cellI] = alg_AV_n_infinity;
            SUPPORT_SoA[FABBRI_BATCH_SUPPORT_gInf_a * N + cellI] = alg_AV_a_infinity;
            SUPPORT_SoA[FABBRI_BATCH_SUPPORT_Iion_cm * N + cellI] = alg_Iion_cm;
        }
    }
}

} // End namespace Foam

// ************************************************************************* //
