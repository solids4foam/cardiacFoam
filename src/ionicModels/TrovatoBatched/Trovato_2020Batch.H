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
    State-major (SoA) batched evaluator for the Trovato 2020 ionic model.
Author
    Simao Nieto de Castro, UCD.
\*---------------------------------------------------------------------------*/

#pragma once

#include "Trovato_2020Names.H"
#include "StimulusProtocolPOD.H"
#include <cmath>

namespace Foam
{

enum TROVATO_BATCH_SUPPORT_INDEX
{
    TROVATO_BATCH_SUPPORT_tau_ICaT_b,
    TROVATO_BATCH_SUPPORT_gInf_ICaT_b,
    TROVATO_BATCH_SUPPORT_tau_ICaT_g,
    TROVATO_BATCH_SUPPORT_gInf_ICaT_g,
    TROVATO_BATCH_SUPPORT_tau_IK1_xk1,
    TROVATO_BATCH_SUPPORT_gInf_IK1_xk1,
    TROVATO_BATCH_SUPPORT_tau_IKr_xrf,
    TROVATO_BATCH_SUPPORT_gInf_IKr_xrf,
    TROVATO_BATCH_SUPPORT_tau_IKr_xrs,
    TROVATO_BATCH_SUPPORT_gInf_IKr_xrs,
    TROVATO_BATCH_SUPPORT_tau_IKs_xs1,
    TROVATO_BATCH_SUPPORT_gInf_IKs_xs1,
    TROVATO_BATCH_SUPPORT_tau_IKs_xs2,
    TROVATO_BATCH_SUPPORT_gInf_IKs_xs2,
    TROVATO_BATCH_SUPPORT_tau_INa_hf,
    TROVATO_BATCH_SUPPORT_gInf_INa_hf,
    TROVATO_BATCH_SUPPORT_tau_INa_hs,
    TROVATO_BATCH_SUPPORT_gInf_INa_hs,
    TROVATO_BATCH_SUPPORT_tau_INa_hsp,
    TROVATO_BATCH_SUPPORT_gInf_INa_hsp,
    TROVATO_BATCH_SUPPORT_tau_INa_m,
    TROVATO_BATCH_SUPPORT_gInf_INa_m,
    TROVATO_BATCH_SUPPORT_tau_INa_j,
    TROVATO_BATCH_SUPPORT_gInf_INa_j,
    TROVATO_BATCH_SUPPORT_tau_INa_jp,
    TROVATO_BATCH_SUPPORT_gInf_INa_jp,
    TROVATO_BATCH_SUPPORT_tau_If_y,
    TROVATO_BATCH_SUPPORT_gInf_If_y,
    TROVATO_BATCH_SUPPORT_tau_Ito_a,
    TROVATO_BATCH_SUPPORT_gInf_Ito_a,
    TROVATO_BATCH_SUPPORT_tau_Ito_i1,
    TROVATO_BATCH_SUPPORT_gInf_Ito_i1,
    TROVATO_BATCH_SUPPORT_tau_Ito_i2,
    TROVATO_BATCH_SUPPORT_gInf_Ito_i2,
    TROVATO_BATCH_SUPPORT_tau_INaL_hL,
    TROVATO_BATCH_SUPPORT_gInf_INaL_hL,
    TROVATO_BATCH_SUPPORT_tau_INaL_mL,
    TROVATO_BATCH_SUPPORT_gInf_INaL_mL,
    TROVATO_BATCH_SUPPORT_tau_INaL_hLp,
    TROVATO_BATCH_SUPPORT_gInf_INaL_hLp,
    TROVATO_BATCH_SUPPORT_tau_ICaL_d,
    TROVATO_BATCH_SUPPORT_gInf_ICaL_d,
    TROVATO_BATCH_SUPPORT_tau_ICaL_ff,
    TROVATO_BATCH_SUPPORT_gInf_ICaL_ff,
    TROVATO_BATCH_SUPPORT_tau_ICaL_fs,
    TROVATO_BATCH_SUPPORT_gInf_ICaL_fs,
    TROVATO_BATCH_SUPPORT_tau_ICaL_fcaf,
    TROVATO_BATCH_SUPPORT_gInf_ICaL_fcaf,
    TROVATO_BATCH_SUPPORT_tau_ICaL_fcafp,
    TROVATO_BATCH_SUPPORT_gInf_ICaL_fcafp,
    TROVATO_BATCH_SUPPORT_tau_ICaL_fcas,
    TROVATO_BATCH_SUPPORT_gInf_ICaL_fcas,
    TROVATO_BATCH_SUPPORT_tau_ICaL_ffp,
    TROVATO_BATCH_SUPPORT_gInf_ICaL_ffp,
    TROVATO_BATCH_SUPPORT_tau_ICaL_jca,
    TROVATO_BATCH_SUPPORT_gInf_ICaL_jca,
    TROVATO_BATCH_SUPPORT_tau_ryr_Jrel1,
    TROVATO_BATCH_SUPPORT_gInf_ryr_Jrel1,
    TROVATO_BATCH_SUPPORT_tau_ryr_Jrel2,
    TROVATO_BATCH_SUPPORT_gInf_ryr_Jrel2,
    TROVATO_BATCH_SUPPORT_Iion_cm,

    NUM_TROVATO_BATCH_SUPPORT
};

static_assert(NUM_TROVATO_BATCH_SUPPORT == 61, "Trovato compact support layout drift");

CARDIAC_HOST_DEVICE inline void TrovatoComputeVariablesBatch
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
        const double state_CaMKt = STATES_SoA[CaMKt * N + cellI];
        const double state_Cacsr = STATES_SoA[Cacsr * N + cellI];
        const double state_Cai = STATES_SoA[Cai * N + cellI];
        const double state_Cajsr = STATES_SoA[Cajsr * N + cellI];
        const double state_Cansr = STATES_SoA[Cansr * N + cellI];
        const double state_Casl = STATES_SoA[Casl * N + cellI];
        const double state_Cass = STATES_SoA[Cass * N + cellI];
        const double state_ICaL_d = STATES_SoA[ICaL_d * N + cellI];
        const double state_ICaL_fcaf = STATES_SoA[ICaL_fcaf * N + cellI];
        const double state_ICaL_fcafp = STATES_SoA[ICaL_fcafp * N + cellI];
        const double state_ICaL_fcas = STATES_SoA[ICaL_fcas * N + cellI];
        const double state_ICaL_ff = STATES_SoA[ICaL_ff * N + cellI];
        const double state_ICaL_ffp = STATES_SoA[ICaL_ffp * N + cellI];
        const double state_ICaL_fs = STATES_SoA[ICaL_fs * N + cellI];
        const double state_ICaL_jca = STATES_SoA[ICaL_jca * N + cellI];
        const double state_ICaL_nca = STATES_SoA[ICaL_nca * N + cellI];
        const double state_ICaT_b = STATES_SoA[ICaT_b * N + cellI];
        const double state_ICaT_g = STATES_SoA[ICaT_g * N + cellI];
        const double state_IK1_xk1 = STATES_SoA[IK1_xk1 * N + cellI];
        const double state_IKr_xrf = STATES_SoA[IKr_xrf * N + cellI];
        const double state_IKr_xrs = STATES_SoA[IKr_xrs * N + cellI];
        const double state_IKs_xs1 = STATES_SoA[IKs_xs1 * N + cellI];
        const double state_IKs_xs2 = STATES_SoA[IKs_xs2 * N + cellI];
        const double state_INaL_hL = STATES_SoA[INaL_hL * N + cellI];
        const double state_INaL_hLp = STATES_SoA[INaL_hLp * N + cellI];
        const double state_INaL_mL = STATES_SoA[INaL_mL * N + cellI];
        const double state_INa_hf = STATES_SoA[INa_hf * N + cellI];
        const double state_INa_hs = STATES_SoA[INa_hs * N + cellI];
        const double state_INa_hsp = STATES_SoA[INa_hsp * N + cellI];
        const double state_INa_j = STATES_SoA[INa_j * N + cellI];
        const double state_INa_jp = STATES_SoA[INa_jp * N + cellI];
        const double state_INa_m = STATES_SoA[INa_m * N + cellI];
        const double state_IP3_u = STATES_SoA[IP3_u * N + cellI];
        const double state_If_y = STATES_SoA[If_y * N + cellI];
        const double state_Ito_a = STATES_SoA[Ito_a * N + cellI];
        const double state_Ito_i1 = STATES_SoA[Ito_i1 * N + cellI];
        const double state_Ito_i2 = STATES_SoA[Ito_i2 * N + cellI];
        const double state_Ki = STATES_SoA[Ki * N + cellI];
        const double state_Ksl = STATES_SoA[Ksl * N + cellI];
        const double state_Kss = STATES_SoA[Kss * N + cellI];
        const double state_Nai = STATES_SoA[Nai * N + cellI];
        const double state_Nasl = STATES_SoA[Nasl * N + cellI];
        const double state_Nass = STATES_SoA[Nass * N + cellI];
        const double state_membrane_v = STATES_SoA[membrane_v * N + cellI];
        const double state_ryr_Jrel1 = STATES_SoA[ryr_Jrel1 * N + cellI];
        const double state_ryr_Jrel2 = STATES_SoA[ryr_Jrel2 * N + cellI];

        double alg_AV_Afcaf = 0.0;
        double alg_AV_Afcas = 0.0;
        double alg_AV_Axrf = 0.0;
        double alg_AV_Axrs = 0.0;
        double alg_AV_Bcacsr = 0.0;
        double alg_AV_Bcai = 0.0;
        double alg_AV_Bcajsr = 0.0;
        double alg_AV_Bcasl = 0.0;
        double alg_AV_Bcass = 0.0;
        double alg_AV_CaMKa = 0.0;
        double alg_AV_CaMKb = 0.0;
        double alg_AV_E1 = 0.0;
        double alg_AV_E1_i = 0.0;
        double alg_AV_E1_ss = 0.0;
        double alg_AV_E2 = 0.0;
        double alg_AV_E2_i = 0.0;
        double alg_AV_E2_ss = 0.0;
        double alg_AV_E3 = 0.0;
        double alg_AV_E3_i = 0.0;
        double alg_AV_E3_ss = 0.0;
        double alg_AV_E4 = 0.0;
        double alg_AV_E4_i = 0.0;
        double alg_AV_E4_ss = 0.0;
        double alg_AV_ECa = 0.0;
        double alg_AV_EK = 0.0;
        double alg_AV_EKs = 0.0;
        double alg_AV_ENa = 0.0;
        double alg_AV_ICaK = 0.0;
        double alg_AV_ICaL_ICaL = 0.0;
        double alg_AV_ICaNa = 0.0;
        double alg_AV_ICaT_ICaT = 0.0;
        double alg_AV_ICab_ICab = 0.0;
        double alg_AV_IK1_IK1 = 0.0;
        double alg_AV_IKr_IKr = 0.0;
        double alg_AV_IKs_IKs = 0.0;
        double alg_AV_INaCa_i_INaCa_i = 0.0;
        double alg_AV_INaCa_ss = 0.0;
        double alg_AV_INaK_INaK = 0.0;
        double alg_AV_INaL_INaL = 0.0;
        double alg_AV_INa_INa = 0.0;
        double alg_AV_INab_INab = 0.0;
        double alg_AV_IfK = 0.0;
        double alg_AV_IfNa = 0.0;
        double alg_AV_If_If = 0.0;
        double alg_AV_IpCa_IpCa = 0.0;
        double alg_AV_Isus_Isus = 0.0;
        double alg_AV_Ito_Ito = 0.0;
        double alg_AV_Jdiff = 0.0;
        double alg_AV_JdiffK = 0.0;
        double alg_AV_JdiffNa = 0.0;
        double alg_AV_Jgap = 0.0;
        double alg_AV_JgapK = 0.0;
        double alg_AV_JgapNa = 0.0;
        double alg_AV_Jip3 = 0.0;
        double alg_AV_JnakK = 0.0;
        double alg_AV_JnakNa = 0.0;
        double alg_AV_JncxCa_i = 0.0;
        double alg_AV_JncxCa_ss = 0.0;
        double alg_AV_JncxNa_i = 0.0;
        double alg_AV_JncxNa_ss = 0.0;
        double alg_AV_Jtr1 = 0.0;
        double alg_AV_Jtr2 = 0.0;
        double alg_AV_Jup1 = 0.0;
        double alg_AV_Jup2 = 0.0;
        double alg_AV_Knai = 0.0;
        double alg_AV_Knao = 0.0;
        double alg_AV_KsCa = 0.0;
        double alg_AV_P = 0.0;
        double alg_AV_POip3 = 0.0;
        double alg_AV_PhiCaK = 0.0;
        double alg_AV_PhiCaL = 0.0;
        double alg_AV_PhiCaNa = 0.0;
        double alg_AV_REL = 0.0;
        double alg_AV_REL2 = 0.0;
        double alg_AV_a1 = 0.0;
        double alg_AV_a3 = 0.0;
        double alg_AV_allo_i = 0.0;
        double alg_AV_allo_ss = 0.0;
        double alg_AV_anca = 0.0;
        double alg_AV_ass = 0.0;
        double alg_AV_asus = 0.0;
        double alg_AV_b2 = 0.0;
        double alg_AV_b3 = 0.0;
        double alg_AV_b4 = 0.0;
        double alg_AV_bss = 0.0;
        double alg_AV_dkmplb = 0.0;
        double alg_AV_dqupcamk = 0.0;
        double alg_AV_dss = 0.0;
        double alg_AV_f = 0.0;
        double alg_AV_fICaLp = 0.0;
        double alg_AV_fINaLp = 0.0;
        double alg_AV_fINap = 0.0;
        double alg_AV_fca = 0.0;
        double alg_AV_fcap = 0.0;
        double alg_AV_fcass = 0.0;
        double alg_AV_fp = 0.0;
        double alg_AV_fss = 0.0;
        double alg_AV_gss = 0.0;
        double alg_AV_h = 0.0;
        double alg_AV_h1_i = 0.0;
        double alg_AV_h1_ss = 0.0;
        double alg_AV_h2_i = 0.0;
        double alg_AV_h2_ss = 0.0;
        double alg_AV_h3_i = 0.0;
        double alg_AV_h3_ss = 0.0;
        double alg_AV_h4_i = 0.0;
        double alg_AV_h4_ss = 0.0;
        double alg_AV_h5_i = 0.0;
        double alg_AV_h5_ss = 0.0;
        double alg_AV_h6_i = 0.0;
        double alg_AV_h6_ss = 0.0;
        double alg_AV_h7_i = 0.0;
        double alg_AV_h7_ss = 0.0;
        double alg_AV_h8_i = 0.0;
        double alg_AV_h8_ss = 0.0;
        double alg_AV_h9_i = 0.0;
        double alg_AV_h9_ss = 0.0;
        double alg_AV_hLss = 0.0;
        double alg_AV_hLssp = 0.0;
        double alg_AV_hca = 0.0;
        double alg_AV_hna = 0.0;
        double alg_AV_hp = 0.0;
        double alg_AV_hss = 0.0;
        double alg_AV_hssp = 0.0;
        double alg_AV_irelss = 0.0;
        double alg_AV_irelss2 = 0.0;
        double alg_AV_ireltau = 0.0;
        double alg_AV_ireltau2 = 0.0;
        double alg_AV_iss = 0.0;
        double alg_AV_jss = 0.0;
        double alg_AV_k3_i = 0.0;
        double alg_AV_k3_ss = 0.0;
        double alg_AV_k3p_i = 0.0;
        double alg_AV_k3p_ss = 0.0;
        double alg_AV_k3pp_i = 0.0;
        double alg_AV_k3pp_ss = 0.0;
        double alg_AV_k4_i = 0.0;
        double alg_AV_k4_ss = 0.0;
        double alg_AV_k4p_i = 0.0;
        double alg_AV_k4p_ss = 0.0;
        double alg_AV_k4pp_i = 0.0;
        double alg_AV_k4pp_ss = 0.0;
        double alg_AV_k6_i = 0.0;
        double alg_AV_k6_ss = 0.0;
        double alg_AV_k7_i = 0.0;
        double alg_AV_k7_ss = 0.0;
        double alg_AV_k8_i = 0.0;
        double alg_AV_k8_ss = 0.0;
        double alg_AV_km2n = 0.0;
        double alg_AV_mLss = 0.0;
        double alg_AV_mss = 0.0;
        double alg_AV_rk1 = 0.0;
        double alg_AV_rkr = 0.0;
        double alg_AV_taua = 0.0;
        double alg_AV_taub = 0.0;
        double alg_AV_taug = 0.0;
        double alg_AV_tauif = 0.0;
        double alg_AV_tauis = 0.0;
        double alg_AV_tauy = 0.0;
        double alg_AV_td = 0.0;
        double alg_AV_tfcaf = 0.0;
        double alg_AV_tfcafp = 0.0;
        double alg_AV_tfcas = 0.0;
        double alg_AV_tff = 0.0;
        double alg_AV_tffp = 0.0;
        double alg_AV_tfs = 0.0;
        double alg_AV_thf = 0.0;
        double alg_AV_ths = 0.0;
        double alg_AV_thsp = 0.0;
        double alg_AV_tj = 0.0;
        double alg_AV_tjp = 0.0;
        double alg_AV_tm = 0.0;
        double alg_AV_tmL = 0.0;
        double alg_AV_txk1 = 0.0;
        double alg_AV_txrf = 0.0;
        double alg_AV_txrs = 0.0;
        double alg_AV_txs1 = 0.0;
        double alg_AV_txs2 = 0.0;
        double alg_AV_vffrt = 0.0;
        double alg_AV_vfrt = 0.0;
        double alg_AV_x1 = 0.0;
        double alg_AV_x1_i = 0.0;
        double alg_AV_x1_ss = 0.0;
        double alg_AV_x2 = 0.0;
        double alg_AV_x2_i = 0.0;
        double alg_AV_x2_ss = 0.0;
        double alg_AV_x3 = 0.0;
        double alg_AV_x3_i = 0.0;
        double alg_AV_x3_ss = 0.0;
        double alg_AV_x4 = 0.0;
        double alg_AV_x4_i = 0.0;
        double alg_AV_x4_ss = 0.0;
        double alg_AV_xk1ss = 0.0;
        double alg_AV_xr = 0.0;
        double alg_AV_xrss = 0.0;
        double alg_AV_xs1ss = 0.0;
        double alg_AV_xs2ss = 0.0;
        double alg_AV_yss = 0.0;
        double alg_Iion_cm = 0.0;
        double alg_Istim = 0.0;

        double rate_CaMKt = 0.0;
        double rate_Cacsr = 0.0;
        double rate_Cai = 0.0;
        double rate_Cajsr = 0.0;
        double rate_Cansr = 0.0;
        double rate_Casl = 0.0;
        double rate_Cass = 0.0;
        double rate_ICaL_d = 0.0;
        double rate_ICaL_fcaf = 0.0;
        double rate_ICaL_fcafp = 0.0;
        double rate_ICaL_fcas = 0.0;
        double rate_ICaL_ff = 0.0;
        double rate_ICaL_ffp = 0.0;
        double rate_ICaL_fs = 0.0;
        double rate_ICaL_jca = 0.0;
        double rate_ICaL_nca = 0.0;
        double rate_ICaT_b = 0.0;
        double rate_ICaT_g = 0.0;
        double rate_IK1_xk1 = 0.0;
        double rate_IKr_xrf = 0.0;
        double rate_IKr_xrs = 0.0;
        double rate_IKs_xs1 = 0.0;
        double rate_IKs_xs2 = 0.0;
        double rate_INaL_hL = 0.0;
        double rate_INaL_hLp = 0.0;
        double rate_INaL_mL = 0.0;
        double rate_INa_hf = 0.0;
        double rate_INa_hs = 0.0;
        double rate_INa_hsp = 0.0;
        double rate_INa_j = 0.0;
        double rate_INa_jp = 0.0;
        double rate_INa_m = 0.0;
        double rate_IP3_u = 0.0;
        double rate_If_y = 0.0;
        double rate_Ito_a = 0.0;
        double rate_Ito_i1 = 0.0;
        double rate_Ito_i2 = 0.0;
        double rate_Ki = 0.0;
        double rate_Ksl = 0.0;
        double rate_Kss = 0.0;
        double rate_Nai = 0.0;
        double rate_Nasl = 0.0;
        double rate_Nass = 0.0;
        double rate_membrane_v = 0.0;
        double rate_ryr_Jrel1 = 0.0;
        double rate_ryr_Jrel2 = 0.0;


/* CaMK */
    alg_AV_CaMKb = CONSTANTS[AC_CaMKo] * (1.0 - state_CaMKt) / (1.0 + CONSTANTS[AC_KmCaM] / state_Cass);
    alg_AV_CaMKa = alg_AV_CaMKb + state_CaMKt;
    rate_CaMKt = CONSTANTS[AC_aCaMK] * alg_AV_CaMKb * (alg_AV_CaMKb + state_CaMKt) - CONSTANTS[AC_bCaMK] * state_CaMKt;

    /* IP3 */
    alg_AV_POip3 = CONSTANTS[AC_tauip3] * CONSTANTS[AC_IP3_IP3] * state_Cass * (1.0 - state_IP3_u) / ((1.0 + CONSTANTS[AC_IP3_IP3] * CONSTANTS[AC_k0] / CONSTANTS[AC_k0a]) * (1.0 + state_Cass * CONSTANTS[AC_k1] / CONSTANTS[AC_k1a]));
    rate_IP3_u = state_Cass * CONSTANTS[AC_k2] * (1.0 - state_IP3_u) - CONSTANTS[AC_k2a] * state_IP3_u;
    alg_AV_Jip3 = 10.92 * alg_AV_POip3 * (state_Cajsr - state_Cass);

    /* IpCa */
    alg_AV_IpCa_IpCa = CONSTANTS[AC_GpCa] * state_Casl / (CONSTANTS[AC_KmCap] + state_Casl);

    /* diff */
    alg_AV_Jdiff = (state_Cass - state_Casl) / CONSTANTS[AC_sstau];
    alg_AV_JdiffK = (state_Kss - state_Ksl) / CONSTANTS[AC_sstau];
    alg_AV_JdiffNa = (state_Nass - state_Nasl) / CONSTANTS[AC_sstau];
    alg_AV_Jgap = (state_Casl - state_Cai) / CONSTANTS[AC_gaptau];
    alg_AV_JgapK = (state_Ksl - state_Ki) / CONSTANTS[AC_gaptau];
    alg_AV_JgapNa = (state_Nasl - state_Nai) / CONSTANTS[AC_gaptau];



    /* trans_flux */
    alg_AV_Jtr1 = (state_Cansr - state_Cajsr) / 120.0;
    alg_AV_Jtr2 = (state_Cansr - state_Cacsr) / 120.0;

    /* INaCa_i */
    alg_AV_allo_i = 1.0 / (1.0 + pow(CONSTANTS[AC_KmCaAct] / state_Casl, 2.0));
    alg_AV_allo_ss = 1.0 / (1.0 + pow(CONSTANTS[AC_KmCaAct] / state_Cass, 2.0));
    alg_AV_h4_i = 1.0 + state_Nasl / CONSTANTS[AC_kna1] * (1.0 + state_Nasl / CONSTANTS[AC_kna2]);
    alg_AV_h4_ss = 1.0 + state_Nass / CONSTANTS[AC_kna1] * (1.0 + state_Nass / CONSTANTS[AC_kna2]);
    alg_AV_hca = exp(CONSTANTS[AC_qca] * state_membrane_v * CONSTANTS[AC_F] / (CONSTANTS[AC_R] * CONSTANTS[AC_T]));
    alg_AV_hna = exp(CONSTANTS[AC_qna] * state_membrane_v * CONSTANTS[AC_F] / (CONSTANTS[AC_R] * CONSTANTS[AC_T]));
    alg_AV_h1_i = 1.0 + state_Nasl / CONSTANTS[AC_kna3] * (1.0 + alg_AV_hna);
    alg_AV_h1_ss = 1.0 + state_Nass / CONSTANTS[AC_kna3] * (1.0 + alg_AV_hna);
    alg_AV_h5_i = state_Nasl * state_Nasl / (alg_AV_h4_i * CONSTANTS[AC_kna1] * CONSTANTS[AC_kna2]);
    alg_AV_h5_ss = state_Nass * state_Nass / (alg_AV_h4_ss * CONSTANTS[AC_kna1] * CONSTANTS[AC_kna2]);
    alg_AV_h6_i = 1.0 / alg_AV_h4_i;
    alg_AV_h6_ss = 1.0 / alg_AV_h4_ss;
    alg_AV_h7_i = 1.0 + CONSTANTS[AC_nao] / CONSTANTS[AC_kna3] * (1.0 + 1.0 / alg_AV_hna);
    alg_AV_h7_ss = 1.0 + CONSTANTS[AC_nao] / CONSTANTS[AC_kna3] * (1.0 + 1.0 / alg_AV_hna);
    alg_AV_h2_i = state_Nasl * alg_AV_hna / (CONSTANTS[AC_kna3] * alg_AV_h1_i);
    alg_AV_h2_ss = state_Nass * alg_AV_hna / (CONSTANTS[AC_kna3] * alg_AV_h1_ss);
    alg_AV_h3_i = 1.0 / alg_AV_h1_i;
    alg_AV_h3_ss = 1.0 / alg_AV_h1_ss;
    alg_AV_h8_i = CONSTANTS[AC_nao] / (CONSTANTS[AC_kna3] * alg_AV_hna * alg_AV_h7_i);
    alg_AV_h8_ss = CONSTANTS[AC_nao] / (CONSTANTS[AC_kna3] * alg_AV_hna * alg_AV_h7_ss);
    alg_AV_h9_i = 1.0 / alg_AV_h7_i;
    alg_AV_h9_ss = 1.0 / alg_AV_h7_ss;
    alg_AV_k6_i = alg_AV_h6_i * state_Casl * CONSTANTS[AC_kcaon];
    alg_AV_k6_ss = alg_AV_h6_ss * state_Cass * CONSTANTS[AC_kcaon];
    alg_AV_k3p_i = alg_AV_h9_i * CONSTANTS[AC_wca];
    alg_AV_k3p_ss = alg_AV_h9_ss * CONSTANTS[AC_wca];
    alg_AV_k3pp_i = alg_AV_h8_i * CONSTANTS[AC_wnaca];
    alg_AV_k3pp_ss = alg_AV_h8_ss * CONSTANTS[AC_wnaca];
    alg_AV_k4p_i = alg_AV_h3_i * CONSTANTS[AC_wca] / alg_AV_hca;
    alg_AV_k4p_ss = alg_AV_h3_ss * CONSTANTS[AC_wca] / alg_AV_hca;
    alg_AV_k4pp_i = alg_AV_h2_i * CONSTANTS[AC_wnaca];
    alg_AV_k4pp_ss = alg_AV_h2_ss * CONSTANTS[AC_wnaca];
    alg_AV_k7_i = alg_AV_h5_i * alg_AV_h2_i * CONSTANTS[AC_wna];
    alg_AV_k7_ss = alg_AV_h5_ss * alg_AV_h2_ss * CONSTANTS[AC_wna];
    alg_AV_k8_i = alg_AV_h8_i * CONSTANTS[AC_h11_i] * CONSTANTS[AC_wna];
    alg_AV_k8_ss = alg_AV_h8_ss * CONSTANTS[AC_h11_ss] * CONSTANTS[AC_wna];
    alg_AV_k3_i = alg_AV_k3p_i + alg_AV_k3pp_i;
    alg_AV_k3_ss = alg_AV_k3p_ss + alg_AV_k3pp_ss;
    alg_AV_k4_i = alg_AV_k4p_i + alg_AV_k4pp_i;
    alg_AV_k4_ss = alg_AV_k4p_ss + alg_AV_k4pp_ss;
    alg_AV_x1_i = CONSTANTS[AC_k2_i] * alg_AV_k4_i * (alg_AV_k7_i + alg_AV_k6_i) + CONSTANTS[AC_k5_i] * alg_AV_k7_i * (CONSTANTS[AC_k2_i] + alg_AV_k3_i);
    alg_AV_x1_ss = CONSTANTS[AC_k2_ss] * alg_AV_k4_ss * (alg_AV_k7_ss + alg_AV_k6_ss) + CONSTANTS[AC_k5_ss] * alg_AV_k7_ss * (CONSTANTS[AC_k2_ss] + alg_AV_k3_ss);
    alg_AV_x2_i = CONSTANTS[AC_k1_i] * alg_AV_k7_i * (alg_AV_k4_i + CONSTANTS[AC_k5_i]) + alg_AV_k4_i * alg_AV_k6_i * (CONSTANTS[AC_k1_i] + alg_AV_k8_i);
    alg_AV_x2_ss = CONSTANTS[AC_k1_ss] * alg_AV_k7_ss * (alg_AV_k4_ss + CONSTANTS[AC_k5_ss]) + alg_AV_k4_ss * alg_AV_k6_ss * (CONSTANTS[AC_k1_ss] + alg_AV_k8_ss);
    alg_AV_x3_i = CONSTANTS[AC_k1_i] * alg_AV_k3_i * (alg_AV_k7_i + alg_AV_k6_i) + alg_AV_k8_i * alg_AV_k6_i * (CONSTANTS[AC_k2_i] + alg_AV_k3_i);
    alg_AV_x3_ss = CONSTANTS[AC_k1_ss] * alg_AV_k3_ss * (alg_AV_k7_ss + alg_AV_k6_ss) + alg_AV_k8_ss * alg_AV_k6_ss * (CONSTANTS[AC_k2_ss] + alg_AV_k3_ss);
    alg_AV_x4_i = CONSTANTS[AC_k2_i] * alg_AV_k8_i * (alg_AV_k4_i + CONSTANTS[AC_k5_i]) + alg_AV_k3_i * CONSTANTS[AC_k5_i] * (CONSTANTS[AC_k1_i] + alg_AV_k8_i);
    alg_AV_x4_ss = CONSTANTS[AC_k2_ss] * alg_AV_k8_ss * (alg_AV_k4_ss + CONSTANTS[AC_k5_ss]) + alg_AV_k3_ss * CONSTANTS[AC_k5_ss] * (CONSTANTS[AC_k1_ss] + alg_AV_k8_ss);
    alg_AV_E1_i = alg_AV_x1_i / (alg_AV_x1_i + alg_AV_x2_i + alg_AV_x3_i + alg_AV_x4_i);
    alg_AV_E1_ss = alg_AV_x1_ss / (alg_AV_x1_ss + alg_AV_x2_ss + alg_AV_x3_ss + alg_AV_x4_ss);
    alg_AV_E2_i = alg_AV_x2_i / (alg_AV_x1_i + alg_AV_x2_i + alg_AV_x3_i + alg_AV_x4_i);
    alg_AV_E2_ss = alg_AV_x2_ss / (alg_AV_x1_ss + alg_AV_x2_ss + alg_AV_x3_ss + alg_AV_x4_ss);
    alg_AV_E3_i = alg_AV_x3_i / (alg_AV_x1_i + alg_AV_x2_i + alg_AV_x3_i + alg_AV_x4_i);
    alg_AV_E3_ss = alg_AV_x3_ss / (alg_AV_x1_ss + alg_AV_x2_ss + alg_AV_x3_ss + alg_AV_x4_ss);
    alg_AV_E4_i = alg_AV_x4_i / (alg_AV_x1_i + alg_AV_x2_i + alg_AV_x3_i + alg_AV_x4_i);
    alg_AV_E4_ss = alg_AV_x4_ss / (alg_AV_x1_ss + alg_AV_x2_ss + alg_AV_x3_ss + alg_AV_x4_ss);
    alg_AV_JncxCa_i = alg_AV_E2_i * CONSTANTS[AC_k2_i] - alg_AV_E1_i * CONSTANTS[AC_k1_i];
    alg_AV_JncxCa_ss = alg_AV_E2_ss * CONSTANTS[AC_k2_ss] - alg_AV_E1_ss * CONSTANTS[AC_k1_ss];
    alg_AV_JncxNa_i = 3.0 * (alg_AV_E4_i * alg_AV_k7_i - alg_AV_E1_i * alg_AV_k8_i) + alg_AV_E3_i * alg_AV_k4pp_i - alg_AV_E2_i * alg_AV_k3pp_i;
    alg_AV_JncxNa_ss = 3.0 * (alg_AV_E4_ss * alg_AV_k7_ss - alg_AV_E1_ss * alg_AV_k8_ss) + alg_AV_E3_ss * alg_AV_k4pp_ss - alg_AV_E2_ss * alg_AV_k3pp_ss;
    alg_AV_INaCa_i_INaCa_i = 0.8 * CONSTANTS[AC_Gncx] * alg_AV_allo_i * (CONSTANTS[AC_zna] * alg_AV_JncxNa_i + CONSTANTS[AC_zca] * alg_AV_JncxCa_i);
    alg_AV_INaCa_ss = 0.2 * CONSTANTS[AC_Gncx] * alg_AV_allo_ss * (CONSTANTS[AC_zna] * alg_AV_JncxNa_ss + CONSTANTS[AC_zca] * alg_AV_JncxCa_ss);

    /* INaK */
    alg_AV_Knai = CONSTANTS[AC_Knai0] * exp(CONSTANTS[AC_delta] * state_membrane_v * CONSTANTS[AC_F] / (3.0 * CONSTANTS[AC_R] * CONSTANTS[AC_T]));
    alg_AV_Knao = CONSTANTS[AC_Knao0] * exp((1.0 - CONSTANTS[AC_delta]) * state_membrane_v * CONSTANTS[AC_F] / (3.0 * CONSTANTS[AC_R] * CONSTANTS[AC_T]));
    alg_AV_P = CONSTANTS[AC_eP] / (1.0 + CONSTANTS[AC_H] / CONSTANTS[AC_Khp] + state_Nasl / CONSTANTS[AC_Knap] + state_Ksl / CONSTANTS[AC_Kxkur]);
    alg_AV_a1 = CONSTANTS[AC_k1p] * pow(state_Nasl / alg_AV_Knai,
                                               3.0) / (pow(1.0 + state_Nasl / alg_AV_Knai, 3.0) + pow(1.0 + state_Ksl / CONSTANTS[AC_Kki], 2.0) - 1.0);
    alg_AV_a3 = CONSTANTS[AC_k3p] * pow(CONSTANTS[AC_ko] / CONSTANTS[AC_Kko], 2.0) / (pow(1.0 + CONSTANTS[AC_nao] / alg_AV_Knao, 3.0) + pow(1.0 + CONSTANTS[AC_ko] / CONSTANTS[AC_Kko], 2.0) - 1.0);
    alg_AV_b2 = CONSTANTS[AC_k2m] * pow(CONSTANTS[AC_nao] / alg_AV_Knao,
                                               3.0) / (pow(1.0 + CONSTANTS[AC_nao] / alg_AV_Knao, 3.0) + pow(1.0 + CONSTANTS[AC_ko] / CONSTANTS[AC_Kko], 2.0) - 1.0);
    alg_AV_b3 = CONSTANTS[AC_k3m] * alg_AV_P * CONSTANTS[AC_H] / (1.0 + CONSTANTS[AC_MgATP] / CONSTANTS[AC_Kmgatp]);
    alg_AV_b4 = CONSTANTS[AC_k4m] * pow(state_Ksl / CONSTANTS[AC_Kki], 2.0) / (pow(1.0 + state_Nasl / alg_AV_Knai, 3.0) + pow(1.0 + state_Ksl / CONSTANTS[AC_Kki], 2.0) - 1.0);
    alg_AV_x1 = CONSTANTS[AC_a4] * alg_AV_a1 * CONSTANTS[AC_a2] + alg_AV_b2 * alg_AV_b4 * alg_AV_b3 + CONSTANTS[AC_a2] * alg_AV_b4 * alg_AV_b3 + alg_AV_b3 * alg_AV_a1 * CONSTANTS[AC_a2];
    alg_AV_x2 = alg_AV_b2 * CONSTANTS[AC_b1] * alg_AV_b4 + alg_AV_a1 * CONSTANTS[AC_a2] * alg_AV_a3 + alg_AV_a3 * CONSTANTS[AC_b1] * alg_AV_b4 + CONSTANTS[AC_a2] * alg_AV_a3 * alg_AV_b4;
    alg_AV_x3 = CONSTANTS[AC_a2] * alg_AV_a3 * CONSTANTS[AC_a4] + alg_AV_b3 * alg_AV_b2 * CONSTANTS[AC_b1] + alg_AV_b2 * CONSTANTS[AC_b1] * CONSTANTS[AC_a4] + alg_AV_a3 * CONSTANTS[AC_a4] * CONSTANTS[AC_b1];
    alg_AV_x4 = alg_AV_b4 * alg_AV_b3 * alg_AV_b2 + alg_AV_a3 * CONSTANTS[AC_a4] * alg_AV_a1 + alg_AV_b2 * CONSTANTS[AC_a4] * alg_AV_a1 + alg_AV_b3 * alg_AV_b2 * alg_AV_a1;
    alg_AV_E1 = alg_AV_x1 / (alg_AV_x1 + alg_AV_x2 + alg_AV_x3 + alg_AV_x4);
    alg_AV_E2 = alg_AV_x2 / (alg_AV_x1 + alg_AV_x2 + alg_AV_x3 + alg_AV_x4);
    alg_AV_E3 = alg_AV_x3 / (alg_AV_x1 + alg_AV_x2 + alg_AV_x3 + alg_AV_x4);
    alg_AV_E4 = alg_AV_x4 / (alg_AV_x1 + alg_AV_x2 + alg_AV_x3 + alg_AV_x4);
    alg_AV_JnakK = 2.0 * (alg_AV_E4 * CONSTANTS[AC_b1] - alg_AV_E3 * alg_AV_a1);
    alg_AV_JnakNa = 3.0 * (alg_AV_E1 * alg_AV_a3 - alg_AV_E2 * alg_AV_b3);
    alg_AV_INaK_INaK = CONSTANTS[AC_Pnak] * (CONSTANTS[AC_zna] * alg_AV_JnakNa + CONSTANTS[AC_zk] * alg_AV_JnakK);

    /* SERCA */
    alg_AV_dkmplb = CONSTANTS[AC_dkmplbbar] * alg_AV_CaMKa / (alg_AV_CaMKa + CONSTANTS[AC_KmCaMK]);
    alg_AV_dqupcamk = CONSTANTS[AC_dqupcamkbar] * alg_AV_CaMKa / (alg_AV_CaMKa + CONSTANTS[AC_KmCaMK]);
    alg_AV_Jup1 = 0.0002 * (alg_AV_dqupcamk + 1.0) / (1.0 + (CONSTANTS[AC_kmup] - alg_AV_dkmplb) / state_Casl) - 0.00105 * state_Cansr / CONSTANTS[AC_nsrbar];
    alg_AV_Jup2 = 0.0026 * (alg_AV_dqupcamk + 1.0) / (1.0 + (CONSTANTS[AC_kmup] - alg_AV_dkmplb) / state_Cai) - 0.0042 * state_Cansr / CONSTANTS[AC_nsrbar];

    /* reversal_potentials */
    alg_AV_ECa = CONSTANTS[AC_R] * CONSTANTS[AC_T] / (2.0 * CONSTANTS[AC_F]) * log(CONSTANTS[AC_cao] / state_Casl);
    alg_AV_EK = CONSTANTS[AC_R] * CONSTANTS[AC_T] / CONSTANTS[AC_F] * log(CONSTANTS[AC_ko] / state_Ksl);
    alg_AV_ENa = CONSTANTS[AC_R] * CONSTANTS[AC_T] / CONSTANTS[AC_F] * log(CONSTANTS[AC_nao] / state_Nasl);
    alg_AV_EKs = CONSTANTS[AC_R] * CONSTANTS[AC_T] / CONSTANTS[AC_F] * log((CONSTANTS[AC_ko] + CONSTANTS[AC_PKNa] * CONSTANTS[AC_nao]) / (state_Ksl + CONSTANTS[AC_PKNa] * state_Nasl));

    /* ICaT */
    alg_AV_bss = 1.0 / (1.0 + exp(-(state_membrane_v + 30.0) / 7.0));
    alg_AV_gss = 1.0 / (1.0 + exp((state_membrane_v + 61.0) / 5.0));
    alg_AV_taub = 1.0 / (1.068 * exp((state_membrane_v + 16.3) / 30.0) + 1.068 * exp(-(state_membrane_v + 16.3) / 30.0));
    alg_AV_taug = 1.0 / (0.015 * exp((state_membrane_v + 71.7) / 15.4) + 0.015 * exp(-(state_membrane_v + 71.7) / 83.3));
    alg_AV_ICaT_ICaT = CONSTANTS[AC_GCaT] * state_ICaT_b * state_ICaT_g * (state_membrane_v - alg_AV_ECa);
    rate_ICaT_b = (alg_AV_bss - state_ICaT_b) / alg_AV_taub;
    rate_ICaT_g = (alg_AV_gss - state_ICaT_g) / alg_AV_taug;

    /* IK1 */
    alg_AV_rk1 = 1.0 / (1.0 + exp((state_membrane_v + 116.0 - 5.5 * CONSTANTS[AC_ko]) / 11.0));
    alg_AV_txk1 = 122.2 / (exp(-(state_membrane_v + 127.2) / 20.36) + exp((state_membrane_v + 236.8) / 69.33));
    alg_AV_xk1ss = 1.0 / (1.0 + exp(-(state_membrane_v + 2.5538 * CONSTANTS[AC_ko] + 144.59) / (1.5692 * CONSTANTS[AC_ko] + 3.8115)));
    alg_AV_IK1_IK1 = CONSTANTS[AC_GK1] * 2.3238 * sqrt(CONSTANTS[AC_ko] / 5.4) * alg_AV_rk1 * state_IK1_xk1 * (state_membrane_v - alg_AV_EK);
    rate_IK1_xk1 = (alg_AV_xk1ss - state_IK1_xk1) / alg_AV_txk1;

    /* IKr */
    alg_AV_Axrf = 1.0 / (1.0 + exp((state_membrane_v + 54.81) / 38.21));
    alg_AV_rkr = 1.0 / (1.0 + exp((state_membrane_v + 55.0) / (0.32 * 75.0))) * 1.0 / (1.0 + exp((state_membrane_v - 10.0) / (0.32 * 30.0)));
    alg_AV_txrf = 12.98 + 1.0 / (0.3652 * exp((state_membrane_v + 17.6 - 31.66) / 3.869) + 4.123e-05 * exp(-(state_membrane_v + 17.6 - 47.78) / 20.38));
    alg_AV_txrs = 1.865 + 1.0 / (0.06629 * exp((state_membrane_v + 17.2 - 34.7) / 7.355) + 1.128e-05 * exp(-(state_membrane_v + 17.2 - 29.74) / 25.94));
    alg_AV_xrss = 1.0 / (1.0 + exp(-(state_membrane_v + 8.337) / 6.789));
    alg_AV_Axrs = 1.0 - alg_AV_Axrf;
    rate_IKr_xrf = (alg_AV_xrss - state_IKr_xrf) / alg_AV_txrf;
    rate_IKr_xrs = (alg_AV_xrss - state_IKr_xrs) / alg_AV_txrs;
    alg_AV_xr = alg_AV_Axrf * state_IKr_xrf + alg_AV_Axrs * state_IKr_xrs;
    alg_AV_IKr_IKr = CONSTANTS[AC_GKr] * sqrt(CONSTANTS[AC_ko] / 5.4) * alg_AV_xr * alg_AV_rkr * (state_membrane_v - alg_AV_EK);

    /* IKs */
    alg_AV_KsCa = 1.0 + 0.6 / (1.0 + pow(3.8e-05 / state_Casl, 1.4));
    alg_AV_txs1 = 817.3 + 1.0 / (0.0002326 * exp((state_membrane_v + 48.28) / 17.8) + 0.001292 * exp(-(state_membrane_v + 210.0) / 230.0));
    alg_AV_txs2 = 1.0 / (0.01 * exp((state_membrane_v - 50.0) / 20.0) + 0.0193 * exp(-(state_membrane_v + 66.54) / 31.0));
    alg_AV_xs1ss = 1.0 / (1.0 + exp(-(state_membrane_v + 11.6) / 8.932));
    alg_AV_IKs_IKs = CONSTANTS[AC_GKs] * alg_AV_KsCa * state_IKs_xs1 * state_IKs_xs2 * (state_membrane_v - alg_AV_EKs);
    alg_AV_xs2ss = alg_AV_xs1ss;
    rate_IKs_xs1 = (alg_AV_xs1ss - state_IKs_xs1) / alg_AV_txs1;
    rate_IKs_xs2 = (alg_AV_xs2ss - state_IKs_xs2) / alg_AV_txs2;

    /* INa */
    alg_AV_fINap = 1.0 / (1.0 + CONSTANTS[AC_KmCaMK] / alg_AV_CaMKa);
    alg_AV_hssp = 1.0 / (1.0 + exp((state_membrane_v + 84.7) / 6.22));
    alg_AV_thf = 1.0 / (3.686e-06 * exp(-(state_membrane_v + 3.8875) / 7.8579) + 16.0 * exp((state_membrane_v - 0.4963) / 9.1843));
    alg_AV_ths = 1.0 / (0.009794 * exp(-(state_membrane_v + 17.95) / 28.05) + 0.3343 * exp((state_membrane_v + 5.73) / 56.66));
    alg_AV_tj = 4.859 + 1.0 / (0.8628 * exp(-(state_membrane_v + 116.7258) / 7.6005) + 1.1096 * exp((state_membrane_v + 6.2719) / 9.0358));
    alg_AV_hss = 1.0 / (1.0 + exp((state_membrane_v + CONSTANTS[AC_hssV1]) / CONSTANTS[AC_hssV2]));
    alg_AV_mss = 1.0 / (1.0 + exp(-(state_membrane_v + CONSTANTS[AC_mssV1]) / CONSTANTS[AC_mssV2]));
    alg_AV_thsp = 3.0 * alg_AV_ths;
    alg_AV_tjp = 1.46 * alg_AV_tj;
    alg_AV_tm = 1.0 / (CONSTANTS[AC_mtD1] * exp((state_membrane_v + CONSTANTS[AC_mtV1]) / CONSTANTS[AC_mtV2]) + CONSTANTS[AC_mtD2] * exp(-(state_membrane_v + CONSTANTS[AC_mtV3]) / CONSTANTS[AC_mtV4]));
    alg_AV_h = CONSTANTS[AC_Ahf] * state_INa_hf + CONSTANTS[AC_Ahs] * state_INa_hs;
    alg_AV_hp = CONSTANTS[AC_Ahf] * state_INa_hf + CONSTANTS[AC_Ahs] * state_INa_hsp;
    alg_AV_jss = alg_AV_hss;
    rate_INa_hf = (alg_AV_hss - state_INa_hf) / alg_AV_thf;
    rate_INa_hs = (alg_AV_hss - state_INa_hs) / alg_AV_ths;
    rate_INa_hsp = (alg_AV_hssp - state_INa_hsp) / alg_AV_thsp;
    rate_INa_m = (alg_AV_mss - state_INa_m) / alg_AV_tm;
    alg_AV_INa_INa = CONSTANTS[AC_GNa] * (state_membrane_v - alg_AV_ENa) * pow(state_INa_m, 3.0) * ((1.0 - alg_AV_fINap) * alg_AV_h * state_INa_j + alg_AV_fINap * alg_AV_hp * state_INa_jp);
    rate_INa_j = (alg_AV_jss - state_INa_j) / alg_AV_tj;
    rate_INa_jp = (alg_AV_jss - state_INa_jp) / alg_AV_tjp;

    /* If */
    alg_AV_tauy = 2000.0 / (exp((state_membrane_v + 57.0) / 60.0) + exp(-(state_membrane_v + 132.0) / 10.0));
    alg_AV_yss = 1.0 / (1.0 + exp((state_membrane_v + 87.0) / 9.5));
    alg_AV_IfK = CONSTANTS[AC_GfK] * state_If_y * state_If_y * (state_membrane_v - alg_AV_EK);
    alg_AV_IfNa = CONSTANTS[AC_GfNa] * state_If_y * state_If_y * (state_membrane_v - alg_AV_ENa);
    rate_If_y = (alg_AV_yss - state_If_y) / alg_AV_tauy;
    alg_AV_If_If = alg_AV_IfNa + alg_AV_IfK;

    /* Isus */
    alg_AV_asus = 1.0 / (1.0 + exp(-(state_membrane_v - 12.0) / 16.0));
    alg_AV_Isus_Isus = CONSTANTS[AC_Gsus] * alg_AV_asus * (state_membrane_v - alg_AV_EK);

    /* Ito */
    alg_AV_ass = 1.0 / (1.0 + exp((20.0 - state_membrane_v) / 13.0));
    alg_AV_iss = 1.0 / (1.0 + exp((27.0 + state_membrane_v) / 13.0));
    alg_AV_taua = 1.0515 / (1.0 / (1.2089 * (1.0 + exp(-(state_membrane_v - 18.4099) / 29.3814))) + 3.5 / (1.0 + exp((state_membrane_v + 100.0) / 29.3814)));
    alg_AV_tauif = 6.162 + 1.0 / (0.3933 * exp(-(state_membrane_v + 100.0) / 100.0) + 0.08004 * exp((state_membrane_v - 8.0) / 8.59));
    alg_AV_tauis = 43.0 + 1.0 / (0.001416 * exp(-(state_membrane_v + 96.52) / 59.05) + 1.78e-08 * exp((state_membrane_v + 114.1) / 8.079));
    alg_AV_Ito_Ito = CONSTANTS[AC_Gto] * state_Ito_a * state_Ito_i1 * state_Ito_i2 * (state_membrane_v - alg_AV_EK);
    rate_Ito_a = (alg_AV_ass - state_Ito_a) / alg_AV_taua;
    rate_Ito_i1 = (alg_AV_iss - state_Ito_i1) / alg_AV_tauis;
    rate_Ito_i2 = (alg_AV_iss - state_Ito_i2) / alg_AV_tauif;

    /* INaL */
    alg_AV_fINaLp = 1.0 / (1.0 + CONSTANTS[AC_KmCaMK] / alg_AV_CaMKa);
    alg_AV_hLss = 1.0 / (1.0 + exp((state_membrane_v + 87.61) / 7.488));
    alg_AV_hLssp = 1.0 / (1.0 + exp((state_membrane_v + 93.81) / 7.488));
    alg_AV_mLss = 1.0 / (1.0 + exp(-(state_membrane_v + 42.85) / 5.264));
    alg_AV_tmL = alg_AV_tm;
    alg_AV_INaL_INaL = CONSTANTS[AC_GNaL] * (state_membrane_v - alg_AV_ENa) * state_INaL_mL * ((1.0 - alg_AV_fINaLp) * state_INaL_hL + alg_AV_fINaLp * state_INaL_hLp);
    rate_INaL_hL = (alg_AV_hLss - state_INaL_hL) / CONSTANTS[AC_thL];
    rate_INaL_mL = (alg_AV_mLss - state_INaL_mL) / alg_AV_tmL;
    rate_INaL_hLp = (alg_AV_hLssp - state_INaL_hLp) / CONSTANTS[AC_thLp];

    /* ICaL */
    alg_AV_Afcaf = 0.3 + 0.6 / (1.0 + exp((state_membrane_v - 10.0) / 10.0));
    alg_AV_dss = 1.0 / (1.0 + exp(-(state_membrane_v + 3.94 + 3.3) / 4.23));
    alg_AV_fICaLp = 1.0 / (1.0 + CONSTANTS[AC_KmCaMK] / alg_AV_CaMKa);
    alg_AV_fss = 1.0 / (1.0 + exp((state_membrane_v + 19.58 + 3.3) / 3.696));
    alg_AV_km2n = state_ICaL_jca * 1.0;
    alg_AV_td = 0.6 + 1.0 / (exp(-0.05 * (state_membrane_v + 6.0)) + exp(0.09 * (state_membrane_v + 14.0)));
    alg_AV_tfcaf = 0.72 * (7.0 + 1.0 / (0.04 * exp(-(state_membrane_v + 15.19 - 4.0) / 7.0) + 0.04 * exp((state_membrane_v + 15.19 - 4.0) / 7.0)));
    alg_AV_tfcas = 0.49 * (100.0 + 1.0 / (0.00012 * exp(-(state_membrane_v + 15.19) / 3.0) + 0.00012 * exp((state_membrane_v + 15.19) / 7.0)));
    alg_AV_tff = 7.0 + 1.0 / (0.0045 * exp(-(state_membrane_v + 20.0 + 15.19) / 10.0) + 0.0045 * exp((state_membrane_v + 20.0 + 15.19) / 10.0));
    alg_AV_tfs = 1000.0 + 1.0 / (3.5e-05 * exp(-(state_membrane_v + 5.0 + 15.19) / 4.0) + 3.5e-05 * exp((state_membrane_v + 5.0 + 15.19) / 6.0));
    alg_AV_Afcas = 1.0 - alg_AV_Afcaf;
    alg_AV_anca = 1.0 / (CONSTANTS[AC_k2n] / alg_AV_km2n + pow(1.0 + CONSTANTS[AC_Kmn] / state_Cass, 4.0));
    alg_AV_fcass = alg_AV_fss;
    alg_AV_tfcafp = 2.5 * alg_AV_tfcaf;
    alg_AV_tffp = 2.5 * alg_AV_tff;
    rate_ICaL_d = (alg_AV_dss - state_ICaL_d) / alg_AV_td;
    rate_ICaL_ff = (alg_AV_fss - state_ICaL_ff) / alg_AV_tff;
    rate_ICaL_fs = (alg_AV_fss - state_ICaL_fs) / alg_AV_tfs;
    alg_AV_f = CONSTANTS[AC_Aff] * state_ICaL_ff + CONSTANTS[AC_Afs] * state_ICaL_fs;
    alg_AV_fca = alg_AV_Afcaf * state_ICaL_fcaf + alg_AV_Afcas * state_ICaL_fcas;
    alg_AV_fcap = alg_AV_Afcaf * state_ICaL_fcafp + alg_AV_Afcas * state_ICaL_fcas;
    alg_AV_fp = CONSTANTS[AC_Aff] * state_ICaL_ffp + CONSTANTS[AC_Afs] * state_ICaL_fs;
    rate_ICaL_fcaf = (alg_AV_fcass - state_ICaL_fcaf) / alg_AV_tfcaf;
    rate_ICaL_fcafp = (alg_AV_fcass - state_ICaL_fcafp) / alg_AV_tfcafp;
    rate_ICaL_fcas = (alg_AV_fcass - state_ICaL_fcas) / alg_AV_tfcas;
    rate_ICaL_ffp = (alg_AV_fss - state_ICaL_ffp) / alg_AV_tffp;
    rate_ICaL_jca = (alg_AV_fcass - state_ICaL_jca) / CONSTANTS[AC_tjca];
    rate_ICaL_nca = alg_AV_anca * CONSTANTS[AC_k2n] - state_ICaL_nca * alg_AV_km2n;

    /* intracellular_ions */
    rate_Cansr = alg_AV_Jup1 + alg_AV_Jup2 - (alg_AV_Jtr1 * CONSTANTS[AC_vjsr] / CONSTANTS[AC_vnsr] + alg_AV_Jtr2 * CONSTANTS[AC_vcsr] / CONSTANTS[AC_vnsr]);
    rate_Ki = alg_AV_JgapK * CONSTANTS[AC_vsl] / CONSTANTS[AC_vmyo];
    rate_Nai = alg_AV_JgapNa * CONSTANTS[AC_vsl] / CONSTANTS[AC_vmyo];
    alg_AV_Bcacsr = 1.0 / (1.0 + CONSTANTS[AC_csqnmax] * CONSTANTS[AC_kmcsqn] / pow(CONSTANTS[AC_kmcsqn] + state_Cacsr, 2.0));
    alg_AV_Bcai = 1.0 / (1.0 + CONSTANTS[AC_cmdnmax] * CONSTANTS[AC_kmcmdn] / pow(CONSTANTS[AC_kmcmdn] + state_Cai, 2.0) + CONSTANTS[AC_trpnmax] * CONSTANTS[AC_kmtrpn] / pow(CONSTANTS[AC_kmtrpn] + state_Cai, 2.0));
    alg_AV_Bcajsr = 1.0 / (1.0 + CONSTANTS[AC_csqnmaxsl] * CONSTANTS[AC_kmcsqn] / pow(CONSTANTS[AC_kmcsqn] + state_Cajsr, 2.0));
    alg_AV_Bcasl = 1.0 / (1.0 + CONSTANTS[AC_cmdnmaxsl] * CONSTANTS[AC_kmcmdn] / pow(CONSTANTS[AC_kmcmdn] + state_Casl, 2.0) + CONSTANTS[AC_trpnmaxsl] * CONSTANTS[AC_kmtrpn] / pow(CONSTANTS[AC_kmtrpn] + state_Casl, 2.0));
    alg_AV_Bcass = 1.0 / (1.0 + CONSTANTS[AC_BSRmax] * CONSTANTS[AC_KmBSR] / pow(CONSTANTS[AC_KmBSR] + state_Cass, 2.0) + CONSTANTS[AC_BSLmax] * CONSTANTS[AC_KmBSL] / pow(CONSTANTS[AC_KmBSL] + state_Cass, 2.0));
    rate_Cacsr = alg_AV_Bcacsr * (alg_AV_Jtr2 - state_ryr_Jrel2);
    rate_Cai = alg_AV_Bcai * (alg_AV_Jgap * CONSTANTS[AC_vsl] / CONSTANTS[AC_vmyo] + state_ryr_Jrel2 * CONSTANTS[AC_vcsr] / CONSTANTS[AC_vmyo] - alg_AV_Jup2 * CONSTANTS[AC_vnsr] / CONSTANTS[AC_vmyo]);
    rate_Cajsr = alg_AV_Bcajsr * (alg_AV_Jtr1 - (state_ryr_Jrel1 + alg_AV_Jip3));

    /* membrane */
    alg_AV_vffrt = state_membrane_v * CONSTANTS[AC_F] * CONSTANTS[AC_F] / (CONSTANTS[AC_R] * CONSTANTS[AC_T]);
    alg_AV_vfrt = state_membrane_v * CONSTANTS[AC_F] / (CONSTANTS[AC_R] * CONSTANTS[AC_T]);


    /* ryr */
    alg_AV_REL2 = alg_AV_Jgap * CONSTANTS[AC_vsl] / CONSTANTS[AC_vmyo] + -alg_AV_Jup2 * CONSTANTS[AC_vnsr] / CONSTANTS[AC_vmyo] + state_ryr_Jrel2 * CONSTANTS[AC_vcsr] / CONSTANTS[AC_vmyo];
    alg_AV_ireltau = 2.0 * (1.0 + 1.0 * 1.0 / (1.0 + pow(0.28 / alg_AV_CaMKa, 8.0))) / (1.0 + 0.0123 / state_Cajsr);
    alg_AV_ireltau2 = 6.0 * (1.0 + 1.0 * 1.0 / (1.0 + pow(0.28 / alg_AV_CaMKa, 8.0))) / (1.0 + 0.0123 / state_Cacsr);
    alg_AV_irelss2 = ((alg_AV_REL2 > 0.0) ? 91.0 * (1.0 + 1.0 * 1.0 / (1.0 + pow(0.28 / alg_AV_CaMKa, 8.0))) * alg_AV_REL2 / (1.0 + pow(1.0 / state_Cacsr, 8.0)) : 0.0);
    rate_ryr_Jrel2 = (alg_AV_irelss2 - state_ryr_Jrel2) / alg_AV_ireltau2;

    /* *remaining* */
    alg_AV_PhiCaK = 1.0 * alg_AV_vffrt * (0.75 * state_Kss * exp(1.0 * alg_AV_vfrt) - 0.75 * CONSTANTS[AC_ko]) / (exp(1.0 * alg_AV_vfrt) - 1.0);
    alg_AV_PhiCaL = 4.0 * alg_AV_vffrt * (state_Cass * exp(2.0 * alg_AV_vfrt) - 0.341 * CONSTANTS[AC_cao]) / (exp(2.0 * alg_AV_vfrt) - 1.0);
    alg_AV_PhiCaNa = 1.0 * alg_AV_vffrt * (0.75 * state_Nass * exp(1.0 * alg_AV_vfrt) - 0.75 * CONSTANTS[AC_nao]) / (exp(1.0 * alg_AV_vfrt) - 1.0);
    alg_AV_ICab_ICab = CONSTANTS[AC_PCab] * 4.0 * alg_AV_vffrt * (state_Casl * exp(2.0 * alg_AV_vfrt) - 0.341 * CONSTANTS[AC_cao]) / (exp(2.0 * alg_AV_vfrt) - 1.0);
    alg_AV_INab_INab = CONSTANTS[AC_PNab] * alg_AV_vffrt * (state_Nasl * exp(alg_AV_vfrt) - CONSTANTS[AC_nao]) / (exp(alg_AV_vfrt) - 1.0);
    rate_Ksl = -(alg_AV_Ito_Ito + alg_AV_Isus_Isus + alg_AV_IKr_IKr + alg_AV_IKs_IKs + alg_AV_IfK + alg_AV_IK1_IK1 + alg_Istim - 2.0 * alg_AV_INaK_INaK) * CONSTANTS[AC_cm] * CONSTANTS[AC_Acap] / (CONSTANTS[AC_F] * CONSTANTS[AC_vsl]) + alg_AV_JdiffK * CONSTANTS[AC_vss] / CONSTANTS[AC_vsl] + -alg_AV_JgapK;
    alg_AV_ICaK = (1.0 - alg_AV_fICaLp) * CONSTANTS[AC_PCaK] * alg_AV_PhiCaK * state_ICaL_d * (alg_AV_f * (1.0 - state_ICaL_nca) + state_ICaL_jca * alg_AV_fca * state_ICaL_nca) + alg_AV_fICaLp * CONSTANTS[AC_PCaKp] * alg_AV_PhiCaK * state_ICaL_d * (alg_AV_fp * (1.0 - state_ICaL_nca) + state_ICaL_jca * alg_AV_fcap * state_ICaL_nca);
    alg_AV_ICaL_ICaL = (1.0 - alg_AV_fICaLp) * CONSTANTS[AC_PCa] * alg_AV_PhiCaL * state_ICaL_d * (alg_AV_f * (1.0 - state_ICaL_nca) + state_ICaL_jca * alg_AV_fca * state_ICaL_nca) + alg_AV_fICaLp * CONSTANTS[AC_PCap] * alg_AV_PhiCaL * state_ICaL_d * (alg_AV_fp * (1.0 - state_ICaL_nca) + state_ICaL_jca * alg_AV_fcap * state_ICaL_nca);
    alg_AV_ICaNa = (1.0 - alg_AV_fICaLp) * CONSTANTS[AC_PCaNa] * alg_AV_PhiCaNa * state_ICaL_d * (alg_AV_f * (1.0 - state_ICaL_nca) + state_ICaL_jca * alg_AV_fca * state_ICaL_nca) + alg_AV_fICaLp * CONSTANTS[AC_PCaNap] * alg_AV_PhiCaNa * state_ICaL_d * (alg_AV_fp * (1.0 - state_ICaL_nca) + state_ICaL_jca * alg_AV_fcap * state_ICaL_nca);
    rate_Casl = alg_AV_Bcasl * (-(alg_AV_IpCa_IpCa + alg_AV_ICab_ICab + alg_AV_ICaT_ICaT - 2.0 * alg_AV_INaCa_i_INaCa_i) * CONSTANTS[AC_cm] * CONSTANTS[AC_Acap] / (2.0 * CONSTANTS[AC_F] * CONSTANTS[AC_vsl]) - alg_AV_Jup1 * CONSTANTS[AC_vnsr] / CONSTANTS[AC_vsl] + -alg_AV_Jgap + alg_AV_Jdiff * CONSTANTS[AC_vss] / CONSTANTS[AC_vsl]);
    rate_Nasl = -(alg_AV_INa_INa + alg_AV_INaL_INaL + 3.0 * alg_AV_INaCa_i_INaCa_i + 3.0 * alg_AV_INaK_INaK + alg_AV_IfNa + alg_AV_INab_INab) * CONSTANTS[AC_Acap] * CONSTANTS[AC_cm] / (CONSTANTS[AC_F] * CONSTANTS[AC_vsl]) + alg_AV_JdiffNa * CONSTANTS[AC_vss] / CONSTANTS[AC_vsl] + -alg_AV_JgapNa;
    rate_Cass = alg_AV_Bcass * (-(alg_AV_ICaL_ICaL - 2.0 * alg_AV_INaCa_ss) * CONSTANTS[AC_cm] * CONSTANTS[AC_Acap] / (2.0 * CONSTANTS[AC_F] * CONSTANTS[AC_vss]) + (state_ryr_Jrel1 + alg_AV_Jip3) * CONSTANTS[AC_vjsr] / CONSTANTS[AC_vss] - alg_AV_Jdiff);
    rate_Kss = -alg_AV_ICaK * CONSTANTS[AC_cm] * CONSTANTS[AC_Acap] / (CONSTANTS[AC_F] * CONSTANTS[AC_vss]) - alg_AV_JdiffK;
    rate_Nass = -(alg_AV_ICaNa + 3.0 * alg_AV_INaCa_ss) * CONSTANTS[AC_cm] * CONSTANTS[AC_Acap] / (CONSTANTS[AC_F] * CONSTANTS[AC_vss]) - alg_AV_JdiffNa;
    alg_AV_REL = -(alg_AV_ICaL_ICaL * CONSTANTS[AC_Acap] / (CONSTANTS[AC_vss] * CONSTANTS[AC_F] * 2.0) + -(state_ryr_Jrel1 + alg_AV_Jip3) * CONSTANTS[AC_vjsr] / CONSTANTS[AC_vss] + alg_AV_Jdiff);
    alg_AV_irelss = ((alg_AV_REL > 0.0) ? 15.0 * (1.0 + 1.0 * 1.0 / (1.0 + pow(0.28 / alg_AV_CaMKa, 8.0))) * alg_AV_REL / (1.0 + pow(1.0 / state_Cajsr, 8.0)) : 0.0);
    rate_ryr_Jrel1 = (alg_AV_irelss - state_ryr_Jrel1) / alg_AV_ireltau;

    alg_Iion_cm = alg_AV_INa_INa + alg_AV_INaL_INaL + alg_AV_Ito_Ito + alg_AV_Isus_Isus + alg_AV_ICaL_ICaL + alg_AV_ICaT_ICaT + alg_AV_ICaNa + alg_AV_ICaK + alg_AV_IKr_IKr + alg_AV_IKs_IKs + alg_AV_If_If + alg_AV_IK1_IK1 + alg_AV_INaCa_i_INaCa_i + alg_AV_INaCa_ss + alg_AV_INaK_INaK + alg_AV_INab_INab + alg_AV_IpCa_IpCa + alg_AV_ICab_ICab;





    alg_Istim = computeStimulusPOD(VOI, stimulus);

    if (solveVm)
    {
        rate_membrane_v = -alg_Iion_cm - alg_Istim;
    }

        RATES_SoA[CaMKt * N + cellI] = rate_CaMKt;
        RATES_SoA[Cacsr * N + cellI] = rate_Cacsr;
        RATES_SoA[Cai * N + cellI] = rate_Cai;
        RATES_SoA[Cajsr * N + cellI] = rate_Cajsr;
        RATES_SoA[Cansr * N + cellI] = rate_Cansr;
        RATES_SoA[Casl * N + cellI] = rate_Casl;
        RATES_SoA[Cass * N + cellI] = rate_Cass;
        RATES_SoA[ICaL_d * N + cellI] = rate_ICaL_d;
        RATES_SoA[ICaL_fcaf * N + cellI] = rate_ICaL_fcaf;
        RATES_SoA[ICaL_fcafp * N + cellI] = rate_ICaL_fcafp;
        RATES_SoA[ICaL_fcas * N + cellI] = rate_ICaL_fcas;
        RATES_SoA[ICaL_ff * N + cellI] = rate_ICaL_ff;
        RATES_SoA[ICaL_ffp * N + cellI] = rate_ICaL_ffp;
        RATES_SoA[ICaL_fs * N + cellI] = rate_ICaL_fs;
        RATES_SoA[ICaL_jca * N + cellI] = rate_ICaL_jca;
        RATES_SoA[ICaL_nca * N + cellI] = rate_ICaL_nca;
        RATES_SoA[ICaT_b * N + cellI] = rate_ICaT_b;
        RATES_SoA[ICaT_g * N + cellI] = rate_ICaT_g;
        RATES_SoA[IK1_xk1 * N + cellI] = rate_IK1_xk1;
        RATES_SoA[IKr_xrf * N + cellI] = rate_IKr_xrf;
        RATES_SoA[IKr_xrs * N + cellI] = rate_IKr_xrs;
        RATES_SoA[IKs_xs1 * N + cellI] = rate_IKs_xs1;
        RATES_SoA[IKs_xs2 * N + cellI] = rate_IKs_xs2;
        RATES_SoA[INaL_hL * N + cellI] = rate_INaL_hL;
        RATES_SoA[INaL_hLp * N + cellI] = rate_INaL_hLp;
        RATES_SoA[INaL_mL * N + cellI] = rate_INaL_mL;
        RATES_SoA[INa_hf * N + cellI] = rate_INa_hf;
        RATES_SoA[INa_hs * N + cellI] = rate_INa_hs;
        RATES_SoA[INa_hsp * N + cellI] = rate_INa_hsp;
        RATES_SoA[INa_j * N + cellI] = rate_INa_j;
        RATES_SoA[INa_jp * N + cellI] = rate_INa_jp;
        RATES_SoA[INa_m * N + cellI] = rate_INa_m;
        RATES_SoA[IP3_u * N + cellI] = rate_IP3_u;
        RATES_SoA[If_y * N + cellI] = rate_If_y;
        RATES_SoA[Ito_a * N + cellI] = rate_Ito_a;
        RATES_SoA[Ito_i1 * N + cellI] = rate_Ito_i1;
        RATES_SoA[Ito_i2 * N + cellI] = rate_Ito_i2;
        RATES_SoA[Ki * N + cellI] = rate_Ki;
        RATES_SoA[Ksl * N + cellI] = rate_Ksl;
        RATES_SoA[Kss * N + cellI] = rate_Kss;
        RATES_SoA[Nai * N + cellI] = rate_Nai;
        RATES_SoA[Nasl * N + cellI] = rate_Nasl;
        RATES_SoA[Nass * N + cellI] = rate_Nass;
        RATES_SoA[membrane_v * N + cellI] = rate_membrane_v;
        RATES_SoA[ryr_Jrel1 * N + cellI] = rate_ryr_Jrel1;
        RATES_SoA[ryr_Jrel2 * N + cellI] = rate_ryr_Jrel2;

        if (SUPPORT_SoA != nullptr)
        {
            SUPPORT_SoA[TROVATO_BATCH_SUPPORT_tau_ICaT_b * N + cellI] = alg_AV_taub;
            SUPPORT_SoA[TROVATO_BATCH_SUPPORT_tau_ICaT_g * N + cellI] = alg_AV_taug;
            SUPPORT_SoA[TROVATO_BATCH_SUPPORT_tau_IK1_xk1 * N + cellI] = alg_AV_txk1;
            SUPPORT_SoA[TROVATO_BATCH_SUPPORT_tau_IKr_xrf * N + cellI] = alg_AV_txrf;
            SUPPORT_SoA[TROVATO_BATCH_SUPPORT_tau_IKr_xrs * N + cellI] = alg_AV_txrs;
            SUPPORT_SoA[TROVATO_BATCH_SUPPORT_tau_IKs_xs1 * N + cellI] = alg_AV_txs1;
            SUPPORT_SoA[TROVATO_BATCH_SUPPORT_tau_IKs_xs2 * N + cellI] = alg_AV_txs2;
            SUPPORT_SoA[TROVATO_BATCH_SUPPORT_tau_INa_hf * N + cellI] = alg_AV_thf;
            SUPPORT_SoA[TROVATO_BATCH_SUPPORT_tau_INa_hs * N + cellI] = alg_AV_ths;
            SUPPORT_SoA[TROVATO_BATCH_SUPPORT_tau_INa_hsp * N + cellI] = alg_AV_thsp;
            SUPPORT_SoA[TROVATO_BATCH_SUPPORT_tau_INa_m * N + cellI] = alg_AV_tm;
            SUPPORT_SoA[TROVATO_BATCH_SUPPORT_tau_INa_j * N + cellI] = alg_AV_tj;
            SUPPORT_SoA[TROVATO_BATCH_SUPPORT_tau_INa_jp * N + cellI] = alg_AV_tjp;
            SUPPORT_SoA[TROVATO_BATCH_SUPPORT_tau_If_y * N + cellI] = alg_AV_tauy;
            SUPPORT_SoA[TROVATO_BATCH_SUPPORT_tau_Ito_a * N + cellI] = alg_AV_taua;
            SUPPORT_SoA[TROVATO_BATCH_SUPPORT_tau_Ito_i1 * N + cellI] = alg_AV_tauis;
            SUPPORT_SoA[TROVATO_BATCH_SUPPORT_tau_Ito_i2 * N + cellI] = alg_AV_tauif;
            SUPPORT_SoA[TROVATO_BATCH_SUPPORT_tau_INaL_hL * N + cellI] = CONSTANTS[AC_thL];
            SUPPORT_SoA[TROVATO_BATCH_SUPPORT_tau_INaL_mL * N + cellI] = alg_AV_tmL;
            SUPPORT_SoA[TROVATO_BATCH_SUPPORT_tau_INaL_hLp * N + cellI] = CONSTANTS[AC_thLp];
            SUPPORT_SoA[TROVATO_BATCH_SUPPORT_tau_ICaL_d * N + cellI] = alg_AV_td;
            SUPPORT_SoA[TROVATO_BATCH_SUPPORT_tau_ICaL_ff * N + cellI] = alg_AV_tff;
            SUPPORT_SoA[TROVATO_BATCH_SUPPORT_tau_ICaL_fs * N + cellI] = alg_AV_tfs;
            SUPPORT_SoA[TROVATO_BATCH_SUPPORT_tau_ICaL_fcaf * N + cellI] = alg_AV_tfcaf;
            SUPPORT_SoA[TROVATO_BATCH_SUPPORT_tau_ICaL_fcafp * N + cellI] = alg_AV_tfcafp;
            SUPPORT_SoA[TROVATO_BATCH_SUPPORT_tau_ICaL_fcas * N + cellI] = alg_AV_tfcas;
            SUPPORT_SoA[TROVATO_BATCH_SUPPORT_tau_ICaL_ffp * N + cellI] = alg_AV_tffp;
            SUPPORT_SoA[TROVATO_BATCH_SUPPORT_tau_ICaL_jca * N + cellI] = CONSTANTS[AC_tjca];
            SUPPORT_SoA[TROVATO_BATCH_SUPPORT_tau_ryr_Jrel1 * N + cellI] = alg_AV_ireltau;
            SUPPORT_SoA[TROVATO_BATCH_SUPPORT_tau_ryr_Jrel2 * N + cellI] = alg_AV_ireltau2;
            SUPPORT_SoA[TROVATO_BATCH_SUPPORT_gInf_ICaT_b * N + cellI] = alg_AV_bss;
            SUPPORT_SoA[TROVATO_BATCH_SUPPORT_gInf_ICaT_g * N + cellI] = alg_AV_gss;
            SUPPORT_SoA[TROVATO_BATCH_SUPPORT_gInf_IK1_xk1 * N + cellI] = alg_AV_xk1ss;
            SUPPORT_SoA[TROVATO_BATCH_SUPPORT_gInf_IKr_xrf * N + cellI] = alg_AV_xrss;
            SUPPORT_SoA[TROVATO_BATCH_SUPPORT_gInf_IKr_xrs * N + cellI] = alg_AV_xrss;
            SUPPORT_SoA[TROVATO_BATCH_SUPPORT_gInf_IKs_xs1 * N + cellI] = alg_AV_xs1ss;
            SUPPORT_SoA[TROVATO_BATCH_SUPPORT_gInf_IKs_xs2 * N + cellI] = alg_AV_xs2ss;
            SUPPORT_SoA[TROVATO_BATCH_SUPPORT_gInf_INa_hf * N + cellI] = alg_AV_hss;
            SUPPORT_SoA[TROVATO_BATCH_SUPPORT_gInf_INa_hs * N + cellI] = alg_AV_hss;
            SUPPORT_SoA[TROVATO_BATCH_SUPPORT_gInf_INa_hsp * N + cellI] = alg_AV_hssp;
            SUPPORT_SoA[TROVATO_BATCH_SUPPORT_gInf_INa_m * N + cellI] = alg_AV_mss;
            SUPPORT_SoA[TROVATO_BATCH_SUPPORT_gInf_INa_j * N + cellI] = alg_AV_jss;
            SUPPORT_SoA[TROVATO_BATCH_SUPPORT_gInf_INa_jp * N + cellI] = alg_AV_jss;
            SUPPORT_SoA[TROVATO_BATCH_SUPPORT_gInf_If_y * N + cellI] = alg_AV_yss;
            SUPPORT_SoA[TROVATO_BATCH_SUPPORT_gInf_Ito_a * N + cellI] = alg_AV_ass;
            SUPPORT_SoA[TROVATO_BATCH_SUPPORT_gInf_Ito_i1 * N + cellI] = alg_AV_iss;
            SUPPORT_SoA[TROVATO_BATCH_SUPPORT_gInf_Ito_i2 * N + cellI] = alg_AV_iss;
            SUPPORT_SoA[TROVATO_BATCH_SUPPORT_gInf_INaL_hL * N + cellI] = alg_AV_hLss;
            SUPPORT_SoA[TROVATO_BATCH_SUPPORT_gInf_INaL_mL * N + cellI] = alg_AV_mLss;
            SUPPORT_SoA[TROVATO_BATCH_SUPPORT_gInf_INaL_hLp * N + cellI] = alg_AV_hLssp;
            SUPPORT_SoA[TROVATO_BATCH_SUPPORT_gInf_ICaL_d * N + cellI] = alg_AV_dss;
            SUPPORT_SoA[TROVATO_BATCH_SUPPORT_gInf_ICaL_ff * N + cellI] = alg_AV_fss;
            SUPPORT_SoA[TROVATO_BATCH_SUPPORT_gInf_ICaL_fs * N + cellI] = alg_AV_fss;
            SUPPORT_SoA[TROVATO_BATCH_SUPPORT_gInf_ICaL_fcaf * N + cellI] = alg_AV_fcass;
            SUPPORT_SoA[TROVATO_BATCH_SUPPORT_gInf_ICaL_fcafp * N + cellI] = alg_AV_fcass;
            SUPPORT_SoA[TROVATO_BATCH_SUPPORT_gInf_ICaL_fcas * N + cellI] = alg_AV_fcass;
            SUPPORT_SoA[TROVATO_BATCH_SUPPORT_gInf_ICaL_ffp * N + cellI] = alg_AV_fss;
            SUPPORT_SoA[TROVATO_BATCH_SUPPORT_gInf_ICaL_jca * N + cellI] = alg_AV_fcass;
            SUPPORT_SoA[TROVATO_BATCH_SUPPORT_gInf_ryr_Jrel1 * N + cellI] = alg_AV_irelss;
            SUPPORT_SoA[TROVATO_BATCH_SUPPORT_gInf_ryr_Jrel2 * N + cellI] = alg_AV_irelss2;
            SUPPORT_SoA[TROVATO_BATCH_SUPPORT_Iion_cm * N + cellI] = alg_Iion_cm;
        }
    }
}

} // End namespace Foam

// ************************************************************************* //
