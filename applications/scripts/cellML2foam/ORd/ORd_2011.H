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

\*---------------------------------------------------------------------------*/

#include "ORd_2011Names.H"
#include "stimulusIO.H"

static const char* ORdSTATES_NAMES[NUM_STATES] = {
    "v",
    "CaMKt",
    "nai",
    "nass",
    "ki",
    "kss",
    "cass",
    "cansr",
    "cajsr",
    "cai",
    "m",
    "hf",
    "hs",
    "j",
    "hsp",
    "jp",
    "mL",
    "hL",
    "hLp",
    "a",
    "iF",
    "iS",
    "ap",
    "iFp",
    "iSp",
    "d",
    "ff",
    "fs",
    "fcaf",
    "fcas",
    "jca",
    "ffp",
    "fcafp",
    "nca",
    "xrf",
    "xrs",
    "xs1",
    "xs2",
    "xk1",
    "Jrelnp",
    "Jrelp",
};

static const char* ORdALGEBRAIC_NAMES[NUM_ALGEBRAIC] = {
    "AV_CaMKa",
    "AV_CaMKb",
    "AV_Afcaf",
    "AV_Afcas",
    "AV_ICaK",
    "AV_ICaL_ICaL",
    "AV_ICaNa",
    "AV_PhiCaK",
    "AV_PhiCaL",
    "AV_PhiCaNa",
    "AV_anca",
    "AV_dss",
    "AV_f",
    "AV_fICaLp",
    "AV_fca",
    "AV_fcap",
    "AV_fcass",
    "AV_fp",
    "AV_fss",
    "AV_km2n",
    "AV_td",
    "AV_tfcaf",
    "AV_tfcafp",
    "AV_tfcas",
    "AV_tff",
    "AV_tffp",
    "AV_tfs",
    "AV_ICab_ICab",
    "AV_IK1_IK1",
    "AV_rk1",
    "AV_txk1",
    "AV_xk1ss",
    "AV_IKb_IKb",
    "AV_xkb",
    "AV_Axrf",
    "AV_Axrs",
    "AV_IKr_IKr",
    "AV_rkr",
    "AV_txrf",
    "AV_txrs",
    "AV_xr",
    "AV_xrss",
    "AV_IKs_IKs",
    "AV_KsCa",
    "AV_txs1",
    "AV_txs2",
    "AV_xs1ss",
    "AV_xs2ss",
    "AV_INa_INa",
    "AV_fINap",
    "AV_h",
    "AV_hp",
    "AV_hss",
    "AV_hssp",
    "AV_jss",
    "AV_mss",
    "AV_thf",
    "AV_ths",
    "AV_thsp",
    "AV_tj",
    "AV_tjp",
    "AV_tm",
    "AV_E1_i",
    "AV_E1_ss",
    "AV_E2_i",
    "AV_E2_ss",
    "AV_E3_i",
    "AV_E3_ss",
    "AV_E4_i",
    "AV_E4_ss",
    "AV_INaCa_i_INaCa_i",
    "AV_INaCa_ss",
    "AV_JncxCa_i",
    "AV_JncxCa_ss",
    "AV_JncxNa_i",
    "AV_JncxNa_ss",
    "AV_allo_i",
    "AV_allo_ss",
    "AV_h1_i",
    "AV_h1_ss",
    "AV_h2_i",
    "AV_h2_ss",
    "AV_h3_i",
    "AV_h3_ss",
    "AV_h4_i",
    "AV_h4_ss",
    "AV_h5_i",
    "AV_h5_ss",
    "AV_h6_i",
    "AV_h6_ss",
    "AV_h7_i",
    "AV_h7_ss",
    "AV_h8_i",
    "AV_h8_ss",
    "AV_h9_i",
    "AV_h9_ss",
    "AV_hca",
    "AV_hna",
    "AV_k3_i",
    "AV_k3_ss",
    "AV_k3p_i",
    "AV_k3p_ss",
    "AV_k3pp_i",
    "AV_k3pp_ss",
    "AV_k4_i",
    "AV_k4_ss",
    "AV_k4p_i",
    "AV_k4p_ss",
    "AV_k4pp_i",
    "AV_k4pp_ss",
    "AV_k6_i",
    "AV_k6_ss",
    "AV_k7_i",
    "AV_k7_ss",
    "AV_k8_i",
    "AV_k8_ss",
    "AV_x1_i",
    "AV_x1_ss",
    "AV_x2_i",
    "AV_x2_ss",
    "AV_x3_i",
    "AV_x3_ss",
    "AV_x4_i",
    "AV_x4_ss",
    "AV_E1",
    "AV_E2",
    "AV_E3",
    "AV_E4",
    "AV_INaK_INaK",
    "AV_JnakK",
    "AV_JnakNa",
    "AV_Knai",
    "AV_Knao",
    "AV_P",
    "AV_a1",
    "AV_a3",
    "AV_b2",
    "AV_b3",
    "AV_b4",
    "AV_x1",
    "AV_x2",
    "AV_x3",
    "AV_x4",
    "AV_INaL_INaL",
    "AV_fINaLp",
    "AV_hLss",
    "AV_hLssp",
    "AV_mLss",
    "AV_tmL",
    "AV_INab_INab",
    "AV_IpCa_IpCa",
    "AV_AiF",
    "AV_AiS",
    "AV_Ito_Ito",
    "AV_ass",
    "AV_assp",
    "AV_delta_epi",
    "AV_dti_develop",
    "AV_dti_recover",
    "AV_fItop",
    "AV_i",
    "AV_ip",
    "AV_iss",
    "AV_ta",
    "AV_tiF",
    "AV_tiF_b",
    "AV_tiFp",
    "AV_tiS",
    "AV_tiS_b",
    "AV_tiSp",
    "AV_Jleak",
    "AV_Jup",
    "AV_Jupnp",
    "AV_Jupp",
    "AV_fJupp",
    "AV_Jdiff",
    "AV_JdiffK",
    "AV_JdiffNa",
    "AV_time",
    "AV_Bcai",
    "AV_Bcajsr",
    "AV_Bcass",
    "AV_Istim",
    "AV_vffrt",
    "AV_vfrt",
    "AV_EK",
    "AV_EKs",
    "AV_ENa",
    "AV_Jrel",
    "AV_Jrel_inf",
    "AV_Jrel_inf_temp",
    "AV_Jrel_infp",
    "AV_Jrel_temp",
    "AV_fJrelp",
    "AV_tau_rel",
    "AV_tau_rel_temp",
    "AV_tau_relp",
    "AV_tau_relp_temp",
    "AV_Jtr",
    "Istim",
    "Iion_cm",
};


static inline const Foam::HashTable<Foam::wordList>&
ORdDependencyMap()
{
   static const Foam::HashTable<Foam::wordList> dep
    (
        {
        }
    );

    return dep;
}

inline Foam::scalar computeIstim(Foam::scalar t, const Foam::StimulusProtocol& stimulus)
{
    return Foam::stimulusIO::computeStimulus(t, stimulus);
}
void
ORdinitConsts(double* CONSTANTS, double* RATES, double* STATES, int tissueFlag, const Foam::dictionary& stimulus)
{
/* CaMK */
    CONSTANTS[AC_CaMKo] = 0.05;
    CONSTANTS[AC_KmCaM] = 0.0015;
    CONSTANTS[AC_KmCaMK] = 0.15;
    CONSTANTS[AC_aCaMK] = 0.05;
    CONSTANTS[AC_bCaMK] = 0.00068;
    
    /* IpCa */
    CONSTANTS[AC_GpCa] = 0.0005;
    CONSTANTS[AC_KmCap] = 0.0005;
    
    /* cell_geometry */
    CONSTANTS[AC_L] = 0.01;
    CONSTANTS[AC_rad] = 0.0011;
    CONSTANTS[AC_Ageo] = 2.0 * 3.14 * CONSTANTS[AC_rad] * CONSTANTS[AC_rad] + 2.0 * 3.14 * CONSTANTS[AC_rad] * CONSTANTS[AC_L];
    CONSTANTS[AC_vcell] = 1000.0 * 3.14 * CONSTANTS[AC_rad] * CONSTANTS[AC_rad] * CONSTANTS[AC_L];
    CONSTANTS[AC_Acap] = 2.0 * CONSTANTS[AC_Ageo];
    CONSTANTS[AC_vjsr] = 0.0048 * CONSTANTS[AC_vcell];
    CONSTANTS[AC_vmyo] = 0.68 * CONSTANTS[AC_vcell];
    CONSTANTS[AC_vnsr] = 0.0552 * CONSTANTS[AC_vcell];
    CONSTANTS[AC_vss] = 0.02 * CONSTANTS[AC_vcell];
    
    /* environment */
    CONSTANTS[AC_celltype] = 0.0;
    
    /* extracellular */
    CONSTANTS[AC_cao] = 1.8;
    CONSTANTS[AC_ko] = 5.4;
    CONSTANTS[AC_nao] = 140.0;
    
    /* physical_constants */
    CONSTANTS[AC_F] = 96485.0;
    CONSTANTS[AC_R] = 8314.0;
    CONSTANTS[AC_T] = 310.0;
    CONSTANTS[AC_zca] = 2.0;
    CONSTANTS[AC_zk] = 1.0;
    CONSTANTS[AC_zna] = 1.0;
    
    /* INaCa_i */
    CONSTANTS[AC_Gncx_b] = 0.0008;
    CONSTANTS[AC_KmCaAct] = 0.00015;
    CONSTANTS[AC_kasymm] = 12.5;
    CONSTANTS[AC_kcaoff] = 5000.0;
    CONSTANTS[AC_kcaon] = 1500000.0;
    CONSTANTS[AC_kna1] = 15.0;
    CONSTANTS[AC_kna2] = 5.0;
    CONSTANTS[AC_kna3] = 88.12;
    CONSTANTS[AC_qca] = 0.167;
    CONSTANTS[AC_qna] = 0.5224;
    CONSTANTS[AC_wca] = 60000.0;
    CONSTANTS[AC_wna] = 60000.0;
    CONSTANTS[AC_wnaca] = 5000.0;
    CONSTANTS[AC_Gncx] = ((CONSTANTS[AC_celltype] == 1.0) ? CONSTANTS[AC_Gncx_b] * 1.1 : ((CONSTANTS[AC_celltype] == 2.0) ? CONSTANTS[AC_Gncx_b] * 1.4 : CONSTANTS[AC_Gncx_b]));
    CONSTANTS[AC_h10_i] = CONSTANTS[AC_kasymm] + 1.0 + CONSTANTS[AC_nao] / CONSTANTS[AC_kna1] * (1.0 + CONSTANTS[AC_nao] / CONSTANTS[AC_kna2]);
    CONSTANTS[AC_h10_ss] = CONSTANTS[AC_kasymm] + 1.0 + CONSTANTS[AC_nao] / CONSTANTS[AC_kna1] * (1.0 + CONSTANTS[AC_nao] / CONSTANTS[AC_kna2]);
    CONSTANTS[AC_k2_i] = CONSTANTS[AC_kcaoff];
    CONSTANTS[AC_k2_ss] = CONSTANTS[AC_kcaoff];
    CONSTANTS[AC_k5_i] = CONSTANTS[AC_kcaoff];
    CONSTANTS[AC_k5_ss] = CONSTANTS[AC_kcaoff];
    CONSTANTS[AC_h11_i] = CONSTANTS[AC_nao] * CONSTANTS[AC_nao] / (CONSTANTS[AC_h10_i] * CONSTANTS[AC_kna1] * CONSTANTS[AC_kna2]);
    CONSTANTS[AC_h11_ss] = CONSTANTS[AC_nao] * CONSTANTS[AC_nao] / (CONSTANTS[AC_h10_ss] * CONSTANTS[AC_kna1] * CONSTANTS[AC_kna2]);
    CONSTANTS[AC_h12_i] = 1.0 / CONSTANTS[AC_h10_i];
    CONSTANTS[AC_h12_ss] = 1.0 / CONSTANTS[AC_h10_ss];
    CONSTANTS[AC_k1_i] = CONSTANTS[AC_h12_i] * CONSTANTS[AC_cao] * CONSTANTS[AC_kcaon];
    CONSTANTS[AC_k1_ss] = CONSTANTS[AC_h12_ss] * CONSTANTS[AC_cao] * CONSTANTS[AC_kcaon];
    
    /* INaK */
    CONSTANTS[AC_H] = 1e-07;
    CONSTANTS[AC_Khp] = 1.698e-07;
    CONSTANTS[AC_Kki] = 0.5;
    CONSTANTS[AC_Kko] = 0.3582;
    CONSTANTS[AC_Kmgatp] = 1.698e-07;
    CONSTANTS[AC_Knai0] = 9.073;
    CONSTANTS[AC_Knao0] = 27.78;
    CONSTANTS[AC_Knap] = 224.0;
    CONSTANTS[AC_Kxkur] = 292.0;
    CONSTANTS[AC_MgADP] = 0.05;
    CONSTANTS[AC_MgATP] = 9.8;
    CONSTANTS[AC_Pnak_b] = 30.0;
    CONSTANTS[AC_delta] = -0.155;
    CONSTANTS[AC_eP] = 4.2;
    CONSTANTS[AC_k1m] = 182.4;
    CONSTANTS[AC_k1p] = 949.5;
    CONSTANTS[AC_k2m] = 39.4;
    CONSTANTS[AC_k2p] = 687.2;
    CONSTANTS[AC_k3m] = 79300.0;
    CONSTANTS[AC_k3p] = 1899.0;
    CONSTANTS[AC_k4m] = 40.0;
    CONSTANTS[AC_k4p] = 639.0;
    CONSTANTS[AC_Pnak] = ((CONSTANTS[AC_celltype] == 1.0) ? CONSTANTS[AC_Pnak_b] * 0.9 : ((CONSTANTS[AC_celltype] == 2.0) ? CONSTANTS[AC_Pnak_b] * 0.7 : CONSTANTS[AC_Pnak_b]));
    CONSTANTS[AC_a2] = CONSTANTS[AC_k2p];
    CONSTANTS[AC_a4] = CONSTANTS[AC_k4p] * CONSTANTS[AC_MgATP] / CONSTANTS[AC_Kmgatp] / (1.0 + CONSTANTS[AC_MgATP] / CONSTANTS[AC_Kmgatp]);
    CONSTANTS[AC_b1] = CONSTANTS[AC_k1m] * CONSTANTS[AC_MgADP];
    
    /* SERCA */
    CONSTANTS[AC_upScale] = ((CONSTANTS[AC_celltype] == 1.0) ? 1.3 : 1.0);
    
    /* reversal_potentials */
    CONSTANTS[AC_PKNa] = 0.01833;
    
    /* IK1 */
    CONSTANTS[AC_GK1_b] = 0.1908;
    CONSTANTS[AC_GK1] = ((CONSTANTS[AC_celltype] == 1.0) ? CONSTANTS[AC_GK1_b] * 1.2 : ((CONSTANTS[AC_celltype] == 2.0) ? CONSTANTS[AC_GK1_b] * 1.3 : CONSTANTS[AC_GK1_b]));
    
    /* IKb */
    CONSTANTS[AC_GKb_b] = 0.003;
    CONSTANTS[AC_GKb] = ((CONSTANTS[AC_celltype] == 1.0) ? CONSTANTS[AC_GKb_b] * 0.6 : CONSTANTS[AC_GKb_b]);
    
    /* IKr */
    CONSTANTS[AC_GKr_b] = 0.046;
    CONSTANTS[AC_GKr] = ((CONSTANTS[AC_celltype] == 1.0) ? CONSTANTS[AC_GKr_b] * 1.3 : ((CONSTANTS[AC_celltype] == 2.0) ? CONSTANTS[AC_GKr_b] * 0.8 : CONSTANTS[AC_GKr_b]));
    
    /* IKs */
    CONSTANTS[AC_GKs_b] = 0.0034;
    CONSTANTS[AC_GKs] = ((CONSTANTS[AC_celltype] == 1.0) ? CONSTANTS[AC_GKs_b] * 1.4 : CONSTANTS[AC_GKs_b]);
    
    /* INa */
    CONSTANTS[AC_Ahf] = 0.99;
    CONSTANTS[AC_GNa] = 75.0;
    CONSTANTS[AC_hssV1] = 82.9;
    CONSTANTS[AC_hssV2] = 6.086;
    CONSTANTS[AC_mssV1] = 39.57;
    CONSTANTS[AC_mssV2] = 9.871;
    CONSTANTS[AC_mtD1] = 6.765;
    CONSTANTS[AC_mtD2] = 8.552;
    CONSTANTS[AC_mtV1] = 11.64;
    CONSTANTS[AC_mtV2] = 34.77;
    CONSTANTS[AC_mtV3] = 77.42;
    CONSTANTS[AC_mtV4] = 5.955;
    CONSTANTS[AC_Ahs] = 1.0 - CONSTANTS[AC_Ahf];
    
    /* Ito */
    CONSTANTS[AC_Gto_b] = 0.02;
    CONSTANTS[AC_Gto] = ((CONSTANTS[AC_celltype] == 1.0) ? CONSTANTS[AC_Gto_b] * 4.0 : ((CONSTANTS[AC_celltype] == 2.0) ? CONSTANTS[AC_Gto_b] * 4.0 : CONSTANTS[AC_Gto_b]));
    
    /* INaL */
    CONSTANTS[AC_GNaL_b] = 0.0075;
    CONSTANTS[AC_thL] = 200.0;
    CONSTANTS[AC_GNaL] = ((CONSTANTS[AC_celltype] == 1.0) ? CONSTANTS[AC_GNaL_b] * 0.6 : CONSTANTS[AC_GNaL_b]);
    CONSTANTS[AC_thLp] = 3.0 * CONSTANTS[AC_thL];
    
    /* ICaL */
    CONSTANTS[AC_Aff] = 0.6;
    CONSTANTS[AC_Kmn] = 0.002;
    CONSTANTS[AC_PCa_b] = 0.0001;
    CONSTANTS[AC_k2n] = 1000.0;
    CONSTANTS[AC_tjca] = 75.0;
    CONSTANTS[AC_Afs] = 1.0 - CONSTANTS[AC_Aff];
    CONSTANTS[AC_PCa] = ((CONSTANTS[AC_celltype] == 1.0) ? CONSTANTS[AC_PCa_b] * 1.2 : ((CONSTANTS[AC_celltype] == 2.0) ? CONSTANTS[AC_PCa_b] * 2.5 : CONSTANTS[AC_PCa_b]));
    CONSTANTS[AC_PCaK] = 0.0003574 * CONSTANTS[AC_PCa];
    CONSTANTS[AC_PCaNa] = 0.00125 * CONSTANTS[AC_PCa];
    CONSTANTS[AC_PCap] = 1.1 * CONSTANTS[AC_PCa];
    CONSTANTS[AC_PCaKp] = 0.0003574 * CONSTANTS[AC_PCap];
    CONSTANTS[AC_PCaNap] = 0.00125 * CONSTANTS[AC_PCap];
    
    /* ICab */
    CONSTANTS[AC_PCab] = 2.5e-08;
    
    /* INab */
    CONSTANTS[AC_PNab] = 3.75e-10;
    
    /* intracellular_ions */
    CONSTANTS[AC_BSLmax] = 1.124;
    CONSTANTS[AC_BSRmax] = 0.047;
    CONSTANTS[AC_KmBSL] = 0.0087;
    CONSTANTS[AC_KmBSR] = 0.00087;
    CONSTANTS[AC_cm] = 1.0;
    CONSTANTS[AC_cmdnmax_b] = 0.05;
    CONSTANTS[AC_csqnmax] = 10.0;
    CONSTANTS[AC_kmcmdn] = 0.00238;
    CONSTANTS[AC_kmcsqn] = 0.8;
    CONSTANTS[AC_kmtrpn] = 0.0005;
    CONSTANTS[AC_trpnmax] = 0.07;
    CONSTANTS[AC_cmdnmax] = ((CONSTANTS[AC_celltype] == 1.0) ? CONSTANTS[AC_cmdnmax_b] * 1.3 : CONSTANTS[AC_cmdnmax_b]);
    
    /* membrane */
    CONSTANTS[AC_amp] = -80.0;
    CONSTANTS[AC_duration] = 0.5;
    CONSTANTS[AC_stimStart] = 50.0;
    
    /* ryr */
    CONSTANTS[AC_bt] = 4.75;
    CONSTANTS[AC_a_rel] = 0.5 * CONSTANTS[AC_bt];
    CONSTANTS[AC_btp] = 1.25 * CONSTANTS[AC_bt];
    CONSTANTS[AC_a_relp] = 0.5 * CONSTANTS[AC_btp];

    /* --- Initial values --- */

STATES[v] = -87.0;
    STATES[CaMKt] = 0.0;
    STATES[nai] = 7.0;
    STATES[nass] = 7.0;
    STATES[ki] = 145.0;
    STATES[kss] = 145.0;
    STATES[cass] = 0.0001;
    STATES[cansr] = 1.2;
    STATES[cajsr] = 1.2;
    STATES[cai] = 0.0001;
    STATES[m] = 0.0;
    STATES[hf] = 1.0;
    STATES[hs] = 1.0;
    STATES[j] = 1.0;
    STATES[hsp] = 1.0;
    STATES[jp] = 1.0;
    STATES[mL] = 0.0;
    STATES[hL] = 1.0;
    STATES[hLp] = 1.0;
    STATES[a] = 0.0;
    STATES[iF] = 1.0;
    STATES[iS] = 1.0;
    STATES[ap] = 0.0;
    STATES[iFp] = 1.0;
    STATES[iSp] = 1.0;
    STATES[d] = 0.0;
    STATES[ff] = 1.0;
    STATES[fs] = 1.0;
    STATES[fcaf] = 1.0;
    STATES[fcas] = 1.0;
    STATES[jca] = 1.0;
    STATES[ffp] = 1.0;
    STATES[fcafp] = 1.0;
    STATES[nca] = 0.0;
    STATES[xrf] = 0.0;
    STATES[xrs] = 0.0;
    STATES[xs1] = 0.0;
    STATES[xs2] = 0.0;
    STATES[xk1] = 1.0;
    STATES[Jrelnp] = 0.0;
    STATES[Jrelp] = 0.0;
}


void
ORdcomputeVariables(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC, int tissueFlag, bool solveVmWithinODESolver, const Foam::StimulusProtocol& stimulus)
{
/* CaMK */
    ALGEBRAIC[AV_CaMKb] = CONSTANTS[AC_CaMKo] * (1.0 - STATES[CaMKt]) / (1.0 + CONSTANTS[AC_KmCaM] / STATES[cass]);
    ALGEBRAIC[AV_CaMKa] = ALGEBRAIC[AV_CaMKb] + STATES[CaMKt];
    RATES[CaMKt] = CONSTANTS[AC_aCaMK] * ALGEBRAIC[AV_CaMKb] * (ALGEBRAIC[AV_CaMKb] + STATES[CaMKt]) - CONSTANTS[AC_bCaMK] * STATES[CaMKt];
    
    /* IpCa */
    ALGEBRAIC[AV_IpCa_IpCa] = CONSTANTS[AC_GpCa] * STATES[cai] / (CONSTANTS[AC_KmCap] + STATES[cai]);
    
    /* diff */
    ALGEBRAIC[AV_Jdiff] = (STATES[cass] - STATES[cai]) / 0.2;
    ALGEBRAIC[AV_JdiffK] = (STATES[kss] - STATES[ki]) / 2.0;
    ALGEBRAIC[AV_JdiffNa] = (STATES[nass] - STATES[nai]) / 2.0;
    
    /* environment */
    ALGEBRAIC[AV_time] = t;
    
    /* trans_flux */
    ALGEBRAIC[AV_Jtr] = (STATES[cansr] - STATES[cajsr]) / 100.0;
    
    /* INaCa_i */
    ALGEBRAIC[AV_allo_i] = 1.0 / (1.0 + pow(CONSTANTS[AC_KmCaAct] / STATES[cai], 2.0));
    ALGEBRAIC[AV_allo_ss] = 1.0 / (1.0 + pow(CONSTANTS[AC_KmCaAct] / STATES[cass], 2.0));
    ALGEBRAIC[AV_h4_i] = 1.0 + STATES[nai] / CONSTANTS[AC_kna1] * (1.0 + STATES[nai] / CONSTANTS[AC_kna2]);
    ALGEBRAIC[AV_h4_ss] = 1.0 + STATES[nass] / CONSTANTS[AC_kna1] * (1.0 + STATES[nass] / CONSTANTS[AC_kna2]);
    ALGEBRAIC[AV_hca] = exp(CONSTANTS[AC_qca] * STATES[v] * CONSTANTS[AC_F] / (CONSTANTS[AC_R] * CONSTANTS[AC_T]));
    ALGEBRAIC[AV_hna] = exp(CONSTANTS[AC_qna] * STATES[v] * CONSTANTS[AC_F] / (CONSTANTS[AC_R] * CONSTANTS[AC_T]));
    ALGEBRAIC[AV_h1_i] = 1.0 + STATES[nai] / CONSTANTS[AC_kna3] * (1.0 + ALGEBRAIC[AV_hna]);
    ALGEBRAIC[AV_h1_ss] = 1.0 + STATES[nass] / CONSTANTS[AC_kna3] * (1.0 + ALGEBRAIC[AV_hna]);
    ALGEBRAIC[AV_h5_i] = STATES[nai] * STATES[nai] / (ALGEBRAIC[AV_h4_i] * CONSTANTS[AC_kna1] * CONSTANTS[AC_kna2]);
    ALGEBRAIC[AV_h5_ss] = STATES[nass] * STATES[nass] / (ALGEBRAIC[AV_h4_ss] * CONSTANTS[AC_kna1] * CONSTANTS[AC_kna2]);
    ALGEBRAIC[AV_h6_i] = 1.0 / ALGEBRAIC[AV_h4_i];
    ALGEBRAIC[AV_h6_ss] = 1.0 / ALGEBRAIC[AV_h4_ss];
    ALGEBRAIC[AV_h7_i] = 1.0 + CONSTANTS[AC_nao] / CONSTANTS[AC_kna3] * (1.0 + 1.0 / ALGEBRAIC[AV_hna]);
    ALGEBRAIC[AV_h7_ss] = 1.0 + CONSTANTS[AC_nao] / CONSTANTS[AC_kna3] * (1.0 + 1.0 / ALGEBRAIC[AV_hna]);
    ALGEBRAIC[AV_h2_i] = STATES[nai] * ALGEBRAIC[AV_hna] / (CONSTANTS[AC_kna3] * ALGEBRAIC[AV_h1_i]);
    ALGEBRAIC[AV_h2_ss] = STATES[nass] * ALGEBRAIC[AV_hna] / (CONSTANTS[AC_kna3] * ALGEBRAIC[AV_h1_ss]);
    ALGEBRAIC[AV_h3_i] = 1.0 / ALGEBRAIC[AV_h1_i];
    ALGEBRAIC[AV_h3_ss] = 1.0 / ALGEBRAIC[AV_h1_ss];
    ALGEBRAIC[AV_h8_i] = CONSTANTS[AC_nao] / (CONSTANTS[AC_kna3] * ALGEBRAIC[AV_hna] * ALGEBRAIC[AV_h7_i]);
    ALGEBRAIC[AV_h8_ss] = CONSTANTS[AC_nao] / (CONSTANTS[AC_kna3] * ALGEBRAIC[AV_hna] * ALGEBRAIC[AV_h7_ss]);
    ALGEBRAIC[AV_h9_i] = 1.0 / ALGEBRAIC[AV_h7_i];
    ALGEBRAIC[AV_h9_ss] = 1.0 / ALGEBRAIC[AV_h7_ss];
    ALGEBRAIC[AV_k6_i] = ALGEBRAIC[AV_h6_i] * STATES[cai] * CONSTANTS[AC_kcaon];
    ALGEBRAIC[AV_k6_ss] = ALGEBRAIC[AV_h6_ss] * STATES[cass] * CONSTANTS[AC_kcaon];
    ALGEBRAIC[AV_k3p_i] = ALGEBRAIC[AV_h9_i] * CONSTANTS[AC_wca];
    ALGEBRAIC[AV_k3p_ss] = ALGEBRAIC[AV_h9_ss] * CONSTANTS[AC_wca];
    ALGEBRAIC[AV_k3pp_i] = ALGEBRAIC[AV_h8_i] * CONSTANTS[AC_wnaca];
    ALGEBRAIC[AV_k3pp_ss] = ALGEBRAIC[AV_h8_ss] * CONSTANTS[AC_wnaca];
    ALGEBRAIC[AV_k4p_i] = ALGEBRAIC[AV_h3_i] * CONSTANTS[AC_wca] / ALGEBRAIC[AV_hca];
    ALGEBRAIC[AV_k4p_ss] = ALGEBRAIC[AV_h3_ss] * CONSTANTS[AC_wca] / ALGEBRAIC[AV_hca];
    ALGEBRAIC[AV_k4pp_i] = ALGEBRAIC[AV_h2_i] * CONSTANTS[AC_wnaca];
    ALGEBRAIC[AV_k4pp_ss] = ALGEBRAIC[AV_h2_ss] * CONSTANTS[AC_wnaca];
    ALGEBRAIC[AV_k7_i] = ALGEBRAIC[AV_h5_i] * ALGEBRAIC[AV_h2_i] * CONSTANTS[AC_wna];
    ALGEBRAIC[AV_k7_ss] = ALGEBRAIC[AV_h5_ss] * ALGEBRAIC[AV_h2_ss] * CONSTANTS[AC_wna];
    ALGEBRAIC[AV_k8_i] = ALGEBRAIC[AV_h8_i] * CONSTANTS[AC_h11_i] * CONSTANTS[AC_wna];
    ALGEBRAIC[AV_k8_ss] = ALGEBRAIC[AV_h8_ss] * CONSTANTS[AC_h11_ss] * CONSTANTS[AC_wna];
    ALGEBRAIC[AV_k3_i] = ALGEBRAIC[AV_k3p_i] + ALGEBRAIC[AV_k3pp_i];
    ALGEBRAIC[AV_k3_ss] = ALGEBRAIC[AV_k3p_ss] + ALGEBRAIC[AV_k3pp_ss];
    ALGEBRAIC[AV_k4_i] = ALGEBRAIC[AV_k4p_i] + ALGEBRAIC[AV_k4pp_i];
    ALGEBRAIC[AV_k4_ss] = ALGEBRAIC[AV_k4p_ss] + ALGEBRAIC[AV_k4pp_ss];
    ALGEBRAIC[AV_x1_i] = CONSTANTS[AC_k2_i] * ALGEBRAIC[AV_k4_i] * (ALGEBRAIC[AV_k7_i] + ALGEBRAIC[AV_k6_i]) + CONSTANTS[AC_k5_i] * ALGEBRAIC[AV_k7_i] * (CONSTANTS[AC_k2_i] + ALGEBRAIC[AV_k3_i]);
    ALGEBRAIC[AV_x1_ss] = CONSTANTS[AC_k2_ss] * ALGEBRAIC[AV_k4_ss] * (ALGEBRAIC[AV_k7_ss] + ALGEBRAIC[AV_k6_ss]) + CONSTANTS[AC_k5_ss] * ALGEBRAIC[AV_k7_ss] * (CONSTANTS[AC_k2_ss] + ALGEBRAIC[AV_k3_ss]);
    ALGEBRAIC[AV_x2_i] = CONSTANTS[AC_k1_i] * ALGEBRAIC[AV_k7_i] * (ALGEBRAIC[AV_k4_i] + CONSTANTS[AC_k5_i]) + ALGEBRAIC[AV_k4_i] * ALGEBRAIC[AV_k6_i] * (CONSTANTS[AC_k1_i] + ALGEBRAIC[AV_k8_i]);
    ALGEBRAIC[AV_x2_ss] = CONSTANTS[AC_k1_ss] * ALGEBRAIC[AV_k7_ss] * (ALGEBRAIC[AV_k4_ss] + CONSTANTS[AC_k5_ss]) + ALGEBRAIC[AV_k4_ss] * ALGEBRAIC[AV_k6_ss] * (CONSTANTS[AC_k1_ss] + ALGEBRAIC[AV_k8_ss]);
    ALGEBRAIC[AV_x3_i] = CONSTANTS[AC_k1_i] * ALGEBRAIC[AV_k3_i] * (ALGEBRAIC[AV_k7_i] + ALGEBRAIC[AV_k6_i]) + ALGEBRAIC[AV_k8_i] * ALGEBRAIC[AV_k6_i] * (CONSTANTS[AC_k2_i] + ALGEBRAIC[AV_k3_i]);
    ALGEBRAIC[AV_x3_ss] = CONSTANTS[AC_k1_ss] * ALGEBRAIC[AV_k3_ss] * (ALGEBRAIC[AV_k7_ss] + ALGEBRAIC[AV_k6_ss]) + ALGEBRAIC[AV_k8_ss] * ALGEBRAIC[AV_k6_ss] * (CONSTANTS[AC_k2_ss] + ALGEBRAIC[AV_k3_ss]);
    ALGEBRAIC[AV_x4_i] = CONSTANTS[AC_k2_i] * ALGEBRAIC[AV_k8_i] * (ALGEBRAIC[AV_k4_i] + CONSTANTS[AC_k5_i]) + ALGEBRAIC[AV_k3_i] * CONSTANTS[AC_k5_i] * (CONSTANTS[AC_k1_i] + ALGEBRAIC[AV_k8_i]);
    ALGEBRAIC[AV_x4_ss] = CONSTANTS[AC_k2_ss] * ALGEBRAIC[AV_k8_ss] * (ALGEBRAIC[AV_k4_ss] + CONSTANTS[AC_k5_ss]) + ALGEBRAIC[AV_k3_ss] * CONSTANTS[AC_k5_ss] * (CONSTANTS[AC_k1_ss] + ALGEBRAIC[AV_k8_ss]);
    ALGEBRAIC[AV_E1_i] = ALGEBRAIC[AV_x1_i] / (ALGEBRAIC[AV_x1_i] + ALGEBRAIC[AV_x2_i] + ALGEBRAIC[AV_x3_i] + ALGEBRAIC[AV_x4_i]);
    ALGEBRAIC[AV_E1_ss] = ALGEBRAIC[AV_x1_ss] / (ALGEBRAIC[AV_x1_ss] + ALGEBRAIC[AV_x2_ss] + ALGEBRAIC[AV_x3_ss] + ALGEBRAIC[AV_x4_ss]);
    ALGEBRAIC[AV_E2_i] = ALGEBRAIC[AV_x2_i] / (ALGEBRAIC[AV_x1_i] + ALGEBRAIC[AV_x2_i] + ALGEBRAIC[AV_x3_i] + ALGEBRAIC[AV_x4_i]);
    ALGEBRAIC[AV_E2_ss] = ALGEBRAIC[AV_x2_ss] / (ALGEBRAIC[AV_x1_ss] + ALGEBRAIC[AV_x2_ss] + ALGEBRAIC[AV_x3_ss] + ALGEBRAIC[AV_x4_ss]);
    ALGEBRAIC[AV_E3_i] = ALGEBRAIC[AV_x3_i] / (ALGEBRAIC[AV_x1_i] + ALGEBRAIC[AV_x2_i] + ALGEBRAIC[AV_x3_i] + ALGEBRAIC[AV_x4_i]);
    ALGEBRAIC[AV_E3_ss] = ALGEBRAIC[AV_x3_ss] / (ALGEBRAIC[AV_x1_ss] + ALGEBRAIC[AV_x2_ss] + ALGEBRAIC[AV_x3_ss] + ALGEBRAIC[AV_x4_ss]);
    ALGEBRAIC[AV_E4_i] = ALGEBRAIC[AV_x4_i] / (ALGEBRAIC[AV_x1_i] + ALGEBRAIC[AV_x2_i] + ALGEBRAIC[AV_x3_i] + ALGEBRAIC[AV_x4_i]);
    ALGEBRAIC[AV_E4_ss] = ALGEBRAIC[AV_x4_ss] / (ALGEBRAIC[AV_x1_ss] + ALGEBRAIC[AV_x2_ss] + ALGEBRAIC[AV_x3_ss] + ALGEBRAIC[AV_x4_ss]);
    ALGEBRAIC[AV_JncxCa_i] = ALGEBRAIC[AV_E2_i] * CONSTANTS[AC_k2_i] - ALGEBRAIC[AV_E1_i] * CONSTANTS[AC_k1_i];
    ALGEBRAIC[AV_JncxCa_ss] = ALGEBRAIC[AV_E2_ss] * CONSTANTS[AC_k2_ss] - ALGEBRAIC[AV_E1_ss] * CONSTANTS[AC_k1_ss];
    ALGEBRAIC[AV_JncxNa_i] = 3.0 * (ALGEBRAIC[AV_E4_i] * ALGEBRAIC[AV_k7_i] - ALGEBRAIC[AV_E1_i] * ALGEBRAIC[AV_k8_i]) + ALGEBRAIC[AV_E3_i] * ALGEBRAIC[AV_k4pp_i] - ALGEBRAIC[AV_E2_i] * ALGEBRAIC[AV_k3pp_i];
    ALGEBRAIC[AV_JncxNa_ss] = 3.0 * (ALGEBRAIC[AV_E4_ss] * ALGEBRAIC[AV_k7_ss] - ALGEBRAIC[AV_E1_ss] * ALGEBRAIC[AV_k8_ss]) + ALGEBRAIC[AV_E3_ss] * ALGEBRAIC[AV_k4pp_ss] - ALGEBRAIC[AV_E2_ss] * ALGEBRAIC[AV_k3pp_ss];
    ALGEBRAIC[AV_INaCa_i_INaCa_i] = 0.8 * CONSTANTS[AC_Gncx] * ALGEBRAIC[AV_allo_i] * (CONSTANTS[AC_zna] * ALGEBRAIC[AV_JncxNa_i] + CONSTANTS[AC_zca] * ALGEBRAIC[AV_JncxCa_i]);
    ALGEBRAIC[AV_INaCa_ss] = 0.2 * CONSTANTS[AC_Gncx] * ALGEBRAIC[AV_allo_ss] * (CONSTANTS[AC_zna] * ALGEBRAIC[AV_JncxNa_ss] + CONSTANTS[AC_zca] * ALGEBRAIC[AV_JncxCa_ss]);
    
    /* INaK */
    ALGEBRAIC[AV_Knai] = CONSTANTS[AC_Knai0] * exp(CONSTANTS[AC_delta] * STATES[v] * CONSTANTS[AC_F] / (3.0 * CONSTANTS[AC_R] * CONSTANTS[AC_T]));
    ALGEBRAIC[AV_Knao] = CONSTANTS[AC_Knao0] * exp((1.0 - CONSTANTS[AC_delta]) * STATES[v] * CONSTANTS[AC_F] / (3.0 * CONSTANTS[AC_R] * CONSTANTS[AC_T]));
    ALGEBRAIC[AV_P] = CONSTANTS[AC_eP] / (1.0 + CONSTANTS[AC_H] / CONSTANTS[AC_Khp] + STATES[nai] / CONSTANTS[AC_Knap] + STATES[ki] / CONSTANTS[AC_Kxkur]);
    ALGEBRAIC[AV_a1] = CONSTANTS[AC_k1p] * pow(STATES[nai] / ALGEBRAIC[AV_Knai], 3.0) / (pow(1.0 + STATES[nai] / ALGEBRAIC[AV_Knai], 3.0) + pow(1.0 + STATES[ki] / CONSTANTS[AC_Kki], 2.0) - 1.0);
    ALGEBRAIC[AV_a3] = CONSTANTS[AC_k3p] * pow(CONSTANTS[AC_ko] / CONSTANTS[AC_Kko], 2.0) / (pow(1.0 + CONSTANTS[AC_nao] / ALGEBRAIC[AV_Knao], 3.0) + pow(1.0 + CONSTANTS[AC_ko] / CONSTANTS[AC_Kko], 2.0) - 1.0);
    ALGEBRAIC[AV_b2] = CONSTANTS[AC_k2m] * pow(CONSTANTS[AC_nao] / ALGEBRAIC[AV_Knao], 3.0) / (pow(1.0 + CONSTANTS[AC_nao] / ALGEBRAIC[AV_Knao], 3.0) + pow(1.0 + CONSTANTS[AC_ko] / CONSTANTS[AC_Kko], 2.0) - 1.0);
    ALGEBRAIC[AV_b3] = CONSTANTS[AC_k3m] * ALGEBRAIC[AV_P] * CONSTANTS[AC_H] / (1.0 + CONSTANTS[AC_MgATP] / CONSTANTS[AC_Kmgatp]);
    ALGEBRAIC[AV_b4] = CONSTANTS[AC_k4m] * pow(STATES[ki] / CONSTANTS[AC_Kki], 2.0) / (pow(1.0 + STATES[nai] / ALGEBRAIC[AV_Knai], 3.0) + pow(1.0 + STATES[ki] / CONSTANTS[AC_Kki], 2.0) - 1.0);
    ALGEBRAIC[AV_x1] = CONSTANTS[AC_a4] * ALGEBRAIC[AV_a1] * CONSTANTS[AC_a2] + ALGEBRAIC[AV_b2] * ALGEBRAIC[AV_b4] * ALGEBRAIC[AV_b3] + CONSTANTS[AC_a2] * ALGEBRAIC[AV_b4] * ALGEBRAIC[AV_b3] + ALGEBRAIC[AV_b3] * ALGEBRAIC[AV_a1] * CONSTANTS[AC_a2];
    ALGEBRAIC[AV_x2] = ALGEBRAIC[AV_b2] * CONSTANTS[AC_b1] * ALGEBRAIC[AV_b4] + ALGEBRAIC[AV_a1] * CONSTANTS[AC_a2] * ALGEBRAIC[AV_a3] + ALGEBRAIC[AV_a3] * CONSTANTS[AC_b1] * ALGEBRAIC[AV_b4] + CONSTANTS[AC_a2] * ALGEBRAIC[AV_a3] * ALGEBRAIC[AV_b4];
    ALGEBRAIC[AV_x3] = CONSTANTS[AC_a2] * ALGEBRAIC[AV_a3] * CONSTANTS[AC_a4] + ALGEBRAIC[AV_b3] * ALGEBRAIC[AV_b2] * CONSTANTS[AC_b1] + ALGEBRAIC[AV_b2] * CONSTANTS[AC_b1] * CONSTANTS[AC_a4] + ALGEBRAIC[AV_a3] * CONSTANTS[AC_a4] * CONSTANTS[AC_b1];
    ALGEBRAIC[AV_x4] = ALGEBRAIC[AV_b4] * ALGEBRAIC[AV_b3] * ALGEBRAIC[AV_b2] + ALGEBRAIC[AV_a3] * CONSTANTS[AC_a4] * ALGEBRAIC[AV_a1] + ALGEBRAIC[AV_b2] * CONSTANTS[AC_a4] * ALGEBRAIC[AV_a1] + ALGEBRAIC[AV_b3] * ALGEBRAIC[AV_b2] * ALGEBRAIC[AV_a1];
    ALGEBRAIC[AV_E1] = ALGEBRAIC[AV_x1] / (ALGEBRAIC[AV_x1] + ALGEBRAIC[AV_x2] + ALGEBRAIC[AV_x3] + ALGEBRAIC[AV_x4]);
    ALGEBRAIC[AV_E2] = ALGEBRAIC[AV_x2] / (ALGEBRAIC[AV_x1] + ALGEBRAIC[AV_x2] + ALGEBRAIC[AV_x3] + ALGEBRAIC[AV_x4]);
    ALGEBRAIC[AV_E3] = ALGEBRAIC[AV_x3] / (ALGEBRAIC[AV_x1] + ALGEBRAIC[AV_x2] + ALGEBRAIC[AV_x3] + ALGEBRAIC[AV_x4]);
    ALGEBRAIC[AV_E4] = ALGEBRAIC[AV_x4] / (ALGEBRAIC[AV_x1] + ALGEBRAIC[AV_x2] + ALGEBRAIC[AV_x3] + ALGEBRAIC[AV_x4]);
    ALGEBRAIC[AV_JnakK] = 2.0 * (ALGEBRAIC[AV_E4] * CONSTANTS[AC_b1] - ALGEBRAIC[AV_E3] * ALGEBRAIC[AV_a1]);
    ALGEBRAIC[AV_JnakNa] = 3.0 * (ALGEBRAIC[AV_E1] * ALGEBRAIC[AV_a3] - ALGEBRAIC[AV_E2] * ALGEBRAIC[AV_b3]);
    ALGEBRAIC[AV_INaK_INaK] = CONSTANTS[AC_Pnak] * (CONSTANTS[AC_zna] * ALGEBRAIC[AV_JnakNa] + CONSTANTS[AC_zk] * ALGEBRAIC[AV_JnakK]);
    
    /* SERCA */
    ALGEBRAIC[AV_Jleak] = 0.0039375 * STATES[cansr] / 15.0;
    ALGEBRAIC[AV_fJupp] = 1.0 / (1.0 + CONSTANTS[AC_KmCaMK] / ALGEBRAIC[AV_CaMKa]);
    ALGEBRAIC[AV_Jupnp] = CONSTANTS[AC_upScale] * 0.004375 * STATES[cai] / (STATES[cai] + 0.00092);
    ALGEBRAIC[AV_Jupp] = CONSTANTS[AC_upScale] * 2.75 * 0.004375 * STATES[cai] / (STATES[cai] + 0.00092 - 0.00017);
    ALGEBRAIC[AV_Jup] = (1.0 - ALGEBRAIC[AV_fJupp]) * ALGEBRAIC[AV_Jupnp] + ALGEBRAIC[AV_fJupp] * ALGEBRAIC[AV_Jupp] - ALGEBRAIC[AV_Jleak];
    
    /* reversal_potentials */
    ALGEBRAIC[AV_EK] = CONSTANTS[AC_R] * CONSTANTS[AC_T] / CONSTANTS[AC_F] * log(CONSTANTS[AC_ko] / STATES[ki]);
    ALGEBRAIC[AV_ENa] = CONSTANTS[AC_R] * CONSTANTS[AC_T] / CONSTANTS[AC_F] * log(CONSTANTS[AC_nao] / STATES[nai]);
    ALGEBRAIC[AV_EKs] = CONSTANTS[AC_R] * CONSTANTS[AC_T] / CONSTANTS[AC_F] * log((CONSTANTS[AC_ko] + CONSTANTS[AC_PKNa] * CONSTANTS[AC_nao]) / (STATES[ki] + CONSTANTS[AC_PKNa] * STATES[nai]));
    
    /* IK1 */
    ALGEBRAIC[AV_rk1] = 1.0 / (1.0 + exp((STATES[v] + 105.8 - 2.6 * CONSTANTS[AC_ko]) / 9.493));
    ALGEBRAIC[AV_txk1] = 122.2 / (exp(-(STATES[v] + 127.2) / 20.36) + exp((STATES[v] + 236.8) / 69.33));
    ALGEBRAIC[AV_xk1ss] = 1.0 / (1.0 + exp(-(STATES[v] + 2.5538 * CONSTANTS[AC_ko] + 144.59) / (1.5692 * CONSTANTS[AC_ko] + 3.8115)));
    RATES[xk1] = (ALGEBRAIC[AV_xk1ss] - STATES[xk1]) / ALGEBRAIC[AV_txk1];
    ALGEBRAIC[AV_IK1_IK1] = CONSTANTS[AC_GK1] * sqrt(CONSTANTS[AC_ko]) * ALGEBRAIC[AV_rk1] * STATES[xk1] * (STATES[v] - ALGEBRAIC[AV_EK]);
    
    /* IKb */
    ALGEBRAIC[AV_xkb] = 1.0 / (1.0 + exp(-(STATES[v] - 14.48) / 18.34));
    ALGEBRAIC[AV_IKb_IKb] = CONSTANTS[AC_GKb] * ALGEBRAIC[AV_xkb] * (STATES[v] - ALGEBRAIC[AV_EK]);
    
    /* IKr */
    ALGEBRAIC[AV_Axrf] = 1.0 / (1.0 + exp((STATES[v] + 54.81) / 38.21));
    ALGEBRAIC[AV_rkr] = 1.0 / (1.0 + exp((STATES[v] + 55.0) / 75.0)) * 1.0 / (1.0 + exp((STATES[v] - 10.0) / 30.0));
    ALGEBRAIC[AV_txrf] = 12.98 + 1.0 / (0.3652 * exp((STATES[v] - 31.66) / 3.869) + 4.123e-05 * exp(-(STATES[v] - 47.78) / 20.38));
    ALGEBRAIC[AV_txrs] = 1.865 + 1.0 / (0.06629 * exp((STATES[v] - 34.7) / 7.355) + 1.128e-05 * exp(-(STATES[v] - 29.74) / 25.94));
    ALGEBRAIC[AV_xrss] = 1.0 / (1.0 + exp(-(STATES[v] + 8.337) / 6.789));
    ALGEBRAIC[AV_Axrs] = 1.0 - ALGEBRAIC[AV_Axrf];
    RATES[xrf] = (ALGEBRAIC[AV_xrss] - STATES[xrf]) / ALGEBRAIC[AV_txrf];
    RATES[xrs] = (ALGEBRAIC[AV_xrss] - STATES[xrs]) / ALGEBRAIC[AV_txrs];
    ALGEBRAIC[AV_xr] = ALGEBRAIC[AV_Axrf] * STATES[xrf] + ALGEBRAIC[AV_Axrs] * STATES[xrs];
    ALGEBRAIC[AV_IKr_IKr] = CONSTANTS[AC_GKr] * sqrt(CONSTANTS[AC_ko] / 5.4) * ALGEBRAIC[AV_xr] * ALGEBRAIC[AV_rkr] * (STATES[v] - ALGEBRAIC[AV_EK]);
    
    /* IKs */
    ALGEBRAIC[AV_KsCa] = 1.0 + 0.6 / (1.0 + pow(3.8e-05 / STATES[cai], 1.4));
    ALGEBRAIC[AV_txs1] = 817.3 + 1.0 / (0.0002326 * exp((STATES[v] + 48.28) / 17.8) + 0.001292 * exp(-(STATES[v] + 210.0) / 230.0));
    ALGEBRAIC[AV_txs2] = 1.0 / (0.01 * exp((STATES[v] - 50.0) / 20.0) + 0.0193 * exp(-(STATES[v] + 66.54) / 31.0));
    ALGEBRAIC[AV_xs1ss] = 1.0 / (1.0 + exp(-(STATES[v] + 11.6) / 8.932));
    ALGEBRAIC[AV_xs2ss] = ALGEBRAIC[AV_xs1ss];
    RATES[xs1] = (ALGEBRAIC[AV_xs1ss] - STATES[xs1]) / ALGEBRAIC[AV_txs1];
    ALGEBRAIC[AV_IKs_IKs] = CONSTANTS[AC_GKs] * ALGEBRAIC[AV_KsCa] * STATES[xs1] * STATES[xs2] * (STATES[v] - ALGEBRAIC[AV_EKs]);
    RATES[xs2] = (ALGEBRAIC[AV_xs2ss] - STATES[xs2]) / ALGEBRAIC[AV_txs2];
    
    /* INa */
    ALGEBRAIC[AV_fINap] = 1.0 / (1.0 + CONSTANTS[AC_KmCaMK] / ALGEBRAIC[AV_CaMKa]);
    ALGEBRAIC[AV_hssp] = 1.0 / (1.0 + exp((STATES[v] + 89.1) / 6.086));
    ALGEBRAIC[AV_thf] = 1.0 / (1.432e-05 * exp(-(STATES[v] + 1.196) / 6.285) + 6.149 * exp((STATES[v] + 0.5096) / 20.27));
    ALGEBRAIC[AV_ths] = 1.0 / (0.009794 * exp(-(STATES[v] + 17.95) / 28.05) + 0.3343 * exp((STATES[v] + 5.73) / 56.66));
    ALGEBRAIC[AV_tj] = 2.038 + 1.0 / (0.02136 * exp(-(STATES[v] + 100.6) / 8.281) + 0.3052 * exp((STATES[v] + 0.9941) / 38.45));
    ALGEBRAIC[AV_hss] = 1.0 / (1.0 + exp((STATES[v] + CONSTANTS[AC_hssV1]) / CONSTANTS[AC_hssV2]));
    ALGEBRAIC[AV_mss] = 1.0 / (1.0 + exp(-(STATES[v] + CONSTANTS[AC_mssV1]) / CONSTANTS[AC_mssV2]));
    ALGEBRAIC[AV_thsp] = 3.0 * ALGEBRAIC[AV_ths];
    ALGEBRAIC[AV_tjp] = 1.46 * ALGEBRAIC[AV_tj];
    ALGEBRAIC[AV_tm] = 1.0 / (CONSTANTS[AC_mtD1] * exp((STATES[v] + CONSTANTS[AC_mtV1]) / CONSTANTS[AC_mtV2]) + CONSTANTS[AC_mtD2] * exp(-(STATES[v] + CONSTANTS[AC_mtV3]) / CONSTANTS[AC_mtV4]));
    ALGEBRAIC[AV_h] = CONSTANTS[AC_Ahf] * STATES[hf] + CONSTANTS[AC_Ahs] * STATES[hs];
    ALGEBRAIC[AV_hp] = CONSTANTS[AC_Ahf] * STATES[hf] + CONSTANTS[AC_Ahs] * STATES[hsp];
    ALGEBRAIC[AV_jss] = ALGEBRAIC[AV_hss];
    RATES[hf] = (ALGEBRAIC[AV_hss] - STATES[hf]) / ALGEBRAIC[AV_thf];
    RATES[hs] = (ALGEBRAIC[AV_hss] - STATES[hs]) / ALGEBRAIC[AV_ths];
    RATES[hsp] = (ALGEBRAIC[AV_hssp] - STATES[hsp]) / ALGEBRAIC[AV_thsp];
    RATES[m] = (ALGEBRAIC[AV_mss] - STATES[m]) / ALGEBRAIC[AV_tm];
    ALGEBRAIC[AV_INa_INa] = CONSTANTS[AC_GNa] * (STATES[v] - ALGEBRAIC[AV_ENa]) * pow(STATES[m], 3.0) * ((1.0 - ALGEBRAIC[AV_fINap]) * ALGEBRAIC[AV_h] * STATES[j] + ALGEBRAIC[AV_fINap] * ALGEBRAIC[AV_hp] * STATES[jp]);
    RATES[j] = (ALGEBRAIC[AV_jss] - STATES[j]) / ALGEBRAIC[AV_tj];
    RATES[jp] = (ALGEBRAIC[AV_jss] - STATES[jp]) / ALGEBRAIC[AV_tjp];
    
    /* Ito */
    ALGEBRAIC[AV_AiF] = 1.0 / (1.0 + exp((STATES[v] - 213.6) / 151.2));
    ALGEBRAIC[AV_ass] = 1.0 / (1.0 + exp(-(STATES[v] - 14.34) / 14.82));
    ALGEBRAIC[AV_assp] = 1.0 / (1.0 + exp(-(STATES[v] - 24.34) / 14.82));
    ALGEBRAIC[AV_delta_epi] = ((CONSTANTS[AC_celltype] == 1.0) ? 1.0 - 0.95 / (1.0 + exp((STATES[v] + 70.0) / 5.0)) : 1.0);
    ALGEBRAIC[AV_dti_develop] = 1.354 + 0.0001 / (exp((STATES[v] - 167.4) / 15.89) + exp(-(STATES[v] - 12.23) / 0.2154));
    ALGEBRAIC[AV_dti_recover] = 1.0 - 0.5 / (1.0 + exp((STATES[v] + 70.0) / 20.0));
    ALGEBRAIC[AV_fItop] = 1.0 / (1.0 + CONSTANTS[AC_KmCaMK] / ALGEBRAIC[AV_CaMKa]);
    ALGEBRAIC[AV_iss] = 1.0 / (1.0 + exp((STATES[v] + 43.94) / 5.711));
    ALGEBRAIC[AV_ta] = 1.0515 / (1.0 / (1.2089 * (1.0 + exp(-(STATES[v] - 18.4099) / 29.3814))) + 3.5 / (1.0 + exp((STATES[v] + 100.0) / 29.3814)));
    ALGEBRAIC[AV_tiF_b] = 4.562 + 1.0 / (0.3933 * exp(-(STATES[v] + 100.0) / 100.0) + 0.08004 * exp((STATES[v] + 50.0) / 16.59));
    ALGEBRAIC[AV_tiS_b] = 23.62 + 1.0 / (0.001416 * exp(-(STATES[v] + 96.52) / 59.05) + 1.78e-08 * exp((STATES[v] + 114.1) / 8.079));
    ALGEBRAIC[AV_AiS] = 1.0 - ALGEBRAIC[AV_AiF];
    ALGEBRAIC[AV_tiF] = ALGEBRAIC[AV_tiF_b] * ALGEBRAIC[AV_delta_epi];
    ALGEBRAIC[AV_tiS] = ALGEBRAIC[AV_tiS_b] * ALGEBRAIC[AV_delta_epi];
    RATES[a] = (ALGEBRAIC[AV_ass] - STATES[a]) / ALGEBRAIC[AV_ta];
    RATES[ap] = (ALGEBRAIC[AV_assp] - STATES[ap]) / ALGEBRAIC[AV_ta];
    ALGEBRAIC[AV_i] = ALGEBRAIC[AV_AiF] * STATES[iF] + ALGEBRAIC[AV_AiS] * STATES[iS];
    ALGEBRAIC[AV_ip] = ALGEBRAIC[AV_AiF] * STATES[iFp] + ALGEBRAIC[AV_AiS] * STATES[iSp];
    ALGEBRAIC[AV_tiFp] = ALGEBRAIC[AV_dti_develop] * ALGEBRAIC[AV_dti_recover] * ALGEBRAIC[AV_tiF];
    ALGEBRAIC[AV_tiSp] = ALGEBRAIC[AV_dti_develop] * ALGEBRAIC[AV_dti_recover] * ALGEBRAIC[AV_tiS];
    RATES[iF] = (ALGEBRAIC[AV_iss] - STATES[iF]) / ALGEBRAIC[AV_tiF];
    RATES[iS] = (ALGEBRAIC[AV_iss] - STATES[iS]) / ALGEBRAIC[AV_tiS];
    ALGEBRAIC[AV_Ito_Ito] = CONSTANTS[AC_Gto] * (STATES[v] - ALGEBRAIC[AV_EK]) * ((1.0 - ALGEBRAIC[AV_fItop]) * STATES[a] * ALGEBRAIC[AV_i] + ALGEBRAIC[AV_fItop] * STATES[ap] * ALGEBRAIC[AV_ip]);
    RATES[iFp] = (ALGEBRAIC[AV_iss] - STATES[iFp]) / ALGEBRAIC[AV_tiFp];
    RATES[iSp] = (ALGEBRAIC[AV_iss] - STATES[iSp]) / ALGEBRAIC[AV_tiSp];
    
    /* INaL */
    ALGEBRAIC[AV_fINaLp] = 1.0 / (1.0 + CONSTANTS[AC_KmCaMK] / ALGEBRAIC[AV_CaMKa]);
    ALGEBRAIC[AV_hLss] = 1.0 / (1.0 + exp((STATES[v] + 87.61) / 7.488));
    ALGEBRAIC[AV_hLssp] = 1.0 / (1.0 + exp((STATES[v] + 93.81) / 7.488));
    ALGEBRAIC[AV_mLss] = 1.0 / (1.0 + exp(-(STATES[v] + 42.85) / 5.264));
    ALGEBRAIC[AV_tmL] = ALGEBRAIC[AV_tm];
    RATES[hL] = (ALGEBRAIC[AV_hLss] - STATES[hL]) / CONSTANTS[AC_thL];
    RATES[mL] = (ALGEBRAIC[AV_mLss] - STATES[mL]) / ALGEBRAIC[AV_tmL];
    ALGEBRAIC[AV_INaL_INaL] = CONSTANTS[AC_GNaL] * (STATES[v] - ALGEBRAIC[AV_ENa]) * STATES[mL] * ((1.0 - ALGEBRAIC[AV_fINaLp]) * STATES[hL] + ALGEBRAIC[AV_fINaLp] * STATES[hLp]);
    RATES[hLp] = (ALGEBRAIC[AV_hLssp] - STATES[hLp]) / CONSTANTS[AC_thLp];
    
    /* ICaL */
    ALGEBRAIC[AV_Afcaf] = 0.3 + 0.6 / (1.0 + exp((STATES[v] - 10.0) / 10.0));
    ALGEBRAIC[AV_dss] = 1.0 / (1.0 + exp(-(STATES[v] + 3.94) / 4.23));
    ALGEBRAIC[AV_fICaLp] = 1.0 / (1.0 + CONSTANTS[AC_KmCaMK] / ALGEBRAIC[AV_CaMKa]);
    ALGEBRAIC[AV_fss] = 1.0 / (1.0 + exp((STATES[v] + 19.58) / 3.696));
    ALGEBRAIC[AV_km2n] = STATES[jca] * 1.0;
    ALGEBRAIC[AV_td] = 0.6 + 1.0 / (exp(-0.05 * (STATES[v] + 6.0)) + exp(0.09 * (STATES[v] + 14.0)));
    ALGEBRAIC[AV_tfcaf] = 7.0 + 1.0 / (0.04 * exp(-(STATES[v] - 4.0) / 7.0) + 0.04 * exp((STATES[v] - 4.0) / 7.0));
    ALGEBRAIC[AV_tfcas] = 100.0 + 1.0 / (0.00012 * exp(-STATES[v] / 3.0) + 0.00012 * exp(STATES[v] / 7.0));
    ALGEBRAIC[AV_tff] = 7.0 + 1.0 / (0.0045 * exp(-(STATES[v] + 20.0) / 10.0) + 0.0045 * exp((STATES[v] + 20.0) / 10.0));
    ALGEBRAIC[AV_tfs] = 1000.0 + 1.0 / (3.5e-05 * exp(-(STATES[v] + 5.0) / 4.0) + 3.5e-05 * exp((STATES[v] + 5.0) / 6.0));
    ALGEBRAIC[AV_Afcas] = 1.0 - ALGEBRAIC[AV_Afcaf];
    ALGEBRAIC[AV_anca] = 1.0 / (CONSTANTS[AC_k2n] / ALGEBRAIC[AV_km2n] + pow(1.0 + CONSTANTS[AC_Kmn] / STATES[cass], 4.0));
    ALGEBRAIC[AV_fcass] = ALGEBRAIC[AV_fss];
    ALGEBRAIC[AV_tfcafp] = 2.5 * ALGEBRAIC[AV_tfcaf];
    ALGEBRAIC[AV_tffp] = 2.5 * ALGEBRAIC[AV_tff];
    RATES[d] = (ALGEBRAIC[AV_dss] - STATES[d]) / ALGEBRAIC[AV_td];
    RATES[ff] = (ALGEBRAIC[AV_fss] - STATES[ff]) / ALGEBRAIC[AV_tff];
    RATES[fs] = (ALGEBRAIC[AV_fss] - STATES[fs]) / ALGEBRAIC[AV_tfs];
    ALGEBRAIC[AV_f] = CONSTANTS[AC_Aff] * STATES[ff] + CONSTANTS[AC_Afs] * STATES[fs];
    ALGEBRAIC[AV_fca] = ALGEBRAIC[AV_Afcaf] * STATES[fcaf] + ALGEBRAIC[AV_Afcas] * STATES[fcas];
    ALGEBRAIC[AV_fcap] = ALGEBRAIC[AV_Afcaf] * STATES[fcafp] + ALGEBRAIC[AV_Afcas] * STATES[fcas];
    ALGEBRAIC[AV_fp] = CONSTANTS[AC_Aff] * STATES[ffp] + CONSTANTS[AC_Afs] * STATES[fs];
    RATES[fcaf] = (ALGEBRAIC[AV_fcass] - STATES[fcaf]) / ALGEBRAIC[AV_tfcaf];
    RATES[fcafp] = (ALGEBRAIC[AV_fcass] - STATES[fcafp]) / ALGEBRAIC[AV_tfcafp];
    RATES[fcas] = (ALGEBRAIC[AV_fcass] - STATES[fcas]) / ALGEBRAIC[AV_tfcas];
    RATES[ffp] = (ALGEBRAIC[AV_fss] - STATES[ffp]) / ALGEBRAIC[AV_tffp];
    RATES[jca] = (ALGEBRAIC[AV_fcass] - STATES[jca]) / CONSTANTS[AC_tjca];
    RATES[nca] = ALGEBRAIC[AV_anca] * CONSTANTS[AC_k2n] - STATES[nca] * ALGEBRAIC[AV_km2n];
    
    /* intracellular_ions */
    RATES[cansr] = ALGEBRAIC[AV_Jup] - ALGEBRAIC[AV_Jtr] * CONSTANTS[AC_vjsr] / CONSTANTS[AC_vnsr];
    ALGEBRAIC[AV_Bcajsr] = 1.0 / (1.0 + CONSTANTS[AC_csqnmax] * CONSTANTS[AC_kmcsqn] / pow(CONSTANTS[AC_kmcsqn] + STATES[cajsr], 2.0));
    ALGEBRAIC[AV_Bcass] = 1.0 / (1.0 + CONSTANTS[AC_BSRmax] * CONSTANTS[AC_KmBSR] / pow(CONSTANTS[AC_KmBSR] + STATES[cass], 2.0) + CONSTANTS[AC_BSLmax] * CONSTANTS[AC_KmBSL] / pow(CONSTANTS[AC_KmBSL] + STATES[cass], 2.0));
    ALGEBRAIC[AV_Bcai] = 1.0 / (1.0 + CONSTANTS[AC_cmdnmax] * CONSTANTS[AC_kmcmdn] / pow(CONSTANTS[AC_kmcmdn] + STATES[cai], 2.0) + CONSTANTS[AC_trpnmax] * CONSTANTS[AC_kmtrpn] / pow(CONSTANTS[AC_kmtrpn] + STATES[cai], 2.0));
    
    /* membrane */
    ALGEBRAIC[AV_vffrt] = STATES[v] * CONSTANTS[AC_F] * CONSTANTS[AC_F] / (CONSTANTS[AC_R] * CONSTANTS[AC_T]);
    ALGEBRAIC[AV_vfrt] = STATES[v] * CONSTANTS[AC_F] / (CONSTANTS[AC_R] * CONSTANTS[AC_T]);
    ALGEBRAIC[AV_Istim] = (((ALGEBRAIC[AV_time] > CONSTANTS[AC_stimStart]) && (ALGEBRAIC[AV_time] <= CONSTANTS[AC_stimStart] + CONSTANTS[AC_duration])) ? CONSTANTS[AC_amp] : 0.0);
    
    /* ryr */
    ALGEBRAIC[AV_fJrelp] = 1.0 / (1.0 + CONSTANTS[AC_KmCaMK] / ALGEBRAIC[AV_CaMKa]);
    ALGEBRAIC[AV_Jrel] = (1.0 - ALGEBRAIC[AV_fJrelp]) * STATES[Jrelnp] + ALGEBRAIC[AV_fJrelp] * STATES[Jrelp];
    ALGEBRAIC[AV_tau_rel_temp] = CONSTANTS[AC_bt] / (1.0 + 0.0123 / STATES[cajsr]);
    ALGEBRAIC[AV_tau_rel] = ((ALGEBRAIC[AV_tau_rel_temp] < 0.001) ? 0.001 : ALGEBRAIC[AV_tau_rel_temp]);
    ALGEBRAIC[AV_tau_relp_temp] = CONSTANTS[AC_btp] / (1.0 + 0.0123 / STATES[cajsr]);
    ALGEBRAIC[AV_tau_relp] = ((ALGEBRAIC[AV_tau_relp_temp] < 0.001) ? 0.001 : ALGEBRAIC[AV_tau_relp_temp]);
    
    /* *remaining* */
    ALGEBRAIC[AV_PhiCaK] = 1.0 * ALGEBRAIC[AV_vffrt] * (0.75 * STATES[kss] * exp(1.0 * ALGEBRAIC[AV_vfrt]) - 0.75 * CONSTANTS[AC_ko]) / (exp(1.0 * ALGEBRAIC[AV_vfrt]) - 1.0);
    ALGEBRAIC[AV_PhiCaL] = 4.0 * ALGEBRAIC[AV_vffrt] * (STATES[cass] * exp(2.0 * ALGEBRAIC[AV_vfrt]) - 0.341 * CONSTANTS[AC_cao]) / (exp(2.0 * ALGEBRAIC[AV_vfrt]) - 1.0);
    ALGEBRAIC[AV_PhiCaNa] = 1.0 * ALGEBRAIC[AV_vffrt] * (0.75 * STATES[nass] * exp(1.0 * ALGEBRAIC[AV_vfrt]) - 0.75 * CONSTANTS[AC_nao]) / (exp(1.0 * ALGEBRAIC[AV_vfrt]) - 1.0);
    ALGEBRAIC[AV_ICab_ICab] = CONSTANTS[AC_PCab] * 4.0 * ALGEBRAIC[AV_vffrt] * (STATES[cai] * exp(2.0 * ALGEBRAIC[AV_vfrt]) - 0.341 * CONSTANTS[AC_cao]) / (exp(2.0 * ALGEBRAIC[AV_vfrt]) - 1.0);
    ALGEBRAIC[AV_INab_INab] = CONSTANTS[AC_PNab] * ALGEBRAIC[AV_vffrt] * (STATES[nai] * exp(ALGEBRAIC[AV_vfrt]) - CONSTANTS[AC_nao]) / (exp(ALGEBRAIC[AV_vfrt]) - 1.0);
    RATES[cajsr] = ALGEBRAIC[AV_Bcajsr] * (ALGEBRAIC[AV_Jtr] - ALGEBRAIC[AV_Jrel]);
    RATES[ki] = -(ALGEBRAIC[AV_Ito_Ito] + ALGEBRAIC[AV_IKr_IKr] + ALGEBRAIC[AV_IKs_IKs] + ALGEBRAIC[AV_IK1_IK1] + ALGEBRAIC[AV_IKb_IKb] + ALGEBRAIC[AV_Istim] - 2.0 * ALGEBRAIC[AV_INaK_INaK]) * CONSTANTS[AC_cm] * CONSTANTS[AC_Acap] / (CONSTANTS[AC_F] * CONSTANTS[AC_vmyo]) + ALGEBRAIC[AV_JdiffK] * CONSTANTS[AC_vss] / CONSTANTS[AC_vmyo];
    ALGEBRAIC[AV_ICaK] = (1.0 - ALGEBRAIC[AV_fICaLp]) * CONSTANTS[AC_PCaK] * ALGEBRAIC[AV_PhiCaK] * STATES[d] * (ALGEBRAIC[AV_f] * (1.0 - STATES[nca]) + STATES[jca] * ALGEBRAIC[AV_fca] * STATES[nca]) + ALGEBRAIC[AV_fICaLp] * CONSTANTS[AC_PCaKp] * ALGEBRAIC[AV_PhiCaK] * STATES[d] * (ALGEBRAIC[AV_fp] * (1.0 - STATES[nca]) + STATES[jca] * ALGEBRAIC[AV_fcap] * STATES[nca]);
    ALGEBRAIC[AV_ICaL_ICaL] = (1.0 - ALGEBRAIC[AV_fICaLp]) * CONSTANTS[AC_PCa] * ALGEBRAIC[AV_PhiCaL] * STATES[d] * (ALGEBRAIC[AV_f] * (1.0 - STATES[nca]) + STATES[jca] * ALGEBRAIC[AV_fca] * STATES[nca]) + ALGEBRAIC[AV_fICaLp] * CONSTANTS[AC_PCap] * ALGEBRAIC[AV_PhiCaL] * STATES[d] * (ALGEBRAIC[AV_fp] * (1.0 - STATES[nca]) + STATES[jca] * ALGEBRAIC[AV_fcap] * STATES[nca]);
    ALGEBRAIC[AV_ICaNa] = (1.0 - ALGEBRAIC[AV_fICaLp]) * CONSTANTS[AC_PCaNa] * ALGEBRAIC[AV_PhiCaNa] * STATES[d] * (ALGEBRAIC[AV_f] * (1.0 - STATES[nca]) + STATES[jca] * ALGEBRAIC[AV_fca] * STATES[nca]) + ALGEBRAIC[AV_fICaLp] * CONSTANTS[AC_PCaNap] * ALGEBRAIC[AV_PhiCaNa] * STATES[d] * (ALGEBRAIC[AV_fp] * (1.0 - STATES[nca]) + STATES[jca] * ALGEBRAIC[AV_fcap] * STATES[nca]);
    RATES[cai] = ALGEBRAIC[AV_Bcai] * (-(ALGEBRAIC[AV_IpCa_IpCa] + ALGEBRAIC[AV_ICab_ICab] - 2.0 * ALGEBRAIC[AV_INaCa_i_INaCa_i]) * CONSTANTS[AC_cm] * CONSTANTS[AC_Acap] / (2.0 * CONSTANTS[AC_F] * CONSTANTS[AC_vmyo]) - ALGEBRAIC[AV_Jup] * CONSTANTS[AC_vnsr] / CONSTANTS[AC_vmyo] + ALGEBRAIC[AV_Jdiff] * CONSTANTS[AC_vss] / CONSTANTS[AC_vmyo]);
    RATES[nai] = -(ALGEBRAIC[AV_INa_INa] + ALGEBRAIC[AV_INaL_INaL] + 3.0 * ALGEBRAIC[AV_INaCa_i_INaCa_i] + 3.0 * ALGEBRAIC[AV_INaK_INaK] + ALGEBRAIC[AV_INab_INab]) * CONSTANTS[AC_Acap] * CONSTANTS[AC_cm] / (CONSTANTS[AC_F] * CONSTANTS[AC_vmyo]) + ALGEBRAIC[AV_JdiffNa] * CONSTANTS[AC_vss] / CONSTANTS[AC_vmyo];
    RATES[cass] = ALGEBRAIC[AV_Bcass] * (-(ALGEBRAIC[AV_ICaL_ICaL] - 2.0 * ALGEBRAIC[AV_INaCa_ss]) * CONSTANTS[AC_cm] * CONSTANTS[AC_Acap] / (2.0 * CONSTANTS[AC_F] * CONSTANTS[AC_vss]) + ALGEBRAIC[AV_Jrel] * CONSTANTS[AC_vjsr] / CONSTANTS[AC_vss] - ALGEBRAIC[AV_Jdiff]);
    RATES[kss] = -ALGEBRAIC[AV_ICaK] * CONSTANTS[AC_cm] * CONSTANTS[AC_Acap] / (CONSTANTS[AC_F] * CONSTANTS[AC_vss]) - ALGEBRAIC[AV_JdiffK];
    RATES[nass] = -(ALGEBRAIC[AV_ICaNa] + 3.0 * ALGEBRAIC[AV_INaCa_ss]) * CONSTANTS[AC_cm] * CONSTANTS[AC_Acap] / (CONSTANTS[AC_F] * CONSTANTS[AC_vss]) - ALGEBRAIC[AV_JdiffNa];
    RATES[v] = -(ALGEBRAIC[AV_INa_INa] + ALGEBRAIC[AV_INaL_INaL] + ALGEBRAIC[AV_Ito_Ito] + ALGEBRAIC[AV_ICaL_ICaL] + ALGEBRAIC[AV_ICaNa] + ALGEBRAIC[AV_ICaK] + ALGEBRAIC[AV_IKr_IKr] + ALGEBRAIC[AV_IKs_IKs] + ALGEBRAIC[AV_IK1_IK1] + ALGEBRAIC[AV_INaCa_i_INaCa_i] + ALGEBRAIC[AV_INaCa_ss] + ALGEBRAIC[AV_INaK_INaK] + ALGEBRAIC[AV_INab_INab] + ALGEBRAIC[AV_IKb_IKb] + ALGEBRAIC[AV_IpCa_IpCa] + ALGEBRAIC[AV_ICab_ICab] + ALGEBRAIC[AV_Istim]);
    ALGEBRAIC[AV_Jrel_inf_temp] = CONSTANTS[AC_a_rel] * -ALGEBRAIC[AV_ICaL_ICaL] / (1.0 + 1.0 * pow(1.5 / STATES[cajsr], 8.0));
    ALGEBRAIC[AV_Jrel_temp] = CONSTANTS[AC_a_relp] * -ALGEBRAIC[AV_ICaL_ICaL] / (1.0 + pow(1.5 / STATES[cajsr], 8.0));
    ALGEBRAIC[AV_Jrel_inf] = ((CONSTANTS[AC_celltype] == 2.0) ? ALGEBRAIC[AV_Jrel_inf_temp] * 1.7 : ALGEBRAIC[AV_Jrel_inf_temp]);
    ALGEBRAIC[AV_Jrel_infp] = ((CONSTANTS[AC_celltype] == 2.0) ? ALGEBRAIC[AV_Jrel_temp] * 1.7 : ALGEBRAIC[AV_Jrel_temp]);
    RATES[Jrelnp] = (ALGEBRAIC[AV_Jrel_inf] - STATES[Jrelnp]) / ALGEBRAIC[AV_tau_rel];
    RATES[Jrelp] = (ALGEBRAIC[AV_Jrel_infp] - STATES[Jrelp]) / ALGEBRAIC[AV_tau_relp];
    

    

    ALGEBRAIC[Istim] = computeIstim(VOI, stimulus);

    if (solveVmWithinODESolver)
    {
        RATES[0] = -ALGEBRAIC[Iion_cm] - ALGEBRAIC[Istim];
    }

}
