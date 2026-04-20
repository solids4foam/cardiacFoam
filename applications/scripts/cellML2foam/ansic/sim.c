/*
Ohara_Rudy_2011
Generated on 2026-04-14 22:16:11

Compiling on GCC:
 $ gcc -Wall -lm -lsundials_nvecserial -lsundials_cvode sim.c

Gnuplot example:
set terminal pngcairo enhanced linewidth 2 size 1200, 800;
set output 'V.png'
set size 1.0, 1.0
set xlabel 'time [ms]';
set grid
plot 'V.txt' using 1:2 with lines ls 1 title 'Vm'

*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_config.h>
#if SUNDIALS_VERSION_MAJOR >= 3
  #include <sunmatrix/sunmatrix_dense.h>
  #include <sunlinsol/sunlinsol_dense.h>
  #include <cvodes/cvodes_direct.h>
#else
  #include <cvode/cvode_dense.h>
#endif

#define N_STATE 41

/* Declare intermediary, temporary and system variables */
static realtype t;
static realtype pace;
static realtype ALGEBRAIC[AV_CaMKa];
static realtype ALGEBRAIC[AV_CaMKb];
static realtype CONSTANTS[AC_CaMKo];
static realtype CONSTANTS[AC_KmCaM];
static realtype CONSTANTS[AC_KmCaMK];
static realtype CONSTANTS[AC_aCaMK];
static realtype CONSTANTS[AC_bCaMK];
static realtype ALGEBRAIC[AV_Afcaf];
static realtype ALGEBRAIC[AV_Afcas];
static realtype CONSTANTS[AC_Aff];
static realtype CONSTANTS[AC_Afs];
static realtype ALGEBRAIC[AV_ICaK];
static realtype ALGEBRAIC[AV_ICaL_ICaL];
static realtype ALGEBRAIC[AV_ICaNa];
static realtype CONSTANTS[AC_Kmn];
static realtype CONSTANTS[AC_PCa];
static realtype CONSTANTS[AC_PCaK];
static realtype CONSTANTS[AC_PCaKp];
static realtype CONSTANTS[AC_PCaNa];
static realtype CONSTANTS[AC_PCaNap];
static realtype CONSTANTS[AC_PCa_b];
static realtype CONSTANTS[AC_PCap];
static realtype ALGEBRAIC[AV_PhiCaK];
static realtype ALGEBRAIC[AV_PhiCaL];
static realtype ALGEBRAIC[AV_PhiCaNa];
static realtype ALGEBRAIC[AV_anca];
static realtype ALGEBRAIC[AV_dss];
static realtype ALGEBRAIC[AV_f];
static realtype ALGEBRAIC[AV_fICaLp];
static realtype ALGEBRAIC[AV_fca];
static realtype ALGEBRAIC[AV_fcap];
static realtype ALGEBRAIC[AV_fcass];
static realtype ALGEBRAIC[AV_fp];
static realtype ALGEBRAIC[AV_fss];
static realtype CONSTANTS[AC_k2n];
static realtype ALGEBRAIC[AV_km2n];
static realtype ALGEBRAIC[AV_td];
static realtype ALGEBRAIC[AV_tfcaf];
static realtype ALGEBRAIC[AV_tfcafp];
static realtype ALGEBRAIC[AV_tfcas];
static realtype ALGEBRAIC[AV_tff];
static realtype ALGEBRAIC[AV_tffp];
static realtype ALGEBRAIC[AV_tfs];
static realtype CONSTANTS[AC_tjca];
static realtype ALGEBRAIC[AV_ICab_ICab];
static realtype CONSTANTS[AC_PCab];
static realtype CONSTANTS[AC_GK1];
static realtype CONSTANTS[AC_GK1_b];
static realtype ALGEBRAIC[AV_IK1_IK1];
static realtype ALGEBRAIC[AV_rk1];
static realtype ALGEBRAIC[AV_txk1];
static realtype ALGEBRAIC[AV_xk1ss];
static realtype CONSTANTS[AC_GKb];
static realtype CONSTANTS[AC_GKb_b];
static realtype ALGEBRAIC[AV_IKb_IKb];
static realtype ALGEBRAIC[AV_xkb];
static realtype ALGEBRAIC[AV_Axrf];
static realtype ALGEBRAIC[AV_Axrs];
static realtype CONSTANTS[AC_GKr];
static realtype CONSTANTS[AC_GKr_b];
static realtype ALGEBRAIC[AV_IKr_IKr];
static realtype ALGEBRAIC[AV_rkr];
static realtype ALGEBRAIC[AV_txrf];
static realtype ALGEBRAIC[AV_txrs];
static realtype ALGEBRAIC[AV_xr];
static realtype ALGEBRAIC[AV_xrss];
static realtype CONSTANTS[AC_GKs];
static realtype CONSTANTS[AC_GKs_b];
static realtype ALGEBRAIC[AV_IKs_IKs];
static realtype ALGEBRAIC[AV_KsCa];
static realtype ALGEBRAIC[AV_txs1];
static realtype ALGEBRAIC[AV_txs2];
static realtype ALGEBRAIC[AV_xs1ss];
static realtype ALGEBRAIC[AV_xs2ss];
static realtype CONSTANTS[AC_Ahf];
static realtype CONSTANTS[AC_Ahs];
static realtype CONSTANTS[AC_GNa];
static realtype ALGEBRAIC[AV_INa_INa];
static realtype ALGEBRAIC[AV_fINap];
static realtype ALGEBRAIC[AV_h];
static realtype ALGEBRAIC[AV_hp];
static realtype ALGEBRAIC[AV_hss];
static realtype CONSTANTS[AC_hssV1];
static realtype CONSTANTS[AC_hssV2];
static realtype ALGEBRAIC[AV_hssp];
static realtype ALGEBRAIC[AV_jss];
static realtype ALGEBRAIC[AV_mss];
static realtype CONSTANTS[AC_mssV1];
static realtype CONSTANTS[AC_mssV2];
static realtype CONSTANTS[AC_mtD1];
static realtype CONSTANTS[AC_mtD2];
static realtype CONSTANTS[AC_mtV1];
static realtype CONSTANTS[AC_mtV2];
static realtype CONSTANTS[AC_mtV3];
static realtype CONSTANTS[AC_mtV4];
static realtype ALGEBRAIC[AV_thf];
static realtype ALGEBRAIC[AV_ths];
static realtype ALGEBRAIC[AV_thsp];
static realtype ALGEBRAIC[AV_tj];
static realtype ALGEBRAIC[AV_tjp];
static realtype ALGEBRAIC[AV_tm];
static realtype ALGEBRAIC[AV_E1_i];
static realtype ALGEBRAIC[AV_E1_ss];
static realtype ALGEBRAIC[AV_E2_i];
static realtype ALGEBRAIC[AV_E2_ss];
static realtype ALGEBRAIC[AV_E3_i];
static realtype ALGEBRAIC[AV_E3_ss];
static realtype ALGEBRAIC[AV_E4_i];
static realtype ALGEBRAIC[AV_E4_ss];
static realtype CONSTANTS[AC_Gncx];
static realtype CONSTANTS[AC_Gncx_b];
static realtype ALGEBRAIC[AV_INaCa_i_INaCa_i];
static realtype ALGEBRAIC[AV_INaCa_ss];
static realtype ALGEBRAIC[AV_JncxCa_i];
static realtype ALGEBRAIC[AV_JncxCa_ss];
static realtype ALGEBRAIC[AV_JncxNa_i];
static realtype ALGEBRAIC[AV_JncxNa_ss];
static realtype CONSTANTS[AC_KmCaAct];
static realtype ALGEBRAIC[AV_allo_i];
static realtype ALGEBRAIC[AV_allo_ss];
static realtype CONSTANTS[AC_h10_i];
static realtype CONSTANTS[AC_h10_ss];
static realtype CONSTANTS[AC_h11_i];
static realtype CONSTANTS[AC_h11_ss];
static realtype CONSTANTS[AC_h12_i];
static realtype CONSTANTS[AC_h12_ss];
static realtype ALGEBRAIC[AV_h1_i];
static realtype ALGEBRAIC[AV_h1_ss];
static realtype ALGEBRAIC[AV_h2_i];
static realtype ALGEBRAIC[AV_h2_ss];
static realtype ALGEBRAIC[AV_h3_i];
static realtype ALGEBRAIC[AV_h3_ss];
static realtype ALGEBRAIC[AV_h4_i];
static realtype ALGEBRAIC[AV_h4_ss];
static realtype ALGEBRAIC[AV_h5_i];
static realtype ALGEBRAIC[AV_h5_ss];
static realtype ALGEBRAIC[AV_h6_i];
static realtype ALGEBRAIC[AV_h6_ss];
static realtype ALGEBRAIC[AV_h7_i];
static realtype ALGEBRAIC[AV_h7_ss];
static realtype ALGEBRAIC[AV_h8_i];
static realtype ALGEBRAIC[AV_h8_ss];
static realtype ALGEBRAIC[AV_h9_i];
static realtype ALGEBRAIC[AV_h9_ss];
static realtype ALGEBRAIC[AV_hca];
static realtype ALGEBRAIC[AV_hna];
static realtype CONSTANTS[AC_k1_i];
static realtype CONSTANTS[AC_k1_ss];
static realtype CONSTANTS[AC_k2_i];
static realtype CONSTANTS[AC_k2_ss];
static realtype ALGEBRAIC[AV_k3_i];
static realtype ALGEBRAIC[AV_k3_ss];
static realtype ALGEBRAIC[AV_k3p_i];
static realtype ALGEBRAIC[AV_k3p_ss];
static realtype ALGEBRAIC[AV_k3pp_i];
static realtype ALGEBRAIC[AV_k3pp_ss];
static realtype ALGEBRAIC[AV_k4_i];
static realtype ALGEBRAIC[AV_k4_ss];
static realtype ALGEBRAIC[AV_k4p_i];
static realtype ALGEBRAIC[AV_k4p_ss];
static realtype ALGEBRAIC[AV_k4pp_i];
static realtype ALGEBRAIC[AV_k4pp_ss];
static realtype CONSTANTS[AC_k5_i];
static realtype CONSTANTS[AC_k5_ss];
static realtype ALGEBRAIC[AV_k6_i];
static realtype ALGEBRAIC[AV_k6_ss];
static realtype ALGEBRAIC[AV_k7_i];
static realtype ALGEBRAIC[AV_k7_ss];
static realtype ALGEBRAIC[AV_k8_i];
static realtype ALGEBRAIC[AV_k8_ss];
static realtype CONSTANTS[AC_kasymm];
static realtype CONSTANTS[AC_kcaoff];
static realtype CONSTANTS[AC_kcaon];
static realtype CONSTANTS[AC_kna1];
static realtype CONSTANTS[AC_kna2];
static realtype CONSTANTS[AC_kna3];
static realtype CONSTANTS[AC_qca];
static realtype CONSTANTS[AC_qna];
static realtype CONSTANTS[AC_wca];
static realtype CONSTANTS[AC_wna];
static realtype CONSTANTS[AC_wnaca];
static realtype ALGEBRAIC[AV_x1_i];
static realtype ALGEBRAIC[AV_x1_ss];
static realtype ALGEBRAIC[AV_x2_i];
static realtype ALGEBRAIC[AV_x2_ss];
static realtype ALGEBRAIC[AV_x3_i];
static realtype ALGEBRAIC[AV_x3_ss];
static realtype ALGEBRAIC[AV_x4_i];
static realtype ALGEBRAIC[AV_x4_ss];
static realtype ALGEBRAIC[AV_E1];
static realtype ALGEBRAIC[AV_E2];
static realtype ALGEBRAIC[AV_E3];
static realtype ALGEBRAIC[AV_E4];
static realtype CONSTANTS[AC_H];
static realtype ALGEBRAIC[AV_INaK_INaK];
static realtype ALGEBRAIC[AV_JnakK];
static realtype ALGEBRAIC[AV_JnakNa];
static realtype CONSTANTS[AC_Khp];
static realtype CONSTANTS[AC_Kki];
static realtype CONSTANTS[AC_Kko];
static realtype CONSTANTS[AC_Kmgatp];
static realtype ALGEBRAIC[AV_Knai];
static realtype CONSTANTS[AC_Knai0];
static realtype ALGEBRAIC[AV_Knao];
static realtype CONSTANTS[AC_Knao0];
static realtype CONSTANTS[AC_Knap];
static realtype CONSTANTS[AC_Kxkur];
static realtype CONSTANTS[AC_MgADP];
static realtype CONSTANTS[AC_MgATP];
static realtype ALGEBRAIC[AV_P];
static realtype CONSTANTS[AC_Pnak];
static realtype CONSTANTS[AC_Pnak_b];
static realtype ALGEBRAIC[AV_a1];
static realtype CONSTANTS[AC_a2];
static realtype ALGEBRAIC[AV_a3];
static realtype CONSTANTS[AC_a4];
static realtype CONSTANTS[AC_b1];
static realtype ALGEBRAIC[AV_b2];
static realtype ALGEBRAIC[AV_b3];
static realtype ALGEBRAIC[AV_b4];
static realtype CONSTANTS[AC_delta];
static realtype CONSTANTS[AC_eP];
static realtype CONSTANTS[AC_k1m];
static realtype CONSTANTS[AC_k1p];
static realtype CONSTANTS[AC_k2m];
static realtype CONSTANTS[AC_k2p];
static realtype CONSTANTS[AC_k3m];
static realtype CONSTANTS[AC_k3p];
static realtype CONSTANTS[AC_k4m];
static realtype CONSTANTS[AC_k4p];
static realtype ALGEBRAIC[AV_x1];
static realtype ALGEBRAIC[AV_x2];
static realtype ALGEBRAIC[AV_x3];
static realtype ALGEBRAIC[AV_x4];
static realtype CONSTANTS[AC_GNaL];
static realtype CONSTANTS[AC_GNaL_b];
static realtype ALGEBRAIC[AV_INaL_INaL];
static realtype ALGEBRAIC[AV_fINaLp];
static realtype ALGEBRAIC[AV_hLss];
static realtype ALGEBRAIC[AV_hLssp];
static realtype ALGEBRAIC[AV_mLss];
static realtype CONSTANTS[AC_thL];
static realtype CONSTANTS[AC_thLp];
static realtype ALGEBRAIC[AV_tmL];
static realtype ALGEBRAIC[AV_INab_INab];
static realtype CONSTANTS[AC_PNab];
static realtype CONSTANTS[AC_GpCa];
static realtype ALGEBRAIC[AV_IpCa_IpCa];
static realtype CONSTANTS[AC_KmCap];
static realtype ALGEBRAIC[AV_AiF];
static realtype ALGEBRAIC[AV_AiS];
static realtype CONSTANTS[AC_Gto];
static realtype CONSTANTS[AC_Gto_b];
static realtype ALGEBRAIC[AV_Ito_Ito];
static realtype ALGEBRAIC[AV_ass];
static realtype ALGEBRAIC[AV_assp];
static realtype ALGEBRAIC[AV_delta_epi];
static realtype ALGEBRAIC[AV_dti_develop];
static realtype ALGEBRAIC[AV_dti_recover];
static realtype ALGEBRAIC[AV_fItop];
static realtype ALGEBRAIC[AV_i];
static realtype ALGEBRAIC[AV_ip];
static realtype ALGEBRAIC[AV_iss];
static realtype ALGEBRAIC[AV_ta];
static realtype ALGEBRAIC[AV_tiF];
static realtype ALGEBRAIC[AV_tiF_b];
static realtype ALGEBRAIC[AV_tiFp];
static realtype ALGEBRAIC[AV_tiS];
static realtype ALGEBRAIC[AV_tiS_b];
static realtype ALGEBRAIC[AV_tiSp];
static realtype ALGEBRAIC[AV_Jleak];
static realtype ALGEBRAIC[AV_Jup];
static realtype ALGEBRAIC[AV_Jupnp];
static realtype ALGEBRAIC[AV_Jupp];
static realtype ALGEBRAIC[AV_fJupp];
static realtype CONSTANTS[AC_upScale];
static realtype CONSTANTS[AC_Acap];
static realtype CONSTANTS[AC_Ageo];
static realtype CONSTANTS[AC_L];
static realtype CONSTANTS[AC_rad];
static realtype CONSTANTS[AC_vcell];
static realtype CONSTANTS[AC_vjsr];
static realtype CONSTANTS[AC_vmyo];
static realtype CONSTANTS[AC_vnsr];
static realtype CONSTANTS[AC_vss];
static realtype ALGEBRAIC[AV_Jdiff];
static realtype ALGEBRAIC[AV_JdiffK];
static realtype ALGEBRAIC[AV_JdiffNa];
static realtype CONSTANTS[AC_celltype];
static realtype ALGEBRAIC[AV_time];
static realtype CONSTANTS[AC_cao];
static realtype CONSTANTS[AC_ko];
static realtype CONSTANTS[AC_nao];
static realtype CONSTANTS[AC_BSLmax];
static realtype CONSTANTS[AC_BSRmax];
static realtype ALGEBRAIC[AV_Bcai];
static realtype ALGEBRAIC[AV_Bcajsr];
static realtype ALGEBRAIC[AV_Bcass];
static realtype CONSTANTS[AC_KmBSL];
static realtype CONSTANTS[AC_KmBSR];
static realtype CONSTANTS[AC_cm];
static realtype CONSTANTS[AC_cmdnmax];
static realtype CONSTANTS[AC_cmdnmax_b];
static realtype CONSTANTS[AC_csqnmax];
static realtype CONSTANTS[AC_kmcmdn];
static realtype CONSTANTS[AC_kmcsqn];
static realtype CONSTANTS[AC_kmtrpn];
static realtype CONSTANTS[AC_trpnmax];
static realtype ALGEBRAIC[AV_Istim];
static realtype CONSTANTS[AC_amp];
static realtype CONSTANTS[AC_duration];
static realtype CONSTANTS[AC_stimStart];
static realtype ALGEBRAIC[AV_vffrt];
static realtype ALGEBRAIC[AV_vfrt];
static realtype CONSTANTS[AC_F];
static realtype CONSTANTS[AC_R];
static realtype CONSTANTS[AC_T];
static realtype CONSTANTS[AC_zca];
static realtype CONSTANTS[AC_zk];
static realtype CONSTANTS[AC_zna];
static realtype ALGEBRAIC[AV_EK];
static realtype ALGEBRAIC[AV_EKs];
static realtype ALGEBRAIC[AV_ENa];
static realtype CONSTANTS[AC_PKNa];
static realtype ALGEBRAIC[AV_Jrel];
static realtype ALGEBRAIC[AV_Jrel_inf];
static realtype ALGEBRAIC[AV_Jrel_inf_temp];
static realtype ALGEBRAIC[AV_Jrel_infp];
static realtype ALGEBRAIC[AV_Jrel_temp];
static realtype CONSTANTS[AC_a_rel];
static realtype CONSTANTS[AC_a_relp];
static realtype CONSTANTS[AC_bt];
static realtype CONSTANTS[AC_btp];
static realtype ALGEBRAIC[AV_fJrelp];
static realtype ALGEBRAIC[AV_tau_rel];
static realtype ALGEBRAIC[AV_tau_rel_temp];
static realtype ALGEBRAIC[AV_tau_relp];
static realtype ALGEBRAIC[AV_tau_relp_temp];
static realtype ALGEBRAIC[AV_Jtr];

/* Set values of constants */
static void
updateConstants(void)
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
    
}

/* Right-hand-side function of the model ODE */
static int rhs(realtype t, N_Vector y, N_Vector ydot, void *f_data)
{
    /* CaMK */
    ALGEBRAIC[AV_CaMKb] = CONSTANTS[AC_CaMKo] * (1.0 - STATES[1]) / (1.0 + CONSTANTS[AC_KmCaM] / STATES[6]);
    ALGEBRAIC[AV_CaMKa] = ALGEBRAIC[AV_CaMKb] + STATES[1];
    RATES[1] = CONSTANTS[AC_aCaMK] * ALGEBRAIC[AV_CaMKb] * (ALGEBRAIC[AV_CaMKb] + STATES[1]) - CONSTANTS[AC_bCaMK] * STATES[1];
    
    /* IpCa */
    ALGEBRAIC[AV_IpCa_IpCa] = CONSTANTS[AC_GpCa] * STATES[9] / (CONSTANTS[AC_KmCap] + STATES[9]);
    
    /* diff */
    ALGEBRAIC[AV_Jdiff] = (STATES[6] - STATES[9]) / 0.2;
    ALGEBRAIC[AV_JdiffK] = (STATES[5] - STATES[4]) / 2.0;
    ALGEBRAIC[AV_JdiffNa] = (STATES[3] - STATES[2]) / 2.0;
    
    /* environment */
    ALGEBRAIC[AV_time] = t;
    
    /* trans_flux */
    ALGEBRAIC[AV_Jtr] = (STATES[7] - STATES[8]) / 100.0;
    
    /* INaCa_i */
    ALGEBRAIC[AV_allo_i] = 1.0 / (1.0 + pow(CONSTANTS[AC_KmCaAct] / STATES[9], 2.0));
    ALGEBRAIC[AV_allo_ss] = 1.0 / (1.0 + pow(CONSTANTS[AC_KmCaAct] / STATES[6], 2.0));
    ALGEBRAIC[AV_h4_i] = 1.0 + STATES[2] / CONSTANTS[AC_kna1] * (1.0 + STATES[2] / CONSTANTS[AC_kna2]);
    ALGEBRAIC[AV_h4_ss] = 1.0 + STATES[3] / CONSTANTS[AC_kna1] * (1.0 + STATES[3] / CONSTANTS[AC_kna2]);
    ALGEBRAIC[AV_hca] = exp(CONSTANTS[AC_qca] * STATES[0] * CONSTANTS[AC_F] / (CONSTANTS[AC_R] * CONSTANTS[AC_T]));
    ALGEBRAIC[AV_hna] = exp(CONSTANTS[AC_qna] * STATES[0] * CONSTANTS[AC_F] / (CONSTANTS[AC_R] * CONSTANTS[AC_T]));
    ALGEBRAIC[AV_h1_i] = 1.0 + STATES[2] / CONSTANTS[AC_kna3] * (1.0 + ALGEBRAIC[AV_hna]);
    ALGEBRAIC[AV_h1_ss] = 1.0 + STATES[3] / CONSTANTS[AC_kna3] * (1.0 + ALGEBRAIC[AV_hna]);
    ALGEBRAIC[AV_h5_i] = STATES[2] * STATES[2] / (ALGEBRAIC[AV_h4_i] * CONSTANTS[AC_kna1] * CONSTANTS[AC_kna2]);
    ALGEBRAIC[AV_h5_ss] = STATES[3] * STATES[3] / (ALGEBRAIC[AV_h4_ss] * CONSTANTS[AC_kna1] * CONSTANTS[AC_kna2]);
    ALGEBRAIC[AV_h6_i] = 1.0 / ALGEBRAIC[AV_h4_i];
    ALGEBRAIC[AV_h6_ss] = 1.0 / ALGEBRAIC[AV_h4_ss];
    ALGEBRAIC[AV_h7_i] = 1.0 + CONSTANTS[AC_nao] / CONSTANTS[AC_kna3] * (1.0 + 1.0 / ALGEBRAIC[AV_hna]);
    ALGEBRAIC[AV_h7_ss] = 1.0 + CONSTANTS[AC_nao] / CONSTANTS[AC_kna3] * (1.0 + 1.0 / ALGEBRAIC[AV_hna]);
    ALGEBRAIC[AV_h2_i] = STATES[2] * ALGEBRAIC[AV_hna] / (CONSTANTS[AC_kna3] * ALGEBRAIC[AV_h1_i]);
    ALGEBRAIC[AV_h2_ss] = STATES[3] * ALGEBRAIC[AV_hna] / (CONSTANTS[AC_kna3] * ALGEBRAIC[AV_h1_ss]);
    ALGEBRAIC[AV_h3_i] = 1.0 / ALGEBRAIC[AV_h1_i];
    ALGEBRAIC[AV_h3_ss] = 1.0 / ALGEBRAIC[AV_h1_ss];
    ALGEBRAIC[AV_h8_i] = CONSTANTS[AC_nao] / (CONSTANTS[AC_kna3] * ALGEBRAIC[AV_hna] * ALGEBRAIC[AV_h7_i]);
    ALGEBRAIC[AV_h8_ss] = CONSTANTS[AC_nao] / (CONSTANTS[AC_kna3] * ALGEBRAIC[AV_hna] * ALGEBRAIC[AV_h7_ss]);
    ALGEBRAIC[AV_h9_i] = 1.0 / ALGEBRAIC[AV_h7_i];
    ALGEBRAIC[AV_h9_ss] = 1.0 / ALGEBRAIC[AV_h7_ss];
    ALGEBRAIC[AV_k6_i] = ALGEBRAIC[AV_h6_i] * STATES[9] * CONSTANTS[AC_kcaon];
    ALGEBRAIC[AV_k6_ss] = ALGEBRAIC[AV_h6_ss] * STATES[6] * CONSTANTS[AC_kcaon];
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
    ALGEBRAIC[AV_Knai] = CONSTANTS[AC_Knai0] * exp(CONSTANTS[AC_delta] * STATES[0] * CONSTANTS[AC_F] / (3.0 * CONSTANTS[AC_R] * CONSTANTS[AC_T]));
    ALGEBRAIC[AV_Knao] = CONSTANTS[AC_Knao0] * exp((1.0 - CONSTANTS[AC_delta]) * STATES[0] * CONSTANTS[AC_F] / (3.0 * CONSTANTS[AC_R] * CONSTANTS[AC_T]));
    ALGEBRAIC[AV_P] = CONSTANTS[AC_eP] / (1.0 + CONSTANTS[AC_H] / CONSTANTS[AC_Khp] + STATES[2] / CONSTANTS[AC_Knap] + STATES[4] / CONSTANTS[AC_Kxkur]);
    ALGEBRAIC[AV_a1] = CONSTANTS[AC_k1p] * pow(STATES[2] / ALGEBRAIC[AV_Knai], 3.0) / (pow(1.0 + STATES[2] / ALGEBRAIC[AV_Knai], 3.0) + pow(1.0 + STATES[4] / CONSTANTS[AC_Kki], 2.0) - 1.0);
    ALGEBRAIC[AV_a3] = CONSTANTS[AC_k3p] * pow(CONSTANTS[AC_ko] / CONSTANTS[AC_Kko], 2.0) / (pow(1.0 + CONSTANTS[AC_nao] / ALGEBRAIC[AV_Knao], 3.0) + pow(1.0 + CONSTANTS[AC_ko] / CONSTANTS[AC_Kko], 2.0) - 1.0);
    ALGEBRAIC[AV_b2] = CONSTANTS[AC_k2m] * pow(CONSTANTS[AC_nao] / ALGEBRAIC[AV_Knao], 3.0) / (pow(1.0 + CONSTANTS[AC_nao] / ALGEBRAIC[AV_Knao], 3.0) + pow(1.0 + CONSTANTS[AC_ko] / CONSTANTS[AC_Kko], 2.0) - 1.0);
    ALGEBRAIC[AV_b3] = CONSTANTS[AC_k3m] * ALGEBRAIC[AV_P] * CONSTANTS[AC_H] / (1.0 + CONSTANTS[AC_MgATP] / CONSTANTS[AC_Kmgatp]);
    ALGEBRAIC[AV_b4] = CONSTANTS[AC_k4m] * pow(STATES[4] / CONSTANTS[AC_Kki], 2.0) / (pow(1.0 + STATES[2] / ALGEBRAIC[AV_Knai], 3.0) + pow(1.0 + STATES[4] / CONSTANTS[AC_Kki], 2.0) - 1.0);
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
    ALGEBRAIC[AV_Jleak] = 0.0039375 * STATES[7] / 15.0;
    ALGEBRAIC[AV_fJupp] = 1.0 / (1.0 + CONSTANTS[AC_KmCaMK] / ALGEBRAIC[AV_CaMKa]);
    ALGEBRAIC[AV_Jupnp] = CONSTANTS[AC_upScale] * 0.004375 * STATES[9] / (STATES[9] + 0.00092);
    ALGEBRAIC[AV_Jupp] = CONSTANTS[AC_upScale] * 2.75 * 0.004375 * STATES[9] / (STATES[9] + 0.00092 - 0.00017);
    ALGEBRAIC[AV_Jup] = (1.0 - ALGEBRAIC[AV_fJupp]) * ALGEBRAIC[AV_Jupnp] + ALGEBRAIC[AV_fJupp] * ALGEBRAIC[AV_Jupp] - ALGEBRAIC[AV_Jleak];
    
    /* reversal_potentials */
    ALGEBRAIC[AV_EK] = CONSTANTS[AC_R] * CONSTANTS[AC_T] / CONSTANTS[AC_F] * log(CONSTANTS[AC_ko] / STATES[4]);
    ALGEBRAIC[AV_ENa] = CONSTANTS[AC_R] * CONSTANTS[AC_T] / CONSTANTS[AC_F] * log(CONSTANTS[AC_nao] / STATES[2]);
    ALGEBRAIC[AV_EKs] = CONSTANTS[AC_R] * CONSTANTS[AC_T] / CONSTANTS[AC_F] * log((CONSTANTS[AC_ko] + CONSTANTS[AC_PKNa] * CONSTANTS[AC_nao]) / (STATES[4] + CONSTANTS[AC_PKNa] * STATES[2]));
    
    /* IK1 */
    ALGEBRAIC[AV_rk1] = 1.0 / (1.0 + exp((STATES[0] + 105.8 - 2.6 * CONSTANTS[AC_ko]) / 9.493));
    ALGEBRAIC[AV_txk1] = 122.2 / (exp(-(STATES[0] + 127.2) / 20.36) + exp((STATES[0] + 236.8) / 69.33));
    ALGEBRAIC[AV_xk1ss] = 1.0 / (1.0 + exp(-(STATES[0] + 2.5538 * CONSTANTS[AC_ko] + 144.59) / (1.5692 * CONSTANTS[AC_ko] + 3.8115)));
    RATES[38] = (ALGEBRAIC[AV_xk1ss] - STATES[38]) / ALGEBRAIC[AV_txk1];
    ALGEBRAIC[AV_IK1_IK1] = CONSTANTS[AC_GK1] * sqrt(CONSTANTS[AC_ko]) * ALGEBRAIC[AV_rk1] * STATES[38] * (STATES[0] - ALGEBRAIC[AV_EK]);
    
    /* IKb */
    ALGEBRAIC[AV_xkb] = 1.0 / (1.0 + exp(-(STATES[0] - 14.48) / 18.34));
    ALGEBRAIC[AV_IKb_IKb] = CONSTANTS[AC_GKb] * ALGEBRAIC[AV_xkb] * (STATES[0] - ALGEBRAIC[AV_EK]);
    
    /* IKr */
    ALGEBRAIC[AV_Axrf] = 1.0 / (1.0 + exp((STATES[0] + 54.81) / 38.21));
    ALGEBRAIC[AV_rkr] = 1.0 / (1.0 + exp((STATES[0] + 55.0) / 75.0)) * 1.0 / (1.0 + exp((STATES[0] - 10.0) / 30.0));
    ALGEBRAIC[AV_txrf] = 12.98 + 1.0 / (0.3652 * exp((STATES[0] - 31.66) / 3.869) + 4.123e-05 * exp(-(STATES[0] - 47.78) / 20.38));
    ALGEBRAIC[AV_txrs] = 1.865 + 1.0 / (0.06629 * exp((STATES[0] - 34.7) / 7.355) + 1.128e-05 * exp(-(STATES[0] - 29.74) / 25.94));
    ALGEBRAIC[AV_xrss] = 1.0 / (1.0 + exp(-(STATES[0] + 8.337) / 6.789));
    ALGEBRAIC[AV_Axrs] = 1.0 - ALGEBRAIC[AV_Axrf];
    RATES[34] = (ALGEBRAIC[AV_xrss] - STATES[34]) / ALGEBRAIC[AV_txrf];
    RATES[35] = (ALGEBRAIC[AV_xrss] - STATES[35]) / ALGEBRAIC[AV_txrs];
    ALGEBRAIC[AV_xr] = ALGEBRAIC[AV_Axrf] * STATES[34] + ALGEBRAIC[AV_Axrs] * STATES[35];
    ALGEBRAIC[AV_IKr_IKr] = CONSTANTS[AC_GKr] * sqrt(CONSTANTS[AC_ko] / 5.4) * ALGEBRAIC[AV_xr] * ALGEBRAIC[AV_rkr] * (STATES[0] - ALGEBRAIC[AV_EK]);
    
    /* IKs */
    ALGEBRAIC[AV_KsCa] = 1.0 + 0.6 / (1.0 + pow(3.8e-05 / STATES[9], 1.4));
    ALGEBRAIC[AV_txs1] = 817.3 + 1.0 / (0.0002326 * exp((STATES[0] + 48.28) / 17.8) + 0.001292 * exp(-(STATES[0] + 210.0) / 230.0));
    ALGEBRAIC[AV_txs2] = 1.0 / (0.01 * exp((STATES[0] - 50.0) / 20.0) + 0.0193 * exp(-(STATES[0] + 66.54) / 31.0));
    ALGEBRAIC[AV_xs1ss] = 1.0 / (1.0 + exp(-(STATES[0] + 11.6) / 8.932));
    ALGEBRAIC[AV_xs2ss] = ALGEBRAIC[AV_xs1ss];
    RATES[36] = (ALGEBRAIC[AV_xs1ss] - STATES[36]) / ALGEBRAIC[AV_txs1];
    ALGEBRAIC[AV_IKs_IKs] = CONSTANTS[AC_GKs] * ALGEBRAIC[AV_KsCa] * STATES[36] * STATES[37] * (STATES[0] - ALGEBRAIC[AV_EKs]);
    RATES[37] = (ALGEBRAIC[AV_xs2ss] - STATES[37]) / ALGEBRAIC[AV_txs2];
    
    /* INa */
    ALGEBRAIC[AV_fINap] = 1.0 / (1.0 + CONSTANTS[AC_KmCaMK] / ALGEBRAIC[AV_CaMKa]);
    ALGEBRAIC[AV_hssp] = 1.0 / (1.0 + exp((STATES[0] + 89.1) / 6.086));
    ALGEBRAIC[AV_thf] = 1.0 / (1.432e-05 * exp(-(STATES[0] + 1.196) / 6.285) + 6.149 * exp((STATES[0] + 0.5096) / 20.27));
    ALGEBRAIC[AV_ths] = 1.0 / (0.009794 * exp(-(STATES[0] + 17.95) / 28.05) + 0.3343 * exp((STATES[0] + 5.73) / 56.66));
    ALGEBRAIC[AV_tj] = 2.038 + 1.0 / (0.02136 * exp(-(STATES[0] + 100.6) / 8.281) + 0.3052 * exp((STATES[0] + 0.9941) / 38.45));
    ALGEBRAIC[AV_hss] = 1.0 / (1.0 + exp((STATES[0] + CONSTANTS[AC_hssV1]) / CONSTANTS[AC_hssV2]));
    ALGEBRAIC[AV_mss] = 1.0 / (1.0 + exp(-(STATES[0] + CONSTANTS[AC_mssV1]) / CONSTANTS[AC_mssV2]));
    ALGEBRAIC[AV_thsp] = 3.0 * ALGEBRAIC[AV_ths];
    ALGEBRAIC[AV_tjp] = 1.46 * ALGEBRAIC[AV_tj];
    ALGEBRAIC[AV_tm] = 1.0 / (CONSTANTS[AC_mtD1] * exp((STATES[0] + CONSTANTS[AC_mtV1]) / CONSTANTS[AC_mtV2]) + CONSTANTS[AC_mtD2] * exp(-(STATES[0] + CONSTANTS[AC_mtV3]) / CONSTANTS[AC_mtV4]));
    ALGEBRAIC[AV_h] = CONSTANTS[AC_Ahf] * STATES[11] + CONSTANTS[AC_Ahs] * STATES[12];
    ALGEBRAIC[AV_hp] = CONSTANTS[AC_Ahf] * STATES[11] + CONSTANTS[AC_Ahs] * STATES[14];
    ALGEBRAIC[AV_jss] = ALGEBRAIC[AV_hss];
    RATES[11] = (ALGEBRAIC[AV_hss] - STATES[11]) / ALGEBRAIC[AV_thf];
    RATES[12] = (ALGEBRAIC[AV_hss] - STATES[12]) / ALGEBRAIC[AV_ths];
    RATES[14] = (ALGEBRAIC[AV_hssp] - STATES[14]) / ALGEBRAIC[AV_thsp];
    RATES[10] = (ALGEBRAIC[AV_mss] - STATES[10]) / ALGEBRAIC[AV_tm];
    ALGEBRAIC[AV_INa_INa] = CONSTANTS[AC_GNa] * (STATES[0] - ALGEBRAIC[AV_ENa]) * pow(STATES[10], 3.0) * ((1.0 - ALGEBRAIC[AV_fINap]) * ALGEBRAIC[AV_h] * STATES[13] + ALGEBRAIC[AV_fINap] * ALGEBRAIC[AV_hp] * STATES[15]);
    RATES[13] = (ALGEBRAIC[AV_jss] - STATES[13]) / ALGEBRAIC[AV_tj];
    RATES[15] = (ALGEBRAIC[AV_jss] - STATES[15]) / ALGEBRAIC[AV_tjp];
    
    /* Ito */
    ALGEBRAIC[AV_AiF] = 1.0 / (1.0 + exp((STATES[0] - 213.6) / 151.2));
    ALGEBRAIC[AV_ass] = 1.0 / (1.0 + exp(-(STATES[0] - 14.34) / 14.82));
    ALGEBRAIC[AV_assp] = 1.0 / (1.0 + exp(-(STATES[0] - 24.34) / 14.82));
    ALGEBRAIC[AV_delta_epi] = ((CONSTANTS[AC_celltype] == 1.0) ? 1.0 - 0.95 / (1.0 + exp((STATES[0] + 70.0) / 5.0)) : 1.0);
    ALGEBRAIC[AV_dti_develop] = 1.354 + 0.0001 / (exp((STATES[0] - 167.4) / 15.89) + exp(-(STATES[0] - 12.23) / 0.2154));
    ALGEBRAIC[AV_dti_recover] = 1.0 - 0.5 / (1.0 + exp((STATES[0] + 70.0) / 20.0));
    ALGEBRAIC[AV_fItop] = 1.0 / (1.0 + CONSTANTS[AC_KmCaMK] / ALGEBRAIC[AV_CaMKa]);
    ALGEBRAIC[AV_iss] = 1.0 / (1.0 + exp((STATES[0] + 43.94) / 5.711));
    ALGEBRAIC[AV_ta] = 1.0515 / (1.0 / (1.2089 * (1.0 + exp(-(STATES[0] - 18.4099) / 29.3814))) + 3.5 / (1.0 + exp((STATES[0] + 100.0) / 29.3814)));
    ALGEBRAIC[AV_tiF_b] = 4.562 + 1.0 / (0.3933 * exp(-(STATES[0] + 100.0) / 100.0) + 0.08004 * exp((STATES[0] + 50.0) / 16.59));
    ALGEBRAIC[AV_tiS_b] = 23.62 + 1.0 / (0.001416 * exp(-(STATES[0] + 96.52) / 59.05) + 1.78e-08 * exp((STATES[0] + 114.1) / 8.079));
    ALGEBRAIC[AV_AiS] = 1.0 - ALGEBRAIC[AV_AiF];
    ALGEBRAIC[AV_tiF] = ALGEBRAIC[AV_tiF_b] * ALGEBRAIC[AV_delta_epi];
    ALGEBRAIC[AV_tiS] = ALGEBRAIC[AV_tiS_b] * ALGEBRAIC[AV_delta_epi];
    RATES[19] = (ALGEBRAIC[AV_ass] - STATES[19]) / ALGEBRAIC[AV_ta];
    RATES[22] = (ALGEBRAIC[AV_assp] - STATES[22]) / ALGEBRAIC[AV_ta];
    ALGEBRAIC[AV_i] = ALGEBRAIC[AV_AiF] * STATES[20] + ALGEBRAIC[AV_AiS] * STATES[21];
    ALGEBRAIC[AV_ip] = ALGEBRAIC[AV_AiF] * STATES[23] + ALGEBRAIC[AV_AiS] * STATES[24];
    ALGEBRAIC[AV_tiFp] = ALGEBRAIC[AV_dti_develop] * ALGEBRAIC[AV_dti_recover] * ALGEBRAIC[AV_tiF];
    ALGEBRAIC[AV_tiSp] = ALGEBRAIC[AV_dti_develop] * ALGEBRAIC[AV_dti_recover] * ALGEBRAIC[AV_tiS];
    RATES[20] = (ALGEBRAIC[AV_iss] - STATES[20]) / ALGEBRAIC[AV_tiF];
    RATES[21] = (ALGEBRAIC[AV_iss] - STATES[21]) / ALGEBRAIC[AV_tiS];
    ALGEBRAIC[AV_Ito_Ito] = CONSTANTS[AC_Gto] * (STATES[0] - ALGEBRAIC[AV_EK]) * ((1.0 - ALGEBRAIC[AV_fItop]) * STATES[19] * ALGEBRAIC[AV_i] + ALGEBRAIC[AV_fItop] * STATES[22] * ALGEBRAIC[AV_ip]);
    RATES[23] = (ALGEBRAIC[AV_iss] - STATES[23]) / ALGEBRAIC[AV_tiFp];
    RATES[24] = (ALGEBRAIC[AV_iss] - STATES[24]) / ALGEBRAIC[AV_tiSp];
    
    /* INaL */
    ALGEBRAIC[AV_fINaLp] = 1.0 / (1.0 + CONSTANTS[AC_KmCaMK] / ALGEBRAIC[AV_CaMKa]);
    ALGEBRAIC[AV_hLss] = 1.0 / (1.0 + exp((STATES[0] + 87.61) / 7.488));
    ALGEBRAIC[AV_hLssp] = 1.0 / (1.0 + exp((STATES[0] + 93.81) / 7.488));
    ALGEBRAIC[AV_mLss] = 1.0 / (1.0 + exp(-(STATES[0] + 42.85) / 5.264));
    ALGEBRAIC[AV_tmL] = ALGEBRAIC[AV_tm];
    RATES[17] = (ALGEBRAIC[AV_hLss] - STATES[17]) / CONSTANTS[AC_thL];
    RATES[16] = (ALGEBRAIC[AV_mLss] - STATES[16]) / ALGEBRAIC[AV_tmL];
    ALGEBRAIC[AV_INaL_INaL] = CONSTANTS[AC_GNaL] * (STATES[0] - ALGEBRAIC[AV_ENa]) * STATES[16] * ((1.0 - ALGEBRAIC[AV_fINaLp]) * STATES[17] + ALGEBRAIC[AV_fINaLp] * STATES[18]);
    RATES[18] = (ALGEBRAIC[AV_hLssp] - STATES[18]) / CONSTANTS[AC_thLp];
    
    /* ICaL */
    ALGEBRAIC[AV_Afcaf] = 0.3 + 0.6 / (1.0 + exp((STATES[0] - 10.0) / 10.0));
    ALGEBRAIC[AV_dss] = 1.0 / (1.0 + exp(-(STATES[0] + 3.94) / 4.23));
    ALGEBRAIC[AV_fICaLp] = 1.0 / (1.0 + CONSTANTS[AC_KmCaMK] / ALGEBRAIC[AV_CaMKa]);
    ALGEBRAIC[AV_fss] = 1.0 / (1.0 + exp((STATES[0] + 19.58) / 3.696));
    ALGEBRAIC[AV_km2n] = STATES[30] * 1.0;
    ALGEBRAIC[AV_td] = 0.6 + 1.0 / (exp(-0.05 * (STATES[0] + 6.0)) + exp(0.09 * (STATES[0] + 14.0)));
    ALGEBRAIC[AV_tfcaf] = 7.0 + 1.0 / (0.04 * exp(-(STATES[0] - 4.0) / 7.0) + 0.04 * exp((STATES[0] - 4.0) / 7.0));
    ALGEBRAIC[AV_tfcas] = 100.0 + 1.0 / (0.00012 * exp(-STATES[0] / 3.0) + 0.00012 * exp(STATES[0] / 7.0));
    ALGEBRAIC[AV_tff] = 7.0 + 1.0 / (0.0045 * exp(-(STATES[0] + 20.0) / 10.0) + 0.0045 * exp((STATES[0] + 20.0) / 10.0));
    ALGEBRAIC[AV_tfs] = 1000.0 + 1.0 / (3.5e-05 * exp(-(STATES[0] + 5.0) / 4.0) + 3.5e-05 * exp((STATES[0] + 5.0) / 6.0));
    ALGEBRAIC[AV_Afcas] = 1.0 - ALGEBRAIC[AV_Afcaf];
    ALGEBRAIC[AV_anca] = 1.0 / (CONSTANTS[AC_k2n] / ALGEBRAIC[AV_km2n] + pow(1.0 + CONSTANTS[AC_Kmn] / STATES[6], 4.0));
    ALGEBRAIC[AV_fcass] = ALGEBRAIC[AV_fss];
    ALGEBRAIC[AV_tfcafp] = 2.5 * ALGEBRAIC[AV_tfcaf];
    ALGEBRAIC[AV_tffp] = 2.5 * ALGEBRAIC[AV_tff];
    RATES[25] = (ALGEBRAIC[AV_dss] - STATES[25]) / ALGEBRAIC[AV_td];
    RATES[26] = (ALGEBRAIC[AV_fss] - STATES[26]) / ALGEBRAIC[AV_tff];
    RATES[27] = (ALGEBRAIC[AV_fss] - STATES[27]) / ALGEBRAIC[AV_tfs];
    ALGEBRAIC[AV_f] = CONSTANTS[AC_Aff] * STATES[26] + CONSTANTS[AC_Afs] * STATES[27];
    ALGEBRAIC[AV_fca] = ALGEBRAIC[AV_Afcaf] * STATES[28] + ALGEBRAIC[AV_Afcas] * STATES[29];
    ALGEBRAIC[AV_fcap] = ALGEBRAIC[AV_Afcaf] * STATES[32] + ALGEBRAIC[AV_Afcas] * STATES[29];
    ALGEBRAIC[AV_fp] = CONSTANTS[AC_Aff] * STATES[31] + CONSTANTS[AC_Afs] * STATES[27];
    RATES[28] = (ALGEBRAIC[AV_fcass] - STATES[28]) / ALGEBRAIC[AV_tfcaf];
    RATES[32] = (ALGEBRAIC[AV_fcass] - STATES[32]) / ALGEBRAIC[AV_tfcafp];
    RATES[29] = (ALGEBRAIC[AV_fcass] - STATES[29]) / ALGEBRAIC[AV_tfcas];
    RATES[31] = (ALGEBRAIC[AV_fss] - STATES[31]) / ALGEBRAIC[AV_tffp];
    RATES[30] = (ALGEBRAIC[AV_fcass] - STATES[30]) / CONSTANTS[AC_tjca];
    RATES[33] = ALGEBRAIC[AV_anca] * CONSTANTS[AC_k2n] - STATES[33] * ALGEBRAIC[AV_km2n];
    
    /* intracellular_ions */
    RATES[7] = ALGEBRAIC[AV_Jup] - ALGEBRAIC[AV_Jtr] * CONSTANTS[AC_vjsr] / CONSTANTS[AC_vnsr];
    ALGEBRAIC[AV_Bcajsr] = 1.0 / (1.0 + CONSTANTS[AC_csqnmax] * CONSTANTS[AC_kmcsqn] / pow(CONSTANTS[AC_kmcsqn] + STATES[8], 2.0));
    ALGEBRAIC[AV_Bcass] = 1.0 / (1.0 + CONSTANTS[AC_BSRmax] * CONSTANTS[AC_KmBSR] / pow(CONSTANTS[AC_KmBSR] + STATES[6], 2.0) + CONSTANTS[AC_BSLmax] * CONSTANTS[AC_KmBSL] / pow(CONSTANTS[AC_KmBSL] + STATES[6], 2.0));
    ALGEBRAIC[AV_Bcai] = 1.0 / (1.0 + CONSTANTS[AC_cmdnmax] * CONSTANTS[AC_kmcmdn] / pow(CONSTANTS[AC_kmcmdn] + STATES[9], 2.0) + CONSTANTS[AC_trpnmax] * CONSTANTS[AC_kmtrpn] / pow(CONSTANTS[AC_kmtrpn] + STATES[9], 2.0));
    
    /* membrane */
    ALGEBRAIC[AV_vffrt] = STATES[0] * CONSTANTS[AC_F] * CONSTANTS[AC_F] / (CONSTANTS[AC_R] * CONSTANTS[AC_T]);
    ALGEBRAIC[AV_vfrt] = STATES[0] * CONSTANTS[AC_F] / (CONSTANTS[AC_R] * CONSTANTS[AC_T]);
    ALGEBRAIC[AV_Istim] = (((ALGEBRAIC[AV_time] > CONSTANTS[AC_stimStart]) && (ALGEBRAIC[AV_time] <= CONSTANTS[AC_stimStart] + CONSTANTS[AC_duration])) ? CONSTANTS[AC_amp] : 0.0);
    
    /* ryr */
    ALGEBRAIC[AV_fJrelp] = 1.0 / (1.0 + CONSTANTS[AC_KmCaMK] / ALGEBRAIC[AV_CaMKa]);
    ALGEBRAIC[AV_Jrel] = (1.0 - ALGEBRAIC[AV_fJrelp]) * STATES[39] + ALGEBRAIC[AV_fJrelp] * STATES[40];
    ALGEBRAIC[AV_tau_rel_temp] = CONSTANTS[AC_bt] / (1.0 + 0.0123 / STATES[8]);
    ALGEBRAIC[AV_tau_rel] = ((ALGEBRAIC[AV_tau_rel_temp] < 0.001) ? 0.001 : ALGEBRAIC[AV_tau_rel_temp]);
    ALGEBRAIC[AV_tau_relp_temp] = CONSTANTS[AC_btp] / (1.0 + 0.0123 / STATES[8]);
    ALGEBRAIC[AV_tau_relp] = ((ALGEBRAIC[AV_tau_relp_temp] < 0.001) ? 0.001 : ALGEBRAIC[AV_tau_relp_temp]);
    
    /* *remaining* */
    ALGEBRAIC[AV_PhiCaK] = 1.0 * ALGEBRAIC[AV_vffrt] * (0.75 * STATES[5] * exp(1.0 * ALGEBRAIC[AV_vfrt]) - 0.75 * CONSTANTS[AC_ko]) / (exp(1.0 * ALGEBRAIC[AV_vfrt]) - 1.0);
    ALGEBRAIC[AV_PhiCaL] = 4.0 * ALGEBRAIC[AV_vffrt] * (STATES[6] * exp(2.0 * ALGEBRAIC[AV_vfrt]) - 0.341 * CONSTANTS[AC_cao]) / (exp(2.0 * ALGEBRAIC[AV_vfrt]) - 1.0);
    ALGEBRAIC[AV_PhiCaNa] = 1.0 * ALGEBRAIC[AV_vffrt] * (0.75 * STATES[3] * exp(1.0 * ALGEBRAIC[AV_vfrt]) - 0.75 * CONSTANTS[AC_nao]) / (exp(1.0 * ALGEBRAIC[AV_vfrt]) - 1.0);
    ALGEBRAIC[AV_ICab_ICab] = CONSTANTS[AC_PCab] * 4.0 * ALGEBRAIC[AV_vffrt] * (STATES[9] * exp(2.0 * ALGEBRAIC[AV_vfrt]) - 0.341 * CONSTANTS[AC_cao]) / (exp(2.0 * ALGEBRAIC[AV_vfrt]) - 1.0);
    ALGEBRAIC[AV_INab_INab] = CONSTANTS[AC_PNab] * ALGEBRAIC[AV_vffrt] * (STATES[2] * exp(ALGEBRAIC[AV_vfrt]) - CONSTANTS[AC_nao]) / (exp(ALGEBRAIC[AV_vfrt]) - 1.0);
    RATES[8] = ALGEBRAIC[AV_Bcajsr] * (ALGEBRAIC[AV_Jtr] - ALGEBRAIC[AV_Jrel]);
    RATES[4] = -(ALGEBRAIC[AV_Ito_Ito] + ALGEBRAIC[AV_IKr_IKr] + ALGEBRAIC[AV_IKs_IKs] + ALGEBRAIC[AV_IK1_IK1] + ALGEBRAIC[AV_IKb_IKb] + ALGEBRAIC[AV_Istim] - 2.0 * ALGEBRAIC[AV_INaK_INaK]) * CONSTANTS[AC_cm] * CONSTANTS[AC_Acap] / (CONSTANTS[AC_F] * CONSTANTS[AC_vmyo]) + ALGEBRAIC[AV_JdiffK] * CONSTANTS[AC_vss] / CONSTANTS[AC_vmyo];
    ALGEBRAIC[AV_ICaK] = (1.0 - ALGEBRAIC[AV_fICaLp]) * CONSTANTS[AC_PCaK] * ALGEBRAIC[AV_PhiCaK] * STATES[25] * (ALGEBRAIC[AV_f] * (1.0 - STATES[33]) + STATES[30] * ALGEBRAIC[AV_fca] * STATES[33]) + ALGEBRAIC[AV_fICaLp] * CONSTANTS[AC_PCaKp] * ALGEBRAIC[AV_PhiCaK] * STATES[25] * (ALGEBRAIC[AV_fp] * (1.0 - STATES[33]) + STATES[30] * ALGEBRAIC[AV_fcap] * STATES[33]);
    ALGEBRAIC[AV_ICaL_ICaL] = (1.0 - ALGEBRAIC[AV_fICaLp]) * CONSTANTS[AC_PCa] * ALGEBRAIC[AV_PhiCaL] * STATES[25] * (ALGEBRAIC[AV_f] * (1.0 - STATES[33]) + STATES[30] * ALGEBRAIC[AV_fca] * STATES[33]) + ALGEBRAIC[AV_fICaLp] * CONSTANTS[AC_PCap] * ALGEBRAIC[AV_PhiCaL] * STATES[25] * (ALGEBRAIC[AV_fp] * (1.0 - STATES[33]) + STATES[30] * ALGEBRAIC[AV_fcap] * STATES[33]);
    ALGEBRAIC[AV_ICaNa] = (1.0 - ALGEBRAIC[AV_fICaLp]) * CONSTANTS[AC_PCaNa] * ALGEBRAIC[AV_PhiCaNa] * STATES[25] * (ALGEBRAIC[AV_f] * (1.0 - STATES[33]) + STATES[30] * ALGEBRAIC[AV_fca] * STATES[33]) + ALGEBRAIC[AV_fICaLp] * CONSTANTS[AC_PCaNap] * ALGEBRAIC[AV_PhiCaNa] * STATES[25] * (ALGEBRAIC[AV_fp] * (1.0 - STATES[33]) + STATES[30] * ALGEBRAIC[AV_fcap] * STATES[33]);
    RATES[9] = ALGEBRAIC[AV_Bcai] * (-(ALGEBRAIC[AV_IpCa_IpCa] + ALGEBRAIC[AV_ICab_ICab] - 2.0 * ALGEBRAIC[AV_INaCa_i_INaCa_i]) * CONSTANTS[AC_cm] * CONSTANTS[AC_Acap] / (2.0 * CONSTANTS[AC_F] * CONSTANTS[AC_vmyo]) - ALGEBRAIC[AV_Jup] * CONSTANTS[AC_vnsr] / CONSTANTS[AC_vmyo] + ALGEBRAIC[AV_Jdiff] * CONSTANTS[AC_vss] / CONSTANTS[AC_vmyo]);
    RATES[2] = -(ALGEBRAIC[AV_INa_INa] + ALGEBRAIC[AV_INaL_INaL] + 3.0 * ALGEBRAIC[AV_INaCa_i_INaCa_i] + 3.0 * ALGEBRAIC[AV_INaK_INaK] + ALGEBRAIC[AV_INab_INab]) * CONSTANTS[AC_Acap] * CONSTANTS[AC_cm] / (CONSTANTS[AC_F] * CONSTANTS[AC_vmyo]) + ALGEBRAIC[AV_JdiffNa] * CONSTANTS[AC_vss] / CONSTANTS[AC_vmyo];
    RATES[6] = ALGEBRAIC[AV_Bcass] * (-(ALGEBRAIC[AV_ICaL_ICaL] - 2.0 * ALGEBRAIC[AV_INaCa_ss]) * CONSTANTS[AC_cm] * CONSTANTS[AC_Acap] / (2.0 * CONSTANTS[AC_F] * CONSTANTS[AC_vss]) + ALGEBRAIC[AV_Jrel] * CONSTANTS[AC_vjsr] / CONSTANTS[AC_vss] - ALGEBRAIC[AV_Jdiff]);
    RATES[5] = -ALGEBRAIC[AV_ICaK] * CONSTANTS[AC_cm] * CONSTANTS[AC_Acap] / (CONSTANTS[AC_F] * CONSTANTS[AC_vss]) - ALGEBRAIC[AV_JdiffK];
    RATES[3] = -(ALGEBRAIC[AV_ICaNa] + 3.0 * ALGEBRAIC[AV_INaCa_ss]) * CONSTANTS[AC_cm] * CONSTANTS[AC_Acap] / (CONSTANTS[AC_F] * CONSTANTS[AC_vss]) - ALGEBRAIC[AV_JdiffNa];
    RATES[0] = -(ALGEBRAIC[AV_INa_INa] + ALGEBRAIC[AV_INaL_INaL] + ALGEBRAIC[AV_Ito_Ito] + ALGEBRAIC[AV_ICaL_ICaL] + ALGEBRAIC[AV_ICaNa] + ALGEBRAIC[AV_ICaK] + ALGEBRAIC[AV_IKr_IKr] + ALGEBRAIC[AV_IKs_IKs] + ALGEBRAIC[AV_IK1_IK1] + ALGEBRAIC[AV_INaCa_i_INaCa_i] + ALGEBRAIC[AV_INaCa_ss] + ALGEBRAIC[AV_INaK_INaK] + ALGEBRAIC[AV_INab_INab] + ALGEBRAIC[AV_IKb_IKb] + ALGEBRAIC[AV_IpCa_IpCa] + ALGEBRAIC[AV_ICab_ICab] + ALGEBRAIC[AV_Istim]);
    ALGEBRAIC[AV_Jrel_inf_temp] = CONSTANTS[AC_a_rel] * -ALGEBRAIC[AV_ICaL_ICaL] / (1.0 + 1.0 * pow(1.5 / STATES[8], 8.0));
    ALGEBRAIC[AV_Jrel_temp] = CONSTANTS[AC_a_relp] * -ALGEBRAIC[AV_ICaL_ICaL] / (1.0 + pow(1.5 / STATES[8], 8.0));
    ALGEBRAIC[AV_Jrel_inf] = ((CONSTANTS[AC_celltype] == 2.0) ? ALGEBRAIC[AV_Jrel_inf_temp] * 1.7 : ALGEBRAIC[AV_Jrel_inf_temp]);
    ALGEBRAIC[AV_Jrel_infp] = ((CONSTANTS[AC_celltype] == 2.0) ? ALGEBRAIC[AV_Jrel_temp] * 1.7 : ALGEBRAIC[AV_Jrel_temp]);
    RATES[39] = (ALGEBRAIC[AV_Jrel_inf] - STATES[39]) / ALGEBRAIC[AV_tau_rel];
    RATES[40] = (ALGEBRAIC[AV_Jrel_infp] - STATES[40]) / ALGEBRAIC[AV_tau_relp];
    

    return 0;
}

/* Set initial values */
static void
default_initial_values(N_Vector y)
{
    STATES[0] = -87.0;
    STATES[1] = 0.0;
    STATES[2] = 7.0;
    STATES[3] = 7.0;
    STATES[4] = 145.0;
    STATES[5] = 145.0;
    STATES[6] = 0.0001;
    STATES[7] = 1.2;
    STATES[8] = 1.2;
    STATES[9] = 0.0001;
    STATES[10] = 0.0;
    STATES[11] = 1.0;
    STATES[12] = 1.0;
    STATES[13] = 1.0;
    STATES[14] = 1.0;
    STATES[15] = 1.0;
    STATES[16] = 0.0;
    STATES[17] = 1.0;
    STATES[18] = 1.0;
    STATES[19] = 0.0;
    STATES[20] = 1.0;
    STATES[21] = 1.0;
    STATES[22] = 0.0;
    STATES[23] = 1.0;
    STATES[24] = 1.0;
    STATES[25] = 0.0;
    STATES[26] = 1.0;
    STATES[27] = 1.0;
    STATES[28] = 1.0;
    STATES[29] = 1.0;
    STATES[30] = 1.0;
    STATES[31] = 1.0;
    STATES[32] = 1.0;
    STATES[33] = 0.0;
    STATES[34] = 0.0;
    STATES[35] = 0.0;
    STATES[36] = 0.0;
    STATES[37] = 0.0;
    STATES[38] = 1.0;
    STATES[39] = 0.0;
    STATES[40] = 0.0;

}

/* Pacing event (non-zero stimulus) */
struct PacingEventS {
    double level;       /* The stimulus level (dimensionless, normal range [0,1]) */
    double start;       /* The time this stimulus starts */
    double duration;    /* The stimulus duration */
    double period;      /* The period with which it repeats (or 0 if it doesn't) */
    double multiplier;  /* The number of times this period occurs (or 0 if it doesn't) */
    struct PacingEventS* next;
};
typedef struct PacingEventS PacingEvent;

/*
 * Schedules a pacing event.
 * @param top The first event in a stack (the stack's head)
 * @param add The event to schedule
 * @return The new pointer to the head of the stack
 */
static PacingEvent*
PacingEvent_Schedule(PacingEvent* top, PacingEvent* add)
{
    add->next = 0;
    if (add == 0) return top;
    if (top == 0) return add;
    if (add->start <= top->start) {
        add->next = top;
        return add;
    }
    PacingEvent* evt = top;
    while(evt->next != 0 && evt->next->start <= add->start) {
        evt = evt->next;
    }
    add->next = evt->next;
    evt->next = add;
    return top;
}

/* CVODE Flags */
static int check_flag(void *flagvalue, char *funcname, int opt)
{
    int *errflag;
    /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
    if (opt == 0 && flagvalue == NULL) {
        fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n", funcname);
        return(1);
    } /* Check if flag < 0 */
    else if (opt == 1) {
        errflag = (int *) flagvalue;
        if (*errflag < 0) {
            fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n", funcname, *errflag);
            return(1);
        }
    } /* Check if function returned NULL pointer - no memory allocated */
    else if (opt == 2 && flagvalue == NULL) {
        fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n", funcname);
        return(1);
    }
    return 0;
}

/* Show output */
static void PrintOutput(realtype t, realtype y)
{
    #if defined(SUNDIALS_EXTENDED_PRECISION)
        printf("%4.1f     %14.6Le\n", t, y);
    #elif defined(SUNDIALS_DOUBLE_PRECISION)
        printf("%4.1f     %14.6le\n", t, y);
    #else
        printf("%4.1f     %14.6e\n", t, y);
    #endif
    return;
}

/* Run a simulation */
int main()
{
    int success = -1;

    /* Sundials error flag */
    int flag;

    /* Declare variables that will need freeing */
    PacingEvent* events = NULL;
    N_Vector y = NULL;
    N_Vector dy = NULL;
    void *cvode_mem = NULL;

    #if SUNDIALS_VERSION_MAJOR >= 3
    SUNMatrix sundials_dense_matrix;
    SUNLinearSolver sundials_linear_solver;
    #endif

    #if SUNDIALS_VERSION_MAJOR >= 6
    /* Create sundials context */
    SUNContext sundials_context;
    flag = SUNContext_Create(NULL, &sundials_context);
    if (check_flag(&flag, "SUNContext_Create", 1)) goto error;

    /* Create state vectors */
    y = N_VNew_Serial(N_STATE, sundials_context);
    dy = N_VNew_Serial(N_STATE, sundials_context);
    #else
    /* Create state vectors */
    y = N_VNew_Serial(N_STATE);
    dy = N_VNew_Serial(N_STATE);
    #endif
    if (check_flag((void*)y, "N_VNew_Serial", 0)) goto error;
    if (check_flag((void*)dy, "N_VNew_Serial", 0)) goto error;

    /* Set calculated constants */
    updateConstants();

    /* Set initial values */
    default_initial_values(y);

    /* Set integration times */
    double tMin = 0;
    double tMax = 1000;
    double tLog = 0;

    /* Create pacing events */

    int nPacing = 1;
    events = (PacingEvent*)malloc(sizeof(PacingEvent)*nPacing);
    if (events == 0) goto error;
    int iPacing = 0;
    events[iPacing].level = 1.0;
    events[iPacing].start = 100.0;
    events[iPacing].duration = 0.5;
    events[iPacing].period = 1000.0;
    events[iPacing].multiplier = 0;
    iPacing++;

    /* Schedule events, make "next" point to the first event */
    PacingEvent* next = events;
    PacingEvent* fire = events + 1;
    for(iPacing=1; iPacing<nPacing; iPacing++) {
        next = PacingEvent_Schedule(next, fire++);
    }

    /* Fast forward events to starting time */
    double tNext = next->start;
    double tDown = 0.0;
    fire = 0;
    while (tNext <= tMin) {
        /* Event over? */
        if (fire != 0 && tNext >= tDown) {
            fire = 0;
        }
        /* New event? */
        if (next != 0 && tNext >= next->start) {
            fire = next;
            next = next->next;
            tDown = fire->start + fire->duration;
            if (fire->period > 0) {
                if (fire->multiplier != 1) {
                    if (fire->multiplier > 1) fire->multiplier--;
                    fire->start += fire->period;
                    next = PacingEvent_Schedule(next, fire);
                } else {
                    fire->period = 0;
                }
            }
        }
        /* Set next time */
        tNext = tMax;
        if (fire != 0 && tDown < tNext) tNext = tDown;
        if (next != 0 && next->start < tNext) tNext = next->start;
    }
    if (fire != 0) {
        pace = fire->level;
    } else {
        pace = 0.0;
    }

    /* Set simulation starting time */
    t = tMin;

    /* Create solver */
    #if SUNDIALS_VERSION_MAJOR >= 6
    cvode_mem = CVodeCreate(CV_BDF, sundials_context);
    #elif SUNDIALS_VERSION_MAJOR >= 4
    cvode_mem = CVodeCreate(CV_BDF);
    #else
    cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
    #endif
    if (check_flag((void*)cvode_mem, "CVodeCreate", 0)) goto error;
    flag = CVodeInit(cvode_mem, rhs, t, y);
    if (check_flag(&flag, "CVodeInit", 1)) goto error;

    #if SUNDIALS_VERSION_MAJOR >= 6
        /* Create dense matrix for use in linear solves */
        sundials_dense_matrix = SUNDenseMatrix(N_STATE, N_STATE, sundials_context);
        if (check_flag((void *)sundials_dense_matrix, "SUNDenseMatrix", 0)) goto error;

        /* Create dense linear solver object with matrix */
        sundials_linear_solver = SUNLinSol_Dense(y, sundials_dense_matrix, sundials_context);
        if (check_flag((void *)sundials_linear_solver, "SUNLinSol_Dense", 0)) goto error;

        /* Attach the matrix and solver to cvode */
        flag = CVodeSetLinearSolver(cvode_mem, sundials_linear_solver, sundials_dense_matrix);
        if (check_flag(&flag, "CVodeSetLinearSolver", 1)) goto error;
    #elif SUNDIALS_VERSION_MAJOR >= 4
        /* Create dense matrix for use in linear solves */
        sundials_dense_matrix = SUNDenseMatrix(N_STATE, N_STATE);
        if (check_flag((void *)sundials_dense_matrix, "SUNDenseMatrix", 0)) goto error;

        /* Create dense linear solver object with matrix */
        sundials_linear_solver = SUNLinSol_Dense(y, sundials_dense_matrix);
        if (check_flag((void *)sundials_linear_solver, "SUNLinSol_Dense", 0)) goto error;

        /* Attach the matrix and solver to cvode */
        flag = CVodeSetLinearSolver(cvode_mem, sundials_linear_solver, sundials_dense_matrix);
        if (check_flag(&flag, "CVodeSetLinearSolver", 1)) goto error;
    #elif SUNDIALS_VERSION_MAJOR >= 3
        /* Create dense matrix for use in linear solves */
        sundials_dense_matrix = SUNDenseMatrix(N_STATE,N_STATE);
        if(check_flag((void *)sundials_dense_matrix, "SUNDenseMatrix", 0)) goto error;

        /* Create dense linear solver object with matrix */
        sundials_linear_solver = SUNDenseLinearSolver(y, sundials_dense_matrix);
        if(check_flag((void *)sundials_linear_solver, "SUNDenseLinearSolver", 0)) goto error;

        /* Attach the matrix and linear solver to cvode */
        flag = CVDlsSetLinearSolver(cvode_mem, sundials_linear_solver, sundials_dense_matrix);
        if(check_flag(&flag, "CVDlsSetLinearSolver", 1)) goto error;
    #else
        /* Create dense matrix for use in linear solves */
        flag = CVDense(cvode_mem, N_STATE);
        if (check_flag(&flag, "CVDense", 1)) goto error;
    #endif

    /* Set tolerances */
    double reltol = RCONST(1.0e-4);
    double abstol = RCONST(1.0e-6);
    flag = CVodeSStolerances(cvode_mem, reltol, abstol);
    if (check_flag(&flag, "CVodeSStolerances", 1)) goto error;

    /* Go! */
    while(1) {
        if (tMax < tNext) tNext = tMax;
        flag = CVode(cvode_mem, tNext, y, &t, CV_ONE_STEP);
        if (check_flag(&flag, "CVode", 1)) break;
        if (flag == CV_SUCCESS) {
            /* Shot past next discontinuity? */
            if (t > tNext) {
                flag = CVodeGetDky(cvode_mem, tNext, 0, y);
                if (check_flag(&flag, "CVodeGetDky", 1)) goto error;
                t = tNext;
                flag = CVodeReInit(cvode_mem, t, y);
                if (check_flag(&flag, "CVodeReInit", 1)) goto error;
                /* Recalculate logging values at this point in time */
                rhs(t, y, dy, 0);
            }
            /* Event over? */
            if (fire != 0 && t >= tDown) {
                pace = 0;
                fire = 0;
                flag = CVodeReInit(cvode_mem, t, y);
                if (check_flag(&flag, "CVodeReInit", 1)) goto error;
            }
            /* New event? */
            if (next != 0 && t >= next->start) {
                fire = next;
                next = next->next;
                pace = fire->level;
                tDown = fire->start + fire->duration;
                if (fire->period > 0) {
                    if (fire->multiplier == 1) {
                        fire->period = 0;
                    } else {
                        if (fire->multiplier > 1) fire->multiplier--;
                        fire->start += fire->period;
                        next = PacingEvent_Schedule(next, fire);
                    }
                }
                flag = CVodeReInit(cvode_mem, t, y);
                if (check_flag(&flag, "CVodeReInit", 1)) goto error;
            }
            /* Set next time */
            tNext = tMax;
            if (fire != 0 && tDown < tNext) tNext = tDown;
            if (next != 0 && next->start < tNext) tNext = next->start;
            /* Log */
            if (t >= tLog) {
                /* Log current position */
                PrintOutput(t, STATES[0]);
            }
        }
        if (t >= tMax) break;
    }

    /* Success! */
    success = 0;

error:
    /* Free allocated space */
    free(events);

    /* Free CVODE space */
    N_VDestroy_Serial(y);
    N_VDestroy_Serial(dy);
    CVodeFree(&cvode_mem);
    #if SUNDIALS_VERSION_MAJOR >= 6
    SUNContext_Free(&sundials_context);
    #endif

    /* Return */
    return success;
}
