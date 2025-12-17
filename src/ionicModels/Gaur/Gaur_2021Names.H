#pragma once

enum CONSTANTS_INDEX{ 
    AC_Acap, AC_Ageo, AC_BSLmax, AC_BSRmax, AC_CaMKo, AC_ECl, AC_F, AC_GK1,
    AC_GKb, AC_GKr, AC_GKs, AC_GNa, AC_GNaL, AC_Gncx, AC_GpCa, AC_Gto,
    AC_H, AC_Khp, AC_Kki, AC_Kko, AC_KmBSL, AC_KmBSR, AC_KmCaAct, AC_KmCaM,
    AC_KmCaMK, AC_Kmgatp, AC_Knai0, AC_Knao0, AC_Knap, AC_Kxkur, AC_L, AC_MgADP,
    AC_MgATP, AC_PCa, AC_PCaK, AC_PCaNa, AC_PCab, AC_PKNa, AC_PNab, AC_Pnak,
    AC_R, AC_SOICR, AC_T, AC_a2, AC_a4, AC_aCaMK, AC_b1, AC_bCaMK,
    AC_cao, AC_cli, AC_clo, AC_cmdnmax, AC_csqnmax, AC_delta, AC_eP, AC_grelbarjsrol,
    AC_h10, AC_h101, AC_h11, AC_h1111, AC_h12, AC_h121, AC_k1, AC_k11,
    AC_k1m, AC_k1p, AC_k2, AC_k21, AC_k2m, AC_k2p, AC_k3m, AC_k3p2,
    AC_k4m, AC_k4p2, AC_k5, AC_k51, AC_kasymm, AC_kcaoff, AC_kcaon, AC_kmcmdn,
    AC_kmcsqn, AC_kmtrpn, AC_kna1, AC_kna2, AC_kna3, AC_ko, AC_nao, AC_pi,
    AC_qca, AC_qna, AC_rad, AC_tau_gap, AC_tau_hl, AC_tauoff, AC_tauon,
    AC_trpnmax, AC_vcell, AC_vcsr, AC_vhalf_d, AC_vhalff, AC_vjsr, AC_vmyo, AC_vmyo1,
    AC_vmyo1frac, AC_vmyo2, AC_vnsr, AC_vnsr1, AC_vnsr2, AC_vss, AC_wca, AC_wna,
    AC_wnaca, AC_zca, AC_zk, AC_zna, 

    stim_start, stim_duration, stim_amplitude,  
    stim_period_S1, nstim1,        
    stim_period_S2, nstim2,

    NUM_CONSTANTS
}; 

enum STATES_INDEX{
    cell_v,
    ionic_concentrations_nai,
    ionic_concentrations_nass,
    ionic_concentrations_ki,
    ionic_concentrations_kss,
    ionic_concentrations_cai,
    ionic_concentrations_cai2,
    ionic_concentrations_cass,
    ionic_concentrations_cansr,
    ionic_concentrations_cajsr,
    ionic_concentrations_cacsr,
    I_Na_m,
    I_Na_h,
    I_Na_j,
    INaL_ml,
    INaL_hl,
    ICaL_d,
    ICaL_fca,
    IKr_xr,
    IKs_xs1,
    IKs_xs2,
    ITo_aa,
    CICR_Jrel2,
    CICR_Jrel1,
    CaMK_CaMKt,
    CICR_tjsrol,
    CICR_A,
    ICaL_fs,
    ICaL_ff,

    NUM_STATES
};


enum ALGEBRAIC_INDEX{  
    Iion_cm,
    Istim,  
    AV_ICaK,    
    AV_ICaL,    
    AV_ICaNa,    
    AV_ICab,    
    AV_IK1,    
    AV_IKb,    
    AV_IKr,    
    AV_IKs,    
    AV_INa,    
    AV_INaCa_i,    
    AV_INaCa_ss, 
    AV_INaCa,   
    AV_INaK,    
    AV_INaL,    
    AV_INab,   
    AV_ITo,    
    AV_IpCa,
    // ---------- INa gating ----------
    AV_m_inf,
    AV_h_inf,
    AV_j_inf,
    AV_tau_m,
    AV_tau_h,
    AV_tau_j,
    // ---------- INaL gating ----------
    AV_ml_inf,
    AV_hl_inf,
    AV_tau_ml,
    // ---------- ICaL gating ----------
    AV_d_inf,
    AV_ff_inf,
    AV_fs_inf,
    AV_fca_inf,
    AV_tau_d,
    AV_tau_ff,
    AV_tau_fs,
    AV_tau_fca,
    // ---------- IKr gating ----------
    AV_xr_inf,
    AV_tau_xr,
    AV_rkr,
    // ---------- IKs gating ----------
    AV_xs1_inf,
    AV_xs2_inf,
    AV_tau_xs1,
    AV_tau_xs2,
    // ---------- Ito gating ----------
    AV_beta_aa,
    AV_alpha_aa,
    AV_rito2,
    // ---------- IKb gating ----------
    AV_xkb,

    NUM_ALGEBRAIC
};

void
GaurinitConsts(double* CONSTANTS, double* RATES, double *STATES, int tissueFlag, const Foam::dictionary& stimulus);
void
GaurcomputeVariables(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC,int tissueFlag, bool solveVmWithinODESolver);

