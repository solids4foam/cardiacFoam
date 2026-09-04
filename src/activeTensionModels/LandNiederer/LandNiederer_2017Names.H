/*---------------------------------------------------------------------------*\
Description
    Index definitions for the original Land et al. human contraction model.

    This is the full 7-state formulation corresponding to the released
    Land2016.m / Land et al. (2017) human model.

    ORIGINAL INPUT CONVENTIONS:
      AV_Cai         : intracellular Ca2+ [uM]
      AV_lambda      : sarcomere stretch ratio [-]
      AV_lambda_rate : d(lambda)/dt [ms^-1]

    ORIGINAL OUTPUTS:
      AV_Ta : active tension  [kPa]
      AV_Tp : passive tension [kPa]
      AV_T  : total tension = Ta + Tp [kPa]

Reference
    Land S, Park-Holohan SJ, Smith NP, dos Remedios CG, Kentish JC,
    Niederer SA (2017).
    A model of cardiac contraction based on novel measurements of tension
    development in human cardiomyocytes.
    Journal of Molecular and Cellular Cardiology 106:68-83.
    DOI: 10.1016/j.yjmcc.2017.03.008
\*---------------------------------------------------------------------------*/

#ifndef LandNiederer_2017Names_H
#define LandNiederer_2017Names_H

enum STATES_INDEX
{
    XS = 0,
    XW,
    TRPN,
    TmBlocked,
    ZETAS,
    ZETAW,
    Cd,
    NUM_STATES
};

enum ALGEBRAIC_INDEX
{
    AV_XU = 0,
    AV_gamma_rate,
    AV_gamma_rate_w,
    AV_xb_uw,
    AV_xb_ws,
    AV_xb_su,
    AV_xb_wu,
    AV_xb_su_gamma,
    AV_xb_wu_gamma,
    AV_ca50,
    AV_Lfac,
    AV_dCd_dt,
    AV_Fd,
    AV_F1,
    AV_Ta,
    AV_Tp,
    AV_T,

    // External inputs
    AV_Cai,
    AV_lambda,
    AV_lambda_rate,

    NUM_ALGEBRAIC
};

enum CONSTANTS_INDEX
{
    AC_perm50 = 0,
    AC_TRPN_n,
    AC_koff,
    AC_dr,
    AC_wfrac,
    AC_TOT_A,
    AC_ktm_unblock,
    AC_beta_1,
    AC_beta_0,
    AC_gamma,
    AC_gamma_wu,
    AC_phi,

    // intact-human values
    AC_nperm,
    AC_ca50_ref,
    AC_Tref,
    AC_nu,
    AC_mu,

    // derived active-model constants
    AC_k_ws,
    AC_k_uw,
    AC_cdw,
    AC_cds,
    AC_k_wu,
    AC_k_su,
    AC_A,
    AC_XSSS,
    AC_XWSS,
    AC_ktm_block,

    // passive model
    AC_par_k,
    AC_b,
    AC_eta_l,
    AC_eta_s,
    AC_passive_a,

    NUM_CONSTANTS
};

inline void LandNiederer2017initConsts
(
    double* CONSTANTS,
    double* RATES,
    double* STATES
);

#ifndef CARDIAC_HOST_DEVICE
#define CARDIAC_HOST_DEVICE
#endif

CARDIAC_HOST_DEVICE inline void LandNiederer2017computeVariables
(
    double VOI,
    double* CONSTANTS,
    double* RATES,
    double* STATES,
    double* ALGEBRAIC
);

#endif
