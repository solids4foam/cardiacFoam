/*---------------------------------------------------------------------------*\
Description
    Original full Land human cardiac contraction model translated into the
    same self-contained C++ style used by the supplied cardiacFoam headers.

    This file intentionally preserves the released Land model conventions
    rather than the TWorld coupling conventions.

    IMPORTANT ORIGINAL UNITS
      Time / ODE rates : ms
      Cai              : uM
      lambda           : dimensionless
      lambda_rate      : ms^-1
      tension          : kPa

    Therefore, if a 3-D mechanics solver provides lambda_rate in s^-1:
        lambda_rate_for_Land = lambda_rate_3D / 1000.0;

    If an ionic model provides Cai in mM:
        Cai_for_Land_uM = Cai_ionic_mM * 1000.0;

    Those conversions belong at the coupling interface if the objective is
    to keep this model itself identical in convention to the original.

    Full state vector:
      XS, XW, TRPN, TmBlocked, ZETAS, ZETAW, Cd

    The seventh state Cd is the original passive viscoelastic distortion
    state. Total tension is T = Ta + Tp.

Reference
    Land S, Park-Holohan SJ, Smith NP, dos Remedios CG, Kentish JC,
    Niederer SA (2017).
    A model of cardiac contraction based on novel measurements of tension
    development in human cardiomyocytes.
    Journal of Molecular and Cellular Cardiology 106:68-83.
    DOI: 10.1016/j.yjmcc.2017.03.008
\*---------------------------------------------------------------------------*/

#ifndef LandNiederer_2017_H
#define LandNiederer_2017_H

#include <cmath>
#include <algorithm>

#include "LandNiederer_2017Names.H"

static const char* LandNiedererSTATES_NAMES[NUM_STATES] =
{
    "XS",
    "XW",
    "TRPN",
    "TmBlocked",
    "ZETAS",
    "ZETAW",
    "Cd"
};

static const char* LandNiedererALGEBRAIC_NAMES[NUM_ALGEBRAIC] =
{
    "AV_XU",
    "AV_gamma_rate",
    "AV_gamma_rate_w",
    "AV_xb_uw",
    "AV_xb_ws",
    "AV_xb_su",
    "AV_xb_wu",
    "AV_xb_su_gamma",
    "AV_xb_wu_gamma",
    "AV_ca50",
    "AV_Lfac",
    "AV_dCd_dt",
    "AV_Fd",
    "AV_F1",
    "AV_Ta",
    "AV_Tp",
    "AV_T",
    "AV_Cai",
    "AV_lambda",
    "AV_lambda_rate"
};

static const char* LandNiedererCONSTANTS_NAMES[NUM_CONSTANTS] =
{
    "AC_perm50",
    "AC_TRPN_n",
    "AC_koff",
    "AC_dr",
    "AC_wfrac",
    "AC_TOT_A",
    "AC_ktm_unblock",
    "AC_beta_1",
    "AC_beta_0",
    "AC_gamma",
    "AC_gamma_wu",
    "AC_phi",
    "AC_nperm",
    "AC_ca50_ref",
    "AC_Tref",
    "AC_nu",
    "AC_mu",
    "AC_k_ws",
    "AC_k_uw",
    "AC_cdw",
    "AC_cds",
    "AC_k_wu",
    "AC_k_su",
    "AC_A",
    "AC_XSSS",
    "AC_XWSS",
    "AC_ktm_block",
    "AC_par_k",
    "AC_b",
    "AC_eta_l",
    "AC_eta_s",
    "AC_passive_a"
};

inline void LandNiederer2017initConsts
(
    double* CONSTANTS,
    double* RATES,
    double* STATES
)
{
    // -----------------------------------------------------------------
    // Original Land active-model parameters
    // -----------------------------------------------------------------
    CONSTANTS[AC_perm50]      = 0.35;
    CONSTANTS[AC_TRPN_n]      = 2.0;
    CONSTANTS[AC_koff]        = 0.1;
    CONSTANTS[AC_dr]          = 0.25;
    CONSTANTS[AC_wfrac]       = 0.5;
    CONSTANTS[AC_TOT_A]       = 25.0;
    CONSTANTS[AC_ktm_unblock] = 1.0;
    CONSTANTS[AC_beta_1]      = -2.4;
    CONSTANTS[AC_beta_0]      = 2.3;
    CONSTANTS[AC_gamma]       = 0.0085;
    CONSTANTS[AC_gamma_wu]    = 0.615;
    CONSTANTS[AC_phi]         = 2.23;

    // -----------------------------------------------------------------
    // Original INTACT human parameter branch
    // -----------------------------------------------------------------
    CONSTANTS[AC_nperm]      = 5.0;
    CONSTANTS[AC_ca50_ref]   = 0.805; // uM
    CONSTANTS[AC_Tref]       = 120.0; // kPa
    CONSTANTS[AC_nu]         = 7.0;
    CONSTANTS[AC_mu]         = 3.0;

    // -----------------------------------------------------------------
    // Derived active-model constants: same algebra as released MATLAB
    // -----------------------------------------------------------------
    CONSTANTS[AC_k_ws] =
        0.004 * CONSTANTS[AC_mu];

    CONSTANTS[AC_k_uw] =
        0.026 * CONSTANTS[AC_nu];

    CONSTANTS[AC_cdw] =
        CONSTANTS[AC_phi]
      * CONSTANTS[AC_k_uw]
      * (1.0 - CONSTANTS[AC_dr])
      * (1.0 - CONSTANTS[AC_wfrac])
      / ((1.0 - CONSTANTS[AC_dr]) * CONSTANTS[AC_wfrac]);

    CONSTANTS[AC_cds] =
        CONSTANTS[AC_phi]
      * CONSTANTS[AC_k_ws]
      * (1.0 - CONSTANTS[AC_dr])
      * CONSTANTS[AC_wfrac]
      / CONSTANTS[AC_dr];

    CONSTANTS[AC_k_wu] =
        CONSTANTS[AC_k_uw]
      * (1.0 / CONSTANTS[AC_wfrac] - 1.0)
      - CONSTANTS[AC_k_ws];

    CONSTANTS[AC_k_su] =
        CONSTANTS[AC_k_ws]
      * (1.0 / CONSTANTS[AC_dr] - 1.0)
      * CONSTANTS[AC_wfrac];

    CONSTANTS[AC_A] =
        (0.25 * CONSTANTS[AC_TOT_A])
      / ((1.0 - CONSTANTS[AC_dr]) * CONSTANTS[AC_wfrac]
         + CONSTANTS[AC_dr])
      * (CONSTANTS[AC_dr] / 0.25);

    CONSTANTS[AC_XSSS] =
        CONSTANTS[AC_dr] * 0.5;

    CONSTANTS[AC_XWSS] =
        (1.0 - CONSTANTS[AC_dr])
      * CONSTANTS[AC_wfrac]
      * 0.5;

    CONSTANTS[AC_ktm_block] =
        CONSTANTS[AC_ktm_unblock]
      * std::pow(CONSTANTS[AC_perm50], CONSTANTS[AC_nperm])
      * 0.5
      / (0.5 - CONSTANTS[AC_XSSS] - CONSTANTS[AC_XWSS]);

    // -----------------------------------------------------------------
    // Original passive-model parameters
    // -----------------------------------------------------------------
    CONSTANTS[AC_par_k]     = 7.0;
    CONSTANTS[AC_b]         = 9.1;
    CONSTANTS[AC_eta_l]     = 200.0;
    CONSTANTS[AC_eta_s]     = 20.0;
    CONSTANTS[AC_passive_a] = 2.1;

    // Original reasonable initial condition:
    // [XS XW TRPN TmBlocked ZETAS ZETAW Cd] = [0 0 0 1 0 0 0]
    STATES[XS]        = 0.0;
    STATES[XW]        = 0.0;
    STATES[TRPN]      = 0.0;
    STATES[TmBlocked] = 1.0;
    STATES[ZETAS]     = 0.0;
    STATES[ZETAW]     = 0.0;
    STATES[Cd]        = 0.0;

    for (int i = 0; i < NUM_STATES; ++i)
    {
        RATES[i] = 0.0;
    }
}

#ifndef CARDIAC_HOST_DEVICE
#define CARDIAC_HOST_DEVICE
#endif

CARDIAC_HOST_DEVICE inline void LandNiederer2017computeVariables
(
    double /*VOI*/,
    double* CONSTANTS,
    double* RATES,
    double* STATES,
    double* ALGEBRAIC
)
{
    // -----------------------------------------------------------------
    // Inputs in ORIGINAL Land units
    // -----------------------------------------------------------------
    const double Cai = ALGEBRAIC[AV_Cai];                 // uM
    double lambda = ALGEBRAIC[AV_lambda];                 // -
    const double dlambda_dt = ALGEBRAIC[AV_lambda_rate];  // ms^-1

    // Original code: lambda = min(1.2,lambda)
    lambda = (lambda <= 1.2) ? lambda : 1.2;

    // Original state handling
    const double xs =
        (STATES[XS] >= 0.0) ? STATES[XS] : 0.0;

    const double xw =
        (STATES[XW] >= 0.0) ? STATES[XW] : 0.0;

    const double trpn =
        (STATES[TRPN] >= 0.0) ? STATES[TRPN] : 0.0;

    const double tmBlocked = STATES[TmBlocked];
    const double zetas = STATES[ZETAS];
    const double zetaw = STATES[ZETAW];

    // -----------------------------------------------------------------
    // Length dependence
    // -----------------------------------------------------------------
    const double lambda_min_087 =
        (0.87 <= lambda) ? 0.87 : lambda;

    ALGEBRAIC[AV_Lfac] =
        1.0
      + CONSTANTS[AC_beta_0]
      * (lambda + lambda_min_087 - 1.87);

    if (ALGEBRAIC[AV_Lfac] < 0.0)
    {
        ALGEBRAIC[AV_Lfac] = 0.0;
    }

    ALGEBRAIC[AV_XU] =
        (1.0 - tmBlocked) - xw - xs;

    // -----------------------------------------------------------------
    // Crossbridge model
    // -----------------------------------------------------------------
    ALGEBRAIC[AV_xb_ws] =
        CONSTANTS[AC_k_ws] * xw;

    ALGEBRAIC[AV_xb_uw] =
        CONSTANTS[AC_k_uw] * ALGEBRAIC[AV_XU];

    ALGEBRAIC[AV_xb_wu] =
        CONSTANTS[AC_k_wu] * xw;

    ALGEBRAIC[AV_xb_su] =
        CONSTANTS[AC_k_su] * xs;

    const double positive_branch =
        (zetas > 0.0) ? zetas : 0.0;

    const double negative_branch =
        (zetas < -1.0) ? (-zetas - 1.0) : 0.0;

    ALGEBRAIC[AV_gamma_rate] =
        CONSTANTS[AC_gamma]
      * ((positive_branch >= negative_branch)
         ? positive_branch
         : negative_branch);

    ALGEBRAIC[AV_xb_su_gamma] =
        ALGEBRAIC[AV_gamma_rate] * xs;

    ALGEBRAIC[AV_gamma_rate_w] =
        CONSTANTS[AC_gamma_wu] * std::fabs(zetaw);

    ALGEBRAIC[AV_xb_wu_gamma] =
        ALGEBRAIC[AV_gamma_rate_w] * xw;

    RATES[XS] =
        ALGEBRAIC[AV_xb_ws]
      - ALGEBRAIC[AV_xb_su]
      - ALGEBRAIC[AV_xb_su_gamma];

    RATES[XW] =
        ALGEBRAIC[AV_xb_uw]
      - ALGEBRAIC[AV_xb_wu]
      - ALGEBRAIC[AV_xb_ws]
      - ALGEBRAIC[AV_xb_wu_gamma];

    // -----------------------------------------------------------------
    // Calcium -> troponin
    // Original intact human reference: ca50 = 0.805 uM
    // -----------------------------------------------------------------
    const double dl =
        lambda - 1.0;

    const double dl_for_ca50 =
        (dl <= 0.2) ? dl : 0.2;

    ALGEBRAIC[AV_ca50] =
        CONSTANTS[AC_ca50_ref]
      + CONSTANTS[AC_beta_1] * dl_for_ca50;

    RATES[TRPN] =
        CONSTANTS[AC_koff]
      * (std::pow(Cai / ALGEBRAIC[AV_ca50],
                  CONSTANTS[AC_TRPN_n])
         * (1.0 - trpn)
         - trpn);

    // -----------------------------------------------------------------
    // Tropomyosin blocking/unblocking
    // -----------------------------------------------------------------
    // Original MATLAB: min(100, TRPN^(-nperm/2)).
    // A positive floor only prevents undefined pow(0,negative); because of
    // the subsequent cap it gives the same limiting factor of 100 at TRPN=0.
    const double trpn_safe =
        (trpn > 1.0e-15) ? trpn : 1.0e-15;

    const double blocking_factor_raw =
        std::pow(
            trpn_safe,
            -CONSTANTS[AC_nperm] / 2.0);

    const double blocking_factor =
        (blocking_factor_raw <= 100.0)
        ? blocking_factor_raw
        : 100.0;

    RATES[TmBlocked] =
        CONSTANTS[AC_ktm_block]
      * blocking_factor
      * ALGEBRAIC[AV_XU]
      - CONSTANTS[AC_ktm_unblock]
      * std::pow(
            trpn_safe,
            CONSTANTS[AC_nperm] / 2.0)
      * tmBlocked;

    // -----------------------------------------------------------------
    // Velocity dependence
    // ORIGINAL: dlambda_dt must already be in ms^-1
    // -----------------------------------------------------------------
    RATES[ZETAS] =
        CONSTANTS[AC_A] * dlambda_dt
      - CONSTANTS[AC_cds] * zetas;

    RATES[ZETAW] =
        CONSTANTS[AC_A] * dlambda_dt
      - CONSTANTS[AC_cdw] * zetaw;

    // -----------------------------------------------------------------
    // Original passive viscoelastic model
    // -----------------------------------------------------------------
    const double C =
        lambda - 1.0;

    const double cd =
        STATES[Cd];

    const double eta =
        ((C - cd) < 0.0)
        ? CONSTANTS[AC_eta_s]
        : CONSTANTS[AC_eta_l];

    ALGEBRAIC[AV_dCd_dt] =
        CONSTANTS[AC_par_k]
      * (C - cd)
      / eta;

    RATES[Cd] =
        ALGEBRAIC[AV_dCd_dt];

    ALGEBRAIC[AV_Fd] =
        eta * ALGEBRAIC[AV_dCd_dt];

    ALGEBRAIC[AV_F1] =
        std::exp(CONSTANTS[AC_b] * C) - 1.0;

    ALGEBRAIC[AV_Tp] =
        CONSTANTS[AC_passive_a]
      * (ALGEBRAIC[AV_F1] + ALGEBRAIC[AV_Fd]);

    // -----------------------------------------------------------------
    // Active and total tension
    // -----------------------------------------------------------------
    ALGEBRAIC[AV_Ta] =
        ALGEBRAIC[AV_Lfac]
      * (CONSTANTS[AC_Tref] / CONSTANTS[AC_dr])
      * ((zetas + 1.0) * xs + zetaw * xw);

    ALGEBRAIC[AV_T] =
        ALGEBRAIC[AV_Ta]
      + ALGEBRAIC[AV_Tp];
}

#endif
