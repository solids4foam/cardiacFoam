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
    Lewalle et al. contraction model translated into the same
    self-contained C++ style used by the supplied cardiacFoam headers.

    The model is the Land et al. (2017) contraction model with two extra
    myosin OFF states (Boff, Uoff). The OFF->ON rate k1 is increased by the
    total tension T = Ta + Tp:

      k1 = k1_0 * (1 + kforce * max(T, 0)),   k2 = k2_0

    The length dependence of the force amplitude (beta0, beta1) is not part
    of this formulation: Ta has no h(lambda) factor and Ca50 does not depend
    on lambda.

    UNITS
      Time / ODE rates : s
      Cai              : uM
      lambda           : dimensionless
      lambda_rate      : s^-1
      tension          : kPa  (Tref, passive a in kPa; kforce in kPa^-1)

    If an ionic model provides Cai in mM:
        Cai_for_Lewalle_uM = Cai_ionic_mM * 1000.0;

    Full state vector:
      CaTRPN, TmBlocked (B), XS (S), XW (W), TmBlockedOff (Boff),
      TmUnblockedOff (Uoff), ZETAS, ZETAW, Cd

    The unbound fraction is U = 1 - B - S - W - Boff - Uoff. The default
    initial state is the steady state at lambda = 1 and Cai = 0.1 uM of the
    default constants. Cd is the
    passive viscoelastic distortion state; total tension is T = Ta + Tp.

    DEFAULT CONSTANTS
      The defaults are the skinned human myocardium calibration of
      Lewalle et al. (2024), Fig. 6 (fitted to Awinda et al. 2020):
      [Ca2+]50 = 5.62 uM and Tref = 23 kPa. They are not an intact-tissue
      calibration. Intact constants are supplied through the "constants"
      sub-dictionary.

      With these defaults the system is stiff at diastolic Cai: the
      tropomyosin blocking rate ktm_block * CaTRPN^(-nperm/2) is of order
      1e7 s^-1 at Cai = 0.1 uM, because the factor CaTRPN^(-nperm/2) is not
      capped. Land et al. (2017) caps the same factor at 100.

    INTACT CALIBRATION
      No intact-tissue constant set is provided. Land et al. (2017) has a
      skinned branch (nperm 2.2, [Ca2+]50 2.5 uM, Tref 40.5 kPa, nu 1, mu 1)
      and an intact branch (nperm 5, [Ca2+]50 0.805 uM, Tref 120 kPa,
      nu 7, mu 3; the LandNiederer defaults). With the steady-state
      solution of the reference implementation, at SL 1.8 and 2.2 um:

      - The paper's intact-emulating values (KOFF0 = k1_0/k2_0 = 0.0158,
        kforce = 1.48e-4 Pa^-1, Tref = 109 kPa, [Ca2+]50 = 1.17 uM) on the
        Land intact branch with beta0 = beta1 = 0 give Fmax = 10 kPa against
        119 kPa for Land intact, and do not reproduce the Land skinned branch
        either.
      - Fitting kforce, Tref and [Ca2+]50 to Land intact reproduces its Fmax
        length dependence (119 -> 174 kPa) but only about a third of its
        [Ca2+]50 shift (0.47 -> 0.39 uM against 0.59 -> 0.24 uM), and needs
        Tref of about 940 kPa at KOFF0 = 1e-3. Larger KOFF0 matches worse.
        Force feedback on the OFF->ON rate alone does not reproduce the
        beta1 calcium-sensitivity shift of Land intact.

      The paper gives only the ratio k1_0/k2_0 for its intact-emulating set,
      not k1_0 itself, so twitch kinetics are not constrained by it.

Reference
    Lewalle A, Milburn GN, Campbell KS, Niederer SA (2024).
    Cardiac length-dependent activation driven by force-dependent
    thick-filament dynamics.
    Biophysical Journal 123(18).
    DOI: 10.1016/j.bpj.2024.05.025

    Land S, Park-Holohan SJ, Smith NP, dos Remedios CG, Kentish JC,
    Niederer SA (2017).
    A model of cardiac contraction based on novel measurements of tension
    development in human cardiomyocytes.
    Journal of Molecular and Cellular Cardiology 106:68-83.
    DOI: 10.1016/j.yjmcc.2017.03.008

Author
    Simao Nieto de Castro, UCD.
\*---------------------------------------------------------------------------*/

#ifndef LewalleNiederer_2024_H
#define LewalleNiederer_2024_H

#include <cmath>
#include <algorithm>

#include "LewalleNiederer_2024Names.H"

static const char* LewalleNiedererSTATES_NAMES[NUM_STATES] =
{
    "CaTRPN",
    "TmBlocked",
    "XS",
    "XW",
    "TmBlockedOff",
    "TmUnblockedOff",
    "ZETAS",
    "ZETAW",
    "Cd"
};

static const char* LewalleNiedererALGEBRAIC_NAMES[NUM_ALGEBRAIC] =
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
    "AV_dCd_dt",
    "AV_F1",
    "AV_F2",
    "AV_Ta",
    "AV_Tp",
    "AV_T",
    "AV_k1",
    "AV_k2",
    "AV_Cai",
    "AV_lambda",
    "AV_lambda_rate"
};

static const char* LewalleNiedererCONSTANTS_NAMES[NUM_CONSTANTS] =
{
    "AC_TRPN_n",
    "AC_k_trpn_on",
    "AC_k_trpn_off",
    "AC_ca50_ref",
    "AC_perm50",
    "AC_nperm",
    "AC_ktm_unblock",
    "AC_dr",
    "AC_wfrac",
    "AC_TOT_A",
    "AC_gamma",
    "AC_gamma_wu",
    "AC_phi",
    "AC_Tref",
    "AC_k_uw",
    "AC_k_ws",
    "AC_k1_0",
    "AC_k2_0",
    "AC_kforce",
    "AC_k_wu",
    "AC_k_su",
    "AC_cdw",
    "AC_cds",
    "AC_A",
    "AC_ktm_block",
    "AC_par_k",
    "AC_b",
    "AC_eta_l",
    "AC_eta_s",
    "AC_passive_a"
};

inline void LewalleNiederer2024initConsts
(
    double* CONSTANTS,
    double* RATES,
    double* STATES
)
{
    // -----------------------------------------------------------------
    // Calcium-troponin
    // -----------------------------------------------------------------
    CONSTANTS[AC_TRPN_n]      = 2.58;
    CONSTANTS[AC_k_trpn_on]   = 100.0;                // s^-1
    CONSTANTS[AC_k_trpn_off]  = 100.0;                // s^-1
    CONSTANTS[AC_ca50_ref]    = std::pow(10.0, -5.25) * 1.0e6; // uM

    // -----------------------------------------------------------------
    // Tropomyosin and crossbridge cycle
    // -----------------------------------------------------------------
    CONSTANTS[AC_perm50]      = 0.35;
    CONSTANTS[AC_nperm]       = 2.2;
    CONSTANTS[AC_ktm_unblock] = 1000.0;               // s^-1
    CONSTANTS[AC_dr]          = 0.25;
    CONSTANTS[AC_wfrac]       = 0.5;
    CONSTANTS[AC_TOT_A]       = 125.0;
    CONSTANTS[AC_gamma]       = 42.1;                 // s^-1
    CONSTANTS[AC_gamma_wu]    = 28.3;                 // s^-1
    CONSTANTS[AC_phi]         = 0.1498;
    CONSTANTS[AC_Tref]        = 23.0;                 // kPa
    CONSTANTS[AC_k_uw]        = 4.98;                 // s^-1
    CONSTANTS[AC_k_ws]        = 19.10;                // s^-1

    // -----------------------------------------------------------------
    // Myosin OFF states
    // -----------------------------------------------------------------
    CONSTANTS[AC_k1_0]        = 0.877;                // s^-1
    CONSTANTS[AC_k2_0]        = 12.6;                 // s^-1
    CONSTANTS[AC_kforce]      = 1.44064;              // kPa^-1

    // -----------------------------------------------------------------
    // Derived active-model constants
    // -----------------------------------------------------------------
    CONSTANTS[AC_k_wu] =
        CONSTANTS[AC_k_uw]
      * (1.0 / CONSTANTS[AC_wfrac] - 1.0)
      - CONSTANTS[AC_k_ws];

    CONSTANTS[AC_k_su] =
        CONSTANTS[AC_k_ws]
      * CONSTANTS[AC_wfrac]
      * (1.0 / CONSTANTS[AC_dr] - 1.0);

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

    CONSTANTS[AC_A] =
        CONSTANTS[AC_TOT_A]
      * CONSTANTS[AC_dr]
      / ((1.0 - CONSTANTS[AC_dr]) * CONSTANTS[AC_wfrac]
         + CONSTANTS[AC_dr]);

    CONSTANTS[AC_ktm_block] =
        CONSTANTS[AC_ktm_unblock]
      * std::pow(CONSTANTS[AC_perm50], CONSTANTS[AC_nperm])
      / ((1.0 - CONSTANTS[AC_dr]) * (1.0 - CONSTANTS[AC_wfrac]));

    // -----------------------------------------------------------------
    // Passive-model parameters
    // -----------------------------------------------------------------
    CONSTANTS[AC_par_k]     = 8.86;
    CONSTANTS[AC_b]         = 9.1;
    CONSTANTS[AC_eta_l]     = 0.2;                    // s
    CONSTANTS[AC_eta_s]     = 0.02;                   // s
    CONSTANTS[AC_passive_a] = 0.241;                  // kPa

    // Initial condition: steady state at lambda = 1, Cai = 0.1 uM
    STATES[CaTRPN]         = 3.054828e-05;
    STATES[TmBlocked]      = 0.065074;
    STATES[XS]             = 0.0;
    STATES[XW]             = 0.0;
    STATES[TmBlockedOff]   = 0.934926;
    STATES[TmUnblockedOff] = 0.0;
    STATES[ZETAS]          = 0.0;
    STATES[ZETAW]          = 0.0;
    STATES[Cd]             = 0.0;

    for (int i = 0; i < NUM_STATES; ++i)
    {
        RATES[i] = 0.0;
    }
}

#ifndef CARDIAC_HOST_DEVICE
#define CARDIAC_HOST_DEVICE
#endif

CARDIAC_HOST_DEVICE inline void LewalleNiederer2024computeVariables
(
    double /*VOI*/,
    double* CONSTANTS,
    double* RATES,
    double* STATES,
    double* ALGEBRAIC
)
{
    // -----------------------------------------------------------------
    // Inputs
    // -----------------------------------------------------------------
    const double Cai = ALGEBRAIC[AV_Cai];                 // uM
    const double lambda = ALGEBRAIC[AV_lambda];           // -
    const double dlambda_dt = ALGEBRAIC[AV_lambda_rate];  // s^-1

    const double trpn =
        (STATES[CaTRPN] >= 0.0) ? STATES[CaTRPN] : 0.0;

    const double b = STATES[TmBlocked];
    const double xs = STATES[XS];
    const double xw = STATES[XW];
    const double bOff = STATES[TmBlockedOff];
    const double uOff = STATES[TmUnblockedOff];
    const double zetas = STATES[ZETAS];
    const double zetaw = STATES[ZETAW];
    const double cd = STATES[Cd];

    ALGEBRAIC[AV_XU] = 1.0 - b - xs - xw - bOff - uOff;

    // -----------------------------------------------------------------
    // Passive viscoelastic model
    // -----------------------------------------------------------------
    const double C = lambda - 1.0;

    const double eta =
        ((C - cd) < 0.0)
        ? CONSTANTS[AC_eta_s]
        : CONSTANTS[AC_eta_l];

    ALGEBRAIC[AV_dCd_dt] =
        CONSTANTS[AC_par_k]
      * (C - cd)
      / eta;

    RATES[Cd] = ALGEBRAIC[AV_dCd_dt];

    ALGEBRAIC[AV_F1] =
        CONSTANTS[AC_passive_a]
      * (std::exp(CONSTANTS[AC_b] * C) - 1.0);

    ALGEBRAIC[AV_F2] =
        CONSTANTS[AC_passive_a]
      * CONSTANTS[AC_par_k]
      * (C - cd);

    ALGEBRAIC[AV_Tp] =
        ALGEBRAIC[AV_F1]
      + ALGEBRAIC[AV_F2];

    // -----------------------------------------------------------------
    // Active and total tension
    // -----------------------------------------------------------------
    ALGEBRAIC[AV_Ta] =
        (CONSTANTS[AC_Tref] / CONSTANTS[AC_dr])
      * ((zetas + 1.0) * xs + zetaw * xw);

    ALGEBRAIC[AV_T] =
        ALGEBRAIC[AV_Ta]
      + ALGEBRAIC[AV_Tp];

    // -----------------------------------------------------------------
    // Force feedback on the OFF->ON rate
    // -----------------------------------------------------------------
    ALGEBRAIC[AV_k1] =
        CONSTANTS[AC_k1_0]
      * (1.0
       + CONSTANTS[AC_kforce]
       * ((ALGEBRAIC[AV_T] > 0.0) ? ALGEBRAIC[AV_T] : 0.0));

    ALGEBRAIC[AV_k2] = CONSTANTS[AC_k2_0];

    const double k1 = ALGEBRAIC[AV_k1];
    const double k2 = ALGEBRAIC[AV_k2];

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
    // -----------------------------------------------------------------
    RATES[CaTRPN] =
        CONSTANTS[AC_k_trpn_on]
      * std::pow(Cai / CONSTANTS[AC_ca50_ref], CONSTANTS[AC_TRPN_n])
      * (1.0 - trpn)
      - CONSTANTS[AC_k_trpn_off] * trpn;

    // -----------------------------------------------------------------
    // Tropomyosin blocking/unblocking, ON and OFF myosin pools
    // -----------------------------------------------------------------
    // A positive floor only prevents undefined pow(0,negative).
    const double trpn_safe =
        (trpn > 1.0e-15) ? trpn : 1.0e-15;

    const double blocking_rate =
        CONSTANTS[AC_ktm_block]
      * std::pow(trpn_safe, -CONSTANTS[AC_nperm] / 2.0);

    const double unblocking_rate =
        CONSTANTS[AC_ktm_unblock]
      * std::pow(trpn_safe, CONSTANTS[AC_nperm] / 2.0);

    RATES[TmBlocked] =
        blocking_rate * ALGEBRAIC[AV_XU]
      - unblocking_rate * b
      - k2 * b
      + k1 * bOff;

    RATES[TmBlockedOff] =
        blocking_rate * uOff
      - unblocking_rate * bOff
      + k2 * b
      - k1 * bOff;

    RATES[TmUnblockedOff] =
      - blocking_rate * uOff
      + unblocking_rate * bOff
      + k2 * ALGEBRAIC[AV_XU]
      - k1 * uOff;

    // -----------------------------------------------------------------
    // Velocity dependence
    // dlambda_dt must be in s^-1
    // -----------------------------------------------------------------
    RATES[ZETAS] =
        CONSTANTS[AC_A] * dlambda_dt
      - CONSTANTS[AC_cds] * zetas;

    RATES[ZETAW] =
        CONSTANTS[AC_A] * dlambda_dt
      - CONSTANTS[AC_cdw] * zetaw;
}

#endif
