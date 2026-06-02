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
    Self-contained ODE equations for the Land-Niederer 2017 cardiac
    myofilament contraction model.

    This file is intentionally free of OpenFOAM headers so that it can be
    included from both the OpenFOAM wrapper and any standalone C++ harness.

    The equations are extracted from the TWorld 2024 ionic model contraction
    subsystem (src/ionicModels/TWorld/TWorld_2024.H, lines 1341-1358).
    The electrophysiology coupling is removed: the only external inputs are
      - ALGEBRAIC[AV_Cai]        : [Ca2+]_i  in mM   (set by caller)
      - ALGEBRAIC[AV_lambda]     : sarcomere stretch  (set by caller)
      - ALGEBRAIC[AV_lambda_rate]: d(lambda)/dt  [s^-1]  (set by caller)

    Reference:
      Land S, Wang YT, Vicente-Steijn R et al. (2017).
      An analysis of deformation-dependent electromechanical coupling in
      the mouse heart. J Physiol, 595, 2547-2562.

Author
    Simao Nieto de Castro, UCD.
\*---------------------------------------------------------------------------*/

#ifndef LandNiederer_2017_H
#define LandNiederer_2017_H

#include <cmath>
#include <algorithm>

#include "LandNiederer_2017Names.H"

// State variable names (for I/O, matched to STATES_INDEX order)
static const char* LandNiedererSTATES_NAMES[NUM_STATES] =
{
    "Ca_TRPN",
    "TmBlocked",
    "XW",
    "XS",
    "ZETAS",
    "ZETAW"
};

// Algebraic variable names (for I/O, matched to ALGEBRAIC_INDEX order)
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
    "AV_Ta",
    "AV_d_Ca_TRPN",
    "AV_Cai",
    "AV_lambda",
    "AV_lambda_rate",
    "AV_Lfac",
    "AV_ca50"
};

// Constant names (for I/O, matched to CONSTANTS_INDEX order)
static const char* LandNiedererCONSTANTS_NAMES[NUM_CONSTANTS] =
{
    "AC_TOT_A",
    "AC_TRPN_n",
    "AC_Tref",
    "AC_beta_0",
    "AC_beta_1",
    "AC_dr",
    "AC_fracTnIpo",
    "AC_contraction_gamma",
    "AC_gamma_wu",
    "AC_koff",
    "AC_ktm_unblock",
    "AC_lambda_max",
    "AC_lambda_min",
    "AC_mu",
    "AC_nperm",
    "AC_nu",
    "AC_perm50",
    "AC_phi",
    "AC_wfrac",
    "AC_fMyBPC_PKA",
    "AC_fTnI_PKA",
    "AC_A",
    "AC_XSSS",
    "AC_XWSS",
    "AC_fPKA_TnI",
    "AC_k_uw",
    "AC_k_ws",
    "AC_k_wu",
    "AC_k_su",
    "AC_cds",
    "AC_cdw",
    "AC_ktm_block",
    "AC_PKAForceMultiplier"
};


// =========================================================================
//  initConsts
//  Set default parameter values and initial state conditions.
//  Called once per model instantiation.
// =========================================================================
void LandNiederer2017initConsts
(
    double* CONSTANTS,
    double* RATES,
    double* STATES
)
{
    // -----------------------------------------------------------------
    // User-facing parameters
    // (All values from the TWorld contraction block, Land 2017 defaults)
    // -----------------------------------------------------------------

    CONSTANTS[AC_TOT_A]         = 25.0;
    CONSTANTS[AC_TRPN_n]        = 1.65;
    CONSTANTS[AC_Tref]          = 80.0;
    CONSTANTS[AC_beta_0]        = 2.3;
    CONSTANTS[AC_beta_1]        = -2.4;
    CONSTANTS[AC_dr]            = 0.25;
    CONSTANTS[AC_fracTnIpo]     = 0.0031;
    CONSTANTS[AC_contraction_gamma] = 0.0085;
    CONSTANTS[AC_gamma_wu]      = 0.615;
    CONSTANTS[AC_koff]          = 0.07854;
    CONSTANTS[AC_ktm_unblock]   = 0.02626;
    CONSTANTS[AC_lambda_max]    = 1.2;
    CONSTANTS[AC_lambda_min]    = 0.87;
    CONSTANTS[AC_mu]            = 3.94046;
    CONSTANTS[AC_nperm]         = 2.036;
    CONSTANTS[AC_nu]            = 10.15996;
    CONSTANTS[AC_perm50]        = 0.35;
    CONSTANTS[AC_phi]           = 2.23;
    CONSTANTS[AC_wfrac]         = 0.5;
    CONSTANTS[AC_fMyBPC_PKA]    = 0.0;
    CONSTANTS[AC_fTnI_PKA]      = 0.0;

    // -----------------------------------------------------------------
    // Derived constants (computed once from user-facing parameters)
    // -----------------------------------------------------------------

    // fPKA_TnI: effective TnI phosphorylation factor
    //   fPKA_TnI = 1.45 - 0.45 * (1 - fTnI_PKA) / (1 - fracTnIpo)
    CONSTANTS[AC_fPKA_TnI] =
        1.45
      - 0.45 * (1.0 - CONSTANTS[AC_fTnI_PKA])
             / (1.0 - CONSTANTS[AC_fracTnIpo]);

    // Crossbridge geometry
    CONSTANTS[AC_XSSS] = CONSTANTS[AC_dr] * 0.5;
    CONSTANTS[AC_XWSS] =
        (1.0 - CONSTANTS[AC_dr]) * CONSTANTS[AC_wfrac] * 0.5;

    // Effective crossbridge amplitude
    CONSTANTS[AC_A] =
        CONSTANTS[AC_TOT_A] * CONSTANTS[AC_dr]
      / (  (1.0 - CONSTANTS[AC_dr]) * CONSTANTS[AC_wfrac]
         + CONSTANTS[AC_dr]);

    // PKA-mediated force enhancement via MyBPC phosphorylation
    CONSTANTS[AC_PKAForceMultiplier] =
        1.0 + 0.26 * CONSTANTS[AC_fMyBPC_PKA];

    // Crossbridge transition rates
    CONSTANTS[AC_k_uw] = 0.026 * CONSTANTS[AC_nu];

    CONSTANTS[AC_k_ws] =
        0.004
      * (1.0 + CONSTANTS[AC_fMyBPC_PKA] / 2.0)
      * CONSTANTS[AC_mu];

    CONSTANTS[AC_k_wu] =
        CONSTANTS[AC_k_uw] * (1.0 / CONSTANTS[AC_wfrac] - 1.0)
      - CONSTANTS[AC_k_ws];

    CONSTANTS[AC_k_su] =
        CONSTANTS[AC_k_ws]
      * (1.0 / CONSTANTS[AC_dr] - 1.0)
      * CONSTANTS[AC_wfrac];

    // Velocity-coupling constants
    CONSTANTS[AC_cds] =
        CONSTANTS[AC_phi] * CONSTANTS[AC_k_ws]
      * CONSTANTS[AC_wfrac] * (1.0 - CONSTANTS[AC_dr])
      / CONSTANTS[AC_dr];

    CONSTANTS[AC_cdw] =
        CONSTANTS[AC_phi] * CONSTANTS[AC_k_uw]
      * (1.0 - CONSTANTS[AC_wfrac]) / CONSTANTS[AC_wfrac];

    // Blocking rate (detailed-balance constraint)
    CONSTANTS[AC_ktm_block] =
        CONSTANTS[AC_ktm_unblock]
      * std::pow(CONSTANTS[AC_perm50], CONSTANTS[AC_nperm])
      * 0.5
      / (0.5 - CONSTANTS[AC_XSSS] - CONSTANTS[AC_XWSS]);

    // -----------------------------------------------------------------
    // Initial state conditions
    //   At rest (zero tension, no crossbridges, Tm fully blocked)
    //   Users may override these from the case dictionary.
    // -----------------------------------------------------------------

    STATES[Ca_TRPN]   = 0.0;
    STATES[TmBlocked] = 1.0;
    STATES[XW]        = 0.0;
    STATES[XS]        = 0.0;
    STATES[ZETAS]     = 0.0;
    STATES[ZETAW]     = 0.0;

    // Initialise rates to zero
    for (int i = 0; i < NUM_STATES; ++i)
    {
        RATES[i] = 0.0;
    }
}


// =========================================================================
//  computeVariables
//
//  Compute all algebraic variables and ODE rates.
//
//  Pre-conditions (caller must set before invoking):
//    ALGEBRAIC[AV_Cai]        = [Ca2+]_i  in mM
//    ALGEBRAIC[AV_lambda]     = sarcomere stretch ratio  [-]
//    ALGEBRAIC[AV_lambda_rate]= d(lambda)/dt  [s^-1]
//
//  VOI is unused (model has no explicit time dependence beyond the inputs).
// =========================================================================
void LandNiederer2017computeVariables
(
    double /*VOI*/,
    double* CONSTANTS,
    double* RATES,
    double* STATES,
    double* ALGEBRAIC
)
{
    // -----------------------------------------------------------------
    // 1.  Length-dependent quantities (updated every timestep)
    // -----------------------------------------------------------------

    // Clamp lambda from above at lambda_max
    const double lambda  = ALGEBRAIC[AV_lambda];
    const double lambda0 = (lambda <= CONSTANTS[AC_lambda_max])
                         ? lambda
                         : CONSTANTS[AC_lambda_max];

    // Effective minimum stretch for this sarcomere length
    const double lambda0_min = (CONSTANTS[AC_lambda_min] <= lambda0)
                             ? CONSTANTS[AC_lambda_min]
                             : lambda0;

    // Frank-Starling length-dependence factor
    const double Lfac_val =
        1.0 + CONSTANTS[AC_beta_0]
            * (lambda0 + lambda0_min - (1.0 + CONSTANTS[AC_lambda_min]));

    ALGEBRAIC[AV_Lfac] = (Lfac_val >= 0.0) ? Lfac_val : 0.0;

    // Effective half-activation concentration [µM] (length-sensitive)
    const double dl = lambda - 1.0;
    const double dl_clamped = (dl <= 0.2) ? dl : 0.2;

    ALGEBRAIC[AV_ca50] =
        CONSTANTS[AC_fPKA_TnI] * 0.7645
      + CONSTANTS[AC_beta_1] * dl_clamped;

    // -----------------------------------------------------------------
    // 2.  Ca-troponin ODE rate
    //     Cai is in mM; ca50 is in µM → multiply Cai by 1000
    // -----------------------------------------------------------------

    const double cai_uM = ALGEBRAIC[AV_Cai] * 1000.0;   // mM → µM

    ALGEBRAIC[AV_d_Ca_TRPN] =
        CONSTANTS[AC_koff]
      * (  std::pow(cai_uM / ALGEBRAIC[AV_ca50], CONSTANTS[AC_TRPN_n])
         * (1.0 - STATES[Ca_TRPN])
         - STATES[Ca_TRPN]);

    RATES[Ca_TRPN] = ALGEBRAIC[AV_d_Ca_TRPN];

    // -----------------------------------------------------------------
    // 3.  Tropomyosin blocking/unblocking ODE rate
    // -----------------------------------------------------------------

    // Unbound regulatory unit fraction
    ALGEBRAIC[AV_XU] =
        1.0 - STATES[TmBlocked] - STATES[XW] - STATES[XS];

    // Blocking flux: capped at 100 to guard against Ca_TRPN → 0
    const double Ca_TRPN_safe = (STATES[Ca_TRPN] > 0.0)
                               ? STATES[Ca_TRPN]
                               : 1.0e-12;

    const double blocking_Ca_factor =
        std::pow(Ca_TRPN_safe, -CONSTANTS[AC_nperm] / 2.0);

    const double blocking_Ca_capped =
        (blocking_Ca_factor <= 100.0) ? blocking_Ca_factor : 100.0;

    RATES[TmBlocked] =
        CONSTANTS[AC_ktm_block]
      * blocking_Ca_capped
      * ALGEBRAIC[AV_XU]
      - CONSTANTS[AC_ktm_unblock]
      * std::pow(Ca_TRPN_safe, CONSTANTS[AC_nperm] / 2.0)
      * STATES[TmBlocked];

    // -----------------------------------------------------------------
    // 4.  Velocity-dependent detachment rates
    // -----------------------------------------------------------------

    // Strong-XB velocity term: gamma * ZETAS  (for ZETAS>0 or ZETAS<-1)
    if (STATES[ZETAS] > 0.0)
    {
        ALGEBRAIC[AV_gamma_rate] =
            CONSTANTS[AC_contraction_gamma] * STATES[ZETAS];
    }
    else if (STATES[ZETAS] < -1.0)
    {
        ALGEBRAIC[AV_gamma_rate] =
            CONSTANTS[AC_contraction_gamma] * (-STATES[ZETAS] - 1.0);
    }
    else
    {
        ALGEBRAIC[AV_gamma_rate] = 0.0;
    }

    // Weak-XB velocity term: gamma_wu * |ZETAW|
    ALGEBRAIC[AV_gamma_rate_w] =
        CONSTANTS[AC_gamma_wu] * std::fabs(STATES[ZETAW]);

    // -----------------------------------------------------------------
    // 5.  Crossbridge fluxes
    // -----------------------------------------------------------------

    ALGEBRAIC[AV_xb_uw] =
        CONSTANTS[AC_k_uw] * ALGEBRAIC[AV_XU]
      * (1.0 + CONSTANTS[AC_fMyBPC_PKA] / 2.0);       // U → W

    ALGEBRAIC[AV_xb_ws]       = CONSTANTS[AC_k_ws] * STATES[XW]; // W → S
    ALGEBRAIC[AV_xb_su]       = CONSTANTS[AC_k_su] * STATES[XS]; // S → U
    ALGEBRAIC[AV_xb_wu]       = CONSTANTS[AC_k_wu] * STATES[XW]; // W → U
    ALGEBRAIC[AV_xb_su_gamma] = ALGEBRAIC[AV_gamma_rate] * STATES[XS];
    ALGEBRAIC[AV_xb_wu_gamma] = ALGEBRAIC[AV_gamma_rate_w] * STATES[XW];

    // -----------------------------------------------------------------
    // 6.  Velocity-distortion ODEs (ZETAS, ZETAW)
    //     lambda_rate must be supplied by the caller in s^-1.
    //     The factors A, cds, cdw are already in ms^-1 units (from
    //     TWorld), but here we are working in SI (s^-1).
    //     The caller is responsible for unit consistency.
    // -----------------------------------------------------------------

    const double lambda_rate = ALGEBRAIC[AV_lambda_rate];

    RATES[ZETAS] =
        CONSTANTS[AC_A] * lambda_rate
      - CONSTANTS[AC_cds] * STATES[ZETAS];

    RATES[ZETAW] =
        CONSTANTS[AC_A] * lambda_rate
      - CONSTANTS[AC_cdw] * STATES[ZETAW];

    // -----------------------------------------------------------------
    // 7.  Crossbridge population ODEs
    // -----------------------------------------------------------------

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
    // 8.  Active tension (algebraic output)
    //     Ta = Lfac * PKAForceMultiplier * Tref/dr
    //          * [(ZETAS+1)*XS + ZETAW*XW]
    // -----------------------------------------------------------------

    ALGEBRAIC[AV_Ta] =
        ALGEBRAIC[AV_Lfac]
      * CONSTANTS[AC_PKAForceMultiplier]
      * CONSTANTS[AC_Tref]
      / CONSTANTS[AC_dr]
      * (  (STATES[ZETAS] + 1.0) * STATES[XS]
         + STATES[ZETAW] * STATES[XW]);
}

#endif

// ************************************************************************* //
