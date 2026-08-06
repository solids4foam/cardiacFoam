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
    Index enumerations and function prototype declarations for the
    Land-Niederer 2017 myofilament contraction model.

    Extracted from the TWorld 2024 ionic model contraction subsystem.
    The electrophysiology is removed: this file only covers the 6 ODE
    states of the crossbridge/tropomyosin system plus their driving
    calcium signal and stretch inputs.

    Reference:
      Land S, Wang YT, Vicente-Steijn R et al. (2017).
      An analysis of deformation-dependent electromechanical coupling in
      the mouse heart. J Physiol, 595, 2547-2562.

Author
    Simao Nieto de Castro, UCD.
\*---------------------------------------------------------------------------*/

#ifndef LandNiederer_2017Names_H
#define LandNiederer_2017Names_H

// -------------------------------------------------------------------------
// State variables  (NUM_STATES = 6)
// -------------------------------------------------------------------------
enum STATES_INDEX
{
    Ca_TRPN = 0,   // Normalised Ca-troponin occupancy          [-]
    TmBlocked,     // Fraction of Tm-blocked regulatory units   [-]
    XW,            // Weak-binding (pre-force) crossbridges      [-]
    XS,            // Strong-binding (force-generating) XBs     [-]
    ZETAS,         // Velocity distortion, strong XBs            [-]
    ZETAW,         // Velocity distortion, weak XBs              [-]
    NUM_STATES     // = 6
};

// -------------------------------------------------------------------------
// Algebraic variables  (NUM_ALGEBRAIC = 16)
// -------------------------------------------------------------------------
enum ALGEBRAIC_INDEX
{
    AV_XU = 0,         // Unbound regulatory unit fraction  = 1-TmBlocked-XW-XS
    AV_gamma_rate,     // Velocity-dependent detachment rate, strong XBs
    AV_gamma_rate_w,   // Velocity-dependent detachment rate, weak XBs
    AV_xb_uw,         // Crossbridge flux U→W
    AV_xb_ws,         // Crossbridge flux W→S
    AV_xb_su,         // Crossbridge flux S→U (passive)
    AV_xb_wu,         // Crossbridge flux W→U (passive)
    AV_xb_su_gamma,   // Crossbridge flux S→U (velocity)
    AV_xb_wu_gamma,   // Crossbridge flux W→U (velocity)
    AV_Ta,            // Active tension  [kPa]
    AV_d_Ca_TRPN,     // dCa_TRPN/dt  (same as RATES[Ca_TRPN])
    AV_Cai,           // Intracellular Ca concentration (drive signal)  [mM]
    AV_lambda,        // Sarcomere stretch ratio  [-]
    AV_lambda_rate,   // Time-derivative of lambda  [s^-1]
    AV_Lfac,          // Length-dependent scaling factor (Frank-Starling) [-]
    AV_ca50,          // Effective half-activation concentration  [µM]
    NUM_ALGEBRAIC     // = 16
};

// -------------------------------------------------------------------------
// Constants  (NUM_CONSTANTS = 33)
//   First 21: user-facing parameters (set in initConsts, overridable by dict)
//   Last 12:  derived quantities (computed once in initConsts)
// -------------------------------------------------------------------------
enum CONSTANTS_INDEX
{
    // User-facing parameters
    AC_TOT_A = 0,          // Total crossbridge amplitude  [kPa]         = 25.0
    AC_TRPN_n,             // Hill coefficient, Ca-troponin               = 1.65
    AC_Tref,               // Reference active tension  [kPa]             = 80.0
    AC_beta_0,             // Length-dependence slope                     = 2.3
    AC_beta_1,             // Ca50 length-sensitivity                     = -2.4
    AC_dr,                 // Duty ratio                                  = 0.25
    AC_fracTnIpo,          // Fraction TnI phosphorylated at baseline     = 0.0031
    AC_contraction_gamma,  // Velocity detachment rate constant (strong)  = 0.0085
    AC_gamma_wu,           // Velocity detachment rate constant (weak)    = 0.615
    AC_koff,               // Ca-troponin off-rate  [ms^-1]               = 0.07854
    AC_ktm_unblock,        // Tropomyosin unblocking rate  [ms^-1]        = 0.02626
    AC_lambda_max,         // Maximum stretch for Lfac computation        = 1.2
    AC_lambda_min,         // Minimum stretch for Lfac computation        = 0.87
    AC_mu,                 // Crossbridge kinetics factor (W→S)           = 3.94046
    AC_nperm,              // Hill coefficient, TmBlocked permissivity    = 2.036
    AC_nu,                 // Crossbridge kinetics factor (U→W)           = 10.15996
    AC_perm50,             // Half-activation permissivity                = 0.35
    AC_phi,                // Ratio of detachment to attachment rates     = 2.23
    AC_wfrac,              // Weak-binding fraction                       = 0.5
    AC_fMyBPC_PKA,         // MyBPC phosphorylation fraction (PKA)        = 0.0
    AC_fTnI_PKA,           // TnI phosphorylation fraction (PKA)          = 0.0

    // Derived quantities (computed once from the above in initConsts)
    AC_A,                  // Effective crossbridge amplitude
    AC_XSSS,               // Steady-state strong XB fraction
    AC_XWSS,               // Steady-state weak XB fraction
    AC_fPKA_TnI,           // Effective TnI PKA factor
    AC_k_uw,               // Transition rate U→W  [ms^-1]
    AC_k_ws,               // Transition rate W→S  [ms^-1]
    AC_k_wu,               // Transition rate W→U  [ms^-1]
    AC_k_su,               // Transition rate S→U  [ms^-1]
    AC_cds,                // Velocity-coupling constant, strong XBs
    AC_cdw,                // Velocity-coupling constant, weak XBs
    AC_ktm_block,          // Tropomyosin blocking rate  [ms^-1]
    AC_PKAForceMultiplier, // PKA-mediated force enhancement factor

    NUM_CONSTANTS          // = 33
};

// -------------------------------------------------------------------------
// Function prototypes (implemented in LandNiederer_2017.H)
// -------------------------------------------------------------------------

//- Initialise constants and state initial conditions.
inline void LandNiederer2017initConsts
(
    double* CONSTANTS,
    double* RATES,
    double* STATES
);

#ifndef CARDIAC_HOST_DEVICE
#define CARDIAC_HOST_DEVICE
#endif

//- Compute algebraic variables and ODE rates.
//  ALGEBRAIC[AV_Cai], ALGEBRAIC[AV_lambda], and ALGEBRAIC[AV_lambda_rate]
//  must be set by the caller before this function is invoked.
CARDIAC_HOST_DEVICE inline void LandNiederer2017computeVariables
(
    double VOI,
    double* CONSTANTS,
    double* RATES,
    double* STATES,
    double* ALGEBRAIC
);

#endif

// ************************************************************************* //
