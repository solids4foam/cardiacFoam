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
    Index definitions for the Lewalle et al. (2024) contraction model: the
    Land et al. (2017) model extended with two myosin OFF states (Boff, Uoff)
    and force feedback on the OFF->ON rate.

    Full state vector (9 states):
      CaTRPN, TmBlocked (B), XS (S), XW (W), TmBlockedOff (Boff),
      TmUnblockedOff (Uoff), ZETAS, ZETAW, Cd

    CONVENTIONS:
      AV_Cai         : intracellular Ca2+ [uM]
      AV_lambda      : sarcomere stretch ratio [-]
      AV_lambda_rate : d(lambda)/dt [s^-1]
      Time / ODE rates : s
      Tension          : kPa (AC_Tref, AC_passive_a in kPa; AC_kforce in
                         kPa^-1)

    OUTPUTS:
      AV_Ta : active tension  [kPa]
      AV_Tp : passive tension [kPa]
      AV_T  : total tension = Ta + Tp [kPa]

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

#ifndef LewalleNiederer_2024Names_H
#define LewalleNiederer_2024Names_H

enum STATES_INDEX
{
    CaTRPN = 0,
    TmBlocked,
    XS,
    XW,
    TmBlockedOff,
    TmUnblockedOff,
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
    AV_dCd_dt,
    AV_F1,
    AV_F2,
    AV_Ta,
    AV_Tp,
    AV_T,
    AV_k1,
    AV_k2,

    // External inputs
    AV_Cai,
    AV_lambda,
    AV_lambda_rate,

    NUM_ALGEBRAIC
};

enum CONSTANTS_INDEX
{
    // calcium-troponin
    AC_TRPN_n = 0,
    AC_k_trpn_on,
    AC_k_trpn_off,
    AC_ca50_ref,

    // tropomyosin and crossbridge cycle
    AC_perm50,
    AC_nperm,
    AC_ktm_unblock,
    AC_dr,
    AC_wfrac,
    AC_TOT_A,
    AC_gamma,
    AC_gamma_wu,
    AC_phi,
    AC_Tref,
    AC_k_uw,
    AC_k_ws,

    // myosin OFF states
    AC_k1_0,
    AC_k2_0,
    AC_kforce,

    // derived active-model constants
    AC_k_wu,
    AC_k_su,
    AC_cdw,
    AC_cds,
    AC_A,
    AC_ktm_block,

    // passive model
    AC_par_k,
    AC_b,
    AC_eta_l,
    AC_eta_s,
    AC_passive_a,

    NUM_CONSTANTS
};

inline void LewalleNiederer2024initConsts
(
    double* CONSTANTS,
    double* RATES,
    double* STATES
);

#ifndef CARDIAC_HOST_DEVICE
#define CARDIAC_HOST_DEVICE
#endif

CARDIAC_HOST_DEVICE inline void LewalleNiederer2024computeVariables
(
    double VOI,
    double* CONSTANTS,
    double* RATES,
    double* STATES,
    double* ALGEBRAIC
);

#endif
