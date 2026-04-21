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

#pragma once
#include "stimulusIO.H"
#include "dictionary.H"

enum CONSTANT_INDEX {
    AC_CMDN_max,
    AC_CSQN_max,
    AC_Ca_up_max,
    AC_Cao,
    AC_Cm,
    AC_ECaL,
    AC_F,
    AC_FRT,
    AC_INaCa_max,
    AC_INaK_max,
    AC_I_diff,
    AC_I_up_max,
    AC_KQ10,
    AC_krel,          // RENAMED from AC_K_rel to match equations
    AC_K_up,
    AC_KmCa,
    AC_KmKo,
    AC_KmNa,
    AC_KmNai,
    AC_Km_CMDN,
    AC_Km_CSQN,
    AC_Km_TRPN,
    AC_Ko,
    AC_Nao,
    AC_Clo,           // NEW: Extracellular chloride
    AC_R,
    AC_RTF,
    AC_T,
    AC_TRPN_max,
    AC_V_cell,
    AC_V_i,
    AC_V_rel,
    AC_V_up,
    AC_c1,
    AC_c2,
    AC_cajsr_u_tau,
    AC_gNa,
    AC_gK1,
    AC_gKr,
    AC_gKs,
    AC_gKur_amp,      // RENAMED from gKur_base
    AC_gClCa,         // NEW: gClCa (gto,2)
    AC_gbCa,
    AC_gbNa,
    AC_gbK,           // NEW: Background K conductance if needed
    AC_ical_fCa_tau,
    AC_gamma,         // RENAMED/Standardized from AC_g
    AC_ksat,
    AC_sigma,
    AC_tau_tr,
    AC_gCaL,          // Moved to end to preserve enum order
    AC_IpCa_max,      // Moved to end to preserve enum order

    NUM_CONSTANTS};

enum STATE_INDEX {
    membrane_V,
    sodium_Nai,
    potassium_Ki,
    chloride_Cli,     // NEW
    calcium_Cai,
    calcium_CaUp,
    calcium_CaRel,
    ina_m,
    ina_h,
    ina_j,
    ikur_ua,
    ikur_uif,         // NEW: Fast inactivation
    ikur_uis,         // NEW: Slow inactivation
    ikr_xr,
    iks_xs,
    ical_d,
    ical_f,
    ical_fCa,
    iclca_qCa,        // NEW: Cl(Ca) gating
    ryr_u,            // RENAMED from cajsr_u
    ryr_v,            // RENAMED from cajsr_v
    ryr_w,            // RENAMED from cajsr_w

    NUM_STATES
};

enum ALGEBRAIC_INDEX
{
    // --- Currents ---
    AV_ICaL,
    AV_IK1,
    AV_IKr,
    AV_IKs,
    AV_IKur,
    AV_INa,
    AV_INaCa,
    AV_INaK,
    AV_IClCa,      // NEW
    AV_IbCa,
    AV_IbNa,
    AV_IbK,        // NEW (if used)
    AV_IpCa,

    // --- INa gating ---
    AV_ina_m_inf,
    AV_ina_m_tau,
    AV_ina_h_inf,
    AV_ina_h_tau,
    AV_ina_j_inf,
    AV_ina_j_tau,

    // --- ICaL gating ---
    AV_ical_d_inf,
    AV_ical_d_tau,
    AV_ical_f_inf,
    AV_ical_f_tau,

    // --- IKur gating ---
    AV_ikur_ua_inf,
    AV_ikur_ua_tau,
    AV_ikur_uif_inf, // NEW
    AV_ikur_uif_tau, // NEW
    AV_ikur_uis_inf, // NEW
    AV_ikur_uis_tau, // NEW
    AV_gKur,

    // --- IKr gating ---
    AV_ikr_xr_inf,
    AV_ikr_xr_tau,

    // --- IKs gating ---
    AV_iks_xs_inf,
    AV_iks_xs_tau,

    // --- IClCa gating ---
    AV_iclca_qCa_inf, // NEW
    AV_iclca_qCa_tau, // NEW

    // --- Other V-dependent ---
    AV_fNaK,

    // --- RyR gating ---
    AV_ryr_w_inf,
    AV_ryr_w_tau,

    // --- Totals / stimulus ---
    Iion_cm,
    Istim,

    NUM_ALGEBRAIC
};

void
PerisYague_2022initConsts(double* CONSTANTS, double* RATES, double *STATES, int tissueFlag, const Foam::dictionary& stimulus);
void
PerisYague_2022computeVariables(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC,int tissueFlag, bool solveVmWithinODESolver, const Foam::StimulusProtocol& stimulus);
