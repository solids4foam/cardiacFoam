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

#ifndef ORd_2011Optimized_H
#define ORd_2011Optimized_H

#include "ORd_2011.H"

inline double ordSqr(const double x)
{
    return x*x;
}

inline double ordCube(const double x)
{
    return x*x*x;
}

inline double ordPow4(const double x)
{
    const double x2 = x*x;
    return x2*x2;
}

inline double ordPow8(const double x)
{
    const double x2 = x*x;
    const double x4 = x2*x2;
    return x4*x4;
}

inline void
ORdcomputeVariablesOptimized
(
    double VOI,
    double* CONSTANTS,
    double* RATES,
    double* STATES,
    double* ALGEBRAIC,
    int tissueFlag,
    bool solveVmWithinODESolver,
    const Foam::StimulusProtocol& stimulus
)
{
    const double Vm = STATES[membrane_V];
    const double cai = STATES[Ca_i];
    const double cass = STATES[Ca_ss];
    const double cajsr = STATES[Ca_jsr];
    const double cansr = STATES[Ca_nsr];
    const double nai = STATES[Na_i];
    const double nass = STATES[Na_ss];
    const double ki = STATES[K_i];
    const double kss = STATES[K_ss];

    const double F = CONSTANTS[AC_F];
    const double R = CONSTANTS[AC_R];
    const double T = CONSTANTS[AC_T];
    const double frt = F/(R*T);
    const double vfrt = Vm*frt;
    const double vffrt = Vm*F*frt;
    const double expVfrt = exp(vfrt);
    const double exp2Vfrt = exp(2.0*vfrt);

/* CaMK */
    ALGEBRAIC[AV_CaMKb] = CONSTANTS[AC_CaMKo] * (1.0 - STATES[CaMK_CaMKt]) / (1.0 + CONSTANTS[AC_KmCaM] / cass);
    ALGEBRAIC[AV_CaMKa] = ALGEBRAIC[AV_CaMKb] + STATES[CaMK_CaMKt];
    RATES[CaMK_CaMKt] = CONSTANTS[AC_aCaMK] * ALGEBRAIC[AV_CaMKb] * ALGEBRAIC[AV_CaMKa] - CONSTANTS[AC_bCaMK] * STATES[CaMK_CaMKt];

    /* IpCa */
    ALGEBRAIC[AV_IpCa_IpCa] = CONSTANTS[AC_GpCa] * cai / (CONSTANTS[AC_KmCap] + cai);

    /* diff */
    ALGEBRAIC[AV_Jdiff] = (cass - cai) / 0.2;
    ALGEBRAIC[AV_JdiffK] = (kss - ki) / 2.0;
    ALGEBRAIC[AV_JdiffNa] = (nass - nai) / 2.0;

    /* trans_flux */
    ALGEBRAIC[AV_Jtr] = (cansr - cajsr) / 100.0;

    /* INaCa_i */
    const double kmCaActOverCai = CONSTANTS[AC_KmCaAct] / cai;
    const double kmCaActOverCass = CONSTANTS[AC_KmCaAct] / cass;
    ALGEBRAIC[AV_allo_i] = 1.0 / (1.0 + ordSqr(kmCaActOverCai));
    ALGEBRAIC[AV_allo_ss] = 1.0 / (1.0 + ordSqr(kmCaActOverCass));
    ALGEBRAIC[AV_h4_i] = 1.0 + nai / CONSTANTS[AC_kna1] * (1.0 + nai / CONSTANTS[AC_kna2]);
    ALGEBRAIC[AV_h4_ss] = 1.0 + nass / CONSTANTS[AC_kna1] * (1.0 + nass / CONSTANTS[AC_kna2]);
    const double expQcaVfrt = exp(CONSTANTS[AC_qca] * vfrt);
    const double expQnaVfrt = exp(CONSTANTS[AC_qna] * vfrt);
    ALGEBRAIC[AV_hca] = expQcaVfrt;
    ALGEBRAIC[AV_hna] = expQnaVfrt;
    ALGEBRAIC[AV_h1_i] = 1.0 + nai / CONSTANTS[AC_kna3] * (1.0 + ALGEBRAIC[AV_hna]);
    ALGEBRAIC[AV_h1_ss] = 1.0 + nass / CONSTANTS[AC_kna3] * (1.0 + ALGEBRAIC[AV_hna]);
    ALGEBRAIC[AV_h5_i] = ordSqr(nai) / (ALGEBRAIC[AV_h4_i] * CONSTANTS[AC_kna1] * CONSTANTS[AC_kna2]);
    ALGEBRAIC[AV_h5_ss] = ordSqr(nass) / (ALGEBRAIC[AV_h4_ss] * CONSTANTS[AC_kna1] * CONSTANTS[AC_kna2]);
    ALGEBRAIC[AV_h6_i] = 1.0 / ALGEBRAIC[AV_h4_i];
    ALGEBRAIC[AV_h6_ss] = 1.0 / ALGEBRAIC[AV_h4_ss];
    const double invHna = 1.0 / ALGEBRAIC[AV_hna];
    ALGEBRAIC[AV_h7_i] = 1.0 + CONSTANTS[AC_nao] / CONSTANTS[AC_kna3] * (1.0 + invHna);
    ALGEBRAIC[AV_h7_ss] = ALGEBRAIC[AV_h7_i];
    ALGEBRAIC[AV_h2_i] = nai * ALGEBRAIC[AV_hna] / (CONSTANTS[AC_kna3] * ALGEBRAIC[AV_h1_i]);
    ALGEBRAIC[AV_h2_ss] = nass * ALGEBRAIC[AV_hna] / (CONSTANTS[AC_kna3] * ALGEBRAIC[AV_h1_ss]);
    ALGEBRAIC[AV_h3_i] = 1.0 / ALGEBRAIC[AV_h1_i];
    ALGEBRAIC[AV_h3_ss] = 1.0 / ALGEBRAIC[AV_h1_ss];
    ALGEBRAIC[AV_h8_i] = CONSTANTS[AC_nao] / (CONSTANTS[AC_kna3] * ALGEBRAIC[AV_hna] * ALGEBRAIC[AV_h7_i]);
    ALGEBRAIC[AV_h8_ss] = CONSTANTS[AC_nao] / (CONSTANTS[AC_kna3] * ALGEBRAIC[AV_hna] * ALGEBRAIC[AV_h7_ss]);
    ALGEBRAIC[AV_h9_i] = 1.0 / ALGEBRAIC[AV_h7_i];
    ALGEBRAIC[AV_h9_ss] = 1.0 / ALGEBRAIC[AV_h7_ss];
    ALGEBRAIC[AV_k6_i] = ALGEBRAIC[AV_h6_i] * cai * CONSTANTS[AC_kcaon];
    ALGEBRAIC[AV_k6_ss] = ALGEBRAIC[AV_h6_ss] * cass * CONSTANTS[AC_kcaon];
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
    const double inacaDenI = ALGEBRAIC[AV_x1_i] + ALGEBRAIC[AV_x2_i] + ALGEBRAIC[AV_x3_i] + ALGEBRAIC[AV_x4_i];
    const double inacaDenSs = ALGEBRAIC[AV_x1_ss] + ALGEBRAIC[AV_x2_ss] + ALGEBRAIC[AV_x3_ss] + ALGEBRAIC[AV_x4_ss];
    ALGEBRAIC[AV_E1_i] = ALGEBRAIC[AV_x1_i] / inacaDenI;
    ALGEBRAIC[AV_E1_ss] = ALGEBRAIC[AV_x1_ss] / inacaDenSs;
    ALGEBRAIC[AV_E2_i] = ALGEBRAIC[AV_x2_i] / inacaDenI;
    ALGEBRAIC[AV_E2_ss] = ALGEBRAIC[AV_x2_ss] / inacaDenSs;
    ALGEBRAIC[AV_E3_i] = ALGEBRAIC[AV_x3_i] / inacaDenI;
    ALGEBRAIC[AV_E3_ss] = ALGEBRAIC[AV_x3_ss] / inacaDenSs;
    ALGEBRAIC[AV_E4_i] = ALGEBRAIC[AV_x4_i] / inacaDenI;
    ALGEBRAIC[AV_E4_ss] = ALGEBRAIC[AV_x4_ss] / inacaDenSs;
    ALGEBRAIC[AV_JncxCa_i] = ALGEBRAIC[AV_E2_i] * CONSTANTS[AC_k2_i] - ALGEBRAIC[AV_E1_i] * CONSTANTS[AC_k1_i];
    ALGEBRAIC[AV_JncxCa_ss] = ALGEBRAIC[AV_E2_ss] * CONSTANTS[AC_k2_ss] - ALGEBRAIC[AV_E1_ss] * CONSTANTS[AC_k1_ss];
    ALGEBRAIC[AV_JncxNa_i] = 3.0 * (ALGEBRAIC[AV_E4_i] * ALGEBRAIC[AV_k7_i] - ALGEBRAIC[AV_E1_i] * ALGEBRAIC[AV_k8_i]) + ALGEBRAIC[AV_E3_i] * ALGEBRAIC[AV_k4pp_i] - ALGEBRAIC[AV_E2_i] * ALGEBRAIC[AV_k3pp_i];
    ALGEBRAIC[AV_JncxNa_ss] = 3.0 * (ALGEBRAIC[AV_E4_ss] * ALGEBRAIC[AV_k7_ss] - ALGEBRAIC[AV_E1_ss] * ALGEBRAIC[AV_k8_ss]) + ALGEBRAIC[AV_E3_ss] * ALGEBRAIC[AV_k4pp_ss] - ALGEBRAIC[AV_E2_ss] * ALGEBRAIC[AV_k3pp_ss];
    ALGEBRAIC[AV_INaCa_i_INaCa_i] = 0.8 * CONSTANTS[AC_Gncx] * ALGEBRAIC[AV_allo_i] * (CONSTANTS[AC_zna] * ALGEBRAIC[AV_JncxNa_i] + CONSTANTS[AC_zca] * ALGEBRAIC[AV_JncxCa_i]);
    ALGEBRAIC[AV_INaCa_ss] = 0.2 * CONSTANTS[AC_Gncx] * ALGEBRAIC[AV_allo_ss] * (CONSTANTS[AC_zna] * ALGEBRAIC[AV_JncxNa_ss] + CONSTANTS[AC_zca] * ALGEBRAIC[AV_JncxCa_ss]);

    /* INaK */
    const double deltaVfrtOver3 = CONSTANTS[AC_delta] * vfrt / 3.0;
    ALGEBRAIC[AV_Knai] = CONSTANTS[AC_Knai0] * exp(deltaVfrtOver3);
    ALGEBRAIC[AV_Knao] = CONSTANTS[AC_Knao0] * exp((1.0 - CONSTANTS[AC_delta]) * vfrt / 3.0);
    ALGEBRAIC[AV_P] = CONSTANTS[AC_eP] / (1.0 + CONSTANTS[AC_H] / CONSTANTS[AC_Khp] + nai / CONSTANTS[AC_Knap] + ki / CONSTANTS[AC_Kxkur]);
    const double naiOverKnai = nai / ALGEBRAIC[AV_Knai];
    const double naoOverKnao = CONSTANTS[AC_nao] / ALGEBRAIC[AV_Knao];
    const double kiOverKki = ki / CONSTANTS[AC_Kki];
    const double koOverKko = CONSTANTS[AC_ko] / CONSTANTS[AC_Kko];
    const double nakDenInner = ordCube(1.0 + naiOverKnai) + ordSqr(1.0 + kiOverKki) - 1.0;
    const double nakDenOuter = ordCube(1.0 + naoOverKnao) + ordSqr(1.0 + koOverKko) - 1.0;
    ALGEBRAIC[AV_a1] = CONSTANTS[AC_k1p] * ordCube(naiOverKnai) / nakDenInner;
    ALGEBRAIC[AV_a3] = CONSTANTS[AC_k3p] * ordSqr(koOverKko) / nakDenOuter;
    ALGEBRAIC[AV_b2] = CONSTANTS[AC_k2m] * ordCube(naoOverKnao) / nakDenOuter;
    ALGEBRAIC[AV_b3] = CONSTANTS[AC_k3m] * ALGEBRAIC[AV_P] * CONSTANTS[AC_H] / (1.0 + CONSTANTS[AC_MgATP] / CONSTANTS[AC_Kmgatp]);
    ALGEBRAIC[AV_b4] = CONSTANTS[AC_k4m] * ordSqr(kiOverKki) / nakDenInner;
    ALGEBRAIC[AV_x1] = CONSTANTS[AC_a4] * ALGEBRAIC[AV_a1] * CONSTANTS[AC_a2] + ALGEBRAIC[AV_b2] * ALGEBRAIC[AV_b4] * ALGEBRAIC[AV_b3] + CONSTANTS[AC_a2] * ALGEBRAIC[AV_b4] * ALGEBRAIC[AV_b3] + ALGEBRAIC[AV_b3] * ALGEBRAIC[AV_a1] * CONSTANTS[AC_a2];
    ALGEBRAIC[AV_x2] = ALGEBRAIC[AV_b2] * CONSTANTS[AC_b1] * ALGEBRAIC[AV_b4] + ALGEBRAIC[AV_a1] * CONSTANTS[AC_a2] * ALGEBRAIC[AV_a3] + ALGEBRAIC[AV_a3] * CONSTANTS[AC_b1] * ALGEBRAIC[AV_b4] + CONSTANTS[AC_a2] * ALGEBRAIC[AV_a3] * ALGEBRAIC[AV_b4];
    ALGEBRAIC[AV_x3] = CONSTANTS[AC_a2] * ALGEBRAIC[AV_a3] * CONSTANTS[AC_a4] + ALGEBRAIC[AV_b3] * ALGEBRAIC[AV_b2] * CONSTANTS[AC_b1] + ALGEBRAIC[AV_b2] * CONSTANTS[AC_b1] * CONSTANTS[AC_a4] + ALGEBRAIC[AV_a3] * CONSTANTS[AC_a4] * CONSTANTS[AC_b1];
    ALGEBRAIC[AV_x4] = ALGEBRAIC[AV_b4] * ALGEBRAIC[AV_b3] * ALGEBRAIC[AV_b2] + ALGEBRAIC[AV_a3] * CONSTANTS[AC_a4] * ALGEBRAIC[AV_a1] + ALGEBRAIC[AV_b2] * CONSTANTS[AC_a4] * ALGEBRAIC[AV_a1] + ALGEBRAIC[AV_b3] * ALGEBRAIC[AV_b2] * ALGEBRAIC[AV_a1];
    const double inakDen = ALGEBRAIC[AV_x1] + ALGEBRAIC[AV_x2] + ALGEBRAIC[AV_x3] + ALGEBRAIC[AV_x4];
    ALGEBRAIC[AV_E1] = ALGEBRAIC[AV_x1] / inakDen;
    ALGEBRAIC[AV_E2] = ALGEBRAIC[AV_x2] / inakDen;
    ALGEBRAIC[AV_E3] = ALGEBRAIC[AV_x3] / inakDen;
    ALGEBRAIC[AV_E4] = ALGEBRAIC[AV_x4] / inakDen;
    ALGEBRAIC[AV_JnakK] = 2.0 * (ALGEBRAIC[AV_E4] * CONSTANTS[AC_b1] - ALGEBRAIC[AV_E3] * ALGEBRAIC[AV_a1]);
    ALGEBRAIC[AV_JnakNa] = 3.0 * (ALGEBRAIC[AV_E1] * ALGEBRAIC[AV_a3] - ALGEBRAIC[AV_E2] * ALGEBRAIC[AV_b3]);
    ALGEBRAIC[AV_INaK_INaK] = CONSTANTS[AC_Pnak] * (CONSTANTS[AC_zna] * ALGEBRAIC[AV_JnakNa] + CONSTANTS[AC_zk] * ALGEBRAIC[AV_JnakK]);

    /* SERCA */
    ALGEBRAIC[AV_Jleak] = 0.0039375 * cansr / 15.0;
    ALGEBRAIC[AV_fJupp] = 1.0 / (1.0 + CONSTANTS[AC_KmCaMK] / ALGEBRAIC[AV_CaMKa]);
    const double oneMinusFCaMK = 1.0 - ALGEBRAIC[AV_fJupp];
    ALGEBRAIC[AV_Jupnp] = CONSTANTS[AC_upScale] * 0.004375 * cai / (cai + 0.00092);
    ALGEBRAIC[AV_Jupp] = CONSTANTS[AC_upScale] * 2.75 * 0.004375 * cai / (cai + 0.00092 - 0.00017);
    ALGEBRAIC[AV_Jup] = oneMinusFCaMK * ALGEBRAIC[AV_Jupnp] + ALGEBRAIC[AV_fJupp] * ALGEBRAIC[AV_Jupp] - ALGEBRAIC[AV_Jleak];

    /* reversal_potentials */
    const double RTonF = R * T / F;
    ALGEBRAIC[AV_EK] = RTonF * log(CONSTANTS[AC_ko] / ki);
    ALGEBRAIC[AV_ENa] = RTonF * log(CONSTANTS[AC_nao] / nai);
    ALGEBRAIC[AV_EKs] = RTonF * log((CONSTANTS[AC_ko] + CONSTANTS[AC_PKNa] * CONSTANTS[AC_nao]) / (ki + CONSTANTS[AC_PKNa] * nai));

    /* IK1 */
    ALGEBRAIC[AV_rk1] = 1.0 / (1.0 + exp((Vm + 105.8 - 2.6 * CONSTANTS[AC_ko]) / 9.493));
    ALGEBRAIC[AV_txk1] = 122.2 / (exp(-(Vm + 127.2) / 20.36) + exp((Vm + 236.8) / 69.33));
    ALGEBRAIC[AV_xk1ss] = 1.0 / (1.0 + exp(-(Vm + 2.5538 * CONSTANTS[AC_ko] + 144.59) / (1.5692 * CONSTANTS[AC_ko] + 3.8115)));
    RATES[IK1_xk1] = (ALGEBRAIC[AV_xk1ss] - STATES[IK1_xk1]) / ALGEBRAIC[AV_txk1];
    ALGEBRAIC[AV_IK1_IK1] = CONSTANTS[AC_GK1] * sqrt(CONSTANTS[AC_ko]) * ALGEBRAIC[AV_rk1] * STATES[IK1_xk1] * (Vm - ALGEBRAIC[AV_EK]);

    /* IKb */
    ALGEBRAIC[AV_xkb] = 1.0 / (1.0 + exp(-(Vm - 14.48) / 18.34));
    ALGEBRAIC[AV_IKb_IKb] = CONSTANTS[AC_GKb] * ALGEBRAIC[AV_xkb] * (Vm - ALGEBRAIC[AV_EK]);

    /* IKr */
    ALGEBRAIC[AV_Axrf] = 1.0 / (1.0 + exp((Vm + 54.81) / 38.21));
    ALGEBRAIC[AV_rkr] = 1.0 / (1.0 + exp((Vm + 55.0) / 75.0)) * 1.0 / (1.0 + exp((Vm - 10.0) / 30.0));
    ALGEBRAIC[AV_txrf] = 12.98 + 1.0 / (0.3652 * exp((Vm - 31.66) / 3.869) + 4.123e-05 * exp(-(Vm - 47.78) / 20.38));
    ALGEBRAIC[AV_txrs] = 1.865 + 1.0 / (0.06629 * exp((Vm - 34.7) / 7.355) + 1.128e-05 * exp(-(Vm - 29.74) / 25.94));
    ALGEBRAIC[AV_xrss] = 1.0 / (1.0 + exp(-(Vm + 8.337) / 6.789));
    ALGEBRAIC[AV_Axrs] = 1.0 - ALGEBRAIC[AV_Axrf];
    RATES[IKr_xrf] = (ALGEBRAIC[AV_xrss] - STATES[IKr_xrf]) / ALGEBRAIC[AV_txrf];
    RATES[IKr_xrs] = (ALGEBRAIC[AV_xrss] - STATES[IKr_xrs]) / ALGEBRAIC[AV_txrs];
    ALGEBRAIC[AV_xr] = ALGEBRAIC[AV_Axrf] * STATES[IKr_xrf] + ALGEBRAIC[AV_Axrs] * STATES[IKr_xrs];
    ALGEBRAIC[AV_IKr_IKr] = CONSTANTS[AC_GKr] * sqrt(CONSTANTS[AC_ko] / 5.4) * ALGEBRAIC[AV_xr] * ALGEBRAIC[AV_rkr] * (Vm - ALGEBRAIC[AV_EK]);

    /* IKs */
    ALGEBRAIC[AV_KsCa] = 1.0 + 0.6 / (1.0 + pow(3.8e-05 / cai, 1.4));
    ALGEBRAIC[AV_txs1] = 817.3 + 1.0 / (0.0002326 * exp((Vm + 48.28) / 17.8) + 0.001292 * exp(-(Vm + 210.0) / 230.0));
    ALGEBRAIC[AV_txs2] = 1.0 / (0.01 * exp((Vm - 50.0) / 20.0) + 0.0193 * exp(-(Vm + 66.54) / 31.0));
    ALGEBRAIC[AV_xs1ss] = 1.0 / (1.0 + exp(-(Vm + 11.6) / 8.932));
    ALGEBRAIC[AV_xs2ss] = ALGEBRAIC[AV_xs1ss];
    RATES[IKs_xs1] = (ALGEBRAIC[AV_xs1ss] - STATES[IKs_xs1]) / ALGEBRAIC[AV_txs1];
    ALGEBRAIC[AV_IKs_IKs] = CONSTANTS[AC_GKs] * ALGEBRAIC[AV_KsCa] * STATES[IKs_xs1] * STATES[IKs_xs2] * (Vm - ALGEBRAIC[AV_EKs]);
    RATES[IKs_xs2] = (ALGEBRAIC[AV_xs2ss] - STATES[IKs_xs2]) / ALGEBRAIC[AV_txs2];

    /* INa */
    ALGEBRAIC[AV_fINap] = ALGEBRAIC[AV_fJupp];
    ALGEBRAIC[AV_hssp] = 1.0 / (1.0 + exp((Vm + 89.1) / 6.086));
    ALGEBRAIC[AV_thf] = 1.0 / (1.432e-05 * exp(-(Vm + 1.196) / 6.285) + 6.149 * exp((Vm + 0.5096) / 20.27));
    ALGEBRAIC[AV_ths] = 1.0 / (0.009794 * exp(-(Vm + 17.95) / 28.05) + 0.3343 * exp((Vm + 5.73) / 56.66));
    ALGEBRAIC[AV_tj] = 2.038 + 1.0 / (0.02136 * exp(-(Vm + 100.6) / 8.281) + 0.3052 * exp((Vm + 0.9941) / 38.45));
    ALGEBRAIC[AV_hss] = 1.0 / (1.0 + exp((Vm + CONSTANTS[AC_hssV1]) / CONSTANTS[AC_hssV2]));
    ALGEBRAIC[AV_mss] = 1.0 / (1.0 + exp(-(Vm + CONSTANTS[AC_mssV1]) / CONSTANTS[AC_mssV2]));
    ALGEBRAIC[AV_thsp] = 3.0 * ALGEBRAIC[AV_ths];
    ALGEBRAIC[AV_tjp] = 1.46 * ALGEBRAIC[AV_tj];
    ALGEBRAIC[AV_tm] = 1.0 / (CONSTANTS[AC_mtD1] * exp((Vm + CONSTANTS[AC_mtV1]) / CONSTANTS[AC_mtV2]) + CONSTANTS[AC_mtD2] * exp(-(Vm + CONSTANTS[AC_mtV3]) / CONSTANTS[AC_mtV4]));
    ALGEBRAIC[AV_h] = CONSTANTS[AC_Ahf] * STATES[INa_hf] + CONSTANTS[AC_Ahs] * STATES[INa_hs];
    ALGEBRAIC[AV_hp] = CONSTANTS[AC_Ahf] * STATES[INa_hf] + CONSTANTS[AC_Ahs] * STATES[INa_hsp];
    ALGEBRAIC[AV_jss] = ALGEBRAIC[AV_hss];
    RATES[INa_hf] = (ALGEBRAIC[AV_hss] - STATES[INa_hf]) / ALGEBRAIC[AV_thf];
    RATES[INa_hs] = (ALGEBRAIC[AV_hss] - STATES[INa_hs]) / ALGEBRAIC[AV_ths];
    RATES[INa_hsp] = (ALGEBRAIC[AV_hssp] - STATES[INa_hsp]) / ALGEBRAIC[AV_thsp];
    RATES[INa_m] = (ALGEBRAIC[AV_mss] - STATES[INa_m]) / ALGEBRAIC[AV_tm];
    ALGEBRAIC[AV_INa_INa] = CONSTANTS[AC_GNa] * (Vm - ALGEBRAIC[AV_ENa]) * ordCube(STATES[INa_m]) * (oneMinusFCaMK * ALGEBRAIC[AV_h] * STATES[INa_j] + ALGEBRAIC[AV_fINap] * ALGEBRAIC[AV_hp] * STATES[INa_jp]);
    RATES[INa_j] = (ALGEBRAIC[AV_jss] - STATES[INa_j]) / ALGEBRAIC[AV_tj];
    RATES[INa_jp] = (ALGEBRAIC[AV_jss] - STATES[INa_jp]) / ALGEBRAIC[AV_tjp];

    /* Ito */
    ALGEBRAIC[AV_AiF] = 1.0 / (1.0 + exp((Vm - 213.6) / 151.2));
    ALGEBRAIC[AV_ass] = 1.0 / (1.0 + exp(-(Vm - 14.34) / 14.82));
    ALGEBRAIC[AV_assp] = 1.0 / (1.0 + exp(-(Vm - 24.34) / 14.82));
    ALGEBRAIC[AV_delta_epi] = ((tissueFlag == 1) ? 1.0 - 0.95 / (1.0 + exp((Vm + 70.0) / 5.0)) : 1.0);
    ALGEBRAIC[AV_dti_develop] = 1.354 + 0.0001 / (exp((Vm - 167.4) / 15.89) + exp(-(Vm - 12.23) / 0.2154));
    ALGEBRAIC[AV_dti_recover] = 1.0 - 0.5 / (1.0 + exp((Vm + 70.0) / 20.0));
    ALGEBRAIC[AV_fItop] = ALGEBRAIC[AV_fJupp];
    ALGEBRAIC[AV_iss] = 1.0 / (1.0 + exp((Vm + 43.94) / 5.711));
    ALGEBRAIC[AV_ta] = 1.0515 / (1.0 / (1.2089 * (1.0 + exp(-(Vm - 18.4099) / 29.3814))) + 3.5 / (1.0 + exp((Vm + 100.0) / 29.3814)));
    ALGEBRAIC[AV_tiF_b] = 4.562 + 1.0 / (0.3933 * exp(-(Vm + 100.0) / 100.0) + 0.08004 * exp((Vm + 50.0) / 16.59));
    ALGEBRAIC[AV_tiS_b] = 23.62 + 1.0 / (0.001416 * exp(-(Vm + 96.52) / 59.05) + 1.78e-08 * exp((Vm + 114.1) / 8.079));
    ALGEBRAIC[AV_AiS] = 1.0 - ALGEBRAIC[AV_AiF];
    ALGEBRAIC[AV_tiF] = ALGEBRAIC[AV_tiF_b] * ALGEBRAIC[AV_delta_epi];
    ALGEBRAIC[AV_tiS] = ALGEBRAIC[AV_tiS_b] * ALGEBRAIC[AV_delta_epi];
    RATES[Ito_a] = (ALGEBRAIC[AV_ass] - STATES[Ito_a]) / ALGEBRAIC[AV_ta];
    RATES[Ito_ap] = (ALGEBRAIC[AV_assp] - STATES[Ito_ap]) / ALGEBRAIC[AV_ta];
    ALGEBRAIC[AV_i] = ALGEBRAIC[AV_AiF] * STATES[Ito_iF] + ALGEBRAIC[AV_AiS] * STATES[Ito_iS];
    ALGEBRAIC[AV_ip] = ALGEBRAIC[AV_AiF] * STATES[Ito_iFp] + ALGEBRAIC[AV_AiS] * STATES[Ito_iSp];
    ALGEBRAIC[AV_tiFp] = ALGEBRAIC[AV_dti_develop] * ALGEBRAIC[AV_dti_recover] * ALGEBRAIC[AV_tiF];
    ALGEBRAIC[AV_tiSp] = ALGEBRAIC[AV_dti_develop] * ALGEBRAIC[AV_dti_recover] * ALGEBRAIC[AV_tiS];
    RATES[Ito_iF] = (ALGEBRAIC[AV_iss] - STATES[Ito_iF]) / ALGEBRAIC[AV_tiF];
    RATES[Ito_iS] = (ALGEBRAIC[AV_iss] - STATES[Ito_iS]) / ALGEBRAIC[AV_tiS];
    ALGEBRAIC[AV_Ito_Ito] = CONSTANTS[AC_Gto] * (Vm - ALGEBRAIC[AV_EK]) * (oneMinusFCaMK * STATES[Ito_a] * ALGEBRAIC[AV_i] + ALGEBRAIC[AV_fItop] * STATES[Ito_ap] * ALGEBRAIC[AV_ip]);
    RATES[Ito_iFp] = (ALGEBRAIC[AV_iss] - STATES[Ito_iFp]) / ALGEBRAIC[AV_tiFp];
    RATES[Ito_iSp] = (ALGEBRAIC[AV_iss] - STATES[Ito_iSp]) / ALGEBRAIC[AV_tiSp];

    /* INaL */
    ALGEBRAIC[AV_fINaLp] = ALGEBRAIC[AV_fJupp];
    ALGEBRAIC[AV_hLss] = 1.0 / (1.0 + exp((Vm + 87.61) / 7.488));
    ALGEBRAIC[AV_hLssp] = 1.0 / (1.0 + exp((Vm + 93.81) / 7.488));
    ALGEBRAIC[AV_mLss] = 1.0 / (1.0 + exp(-(Vm + 42.85) / 5.264));
    ALGEBRAIC[AV_tmL] = ALGEBRAIC[AV_tm];
    RATES[INaL_hL] = (ALGEBRAIC[AV_hLss] - STATES[INaL_hL]) / CONSTANTS[AC_thL];
    RATES[INaL_mL] = (ALGEBRAIC[AV_mLss] - STATES[INaL_mL]) / ALGEBRAIC[AV_tmL];
    ALGEBRAIC[AV_INaL_INaL] = CONSTANTS[AC_GNaL] * (Vm - ALGEBRAIC[AV_ENa]) * STATES[INaL_mL] * (oneMinusFCaMK * STATES[INaL_hL] + ALGEBRAIC[AV_fINaLp] * STATES[INaL_hLp]);
    RATES[INaL_hLp] = (ALGEBRAIC[AV_hLssp] - STATES[INaL_hLp]) / CONSTANTS[AC_thLp];

    /* ICaL */
    ALGEBRAIC[AV_Afcaf] = 0.3 + 0.6 / (1.0 + exp((Vm - 10.0) / 10.0));
    ALGEBRAIC[AV_dss] = 1.0 / (1.0 + exp(-(Vm + 3.94) / 4.23));
    ALGEBRAIC[AV_fICaLp] = ALGEBRAIC[AV_fJupp];
    ALGEBRAIC[AV_fss] = 1.0 / (1.0 + exp((Vm + 19.58) / 3.696));
    ALGEBRAIC[AV_km2n] = STATES[ICaL_jca];
    ALGEBRAIC[AV_td] = 0.6 + 1.0 / (exp(-0.05 * (Vm + 6.0)) + exp(0.09 * (Vm + 14.0)));
    ALGEBRAIC[AV_tfcaf] = 7.0 + 1.0 / (0.04 * exp(-(Vm - 4.0) / 7.0) + 0.04 * exp((Vm - 4.0) / 7.0));
    ALGEBRAIC[AV_tfcas] = 100.0 + 1.0 / (0.00012 * exp(-Vm / 3.0) + 0.00012 * exp(Vm / 7.0));
    ALGEBRAIC[AV_tff] = 7.0 + 1.0 / (0.0045 * exp(-(Vm + 20.0) / 10.0) + 0.0045 * exp((Vm + 20.0) / 10.0));
    ALGEBRAIC[AV_tfs] = 1000.0 + 1.0 / (3.5e-05 * exp(-(Vm + 5.0) / 4.0) + 3.5e-05 * exp((Vm + 5.0) / 6.0));
    ALGEBRAIC[AV_Afcas] = 1.0 - ALGEBRAIC[AV_Afcaf];
    ALGEBRAIC[AV_anca] = 1.0 / (CONSTANTS[AC_k2n] / ALGEBRAIC[AV_km2n] + ordPow4(1.0 + CONSTANTS[AC_Kmn] / cass));
    ALGEBRAIC[AV_fcass] = ALGEBRAIC[AV_fss];
    ALGEBRAIC[AV_tfcafp] = 2.5 * ALGEBRAIC[AV_tfcaf];
    ALGEBRAIC[AV_tffp] = 2.5 * ALGEBRAIC[AV_tff];
    RATES[ICaL_d] = (ALGEBRAIC[AV_dss] - STATES[ICaL_d]) / ALGEBRAIC[AV_td];
    RATES[ICaL_ff] = (ALGEBRAIC[AV_fss] - STATES[ICaL_ff]) / ALGEBRAIC[AV_tff];
    RATES[ICaL_fs] = (ALGEBRAIC[AV_fss] - STATES[ICaL_fs]) / ALGEBRAIC[AV_tfs];
    ALGEBRAIC[AV_f] = CONSTANTS[AC_Aff] * STATES[ICaL_ff] + CONSTANTS[AC_Afs] * STATES[ICaL_fs];
    ALGEBRAIC[AV_fca] = ALGEBRAIC[AV_Afcaf] * STATES[ICaL_fcaf] + ALGEBRAIC[AV_Afcas] * STATES[ICaL_fcas];
    ALGEBRAIC[AV_fcap] = ALGEBRAIC[AV_Afcaf] * STATES[ICaL_fcafp] + ALGEBRAIC[AV_Afcas] * STATES[ICaL_fcas];
    ALGEBRAIC[AV_fp] = CONSTANTS[AC_Aff] * STATES[ICaL_ffp] + CONSTANTS[AC_Afs] * STATES[ICaL_fs];
    RATES[ICaL_fcaf] = (ALGEBRAIC[AV_fcass] - STATES[ICaL_fcaf]) / ALGEBRAIC[AV_tfcaf];
    RATES[ICaL_fcafp] = (ALGEBRAIC[AV_fcass] - STATES[ICaL_fcafp]) / ALGEBRAIC[AV_tfcafp];
    RATES[ICaL_fcas] = (ALGEBRAIC[AV_fcass] - STATES[ICaL_fcas]) / ALGEBRAIC[AV_tfcas];
    RATES[ICaL_ffp] = (ALGEBRAIC[AV_fss] - STATES[ICaL_ffp]) / ALGEBRAIC[AV_tffp];
    RATES[ICaL_jca] = (ALGEBRAIC[AV_fcass] - STATES[ICaL_jca]) / CONSTANTS[AC_tjca];
    RATES[ICaL_nca] = ALGEBRAIC[AV_anca] * CONSTANTS[AC_k2n] - STATES[ICaL_nca] * ALGEBRAIC[AV_km2n];

    /* intracellular_ions */
    RATES[Ca_nsr] = ALGEBRAIC[AV_Jup] - ALGEBRAIC[AV_Jtr] * CONSTANTS[AC_vjsr] / CONSTANTS[AC_vnsr];
    ALGEBRAIC[AV_Bcajsr] = 1.0 / (1.0 + CONSTANTS[AC_csqnmax] * CONSTANTS[AC_kmcsqn] / ordSqr(CONSTANTS[AC_kmcsqn] + cajsr));
    ALGEBRAIC[AV_Bcass] = 1.0 / (1.0 + CONSTANTS[AC_BSRmax] * CONSTANTS[AC_KmBSR] / ordSqr(CONSTANTS[AC_KmBSR] + cass) + CONSTANTS[AC_BSLmax] * CONSTANTS[AC_KmBSL] / ordSqr(CONSTANTS[AC_KmBSL] + cass));
    ALGEBRAIC[AV_Bcai] = 1.0 / (1.0 + CONSTANTS[AC_cmdnmax] * CONSTANTS[AC_kmcmdn] / ordSqr(CONSTANTS[AC_kmcmdn] + cai) + CONSTANTS[AC_trpnmax] * CONSTANTS[AC_kmtrpn] / ordSqr(CONSTANTS[AC_kmtrpn] + cai));

    /* membrane */
    ALGEBRAIC[AV_vffrt] = vffrt;
    ALGEBRAIC[AV_vfrt] = vfrt;

    /* ryr */
    ALGEBRAIC[AV_fJrelp] = ALGEBRAIC[AV_fJupp];
    ALGEBRAIC[AV_Jrel] = oneMinusFCaMK * STATES[ryr_Jrelnp] + ALGEBRAIC[AV_fJrelp] * STATES[ryr_Jrelp];
    ALGEBRAIC[AV_tau_rel_temp] = CONSTANTS[AC_bt] / (1.0 + 0.0123 / cajsr);
    ALGEBRAIC[AV_tau_rel] = ((ALGEBRAIC[AV_tau_rel_temp] < 0.001) ? 0.001 : ALGEBRAIC[AV_tau_rel_temp]);
    ALGEBRAIC[AV_tau_relp_temp] = CONSTANTS[AC_btp] / (1.0 + 0.0123 / cajsr);
    ALGEBRAIC[AV_tau_relp] = ((ALGEBRAIC[AV_tau_relp_temp] < 0.001) ? 0.001 : ALGEBRAIC[AV_tau_relp_temp]);

    /* *remaining* */
    const double expVfrtMinus1 = expVfrt - 1.0;
    const double exp2VfrtMinus1 = exp2Vfrt - 1.0;
    ALGEBRAIC[AV_PhiCaK] = vffrt * (0.75 * kss * expVfrt - 0.75 * CONSTANTS[AC_ko]) / expVfrtMinus1;
    ALGEBRAIC[AV_PhiCaL] = 4.0 * vffrt * (cass * exp2Vfrt - 0.341 * CONSTANTS[AC_cao]) / exp2VfrtMinus1;
    ALGEBRAIC[AV_PhiCaNa] = vffrt * (0.75 * nass * expVfrt - 0.75 * CONSTANTS[AC_nao]) / expVfrtMinus1;
    ALGEBRAIC[AV_ICab_ICab] = CONSTANTS[AC_PCab] * 4.0 * vffrt * (cai * exp2Vfrt - 0.341 * CONSTANTS[AC_cao]) / exp2VfrtMinus1;
    ALGEBRAIC[AV_INab_INab] = CONSTANTS[AC_PNab] * vffrt * (nai * expVfrt - CONSTANTS[AC_nao]) / expVfrtMinus1;
    RATES[Ca_jsr] = ALGEBRAIC[AV_Bcajsr] * (ALGEBRAIC[AV_Jtr] - ALGEBRAIC[AV_Jrel]);
    RATES[K_i] = -(ALGEBRAIC[AV_Ito_Ito] + ALGEBRAIC[AV_IKr_IKr] + ALGEBRAIC[AV_IKs_IKs] + ALGEBRAIC[AV_IK1_IK1] + ALGEBRAIC[AV_IKb_IKb] + ALGEBRAIC[Istim] - 2.0 * ALGEBRAIC[AV_INaK_INaK]) * CONSTANTS[AC_cm] * CONSTANTS[AC_Acap] / (F * CONSTANTS[AC_vmyo]) + ALGEBRAIC[AV_JdiffK] * CONSTANTS[AC_vss] / CONSTANTS[AC_vmyo];
    ALGEBRAIC[AV_ICaK] = oneMinusFCaMK * CONSTANTS[AC_PCaK] * ALGEBRAIC[AV_PhiCaK] * STATES[ICaL_d] * (ALGEBRAIC[AV_f] * (1.0 - STATES[ICaL_nca]) + STATES[ICaL_jca] * ALGEBRAIC[AV_fca] * STATES[ICaL_nca]) + ALGEBRAIC[AV_fICaLp] * CONSTANTS[AC_PCaKp] * ALGEBRAIC[AV_PhiCaK] * STATES[ICaL_d] * (ALGEBRAIC[AV_fp] * (1.0 - STATES[ICaL_nca]) + STATES[ICaL_jca] * ALGEBRAIC[AV_fcap] * STATES[ICaL_nca]);
    ALGEBRAIC[AV_ICaL_ICaL] = oneMinusFCaMK * CONSTANTS[AC_PCa] * ALGEBRAIC[AV_PhiCaL] * STATES[ICaL_d] * (ALGEBRAIC[AV_f] * (1.0 - STATES[ICaL_nca]) + STATES[ICaL_jca] * ALGEBRAIC[AV_fca] * STATES[ICaL_nca]) + ALGEBRAIC[AV_fICaLp] * CONSTANTS[AC_PCap] * ALGEBRAIC[AV_PhiCaL] * STATES[ICaL_d] * (ALGEBRAIC[AV_fp] * (1.0 - STATES[ICaL_nca]) + STATES[ICaL_jca] * ALGEBRAIC[AV_fcap] * STATES[ICaL_nca]);
    ALGEBRAIC[AV_ICaNa] = oneMinusFCaMK * CONSTANTS[AC_PCaNa] * ALGEBRAIC[AV_PhiCaNa] * STATES[ICaL_d] * (ALGEBRAIC[AV_f] * (1.0 - STATES[ICaL_nca]) + STATES[ICaL_jca] * ALGEBRAIC[AV_fca] * STATES[ICaL_nca]) + ALGEBRAIC[AV_fICaLp] * CONSTANTS[AC_PCaNap] * ALGEBRAIC[AV_PhiCaNa] * STATES[ICaL_d] * (ALGEBRAIC[AV_fp] * (1.0 - STATES[ICaL_nca]) + STATES[ICaL_jca] * ALGEBRAIC[AV_fcap] * STATES[ICaL_nca]);
    RATES[Ca_i] = ALGEBRAIC[AV_Bcai] * (-(ALGEBRAIC[AV_IpCa_IpCa] + ALGEBRAIC[AV_ICab_ICab] - 2.0 * ALGEBRAIC[AV_INaCa_i_INaCa_i]) * CONSTANTS[AC_cm] * CONSTANTS[AC_Acap] / (2.0 * F * CONSTANTS[AC_vmyo]) - ALGEBRAIC[AV_Jup] * CONSTANTS[AC_vnsr] / CONSTANTS[AC_vmyo] + ALGEBRAIC[AV_Jdiff] * CONSTANTS[AC_vss] / CONSTANTS[AC_vmyo]);
    RATES[Na_i] = -(ALGEBRAIC[AV_INa_INa] + ALGEBRAIC[AV_INaL_INaL] + 3.0 * ALGEBRAIC[AV_INaCa_i_INaCa_i] + 3.0 * ALGEBRAIC[AV_INaK_INaK] + ALGEBRAIC[AV_INab_INab]) * CONSTANTS[AC_Acap] * CONSTANTS[AC_cm] / (F * CONSTANTS[AC_vmyo]) + ALGEBRAIC[AV_JdiffNa] * CONSTANTS[AC_vss] / CONSTANTS[AC_vmyo];
    RATES[Ca_ss] = ALGEBRAIC[AV_Bcass] * (-(ALGEBRAIC[AV_ICaL_ICaL] - 2.0 * ALGEBRAIC[AV_INaCa_ss]) * CONSTANTS[AC_cm] * CONSTANTS[AC_Acap] / (2.0 * F * CONSTANTS[AC_vss]) + ALGEBRAIC[AV_Jrel] * CONSTANTS[AC_vjsr] / CONSTANTS[AC_vss] - ALGEBRAIC[AV_Jdiff]);
    RATES[K_ss] = -ALGEBRAIC[AV_ICaK] * CONSTANTS[AC_cm] * CONSTANTS[AC_Acap] / (F * CONSTANTS[AC_vss]) - ALGEBRAIC[AV_JdiffK];
    RATES[Na_ss] = -(ALGEBRAIC[AV_ICaNa] + 3.0 * ALGEBRAIC[AV_INaCa_ss]) * CONSTANTS[AC_cm] * CONSTANTS[AC_Acap] / (F * CONSTANTS[AC_vss]) - ALGEBRAIC[AV_JdiffNa];
    const double jrelPow = ordPow8(1.5 / cajsr);
    ALGEBRAIC[AV_Jrel_inf_temp] = CONSTANTS[AC_a_rel] * -ALGEBRAIC[AV_ICaL_ICaL] / (1.0 + jrelPow);
    ALGEBRAIC[AV_Jrel_temp] = CONSTANTS[AC_a_relp] * -ALGEBRAIC[AV_ICaL_ICaL] / (1.0 + jrelPow);
    ALGEBRAIC[AV_Jrel_inf] = ((tissueFlag == 2) ? ALGEBRAIC[AV_Jrel_inf_temp] * 1.7 : ALGEBRAIC[AV_Jrel_inf_temp]);
    ALGEBRAIC[AV_Jrel_infp] = ((tissueFlag == 2) ? ALGEBRAIC[AV_Jrel_temp] * 1.7 : ALGEBRAIC[AV_Jrel_temp]);
    RATES[ryr_Jrelnp] = (ALGEBRAIC[AV_Jrel_inf] - STATES[ryr_Jrelnp]) / ALGEBRAIC[AV_tau_rel];
    RATES[ryr_Jrelp] = (ALGEBRAIC[AV_Jrel_infp] - STATES[ryr_Jrelp]) / ALGEBRAIC[AV_tau_relp];

    ALGEBRAIC[Iion_cm] = ALGEBRAIC[AV_INa_INa] + ALGEBRAIC[AV_INaL_INaL] + ALGEBRAIC[AV_Ito_Ito] + ALGEBRAIC[AV_ICaL_ICaL] + ALGEBRAIC[AV_ICaNa] + ALGEBRAIC[AV_ICaK] + ALGEBRAIC[AV_IKr_IKr] + ALGEBRAIC[AV_IKs_IKs] + ALGEBRAIC[AV_IK1_IK1] + ALGEBRAIC[AV_INaCa_i_INaCa_i] + ALGEBRAIC[AV_INaCa_ss] + ALGEBRAIC[AV_INaK_INaK] + ALGEBRAIC[AV_INab_INab] + ALGEBRAIC[AV_IKb_IKb] + ALGEBRAIC[AV_IpCa_IpCa] + ALGEBRAIC[AV_ICab_ICab];

    ALGEBRAIC[Istim] = computeIstim(VOI, stimulus);

    if (solveVmWithinODESolver)
    {
        RATES[membrane_V] = -ALGEBRAIC[Iion_cm] - ALGEBRAIC[Istim];
    }
}

#endif

// ************************************************************************* //
