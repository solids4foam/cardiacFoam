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

#ifndef ORd_2011Compact_H
#define ORd_2011Compact_H

#include "ORd_2011Optimized.H"
#include "../ORd/ORdCompactSupport.H"
#include <cmath>

inline double ordSafeExpm1(const double x)
{
    if (std::abs(x) < 1.0e-7)
    {
        const double x2 = x*x;
        return x + 0.5*x2 + x2*x/6.0;
    }

    return std::exp(x) - 1.0;
}

template<bool fullSupport>
inline void
ORdcomputeVariablesCompactImpl
(
    double VOI,
    double* CONSTANTS,
    double* RATES,
    double* STATES,
    double* SUPPORT,
    int tissueFlag,
    bool solveVmWithinODESolver,
    const Foam::StimulusProtocol& stimulus
)
{
#define ORD_SET_SUPPORT(index, value)                                          \
    do                                                                         \
    {                                                                          \
        if constexpr (fullSupport)                                             \
        {                                                                      \
            SUPPORT[index] = value;                                            \
        }                                                                      \
        else if constexpr ((index) == ORdCS_Iion_cm)                           \
        {                                                                      \
            SUPPORT[0] = value;                                                \
        }                                                                      \
    } while (false)

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
    const double expVfrt = std::exp(vfrt);
    const double exp2Vfrt = std::exp(2.0*vfrt);

    const double CaMKb =
        CONSTANTS[AC_CaMKo]*(1.0 - STATES[CaMK_CaMKt])
       /(1.0 + CONSTANTS[AC_KmCaM]/cass);
    const double CaMKa = CaMKb + STATES[CaMK_CaMKt];
    RATES[CaMK_CaMKt] =
        CONSTANTS[AC_aCaMK]*CaMKb*CaMKa
      - CONSTANTS[AC_bCaMK]*STATES[CaMK_CaMKt];

    const double IpCa =
        CONSTANTS[AC_GpCa]*cai/(CONSTANTS[AC_KmCap] + cai);
    const double Jdiff = (cass - cai)/0.2;
    const double JdiffK = (kss - ki)/2.0;
    const double JdiffNa = (nass - nai)/2.0;
    const double Jtr = (cansr - cajsr)/100.0;

    const double kmCaActOverCai = CONSTANTS[AC_KmCaAct]/cai;
    const double kmCaActOverCass = CONSTANTS[AC_KmCaAct]/cass;
    const double allo_i = 1.0/(1.0 + ordSqr(kmCaActOverCai));
    const double allo_ss = 1.0/(1.0 + ordSqr(kmCaActOverCass));
    const double h4_i = 1.0 + nai/CONSTANTS[AC_kna1]*(1.0 + nai/CONSTANTS[AC_kna2]);
    const double h4_ss = 1.0 + nass/CONSTANTS[AC_kna1]*(1.0 + nass/CONSTANTS[AC_kna2]);
    const double hca = std::exp(CONSTANTS[AC_qca]*vfrt);
    const double hna = std::exp(CONSTANTS[AC_qna]*vfrt);
    const double h1_i = 1.0 + nai/CONSTANTS[AC_kna3]*(1.0 + hna);
    const double h1_ss = 1.0 + nass/CONSTANTS[AC_kna3]*(1.0 + hna);
    const double h5_i = ordSqr(nai)/(h4_i*CONSTANTS[AC_kna1]*CONSTANTS[AC_kna2]);
    const double h5_ss = ordSqr(nass)/(h4_ss*CONSTANTS[AC_kna1]*CONSTANTS[AC_kna2]);
    const double h6_i = 1.0/h4_i;
    const double h6_ss = 1.0/h4_ss;
    const double invHna = 1.0/hna;
    const double h7_i = 1.0 + CONSTANTS[AC_nao]/CONSTANTS[AC_kna3]*(1.0 + invHna);
    const double h7_ss = h7_i;
    const double h2_i = nai*hna/(CONSTANTS[AC_kna3]*h1_i);
    const double h2_ss = nass*hna/(CONSTANTS[AC_kna3]*h1_ss);
    const double h3_i = 1.0/h1_i;
    const double h3_ss = 1.0/h1_ss;
    const double h8_i = CONSTANTS[AC_nao]/(CONSTANTS[AC_kna3]*hna*h7_i);
    const double h8_ss = CONSTANTS[AC_nao]/(CONSTANTS[AC_kna3]*hna*h7_ss);
    const double h9_i = 1.0/h7_i;
    const double h9_ss = 1.0/h7_ss;
    const double k6_i = h6_i*cai*CONSTANTS[AC_kcaon];
    const double k6_ss = h6_ss*cass*CONSTANTS[AC_kcaon];
    const double k3p_i = h9_i*CONSTANTS[AC_wca];
    const double k3p_ss = h9_ss*CONSTANTS[AC_wca];
    const double k3pp_i = h8_i*CONSTANTS[AC_wnaca];
    const double k3pp_ss = h8_ss*CONSTANTS[AC_wnaca];
    const double k4p_i = h3_i*CONSTANTS[AC_wca]/hca;
    const double k4p_ss = h3_ss*CONSTANTS[AC_wca]/hca;
    const double k4pp_i = h2_i*CONSTANTS[AC_wnaca];
    const double k4pp_ss = h2_ss*CONSTANTS[AC_wnaca];
    const double k7_i = h5_i*h2_i*CONSTANTS[AC_wna];
    const double k7_ss = h5_ss*h2_ss*CONSTANTS[AC_wna];
    const double k8_i = h8_i*CONSTANTS[AC_h11_i]*CONSTANTS[AC_wna];
    const double k8_ss = h8_ss*CONSTANTS[AC_h11_ss]*CONSTANTS[AC_wna];
    const double k3_i = k3p_i + k3pp_i;
    const double k3_ss = k3p_ss + k3pp_ss;
    const double k4_i = k4p_i + k4pp_i;
    const double k4_ss = k4p_ss + k4pp_ss;
    const double x1_i =
        CONSTANTS[AC_k2_i]*k4_i*(k7_i + k6_i)
      + CONSTANTS[AC_k5_i]*k7_i*(CONSTANTS[AC_k2_i] + k3_i);
    const double x1_ss =
        CONSTANTS[AC_k2_ss]*k4_ss*(k7_ss + k6_ss)
      + CONSTANTS[AC_k5_ss]*k7_ss*(CONSTANTS[AC_k2_ss] + k3_ss);
    const double x2_i =
        CONSTANTS[AC_k1_i]*k7_i*(k4_i + CONSTANTS[AC_k5_i])
      + k4_i*k6_i*(CONSTANTS[AC_k1_i] + k8_i);
    const double x2_ss =
        CONSTANTS[AC_k1_ss]*k7_ss*(k4_ss + CONSTANTS[AC_k5_ss])
      + k4_ss*k6_ss*(CONSTANTS[AC_k1_ss] + k8_ss);
    const double x3_i =
        CONSTANTS[AC_k1_i]*k3_i*(k7_i + k6_i)
      + k8_i*k6_i*(CONSTANTS[AC_k2_i] + k3_i);
    const double x3_ss =
        CONSTANTS[AC_k1_ss]*k3_ss*(k7_ss + k6_ss)
      + k8_ss*k6_ss*(CONSTANTS[AC_k2_ss] + k3_ss);
    const double x4_i =
        CONSTANTS[AC_k2_i]*k8_i*(k4_i + CONSTANTS[AC_k5_i])
      + k3_i*CONSTANTS[AC_k5_i]*(CONSTANTS[AC_k1_i] + k8_i);
    const double x4_ss =
        CONSTANTS[AC_k2_ss]*k8_ss*(k4_ss + CONSTANTS[AC_k5_ss])
      + k3_ss*CONSTANTS[AC_k5_ss]*(CONSTANTS[AC_k1_ss] + k8_ss);
    const double inacaDenI = x1_i + x2_i + x3_i + x4_i;
    const double inacaDenSs = x1_ss + x2_ss + x3_ss + x4_ss;
    const double invInacaDenI = 1.0/inacaDenI;
    const double invInacaDenSs = 1.0/inacaDenSs;
    const double JncxCa_i =
        (x2_i*CONSTANTS[AC_k2_i] - x1_i*CONSTANTS[AC_k1_i])*invInacaDenI;
    const double JncxCa_ss =
        (x2_ss*CONSTANTS[AC_k2_ss] - x1_ss*CONSTANTS[AC_k1_ss])*invInacaDenSs;
    const double JncxNa_i =
        (
            3.0*(x4_i*k7_i - x1_i*k8_i)
          + x3_i*k4pp_i
          - x2_i*k3pp_i
        )*invInacaDenI;
    const double JncxNa_ss =
        (
            3.0*(x4_ss*k7_ss - x1_ss*k8_ss)
          + x3_ss*k4pp_ss
          - x2_ss*k3pp_ss
        )*invInacaDenSs;
    const double INaCa_i =
        0.8*CONSTANTS[AC_Gncx]*allo_i
       *(CONSTANTS[AC_zna]*JncxNa_i + CONSTANTS[AC_zca]*JncxCa_i);
    const double INaCa_ss =
        0.2*CONSTANTS[AC_Gncx]*allo_ss
       *(CONSTANTS[AC_zna]*JncxNa_ss + CONSTANTS[AC_zca]*JncxCa_ss);

    const double deltaVfrtOver3 = CONSTANTS[AC_delta]*vfrt/3.0;
    const double Knai = CONSTANTS[AC_Knai0]*std::exp(deltaVfrtOver3);
    const double Knao =
        CONSTANTS[AC_Knao0]*std::exp((1.0 - CONSTANTS[AC_delta])*vfrt/3.0);
    const double P =
        CONSTANTS[AC_eP]
       /(1.0 + CONSTANTS[AC_H]/CONSTANTS[AC_Khp]
             + nai/CONSTANTS[AC_Knap]
             + ki/CONSTANTS[AC_Kxkur]);
    const double naiOverKnai = nai/Knai;
    const double naoOverKnao = CONSTANTS[AC_nao]/Knao;
    const double kiOverKki = ki/CONSTANTS[AC_Kki];
    const double koOverKko = CONSTANTS[AC_ko]/CONSTANTS[AC_Kko];
    const double nakDenInner =
        ordCube(1.0 + naiOverKnai) + ordSqr(1.0 + kiOverKki) - 1.0;
    const double nakDenOuter =
        ordCube(1.0 + naoOverKnao) + ordSqr(1.0 + koOverKko) - 1.0;
    const double a1 = CONSTANTS[AC_k1p]*ordCube(naiOverKnai)/nakDenInner;
    const double a3 = CONSTANTS[AC_k3p]*ordSqr(koOverKko)/nakDenOuter;
    const double b2 = CONSTANTS[AC_k2m]*ordCube(naoOverKnao)/nakDenOuter;
    const double b3 =
        CONSTANTS[AC_k3m]*P*CONSTANTS[AC_H]
       /(1.0 + CONSTANTS[AC_MgATP]/CONSTANTS[AC_Kmgatp]);
    const double b4 = CONSTANTS[AC_k4m]*ordSqr(kiOverKki)/nakDenInner;
    const double x1 =
        CONSTANTS[AC_a4]*a1*CONSTANTS[AC_a2]
      + b2*b4*b3
      + CONSTANTS[AC_a2]*b4*b3
      + b3*a1*CONSTANTS[AC_a2];
    const double x2 =
        b2*CONSTANTS[AC_b1]*b4
      + a1*CONSTANTS[AC_a2]*a3
      + a3*CONSTANTS[AC_b1]*b4
      + CONSTANTS[AC_a2]*a3*b4;
    const double x3 =
        CONSTANTS[AC_a2]*a3*CONSTANTS[AC_a4]
      + b3*b2*CONSTANTS[AC_b1]
      + b2*CONSTANTS[AC_b1]*CONSTANTS[AC_a4]
      + a3*CONSTANTS[AC_a4]*CONSTANTS[AC_b1];
    const double x4 =
        b4*b3*b2
      + a3*CONSTANTS[AC_a4]*a1
      + b2*CONSTANTS[AC_a4]*a1
      + b3*b2*a1;
    const double inakDen = x1 + x2 + x3 + x4;
    const double invInakDen = 1.0/inakDen;
    const double JnakK =
        2.0*(x4*CONSTANTS[AC_b1] - x3*a1)*invInakDen;
    const double JnakNa =
        3.0*(x1*a3 - x2*b3)*invInakDen;
    const double INaK =
        CONSTANTS[AC_Pnak]
       *(CONSTANTS[AC_zna]*JnakNa + CONSTANTS[AC_zk]*JnakK);

    const double Jleak = 0.0039375*cansr/15.0;
    const double fCaMK = 1.0/(1.0 + CONSTANTS[AC_KmCaMK]/CaMKa);
    const double oneMinusFCaMK = 1.0 - fCaMK;
    const double Jupnp = CONSTANTS[AC_upScale]*0.004375*cai/(cai + 0.00092);
    const double Jupp =
        CONSTANTS[AC_upScale]*2.75*0.004375*cai/(cai + 0.00092 - 0.00017);
    const double Jup =
        oneMinusFCaMK*Jupnp + fCaMK*Jupp - Jleak;

    const double RTonF = R*T/F;
    const double EK = RTonF*std::log(CONSTANTS[AC_ko]/ki);
    const double ENa = RTonF*std::log(CONSTANTS[AC_nao]/nai);
    const double EKs =
        RTonF*std::log
        (
            (CONSTANTS[AC_ko] + CONSTANTS[AC_PKNa]*CONSTANTS[AC_nao])
           /(ki + CONSTANTS[AC_PKNa]*nai)
        );

    const double rk1 =
        1.0/(1.0 + std::exp((Vm + 105.8 - 2.6*CONSTANTS[AC_ko])/9.493));
    const double txk1 =
        122.2/(std::exp(-(Vm + 127.2)/20.36) + std::exp((Vm + 236.8)/69.33));
    const double xk1ss =
        1.0
       /(1.0
         + std::exp
           (
               -(Vm + 2.5538*CONSTANTS[AC_ko] + 144.59)
              /(1.5692*CONSTANTS[AC_ko] + 3.8115)
           ));
    ORD_SET_SUPPORT(ORdCS_txk1, txk1);
    ORD_SET_SUPPORT(ORdCS_xk1ss, xk1ss);
    RATES[IK1_xk1] = (xk1ss - STATES[IK1_xk1])/txk1;
    const double sqrtKo = std::sqrt(CONSTANTS[AC_ko]);
    const double IK1 =
        CONSTANTS[AC_GK1]*sqrtKo*rk1
       *STATES[IK1_xk1]*(Vm - EK);

    const double xkb = 1.0/(1.0 + std::exp(-(Vm - 14.48)/18.34));
    const double IKb = CONSTANTS[AC_GKb]*xkb*(Vm - EK);

    const double Axrf = 1.0/(1.0 + std::exp((Vm + 54.81)/38.21));
    const double rkr =
        1.0/(1.0 + std::exp((Vm + 55.0)/75.0))
       *1.0/(1.0 + std::exp((Vm - 10.0)/30.0));
    const double txrf =
        12.98 + 1.0
       /(0.3652*std::exp((Vm - 31.66)/3.869)
         + 4.123e-05*std::exp(-(Vm - 47.78)/20.38));
    const double txrs =
        1.865 + 1.0
       /(0.06629*std::exp((Vm - 34.7)/7.355)
         + 1.128e-05*std::exp(-(Vm - 29.74)/25.94));
    const double xrss = 1.0/(1.0 + std::exp(-(Vm + 8.337)/6.789));
    ORD_SET_SUPPORT(ORdCS_txrf, txrf);
    ORD_SET_SUPPORT(ORdCS_txrs, txrs);
    ORD_SET_SUPPORT(ORdCS_xrss, xrss);
    const double Axrs = 1.0 - Axrf;
    RATES[IKr_xrf] = (xrss - STATES[IKr_xrf])/txrf;
    RATES[IKr_xrs] = (xrss - STATES[IKr_xrs])/txrs;
    const double xr = Axrf*STATES[IKr_xrf] + Axrs*STATES[IKr_xrs];
    const double IKr =
        CONSTANTS[AC_GKr]*std::sqrt(CONSTANTS[AC_ko]/5.4)*xr*rkr*(Vm - EK);

    const double KsCa = 1.0 + 0.6/(1.0 + std::pow(3.8e-05/cai, 1.4));
    const double txs1 =
        817.3 + 1.0
       /(0.0002326*std::exp((Vm + 48.28)/17.8)
         + 0.001292*std::exp(-(Vm + 210.0)/230.0));
    const double txs2 =
        1.0
       /(0.01*std::exp((Vm - 50.0)/20.0)
         + 0.0193*std::exp(-(Vm + 66.54)/31.0));
    const double xs1ss = 1.0/(1.0 + std::exp(-(Vm + 11.6)/8.932));
    const double xs2ss = xs1ss;
    ORD_SET_SUPPORT(ORdCS_txs1, txs1);
    ORD_SET_SUPPORT(ORdCS_txs2, txs2);
    ORD_SET_SUPPORT(ORdCS_xs1ss, xs1ss);
    ORD_SET_SUPPORT(ORdCS_xs2ss, xs2ss);
    RATES[IKs_xs1] = (xs1ss - STATES[IKs_xs1])/txs1;
    RATES[IKs_xs2] = (xs2ss - STATES[IKs_xs2])/txs2;
    const double IKs =
        CONSTANTS[AC_GKs]*KsCa*STATES[IKs_xs1]*STATES[IKs_xs2]*(Vm - EKs);

    const double hssp = 1.0/(1.0 + std::exp((Vm + 89.1)/6.086));
    const double thf =
        1.0
       /(1.432e-05*std::exp(-(Vm + 1.196)/6.285)
         + 6.149*std::exp((Vm + 0.5096)/20.27));
    const double ths =
        1.0
       /(0.009794*std::exp(-(Vm + 17.95)/28.05)
         + 0.3343*std::exp((Vm + 5.73)/56.66));
    const double tj =
        2.038 + 1.0
       /(0.02136*std::exp(-(Vm + 100.6)/8.281)
         + 0.3052*std::exp((Vm + 0.9941)/38.45));
    const double hss =
        1.0/(1.0 + std::exp((Vm + CONSTANTS[AC_hssV1])/CONSTANTS[AC_hssV2]));
    const double mss =
        1.0/(1.0 + std::exp(-(Vm + CONSTANTS[AC_mssV1])/CONSTANTS[AC_mssV2]));
    const double thsp = 3.0*ths;
    const double tjp = 1.46*tj;
    const double tm =
        1.0
       /(CONSTANTS[AC_mtD1]*std::exp((Vm + CONSTANTS[AC_mtV1])/CONSTANTS[AC_mtV2])
         + CONSTANTS[AC_mtD2]*std::exp(-(Vm + CONSTANTS[AC_mtV3])/CONSTANTS[AC_mtV4]));
    const double h =
        CONSTANTS[AC_Ahf]*STATES[INa_hf] + CONSTANTS[AC_Ahs]*STATES[INa_hs];
    const double hp =
        CONSTANTS[AC_Ahf]*STATES[INa_hf] + CONSTANTS[AC_Ahs]*STATES[INa_hsp];
    const double jss = hss;
    ORD_SET_SUPPORT(ORdCS_hssp, hssp);
    ORD_SET_SUPPORT(ORdCS_thf, thf);
    ORD_SET_SUPPORT(ORdCS_ths, ths);
    ORD_SET_SUPPORT(ORdCS_tj, tj);
    ORD_SET_SUPPORT(ORdCS_hss, hss);
    ORD_SET_SUPPORT(ORdCS_mss, mss);
    ORD_SET_SUPPORT(ORdCS_thsp, thsp);
    ORD_SET_SUPPORT(ORdCS_tjp, tjp);
    ORD_SET_SUPPORT(ORdCS_tm, tm);
    ORD_SET_SUPPORT(ORdCS_jss, jss);
    RATES[INa_hf] = (hss - STATES[INa_hf])/thf;
    RATES[INa_hs] = (hss - STATES[INa_hs])/ths;
    RATES[INa_hsp] = (hssp - STATES[INa_hsp])/thsp;
    RATES[INa_m] = (mss - STATES[INa_m])/tm;
    const double INa =
        CONSTANTS[AC_GNa]*(Vm - ENa)*ordCube(STATES[INa_m])
       *(oneMinusFCaMK*h*STATES[INa_j] + fCaMK*hp*STATES[INa_jp]);
    RATES[INa_j] = (jss - STATES[INa_j])/tj;
    RATES[INa_jp] = (jss - STATES[INa_jp])/tjp;

    const double AiF = 1.0/(1.0 + std::exp((Vm - 213.6)/151.2));
    const double ass = 1.0/(1.0 + std::exp(-(Vm - 14.34)/14.82));
    const double assp = 1.0/(1.0 + std::exp(-(Vm - 24.34)/14.82));
    const double delta_epi =
        (tissueFlag == 1)
      ? 1.0 - 0.95/(1.0 + std::exp((Vm + 70.0)/5.0))
      : 1.0;
    const double dti_develop =
        1.354
      + 0.0001
       /(std::exp((Vm - 167.4)/15.89) + std::exp(-(Vm - 12.23)/0.2154));
    const double dti_recover =
        1.0 - 0.5/(1.0 + std::exp((Vm + 70.0)/20.0));
    const double iss = 1.0/(1.0 + std::exp((Vm + 43.94)/5.711));
    const double ta =
        1.0515
       /(1.0/(1.2089*(1.0 + std::exp(-(Vm - 18.4099)/29.3814)))
         + 3.5/(1.0 + std::exp((Vm + 100.0)/29.3814)));
    const double tiF_b =
        4.562 + 1.0
       /(0.3933*std::exp(-(Vm + 100.0)/100.0)
         + 0.08004*std::exp((Vm + 50.0)/16.59));
    const double tiS_b =
        23.62 + 1.0
       /(0.001416*std::exp(-(Vm + 96.52)/59.05)
         + 1.78e-08*std::exp((Vm + 114.1)/8.079));
    const double AiS = 1.0 - AiF;
    const double tiF = tiF_b*delta_epi;
    const double tiS = tiS_b*delta_epi;
    ORD_SET_SUPPORT(ORdCS_ass, ass);
    ORD_SET_SUPPORT(ORdCS_assp, assp);
    ORD_SET_SUPPORT(ORdCS_iss, iss);
    ORD_SET_SUPPORT(ORdCS_ta, ta);
    ORD_SET_SUPPORT(ORdCS_tiF, tiF);
    ORD_SET_SUPPORT(ORdCS_tiS, tiS);
    RATES[Ito_a] = (ass - STATES[Ito_a])/ta;
    RATES[Ito_ap] = (assp - STATES[Ito_ap])/ta;
    const double i = AiF*STATES[Ito_iF] + AiS*STATES[Ito_iS];
    const double ip = AiF*STATES[Ito_iFp] + AiS*STATES[Ito_iSp];
    const double tiFp = dti_develop*dti_recover*tiF;
    const double tiSp = dti_develop*dti_recover*tiS;
    ORD_SET_SUPPORT(ORdCS_tiFp, tiFp);
    ORD_SET_SUPPORT(ORdCS_tiSp, tiSp);
    RATES[Ito_iF] = (iss - STATES[Ito_iF])/tiF;
    RATES[Ito_iS] = (iss - STATES[Ito_iS])/tiS;
    const double Ito =
        CONSTANTS[AC_Gto]*(Vm - EK)
       *(oneMinusFCaMK*STATES[Ito_a]*i + fCaMK*STATES[Ito_ap]*ip);
    RATES[Ito_iFp] = (iss - STATES[Ito_iFp])/tiFp;
    RATES[Ito_iSp] = (iss - STATES[Ito_iSp])/tiSp;

    const double hLss = 1.0/(1.0 + std::exp((Vm + 87.61)/7.488));
    const double hLssp = 1.0/(1.0 + std::exp((Vm + 93.81)/7.488));
    const double mLss = 1.0/(1.0 + std::exp(-(Vm + 42.85)/5.264));
    const double tmL = tm;
    ORD_SET_SUPPORT(ORdCS_hLss, hLss);
    ORD_SET_SUPPORT(ORdCS_hLssp, hLssp);
    ORD_SET_SUPPORT(ORdCS_mLss, mLss);
    ORD_SET_SUPPORT(ORdCS_tmL, tmL);
    RATES[INaL_hL] = (hLss - STATES[INaL_hL])/CONSTANTS[AC_thL];
    RATES[INaL_mL] = (mLss - STATES[INaL_mL])/tmL;
    const double INaL =
        CONSTANTS[AC_GNaL]*(Vm - ENa)*STATES[INaL_mL]
       *(oneMinusFCaMK*STATES[INaL_hL] + fCaMK*STATES[INaL_hLp]);
    RATES[INaL_hLp] = (hLssp - STATES[INaL_hLp])/CONSTANTS[AC_thLp];

    const double Afcaf = 0.3 + 0.6/(1.0 + std::exp((Vm - 10.0)/10.0));
    const double dss = 1.0/(1.0 + std::exp(-(Vm + 3.94)/4.23));
    const double fss = 1.0/(1.0 + std::exp((Vm + 19.58)/3.696));
    const double km2n = STATES[ICaL_jca];
    const double td =
        0.6 + 1.0/(std::exp(-0.05*(Vm + 6.0)) + std::exp(0.09*(Vm + 14.0)));
    const double tfcaf =
        7.0
      + 1.0
       /(0.04*std::exp(-(Vm - 4.0)/7.0) + 0.04*std::exp((Vm - 4.0)/7.0));
    const double tfcas =
        100.0
      + 1.0
       /(0.00012*std::exp(-Vm/3.0) + 0.00012*std::exp(Vm/7.0));
    const double tff =
        7.0
      + 1.0
       /(0.0045*std::exp(-(Vm + 20.0)/10.0)
         + 0.0045*std::exp((Vm + 20.0)/10.0));
    const double tfs =
        1000.0
      + 1.0
       /(3.5e-05*std::exp(-(Vm + 5.0)/4.0)
         + 3.5e-05*std::exp((Vm + 5.0)/6.0));
    const double Afcas = 1.0 - Afcaf;
    const double anca =
        1.0/(CONSTANTS[AC_k2n]/km2n
            + ordPow4(1.0 + CONSTANTS[AC_Kmn]/cass));
    const double fcass = fss;
    const double tfcafp = 2.5*tfcaf;
    const double tffp = 2.5*tff;
    ORD_SET_SUPPORT(ORdCS_dss, dss);
    ORD_SET_SUPPORT(ORdCS_fss, fss);
    ORD_SET_SUPPORT(ORdCS_km2n, km2n);
    ORD_SET_SUPPORT(ORdCS_td, td);
    ORD_SET_SUPPORT(ORdCS_tfcaf, tfcaf);
    ORD_SET_SUPPORT(ORdCS_tfcas, tfcas);
    ORD_SET_SUPPORT(ORdCS_tff, tff);
    ORD_SET_SUPPORT(ORdCS_tfs, tfs);
    ORD_SET_SUPPORT(ORdCS_fcass, fcass);
    ORD_SET_SUPPORT(ORdCS_tfcafp, tfcafp);
    ORD_SET_SUPPORT(ORdCS_tffp, tffp);
    RATES[ICaL_d] = (dss - STATES[ICaL_d])/td;
    RATES[ICaL_ff] = (fss - STATES[ICaL_ff])/tff;
    RATES[ICaL_fs] = (fss - STATES[ICaL_fs])/tfs;
    const double f = CONSTANTS[AC_Aff]*STATES[ICaL_ff]
                   + CONSTANTS[AC_Afs]*STATES[ICaL_fs];
    const double fca = Afcaf*STATES[ICaL_fcaf] + Afcas*STATES[ICaL_fcas];
    const double fcap = Afcaf*STATES[ICaL_fcafp] + Afcas*STATES[ICaL_fcas];
    const double fp = CONSTANTS[AC_Aff]*STATES[ICaL_ffp]
                    + CONSTANTS[AC_Afs]*STATES[ICaL_fs];
    RATES[ICaL_fcaf] = (fcass - STATES[ICaL_fcaf])/tfcaf;
    RATES[ICaL_fcafp] = (fcass - STATES[ICaL_fcafp])/tfcafp;
    RATES[ICaL_fcas] = (fcass - STATES[ICaL_fcas])/tfcas;
    RATES[ICaL_ffp] = (fss - STATES[ICaL_ffp])/tffp;
    RATES[ICaL_jca] = (fcass - STATES[ICaL_jca])/CONSTANTS[AC_tjca];
    RATES[ICaL_nca] = anca*CONSTANTS[AC_k2n] - STATES[ICaL_nca]*km2n;

    RATES[Ca_nsr] = Jup - Jtr*CONSTANTS[AC_vjsr]/CONSTANTS[AC_vnsr];
    const double Bcajsr =
        1.0
       /(1.0
         + CONSTANTS[AC_csqnmax]*CONSTANTS[AC_kmcsqn]
          /ordSqr(CONSTANTS[AC_kmcsqn] + cajsr));
    const double Bcass =
        1.0
       /(1.0
         + CONSTANTS[AC_BSRmax]*CONSTANTS[AC_KmBSR]
          /ordSqr(CONSTANTS[AC_KmBSR] + cass)
         + CONSTANTS[AC_BSLmax]*CONSTANTS[AC_KmBSL]
          /ordSqr(CONSTANTS[AC_KmBSL] + cass));
    const double Bcai =
        1.0
       /(1.0
         + CONSTANTS[AC_cmdnmax]*CONSTANTS[AC_kmcmdn]
          /ordSqr(CONSTANTS[AC_kmcmdn] + cai)
         + CONSTANTS[AC_trpnmax]*CONSTANTS[AC_kmtrpn]
          /ordSqr(CONSTANTS[AC_kmtrpn] + cai));

    const double Jrel =
        oneMinusFCaMK*STATES[ryr_Jrelnp] + fCaMK*STATES[ryr_Jrelp];
    const double tau_rel_temp = CONSTANTS[AC_bt]/(1.0 + 0.0123/cajsr);
    const double tau_relp_temp = CONSTANTS[AC_btp]/(1.0 + 0.0123/cajsr);
    const double tau_rel = tau_rel_temp < 0.001 ? 0.001 : tau_rel_temp;
    const double tau_relp = tau_relp_temp < 0.001 ? 0.001 : tau_relp_temp;
    ORD_SET_SUPPORT(ORdCS_tau_rel, tau_rel);
    ORD_SET_SUPPORT(ORdCS_tau_relp, tau_relp);

    const double expVfrtMinus1 = ordSafeExpm1(vfrt);
    const double exp2VfrtMinus1 = ordSafeExpm1(2.0*vfrt);
    const double PhiCaK =
        vffrt*(0.75*kss*expVfrt - 0.75*CONSTANTS[AC_ko])/expVfrtMinus1;
    const double PhiCaL =
        4.0*vffrt*(cass*exp2Vfrt - 0.341*CONSTANTS[AC_cao])/exp2VfrtMinus1;
    const double PhiCaNa =
        vffrt*(0.75*nass*expVfrt - 0.75*CONSTANTS[AC_nao])/expVfrtMinus1;
    const double ICab =
        CONSTANTS[AC_PCab]*4.0*vffrt
       *(cai*exp2Vfrt - 0.341*CONSTANTS[AC_cao])/exp2VfrtMinus1;
    const double INab =
        CONSTANTS[AC_PNab]*vffrt*(nai*expVfrt - CONSTANTS[AC_nao])/expVfrtMinus1;
    const double nca = STATES[ICaL_nca];
    const double jca = STATES[ICaL_jca];
    const double dGate = STATES[ICaL_d];
    const double icaGate =
        dGate*(f*(1.0 - nca) + jca*fca*nca);
    const double icaGateP =
        dGate*(fp*(1.0 - nca) + jca*fcap*nca);
    const double ICaK =
        PhiCaK
       *(
            oneMinusFCaMK*CONSTANTS[AC_PCaK]*icaGate
          + fCaMK*CONSTANTS[AC_PCaKp]*icaGateP
        );
    const double ICaL =
        PhiCaL
       *(
            oneMinusFCaMK*CONSTANTS[AC_PCa]*icaGate
          + fCaMK*CONSTANTS[AC_PCap]*icaGateP
        );
    const double ICaNa =
        PhiCaNa
       *(
            oneMinusFCaMK*CONSTANTS[AC_PCaNa]*icaGate
          + fCaMK*CONSTANTS[AC_PCaNap]*icaGateP
        );

    const double iStim = computeIstim(VOI, stimulus);
    const double cmAcap = CONSTANTS[AC_cm]*CONSTANTS[AC_Acap];
    const double cmAcapOverFVmyo = cmAcap/(F*CONSTANTS[AC_vmyo]);
    const double cmAcapOverFVss = cmAcap/(F*CONSTANTS[AC_vss]);
    const double cmAcapOver2FVmyo = 0.5*cmAcapOverFVmyo;
    const double cmAcapOver2FVss = 0.5*cmAcapOverFVss;

    RATES[Ca_jsr] = Bcajsr*(Jtr - Jrel);
    RATES[K_i] =
        -(Ito + IKr + IKs + IK1 + IKb + iStim - 2.0*INaK)
        *cmAcapOverFVmyo
      + JdiffK*CONSTANTS[AC_vss]/CONSTANTS[AC_vmyo];
    RATES[Ca_i] =
        Bcai
       *(
            -(IpCa + ICab - 2.0*INaCa_i)
             *cmAcapOver2FVmyo
          - Jup*CONSTANTS[AC_vnsr]/CONSTANTS[AC_vmyo]
          + Jdiff*CONSTANTS[AC_vss]/CONSTANTS[AC_vmyo]
        );
    RATES[Na_i] =
        -(INa + INaL + 3.0*INaCa_i + 3.0*INaK + INab)
        *cmAcapOverFVmyo
      + JdiffNa*CONSTANTS[AC_vss]/CONSTANTS[AC_vmyo];
    RATES[Ca_ss] =
        Bcass
       *(
            -(ICaL - 2.0*INaCa_ss)
             *cmAcapOver2FVss
          + Jrel*CONSTANTS[AC_vjsr]/CONSTANTS[AC_vss]
          - Jdiff
        );
    RATES[K_ss] =
        -ICaK*cmAcapOverFVss - JdiffK;
    RATES[Na_ss] =
        -(ICaNa + 3.0*INaCa_ss)
        *cmAcapOverFVss - JdiffNa;

    const double jrelPow = ordPow8(1.5/cajsr);
    const double Jrel_inf_temp = CONSTANTS[AC_a_rel]*-ICaL/(1.0 + jrelPow);
    const double Jrel_infp_temp = CONSTANTS[AC_a_relp]*-ICaL/(1.0 + jrelPow);
    const double Jrel_inf = tissueFlag == 2 ? Jrel_inf_temp*1.7 : Jrel_inf_temp;
    const double Jrel_infp =
        tissueFlag == 2 ? Jrel_infp_temp*1.7 : Jrel_infp_temp;
    ORD_SET_SUPPORT(ORdCS_Jrel_inf, Jrel_inf);
    ORD_SET_SUPPORT(ORdCS_Jrel_infp, Jrel_infp);
    RATES[ryr_Jrelnp] = (Jrel_inf - STATES[ryr_Jrelnp])/tau_rel;
    RATES[ryr_Jrelp] = (Jrel_infp - STATES[ryr_Jrelp])/tau_relp;

    const double Iion_cm =
        INa + INaL + Ito + ICaL + ICaNa + ICaK + IKr + IKs + IK1 + INaCa_i
      + INaCa_ss + INaK + INab + IKb + IpCa + ICab;
    ORD_SET_SUPPORT(ORdCS_Iion_cm, Iion_cm);

    if (solveVmWithinODESolver)
    {
        RATES[membrane_V] = -Iion_cm - iStim;
    }

#undef ORD_SET_SUPPORT
}

inline void
ORdcomputeVariablesCompact
(
    double VOI,
    double* CONSTANTS,
    double* RATES,
    double* STATES,
    double* SUPPORT,
    int supportCount,
    int tissueFlag,
    bool solveVmWithinODESolver,
    const Foam::StimulusProtocol& stimulus
)
{
    if (supportCount >= NUM_ORD_COMPACT_SUPPORT)
    {
        ORdcomputeVariablesCompactImpl<true>
        (
            VOI,
            CONSTANTS,
            RATES,
            STATES,
            SUPPORT,
            tissueFlag,
            solveVmWithinODESolver,
            stimulus
        );
    }
    else
    {
        ORdcomputeVariablesCompactImpl<false>
        (
            VOI,
            CONSTANTS,
            RATES,
            STATES,
            SUPPORT,
            tissueFlag,
            solveVmWithinODESolver,
            stimulus
        );
    }
}

#endif

// ************************************************************************* //
