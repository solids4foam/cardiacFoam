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
    Exact activation time psi, its gradient, the manufactured eikonal source
    and the conductivity sigma of the heart-only eikonal manufactured
    solution on the spherical shell, sealed on both spheres, with the
    stimulus point xs, the pinned radius rPin and self-test samples.

    Regenerate from generator/:
        python shell_eikonal_mms.py --write ../eikonalVerification/manufacturedShellEikonalReference.H

    sympy 1.14.0
\*---------------------------------------------------------------------------*/

#ifndef manufacturedShellEikonalReference_H
#define manufacturedShellEikonalReference_H

#include "error.H"
#include "label.H"
#include "mathematicalConstants.H"
#include "point.H"
#include "scalar.H"
#include "symmTensor.H"
#include "unitConversion.H"
#include "vector.H"

namespace Foam
{
namespace shellEikonalMMS
{

struct Params
{
    scalar r1, r2, gamma1, gamma2;
    scalar gil, git, gin, theta0, theta1;
    scalar chi, cm, c0, delta, rPin;
    scalar xs[3];
    scalar Msi[6];
    label solution;
};

inline Params defaultParams()
{
    Params p;
    p.r1 = 0.5;
    p.r2 = 1.0;
    p.gamma1 = degToRad(20.0);
    p.gamma2 = degToRad(50.0);
    p.gil = 0.17;
    p.git = 0.019;
    p.gin = 0.01;
    p.theta0 = degToRad(-20.0);
    p.theta1 = degToRad(90.0);
    p.chi = 3.0;
    p.cm = 2.0;
    p.c0 = 1.0;
    p.delta = 0.05;
    p.rPin = 0.15;
    p.xs[0] = 0.7055406512876979;
    p.xs[1] = 0.21166219538630934;
    p.xs[2] = 0.1411081302575396;
    p.Msi[0] = 142.02325256766676;
    p.Msi[1] = -171.84261697625686;
    p.Msi[2] = -5.770430748139214;
    p.Msi[3] = 390.72583812334926;
    p.Msi[4] = -164.3186863942699;
    p.Msi[5] = 418.33450064025334;
    p.solution = 0;
    return p;
}

inline Params alternateParams()
{
    Params p;
    p.r1 = 0.5;
    p.r2 = 1.0;
    p.gamma1 = degToRad(20.0);
    p.gamma2 = degToRad(50.0);
    p.gil = 0.17;
    p.git = 0.019;
    p.gin = 0.01;
    p.theta0 = degToRad(-20.0);
    p.theta1 = degToRad(70.0);
    p.chi = 3.0;
    p.cm = 2.0;
    p.c0 = 1.0;
    p.delta = 0.06;
    p.rPin = 0.12;
    p.xs[0] = 0.6315809309675278;
    p.xs[1] = 0.24561480648737188;
    p.xs[2] = 0.1754391474909799;
    p.Msi[0] = 92.87501466313385;
    p.Msi[1] = -149.30730396406722;
    p.Msi[2] = 39.49946725392798;
    p.Msi[3] = 425.2788000046807;
    p.Msi[4] = -134.8089221185726;
    p.Msi[5] = 432.9297766634548;
    p.solution = 0;
    return p;
}

inline scalar psi(const point& pt, const Params& p)
{
    const scalar r1 = p.r1;
    const scalar r2 = p.r2;
    const scalar gamma1 = p.gamma1;
    const scalar gamma2 = p.gamma2;
    const scalar gil = p.gil;
    const scalar git = p.git;
    const scalar gin = p.gin;
    const scalar theta0 = p.theta0;
    const scalar theta1 = p.theta1;
    const scalar chi = p.chi;
    const scalar cm = p.cm;
    const scalar c0 = p.c0;
    const scalar delta = p.delta;
    const scalar xs0 = p.xs[0];
    const scalar xs1 = p.xs[1];
    const scalar xs2 = p.xs[2];
    const scalar Msi0 = p.Msi[0];
    const scalar Msi1 = p.Msi[1];
    const scalar Msi2 = p.Msi[2];
    const scalar Msi3 = p.Msi[3];
    const scalar Msi4 = p.Msi[4];
    const scalar Msi5 = p.Msi[5];
    const scalar x = pt.x();
    const scalar y = pt.y();
    const scalar z = pt.z();

    const scalar w0 = (1.0/(c0));
    const scalar w1 = x - xs0;
    const scalar w2 = -xs1 + y;
    const scalar w3 = -xs2 + z;
    const scalar w4 = Msi0*w1 + Msi1*w2 + Msi2*w3;
    const scalar w5 = Msi1*w1 + Msi3*w2 + Msi4*w3;
    const scalar w6 = Msi2*w1 + Msi4*w2 + Msi5*w3;
    const scalar w7 = sqrt((delta*delta) + w1*w4 + w2*w5 + w3*w6);
    const scalar w8 = r1 + r2;
    const scalar w9 = (1.0/(r2));
    const scalar w10 = cos(gamma2);
    const scalar w11 = sin(gamma2);
    const scalar w12 = (w11*w11);
    const scalar w13 = (r1*r1);
    const scalar w14 = -(r2*r2);
    const scalar w15 = (x*x) + (y*y) + (z*z);
    const scalar w16 = (-w13 + w15)/(-w13 - w14);
    const scalar w17 = theta0 + theta1*w16;
    const scalar w18 = sin(w17);
    const scalar w19 = cos(w17);
    const scalar w20 = gil*w12*(w18*w18) + gin*(w10*w10) + git*w12*(w19*w19);
    const scalar w21 = w9*z;
    const scalar w22 = 1/(chi*cm);
    const scalar w23 = w21*w22;
    const scalar w24 = cos(gamma1);
    const scalar w25 = gin*w10*w11;
    const scalar w26 = sin(gamma1);
    const scalar w27 = w19*w26;
    const scalar w28 = w18*w24;
    const scalar w29 = w10*w28 + w27;
    const scalar w30 = w11*w18;
    const scalar w31 = w18*w26;
    const scalar w32 = w10*w19*w24 - w31;
    const scalar w33 = w11*w19;
    const scalar w34 = gil*w29*w30 + git*w32*w33 - w24*w25;
    const scalar w35 = w9*y;
    const scalar w36 = w22*w35;
    const scalar w37 = -w10*w31 + w19*w24;
    const scalar w38 = gil*w37;
    const scalar w39 = -w10*w27 - w28;
    const scalar w40 = git*w39;
    const scalar w41 = w25*w26 + w30*w38 + w33*w40;
    const scalar w42 = w9*x;
    const scalar w43 = w22*w42;
    const scalar w44 = gin*w12;
    const scalar w45 = gil*(w29*w29) + git*(w32*w32) + (w24*w24)*w44;
    const scalar w46 = -w24*w26*w44 + w29*w38 + w32*w40;
    const scalar w47 = gil*(w37*w37) + git*(w39*w39) + (w26*w26)*w44;
    const scalar w48 = w0*w22/w7;
    const scalar w49 = w48*w6;
    const scalar w50 = w48*w5;
    const scalar w51 = w4*w48;
    const scalar w52 = w20*w49 + w34*w50 + w41*w51;
    const scalar w53 = w34*w49 + w45*w50 + w46*w51;
    const scalar w54 = w41*w49 + w46*w50 + w47*w51;
    const scalar w55 = (1.0/2.0)*w8;
    const scalar w56 = (1.0/(r1));
    const scalar w57 = w56*z;
    const scalar w58 = w22*w57;
    const scalar w59 = w56*y;
    const scalar w60 = w22*w59;
    const scalar w61 = w56*x;
    const scalar w62 = w22*w61;
    return w0*w7 + w16*(w14 + w15)*(w16*w55*w9*(-w21*w52 - w35*w53 - w42*w54)/(w21*(w20*w23 + w34*w36 + w41*w43) + w35*(w23*w34 + w36*w45 + w43*w46) + w42*(w23*w41 + w36*w46 + w43*w47)) + w55*w56*(1 - w16)*(w52*w57 + w53*w59 + w54*w61)/(-w57*(-w20*w58 - w34*w60 - w41*w62) - w59*(-w34*w58 - w45*w60 - w46*w62) - w61*(-w41*w58 - w46*w60 - w47*w62)))/w8;
}

inline vector gradPsi(const point& pt, const Params& p)
{
    const scalar r1 = p.r1;
    const scalar r2 = p.r2;
    const scalar gamma1 = p.gamma1;
    const scalar gamma2 = p.gamma2;
    const scalar gil = p.gil;
    const scalar git = p.git;
    const scalar gin = p.gin;
    const scalar theta0 = p.theta0;
    const scalar theta1 = p.theta1;
    const scalar chi = p.chi;
    const scalar cm = p.cm;
    const scalar c0 = p.c0;
    const scalar delta = p.delta;
    const scalar xs0 = p.xs[0];
    const scalar xs1 = p.xs[1];
    const scalar xs2 = p.xs[2];
    const scalar Msi0 = p.Msi[0];
    const scalar Msi1 = p.Msi[1];
    const scalar Msi2 = p.Msi[2];
    const scalar Msi3 = p.Msi[3];
    const scalar Msi4 = p.Msi[4];
    const scalar Msi5 = p.Msi[5];
    const scalar x = pt.x();
    const scalar y = pt.y();
    const scalar z = pt.z();

    const scalar w0 = x - xs0;
    const scalar w1 = -xs1 + y;
    const scalar w2 = -xs2 + z;
    const scalar w3 = Msi0*w0 + Msi1*w1 + Msi2*w2;
    const scalar w4 = (1.0/(c0));
    const scalar w5 = Msi1*w0 + Msi3*w1 + Msi4*w2;
    const scalar w6 = Msi2*w0 + Msi4*w1 + Msi5*w2;
    const scalar w7 = (delta*delta) + w0*w3 + w1*w5 + w2*w6;
    const scalar w8 = w4/sqrt(w7);
    const scalar w9 = w3*w8;
    const scalar w10 = r1 + r2;
    const scalar w11 = (1.0/(w10));
    const scalar w12 = (r1*r1);
    const scalar w13 = -(r2*r2);
    const scalar w14 = (1.0/((-w12 - w13)));
    const scalar w15 = (x*x) + (y*y) + (z*z);
    const scalar w16 = w14*(-w12 + w15);
    const scalar w17 = w11*w16;
    const scalar w18 = (1.0/(r2));
    const scalar w19 = cos(gamma2);
    const scalar w20 = theta0 + theta1*w16;
    const scalar w21 = sin(w20);
    const scalar w22 = sin(gamma2);
    const scalar w23 = (w22*w22);
    const scalar w24 = gil*w23;
    const scalar w25 = cos(w20);
    const scalar w26 = git*w23;
    const scalar w27 = gin*(w19*w19) + (w21*w21)*w24 + (w25*w25)*w26;
    const scalar w28 = w6*w8;
    const scalar w29 = 1/(chi*cm);
    const scalar w30 = w28*w29;
    const scalar w31 = cos(gamma1);
    const scalar w32 = gin*w19*w22;
    const scalar w33 = sin(gamma1);
    const scalar w34 = w25*w33;
    const scalar w35 = w21*w31;
    const scalar w36 = w19*w35;
    const scalar w37 = w34 + w36;
    const scalar w38 = gil*w37;
    const scalar w39 = w21*w22;
    const scalar w40 = w21*w33;
    const scalar w41 = w19*w25*w31 - w40;
    const scalar w42 = git*w41;
    const scalar w43 = w22*w25;
    const scalar w44 = -w31*w32 + w38*w39 + w42*w43;
    const scalar w45 = w5*w8;
    const scalar w46 = w29*w45;
    const scalar w47 = w19*w40;
    const scalar w48 = w25*w31 - w47;
    const scalar w49 = gil*w48;
    const scalar w50 = w19*w34;
    const scalar w51 = -w35 - w50;
    const scalar w52 = git*w51;
    const scalar w53 = w32*w33 + w39*w49 + w43*w52;
    const scalar w54 = w29*w53;
    const scalar w55 = w27*w30 + w44*w46 + w54*w9;
    const scalar w56 = w18*w55;
    const scalar w57 = gin*w23;
    const scalar w58 = gil*(w37*w37) + git*(w41*w41) + (w31*w31)*w57;
    const scalar w59 = -w31*w33*w57 + w37*w49 + w41*w52;
    const scalar w60 = w29*w59;
    const scalar w61 = w30*w44 + w46*w58 + w60*w9;
    const scalar w62 = w18*w61;
    const scalar w63 = gil*(w48*w48) + git*(w51*w51) + (w33*w33)*w57;
    const scalar w64 = w29*w63;
    const scalar w65 = w28*w54 + w45*w60 + w64*w9;
    const scalar w66 = w18*w65;
    const scalar w67 = -w56*z - w62*y - w66*x;
    const scalar w68 = w18*z;
    const scalar w69 = w29*w68;
    const scalar w70 = w18*y;
    const scalar w71 = w29*w70;
    const scalar w72 = w18*x;
    const scalar w73 = w29*w72;
    const scalar w74 = w18*(w27*w69 + w44*w71 + w53*w73);
    const scalar w75 = w18*(w44*w69 + w58*w71 + w59*w73);
    const scalar w76 = w18*(w53*w69 + w59*w71 + w63*w73);
    const scalar w77 = w74*z + w75*y + w76*x;
    const scalar w78 = w10/w77;
    const scalar w79 = w67*w78;
    const scalar w80 = (1.0/2.0)*w16*w18;
    const scalar w81 = (1.0/(r1));
    const scalar w82 = w55*w81;
    const scalar w83 = w61*w81;
    const scalar w84 = w65*w81;
    const scalar w85 = w82*z + w83*y + w84*x;
    const scalar w86 = w81*z;
    const scalar w87 = w29*w86;
    const scalar w88 = w81*y;
    const scalar w89 = w29*w88;
    const scalar w90 = w81*x;
    const scalar w91 = w29*w90;
    const scalar w92 = w81*(-w27*w87 - w44*w89 - w53*w91);
    const scalar w93 = w81*(-w44*w87 - w58*w89 - w59*w91);
    const scalar w94 = w81*(-w53*w87 - w59*w89 - w63*w91);
    const scalar w95 = -w92*z - w93*y - w94*x;
    const scalar w96 = w10/w95;
    const scalar w97 = (1.0/2.0)*w81*(1 - w16);
    const scalar w98 = w96*w97;
    const scalar w99 = 2*w79*w80 + 2*w85*w98;
    const scalar w100 = w17*w99;
    const scalar w101 = w14*x;
    const scalar w102 = w13 + w15;
    const scalar w103 = w102*w11*w99;
    const scalar w104 = w14*w79;
    const scalar w105 = w14*w85*w96;
    const scalar w106 = theta1*w101;
    const scalar w107 = 4*w106;
    const scalar w108 = w21*w25;
    const scalar w109 = w107*w108;
    const scalar w110 = w109*w24 - w109*w26;
    const scalar w111 = w18*w29;
    const scalar w112 = w111*w53;
    const scalar w113 = 2*w106;
    const scalar w114 = w113*w43;
    const scalar w115 = w113*w39;
    const scalar w116 = 2*theta1*w14*w19*w25*w31*x - w113*w40;
    const scalar w117 = gil*w39;
    const scalar w118 = -w113*w34 - w113*w36;
    const scalar w119 = git*w43;
    const scalar w120 = w114*w38 - w115*w42 + w116*w117 + w118*w119;
    const scalar w121 = -w113*w35 - w113*w50;
    const scalar w122 = w25*w31;
    const scalar w123 = -w113*w122 + w113*w47;
    const scalar w124 = w114*w49 - w115*w52 + w117*w121 + w119*w123;
    const scalar w125 = w49*(-w107*w35 - w107*w50) + w52*(-w107*w122 + w107*w47);
    const scalar w126 = w116*w49 + w118*w52 + w121*w38 + w123*w42;
    const scalar w127 = w111*w59;
    const scalar w128 = w38*(4*theta1*w14*w19*w25*w31*x - w107*w40) + w42*(-w107*w34 - w107*w36);
    const scalar w129 = w10*w67*w80/(w77*w77);
    const scalar w130 = w29*w81;
    const scalar w131 = w130*w53;
    const scalar w132 = w130*w59;
    const scalar w133 = w10*w85*w97/(w95*w95);
    const scalar w134 = Msi2*w8;
    const scalar w135 = w134*w29;
    const scalar w136 = w4/pow(w7, 1.5);
    const scalar w137 = -w136*w3;
    const scalar w138 = w137*w6;
    const scalar w139 = w138*w29;
    const scalar w140 = Msi1*w8;
    const scalar w141 = w140*w29;
    const scalar w142 = Msi0*w8;
    const scalar w143 = w137*w5;
    const scalar w144 = w143*w29;
    const scalar w145 = w137*w3;
    const scalar w146 = w29*w9;
    const scalar w147 = w110*w30 + w120*w46 + w124*w146 + w135*w27 + w139*w27 + w141*w44 + w142*w54 + w144*w44 + w145*w54;
    const scalar w148 = w120*w30 + w126*w146 + w128*w46 + w135*w44 + w139*w44 + w141*w58 + w142*w60 + w144*w58 + w145*w60;
    const scalar w149 = w134*w54;
    const scalar w150 = w140*w60;
    const scalar w151 = w124*w30 + w125*w146 + w126*w46 + w138*w54 + w142*w64 + w143*w60 + w145*w64 + w149 + w150;
    const scalar w152 = w78*w80;
    const scalar w153 = w102*w17;
    const scalar w154 = w14*y;
    const scalar w155 = theta1*w154;
    const scalar w156 = 4*w155;
    const scalar w157 = w108*w156;
    const scalar w158 = w157*w24 - w157*w26;
    const scalar w159 = w111*w44;
    const scalar w160 = 2*w155;
    const scalar w161 = w160*w43;
    const scalar w162 = w160*w39;
    const scalar w163 = 2*theta1*w14*w19*w25*w31*y - w160*w40;
    const scalar w164 = -w160*w34 - w160*w36;
    const scalar w165 = w117*w163 + w119*w164 + w161*w38 - w162*w42;
    const scalar w166 = -w160*w35 - w160*w50;
    const scalar w167 = -w122*w160 + w160*w47;
    const scalar w168 = w117*w166 + w119*w167 + w161*w49 - w162*w52;
    const scalar w169 = w38*(4*theta1*w14*w19*w25*w31*y - w156*w40) + w42*(-w156*w34 - w156*w36);
    const scalar w170 = w163*w49 + w164*w52 + w166*w38 + w167*w42;
    const scalar w171 = w49*(-w156*w35 - w156*w50) + w52*(-w122*w156 + w156*w47);
    const scalar w172 = w130*w44;
    const scalar w173 = Msi4*w8;
    const scalar w174 = w173*w29;
    const scalar w175 = -w136*w5;
    const scalar w176 = w175*w6;
    const scalar w177 = w176*w29;
    const scalar w178 = Msi3*w8;
    const scalar w179 = w178*w29;
    const scalar w180 = w175*w5;
    const scalar w181 = w180*w29;
    const scalar w182 = w175*w3;
    const scalar w183 = w140*w54 + w146*w168 + w158*w30 + w165*w46 + w174*w27 + w177*w27 + w179*w44 + w181*w44 + w182*w54;
    const scalar w184 = w174*w44;
    const scalar w185 = w146*w170 + w150 + w165*w30 + w169*w46 + w177*w44 + w179*w58 + w181*w58 + w182*w60 + w184;
    const scalar w186 = w140*w64 + w146*w171 + w168*w30 + w170*w46 + w173*w54 + w176*w54 + w178*w60 + w180*w60 + w182*w64;
    const scalar w187 = w14*z;
    const scalar w188 = theta1*w187;
    const scalar w189 = 4*w188;
    const scalar w190 = w108*w189;
    const scalar w191 = w190*w24 - w190*w26;
    const scalar w192 = 2*w188;
    const scalar w193 = w192*w43;
    const scalar w194 = w192*w39;
    const scalar w195 = 2*theta1*w14*w19*w25*w31*z - w192*w40;
    const scalar w196 = -w192*w34 - w192*w36;
    const scalar w197 = w117*w195 + w119*w196 + w193*w38 - w194*w42;
    const scalar w198 = -w192*w35 - w192*w50;
    const scalar w199 = -w122*w192 + w192*w47;
    const scalar w200 = w117*w198 + w119*w199 + w193*w49 - w194*w52;
    const scalar w201 = w38*(4*theta1*w14*w19*w25*w31*z - w189*w40) + w42*(-w189*w34 - w189*w36);
    const scalar w202 = w195*w49 + w196*w52 + w198*w38 + w199*w42;
    const scalar w203 = w49*(-w189*w35 - w189*w50) + w52*(-w122*w189 + w189*w47);
    const scalar w204 = Msi5*w8;
    const scalar w205 = w204*w29;
    const scalar w206 = -w136*w6;
    const scalar w207 = w206*w6;
    const scalar w208 = w207*w29;
    const scalar w209 = w206*w5;
    const scalar w210 = w209*w29;
    const scalar w211 = w206*w3;
    const scalar w212 = w146*w200 + w149 + w184 + w191*w30 + w197*w46 + w205*w27 + w208*w27 + w210*w44 + w211*w54;
    const scalar w213 = w134*w60 + w146*w202 + w174*w58 + w197*w30 + w201*w46 + w205*w44 + w208*w44 + w210*w58 + w211*w60;
    const scalar w214 = w134*w64 + w146*w203 + w173*w60 + w200*w30 + w202*w46 + w204*w54 + w207*w54 + w209*w60 + w211*w64;
    return vector(w100*x + w101*w103 + w153*(w104*w72 - w105*w90 + w129*(-w68*(w110*w69 + w112 + w120*w71 + w124*w73) - w70*(w120*w69 + w126*w73 + w127 + w128*w71) - w72*(w111*w63 + w124*w69 + w125*w73 + w126*w71) - w76) + w133*(w86*(-w110*w87 - w120*w89 - w124*w91 - w131) + w88*(-w120*w87 - w126*w91 - w128*w89 - w132) + w90*(-w124*w87 - w125*w91 - w126*w89 - w130*w63) + w94) + w152*(-w147*w68 - w148*w70 - w151*w72 - w66) + w98*(w147*w86 + w148*w88 + w151*w90 + w84)) + w9, w100*y + w103*w154 + w153*(w104*w70 - w105*w88 + w129*(-w68*(w158*w69 + w159 + w165*w71 + w168*w73) - w70*(w111*w58 + w165*w69 + w169*w71 + w170*w73) - w72*(w127 + w168*w69 + w170*w71 + w171*w73) - w75) + w133*(w86*(-w158*w87 - w165*w89 - w168*w91 - w172) + w88*(-w130*w58 - w165*w87 - w169*w89 - w170*w91) + w90*(-w132 - w168*w87 - w170*w89 - w171*w91) + w93) + w152*(-w183*w68 - w185*w70 - w186*w72 - w62) + w98*(w183*w86 + w185*w88 + w186*w90 + w83)) + w45, w100*z + w103*w187 + w153*(w104*w68 - w105*w86 + w129*(-w68*(w111*w27 + w191*w69 + w197*w71 + w200*w73) - w70*(w159 + w197*w69 + w201*w71 + w202*w73) - w72*(w112 + w200*w69 + w202*w71 + w203*w73) - w74) + w133*(w86*(-w130*w27 - w191*w87 - w197*w89 - w200*w91) + w88*(-w172 - w197*w87 - w201*w89 - w202*w91) + w90*(-w131 - w200*w87 - w202*w89 - w203*w91) + w92) + w152*(-w212*w68 - w213*w70 - w214*w72 - w56) + w98*(w212*w86 + w213*w88 + w214*w90 + w82)) + w28);
}

inline scalar source(const point& pt, const Params& p)
{
    const scalar r1 = p.r1;
    const scalar r2 = p.r2;
    const scalar gamma1 = p.gamma1;
    const scalar gamma2 = p.gamma2;
    const scalar gil = p.gil;
    const scalar git = p.git;
    const scalar gin = p.gin;
    const scalar theta0 = p.theta0;
    const scalar theta1 = p.theta1;
    const scalar chi = p.chi;
    const scalar cm = p.cm;
    const scalar c0 = p.c0;
    const scalar delta = p.delta;
    const scalar xs0 = p.xs[0];
    const scalar xs1 = p.xs[1];
    const scalar xs2 = p.xs[2];
    const scalar Msi0 = p.Msi[0];
    const scalar Msi1 = p.Msi[1];
    const scalar Msi2 = p.Msi[2];
    const scalar Msi3 = p.Msi[3];
    const scalar Msi4 = p.Msi[4];
    const scalar Msi5 = p.Msi[5];
    const scalar x = pt.x();
    const scalar y = pt.y();
    const scalar z = pt.z();

    const scalar w0 = sin(gamma2);
    const scalar w1 = (w0*w0);
    const scalar w2 = gil*w1;
    const scalar w3 = (r1*r1);
    const scalar w4 = -(r2*r2);
    const scalar w5 = -w3 - w4;
    const scalar w6 = (1.0/(w5));
    const scalar w7 = (x*x);
    const scalar w8 = (y*y);
    const scalar w9 = (z*z);
    const scalar w10 = w7 + w8 + w9;
    const scalar w11 = w6*(w10 - w3);
    const scalar w12 = theta0 + theta1*w11;
    const scalar w13 = cos(w12);
    const scalar w14 = sin(w12);
    const scalar w15 = theta1*w6;
    const scalar w16 = 4*w15;
    const scalar w17 = w13*w14*w16;
    const scalar w18 = w17*w2;
    const scalar w19 = git*w1;
    const scalar w20 = w17*w19;
    const scalar w21 = w18*z - w20*z;
    const scalar w22 = x - xs0;
    const scalar w23 = Msi2*w22;
    const scalar w24 = -xs1 + y;
    const scalar w25 = Msi4*w24;
    const scalar w26 = -xs2 + z;
    const scalar w27 = Msi5*w26;
    const scalar w28 = w23 + w25 + w27;
    const scalar w29 = (1.0/(c0));
    const scalar w30 = Msi0*w22;
    const scalar w31 = Msi1*w24;
    const scalar w32 = Msi2*w26;
    const scalar w33 = w30 + w31 + w32;
    const scalar w34 = Msi1*w22;
    const scalar w35 = Msi3*w24;
    const scalar w36 = Msi4*w26;
    const scalar w37 = w34 + w35 + w36;
    const scalar w38 = (delta*delta) + w22*w33 + w24*w37 + w26*w28;
    const scalar w39 = w29/sqrt(w38);
    const scalar w40 = w28*w39;
    const scalar w41 = r1 + r2;
    const scalar w42 = (1.0/(w41));
    const scalar w43 = w11*w42;
    const scalar w44 = (1.0/(r2));
    const scalar w45 = w44*z;
    const scalar w46 = cos(gamma2);
    const scalar w47 = (w14*w14);
    const scalar w48 = w2*w47;
    const scalar w49 = (w13*w13);
    const scalar w50 = w19*w49;
    const scalar w51 = gin*(w46*w46) + w48 + w50;
    const scalar w52 = 1/(chi*cm);
    const scalar w53 = w51*w52;
    const scalar w54 = w44*y;
    const scalar w55 = cos(gamma1);
    const scalar w56 = gin*w0*w46;
    const scalar w57 = sin(gamma1);
    const scalar w58 = w13*w57;
    const scalar w59 = w14*w55;
    const scalar w60 = w46*w59;
    const scalar w61 = w58 + w60;
    const scalar w62 = gil*w61;
    const scalar w63 = w0*w14;
    const scalar w64 = w62*w63;
    const scalar w65 = w14*w57;
    const scalar w66 = w13*w46*w55 - w65;
    const scalar w67 = git*w66;
    const scalar w68 = w0*w13;
    const scalar w69 = w67*w68;
    const scalar w70 = -w55*w56 + w64 + w69;
    const scalar w71 = w52*w70;
    const scalar w72 = w44*x;
    const scalar w73 = w46*w65;
    const scalar w74 = w13*w55 - w73;
    const scalar w75 = gil*w74;
    const scalar w76 = w63*w75;
    const scalar w77 = w46*w58;
    const scalar w78 = -w59 - w77;
    const scalar w79 = git*w78;
    const scalar w80 = w68*w79;
    const scalar w81 = w56*w57 + w76 + w80;
    const scalar w82 = w52*w81;
    const scalar w83 = w44*(w45*w53 + w54*w71 + w72*w82);
    const scalar w84 = gin*w1;
    const scalar w85 = gil*(w61*w61) + git*(w66*w66) + (w55*w55)*w84;
    const scalar w86 = w52*w85;
    const scalar w87 = -w55*w57*w84 + w61*w75 + w66*w79;
    const scalar w88 = w52*w87;
    const scalar w89 = w44*(w45*w71 + w54*w86 + w72*w88);
    const scalar w90 = gil*(w74*w74) + git*(w78*w78) + (w57*w57)*w84;
    const scalar w91 = w52*w90;
    const scalar w92 = w44*(w45*w82 + w54*w88 + w72*w91);
    const scalar w93 = w83*z + w89*y + w92*x;
    const scalar w94 = (1.0/(w93));
    const scalar w95 = (1.0/2.0)*w94;
    const scalar w96 = w37*w39;
    const scalar w97 = w33*w39;
    const scalar w98 = w40*w53 + w71*w96 + w82*w97;
    const scalar w99 = w44*w98;
    const scalar w100 = w40*w71 + w86*w96 + w88*w97;
    const scalar w101 = w100*w44;
    const scalar w102 = w40*w82 + w88*w96 + w91*w97;
    const scalar w103 = w102*w44;
    const scalar w104 = -w101*y - w103*x - w99*z;
    const scalar w105 = w104*w44;
    const scalar w106 = w105*w11*w41;
    const scalar w107 = (1.0/(r1));
    const scalar w108 = w107*w98;
    const scalar w109 = w100*w107;
    const scalar w110 = w102*w107;
    const scalar w111 = w108*z + w109*y + w110*x;
    const scalar w112 = w41*(1 - w11);
    const scalar w113 = w107*z;
    const scalar w114 = w107*y;
    const scalar w115 = w107*x;
    const scalar w116 = w107*(-w113*w53 - w114*w71 - w115*w82);
    const scalar w117 = w107*(-w113*w71 - w114*w86 - w115*w88);
    const scalar w118 = w107*(-w113*w82 - w114*w88 - w115*w91);
    const scalar w119 = -w116*z - w117*y - w118*x;
    const scalar w120 = (1.0/(w119));
    const scalar w121 = w107*w120;
    const scalar w122 = (1.0/2.0)*w121;
    const scalar w123 = w112*w122;
    const scalar w124 = w106*w95 + w111*w123;
    const scalar w125 = 2*w124;
    const scalar w126 = w125*w43;
    const scalar w127 = w10 + w4;
    const scalar w128 = w42*w6;
    const scalar w129 = w127*w128;
    const scalar w130 = w125*w129;
    const scalar w131 = w41*w6;
    const scalar w132 = w131*w94;
    const scalar w133 = w104*w132;
    const scalar w134 = w111*w131;
    const scalar w135 = w120*w134;
    const scalar w136 = (1.0/(w93*w93));
    const scalar w137 = w44*w52;
    const scalar w138 = w45*w52;
    const scalar w139 = 2*w15;
    const scalar w140 = w139*w68;
    const scalar w141 = w140*w62;
    const scalar w142 = w139*w63;
    const scalar w143 = w142*w67;
    const scalar w144 = w139*w65;
    const scalar w145 = 2*theta1*w13*w46*w55*w6*z - w144*z;
    const scalar w146 = gil*w145;
    const scalar w147 = w139*w58;
    const scalar w148 = w139*w60;
    const scalar w149 = -w147*z - w148*z;
    const scalar w150 = git*w149;
    const scalar w151 = w141*z - w143*z + w146*w63 + w150*w68;
    const scalar w152 = w52*w54;
    const scalar w153 = w140*w75;
    const scalar w154 = w142*w79;
    const scalar w155 = w139*w59;
    const scalar w156 = w139*w77;
    const scalar w157 = -w155*z - w156*z;
    const scalar w158 = gil*w63;
    const scalar w159 = w13*w55;
    const scalar w160 = w15*w159;
    const scalar w161 = 2*w160;
    const scalar w162 = w139*w73;
    const scalar w163 = -w161*z + w162*z;
    const scalar w164 = git*w68;
    const scalar w165 = w153*z - w154*z + w157*w158 + w163*w164;
    const scalar w166 = w52*w72;
    const scalar w167 = w137*w51 + w138*w21 + w151*w152 + w165*w166;
    const scalar w168 = w167*w45;
    const scalar w169 = w137*w70;
    const scalar w170 = w16*w65;
    const scalar w171 = 4*theta1*w13*w46*w55*w6*z - w170*z;
    const scalar w172 = w16*w58;
    const scalar w173 = w16*w60;
    const scalar w174 = -w172*z - w173*z;
    const scalar w175 = w171*w62 + w174*w67;
    const scalar w176 = w145*w75 + w149*w79 + w157*w62 + w163*w67;
    const scalar w177 = w44*(w138*w151 + w152*w175 + w166*w176 + w169);
    const scalar w178 = w177*y;
    const scalar w179 = w137*w81;
    const scalar w180 = 4*w59;
    const scalar w181 = w15*w180;
    const scalar w182 = w16*w77;
    const scalar w183 = -w181*z - w182*z;
    const scalar w184 = 4*w160;
    const scalar w185 = w16*w73;
    const scalar w186 = -w184*z + w185*z;
    const scalar w187 = w183*w75 + w186*w79;
    const scalar w188 = w44*(w138*w165 + w152*w176 + w166*w187 + w179);
    const scalar w189 = w188*x;
    const scalar w190 = -w168 - w178 - w189 - w83;
    const scalar w191 = w136*w190;
    const scalar w192 = (1.0/2.0)*w106;
    const scalar w193 = w107*w52;
    const scalar w194 = w113*w52;
    const scalar w195 = w114*w52;
    const scalar w196 = w115*w52;
    const scalar w197 = -w151*w195 - w165*w196 - w193*w51 - w194*w21;
    const scalar w198 = w113*w197;
    const scalar w199 = w193*w70;
    const scalar w200 = w107*(-w151*w194 - w175*w195 - w176*w196 - w199);
    const scalar w201 = w200*y;
    const scalar w202 = w193*w81;
    const scalar w203 = w107*(-w165*w194 - w176*w195 - w187*w196 - w202);
    const scalar w204 = w203*x;
    const scalar w205 = w116 + w198 + w201 + w204;
    const scalar w206 = (1.0/(w119*w119));
    const scalar w207 = w107*w206;
    const scalar w208 = w112*w207;
    const scalar w209 = w205*w208;
    const scalar w210 = (1.0/2.0)*w111;
    const scalar w211 = Msi5*w39;
    const scalar w212 = -w28;
    const scalar w213 = w29/pow(w38, 1.5);
    const scalar w214 = w213*w28;
    const scalar w215 = w212*w214;
    const scalar w216 = w40*w52;
    const scalar w217 = Msi4*w39;
    const scalar w218 = w217*w71;
    const scalar w219 = Msi2*w39;
    const scalar w220 = w219*w82;
    const scalar w221 = w213*w37;
    const scalar w222 = w212*w221;
    const scalar w223 = w213*w33;
    const scalar w224 = w212*w223;
    const scalar w225 = w52*w96;
    const scalar w226 = w52*w97;
    const scalar w227 = w151*w225 + w165*w226 + w21*w216 + w211*w53 + w215*w53 + w218 + w220 + w222*w71 + w224*w82;
    const scalar w228 = w151*w216 + w175*w225 + w176*w226 + w211*w71 + w215*w71 + w217*w86 + w219*w88 + w222*w86 + w224*w88;
    const scalar w229 = w228*w44;
    const scalar w230 = w165*w216 + w176*w225 + w187*w226 + w211*w82 + w215*w82 + w217*w88 + w219*w91 + w222*w88 + w224*w91;
    const scalar w231 = w230*w44;
    const scalar w232 = -w227*w45 - w229*y - w231*x - w99;
    const scalar w233 = w11*w41*w44;
    const scalar w234 = w233*w95;
    const scalar w235 = w107*w228;
    const scalar w236 = w107*w230;
    const scalar w237 = w108 + w113*w227 + w235*y + w236*x;
    const scalar w238 = -w113*w135 + w123*w237 + w133*w45 + w191*w192 + w209*w210 + w232*w234;
    const scalar w239 = w127*w43;
    const scalar w240 = w126*z + w130*z + w238*w239 + w40;
    const scalar w241 = w240*w52;
    const scalar w242 = 2*theta1*w13*w46*w55*w6*y - w144*y;
    const scalar w243 = gil*w242;
    const scalar w244 = -w147*y - w148*y;
    const scalar w245 = git*w244;
    const scalar w246 = w141*y - w143*y + w243*w63 + w245*w68;
    const scalar w247 = -w155*x - w156*x;
    const scalar w248 = -w161*x + w162*x;
    const scalar w249 = w153*x - w154*x + w158*w247 + w164*w248;
    const scalar w250 = 4*theta1*w13*w46*w55*w6*y - w170*y;
    const scalar w251 = -w172*y - w173*y;
    const scalar w252 = w250*w62 + w251*w67;
    const scalar w253 = w18*y - w20*y;
    const scalar w254 = -w155*y - w156*y;
    const scalar w255 = -w161*y + w162*y;
    const scalar w256 = w153*y - w154*y + w158*w254 + w164*w255;
    const scalar w257 = w44*(w138*w253 + w152*w246 + w166*w256 + w169);
    const scalar w258 = w257*z;
    const scalar w259 = w242*w75 + w244*w79 + w254*w62 + w255*w67;
    const scalar w260 = w137*w85 + w138*w246 + w152*w252 + w166*w259;
    const scalar w261 = w260*w54;
    const scalar w262 = w137*w87;
    const scalar w263 = -w181*y - w182*y;
    const scalar w264 = -w184*y + w185*y;
    const scalar w265 = w263*w75 + w264*w79;
    const scalar w266 = w44*(w138*w256 + w152*w259 + w166*w265 + w262);
    const scalar w267 = w266*x;
    const scalar w268 = -w258 - w261 - w267 - w89;
    const scalar w269 = w136*w268;
    const scalar w270 = w107*(-w194*w253 - w195*w246 - w196*w256 - w199);
    const scalar w271 = w270*z;
    const scalar w272 = -w193*w85 - w194*w246 - w195*w252 - w196*w259;
    const scalar w273 = w114*w272;
    const scalar w274 = w193*w87;
    const scalar w275 = w107*(-w194*w256 - w195*w259 - w196*w265 - w274);
    const scalar w276 = w275*x;
    const scalar w277 = w117 + w271 + w273 + w276;
    const scalar w278 = w208*w277;
    const scalar w279 = -w37;
    const scalar w280 = w214*w279;
    const scalar w281 = Msi3*w39;
    const scalar w282 = Msi1*w39;
    const scalar w283 = w221*w279;
    const scalar w284 = w223*w279;
    const scalar w285 = w216*w253 + w217*w53 + w225*w246 + w226*w256 + w280*w53 + w281*w71 + w282*w82 + w283*w71 + w284*w82;
    const scalar w286 = w285*w44;
    const scalar w287 = w282*w88;
    const scalar w288 = w216*w246 + w218 + w225*w252 + w226*w259 + w280*w71 + w281*w86 + w283*w86 + w284*w88 + w287;
    const scalar w289 = w216*w256 + w217*w82 + w225*w259 + w226*w265 + w280*w82 + w281*w88 + w282*w91 + w283*w88 + w284*w91;
    const scalar w290 = w289*w44;
    const scalar w291 = -w101 - w286*z - w288*w54 - w290*x;
    const scalar w292 = w107*w285;
    const scalar w293 = w107*w289;
    const scalar w294 = w109 + w114*w288 + w292*z + w293*x;
    const scalar w295 = -w114*w135 + w123*w294 + w133*w54 + w192*w269 + w210*w278 + w234*w291;
    const scalar w296 = w126*y + w130*y + w239*w295 + w96;
    const scalar w297 = w296*w52;
    const scalar w298 = -w181*x - w182*x;
    const scalar w299 = -w184*x + w185*x;
    const scalar w300 = w298*w75 + w299*w79;
    const scalar w301 = w18*x - w20*x;
    const scalar w302 = 2*theta1*w13*w46*w55*w6*x - w144*x;
    const scalar w303 = gil*w302;
    const scalar w304 = -w147*x - w148*x;
    const scalar w305 = git*w304;
    const scalar w306 = w141*x - w143*x + w303*w63 + w305*w68;
    const scalar w307 = w44*(w138*w301 + w152*w306 + w166*w249 + w179);
    const scalar w308 = w307*z;
    const scalar w309 = w247*w62 + w248*w67 + w302*w75 + w304*w79;
    const scalar w310 = w137*w90 + w138*w249 + w152*w309 + w166*w300;
    const scalar w311 = w310*w72;
    const scalar w312 = 4*theta1*w13*w46*w55*w6*x - w170*x;
    const scalar w313 = -w172*x - w173*x;
    const scalar w314 = w312*w62 + w313*w67;
    const scalar w315 = w44*(w138*w306 + w152*w314 + w166*w309 + w262);
    const scalar w316 = w315*y;
    const scalar w317 = -w308 - w311 - w316 - w92;
    const scalar w318 = w136*w317;
    const scalar w319 = w107*(-w194*w301 - w195*w306 - w196*w249 - w202);
    const scalar w320 = w319*z;
    const scalar w321 = -w193*w90 - w194*w249 - w195*w309 - w196*w300;
    const scalar w322 = w115*w321;
    const scalar w323 = w107*(-w194*w306 - w195*w314 - w196*w309 - w274);
    const scalar w324 = w323*y;
    const scalar w325 = w118 + w320 + w322 + w324;
    const scalar w326 = w207*w325;
    const scalar w327 = w112*w326;
    const scalar w328 = -w33;
    const scalar w329 = w214*w328;
    const scalar w330 = Msi0*w39;
    const scalar w331 = w221*w328;
    const scalar w332 = w223*w328;
    const scalar w333 = w216*w301 + w219*w53 + w225*w306 + w226*w249 + w282*w71 + w329*w53 + w330*w82 + w331*w71 + w332*w82;
    const scalar w334 = w333*w44;
    const scalar w335 = w216*w306 + w219*w71 + w225*w314 + w226*w309 + w282*w86 + w329*w71 + w330*w88 + w331*w86 + w332*w88;
    const scalar w336 = w335*w44;
    const scalar w337 = w216*w249 + w220 + w225*w309 + w226*w300 + w287 + w329*w82 + w330*w91 + w331*w88 + w332*w91;
    const scalar w338 = -w103 - w334*z - w336*y - w337*w72;
    const scalar w339 = w233*w338;
    const scalar w340 = w107*w333;
    const scalar w341 = w107*w335;
    const scalar w342 = w110 + w115*w337 + w340*z + w341*y;
    const scalar w343 = w112*w342;
    const scalar w344 = -w115*w135 + w122*w343 + w133*w72 + w192*w318 + w210*w327 + w339*w95;
    const scalar w345 = w126*x + w130*x + w239*w344 + w97;
    const scalar w346 = w345*w52;
    const scalar w347 = w124*w128;
    const scalar w348 = 8*w347;
    const scalar w349 = 2*w107;
    const scalar w350 = Msi5*w214;
    const scalar w351 = 2*w53;
    const scalar w352 = w212*w213;
    const scalar w353 = Msi5*w352;
    const scalar w354 = w29/pow(w38, 2.5);
    const scalar w355 = -3*w23 - 3*w25 - 3*w27;
    const scalar w356 = w212*w355;
    const scalar w357 = w354*w356;
    const scalar w358 = w28*w357;
    const scalar w359 = w211*w52;
    const scalar w360 = 2*w359;
    const scalar w361 = w215*w52;
    const scalar w362 = 2*w361;
    const scalar w363 = Msi5*w221;
    const scalar w364 = Msi5*w223;
    const scalar w365 = Msi4*w352;
    const scalar w366 = w365*w71;
    const scalar w367 = w213*w82;
    const scalar w368 = w212*w367;
    const scalar w369 = Msi2*w368;
    const scalar w370 = w357*w37;
    const scalar w371 = w33*w354;
    const scalar w372 = w356*w371;
    const scalar w373 = (theta1*theta1);
    const scalar w374 = (1.0/(w5*w5));
    const scalar w375 = w373*w374;
    const scalar w376 = w375*w9;
    const scalar w377 = 8*w376;
    const scalar w378 = -4*gil*theta1*w1*w13*w14*w6 + w20;
    const scalar w379 = 8*gil*w1*w373*w374*w49*w9 + 8*git*w1*w373*w374*w47*w9 - w377*w48 - w377*w50 - w378;
    const scalar w380 = w217*w52;
    const scalar w381 = w151*w380;
    const scalar w382 = w219*w52;
    const scalar w383 = w165*w382;
    const scalar w384 = w222*w52;
    const scalar w385 = 2*w384;
    const scalar w386 = w224*w52;
    const scalar w387 = 2*w386;
    const scalar w388 = 4*w376;
    const scalar w389 = -2*theta1*w13*w46*w55*w6 + w144;
    const scalar w390 = -w388*w58 - w388*w60 - w389;
    const scalar w391 = w159*w46;
    const scalar w392 = w147 + w148;
    const scalar w393 = 4*w14*w373*w374*w57*w9 - w388*w391 - w392;
    const scalar w394 = w16*z;
    const scalar w395 = -2*gil*theta1*w0*w13*w6*w61 + w143;
    const scalar w396 = 4*gil*theta1*w0*w13*w145*w6*z + gil*w0*w14*w390 + git*w0*w13*w393 - w150*w394*w63 - w388*w64 - w388*w69 - w395;
    const scalar w397 = w155 + w156;
    const scalar w398 = 4*w14*w373*w374*w46*w57*w9 - w159*w388 - w397;
    const scalar w399 = -w161 + w162;
    const scalar w400 = w180*w376 + w388*w77 + w399;
    const scalar w401 = git*w63;
    const scalar w402 = -2*gil*theta1*w0*w13*w6*w74 + w154;
    const scalar w403 = 4*gil*theta1*w0*w13*w157*w6*z + gil*w0*w14*w398 + git*w0*w13*w400 - w163*w394*w401 - w388*w76 - w388*w80 - w402;
    const scalar w404 = w151*w385 + w165*w387 + w21*w360 + w21*w362 + w216*w379 + w225*w396 + w226*w403 - w350*w53 + w351*w353 + w358*w53 - w363*w71 - w364*w82 + 2*w366 + 2*w369 + w370*w71 + w372*w82 + 2*w381 + 2*w383;
    const scalar w405 = Msi4*w86;
    const scalar w406 = 2*w71;
    const scalar w407 = w213*w88;
    const scalar w408 = 2*w212;
    const scalar w409 = w407*w408;
    const scalar w410 = 2*w380;
    const scalar w411 = 2*w382;
    const scalar w412 = -4*theta1*w13*w46*w55*w6 + w170;
    const scalar w413 = w172 + w173;
    const scalar w414 = w146*w171 + w150*w174 + w62*(-w377*w58 - w377*w60 - w412) + w67*(8*w14*w373*w374*w57*w9 - w377*w391 - w413);
    const scalar w415 = 2*w163;
    const scalar w416 = 2*w146*w157 + w150*w415 + w390*w75 + w393*w79 + w398*w62 + w400*w67;
    const scalar w417 = Msi2*w409 + w151*w360 + w151*w362 + w175*w385 + w175*w410 + w176*w387 + w176*w411 + w216*w396 + w225*w414 + w226*w416 - w350*w71 + 2*w352*w405 + w353*w406 + w358*w71 - w363*w86 - w364*w88 + w370*w86 + w372*w88;
    const scalar w418 = w213*w91;
    const scalar w419 = Msi2*w418;
    const scalar w420 = w371*w91;
    const scalar w421 = gil*w183;
    const scalar w422 = git*w186;
    const scalar w423 = w181 + w182;
    const scalar w424 = -w184 + w185;
    const scalar w425 = w157*w421 + w163*w422 + w75*(8*w14*w373*w374*w46*w57*w9 - w159*w377 - w423) + w79*(w377*w59 + w377*w77 + w424);
    const scalar w426 = Msi4*w409 + 2*Msi5*w368 + w165*w360 + w165*w362 + w176*w385 + w176*w410 + w187*w387 + w187*w411 + w216*w403 + w225*w416 + w226*w425 - w350*w82 + w356*w420 + w358*w82 - w363*w88 - w364*w91 + w370*w88 + w408*w419;
    const scalar w427 = w349*w52;
    const scalar w428 = w208*w210;
    const scalar w429 = 2*w44;
    const scalar w430 = w120*w131;
    const scalar w431 = w113*w430;
    const scalar w432 = w132*w45;
    const scalar w433 = w136*w233;
    const scalar w434 = w190*w433;
    const scalar w435 = 2*w116 + 2*w198 + 2*w201 + 2*w204;
    const scalar w436 = w107*w112*w210/(w119*w119*w119);
    const scalar w437 = w205*w436;
    const scalar w438 = w429*w52;
    const scalar w439 = w136*w192;
    const scalar w440 = w134*w206;
    const scalar w441 = w113*w440;
    const scalar w442 = w104*w131;
    const scalar w443 = w442*w45;
    const scalar w444 = -2*w168 - 2*w178 - 2*w189 - 2*w83;
    const scalar w445 = w192/(w93*w93*w93);
    const scalar w446 = w190*w445;
    const scalar w447 = w105*w132 - w121*w134;
    const scalar w448 = 4*w43;
    const scalar w449 = w238*z;
    const scalar w450 = 4*w129;
    const scalar w451 = w126 + w130;
    const scalar w452 = Msi3*w214;
    const scalar w453 = w213*w279;
    const scalar w454 = -3*w34 - 3*w35 - 3*w36;
    const scalar w455 = w279*w454;
    const scalar w456 = w354*w455;
    const scalar w457 = w28*w456;
    const scalar w458 = w280*w52;
    const scalar w459 = 2*w458;
    const scalar w460 = Msi3*w221;
    const scalar w461 = Msi3*w223;
    const scalar w462 = Msi3*w453;
    const scalar w463 = 2*w453;
    const scalar w464 = w463*w82;
    const scalar w465 = w37*w456;
    const scalar w466 = w371*w455;
    const scalar w467 = w375*w8;
    const scalar w468 = 8*w467;
    const scalar w469 = 8*gil*w1*w373*w374*w49*w8 + 8*git*w1*w373*w374*w47*w8 - w378 - w468*w48 - w468*w50;
    const scalar w470 = w281*w52;
    const scalar w471 = 2*w470;
    const scalar w472 = w282*w52;
    const scalar w473 = 2*w472;
    const scalar w474 = w283*w52;
    const scalar w475 = 2*w474;
    const scalar w476 = w284*w52;
    const scalar w477 = 2*w476;
    const scalar w478 = 4*w467;
    const scalar w479 = -w389 - w478*w58 - w478*w60;
    const scalar w480 = 4*w14*w373*w374*w57*w8 - w391*w478 - w392;
    const scalar w481 = w16*y;
    const scalar w482 = 4*gil*theta1*w0*w13*w242*w6*y + gil*w0*w14*w479 + git*w0*w13*w480 - w245*w481*w63 - w395 - w478*w64 - w478*w69;
    const scalar w483 = 4*w14*w373*w374*w46*w57*w8 - w159*w478 - w397;
    const scalar w484 = w180*w467 + w399 + w478*w77;
    const scalar w485 = 4*gil*theta1*w0*w13*w254*w6*y + gil*w0*w14*w483 + git*w0*w13*w484 - w255*w401*w481 - w402 - w478*w76 - w478*w80;
    const scalar w486 = Msi1*w464 + Msi4*w351*w453 + w216*w469 + w225*w482 + w226*w485 + w246*w471 + w246*w475 + w253*w410 + w253*w459 + w256*w473 + w256*w477 + w406*w462 - w452*w53 + w457*w53 - w460*w71 - w461*w82 + w465*w71 + w466*w82;
    const scalar w487 = 2*w86;
    const scalar w488 = Msi4*w71;
    const scalar w489 = w453*w488;
    const scalar w490 = Msi1*w88;
    const scalar w491 = w453*w490;
    const scalar w492 = w246*w380;
    const scalar w493 = w259*w472;
    const scalar w494 = w243*w250 + w245*w251 + w62*(-w412 - w468*w58 - w468*w60) + w67*(8*w14*w373*w374*w57*w8 - w391*w468 - w413);
    const scalar w495 = 2*w255;
    const scalar w496 = 2*w243*w254 + w245*w495 + w479*w75 + w480*w79 + w483*w62 + w484*w67;
    const scalar w497 = w216*w482 + w225*w494 + w226*w496 + w246*w459 + w252*w471 + w252*w475 + w259*w477 - w452*w71 + w457*w71 - w460*w86 - w461*w88 + w462*w487 + w465*w86 + w466*w88 + 2*w489 + 2*w491 + 2*w492 + 2*w493;
    const scalar w498 = Msi1*w418;
    const scalar w499 = Msi3*w88;
    const scalar w500 = gil*w263;
    const scalar w501 = git*w264;
    const scalar w502 = w254*w500 + w255*w501 + w75*(8*w14*w373*w374*w46*w57*w8 - w159*w468 - w423) + w79*(w424 + w468*w59 + w468*w77);
    const scalar w503 = Msi4*w464 + w216*w485 + w225*w496 + w226*w502 + w256*w410 + w256*w459 + w259*w471 + w259*w475 + w265*w473 + w265*w477 + 2*w279*w498 + w420*w455 - w452*w82 + w457*w82 - w460*w88 - w461*w91 + w463*w499 + w465*w88;
    const scalar w504 = w114*w430;
    const scalar w505 = w132*w54;
    const scalar w506 = w268*w433;
    const scalar w507 = 2*w117 + 2*w271 + 2*w273 + 2*w276;
    const scalar w508 = w277*w436;
    const scalar w509 = w114*w440;
    const scalar w510 = w442*w54;
    const scalar w511 = -2*w258 - 2*w261 - 2*w267 - 2*w89;
    const scalar w512 = w268*w445;
    const scalar w513 = w295*y;
    const scalar w514 = 8*w7;
    const scalar w515 = Msi0*w214;
    const scalar w516 = w213*w328;
    const scalar w517 = Msi2*w516;
    const scalar w518 = -3*w30 - 3*w31 - 3*w32;
    const scalar w519 = w328*w518;
    const scalar w520 = w354*w519;
    const scalar w521 = w28*w520;
    const scalar w522 = w329*w52;
    const scalar w523 = 2*w522;
    const scalar w524 = Msi0*w221;
    const scalar w525 = Msi0*w223;
    const scalar w526 = Msi1*w328;
    const scalar w527 = w213*w526;
    const scalar w528 = 2*Msi0;
    const scalar w529 = w516*w528;
    const scalar w530 = w37*w520;
    const scalar w531 = w371*w519;
    const scalar w532 = w375*w514;
    const scalar w533 = 8*gil*w1*w373*w374*w49*w7 + 8*git*w1*w373*w374*w47*w7 - w378 - w48*w532 - w50*w532;
    const scalar w534 = w330*w52;
    const scalar w535 = 2*w534;
    const scalar w536 = w331*w52;
    const scalar w537 = 2*w536;
    const scalar w538 = w332*w52;
    const scalar w539 = 2*w538;
    const scalar w540 = w375*w7;
    const scalar w541 = 4*w540;
    const scalar w542 = -w389 - w541*w58 - w541*w60;
    const scalar w543 = 4*w14*w373*w374*w57*w7 - w391*w541 - w392;
    const scalar w544 = w16*x;
    const scalar w545 = 4*gil*theta1*w0*w13*w302*w6*x + gil*w0*w14*w542 + git*w0*w13*w543 - w305*w544*w63 - w395 - w541*w64 - w541*w69;
    const scalar w546 = 4*w14*w373*w374*w46*w57*w7 - w159*w541 - w397;
    const scalar w547 = w180*w540 + w399 + w541*w77;
    const scalar w548 = 4*gil*theta1*w0*w13*w247*w6*x + gil*w0*w14*w546 + git*w0*w13*w547 - w248*w401*w544 - w402 - w541*w76 - w541*w80;
    const scalar w549 = w216*w533 + w225*w545 + w226*w548 + w249*w535 + w249*w539 + w301*w411 + w301*w523 + w306*w473 + w306*w537 + w351*w517 + w406*w527 - w515*w53 + w521*w53 - w524*w71 - w525*w82 + w529*w82 + w530*w71 + w531*w82;
    const scalar w550 = w303*w312 + w305*w313 + w62*(-w412 - w532*w58 - w532*w60) + w67*(8*w14*w373*w374*w57*w7 - w391*w532 - w413);
    const scalar w551 = 2*w248;
    const scalar w552 = 2*w247*w303 + w305*w551 + w542*w75 + w543*w79 + w546*w62 + w547*w67;
    const scalar w553 = w216*w545 + w225*w550 + w226*w552 + w306*w411 + w306*w523 + w309*w535 + w309*w539 + w314*w473 + w314*w537 + w406*w517 + w487*w527 - w515*w71 + w521*w71 - w524*w86 - w525*w88 + w529*w88 + w530*w86 + w531*w88;
    const scalar w554 = w517*w82;
    const scalar w555 = w407*w526;
    const scalar w556 = w249*w382;
    const scalar w557 = w309*w472;
    const scalar w558 = gil*w298;
    const scalar w559 = git*w299;
    const scalar w560 = w247*w558 + w248*w559 + w75*(8*w14*w373*w374*w46*w57*w7 - w159*w532 - w423) + w79*(w424 + w532*w59 + w532*w77);
    const scalar w561 = w216*w548 + w225*w552 + w226*w560 + w249*w523 + w300*w535 + w300*w539 + w309*w537 + w328*w418*w528 + w420*w519 - w515*w82 + w521*w82 - w524*w88 - w525*w91 + w530*w88 + 2*w554 + 2*w555 + 2*w556 + 2*w557;
    const scalar w562 = w115*w430;
    const scalar w563 = w132*w72;
    const scalar w564 = 2*w118 + 2*w320 + 2*w322 + 2*w324;
    const scalar w565 = w325*w436;
    const scalar w566 = w115*w440;
    const scalar w567 = w442*w72;
    const scalar w568 = -2*w308 - 2*w311 - 2*w316 - 2*w92;
    const scalar w569 = w317*w445;
    const scalar w570 = w344*x;
    const scalar w571 = 8*w159;
    const scalar w572 = y*z;
    const scalar w573 = w375*w572;
    const scalar w574 = 8*w573;
    const scalar w575 = w75*(-w571*w573 + w574*w73) + w79*(w574*w59 + w574*w77);
    const scalar w576 = w254*w421 + w255*w422 + w575;
    const scalar w577 = w212*w454;
    const scalar w578 = w354*w577;
    const scalar w579 = w37*w578;
    const scalar w580 = w28*w578;
    const scalar w581 = 4*w573;
    const scalar w582 = -w159*w581 + w581*w73;
    const scalar w583 = -w58*w581 - w581*w60;
    const scalar w584 = w180*w573 + w581*w77;
    const scalar w585 = -w391*w581 + w581*w65;
    const scalar w586 = w146*w254 + w150*w255 + w157*w243 + w163*w245 + w582*w62 + w583*w75 + w584*w67 + w585*w79;
    const scalar w587 = w15*w401;
    const scalar w588 = 2*gil*theta1*w0*w13*w157*w6*y + 2*gil*theta1*w0*w13*w254*w6*z + gil*w0*w14*w582 + git*w0*w13*w584 - w415*w587*y - w495*w587*z - w581*w76 - w581*w80;
    const scalar w589 = Msi4*w88;
    const scalar w590 = Msi5*w453;
    const scalar w591 = Msi4*w223;
    const scalar w592 = Msi4*w214;
    const scalar w593 = Msi4*w368 + w165*w380 + w165*w458 + w176*w470 + w176*w474 + w187*w472 + w187*w476 + w212*w498 + w216*w588 - w221*w589 + w225*w586 + w256*w359 + w256*w361 + w259*w380 + w259*w384 + w265*w382 + w265*w386 + w279*w419 + w352*w499 + w453*w589 + w590*w82 - w591*w91 - w592*w82;
    const scalar w594 = w226*w576 + w420*w577 + w579*w88 + w580*w82 + w593;
    const scalar w595 = w62*(-w574*w58 - w574*w60) + w67*(-w391*w574 + w574*w65);
    const scalar w596 = w171*w243 + w174*w245 + w595;
    const scalar w597 = w371*w577;
    const scalar w598 = 2*gil*theta1*w0*w13*w145*w6*y + 2*gil*theta1*w0*w13*w242*w6*z + gil*w0*w14*w583 + git*w0*w13*w585 - w142*w150*y - w142*w245*z - w581*w64 - w581*w69;
    const scalar w599 = Msi2*w453;
    const scalar w600 = Msi3*w352;
    const scalar w601 = w176*w472 + w352*w490;
    const scalar w602 = w151*w458 + w175*w470 + w175*w474 + w176*w476 + w216*w598 - w221*w405 - w223*w589 + w226*w586 + w246*w359 + w246*w361 + w252*w380 + w252*w384 + w259*w382 + w259*w386 + w366 + w381 + w405*w453 + w590*w71 - w592*w71 + w599*w88 + w600*w86 + w601;
    const scalar w603 = w225*w596 + w579*w86 + w580*w71 + w597*w88 + w602;
    const scalar w604 = 8*gil*w1*w373*w374*w49*y*z + 8*git*w1*w373*w374*w47*y*z - w48*w574 - w50*w574;
    const scalar w605 = w256*w382 + w599*w82;
    const scalar w606 = Msi1*w368 + w151*w470 + w151*w474 + w165*w472 + w165*w476 + w21*w380 + w21*w458 + w216*w604 - w221*w488 + w225*w598 + w226*w588 + w246*w384 + w253*w359 + w253*w361 + w256*w386 + w365*w53 + w489 + w492 + w53*w590 - w53*w592 - w591*w82 + w600*w71 + w605;
    const scalar w607 = w53*w580 + w579*w71 + w597*w82 + w606;
    const scalar w608 = w235 + w292;
    const scalar w609 = w193*w256;
    const scalar w610 = w176*w193;
    const scalar w611 = w194*w588 + w195*w586 + w609 + w610;
    const scalar w612 = w175*w193 + w193*w246 + w194*w598 + w196*w586;
    const scalar w613 = w113*(-w151*w193 - w193*w253 - w194*w604 - w195*w598 - w196*w588) + w200 + w270;
    const scalar w614 = w229 + w286;
    const scalar w615 = w137*w256;
    const scalar w616 = w137*w176;
    const scalar w617 = w138*w588 + w152*w586 + w615 + w616;
    const scalar w618 = w137*w175 + w137*w246 + w138*w598 + w166*w586;
    const scalar w619 = w177 + w257 + w45*(w137*w151 + w137*w253 + w138*w604 + w152*w598 + w166*w588);
    const scalar w620 = (1.0/2.0)*w232;
    const scalar w621 = (1.0/2.0)*w291;
    const scalar w622 = w191*w510 - w205*w509 + (1.0/2.0)*w209*w294 + w232*w505 + (1.0/2.0)*w237*w278 - w237*w504 + w269*w443 - w277*w441 + w291*w432 - w294*w431 + w434*w621 + w506*w620;
    const scalar w623 = 2*w43;
    const scalar w624 = w238*y;
    const scalar w625 = 2*w129;
    const scalar w626 = w295*z;
    const scalar w627 = w217 + w348*w572 + w623*w624 + w623*w626 + w624*w625 + w625*w626;
    const scalar w628 = w157*w500 + w163*w501 + w575;
    const scalar w629 = w279*w355;
    const scalar w630 = w354*w629;
    const scalar w631 = w37*w630;
    const scalar w632 = w28*w630;
    const scalar w633 = w226*w628 + w420*w629 + w593 + w631*w88 + w632*w82;
    const scalar w634 = w146*w250 + w150*w251 + w595;
    const scalar w635 = w371*w629;
    const scalar w636 = w225*w634 + w602 + w631*w86 + w632*w71 + w635*w88;
    const scalar w637 = w53*w632 + w606 + w631*w71 + w635*w82;
    const scalar w638 = w375*x;
    const scalar w639 = w638*z;
    const scalar w640 = 8*w639;
    const scalar w641 = w75*(-w571*w639 + w640*w73) + w79*(w59*w640 + w640*w77);
    const scalar w642 = w247*w421 + w248*w422 + w641;
    const scalar w643 = w212*w518;
    const scalar w644 = w354*w643;
    const scalar w645 = w37*w644;
    const scalar w646 = w28*w644;
    const scalar w647 = 4*w639;
    const scalar w648 = -w159*w647 + w647*w73;
    const scalar w649 = -w58*w647 - w60*w647;
    const scalar w650 = w180*w639 + w647*w77;
    const scalar w651 = -w391*w647 + w647*w65;
    const scalar w652 = w146*w247 + w150*w248 + w157*w303 + w163*w305 + w62*w648 + w649*w75 + w650*w67 + w651*w79;
    const scalar w653 = w587*x;
    const scalar w654 = w551*w587;
    const scalar w655 = 2*gil*theta1*w0*w13*w157*w6*x + 2*gil*theta1*w0*w13*w247*w6*z + gil*w0*w14*w648 + git*w0*w13*w650 - w415*w653 - w647*w76 - w647*w80 - w654*z;
    const scalar w656 = Msi0*w212;
    const scalar w657 = Msi4*w516;
    const scalar w658 = Msi5*w516;
    const scalar w659 = Msi2*w223;
    const scalar w660 = Msi2*w88;
    const scalar w661 = Msi2*w214;
    const scalar w662 = w165*w522 + w176*w536 + w187*w534 + w187*w538 + w216*w655 - w221*w660 + w225*w652 + w249*w359 + w249*w361 + w300*w382 + w300*w386 + w309*w380 + w309*w384 + w328*w419 + w369 + w383 + w418*w656 + w601 + w657*w88 + w658*w82 - w659*w91 - w661*w82;
    const scalar w663 = w226*w642 + w420*w643 + w645*w88 + w646*w82 + w662;
    const scalar w664 = w62*(-w58*w640 - w60*w640) + w67*(-w391*w640 + w640*w65);
    const scalar w665 = w171*w303 + w174*w305 + w664;
    const scalar w666 = w371*w643;
    const scalar w667 = w142*x;
    const scalar w668 = w142*w305;
    const scalar w669 = 2*gil*theta1*w0*w13*w145*w6*x + 2*gil*theta1*w0*w13*w302*w6*z + gil*w0*w14*w649 + git*w0*w13*w651 - w150*w667 - w64*w647 - w647*w69 - w668*z;
    const scalar w670 = Msi1*w86;
    const scalar w671 = Msi2*w352;
    const scalar w672 = Msi2*w221;
    const scalar w673 = w151*w382 + w151*w522 + w175*w472 + w175*w536 + w176*w534 + w176*w538 + w216*w669 + w226*w652 + w306*w359 + w306*w361 + w309*w382 + w309*w386 + w314*w380 + w314*w384 + w352*w670 + w407*w656 + w516*w660 + w657*w86 + w658*w71 - w659*w88 - w661*w71 + w671*w71 - w672*w86;
    const scalar w674 = w225*w665 + w645*w86 + w646*w71 + w666*w88 + w673;
    const scalar w675 = 8*gil*w1*w373*w374*w49*x*z + 8*git*w1*w373*w374*w47*x*z - w48*w640 - w50*w640;
    const scalar w676 = Msi1*w71;
    const scalar w677 = w306*w380 + w657*w71;
    const scalar w678 = w151*w472 + w151*w536 + w165*w534 + w165*w538 + w21*w382 + w21*w522 + w216*w675 + w225*w669 + w226*w655 + w249*w386 + w301*w359 + w301*w361 + w306*w384 + w352*w676 + w367*w656 + w53*w658 - w53*w661 + w53*w671 + w554 + w556 - w659*w82 - w672*w71 + w677;
    const scalar w679 = w53*w646 + w645*w71 + w666*w82 + w678;
    const scalar w680 = w236 + w340;
    const scalar w681 = w187*w193 + w193*w249 + w194*w655 + w195*w652;
    const scalar w682 = w193*w306;
    const scalar w683 = w194*w669 + w196*w652 + w610 + w682;
    const scalar w684 = w113*(-w165*w193 - w193*w301 - w194*w675 - w195*w669 - w196*w655) + w203 + w319;
    const scalar w685 = w231 + w334;
    const scalar w686 = w137*w187 + w137*w249 + w138*w655 + w152*w652;
    const scalar w687 = w137*w306;
    const scalar w688 = w138*w669 + w166*w652 + w616 + w687;
    const scalar w689 = w188 + w307 + w45*(w137*w165 + w137*w301 + w138*w675 + w152*w669 + w166*w655);
    const scalar w690 = w132*w338;
    const scalar w691 = (1.0/2.0)*w327;
    const scalar w692 = (1.0/2.0)*w207*w343;
    const scalar w693 = w342*w430;
    const scalar w694 = w318*w442;
    const scalar w695 = w233*w318;
    const scalar w696 = (1.0/2.0)*w339;
    const scalar w697 = w325*w440;
    const scalar w698 = -w113*w693 - w113*w697 + w191*w567 + w191*w696 - w205*w566 + w205*w692 + w232*w563 - w237*w562 + w237*w691 + w45*w690 + w45*w694 + w620*w695;
    const scalar w699 = w238*x;
    const scalar w700 = w344*z;
    const scalar w701 = w348*x;
    const scalar w702 = w219 + w623*w699 + w623*w700 + w625*w699 + w625*w700 + w701*z;
    const scalar w703 = w157*w558 + w163*w559 + w641;
    const scalar w704 = w328*w355;
    const scalar w705 = w354*w704;
    const scalar w706 = w37*w705;
    const scalar w707 = w28*w705;
    const scalar w708 = w226*w703 + w420*w704 + w662 + w706*w88 + w707*w82;
    const scalar w709 = w146*w312 + w150*w313 + w664;
    const scalar w710 = w371*w704;
    const scalar w711 = w225*w709 + w673 + w706*w86 + w707*w71 + w710*w88;
    const scalar w712 = w53*w707 + w678 + w706*w71 + w710*w82;
    const scalar w713 = w638*y;
    const scalar w714 = 8*w713;
    const scalar w715 = w75*(-w571*w713 + w714*w73) + w79*(w59*w714 + w714*w77);
    const scalar w716 = w247*w500 + w248*w501 + w715;
    const scalar w717 = w279*w518;
    const scalar w718 = w354*w717;
    const scalar w719 = w37*w718;
    const scalar w720 = w28*w718;
    const scalar w721 = 4*w713;
    const scalar w722 = -w159*w721 + w721*w73;
    const scalar w723 = -w58*w721 - w60*w721;
    const scalar w724 = w180*w713 + w721*w77;
    const scalar w725 = -w391*w721 + w65*w721;
    const scalar w726 = w243*w247 + w245*w248 + w254*w303 + w255*w305 + w62*w722 + w67*w724 + w723*w75 + w725*w79;
    const scalar w727 = 2*gil*theta1*w0*w13*w247*w6*y + 2*gil*theta1*w0*w13*w254*w6*x + gil*w0*w14*w722 + git*w0*w13*w724 - w495*w653 - w654*y - w721*w76 - w721*w80;
    const scalar w728 = Msi0*w279;
    const scalar w729 = Msi1*w223;
    const scalar w730 = Msi1*w214;
    const scalar w731 = w216*w727 - w221*w490 + w225*w726 + w249*w380 + w249*w458 + w256*w522 + w259*w536 + w265*w534 + w265*w538 + w300*w472 + w300*w476 + w309*w470 + w309*w474 + w418*w526 + w418*w728 + w491 + w493 + w499*w516 + w605 + w657*w82 - w729*w91 - w730*w82;
    const scalar w732 = w226*w716 + w420*w717 + w719*w88 + w720*w82 + w731;
    const scalar w733 = w62*(-w58*w714 - w60*w714) + w67*(-w391*w714 + w65*w714);
    const scalar w734 = w250*w303 + w251*w305 + w733;
    const scalar w735 = w371*w717;
    const scalar w736 = 2*gil*theta1*w0*w13*w242*w6*x + 2*gil*theta1*w0*w13*w302*w6*y + gil*w0*w14*w723 + git*w0*w13*w725 - w245*w667 - w64*w721 - w668*y - w69*w721;
    const scalar w737 = Msi3*w516;
    const scalar w738 = w216*w736 - w221*w670 - w223*w490 + w226*w726 + w246*w382 + w246*w522 + w252*w472 + w252*w536 + w259*w534 + w259*w538 + w306*w458 + w309*w476 + w314*w470 + w314*w474 + w407*w728 + w453*w670 + w555 + w557 + w599*w71 + w677 - w71*w730 + w737*w86;
    const scalar w739 = w225*w734 + w71*w720 + w719*w86 + w735*w88 + w738;
    const scalar w740 = 8*gil*w1*w373*w374*w49*x*y + 8*git*w1*w373*w374*w47*x*y - w48*w714 - w50*w714;
    const scalar w741 = w216*w740 - w221*w676 + w225*w736 + w226*w727 + w246*w472 + w246*w536 + w249*w472 + w249*w476 + w253*w382 + w253*w522 + w256*w534 + w256*w538 + w301*w380 + w301*w458 + w306*w470 + w306*w474 + w367*w526 + w367*w728 + w453*w676 + w53*w599 + w53*w657 - w53*w730 + w71*w737 - w729*w82;
    const scalar w742 = w53*w720 + w71*w719 + w735*w82 + w741;
    const scalar w743 = w293 + w341;
    const scalar w744 = w193*w265 + w193*w309 + w194*w727 + w195*w726;
    const scalar w745 = w193*w259 + w193*w314 + w194*w736 + w196*w726;
    const scalar w746 = w113*(-w194*w740 - w195*w736 - w196*w727 - w609 - w682) + w275 + w323;
    const scalar w747 = w290 + w336;
    const scalar w748 = w137*w265 + w137*w309 + w138*w727 + w152*w726;
    const scalar w749 = w137*w259 + w137*w314 + w138*w736 + w166*w726;
    const scalar w750 = w266 + w315 + w45*(w138*w740 + w152*w736 + w166*w727 + w615 + w687);
    const scalar w751 = -w114*w693 - w114*w697 + w269*w567 + w269*w696 - w277*w566 + w277*w692 + w291*w563 - w294*w562 + w294*w691 + w54*w690 + w54*w694 + w621*w695;
    const scalar w752 = w295*x;
    const scalar w753 = w344*y;
    const scalar w754 = w282 + w623*w752 + w623*w753 + w625*w752 + w625*w753 + w701*y;
    const scalar w755 = w254*w558 + w255*w559 + w715;
    const scalar w756 = w328*w454;
    const scalar w757 = w354*w756;
    const scalar w758 = w37*w757;
    const scalar w759 = w28*w757;
    const scalar w760 = w226*w755 + w420*w756 + w731 + w758*w88 + w759*w82;
    const scalar w761 = w243*w312 + w245*w313 + w733;
    const scalar w762 = w371*w756;
    const scalar w763 = w225*w761 + w71*w759 + w738 + w758*w86 + w762*w88;
    const scalar w764 = w53*w759 + w71*w758 + w741 + w762*w82;
    return c0*sqrt(w240*(w241*w51 + w297*w70 + w346*w81) + w296*(w241*w70 + w297*w85 + w346*w87) + w345*(w241*w81 + w297*w87 + w346*w90)) - w151*w297 - w165*w346 - w21*w241 - w241*w246 - w241*w249 - w252*w297 - w259*w346 - w297*w309 - w300*w346 - w53*(w211 + w215 + w239*(w123*(w113*w404 + w114*w417 + w115*w426 + w227*w349) + 2*w191*w443 - 2*w205*w441 + w209*w237 + 2*w232*w432 + w232*w434 + w234*(-w227*w429 - w404*w45 - w417*w54 - w426*w72) - 2*w237*w431 + w428*(w113*(-w194*w379 - w195*w396 - w196*w403 - w21*w427) + w114*(-w151*w427 - w194*w396 - w195*w414 - w196*w416) + w115*(-w165*w427 - w194*w403 - w195*w416 - w196*w425) + w197*w349) + w435*w437 + w439*(-w167*w429 - w45*(w138*w379 + w152*w396 + w166*w403 + w21*w438) - w54*(w138*w396 + w151*w438 + w152*w414 + w166*w416) - w72*(w138*w403 + w152*w416 + w165*w438 + w166*w425)) + w444*w446 + w447) + w348*w9 + w448*w449 + w449*w450 + w451) - w71*(w222 + w239*(w123*(w113*w637 + w114*w636 + w115*w633 + w608) + w234*(-w45*w637 - w54*w636 - w614 - w633*w72) + w428*(w114*(-w195*w634 - w612) + w115*(-w196*w628 - w611) + w613) + w435*w508 + w439*(-w54*(w152*w634 + w618) - w619 - w72*(w166*w628 + w617)) + w444*w512 + w622) + w627) - w71*(w239*(w123*(w113*w607 + w114*w603 + w115*w594 + w608) + w234*(-w45*w607 - w54*w603 - w594*w72 - w614) + w428*(w114*(-w195*w596 - w612) + w115*(-w196*w576 - w611) + w613) + w437*w507 + w439*(-w54*(w152*w596 + w618) - w619 - w72*(w166*w576 + w617)) + w446*w511 + w622) + w280 + w627) - w82*(w224 + w239*(w123*(w113*w712 + w114*w711 + w115*w708 + w680) + w234*(-w45*w712 - w54*w711 - w685 - w708*w72) + w428*(w114*(-w195*w709 - w683) + w115*(-w196*w703 - w681) + w684) + w435*w565 + w439*(-w54*(w152*w709 + w688) - w689 - w72*(w166*w703 + w686)) + w444*w569 + w698) + w702) - w82*(w239*(w123*(w113*w679 + w114*w674 + w115*w663 + w680) + w234*(-w45*w679 - w54*w674 - w663*w72 - w685) + w428*(w114*(-w195*w665 - w683) + w115*(-w196*w642 - w681) + w684) + w437*w564 + w439*(-w54*(w152*w665 + w688) - w689 - w72*(w166*w642 + w686)) + w446*w568 + w698) + w329 + w702) - w86*(w239*(w123*(w113*w486 + w114*w497 + w115*w503 + w288*w349) + w234*(-w288*w429 - w45*w486 - w497*w54 - w503*w72) + 2*w269*w510 - 2*w277*w509 + w278*w294 + 2*w291*w505 + w291*w506 - 2*w294*w504 + w428*(w113*(-w194*w469 - w195*w482 - w196*w485 - w246*w427) + w114*(-w194*w482 - w195*w494 - w196*w496 - w252*w427) + w115*(-w194*w485 - w195*w496 - w196*w502 - w259*w427) + w272*w349) + w439*(-w260*w429 - w45*(w138*w469 + w152*w482 + w166*w485 + w246*w438) - w54*(w138*w482 + w152*w494 + w166*w496 + w252*w438) - w72*(w138*w485 + w152*w496 + w166*w502 + w259*w438)) + w447 + w507*w508 + w511*w512) + w281 + w283 + w348*w8 + w448*w513 + w450*w513 + w451) - w88*(w239*(w123*(w113*w742 + w114*w739 + w115*w732 + w743) + w234*(-w45*w742 - w54*w739 - w72*w732 - w747) + w428*(w114*(-w195*w734 - w745) + w115*(-w196*w716 - w744) + w746) + w439*(-w54*(w152*w734 + w749) - w72*(w166*w716 + w748) - w750) + w508*w564 + w512*w568 + w751) + w331 + w754) - w88*(w239*(w123*(w113*w764 + w114*w763 + w115*w760 + w743) + w234*(-w45*w764 - w54*w763 - w72*w760 - w747) + w428*(w114*(-w195*w761 - w745) + w115*(-w196*w755 - w744) + w746) + w439*(-w54*(w152*w761 + w749) - w72*(w166*w755 + w748) - w750) + w507*w565 + w511*w569 + w751) + w284 + w754) - w91*(w239*(w123*(w113*w549 + w114*w553 + w115*w561 + w337*w349) + w234*(-w337*w429 - w45*w549 - w54*w553 - w561*w72) + w318*w339 + 2*w318*w567 - 2*w325*w566 + w326*w343 + 2*w338*w563 - 2*w342*w562 + w428*(w113*(-w194*w533 - w195*w545 - w196*w548 - w249*w427) + w114*(-w194*w545 - w195*w550 - w196*w552 - w309*w427) + w115*(-w194*w548 - w195*w552 - w196*w560 - w300*w427) + w321*w349) + w439*(-w310*w429 - w45*(w138*w533 + w152*w545 + w166*w548 + w249*w438) - w54*(w138*w545 + w152*w550 + w166*w552 + w309*w438) - w72*(w138*w548 + w152*w552 + w166*w560 + w300*w438)) + w447 + w564*w565 + w568*w569) + w330 + w332 + w347*w514 + w448*w570 + w450*w570 + w451) - 1;
}

inline symmTensor sigma(const point& pt, const Params& p)
{
    const scalar r1 = p.r1;
    const scalar r2 = p.r2;
    const scalar gamma1 = p.gamma1;
    const scalar gamma2 = p.gamma2;
    const scalar gil = p.gil;
    const scalar git = p.git;
    const scalar gin = p.gin;
    const scalar theta0 = p.theta0;
    const scalar theta1 = p.theta1;
    const scalar x = pt.x();
    const scalar y = pt.y();
    const scalar z = pt.z();

    const scalar w0 = sin(gamma1);
    const scalar w1 = sin(gamma2);
    const scalar w2 = (w1*w1);
    const scalar w3 = gin*w2;
    const scalar w4 = cos(gamma1);
    const scalar w5 = (r1*r1);
    const scalar w6 = theta0 + theta1*(-w5 + (x*x) + (y*y) + (z*z))/((r2*r2) - w5);
    const scalar w7 = cos(w6);
    const scalar w8 = cos(gamma2);
    const scalar w9 = sin(w6);
    const scalar w10 = w0*w9;
    const scalar w11 = -w10*w8 + w4*w7;
    const scalar w12 = w4*w9;
    const scalar w13 = w0*w7;
    const scalar w14 = -w12 - w13*w8;
    const scalar w15 = w12*w8 + w13;
    const scalar w16 = gil*w11;
    const scalar w17 = -w10 + w4*w7*w8;
    const scalar w18 = git*w14;
    const scalar w19 = gin*w1*w8;
    const scalar w20 = w1*w9;
    const scalar w21 = w1*w7;
    return symmTensor(gil*(w11*w11) + git*(w14*w14) + (w0*w0)*w3, -w0*w3*w4 + w15*w16 + w17*w18, w0*w19 + w16*w20 + w18*w21, gil*(w15*w15) + git*(w17*w17) + w3*(w4*w4), gil*w15*w20 + git*w17*w21 - w19*w4, gil*w2*(w9*w9) + gin*(w8*w8) + git*w2*(w7*w7));
}

struct SelfTestSample
{
    label fn;
    label param;
    label solution;
    label zone;
    label component;
    scalar x, y, z, t, value;
};

inline const SelfTestSample selfTestSamples[] =
{
    {0, 0, 0, 1, 0, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, 9.3985310143272542},
    {0, 0, 0, 1, 0, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 14.814692639850362},
    {0, 0, 0, 1, 0, -0.12, 0.64000000000000001, 0.70999999999999996, 0, 19.078897677070074},
    {0, 1, 0, 1, 0, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, 9.5565516219834965},
    {0, 1, 0, 1, 0, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 15.842181149121187},
    {0, 1, 0, 1, 0, -0.12, 0.64000000000000001, 0.70999999999999996, 0, 15.945596600352482},
    {1, 0, 0, 1, 0, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, -7.0036365056872976},
    {1, 0, 0, 1, 1, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, 14.584055108628839},
    {1, 0, 0, 1, 2, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, -16.650790406477849},
    {1, 0, 0, 1, 0, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 8.1826297038048619},
    {1, 0, 0, 1, 1, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, -20.621836632627573},
    {1, 0, 0, 1, 2, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 15.163531441268217},
    {1, 0, 0, 1, 0, -0.12, 0.64000000000000001, 0.70999999999999996, 0, -9.7634037636906363},
    {1, 0, 0, 1, 1, -0.12, 0.64000000000000001, 0.70999999999999996, 0, 7.6610018251784302},
    {1, 0, 0, 1, 2, -0.12, 0.64000000000000001, 0.70999999999999996, 0, 5.078690821483324},
    {1, 1, 0, 1, 0, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, -6.5882756084170682},
    {1, 1, 0, 1, 1, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, 9.6473290512424512},
    {1, 1, 0, 1, 2, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, -18.54869664061588},
    {1, 1, 0, 1, 0, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 7.5760927567767435},
    {1, 1, 0, 1, 1, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, -21.038959788705938},
    {1, 1, 0, 1, 2, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 11.643796930074252},
    {1, 1, 0, 1, 0, -0.12, 0.64000000000000001, 0.70999999999999996, 0, -6.2862615219541329},
    {1, 1, 0, 1, 1, -0.12, 0.64000000000000001, 0.70999999999999996, 0, 8.8133451429891956},
    {1, 1, 0, 1, 2, -0.12, 0.64000000000000001, 0.70999999999999996, 0, 4.7468115375382975},
    {2, 0, 0, 1, 0, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, -0.29956970838805919},
    {2, 0, 0, 1, 0, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 0.035359330129071616},
    {2, 0, 0, 1, 0, -0.12, 0.64000000000000001, 0.70999999999999996, 0, 3.4394229227260102},
    {2, 1, 0, 1, 0, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, -0.92946209604779273},
    {2, 1, 0, 1, 0, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 0.1864695848774669},
    {2, 1, 0, 1, 0, -0.12, 0.64000000000000001, 0.70999999999999996, 0, 2.2641900722337445},
    {3, 0, 0, 1, 0, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, 0.15738472957834365},
    {3, 0, 0, 1, 1, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, 0.019450757721302565},
    {3, 0, 0, 1, 2, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, -0.038292011629080386},
    {3, 0, 0, 1, 3, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, 0.016603840783215282},
    {3, 0, 0, 1, 4, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, -0.00053267868819488625},
    {3, 0, 0, 1, 5, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, 0.025011429638441073},
    {3, 0, 0, 1, 0, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 0.097201186639712764},
    {3, 0, 0, 1, 1, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 0.065832905238739109},
    {3, 0, 0, 1, 2, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 0.038181559989613846},
    {3, 0, 0, 1, 3, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 0.066523826271050968},
    {3, 0, 0, 1, 4, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 0.036466280788684788},
    {3, 0, 0, 1, 5, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 0.035274987089236293},
    {3, 0, 0, 1, 0, -0.12, 0.64000000000000001, 0.70999999999999996, 0, 0.02837550659096431},
    {3, 0, 0, 1, 1, -0.12, 0.64000000000000001, 0.70999999999999996, 0, 0.028658266129057097},
    {3, 0, 0, 1, 2, -0.12, 0.64000000000000001, 0.70999999999999996, 0, 0.024603952904978122},
    {3, 0, 0, 1, 3, -0.12, 0.64000000000000001, 0.70999999999999996, 0, 0.087073767876423389},
    {3, 0, 0, 1, 4, -0.12, 0.64000000000000001, 0.70999999999999996, 0, 0.074632313391501143},
    {3, 0, 0, 1, 5, -0.12, 0.64000000000000001, 0.70999999999999996, 0, 0.083550725532612308},
    {3, 1, 0, 1, 0, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, 0.15730815404489815},
    {3, 1, 0, 1, 1, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, 0.019197104265080218},
    {3, 1, 0, 1, 2, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, -0.038544861109930509},
    {3, 1, 0, 1, 3, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, 0.016540724988994872},
    {3, 1, 0, 1, 4, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, -0.00049997083760299551},
    {3, 1, 0, 1, 5, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, 0.025151120966107002},
    {3, 1, 0, 1, 0, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 0.12218473194461883},
    {3, 1, 0, 1, 1, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 0.065392252779922957},
    {3, 1, 0, 1, 2, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 0.0275046948391926},
    {3, 1, 0, 1, 3, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 0.053420525692549958},
    {3, 1, 0, 1, 4, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 0.021971741174393843},
    {3, 1, 0, 1, 5, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 0.023394742362831223},
    {3, 1, 0, 1, 0, -0.12, 0.64000000000000001, 0.70999999999999996, 0, 0.061321139387836475},
    {3, 1, 0, 1, 1, -0.12, 0.64000000000000001, 0.70999999999999996, 0, 0.055093492140707412},
    {3, 1, 0, 1, 2, -0.12, 0.64000000000000001, 0.70999999999999996, 0, 0.040779525216660888},
    {3, 1, 0, 1, 3, -0.12, 0.64000000000000001, 0.70999999999999996, 0, 0.080736274206694686},
    {3, 1, 0, 1, 4, -0.12, 0.64000000000000001, 0.70999999999999996, 0, 0.056759971100000188},
    {3, 1, 0, 1, 5, -0.12, 0.64000000000000001, 0.70999999999999996, 0, 0.056942586405468849},
};

inline const label nSelfTestSamples = 66;

} // End namespace shellEikonalMMS
} // End namespace Foam

#endif
