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
    Exact field F, its gradient, div(G grad F), the conductivity G and the
    positive ionic function g of the heart-only monodomain manufactured
    solution on the spherical shell, sealed on both spheres, with self-test
    samples.

    Regenerate from bathBidomainVerification/generator/:
        python shell_monodomain_mms.py --write ../../monodomainVerification/manufacturedShellMonodomainReference.H

    sympy 1.14.0
\*---------------------------------------------------------------------------*/

#ifndef manufacturedShellMonodomainReference_H
#define manufacturedShellMonodomainReference_H

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
namespace shellMonodomainMMS
{

struct Params
{
    scalar r1, r2, gamma1, gamma2, k;
    scalar gil, git, gin, theta0, theta1;
    scalar hA[3];
    label solution;
};

inline Params defaultParams()
{
    Params p;
    p.r1 = 0.5;
    p.r2 = 1.0;
    p.gamma1 = degToRad(20.0);
    p.gamma2 = degToRad(50.0);
    p.k = 3.141592653589793;
    p.gil = 0.17;
    p.git = 0.019;
    p.gin = 0.01;
    p.theta0 = degToRad(-20.0);
    p.theta1 = degToRad(90.0);
    p.hA[0] = -0.1;
    p.hA[1] = 0.3;
    p.hA[2] = 0.25;
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
    p.k = 3.141592653589793;
    p.gil = 0.17;
    p.git = 0.019;
    p.gin = 0.01;
    p.theta0 = degToRad(-20.0);
    p.theta1 = degToRad(70.0);
    p.hA[0] = -0.15;
    p.hA[1] = 0.35;
    p.hA[2] = 0.15;
    p.solution = 0;
    return p;
}

inline scalar F(const point& pt, const Params& p)
{
    const scalar r1 = p.r1;
    const scalar r2 = p.r2;
    const scalar gamma1 = p.gamma1;
    const scalar gamma2 = p.gamma2;
    const scalar k = p.k;
    const scalar gil = p.gil;
    const scalar git = p.git;
    const scalar gin = p.gin;
    const scalar theta0 = p.theta0;
    const scalar theta1 = p.theta1;
    const scalar hA0 = p.hA[0];
    const scalar hA1 = p.hA[1];
    const scalar hA2 = p.hA[2];
    const scalar x = pt.x();
    const scalar y = pt.y();
    const scalar z = pt.z();

    const scalar w0 = k*z;
    const scalar w1 = sin(w0);
    const scalar w2 = (x*x);
    const scalar w3 = 1 - 3.0/10.0*w2;
    const scalar w4 = k*y;
    const scalar w5 = hA1*cos(w4);
    const scalar w6 = r1 + r2;
    const scalar w7 = (1.0/(r2));
    const scalar w8 = cos(gamma2);
    const scalar w9 = sin(gamma2);
    const scalar w10 = (w9*w9);
    const scalar w11 = (r1*r1);
    const scalar w12 = -(r2*r2);
    const scalar w13 = w2 + (y*y) + (z*z);
    const scalar w14 = (-w11 + w13)/(-w11 - w12);
    const scalar w15 = theta0 + theta1*w14;
    const scalar w16 = sin(w15);
    const scalar w17 = cos(w15);
    const scalar w18 = gil*w10*(w16*w16) + gin*(w8*w8) + git*w10*(w17*w17);
    const scalar w19 = w7*z;
    const scalar w20 = cos(gamma1);
    const scalar w21 = gin*w8*w9;
    const scalar w22 = sin(gamma1);
    const scalar w23 = w17*w22;
    const scalar w24 = w16*w20;
    const scalar w25 = w23 + w24*w8;
    const scalar w26 = w16*w9;
    const scalar w27 = w16*w22;
    const scalar w28 = w17*w20*w8 - w27;
    const scalar w29 = w17*w9;
    const scalar w30 = gil*w25*w26 + git*w28*w29 - w20*w21;
    const scalar w31 = w7*y;
    const scalar w32 = w17*w20 - w27*w8;
    const scalar w33 = gil*w32;
    const scalar w34 = -w23*w8 - w24;
    const scalar w35 = git*w34;
    const scalar w36 = w21*w22 + w26*w33 + w29*w35;
    const scalar w37 = w7*x;
    const scalar w38 = gin*w10;
    const scalar w39 = gil*(w25*w25) + git*(w28*w28) + (w20*w20)*w38;
    const scalar w40 = -w20*w22*w38 + w25*w33 + w28*w35;
    const scalar w41 = gil*(w32*w32) + git*(w34*w34) + (w22*w22)*w38;
    const scalar w42 = hA0 + hA2*k*x*cos(w0);
    const scalar w43 = hA1*k*w3*sin(w4);
    const scalar w44 = hA2*w1 - 3.0/5.0*w5*x;
    const scalar w45 = w18*w42 - w30*w43 + w36*w44;
    const scalar w46 = w30*w42 - w39*w43 + w40*w44;
    const scalar w47 = w36*w42 - w40*w43 + w41*w44;
    const scalar w48 = (1.0/2.0)*w6;
    const scalar w49 = (1.0/(r1));
    const scalar w50 = w49*z;
    const scalar w51 = w49*y;
    const scalar w52 = w49*x;
    return hA0*z + hA2*w1*x + w14*(w12 + w13)*(w14*w48*w7*(-w19*w45 - w31*w46 - w37*w47)/(w19*(w18*w19 + w30*w31 + w36*w37) + w31*(w19*w30 + w31*w39 + w37*w40) + w37*(w19*w36 + w31*w40 + w37*w41)) + w48*w49*(1 - w14)*(w45*w50 + w46*w51 + w47*w52)/(-w50*(-w18*w50 - w30*w51 - w36*w52) - w51*(-w30*w50 - w39*w51 - w40*w52) - w52*(-w36*w50 - w40*w51 - w41*w52)))/w6 + w3*w5;
}

inline vector gradF(const point& pt, const Params& p)
{
    const scalar r1 = p.r1;
    const scalar r2 = p.r2;
    const scalar gamma1 = p.gamma1;
    const scalar gamma2 = p.gamma2;
    const scalar k = p.k;
    const scalar gil = p.gil;
    const scalar git = p.git;
    const scalar gin = p.gin;
    const scalar theta0 = p.theta0;
    const scalar theta1 = p.theta1;
    const scalar hA0 = p.hA[0];
    const scalar hA1 = p.hA[1];
    const scalar hA2 = p.hA[2];
    const scalar x = pt.x();
    const scalar y = pt.y();
    const scalar z = pt.z();

    const scalar w0 = k*z;
    const scalar w1 = sin(w0);
    const scalar w2 = hA2*w1;
    const scalar w3 = k*y;
    const scalar w4 = hA1*cos(w3);
    const scalar w5 = (3.0/5.0)*w4;
    const scalar w6 = w5*x;
    const scalar w7 = r1 + r2;
    const scalar w8 = (1.0/(w7));
    const scalar w9 = (r1*r1);
    const scalar w10 = -(r2*r2);
    const scalar w11 = (1.0/((-w10 - w9)));
    const scalar w12 = (x*x);
    const scalar w13 = w12 + (y*y) + (z*z);
    const scalar w14 = w11*(w13 - w9);
    const scalar w15 = w14*w8;
    const scalar w16 = (1.0/(r2));
    const scalar w17 = k*x;
    const scalar w18 = hA2*cos(w0);
    const scalar w19 = hA0 + w17*w18;
    const scalar w20 = cos(gamma2);
    const scalar w21 = theta0 + theta1*w14;
    const scalar w22 = sin(w21);
    const scalar w23 = sin(gamma2);
    const scalar w24 = (w23*w23);
    const scalar w25 = gil*w24;
    const scalar w26 = cos(w21);
    const scalar w27 = git*w24;
    const scalar w28 = gin*(w20*w20) + (w22*w22)*w25 + (w26*w26)*w27;
    const scalar w29 = cos(gamma1);
    const scalar w30 = gin*w20*w23;
    const scalar w31 = sin(gamma1);
    const scalar w32 = w26*w31;
    const scalar w33 = w22*w29;
    const scalar w34 = w20*w33;
    const scalar w35 = w32 + w34;
    const scalar w36 = gil*w35;
    const scalar w37 = w22*w23;
    const scalar w38 = w22*w31;
    const scalar w39 = w20*w26*w29 - w38;
    const scalar w40 = git*w39;
    const scalar w41 = w23*w26;
    const scalar w42 = -w29*w30 + w36*w37 + w40*w41;
    const scalar w43 = 1 - 3.0/10.0*w12;
    const scalar w44 = hA1*sin(w3);
    const scalar w45 = k*w43*w44;
    const scalar w46 = hA2*w1 - w6;
    const scalar w47 = w20*w38;
    const scalar w48 = w26*w29 - w47;
    const scalar w49 = gil*w48;
    const scalar w50 = w20*w32;
    const scalar w51 = -w33 - w50;
    const scalar w52 = git*w51;
    const scalar w53 = w30*w31 + w37*w49 + w41*w52;
    const scalar w54 = w19*w28 - w42*w45 + w46*w53;
    const scalar w55 = w16*w54;
    const scalar w56 = gin*w24;
    const scalar w57 = gil*(w35*w35) + git*(w39*w39) + (w29*w29)*w56;
    const scalar w58 = -w29*w31*w56 + w35*w49 + w39*w52;
    const scalar w59 = w19*w42 - w45*w57 + w46*w58;
    const scalar w60 = w16*w59;
    const scalar w61 = gil*(w48*w48) + git*(w51*w51) + (w31*w31)*w56;
    const scalar w62 = w19*w53 - w45*w58 + w46*w61;
    const scalar w63 = w16*w62;
    const scalar w64 = -w55*z - w60*y - w63*x;
    const scalar w65 = w16*w28;
    const scalar w66 = w16*w42;
    const scalar w67 = w16*w53;
    const scalar w68 = w16*(w65*z + w66*y + w67*x);
    const scalar w69 = w16*w57;
    const scalar w70 = w16*w58;
    const scalar w71 = w16*(w66*z + w69*y + w70*x);
    const scalar w72 = w16*w61;
    const scalar w73 = w16*(w67*z + w70*y + w72*x);
    const scalar w74 = w68*z + w71*y + w73*x;
    const scalar w75 = w7/w74;
    const scalar w76 = w64*w75;
    const scalar w77 = (1.0/2.0)*w14*w16;
    const scalar w78 = (1.0/(r1));
    const scalar w79 = w54*w78;
    const scalar w80 = w59*w78;
    const scalar w81 = w62*w78;
    const scalar w82 = w79*z + w80*y + w81*x;
    const scalar w83 = w28*w78;
    const scalar w84 = w42*w78;
    const scalar w85 = w53*w78;
    const scalar w86 = w78*(-w83*z - w84*y - w85*x);
    const scalar w87 = w57*w78;
    const scalar w88 = w58*w78;
    const scalar w89 = w78*(-w84*z - w87*y - w88*x);
    const scalar w90 = w61*w78;
    const scalar w91 = w78*(-w85*z - w88*y - w90*x);
    const scalar w92 = -w86*z - w89*y - w91*x;
    const scalar w93 = w7/w92;
    const scalar w94 = (1.0/2.0)*w78*(1 - w14);
    const scalar w95 = w93*w94;
    const scalar w96 = w76*w77 + w82*w95;
    const scalar w97 = 2*x;
    const scalar w98 = w96*w97;
    const scalar w99 = w10 + w13;
    const scalar w100 = w11*w8*w99;
    const scalar w101 = w16*x;
    const scalar w102 = w11*w76;
    const scalar w103 = w78*x;
    const scalar w104 = w11*w82*w93;
    const scalar w105 = k*w18;
    const scalar w106 = theta1*w11;
    const scalar w107 = 4*x;
    const scalar w108 = w106*w107;
    const scalar w109 = w22*w26;
    const scalar w110 = w108*w109;
    const scalar w111 = w110*w25 - w110*w27;
    const scalar w112 = (3.0/5.0)*w17*w44;
    const scalar w113 = w106*w97;
    const scalar w114 = w113*w41;
    const scalar w115 = w113*w37;
    const scalar w116 = 2*theta1*w11*w20*w26*w29*x - w113*w38;
    const scalar w117 = gil*w37;
    const scalar w118 = -w113*w32 - w113*w34;
    const scalar w119 = git*w41;
    const scalar w120 = w114*w36 - w115*w40 + w116*w117 + w118*w119;
    const scalar w121 = w106*w33;
    const scalar w122 = -w113*w50 - w121*w97;
    const scalar w123 = w106*w26*w29;
    const scalar w124 = w113*w47 - w123*w97;
    const scalar w125 = w114*w49 - w115*w52 + w117*w122 + w119*w124;
    const scalar w126 = w105*w28 + w111*w19 + w112*w42 - w120*w45 + w125*w46 - w5*w53;
    const scalar w127 = w16*z;
    const scalar w128 = w36*(4*theta1*w11*w20*w26*w29*x - w108*w38) + w40*(-w108*w32 - w108*w34);
    const scalar w129 = w116*w49 + w118*w52 + w122*w36 + w124*w40;
    const scalar w130 = w105*w42 + w112*w57 + w120*w19 - w128*w45 + w129*w46 - w5*w58;
    const scalar w131 = w16*y;
    const scalar w132 = w105*w53;
    const scalar w133 = w112*w58;
    const scalar w134 = w49*(-w107*w121 - w108*w50) + w52*(-w107*w123 + w108*w47);
    const scalar w135 = w125*w19 - w129*w45 + w132 + w133 + w134*w46 - w5*w61;
    const scalar w136 = w75*w77;
    const scalar w137 = w78*z;
    const scalar w138 = w78*y;
    const scalar w139 = w64*w7*w77/(w74*w74);
    const scalar w140 = w7*w82*w94/(w92*w92);
    const scalar w141 = w15*w99;
    const scalar w142 = 2*y;
    const scalar w143 = w142*w96;
    const scalar w144 = 4*y;
    const scalar w145 = w106*w144;
    const scalar w146 = w109*w145;
    const scalar w147 = w146*w25 - w146*w27;
    const scalar w148 = (k*k);
    const scalar w149 = w148*w4*w43;
    const scalar w150 = w106*w142;
    const scalar w151 = w150*w41;
    const scalar w152 = w150*w37;
    const scalar w153 = 2*theta1*w11*w20*w26*w29*y - w150*w38;
    const scalar w154 = -w150*w32 - w150*w34;
    const scalar w155 = w117*w153 + w119*w154 + w151*w36 - w152*w40;
    const scalar w156 = -w121*w142 - w150*w50;
    const scalar w157 = -w123*w142 + w150*w47;
    const scalar w158 = w117*w156 + w119*w157 + w151*w49 - w152*w52;
    const scalar w159 = w112*w53 + w147*w19 - w149*w42 - w155*w45 + w158*w46;
    const scalar w160 = w36*(4*theta1*w11*w20*w26*w29*y - w145*w38) + w40*(-w145*w32 - w145*w34);
    const scalar w161 = w153*w49 + w154*w52 + w156*w36 + w157*w40;
    const scalar w162 = w133 - w149*w57 + w155*w19 - w160*w45 + w161*w46;
    const scalar w163 = w49*(-w121*w144 - w145*w50) + w52*(-w123*w144 + w145*w47);
    const scalar w164 = w112*w61 - w149*w58 + w158*w19 - w161*w45 + w163*w46;
    const scalar w165 = w148*w2*x;
    const scalar w166 = 4*z;
    const scalar w167 = w106*w166;
    const scalar w168 = w109*w167;
    const scalar w169 = w168*w25 - w168*w27;
    const scalar w170 = 2*z;
    const scalar w171 = w106*w170;
    const scalar w172 = w171*w41;
    const scalar w173 = w171*w37;
    const scalar w174 = 2*theta1*w11*w20*w26*w29*z - w171*w38;
    const scalar w175 = -w171*w32 - w171*w34;
    const scalar w176 = w117*w174 + w119*w175 + w172*w36 - w173*w40;
    const scalar w177 = -w121*w170 - w171*w50;
    const scalar w178 = -w123*w170 + w171*w47;
    const scalar w179 = w117*w177 + w119*w178 + w172*w49 - w173*w52;
    const scalar w180 = w132 - w165*w28 + w169*w19 - w176*w45 + w179*w46;
    const scalar w181 = w49*(-w121*w166 - w167*w50) + w52*(-w123*w166 + w167*w47);
    const scalar w182 = w174*w49 + w175*w52 + w177*w36 + w178*w40;
    const scalar w183 = w105*w61 - w165*w53 + w179*w19 + w181*w46 - w182*w45;
    const scalar w184 = w36*(4*theta1*w11*w20*w26*w29*z - w167*w38) + w40*(-w167*w32 - w167*w34);
    const scalar w185 = w105*w58 - w165*w42 + w176*w19 + w182*w46 - w184*w45;
    const scalar w186 = w170*w96;
    return vector(w100*w98 + w141*(w101*w102 - w103*w104 + w136*(-w101*w135 - w126*w127 - w130*w131 - w63) + w139*(-w101*(w101*w134 + w125*w127 + w129*w131 + w72) - w127*(w101*w125 + w111*w127 + w120*w131 + w67) - w131*(w101*w129 + w120*w127 + w128*w131 + w70) - w73) + w140*(w103*(-w103*w134 - w125*w137 - w129*w138 - w90) + w137*(-w103*w125 - w111*w137 - w120*w138 - w85) + w138*(-w103*w129 - w120*w137 - w128*w138 - w88) + w91) + w95*(w103*w135 + w126*w137 + w130*w138 + w81)) + w15*w98 + w2 - w6, w100*w143 + w141*(w102*w131 - w104*w138 + w136*(-w101*w164 - w127*w159 - w131*w162 - w60) + w139*(-w101*(w101*w163 + w127*w158 + w131*w161 + w70) - w127*(w101*w158 + w127*w147 + w131*w155 + w66) - w131*(w101*w161 + w127*w155 + w131*w160 + w69) - w71) + w140*(w103*(-w103*w163 - w137*w158 - w138*w161 - w88) + w137*(-w103*w158 - w137*w147 - w138*w155 - w84) + w138*(-w103*w161 - w137*w155 - w138*w160 - w87) + w89) + w95*(w103*w164 + w137*w159 + w138*w162 + w80)) + w143*w15 - w45, w100*w186 + w141*(w102*w127 - w104*w137 + w136*(-w101*w183 - w127*w180 - w131*w185 - w55) + w139*(-w101*(w101*w181 + w127*w179 + w131*w182 + w67) - w127*(w101*w179 + w127*w169 + w131*w176 + w65) - w131*(w101*w182 + w127*w176 + w131*w184 + w66) - w68) + w140*(w103*(-w103*w181 - w137*w179 - w138*w182 - w85) + w137*(-w103*w179 - w137*w169 - w138*w176 - w83) + w138*(-w103*w182 - w137*w176 - w138*w184 - w84) + w86) + w95*(w103*w183 + w137*w180 + w138*w185 + w79)) + w15*w186 + w19);
}

inline scalar divGGradF(const point& pt, const Params& p)
{
    const scalar r1 = p.r1;
    const scalar r2 = p.r2;
    const scalar gamma1 = p.gamma1;
    const scalar gamma2 = p.gamma2;
    const scalar k = p.k;
    const scalar gil = p.gil;
    const scalar git = p.git;
    const scalar gin = p.gin;
    const scalar theta0 = p.theta0;
    const scalar theta1 = p.theta1;
    const scalar hA0 = p.hA[0];
    const scalar hA1 = p.hA[1];
    const scalar hA2 = p.hA[2];
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
    const scalar w22 = (1.0/(r2));
    const scalar w23 = w22*z;
    const scalar w24 = k*z;
    const scalar w25 = cos(w24);
    const scalar w26 = hA2*w25;
    const scalar w27 = k*w26;
    const scalar w28 = hA0 + w27*x;
    const scalar w29 = cos(gamma2);
    const scalar w30 = (w14*w14);
    const scalar w31 = w2*w30;
    const scalar w32 = (w13*w13);
    const scalar w33 = w19*w32;
    const scalar w34 = gin*(w29*w29) + w31 + w33;
    const scalar w35 = cos(gamma1);
    const scalar w36 = gin*w0*w29;
    const scalar w37 = sin(gamma1);
    const scalar w38 = w13*w37;
    const scalar w39 = w14*w35;
    const scalar w40 = w29*w39;
    const scalar w41 = w38 + w40;
    const scalar w42 = gil*w41;
    const scalar w43 = w0*w14;
    const scalar w44 = w42*w43;
    const scalar w45 = w14*w37;
    const scalar w46 = w13*w29*w35 - w45;
    const scalar w47 = git*w46;
    const scalar w48 = w0*w13;
    const scalar w49 = w47*w48;
    const scalar w50 = -w35*w36 + w44 + w49;
    const scalar w51 = 1 - 3.0/10.0*w7;
    const scalar w52 = k*y;
    const scalar w53 = sin(w52);
    const scalar w54 = hA1*k*w53;
    const scalar w55 = w51*w54;
    const scalar w56 = sin(w24);
    const scalar w57 = hA1*cos(w52);
    const scalar w58 = (3.0/5.0)*w57;
    const scalar w59 = w58*x;
    const scalar w60 = hA2*w56 - w59;
    const scalar w61 = w29*w45;
    const scalar w62 = w13*w35 - w61;
    const scalar w63 = gil*w62;
    const scalar w64 = w43*w63;
    const scalar w65 = w29*w38;
    const scalar w66 = -w39 - w65;
    const scalar w67 = git*w66;
    const scalar w68 = w48*w67;
    const scalar w69 = w36*w37 + w64 + w68;
    const scalar w70 = w28*w34 - w50*w55 + w60*w69;
    const scalar w71 = w22*w70;
    const scalar w72 = gin*w1;
    const scalar w73 = gil*(w41*w41) + git*(w46*w46) + (w35*w35)*w72;
    const scalar w74 = -w35*w37*w72 + w41*w63 + w46*w67;
    const scalar w75 = w28*w50 - w55*w73 + w60*w74;
    const scalar w76 = w22*w75;
    const scalar w77 = gil*(w62*w62) + git*(w66*w66) + (w37*w37)*w72;
    const scalar w78 = w28*w69 - w55*w74 + w60*w77;
    const scalar w79 = w22*w78;
    const scalar w80 = -w71*z - w76*y - w79*x;
    const scalar w81 = w22*w34;
    const scalar w82 = w22*w50;
    const scalar w83 = w22*w69;
    const scalar w84 = w22*(w81*z + w82*y + w83*x);
    const scalar w85 = w22*w73;
    const scalar w86 = w22*w74;
    const scalar w87 = w22*(w82*z + w85*y + w86*x);
    const scalar w88 = w22*w77;
    const scalar w89 = w22*(w83*z + w86*y + w88*x);
    const scalar w90 = w84*z + w87*y + w89*x;
    const scalar w91 = (1.0/(w90));
    const scalar w92 = r1 + r2;
    const scalar w93 = w6*w92;
    const scalar w94 = w91*w93;
    const scalar w95 = w80*w94;
    const scalar w96 = (1.0/(r1));
    const scalar w97 = w96*z;
    const scalar w98 = w34*w96;
    const scalar w99 = w50*w96;
    const scalar w100 = w69*w96;
    const scalar w101 = w96*(-w100*x - w98*z - w99*y);
    const scalar w102 = w73*w96;
    const scalar w103 = w74*w96;
    const scalar w104 = w96*(-w102*y - w103*x - w99*z);
    const scalar w105 = w77*w96;
    const scalar w106 = w96*(-w100*z - w103*y - w105*x);
    const scalar w107 = -w101*z - w104*y - w106*x;
    const scalar w108 = (1.0/(w107));
    const scalar w109 = w70*w96;
    const scalar w110 = w75*w96;
    const scalar w111 = w78*w96;
    const scalar w112 = w109*z + w110*y + w111*x;
    const scalar w113 = w112*w93;
    const scalar w114 = w108*w113;
    const scalar w115 = (k*k);
    const scalar w116 = hA2*w56;
    const scalar w117 = w115*w116;
    const scalar w118 = w117*x;
    const scalar w119 = w27*w69;
    const scalar w120 = 2*w15;
    const scalar w121 = w120*w48;
    const scalar w122 = w121*w42;
    const scalar w123 = w120*w43;
    const scalar w124 = w123*w47;
    const scalar w125 = w120*w45;
    const scalar w126 = 2*theta1*w13*w29*w35*w6*z - w125*z;
    const scalar w127 = gil*w126;
    const scalar w128 = w120*w38;
    const scalar w129 = w120*w40;
    const scalar w130 = -w128*z - w129*z;
    const scalar w131 = git*w130;
    const scalar w132 = w122*z - w124*z + w127*w43 + w131*w48;
    const scalar w133 = w121*w63;
    const scalar w134 = w123*w67;
    const scalar w135 = w120*w39;
    const scalar w136 = w120*w65;
    const scalar w137 = -w135*z - w136*z;
    const scalar w138 = gil*w43;
    const scalar w139 = w13*w35;
    const scalar w140 = w120*w139;
    const scalar w141 = w120*w61;
    const scalar w142 = -w140*z + w141*z;
    const scalar w143 = git*w48;
    const scalar w144 = w133*z - w134*z + w137*w138 + w142*w143;
    const scalar w145 = -w118*w34 + w119 - w132*w55 + w144*w60 + w21*w28;
    const scalar w146 = w16*w39;
    const scalar w147 = w16*w65;
    const scalar w148 = -w146*z - w147*z;
    const scalar w149 = w139*w16;
    const scalar w150 = w16*w61;
    const scalar w151 = -w149*z + w150*z;
    const scalar w152 = w148*w63 + w151*w67;
    const scalar w153 = w126*w63 + w130*w67 + w137*w42 + w142*w47;
    const scalar w154 = -w118*w69 + w144*w28 + w152*w60 - w153*w55 + w27*w77;
    const scalar w155 = w154*w22;
    const scalar w156 = w16*w45;
    const scalar w157 = 4*theta1*w13*w29*w35*w6*z - w156*z;
    const scalar w158 = w16*w38;
    const scalar w159 = w16*w40;
    const scalar w160 = -w158*z - w159*z;
    const scalar w161 = w157*w42 + w160*w47;
    const scalar w162 = -w118*w50 + w132*w28 + w153*w60 - w161*w55 + w27*w74;
    const scalar w163 = w162*w22;
    const scalar w164 = -w145*w23 - w155*x - w163*y - w71;
    const scalar w165 = w11*w22*w92;
    const scalar w166 = (1.0/2.0)*w91;
    const scalar w167 = w165*w166;
    const scalar w168 = w154*w96;
    const scalar w169 = w162*w96;
    const scalar w170 = w109 + w145*w97 + w168*x + w169*y;
    const scalar w171 = w92*(1 - w11);
    const scalar w172 = w108*w96;
    const scalar w173 = (1.0/2.0)*w172;
    const scalar w174 = w171*w173;
    const scalar w175 = (1.0/(w90*w90));
    const scalar w176 = w132*w22;
    const scalar w177 = w144*w22;
    const scalar w178 = w176*y + w177*x + w21*w23 + w81;
    const scalar w179 = w178*w23;
    const scalar w180 = w161*w22;
    const scalar w181 = w153*w22;
    const scalar w182 = w22*(w176*z + w180*y + w181*x + w82);
    const scalar w183 = w182*y;
    const scalar w184 = w152*w22;
    const scalar w185 = w22*(w177*z + w181*y + w184*x + w83);
    const scalar w186 = w185*x;
    const scalar w187 = -w179 - w183 - w186 - w84;
    const scalar w188 = w175*w187;
    const scalar w189 = w22*w80;
    const scalar w190 = w11*w189*w92;
    const scalar w191 = (1.0/2.0)*w190;
    const scalar w192 = w132*w96;
    const scalar w193 = w144*w96;
    const scalar w194 = -w192*y - w193*x - w21*w97 - w98;
    const scalar w195 = w194*w97;
    const scalar w196 = w161*w96;
    const scalar w197 = w153*w96;
    const scalar w198 = w96*(-w192*z - w196*y - w197*x - w99);
    const scalar w199 = w198*y;
    const scalar w200 = w152*w96;
    const scalar w201 = w96*(-w100 - w193*z - w197*y - w200*x);
    const scalar w202 = w201*x;
    const scalar w203 = w101 + w195 + w199 + w202;
    const scalar w204 = (1.0/(w107*w107));
    const scalar w205 = w204*w96;
    const scalar w206 = w171*w205;
    const scalar w207 = w203*w206;
    const scalar w208 = (1.0/2.0)*w112;
    const scalar w209 = -w114*w97 + w164*w167 + w170*w174 + w188*w191 + w207*w208 + w23*w95;
    const scalar w210 = w10 + w4;
    const scalar w211 = (1.0/(w92));
    const scalar w212 = w11*w211;
    const scalar w213 = w210*w212;
    const scalar w214 = w112*w174 + w166*w190;
    const scalar w215 = 2*w214;
    const scalar w216 = w212*w215;
    const scalar w217 = w211*w6;
    const scalar w218 = w210*w217;
    const scalar w219 = w215*w218;
    const scalar w220 = w209*w213 + w216*z + w219*z + w28;
    const scalar w221 = 2*theta1*w13*w29*w35*w6*y - w125*y;
    const scalar w222 = gil*w221;
    const scalar w223 = -w128*y - w129*y;
    const scalar w224 = git*w223;
    const scalar w225 = w122*y - w124*y + w222*w43 + w224*w48;
    const scalar w226 = -w135*x - w136*x;
    const scalar w227 = -w140*x + w141*x;
    const scalar w228 = w133*x - w134*x + w138*w226 + w143*w227;
    const scalar w229 = 4*theta1*w13*w29*w35*w6*y - w156*y;
    const scalar w230 = -w158*y - w159*y;
    const scalar w231 = w229*w42 + w230*w47;
    const scalar w232 = w22*y;
    const scalar w233 = w96*y;
    const scalar w234 = w18*y - w20*y;
    const scalar w235 = (3.0/5.0)*w54;
    const scalar w236 = w69*x;
    const scalar w237 = w115*w51*w57;
    const scalar w238 = -w135*y - w136*y;
    const scalar w239 = -w140*y + w141*y;
    const scalar w240 = w133*y - w134*y + w138*w238 + w143*w239;
    const scalar w241 = -w225*w55 + w234*w28 + w235*w236 - w237*w50 + w240*w60;
    const scalar w242 = w22*w241;
    const scalar w243 = w235*w74;
    const scalar w244 = w243*x;
    const scalar w245 = w221*w63 + w223*w67 + w238*w42 + w239*w47;
    const scalar w246 = w225*w28 - w231*w55 - w237*w73 + w244 + w245*w60;
    const scalar w247 = w235*x;
    const scalar w248 = -w146*y - w147*y;
    const scalar w249 = -w149*y + w150*y;
    const scalar w250 = w248*w63 + w249*w67;
    const scalar w251 = -w237*w74 + w240*w28 - w245*w55 + w247*w77 + w250*w60;
    const scalar w252 = w22*w251;
    const scalar w253 = -w232*w246 - w242*z - w252*x - w76;
    const scalar w254 = w241*w96;
    const scalar w255 = w251*w96;
    const scalar w256 = w110 + w233*w246 + w254*z + w255*x;
    const scalar w257 = w22*w234;
    const scalar w258 = w22*w225;
    const scalar w259 = w22*w240;
    const scalar w260 = w22*(w257*z + w258*y + w259*x + w82);
    const scalar w261 = w260*z;
    const scalar w262 = w22*w245;
    const scalar w263 = w231*w232 + w258*z + w262*x + w85;
    const scalar w264 = w232*w263;
    const scalar w265 = w22*w250;
    const scalar w266 = w22*(w259*z + w262*y + w265*x + w86);
    const scalar w267 = w266*x;
    const scalar w268 = -w261 - w264 - w267 - w87;
    const scalar w269 = w175*w268;
    const scalar w270 = w234*w96;
    const scalar w271 = w225*w96;
    const scalar w272 = w240*w96;
    const scalar w273 = w96*(-w270*z - w271*y - w272*x - w99);
    const scalar w274 = w273*z;
    const scalar w275 = w245*w96;
    const scalar w276 = -w102 - w231*w233 - w271*z - w275*x;
    const scalar w277 = w233*w276;
    const scalar w278 = w250*w96;
    const scalar w279 = w96*(-w103 - w272*z - w275*y - w278*x);
    const scalar w280 = w279*x;
    const scalar w281 = w104 + w274 + w277 + w280;
    const scalar w282 = w206*w281;
    const scalar w283 = -w114*w233 + w167*w253 + w174*w256 + w191*w269 + w208*w282 + w232*w95;
    const scalar w284 = w213*w283 + w216*y + w219*y - w55;
    const scalar w285 = 2*theta1*w13*w29*w35*w6*x - w125*x;
    const scalar w286 = -w128*x - w129*x;
    const scalar w287 = w226*w42 + w227*w47 + w285*w63 + w286*w67;
    const scalar w288 = -w146*x - w147*x;
    const scalar w289 = -w149*x + w150*x;
    const scalar w290 = w288*w63 + w289*w67;
    const scalar w291 = w22*x;
    const scalar w292 = w96*x;
    const scalar w293 = w18*x - w20*x;
    const scalar w294 = w235*w50;
    const scalar w295 = gil*w285;
    const scalar w296 = git*w286;
    const scalar w297 = w122*x - w124*x + w295*w43 + w296*w48;
    const scalar w298 = w228*w60 + w27*w34 + w28*w293 + w294*x - w297*w55 - w58*w69;
    const scalar w299 = w22*w298;
    const scalar w300 = w235*w73;
    const scalar w301 = 4*theta1*w13*w29*w35*w6*x - w156*x;
    const scalar w302 = -w158*x - w159*x;
    const scalar w303 = w301*w42 + w302*w47;
    const scalar w304 = w27*w50 + w28*w297 + w287*w60 + w300*x - w303*w55 - w58*w74;
    const scalar w305 = w22*w304;
    const scalar w306 = w119 + w228*w28 + w244 - w287*w55 + w290*w60 - w58*w77;
    const scalar w307 = -w291*w306 - w299*z - w305*y - w79;
    const scalar w308 = w165*w307;
    const scalar w309 = w298*w96;
    const scalar w310 = w304*w96;
    const scalar w311 = w111 + w292*w306 + w309*z + w310*y;
    const scalar w312 = w171*w311;
    const scalar w313 = w22*w293;
    const scalar w314 = w22*w297;
    const scalar w315 = w22*w228;
    const scalar w316 = w22*(w313*z + w314*y + w315*x + w83);
    const scalar w317 = w316*z;
    const scalar w318 = w22*w287;
    const scalar w319 = w290*w291 + w315*z + w318*y + w88;
    const scalar w320 = w291*w319;
    const scalar w321 = w22*w303;
    const scalar w322 = w22*(w314*z + w318*x + w321*y + w86);
    const scalar w323 = w322*y;
    const scalar w324 = -w317 - w320 - w323 - w89;
    const scalar w325 = w175*w324;
    const scalar w326 = w293*w96;
    const scalar w327 = w297*w96;
    const scalar w328 = w228*w96;
    const scalar w329 = w96*(-w100 - w326*z - w327*y - w328*x);
    const scalar w330 = w329*z;
    const scalar w331 = w287*w96;
    const scalar w332 = -w105 - w290*w292 - w328*z - w331*y;
    const scalar w333 = w292*w332;
    const scalar w334 = w303*w96;
    const scalar w335 = w96*(-w103 - w327*z - w331*x - w334*y);
    const scalar w336 = w335*y;
    const scalar w337 = w106 + w330 + w333 + w336;
    const scalar w338 = w205*w337;
    const scalar w339 = w171*w338;
    const scalar w340 = -w114*w292 + w166*w308 + w173*w312 + w191*w325 + w208*w339 + w291*w95;
    const scalar w341 = w116 + w213*w340 + w216*x + w219*x - w59;
    const scalar w342 = w214*w217;
    const scalar w343 = 8*w342;
    const scalar w344 = 2*w96;
    const scalar w345 = (k*k*k);
    const scalar w346 = w26*w345;
    const scalar w347 = w346*x;
    const scalar w348 = 2*w118;
    const scalar w349 = w117*w69;
    const scalar w350 = (theta1*theta1);
    const scalar w351 = (1.0/(w5*w5));
    const scalar w352 = w350*w351;
    const scalar w353 = 8*w9;
    const scalar w354 = w352*w353;
    const scalar w355 = -4*gil*theta1*w1*w13*w14*w6 + w20;
    const scalar w356 = 8*gil*w1*w32*w350*w351*w9 + 8*git*w1*w30*w350*w351*w9 - w31*w354 - w33*w354 - w355;
    const scalar w357 = 4*w9;
    const scalar w358 = w352*w357;
    const scalar w359 = -2*theta1*w13*w29*w35*w6 + w125;
    const scalar w360 = -w358*w38 - w358*w40 - w359;
    const scalar w361 = w139*w29;
    const scalar w362 = w128 + w129;
    const scalar w363 = 4*w14*w350*w351*w37*w9 - w358*w361 - w362;
    const scalar w364 = w16*z;
    const scalar w365 = -2*gil*theta1*w0*w13*w41*w6 + w124;
    const scalar w366 = 4*gil*theta1*w0*w126*w13*w6*z + gil*w0*w14*w360 + git*w0*w13*w363 - w131*w364*w43 - w358*w44 - w358*w49 - w365;
    const scalar w367 = w139*w352;
    const scalar w368 = w135 + w136;
    const scalar w369 = 4*w14*w29*w350*w351*w37*w9 - w357*w367 - w368;
    const scalar w370 = w352*w39;
    const scalar w371 = -w140 + w141;
    const scalar w372 = w357*w370 + w358*w65 + w371;
    const scalar w373 = git*w43;
    const scalar w374 = -2*gil*theta1*w0*w13*w6*w62 + w134;
    const scalar w375 = 4*gil*theta1*w0*w13*w137*w6*z + gil*w0*w14*w369 + git*w0*w13*w372 - w142*w364*w373 - w358*w64 - w358*w68 - w374;
    const scalar w376 = 2*hA2*k*w144*w25 - w21*w348 + w28*w356 - w34*w347 - w349 - w366*w55 + w375*w60;
    const scalar w377 = gil*w148;
    const scalar w378 = git*w151;
    const scalar w379 = w146 + w147;
    const scalar w380 = -w149 + w150;
    const scalar w381 = w137*w377 + w142*w378 + w63*(8*w14*w29*w350*w351*w37*w9 - w353*w367 - w379) + w67*(w353*w370 + w354*w65 + w380);
    const scalar w382 = 2*w142;
    const scalar w383 = 2*w127*w137 + w131*w382 + w360*w63 + w363*w67 + w369*w42 + w372*w47;
    const scalar w384 = 2*hA2*k*w152*w25 - w117*w77 - w144*w348 - w236*w346 + w28*w375 + w381*w60 - w383*w55;
    const scalar w385 = -4*theta1*w13*w29*w35*w6 + w156;
    const scalar w386 = w158 + w159;
    const scalar w387 = w127*w157 + w131*w160 + w42*(-w354*w38 - w354*w40 - w385) + w47*(8*w14*w350*w351*w37*w9 - w354*w361 - w386);
    const scalar w388 = 2*hA2*k*w153*w25 - w117*w74 - w132*w348 + w28*w366 - w347*w50 + w383*w60 - w387*w55;
    const scalar w389 = w206*w208;
    const scalar w390 = 2*w22;
    const scalar w391 = w108*w93;
    const scalar w392 = w391*w97;
    const scalar w393 = w23*w94;
    const scalar w394 = w165*w175;
    const scalar w395 = w187*w394;
    const scalar w396 = 2*w101 + 2*w195 + 2*w199 + 2*w202;
    const scalar w397 = w171*w208*w96/(w107*w107*w107);
    const scalar w398 = w203*w397;
    const scalar w399 = w175*w191;
    const scalar w400 = w113*w204;
    const scalar w401 = w400*w97;
    const scalar w402 = w80*w93;
    const scalar w403 = w23*w402;
    const scalar w404 = -2*w179 - 2*w183 - 2*w186 - 2*w84;
    const scalar w405 = w191/(w90*w90*w90);
    const scalar w406 = w187*w405;
    const scalar w407 = -w113*w172 + w189*w94;
    const scalar w408 = 4*w212;
    const scalar w409 = w209*z;
    const scalar w410 = 4*w218;
    const scalar w411 = w216 + w219;
    const scalar w412 = w115*w59;
    const scalar w413 = hA1*w345*w51*w53;
    const scalar w414 = 8*w8;
    const scalar w415 = w352*w414;
    const scalar w416 = 8*gil*w1*w32*w350*w351*w8 + 8*git*w1*w30*w350*w351*w8 - w31*w415 - w33*w415 - w355;
    const scalar w417 = (6.0/5.0)*w54*x;
    const scalar w418 = 2*w237;
    const scalar w419 = 4*w8;
    const scalar w420 = w352*w419;
    const scalar w421 = -w359 - w38*w420 - w40*w420;
    const scalar w422 = 4*w14*w350*w351*w37*w8 - w361*w420 - w362;
    const scalar w423 = w16*y;
    const scalar w424 = 4*gil*theta1*w0*w13*w221*w6*y + gil*w0*w14*w421 + git*w0*w13*w422 - w224*w423*w43 - w365 - w420*w44 - w420*w49;
    const scalar w425 = 4*w14*w29*w350*w351*w37*w8 - w367*w419 - w368;
    const scalar w426 = w370*w419 + w371 + w420*w65;
    const scalar w427 = 4*gil*theta1*w0*w13*w238*w6*y + gil*w0*w14*w425 + git*w0*w13*w426 - w239*w373*w423 - w374 - w420*w64 - w420*w68;
    const scalar w428 = -w225*w418 + w240*w417 + w28*w416 + w412*w69 + w413*w50 - w424*w55 + w427*w60;
    const scalar w429 = w412*w74;
    const scalar w430 = w222*w229 + w224*w230 + w42*(-w38*w415 - w385 - w40*w415) + w47*(8*w14*w350*w351*w37*w8 - w361*w415 - w386);
    const scalar w431 = 2*w239;
    const scalar w432 = 2*w222*w238 + w224*w431 + w42*w425 + w421*w63 + w422*w67 + w426*w47;
    const scalar w433 = -w231*w418 + w245*w417 + w28*w424 + w413*w73 + w429 - w430*w55 + w432*w60;
    const scalar w434 = gil*w248;
    const scalar w435 = git*w249;
    const scalar w436 = w238*w434 + w239*w435 + w63*(8*w14*w29*w350*w351*w37*w8 - w367*w414 - w379) + w67*(w370*w414 + w380 + w415*w65);
    const scalar w437 = -w245*w418 + w250*w417 + w28*w427 + w412*w77 + w413*w74 - w432*w55 + w436*w60;
    const scalar w438 = w233*w391;
    const scalar w439 = w232*w94;
    const scalar w440 = w268*w394;
    const scalar w441 = 2*w104 + 2*w274 + 2*w277 + 2*w280;
    const scalar w442 = w281*w397;
    const scalar w443 = w233*w400;
    const scalar w444 = w232*w402;
    const scalar w445 = -2*w261 - 2*w264 - 2*w267 - 2*w87;
    const scalar w446 = w268*w405;
    const scalar w447 = w283*y;
    const scalar w448 = 8*w7;
    const scalar w449 = 2*w27;
    const scalar w450 = w352*w448;
    const scalar w451 = 8*gil*w1*w32*w350*w351*w7 + 8*git*w1*w30*w350*w351*w7 - w31*w450 - w33*w450 - w355;
    const scalar w452 = (6.0/5.0)*w57;
    const scalar w453 = 4*w7;
    const scalar w454 = w352*w453;
    const scalar w455 = -w359 - w38*w454 - w40*w454;
    const scalar w456 = 4*w14*w350*w351*w37*w7 - w361*w454 - w362;
    const scalar w457 = w16*x;
    const scalar w458 = 4*gil*theta1*w0*w13*w285*w6*x + gil*w0*w14*w455 + git*w0*w13*w456 - w296*w43*w457 - w365 - w44*w454 - w454*w49;
    const scalar w459 = 4*w14*w29*w350*w351*w37*w7 - w367*w453 - w368;
    const scalar w460 = w370*w453 + w371 + w454*w65;
    const scalar w461 = 4*gil*theta1*w0*w13*w226*w6*x + gil*w0*w14*w459 + git*w0*w13*w460 - w227*w373*w457 - w374 - w454*w64 - w454*w68;
    const scalar w462 = -w228*w452 + w28*w451 + w293*w449 + w294 + w297*w417 - w458*w55 + w461*w60;
    const scalar w463 = w295*w301 + w296*w302 + w42*(-w38*w450 - w385 - w40*w450) + w47*(8*w14*w350*w351*w37*w7 - w361*w450 - w386);
    const scalar w464 = 2*w227;
    const scalar w465 = 2*w226*w295 + w296*w464 + w42*w459 + w455*w63 + w456*w67 + w460*w47;
    const scalar w466 = w28*w458 - w287*w452 + w297*w449 + w300 + w303*w417 - w463*w55 + w465*w60;
    const scalar w467 = w228*w27;
    const scalar w468 = gil*w288;
    const scalar w469 = git*w289;
    const scalar w470 = w226*w468 + w227*w469 + w63*(8*w14*w29*w350*w351*w37*w7 - w367*w448 - w379) + w67*(w370*w448 + w380 + w450*w65);
    const scalar w471 = w243 + w28*w461 + w287*w417 - w290*w452 - w465*w55 + 2*w467 + w470*w60;
    const scalar w472 = w292*w391;
    const scalar w473 = w291*w94;
    const scalar w474 = 2*w106 + 2*w330 + 2*w333 + 2*w336;
    const scalar w475 = w337*w397;
    const scalar w476 = w292*w400;
    const scalar w477 = w291*w402;
    const scalar w478 = -2*w317 - 2*w320 - 2*w323 - 2*w89;
    const scalar w479 = w324*w405;
    const scalar w480 = w340*x;
    const scalar w481 = y*z;
    const scalar w482 = 8*w481;
    const scalar w483 = w352*w482;
    const scalar w484 = w370*w481;
    const scalar w485 = w63*(-w367*w482 + w483*w61) + w67*(w483*w65 + 8*w484);
    const scalar w486 = w137*w434 + w142*w435 + w485;
    const scalar w487 = 4*w481;
    const scalar w488 = w352*w487;
    const scalar w489 = 4*w367;
    const scalar w490 = -w481*w489 + w488*w61;
    const scalar w491 = 4*w484 + w488*w65;
    const scalar w492 = w15*w373;
    const scalar w493 = 2*gil*theta1*w0*w13*w137*w6*y + 2*gil*theta1*w0*w13*w238*w6*z + gil*w0*w14*w490 + git*w0*w13*w491 - w382*w492*y - w431*w492*z - w488*w64 - w488*w68;
    const scalar w494 = w352*w38;
    const scalar w495 = -w40*w488 - w487*w494;
    const scalar w496 = w352*w45;
    const scalar w497 = -w361*w488 + w487*w496;
    const scalar w498 = w127*w238 + w131*w239 + w137*w222 + w142*w224 + w42*w490 + w47*w491 + w495*w63 + w497*w67;
    const scalar w499 = -w118*w240 + w152*w247 - w153*w237 + w250*w27 + w28*w493 - w498*w55;
    const scalar w500 = w486*w60 + w499;
    const scalar w501 = w42*(-w40*w483 - w482*w494) + w47*(-w361*w483 + w482*w496);
    const scalar w502 = w127*w229 + w131*w230 + w501;
    const scalar w503 = 2*gil*theta1*w0*w126*w13*w6*y + 2*gil*theta1*w0*w13*w221*w6*z + gil*w0*w14*w495 + git*w0*w13*w497 - w123*w131*y - w123*w224*z - w44*w488 - w488*w49;
    const scalar w504 = w153*w247;
    const scalar w505 = -w118*w225 - w161*w237 + w245*w27 + w28*w503 + w498*w60 + w504;
    const scalar w506 = -w502*w55 + w505;
    const scalar w507 = 8*gil*w1*w32*w350*w351*y*z + 8*git*w1*w30*w350*w351*y*z - w31*w483 - w33*w483;
    const scalar w508 = w240*w27;
    const scalar w509 = -w118*w234 - w132*w237 + w144*w247 + w28*w507 + w493*w60 - w503*w55 + w508;
    const scalar w510 = w169 + w254 + w509*w97;
    const scalar w511 = w197 + w233*w498 + w272 + w493*w97;
    const scalar w512 = w196 + w271 + w292*w498 + w503*w97;
    const scalar w513 = w198 + w273 + w97*(-w192 - w233*w503 - w270 - w292*w493 - w507*w97);
    const scalar w514 = w163 + w23*w509 + w242;
    const scalar w515 = w181 + w23*w493 + w232*w498 + w259;
    const scalar w516 = w180 + w23*w503 + w258 + w291*w498;
    const scalar w517 = w182 + w23*(w176 + w23*w507 + w232*w503 + w257 + w291*w493) + w260;
    const scalar w518 = (1.0/2.0)*w164;
    const scalar w519 = (1.0/2.0)*w253;
    const scalar w520 = w164*w439 + (1.0/2.0)*w170*w282 - w170*w438 + w188*w444 - w203*w443 + (1.0/2.0)*w207*w256 + w253*w393 - w256*w392 + w269*w403 - w281*w401 + w395*w519 + w440*w518;
    const scalar w521 = 2*w212;
    const scalar w522 = w209*y;
    const scalar w523 = 2*w218;
    const scalar w524 = w283*z;
    const scalar w525 = w343*w481 + w521*w522 + w521*w524 + w522*w523 + w523*w524;
    const scalar w526 = w238*w377 + w239*w378 + w485;
    const scalar w527 = w499 + w526*w60;
    const scalar w528 = w157*w222 + w160*w224 + w501;
    const scalar w529 = w505 - w528*w55;
    const scalar w530 = x*z;
    const scalar w531 = 8*w530;
    const scalar w532 = w352*w531;
    const scalar w533 = w370*w530;
    const scalar w534 = w63*(-w367*w531 + w532*w61) + w67*(w532*w65 + 8*w533);
    const scalar w535 = w137*w468 + w142*w469 + w534;
    const scalar w536 = 4*w530;
    const scalar w537 = w352*w536;
    const scalar w538 = -w489*w530 + w537*w61;
    const scalar w539 = 4*w533 + w537*w65;
    const scalar w540 = w492*x;
    const scalar w541 = w464*w492;
    const scalar w542 = 2*gil*theta1*w0*w13*w137*w6*x + 2*gil*theta1*w0*w13*w226*w6*z + gil*w0*w14*w538 + git*w0*w13*w539 - w382*w540 - w537*w64 - w537*w68 - w541*z;
    const scalar w543 = -w40*w537 - w494*w536;
    const scalar w544 = -w361*w537 + w496*w536;
    const scalar w545 = w127*w226 + w131*w227 + w137*w295 + w142*w296 + w42*w538 + w47*w539 + w543*w63 + w544*w67;
    const scalar w546 = -w118*w228 + w144*w27 - w152*w58 + w27*w290 + w28*w542 - w349 + w504 - w545*w55;
    const scalar w547 = w535*w60 + w546;
    const scalar w548 = w42*(-w40*w532 - w494*w531) + w47*(-w361*w532 + w496*w531);
    const scalar w549 = w127*w301 + w131*w302 + w548;
    const scalar w550 = w123*x;
    const scalar w551 = w123*w296;
    const scalar w552 = 2*gil*theta1*w0*w126*w13*w6*x + 2*gil*theta1*w0*w13*w285*w6*z + gil*w0*w14*w543 + git*w0*w13*w544 - w131*w550 - w44*w537 - w49*w537 - w551*z;
    const scalar w553 = -w117*w50 - w118*w297 + w132*w27 - w153*w58 + w161*w247 + w27*w287 + w28*w552 + w545*w60;
    const scalar w554 = -w549*w55 + w553;
    const scalar w555 = 8*gil*w1*w32*w350*w351*x*z + 8*git*w1*w30*w350*w351*x*z - w31*w532 - w33*w532;
    const scalar w556 = -w117*w34 - w118*w293 + w132*w247 - w144*w58 + w21*w27 + w28*w555 + w467 + w542*w60 - w55*w552;
    const scalar w557 = w168 + w309 + w556*w97;
    const scalar w558 = w200 + w233*w545 + w328 + w542*w97;
    const scalar w559 = w197 + w292*w545 + w327 + w552*w97;
    const scalar w560 = w201 + w329 + w97*(-w193 - w233*w552 - w292*w542 - w326 - w555*w97);
    const scalar w561 = w155 + w23*w556 + w299;
    const scalar w562 = w184 + w23*w542 + w232*w545 + w315;
    const scalar w563 = w181 + w23*w552 + w291*w545 + w314;
    const scalar w564 = w185 + w23*(w177 + w23*w555 + w232*w552 + w291*w542 + w313) + w316;
    const scalar w565 = w307*w94;
    const scalar w566 = (1.0/2.0)*w339;
    const scalar w567 = (1.0/2.0)*w205*w312;
    const scalar w568 = w311*w391;
    const scalar w569 = w325*w402;
    const scalar w570 = w165*w325;
    const scalar w571 = (1.0/2.0)*w308;
    const scalar w572 = w337*w400;
    const scalar w573 = w164*w473 - w170*w472 + w170*w566 + w188*w477 + w188*w571 - w203*w476 + w203*w567 + w23*w565 + w23*w569 + w518*w570 - w568*w97 - w572*w97;
    const scalar w574 = w209*x;
    const scalar w575 = w340*z;
    const scalar w576 = w343*x;
    const scalar w577 = w27 + w521*w574 + w521*w575 + w523*w574 + w523*w575 + w576*z;
    const scalar w578 = w226*w377 + w227*w378 + w534;
    const scalar w579 = w546 + w578*w60;
    const scalar w580 = w157*w295 + w160*w296 + w548;
    const scalar w581 = -w55*w580 + w553;
    const scalar w582 = x*y;
    const scalar w583 = 8*w582;
    const scalar w584 = w352*w583;
    const scalar w585 = w370*w582;
    const scalar w586 = w63*(-w367*w583 + w584*w61) + w67*(w584*w65 + 8*w585);
    const scalar w587 = w238*w468 + w239*w469 + w586;
    const scalar w588 = 4*w582;
    const scalar w589 = w352*w588;
    const scalar w590 = -w489*w582 + w589*w61;
    const scalar w591 = 4*w585 + w589*w65;
    const scalar w592 = 2*gil*theta1*w0*w13*w226*w6*y + 2*gil*theta1*w0*w13*w238*w6*x + gil*w0*w14*w590 + git*w0*w13*w591 - w431*w540 - w541*y - w589*w64 - w589*w68;
    const scalar w593 = -w40*w589 - w494*w588;
    const scalar w594 = -w361*w589 + w496*w588;
    const scalar w595 = w222*w226 + w224*w227 + w238*w295 + w239*w296 + w42*w590 + w47*w591 + w593*w63 + w594*w67;
    const scalar w596 = w235*w77 - w237*w287 + w245*w247 + w247*w290 - w250*w58 + w28*w592 + w429 + w508 - w55*w595;
    const scalar w597 = w587*w60 + w596;
    const scalar w598 = w42*(-w40*w584 - w494*w583) + w47*(-w361*w584 + w496*w583);
    const scalar w599 = w222*w301 + w224*w302 + w598;
    const scalar w600 = 2*gil*theta1*w0*w13*w221*w6*x + 2*gil*theta1*w0*w13*w285*w6*y + gil*w0*w14*w593 + git*w0*w13*w594 - w224*w550 - w44*w589 - w49*w589 - w551*y;
    const scalar w601 = w225*w27 + w231*w247 - w237*w303 + w243 - w245*w58 + w247*w287 + w28*w600 + w412*w73 + w595*w60;
    const scalar w602 = -w55*w599 + w601;
    const scalar w603 = 8*gil*w1*w32*w350*w351*x*y + 8*git*w1*w30*w350*w351*x*y - w31*w584 - w33*w584;
    const scalar w604 = w225*w247 + w228*w247 + w234*w27 + w235*w69 - w237*w297 - w240*w58 + w28*w603 + w412*w50 - w55*w600 + w592*w60;
    const scalar w605 = w255 + w310 + w604*w97;
    const scalar w606 = w233*w595 + w278 + w331 + w592*w97;
    const scalar w607 = w275 + w292*w595 + w334 + w600*w97;
    const scalar w608 = w279 + w335 + w97*(-w233*w600 - w272 - w292*w592 - w327 - w603*w97);
    const scalar w609 = w23*w604 + w252 + w305;
    const scalar w610 = w23*w592 + w232*w595 + w265 + w318;
    const scalar w611 = w23*w600 + w262 + w291*w595 + w321;
    const scalar w612 = w23*(w23*w603 + w232*w600 + w259 + w291*w592 + w314) + w266 + w322;
    const scalar w613 = w232*w565 + w232*w569 - w233*w568 - w233*w572 + w253*w473 - w256*w472 + w256*w566 + w269*w477 + w269*w571 - w281*w476 + w281*w567 + w519*w570;
    const scalar w614 = w283*x;
    const scalar w615 = w340*y;
    const scalar w616 = w247 + w521*w614 + w521*w615 + w523*w614 + w523*w615 + w576*y;
    const scalar w617 = w226*w434 + w227*w435 + w586;
    const scalar w618 = w596 + w60*w617;
    const scalar w619 = w229*w295 + w230*w296 + w598;
    const scalar w620 = -w55*w619 + w601;
    return w132*w284 + w144*w341 + w21*w220 + w220*w225 + w220*w228 + w231*w284 + w245*w341 + w284*w287 + w290*w341 + w34*(-w118 + w213*(2*w164*w393 + w164*w395 + w167*(-w145*w390 - w23*w376 - w232*w388 - w291*w384) + w170*w207 - 2*w170*w392 + w174*(w145*w344 + w233*w388 + w292*w384 + w376*w97) + 2*w188*w403 - 2*w203*w401 + w389*(w194*w344 + w233*(-2*w192 - w233*w387 - w292*w383 - w366*w97) + w292*(-2*w193 - w233*w383 - w292*w381 - w375*w97) + w97*(-w21*w344 - w233*w366 - w292*w375 - w356*w97)) + w396*w398 + w399*(-w178*w390 - w23*(w21*w390 + w23*w356 + w232*w366 + w291*w375) - w232*(2*w176 + w23*w366 + w232*w387 + w291*w383) - w291*(2*w177 + w23*w375 + w232*w383 + w291*w381)) + w404*w406 + w407) + w343*w9 + w408*w409 + w409*w410 + w411) + w50*(w213*(w167*(-w232*w506 - w291*w500 - w514) + w174*(w233*w506 + w292*w500 + w510) + w389*(w233*(-w233*w502 - w512) + w292*(-w292*w486 - w511) + w513) + w396*w442 + w399*(-w232*(w232*w502 + w516) - w291*(w291*w486 + w515) - w517) + w404*w446 + w520) + w525) + w50*(w213*(w167*(-w232*w529 - w291*w527 - w514) + w174*(w233*w529 + w292*w527 + w510) + w389*(w233*(-w233*w528 - w512) + w292*(-w292*w526 - w511) + w513) + w398*w441 + w399*(-w232*(w232*w528 + w516) - w291*(w291*w526 + w515) - w517) + w406*w445 + w520) + w525) + w69*(w213*(w167*(-w232*w554 - w291*w547 - w561) + w174*(w233*w554 + w292*w547 + w557) + w389*(w233*(-w233*w549 - w559) + w292*(-w292*w535 - w558) + w560) + w396*w475 + w399*(-w232*(w232*w549 + w563) - w291*(w291*w535 + w562) - w564) + w404*w479 + w573) + w577) + w69*(w213*(w167*(-w232*w581 - w291*w579 - w561) + w174*(w233*w581 + w292*w579 + w557) + w389*(w233*(-w233*w580 - w559) + w292*(-w292*w578 - w558) + w560) + w398*w474 + w399*(-w232*(w232*w580 + w563) - w291*(w291*w578 + w562) - w564) + w406*w478 + w573) + w577) + w73*(w213*(w167*(-w23*w428 - w232*w433 - w246*w390 - w291*w437) + w174*(w233*w433 + w246*w344 + w292*w437 + w428*w97) + 2*w253*w439 + w253*w440 + w256*w282 - 2*w256*w438 + 2*w269*w444 - 2*w281*w443 + w389*(w233*(-w231*w344 - w233*w430 - w292*w432 - w424*w97) + w276*w344 + w292*(-w233*w432 - 2*w275 - w292*w436 - w427*w97) + w97*(-w233*w424 - 2*w271 - w292*w427 - w416*w97)) + w399*(-w23*(w23*w416 + w232*w424 + 2*w258 + w291*w427) - w232*(w23*w424 + w231*w390 + w232*w430 + w291*w432) - w263*w390 - w291*(w23*w427 + w232*w432 + 2*w262 + w291*w436)) + w407 + w441*w442 + w445*w446) - w237 + w343*w8 + w408*w447 + w410*w447 + w411) + w74*(w213*(w167*(-w232*w602 - w291*w597 - w609) + w174*(w233*w602 + w292*w597 + w605) + w389*(w233*(-w233*w599 - w607) + w292*(-w292*w587 - w606) + w608) + w399*(-w232*(w232*w599 + w611) - w291*(w291*w587 + w610) - w612) + w441*w475 + w445*w479 + w613) + w616) + w74*(w213*(w167*(-w232*w620 - w291*w618 - w609) + w174*(w233*w620 + w292*w618 + w605) + w389*(w233*(-w233*w619 - w607) + w292*(-w292*w617 - w606) + w608) + w399*(-w232*(w232*w619 + w611) - w291*(w291*w617 + w610) - w612) + w442*w474 + w446*w478 + w613) + w616) + w77*(w213*(w167*(-w23*w462 - w232*w466 - w291*w471 - w306*w390) + w174*(w233*w466 + w292*w471 + w306*w344 + w462*w97) + 2*w307*w473 + w308*w325 - 2*w311*w472 + w312*w338 + 2*w325*w477 - 2*w337*w476 + w389*(w233*(-w233*w463 - w292*w465 - 2*w331 - w458*w97) + w292*(-w233*w465 - w290*w344 - w292*w470 - w461*w97) + w332*w344 + w97*(-w233*w458 - w292*w461 - 2*w328 - w451*w97)) + w399*(-w23*(w23*w451 + w232*w458 + w291*w461 + 2*w315) - w232*(w23*w458 + w232*w463 + w291*w465 + 2*w318) - w291*(w23*w461 + w232*w465 + w290*w390 + w291*w470) - w319*w390) + w407 + w474*w475 + w478*w479) + w342*w448 + w408*w480 + w410*w480 + w411 - w58);
}

inline symmTensor G(const point& pt, const Params& p)
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

inline scalar g(const point& pt, const Params& p)
{
    const scalar k = p.k;
    const scalar x = pt.x();
    const scalar y = pt.y();
    const scalar z = pt.z();

    return (1.0/4.0)*sin(k*x + 1.0/5.0)*cos(k*z)*cos(k*y - 1.0/3.0) + 1;
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
    {0, 0, 0, 1, 0, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, 0.15614150504445989},
    {0, 0, 0, 1, 0, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 0.12394595957227672},
    {0, 0, 0, 1, 0, -0.12, 0.64000000000000001, 0.70999999999999996, 0, -0.23864962461796799},
    {0, 1, 0, 1, 0, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, 0.2227820083505696},
    {0, 1, 0, 1, 0, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 0.058370650948340416},
    {0, 1, 0, 1, 0, -0.12, 0.64000000000000001, 0.70999999999999996, 0, -0.29215379627359872},
    {1, 0, 0, 1, 0, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, 0.049866789927802149},
    {1, 0, 0, 1, 1, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, -0.48297797547781068},
    {1, 0, 0, 1, 2, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, -0.055620940220825317},
    {1, 0, 0, 1, 0, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 0.59994513164926888},
    {1, 0, 0, 1, 1, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 0.52079149391880963},
    {1, 0, 0, 1, 2, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 0.33555481058203529},
    {1, 0, 0, 1, 0, -0.12, 0.64000000000000001, 0.70999999999999996, 0, 0.14466487265817712},
    {1, 0, 0, 1, 1, -0.12, 0.64000000000000001, 0.70999999999999996, 0, -0.61271501746147194},
    {1, 0, 0, 1, 2, -0.12, 0.64000000000000001, 0.70999999999999996, 0, 0.19428310152250294},
    {1, 1, 0, 1, 0, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, 0.051692143322643727},
    {1, 1, 0, 1, 1, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, -0.6390782622162603},
    {1, 1, 0, 1, 2, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, -0.1526392153754175},
    {1, 1, 0, 1, 0, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 0.42872279466482383},
    {1, 1, 0, 1, 1, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 0.72630487658745557},
    {1, 1, 0, 1, 2, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 0.16367075524364474},
    {1, 1, 0, 1, 0, -0.12, 0.64000000000000001, 0.70999999999999996, 0, 0.061928931092193919},
    {1, 1, 0, 1, 1, -0.12, 0.64000000000000001, 0.70999999999999996, 0, -0.65955981348154447},
    {1, 1, 0, 1, 2, -0.12, 0.64000000000000001, 0.70999999999999996, 0, 0.21216184559697684},
    {2, 0, 0, 1, 0, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, -0.36682944027802533},
    {2, 0, 0, 1, 0, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, -0.15680035073312398},
    {2, 0, 0, 1, 0, -0.12, 0.64000000000000001, 0.70999999999999996, 0, 0.98040474129199973},
    {2, 1, 0, 1, 0, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, -0.26624417384532634},
    {2, 1, 0, 1, 0, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, -0.062437208355870133},
    {2, 1, 0, 1, 0, -0.12, 0.64000000000000001, 0.70999999999999996, 0, 0.89310210920615674},
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
    {4, 0, 0, 1, 0, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, 1.1585424265442228},
    {4, 0, 0, 1, 0, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 0.97638240328397263},
    {4, 0, 0, 1, 0, -0.12, 0.64000000000000001, 0.70999999999999996, 0, 0.99713250605483494},
    {4, 1, 0, 1, 0, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, 1.1585424265442228},
    {4, 1, 0, 1, 0, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 0.97638240328397263},
    {4, 1, 0, 1, 0, -0.12, 0.64000000000000001, 0.70999999999999996, 0, 0.99713250605483494},
};

inline const label nSelfTestSamples = 72;

} // End namespace shellMonodomainMMS
} // End namespace Foam

#endif
