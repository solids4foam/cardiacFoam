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
    Exact fields, conductivities, sources and self-test samples of the
    bath-bidomain transmission manufactured solution on the spherical shell.

    Regenerate from generator/:
        python shell_mms.py --write ../../../src/verificationModels/bathBidomainVerification/manufacturedShellTransmissionBathBidomainReference.H

    sympy 1.14.0
\*---------------------------------------------------------------------------*/

#ifndef manufacturedShellTransmissionBathBidomainReference_H
#define manufacturedShellTransmissionBathBidomainReference_H

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
namespace shellTransmissionMMS
{

struct Params
{
    scalar r1, r2, gamma1, gamma2, k;
    scalar gil, git, gin, gel, get, gen, theta0, theta1;
    scalar sigma0In, sigma0Out, sigmaEps;
    scalar fA[4], hA[3];
    scalar cmAmplitude, cmOmega;
    label solution;
};

enum Zone : label { bloodPool = 0, myocardium = 1, bath = 2 };

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
    p.gel = 0.62;
    p.get = 0.12;
    p.gen = 0.08;
    p.theta0 = degToRad(-20.0);
    p.theta1 = degToRad(90.0);
    p.sigma0In = 0.6;
    p.sigma0Out = 0.25;
    p.sigmaEps = 0.2;
    p.cmAmplitude = 0.0;
    p.cmOmega = 0.0;
    p.fA[0] = 0.3;
    p.fA[1] = 0.2;
    p.fA[2] = 0.4;
    p.fA[3] = 0.15;
    p.hA[0] = -0.1;
    p.hA[1] = 0.3;
    p.hA[2] = 0.25;
    p.solution = 1;
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
    p.gel = 0.62;
    p.get = 0.12;
    p.gen = 0.08;
    p.theta0 = degToRad(-20.0);
    p.theta1 = degToRad(90.0);
    p.sigma0In = 0.6;
    p.sigma0Out = 0.25;
    p.sigmaEps = 0.2;
    p.cmAmplitude = 0.3;
    p.cmOmega = 20.0;
    p.fA[0] = 0.25;
    p.fA[1] = 0.25;
    p.fA[2] = 0.35;
    p.fA[3] = 0.2;
    p.hA[0] = -0.15;
    p.hA[1] = 0.35;
    p.hA[2] = 0.15;
    p.solution = 1;
    return p;
}

inline label zoneOf(const point& pt, const Params& p)
{
    const scalar rho = magSqr(pt);
    if (rho < p.r1*p.r1)
    {
        return bloodPool;
    }
    if (rho > p.r2*p.r2)
    {
        return bath;
    }
    return myocardium;
}

inline void checkZone(const label zone, const bool heartAllowed)
{
    if (zone < bloodPool || zone > bath || (!heartAllowed && zone == myocardium))
    {
        FatalErrorInFunction
            << "Zone " << zone << " is not valid here" << exit(FatalError);
    }
}

inline scalar Vm(const point& pt, const scalar t, const Params& p)
{
    const scalar r1 = p.r1;
    const scalar r2 = p.r2;
    const scalar gamma1 = p.gamma1;
    const scalar gamma2 = p.gamma2;
    const scalar k = p.k;
    const scalar gil = p.gil;
    const scalar git = p.git;
    const scalar gin = p.gin;
    const scalar gel = p.gel;
    const scalar get = p.get;
    const scalar gen = p.gen;
    const scalar theta0 = p.theta0;
    const scalar theta1 = p.theta1;
    const scalar sigma0In = p.sigma0In;
    const scalar sigma0Out = p.sigma0Out;
    const scalar sigmaEps = p.sigmaEps;
    const scalar fA0 = p.fA[0];
    const scalar fA1 = p.fA[1];
    const scalar fA2 = p.fA[2];
    const scalar fA3 = p.fA[3];
    const scalar hA0 = p.hA[0];
    const scalar hA1 = p.hA[1];
    const scalar hA2 = p.hA[2];
    const scalar x = pt.x();
    const scalar y = pt.y();
    const scalar z = pt.z();
    const scalar St = sqrt(1 + t);

    const scalar w0 = k*z;
    const scalar w1 = sin(w0);
    const scalar w2 = (x*x);
    const scalar w3 = 1 - 3.0/10.0*w2;
    const scalar w4 = k*y;
    const scalar w5 = cos(w4);
    const scalar w6 = hA1*w5;
    const scalar w7 = fA1*y;
    const scalar w8 = (1.0/2.0)*y + 1;
    const scalar w9 = k*x;
    const scalar w10 = fA2*sin(w9);
    const scalar w11 = w4 + 3.0/10.0;
    const scalar w12 = fA3*exp((2.0/5.0)*z);
    const scalar w13 = w12*cos(w11);
    const scalar w14 = cos(gamma2);
    const scalar w15 = (w14*w14);
    const scalar w16 = sin(gamma2);
    const scalar w17 = (w16*w16);
    const scalar w18 = (r1*r1);
    const scalar w19 = -(r2*r2);
    const scalar w20 = w2 + (y*y) + (z*z);
    const scalar w21 = (-w18 + w20)/(-w18 - w19);
    const scalar w22 = theta0 + theta1*w21;
    const scalar w23 = sin(w22);
    const scalar w24 = w17*(w23*w23);
    const scalar w25 = cos(w22);
    const scalar w26 = w17*(w25*w25);
    const scalar w27 = gil*w24 + gin*w15 + git*w26;
    const scalar w28 = (1.0/(r2));
    const scalar w29 = w28*z;
    const scalar w30 = cos(gamma1);
    const scalar w31 = w14*w16;
    const scalar w32 = gin*w31;
    const scalar w33 = sin(gamma1);
    const scalar w34 = w25*w33;
    const scalar w35 = w23*w30;
    const scalar w36 = w14*w35 + w34;
    const scalar w37 = w16*w23;
    const scalar w38 = gil*w37;
    const scalar w39 = w23*w33;
    const scalar w40 = w14*w25*w30 - w39;
    const scalar w41 = w16*w25;
    const scalar w42 = git*w41;
    const scalar w43 = -w30*w32 + w36*w38 + w40*w42;
    const scalar w44 = w28*y;
    const scalar w45 = -w14*w39 + w25*w30;
    const scalar w46 = -w14*w34 - w35;
    const scalar w47 = w32*w33 + w38*w45 + w42*w46;
    const scalar w48 = w28*x;
    const scalar w49 = (w30*w30);
    const scalar w50 = gin*w17;
    const scalar w51 = (w36*w36);
    const scalar w52 = (w40*w40);
    const scalar w53 = gil*w51 + git*w52 + w49*w50;
    const scalar w54 = w30*w33;
    const scalar w55 = w36*w45;
    const scalar w56 = w40*w46;
    const scalar w57 = gil*w55 + git*w56 - w50*w54;
    const scalar w58 = (w33*w33);
    const scalar w59 = (w45*w45);
    const scalar w60 = (w46*w46);
    const scalar w61 = gil*w59 + git*w60 + w50*w58;
    const scalar w62 = St*(hA0 + hA2*w9*cos(w0));
    const scalar w63 = St*hA1*k*w3*sin(w4);
    const scalar w64 = St*(hA2*w1 - 3.0/5.0*w6*x);
    const scalar w65 = w27*w62 - w43*w63 + w47*w64;
    const scalar w66 = w43*w62 - w53*w63 + w57*w64;
    const scalar w67 = w47*w62 - w57*w63 + w61*w64;
    const scalar w68 = r1 + r2;
    const scalar w69 = (1.0/2.0)*w68;
    const scalar w70 = w21*w28*w69;
    const scalar w71 = (1.0/(r1));
    const scalar w72 = w71*z;
    const scalar w73 = w71*y;
    const scalar w74 = w71*x;
    const scalar w75 = w69*w71*(1 - w21);
    const scalar w76 = w21*(w19 + w20)/w68;
    const scalar w77 = gel*w24 + gen*w15 + get*w26;
    const scalar w78 = gen*w31;
    const scalar w79 = gel*w37;
    const scalar w80 = get*w41;
    const scalar w81 = -w30*w78 + w36*w79 + w40*w80;
    const scalar w82 = w33*w78 + w45*w79 + w46*w80;
    const scalar w83 = gen*w17;
    const scalar w84 = gel*w51 + get*w52 + w49*w83;
    const scalar w85 = gel*w55 + get*w56 - w54*w83;
    const scalar w86 = gel*w59 + get*w60 + w58*w83;
    const scalar w87 = sigmaEps*w5*((1.0/2.0)*w1 + 1)*sin(w9 + 1.0/5.0) + 1;
    const scalar w88 = -sigma0Out*w87;
    const scalar w89 = St*((2.0/5.0)*w13 + w7);
    const scalar w90 = -w82;
    const scalar w91 = St*(fA0 + fA2*k*w8*cos(w9));
    const scalar w92 = -w81;
    const scalar w93 = St*(fA1*z - k*w12*sin(w11) + (1.0/2.0)*w10);
    const scalar w94 = w90*w91 + w92*w93;
    const scalar w95 = -w85;
    const scalar w96 = w89*w92 + w91*w95;
    const scalar w97 = w89*w90 + w93*w95;
    const scalar w98 = -sigma0In*w87;
    return St*(hA0*z + hA2*w1*x + w3*w6) - St*(fA0*x + w10*w8 + w13 + w7*z) + w76*(w70*(-w29*w65 - w44*w66 - w48*w67)/(w29*(w27*w29 + w43*w44 + w47*w48) + w44*(w29*w43 + w44*w53 + w48*w57) + w48*(w29*w47 + w44*w57 + w48*w61)) + w75*(w65*w72 + w66*w73 + w67*w74)/(-w72*(-w27*w72 - w43*w73 - w47*w74) - w73*(-w43*w72 - w53*w73 - w57*w74) - w74*(-w47*w72 - w57*w73 - w61*w74))) - w76*(w70*(w29*(w89*(-w77 - w88) + w94) + w44*(w93*(-w84 - w88) + w96) + w48*(w91*(-w86 - w88) + w97))/(w29*(w29*w77 + w44*w81 + w48*w82) + w44*(w29*w81 + w44*w84 + w48*w85) + w48*(w29*w82 + w44*w85 + w48*w86)) + w75*(-w72*(w89*(-w77 - w98) + w94) - w73*(w93*(-w84 - w98) + w96) - w74*(w91*(-w86 - w98) + w97))/(-w72*(-w72*w77 - w73*w81 - w74*w82) - w73*(-w72*w81 - w73*w84 - w74*w85) - w74*(-w72*w82 - w73*w85 - w74*w86)));
}

inline scalar phiE(const point& pt, const scalar t, const label zone, const Params& p)
{
    const scalar r1 = p.r1;
    const scalar r2 = p.r2;
    const scalar gamma1 = p.gamma1;
    const scalar gamma2 = p.gamma2;
    const scalar k = p.k;
    const scalar gil = p.gil;
    const scalar git = p.git;
    const scalar gin = p.gin;
    const scalar gel = p.gel;
    const scalar get = p.get;
    const scalar gen = p.gen;
    const scalar theta0 = p.theta0;
    const scalar theta1 = p.theta1;
    const scalar sigma0In = p.sigma0In;
    const scalar sigma0Out = p.sigma0Out;
    const scalar sigmaEps = p.sigmaEps;
    const scalar fA0 = p.fA[0];
    const scalar fA1 = p.fA[1];
    const scalar fA2 = p.fA[2];
    const scalar fA3 = p.fA[3];
    const scalar cmAmplitude = p.cmAmplitude;
    const scalar cmOmega = p.cmOmega;
    const scalar x = pt.x();
    const scalar y = pt.y();
    const scalar z = pt.z();
    const scalar St = sqrt(1 + t);
    const scalar Tt = cmAmplitude*sin(cmOmega*t);
    checkZone(zone, true);

    if (zone == bloodPool)
    {
        const scalar w0 = k*x;
        const scalar w1 = sin(w0);
        const scalar w2 = k*y;
        const scalar w3 = (1.0/10.0)*y;
        const scalar w4 = (z*z);
        const scalar w5 = (1.0/(r1));
        const scalar w6 = (r1*r1);
        const scalar w7 = -w4 + w6 - (x*x) - (y*y);
        const scalar w8 = sigmaEps*((1.0/2.0)*sin(k*z) + 1)*sin(w0 + 1.0/5.0)*cos(w2) + 1;
        const scalar w9 = cos(gamma1);
        const scalar w10 = cos(gamma2);
        const scalar w11 = sin(gamma2);
        const scalar w12 = w10*w11;
        const scalar w13 = w12*w9;
        const scalar w14 = gen*w13;
        const scalar w15 = sin(gamma1);
        const scalar w16 = cos(theta0);
        const scalar w17 = w15*w16;
        const scalar w18 = sin(theta0);
        const scalar w19 = w18*w9;
        const scalar w20 = w10*w19 + w17;
        const scalar w21 = w11*w18;
        const scalar w22 = gel*w21;
        const scalar w23 = w20*w22;
        const scalar w24 = w15*w18;
        const scalar w25 = w10*w16*w9 - w24;
        const scalar w26 = w11*w16;
        const scalar w27 = get*w26;
        const scalar w28 = w25*w27;
        const scalar w29 = (w10*w10);
        const scalar w30 = (w11*w11);
        const scalar w31 = (w18*w18)*w30;
        const scalar w32 = (w16*w16)*w30;
        const scalar w33 = gel*w31 + gen*w29 + get*w32;
        const scalar w34 = (1.0/5.0)*z;
        const scalar w35 = gen*w15;
        const scalar w36 = -w10*w24 + w16*w9;
        const scalar w37 = -w10*w17 - w19;
        const scalar w38 = w12*w35 + w22*w36 + w27*w37;
        const scalar w39 = (1.0/5.0)*k*w1;
        const scalar w40 = w5*z;
        const scalar w41 = w30*(w9*w9);
        const scalar w42 = gen*w41;
        const scalar w43 = (w20*w20);
        const scalar w44 = gel*w43;
        const scalar w45 = (w25*w25);
        const scalar w46 = get*w45;
        const scalar w47 = -w14 + w23 + w28;
        const scalar w48 = w30*w9;
        const scalar w49 = w35*w48;
        const scalar w50 = w20*w36;
        const scalar w51 = gel*w50;
        const scalar w52 = w25*w37;
        const scalar w53 = get*w52;
        const scalar w54 = -w49 + w51 + w53;
        const scalar w55 = w5*y;
        const scalar w56 = (w15*w15)*w30;
        const scalar w57 = (w36*w36);
        const scalar w58 = (w37*w37);
        const scalar w59 = gel*w57 + gen*w56 + get*w58;
        const scalar w60 = w5*x;
        const scalar w61 = gin*w13;
        const scalar w62 = gil*w21;
        const scalar w63 = w20*w62;
        const scalar w64 = git*w26;
        const scalar w65 = w25*w64;
        const scalar w66 = z*(gil*w31 + gin*w29 + git*w32);
        const scalar w67 = gin*w15;
        const scalar w68 = w12*w67 + w36*w62 + w37*w64;
        const scalar w69 = gin*w41;
        const scalar w70 = gil*w43;
        const scalar w71 = git*w45;
        const scalar w72 = -w61 + w63 + w65;
        const scalar w73 = w72*z;
        const scalar w74 = w48*w67;
        const scalar w75 = gil*w50;
        const scalar w76 = git*w52;
        const scalar w77 = -w74 + w75 + w76;
        const scalar w78 = w68*z;
        const scalar w79 = gil*w57 + gin*w56 + git*w58;
        return St*(fA0*x + fA1*y*z + fA2*w1*((1.0/2.0)*y + 1) + fA3*exp((2.0/5.0)*z)*cos(w2 + 3.0/10.0)) + Tt*(w3 + (1.0/10.0)*w4 + (1.0/5.0)*cos(w0) + (1.0/2.0)*w5*w7*(-sigma0In*w8*((1.0/5.0)*k*w1*w5*x - w3*w5 - 1.0/5.0*w4*w5) - w40*(-1.0/10.0*w14 + (1.0/10.0)*w23 + (1.0/10.0)*w28 + w33*w34 - w38*w39) - w55*(w34*w47 - w39*w54 + (1.0/10.0)*w42 + (1.0/10.0)*w44 + (1.0/10.0)*w46) - w6*(-w40*(-w33*w40 - w38*w60 - w47*w55) - w55*(-w40*w47 - w54*w60 - w55*(w42 + w44 + w46)) - w60*(-w38*w40 - w54*w55 - w59*w60))*(-w40*(-w39*w68 - 1.0/10.0*w61 + (1.0/10.0)*w63 + (1.0/10.0)*w65 + (1.0/5.0)*w66) - w55*(-w39*w77 + (1.0/10.0)*w69 + (1.0/10.0)*w70 + (1.0/10.0)*w71 + (1.0/5.0)*w73) - w60*(-w39*w79 - 1.0/10.0*w74 + (1.0/10.0)*w75 + (1.0/10.0)*w76 + (1.0/5.0)*w78))/(w7/((1.0/(git)) + (1.0/(gin)) + (1.0/(gil))) + x*(w77*y + w78 + w79*x) + y*(w73 + w77*x + y*(w69 + w70 + w71)) + z*(w66 + w68*x + w72*y)) - w60*(w34*w38 - w39*w59 - 1.0/10.0*w49 + (1.0/10.0)*w51 + (1.0/10.0)*w53))/(sigma0In*w8));
    }
    if (zone == bath)
    {
        const scalar w0 = k*x;
        const scalar w1 = sin(w0);
        const scalar w2 = k*y;
        const scalar w3 = (1.0/10.0)*y;
        const scalar w4 = (z*z);
        const scalar w5 = (1.0/(r2));
        const scalar w6 = (r2*r2);
        const scalar w7 = sigmaEps*((1.0/2.0)*sin(k*z) + 1)*sin(w0 + 1.0/5.0)*cos(w2) + 1;
        const scalar w8 = (1.0/5.0)*w5;
        const scalar w9 = cos(gamma1);
        const scalar w10 = cos(gamma2);
        const scalar w11 = sin(gamma2);
        const scalar w12 = w10*w11;
        const scalar w13 = w12*w9;
        const scalar w14 = gen*w13;
        const scalar w15 = sin(gamma1);
        const scalar w16 = theta0 + theta1;
        const scalar w17 = cos(w16);
        const scalar w18 = w15*w17;
        const scalar w19 = sin(w16);
        const scalar w20 = w19*w9;
        const scalar w21 = w10*w20 + w18;
        const scalar w22 = w11*w19;
        const scalar w23 = gel*w22;
        const scalar w24 = w21*w23;
        const scalar w25 = w15*w19;
        const scalar w26 = w10*w17*w9 - w25;
        const scalar w27 = w11*w17;
        const scalar w28 = get*w27;
        const scalar w29 = w26*w28;
        const scalar w30 = (w10*w10);
        const scalar w31 = (w11*w11);
        const scalar w32 = (w19*w19)*w31;
        const scalar w33 = (w17*w17)*w31;
        const scalar w34 = gel*w32 + gen*w30 + get*w33;
        const scalar w35 = (1.0/5.0)*z;
        const scalar w36 = gen*w15;
        const scalar w37 = -w10*w25 + w17*w9;
        const scalar w38 = -w10*w18 - w20;
        const scalar w39 = w12*w36 + w23*w37 + w28*w38;
        const scalar w40 = (1.0/5.0)*k*w1;
        const scalar w41 = w5*z;
        const scalar w42 = w31*(w9*w9);
        const scalar w43 = gen*w42;
        const scalar w44 = (w21*w21);
        const scalar w45 = gel*w44;
        const scalar w46 = (w26*w26);
        const scalar w47 = get*w46;
        const scalar w48 = -w14 + w24 + w29;
        const scalar w49 = w31*w9;
        const scalar w50 = w36*w49;
        const scalar w51 = w21*w37;
        const scalar w52 = gel*w51;
        const scalar w53 = w26*w38;
        const scalar w54 = get*w53;
        const scalar w55 = -w50 + w52 + w54;
        const scalar w56 = w5*y;
        const scalar w57 = (w15*w15)*w31;
        const scalar w58 = (w37*w37);
        const scalar w59 = (w38*w38);
        const scalar w60 = gel*w58 + gen*w57 + get*w59;
        const scalar w61 = w5*x;
        const scalar w62 = z*(gil*w32 + gin*w30 + git*w33);
        const scalar w63 = gin*w13;
        const scalar w64 = gil*w22;
        const scalar w65 = w21*w64;
        const scalar w66 = git*w27;
        const scalar w67 = w26*w66;
        const scalar w68 = -w63 + w65 + w67;
        const scalar w69 = gin*w15;
        const scalar w70 = w12*w69 + w37*w64 + w38*w66;
        const scalar w71 = gin*w42;
        const scalar w72 = gil*w44;
        const scalar w73 = git*w46;
        const scalar w74 = w68*z;
        const scalar w75 = w49*w69;
        const scalar w76 = gil*w51;
        const scalar w77 = git*w53;
        const scalar w78 = -w75 + w76 + w77;
        const scalar w79 = gil*w58 + gin*w57 + git*w59;
        const scalar w80 = w70*z;
        return St*(fA0*x + fA1*y*z + fA2*w1*((1.0/2.0)*y + 1) + fA3*exp((2.0/5.0)*z)*cos(w2 + 3.0/10.0)) + Tt*(w3 + (1.0/10.0)*w4 + (1.0/5.0)*cos(w0) + (1.0/2.0)*w5*(w4 - w6 + (x*x) + (y*y))*(-sigma0Out*w7*(-w0*w1*w8 + w3*w5 + w4*w8) + w41*(-1.0/10.0*w14 + (1.0/10.0)*w24 + (1.0/10.0)*w29 + w34*w35 - w39*w40) + w56*(w35*w48 - w40*w55 + (1.0/10.0)*w43 + (1.0/10.0)*w45 + (1.0/10.0)*w47) - w6*(w41*(w34*w41 + w39*w61 + w48*w56) + w56*(w41*w48 + w55*w61 + w56*(w43 + w45 + w47)) + w61*(w39*w41 + w55*w56 + w60*w61))*(w41*(-w40*w70 + (1.0/5.0)*w62 - 1.0/10.0*w63 + (1.0/10.0)*w65 + (1.0/10.0)*w67) + w56*(-w40*w78 + (1.0/10.0)*w71 + (1.0/10.0)*w72 + (1.0/10.0)*w73 + (1.0/5.0)*w74) + w61*(-w40*w79 - 1.0/10.0*w75 + (1.0/10.0)*w76 + (1.0/10.0)*w77 + (1.0/5.0)*w80))/(x*(w78*y + w79*x + w80) + y*(w74 + w78*x + y*(w71 + w72 + w73)) + z*(w62 + w68*y + w70*x)) + w61*(w35*w39 - w40*w60 - 1.0/10.0*w50 + (1.0/10.0)*w52 + (1.0/10.0)*w54))/(sigma0Out*w7));
    }
    const scalar w0 = fA1*y;
    const scalar w1 = (1.0/2.0)*y + 1;
    const scalar w2 = k*x;
    const scalar w3 = sin(w2);
    const scalar w4 = fA2*w3;
    const scalar w5 = k*y;
    const scalar w6 = w5 + 3.0/10.0;
    const scalar w7 = fA3*exp((2.0/5.0)*z);
    const scalar w8 = w7*cos(w6);
    const scalar w9 = (z*z);
    const scalar w10 = cos(w2);
    const scalar w11 = cos(gamma2);
    const scalar w12 = (w11*w11);
    const scalar w13 = sin(gamma2);
    const scalar w14 = (w13*w13);
    const scalar w15 = (r1*r1);
    const scalar w16 = -(r2*r2);
    const scalar w17 = w9 + (x*x) + (y*y);
    const scalar w18 = (-w15 + w17)/(-w15 - w16);
    const scalar w19 = theta0 + theta1*w18;
    const scalar w20 = sin(w19);
    const scalar w21 = w14*(w20*w20);
    const scalar w22 = cos(w19);
    const scalar w23 = w14*(w22*w22);
    const scalar w24 = gil*w21 + gin*w12 + git*w23;
    const scalar w25 = (1.0/(r2));
    const scalar w26 = w25*z;
    const scalar w27 = cos(gamma1);
    const scalar w28 = w11*w13;
    const scalar w29 = gin*w28;
    const scalar w30 = w27*w29;
    const scalar w31 = sin(gamma1);
    const scalar w32 = w22*w31;
    const scalar w33 = w20*w27;
    const scalar w34 = w11*w33 + w32;
    const scalar w35 = w13*w20;
    const scalar w36 = gil*w35;
    const scalar w37 = w34*w36;
    const scalar w38 = w20*w31;
    const scalar w39 = w11*w22*w27 - w38;
    const scalar w40 = w13*w22;
    const scalar w41 = git*w40;
    const scalar w42 = w39*w41;
    const scalar w43 = -w30 + w37 + w42;
    const scalar w44 = w25*y;
    const scalar w45 = -w11*w38 + w22*w27;
    const scalar w46 = -w11*w32 - w33;
    const scalar w47 = w29*w31 + w36*w45 + w41*w46;
    const scalar w48 = w25*x;
    const scalar w49 = (w27*w27);
    const scalar w50 = gin*w14;
    const scalar w51 = w49*w50;
    const scalar w52 = (w34*w34);
    const scalar w53 = gil*w52;
    const scalar w54 = (w39*w39);
    const scalar w55 = git*w54;
    const scalar w56 = w51 + w53 + w55;
    const scalar w57 = w27*w31;
    const scalar w58 = w50*w57;
    const scalar w59 = w34*w45;
    const scalar w60 = gil*w59;
    const scalar w61 = w39*w46;
    const scalar w62 = git*w61;
    const scalar w63 = -w58 + w60 + w62;
    const scalar w64 = (w31*w31);
    const scalar w65 = (w45*w45);
    const scalar w66 = (w46*w46);
    const scalar w67 = gil*w65 + git*w66 + w50*w64;
    const scalar w68 = (1.0/5.0)*z;
    const scalar w69 = (1.0/5.0)*k*w3;
    const scalar w70 = w24*w68 - 1.0/10.0*w30 + (1.0/10.0)*w37 + (1.0/10.0)*w42 - w47*w69;
    const scalar w71 = w43*w68 + (1.0/10.0)*w51 + (1.0/10.0)*w53 + (1.0/10.0)*w55 - w63*w69;
    const scalar w72 = w47*w68 - 1.0/10.0*w58 + (1.0/10.0)*w60 + (1.0/10.0)*w62 - w67*w69;
    const scalar w73 = r1 + r2;
    const scalar w74 = (1.0/2.0)*w73;
    const scalar w75 = w18*w25*w74;
    const scalar w76 = (1.0/(r1));
    const scalar w77 = w76*z;
    const scalar w78 = w76*y;
    const scalar w79 = w76*x;
    const scalar w80 = w74*w76*(1 - w18);
    const scalar w81 = w18*(w16 + w17)/w73;
    const scalar w82 = gel*w21 + gen*w12 + get*w23;
    const scalar w83 = gen*w28;
    const scalar w84 = gel*w35;
    const scalar w85 = get*w40;
    const scalar w86 = -w27*w83 + w34*w84 + w39*w85;
    const scalar w87 = w31*w83 + w45*w84 + w46*w85;
    const scalar w88 = gen*w14;
    const scalar w89 = gel*w52 + get*w54 + w49*w88;
    const scalar w90 = gel*w59 + get*w61 - w57*w88;
    const scalar w91 = gel*w65 + get*w66 + w64*w88;
    const scalar w92 = sigmaEps*((1.0/2.0)*sin(k*z) + 1)*sin(w2 + 1.0/5.0)*cos(w5) + 1;
    const scalar w93 = -sigma0Out*w92;
    const scalar w94 = St*(w0 + (2.0/5.0)*w8);
    const scalar w95 = -w87;
    const scalar w96 = St*(fA0 + fA2*k*w1*w10);
    const scalar w97 = -w86;
    const scalar w98 = St*(fA1*z - k*w7*sin(w6) + (1.0/2.0)*w4);
    const scalar w99 = w95*w96 + w97*w98;
    const scalar w100 = -w90;
    const scalar w101 = w100*w96 + w94*w97;
    const scalar w102 = w100*w98 + w94*w95;
    const scalar w103 = -sigma0In*w92;
    return St*(fA0*x + w0*z + w1*w4 + w8) + Tt*((1.0/5.0)*w10 + w81*(w75*(-w26*w70 - w44*w71 - w48*w72)/(w26*(w24*w26 + w43*w44 + w47*w48) + w44*(w26*w43 + w44*w56 + w48*w63) + w48*(w26*w47 + w44*w63 + w48*w67)) + w80*(w70*w77 + w71*w78 + w72*w79)/(-w77*(-w24*w77 - w43*w78 - w47*w79) - w78*(-w43*w77 - w56*w78 - w63*w79) - w79*(-w47*w77 - w63*w78 - w67*w79))) + (1.0/10.0)*w9 + (1.0/10.0)*y) + w81*(w75*(w26*(w94*(-w82 - w93) + w99) + w44*(w101 + w98*(-w89 - w93)) + w48*(w102 + w96*(-w91 - w93)))/(w26*(w26*w82 + w44*w86 + w48*w87) + w44*(w26*w86 + w44*w89 + w48*w90) + w48*(w26*w87 + w44*w90 + w48*w91)) + w80*(-w77*(w94*(-w103 - w82) + w99) - w78*(w101 + w98*(-w103 - w89)) - w79*(w102 + w96*(-w103 - w91)))/(-w77*(-w77*w82 - w78*w86 - w79*w87) - w78*(-w77*w86 - w78*w89 - w79*w90) - w79*(-w77*w87 - w78*w90 - w79*w91)));
}

inline scalar phiI(const point& pt, const scalar t, const Params& p)
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
    const scalar cmAmplitude = p.cmAmplitude;
    const scalar cmOmega = p.cmOmega;
    const scalar x = pt.x();
    const scalar y = pt.y();
    const scalar z = pt.z();
    const scalar St = sqrt(1 + t);
    const scalar Tt = cmAmplitude*sin(cmOmega*t);

    const scalar w0 = k*z;
    const scalar w1 = sin(w0);
    const scalar w2 = (x*x);
    const scalar w3 = 1 - 3.0/10.0*w2;
    const scalar w4 = k*y;
    const scalar w5 = hA1*cos(w4);
    const scalar w6 = (z*z);
    const scalar w7 = k*x;
    const scalar w8 = cos(gamma1);
    const scalar w9 = cos(gamma2);
    const scalar w10 = sin(gamma2);
    const scalar w11 = gin*w10*w9;
    const scalar w12 = w11*w8;
    const scalar w13 = (w10*w10);
    const scalar w14 = (r1*r1);
    const scalar w15 = -(r2*r2);
    const scalar w16 = w2 + w6 + (y*y);
    const scalar w17 = (-w14 + w16)/(-w14 - w15);
    const scalar w18 = theta0 + theta1*w17;
    const scalar w19 = sin(w18);
    const scalar w20 = cos(w18);
    const scalar w21 = gil*w13*(w19*w19) + gin*(w9*w9) + git*w13*(w20*w20);
    const scalar w22 = (1.0/5.0)*z;
    const scalar w23 = sin(gamma1);
    const scalar w24 = w20*w23;
    const scalar w25 = w19*w8;
    const scalar w26 = w24 + w25*w9;
    const scalar w27 = w10*w19;
    const scalar w28 = gil*w26*w27;
    const scalar w29 = w19*w23;
    const scalar w30 = w20*w8*w9 - w29;
    const scalar w31 = w10*w20;
    const scalar w32 = git*w30*w31;
    const scalar w33 = w20*w8 - w29*w9;
    const scalar w34 = gil*w33;
    const scalar w35 = -w24*w9 - w25;
    const scalar w36 = git*w35;
    const scalar w37 = w11*w23 + w27*w34 + w31*w36;
    const scalar w38 = (1.0/5.0)*k*sin(w7);
    const scalar w39 = -1.0/10.0*w12 + w21*w22 + (1.0/10.0)*w28 + (1.0/10.0)*w32 - w37*w38;
    const scalar w40 = (1.0/(r2));
    const scalar w41 = w40*z;
    const scalar w42 = gin*w13;
    const scalar w43 = w42*(w8*w8);
    const scalar w44 = gil*(w26*w26);
    const scalar w45 = git*(w30*w30);
    const scalar w46 = -w12 + w28 + w32;
    const scalar w47 = w23*w42*w8;
    const scalar w48 = w26*w34;
    const scalar w49 = w30*w36;
    const scalar w50 = -w47 + w48 + w49;
    const scalar w51 = w22*w46 - w38*w50 + (1.0/10.0)*w43 + (1.0/10.0)*w44 + (1.0/10.0)*w45;
    const scalar w52 = w40*y;
    const scalar w53 = gil*(w33*w33) + git*(w35*w35) + (w23*w23)*w42;
    const scalar w54 = w22*w37 - w38*w53 - 1.0/10.0*w47 + (1.0/10.0)*w48 + (1.0/10.0)*w49;
    const scalar w55 = w40*x;
    const scalar w56 = w43 + w44 + w45;
    const scalar w57 = r1 + r2;
    const scalar w58 = (1.0/2.0)*w57;
    const scalar w59 = w17*w40*w58/(w41*(w21*w41 + w37*w55 + w46*w52) + w52*(w41*w46 + w50*w55 + w52*w56) + w55*(w37*w41 + w50*w52 + w53*w55));
    const scalar w60 = (1.0/(r1));
    const scalar w61 = w60*z;
    const scalar w62 = w60*y;
    const scalar w63 = w60*x;
    const scalar w64 = w58*w60*(1 - w17)/(-w61*(-w21*w61 - w37*w63 - w46*w62) - w62*(-w46*w61 - w50*w63 - w56*w62) - w63*(-w37*w61 - w50*w62 - w53*w63));
    const scalar w65 = w17*(w15 + w16)/w57;
    const scalar w66 = St*(hA0 + hA2*w7*cos(w0));
    const scalar w67 = St*hA1*k*w3*sin(w4);
    const scalar w68 = St*(hA2*w1 - 3.0/5.0*w5*x);
    const scalar w69 = w21*w66 + w37*w68 - w46*w67;
    const scalar w70 = w46*w66 + w50*w68 - w56*w67;
    const scalar w71 = w37*w66 - w50*w67 + w53*w68;
    return St*(hA0*z + hA2*w1*x + w3*w5) + Tt*((1.0/10.0)*w6 + w65*(w59*(-w39*w41 - w51*w52 - w54*w55) + w64*(w39*w61 + w51*w62 + w54*w63)) + (1.0/10.0)*y + (1.0/5.0)*cos(w7)) + w65*(w59*(-w41*w69 - w52*w70 - w55*w71) + w64*(w61*w69 + w62*w70 + w63*w71));
}

inline vector gradPhiE(const point& pt, const scalar t, const label zone, const Params& p)
{
    const scalar r1 = p.r1;
    const scalar r2 = p.r2;
    const scalar gamma1 = p.gamma1;
    const scalar gamma2 = p.gamma2;
    const scalar k = p.k;
    const scalar gil = p.gil;
    const scalar git = p.git;
    const scalar gin = p.gin;
    const scalar gel = p.gel;
    const scalar get = p.get;
    const scalar gen = p.gen;
    const scalar theta0 = p.theta0;
    const scalar theta1 = p.theta1;
    const scalar sigma0In = p.sigma0In;
    const scalar sigma0Out = p.sigma0Out;
    const scalar sigmaEps = p.sigmaEps;
    const scalar fA0 = p.fA[0];
    const scalar fA1 = p.fA[1];
    const scalar fA2 = p.fA[2];
    const scalar fA3 = p.fA[3];
    const scalar cmAmplitude = p.cmAmplitude;
    const scalar cmOmega = p.cmOmega;
    const scalar x = pt.x();
    const scalar y = pt.y();
    const scalar z = pt.z();
    const scalar St = sqrt(1 + t);
    const scalar Tt = cmAmplitude*sin(cmOmega*t);
    checkZone(zone, true);

    if (zone == bloodPool)
    {
        const scalar w0 = k*x;
        const scalar w1 = cos(w0);
        const scalar w2 = sin(w0);
        const scalar w3 = (1.0/5.0)*k*w2;
        const scalar w4 = (1.0/(r1));
        const scalar w5 = w4*x;
        const scalar w6 = w4*y;
        const scalar w7 = (z*z);
        const scalar w8 = (1.0/5.0)*w4;
        const scalar w9 = (1.0/5.0)*k*w2*w4*x - 1.0/10.0*w6 - w7*w8;
        const scalar w10 = w0 + 1.0/5.0;
        const scalar w11 = sin(w10);
        const scalar w12 = k*z;
        const scalar w13 = (1.0/2.0)*sin(w12) + 1;
        const scalar w14 = k*y;
        const scalar w15 = cos(w14);
        const scalar w16 = sigmaEps*w13*w15;
        const scalar w17 = w11*w16 + 1;
        const scalar w18 = sigma0In*w17;
        const scalar w19 = cos(gamma1);
        const scalar w20 = cos(gamma2);
        const scalar w21 = sin(gamma2);
        const scalar w22 = w20*w21;
        const scalar w23 = w19*w22;
        const scalar w24 = gen*w23;
        const scalar w25 = sin(gamma1);
        const scalar w26 = cos(theta0);
        const scalar w27 = w25*w26;
        const scalar w28 = sin(theta0);
        const scalar w29 = w19*w28;
        const scalar w30 = w20*w29 + w27;
        const scalar w31 = w21*w28;
        const scalar w32 = w30*w31;
        const scalar w33 = gel*w32;
        const scalar w34 = w25*w28;
        const scalar w35 = w19*w20*w26 - w34;
        const scalar w36 = w21*w26;
        const scalar w37 = w35*w36;
        const scalar w38 = get*w37;
        const scalar w39 = (w20*w20);
        const scalar w40 = gen*w39;
        const scalar w41 = (w21*w21);
        const scalar w42 = (w28*w28)*w41;
        const scalar w43 = gel*w42;
        const scalar w44 = (w26*w26)*w41;
        const scalar w45 = get*w44;
        const scalar w46 = w40 + w43 + w45;
        const scalar w47 = (1.0/5.0)*z;
        const scalar w48 = gen*w25;
        const scalar w49 = w22*w48;
        const scalar w50 = w19*w26 - w20*w34;
        const scalar w51 = w31*w50;
        const scalar w52 = gel*w51;
        const scalar w53 = -w20*w27 - w29;
        const scalar w54 = w36*w53;
        const scalar w55 = get*w54;
        const scalar w56 = w49 + w52 + w55;
        const scalar w57 = w4*(-1.0/10.0*w24 - w3*w56 + (1.0/10.0)*w33 + (1.0/10.0)*w38 + w46*w47);
        const scalar w58 = (w19*w19)*w41;
        const scalar w59 = gen*w58;
        const scalar w60 = (w30*w30);
        const scalar w61 = gel*w60;
        const scalar w62 = (w35*w35);
        const scalar w63 = get*w62;
        const scalar w64 = -w24 + w33 + w38;
        const scalar w65 = w19*w41;
        const scalar w66 = w48*w65;
        const scalar w67 = w30*w50;
        const scalar w68 = gel*w67;
        const scalar w69 = w35*w53;
        const scalar w70 = get*w69;
        const scalar w71 = -w66 + w68 + w70;
        const scalar w72 = w4*(-w3*w71 + w47*w64 + (1.0/10.0)*w59 + (1.0/10.0)*w61 + (1.0/10.0)*w63);
        const scalar w73 = (w25*w25)*w41;
        const scalar w74 = (w50*w50);
        const scalar w75 = (w53*w53);
        const scalar w76 = gel*w74 + gen*w73 + get*w75;
        const scalar w77 = w4*(-w3*w76 + w47*w56 - 1.0/10.0*w66 + (1.0/10.0)*w68 + (1.0/10.0)*w70);
        const scalar w78 = gin*w23;
        const scalar w79 = gil*w32;
        const scalar w80 = git*w37;
        const scalar w81 = gin*w39;
        const scalar w82 = gil*w42;
        const scalar w83 = git*w44;
        const scalar w84 = z*(w81 + w82 + w83);
        const scalar w85 = gin*w25;
        const scalar w86 = w22*w85;
        const scalar w87 = gil*w51;
        const scalar w88 = git*w54;
        const scalar w89 = w86 + w87 + w88;
        const scalar w90 = w4*(-w3*w89 - 1.0/10.0*w78 + (1.0/10.0)*w79 + (1.0/10.0)*w80 + (1.0/5.0)*w84);
        const scalar w91 = gin*w58;
        const scalar w92 = gil*w60;
        const scalar w93 = git*w62;
        const scalar w94 = -w78 + w79 + w80;
        const scalar w95 = w94*z;
        const scalar w96 = w65*w85;
        const scalar w97 = gil*w67;
        const scalar w98 = git*w69;
        const scalar w99 = -w96 + w97 + w98;
        const scalar w100 = -w3*w99 + (1.0/10.0)*w91 + (1.0/10.0)*w92 + (1.0/10.0)*w93 + (1.0/5.0)*w95;
        const scalar w101 = w89*z;
        const scalar w102 = (1.0/5.0)*w101;
        const scalar w103 = gil*w74 + gin*w73 + git*w75;
        const scalar w104 = w4*(w102 - w103*w3 - 1.0/10.0*w96 + (1.0/10.0)*w97 + (1.0/10.0)*w98);
        const scalar w105 = -w100*w6 - w104*x - w90*z;
        const scalar w106 = w4*z;
        const scalar w107 = w4*(-w106*w46 - w5*w56 - w6*w64);
        const scalar w108 = w59 + w61 + w63;
        const scalar w109 = w4*(-w106*w64 - w108*w6 - w5*w71);
        const scalar w110 = w4*(-w106*w56 - w5*w76 - w6*w71);
        const scalar w111 = -w107*z - w109*y - w110*x;
        const scalar w112 = (r1*r1);
        const scalar w113 = (1.0/(((1.0/(git)) + (1.0/(gin)) + (1.0/(gil)))));
        const scalar w114 = w112 - w7 - (x*x) - (y*y);
        const scalar w115 = w94*y;
        const scalar w116 = w89*x;
        const scalar w117 = y*(w91 + w92 + w93);
        const scalar w118 = w99*x;
        const scalar w119 = w103*x;
        const scalar w120 = w99*y;
        const scalar w121 = w113*w114 + x*(w101 + w119 + w120) + y*(w117 + w118 + w95) + z*(w115 + w116 + w84);
        const scalar w122 = (1.0/(w121));
        const scalar w123 = w112*w122;
        const scalar w124 = w111*w123;
        const scalar w125 = -w105*w124 - w18*w9 - w57*z - w72*y - w77*x;
        const scalar w126 = (1.0/(sigma0In));
        const scalar w127 = (1.0/(w17));
        const scalar w128 = w126*w127;
        const scalar w129 = w125*w128;
        const scalar w130 = w16*cos(w10);
        const scalar w131 = w114*w4;
        const scalar w132 = (1.0/2.0)*w131;
        const scalar w133 = k*w125*w126/(w17*w17);
        const scalar w134 = w132*w133;
        const scalar w135 = (k*k);
        const scalar w136 = w1*w135;
        const scalar w137 = k*sigma0In*w9;
        const scalar w138 = (1.0/(w112));
        const scalar w139 = w138*x;
        const scalar w140 = w138*z;
        const scalar w141 = w138*y;
        const scalar w142 = w105*w123;
        const scalar w143 = w136*w8;
        const scalar w144 = w105*w111*w112/(w121*w121);
        const scalar w145 = w14 + 3.0/10.0;
        const scalar w146 = fA3*exp((2.0/5.0)*z);
        const scalar w147 = sigmaEps*w11;
        const scalar w148 = w13*sin(w14);
        const scalar w149 = w137*w147;
        const scalar w150 = w15*cos(w12);
        return vector(St*(fA0 + fA2*k*w1*((1.0/2.0)*y + 1)) + Tt*((1.0/2.0)*w114*w126*w127*w4*((1.0/5.0)*w1*w135*w4*w56*z + (1.0/5.0)*w1*w135*w4*w71*y + (1.0/5.0)*w1*w135*w4*w76*x - w124*(w102*w136*w4 - w104 + w119*w143 + w120*w143) - w130*w137 - w142*(-w110 + w139*w76 + w140*w56 + w141*w71) - w144*(-2*w101 + 2*w113*x - 2*w119 - 2*w120) - w18*((1.0/5.0)*w136*w5 + w3*w4) - w77) - w129*w5 - w130*w134 - w3), St*(fA1*z + (1.0/2.0)*fA2*w2 - k*w146*sin(w145)) + Tt*(w128*w132*(r1*w100*w111*w122 - w142*(w108*w141 - w109 + w139*w71 + w140*w64) - w144*(2*w113*y - 2*w117 - 2*w118 - 2*w95) + w148*w149 + (1.0/10.0)*w18*w4 - w72) - w129*w6 + w134*w147*w148 + 1.0/10.0), St*(fA1*y + (2.0/5.0)*w146*cos(w145)) + Tt*(-w106*w129 + (1.0/2.0)*w114*w126*w127*w4*((2.0/5.0)*sigma0In*w17*w4*z - w106*((1.0/5.0)*w40 + (1.0/5.0)*w43 + (1.0/5.0)*w45) - w124*(-w106*((1.0/5.0)*w81 + (1.0/5.0)*w82 + (1.0/5.0)*w83) - w5*((1.0/5.0)*w86 + (1.0/5.0)*w87 + (1.0/5.0)*w88) - w6*(-1.0/5.0*w78 + (1.0/5.0)*w79 + (1.0/5.0)*w80) - w90) - w142*(-w107 + w139*w56 + w140*w46 + w141*w64) - w144*(2*w113*z - 2*w115 - 2*w116 - 2*w84) - 1.0/2.0*w149*w150 - w5*((1.0/5.0)*w49 + (1.0/5.0)*w52 + (1.0/5.0)*w55) - w57 - w6*(-1.0/5.0*w24 + (1.0/5.0)*w33 + (1.0/5.0)*w38)) - 1.0/4.0*w131*w133*w147*w150 + (1.0/5.0)*z));
    }
    if (zone == bath)
    {
        const scalar w0 = k*x;
        const scalar w1 = cos(w0);
        const scalar w2 = sin(w0);
        const scalar w3 = (1.0/5.0)*w2;
        const scalar w4 = k*w3;
        const scalar w5 = (1.0/(r2));
        const scalar w6 = (1.0/(sigma0Out));
        const scalar w7 = w0 + 1.0/5.0;
        const scalar w8 = sin(w7);
        const scalar w9 = k*z;
        const scalar w10 = (1.0/2.0)*sin(w9) + 1;
        const scalar w11 = k*y;
        const scalar w12 = cos(w11);
        const scalar w13 = sigmaEps*w10*w12;
        const scalar w14 = w13*w8 + 1;
        const scalar w15 = (1.0/(w14));
        const scalar w16 = w5*y;
        const scalar w17 = (z*z);
        const scalar w18 = (1.0/5.0)*w5;
        const scalar w19 = -w0*w3*w5 + (1.0/10.0)*w16 + w17*w18;
        const scalar w20 = sigma0Out*w14;
        const scalar w21 = cos(gamma1);
        const scalar w22 = cos(gamma2);
        const scalar w23 = sin(gamma2);
        const scalar w24 = w22*w23;
        const scalar w25 = w21*w24;
        const scalar w26 = gen*w25;
        const scalar w27 = sin(gamma1);
        const scalar w28 = theta0 + theta1;
        const scalar w29 = cos(w28);
        const scalar w30 = w27*w29;
        const scalar w31 = sin(w28);
        const scalar w32 = w21*w31;
        const scalar w33 = w22*w32 + w30;
        const scalar w34 = w23*w31;
        const scalar w35 = gel*w34;
        const scalar w36 = w33*w35;
        const scalar w37 = w27*w31;
        const scalar w38 = w21*w22*w29 - w37;
        const scalar w39 = w23*w29;
        const scalar w40 = get*w39;
        const scalar w41 = w38*w40;
        const scalar w42 = (w22*w22);
        const scalar w43 = gen*w42;
        const scalar w44 = (w23*w23);
        const scalar w45 = (w31*w31)*w44;
        const scalar w46 = gel*w45;
        const scalar w47 = (w29*w29)*w44;
        const scalar w48 = get*w47;
        const scalar w49 = w43 + w46 + w48;
        const scalar w50 = (1.0/5.0)*z;
        const scalar w51 = gen*w27;
        const scalar w52 = w24*w51;
        const scalar w53 = w21*w29 - w22*w37;
        const scalar w54 = w35*w53;
        const scalar w55 = -w22*w30 - w32;
        const scalar w56 = w40*w55;
        const scalar w57 = w52 + w54 + w56;
        const scalar w58 = -1.0/10.0*w26 + (1.0/10.0)*w36 - w4*w57 + (1.0/10.0)*w41 + w49*w50;
        const scalar w59 = (w21*w21)*w44;
        const scalar w60 = gen*w59;
        const scalar w61 = (w33*w33);
        const scalar w62 = gel*w61;
        const scalar w63 = (w38*w38);
        const scalar w64 = get*w63;
        const scalar w65 = -w26 + w36 + w41;
        const scalar w66 = w21*w44;
        const scalar w67 = w51*w66;
        const scalar w68 = w33*w53;
        const scalar w69 = gel*w68;
        const scalar w70 = w38*w55;
        const scalar w71 = get*w70;
        const scalar w72 = -w67 + w69 + w71;
        const scalar w73 = -w4*w72 + w50*w65 + (1.0/10.0)*w60 + (1.0/10.0)*w62 + (1.0/10.0)*w64;
        const scalar w74 = (w27*w27)*w44;
        const scalar w75 = (w53*w53);
        const scalar w76 = (w55*w55);
        const scalar w77 = gel*w75 + gen*w74 + get*w76;
        const scalar w78 = w50*w57;
        const scalar w79 = -w4*w77 - 1.0/10.0*w67 + (1.0/10.0)*w69 + (1.0/10.0)*w71 + w78;
        const scalar w80 = gin*w25;
        const scalar w81 = gil*w34;
        const scalar w82 = w33*w81;
        const scalar w83 = git*w39;
        const scalar w84 = w38*w83;
        const scalar w85 = gin*w42;
        const scalar w86 = gil*w45;
        const scalar w87 = git*w47;
        const scalar w88 = z*(w85 + w86 + w87);
        const scalar w89 = gin*w27;
        const scalar w90 = w24*w89;
        const scalar w91 = w53*w81;
        const scalar w92 = w55*w83;
        const scalar w93 = w90 + w91 + w92;
        const scalar w94 = w5*(-w4*w93 - 1.0/10.0*w80 + (1.0/10.0)*w82 + (1.0/10.0)*w84 + (1.0/5.0)*w88);
        const scalar w95 = gin*w59;
        const scalar w96 = gil*w61;
        const scalar w97 = git*w63;
        const scalar w98 = -w80 + w82 + w84;
        const scalar w99 = w98*z;
        const scalar w100 = w66*w89;
        const scalar w101 = gil*w68;
        const scalar w102 = git*w70;
        const scalar w103 = -w100 + w101 + w102;
        const scalar w104 = -w103*w4 + (1.0/10.0)*w95 + (1.0/10.0)*w96 + (1.0/10.0)*w97 + (1.0/5.0)*w99;
        const scalar w105 = gil*w75 + gin*w74 + git*w76;
        const scalar w106 = w93*z;
        const scalar w107 = -1.0/10.0*w100 + (1.0/10.0)*w101 + (1.0/10.0)*w102 - w105*w4 + (1.0/5.0)*w106;
        const scalar w108 = w104*w16 + w107*w5*x + w94*z;
        const scalar w109 = (r2*r2);
        const scalar w110 = w98*y;
        const scalar w111 = w93*x;
        const scalar w112 = y*(w95 + w96 + w97);
        const scalar w113 = w103*x;
        const scalar w114 = w105*x;
        const scalar w115 = w103*y;
        const scalar w116 = x*(w106 + w114 + w115) + y*(w112 + w113 + w99) + z*(w110 + w111 + w88);
        const scalar w117 = (1.0/(w116));
        const scalar w118 = w5*z;
        const scalar w119 = w5*x;
        const scalar w120 = w5*(w118*w49 + w119*w57 + w16*w65);
        const scalar w121 = w60 + w62 + w64;
        const scalar w122 = w5*(w118*w65 + w119*w72 + w121*w16);
        const scalar w123 = w5*(w118*w57 + w119*w77 + w16*w72);
        const scalar w124 = w120*z + w122*y + w123*x;
        const scalar w125 = w117*w124;
        const scalar w126 = w109*w125;
        const scalar w127 = -w108*w126 - w19*w20 + w5*w58*z + w5*w73*y + w5*w79*x;
        const scalar w128 = w13*cos(w7);
        const scalar w129 = -w109 + w17 + (x*x) + (y*y);
        const scalar w130 = w129*w5;
        const scalar w131 = (1.0/2.0)*w130;
        const scalar w132 = k*w127*w6/(w14*w14);
        const scalar w133 = w131*w132;
        const scalar w134 = (k*k)*w1;
        const scalar w135 = w134*w18;
        const scalar w136 = w135*x;
        const scalar w137 = k*sigma0Out*w19;
        const scalar w138 = (1.0/(w109));
        const scalar w139 = w138*x;
        const scalar w140 = w138*z;
        const scalar w141 = w138*y;
        const scalar w142 = w108*w109;
        const scalar w143 = w117*w142;
        const scalar w144 = w124*w142/(w116*w116);
        const scalar w145 = w11 + 3.0/10.0;
        const scalar w146 = (2.0/5.0)*z;
        const scalar w147 = fA3*exp(w146);
        const scalar w148 = w15*w6;
        const scalar w149 = w127*w148;
        const scalar w150 = sigmaEps*w8;
        const scalar w151 = sin(w11);
        const scalar w152 = w20*w5;
        const scalar w153 = w131*w148;
        const scalar w154 = w12*cos(w9);
        return vector(St*(fA0 + fA2*k*w1*((1.0/2.0)*y + 1)) + Tt*(w127*w15*w5*w6*x - w128*w133 + (1.0/2.0)*w129*w15*w5*w6*(-w126*(-w106*w135 + w107*w5 - w114*w135 - w115*w135) - w128*w137 - w134*w5*w78 - w135*w72*y - w136*w77 - w143*(w123 + w139*w77 + w140*w57 + w141*w72) - w144*(-2*w106 - 2*w114 - 2*w115) - w20*(-w136 - w4*w5) + w5*w79) - w4), St*(fA1*z + (1.0/2.0)*fA2*w2 - k*w147*sin(w145)) + Tt*(w10*w133*w150*w151 + w149*w16 + w153*(k*sigma0Out*sigmaEps*w10*w151*w19*w8 - r2*w104*w125 - w143*(w121*w141 + w122 + w139*w72 + w140*w65) - w144*(-2*w112 - 2*w113 - 2*w99) - 1.0/10.0*w152 + w5*w73) + 1.0/10.0), St*(fA1*y + (2.0/5.0)*w147*cos(w145)) + Tt*(w118*w149 - 1.0/4.0*w130*w132*w150*w154 + w153*(-w126*(w118*((1.0/5.0)*w85 + (1.0/5.0)*w86 + (1.0/5.0)*w87) + w119*((1.0/5.0)*w90 + (1.0/5.0)*w91 + (1.0/5.0)*w92) + w16*(-1.0/5.0*w80 + (1.0/5.0)*w82 + (1.0/5.0)*w84) + w94) - 1.0/2.0*w137*w150*w154 - w143*(w120 + w139*w57 + w140*w49 + w141*w65) - w144*(-2*w110 - 2*w111 - 2*w88) - w146*w152 + w5*w58 + w5*x*((1.0/5.0)*w52 + (1.0/5.0)*w54 + (1.0/5.0)*w56) + w5*y*(-1.0/5.0*w26 + (1.0/5.0)*w36 + (1.0/5.0)*w41) + w5*z*((1.0/5.0)*w43 + (1.0/5.0)*w46 + (1.0/5.0)*w48)) + w50));
    }
    const scalar w0 = k*x;
    const scalar w1 = cos(w0);
    const scalar w2 = k*w1;
    const scalar w3 = fA2*((1.0/2.0)*y + 1);
    const scalar w4 = St*(fA0 + w2*w3);
    const scalar w5 = r1 + r2;
    const scalar w6 = (1.0/(w5));
    const scalar w7 = (r1*r1);
    const scalar w8 = -(r2*r2);
    const scalar w9 = (1.0/((-w7 - w8)));
    const scalar w10 = (x*x) + (y*y) + (z*z);
    const scalar w11 = w9*(w10 - w7);
    const scalar w12 = w11*w6;
    const scalar w13 = (1.0/(r2));
    const scalar w14 = cos(gamma2);
    const scalar w15 = (w14*w14);
    const scalar w16 = sin(gamma2);
    const scalar w17 = (w16*w16);
    const scalar w18 = theta0 + theta1*w11;
    const scalar w19 = sin(w18);
    const scalar w20 = w17*(w19*w19);
    const scalar w21 = cos(w18);
    const scalar w22 = w17*(w21*w21);
    const scalar w23 = gel*w20 + gen*w15 + get*w22;
    const scalar w24 = w13*w23;
    const scalar w25 = cos(gamma1);
    const scalar w26 = w14*w16;
    const scalar w27 = gen*w26;
    const scalar w28 = sin(gamma1);
    const scalar w29 = w21*w28;
    const scalar w30 = w19*w25;
    const scalar w31 = w14*w30;
    const scalar w32 = w29 + w31;
    const scalar w33 = gel*w32;
    const scalar w34 = w16*w19;
    const scalar w35 = w19*w28;
    const scalar w36 = w14*w21*w25 - w35;
    const scalar w37 = get*w36;
    const scalar w38 = w16*w21;
    const scalar w39 = -w25*w27 + w33*w34 + w37*w38;
    const scalar w40 = w13*w39;
    const scalar w41 = w14*w35;
    const scalar w42 = w21*w25 - w41;
    const scalar w43 = gel*w42;
    const scalar w44 = w14*w29;
    const scalar w45 = -w30 - w44;
    const scalar w46 = get*w45;
    const scalar w47 = w27*w28 + w34*w43 + w38*w46;
    const scalar w48 = w13*w47;
    const scalar w49 = w13*(w24*z + w40*y + w48*x);
    const scalar w50 = (w25*w25);
    const scalar w51 = gen*w17;
    const scalar w52 = (w32*w32);
    const scalar w53 = (w36*w36);
    const scalar w54 = gel*w52 + get*w53 + w50*w51;
    const scalar w55 = w13*w54;
    const scalar w56 = w25*w28;
    const scalar w57 = w32*w43 + w36*w46 - w51*w56;
    const scalar w58 = w13*w57;
    const scalar w59 = w13*(w40*z + w55*y + w58*x);
    const scalar w60 = (w28*w28);
    const scalar w61 = (w42*w42);
    const scalar w62 = (w45*w45);
    const scalar w63 = gel*w61 + get*w62 + w51*w60;
    const scalar w64 = w13*w63;
    const scalar w65 = w13*(w48*z + w58*y + w64*x);
    const scalar w66 = w49*z + w59*y + w65*x;
    const scalar w67 = (1.0/(w66));
    const scalar w68 = w0 + 1.0/5.0;
    const scalar w69 = sin(w68);
    const scalar w70 = k*z;
    const scalar w71 = (1.0/2.0)*sin(w70) + 1;
    const scalar w72 = k*y;
    const scalar w73 = cos(w72);
    const scalar w74 = sigmaEps*w71*w73;
    const scalar w75 = w69*w74 + 1;
    const scalar w76 = -sigma0Out*w75;
    const scalar w77 = -w23 - w76;
    const scalar w78 = w72 + 3.0/10.0;
    const scalar w79 = cos(w78);
    const scalar w80 = fA3*exp((2.0/5.0)*z);
    const scalar w81 = (2.0/5.0)*w80;
    const scalar w82 = St*(fA1*y + w79*w81);
    const scalar w83 = -w47;
    const scalar w84 = -w39;
    const scalar w85 = sin(w0);
    const scalar w86 = (1.0/2.0)*fA2;
    const scalar w87 = k*sin(w78);
    const scalar w88 = St*(fA1*z - w80*w87 + w85*w86);
    const scalar w89 = w4*w83 + w84*w88;
    const scalar w90 = w13*(w77*w82 + w89);
    const scalar w91 = -w54 - w76;
    const scalar w92 = -w57;
    const scalar w93 = w4*w92 + w82*w84;
    const scalar w94 = w13*(w88*w91 + w93);
    const scalar w95 = -w63 - w76;
    const scalar w96 = w82*w83 + w88*w92;
    const scalar w97 = w13*(w4*w95 + w96);
    const scalar w98 = w90*z + w94*y + w97*x;
    const scalar w99 = w67*w98;
    const scalar w100 = (1.0/2.0)*w11*w13*w5;
    const scalar w101 = (1.0/(r1));
    const scalar w102 = -sigma0In*w75;
    const scalar w103 = -w102 - w23;
    const scalar w104 = w101*(w103*w82 + w89);
    const scalar w105 = -w102 - w54;
    const scalar w106 = w101*(w105*w88 + w93);
    const scalar w107 = -w102 - w63;
    const scalar w108 = w101*(w107*w4 + w96);
    const scalar w109 = -w104*z - w106*y - w108*x;
    const scalar w110 = w101*w23;
    const scalar w111 = w101*w39;
    const scalar w112 = w101*w47;
    const scalar w113 = w101*(-w110*z - w111*y - w112*x);
    const scalar w114 = w101*w54;
    const scalar w115 = w101*w57;
    const scalar w116 = w101*(-w111*z - w114*y - w115*x);
    const scalar w117 = w101*w63;
    const scalar w118 = w101*(-w112*z - w115*y - w117*x);
    const scalar w119 = -w113*z - w116*y - w118*x;
    const scalar w120 = (1.0/(w119));
    const scalar w121 = (1.0/2.0)*w101*(1 - w11);
    const scalar w122 = w121*w5;
    const scalar w123 = w120*w122;
    const scalar w124 = w100*w99 + w109*w123;
    const scalar w125 = 2*x;
    const scalar w126 = w124*w125;
    const scalar w127 = w10 + w8;
    const scalar w128 = w127*w6*w9;
    const scalar w129 = w13*x;
    const scalar w130 = w5*w9;
    const scalar w131 = w129*w130;
    const scalar w132 = w101*x;
    const scalar w133 = w109*w120*w130;
    const scalar w134 = cos(w68);
    const scalar w135 = w134*w74;
    const scalar w136 = k*sigma0Out;
    const scalar w137 = theta1*w9;
    const scalar w138 = 4*x;
    const scalar w139 = w137*w138;
    const scalar w140 = w17*w19*w21;
    const scalar w141 = w139*w140;
    const scalar w142 = gel*w141;
    const scalar w143 = get*w141;
    const scalar w144 = -w142 + w143;
    const scalar w145 = w125*w137;
    const scalar w146 = w145*w38;
    const scalar w147 = w137*w30;
    const scalar w148 = -w125*w147 - w145*w44;
    const scalar w149 = gel*w34;
    const scalar w150 = w137*w21*w25;
    const scalar w151 = -w125*w150 + w145*w41;
    const scalar w152 = get*w38;
    const scalar w153 = -2*get*theta1*w16*w19*w45*w9*x + w146*w43 + w148*w149 + w151*w152;
    const scalar w154 = -w153;
    const scalar w155 = 2*theta1*w14*w21*w25*w9*x - w145*w35;
    const scalar w156 = -w145*w29 - w145*w31;
    const scalar w157 = -2*get*theta1*w16*w19*w36*w9*x + w146*w33 + w149*w155 + w152*w156;
    const scalar w158 = -w157;
    const scalar w159 = w2*w86;
    const scalar w160 = St*w159;
    const scalar w161 = (k*k);
    const scalar w162 = w161*w3*w85;
    const scalar w163 = St*w162;
    const scalar w164 = w154*w4 + w158*w88 + w160*w84 - w163*w83;
    const scalar w165 = w13*z;
    const scalar w166 = -k*sigma0Out*sigmaEps*w134*w71*w73;
    const scalar w167 = 4*theta1*w14*w21*w25*w9*x - w139*w35;
    const scalar w168 = -w139*w29 - w139*w31;
    const scalar w169 = w167*w33 + w168*w37;
    const scalar w170 = w148*w33 + w151*w37 + w155*w43 + w156*w46;
    const scalar w171 = -w170;
    const scalar w172 = St*w92;
    const scalar w173 = w158*w82 - w162*w172 + w171*w4;
    const scalar w174 = w13*y;
    const scalar w175 = -w138*w147 - w139*w44;
    const scalar w176 = -w138*w150 + w139*w41;
    const scalar w177 = w175*w43 + w176*w46;
    const scalar w178 = w159*w172;
    const scalar w179 = w154*w82 + w171*w88 + w178;
    const scalar w180 = w100*w67;
    const scalar w181 = k*sigma0In;
    const scalar w182 = w101*z;
    const scalar w183 = -k*sigma0In*sigmaEps*w134*w71*w73;
    const scalar w184 = w101*y;
    const scalar w185 = St*w107;
    const scalar w186 = w142 - w143;
    const scalar w187 = w100*w98/(w66*w66);
    const scalar w188 = w109*w122/(w119*w119);
    const scalar w189 = w12*w127;
    const scalar w190 = (1.0/5.0)*k*w85;
    const scalar w191 = gin*w15;
    const scalar w192 = gil*w20;
    const scalar w193 = git*w22;
    const scalar w194 = w191 + w192 + w193;
    const scalar w195 = w13*w194;
    const scalar w196 = gin*w26;
    const scalar w197 = w196*w25;
    const scalar w198 = gil*w32;
    const scalar w199 = w198*w34;
    const scalar w200 = git*w36;
    const scalar w201 = w200*w38;
    const scalar w202 = -w197 + w199 + w201;
    const scalar w203 = w13*w202;
    const scalar w204 = w196*w28;
    const scalar w205 = gil*w42;
    const scalar w206 = w205*w34;
    const scalar w207 = git*w45;
    const scalar w208 = w207*w38;
    const scalar w209 = w204 + w206 + w208;
    const scalar w210 = w13*w209;
    const scalar w211 = w13*(w195*z + w203*y + w210*x);
    const scalar w212 = gin*w17;
    const scalar w213 = w212*w50;
    const scalar w214 = gil*w52;
    const scalar w215 = git*w53;
    const scalar w216 = w213 + w214 + w215;
    const scalar w217 = w13*w216;
    const scalar w218 = w212*w56;
    const scalar w219 = w205*w32;
    const scalar w220 = w207*w36;
    const scalar w221 = -w218 + w219 + w220;
    const scalar w222 = w13*w221;
    const scalar w223 = w13*(w203*z + w217*y + w222*x);
    const scalar w224 = gil*w61 + git*w62 + w212*w60;
    const scalar w225 = w13*w224;
    const scalar w226 = w13*(w210*z + w222*y + w225*x);
    const scalar w227 = w211*z + w223*y + w226*x;
    const scalar w228 = (1.0/(w227));
    const scalar w229 = (1.0/5.0)*z;
    const scalar w230 = -w190*w209 + w194*w229 - 1.0/10.0*w197 + (1.0/10.0)*w199 + (1.0/10.0)*w201;
    const scalar w231 = w13*w230;
    const scalar w232 = -w190*w221 + w202*w229 + (1.0/10.0)*w213 + (1.0/10.0)*w214 + (1.0/10.0)*w215;
    const scalar w233 = w13*w232;
    const scalar w234 = -w190*w224 + w209*w229 - 1.0/10.0*w218 + (1.0/10.0)*w219 + (1.0/10.0)*w220;
    const scalar w235 = w13*w234;
    const scalar w236 = -w231*z - w233*y - w235*x;
    const scalar w237 = w228*w236;
    const scalar w238 = w101*w230;
    const scalar w239 = w101*w232;
    const scalar w240 = w101*w234;
    const scalar w241 = w238*z + w239*y + w240*x;
    const scalar w242 = w101*w194;
    const scalar w243 = w101*w202;
    const scalar w244 = w101*w209;
    const scalar w245 = w101*(-w242*z - w243*y - w244*x);
    const scalar w246 = w101*w216;
    const scalar w247 = w101*w221;
    const scalar w248 = w101*(-w243*z - w246*y - w247*x);
    const scalar w249 = w101*w224;
    const scalar w250 = w101*(-w244*z - w247*y - w249*x);
    const scalar w251 = -w245*z - w248*y - w250*x;
    const scalar w252 = w5/w251;
    const scalar w253 = w121*w252;
    const scalar w254 = w100*w237 + w241*w253;
    const scalar w255 = w125*w254;
    const scalar w256 = w241*w252*w9;
    const scalar w257 = (1.0/5.0)*w137;
    const scalar w258 = w257*x;
    const scalar w259 = w198*w38;
    const scalar w260 = w200*w34;
    const scalar w261 = gil*w34;
    const scalar w262 = w155*w261;
    const scalar w263 = git*w38;
    const scalar w264 = w156*w263;
    const scalar w265 = gil*w141 - git*w141;
    const scalar w266 = (1.0/5.0)*w1*w161;
    const scalar w267 = w145*w34;
    const scalar w268 = w146*w205 + w148*w261 + w151*w263 - w207*w267;
    const scalar w269 = -w190*w268 - w209*w266 + w229*w265 + w258*w259 - w258*w260 + (1.0/10.0)*w262 + (1.0/10.0)*w264;
    const scalar w270 = w148*w198;
    const scalar w271 = w155*w205;
    const scalar w272 = w151*w200;
    const scalar w273 = w156*w207;
    const scalar w274 = w175*w205 + w176*w207;
    const scalar w275 = -w190*w274 - w224*w266 + w229*w268 + (1.0/10.0)*w270 + (1.0/10.0)*w271 + (1.0/10.0)*w272 + (1.0/10.0)*w273;
    const scalar w276 = w167*w198;
    const scalar w277 = w168*w200;
    const scalar w278 = w146*w198 - w200*w267 + w262 + w264;
    const scalar w279 = w270 + w271 + w272 + w273;
    const scalar w280 = -w190*w279 - w221*w266 + w229*w278 + (1.0/10.0)*w276 + (1.0/10.0)*w277;
    const scalar w281 = w100*w228;
    const scalar w282 = w276 + w277;
    const scalar w283 = w100*w236/(w227*w227);
    const scalar w284 = w122*w241/(w251*w251);
    const scalar w285 = 2*y;
    const scalar w286 = w124*w285;
    const scalar w287 = w130*w174;
    const scalar w288 = St*(fA1 - w81*w87);
    const scalar w289 = w71*sin(w72);
    const scalar w290 = sigmaEps*w69;
    const scalar w291 = w136*w290;
    const scalar w292 = w289*w291;
    const scalar w293 = 4*y;
    const scalar w294 = w137*w293;
    const scalar w295 = w140*w294;
    const scalar w296 = gel*w295 - get*w295;
    const scalar w297 = w137*w285;
    const scalar w298 = w297*w38;
    const scalar w299 = -w147*w285 - w297*w44;
    const scalar w300 = -w150*w285 + w297*w41;
    const scalar w301 = -2*get*theta1*w16*w19*w45*w9*y + w149*w299 + w152*w300 + w298*w43;
    const scalar w302 = -w301;
    const scalar w303 = 2*theta1*w14*w21*w25*w9*y - w297*w35;
    const scalar w304 = -w29*w297 - w297*w31;
    const scalar w305 = -2*get*theta1*w16*w19*w36*w9*y + w149*w303 + w152*w304 + w298*w33;
    const scalar w306 = -w305;
    const scalar w307 = w79*w80;
    const scalar w308 = w161*w307;
    const scalar w309 = St*w308;
    const scalar w310 = w160*w83 + w302*w4 + w306*w88 - w309*w84;
    const scalar w311 = 4*theta1*w14*w21*w25*w9*y - w294*w35;
    const scalar w312 = -w29*w294 - w294*w31;
    const scalar w313 = w311*w33 + w312*w37;
    const scalar w314 = w288*w84;
    const scalar w315 = w299*w33 + w300*w37 + w303*w43 + w304*w46;
    const scalar w316 = -w315;
    const scalar w317 = w178 + w306*w82 + w314 + w316*w4;
    const scalar w318 = -w147*w293 - w294*w44;
    const scalar w319 = -w150*w293 + w294*w41;
    const scalar w320 = w318*w43 + w319*w46;
    const scalar w321 = -w172*w308 + w288*w83 + w302*w82 + w316*w88;
    const scalar w322 = w181*w290;
    const scalar w323 = w289*w322;
    const scalar w324 = w254*w285;
    const scalar w325 = w257*y;
    const scalar w326 = w261*w303;
    const scalar w327 = w263*w304;
    const scalar w328 = gil*w295 - git*w295;
    const scalar w329 = w205*w298 - w207*w297*w34 + w261*w299 + w263*w300;
    const scalar w330 = -w190*w329 + w229*w328 + w259*w325 - w260*w325 + (1.0/10.0)*w326 + (1.0/10.0)*w327;
    const scalar w331 = w198*w311;
    const scalar w332 = w200*w312;
    const scalar w333 = w259*w297 - w260*w297 + w326 + w327;
    const scalar w334 = w198*w299;
    const scalar w335 = w205*w303;
    const scalar w336 = w200*w300;
    const scalar w337 = w207*w304;
    const scalar w338 = w334 + w335 + w336 + w337;
    const scalar w339 = -w190*w338 + w229*w333 + (1.0/10.0)*w331 + (1.0/10.0)*w332;
    const scalar w340 = w205*w318 + w207*w319;
    const scalar w341 = -w190*w340 + w229*w329 + (1.0/10.0)*w334 + (1.0/10.0)*w335 + (1.0/10.0)*w336 + (1.0/10.0)*w337;
    const scalar w342 = w331 + w332;
    const scalar w343 = 2*z;
    const scalar w344 = w124*w343;
    const scalar w345 = w130*w165;
    const scalar w346 = cos(w70);
    const scalar w347 = (1.0/2.0)*w346*w73;
    const scalar w348 = 4*z;
    const scalar w349 = w137*w348;
    const scalar w350 = w140*w349;
    const scalar w351 = gel*w350;
    const scalar w352 = get*w350;
    const scalar w353 = -w351 + w352;
    const scalar w354 = (4.0/25.0)*St*w307;
    const scalar w355 = w137*w343;
    const scalar w356 = w355*w38;
    const scalar w357 = -w147*w343 - w355*w44;
    const scalar w358 = -w150*w343 + w355*w41;
    const scalar w359 = -2*get*theta1*w16*w19*w45*w9*z + w149*w357 + w152*w358 + w356*w43;
    const scalar w360 = -w359;
    const scalar w361 = 2*theta1*w14*w21*w25*w9*z - w35*w355;
    const scalar w362 = -w29*w355 - w31*w355;
    const scalar w363 = -2*get*theta1*w16*w19*w36*w9*z + w149*w361 + w152*w362 + w33*w356;
    const scalar w364 = -w363;
    const scalar w365 = w314 + w360*w4 + w364*w88;
    const scalar w366 = -1.0/2.0*k*sigma0Out*sigmaEps*w346*w69*w73;
    const scalar w367 = 4*theta1*w14*w21*w25*w9*z - w349*w35;
    const scalar w368 = -w29*w349 - w31*w349;
    const scalar w369 = w33*w367 + w368*w37;
    const scalar w370 = w33*w357 + w358*w37 + w361*w43 + w362*w46;
    const scalar w371 = -w370;
    const scalar w372 = w354*w84 + w364*w82 + w371*w4;
    const scalar w373 = -w147*w348 - w349*w44;
    const scalar w374 = -w150*w348 + w349*w41;
    const scalar w375 = w373*w43 + w374*w46;
    const scalar w376 = w288*w92 + w354*w83 + w360*w82 + w371*w88;
    const scalar w377 = -1.0/2.0*k*sigma0In*sigmaEps*w346*w69*w73;
    const scalar w378 = w351 - w352;
    const scalar w379 = w254*w343;
    const scalar w380 = w137*w229;
    const scalar w381 = w261*w361;
    const scalar w382 = w263*w362;
    const scalar w383 = gil*w350 - git*w350;
    const scalar w384 = w205*w356 - w207*w34*w355 + w261*w357 + w263*w358;
    const scalar w385 = -w190*w384 + (1.0/5.0)*w191 + (1.0/5.0)*w192 + (1.0/5.0)*w193 + w229*w383 + w259*w380 - w260*w380 + (1.0/10.0)*w381 + (1.0/10.0)*w382;
    const scalar w386 = w198*w367;
    const scalar w387 = w200*w368;
    const scalar w388 = w259*w355 - w260*w355 + w381 + w382;
    const scalar w389 = w198*w357;
    const scalar w390 = w205*w361;
    const scalar w391 = w200*w358;
    const scalar w392 = w207*w362;
    const scalar w393 = w389 + w390 + w391 + w392;
    const scalar w394 = -w190*w393 - 1.0/5.0*w197 + (1.0/5.0)*w199 + (1.0/5.0)*w201 + w229*w388 + (1.0/10.0)*w386 + (1.0/10.0)*w387;
    const scalar w395 = w205*w373 + w207*w374;
    const scalar w396 = -w190*w395 + (1.0/5.0)*w204 + (1.0/5.0)*w206 + (1.0/5.0)*w208 + w229*w384 + (1.0/10.0)*w389 + (1.0/10.0)*w390 + (1.0/10.0)*w391 + (1.0/10.0)*w392;
    const scalar w397 = w386 + w387;
    return vector(Tt*(w12*w255 + w128*w255 + w189*(w131*w237 - w132*w256 + w253*(w132*w275 + w182*w269 + w184*w280 + w240) + w281*(-w129*w275 - w165*w269 - w174*w280 - w235) + w283*(-w129*(w129*w274 + w165*w268 + w174*w279 + w225) - w165*(w129*w268 + w165*w265 + w174*w278 + w210) - w174*(w129*w279 + w165*w278 + w174*w282 + w222) - w226) + w284*(w132*(-w132*w274 - w182*w268 - w184*w279 - w249) + w182*(-w132*w268 - w182*w265 - w184*w278 - w244) + w184*(-w132*w279 - w182*w278 - w184*w282 - w247) + w250)) - w190) + w12*w126 + w126*w128 + w189*(w123*(-w108 - w132*(-w162*w185 + w179 + w4*(-w177 - w183)) - w182*(w164 + w82*(w135*w181 + w144)) - w184*(w105*w160 + w173 + w88*(-w169 - w183))) + w131*w99 - w132*w133 + w180*(w129*(-w163*w95 + w179 + w4*(-w166 - w177)) + w165*(w164 + w82*(w135*w136 + w144)) + w174*(w160*w91 + w173 + w88*(-w166 - w169)) + w97) + w187*(-w129*(w129*w177 + w153*w165 + w170*w174 + w64) - w165*(w129*w153 + w157*w174 + w165*w186 + w48) - w174*(w129*w170 + w157*w165 + w169*w174 + w58) - w65) + w188*(w118 + w132*(-w117 - w132*w177 - w153*w182 - w170*w184) + w182*(-w112 - w132*w153 - w157*w184 - w182*w186) + w184*(-w115 - w132*w170 - w157*w182 - w169*w184))) + w4, Tt*(w12*w324 + w128*w324 + w189*(-w184*w256 + w237*w287 + w253*(w132*w341 + w182*w330 + w184*w339 + w239) + w281*(-w129*w341 - w165*w330 - w174*w339 - w233) + w283*(-w129*(w129*w340 + w165*w329 + w174*w338 + w222) - w165*(w129*w329 + w165*w328 + w174*w333 + w203) - w174*(w129*w338 + w165*w333 + w174*w342 + w217) - w223) + w284*(w132*(-w132*w340 - w182*w329 - w184*w338 - w247) + w182*(-w132*w329 - w182*w328 - w184*w333 - w243) + w184*(-w132*w338 - w182*w333 - w184*w342 - w246) + w248)) + 1.0/10.0) + w12*w286 + w128*w286 + w189*(w123*(-w106 - w132*(w159*w185 + w321 + w4*(-w320 - w323)) - w182*(w103*w288 + w310 + w82*(-w296 - w323)) - w184*(-w105*w309 + w317 + w88*(-w313 - w323))) - w133*w184 + w180*(w129*(w160*w95 + w321 + w4*(-w292 - w320)) + w165*(w288*w77 + w310 + w82*(-w292 - w296)) + w174*(-w309*w91 + w317 + w88*(-w292 - w313)) + w94) + w187*(-w129*(w129*w320 + w165*w301 + w174*w315 + w58) - w165*(w129*w301 + w165*w296 + w174*w305 + w40) - w174*(w129*w315 + w165*w305 + w174*w313 + w55) - w59) + w188*(w116 + w132*(-w115 - w132*w320 - w182*w301 - w184*w315) + w182*(-w111 - w132*w301 - w182*w296 - w184*w305) + w184*(-w114 - w132*w315 - w182*w305 - w184*w313)) + w287*w99) + w88, Tt*(w12*w379 + w128*w379 + w189*(-w182*w256 + w237*w345 + w253*(w132*w396 + w182*w385 + w184*w394 + w238) + w281*(-w129*w396 - w165*w385 - w174*w394 - w231) + w283*(-w129*(w129*w395 + w165*w384 + w174*w393 + w210) - w165*(w129*w384 + w165*w383 + w174*w388 + w195) - w174*(w129*w393 + w165*w388 + w174*w397 + w203) - w211) + w284*(w132*(-w132*w395 - w182*w384 - w184*w393 - w244) + w182*(-w132*w384 - w182*w383 - w184*w388 - w242) + w184*(-w132*w393 - w182*w388 - w184*w397 - w243) + w245)) + w229) + w12*w344 + w128*w344 + w189*(w123*(-w104 - w132*(w376 + w4*(-w375 - w377)) - w182*(w103*w354 + w365 + w82*(w322*w347 + w353)) - w184*(w105*w288 + w372 + w88*(-w369 - w377))) - w133*w182 + w180*(w129*(w376 + w4*(-w366 - w375)) + w165*(w354*w77 + w365 + w82*(w291*w347 + w353)) + w174*(w288*w91 + w372 + w88*(-w366 - w369)) + w90) + w187*(-w129*(w129*w375 + w165*w359 + w174*w370 + w48) - w165*(w129*w359 + w165*w378 + w174*w363 + w24) - w174*(w129*w370 + w165*w363 + w174*w369 + w40) - w49) + w188*(w113 + w132*(-w112 - w132*w375 - w182*w359 - w184*w370) + w182*(-w110 - w132*w359 - w182*w378 - w184*w363) + w184*(-w111 - w132*w370 - w182*w363 - w184*w369)) + w345*w99) + w82);
}

inline vector gradVm(const point& pt, const scalar t, const Params& p)
{
    const scalar r1 = p.r1;
    const scalar r2 = p.r2;
    const scalar gamma1 = p.gamma1;
    const scalar gamma2 = p.gamma2;
    const scalar k = p.k;
    const scalar gil = p.gil;
    const scalar git = p.git;
    const scalar gin = p.gin;
    const scalar gel = p.gel;
    const scalar get = p.get;
    const scalar gen = p.gen;
    const scalar theta0 = p.theta0;
    const scalar theta1 = p.theta1;
    const scalar sigma0In = p.sigma0In;
    const scalar sigma0Out = p.sigma0Out;
    const scalar sigmaEps = p.sigmaEps;
    const scalar fA0 = p.fA[0];
    const scalar fA1 = p.fA[1];
    const scalar fA2 = p.fA[2];
    const scalar fA3 = p.fA[3];
    const scalar hA0 = p.hA[0];
    const scalar hA1 = p.hA[1];
    const scalar hA2 = p.hA[2];
    const scalar x = pt.x();
    const scalar y = pt.y();
    const scalar z = pt.z();
    const scalar St = sqrt(1 + t);

    const scalar w0 = k*x;
    const scalar w1 = k*cos(w0);
    const scalar w2 = fA2*((1.0/2.0)*y + 1);
    const scalar w3 = St*(fA0 + w1*w2);
    const scalar w4 = k*z;
    const scalar w5 = sin(w4);
    const scalar w6 = k*y;
    const scalar w7 = cos(w6);
    const scalar w8 = hA1*w7;
    const scalar w9 = (3.0/5.0)*w8;
    const scalar w10 = hA2*w5 - w9*x;
    const scalar w11 = r1 + r2;
    const scalar w12 = (1.0/(w11));
    const scalar w13 = (r1*r1);
    const scalar w14 = -(r2*r2);
    const scalar w15 = (1.0/((-w13 - w14)));
    const scalar w16 = (x*x);
    const scalar w17 = w16 + (y*y) + (z*z);
    const scalar w18 = -w13 + w17;
    const scalar w19 = (1.0/(r2));
    const scalar w20 = cos(gamma2);
    const scalar w21 = (w20*w20);
    const scalar w22 = sin(gamma2);
    const scalar w23 = (w22*w22);
    const scalar w24 = w15*w18;
    const scalar w25 = theta0 + theta1*w24;
    const scalar w26 = sin(w25);
    const scalar w27 = w23*(w26*w26);
    const scalar w28 = cos(w25);
    const scalar w29 = w23*(w28*w28);
    const scalar w30 = gil*w27 + gin*w21 + git*w29;
    const scalar w31 = w19*w30;
    const scalar w32 = cos(gamma1);
    const scalar w33 = w20*w22;
    const scalar w34 = gin*w33;
    const scalar w35 = sin(gamma1);
    const scalar w36 = w28*w35;
    const scalar w37 = w26*w32;
    const scalar w38 = w20*w37;
    const scalar w39 = w36 + w38;
    const scalar w40 = gil*w39;
    const scalar w41 = w22*w26;
    const scalar w42 = w26*w35;
    const scalar w43 = w20*w28*w32 - w42;
    const scalar w44 = git*w43;
    const scalar w45 = w22*w28;
    const scalar w46 = -w32*w34 + w40*w41 + w44*w45;
    const scalar w47 = w19*w46;
    const scalar w48 = w20*w42;
    const scalar w49 = w28*w32 - w48;
    const scalar w50 = gil*w49;
    const scalar w51 = w20*w36;
    const scalar w52 = -w37 - w51;
    const scalar w53 = git*w52;
    const scalar w54 = w34*w35 + w41*w50 + w45*w53;
    const scalar w55 = w19*w54;
    const scalar w56 = w19*(w31*z + w47*y + w55*x);
    const scalar w57 = (w32*w32);
    const scalar w58 = gin*w23;
    const scalar w59 = (w39*w39);
    const scalar w60 = (w43*w43);
    const scalar w61 = gil*w59 + git*w60 + w57*w58;
    const scalar w62 = w19*w61;
    const scalar w63 = w32*w35;
    const scalar w64 = w39*w50 + w43*w53 - w58*w63;
    const scalar w65 = w19*w64;
    const scalar w66 = w19*(w47*z + w62*y + w65*x);
    const scalar w67 = (w35*w35);
    const scalar w68 = (w49*w49);
    const scalar w69 = (w52*w52);
    const scalar w70 = gil*w68 + git*w69 + w58*w67;
    const scalar w71 = w19*w70;
    const scalar w72 = w19*(w55*z + w65*y + w71*x);
    const scalar w73 = w56*z + w66*y + w72*x;
    const scalar w74 = (1.0/(w73));
    const scalar w75 = cos(w4);
    const scalar w76 = hA2*w75;
    const scalar w77 = St*(hA0 + w0*w76);
    const scalar w78 = 1 - 3.0/10.0*w16;
    const scalar w79 = St*k;
    const scalar w80 = sin(w6);
    const scalar w81 = hA1*w80;
    const scalar w82 = w78*w79*w81;
    const scalar w83 = St*w10;
    const scalar w84 = w30*w77 - w46*w82 + w54*w83;
    const scalar w85 = w19*w84;
    const scalar w86 = w46*w77 - w61*w82 + w64*w83;
    const scalar w87 = w19*w86;
    const scalar w88 = w54*w77 - w64*w82 + w70*w83;
    const scalar w89 = w19*w88;
    const scalar w90 = -w85*z - w87*y - w89*x;
    const scalar w91 = w74*w90;
    const scalar w92 = (1.0/2.0)*w11*w19*w24;
    const scalar w93 = (1.0/(r1));
    const scalar w94 = w84*w93;
    const scalar w95 = w86*w93;
    const scalar w96 = w88*w93;
    const scalar w97 = w94*z + w95*y + w96*x;
    const scalar w98 = w30*w93;
    const scalar w99 = w46*w93;
    const scalar w100 = w54*w93;
    const scalar w101 = w93*(-w100*x - w98*z - w99*y);
    const scalar w102 = w61*w93;
    const scalar w103 = w64*w93;
    const scalar w104 = w93*(-w102*y - w103*x - w99*z);
    const scalar w105 = w70*w93;
    const scalar w106 = w93*(-w100*z - w103*y - w105*x);
    const scalar w107 = -w101*z - w104*y - w106*x;
    const scalar w108 = w11/w107;
    const scalar w109 = (1.0/2.0)*w93*(1 - w24);
    const scalar w110 = w108*w109;
    const scalar w111 = w110*w97 + w91*w92;
    const scalar w112 = w14 + w17;
    const scalar w113 = w12*w24;
    const scalar w114 = gel*w27 + gen*w21 + get*w29;
    const scalar w115 = w114*w19;
    const scalar w116 = gen*w33;
    const scalar w117 = gel*w39;
    const scalar w118 = get*w43;
    const scalar w119 = -w116*w32 + w117*w41 + w118*w45;
    const scalar w120 = w119*w19;
    const scalar w121 = gel*w49;
    const scalar w122 = get*w52;
    const scalar w123 = w116*w35 + w121*w41 + w122*w45;
    const scalar w124 = w123*w19;
    const scalar w125 = w19*(w115*z + w120*y + w124*x);
    const scalar w126 = gen*w23;
    const scalar w127 = gel*w59 + get*w60 + w126*w57;
    const scalar w128 = w127*w19;
    const scalar w129 = w121*w39 + w122*w43 - w126*w63;
    const scalar w130 = w129*w19;
    const scalar w131 = w19*(w120*z + w128*y + w130*x);
    const scalar w132 = gel*w68 + get*w69 + w126*w67;
    const scalar w133 = w132*w19;
    const scalar w134 = w19*(w124*z + w130*y + w133*x);
    const scalar w135 = w125*z + w131*y + w134*x;
    const scalar w136 = (1.0/(w135));
    const scalar w137 = w0 + 1.0/5.0;
    const scalar w138 = sin(w137);
    const scalar w139 = (1.0/2.0)*w5 + 1;
    const scalar w140 = sigmaEps*w139*w7;
    const scalar w141 = w138*w140 + 1;
    const scalar w142 = -sigma0Out*w141;
    const scalar w143 = -w114 - w142;
    const scalar w144 = w6 + 3.0/10.0;
    const scalar w145 = cos(w144);
    const scalar w146 = fA3*exp((2.0/5.0)*z);
    const scalar w147 = (2.0/5.0)*w146;
    const scalar w148 = St*(fA1*y + w145*w147);
    const scalar w149 = -w123;
    const scalar w150 = -w119;
    const scalar w151 = sin(w0);
    const scalar w152 = (1.0/2.0)*fA2;
    const scalar w153 = k*sin(w144);
    const scalar w154 = St*(fA1*z - w146*w153 + w151*w152);
    const scalar w155 = w149*w3 + w150*w154;
    const scalar w156 = w19*(w143*w148 + w155);
    const scalar w157 = -w127 - w142;
    const scalar w158 = -w129;
    const scalar w159 = w148*w150 + w158*w3;
    const scalar w160 = w19*(w154*w157 + w159);
    const scalar w161 = -w132 - w142;
    const scalar w162 = w148*w149 + w154*w158;
    const scalar w163 = w19*(w161*w3 + w162);
    const scalar w164 = w156*z + w160*y + w163*x;
    const scalar w165 = w136*w164;
    const scalar w166 = -sigma0In*w141;
    const scalar w167 = -w114 - w166;
    const scalar w168 = w93*(w148*w167 + w155);
    const scalar w169 = -w127 - w166;
    const scalar w170 = w93*(w154*w169 + w159);
    const scalar w171 = -w132 - w166;
    const scalar w172 = w93*(w162 + w171*w3);
    const scalar w173 = -w168*z - w170*y - w172*x;
    const scalar w174 = w114*w93;
    const scalar w175 = w119*w93;
    const scalar w176 = w123*w93;
    const scalar w177 = w93*(-w174*z - w175*y - w176*x);
    const scalar w178 = w127*w93;
    const scalar w179 = w129*w93;
    const scalar w180 = w93*(-w175*z - w178*y - w179*x);
    const scalar w181 = w132*w93;
    const scalar w182 = w93*(-w176*z - w179*y - w181*x);
    const scalar w183 = -w177*z - w180*y - w182*x;
    const scalar w184 = (1.0/(w183));
    const scalar w185 = w109*w11;
    const scalar w186 = w184*w185;
    const scalar w187 = w165*w92 + w173*w186;
    const scalar w188 = 2*x;
    const scalar w189 = w187*w188;
    const scalar w190 = w112*w12*w15;
    const scalar w191 = w19*x;
    const scalar w192 = w11*w15;
    const scalar w193 = w191*w192;
    const scalar w194 = w93*x;
    const scalar w195 = w173*w184*w192;
    const scalar w196 = cos(w137);
    const scalar w197 = w140*w196;
    const scalar w198 = k*sigma0Out;
    const scalar w199 = theta1*w15;
    const scalar w200 = 4*x;
    const scalar w201 = w199*w200;
    const scalar w202 = w23*w26*w28;
    const scalar w203 = w201*w202;
    const scalar w204 = gel*w203;
    const scalar w205 = get*w203;
    const scalar w206 = -w204 + w205;
    const scalar w207 = w188*w199;
    const scalar w208 = w207*w45;
    const scalar w209 = w199*w37;
    const scalar w210 = -w188*w209 - w207*w51;
    const scalar w211 = gel*w41;
    const scalar w212 = w199*w28*w32;
    const scalar w213 = -w188*w212 + w207*w48;
    const scalar w214 = get*w45;
    const scalar w215 = -2*get*theta1*w15*w22*w26*w52*x + w121*w208 + w210*w211 + w213*w214;
    const scalar w216 = -w215;
    const scalar w217 = 2*theta1*w15*w20*w28*w32*x - w207*w42;
    const scalar w218 = -w207*w36 - w207*w38;
    const scalar w219 = -2*get*theta1*w15*w22*w26*w43*x + w117*w208 + w211*w217 + w214*w218;
    const scalar w220 = -w219;
    const scalar w221 = w1*w152;
    const scalar w222 = St*w221;
    const scalar w223 = (k*k);
    const scalar w224 = St*w223;
    const scalar w225 = w151*w2;
    const scalar w226 = w224*w225;
    const scalar w227 = -w149*w226 + w150*w222 + w154*w220 + w216*w3;
    const scalar w228 = w19*z;
    const scalar w229 = -k*sigma0Out*sigmaEps*w139*w196*w7;
    const scalar w230 = 4*theta1*w15*w20*w28*w32*x - w201*w42;
    const scalar w231 = -w201*w36 - w201*w38;
    const scalar w232 = w117*w230 + w118*w231;
    const scalar w233 = w117*w210 + w118*w213 + w121*w217 + w122*w218;
    const scalar w234 = -w233;
    const scalar w235 = St*w158;
    const scalar w236 = w223*w235;
    const scalar w237 = w148*w220 - w225*w236 + w234*w3;
    const scalar w238 = w19*y;
    const scalar w239 = -w200*w209 - w201*w51;
    const scalar w240 = -w200*w212 + w201*w48;
    const scalar w241 = w121*w239 + w122*w240;
    const scalar w242 = w221*w235;
    const scalar w243 = w148*w216 + w154*w234 + w242;
    const scalar w244 = w136*w92;
    const scalar w245 = k*sigma0In;
    const scalar w246 = w93*z;
    const scalar w247 = -k*sigma0In*sigmaEps*w139*w196*w7;
    const scalar w248 = w93*y;
    const scalar w249 = w204 - w205;
    const scalar w250 = w164*w92/(w135*w135);
    const scalar w251 = w173*w185/(w183*w183);
    const scalar w252 = w112*w113;
    const scalar w253 = w108*w15*w97;
    const scalar w254 = w76*w79;
    const scalar w255 = gil*w203 - git*w203;
    const scalar w256 = St*w9;
    const scalar w257 = (3.0/5.0)*St*w0*w81;
    const scalar w258 = w207*w41;
    const scalar w259 = gil*w41;
    const scalar w260 = git*w45;
    const scalar w261 = w208*w40 + w217*w259 + w218*w260 - w258*w44;
    const scalar w262 = w208*w50 + w210*w259 + w213*w260 - w258*w53;
    const scalar w263 = w254*w30 + w255*w77 - w256*w54 + w257*w46 - w261*w82 + w262*w83;
    const scalar w264 = w230*w40 + w231*w44;
    const scalar w265 = w210*w40 + w213*w44 + w217*w50 + w218*w53;
    const scalar w266 = w254*w46 - w256*w64 + w257*w61 + w261*w77 - w264*w82 + w265*w83;
    const scalar w267 = w254*w54;
    const scalar w268 = w257*w64;
    const scalar w269 = w239*w50 + w240*w53;
    const scalar w270 = -w256*w70 + w262*w77 - w265*w82 + w267 + w268 + w269*w83;
    const scalar w271 = w74*w92;
    const scalar w272 = w90*w92/(w73*w73);
    const scalar w273 = w185*w97/(w107*w107);
    const scalar w274 = 2*y;
    const scalar w275 = w187*w274;
    const scalar w276 = w192*w238;
    const scalar w277 = 4*y;
    const scalar w278 = w199*w277;
    const scalar w279 = w202*w278;
    const scalar w280 = gil*w279 - git*w279;
    const scalar w281 = w224*w78*w8;
    const scalar w282 = w199*w274;
    const scalar w283 = w282*w45;
    const scalar w284 = w282*w41;
    const scalar w285 = 2*theta1*w15*w20*w28*w32*y - w282*w42;
    const scalar w286 = -w282*w36 - w282*w38;
    const scalar w287 = w259*w285 + w260*w286 + w283*w40 - w284*w44;
    const scalar w288 = -w209*w274 - w282*w51;
    const scalar w289 = -w212*w274 + w282*w48;
    const scalar w290 = w259*w288 + w260*w289 + w283*w50 - w284*w53;
    const scalar w291 = w257*w54 + w280*w77 - w281*w46 - w287*w82 + w290*w83;
    const scalar w292 = 4*theta1*w15*w20*w28*w32*y - w278*w42;
    const scalar w293 = -w278*w36 - w278*w38;
    const scalar w294 = w292*w40 + w293*w44;
    const scalar w295 = w285*w50 + w286*w53 + w288*w40 + w289*w44;
    const scalar w296 = w268 - w281*w61 + w287*w77 - w294*w82 + w295*w83;
    const scalar w297 = -w209*w277 - w278*w51;
    const scalar w298 = -w212*w277 + w278*w48;
    const scalar w299 = w297*w50 + w298*w53;
    const scalar w300 = w257*w70 - w281*w64 + w290*w77 - w295*w82 + w299*w83;
    const scalar w301 = St*(fA1 - w147*w153);
    const scalar w302 = w139*w80;
    const scalar w303 = sigmaEps*w138;
    const scalar w304 = w198*w303;
    const scalar w305 = w302*w304;
    const scalar w306 = gel*w279 - get*w279;
    const scalar w307 = -2*get*theta1*w15*w22*w26*w52*y + w121*w283 + w211*w288 + w214*w289;
    const scalar w308 = -w307;
    const scalar w309 = -2*get*theta1*w15*w22*w26*w43*y + w117*w283 + w211*w285 + w214*w286;
    const scalar w310 = -w309;
    const scalar w311 = w145*w146;
    const scalar w312 = w224*w311;
    const scalar w313 = w149*w222 - w150*w312 + w154*w310 + w3*w308;
    const scalar w314 = w117*w292 + w118*w293;
    const scalar w315 = w150*w301;
    const scalar w316 = w117*w288 + w118*w289 + w121*w285 + w122*w286;
    const scalar w317 = -w316;
    const scalar w318 = w148*w310 + w242 + w3*w317 + w315;
    const scalar w319 = w121*w297 + w122*w298;
    const scalar w320 = w148*w308 + w149*w301 + w154*w317 - w236*w311;
    const scalar w321 = w245*w303;
    const scalar w322 = w302*w321;
    const scalar w323 = 2*z;
    const scalar w324 = w111*w323;
    const scalar w325 = w187*w323;
    const scalar w326 = w192*w228;
    const scalar w327 = hA2*w224*w5*x;
    const scalar w328 = 4*z;
    const scalar w329 = w199*w328;
    const scalar w330 = w202*w329;
    const scalar w331 = gil*w330 - git*w330;
    const scalar w332 = w199*w323;
    const scalar w333 = w332*w45;
    const scalar w334 = w332*w41;
    const scalar w335 = 2*theta1*w15*w20*w28*w32*z - w332*w42;
    const scalar w336 = -w332*w36 - w332*w38;
    const scalar w337 = w259*w335 + w260*w336 + w333*w40 - w334*w44;
    const scalar w338 = -w209*w323 - w332*w51;
    const scalar w339 = -w212*w323 + w332*w48;
    const scalar w340 = w259*w338 + w260*w339 + w333*w50 - w334*w53;
    const scalar w341 = w267 - w30*w327 + w331*w77 - w337*w82 + w340*w83;
    const scalar w342 = -w209*w328 - w329*w51;
    const scalar w343 = -w212*w328 + w329*w48;
    const scalar w344 = w342*w50 + w343*w53;
    const scalar w345 = w335*w50 + w336*w53 + w338*w40 + w339*w44;
    const scalar w346 = w254*w70 - w327*w54 + w340*w77 + w344*w83 - w345*w82;
    const scalar w347 = 4*theta1*w15*w20*w28*w32*z - w329*w42;
    const scalar w348 = -w329*w36 - w329*w38;
    const scalar w349 = w347*w40 + w348*w44;
    const scalar w350 = w254*w64 - w327*w46 + w337*w77 + w345*w83 - w349*w82;
    const scalar w351 = (1.0/2.0)*w7*w75;
    const scalar w352 = gel*w330;
    const scalar w353 = get*w330;
    const scalar w354 = -w352 + w353;
    const scalar w355 = (4.0/25.0)*St*w311;
    const scalar w356 = -2*get*theta1*w15*w22*w26*w52*z + w121*w333 + w211*w338 + w214*w339;
    const scalar w357 = -w356;
    const scalar w358 = -2*get*theta1*w15*w22*w26*w43*z + w117*w333 + w211*w335 + w214*w336;
    const scalar w359 = -w358;
    const scalar w360 = w154*w359 + w3*w357 + w315;
    const scalar w361 = -1.0/2.0*k*sigma0Out*sigmaEps*w138*w7*w75;
    const scalar w362 = w117*w347 + w118*w348;
    const scalar w363 = w117*w338 + w118*w339 + w121*w335 + w122*w336;
    const scalar w364 = -w363;
    const scalar w365 = w148*w359 + w150*w355 + w3*w364;
    const scalar w366 = w121*w342 + w122*w343;
    const scalar w367 = w148*w357 + w149*w355 + w154*w364 + w158*w301;
    const scalar w368 = -1.0/2.0*k*sigma0In*sigmaEps*w138*w7*w75;
    const scalar w369 = w352 - w353;
    return vector(St*w10 + 2*w111*w112*w12*w15*x + 2*w111*w12*w15*w18*x + w112*w12*w15*w18*(w110*(w194*w270 + w246*w263 + w248*w266 + w96) + w193*w91 - w194*w253 + w271*(-w191*w270 - w228*w263 - w238*w266 - w89) + w272*(-w191*(w191*w269 + w228*w262 + w238*w265 + w71) - w228*(w191*w262 + w228*w255 + w238*w261 + w55) - w238*(w191*w265 + w228*w261 + w238*w264 + w65) - w72) + w273*(w106 + w194*(-w105 - w194*w269 - w246*w262 - w248*w265) + w246*(-w100 - w194*w262 - w246*w255 - w248*w261) + w248*(-w103 - w194*w265 - w246*w261 - w248*w264))) - w113*w189 - w189*w190 - w252*(w165*w193 + w186*(-w172 - w194*(-w171*w226 + w243 + w3*(-w241 - w247)) - w246*(w148*(w197*w245 + w206) + w227) - w248*(w154*(-w232 - w247) + w169*w222 + w237)) - w194*w195 + w244*(w163 + w191*(-w161*w226 + w243 + w3*(-w229 - w241)) + w228*(w148*(w197*w198 + w206) + w227) + w238*(w154*(-w229 - w232) + w157*w222 + w237)) + w250*(-w134 - w191*(w133 + w191*w241 + w215*w228 + w233*w238) - w228*(w124 + w191*w215 + w219*w238 + w228*w249) - w238*(w130 + w191*w233 + w219*w228 + w232*w238)) + w251*(w182 + w194*(-w181 - w194*w241 - w215*w246 - w233*w248) + w246*(-w176 - w194*w215 - w219*w248 - w246*w249) + w248*(-w179 - w194*w233 - w219*w246 - w232*w248))) - w3, 2*w111*w112*w12*w15*y + 2*w111*w12*w15*w18*y + w112*w12*w15*w18*(w110*(w194*w300 + w246*w291 + w248*w296 + w95) - w248*w253 + w271*(-w191*w300 - w228*w291 - w238*w296 - w87) + w272*(-w191*(w191*w299 + w228*w290 + w238*w295 + w65) - w228*(w191*w290 + w228*w280 + w238*w287 + w47) - w238*(w191*w295 + w228*w287 + w238*w294 + w62) - w66) + w273*(w104 + w194*(-w103 - w194*w299 - w246*w290 - w248*w295) + w246*(-w194*w290 - w246*w280 - w248*w287 - w99) + w248*(-w102 - w194*w295 - w246*w287 - w248*w294)) + w276*w91) - w113*w275 - w154 - w190*w275 - w252*(w165*w276 + w186*(-w170 - w194*(w171*w222 + w3*(-w319 - w322) + w320) - w246*(w148*(-w306 - w322) + w167*w301 + w313) - w248*(w154*(-w314 - w322) - w169*w312 + w318)) - w195*w248 + w244*(w160 + w191*(w161*w222 + w3*(-w305 - w319) + w320) + w228*(w143*w301 + w148*(-w305 - w306) + w313) + w238*(w154*(-w305 - w314) - w157*w312 + w318)) + w250*(-w131 - w191*(w130 + w191*w319 + w228*w307 + w238*w316) - w228*(w120 + w191*w307 + w228*w306 + w238*w309) - w238*(w128 + w191*w316 + w228*w309 + w238*w314)) + w251*(w180 + w194*(-w179 - w194*w319 - w246*w307 - w248*w316) + w246*(-w175 - w194*w307 - w246*w306 - w248*w309) + w248*(-w178 - w194*w316 - w246*w309 - w248*w314))) - w82, w113*w324 - w113*w325 - w148 + w190*w324 - w190*w325 + w252*(w110*(w194*w346 + w246*w341 + w248*w350 + w94) - w246*w253 + w271*(-w191*w346 - w228*w341 - w238*w350 - w85) + w272*(-w191*(w191*w344 + w228*w340 + w238*w345 + w55) - w228*(w191*w340 + w228*w331 + w238*w337 + w31) - w238*(w191*w345 + w228*w337 + w238*w349 + w47) - w56) + w273*(w101 + w194*(-w100 - w194*w344 - w246*w340 - w248*w345) + w246*(-w194*w340 - w246*w331 - w248*w337 - w98) + w248*(-w194*w345 - w246*w337 - w248*w349 - w99)) + w326*w91) - w252*(w165*w326 + w186*(-w168 - w194*(w3*(-w366 - w368) + w367) - w246*(w148*(w321*w351 + w354) + w167*w355 + w360) - w248*(w154*(-w362 - w368) + w169*w301 + w365)) - w195*w246 + w244*(w156 + w191*(w3*(-w361 - w366) + w367) + w228*(w143*w355 + w148*(w304*w351 + w354) + w360) + w238*(w154*(-w361 - w362) + w157*w301 + w365)) + w250*(-w125 - w191*(w124 + w191*w366 + w228*w356 + w238*w363) - w228*(w115 + w191*w356 + w228*w369 + w238*w358) - w238*(w120 + w191*w363 + w228*w358 + w238*w362)) + w251*(w177 + w194*(-w176 - w194*w366 - w246*w356 - w248*w363) + w246*(-w174 - w194*w356 - w246*w369 - w248*w358) + w248*(-w175 - w194*w363 - w246*w358 - w248*w362))) + w77);
}

inline scalar u1(const point& pt, const scalar t, const Params& p)
{
    const scalar r1 = p.r1;
    const scalar r2 = p.r2;
    const scalar gamma1 = p.gamma1;
    const scalar gamma2 = p.gamma2;
    const scalar k = p.k;
    const scalar gil = p.gil;
    const scalar git = p.git;
    const scalar gin = p.gin;
    const scalar gel = p.gel;
    const scalar get = p.get;
    const scalar gen = p.gen;
    const scalar theta0 = p.theta0;
    const scalar theta1 = p.theta1;
    const scalar sigma0In = p.sigma0In;
    const scalar sigma0Out = p.sigma0Out;
    const scalar sigmaEps = p.sigmaEps;
    const scalar fA0 = p.fA[0];
    const scalar fA1 = p.fA[1];
    const scalar fA2 = p.fA[2];
    const scalar fA3 = p.fA[3];
    const scalar hA0 = p.hA[0];
    const scalar hA1 = p.hA[1];
    const scalar hA2 = p.hA[2];
    const scalar x = pt.x();
    const scalar y = pt.y();
    const scalar z = pt.z();
    const scalar St = sqrt(1 + t);

    const scalar w0 = k*z;
    const scalar w1 = cos(w0);
    const scalar w2 = k*x;
    const scalar w3 = sin(w2 + 1.0/5.0);
    const scalar w4 = k*y;
    const scalar w5 = sin(w0);
    const scalar w6 = (x*x);
    const scalar w7 = 1 - 3.0/10.0*w6;
    const scalar w8 = cos(w4);
    const scalar w9 = hA1*w8;
    const scalar w10 = fA1*y;
    const scalar w11 = (1.0/2.0)*y + 1;
    const scalar w12 = fA2*sin(w2);
    const scalar w13 = w4 + 3.0/10.0;
    const scalar w14 = fA3*exp((2.0/5.0)*z);
    const scalar w15 = w14*cos(w13);
    const scalar w16 = cos(gamma2);
    const scalar w17 = (w16*w16);
    const scalar w18 = sin(gamma2);
    const scalar w19 = (w18*w18);
    const scalar w20 = (r1*r1);
    const scalar w21 = -(r2*r2);
    const scalar w22 = w6 + (y*y) + (z*z);
    const scalar w23 = (-w20 + w22)/(-w20 - w21);
    const scalar w24 = theta0 + theta1*w23;
    const scalar w25 = sin(w24);
    const scalar w26 = w19*(w25*w25);
    const scalar w27 = cos(w24);
    const scalar w28 = w19*(w27*w27);
    const scalar w29 = gil*w26 + gin*w17 + git*w28;
    const scalar w30 = (1.0/(r2));
    const scalar w31 = w30*z;
    const scalar w32 = cos(gamma1);
    const scalar w33 = w16*w18;
    const scalar w34 = gin*w33;
    const scalar w35 = sin(gamma1);
    const scalar w36 = w27*w35;
    const scalar w37 = w25*w32;
    const scalar w38 = w16*w37 + w36;
    const scalar w39 = w18*w25;
    const scalar w40 = gil*w39;
    const scalar w41 = w25*w35;
    const scalar w42 = w16*w27*w32 - w41;
    const scalar w43 = w18*w27;
    const scalar w44 = git*w43;
    const scalar w45 = -w32*w34 + w38*w40 + w42*w44;
    const scalar w46 = w30*y;
    const scalar w47 = -w16*w41 + w27*w32;
    const scalar w48 = -w16*w36 - w37;
    const scalar w49 = w34*w35 + w40*w47 + w44*w48;
    const scalar w50 = w30*x;
    const scalar w51 = (w32*w32);
    const scalar w52 = gin*w19;
    const scalar w53 = (w38*w38);
    const scalar w54 = (w42*w42);
    const scalar w55 = gil*w53 + git*w54 + w51*w52;
    const scalar w56 = w32*w35;
    const scalar w57 = w38*w47;
    const scalar w58 = w42*w48;
    const scalar w59 = gil*w57 + git*w58 - w52*w56;
    const scalar w60 = (w35*w35);
    const scalar w61 = (w47*w47);
    const scalar w62 = (w48*w48);
    const scalar w63 = gil*w61 + git*w62 + w52*w60;
    const scalar w64 = St*(hA0 + hA2*w1*w2);
    const scalar w65 = St*hA1*k*w7*sin(w4);
    const scalar w66 = St*(hA2*w5 - 3.0/5.0*w9*x);
    const scalar w67 = w29*w64 - w45*w65 + w49*w66;
    const scalar w68 = w45*w64 - w55*w65 + w59*w66;
    const scalar w69 = w49*w64 - w59*w65 + w63*w66;
    const scalar w70 = r1 + r2;
    const scalar w71 = (1.0/2.0)*w70;
    const scalar w72 = w23*w30*w71;
    const scalar w73 = (1.0/(r1));
    const scalar w74 = w73*z;
    const scalar w75 = w73*y;
    const scalar w76 = w73*x;
    const scalar w77 = w71*w73*(1 - w23);
    const scalar w78 = w23*(w21 + w22)/w70;
    const scalar w79 = gel*w26 + gen*w17 + get*w28;
    const scalar w80 = gen*w33;
    const scalar w81 = gel*w39;
    const scalar w82 = get*w43;
    const scalar w83 = -w32*w80 + w38*w81 + w42*w82;
    const scalar w84 = w35*w80 + w47*w81 + w48*w82;
    const scalar w85 = gen*w19;
    const scalar w86 = gel*w53 + get*w54 + w51*w85;
    const scalar w87 = gel*w57 + get*w58 - w56*w85;
    const scalar w88 = gel*w61 + get*w62 + w60*w85;
    const scalar w89 = sigmaEps*w3*w8*((1.0/2.0)*w5 + 1) + 1;
    const scalar w90 = -sigma0Out*w89;
    const scalar w91 = St*(w10 + (2.0/5.0)*w15);
    const scalar w92 = -w84;
    const scalar w93 = St*(fA0 + fA2*k*w11*cos(w2));
    const scalar w94 = -w83;
    const scalar w95 = St*(fA1*z - k*w14*sin(w13) + (1.0/2.0)*w12);
    const scalar w96 = w92*w93 + w94*w95;
    const scalar w97 = -w87;
    const scalar w98 = w91*w94 + w93*w97;
    const scalar w99 = w91*w92 + w95*w97;
    const scalar w100 = -sigma0In*w89;
    return St*(hA0*z + hA2*w5*x + w7*w9) - St*(fA0*x + w10*z + w11*w12 + w15) + w78*(w72*(-w31*w67 - w46*w68 - w50*w69)/(w31*(w29*w31 + w45*w46 + w49*w50) + w46*(w31*w45 + w46*w55 + w50*w59) + w50*(w31*w49 + w46*w59 + w50*w63)) + w77*(w67*w74 + w68*w75 + w69*w76)/(-w74*(-w29*w74 - w45*w75 - w49*w76) - w75*(-w45*w74 - w55*w75 - w59*w76) - w76*(-w49*w74 - w59*w75 - w63*w76))) - w78*(w72*(w31*(w91*(-w79 - w90) + w96) + w46*(w95*(-w86 - w90) + w98) + w50*(w93*(-w88 - w90) + w99))/(w31*(w31*w79 + w46*w83 + w50*w84) + w46*(w31*w83 + w46*w86 + w50*w87) + w50*(w31*w84 + w46*w87 + w50*w88)) + w77*(-w74*(w91*(-w100 - w79) + w96) - w75*(w95*(-w100 - w86) + w98) - w76*(w93*(-w100 - w88) + w99))/(-w74*(-w74*w79 - w75*w83 - w76*w84) - w75*(-w74*w83 - w75*w86 - w76*w87) - w76*(-w74*w84 - w75*w87 - w76*w88))) + (t + 1)*((1.0/4.0)*w1*w3*cos(w4 - 1.0/3.0) + 1);
}

inline scalar u2(const point& pt, const scalar t, const Params& p)
{
    const scalar k = p.k;
    const scalar x = pt.x();
    const scalar y = pt.y();
    const scalar z = pt.z();

    return 1/((t + 1)*sqrt((1.0/4.0)*sin(k*x + 1.0/5.0)*cos(k*z)*cos(k*y - 1.0/3.0) + 1));
}

inline symmTensor Gi(const point& pt, const Params& p)
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

inline symmTensor Ge(const point& pt, const Params& p)
{
    const scalar r1 = p.r1;
    const scalar r2 = p.r2;
    const scalar gamma1 = p.gamma1;
    const scalar gamma2 = p.gamma2;
    const scalar gel = p.gel;
    const scalar get = p.get;
    const scalar gen = p.gen;
    const scalar theta0 = p.theta0;
    const scalar theta1 = p.theta1;
    const scalar x = pt.x();
    const scalar y = pt.y();
    const scalar z = pt.z();

    const scalar w0 = sin(gamma1);
    const scalar w1 = sin(gamma2);
    const scalar w2 = (w1*w1);
    const scalar w3 = gen*w2;
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
    const scalar w16 = gel*w11;
    const scalar w17 = -w10 + w4*w7*w8;
    const scalar w18 = get*w14;
    const scalar w19 = gen*w1*w8;
    const scalar w20 = w1*w9;
    const scalar w21 = w1*w7;
    return symmTensor(gel*(w11*w11) + get*(w14*w14) + (w0*w0)*w3, -w0*w3*w4 + w15*w16 + w17*w18, w0*w19 + w16*w20 + w18*w21, gel*(w15*w15) + get*(w17*w17) + w3*(w4*w4), gel*w15*w20 + get*w17*w21 - w19*w4, gel*w2*(w9*w9) + gen*(w8*w8) + get*w2*(w7*w7));
}

inline scalar sigmaB(const point& pt, const label zone, const Params& p)
{
    const scalar k = p.k;
    const scalar sigma0In = p.sigma0In;
    const scalar sigma0Out = p.sigma0Out;
    const scalar sigmaEps = p.sigmaEps;
    const scalar x = pt.x();
    const scalar y = pt.y();
    const scalar z = pt.z();
    checkZone(zone, false);

    if (zone == bloodPool)
    {
        return sigma0In*(sigmaEps*((1.0/2.0)*sin(k*z) + 1)*sin(k*x + 1.0/5.0)*cos(k*y) + 1);
    }
    if (zone == bath)
    {
        return sigma0Out*(sigmaEps*((1.0/2.0)*sin(k*z) + 1)*sin(k*x + 1.0/5.0)*cos(k*y) + 1);
    }
    return 0;
}

inline scalar divGiGradPhiI(const point& pt, const scalar t, const Params& p)
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
    const scalar cmAmplitude = p.cmAmplitude;
    const scalar cmOmega = p.cmOmega;
    const scalar x = pt.x();
    const scalar y = pt.y();
    const scalar z = pt.z();
    const scalar St = sqrt(1 + t);
    const scalar Tt = cmAmplitude*sin(cmOmega*t);

    const scalar w0 = sin(gamma2);
    const scalar w1 = (w0*w0);
    const scalar w2 = gil*w1;
    const scalar w3 = (r1*r1);
    const scalar w4 = -(r2*r2);
    const scalar w5 = -w3 - w4;
    const scalar w6 = (1.0/(w5));
    const scalar w7 = theta1*w6;
    const scalar w8 = 4*w7;
    const scalar w9 = (x*x);
    const scalar w10 = (y*y);
    const scalar w11 = (z*z);
    const scalar w12 = w10 + w11 + w9;
    const scalar w13 = w6*(w12 - w3);
    const scalar w14 = theta0 + theta1*w13;
    const scalar w15 = cos(w14);
    const scalar w16 = sin(w14);
    const scalar w17 = w15*w16;
    const scalar w18 = w17*w8;
    const scalar w19 = w18*w2;
    const scalar w20 = git*w1;
    const scalar w21 = w18*w20;
    const scalar w22 = w19*z - w21*z;
    const scalar w23 = k*x;
    const scalar w24 = k*z;
    const scalar w25 = cos(w24);
    const scalar w26 = hA2*w25;
    const scalar w27 = hA0 + w23*w26;
    const scalar w28 = St*w27;
    const scalar w29 = r1 + r2;
    const scalar w30 = (1.0/(w29));
    const scalar w31 = w13*w30;
    const scalar w32 = (1.0/(r2));
    const scalar w33 = cos(gamma2);
    const scalar w34 = gin*(w33*w33);
    const scalar w35 = (w16*w16);
    const scalar w36 = w2*w35;
    const scalar w37 = (w15*w15);
    const scalar w38 = w20*w37;
    const scalar w39 = w34 + w36 + w38;
    const scalar w40 = cos(gamma1);
    const scalar w41 = gin*w0*w33;
    const scalar w42 = w40*w41;
    const scalar w43 = sin(gamma1);
    const scalar w44 = w15*w43;
    const scalar w45 = w16*w40;
    const scalar w46 = w33*w45;
    const scalar w47 = w44 + w46;
    const scalar w48 = gil*w47;
    const scalar w49 = w0*w16;
    const scalar w50 = w48*w49;
    const scalar w51 = w16*w43;
    const scalar w52 = w15*w33*w40 - w51;
    const scalar w53 = git*w52;
    const scalar w54 = w0*w15;
    const scalar w55 = w53*w54;
    const scalar w56 = -w42 + w50 + w55;
    const scalar w57 = 1 - 3.0/10.0*w9;
    const scalar w58 = St*k;
    const scalar w59 = k*y;
    const scalar w60 = hA1*sin(w59);
    const scalar w61 = w58*w60;
    const scalar w62 = w57*w61;
    const scalar w63 = w41*w43;
    const scalar w64 = w33*w51;
    const scalar w65 = w15*w40 - w64;
    const scalar w66 = gil*w65;
    const scalar w67 = w49*w66;
    const scalar w68 = w33*w44;
    const scalar w69 = -w45 - w68;
    const scalar w70 = git*w69;
    const scalar w71 = w54*w70;
    const scalar w72 = w63 + w67 + w71;
    const scalar w73 = sin(w24);
    const scalar w74 = hA1*cos(w59);
    const scalar w75 = (3.0/5.0)*w74;
    const scalar w76 = w75*x;
    const scalar w77 = hA2*w73 - w76;
    const scalar w78 = St*w77;
    const scalar w79 = w28*w39 - w56*w62 + w72*w78;
    const scalar w80 = w32*w79;
    const scalar w81 = gin*w1;
    const scalar w82 = (w40*w40)*w81;
    const scalar w83 = gil*(w47*w47);
    const scalar w84 = git*(w52*w52);
    const scalar w85 = w82 + w83 + w84;
    const scalar w86 = w40*w43*w81;
    const scalar w87 = w47*w66;
    const scalar w88 = w52*w70;
    const scalar w89 = -w86 + w87 + w88;
    const scalar w90 = w28*w56 - w62*w85 + w78*w89;
    const scalar w91 = w32*w90;
    const scalar w92 = gil*(w65*w65) + git*(w69*w69) + (w43*w43)*w81;
    const scalar w93 = w28*w72 - w62*w89 + w78*w92;
    const scalar w94 = w32*w93;
    const scalar w95 = -w80*z - w91*y - w94*x;
    const scalar w96 = w13*w29*w32;
    const scalar w97 = w32*w39;
    const scalar w98 = w32*w56;
    const scalar w99 = w32*w72;
    const scalar w100 = w32*(w97*z + w98*y + w99*x);
    const scalar w101 = w32*w85;
    const scalar w102 = w32*w89;
    const scalar w103 = w32*(w101*y + w102*x + w98*z);
    const scalar w104 = w32*w92;
    const scalar w105 = w32*(w102*y + w104*x + w99*z);
    const scalar w106 = w100*z + w103*y + w105*x;
    const scalar w107 = (1.0/(w106));
    const scalar w108 = (1.0/2.0)*w107;
    const scalar w109 = w108*w96;
    const scalar w110 = (1.0/(r1));
    const scalar w111 = w110*w79;
    const scalar w112 = w110*w90;
    const scalar w113 = w110*w93;
    const scalar w114 = w111*z + w112*y + w113*x;
    const scalar w115 = w29*(1 - w13);
    const scalar w116 = w110*w39;
    const scalar w117 = w110*w56;
    const scalar w118 = w110*w72;
    const scalar w119 = w110*(-w116*z - w117*y - w118*x);
    const scalar w120 = w110*w85;
    const scalar w121 = w110*w89;
    const scalar w122 = w110*(-w117*z - w120*y - w121*x);
    const scalar w123 = w110*w92;
    const scalar w124 = w110*(-w118*z - w121*y - w123*x);
    const scalar w125 = -w119*z - w122*y - w124*x;
    const scalar w126 = (1.0/(w125));
    const scalar w127 = w110*w126;
    const scalar w128 = (1.0/2.0)*w127;
    const scalar w129 = w115*w128;
    const scalar w130 = w109*w95 + w114*w129;
    const scalar w131 = 2*w130;
    const scalar w132 = w131*w31;
    const scalar w133 = w12 + w4;
    const scalar w134 = w30*w6;
    const scalar w135 = w133*w134;
    const scalar w136 = w131*w135;
    const scalar w137 = w32*z;
    const scalar w138 = w29*w6;
    const scalar w139 = w107*w138;
    const scalar w140 = w139*w95;
    const scalar w141 = w114*w138;
    const scalar w142 = w110*z;
    const scalar w143 = w126*w142;
    const scalar w144 = (k*k);
    const scalar w145 = hA2*w144*w73;
    const scalar w146 = St*w145;
    const scalar w147 = w146*x;
    const scalar w148 = w26*w58;
    const scalar w149 = w148*w72;
    const scalar w150 = 2*w7;
    const scalar w151 = w150*w54;
    const scalar w152 = w151*w48;
    const scalar w153 = w150*w49;
    const scalar w154 = w153*w53;
    const scalar w155 = w150*w51;
    const scalar w156 = 2*theta1*w15*w33*w40*w6*z - w155*z;
    const scalar w157 = gil*w156;
    const scalar w158 = w157*w49;
    const scalar w159 = w150*w44;
    const scalar w160 = w150*w46;
    const scalar w161 = -w159*z - w160*z;
    const scalar w162 = git*w161;
    const scalar w163 = w162*w54;
    const scalar w164 = w152*z - w154*z + w158 + w163;
    const scalar w165 = w151*w66;
    const scalar w166 = w153*w70;
    const scalar w167 = w150*w45;
    const scalar w168 = w150*w68;
    const scalar w169 = -w167*z - w168*z;
    const scalar w170 = gil*w49;
    const scalar w171 = w169*w170;
    const scalar w172 = w15*w40;
    const scalar w173 = w150*w172;
    const scalar w174 = w150*w64;
    const scalar w175 = -w173*z + w174*z;
    const scalar w176 = git*w54;
    const scalar w177 = w175*w176;
    const scalar w178 = w165*z - w166*z + w171 + w177;
    const scalar w179 = -w147*w39 + w149 - w164*w62 + w178*w78 + w22*w28;
    const scalar w180 = w146*w72;
    const scalar w181 = w45*w8;
    const scalar w182 = w68*w8;
    const scalar w183 = -w181*z - w182*z;
    const scalar w184 = w172*w8;
    const scalar w185 = w64*w8;
    const scalar w186 = -w184*z + w185*z;
    const scalar w187 = w183*w66 + w186*w70;
    const scalar w188 = w169*w48;
    const scalar w189 = w156*w66;
    const scalar w190 = w175*w53;
    const scalar w191 = w161*w70;
    const scalar w192 = w188 + w189 + w190 + w191;
    const scalar w193 = w148*w92 + w178*w28 - w180*x + w187*w78 - w192*w62;
    const scalar w194 = w193*w32;
    const scalar w195 = w51*w8;
    const scalar w196 = 4*theta1*w15*w33*w40*w6*z - w195*z;
    const scalar w197 = w196*w48;
    const scalar w198 = w44*w8;
    const scalar w199 = w46*w8;
    const scalar w200 = -w198*z - w199*z;
    const scalar w201 = w200*w53;
    const scalar w202 = w197 + w201;
    const scalar w203 = -w147*w56 + w148*w89 + w164*w28 + w192*w78 - w202*w62;
    const scalar w204 = w203*w32;
    const scalar w205 = -w137*w179 - w194*x - w204*y - w80;
    const scalar w206 = w110*w193;
    const scalar w207 = w110*w203;
    const scalar w208 = w111 + w142*w179 + w206*x + w207*y;
    const scalar w209 = w164*w32;
    const scalar w210 = w178*w32;
    const scalar w211 = w137*w22 + w209*y + w210*x + w97;
    const scalar w212 = w137*w211;
    const scalar w213 = w202*w32;
    const scalar w214 = w192*w32;
    const scalar w215 = w32*(w209*z + w213*y + w214*x + w98);
    const scalar w216 = w215*y;
    const scalar w217 = w187*w32;
    const scalar w218 = w32*(w210*z + w214*y + w217*x + w99);
    const scalar w219 = w218*x;
    const scalar w220 = -w100 - w212 - w216 - w219;
    const scalar w221 = (1.0/(w106*w106));
    const scalar w222 = w221*w96;
    const scalar w223 = w220*w222;
    const scalar w224 = (1.0/2.0)*w95;
    const scalar w225 = w110*w164;
    const scalar w226 = w110*w178;
    const scalar w227 = -w116 - w142*w22 - w225*y - w226*x;
    const scalar w228 = w142*w227;
    const scalar w229 = w110*w202;
    const scalar w230 = w110*w192;
    const scalar w231 = w110*(-w117 - w225*z - w229*y - w230*x);
    const scalar w232 = w231*y;
    const scalar w233 = w110*w187;
    const scalar w234 = w110*(-w118 - w226*z - w230*y - w233*x);
    const scalar w235 = w234*x;
    const scalar w236 = w119 + w228 + w232 + w235;
    const scalar w237 = (1.0/(w125*w125));
    const scalar w238 = w110*w237;
    const scalar w239 = w115*w238;
    const scalar w240 = w236*w239;
    const scalar w241 = (1.0/2.0)*w240;
    const scalar w242 = w109*w205 + w114*w241 + w129*w208 + w137*w140 - w141*w143 + w223*w224;
    const scalar w243 = w133*w31;
    const scalar w244 = (1.0/5.0)*z;
    const scalar w245 = sin(w23);
    const scalar w246 = (1.0/5.0)*w245;
    const scalar w247 = k*w246;
    const scalar w248 = w244*w39 - w247*w72 - 1.0/10.0*w42 + (1.0/10.0)*w50 + (1.0/10.0)*w55;
    const scalar w249 = w248*w32;
    const scalar w250 = w244*w56 - w247*w89 + (1.0/10.0)*w82 + (1.0/10.0)*w83 + (1.0/10.0)*w84;
    const scalar w251 = w250*w32;
    const scalar w252 = w244*w72 - w247*w92 - 1.0/10.0*w86 + (1.0/10.0)*w87 + (1.0/10.0)*w88;
    const scalar w253 = w252*w32;
    const scalar w254 = -w249*z - w251*y - w253*x;
    const scalar w255 = w254*w32;
    const scalar w256 = w13*w255*w29;
    const scalar w257 = w110*w248;
    const scalar w258 = w110*w250;
    const scalar w259 = w110*w252;
    const scalar w260 = w257*z + w258*y + w259*x;
    const scalar w261 = w108*w256 + w129*w260;
    const scalar w262 = 2*w261;
    const scalar w263 = w262*w31;
    const scalar w264 = w135*w262;
    const scalar w265 = w139*w254;
    const scalar w266 = w138*w260;
    const scalar w267 = w48*w54;
    const scalar w268 = w244*w7;
    const scalar w269 = w49*w53;
    const scalar w270 = (1.0/10.0)*w158 + (1.0/10.0)*w163 - w178*w247 + w22*w244 + w267*w268 - w268*w269 + (1.0/5.0)*w34 + (1.0/5.0)*w36 + (1.0/5.0)*w38;
    const scalar w271 = w164*w244 - w192*w247 + (1.0/10.0)*w197 + (1.0/10.0)*w201 - 1.0/5.0*w42 + (1.0/5.0)*w50 + (1.0/5.0)*w55;
    const scalar w272 = w271*w32;
    const scalar w273 = w178*w244 - w187*w247 + (1.0/10.0)*w188 + (1.0/10.0)*w189 + (1.0/10.0)*w190 + (1.0/10.0)*w191 + (1.0/5.0)*w63 + (1.0/5.0)*w67 + (1.0/5.0)*w71;
    const scalar w274 = w273*w32;
    const scalar w275 = -w137*w270 - w249 - w272*y - w274*x;
    const scalar w276 = w110*w271;
    const scalar w277 = w110*w273;
    const scalar w278 = w142*w270 + w257 + w276*y + w277*x;
    const scalar w279 = w220*w221;
    const scalar w280 = (1.0/2.0)*w256;
    const scalar w281 = w109*w275 + w129*w278 + w137*w265 - w143*w266 + w241*w260 + w279*w280;
    const scalar w282 = Tt*(w243*w281 + w244 + w263*z + w264*z) + w132*z + w136*z + w242*w243 + w28;
    const scalar w283 = 2*theta1*w15*w33*w40*w6*y - w155*y;
    const scalar w284 = gil*w283;
    const scalar w285 = w284*w49;
    const scalar w286 = -w159*y - w160*y;
    const scalar w287 = git*w286;
    const scalar w288 = w287*w54;
    const scalar w289 = w152*y - w154*y + w285 + w288;
    const scalar w290 = -w167*x - w168*x;
    const scalar w291 = w170*w290;
    const scalar w292 = -w173*x + w174*x;
    const scalar w293 = w176*w292;
    const scalar w294 = w165*x - w166*x + w291 + w293;
    const scalar w295 = 4*theta1*w15*w33*w40*w6*y - w195*y;
    const scalar w296 = w295*w48;
    const scalar w297 = -w198*y - w199*y;
    const scalar w298 = w297*w53;
    const scalar w299 = w296 + w298;
    const scalar w300 = w32*y;
    const scalar w301 = w110*y;
    const scalar w302 = w126*w301;
    const scalar w303 = w19*y - w21*y;
    const scalar w304 = (3.0/5.0)*St*w23;
    const scalar w305 = w304*w60;
    const scalar w306 = St*w74;
    const scalar w307 = w144*w306*w57;
    const scalar w308 = -w167*y - w168*y;
    const scalar w309 = w170*w308;
    const scalar w310 = -w173*y + w174*y;
    const scalar w311 = w176*w310;
    const scalar w312 = w165*y - w166*y + w309 + w311;
    const scalar w313 = w28*w303 - w289*w62 + w305*w72 - w307*w56 + w312*w78;
    const scalar w314 = w313*w32;
    const scalar w315 = w60*w89;
    const scalar w316 = w304*w315;
    const scalar w317 = w308*w48;
    const scalar w318 = w283*w66;
    const scalar w319 = w310*w53;
    const scalar w320 = w286*w70;
    const scalar w321 = w317 + w318 + w319 + w320;
    const scalar w322 = w28*w289 - w299*w62 - w307*w85 + w316 + w321*w78;
    const scalar w323 = St*w92;
    const scalar w324 = w23*w60;
    const scalar w325 = -w181*y - w182*y;
    const scalar w326 = -w184*y + w185*y;
    const scalar w327 = w325*w66 + w326*w70;
    const scalar w328 = w28*w312 - w307*w89 - w321*w62 + (3.0/5.0)*w323*w324 + w327*w78;
    const scalar w329 = w32*w328;
    const scalar w330 = -w300*w322 - w314*z - w329*x - w91;
    const scalar w331 = w110*w313;
    const scalar w332 = w110*w328;
    const scalar w333 = w112 + w301*w322 + w331*z + w332*x;
    const scalar w334 = w303*w32;
    const scalar w335 = w289*w32;
    const scalar w336 = w312*w32;
    const scalar w337 = w32*(w334*z + w335*y + w336*x + w98);
    const scalar w338 = w337*z;
    const scalar w339 = w32*w321;
    const scalar w340 = w101 + w299*w300 + w335*z + w339*x;
    const scalar w341 = w300*w340;
    const scalar w342 = w32*w327;
    const scalar w343 = w32*(w102 + w336*z + w339*y + w342*x);
    const scalar w344 = w343*x;
    const scalar w345 = -w103 - w338 - w341 - w344;
    const scalar w346 = w222*w345;
    const scalar w347 = w110*w303;
    const scalar w348 = w110*w289;
    const scalar w349 = w110*w312;
    const scalar w350 = w110*(-w117 - w347*z - w348*y - w349*x);
    const scalar w351 = w350*z;
    const scalar w352 = w110*w321;
    const scalar w353 = -w120 - w299*w301 - w348*z - w352*x;
    const scalar w354 = w301*w353;
    const scalar w355 = w110*w327;
    const scalar w356 = w110*(-w121 - w349*z - w352*y - w355*x);
    const scalar w357 = w356*x;
    const scalar w358 = w122 + w351 + w354 + w357;
    const scalar w359 = w239*w358;
    const scalar w360 = (1.0/2.0)*w359;
    const scalar w361 = w109*w330 + w114*w360 + w129*w333 + w140*w300 - w141*w302 + w224*w346;
    const scalar w362 = (1.0/5.0)*w7;
    const scalar w363 = w267*w362;
    const scalar w364 = w269*w362;
    const scalar w365 = w244*w303 - w247*w312 + (1.0/10.0)*w285 + (1.0/10.0)*w288 + w363*y - w364*y;
    const scalar w366 = w32*w365;
    const scalar w367 = w244*w289 - w247*w321 + (1.0/10.0)*w296 + (1.0/10.0)*w298;
    const scalar w368 = w244*w312 - w247*w327 + (1.0/10.0)*w317 + (1.0/10.0)*w318 + (1.0/10.0)*w319 + (1.0/10.0)*w320;
    const scalar w369 = w32*w368;
    const scalar w370 = -w251 - w300*w367 - w366*z - w369*x;
    const scalar w371 = w110*w365;
    const scalar w372 = w110*w368;
    const scalar w373 = w258 + w301*w367 + w371*z + w372*x;
    const scalar w374 = w221*w345;
    const scalar w375 = w109*w370 + w129*w373 + w260*w360 + w265*w300 - w266*w302 + w280*w374;
    const scalar w376 = Tt*(w243*w375 + w263*y + w264*y + 1.0/10.0) + w132*y + w136*y + w243*w361 - w62;
    const scalar w377 = w290*w48;
    const scalar w378 = 2*theta1*w15*w33*w40*w6*x - w155*x;
    const scalar w379 = w378*w66;
    const scalar w380 = w292*w53;
    const scalar w381 = -w159*x - w160*x;
    const scalar w382 = w381*w70;
    const scalar w383 = w377 + w379 + w380 + w382;
    const scalar w384 = -w181*x - w182*x;
    const scalar w385 = -w184*x + w185*x;
    const scalar w386 = w384*w66 + w385*w70;
    const scalar w387 = w32*x;
    const scalar w388 = w110*x;
    const scalar w389 = w126*w388;
    const scalar w390 = w19*x - w21*x;
    const scalar w391 = St*w75;
    const scalar w392 = gil*w378;
    const scalar w393 = w392*w49;
    const scalar w394 = git*w381;
    const scalar w395 = w394*w54;
    const scalar w396 = w152*x - w154*x + w393 + w395;
    const scalar w397 = w148*w39 + w28*w390 + w294*w78 + w305*w56 - w391*w72 - w396*w62;
    const scalar w398 = w32*w397;
    const scalar w399 = 4*theta1*w15*w33*w40*w6*x - w195*x;
    const scalar w400 = w399*w48;
    const scalar w401 = -w198*x - w199*x;
    const scalar w402 = w401*w53;
    const scalar w403 = w400 + w402;
    const scalar w404 = w148*w56 + w28*w396 + w305*w85 + w383*w78 - w391*w89 - w403*w62;
    const scalar w405 = w32*w404;
    const scalar w406 = w149 + w28*w294 + w316 - w383*w62 + w386*w78 - w391*w92;
    const scalar w407 = -w387*w406 - w398*z - w405*y - w94;
    const scalar w408 = w110*w397;
    const scalar w409 = w110*w404;
    const scalar w410 = w113 + w388*w406 + w408*z + w409*y;
    const scalar w411 = w32*w390;
    const scalar w412 = w32*w396;
    const scalar w413 = w294*w32;
    const scalar w414 = w32*(w411*z + w412*y + w413*x + w99);
    const scalar w415 = w414*z;
    const scalar w416 = w32*w383;
    const scalar w417 = w104 + w386*w387 + w413*z + w416*y;
    const scalar w418 = w387*w417;
    const scalar w419 = w32*w403;
    const scalar w420 = w32*(w102 + w412*z + w416*x + w419*y);
    const scalar w421 = w420*y;
    const scalar w422 = -w105 - w415 - w418 - w421;
    const scalar w423 = w221*w422;
    const scalar w424 = w423*w96;
    const scalar w425 = w110*w390;
    const scalar w426 = w110*w396;
    const scalar w427 = w110*w294;
    const scalar w428 = w110*(-w118 - w425*z - w426*y - w427*x);
    const scalar w429 = w428*z;
    const scalar w430 = w110*w383;
    const scalar w431 = -w123 - w386*w388 - w427*z - w430*y;
    const scalar w432 = w388*w431;
    const scalar w433 = w110*w403;
    const scalar w434 = w110*(-w121 - w426*z - w430*x - w433*y);
    const scalar w435 = w434*y;
    const scalar w436 = w124 + w429 + w432 + w435;
    const scalar w437 = w238*w436;
    const scalar w438 = w115*w437;
    const scalar w439 = (1.0/2.0)*w438;
    const scalar w440 = w109*w407 + w114*w439 + w129*w410 + w140*w387 - w141*w389 + w224*w424;
    const scalar w441 = w144*cos(w23);
    const scalar w442 = (1.0/5.0)*w441;
    const scalar w443 = w244*w390 - w247*w294 + w363*x - w364*x + (1.0/10.0)*w393 + (1.0/10.0)*w395 - w442*w72;
    const scalar w444 = w32*w443;
    const scalar w445 = w244*w294 - w247*w386 + (1.0/10.0)*w377 + (1.0/10.0)*w379 + (1.0/10.0)*w380 + (1.0/10.0)*w382 - w442*w92;
    const scalar w446 = w244*w396 - w247*w383 + (1.0/10.0)*w400 + (1.0/10.0)*w402 - w442*w89;
    const scalar w447 = w32*w446;
    const scalar w448 = -w253 - w387*w445 - w444*z - w447*y;
    const scalar w449 = w448*w96;
    const scalar w450 = w110*w443;
    const scalar w451 = w110*w446;
    const scalar w452 = w259 + w388*w445 + w450*z + w451*y;
    const scalar w453 = w115*w452;
    const scalar w454 = w108*w449 + w128*w453 + w260*w439 + w265*w387 - w266*w389 + w280*w423;
    const scalar w455 = Tt*(w243*w454 - w247 + w263*x + w264*x) + w132*x + w136*x + w243*w440 + w78;
    const scalar w456 = 8*w134;
    const scalar w457 = w11*w456;
    const scalar w458 = 2*w110;
    const scalar w459 = (theta1*theta1);
    const scalar w460 = (1.0/(w5*w5));
    const scalar w461 = w459*w460;
    const scalar w462 = w11*w461;
    const scalar w463 = 8*w462;
    const scalar w464 = -4*gil*theta1*w1*w15*w16*w6 + w21;
    const scalar w465 = 8*gil*w1*w11*w37*w459*w460 + 8*git*w1*w11*w35*w459*w460 - w36*w463 - w38*w463 - w464;
    const scalar w466 = 4*w462;
    const scalar w467 = w167 + w168;
    const scalar w468 = 4*w11*w16*w33*w43*w459*w460 - w172*w466 - w467;
    const scalar w469 = -w173 + w174;
    const scalar w470 = w45*w466 + w466*w68 + w469;
    const scalar w471 = git*w175;
    const scalar w472 = w49*w8*z;
    const scalar w473 = -2*gil*theta1*w0*w15*w6*w65 + w166;
    const scalar w474 = 4*gil*theta1*w0*w15*w169*w6*z + gil*w0*w16*w468 + git*w0*w15*w470 - w466*w67 - w466*w71 - w471*w472 - w473;
    const scalar w475 = -2*theta1*w15*w33*w40*w6 + w155;
    const scalar w476 = -w44*w466 - w46*w466 - w475;
    const scalar w477 = w172*w33;
    const scalar w478 = w159 + w160;
    const scalar w479 = 4*w11*w16*w43*w459*w460 - w466*w477 - w478;
    const scalar w480 = (2.0/5.0)*w462;
    const scalar w481 = (2.0/5.0)*w7;
    const scalar w482 = w481*w49;
    const scalar w483 = w17*w20;
    const scalar w484 = w7*z;
    const scalar w485 = -1.0/5.0*gil*theta1*w0*w15*w47*w6 + w364;
    const scalar w486 = (2.0/5.0)*gil*theta1*w0*w15*w156*w6*z + (8.0/5.0)*gil*theta1*w1*w15*w16*w6*z + (1.0/10.0)*gil*w0*w16*w476 + (1.0/10.0)*git*w0*w15*w479 - w162*w482*z - w247*w474 + (1.0/5.0)*w465*z - w480*w50 - w480*w55 - 8.0/5.0*w483*w484 - w485;
    const scalar w487 = (4.0/5.0)*w484;
    const scalar w488 = w157*w196;
    const scalar w489 = w162*w200;
    const scalar w490 = -4*theta1*w15*w33*w40*w6 + w195;
    const scalar w491 = w48*(-w44*w463 - w46*w463 - w490);
    const scalar w492 = w198 + w199;
    const scalar w493 = w53*(8*w11*w16*w43*w459*w460 - w463*w477 - w492);
    const scalar w494 = -2*gil*theta1*w0*w15*w47*w6 + w154;
    const scalar w495 = 4*gil*theta1*w0*w15*w156*w6*z + gil*w0*w16*w476 + git*w0*w15*w479 - w162*w472 - w466*w50 - w466*w55 - w494;
    const scalar w496 = w157*w169;
    const scalar w497 = w162*w175;
    const scalar w498 = w468*w48;
    const scalar w499 = w476*w66;
    const scalar w500 = w470*w53;
    const scalar w501 = w479*w70;
    const scalar w502 = 2*w496 + 2*w497 + w498 + w499 + w500 + w501;
    const scalar w503 = (2.0/5.0)*w158 + (2.0/5.0)*w163 + w244*w495 - w247*w502 + w267*w487 - w269*w487 + (1.0/10.0)*w488 + (1.0/10.0)*w489 + (1.0/10.0)*w491 + (1.0/10.0)*w493;
    const scalar w504 = gil*w183;
    const scalar w505 = git*w186;
    const scalar w506 = w181 + w182;
    const scalar w507 = -w184 + w185;
    const scalar w508 = w169*w504 + w175*w505 + w66*(8*w11*w16*w33*w43*w459*w460 - w172*w463 - w506) + w70*(w45*w463 + w463*w68 + w507);
    const scalar w509 = (2.0/5.0)*w171 + (2.0/5.0)*w177 + w244*w474 - w247*w508 - w487*w49*w70 + w487*w54*w66 + (1.0/5.0)*w496 + (1.0/5.0)*w497 + (1.0/10.0)*w498 + (1.0/10.0)*w499 + (1.0/10.0)*w500 + (1.0/10.0)*w501;
    const scalar w510 = w488 + w489 + w491 + w493;
    const scalar w511 = (1.0/2.0)*w239;
    const scalar w512 = w511*(w142*(-w142*w465 - w22*w458 - w301*w495 - w388*w474) + w227*w458 + w301*(-w142*w495 - 2*w225 - w301*w510 - w388*w502) + w388*(-w142*w474 - 2*w226 - w301*w502 - w388*w508));
    const scalar w513 = 2*w32;
    const scalar w514 = w138*w143;
    const scalar w515 = 2*w514;
    const scalar w516 = w137*w139;
    const scalar w517 = 2*w516;
    const scalar w518 = w110*w115/(w125*w125*w125);
    const scalar w519 = w260*w518;
    const scalar w520 = 2*w119 + 2*w228 + 2*w232 + 2*w235;
    const scalar w521 = (1.0/2.0)*w236;
    const scalar w522 = w520*w521;
    const scalar w523 = -w137*(w137*w465 + w22*w513 + w300*w495 + w387*w474) - w211*w513 - w300*(w137*w495 + 2*w209 + w300*w510 + w387*w502) - w387*(w137*w474 + 2*w210 + w300*w502 + w387*w508);
    const scalar w524 = w221*w280;
    const scalar w525 = w237*w266;
    const scalar w526 = 2*w525;
    const scalar w527 = w142*w236;
    const scalar w528 = w138*w254;
    const scalar w529 = 2*w528;
    const scalar w530 = w137*w279;
    const scalar w531 = (1.0/(w106*w106*w106));
    const scalar w532 = w280*w531;
    const scalar w533 = -2*w100 - 2*w212 - 2*w216 - 2*w219;
    const scalar w534 = w220*w533;
    const scalar w535 = -w127*w266 + w139*w255;
    const scalar w536 = 4*z;
    const scalar w537 = w281*w31;
    const scalar w538 = w135*w536;
    const scalar w539 = w263 + w264;
    const scalar w540 = (k*k*k);
    const scalar w541 = St*w540;
    const scalar w542 = w26*w541*x;
    const scalar w543 = 2*w147;
    const scalar w544 = 2*St*hA2*k*w178*w25 + St*w27*w465 + St*w474*w77 - w180 - w22*w543 - w39*w542 - w495*w62;
    const scalar w545 = 2*St*hA2*k*w187*w25 + St*w27*w474 + St*w508*w77 - w145*w323 - w178*w543 - w502*w62 - w542*w72;
    const scalar w546 = 2*St*hA2*k*w192*w25 + St*w27*w495 + St*w502*w77 - w146*w89 - w164*w543 - w510*w62 - w542*w56;
    const scalar w547 = w114*w518;
    const scalar w548 = w222*w224;
    const scalar w549 = w141*w237;
    const scalar w550 = 2*w549;
    const scalar w551 = w138*w95;
    const scalar w552 = 2*w551;
    const scalar w553 = w224*w531*w96;
    const scalar w554 = -w127*w141 + w140*w32;
    const scalar w555 = w242*w31;
    const scalar w556 = w132 + w136;
    const scalar w557 = w10*w456;
    const scalar w558 = w10*w461;
    const scalar w559 = 8*w558;
    const scalar w560 = 8*gil*w1*w10*w37*w459*w460 + 8*git*w1*w10*w35*w459*w460 - w36*w559 - w38*w559 - w464;
    const scalar w561 = 4*w558;
    const scalar w562 = 4*w10*w16*w33*w43*w459*w460 - w172*w561 - w467;
    const scalar w563 = w45*w561 + w469 + w561*w68;
    const scalar w564 = git*w310;
    const scalar w565 = w49*w8*y;
    const scalar w566 = 4*gil*theta1*w0*w15*w308*w6*y + gil*w0*w16*w562 + git*w0*w15*w563 - w473 - w561*w67 - w561*w71 - w564*w565;
    const scalar w567 = -w44*w561 - w46*w561 - w475;
    const scalar w568 = 4*w10*w16*w43*w459*w460 - w477*w561 - w478;
    const scalar w569 = (2.0/5.0)*w558;
    const scalar w570 = w482*y;
    const scalar w571 = (2.0/5.0)*gil*theta1*w0*w15*w283*w6*y + (1.0/10.0)*gil*w0*w16*w567 + (1.0/10.0)*git*w0*w15*w568 - w247*w566 - w287*w570 - w485 - w50*w569 - w55*w569 + (1.0/5.0)*w560*z;
    const scalar w572 = w284*w295;
    const scalar w573 = w287*w297;
    const scalar w574 = w48*(-w44*w559 - w46*w559 - w490);
    const scalar w575 = w53*(8*w10*w16*w43*w459*w460 - w477*w559 - w492);
    const scalar w576 = 4*gil*theta1*w0*w15*w283*w6*y + gil*w0*w16*w567 + git*w0*w15*w568 - w287*w565 - w494 - w50*w561 - w55*w561;
    const scalar w577 = w284*w308;
    const scalar w578 = w287*w310;
    const scalar w579 = w48*w562;
    const scalar w580 = w567*w66;
    const scalar w581 = w53*w563;
    const scalar w582 = w568*w70;
    const scalar w583 = 2*w577 + 2*w578 + w579 + w580 + w581 + w582;
    const scalar w584 = w244*w576 - w247*w583 + (1.0/10.0)*w572 + (1.0/10.0)*w573 + (1.0/10.0)*w574 + (1.0/10.0)*w575;
    const scalar w585 = gil*w325;
    const scalar w586 = git*w326;
    const scalar w587 = w308*w585 + w310*w586 + w66*(8*w10*w16*w33*w43*w459*w460 - w172*w559 - w506) + w70*(w45*w559 + w507 + w559*w68);
    const scalar w588 = w244*w566 - w247*w587 + (1.0/5.0)*w577 + (1.0/5.0)*w578 + (1.0/10.0)*w579 + (1.0/10.0)*w580 + (1.0/10.0)*w581 + (1.0/10.0)*w582;
    const scalar w589 = w572 + w573 + w574 + w575;
    const scalar w590 = w511*(w142*(-w142*w560 - w301*w576 - 2*w348 - w388*w566) + w301*(-w142*w576 - w299*w458 - w301*w589 - w388*w583) + w353*w458 + w388*(-w142*w566 - w301*w583 - 2*w352 - w388*w587));
    const scalar w591 = w138*w302;
    const scalar w592 = 2*w591;
    const scalar w593 = w139*w300;
    const scalar w594 = 2*w593;
    const scalar w595 = 2*w122 + 2*w351 + 2*w354 + 2*w357;
    const scalar w596 = (1.0/2.0)*w358;
    const scalar w597 = w595*w596;
    const scalar w598 = -w137*(w137*w560 + w300*w576 + 2*w335 + w387*w566) - w300*(w137*w576 + w299*w513 + w300*w589 + w387*w583) - w340*w513 - w387*(w137*w566 + w300*w583 + 2*w339 + w387*w587);
    const scalar w599 = w301*w358;
    const scalar w600 = w300*w374;
    const scalar w601 = -2*w103 - 2*w338 - 2*w341 - 2*w344;
    const scalar w602 = w345*w601;
    const scalar w603 = 4*y;
    const scalar w604 = w31*w375;
    const scalar w605 = w135*w603;
    const scalar w606 = w144*w76;
    const scalar w607 = St*w606;
    const scalar w608 = w541*w57;
    const scalar w609 = w60*w608;
    const scalar w610 = (6.0/5.0)*St*w324;
    const scalar w611 = 2*w307;
    const scalar w612 = w28*w560 - w289*w611 + w312*w610 + w56*w609 + w566*w78 - w576*w62 + w607*w72;
    const scalar w613 = w607*w89;
    const scalar w614 = w28*w576 - w299*w611 + w321*w610 + w583*w78 - w589*w62 + w609*w85 + w613;
    const scalar w615 = w28*w566 + w315*w608 - w321*w611 + w323*w606 + w327*w610 - w583*w62 + w587*w78;
    const scalar w616 = w31*w361;
    const scalar w617 = w456*w9;
    const scalar w618 = w461*w9;
    const scalar w619 = 8*w618;
    const scalar w620 = 8*gil*w1*w37*w459*w460*w9 + 8*git*w1*w35*w459*w460*w9 - w36*w619 - w38*w619 - w464;
    const scalar w621 = 4*w618;
    const scalar w622 = 4*w16*w33*w43*w459*w460*w9 - w172*w621 - w467;
    const scalar w623 = w45*w621 + w469 + w621*w68;
    const scalar w624 = git*w292;
    const scalar w625 = w49*w8*x;
    const scalar w626 = 4*gil*theta1*w0*w15*w290*w6*x + gil*w0*w16*w622 + git*w0*w15*w623 - w473 - w621*w67 - w621*w71 - w624*w625;
    const scalar w627 = (2.0/5.0)*w441;
    const scalar w628 = -w44*w621 - w46*w621 - w475;
    const scalar w629 = 4*w16*w43*w459*w460*w9 - w477*w621 - w478;
    const scalar w630 = (2.0/5.0)*w618;
    const scalar w631 = w49*x;
    const scalar w632 = w481*w631;
    const scalar w633 = (2.0/5.0)*gil*theta1*w0*w15*w378*w6*x + (1.0/10.0)*gil*w0*w16*w628 + (1.0/10.0)*git*w0*w15*w629 + (1.0/5.0)*w245*w540*w72 - w247*w626 - w294*w627 - w394*w632 - w485 - w50*w630 - w55*w630 + (1.0/5.0)*w620*z;
    const scalar w634 = w246*w540;
    const scalar w635 = w290*w392;
    const scalar w636 = w292*w394;
    const scalar w637 = w48*w622;
    const scalar w638 = w628*w66;
    const scalar w639 = w53*w623;
    const scalar w640 = w629*w70;
    const scalar w641 = gil*w384;
    const scalar w642 = git*w385;
    const scalar w643 = w290*w641 + w292*w642 + w66*(8*w16*w33*w43*w459*w460*w9 - w172*w619 - w506) + w70*(w45*w619 + w507 + w619*w68);
    const scalar w644 = w244*w626 - w247*w643 - w386*w627 + w634*w92 + (1.0/5.0)*w635 + (1.0/5.0)*w636 + (1.0/10.0)*w637 + (1.0/10.0)*w638 + (1.0/10.0)*w639 + (1.0/10.0)*w640;
    const scalar w645 = w392*w399;
    const scalar w646 = w394*w401;
    const scalar w647 = w48*(-w44*w619 - w46*w619 - w490);
    const scalar w648 = w53*(8*w16*w43*w459*w460*w9 - w477*w619 - w492);
    const scalar w649 = 4*gil*theta1*w0*w15*w378*w6*x + gil*w0*w16*w628 + git*w0*w15*w629 - w394*w625 - w494 - w50*w621 - w55*w621;
    const scalar w650 = 2*w635 + 2*w636 + w637 + w638 + w639 + w640;
    const scalar w651 = w244*w649 - w247*w650 - w383*w627 + w634*w89 + (1.0/10.0)*w645 + (1.0/10.0)*w646 + (1.0/10.0)*w647 + (1.0/10.0)*w648;
    const scalar w652 = w645 + w646 + w647 + w648;
    const scalar w653 = w511*(w142*(-w142*w620 - w301*w649 - w388*w626 - 2*w427) + w301*(-w142*w649 - w301*w652 - w388*w650 - 2*w430) + w388*(-w142*w626 - w301*w650 - w386*w458 - w388*w643) + w431*w458);
    const scalar w654 = w138*w389;
    const scalar w655 = 2*w654;
    const scalar w656 = w139*w387;
    const scalar w657 = 2*w656;
    const scalar w658 = (1.0/2.0)*w436;
    const scalar w659 = 2*w124 + 2*w429 + 2*w432 + 2*w435;
    const scalar w660 = w519*w659;
    const scalar w661 = -w137*(w137*w620 + w300*w649 + w387*w626 + 2*w413) - w300*(w137*w649 + w300*w652 + w387*w650 + 2*w416) - w387*(w137*w626 + w300*w650 + w386*w513 + w387*w643) - w417*w513;
    const scalar w662 = w388*w436;
    const scalar w663 = w387*w528;
    const scalar w664 = 2*w423;
    const scalar w665 = -2*w105 - 2*w415 - 2*w418 - 2*w421;
    const scalar w666 = w422*w665;
    const scalar w667 = 4*x;
    const scalar w668 = w31*w454;
    const scalar w669 = w135*w667;
    const scalar w670 = 2*w148;
    const scalar w671 = (3.0/5.0)*w61;
    const scalar w672 = (6.0/5.0)*w306;
    const scalar w673 = w28*w620 - w294*w672 + w390*w670 + w396*w610 + w56*w671 - w62*w649 + w626*w78;
    const scalar w674 = w28*w649 - w383*w672 + w396*w670 + w403*w610 - w62*w652 + w650*w78 + w671*w85;
    const scalar w675 = w671*w89;
    const scalar w676 = w148*w294;
    const scalar w677 = w28*w626 + w383*w610 - w386*w672 - w62*w650 + w643*w78 + w675 + 2*w676;
    const scalar w678 = w547*w658;
    const scalar w679 = w387*w551;
    const scalar w680 = w31*w440;
    const scalar w681 = y*z;
    const scalar w682 = w461*w681;
    const scalar w683 = 8*w682;
    const scalar w684 = w66*(-w172*w683 + w64*w683) + w70*(w45*w683 + w68*w683);
    const scalar w685 = w169*w585 + w175*w586 + w684;
    const scalar w686 = w461*z;
    const scalar w687 = w603*w686;
    const scalar w688 = -w172*w687 + w64*w687;
    const scalar w689 = w45*w687 + w68*w687;
    const scalar w690 = w153*y;
    const scalar w691 = w153*z;
    const scalar w692 = 2*gil*theta1*w0*w15*w169*w6*y + 2*gil*theta1*w0*w15*w308*w6*z + gil*w0*w16*w688 + git*w0*w15*w689 - w471*w690 - w564*w691 - w67*w687 - w687*w71;
    const scalar w693 = -w44*w687 - w46*w687;
    const scalar w694 = w66*w693;
    const scalar w695 = w48*w688;
    const scalar w696 = w157*w308;
    const scalar w697 = w169*w284;
    const scalar w698 = -w477*w687 + w51*w687;
    const scalar w699 = w698*w70;
    const scalar w700 = w53*w689;
    const scalar w701 = w162*w310;
    const scalar w702 = w175*w287;
    const scalar w703 = w481*w54;
    const scalar w704 = w244*w692 + (1.0/5.0)*w309 + (1.0/5.0)*w311 - w570*w70 + w66*w703*y + (1.0/10.0)*w694 + (1.0/10.0)*w695 + (1.0/10.0)*w696 + (1.0/10.0)*w697 + (1.0/10.0)*w699 + (1.0/10.0)*w700 + (1.0/10.0)*w701 + (1.0/10.0)*w702;
    const scalar w705 = -w247*w685 + w704;
    const scalar w706 = w157*w295;
    const scalar w707 = w162*w297;
    const scalar w708 = 2*gil*theta1*w0*w15*w156*w6*y + 2*gil*theta1*w0*w15*w283*w6*z + gil*w0*w16*w693 + git*w0*w15*w698 - w162*w690 - w287*w691 - w50*w687 - w55*w687;
    const scalar w709 = w694 + w695 + w696 + w697 + w699 + w700 + w701 + w702;
    const scalar w710 = w48*(-w44*w683 - w46*w683);
    const scalar w711 = w53*(-w477*w683 + w51*w683);
    const scalar w712 = w481*y;
    const scalar w713 = w244*w708 - w247*w709 + w267*w712 - w269*w712 + (1.0/5.0)*w285 + (1.0/5.0)*w288 + (1.0/10.0)*w710 + (1.0/10.0)*w711;
    const scalar w714 = (1.0/10.0)*w706 + (1.0/10.0)*w707 + w713;
    const scalar w715 = (4.0/5.0)*w7;
    const scalar w716 = (2.0/5.0)*w682;
    const scalar w717 = w362*y;
    const scalar w718 = 8*gil*w1*w37*w459*w460*y*z + 8*git*w1*w35*w459*w460*y*z - w36*w683 - w38*w683;
    const scalar w719 = (1.0/5.0)*gil*theta1*w0*w15*w156*w6*y + (1.0/5.0)*gil*theta1*w0*w15*w283*w6*z + (4.0/5.0)*gil*theta1*w1*w15*w16*w6*y + (1.0/10.0)*gil*w0*w16*w693 + (1.0/10.0)*git*w0*w15*w698 - w162*w49*w717 - w247*w692 - w268*w287*w49 - w483*w715*y - w50*w716 - w55*w716 + (1.0/5.0)*w718*z;
    const scalar w720 = w142*w719 + w276 + w371;
    const scalar w721 = w142*w692 + w230 + w301*w709 + w349;
    const scalar w722 = w710 + w711;
    const scalar w723 = w706 + w707 + w722;
    const scalar w724 = w142*w708 + w229 + w348 + w388*w709;
    const scalar w725 = w142*(-w142*w718 - w225 - w301*w708 - w347 - w388*w692) + w231 + w350;
    const scalar w726 = w511*(w301*(-w301*w723 - w724) + w388*(-w388*w685 - w721) + w725);
    const scalar w727 = w137*w719 + w272 + w366;
    const scalar w728 = w520*w596;
    const scalar w729 = w137*w692 + w214 + w300*w709 + w336;
    const scalar w730 = w137*w708 + w213 + w335 + w387*w709;
    const scalar w731 = w137*(w137*w718 + w209 + w300*w708 + w334 + w387*w692) + w215 + w337;
    const scalar w732 = -w300*(w300*w723 + w730) - w387*(w387*w685 + w729) - w731;
    const scalar w733 = w345*w533;
    const scalar w734 = w279*w300;
    const scalar w735 = w137*w374;
    const scalar w736 = (1.0/2.0)*w346;
    const scalar w737 = (1.0/2.0)*w223;
    const scalar w738 = w236*w301;
    const scalar w739 = w142*w358;
    const scalar w740 = w241*w373 + w275*w593 + w275*w736 + w278*w360 - w278*w591 + w370*w516 + w370*w737 - w373*w514 - w525*w738 - w525*w739 + w528*w734 + w528*w735;
    const scalar w741 = 2*y;
    const scalar w742 = w135*w741;
    const scalar w743 = 2*z;
    const scalar w744 = w135*w743;
    const scalar w745 = w456*w681;
    const scalar w746 = w261*w745 + w281*w742 + w375*w744 + w537*w741 + w604*w743;
    const scalar w747 = -w147*w312 + w148*w327 + w187*w305 - w192*w307 + w28*w692 - w62*w709;
    const scalar w748 = w685*w78 + w747;
    const scalar w749 = w192*w305;
    const scalar w750 = -w147*w289 + w148*w321 - w202*w307 + w28*w708 + w709*w78 + w749;
    const scalar w751 = -w62*w723 + w750;
    const scalar w752 = w148*w312;
    const scalar w753 = -w147*w303 - w164*w307 + w178*w305 + w28*w718 - w62*w708 + w692*w78 + w752;
    const scalar w754 = w142*w753 + w207 + w331;
    const scalar w755 = w137*w753 + w204 + w314;
    const scalar w756 = w205*w593 + w205*w736 + w208*w360 - w208*w591 + w241*w333 + w330*w516 + w330*w737 - w333*w514 - w549*w738 - w549*w739 + w551*w734 + w551*w735;
    const scalar w757 = w130*w745 + w242*w742 + w361*w744 + w555*w741 + w616*w743;
    const scalar w758 = w308*w504 + w310*w505 + w684;
    const scalar w759 = -w247*w758 + w704;
    const scalar w760 = w196*w284;
    const scalar w761 = w200*w287;
    const scalar w762 = w713 + (1.0/10.0)*w760 + (1.0/10.0)*w761;
    const scalar w763 = w722 + w760 + w761;
    const scalar w764 = w511*(w301*(-w301*w763 - w724) + w388*(-w388*w758 - w721) + w725);
    const scalar w765 = w521*w595;
    const scalar w766 = -w300*(w300*w763 + w730) - w387*(w387*w758 + w729) - w731;
    const scalar w767 = w220*w601;
    const scalar w768 = w747 + w758*w78;
    const scalar w769 = -w62*w763 + w750;
    const scalar w770 = 8*x;
    const scalar w771 = w686*w770;
    const scalar w772 = w66*(-w172*w771 + w64*w771) + w70*(w45*w771 + w68*w771);
    const scalar w773 = w169*w641 + w175*w642 + w772;
    const scalar w774 = w667*w686;
    const scalar w775 = -w172*w774 + w64*w774;
    const scalar w776 = w45*w774 + w68*w774;
    const scalar w777 = w153*x;
    const scalar w778 = 2*gil*theta1*w0*w15*w169*w6*x + 2*gil*theta1*w0*w15*w290*w6*z + gil*w0*w16*w775 + git*w0*w15*w776 - w471*w777 - w624*w691 - w67*w774 - w71*w774;
    const scalar w779 = -w44*w774 - w46*w774;
    const scalar w780 = w66*w779;
    const scalar w781 = w48*w775;
    const scalar w782 = w157*w290;
    const scalar w783 = w169*w392;
    const scalar w784 = -w477*w774 + w51*w774;
    const scalar w785 = w70*w784;
    const scalar w786 = w53*w776;
    const scalar w787 = w162*w292;
    const scalar w788 = w175*w394;
    const scalar w789 = -w187*w442 + w244*w778 + (1.0/5.0)*w291 + (1.0/5.0)*w293 - w632*w70 + w66*w703*x + (1.0/10.0)*w780 + (1.0/10.0)*w781 + (1.0/10.0)*w782 + (1.0/10.0)*w783 + (1.0/10.0)*w785 + (1.0/10.0)*w786 + (1.0/10.0)*w787 + (1.0/10.0)*w788;
    const scalar w790 = -w247*w773 + w789;
    const scalar w791 = w157*w399;
    const scalar w792 = w162*w401;
    const scalar w793 = 2*gil*theta1*w0*w15*w156*w6*x + 2*gil*theta1*w0*w15*w378*w6*z + gil*w0*w16*w779 + git*w0*w15*w784 - w162*w777 - w394*w691 - w50*w774 - w55*w774;
    const scalar w794 = w780 + w781 + w782 + w783 + w785 + w786 + w787 + w788;
    const scalar w795 = w48*(-w44*w771 - w46*w771);
    const scalar w796 = w53*(-w477*w771 + w51*w771);
    const scalar w797 = w481*x;
    const scalar w798 = -w192*w442 + w244*w793 - w247*w794 + w267*w797 - w269*w797 + (1.0/5.0)*w393 + (1.0/5.0)*w395 + (1.0/10.0)*w795 + (1.0/10.0)*w796;
    const scalar w799 = (1.0/10.0)*w791 + (1.0/10.0)*w792 + w798;
    const scalar w800 = (2.0/5.0)*x;
    const scalar w801 = w686*w800;
    const scalar w802 = w362*w631;
    const scalar w803 = w394*w49;
    const scalar w804 = 8*gil*w1*w37*w459*w460*x*z + 8*git*w1*w35*w459*w460*x*z - w36*w771 - w38*w771;
    const scalar w805 = (1.0/5.0)*gil*theta1*w0*w15*w156*w6*x + (1.0/5.0)*gil*theta1*w0*w15*w378*w6*z + (4.0/5.0)*gil*theta1*w1*w15*w16*w6*x + (1.0/10.0)*gil*w0*w16*w779 + (1.0/10.0)*git*w0*w15*w784 - w162*w802 - w178*w442 - w247*w778 - w268*w803 - w483*w715*x - w50*w801 - w55*w801 + (1.0/5.0)*w804*z;
    const scalar w806 = w142*w805 + w277 + w450;
    const scalar w807 = w142*w778 + w233 + w301*w794 + w427;
    const scalar w808 = w795 + w796;
    const scalar w809 = w791 + w792 + w808;
    const scalar w810 = w142*w793 + w230 + w388*w794 + w426;
    const scalar w811 = w142*(-w142*w804 - w226 - w301*w793 - w388*w778 - w425) + w234 + w428;
    const scalar w812 = w511*(w301*(-w301*w809 - w810) + w388*(-w388*w773 - w807) + w811);
    const scalar w813 = w137*w805 + w274 + w444;
    const scalar w814 = w519*w658;
    const scalar w815 = w137*w778 + w217 + w300*w794 + w413;
    const scalar w816 = w137*w793 + w214 + w387*w794 + w412;
    const scalar w817 = w137*(w137*w804 + w210 + w300*w793 + w387*w778 + w411) + w218 + w414;
    const scalar w818 = -w300*(w300*w809 + w816) - w387*(w387*w773 + w815) - w817;
    const scalar w819 = w422*w533;
    const scalar w820 = w238*w453;
    const scalar w821 = w137*w423;
    const scalar w822 = (1.0/2.0)*w424;
    const scalar w823 = (1.0/2.0)*w449;
    const scalar w824 = w236*w388;
    const scalar w825 = w142*w436;
    const scalar w826 = w275*w656 + w275*w822 + w278*w439 - w278*w654 + w279*w663 + w279*w823 + w448*w516 - w452*w514 + w521*w820 - w525*w824 - w525*w825 + w528*w821;
    const scalar w827 = 2*x;
    const scalar w828 = w135*w827;
    const scalar w829 = w456*x;
    const scalar w830 = w829*z;
    const scalar w831 = w261*w830 + w281*w828 + w454*w744 + w537*w827 + w668*w743;
    const scalar w832 = -w147*w294 + w148*w178 + w148*w386 - w180 - w187*w391 + w28*w778 - w62*w794 + w749;
    const scalar w833 = w773*w78 + w832;
    const scalar w834 = -w146*w56 - w147*w396 + w148*w164 + w148*w383 - w192*w391 + w202*w305 + w28*w793 + w78*w794;
    const scalar w835 = -w62*w809 + w834;
    const scalar w836 = -w146*w39 - w147*w390 + w148*w22 + w164*w305 - w178*w391 + w28*w804 - w62*w793 + w676 + w778*w78;
    const scalar w837 = w142*w836 + w206 + w408;
    const scalar w838 = w137*w836 + w194 + w398;
    const scalar w839 = (1.0/2.0)*w407;
    const scalar w840 = w205*w656 + w205*w822 + w208*w439 - w208*w654 + w223*w839 + w241*w410 + w279*w679 + w407*w516 - w410*w514 - w549*w824 - w549*w825 + w551*w821;
    const scalar w841 = w130*w830 + w148 + w242*w828 + w440*w744 + w555*w827 + w680*w743;
    const scalar w842 = w290*w504 + w292*w505 + w772;
    const scalar w843 = -w247*w842 + w789;
    const scalar w844 = w196*w392;
    const scalar w845 = w200*w394;
    const scalar w846 = w798 + (1.0/10.0)*w844 + (1.0/10.0)*w845;
    const scalar w847 = w808 + w844 + w845;
    const scalar w848 = w511*(w301*(-w301*w847 - w810) + w388*(-w388*w842 - w807) + w811);
    const scalar w849 = -w300*(w300*w847 + w816) - w387*(w387*w842 + w815) - w817;
    const scalar w850 = w220*w665;
    const scalar w851 = w78*w842 + w832;
    const scalar w852 = -w62*w847 + w834;
    const scalar w853 = w547*w659;
    const scalar w854 = w461*y;
    const scalar w855 = w770*w854;
    const scalar w856 = w66*(-w172*w855 + w64*w855) + w70*(w45*w855 + w68*w855);
    const scalar w857 = w308*w641 + w310*w642 + w856;
    const scalar w858 = w667*w854;
    const scalar w859 = -w172*w858 + w64*w858;
    const scalar w860 = w45*w858 + w68*w858;
    const scalar w861 = 2*gil*theta1*w0*w15*w290*w6*y + 2*gil*theta1*w0*w15*w308*w6*x + gil*w0*w16*w859 + git*w0*w15*w860 - w564*w777 - w624*w690 - w67*w858 - w71*w858;
    const scalar w862 = -w44*w858 - w46*w858;
    const scalar w863 = w66*w862;
    const scalar w864 = w48*w859;
    const scalar w865 = w284*w290;
    const scalar w866 = w308*w392;
    const scalar w867 = -w477*w858 + w51*w858;
    const scalar w868 = w70*w867;
    const scalar w869 = w53*w860;
    const scalar w870 = w287*w292;
    const scalar w871 = w310*w394;
    const scalar w872 = w244*w861 - w327*w442 + (1.0/10.0)*w863 + (1.0/10.0)*w864 + (1.0/10.0)*w865 + (1.0/10.0)*w866 + (1.0/10.0)*w868 + (1.0/10.0)*w869 + (1.0/10.0)*w870 + (1.0/10.0)*w871;
    const scalar w873 = -w247*w857 + w872;
    const scalar w874 = w284*w399;
    const scalar w875 = w287*w401;
    const scalar w876 = 2*gil*theta1*w0*w15*w283*w6*x + 2*gil*theta1*w0*w15*w378*w6*y + gil*w0*w16*w862 + git*w0*w15*w867 - w287*w777 - w394*w690 - w50*w858 - w55*w858;
    const scalar w877 = w863 + w864 + w865 + w866 + w868 + w869 + w870 + w871;
    const scalar w878 = w48*(-w44*w855 - w46*w855);
    const scalar w879 = w53*(-w477*w855 + w51*w855);
    const scalar w880 = w244*w876 - w247*w877 - w321*w442 + (1.0/10.0)*w878 + (1.0/10.0)*w879;
    const scalar w881 = (1.0/10.0)*w874 + (1.0/10.0)*w875 + w880;
    const scalar w882 = w800*w854;
    const scalar w883 = 8*gil*w1*w37*w459*w460*x*y + 8*git*w1*w35*w459*w460*x*y - w36*w855 - w38*w855;
    const scalar w884 = (1.0/5.0)*gil*theta1*w0*w15*w283*w6*x + (1.0/5.0)*gil*theta1*w0*w15*w378*w6*y + (1.0/10.0)*gil*w0*w16*w862 + (1.0/10.0)*git*w0*w15*w867 - w247*w861 - w287*w802 - w312*w442 - w50*w882 - w55*w882 - w717*w803 + (1.0/5.0)*w883*z;
    const scalar w885 = w142*w884 + w372 + w451;
    const scalar w886 = w142*w861 + w301*w877 + w355 + w430;
    const scalar w887 = w878 + w879;
    const scalar w888 = w874 + w875 + w887;
    const scalar w889 = w142*w876 + w352 + w388*w877 + w433;
    const scalar w890 = w142*(-w142*w883 - w301*w876 - w349 - w388*w861 - w426) + w356 + w434;
    const scalar w891 = w511*(w301*(-w301*w888 - w889) + w388*(-w388*w857 - w886) + w890);
    const scalar w892 = w137*w884 + w369 + w447;
    const scalar w893 = w137*w861 + w300*w877 + w342 + w416;
    const scalar w894 = w137*w876 + w339 + w387*w877 + w419;
    const scalar w895 = w137*(w137*w883 + w300*w876 + w336 + w387*w861 + w412) + w343 + w420;
    const scalar w896 = -w300*(w300*w888 + w894) - w387*(w387*w857 + w893) - w895;
    const scalar w897 = w422*w601;
    const scalar w898 = w300*w423;
    const scalar w899 = w358*w388;
    const scalar w900 = w301*w436;
    const scalar w901 = w370*w656 + w370*w822 + w373*w439 - w373*w654 + w374*w663 + w374*w823 + w448*w593 - w452*w591 - w525*w899 - w525*w900 + w528*w898 + w596*w820;
    const scalar w902 = w829*y;
    const scalar w903 = w261*w902 + w375*w828 + w454*w742 + w604*w827 + w668*w741;
    const scalar w904 = w28*w861 + w305*w321 + w305*w386 - w307*w383 - w327*w391 + w613 - w62*w877 + w671*w92 + w752;
    const scalar w905 = w78*w857 + w904;
    const scalar w906 = w148*w289 + w28*w876 + w299*w305 + w305*w383 - w307*w403 - w321*w391 + w607*w85 + w675 + w78*w877;
    const scalar w907 = -w62*w888 + w906;
    const scalar w908 = w148*w303 + w28*w883 + w289*w305 + w294*w305 - w307*w396 - w312*w391 + w56*w607 - w62*w876 + w671*w72 + w78*w861;
    const scalar w909 = w142*w908 + w332 + w409;
    const scalar w910 = w137*w908 + w329 + w405;
    const scalar w911 = w330*w656 + w330*w822 + w333*w439 - w333*w654 + w346*w839 + w360*w410 + w374*w679 + w407*w593 - w410*w591 - w549*w899 - w549*w900 + w551*w898;
    const scalar w912 = w130*w902 + w305 + w361*w828 + w440*w742 + w616*w827 + w680*w741;
    const scalar w913 = w290*w585 + w292*w586 + w856;
    const scalar w914 = -w247*w913 + w872;
    const scalar w915 = w295*w392;
    const scalar w916 = w297*w394;
    const scalar w917 = w880 + (1.0/10.0)*w915 + (1.0/10.0)*w916;
    const scalar w918 = w887 + w915 + w916;
    const scalar w919 = w511*(w301*(-w301*w918 - w889) + w388*(-w388*w913 - w886) + w890);
    const scalar w920 = -w300*(w300*w918 + w894) - w387*(w387*w913 + w893) - w895;
    const scalar w921 = w345*w665;
    const scalar w922 = w78*w913 + w904;
    const scalar w923 = -w62*w918 + w906;
    return w164*w376 + w178*w455 + w22*w282 + w282*w289 + w282*w294 + w299*w376 + w321*w455 + w376*w383 + w386*w455 + w39*(Tt*(w243*(w109*(-w137*w486 - w270*w513 - w300*w503 - w387*w509) + w129*(w142*w486 + w270*w458 + w301*w503 + w388*w509) + w223*w275 + w240*w278 + w260*w512 + w275*w517 - w278*w515 + w519*w522 + w523*w524 - w526*w527 + w529*w530 + w532*w534 + w535) + w261*w457 + w281*w538 + w536*w537 + w539 + 1.0/5.0) + w130*w457 - w147 + w242*w538 + w243*(w109*(-w137*w544 - w179*w513 - w300*w546 - w387*w545) + w114*w512 + w129*(w142*w544 + w179*w458 + w301*w546 + w388*w545) + w205*w223 + w205*w517 + w208*w240 - w208*w515 + w522*w547 + w523*w548 - w527*w550 + w530*w552 + w534*w553 + w554) + w536*w555 + w556) + w56*(Tt*(w243*(w109*(-w300*w714 - w387*w705 - w727) + w129*(w301*w714 + w388*w705 + w720) + w260*w726 + w519*w728 + w524*w732 + w532*w733 + w740) + w746) + w243*(w109*(-w300*w751 - w387*w748 - w755) + w114*w726 + w129*(w301*w751 + w388*w748 + w754) + w547*w728 + w548*w732 + w553*w733 + w756) + w757) + w56*(Tt*(w243*(w109*(-w300*w762 - w387*w759 - w727) + w129*(w301*w762 + w388*w759 + w720) + w260*w764 + w519*w765 + w524*w766 + w532*w767 + w740) + w746) + w243*(w109*(-w300*w769 - w387*w768 - w755) + w114*w764 + w129*(w301*w769 + w388*w768 + w754) + w547*w765 + w548*w766 + w553*w767 + w756) + w757) + w72*(Tt*(w243*(w109*(-w300*w799 - w387*w790 - w813) + w129*(w301*w799 + w388*w790 + w806) + w260*w812 + w520*w814 + w524*w818 + w532*w819 + w826) + w831) + w243*(w109*(-w300*w835 - w387*w833 - w838) + w114*w812 + w129*(w301*w835 + w388*w833 + w837) + w520*w678 + w548*w818 + w553*w819 + w840) + w841) + w72*(Tt*(w243*(w109*(-w300*w846 - w387*w843 - w813) + w129*(w301*w846 + w388*w843 + w806) + w260*w848 + w521*w660 + w524*w849 + w532*w850 + w826) + w831) + w243*(w109*(-w300*w852 - w387*w851 - w838) + w114*w848 + w129*(w301*w852 + w388*w851 + w837) + w521*w853 + w548*w849 + w553*w850 + w840) + w841) + w85*(Tt*(w243*(w109*(-w137*w571 - w300*w584 - w367*w513 - w387*w588) + w129*(w142*w571 + w301*w584 + w367*w458 + w388*w588) + w260*w590 + w346*w370 + w359*w373 + w370*w594 - w373*w592 + w519*w597 + w524*w598 - w526*w599 + w529*w600 + w532*w602 + w535) + w261*w557 + w375*w605 + w539 + w603*w604) + w130*w557 + w243*(w109*(-w137*w612 - w300*w614 - w322*w513 - w387*w615) + w114*w590 + w129*(w142*w612 + w301*w614 + w322*w458 + w388*w615) + w330*w346 + w330*w594 + w333*w359 - w333*w592 + w547*w597 + w548*w598 - w550*w599 + w552*w600 + w553*w602 + w554) - w307 + w361*w605 + w556 + w603*w616) + w89*(Tt*(w243*(w109*(-w300*w881 - w387*w873 - w892) + w129*(w301*w881 + w388*w873 + w885) + w260*w891 + w524*w896 + w532*w897 + w595*w814 + w901) + w903) + w243*(w109*(-w300*w907 - w387*w905 - w910) + w114*w891 + w129*(w301*w907 + w388*w905 + w909) + w548*w896 + w553*w897 + w595*w678 + w911) + w912) + w89*(Tt*(w243*(w109*(-w300*w917 - w387*w914 - w892) + w129*(w301*w917 + w388*w914 + w885) + w260*w919 + w524*w920 + w532*w921 + w596*w660 + w901) + w903) + w243*(w109*(-w300*w923 - w387*w922 - w910) + w114*w919 + w129*(w301*w923 + w388*w922 + w909) + w548*w920 + w553*w921 + w596*w853 + w911) + w912) + w92*(Tt*(w243*(w109*(-w137*w633 - w300*w651 - w387*w644 - w445*w513) + w129*(w142*w633 + w301*w651 + w388*w644 + w445*w458) + w260*w653 + w423*w449 + w437*w453 + w448*w657 - w452*w655 + w524*w661 - w526*w662 + w532*w666 + w535 + w658*w660 + w663*w664) + w261*w617 - w442 + w454*w669 + w539 + w667*w668) + w130*w617 + w243*(w109*(-w137*w673 - w300*w674 - w387*w677 - w406*w513) + w114*w653 + w129*(w142*w673 + w301*w674 + w388*w677 + w406*w458) + w407*w424 + w407*w657 + w410*w438 - w410*w655 + w548*w661 - w550*w662 + w553*w666 + w554 + w659*w678 + w664*w679) - w391 + w440*w669 + w556 + w667*w680);
}

inline scalar phiESource(const point& pt, const scalar t, const label zone, const Params& p)
{
    const scalar r1 = p.r1;
    const scalar r2 = p.r2;
    const scalar gamma1 = p.gamma1;
    const scalar gamma2 = p.gamma2;
    const scalar k = p.k;
    const scalar gil = p.gil;
    const scalar git = p.git;
    const scalar gin = p.gin;
    const scalar gel = p.gel;
    const scalar get = p.get;
    const scalar gen = p.gen;
    const scalar theta0 = p.theta0;
    const scalar theta1 = p.theta1;
    const scalar sigma0In = p.sigma0In;
    const scalar sigma0Out = p.sigma0Out;
    const scalar sigmaEps = p.sigmaEps;
    const scalar fA0 = p.fA[0];
    const scalar fA1 = p.fA[1];
    const scalar fA2 = p.fA[2];
    const scalar fA3 = p.fA[3];
    const scalar hA0 = p.hA[0];
    const scalar hA1 = p.hA[1];
    const scalar hA2 = p.hA[2];
    const scalar cmAmplitude = p.cmAmplitude;
    const scalar cmOmega = p.cmOmega;
    const scalar x = pt.x();
    const scalar y = pt.y();
    const scalar z = pt.z();
    const scalar St = sqrt(1 + t);
    const scalar Tt = cmAmplitude*sin(cmOmega*t);
    checkZone(zone, true);

    if (zone == bloodPool)
    {
        const scalar w0 = k*x;
        const scalar w1 = sin(w0);
        const scalar w2 = fA2*w1;
        const scalar w3 = k*y;
        const scalar w4 = w3 + 3.0/10.0;
        const scalar w5 = fA3*exp((2.0/5.0)*z);
        const scalar w6 = (1.0/(r1));
        const scalar w7 = w6*y;
        const scalar w8 = (z*z);
        const scalar w9 = (1.0/5.0)*w6;
        const scalar w10 = (1.0/5.0)*k*w1*w6*x - 1.0/10.0*w7 - w8*w9;
        const scalar w11 = cos(w3);
        const scalar w12 = k*z;
        const scalar w13 = sin(w12);
        const scalar w14 = (1.0/2.0)*w13 + 1;
        const scalar w15 = w0 + 1.0/5.0;
        const scalar w16 = sin(w15);
        const scalar w17 = sigmaEps*w14*w16;
        const scalar w18 = w11*w17;
        const scalar w19 = w18 + 1;
        const scalar w20 = sigma0In*w19;
        const scalar w21 = cos(gamma1);
        const scalar w22 = cos(gamma2);
        const scalar w23 = sin(gamma2);
        const scalar w24 = w22*w23;
        const scalar w25 = w21*w24;
        const scalar w26 = gen*w25;
        const scalar w27 = sin(gamma1);
        const scalar w28 = cos(theta0);
        const scalar w29 = w27*w28;
        const scalar w30 = sin(theta0);
        const scalar w31 = w21*w30;
        const scalar w32 = w22*w31 + w29;
        const scalar w33 = w23*w30;
        const scalar w34 = gel*w33;
        const scalar w35 = w32*w34;
        const scalar w36 = w27*w30;
        const scalar w37 = w21*w22*w28 - w36;
        const scalar w38 = w23*w28;
        const scalar w39 = get*w38;
        const scalar w40 = w37*w39;
        const scalar w41 = (w22*w22);
        const scalar w42 = gen*w41;
        const scalar w43 = (w23*w23);
        const scalar w44 = (w30*w30)*w43;
        const scalar w45 = gel*w44;
        const scalar w46 = (w28*w28)*w43;
        const scalar w47 = get*w46;
        const scalar w48 = w42 + w45 + w47;
        const scalar w49 = (1.0/5.0)*z;
        const scalar w50 = gen*w27;
        const scalar w51 = w24*w50;
        const scalar w52 = w21*w28 - w22*w36;
        const scalar w53 = w34*w52;
        const scalar w54 = -w22*w29 - w31;
        const scalar w55 = w39*w54;
        const scalar w56 = w51 + w53 + w55;
        const scalar w57 = k*w1;
        const scalar w58 = (1.0/5.0)*w57;
        const scalar w59 = w6*(-1.0/10.0*w26 + (1.0/10.0)*w35 + (1.0/10.0)*w40 + w48*w49 - w56*w58);
        const scalar w60 = (w21*w21)*w43;
        const scalar w61 = gen*w60;
        const scalar w62 = (w32*w32);
        const scalar w63 = gel*w62;
        const scalar w64 = (w37*w37);
        const scalar w65 = get*w64;
        const scalar w66 = -w26 + w35 + w40;
        const scalar w67 = w21*w43;
        const scalar w68 = w50*w67;
        const scalar w69 = w32*w52;
        const scalar w70 = gel*w69;
        const scalar w71 = w37*w54;
        const scalar w72 = get*w71;
        const scalar w73 = -w68 + w70 + w72;
        const scalar w74 = w6*(w49*w66 - w58*w73 + (1.0/10.0)*w61 + (1.0/10.0)*w63 + (1.0/10.0)*w65);
        const scalar w75 = w49*w56;
        const scalar w76 = (w27*w27)*w43;
        const scalar w77 = (w52*w52);
        const scalar w78 = (w54*w54);
        const scalar w79 = gel*w77 + gen*w76 + get*w78;
        const scalar w80 = w6*(-w58*w79 - 1.0/10.0*w68 + (1.0/10.0)*w70 + (1.0/10.0)*w72 + w75);
        const scalar w81 = (1.0/(((1.0/(git)) + (1.0/(gin)) + (1.0/(gil)))));
        const scalar w82 = (r1*r1);
        const scalar w83 = -w8 + w82 - (x*x) - (y*y);
        const scalar w84 = gin*w41;
        const scalar w85 = gil*w44;
        const scalar w86 = git*w46;
        const scalar w87 = z*(w84 + w85 + w86);
        const scalar w88 = gin*w25;
        const scalar w89 = gil*w33;
        const scalar w90 = w32*w89;
        const scalar w91 = git*w38;
        const scalar w92 = w37*w91;
        const scalar w93 = -w88 + w90 + w92;
        const scalar w94 = w93*y;
        const scalar w95 = gin*w27;
        const scalar w96 = w24*w95;
        const scalar w97 = w52*w89;
        const scalar w98 = w54*w91;
        const scalar w99 = w96 + w97 + w98;
        const scalar w100 = w99*x;
        const scalar w101 = gin*w60;
        const scalar w102 = gil*w62;
        const scalar w103 = git*w64;
        const scalar w104 = y*(w101 + w102 + w103);
        const scalar w105 = w93*z;
        const scalar w106 = w67*w95;
        const scalar w107 = gil*w69;
        const scalar w108 = git*w71;
        const scalar w109 = -w106 + w107 + w108;
        const scalar w110 = w109*x;
        const scalar w111 = gin*w76;
        const scalar w112 = gil*w77;
        const scalar w113 = git*w78;
        const scalar w114 = w111 + w112 + w113;
        const scalar w115 = w114*x;
        const scalar w116 = w99*z;
        const scalar w117 = w109*y;
        const scalar w118 = w81*w83 + x*(w115 + w116 + w117) + y*(w104 + w105 + w110) + z*(w100 + w87 + w94);
        const scalar w119 = (1.0/(w118));
        const scalar w120 = w6*(-w58*w99 + (1.0/5.0)*w87 - 1.0/10.0*w88 + (1.0/10.0)*w90 + (1.0/10.0)*w92);
        const scalar w121 = (1.0/10.0)*w101 + (1.0/10.0)*w102 + (1.0/10.0)*w103 + (1.0/5.0)*w105 - w109*w58;
        const scalar w122 = w6*(-1.0/10.0*w106 + (1.0/10.0)*w107 + (1.0/10.0)*w108 - w114*w58 + (1.0/5.0)*w116);
        const scalar w123 = -w120*z - w121*w7 - w122*x;
        const scalar w124 = w119*w123;
        const scalar w125 = w6*z;
        const scalar w126 = w6*x;
        const scalar w127 = w6*(-w125*w48 - w126*w56 - w66*w7);
        const scalar w128 = w61 + w63 + w65;
        const scalar w129 = w6*(-w125*w66 - w126*w73 - w128*w7);
        const scalar w130 = w6*(-w125*w56 - w126*w79 - w7*w73);
        const scalar w131 = -w127*z - w129*y - w130*x;
        const scalar w132 = w131*w82;
        const scalar w133 = -w10*w20 - w124*w132 - w59*z - w74*y - w80*x;
        const scalar w134 = (1.0/(sigma0In));
        const scalar w135 = (1.0/(w19));
        const scalar w136 = w134*w135;
        const scalar w137 = w133*w136;
        const scalar w138 = (1.0/(w19*w19));
        const scalar w139 = sin(w3);
        const scalar w140 = w139*w17;
        const scalar w141 = w138*w140;
        const scalar w142 = w133*w134;
        const scalar w143 = w141*w142;
        const scalar w144 = w6*w83;
        const scalar w145 = (1.0/2.0)*w144;
        const scalar w146 = k*w145;
        const scalar w147 = k*sigma0In*w10;
        const scalar w148 = (1.0/(w82));
        const scalar w149 = w148*y;
        const scalar w150 = w148*z;
        const scalar w151 = w148*x;
        const scalar w152 = w128*w149 - w129 + w150*w66 + w151*w73;
        const scalar w153 = w124*w82;
        const scalar w154 = -2*w104 - 2*w105 - 2*w110 + 2*w81*y;
        const scalar w155 = (1.0/(w118*w118));
        const scalar w156 = w123*w155;
        const scalar w157 = w132*w156;
        const scalar w158 = r1*w119*w121*w131 + w140*w147 - w152*w153 - w154*w157 + (1.0/10.0)*w20*w6 - w74;
        const scalar w159 = w136*w145;
        const scalar w160 = k*sigma0In;
        const scalar w161 = w140*w160;
        const scalar w162 = w5*cos(w4);
        const scalar w163 = cos(w12);
        const scalar w164 = sigmaEps*w11;
        const scalar w165 = w16*w164;
        const scalar w166 = w163*w165;
        const scalar w167 = w134*w144;
        const scalar w168 = k*w167;
        const scalar w169 = w138*w168;
        const scalar w170 = (1.0/5.0)*w42 + (1.0/5.0)*w45 + (1.0/5.0)*w47;
        const scalar w171 = w147*w164;
        const scalar w172 = (1.0/2.0)*w163;
        const scalar w173 = -w127 + w149*w66 + w150*w48 + w151*w56;
        const scalar w174 = (1.0/5.0)*w84 + (1.0/5.0)*w85 + (1.0/5.0)*w86;
        const scalar w175 = -w120 - w125*w174 - w126*((1.0/5.0)*w96 + (1.0/5.0)*w97 + (1.0/5.0)*w98) - w7*(-1.0/5.0*w88 + (1.0/5.0)*w90 + (1.0/5.0)*w92);
        const scalar w176 = w119*w132;
        const scalar w177 = -2*w100 + 2*w81*z - 2*w87 - 2*w94;
        const scalar w178 = (2.0/5.0)*sigma0In*w19*w6*z - w125*w170 - w126*((1.0/5.0)*w51 + (1.0/5.0)*w53 + (1.0/5.0)*w55) - w153*w173 - w157*w177 - w16*w171*w172 - w175*w176 - w59 - w7*(-1.0/5.0*w26 + (1.0/5.0)*w35 + (1.0/5.0)*w40);
        const scalar w179 = cos(w0);
        const scalar w180 = (1.0/2.0)*y + 1;
        const scalar w181 = cos(w15);
        const scalar w182 = w14*w181;
        const scalar w183 = w164*w182;
        const scalar w184 = w138*w142;
        const scalar w185 = w145*w184;
        const scalar w186 = (k*k);
        const scalar w187 = w179*w186;
        const scalar w188 = w187*w9;
        const scalar w189 = w188*x + w57*w9;
        const scalar w190 = -w130 + w149*w73 + w150*w56 + w151*w79;
        const scalar w191 = w115*w188 + w116*w188 + w117*w188 - w122;
        const scalar w192 = -2*w115 - 2*w116 - 2*w117 + 2*w81*x;
        const scalar w193 = -w153*w190 - w157*w192 - w171*w182 - w176*w191 + (1.0/5.0)*w179*w186*w56*w6*z + (1.0/5.0)*w179*w186*w6*w73*y + (1.0/5.0)*w179*w186*w6*w79*x - w189*w20 - w80;
        const scalar w194 = w160*w183;
        const scalar w195 = St*w162;
        const scalar w196 = w137*w6;
        const scalar w197 = -w196;
        const scalar w198 = 2*w6;
        const scalar w199 = (w16*w16);
        const scalar w200 = (sigmaEps*sigmaEps);
        const scalar w201 = (w14*w14);
        const scalar w202 = (1.0/(w19*w19*w19));
        const scalar w203 = w133*w167;
        const scalar w204 = 2*w136;
        const scalar w205 = -sigma0In*sigmaEps*w10*w11*w14*w16*w186;
        const scalar w206 = 2*w124;
        const scalar w207 = w156*w82;
        const scalar w208 = -2*w81;
        const scalar w209 = w123*w132/(w118*w118*w118);
        const scalar w210 = (1.0/4.0)*w203;
        const scalar w211 = (w11*w11);
        const scalar w212 = 2*w119*w82;
        const scalar w213 = 2*w207;
        const scalar w214 = 2*w132*w155;
        const scalar w215 = (k*k*k);
        const scalar w216 = w1*w215*w9;
        const scalar w217 = w216*x;
        return w160*w165*w172*(St*(fA1*y + (2.0/5.0)*w162) + Tt*(-w125*w137 - 1.0/4.0*w133*w166*w169 + (1.0/2.0)*w134*w135*w178*w6*w83 + (1.0/5.0)*z)) - w161*(St*(fA1*z - k*w5*sin(w4) + (1.0/2.0)*w2) + Tt*(-w137*w7 + w143*w146 + w158*w159 + 1.0/10.0)) + w194*(St*(fA0 + fA2*k*w179*w180) + Tt*(-k*w183*w185 - w126*w137 + (1.0/2.0)*w134*w135*w193*w6*w83 - w58)) + w20*(Tt*((w139*w139)*w186*w199*w200*w201*w202*w203 + w141*w158*w168 - w143*w198*w3 - w158*w204*w7 + w159*(2*r1*w119*w121*w152 + 2*r1*w121*w131*w154*w155 - w128*w206 - 2*w152*w154*w207 - w154*w209*(-4*w104 - 4*w105 - 4*w110 + 4*w81*y) - w157*(-2*w101 - 2*w102 - 2*w103 - w208) - w161*w9 - w205) + w18*w185*w186 + w197) - w186*w195) + w20*(Tt*(w12*w166*w184*w6 - w125*w178*w204 + w13*w138*w165*w186*w210 - w134*w138*w146*w166*w178 + w159*((2.0/5.0)*k*sigma0In*sigmaEps*w11*w16*w163*w6*z + 2*r1*w119*w131*w174 + (1.0/2.0)*sigma0In*sigmaEps*w10*w11*w13*w16*w186 + (2.0/5.0)*sigma0In*w19*w6 - w157*(-w208 - 2*w84 - 2*w85 - 2*w86) - w170*w198 - w173*w175*w212 - w173*w177*w213 - w175*w177*w214 - w177*w209*(-4*w100 + 4*w81*z - 4*w87 - 4*w94) - w206*w48) + (w163*w163)*w186*w199*w200*w202*w210*w211 + w197 + 1.0/5.0) + (4.0/25.0)*w195) + w20*(-St*w180*w186*w2 + Tt*(2*k*sigmaEps*w11*w133*w134*w138*w14*w181*w6*x + (1.0/2.0)*sigmaEps*w11*w133*w134*w138*w14*w16*w186*w6*w83 - w126*w193*w204 + w133*w134*(w181*w181)*w186*w200*w201*w202*w211*w6*w83 + (1.0/2.0)*w134*w135*w6*w83*(-w1*w215*w6*w75 - w157*(-2*w111 - 2*w112 - 2*w113 - w208) - w176*((2.0/5.0)*w114*w179*w186*w6 - w115*w216 - w116*w216 - w117*w216) + (2.0/5.0)*w179*w186*w6*w79 - 2*w189*w194 - w190*w191*w212 - w190*w192*w213 - w191*w192*w214 - w192*w209*(-4*w115 - 4*w116 - 4*w117 + 4*w81*x) - w20*((2.0/5.0)*w179*w186*w6 - w217) - w205 - w206*w79 - w216*w73*y - w217*w79) - w169*w183*w193 - 1.0/5.0*w187 - w196));
    }
    if (zone == bath)
    {
        const scalar w0 = k*x;
        const scalar w1 = sin(w0);
        const scalar w2 = fA2*w1;
        const scalar w3 = k*y;
        const scalar w4 = w3 + 3.0/10.0;
        const scalar w5 = (2.0/5.0)*z;
        const scalar w6 = fA3*exp(w5);
        const scalar w7 = (1.0/(r2));
        const scalar w8 = w7*y;
        const scalar w9 = (z*z);
        const scalar w10 = (1.0/5.0)*w7;
        const scalar w11 = w1*w10;
        const scalar w12 = -w0*w11 + w10*w9 + (1.0/10.0)*w8;
        const scalar w13 = cos(w3);
        const scalar w14 = k*z;
        const scalar w15 = sin(w14);
        const scalar w16 = (1.0/2.0)*w15 + 1;
        const scalar w17 = w0 + 1.0/5.0;
        const scalar w18 = sin(w17);
        const scalar w19 = sigmaEps*w16*w18;
        const scalar w20 = w13*w19;
        const scalar w21 = w20 + 1;
        const scalar w22 = sigma0Out*w21;
        const scalar w23 = cos(gamma1);
        const scalar w24 = cos(gamma2);
        const scalar w25 = sin(gamma2);
        const scalar w26 = w24*w25;
        const scalar w27 = w23*w26;
        const scalar w28 = gen*w27;
        const scalar w29 = sin(gamma1);
        const scalar w30 = theta0 + theta1;
        const scalar w31 = cos(w30);
        const scalar w32 = w29*w31;
        const scalar w33 = sin(w30);
        const scalar w34 = w23*w33;
        const scalar w35 = w24*w34 + w32;
        const scalar w36 = w25*w33;
        const scalar w37 = gel*w36;
        const scalar w38 = w35*w37;
        const scalar w39 = w29*w33;
        const scalar w40 = w23*w24*w31 - w39;
        const scalar w41 = w25*w31;
        const scalar w42 = get*w41;
        const scalar w43 = w40*w42;
        const scalar w44 = (w24*w24);
        const scalar w45 = gen*w44;
        const scalar w46 = (w25*w25);
        const scalar w47 = (w33*w33)*w46;
        const scalar w48 = gel*w47;
        const scalar w49 = (w31*w31)*w46;
        const scalar w50 = get*w49;
        const scalar w51 = w45 + w48 + w50;
        const scalar w52 = (1.0/5.0)*z;
        const scalar w53 = gen*w29;
        const scalar w54 = w26*w53;
        const scalar w55 = w23*w31 - w24*w39;
        const scalar w56 = w37*w55;
        const scalar w57 = -w24*w32 - w34;
        const scalar w58 = w42*w57;
        const scalar w59 = w54 + w56 + w58;
        const scalar w60 = k*w1;
        const scalar w61 = (1.0/5.0)*w60;
        const scalar w62 = -1.0/10.0*w28 + (1.0/10.0)*w38 + (1.0/10.0)*w43 + w51*w52 - w59*w61;
        const scalar w63 = (w23*w23)*w46;
        const scalar w64 = gen*w63;
        const scalar w65 = (w35*w35);
        const scalar w66 = gel*w65;
        const scalar w67 = (w40*w40);
        const scalar w68 = get*w67;
        const scalar w69 = -w28 + w38 + w43;
        const scalar w70 = w23*w46;
        const scalar w71 = w53*w70;
        const scalar w72 = w35*w55;
        const scalar w73 = gel*w72;
        const scalar w74 = w40*w57;
        const scalar w75 = get*w74;
        const scalar w76 = -w71 + w73 + w75;
        const scalar w77 = w52*w69 - w61*w76 + (1.0/10.0)*w64 + (1.0/10.0)*w66 + (1.0/10.0)*w68;
        const scalar w78 = (w29*w29)*w46;
        const scalar w79 = (w55*w55);
        const scalar w80 = (w57*w57);
        const scalar w81 = gel*w79 + gen*w78 + get*w80;
        const scalar w82 = w52*w59;
        const scalar w83 = -w61*w81 - 1.0/10.0*w71 + (1.0/10.0)*w73 + (1.0/10.0)*w75 + w82;
        const scalar w84 = gin*w44;
        const scalar w85 = gil*w47;
        const scalar w86 = git*w49;
        const scalar w87 = z*(w84 + w85 + w86);
        const scalar w88 = gin*w27;
        const scalar w89 = gil*w36;
        const scalar w90 = w35*w89;
        const scalar w91 = git*w41;
        const scalar w92 = w40*w91;
        const scalar w93 = -w88 + w90 + w92;
        const scalar w94 = w93*y;
        const scalar w95 = gin*w29;
        const scalar w96 = w26*w95;
        const scalar w97 = w55*w89;
        const scalar w98 = w57*w91;
        const scalar w99 = w96 + w97 + w98;
        const scalar w100 = w99*x;
        const scalar w101 = gin*w63;
        const scalar w102 = gil*w65;
        const scalar w103 = git*w67;
        const scalar w104 = y*(w101 + w102 + w103);
        const scalar w105 = w93*z;
        const scalar w106 = w70*w95;
        const scalar w107 = gil*w72;
        const scalar w108 = git*w74;
        const scalar w109 = -w106 + w107 + w108;
        const scalar w110 = w109*x;
        const scalar w111 = gin*w78;
        const scalar w112 = gil*w79;
        const scalar w113 = git*w80;
        const scalar w114 = w111 + w112 + w113;
        const scalar w115 = w114*x;
        const scalar w116 = w99*z;
        const scalar w117 = w109*y;
        const scalar w118 = x*(w115 + w116 + w117) + y*(w104 + w105 + w110) + z*(w100 + w87 + w94);
        const scalar w119 = (1.0/(w118));
        const scalar w120 = w7*(-w61*w99 + (1.0/5.0)*w87 - 1.0/10.0*w88 + (1.0/10.0)*w90 + (1.0/10.0)*w92);
        const scalar w121 = (1.0/10.0)*w101 + (1.0/10.0)*w102 + (1.0/10.0)*w103 + (1.0/5.0)*w105 - w109*w61;
        const scalar w122 = -1.0/10.0*w106 + (1.0/10.0)*w107 + (1.0/10.0)*w108 - w114*w61 + (1.0/5.0)*w116;
        const scalar w123 = w120*z + w121*w8 + w122*w7*x;
        const scalar w124 = w119*w123;
        const scalar w125 = (r2*r2);
        const scalar w126 = w7*z;
        const scalar w127 = w7*x;
        const scalar w128 = w7*(w126*w51 + w127*w59 + w69*w8);
        const scalar w129 = w64 + w66 + w68;
        const scalar w130 = w7*(w126*w69 + w127*w76 + w129*w8);
        const scalar w131 = w7*(w126*w59 + w127*w81 + w76*w8);
        const scalar w132 = w128*z + w130*y + w131*x;
        const scalar w133 = w125*w132;
        const scalar w134 = -w12*w22 - w124*w133 + w62*w7*z + w7*w77*y + w7*w83*x;
        const scalar w135 = (1.0/(sigma0Out));
        const scalar w136 = (1.0/(w21));
        const scalar w137 = w135*w136;
        const scalar w138 = w134*w137;
        const scalar w139 = (1.0/(w21*w21));
        const scalar w140 = sin(w3);
        const scalar w141 = w140*w19;
        const scalar w142 = w139*w141;
        const scalar w143 = w134*w135;
        const scalar w144 = w142*w143;
        const scalar w145 = -w125 + w9 + (x*x) + (y*y);
        const scalar w146 = w145*w7;
        const scalar w147 = (1.0/2.0)*w146;
        const scalar w148 = k*w147;
        const scalar w149 = w22*w7;
        const scalar w150 = r2*w119*w121;
        const scalar w151 = (1.0/(w125));
        const scalar w152 = w151*y;
        const scalar w153 = w151*z;
        const scalar w154 = w151*x;
        const scalar w155 = w129*w152 + w130 + w153*w69 + w154*w76;
        const scalar w156 = w124*w125;
        const scalar w157 = -2*w104 - 2*w105 - 2*w110;
        const scalar w158 = (1.0/(w118*w118));
        const scalar w159 = w123*w158;
        const scalar w160 = w133*w159;
        const scalar w161 = k*sigma0Out*sigmaEps*w12*w140*w16*w18 - w132*w150 - 1.0/10.0*w149 - w155*w156 - w157*w160 + w7*w77;
        const scalar w162 = w137*w147;
        const scalar w163 = k*sigma0Out;
        const scalar w164 = w6*cos(w4);
        const scalar w165 = cos(w14);
        const scalar w166 = sigmaEps*w13;
        const scalar w167 = w166*w18;
        const scalar w168 = w165*w167;
        const scalar w169 = w135*w146;
        const scalar w170 = k*w169;
        const scalar w171 = w139*w170;
        const scalar w172 = (1.0/5.0)*w45 + (1.0/5.0)*w48 + (1.0/5.0)*w50;
        const scalar w173 = k*sigma0Out*w12*w166;
        const scalar w174 = (1.0/2.0)*w165;
        const scalar w175 = w128 + w152*w69 + w153*w51 + w154*w59;
        const scalar w176 = (1.0/5.0)*w84 + (1.0/5.0)*w85 + (1.0/5.0)*w86;
        const scalar w177 = w120 + w126*w176 + w127*((1.0/5.0)*w96 + (1.0/5.0)*w97 + (1.0/5.0)*w98) + w8*(-1.0/5.0*w88 + (1.0/5.0)*w90 + (1.0/5.0)*w92);
        const scalar w178 = w119*w133;
        const scalar w179 = -2*w100 - 2*w87 - 2*w94;
        const scalar w180 = -w149*w5 - w156*w175 - w160*w179 + w172*w7*z - w173*w174*w18 - w177*w178 + w62*w7 + w7*x*((1.0/5.0)*w54 + (1.0/5.0)*w56 + (1.0/5.0)*w58) + w7*y*(-1.0/5.0*w28 + (1.0/5.0)*w38 + (1.0/5.0)*w43);
        const scalar w181 = cos(w0);
        const scalar w182 = (1.0/2.0)*y + 1;
        const scalar w183 = cos(w17);
        const scalar w184 = w16*w183;
        const scalar w185 = w166*w184;
        const scalar w186 = w139*w143;
        const scalar w187 = w147*w186;
        const scalar w188 = (k*k);
        const scalar w189 = w181*w188;
        const scalar w190 = w10*w189;
        const scalar w191 = w190*x;
        const scalar w192 = -w10*w60 - w191;
        const scalar w193 = w189*w7;
        const scalar w194 = w131 + w152*w76 + w153*w59 + w154*w81;
        const scalar w195 = -w115*w190 - w116*w190 - w117*w190 + w122*w7;
        const scalar w196 = -2*w115 - 2*w116 - 2*w117;
        const scalar w197 = -w156*w194 - w160*w196 - w173*w184 - w178*w195 - w190*w76*y - w191*w81 - w192*w22 - w193*w82 + w7*w83;
        const scalar w198 = w163*w185;
        const scalar w199 = St*w164;
        const scalar w200 = -sigma0Out*sigmaEps*w12*w13*w16*w18*w188;
        const scalar w201 = 2*w155;
        const scalar w202 = 2*w124;
        const scalar w203 = r2*w132;
        const scalar w204 = 2*w158;
        const scalar w205 = w125*w159;
        const scalar w206 = w123*w133/(w118*w118*w118);
        const scalar w207 = 2*w137;
        const scalar w208 = (w18*w18);
        const scalar w209 = (sigmaEps*sigmaEps);
        const scalar w210 = (1.0/(w21*w21*w21));
        const scalar w211 = w134*w169;
        const scalar w212 = (w16*w16)*w188*w209*w210*w211;
        const scalar w213 = 2*w7;
        const scalar w214 = w138*w7;
        const scalar w215 = w187*w188*w20 + w214;
        const scalar w216 = w14*w168*w7;
        const scalar w217 = (1.0/4.0)*w211;
        const scalar w218 = (w13*w13);
        const scalar w219 = 2*w119;
        const scalar w220 = w125*w219;
        const scalar w221 = 2*w205;
        const scalar w222 = w133*w204;
        const scalar w223 = (2.0/5.0)*w193;
        const scalar w224 = (k*k*k);
        const scalar w225 = w11*w224;
        return -w141*w163*(St*(fA1*z - k*w6*sin(w4) + (1.0/2.0)*w2) + Tt*(w138*w8 + w144*w148 + w161*w162 + 1.0/10.0)) + w163*w167*w174*(St*(fA1*y + (2.0/5.0)*w164) + Tt*(w126*w138 - 1.0/4.0*w134*w168*w171 + w162*w180 + w52)) + w198*(St*(fA0 + fA2*k*w181*w182) + Tt*(-k*w185*w187 + w134*w135*w136*w7*x + (1.0/2.0)*w135*w136*w145*w197*w7 - w61)) + w22*(Tt*((w140*w140)*w208*w212 + w142*w161*w170 + w144*w213*w3 + w161*w207*w8 + w162*((1.0/5.0)*k*sigma0Out*sigmaEps*w140*w16*w18*w7 - w121*w157*w203*w204 - w129*w202 - w150*w201 - w157*w201*w205 - w157*w206*(-4*w104 - 4*w105 - 4*w110) - w160*(-2*w101 - 2*w102 - 2*w103) - w200) + w215) - w188*w199) + w22*(Tt*(w126*w180*w207 - w135*w139*w148*w168*w180 + w139*w15*w167*w188*w217 + w162*((1.0/2.0)*sigma0Out*sigmaEps*w12*w13*w15*w18*w188 - 2.0/5.0*sigma0Out*w216 - 2.0/5.0*w149 - w160*(-2*w84 - 2*w85 - 2*w86) + 2*w172*w7 - w175*w177*w220 - w175*w179*w221 - w176*w203*w219 - w177*w179*w222 - w179*w206*(-4*w100 - 4*w87 - 4*w94) - w202*w51) + (w165*w165)*w188*w208*w209*w210*w217*w218 - w186*w216 + w214 + 1.0/5.0) + (4.0/25.0)*w199) + w22*(-St*w182*w188*w2 + Tt*(-w0*w185*w186*w213 + w127*w197*w207 + w162*((1.0/5.0)*w1*w224*w59*w7*z + (1.0/5.0)*w1*w224*w7*w76*y + (1.0/5.0)*w1*w224*w7*w81*x - w160*(-2*w111 - 2*w112 - 2*w113) - w178*(-w114*w223 + w115*w225 + w116*w225 + w117*w225) - 2*w192*w198 - w194*w195*w220 - w194*w196*w221 - w195*w196*w222 - w196*w206*(-4*w115 - 4*w116 - 4*w117) - w200 - w202*w81 - w22*(-w223 + w225*x) - w223*w81) - w171*w185*w197 + (w183*w183)*w212*w218 - 1.0/5.0*w189 + w215));
    }
    const scalar w0 = (r1*r1);
    const scalar w1 = -(r2*r2);
    const scalar w2 = -w0 - w1;
    const scalar w3 = (1.0/(w2));
    const scalar w4 = theta1*w3;
    const scalar w5 = 4*w4;
    const scalar w6 = (x*x);
    const scalar w7 = (y*y);
    const scalar w8 = (z*z);
    const scalar w9 = w6 + w7 + w8;
    const scalar w10 = w3*(-w0 + w9);
    const scalar w11 = theta0 + theta1*w10;
    const scalar w12 = cos(w11);
    const scalar w13 = sin(w11);
    const scalar w14 = sin(gamma2);
    const scalar w15 = (w14*w14);
    const scalar w16 = w12*w13*w15;
    const scalar w17 = w16*w5;
    const scalar w18 = gil*w17;
    const scalar w19 = git*w5;
    const scalar w20 = w16*w19;
    const scalar w21 = w18*z - w20*z;
    const scalar w22 = k*x;
    const scalar w23 = k*z;
    const scalar w24 = cos(w23);
    const scalar w25 = hA2*w24;
    const scalar w26 = hA0 + w22*w25;
    const scalar w27 = St*w26;
    const scalar w28 = r1 + r2;
    const scalar w29 = (1.0/(w28));
    const scalar w30 = w10*w29;
    const scalar w31 = (1.0/(r2));
    const scalar w32 = cos(gamma2);
    const scalar w33 = (w32*w32);
    const scalar w34 = gin*w33;
    const scalar w35 = (w13*w13);
    const scalar w36 = w15*w35;
    const scalar w37 = gil*w36;
    const scalar w38 = (w12*w12);
    const scalar w39 = w15*w38;
    const scalar w40 = git*w39;
    const scalar w41 = w34 + w37 + w40;
    const scalar w42 = w31*w41;
    const scalar w43 = cos(gamma1);
    const scalar w44 = w14*w32;
    const scalar w45 = gin*w44;
    const scalar w46 = w43*w45;
    const scalar w47 = sin(gamma1);
    const scalar w48 = w12*w47;
    const scalar w49 = w13*w43;
    const scalar w50 = w32*w49;
    const scalar w51 = w48 + w50;
    const scalar w52 = gil*w51;
    const scalar w53 = w13*w14;
    const scalar w54 = w52*w53;
    const scalar w55 = w13*w47;
    const scalar w56 = w12*w32*w43 - w55;
    const scalar w57 = git*w56;
    const scalar w58 = w12*w14;
    const scalar w59 = w57*w58;
    const scalar w60 = -w46 + w54 + w59;
    const scalar w61 = w31*w60;
    const scalar w62 = w45*w47;
    const scalar w63 = w32*w55;
    const scalar w64 = w12*w43 - w63;
    const scalar w65 = gil*w64;
    const scalar w66 = w53*w65;
    const scalar w67 = w32*w48;
    const scalar w68 = -w49 - w67;
    const scalar w69 = git*w68;
    const scalar w70 = w58*w69;
    const scalar w71 = w62 + w66 + w70;
    const scalar w72 = w31*w71;
    const scalar w73 = w31*(w42*z + w61*y + w72*x);
    const scalar w74 = (w43*w43);
    const scalar w75 = gin*w15;
    const scalar w76 = w74*w75;
    const scalar w77 = (w51*w51);
    const scalar w78 = gil*w77;
    const scalar w79 = (w56*w56);
    const scalar w80 = git*w79;
    const scalar w81 = w76 + w78 + w80;
    const scalar w82 = w31*w81;
    const scalar w83 = w43*w47;
    const scalar w84 = w75*w83;
    const scalar w85 = w51*w65;
    const scalar w86 = w56*w69;
    const scalar w87 = -w84 + w85 + w86;
    const scalar w88 = w31*w87;
    const scalar w89 = w31*(w61*z + w82*y + w88*x);
    const scalar w90 = (w47*w47);
    const scalar w91 = (w64*w64);
    const scalar w92 = (w68*w68);
    const scalar w93 = gil*w91 + git*w92 + w75*w90;
    const scalar w94 = w31*w93;
    const scalar w95 = w31*(w72*z + w88*y + w94*x);
    const scalar w96 = w73*z + w89*y + w95*x;
    const scalar w97 = (1.0/(w96));
    const scalar w98 = 1 - 3.0/10.0*w6;
    const scalar w99 = St*k;
    const scalar w100 = k*y;
    const scalar w101 = sin(w100);
    const scalar w102 = hA1*w101;
    const scalar w103 = w102*w99;
    const scalar w104 = w103*w98;
    const scalar w105 = sin(w23);
    const scalar w106 = cos(w100);
    const scalar w107 = hA1*w106;
    const scalar w108 = (3.0/5.0)*w107;
    const scalar w109 = w108*x;
    const scalar w110 = hA2*w105 - w109;
    const scalar w111 = St*w110;
    const scalar w112 = -w104*w60 + w111*w71 + w27*w41;
    const scalar w113 = w112*w31;
    const scalar w114 = -w104*w81 + w111*w87 + w27*w60;
    const scalar w115 = w114*w31;
    const scalar w116 = -w104*w87 + w111*w93 + w27*w71;
    const scalar w117 = w116*w31;
    const scalar w118 = -w113*z - w115*y - w117*x;
    const scalar w119 = w118*w97;
    const scalar w120 = (1.0/2.0)*w28;
    const scalar w121 = w10*w31;
    const scalar w122 = w120*w121;
    const scalar w123 = (1.0/(r1));
    const scalar w124 = w112*w123;
    const scalar w125 = w114*w123;
    const scalar w126 = w116*w123;
    const scalar w127 = w124*z + w125*y + w126*x;
    const scalar w128 = w123*w41;
    const scalar w129 = w123*w60;
    const scalar w130 = w123*w71;
    const scalar w131 = w123*(-w128*z - w129*y - w130*x);
    const scalar w132 = w123*w81;
    const scalar w133 = w123*w87;
    const scalar w134 = w123*(-w129*z - w132*y - w133*x);
    const scalar w135 = w123*w93;
    const scalar w136 = w123*(-w130*z - w133*y - w135*x);
    const scalar w137 = -w131*z - w134*y - w136*x;
    const scalar w138 = (1.0/(w137));
    const scalar w139 = w123*(1 - w10);
    const scalar w140 = w120*w139;
    const scalar w141 = w138*w140;
    const scalar w142 = w119*w122 + w127*w141;
    const scalar w143 = 2*w142;
    const scalar w144 = w143*w30;
    const scalar w145 = w1 + w9;
    const scalar w146 = w29*w3;
    const scalar w147 = w145*w146;
    const scalar w148 = w143*w147;
    const scalar w149 = w31*z;
    const scalar w150 = w28*w3;
    const scalar w151 = w149*w150;
    const scalar w152 = w127*w138;
    const scalar w153 = w123*z;
    const scalar w154 = w150*w153;
    const scalar w155 = (k*k);
    const scalar w156 = St*w155;
    const scalar w157 = hA2*w105*w156;
    const scalar w158 = w157*x;
    const scalar w159 = w25*w99;
    const scalar w160 = w159*w71;
    const scalar w161 = 2*w4;
    const scalar w162 = w161*w58;
    const scalar w163 = w162*w52;
    const scalar w164 = w161*w53;
    const scalar w165 = w164*w57;
    const scalar w166 = w161*w55;
    const scalar w167 = 2*theta1*w12*w3*w32*w43*z - w166*z;
    const scalar w168 = gil*w167;
    const scalar w169 = w168*w53;
    const scalar w170 = w161*w48;
    const scalar w171 = w161*w50;
    const scalar w172 = -w170*z - w171*z;
    const scalar w173 = git*w172;
    const scalar w174 = w173*w58;
    const scalar w175 = w163*z - w165*z + w169 + w174;
    const scalar w176 = w162*w65;
    const scalar w177 = w164*w69;
    const scalar w178 = w161*w49;
    const scalar w179 = w161*w67;
    const scalar w180 = -w178*z - w179*z;
    const scalar w181 = gil*w180;
    const scalar w182 = w181*w53;
    const scalar w183 = w12*w43;
    const scalar w184 = w161*w183;
    const scalar w185 = w161*w63;
    const scalar w186 = -w184*z + w185*z;
    const scalar w187 = git*w186;
    const scalar w188 = w187*w58;
    const scalar w189 = w176*z - w177*z + w182 + w188;
    const scalar w190 = -w104*w175 + w111*w189 - w158*w41 + w160 + w21*w27;
    const scalar w191 = w157*w71;
    const scalar w192 = w49*w5;
    const scalar w193 = w5*w67;
    const scalar w194 = -w192*z - w193*z;
    const scalar w195 = w183*w5;
    const scalar w196 = w5*w63;
    const scalar w197 = -w195*z + w196*z;
    const scalar w198 = w194*w65 + w197*w69;
    const scalar w199 = w180*w52;
    const scalar w200 = w167*w65;
    const scalar w201 = w186*w57;
    const scalar w202 = w172*w69;
    const scalar w203 = w199 + w200 + w201 + w202;
    const scalar w204 = -w104*w203 + w111*w198 + w159*w93 + w189*w27 - w191*x;
    const scalar w205 = w204*w31;
    const scalar w206 = w5*w55;
    const scalar w207 = 4*theta1*w12*w3*w32*w43*z - w206*z;
    const scalar w208 = w207*w52;
    const scalar w209 = w48*w5;
    const scalar w210 = w5*w50;
    const scalar w211 = -w209*z - w210*z;
    const scalar w212 = w211*w57;
    const scalar w213 = w208 + w212;
    const scalar w214 = -w104*w213 + w111*w203 - w158*w60 + w159*w87 + w175*w27;
    const scalar w215 = w214*w31;
    const scalar w216 = -w113 - w149*w190 - w205*x - w215*y;
    const scalar w217 = w122*w97;
    const scalar w218 = w123*w204;
    const scalar w219 = w123*w214;
    const scalar w220 = w124 + w153*w190 + w218*x + w219*y;
    const scalar w221 = w175*w31;
    const scalar w222 = w189*w31;
    const scalar w223 = w149*w21 + w221*y + w222*x + w42;
    const scalar w224 = w149*w223;
    const scalar w225 = w213*w31;
    const scalar w226 = w203*w31;
    const scalar w227 = w31*(w221*z + w225*y + w226*x + w61);
    const scalar w228 = w227*y;
    const scalar w229 = w198*w31;
    const scalar w230 = w31*(w222*z + w226*y + w229*x + w72);
    const scalar w231 = w230*x;
    const scalar w232 = -w224 - w228 - w231 - w73;
    const scalar w233 = (1.0/(w96*w96));
    const scalar w234 = w118*w233;
    const scalar w235 = w232*w234;
    const scalar w236 = w127*w140;
    const scalar w237 = (1.0/(w137*w137));
    const scalar w238 = w123*w175;
    const scalar w239 = w123*w189;
    const scalar w240 = -w128 - w153*w21 - w238*y - w239*x;
    const scalar w241 = w153*w240;
    const scalar w242 = w123*w213;
    const scalar w243 = w123*w203;
    const scalar w244 = w123*(-w129 - w238*z - w242*y - w243*x);
    const scalar w245 = w244*y;
    const scalar w246 = w123*w198;
    const scalar w247 = w123*(-w130 - w239*z - w243*y - w246*x);
    const scalar w248 = w247*x;
    const scalar w249 = w131 + w241 + w245 + w248;
    const scalar w250 = w237*w249;
    const scalar w251 = w119*w151 + w122*w235 + w141*w220 - w152*w154 + w216*w217 + w236*w250;
    const scalar w252 = w145*w30;
    const scalar w253 = (1.0/5.0)*z;
    const scalar w254 = sin(w22);
    const scalar w255 = (1.0/5.0)*w254;
    const scalar w256 = k*w255;
    const scalar w257 = w253*w41 - w256*w71 - 1.0/10.0*w46 + (1.0/10.0)*w54 + (1.0/10.0)*w59;
    const scalar w258 = w257*w31;
    const scalar w259 = w253*w60 - w256*w87 + (1.0/10.0)*w76 + (1.0/10.0)*w78 + (1.0/10.0)*w80;
    const scalar w260 = w259*w31;
    const scalar w261 = w253*w71 - w256*w93 - 1.0/10.0*w84 + (1.0/10.0)*w85 + (1.0/10.0)*w86;
    const scalar w262 = w261*w31;
    const scalar w263 = -w258*z - w260*y - w262*x;
    const scalar w264 = w263*w97;
    const scalar w265 = w123*w257;
    const scalar w266 = w123*w259;
    const scalar w267 = w123*w261;
    const scalar w268 = w265*z + w266*y + w267*x;
    const scalar w269 = w138*w268;
    const scalar w270 = w122*w264 + w140*w269;
    const scalar w271 = 2*w270;
    const scalar w272 = w271*w30;
    const scalar w273 = w147*w271;
    const scalar w274 = w52*w58;
    const scalar w275 = w253*w4;
    const scalar w276 = w53*w57;
    const scalar w277 = (1.0/10.0)*w169 + (1.0/10.0)*w174 - w189*w256 + w21*w253 + w274*w275 - w275*w276 + (1.0/5.0)*w34 + (1.0/5.0)*w37 + (1.0/5.0)*w40;
    const scalar w278 = w175*w253 - w203*w256 + (1.0/10.0)*w208 + (1.0/10.0)*w212 - 1.0/5.0*w46 + (1.0/5.0)*w54 + (1.0/5.0)*w59;
    const scalar w279 = w278*w31;
    const scalar w280 = w189*w253 - w198*w256 + (1.0/10.0)*w199 + (1.0/10.0)*w200 + (1.0/10.0)*w201 + (1.0/10.0)*w202 + (1.0/5.0)*w62 + (1.0/5.0)*w66 + (1.0/5.0)*w70;
    const scalar w281 = w280*w31;
    const scalar w282 = -w149*w277 - w258 - w279*y - w281*x;
    const scalar w283 = w123*w278;
    const scalar w284 = w123*w280;
    const scalar w285 = w153*w277 + w265 + w283*y + w284*x;
    const scalar w286 = w233*w263;
    const scalar w287 = w232*w286;
    const scalar w288 = w140*w268;
    const scalar w289 = w122*w287 + w141*w285 + w151*w264 - w154*w269 + w217*w282 + w250*w288;
    const scalar w290 = Tt*(w252*w289 + w253 + w272*z + w273*z);
    const scalar w291 = w144*z + w148*z + w251*w252 + w27 + w290;
    const scalar w292 = 2*theta1*w12*w3*w32*w43*y - w166*y;
    const scalar w293 = gil*w292;
    const scalar w294 = w293*w53;
    const scalar w295 = -w170*y - w171*y;
    const scalar w296 = git*w295;
    const scalar w297 = w296*w58;
    const scalar w298 = w163*y - w165*y + w294 + w297;
    const scalar w299 = -w178*x - w179*x;
    const scalar w300 = gil*w53;
    const scalar w301 = w299*w300;
    const scalar w302 = -w184*x + w185*x;
    const scalar w303 = git*w58;
    const scalar w304 = w302*w303;
    const scalar w305 = w176*x - w177*x + w301 + w304;
    const scalar w306 = 4*theta1*w12*w3*w32*w43*y - w206*y;
    const scalar w307 = w306*w52;
    const scalar w308 = -w209*y - w210*y;
    const scalar w309 = w308*w57;
    const scalar w310 = w307 + w309;
    const scalar w311 = w31*y;
    const scalar w312 = w150*w311;
    const scalar w313 = w123*y;
    const scalar w314 = w150*w313;
    const scalar w315 = w18*y - w20*y;
    const scalar w316 = w102*w22;
    const scalar w317 = (3.0/5.0)*St;
    const scalar w318 = w316*w317;
    const scalar w319 = w107*w156*w98;
    const scalar w320 = -w178*y - w179*y;
    const scalar w321 = gil*w320;
    const scalar w322 = w321*w53;
    const scalar w323 = -w184*y + w185*y;
    const scalar w324 = git*w323;
    const scalar w325 = w324*w58;
    const scalar w326 = w176*y - w177*y + w322 + w325;
    const scalar w327 = -w104*w298 + w111*w326 + w27*w315 + w318*w71 - w319*w60;
    const scalar w328 = w31*w327;
    const scalar w329 = w102*w87;
    const scalar w330 = w22*w317*w329;
    const scalar w331 = w320*w52;
    const scalar w332 = w292*w65;
    const scalar w333 = w323*w57;
    const scalar w334 = w295*w69;
    const scalar w335 = w331 + w332 + w333 + w334;
    const scalar w336 = -w104*w310 + w111*w335 + w27*w298 - w319*w81 + w330;
    const scalar w337 = -w192*y - w193*y;
    const scalar w338 = -w195*y + w196*y;
    const scalar w339 = w337*w65 + w338*w69;
    const scalar w340 = -w104*w335 + w111*w339 + w27*w326 + w318*w93 - w319*w87;
    const scalar w341 = w31*w340;
    const scalar w342 = -w115 - w311*w336 - w328*z - w341*x;
    const scalar w343 = w123*w327;
    const scalar w344 = w123*w340;
    const scalar w345 = w125 + w313*w336 + w343*z + w344*x;
    const scalar w346 = w31*w315;
    const scalar w347 = w298*w31;
    const scalar w348 = w31*w326;
    const scalar w349 = w31*(w346*z + w347*y + w348*x + w61);
    const scalar w350 = w349*z;
    const scalar w351 = w31*w335;
    const scalar w352 = w310*w311 + w347*z + w351*x + w82;
    const scalar w353 = w311*w352;
    const scalar w354 = w31*w339;
    const scalar w355 = w31*(w348*z + w351*y + w354*x + w88);
    const scalar w356 = w355*x;
    const scalar w357 = -w350 - w353 - w356 - w89;
    const scalar w358 = w122*w357;
    const scalar w359 = w123*w315;
    const scalar w360 = w123*w298;
    const scalar w361 = w123*w326;
    const scalar w362 = w123*(-w129 - w359*z - w360*y - w361*x);
    const scalar w363 = w362*z;
    const scalar w364 = w123*w335;
    const scalar w365 = -w132 - w310*w313 - w360*z - w364*x;
    const scalar w366 = w313*w365;
    const scalar w367 = w123*w339;
    const scalar w368 = w123*(-w133 - w361*z - w364*y - w367*x);
    const scalar w369 = w368*x;
    const scalar w370 = w134 + w363 + w366 + w369;
    const scalar w371 = w237*w370;
    const scalar w372 = w119*w312 + w141*w345 - w152*w314 + w217*w342 + w234*w358 + w236*w371;
    const scalar w373 = (1.0/5.0)*w4;
    const scalar w374 = w274*w373;
    const scalar w375 = w276*w373;
    const scalar w376 = w253*w315 - w256*w326 + (1.0/10.0)*w294 + (1.0/10.0)*w297 + w374*y - w375*y;
    const scalar w377 = w31*w376;
    const scalar w378 = w253*w298 - w256*w335 + (1.0/10.0)*w307 + (1.0/10.0)*w309;
    const scalar w379 = w253*w326 - w256*w339 + (1.0/10.0)*w331 + (1.0/10.0)*w332 + (1.0/10.0)*w333 + (1.0/10.0)*w334;
    const scalar w380 = w31*w379;
    const scalar w381 = -w260 - w311*w378 - w377*z - w380*x;
    const scalar w382 = w123*w376;
    const scalar w383 = w123*w379;
    const scalar w384 = w266 + w313*w378 + w382*z + w383*x;
    const scalar w385 = w141*w384 + w217*w381 + w264*w312 - w269*w314 + w286*w358 + w288*w371;
    const scalar w386 = Tt*(w252*w385 + w272*y + w273*y + 1.0/10.0);
    const scalar w387 = -w104 + w144*y + w148*y + w252*w372 + w386;
    const scalar w388 = gel*w17;
    const scalar w389 = w388*z;
    const scalar w390 = get*w17;
    const scalar w391 = w390*z;
    const scalar w392 = w389 - w391;
    const scalar w393 = w100 + 3.0/10.0;
    const scalar w394 = cos(w393);
    const scalar w395 = (2.0/5.0)*z;
    const scalar w396 = fA3*exp(w395);
    const scalar w397 = (2.0/5.0)*w396;
    const scalar w398 = w394*w397;
    const scalar w399 = St*(fA1*y + w398);
    const scalar w400 = gel*w36;
    const scalar w401 = get*w39;
    const scalar w402 = gen*w33 + w400 + w401;
    const scalar w403 = w31*w402;
    const scalar w404 = gen*w44;
    const scalar w405 = gel*w51;
    const scalar w406 = w405*w53;
    const scalar w407 = get*w56;
    const scalar w408 = w407*w58;
    const scalar w409 = -w404*w43 + w406 + w408;
    const scalar w410 = w31*w409;
    const scalar w411 = gel*w64;
    const scalar w412 = w411*w53;
    const scalar w413 = get*w68;
    const scalar w414 = w413*w58;
    const scalar w415 = w404*w47 + w412 + w414;
    const scalar w416 = w31*w415;
    const scalar w417 = w31*(w403*z + w410*y + w416*x);
    const scalar w418 = gen*w15;
    const scalar w419 = gel*w77 + get*w79 + w418*w74;
    const scalar w420 = w31*w419;
    const scalar w421 = w411*w51 + w413*w56 - w418*w83;
    const scalar w422 = w31*w421;
    const scalar w423 = w31*(w410*z + w420*y + w422*x);
    const scalar w424 = gel*w91 + get*w92 + w418*w90;
    const scalar w425 = w31*w424;
    const scalar w426 = w31*(w416*z + w422*y + w425*x);
    const scalar w427 = w417*z + w423*y + w426*x;
    const scalar w428 = (1.0/(w427));
    const scalar w429 = w22 + 1.0/5.0;
    const scalar w430 = sin(w429);
    const scalar w431 = (1.0/2.0)*w105;
    const scalar w432 = w431 + 1;
    const scalar w433 = sigmaEps*w432;
    const scalar w434 = w106*w433;
    const scalar w435 = w430*w434;
    const scalar w436 = w435 + 1;
    const scalar w437 = -sigma0Out*w436;
    const scalar w438 = -w402 - w437;
    const scalar w439 = -w415;
    const scalar w440 = (1.0/2.0)*y + 1;
    const scalar w441 = cos(w22);
    const scalar w442 = fA2*k*w441;
    const scalar w443 = St*(fA0 + w440*w442);
    const scalar w444 = -w409;
    const scalar w445 = fA2*w254;
    const scalar w446 = (1.0/2.0)*w445;
    const scalar w447 = sin(w393);
    const scalar w448 = k*w447;
    const scalar w449 = w396*w448;
    const scalar w450 = St*(fA1*z + w446 - w449);
    const scalar w451 = w439*w443 + w444*w450;
    const scalar w452 = w31*(w399*w438 + w451);
    const scalar w453 = -w419 - w437;
    const scalar w454 = -w421;
    const scalar w455 = w399*w444 + w443*w454;
    const scalar w456 = w31*(w450*w453 + w455);
    const scalar w457 = -w424 - w437;
    const scalar w458 = w399*w439 + w450*w454;
    const scalar w459 = w31*(w443*w457 + w458);
    const scalar w460 = w452*z + w456*y + w459*x;
    const scalar w461 = w428*w460;
    const scalar w462 = w123*w402;
    const scalar w463 = w123*w409;
    const scalar w464 = w123*w415;
    const scalar w465 = w123*(-w462*z - w463*y - w464*x);
    const scalar w466 = w123*w419;
    const scalar w467 = w123*w421;
    const scalar w468 = w123*(-w463*z - w466*y - w467*x);
    const scalar w469 = w123*w424;
    const scalar w470 = w123*(-w464*z - w467*y - w469*x);
    const scalar w471 = -w465*z - w468*y - w470*x;
    const scalar w472 = (1.0/(w471));
    const scalar w473 = -sigma0In*w436;
    const scalar w474 = -w402 - w473;
    const scalar w475 = w123*(w399*w474 + w451);
    const scalar w476 = -w419 - w473;
    const scalar w477 = w123*(w450*w476 + w455);
    const scalar w478 = -w424 - w473;
    const scalar w479 = w123*(w443*w478 + w458);
    const scalar w480 = -w475*z - w477*y - w479*x;
    const scalar w481 = w472*w480;
    const scalar w482 = w122*w461 + w140*w481;
    const scalar w483 = 2*w482;
    const scalar w484 = w30*w483;
    const scalar w485 = w147*w483;
    const scalar w486 = k*sigma0Out;
    const scalar w487 = w430*w486;
    const scalar w488 = sigmaEps*w106;
    const scalar w489 = (1.0/2.0)*w24;
    const scalar w490 = w488*w489;
    const scalar w491 = -w389 + w391;
    const scalar w492 = w487*w490 + w491;
    const scalar w493 = w394*w396;
    const scalar w494 = St*w493;
    const scalar w495 = (4.0/25.0)*w494;
    const scalar w496 = w162*w411;
    const scalar w497 = gel*w180;
    const scalar w498 = get*w186;
    const scalar w499 = -2*get*theta1*w13*w14*w3*w68*z + w496*z + w497*w53 + w498*w58;
    const scalar w500 = -w499;
    const scalar w501 = St*(fA1 - w397*w448);
    const scalar w502 = w444*w501;
    const scalar w503 = w162*w405;
    const scalar w504 = gel*w167;
    const scalar w505 = get*w172;
    const scalar w506 = -2*get*theta1*w13*w14*w3*w56*z + w503*z + w504*w53 + w505*w58;
    const scalar w507 = -w506;
    const scalar w508 = w443*w500 + w450*w507 + w502;
    const scalar w509 = w399*w492 + w438*w495 + w508;
    const scalar w510 = -1.0/2.0*k*sigma0Out*sigmaEps*w106*w24*w430;
    const scalar w511 = w207*w405 + w211*w407;
    const scalar w512 = -w510 - w511;
    const scalar w513 = w167*w411 + w172*w413 + w180*w405 + w186*w407;
    const scalar w514 = -w513;
    const scalar w515 = w399*w507 + w443*w514 + w444*w495;
    const scalar w516 = w31*(w450*w512 + w453*w501 + w515);
    const scalar w517 = w194*w411 + w197*w413;
    const scalar w518 = -w510 - w517;
    const scalar w519 = w399*w500 + w439*w495 + w450*w514 + w454*w501;
    const scalar w520 = w31*(w443*w518 + w519);
    const scalar w521 = w149*w509 + w452 + w516*y + w520*x;
    const scalar w522 = w428*w521;
    const scalar w523 = k*sigma0In;
    const scalar w524 = w430*w523;
    const scalar w525 = w490*w524 + w491;
    const scalar w526 = w399*w525 + w474*w495 + w508;
    const scalar w527 = -1.0/2.0*k*sigma0In*sigmaEps*w106*w24*w430;
    const scalar w528 = -w511 - w527;
    const scalar w529 = w123*(w450*w528 + w476*w501 + w515);
    const scalar w530 = -w517 - w527;
    const scalar w531 = w123*(w443*w530 + w519);
    const scalar w532 = -w153*w526 - w475 - w529*y - w531*x;
    const scalar w533 = w140*w472;
    const scalar w534 = w31*w506;
    const scalar w535 = w31*w499;
    const scalar w536 = w149*w392 + w403 + w534*y + w535*x;
    const scalar w537 = w149*w536;
    const scalar w538 = w31*w511;
    const scalar w539 = w31*w513;
    const scalar w540 = w31*(w410 + w534*z + w538*y + w539*x);
    const scalar w541 = w540*y;
    const scalar w542 = w31*w517;
    const scalar w543 = w31*(w416 + w535*z + w539*y + w542*x);
    const scalar w544 = w543*x;
    const scalar w545 = -w417 - w537 - w541 - w544;
    const scalar w546 = (1.0/(w427*w427));
    const scalar w547 = w460*w546;
    const scalar w548 = w545*w547;
    const scalar w549 = w123*w506;
    const scalar w550 = w123*w499;
    const scalar w551 = -w153*w392 - w462 - w549*y - w550*x;
    const scalar w552 = w153*w551;
    const scalar w553 = w123*w511;
    const scalar w554 = w123*w513;
    const scalar w555 = w123*(-w463 - w549*z - w553*y - w554*x);
    const scalar w556 = w555*y;
    const scalar w557 = w123*w517;
    const scalar w558 = w123*(-w464 - w550*z - w554*y - w557*x);
    const scalar w559 = w558*x;
    const scalar w560 = w465 + w552 + w556 + w559;
    const scalar w561 = (1.0/(w471*w471));
    const scalar w562 = w140*w561;
    const scalar w563 = w480*w562;
    const scalar w564 = w122*w522 + w122*w548 + w151*w461 - w154*w481 + w532*w533 + w560*w563;
    const scalar w565 = w252*w564 + w290 + w399 + w484*z + w485*z;
    const scalar w566 = w299*w52;
    const scalar w567 = 2*theta1*w12*w3*w32*w43*x - w166*x;
    const scalar w568 = w567*w65;
    const scalar w569 = w302*w57;
    const scalar w570 = -w170*x - w171*x;
    const scalar w571 = w570*w69;
    const scalar w572 = w566 + w568 + w569 + w571;
    const scalar w573 = gel*w292;
    const scalar w574 = get*w295;
    const scalar w575 = -2*get*theta1*w13*w14*w3*w56*y + w503*y + w53*w573 + w574*w58;
    const scalar w576 = gel*w53;
    const scalar w577 = get*w58;
    const scalar w578 = -2*get*theta1*w13*w14*w3*w68*x + w299*w576 + w302*w577 + w496*x;
    const scalar w579 = -w192*x - w193*x;
    const scalar w580 = w579*w64;
    const scalar w581 = -w195*x + w196*x;
    const scalar w582 = w581*w68;
    const scalar w583 = gil*w580 + git*w582;
    const scalar w584 = w31*x;
    const scalar w585 = w150*w584;
    const scalar w586 = w123*x;
    const scalar w587 = w150*w586;
    const scalar w588 = w18*x - w20*x;
    const scalar w589 = St*w108;
    const scalar w590 = w300*w567;
    const scalar w591 = w303*w570;
    const scalar w592 = w163*x - w165*x + w590 + w591;
    const scalar w593 = -w104*w592 + w111*w305 + w159*w41 + w27*w588 + w318*w60 - w589*w71;
    const scalar w594 = w31*w593;
    const scalar w595 = 4*theta1*w12*w3*w32*w43*x - w206*x;
    const scalar w596 = w52*w595;
    const scalar w597 = -w209*x - w210*x;
    const scalar w598 = w57*w597;
    const scalar w599 = w596 + w598;
    const scalar w600 = -w104*w599 + w111*w572 + w159*w60 + w27*w592 + w318*w81 - w589*w87;
    const scalar w601 = w31*w600;
    const scalar w602 = -w104*w572 + w111*w583 + w160 + w27*w305 + w330 - w589*w93;
    const scalar w603 = -w117 - w584*w602 - w594*z - w601*y;
    const scalar w604 = w123*w593;
    const scalar w605 = w123*w600;
    const scalar w606 = w126 + w586*w602 + w604*z + w605*y;
    const scalar w607 = w31*w588;
    const scalar w608 = w31*w592;
    const scalar w609 = w305*w31;
    const scalar w610 = w31*(w607*z + w608*y + w609*x + w72);
    const scalar w611 = w610*z;
    const scalar w612 = w31*w572;
    const scalar w613 = w583*w584 + w609*z + w612*y + w94;
    const scalar w614 = w584*w613;
    const scalar w615 = w31*w599;
    const scalar w616 = w31*(w608*z + w612*x + w615*y + w88);
    const scalar w617 = w616*y;
    const scalar w618 = -w611 - w614 - w617 - w95;
    const scalar w619 = w233*w618;
    const scalar w620 = w122*w619;
    const scalar w621 = w123*w588;
    const scalar w622 = w123*w592;
    const scalar w623 = w123*w305;
    const scalar w624 = w123*(-w130 - w621*z - w622*y - w623*x);
    const scalar w625 = w624*z;
    const scalar w626 = w123*w572;
    const scalar w627 = -w135 - w583*w586 - w623*z - w626*y;
    const scalar w628 = w586*w627;
    const scalar w629 = w123*w599;
    const scalar w630 = w123*(-w133 - w622*z - w626*x - w629*y);
    const scalar w631 = w630*y;
    const scalar w632 = w136 + w625 + w628 + w631;
    const scalar w633 = w237*w632;
    const scalar w634 = w118*w620 + w119*w585 + w141*w606 - w152*w587 + w217*w603 + w236*w633;
    const scalar w635 = w155*w441;
    const scalar w636 = (1.0/5.0)*w635;
    const scalar w637 = w253*w588 - w256*w305 + w374*x - w375*x + (1.0/10.0)*w590 + (1.0/10.0)*w591 - w636*w71;
    const scalar w638 = w31*w637;
    const scalar w639 = w253*w305 - w256*w583 + (1.0/10.0)*w566 + (1.0/10.0)*w568 + (1.0/10.0)*w569 + (1.0/10.0)*w571 - w636*w93;
    const scalar w640 = w253*w592 - w256*w572 + (1.0/10.0)*w596 + (1.0/10.0)*w598 - w636*w87;
    const scalar w641 = w31*w640;
    const scalar w642 = -w262 - w584*w639 - w638*z - w641*y;
    const scalar w643 = w123*w637;
    const scalar w644 = w123*w640;
    const scalar w645 = w267 + w586*w639 + w643*z + w644*y;
    const scalar w646 = w138*w645;
    const scalar w647 = w140*w646 + w217*w642 + w263*w620 + w264*w585 - w269*w587 + w288*w633;
    const scalar w648 = Tt*(w252*w647 - w256 + w272*x + w273*x);
    const scalar w649 = w111 + w144*x + w148*x + w252*w634 + w648;
    const scalar w650 = w306*w405 + w308*w407;
    const scalar w651 = w101*w433;
    const scalar w652 = w487*w651;
    const scalar w653 = w388*y - w390*y;
    const scalar w654 = -w652 - w653;
    const scalar w655 = gel*w320;
    const scalar w656 = get*w323;
    const scalar w657 = -2*get*theta1*w13*w14*w3*w68*y + w496*y + w53*w655 + w58*w656;
    const scalar w658 = -w657;
    const scalar w659 = -w575;
    const scalar w660 = St*w442;
    const scalar w661 = (1.0/2.0)*w660;
    const scalar w662 = w156*w493;
    const scalar w663 = w439*w661 + w443*w658 - w444*w662 + w450*w659;
    const scalar w664 = w31*(w399*w654 + w438*w501 + w663);
    const scalar w665 = -w650 - w652;
    const scalar w666 = w454*w661;
    const scalar w667 = w292*w411 + w295*w413 + w320*w405 + w323*w407;
    const scalar w668 = -w667;
    const scalar w669 = w399*w659 + w443*w668 + w502 + w666;
    const scalar w670 = w450*w665 - w453*w662 + w669;
    const scalar w671 = w337*w411 + w338*w413;
    const scalar w672 = -w652 - w671;
    const scalar w673 = w399*w658 + w439*w501 + w450*w668 - w454*w662;
    const scalar w674 = w31*(w443*w672 + w457*w661 + w673);
    const scalar w675 = w311*w670 + w456 + w664*z + w674*x;
    const scalar w676 = w428*w675;
    const scalar w677 = w524*w651;
    const scalar w678 = -w653 - w677;
    const scalar w679 = w123*(w399*w678 + w474*w501 + w663);
    const scalar w680 = -w650 - w677;
    const scalar w681 = w450*w680 - w476*w662 + w669;
    const scalar w682 = -w671 - w677;
    const scalar w683 = w123*(w443*w682 + w478*w661 + w673);
    const scalar w684 = -w313*w681 - w477 - w679*z - w683*x;
    const scalar w685 = w31*w653;
    const scalar w686 = w31*w575;
    const scalar w687 = w31*w657;
    const scalar w688 = w31*(w410 + w685*z + w686*y + w687*x);
    const scalar w689 = w688*z;
    const scalar w690 = w31*w667;
    const scalar w691 = w311*w650 + w420 + w686*z + w690*x;
    const scalar w692 = w311*w691;
    const scalar w693 = w31*w671;
    const scalar w694 = w31*(w422 + w687*z + w690*y + w693*x);
    const scalar w695 = w694*x;
    const scalar w696 = -w423 - w689 - w692 - w695;
    const scalar w697 = w547*w696;
    const scalar w698 = w123*w653;
    const scalar w699 = w123*w575;
    const scalar w700 = w123*w657;
    const scalar w701 = w123*(-w463 - w698*z - w699*y - w700*x);
    const scalar w702 = w701*z;
    const scalar w703 = w123*w667;
    const scalar w704 = -w313*w650 - w466 - w699*z - w703*x;
    const scalar w705 = w313*w704;
    const scalar w706 = w123*w671;
    const scalar w707 = w123*(-w467 - w700*z - w703*y - w706*x);
    const scalar w708 = w707*x;
    const scalar w709 = w468 + w702 + w705 + w708;
    const scalar w710 = w122*w676 + w122*w697 + w312*w461 - w314*w481 + w533*w684 + w563*w709;
    const scalar w711 = w252*w710 + w386 + w450 + w484*y + w485*y;
    const scalar w712 = gel*w580 + get*w582;
    const scalar w713 = cos(w429);
    const scalar w714 = w434*w713;
    const scalar w715 = w388*x;
    const scalar w716 = w390*x;
    const scalar w717 = -w715 + w716;
    const scalar w718 = w486*w714 + w717;
    const scalar w719 = -w578;
    const scalar w720 = -2*get*theta1*w13*w14*w3*w56*x + w503*x + w567*w576 + w570*w577;
    const scalar w721 = -w720;
    const scalar w722 = w156*w440*w445;
    const scalar w723 = -w439*w722 + w443*w719 + w444*w661 + w450*w721;
    const scalar w724 = w31*(w399*w718 + w723);
    const scalar w725 = -k*sigma0Out*sigmaEps*w106*w432*w713;
    const scalar w726 = w405*w595 + w407*w597;
    const scalar w727 = -w725 - w726;
    const scalar w728 = w299*w405 + w302*w407 + w411*w567 + w413*w570;
    const scalar w729 = -w728;
    const scalar w730 = w399*w721 + w443*w729 - w454*w722;
    const scalar w731 = w31*(w450*w727 + w453*w661 + w730);
    const scalar w732 = -w712 - w725;
    const scalar w733 = w399*w719 + w450*w729 + w666;
    const scalar w734 = w443*w732 - w457*w722 + w733;
    const scalar w735 = w459 + w584*w734 + w724*z + w731*y;
    const scalar w736 = w428*w735;
    const scalar w737 = w523*w714 + w717;
    const scalar w738 = w123*(w399*w737 + w723);
    const scalar w739 = -k*sigma0In*sigmaEps*w106*w432*w713;
    const scalar w740 = -w726 - w739;
    const scalar w741 = w123*(w450*w740 + w476*w661 + w730);
    const scalar w742 = -w712 - w739;
    const scalar w743 = w443*w742 - w478*w722 + w733;
    const scalar w744 = -w479 - w586*w743 - w738*z - w741*y;
    const scalar w745 = w715 - w716;
    const scalar w746 = w31*w745;
    const scalar w747 = w31*w720;
    const scalar w748 = w31*w578;
    const scalar w749 = w31*(w416 + w746*z + w747*y + w748*x);
    const scalar w750 = w749*z;
    const scalar w751 = w31*w728;
    const scalar w752 = w425 + w584*w712 + w748*z + w751*y;
    const scalar w753 = w584*w752;
    const scalar w754 = w31*w726;
    const scalar w755 = w31*(w422 + w747*z + w751*x + w754*y);
    const scalar w756 = w755*y;
    const scalar w757 = -w426 - w750 - w753 - w756;
    const scalar w758 = w547*w757;
    const scalar w759 = w123*w745;
    const scalar w760 = w123*w720;
    const scalar w761 = w123*w578;
    const scalar w762 = w123*(-w464 - w759*z - w760*y - w761*x);
    const scalar w763 = w762*z;
    const scalar w764 = w123*w728;
    const scalar w765 = -w469 - w586*w712 - w761*z - w764*y;
    const scalar w766 = w586*w765;
    const scalar w767 = w123*w726;
    const scalar w768 = w123*(-w467 - w760*z - w764*x - w767*y);
    const scalar w769 = w768*y;
    const scalar w770 = w470 + w763 + w766 + w769;
    const scalar w771 = w122*w736 + w122*w758 + w461*w585 - w481*w587 + w533*w744 + w563*w770;
    const scalar w772 = w252*w771 + w443 + w484*x + w485*x + w648;
    const scalar w773 = 8*w146;
    const scalar w774 = w773*w8;
    const scalar w775 = 2*w123;
    const scalar w776 = (theta1*theta1);
    const scalar w777 = (1.0/(w2*w2));
    const scalar w778 = w776*w777;
    const scalar w779 = w778*w8;
    const scalar w780 = 8*w779;
    const scalar w781 = -4*gil*theta1*w12*w13*w15*w3 + w20;
    const scalar w782 = 8*gil*w15*w38*w776*w777*w8 + 8*git*w15*w35*w776*w777*w8 - w37*w780 - w40*w780 - w781;
    const scalar w783 = 4*w779;
    const scalar w784 = w178 + w179;
    const scalar w785 = 4*w13*w32*w47*w776*w777*w8 - w183*w783 - w784;
    const scalar w786 = -w184 + w185;
    const scalar w787 = w49*w783 + w67*w783 + w786;
    const scalar w788 = w53*z;
    const scalar w789 = w5*w788;
    const scalar w790 = -2*gil*theta1*w12*w14*w3*w64 + w177;
    const scalar w791 = 4*gil*theta1*w12*w14*w180*w3*z + gil*w13*w14*w785 + git*w12*w14*w787 - w187*w789 - w66*w783 - w70*w783 - w790;
    const scalar w792 = -2*theta1*w12*w3*w32*w43 + w166;
    const scalar w793 = -w48*w783 - w50*w783 - w792;
    const scalar w794 = w183*w32;
    const scalar w795 = w170 + w171;
    const scalar w796 = 4*w13*w47*w776*w777*w8 - w783*w794 - w795;
    const scalar w797 = (2.0/5.0)*w779;
    const scalar w798 = git*z;
    const scalar w799 = -1.0/5.0*gil*theta1*w12*w14*w3*w51 + w375;
    const scalar w800 = (8.0/5.0)*gil*theta1*w12*w13*w15*w3*z + (2.0/5.0)*gil*theta1*w12*w14*w167*w3*z + (1.0/10.0)*gil*w13*w14*w793 + (1.0/10.0)*git*w12*w14*w796 - 8.0/5.0*w16*w4*w798 - w173*w395*w4*w53 - w256*w791 - w54*w797 - w59*w797 + (1.0/5.0)*w782*z - w799;
    const scalar w801 = (4.0/5.0)*w4;
    const scalar w802 = w801*z;
    const scalar w803 = w167*w207;
    const scalar w804 = gil*w803;
    const scalar w805 = w172*w211;
    const scalar w806 = git*w805;
    const scalar w807 = -4*theta1*w12*w3*w32*w43 + w206;
    const scalar w808 = -w48*w780 - w50*w780 - w807;
    const scalar w809 = w52*w808;
    const scalar w810 = w209 + w210;
    const scalar w811 = 8*w13*w47*w776*w777*w8 - w780*w794 - w810;
    const scalar w812 = w57*w811;
    const scalar w813 = -2*gil*theta1*w12*w14*w3*w51 + w165;
    const scalar w814 = 4*gil*theta1*w12*w14*w167*w3*z + gil*w13*w14*w793 + git*w12*w14*w796 - w173*w789 - w54*w783 - w59*w783 - w813;
    const scalar w815 = w167*w181;
    const scalar w816 = w172*w187;
    const scalar w817 = w52*w785;
    const scalar w818 = w65*w793;
    const scalar w819 = w57*w787;
    const scalar w820 = w69*w796;
    const scalar w821 = 2*w815 + 2*w816 + w817 + w818 + w819 + w820;
    const scalar w822 = (2.0/5.0)*w169 + (2.0/5.0)*w174 + w253*w814 - w256*w821 + w274*w802 - w276*w802 + (1.0/10.0)*w804 + (1.0/10.0)*w806 + (1.0/10.0)*w809 + (1.0/10.0)*w812;
    const scalar w823 = w180*w194;
    const scalar w824 = w186*w197;
    const scalar w825 = w192 + w193;
    const scalar w826 = 8*w13*w32*w47*w776*w777*w8 - w183*w780 - w825;
    const scalar w827 = -w195 + w196;
    const scalar w828 = w49*w780 + w67*w780 + w827;
    const scalar w829 = gil*w823 + git*w824 + w65*w826 + w69*w828;
    const scalar w830 = (2.0/5.0)*w182 + (2.0/5.0)*w188 + w253*w791 - w256*w829 + w58*w65*w802 - w69*w788*w801 + (1.0/5.0)*w815 + (1.0/5.0)*w816 + (1.0/10.0)*w817 + (1.0/10.0)*w818 + (1.0/10.0)*w819 + (1.0/10.0)*w820;
    const scalar w831 = w139*w28;
    const scalar w832 = w237*w831;
    const scalar w833 = w249*w832;
    const scalar w834 = w804 + w806 + w809 + w812;
    const scalar w835 = w153*(-w153*w782 - w21*w775 - w313*w814 - w586*w791) + w240*w775 + w313*(-w153*w814 - 2*w238 - w313*w834 - w586*w821) + w586*(-w153*w791 - 2*w239 - w313*w821 - w586*w829);
    const scalar w836 = w237*w288;
    const scalar w837 = 2*w31;
    const scalar w838 = w138*w154;
    const scalar w839 = 2*w838;
    const scalar w840 = w151*w97;
    const scalar w841 = 2*w840;
    const scalar w842 = w121*w28;
    const scalar w843 = w233*w842;
    const scalar w844 = w232*w843;
    const scalar w845 = (1.0/(w137*w137*w137));
    const scalar w846 = w288*w845;
    const scalar w847 = 2*w131 + 2*w241 + 2*w245 + 2*w248;
    const scalar w848 = w249*w847;
    const scalar w849 = -w149*(w149*w782 + w21*w837 + w311*w814 + w584*w791) - w223*w837 - w311*(w149*w814 + 2*w221 + w311*w834 + w584*w821) - w584*(w149*w791 + 2*w222 + w311*w821 + w584*w829);
    const scalar w850 = w122*w286;
    const scalar w851 = w250*w268;
    const scalar w852 = 2*w154;
    const scalar w853 = 2*w151;
    const scalar w854 = (1.0/(w96*w96*w96));
    const scalar w855 = w263*w854;
    const scalar w856 = -2*w224 - 2*w228 - 2*w231 - 2*w73;
    const scalar w857 = w122*w232;
    const scalar w858 = w856*w857;
    const scalar w859 = w150*w31;
    const scalar w860 = w859*w97;
    const scalar w861 = w123*w150;
    const scalar w862 = w138*w861;
    const scalar w863 = w263*w860 - w268*w862;
    const scalar w864 = 4*z;
    const scalar w865 = w289*w30;
    const scalar w866 = w147*w864;
    const scalar w867 = w272 + w273;
    const scalar w868 = Tt*(w252*(w141*(w153*w800 + w277*w775 + w313*w822 + w586*w830) + w217*(-w149*w800 - w277*w837 - w311*w822 - w584*w830) + w282*w841 + w282*w844 + w285*w833 - w285*w839 + w287*w853 + w835*w836 + w846*w848 + w849*w850 - w851*w852 + w855*w858 + w863) + w270*w774 + w289*w866 + w864*w865 + w867 + 1.0/5.0);
    const scalar w869 = (k*k*k);
    const scalar w870 = St*w869;
    const scalar w871 = w25*w870*x;
    const scalar w872 = 2*w158;
    const scalar w873 = 2*St*hA2*k*w189*w24 + St*w110*w791 + St*w26*w782 - w104*w814 - w191 - w21*w872 - w41*w871;
    const scalar w874 = 2*St*hA2*k*w198*w24 + St*w110*w829 + St*w26*w791 - w104*w821 - w157*w93 - w189*w872 - w71*w871;
    const scalar w875 = 2*St*hA2*k*w203*w24 + St*w110*w821 + St*w26*w814 - w104*w834 - w157*w87 - w175*w872 - w60*w871;
    const scalar w876 = w236*w237;
    const scalar w877 = w236*w845;
    const scalar w878 = w122*w234;
    const scalar w879 = w127*w250;
    const scalar w880 = w118*w854;
    const scalar w881 = w118*w860 - w127*w862;
    const scalar w882 = w30*w864;
    const scalar w883 = w144 + w148;
    const scalar w884 = sigma0In*w155;
    const scalar w885 = w430*w884;
    const scalar w886 = w431*w488;
    const scalar w887 = w885*w886;
    const scalar w888 = w400*w780;
    const scalar w889 = w401*w780;
    const scalar w890 = gel*w39*w780;
    const scalar w891 = get*w36*w780;
    const scalar w892 = w388 - w390;
    const scalar w893 = -w888 - w889 + w890 + w891 + w892;
    const scalar w894 = (8.0/25.0)*w494;
    const scalar w895 = (8.0/125.0)*w494;
    const scalar w896 = w5*w58;
    const scalar w897 = w896*z;
    const scalar w898 = w164*w413 - w496;
    const scalar w899 = w412*w783 + w414*w783 - w497*w897 + w498*w789 - w576*w785 - w577*w787 + w898;
    const scalar w900 = w164*w407 - w503;
    const scalar w901 = w406*w783 + w408*w783 - w504*w897 + w505*w789 - w576*w793 - w577*w796 + w900;
    const scalar w902 = w501*w507;
    const scalar w903 = (4.0/25.0)*St*w449;
    const scalar w904 = -w444*w903;
    const scalar w905 = w443*w899 + w450*w901 + 2*w902 + w904;
    const scalar w906 = gel*w803 + get*w805 + w405*w808 + w407*w811;
    const scalar w907 = 2*w501;
    const scalar w908 = 2*w497;
    const scalar w909 = 2*w498;
    const scalar w910 = w167*w908 + w172*w909 + w405*w785 + w407*w787 + w411*w793 + w413*w796;
    const scalar w911 = -w910;
    const scalar w912 = w399*w901 + w443*w911 + w444*w895 + w507*w894;
    const scalar w913 = gel*w823 + get*w824 + w411*w826 + w413*w828;
    const scalar w914 = w399*w899 + w439*w895 + w450*w911 - w454*w903 + w500*w894 + w514*w907;
    const scalar w915 = w561*w831;
    const scalar w916 = -w388 + w390;
    const scalar w917 = -w888 - w889 + w890 + w891 - w916;
    const scalar w918 = -w901;
    const scalar w919 = -w899;
    const scalar w920 = sigma0Out*w155;
    const scalar w921 = w430*w920;
    const scalar w922 = w886*w921;
    const scalar w923 = w122*w428;
    const scalar w924 = w154*w472;
    const scalar w925 = w546*w842;
    const scalar w926 = 2*w465 + 2*w552 + 2*w556 + 2*w559;
    const scalar w927 = w140*w480/(w471*w471*w471);
    const scalar w928 = w560*w927;
    const scalar w929 = w122*w547;
    const scalar w930 = w480*w561;
    const scalar w931 = w154*w930;
    const scalar w932 = -2*w417 - 2*w537 - 2*w541 - 2*w544;
    const scalar w933 = w122*w460/(w427*w427*w427);
    const scalar w934 = w545*w933;
    const scalar w935 = w461*w859 - w481*w861;
    const scalar w936 = w484 + w485;
    const scalar w937 = w7*w773;
    const scalar w938 = w7*w778;
    const scalar w939 = 8*w938;
    const scalar w940 = 8*gil*w15*w38*w7*w776*w777 + 8*git*w15*w35*w7*w776*w777 - w37*w939 - w40*w939 - w781;
    const scalar w941 = 4*w938;
    const scalar w942 = 4*w13*w32*w47*w7*w776*w777 - w183*w941 - w784;
    const scalar w943 = w49*w941 + w67*w941 + w786;
    const scalar w944 = w53*y;
    const scalar w945 = w5*w944;
    const scalar w946 = 4*gil*theta1*w12*w14*w3*w320*y + gil*w13*w14*w942 + git*w12*w14*w943 - w324*w945 - w66*w941 - w70*w941 - w790;
    const scalar w947 = -w48*w941 - w50*w941 - w792;
    const scalar w948 = 4*w13*w47*w7*w776*w777 - w794*w941 - w795;
    const scalar w949 = (2.0/5.0)*w938;
    const scalar w950 = (2.0/5.0)*w4;
    const scalar w951 = w944*w950;
    const scalar w952 = (2.0/5.0)*gil*theta1*w12*w14*w292*w3*y + (1.0/10.0)*gil*w13*w14*w947 + (1.0/10.0)*git*w12*w14*w948 - w256*w946 - w296*w951 - w54*w949 - w59*w949 - w799 + (1.0/5.0)*w940*z;
    const scalar w953 = w292*w306;
    const scalar w954 = gil*w953;
    const scalar w955 = w295*w308;
    const scalar w956 = git*w955;
    const scalar w957 = -w48*w939 - w50*w939 - w807;
    const scalar w958 = w52*w957;
    const scalar w959 = 8*w13*w47*w7*w776*w777 - w794*w939 - w810;
    const scalar w960 = w57*w959;
    const scalar w961 = 4*gil*theta1*w12*w14*w292*w3*y + gil*w13*w14*w947 + git*w12*w14*w948 - w296*w945 - w54*w941 - w59*w941 - w813;
    const scalar w962 = w292*w321;
    const scalar w963 = w295*w324;
    const scalar w964 = w52*w942;
    const scalar w965 = w65*w947;
    const scalar w966 = w57*w943;
    const scalar w967 = w69*w948;
    const scalar w968 = 2*w962 + 2*w963 + w964 + w965 + w966 + w967;
    const scalar w969 = w253*w961 - w256*w968 + (1.0/10.0)*w954 + (1.0/10.0)*w956 + (1.0/10.0)*w958 + (1.0/10.0)*w960;
    const scalar w970 = w320*w337;
    const scalar w971 = w323*w338;
    const scalar w972 = 8*w13*w32*w47*w7*w776*w777 - w183*w939 - w825;
    const scalar w973 = w49*w939 + w67*w939 + w827;
    const scalar w974 = gil*w970 + git*w971 + w65*w972 + w69*w973;
    const scalar w975 = w253*w946 - w256*w974 + (1.0/5.0)*w962 + (1.0/5.0)*w963 + (1.0/10.0)*w964 + (1.0/10.0)*w965 + (1.0/10.0)*w966 + (1.0/10.0)*w967;
    const scalar w976 = w370*w832;
    const scalar w977 = w954 + w956 + w958 + w960;
    const scalar w978 = w153*(-w153*w940 - w313*w961 - 2*w360 - w586*w946) + w313*(-w153*w961 - w310*w775 - w313*w977 - w586*w968) + w365*w775 + w586*(-w153*w946 - w313*w968 - 2*w364 - w586*w974);
    const scalar w979 = w138*w314;
    const scalar w980 = 2*w979;
    const scalar w981 = w312*w97;
    const scalar w982 = 2*w981;
    const scalar w983 = w357*w843;
    const scalar w984 = 2*w134 + 2*w363 + 2*w366 + 2*w369;
    const scalar w985 = w370*w984;
    const scalar w986 = -w149*(w149*w940 + w311*w961 + 2*w347 + w584*w946) - w311*(w149*w961 + w310*w837 + w311*w977 + w584*w968) - w352*w837 - w584*(w149*w946 + w311*w968 + 2*w351 + w584*w974);
    const scalar w987 = 2*w314*w371;
    const scalar w988 = 2*w312;
    const scalar w989 = w357*w988;
    const scalar w990 = -2*w350 - 2*w353 - 2*w356 - 2*w89;
    const scalar w991 = w358*w990;
    const scalar w992 = 4*y;
    const scalar w993 = w30*w385;
    const scalar w994 = w147*w992;
    const scalar w995 = Tt*(w252*(w141*(w153*w952 + w313*w969 + w378*w775 + w586*w975) + w217*(-w149*w952 - w311*w969 - w378*w837 - w584*w975) - w268*w987 + w286*w989 + w381*w982 + w381*w983 + w384*w976 - w384*w980 + w836*w978 + w846*w985 + w850*w986 + w855*w991 + w863) + w270*w937 + w385*w994 + w867 + w992*w993);
    const scalar w996 = w109*w156;
    const scalar w997 = w870*w98;
    const scalar w998 = w102*w997;
    const scalar w999 = (6.0/5.0)*St;
    const scalar w1000 = w316*w999;
    const scalar w1001 = 2*w319;
    const scalar w1002 = w1000*w326 - w1001*w298 - w104*w961 + w111*w946 + w27*w940 + w60*w998 + w71*w996;
    const scalar w1003 = w87*w996;
    const scalar w1004 = w1000*w335 - w1001*w310 + w1003 - w104*w977 + w111*w968 + w27*w961 + w81*w998;
    const scalar w1005 = w1000*w339 - w1001*w335 - w104*w968 + w111*w974 + w27*w946 + w329*w997 + w93*w996;
    const scalar w1006 = w30*w992;
    const scalar w1007 = w155*w435;
    const scalar w1008 = sigma0In*w1007;
    const scalar w1009 = w1008 + w892;
    const scalar w1010 = w400*w939;
    const scalar w1011 = w401*w939;
    const scalar w1012 = gel*w39*w939;
    const scalar w1013 = get*w36*w939;
    const scalar w1014 = -w1010 - w1011 + w1012 + w1013;
    const scalar w1015 = w156*w398;
    const scalar w1016 = w896*y;
    const scalar w1017 = -w1016*w655 + w412*w941 + w414*w941 - w576*w942 - w577*w943 + w656*w945 + w898;
    const scalar w1018 = -w1016*w573 + w406*w941 + w408*w941 + w574*w945 - w576*w947 - w577*w948 + w900;
    const scalar w1019 = w396*w447*w870;
    const scalar w1020 = 2*w662;
    const scalar w1021 = w1017*w443 + w1018*w450 + w1019*w444 - w1020*w659 + w658*w660;
    const scalar w1022 = gel*w953 + get*w955 + w405*w957 + w407*w959;
    const scalar w1023 = 2*w655;
    const scalar w1024 = 2*w656;
    const scalar w1025 = w1023*w292 + w1024*w295 + w405*w942 + w407*w943 + w411*w947 + w413*w948;
    const scalar w1026 = -w1025;
    const scalar w1027 = w501*w659;
    const scalar w1028 = w660*w668;
    const scalar w1029 = -w1015*w444;
    const scalar w1030 = w1018*w399 + w1026*w443 + 2*w1027 + w1028 + w1029;
    const scalar w1031 = gel*w970 + get*w971 + w411*w972 + w413*w973;
    const scalar w1032 = -w1015*w439 + w1017*w399 + w1019*w454 - w1020*w668 + w1026*w450 + w658*w907;
    const scalar w1033 = -w1010 - w1011 + w1012 + w1013 - w916;
    const scalar w1034 = -w1018;
    const scalar w1035 = -w1017;
    const scalar w1036 = sigma0Out*w1007;
    const scalar w1037 = w1036 + w892;
    const scalar w1038 = w314*w472;
    const scalar w1039 = 2*w468 + 2*w702 + 2*w705 + 2*w708;
    const scalar w1040 = w709*w927;
    const scalar w1041 = w314*w930;
    const scalar w1042 = -2*w423 - 2*w689 - 2*w692 - 2*w695;
    const scalar w1043 = w696*w933;
    const scalar w1044 = w6*w773;
    const scalar w1045 = w6*w778;
    const scalar w1046 = 8*w1045;
    const scalar w1047 = 8*gil*w15*w38*w6*w776*w777 + 8*git*w15*w35*w6*w776*w777 - w1046*w37 - w1046*w40 - w781;
    const scalar w1048 = 4*w1045;
    const scalar w1049 = -w1048*w183 + 4*w13*w32*w47*w6*w776*w777 - w784;
    const scalar w1050 = w1048*w49 + w1048*w67 + w786;
    const scalar w1051 = w53*x;
    const scalar w1052 = w1051*w19;
    const scalar w1053 = 4*gil*theta1*w12*w14*w299*w3*x + gil*w1049*w13*w14 + git*w1050*w12*w14 - w1048*w66 - w1048*w70 - w1052*w302 - w790;
    const scalar w1054 = (2.0/5.0)*w635;
    const scalar w1055 = -w1048*w48 - w1048*w50 - w792;
    const scalar w1056 = -w1048*w794 + 4*w13*w47*w6*w776*w777 - w795;
    const scalar w1057 = (2.0/5.0)*w1045;
    const scalar w1058 = w1051*w570;
    const scalar w1059 = (2.0/5.0)*gil*theta1*w12*w14*w3*w567*x + (1.0/10.0)*gil*w1055*w13*w14 + (1.0/10.0)*git*w1056*w12*w14 - git*w1058*w950 + (1.0/5.0)*w1047*z - w1053*w256 - w1054*w305 - w1057*w54 - w1057*w59 + (1.0/5.0)*w254*w71*w869 - w799;
    const scalar w1060 = w255*w869;
    const scalar w1061 = w299*w567;
    const scalar w1062 = gil*w1061;
    const scalar w1063 = w302*w570;
    const scalar w1064 = git*w1063;
    const scalar w1065 = w1049*w52;
    const scalar w1066 = w1055*w65;
    const scalar w1067 = w1050*w57;
    const scalar w1068 = w1056*w69;
    const scalar w1069 = w299*w579;
    const scalar w1070 = w302*w581;
    const scalar w1071 = -w1046*w183 + 8*w13*w32*w47*w6*w776*w777 - w825;
    const scalar w1072 = w1046*w49 + w1046*w67 + w827;
    const scalar w1073 = gil*w1069 + git*w1070 + w1071*w65 + w1072*w69;
    const scalar w1074 = w1053*w253 - w1054*w583 + w1060*w93 + (1.0/5.0)*w1062 + (1.0/5.0)*w1064 + (1.0/10.0)*w1065 + (1.0/10.0)*w1066 + (1.0/10.0)*w1067 + (1.0/10.0)*w1068 - w1073*w256;
    const scalar w1075 = w567*w595;
    const scalar w1076 = gil*w1075;
    const scalar w1077 = w570*w597;
    const scalar w1078 = git*w1077;
    const scalar w1079 = -w1046*w48 - w1046*w50 - w807;
    const scalar w1080 = w1079*w52;
    const scalar w1081 = -w1046*w794 + 8*w13*w47*w6*w776*w777 - w810;
    const scalar w1082 = w1081*w57;
    const scalar w1083 = 4*gil*theta1*w12*w14*w3*w567*x + gil*w1055*w13*w14 + git*w1056*w12*w14 - w1048*w54 - w1048*w59 - w1052*w570 - w813;
    const scalar w1084 = 2*w1062 + 2*w1064 + w1065 + w1066 + w1067 + w1068;
    const scalar w1085 = -w1054*w572 + w1060*w87 + (1.0/10.0)*w1076 + (1.0/10.0)*w1078 + (1.0/10.0)*w1080 + (1.0/10.0)*w1082 + w1083*w253 - w1084*w256;
    const scalar w1086 = w633*w831;
    const scalar w1087 = w1076 + w1078 + w1080 + w1082;
    const scalar w1088 = w153*(-w1047*w153 - w1053*w586 - w1083*w313 - 2*w623) + w313*(-w1083*w153 - w1084*w586 - w1087*w313 - 2*w626) + w586*(-w1053*w153 - w1073*w586 - w1084*w313 - w583*w775) + w627*w775;
    const scalar w1089 = 2*w587;
    const scalar w1090 = w585*w97;
    const scalar w1091 = 2*w1090;
    const scalar w1092 = w619*w842;
    const scalar w1093 = 2*w136 + 2*w625 + 2*w628 + 2*w631;
    const scalar w1094 = w1093*w632;
    const scalar w1095 = -w149*(w1047*w149 + w1053*w584 + w1083*w311 + 2*w609) - w311*(w1083*w149 + w1084*w584 + w1087*w311 + 2*w612) - w584*(w1053*w149 + w1073*w584 + w1084*w311 + w583*w837) - w613*w837;
    const scalar w1096 = w1089*w633;
    const scalar w1097 = 2*w585*w619;
    const scalar w1098 = -2*w611 - 2*w614 - 2*w617 - 2*w95;
    const scalar w1099 = w1098*w855;
    const scalar w1100 = w122*w618;
    const scalar w1101 = 4*x;
    const scalar w1102 = w30*w647;
    const scalar w1103 = w1101*w147;
    const scalar w1104 = Tt*(w1044*w270 + w1101*w1102 + w1103*w647 + w252*(w1086*w645 + w1088*w836 - w1089*w646 + w1091*w642 + w1092*w642 + w1094*w846 + w1095*w850 - w1096*w268 + w1097*w263 + w1099*w1100 + w141*(w1059*w153 + w1074*w586 + w1085*w313 + w639*w775) + w217*(-w1059*w149 - w1074*w584 - w1085*w311 - w639*w837) + w863) - w636 + w867);
    const scalar w1105 = 2*w159;
    const scalar w1106 = (3.0/5.0)*w103;
    const scalar w1107 = w107*w999;
    const scalar w1108 = w1000*w592 - w104*w1083 + w1047*w27 + w1053*w111 + w1105*w588 + w1106*w60 - w1107*w305;
    const scalar w1109 = w1000*w599 - w104*w1087 + w1083*w27 + w1084*w111 + w1105*w592 + w1106*w81 - w1107*w572;
    const scalar w1110 = w1106*w87;
    const scalar w1111 = w159*w305;
    const scalar w1112 = w1000*w572 - w104*w1084 + w1053*w27 + w1073*w111 - w1107*w583 + w1110 + 2*w1111;
    const scalar w1113 = w138*w587;
    const scalar w1114 = w1100*w880;
    const scalar w1115 = w1101*w30;
    const scalar w1116 = w1046*w400;
    const scalar w1117 = w1046*w401;
    const scalar w1118 = gel*w1046*w39;
    const scalar w1119 = get*w1046*w36;
    const scalar w1120 = -w1116 - w1117 + w1118 + w1119;
    const scalar w1121 = gel*w896*x;
    const scalar w1122 = get*w5;
    const scalar w1123 = w1048*w412 + w1048*w414 - w1049*w576 - w1050*w577 + w1051*w1122*w302 - w1121*w299 + w898;
    const scalar w1124 = w1048*w406 + w1048*w408 - w1055*w576 - w1056*w577 + w1058*w1122 - w1121*w567 + w900;
    const scalar w1125 = w156*w446;
    const scalar w1126 = fA2*w440*w441*w870;
    const scalar w1127 = 2*w722;
    const scalar w1128 = w1123*w443 + w1124*w450 - w1125*w444 - w1126*w439 - w1127*w719 + w660*w721;
    const scalar w1129 = gel*w1075 + get*w1077 + w1079*w405 + w1081*w407;
    const scalar w1130 = 2*gel;
    const scalar w1131 = 2*get;
    const scalar w1132 = w1049*w405 + w1050*w407 + w1055*w411 + w1056*w413 + w1061*w1130 + w1063*w1131;
    const scalar w1133 = -w1132;
    const scalar w1134 = w1124*w399 - w1126*w454 - w1127*w729 + w1133*w443;
    const scalar w1135 = gel*w1069 + get*w1070 + w1071*w411 + w1072*w413;
    const scalar w1136 = w660*w729;
    const scalar w1137 = -w1125*w454;
    const scalar w1138 = w1123*w399 + w1133*w450 + w1136 + w1137;
    const scalar w1139 = -w1116 - w1117 + w1118 + w1119 - w916;
    const scalar w1140 = -w1124;
    const scalar w1141 = -w1123;
    const scalar w1142 = w472*w587;
    const scalar w1143 = w428*w585;
    const scalar w1144 = 2*w470 + 2*w763 + 2*w766 + 2*w769;
    const scalar w1145 = w770*w927;
    const scalar w1146 = w770*w930;
    const scalar w1147 = w547*w585;
    const scalar w1148 = -2*w426 - 2*w750 - 2*w753 - 2*w756;
    const scalar w1149 = w757*w933;
    const scalar w1150 = y*z;
    const scalar w1151 = 8*w1150*w778;
    const scalar w1152 = -w1151*w183 + w1151*w63;
    const scalar w1153 = w1151*w49 + w1151*w67;
    const scalar w1154 = w1152*w65 + w1153*w69;
    const scalar w1155 = w1154 + w181*w337 + w187*w338;
    const scalar w1156 = w778*z;
    const scalar w1157 = w1156*w992;
    const scalar w1158 = -w1157*w183 + w1157*w63;
    const scalar w1159 = w1157*w49 + w1157*w67;
    const scalar w1160 = w164*y;
    const scalar w1161 = w164*z;
    const scalar w1162 = 2*gil*theta1*w12*w14*w180*w3*y + 2*gil*theta1*w12*w14*w3*w320*z + gil*w1158*w13*w14 + git*w1159*w12*w14 - w1157*w66 - w1157*w70 - w1160*w187 - w1161*w324;
    const scalar w1163 = -w1157*w48 - w1157*w50;
    const scalar w1164 = w1163*w65;
    const scalar w1165 = w1158*w52;
    const scalar w1166 = w167*w321;
    const scalar w1167 = w180*w293;
    const scalar w1168 = w1157*w55 - w1157*w794;
    const scalar w1169 = w1168*w69;
    const scalar w1170 = w1159*w57;
    const scalar w1171 = w172*w324;
    const scalar w1172 = w186*w296;
    const scalar w1173 = w58*w950;
    const scalar w1174 = w1162*w253 + (1.0/10.0)*w1164 + (1.0/10.0)*w1165 + (1.0/10.0)*w1166 + (1.0/10.0)*w1167 + (1.0/10.0)*w1169 + (1.0/10.0)*w1170 + (1.0/10.0)*w1171 + (1.0/10.0)*w1172 + w1173*w65*y + (1.0/5.0)*w322 + (1.0/5.0)*w325 - w69*w951;
    const scalar w1175 = -w1155*w256 + w1174;
    const scalar w1176 = w168*w306;
    const scalar w1177 = w173*w308;
    const scalar w1178 = 2*gil*theta1*w12*w14*w167*w3*y + 2*gil*theta1*w12*w14*w292*w3*z + gil*w1163*w13*w14 + git*w1168*w12*w14 - w1157*w54 - w1157*w59 - w1160*w173 - w1161*w296;
    const scalar w1179 = w1164 + w1165 + w1166 + w1167 + w1169 + w1170 + w1171 + w1172;
    const scalar w1180 = -w1151*w48 - w1151*w50;
    const scalar w1181 = w1180*w52;
    const scalar w1182 = w1151*w55 - w1151*w794;
    const scalar w1183 = w1182*w57;
    const scalar w1184 = w950*y;
    const scalar w1185 = w1178*w253 - w1179*w256 + (1.0/10.0)*w1181 + (1.0/10.0)*w1183 + w1184*w274 - w1184*w276 + (1.0/5.0)*w294 + (1.0/5.0)*w297;
    const scalar w1186 = (1.0/10.0)*w1176 + (1.0/10.0)*w1177 + w1185;
    const scalar w1187 = git*y;
    const scalar w1188 = w16*w801;
    const scalar w1189 = w778*y;
    const scalar w1190 = w1189*w395;
    const scalar w1191 = w275*w53;
    const scalar w1192 = 8*gil*w15*w38*w776*w777*y*z + 8*git*w15*w35*w776*w777*y*z - w1151*w37 - w1151*w40;
    const scalar w1193 = (4.0/5.0)*gil*theta1*w12*w13*w15*w3*y + (1.0/5.0)*gil*theta1*w12*w14*w167*w3*y + (1.0/5.0)*gil*theta1*w12*w14*w292*w3*z + (1.0/10.0)*gil*w1163*w13*w14 + (1.0/10.0)*git*w1168*w12*w14 - w1162*w256 - w1187*w1188 - w1190*w54 - w1190*w59 - w1191*w296 + (1.0/5.0)*w1192*z - w173*w373*w944;
    const scalar w1194 = w1193*w153 + w283 + w382;
    const scalar w1195 = w1162*w153 + w1179*w313 + w243 + w361;
    const scalar w1196 = w1181 + w1183;
    const scalar w1197 = w1176 + w1177 + w1196;
    const scalar w1198 = w1178*w153 + w1179*w586 + w242 + w360;
    const scalar w1199 = w153*(-w1162*w586 - w1178*w313 - w1192*w153 - w238 - w359) + w244 + w362;
    const scalar w1200 = w1199 + w313*(-w1197*w313 - w1198) + w586*(-w1155*w586 - w1195);
    const scalar w1201 = w1193*w149 + w279 + w377;
    const scalar w1202 = w846*w847;
    const scalar w1203 = w1162*w149 + w1179*w311 + w226 + w348;
    const scalar w1204 = w1178*w149 + w1179*w584 + w225 + w347;
    const scalar w1205 = w149*(w1162*w584 + w1178*w311 + w1192*w149 + w221 + w346) + w227 + w349;
    const scalar w1206 = -w1205 - w311*(w1197*w311 + w1204) - w584*(w1155*w584 + w1203);
    const scalar w1207 = w358*w856;
    const scalar w1208 = w140*w371;
    const scalar w1209 = w140*w250;
    const scalar w1210 = w151*w357;
    const scalar w1211 = w233*w358;
    const scalar w1212 = w233*w857;
    const scalar w1213 = w154*w268;
    const scalar w1214 = w1208*w285 + w1209*w384 + w1210*w286 + w1211*w282 + w1212*w381 - w1213*w371 + w282*w981 - w285*w979 + w287*w312 - w314*w851 + w381*w840 - w384*w838;
    const scalar w1215 = 2*y;
    const scalar w1216 = w1215*w147;
    const scalar w1217 = 2*z;
    const scalar w1218 = w1217*w147;
    const scalar w1219 = w1150*w773;
    const scalar w1220 = w1215*w865 + w1216*w289 + w1217*w993 + w1218*w385 + w1219*w270;
    const scalar w1221 = Tt*(w1220 + w252*(w1200*w836 + w1202*w370 + w1206*w850 + w1207*w855 + w1214 + w141*(w1175*w586 + w1186*w313 + w1194) + w217*(-w1175*w584 - w1186*w311 - w1201)));
    const scalar w1222 = -w104*w1179 + w1162*w27 - w158*w326 + w159*w339 + w198*w318 - w203*w319;
    const scalar w1223 = w111*w1155 + w1222;
    const scalar w1224 = w203*w318;
    const scalar w1225 = w111*w1179 + w1178*w27 + w1224 - w158*w298 + w159*w335 - w213*w319;
    const scalar w1226 = -w104*w1197 + w1225;
    const scalar w1227 = w159*w326;
    const scalar w1228 = -w104*w1178 + w111*w1162 + w1192*w27 + w1227 - w158*w315 - w175*w319 + w189*w318;
    const scalar w1229 = w1228*w153 + w219 + w343;
    const scalar w1230 = w1228*w149 + w215 + w328;
    const scalar w1231 = w847*w877;
    const scalar w1232 = w127*w154;
    const scalar w1233 = w1208*w220 + w1209*w345 + w1210*w234 + w1211*w216 + w1212*w342 - w1232*w371 + w216*w981 - w220*w979 + w235*w312 - w314*w879 + w342*w840 - w345*w838;
    const scalar w1234 = w1215*w30;
    const scalar w1235 = w1217*w30;
    const scalar w1236 = w1216*w251 + w1218*w372 + w1219*w142 + w1234*w251 + w1235*w372;
    const scalar w1237 = w1154 + w194*w321 + w197*w324;
    const scalar w1238 = w1174 - w1237*w256;
    const scalar w1239 = w207*w293;
    const scalar w1240 = w211*w296;
    const scalar w1241 = w1185 + (1.0/10.0)*w1239 + (1.0/10.0)*w1240;
    const scalar w1242 = w1196 + w1239 + w1240;
    const scalar w1243 = w1199 + w313*(-w1198 - w1242*w313) + w586*(-w1195 - w1237*w586);
    const scalar w1244 = w249*w846;
    const scalar w1245 = -w1205 - w311*(w1204 + w1242*w311) - w584*(w1203 + w1237*w584);
    const scalar w1246 = w857*w990;
    const scalar w1247 = Tt*(w1220 + w252*(w1214 + w1243*w836 + w1244*w984 + w1245*w850 + w1246*w855 + w141*(w1194 + w1238*w586 + w1241*w313) + w217*(-w1201 - w1238*w584 - w1241*w311)));
    const scalar w1248 = w111*w1237 + w1222;
    const scalar w1249 = -w104*w1242 + w1225;
    const scalar w1250 = w249*w877;
    const scalar w1251 = 8*x;
    const scalar w1252 = w1156*w1251;
    const scalar w1253 = -w1252*w183 + w1252*w63;
    const scalar w1254 = w1252*w49 + w1252*w67;
    const scalar w1255 = w1253*w65 + w1254*w69;
    const scalar w1256 = w1255 + w181*w579 + w187*w581;
    const scalar w1257 = w1101*w1156;
    const scalar w1258 = -w1257*w183 + w1257*w63;
    const scalar w1259 = w1257*w49 + w1257*w67;
    const scalar w1260 = w164*x;
    const scalar w1261 = w164*w302;
    const scalar w1262 = 2*gil*theta1*w12*w14*w180*w3*x + 2*gil*theta1*w12*w14*w299*w3*z + gil*w1258*w13*w14 + git*w12*w1259*w14 - w1257*w66 - w1257*w70 - w1260*w187 - w1261*w798;
    const scalar w1263 = -w1257*w48 - w1257*w50;
    const scalar w1264 = w1263*w65;
    const scalar w1265 = w1258*w52;
    const scalar w1266 = w168*w299;
    const scalar w1267 = w181*w567;
    const scalar w1268 = w1257*w55 - w1257*w794;
    const scalar w1269 = w1268*w69;
    const scalar w1270 = w1259*w57;
    const scalar w1271 = w173*w302;
    const scalar w1272 = w187*w570;
    const scalar w1273 = -w1051*w69*w950 + w1173*w65*x + w1262*w253 + (1.0/10.0)*w1264 + (1.0/10.0)*w1265 + (1.0/10.0)*w1266 + (1.0/10.0)*w1267 + (1.0/10.0)*w1269 + (1.0/10.0)*w1270 + (1.0/10.0)*w1271 + (1.0/10.0)*w1272 - w198*w636 + (1.0/5.0)*w301 + (1.0/5.0)*w304;
    const scalar w1274 = -w1256*w256 + w1273;
    const scalar w1275 = w168*w595;
    const scalar w1276 = w173*w597;
    const scalar w1277 = w164*w570;
    const scalar w1278 = 2*gil*theta1*w12*w14*w167*w3*x + 2*gil*theta1*w12*w14*w3*w567*z + gil*w1263*w13*w14 + git*w12*w1268*w14 - w1257*w54 - w1257*w59 - w1260*w173 - w1277*w798;
    const scalar w1279 = w1264 + w1265 + w1266 + w1267 + w1269 + w1270 + w1271 + w1272;
    const scalar w1280 = -w1252*w48 - w1252*w50;
    const scalar w1281 = w1280*w52;
    const scalar w1282 = w1252*w55 - w1252*w794;
    const scalar w1283 = w1282*w57;
    const scalar w1284 = w950*x;
    const scalar w1285 = w1278*w253 - w1279*w256 + (1.0/10.0)*w1281 + (1.0/10.0)*w1283 + w1284*w274 - w1284*w276 - w203*w636 + (1.0/5.0)*w590 + (1.0/5.0)*w591;
    const scalar w1286 = (1.0/10.0)*w1275 + (1.0/10.0)*w1276 + w1285;
    const scalar w1287 = w54*x;
    const scalar w1288 = w395*w778;
    const scalar w1289 = w59*x;
    const scalar w1290 = w1051*w373;
    const scalar w1291 = 8*gil*w15*w38*w776*w777*x*z + 8*git*w15*w35*w776*w777*x*z - w1252*w37 - w1252*w40;
    const scalar w1292 = (4.0/5.0)*gil*theta1*w12*w13*w15*w3*x + (1.0/5.0)*gil*theta1*w12*w14*w167*w3*x + (1.0/5.0)*gil*theta1*w12*w14*w3*w567*z + (1.0/10.0)*gil*w1263*w13*w14 - git*w1188*x - git*w1191*w570 + (1.0/10.0)*git*w12*w1268*w14 - w1262*w256 - w1287*w1288 - w1288*w1289 - w1290*w173 + (1.0/5.0)*w1291*z - w189*w636;
    const scalar w1293 = w1292*w153 + w284 + w643;
    const scalar w1294 = w1262*w153 + w1279*w313 + w246 + w623;
    const scalar w1295 = w1281 + w1283;
    const scalar w1296 = w1275 + w1276 + w1295;
    const scalar w1297 = w1278*w153 + w1279*w586 + w243 + w622;
    const scalar w1298 = w153*(-w1262*w586 - w1278*w313 - w1291*w153 - w239 - w621) + w247 + w624;
    const scalar w1299 = w1298 + w313*(-w1296*w313 - w1297) + w586*(-w1256*w586 - w1294);
    const scalar w1300 = w1292*w149 + w281 + w638;
    const scalar w1301 = w1262*w149 + w1279*w311 + w229 + w609;
    const scalar w1302 = w1278*w149 + w1279*w584 + w226 + w608;
    const scalar w1303 = w149*(w1262*w584 + w1278*w311 + w1291*w149 + w222 + w607) + w230 + w610;
    const scalar w1304 = -w1303 - w311*(w1296*w311 + w1302) - w584*(w1256*w584 + w1301);
    const scalar w1305 = w1100*w855;
    const scalar w1306 = w140*w633;
    const scalar w1307 = w151*w619;
    const scalar w1308 = w1090*w282 - w1113*w285 + w1209*w645 + w1212*w642 - w1213*w633 + w1306*w285 + w1307*w263 - w154*w646 + w282*w620 + w287*w585 - w587*w851 + w642*w840;
    const scalar w1309 = 2*x;
    const scalar w1310 = w1309*w147;
    const scalar w1311 = w773*x;
    const scalar w1312 = w1311*z;
    const scalar w1313 = w1102*w1217 + w1218*w647 + w1309*w865 + w1310*w289 + w1312*w270;
    const scalar w1314 = Tt*(w1313 + w252*(w1202*w632 + w1299*w836 + w1304*w850 + w1305*w856 + w1308 + w141*(w1274*w586 + w1286*w313 + w1293) + w217*(-w1274*w584 - w1286*w311 - w1300)));
    const scalar w1315 = -w104*w1279 + w1224 + w1262*w27 - w158*w305 + w159*w189 + w159*w583 - w191 - w198*w589;
    const scalar w1316 = w111*w1256 + w1315;
    const scalar w1317 = w111*w1279 + w1278*w27 - w157*w60 - w158*w592 + w159*w175 + w159*w572 - w203*w589 + w213*w318;
    const scalar w1318 = -w104*w1296 + w1317;
    const scalar w1319 = -w104*w1278 + w111*w1262 + w1111 + w1291*w27 - w157*w41 - w158*w588 + w159*w21 + w175*w318 - w189*w589;
    const scalar w1320 = w1319*w153 + w218 + w604;
    const scalar w1321 = w1319*w149 + w205 + w594;
    const scalar w1322 = w1090*w216 - w1113*w220 + w118*w1307 + w1209*w606 + w1212*w603 - w1232*w633 + w1306*w220 + w216*w620 + w235*w585 - w587*w879 + w603*w840 - w606*w838;
    const scalar w1323 = w1309*w30;
    const scalar w1324 = w1218*w634 + w1235*w634 + w1310*w251 + w1312*w142 + w1323*w251 + w159;
    const scalar w1325 = w194*w299;
    const scalar w1326 = w197*w302;
    const scalar w1327 = gil*w1325 + git*w1326 + w1255;
    const scalar w1328 = w1273 - w1327*w256;
    const scalar w1329 = w207*w567;
    const scalar w1330 = gil*w1329;
    const scalar w1331 = w211*w570;
    const scalar w1332 = git*w1331;
    const scalar w1333 = w1285 + (1.0/10.0)*w1330 + (1.0/10.0)*w1332;
    const scalar w1334 = w1295 + w1330 + w1332;
    const scalar w1335 = w1298 + w313*(-w1297 - w1334*w313) + w586*(-w1294 - w1327*w586);
    const scalar w1336 = -w1303 - w311*(w1302 + w1334*w311) - w584*(w1301 + w1327*w584);
    const scalar w1337 = Tt*(w1313 + w252*(w1093*w1244 + w1099*w857 + w1308 + w1335*w836 + w1336*w850 + w141*(w1293 + w1328*w586 + w1333*w313) + w217*(-w1300 - w1328*w584 - w1333*w311)));
    const scalar w1338 = w111*w1327 + w1315;
    const scalar w1339 = -w104*w1334 + w1317;
    const scalar w1340 = w1098*w880;
    const scalar w1341 = -1.0/2.0*sigma0In*sigmaEps*w106*w155*w24*w713;
    const scalar w1342 = w1253*w411 + w1254*w413;
    const scalar w1343 = w1342 + w497*w579 + w498*w581;
    const scalar w1344 = w4*w58;
    const scalar w1345 = w1344*x;
    const scalar w1346 = w1130*w1344;
    const scalar w1347 = w1346*z;
    const scalar w1348 = w1051*w4;
    const scalar w1349 = w1131*w4;
    const scalar w1350 = w1349*w302;
    const scalar w1351 = w1257*w412 + w1257*w414 - w1258*w576 - w1259*w577 - w1345*w908 - w1347*w299 + w1348*w909 + w1350*w788;
    const scalar w1352 = w1258*w405 + w1259*w407 + w1263*w411 + w1268*w413 + w299*w504 + w302*w505 + w497*w567 + w498*w570;
    const scalar w1353 = -w1352;
    const scalar w1354 = w514*w661;
    const scalar w1355 = w1351*w399 + w1353*w450 + w1354 + w495*w719 + w501*w729;
    const scalar w1356 = w1355 - w530*w722;
    const scalar w1357 = w1280*w405 + w1282*w407;
    const scalar w1358 = w1357 + w504*w595 + w505*w597;
    const scalar w1359 = w162*x;
    const scalar w1360 = w1257*w406 + w1257*w408 + w1260*w505 - w1263*w576 - w1268*w577 - w1347*w567 + w1349*w570*w788 - w1359*w504;
    const scalar w1361 = w1353*w443 + w1360*w399 + w495*w721 - w514*w722;
    const scalar w1362 = w1361 + w501*w740 + w528*w661;
    const scalar w1363 = w490*w713;
    const scalar w1364 = -gel*w1252*w39 - get*w1252*w36 + w1252*w400 + w1252*w401;
    const scalar w1365 = w501*w721;
    const scalar w1366 = w1351*w443 + w1360*w450 + w1365 - w500*w722 + w507*w661;
    const scalar w1367 = w153*(w1366 + w399*(w1363*w884 + w1364) + w495*w737) + w531 + w738;
    const scalar w1368 = -w1351;
    const scalar w1369 = w1352*w313 + w1368*w153 + w557 + w761;
    const scalar w1370 = -w1360;
    const scalar w1371 = w1352*w586 + w1370*w153 + w554 + w760;
    const scalar w1372 = -w1364;
    const scalar w1373 = w153*(-w1368*w586 - w1370*w313 - w1372*w153 - w550 - w759) + w558 + w762;
    const scalar w1374 = -1.0/2.0*sigma0Out*sigmaEps*w106*w155*w24*w713;
    const scalar w1375 = w1355 - w518*w722;
    const scalar w1376 = w1361 + w501*w727 + w512*w661;
    const scalar w1377 = w149*(w1366 + w399*(w1363*w920 + w1364) + w495*w718) + w520 + w724;
    const scalar w1378 = w1352*w311 + w1368*w149 + w542 + w748;
    const scalar w1379 = w1352*w584 + w1370*w149 + w539 + w747;
    const scalar w1380 = w149*(w1368*w584 + w1370*w311 + w1372*w149 + w535 + w746) + w543 + w749;
    const scalar w1381 = w562*w770;
    const scalar w1382 = w562*w744;
    const scalar w1383 = w472*w744;
    const scalar w1384 = w122*w546;
    const scalar w1385 = w1384*w757;
    const scalar w1386 = w1384*w735;
    const scalar w1387 = w587*w930;
    const scalar w1388 = -w1142*w532 + w1143*w521 - w1146*w154 + w1147*w545 + w1381*w532 + w1382*w560 - w1383*w154 + w1385*w521 + w1386*w545 - w1387*w560 + w151*w736 + w151*w758;
    const scalar w1389 = w1218*w771 + w1235*w771 + w1310*w564 + w1312*w482 + w1323*w564;
    const scalar w1390 = gel*w1325 + get*w1326 + w1342;
    const scalar w1391 = gel*w1329 + get*w1331 + w1357;
    const scalar w1392 = sigmaEps*w101*w489;
    const scalar w1393 = w1392*w885;
    const scalar w1394 = w1152*w411 + w1153*w413;
    const scalar w1395 = w1394 + w337*w497 + w338*w498;
    const scalar w1396 = -w1023*w1344*z + w1024*w4*w788 + w1157*w412 + w1157*w414 - w1158*w576 - w1159*w577 - w1344*w908*y + w4*w909*w944;
    const scalar w1397 = w1158*w405 + w1159*w407 + w1163*w411 + w1168*w413 + w167*w655 + w172*w656 + w180*w573 + w186*w574;
    const scalar w1398 = -w1397;
    const scalar w1399 = -w1015*w454 + w1396*w399 + w1398*w450 - w439*w903 + w495*w658 + w500*w501 + w501*w668 - w514*w662;
    const scalar w1400 = w1399 + w530*w661;
    const scalar w1401 = w1180*w405 + w1182*w407;
    const scalar w1402 = w1401 + w306*w504 + w308*w505;
    const scalar w1403 = w1157*w406 + w1157*w408 + w1160*w505 + w1161*w574 - w1163*w576 - w1168*w577 - w162*w504*y - w162*w573*z;
    const scalar w1404 = w1354 + w1398*w443 + w1403*w399 + w495*w659 + w902 + w904;
    const scalar w1405 = -w1015*w476 + w1404 + w501*w680 - w528*w662;
    const scalar w1406 = w1151*w400;
    const scalar w1407 = w1151*w401;
    const scalar w1408 = gel*w1151*w39;
    const scalar w1409 = get*w1151*w36;
    const scalar w1410 = -w1406 - w1407 + w1408 + w1409;
    const scalar w1411 = w1027 + w1029 + w1396*w443 + w1403*w450 + w500*w661 - w507*w662;
    const scalar w1412 = w153*(w1411 + w399*(-w1393 - w1410) - w474*w903 + w495*w678 + w501*w525) + w529 + w679;
    const scalar w1413 = -w1396;
    const scalar w1414 = w1397*w313 + w1413*w153 + w554 + w700;
    const scalar w1415 = -w1403;
    const scalar w1416 = w1397*w586 + w1415*w153 + w553 + w699;
    const scalar w1417 = -w1406 - w1407 + w1408 + w1409;
    const scalar w1418 = w153*(-w1413*w586 - w1415*w313 - w1417*w153 - w549 - w698) + w555 + w701;
    const scalar w1419 = w1392*w921;
    const scalar w1420 = w1399 + w518*w661;
    const scalar w1421 = -w1015*w453 + w1404 + w501*w665 - w512*w662;
    const scalar w1422 = w149*(w1411 + w399*(-w1410 - w1419) - w438*w903 + w492*w501 + w495*w654) + w516 + w664;
    const scalar w1423 = w1397*w311 + w1413*w149 + w539 + w687;
    const scalar w1424 = w1397*w584 + w1415*w149 + w538 + w686;
    const scalar w1425 = w149*(w1413*w584 + w1415*w311 + w1417*w149 + w534 + w685) + w540 + w688;
    const scalar w1426 = -w1038*w532 - w1041*w560 + w1384*w521*w696 + w1384*w545*w675 + w151*w676 + w151*w697 + w312*w522 + w312*w548 + w532*w562*w709 + w560*w562*w684 - w684*w924 - w709*w931;
    const scalar w1427 = w1216*w564 + w1218*w710 + w1219*w482 + w1234*w564 + w1235*w710 + w501;
    const scalar w1428 = w1394 + w194*w655 + w197*w656;
    const scalar w1429 = w1401 + w207*w573 + w211*w574;
    const scalar w1430 = w1189*w1251;
    const scalar w1431 = -w1430*w183 + w1430*w63;
    const scalar w1432 = w1430*w49 + w1430*w67;
    const scalar w1433 = w1431*w65 + w1432*w69;
    const scalar w1434 = w1433 + w321*w579 + w324*w581;
    const scalar w1435 = w1101*w1189;
    const scalar w1436 = -w1435*w183 + w1435*w63;
    const scalar w1437 = w1435*w49 + w1435*w67;
    const scalar w1438 = 2*gil*theta1*w12*w14*w299*w3*y + 2*gil*theta1*w12*w14*w3*w320*x + gil*w13*w14*w1436 + git*w12*w14*w1437 - w1187*w1261 - w1260*w324 - w1435*w66 - w1435*w70;
    const scalar w1439 = -w1435*w48 - w1435*w50;
    const scalar w1440 = w1439*w65;
    const scalar w1441 = w1436*w52;
    const scalar w1442 = w293*w299;
    const scalar w1443 = w321*w567;
    const scalar w1444 = w1435*w55 - w1435*w794;
    const scalar w1445 = w1444*w69;
    const scalar w1446 = w1437*w57;
    const scalar w1447 = w296*w302;
    const scalar w1448 = w324*w570;
    const scalar w1449 = w1438*w253 + (1.0/10.0)*w1440 + (1.0/10.0)*w1441 + (1.0/10.0)*w1442 + (1.0/10.0)*w1443 + (1.0/10.0)*w1445 + (1.0/10.0)*w1446 + (1.0/10.0)*w1447 + (1.0/10.0)*w1448 - w339*w636;
    const scalar w1450 = -w1434*w256 + w1449;
    const scalar w1451 = w293*w595;
    const scalar w1452 = w296*w597;
    const scalar w1453 = 2*gil*theta1*w12*w14*w292*w3*x + 2*gil*theta1*w12*w14*w3*w567*y + gil*w13*w14*w1439 + git*w12*w14*w1444 - w1187*w1277 - w1260*w296 - w1435*w54 - w1435*w59;
    const scalar w1454 = w1440 + w1441 + w1442 + w1443 + w1445 + w1446 + w1447 + w1448;
    const scalar w1455 = -w1430*w48 - w1430*w50;
    const scalar w1456 = w1455*w52;
    const scalar w1457 = w1430*w55 - w1430*w794;
    const scalar w1458 = w1457*w57;
    const scalar w1459 = w1453*w253 - w1454*w256 + (1.0/10.0)*w1456 + (1.0/10.0)*w1458 - w335*w636;
    const scalar w1460 = (1.0/10.0)*w1451 + (1.0/10.0)*w1452 + w1459;
    const scalar w1461 = (2.0/5.0)*w1189;
    const scalar w1462 = w570*w944;
    const scalar w1463 = 8*gil*w15*w38*w776*w777*x*y + 8*git*w15*w35*w776*w777*x*y - w1430*w37 - w1430*w40;
    const scalar w1464 = (1.0/5.0)*gil*theta1*w12*w14*w292*w3*x + (1.0/5.0)*gil*theta1*w12*w14*w3*w567*y + (1.0/10.0)*gil*w13*w14*w1439 + (1.0/10.0)*git*w12*w14*w1444 - git*w1462*w373 - w1287*w1461 - w1289*w1461 - w1290*w296 - w1438*w256 + (1.0/5.0)*w1463*z - w326*w636;
    const scalar w1465 = w1464*w153 + w383 + w644;
    const scalar w1466 = w1438*w153 + w1454*w313 + w367 + w626;
    const scalar w1467 = w1456 + w1458;
    const scalar w1468 = w1451 + w1452 + w1467;
    const scalar w1469 = w1453*w153 + w1454*w586 + w364 + w629;
    const scalar w1470 = w153*(-w1438*w586 - w1453*w313 - w1463*w153 - w361 - w622) + w368 + w630;
    const scalar w1471 = w1470 + w313*(-w1468*w313 - w1469) + w586*(-w1434*w586 - w1466);
    const scalar w1472 = w1464*w149 + w380 + w641;
    const scalar w1473 = w632*w984;
    const scalar w1474 = w1438*w149 + w1454*w311 + w354 + w612;
    const scalar w1475 = w1453*w149 + w1454*w584 + w351 + w615;
    const scalar w1476 = w149*(w1438*w584 + w1453*w311 + w1463*w149 + w348 + w608) + w355 + w616;
    const scalar w1477 = -w1476 - w311*(w1468*w311 + w1475) - w584*(w1434*w584 + w1474);
    const scalar w1478 = w357*w585;
    const scalar w1479 = w312*w619;
    const scalar w1480 = w371*w587;
    const scalar w1481 = w314*w633;
    const scalar w1482 = w1090*w381 - w1113*w384 + w1208*w645 + w1211*w642 + w1306*w384 + w1478*w286 + w1479*w263 - w1480*w268 - w1481*w268 - w314*w646 + w381*w620 + w642*w981;
    const scalar w1483 = w1311*y;
    const scalar w1484 = w1102*w1215 + w1216*w647 + w1309*w993 + w1310*w385 + w1483*w270;
    const scalar w1485 = Tt*(w1484 + w252*(w1305*w990 + w141*(w1450*w586 + w1460*w313 + w1465) + w1471*w836 + w1473*w846 + w1477*w850 + w1482 + w217*(-w1450*w584 - w1460*w311 - w1472)));
    const scalar w1486 = w1003 - w104*w1454 + w1106*w93 + w1227 + w1438*w27 + w318*w335 + w318*w583 - w319*w572 - w339*w589;
    const scalar w1487 = w111*w1434 + w1486;
    const scalar w1488 = w111*w1454 + w1110 + w1453*w27 + w159*w298 + w310*w318 + w318*w572 - w319*w599 - w335*w589 + w81*w996;
    const scalar w1489 = -w104*w1468 + w1488;
    const scalar w1490 = -w104*w1453 + w1106*w71 + w111*w1438 + w1463*w27 + w159*w315 + w298*w318 + w305*w318 - w319*w592 - w326*w589 + w60*w996;
    const scalar w1491 = w1490*w153 + w344 + w605;
    const scalar w1492 = w149*w1490 + w341 + w601;
    const scalar w1493 = w1090*w342 - w1113*w345 + w118*w1479 + w1208*w606 + w1211*w603 - w127*w1480 - w127*w1481 + w1306*w345 + w1478*w234 + w342*w620 + w603*w981 - w606*w979;
    const scalar w1494 = w1216*w634 + w1234*w634 + w1310*w372 + w1323*w372 + w142*w1483 + w318;
    const scalar w1495 = w299*w337;
    const scalar w1496 = w302*w338;
    const scalar w1497 = gil*w1495 + git*w1496 + w1433;
    const scalar w1498 = w1449 - w1497*w256;
    const scalar w1499 = w306*w567;
    const scalar w1500 = gil*w1499;
    const scalar w1501 = w308*w570;
    const scalar w1502 = git*w1501;
    const scalar w1503 = w1459 + (1.0/10.0)*w1500 + (1.0/10.0)*w1502;
    const scalar w1504 = w1467 + w1500 + w1502;
    const scalar w1505 = w1470 + w313*(-w1469 - w1504*w313) + w586*(-w1466 - w1497*w586);
    const scalar w1506 = w1093*w370;
    const scalar w1507 = -w1476 - w311*(w1475 + w1504*w311) - w584*(w1474 + w1497*w584);
    const scalar w1508 = Tt*(w1484 + w252*(w1099*w358 + w141*(w1465 + w1498*w586 + w1503*w313) + w1482 + w1505*w836 + w1506*w846 + w1507*w850 + w217*(-w1472 - w1498*w584 - w1503*w311)));
    const scalar w1509 = w111*w1497 + w1486;
    const scalar w1510 = -w104*w1504 + w1488;
    const scalar w1511 = w651*w713;
    const scalar w1512 = w1511*w884;
    const scalar w1513 = w1431*w411 + w1432*w413;
    const scalar w1514 = w1513 + w579*w655 + w581*w656;
    const scalar w1515 = w1346*y;
    const scalar w1516 = -w1023*w1345 + w1024*w1348 + w1350*w944 + w1435*w412 + w1435*w414 - w1436*w576 - w1437*w577 - w1515*w299;
    const scalar w1517 = w1436*w405 + w1437*w407 + w1439*w411 + w1444*w413 + w299*w573 + w302*w574 + w567*w655 + w570*w656;
    const scalar w1518 = -w1517;
    const scalar w1519 = (1.0/2.0)*w1028 + w1516*w399 + w1518*w450 + w501*w719 - w662*w729;
    const scalar w1520 = -w1125*w478 + w1519 + w661*w742 - w682*w722;
    const scalar w1521 = w1455*w405 + w1457*w407;
    const scalar w1522 = w1521 + w573*w595 + w574*w597;
    const scalar w1523 = w1260*w574 + w1349*w1462 - w1359*w573 + w1435*w406 + w1435*w408 - w1439*w576 - w1444*w577 - w1515*w567;
    const scalar w1524 = (1.0/2.0)*w1136 + w1137 + w1365 + w1518*w443 + w1523*w399 - w668*w722;
    const scalar w1525 = w1524 + w661*w680 - w662*w740;
    const scalar w1526 = w1430*w400;
    const scalar w1527 = w1430*w401;
    const scalar w1528 = gel*w1430*w39;
    const scalar w1529 = get*w1430*w36;
    const scalar w1530 = -w1526 - w1527 + w1528 + w1529;
    const scalar w1531 = -w1125*w439 + w1516*w443 + w1523*w450 - w658*w722 + w659*w661 + w661*w719 - w662*w721;
    const scalar w1532 = w153*(w1531 + w399*(-w1512 - w1530) + w501*w737) + w683 + w741;
    const scalar w1533 = -w1516;
    const scalar w1534 = w1517*w313 + w153*w1533 + w706 + w764;
    const scalar w1535 = -w1523;
    const scalar w1536 = w1517*w586 + w153*w1535 + w703 + w767;
    const scalar w1537 = -w1526 - w1527 + w1528 + w1529;
    const scalar w1538 = w153*(-w153*w1537 - w1533*w586 - w1535*w313 - w700 - w760) + w707 + w768;
    const scalar w1539 = w1511*w920;
    const scalar w1540 = -w1125*w457 + w1519 + w661*w732 - w672*w722;
    const scalar w1541 = w1524 + w661*w665 - w662*w727;
    const scalar w1542 = w149*(w1531 + w399*(-w1530 - w1539) + w501*w718) + w674 + w731;
    const scalar w1543 = w149*w1533 + w1517*w311 + w693 + w751;
    const scalar w1544 = w149*w1535 + w1517*w584 + w690 + w754;
    const scalar w1545 = w149*(w149*w1537 + w1533*w584 + w1535*w311 + w687 + w747) + w694 + w755;
    const scalar w1546 = -w1142*w684 + w1143*w675 - w1146*w314 + w1147*w696 + w1381*w684 + w1382*w709 - w1383*w314 + w1385*w675 + w1386*w696 - w1387*w709 + w312*w736 + w312*w758;
    const scalar w1547 = w1216*w771 + w1234*w771 + w1310*w710 + w1323*w710 + w1483*w482 + w661;
    const scalar w1548 = gel*w1495 + get*w1496 + w1513;
    const scalar w1549 = gel*w1499 + get*w1501 + w1521;
    return w175*w387 + w189*w649 + w21*w291 + w291*w298 + w291*w305 + w310*w387 + w335*w649 + w387*w572 + w392*w565 + w402*(w252*(w521*w545*w925 + w522*w853 + w532*w560*w915 - 2*w532*w924 + w533*(-w153*(w399*(-w887 - w893) + w474*w895 + w525*w894 + w905) - w313*(w450*(-w887 - w906) - w476*w903 + w528*w907 + w912) - w526*w775 - w586*(w443*(-w887 - w913) + w914)) + w548*w853 - 2*w560*w931 + w563*(w153*(-w153*w917 - w313*w918 - w392*w775 - w586*w919) + w313*(-w153*w918 - w313*w906 - 2*w549 - w586*w910) + w551*w775 + w586*(-w153*w919 - w313*w910 - 2*w550 - w586*w913)) + w923*(w149*(w399*(-w893 - w922) + w438*w895 + w492*w894 + w905) + w311*(w450*(-w906 - w922) - w453*w903 + w512*w907 + w912) + w509*w837 + w584*(w443*(-w913 - w922) + w914)) + w926*w928 + w929*(-w149*(w149*w917 + w311*w918 + w392*w837 + w584*w919) - w311*(w149*w918 + w311*w906 + 2*w534 + w584*w910) - w536*w837 - w584*(w149*w919 + w311*w910 + 2*w535 + w584*w913)) + w932*w934 + w935) + w482*w774 + w495 + w564*w866 + w564*w882 + w868 + w936) + w409*(w1221 + w1427 + w252*(w1040*w926 + w1043*w932 + w1426 + w533*(-w1412 - w313*(w1405 + w450*(-w1393 - w1402)) - w586*(w1400 + w443*(-w1393 - w1395))) + w563*(w1418 + w313*(-w1402*w313 - w1416) + w586*(-w1395*w586 - w1414)) + w923*(w1422 + w311*(w1421 + w450*(-w1402 - w1419)) + w584*(w1420 + w443*(-w1395 - w1419))) + w929*(-w1425 - w311*(w1402*w311 + w1424) - w584*(w1395*w584 + w1423)))) + w409*(w1247 + w1427 + w252*(w1039*w928 + w1042*w934 + w1426 + w533*(-w1412 - w313*(w1405 + w450*(-w1393 - w1429)) - w586*(w1400 + w443*(-w1393 - w1428))) + w563*(w1418 + w313*(-w1416 - w1429*w313) + w586*(-w1414 - w1428*w586)) + w923*(w1422 + w311*(w1421 + w450*(-w1419 - w1429)) + w584*(w1420 + w443*(-w1419 - w1428))) + w929*(-w1425 - w311*(w1424 + w1429*w311) - w584*(w1423 + w1428*w584)))) + w41*(w142*w774 - w158 + w251*w866 + w251*w882 + w252*(w141*(w153*w873 + w190*w775 + w313*w875 + w586*w874) + w216*w841 + w216*w844 + w217*(-w149*w873 - w190*w837 - w311*w875 - w584*w874) + w220*w833 - w220*w839 + w235*w853 + w835*w876 + w848*w877 + w849*w878 - w852*w879 + w858*w880 + w881) + w868 + w883) + w415*(w1314 + w1389 + w252*(w1145*w926 + w1149*w932 + w1388 + w533*(-w1367 - w313*(w1362 + w450*(-w1341 - w1358)) - w586*(w1356 + w443*(-w1341 - w1343))) + w563*(w1373 + w313*(-w1358*w313 - w1371) + w586*(-w1343*w586 - w1369)) + w923*(w1377 + w311*(w1376 + w450*(-w1358 - w1374)) + w584*(w1375 + w443*(-w1343 - w1374))) + w929*(-w1380 - w311*(w1358*w311 + w1379) - w584*(w1343*w584 + w1378)))) + w415*(w1337 + w1389 + w252*(w1144*w928 + w1148*w934 + w1388 + w533*(-w1367 - w313*(w1362 + w450*(-w1341 - w1391)) - w586*(w1356 + w443*(-w1341 - w1390))) + w563*(w1373 + w313*(-w1371 - w1391*w313) + w586*(-w1369 - w1390*w586)) + w923*(w1377 + w311*(w1376 + w450*(-w1374 - w1391)) + w584*(w1375 + w443*(-w1374 - w1390))) + w929*(-w1380 - w311*(w1379 + w1391*w311) - w584*(w1378 + w1390*w584)))) + w419*(w1006*w710 + w252*(-2*w1038*w684 + w1039*w1040 - 2*w1041*w709 + w1042*w1043 + w533*(-w153*(-w1015*w474 + w1021 + w399*(-w1009 - w1014) + w678*w907) - w313*(w1019*w476 - w1020*w680 + w1030 + w450*(-w1008 - w1022)) - w586*(w1032 + w443*(-w1008 - w1031) + w660*w682) - w681*w775) + w563*(w153*(-w1033*w153 - w1034*w313 - w1035*w586 - 2*w699) + w313*(-w1022*w313 - w1025*w586 - w1034*w153 - w650*w775) + w586*(-w1025*w313 - w1031*w586 - w1035*w153 - 2*w703) + w704*w775) + w675*w696*w925 + w676*w988 + w684*w709*w915 + w697*w988 + w923*(w149*(-w1015*w438 + w1021 + w399*(-w1014 - w1037) + w654*w907) + w311*(w1019*w453 - w1020*w665 + w1030 + w450*(-w1022 - w1036)) + w584*(w1032 + w443*(-w1031 - w1036) + w660*w672) + w670*w837) + w929*(-w149*(w1033*w149 + w1034*w311 + w1035*w584 + 2*w686) - w311*(w1022*w311 + w1025*w584 + w1034*w149 + w650*w837) - w584*(w1025*w311 + w1031*w584 + w1035*w149 + 2*w690) - w691*w837) + w935) + w482*w937 - w662 + w710*w994 + w936 + w995) + w421*(w1485 + w1547 + w252*(w1039*w1145 + w1042*w1149 + w1546 + w533*(-w1532 - w313*(w1525 + w450*(-w1512 - w1522)) - w586*(w1520 + w443*(-w1512 - w1514))) + w563*(w1538 + w313*(-w1522*w313 - w1536) + w586*(-w1514*w586 - w1534)) + w923*(w1542 + w311*(w1541 + w450*(-w1522 - w1539)) + w584*(w1540 + w443*(-w1514 - w1539))) + w929*(-w1545 - w311*(w1522*w311 + w1544) - w584*(w1514*w584 + w1543)))) + w421*(w1508 + w1547 + w252*(w1040*w1144 + w1043*w1148 + w1546 + w533*(-w1532 - w313*(w1525 + w450*(-w1512 - w1549)) - w586*(w1520 + w443*(-w1512 - w1548))) + w563*(w1538 + w313*(-w1536 - w1549*w313) + w586*(-w1534 - w1548*w586)) + w923*(w1542 + w311*(w1541 + w450*(-w1539 - w1549)) + w584*(w1540 + w443*(-w1539 - w1548))) + w929*(-w1545 - w311*(w1544 + w1549*w311) - w584*(w1543 + w1548*w584)))) + w424*(w1044*w482 + w1103*w771 + w1104 + w1115*w771 + w252*(-w1089*w1146 - 2*w1142*w744 + 2*w1143*w735 + w1144*w1145 + 2*w1147*w757 + w1148*w1149 + w533*(-w153*(w1128 + w399*(-w1009 - w1120)) - w313*(-w1125*w476 + w1134 + w450*(-w1008 - w1129) + w660*w740) - w586*(-w1126*w478 - w1127*w742 + w1138 + w443*(-w1008 - w1135)) - w743*w775) + w563*(w153*(-w1139*w153 - w1140*w313 - w1141*w586 - 2*w761) + w313*(-w1129*w313 - w1132*w586 - w1140*w153 - 2*w764) + w586*(-w1132*w313 - w1135*w586 - w1141*w153 - w712*w775) + w765*w775) + w735*w757*w925 + w744*w770*w915 + w923*(w149*(w1128 + w399*(-w1037 - w1120)) + w311*(-w1125*w453 + w1134 + w450*(-w1036 - w1129) + w660*w727) + w584*(-w1126*w457 - w1127*w732 + w1138 + w443*(-w1036 - w1135)) + w734*w837) + w929*(-w149*(w1139*w149 + w1140*w311 + w1141*w584 + 2*w748) - w311*(w1129*w311 + w1132*w584 + w1140*w149 + 2*w751) - w584*(w1132*w311 + w1135*w584 + w1141*w149 + w712*w837) - w752*w837) + w935) - w722 + w936) + w499*w772 + w506*w711 + w565*w575 + w565*w578 + w583*w649 + w60*(w1221 + w1236 + w252*(w1200*w876 + w1206*w878 + w1207*w880 + w1231*w370 + w1233 + w141*(w1223*w586 + w1226*w313 + w1229) + w217*(-w1223*w584 - w1226*w311 - w1230))) + w60*(w1236 + w1247 + w252*(w1233 + w1243*w876 + w1245*w878 + w1246*w880 + w1250*w984 + w141*(w1229 + w1248*w586 + w1249*w313) + w217*(-w1230 - w1248*w584 - w1249*w311))) + w650*w711 + w667*w772 + w71*(w1314 + w1324 + w252*(w1114*w856 + w1231*w632 + w1299*w876 + w1304*w878 + w1322 + w141*(w1316*w586 + w1318*w313 + w1320) + w217*(-w1316*w584 - w1318*w311 - w1321))) + w71*(w1324 + w1337 + w252*(w1093*w1250 + w1322 + w1335*w876 + w1336*w878 + w1340*w857 + w141*(w1320 + w1338*w586 + w1339*w313) + w217*(-w1321 - w1338*w584 - w1339*w311))) + w711*w728 + w712*w772 + w81*(w1006*w372 + w142*w937 + w252*(-w127*w987 + w141*(w1002*w153 + w1004*w313 + w1005*w586 + w336*w775) + w217*(-w1002*w149 - w1004*w311 - w1005*w584 - w336*w837) + w234*w989 + w342*w982 + w342*w983 + w345*w976 - w345*w980 + w876*w978 + w877*w985 + w878*w986 + w880*w991 + w881) - w319 + w372*w994 + w883 + w995) + w87*(w1485 + w1494 + w252*(w1114*w990 + w141*(w1487*w586 + w1489*w313 + w1491) + w1471*w876 + w1473*w877 + w1477*w878 + w1493 + w217*(-w1487*w584 - w1489*w311 - w1492))) + w87*(w1494 + w1508 + w252*(w1340*w358 + w141*(w1491 + w1509*w586 + w1510*w313) + w1493 + w1505*w876 + w1506*w877 + w1507*w878 + w217*(-w1492 - w1509*w584 - w1510*w311))) + w93*(w1044*w142 + w1103*w634 + w1104 + w1115*w634 + w252*(w1086*w606 + w1088*w876 + w1091*w603 + w1092*w603 + w1094*w877 + w1095*w878 - w1096*w127 + w1097*w118 + w1098*w1114 - 2*w1113*w606 + w141*(w1108*w153 + w1109*w313 + w1112*w586 + w602*w775) + w217*(-w1108*w149 - w1109*w311 - w1112*w584 - w602*w837) + w881) - w589 + w883);
}

inline scalar interfaceFlux(const point& pt, const vector& n, const scalar t, const Params& p)
{
    const scalar r1 = p.r1;
    const scalar r2 = p.r2;
    const scalar gamma1 = p.gamma1;
    const scalar gamma2 = p.gamma2;
    const scalar k = p.k;
    const scalar gil = p.gil;
    const scalar git = p.git;
    const scalar gin = p.gin;
    const scalar gel = p.gel;
    const scalar get = p.get;
    const scalar gen = p.gen;
    const scalar theta0 = p.theta0;
    const scalar theta1 = p.theta1;
    const scalar sigma0In = p.sigma0In;
    const scalar sigma0Out = p.sigma0Out;
    const scalar sigmaEps = p.sigmaEps;
    const scalar fA0 = p.fA[0];
    const scalar fA1 = p.fA[1];
    const scalar fA2 = p.fA[2];
    const scalar fA3 = p.fA[3];
    const scalar cmAmplitude = p.cmAmplitude;
    const scalar cmOmega = p.cmOmega;
    const scalar x = pt.x();
    const scalar y = pt.y();
    const scalar z = pt.z();
    const scalar St = sqrt(1 + t);
    const scalar Tt = cmAmplitude*sin(cmOmega*t);

    const scalar w0 = cos(gamma2);
    const scalar w1 = sin(gamma2);
    const scalar w2 = w0*w1;
    const scalar w3 = sin(gamma1);
    const scalar w4 = gen*w3;
    const scalar w5 = cos(gamma1);
    const scalar w6 = (r1*r1);
    const scalar w7 = -(r2*r2);
    const scalar w8 = (1.0/((-w6 - w7)));
    const scalar w9 = (x*x) + (y*y) + (z*z);
    const scalar w10 = w8*(-w6 + w9);
    const scalar w11 = theta0 + theta1*w10;
    const scalar w12 = cos(w11);
    const scalar w13 = sin(w11);
    const scalar w14 = w13*w3;
    const scalar w15 = w0*w14;
    const scalar w16 = w12*w5 - w15;
    const scalar w17 = gel*w16;
    const scalar w18 = w1*w13;
    const scalar w19 = w13*w5;
    const scalar w20 = w12*w3;
    const scalar w21 = w0*w20;
    const scalar w22 = -w19 - w21;
    const scalar w23 = get*w22;
    const scalar w24 = w1*w12;
    const scalar w25 = w17*w18 + w2*w4 + w23*w24;
    const scalar w26 = k*y;
    const scalar w27 = w26 + 3.0/10.0;
    const scalar w28 = cos(w27);
    const scalar w29 = fA3*exp((2.0/5.0)*z);
    const scalar w30 = (2.0/5.0)*w29;
    const scalar w31 = St*(fA1*y + w28*w30);
    const scalar w32 = r1 + r2;
    const scalar w33 = (1.0/(w32));
    const scalar w34 = w10*w33;
    const scalar w35 = (1.0/(r2));
    const scalar w36 = (w0*w0);
    const scalar w37 = (w1*w1);
    const scalar w38 = (w13*w13)*w37;
    const scalar w39 = (w12*w12)*w37;
    const scalar w40 = gel*w38 + gen*w36 + get*w39;
    const scalar w41 = w35*w40;
    const scalar w42 = w2*w5;
    const scalar w43 = w0*w19;
    const scalar w44 = w20 + w43;
    const scalar w45 = gel*w44;
    const scalar w46 = w0*w12*w5 - w14;
    const scalar w47 = get*w46;
    const scalar w48 = -gen*w42 + w18*w45 + w24*w47;
    const scalar w49 = w35*w48;
    const scalar w50 = w25*w35;
    const scalar w51 = w35*(w41*z + w49*y + w50*x);
    const scalar w52 = w37*(w5*w5);
    const scalar w53 = (w44*w44);
    const scalar w54 = (w46*w46);
    const scalar w55 = gel*w53 + gen*w52 + get*w54;
    const scalar w56 = w35*w55;
    const scalar w57 = w37*w5;
    const scalar w58 = w17*w44 + w23*w46 - w4*w57;
    const scalar w59 = w35*w58;
    const scalar w60 = w35*(w49*z + w56*y + w59*x);
    const scalar w61 = (w3*w3)*w37;
    const scalar w62 = (w16*w16);
    const scalar w63 = (w22*w22);
    const scalar w64 = gel*w62 + gen*w61 + get*w63;
    const scalar w65 = w35*w64;
    const scalar w66 = w35*(w50*z + w59*y + w65*x);
    const scalar w67 = w51*z + w60*y + w66*x;
    const scalar w68 = (1.0/(w67));
    const scalar w69 = k*x;
    const scalar w70 = w69 + 1.0/5.0;
    const scalar w71 = sin(w70);
    const scalar w72 = k*z;
    const scalar w73 = (1.0/2.0)*sin(w72) + 1;
    const scalar w74 = cos(w26);
    const scalar w75 = sigmaEps*w73*w74;
    const scalar w76 = w71*w75 + 1;
    const scalar w77 = -sigma0Out*w76;
    const scalar w78 = -w40 - w77;
    const scalar w79 = -w25;
    const scalar w80 = cos(w69);
    const scalar w81 = k*w80;
    const scalar w82 = fA2*((1.0/2.0)*y + 1);
    const scalar w83 = St*(fA0 + w81*w82);
    const scalar w84 = -w48;
    const scalar w85 = sin(w69);
    const scalar w86 = (1.0/2.0)*fA2;
    const scalar w87 = k*sin(w27);
    const scalar w88 = St*(fA1*z - w29*w87 + w85*w86);
    const scalar w89 = w79*w83 + w84*w88;
    const scalar w90 = w35*(w31*w78 + w89);
    const scalar w91 = -w55 - w77;
    const scalar w92 = -w58;
    const scalar w93 = w31*w84 + w83*w92;
    const scalar w94 = w35*(w88*w91 + w93);
    const scalar w95 = -w64 - w77;
    const scalar w96 = w31*w79 + w88*w92;
    const scalar w97 = w35*(w83*w95 + w96);
    const scalar w98 = w90*z + w94*y + w97*x;
    const scalar w99 = w68*w98;
    const scalar w100 = (1.0/2.0)*w10*w32*w35;
    const scalar w101 = (1.0/(r1));
    const scalar w102 = -sigma0In*w76;
    const scalar w103 = -w102 - w40;
    const scalar w104 = w101*(w103*w31 + w89);
    const scalar w105 = -w102 - w55;
    const scalar w106 = w101*(w105*w88 + w93);
    const scalar w107 = -w102 - w64;
    const scalar w108 = w101*(w107*w83 + w96);
    const scalar w109 = -w104*z - w106*y - w108*x;
    const scalar w110 = w101*w40;
    const scalar w111 = w101*w48;
    const scalar w112 = w101*w25;
    const scalar w113 = w101*(-w110*z - w111*y - w112*x);
    const scalar w114 = w101*w55;
    const scalar w115 = w101*w58;
    const scalar w116 = w101*(-w111*z - w114*y - w115*x);
    const scalar w117 = w101*w64;
    const scalar w118 = w101*(-w112*z - w115*y - w117*x);
    const scalar w119 = -w113*z - w116*y - w118*x;
    const scalar w120 = (1.0/(w119));
    const scalar w121 = (1.0/2.0)*w101*(1 - w10);
    const scalar w122 = w121*w32;
    const scalar w123 = w120*w122;
    const scalar w124 = w100*w99 + w109*w123;
    const scalar w125 = 2*z;
    const scalar w126 = w124*w125;
    const scalar w127 = w7 + w9;
    const scalar w128 = w127*w33*w8;
    const scalar w129 = w35*z;
    const scalar w130 = w32*w8;
    const scalar w131 = w129*w130;
    const scalar w132 = w101*z;
    const scalar w133 = w109*w120*w130;
    const scalar w134 = k*sigma0Out;
    const scalar w135 = sigmaEps*w71;
    const scalar w136 = w134*w135;
    const scalar w137 = cos(w72);
    const scalar w138 = (1.0/2.0)*w137*w74;
    const scalar w139 = theta1*w8;
    const scalar w140 = 4*z;
    const scalar w141 = w139*w140;
    const scalar w142 = w12*w13*w37;
    const scalar w143 = w141*w142;
    const scalar w144 = gel*w143;
    const scalar w145 = get*w143;
    const scalar w146 = -w144 + w145;
    const scalar w147 = w28*w29;
    const scalar w148 = (4.0/25.0)*St*w147;
    const scalar w149 = w125*w139;
    const scalar w150 = w149*w24;
    const scalar w151 = w139*w19;
    const scalar w152 = -w125*w151 - w149*w21;
    const scalar w153 = gel*w18;
    const scalar w154 = w12*w139*w5;
    const scalar w155 = -w125*w154 + w149*w15;
    const scalar w156 = get*w24;
    const scalar w157 = -2*get*theta1*w1*w13*w22*w8*z + w150*w17 + w152*w153 + w155*w156;
    const scalar w158 = -w157;
    const scalar w159 = St*(fA1 - w30*w87);
    const scalar w160 = w159*w84;
    const scalar w161 = 2*theta1*w0*w12*w5*w8*z - w14*w149;
    const scalar w162 = -w149*w20 - w149*w43;
    const scalar w163 = -2*get*theta1*w1*w13*w46*w8*z + w150*w45 + w153*w161 + w156*w162;
    const scalar w164 = -w163;
    const scalar w165 = w158*w83 + w160 + w164*w88;
    const scalar w166 = -1.0/2.0*k*sigma0Out*sigmaEps*w137*w71*w74;
    const scalar w167 = 4*theta1*w0*w12*w5*w8*z - w14*w141;
    const scalar w168 = -w141*w20 - w141*w43;
    const scalar w169 = w167*w45 + w168*w47;
    const scalar w170 = w152*w45 + w155*w47 + w161*w17 + w162*w23;
    const scalar w171 = -w170;
    const scalar w172 = w148*w84 + w164*w31 + w171*w83;
    const scalar w173 = w35*y;
    const scalar w174 = -w140*w151 - w141*w21;
    const scalar w175 = -w140*w154 + w141*w15;
    const scalar w176 = w17*w174 + w175*w23;
    const scalar w177 = w148*w79 + w158*w31 + w159*w92 + w171*w88;
    const scalar w178 = w35*x;
    const scalar w179 = w100*w68;
    const scalar w180 = k*sigma0In;
    const scalar w181 = w135*w180;
    const scalar w182 = -1.0/2.0*k*sigma0In*sigmaEps*w137*w71*w74;
    const scalar w183 = w101*y;
    const scalar w184 = w101*x;
    const scalar w185 = w144 - w145;
    const scalar w186 = w100*w98/(w67*w67);
    const scalar w187 = w109*w122/(w119*w119);
    const scalar w188 = w127*w34;
    const scalar w189 = (1.0/5.0)*z;
    const scalar w190 = gin*w36;
    const scalar w191 = gil*w38;
    const scalar w192 = git*w39;
    const scalar w193 = w190 + w191 + w192;
    const scalar w194 = w193*w35;
    const scalar w195 = gin*w42;
    const scalar w196 = gil*w44;
    const scalar w197 = w18*w196;
    const scalar w198 = git*w46;
    const scalar w199 = w198*w24;
    const scalar w200 = -w195 + w197 + w199;
    const scalar w201 = w200*w35;
    const scalar w202 = gin*w3;
    const scalar w203 = w2*w202;
    const scalar w204 = gil*w16;
    const scalar w205 = w18*w204;
    const scalar w206 = git*w22;
    const scalar w207 = w206*w24;
    const scalar w208 = w203 + w205 + w207;
    const scalar w209 = w208*w35;
    const scalar w210 = w35*(w194*z + w201*y + w209*x);
    const scalar w211 = gin*w52;
    const scalar w212 = gil*w53;
    const scalar w213 = git*w54;
    const scalar w214 = w211 + w212 + w213;
    const scalar w215 = w214*w35;
    const scalar w216 = w202*w57;
    const scalar w217 = w204*w44;
    const scalar w218 = w206*w46;
    const scalar w219 = -w216 + w217 + w218;
    const scalar w220 = w219*w35;
    const scalar w221 = w35*(w201*z + w215*y + w220*x);
    const scalar w222 = gil*w62 + gin*w61 + git*w63;
    const scalar w223 = w222*w35;
    const scalar w224 = w35*(w209*z + w220*y + w223*x);
    const scalar w225 = w210*z + w221*y + w224*x;
    const scalar w226 = (1.0/(w225));
    const scalar w227 = (1.0/5.0)*k*w85;
    const scalar w228 = w189*w193 - 1.0/10.0*w195 + (1.0/10.0)*w197 + (1.0/10.0)*w199 - w208*w227;
    const scalar w229 = w228*w35;
    const scalar w230 = w189*w200 + (1.0/10.0)*w211 + (1.0/10.0)*w212 + (1.0/10.0)*w213 - w219*w227;
    const scalar w231 = w230*w35;
    const scalar w232 = w189*w208 - 1.0/10.0*w216 + (1.0/10.0)*w217 + (1.0/10.0)*w218 - w222*w227;
    const scalar w233 = w232*w35;
    const scalar w234 = -w229*z - w231*y - w233*x;
    const scalar w235 = w226*w234;
    const scalar w236 = w101*w228;
    const scalar w237 = w101*w230;
    const scalar w238 = w101*w232;
    const scalar w239 = w236*z + w237*y + w238*x;
    const scalar w240 = w101*w193;
    const scalar w241 = w101*w200;
    const scalar w242 = w101*w208;
    const scalar w243 = w101*(-w240*z - w241*y - w242*x);
    const scalar w244 = w101*w214;
    const scalar w245 = w101*w219;
    const scalar w246 = w101*(-w241*z - w244*y - w245*x);
    const scalar w247 = w101*w222;
    const scalar w248 = w101*(-w242*z - w245*y - w247*x);
    const scalar w249 = -w243*z - w246*y - w248*x;
    const scalar w250 = w32/w249;
    const scalar w251 = w121*w250;
    const scalar w252 = w100*w235 + w239*w251;
    const scalar w253 = w125*w252;
    const scalar w254 = w239*w250*w8;
    const scalar w255 = w196*w24;
    const scalar w256 = w139*w189;
    const scalar w257 = w18*w198;
    const scalar w258 = gil*w18;
    const scalar w259 = w161*w258;
    const scalar w260 = git*w24;
    const scalar w261 = w162*w260;
    const scalar w262 = gil*w143 - git*w143;
    const scalar w263 = -w149*w18*w206 + w150*w204 + w152*w258 + w155*w260;
    const scalar w264 = w189*w262 + (1.0/5.0)*w190 + (1.0/5.0)*w191 + (1.0/5.0)*w192 - w227*w263 + w255*w256 - w256*w257 + (1.0/10.0)*w259 + (1.0/10.0)*w261;
    const scalar w265 = w167*w196;
    const scalar w266 = w168*w198;
    const scalar w267 = w149*w255 - w149*w257 + w259 + w261;
    const scalar w268 = w152*w196;
    const scalar w269 = w161*w204;
    const scalar w270 = w155*w198;
    const scalar w271 = w162*w206;
    const scalar w272 = w268 + w269 + w270 + w271;
    const scalar w273 = w189*w267 - 1.0/5.0*w195 + (1.0/5.0)*w197 + (1.0/5.0)*w199 - w227*w272 + (1.0/10.0)*w265 + (1.0/10.0)*w266;
    const scalar w274 = w174*w204 + w175*w206;
    const scalar w275 = w189*w263 + (1.0/5.0)*w203 + (1.0/5.0)*w205 + (1.0/5.0)*w207 - w227*w274 + (1.0/10.0)*w268 + (1.0/10.0)*w269 + (1.0/10.0)*w270 + (1.0/10.0)*w271;
    const scalar w276 = w100*w226;
    const scalar w277 = w265 + w266;
    const scalar w278 = w100*w234/(w225*w225);
    const scalar w279 = w122*w239/(w249*w249);
    const scalar w280 = Tt*(w128*w253 + w188*(w131*w235 - w132*w254 + w251*(w132*w264 + w183*w273 + w184*w275 + w236) + w276*(-w129*w264 - w173*w273 - w178*w275 - w229) + w278*(-w129*(w129*w262 + w173*w267 + w178*w263 + w194) - w173*(w129*w267 + w173*w277 + w178*w272 + w201) - w178*(w129*w263 + w173*w272 + w178*w274 + w209) - w210) + w279*(w132*(-w132*w262 - w183*w267 - w184*w263 - w240) + w183*(-w132*w267 - w183*w277 - w184*w272 - w241) + w184*(-w132*w263 - w183*w272 - w184*w274 - w242) + w243)) + w189 + w253*w34) + w126*w128 + w126*w34 + w188*(w123*(-w104 - w132*(w103*w148 + w165 + w31*(w138*w181 + w146)) - w183*(w105*w159 + w172 + w88*(-w169 - w182)) - w184*(w177 + w83*(-w176 - w182))) + w131*w99 - w132*w133 + w179*(w129*(w148*w78 + w165 + w31*(w136*w138 + w146)) + w173*(w159*w91 + w172 + w88*(-w166 - w169)) + w178*(w177 + w83*(-w166 - w176)) + w90) + w186*(-w129*(w129*w185 + w157*w178 + w163*w173 + w41) - w173*(w129*w163 + w169*w173 + w170*w178 + w49) - w178*(w129*w157 + w170*w173 + w176*w178 + w50) - w51) + w187*(w113 + w132*(-w110 - w132*w185 - w157*w184 - w163*w183) + w183*(-w111 - w132*w163 - w169*w183 - w170*w184) + w184*(-w112 - w132*w157 - w170*w183 - w176*w184))) + w31;
    const scalar w281 = 2*y;
    const scalar w282 = w124*w281;
    const scalar w283 = w130*w173;
    const scalar w284 = w73*sin(w26);
    const scalar w285 = w136*w284;
    const scalar w286 = 4*y;
    const scalar w287 = w139*w286;
    const scalar w288 = w142*w287;
    const scalar w289 = gel*w288 - get*w288;
    const scalar w290 = w139*w281;
    const scalar w291 = w24*w290;
    const scalar w292 = -w151*w281 - w21*w290;
    const scalar w293 = w15*w290 - w154*w281;
    const scalar w294 = -2*get*theta1*w1*w13*w22*w8*y + w153*w292 + w156*w293 + w17*w291;
    const scalar w295 = -w294;
    const scalar w296 = 2*theta1*w0*w12*w5*w8*y - w14*w290;
    const scalar w297 = -w20*w290 - w290*w43;
    const scalar w298 = -2*get*theta1*w1*w13*w46*w8*y + w153*w296 + w156*w297 + w291*w45;
    const scalar w299 = -w298;
    const scalar w300 = w81*w86;
    const scalar w301 = St*w300;
    const scalar w302 = (k*k);
    const scalar w303 = w147*w302;
    const scalar w304 = St*w303;
    const scalar w305 = w295*w83 + w299*w88 + w301*w79 - w304*w84;
    const scalar w306 = 4*theta1*w0*w12*w5*w8*y - w14*w287;
    const scalar w307 = -w20*w287 - w287*w43;
    const scalar w308 = w306*w45 + w307*w47;
    const scalar w309 = St*w92;
    const scalar w310 = w300*w309;
    const scalar w311 = w17*w296 + w23*w297 + w292*w45 + w293*w47;
    const scalar w312 = -w311;
    const scalar w313 = w160 + w299*w31 + w310 + w312*w83;
    const scalar w314 = -w151*w286 - w21*w287;
    const scalar w315 = w15*w287 - w154*w286;
    const scalar w316 = w17*w314 + w23*w315;
    const scalar w317 = w159*w79 + w295*w31 - w303*w309 + w312*w88;
    const scalar w318 = w181*w284;
    const scalar w319 = St*w107;
    const scalar w320 = w252*w281;
    const scalar w321 = (1.0/5.0)*w139;
    const scalar w322 = w321*y;
    const scalar w323 = w258*w296;
    const scalar w324 = w260*w297;
    const scalar w325 = gil*w288 - git*w288;
    const scalar w326 = -w18*w206*w290 + w204*w291 + w258*w292 + w260*w293;
    const scalar w327 = w189*w325 - w227*w326 + w255*w322 - w257*w322 + (1.0/10.0)*w323 + (1.0/10.0)*w324;
    const scalar w328 = w196*w306;
    const scalar w329 = w198*w307;
    const scalar w330 = w255*w290 - w257*w290 + w323 + w324;
    const scalar w331 = w196*w292;
    const scalar w332 = w204*w296;
    const scalar w333 = w198*w293;
    const scalar w334 = w206*w297;
    const scalar w335 = w331 + w332 + w333 + w334;
    const scalar w336 = w189*w330 - w227*w335 + (1.0/10.0)*w328 + (1.0/10.0)*w329;
    const scalar w337 = w204*w314 + w206*w315;
    const scalar w338 = w189*w326 - w227*w337 + (1.0/10.0)*w331 + (1.0/10.0)*w332 + (1.0/10.0)*w333 + (1.0/10.0)*w334;
    const scalar w339 = w328 + w329;
    const scalar w340 = Tt*(w128*w320 + w188*(-w183*w254 + w235*w283 + w251*(w132*w327 + w183*w336 + w184*w338 + w237) + w276*(-w129*w327 - w173*w336 - w178*w338 - w231) + w278*(-w129*(w129*w325 + w173*w330 + w178*w326 + w201) - w173*(w129*w330 + w173*w339 + w178*w335 + w215) - w178*(w129*w326 + w173*w335 + w178*w337 + w220) - w221) + w279*(w132*(-w132*w325 - w183*w330 - w184*w326 - w241) + w183*(-w132*w330 - w183*w339 - w184*w335 - w244) + w184*(-w132*w326 - w183*w335 - w184*w337 - w245) + w246)) + w320*w34 + 1.0/10.0) + w128*w282 + w188*(w123*(-w106 - w132*(w103*w159 + w305 + w31*(-w289 - w318)) - w183*(-w105*w304 + w313 + w88*(-w308 - w318)) - w184*(w300*w319 + w317 + w83*(-w316 - w318))) - w133*w183 + w179*(w129*(w159*w78 + w305 + w31*(-w285 - w289)) + w173*(-w304*w91 + w313 + w88*(-w285 - w308)) + w178*(w301*w95 + w317 + w83*(-w285 - w316)) + w94) + w186*(-w129*(w129*w289 + w173*w298 + w178*w294 + w49) - w173*(w129*w298 + w173*w308 + w178*w311 + w56) - w178*(w129*w294 + w173*w311 + w178*w316 + w59) - w60) + w187*(w116 + w132*(-w111 - w132*w289 - w183*w298 - w184*w294) + w183*(-w114 - w132*w298 - w183*w308 - w184*w311) + w184*(-w115 - w132*w294 - w183*w311 - w184*w316)) + w283*w99) + w282*w34 + w88;
    const scalar w341 = 2*x;
    const scalar w342 = w124*w341;
    const scalar w343 = w130*w178;
    const scalar w344 = cos(w70);
    const scalar w345 = w344*w75;
    const scalar w346 = 4*x;
    const scalar w347 = w139*w346;
    const scalar w348 = w142*w347;
    const scalar w349 = gel*w348;
    const scalar w350 = get*w348;
    const scalar w351 = -w349 + w350;
    const scalar w352 = w139*w341;
    const scalar w353 = w24*w352;
    const scalar w354 = -w151*w341 - w21*w352;
    const scalar w355 = w15*w352 - w154*w341;
    const scalar w356 = -2*get*theta1*w1*w13*w22*w8*x + w153*w354 + w156*w355 + w17*w353;
    const scalar w357 = -w356;
    const scalar w358 = 2*theta1*w0*w12*w5*w8*x - w14*w352;
    const scalar w359 = -w20*w352 - w352*w43;
    const scalar w360 = -2*get*theta1*w1*w13*w46*w8*x + w153*w358 + w156*w359 + w353*w45;
    const scalar w361 = -w360;
    const scalar w362 = w302*w82*w85;
    const scalar w363 = St*w362;
    const scalar w364 = w301*w84 + w357*w83 + w361*w88 - w363*w79;
    const scalar w365 = -k*sigma0Out*sigmaEps*w344*w73*w74;
    const scalar w366 = 4*theta1*w0*w12*w5*w8*x - w14*w347;
    const scalar w367 = -w20*w347 - w347*w43;
    const scalar w368 = w366*w45 + w367*w47;
    const scalar w369 = w17*w358 + w23*w359 + w354*w45 + w355*w47;
    const scalar w370 = -w369;
    const scalar w371 = -w309*w362 + w31*w361 + w370*w83;
    const scalar w372 = -w151*w346 - w21*w347;
    const scalar w373 = w15*w347 - w154*w346;
    const scalar w374 = w17*w372 + w23*w373;
    const scalar w375 = w31*w357 + w310 + w370*w88;
    const scalar w376 = -k*sigma0In*sigmaEps*w344*w73*w74;
    const scalar w377 = w349 - w350;
    const scalar w378 = w252*w341;
    const scalar w379 = w321*x;
    const scalar w380 = w258*w358;
    const scalar w381 = w260*w359;
    const scalar w382 = gil*w348 - git*w348;
    const scalar w383 = (1.0/5.0)*w302*w80;
    const scalar w384 = w18*w352;
    const scalar w385 = w204*w353 - w206*w384 + w258*w354 + w260*w355;
    const scalar w386 = w189*w382 - w208*w383 - w227*w385 + w255*w379 - w257*w379 + (1.0/10.0)*w380 + (1.0/10.0)*w381;
    const scalar w387 = w196*w354;
    const scalar w388 = w204*w358;
    const scalar w389 = w198*w355;
    const scalar w390 = w206*w359;
    const scalar w391 = w204*w372 + w206*w373;
    const scalar w392 = w189*w385 - w222*w383 - w227*w391 + (1.0/10.0)*w387 + (1.0/10.0)*w388 + (1.0/10.0)*w389 + (1.0/10.0)*w390;
    const scalar w393 = w196*w366;
    const scalar w394 = w198*w367;
    const scalar w395 = w196*w353 - w198*w384 + w380 + w381;
    const scalar w396 = w387 + w388 + w389 + w390;
    const scalar w397 = w189*w395 - w219*w383 - w227*w396 + (1.0/10.0)*w393 + (1.0/10.0)*w394;
    const scalar w398 = w393 + w394;
    const scalar w399 = Tt*(w128*w378 + w188*(-w184*w254 + w235*w343 + w251*(w132*w386 + w183*w397 + w184*w392 + w238) + w276*(-w129*w386 - w173*w397 - w178*w392 - w233) + w278*(-w129*(w129*w382 + w173*w395 + w178*w385 + w209) - w173*(w129*w395 + w173*w398 + w178*w396 + w220) - w178*(w129*w385 + w173*w396 + w178*w391 + w223) - w224) + w279*(w132*(-w132*w382 - w183*w395 - w184*w385 - w242) + w183*(-w132*w395 - w183*w398 - w184*w396 - w245) + w184*(-w132*w385 - w183*w396 - w184*w391 - w247) + w248)) - w227 + w34*w378) + w128*w342 + w188*(w123*(-w108 - w132*(w31*(w180*w345 + w351) + w364) - w183*(w105*w301 + w371 + w88*(-w368 - w376)) - w184*(-w319*w362 + w375 + w83*(-w374 - w376))) - w133*w184 + w179*(w129*(w31*(w134*w345 + w351) + w364) + w173*(w301*w91 + w371 + w88*(-w365 - w368)) + w178*(-w363*w95 + w375 + w83*(-w365 - w374)) + w97) + w186*(-w129*(w129*w377 + w173*w360 + w178*w356 + w50) - w173*(w129*w360 + w173*w368 + w178*w369 + w59) - w178*(w129*w356 + w173*w369 + w178*w374 + w65) - w66) + w187*(w118 + w132*(-w112 - w132*w377 - w183*w360 - w184*w356) + w183*(-w115 - w132*w360 - w183*w368 - w184*w369) + w184*(-w117 - w132*w356 - w183*w369 - w184*w374)) + w343*w99) + w34*w342 + w83;
    return n.x()*(w25*w280 + w340*w58 + w399*w64) + n.y()*(w280*w48 + w340*w55 + w399*w58) + n.z()*(w25*w399 + w280*w40 + w340*w48);
}

inline scalar G(const point& pt, const Params& p)
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
    {0, 0, 1, 1, 0, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, -0.40514908577374803},
    {0, 0, 1, 1, 0, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0.13700000000000001, -0.4320112889689996},
    {0, 0, 1, 1, 0, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, -0.33311674555542337},
    {0, 0, 1, 1, 0, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0.13700000000000001, -0.35520305901645821},
    {0, 0, 1, 1, 0, -0.12, 0.64000000000000001, 0.70999999999999996, 0, 0.027329045308920583},
    {0, 0, 1, 1, 0, -0.12, 0.64000000000000001, 0.70999999999999996, 0.13700000000000001, 0.029141016245047596},
    {0, 1, 1, 1, 0, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, -0.28433152194336586},
    {0, 1, 1, 1, 0, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0.13700000000000001, -0.30318327648372556},
    {0, 1, 1, 1, 0, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, -0.33124571675229753},
    {0, 1, 1, 1, 0, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0.13700000000000001, -0.353207977222326},
    {0, 1, 1, 1, 0, -0.12, 0.64000000000000001, 0.70999999999999996, 0, -0.030106933704326478},
    {0, 1, 1, 1, 0, -0.12, 0.64000000000000001, 0.70999999999999996, 0.13700000000000001, -0.032103084255196083},
    {1, 0, 1, 0, 0, 0, 0, 0, 0, 0.14330047336884089},
    {1, 0, 1, 0, 0, 0, 0, 0, 0.13700000000000001, 0.15280158436421207},
    {1, 0, 1, 0, 0, 0.070000000000000007, 0.13, -0.050000000000000003, 0, 0.22428334221294519},
    {1, 0, 1, 0, 0, 0.070000000000000007, 0.13, -0.050000000000000003, 0.13700000000000001, 0.23915378107948812},
    {1, 0, 1, 0, 0, 0.20999999999999999, -0.29999999999999999, 0.25, 0, 0.38911062314617195},
    {1, 0, 1, 0, 0, 0.20999999999999999, -0.29999999999999999, 0.25, 0.13700000000000001, 0.41490944385539713},
    {1, 0, 1, 0, 0, -0.28000000000000003, 0.23999999999999999, -0.33000000000000002, 0, -0.3800782337906668},
    {1, 0, 1, 0, 0, -0.28000000000000003, 0.23999999999999999, -0.33000000000000002, 0.13700000000000001, -0.40527818883111516},
    {1, 0, 1, 1, 0, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, 0.56129059081820798},
    {1, 0, 1, 1, 0, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0.13700000000000001, 0.5985052913607174},
    {1, 0, 1, 1, 0, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 0.45706270512770009},
    {1, 0, 1, 1, 0, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0.13700000000000001, 0.48736688620381879},
    {1, 0, 1, 1, 0, -0.12, 0.64000000000000001, 0.70999999999999996, 0, -0.26597866992688857},
    {1, 0, 1, 1, 0, -0.12, 0.64000000000000001, 0.70999999999999996, 0.13700000000000001, -0.28361359328734448},
    {1, 0, 1, 2, 0, 0.62, 0.57999999999999996, -0.58999999999999997, 0, 0.53526941704720399},
    {1, 0, 1, 2, 0, 0.62, 0.57999999999999996, -0.58999999999999997, 0.13700000000000001, 0.5707588611797656},
    {1, 0, 1, 2, 0, -1.1000000000000001, 0.40000000000000002, 0.82999999999999996, 0, -0.11231176414675761},
    {1, 0, 1, 2, 0, -1.1000000000000001, 0.40000000000000002, 0.82999999999999996, 0.13700000000000001, -0.11975826109235892},
    {1, 0, 1, 2, 0, 1.47, -1.4099999999999999, 1.3799999999999999, 0, -0.2089984284522364},
    {1, 0, 1, 2, 0, 1.47, -1.4099999999999999, 1.3799999999999999, 0.13700000000000001, -0.22285544664555251},
    {1, 1, 1, 0, 0, 0, 0, 0, 0, 0.19106729782512122},
    {1, 1, 1, 0, 0, 0, 0, 0, 0.13700000000000001, 0.22718853255285659},
    {1, 1, 1, 0, 0, 0.070000000000000007, 0.13, -0.050000000000000003, 0, 0.24606032173458534},
    {1, 1, 1, 0, 0, 0.070000000000000007, 0.13, -0.050000000000000003, 0.13700000000000001, 0.28717067343462233},
    {1, 1, 1, 0, 0, 0.20999999999999999, -0.29999999999999999, 0.25, 0, 0.39305281500064043},
    {1, 1, 1, 0, 0, 0.20999999999999999, -0.29999999999999999, 0.25, 0.13700000000000001, 0.43344964078309273},
    {1, 1, 1, 0, 0, -0.28000000000000003, 0.23999999999999999, -0.33000000000000002, 0, -0.30523892587092694},
    {1, 1, 1, 0, 0, -0.28000000000000003, 0.23999999999999999, -0.33000000000000002, 0.13700000000000001, -0.30658337169044825},
    {1, 1, 1, 1, 0, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, 0.50712065405564621},
    {1, 1, 1, 1, 0, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0.13700000000000001, 0.55487287779118766},
    {1, 1, 1, 1, 0, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 0.39057953103526832},
    {1, 1, 1, 1, 0, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0.13700000000000001, 0.40790082896306318},
    {1, 1, 1, 1, 0, -0.12, 0.64000000000000001, 0.70999999999999996, 0, -0.26056033526850664},
    {1, 1, 1, 1, 0, -0.12, 0.64000000000000001, 0.70999999999999996, 0.13700000000000001, -0.24187801523843336},
    {1, 1, 1, 2, 0, 0.62, 0.57999999999999996, -0.58999999999999997, 0, 0.40650375048370047},
    {1, 1, 1, 2, 0, 0.62, 0.57999999999999996, -0.58999999999999997, 0.13700000000000001, 0.43677223904846413},
    {1, 1, 1, 2, 0, -1.1000000000000001, 0.40000000000000002, 0.82999999999999996, 0, -0.058266090958165311},
    {1, 1, 1, 2, 0, -1.1000000000000001, 0.40000000000000002, 0.82999999999999996, 0.13700000000000001, -0.096600158881714221},
    {1, 1, 1, 2, 0, 1.47, -1.4099999999999999, 1.3799999999999999, 0, -0.4128912616840319},
    {1, 1, 1, 2, 0, 1.47, -1.4099999999999999, 1.3799999999999999, 0.13700000000000001, -0.83233060761431177},
    {2, 0, 1, 1, 0, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, 0.15614150504445989},
    {2, 0, 1, 1, 0, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0.13700000000000001, 0.16649400239171777},
    {2, 0, 1, 1, 0, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 0.12394595957227672},
    {2, 0, 1, 1, 0, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0.13700000000000001, 0.13216382718736061},
    {2, 0, 1, 1, 0, -0.12, 0.64000000000000001, 0.70999999999999996, 0, -0.23864962461796799},
    {2, 0, 1, 1, 0, -0.12, 0.64000000000000001, 0.70999999999999996, 0.13700000000000001, -0.25447257704229687},
    {2, 1, 1, 1, 0, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, 0.22278913211228032},
    {2, 1, 1, 1, 0, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0.13700000000000001, 0.2516896013074621},
    {2, 1, 1, 1, 0, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 0.059333814282970763},
    {2, 1, 1, 1, 0, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0.13700000000000001, 0.054692851740737181},
    {2, 1, 1, 1, 0, -0.12, 0.64000000000000001, 0.70999999999999996, 0, -0.29066726897283307},
    {2, 1, 1, 1, 0, -0.12, 0.64000000000000001, 0.70999999999999996, 0.13700000000000001, -0.27398109949362942},
    {3, 0, 1, 0, 0, 0, 0, 0, 0, 1.5566370614359173},
    {3, 0, 1, 0, 1, 0, 0, 0, 0, -0.13926061653519028},
    {3, 0, 1, 0, 2, 0, 0, 0, 0, 0.057320189347536361},
    {3, 0, 1, 0, 0, 0, 0, 0, 0.13700000000000001, 1.6598452445808793},
    {3, 0, 1, 0, 1, 0, 0, 0, 0.13700000000000001, -0.14849387685791832},
    {3, 0, 1, 0, 2, 0, 0, 0, 0.13700000000000001, 0.061120633745684834},
    {3, 0, 1, 0, 0, 0.070000000000000007, 0.13, -0.050000000000000003, 0, 1.6060874281041329},
    {3, 0, 1, 0, 1, 0.070000000000000007, 0.13, -0.050000000000000003, 0, -0.26690003436736309},
    {3, 0, 1, 0, 2, 0.070000000000000007, 0.13, -0.050000000000000003, 0, 0.070661728551207226},
    {3, 0, 1, 0, 0, 0.070000000000000007, 0.13, -0.050000000000000003, 0.13700000000000001, 1.7125742704986509},
    {3, 0, 1, 0, 1, 0.070000000000000007, 0.13, -0.050000000000000003, 0.13700000000000001, -0.28459604605230487},
    {3, 0, 1, 0, 2, 0.070000000000000007, 0.13, -0.050000000000000003, 0.13700000000000001, 0.075346743961880597},
    {3, 0, 1, 0, 0, 0.20999999999999999, -0.29999999999999999, 0.25, 0, 1.1439973619060513},
    {3, 0, 1, 0, 1, 0.20999999999999999, -0.29999999999999999, 0.25, 0, 0.48463460992062346},
    {3, 0, 1, 0, 2, 0.20999999999999999, -0.29999999999999999, 0.25, 0, -0.0069111100383360225},
    {3, 0, 1, 0, 0, 0.20999999999999999, -0.29999999999999999, 0.25, 0.13700000000000001, 1.2198466990251609},
    {3, 0, 1, 0, 1, 0.20999999999999999, -0.29999999999999999, 0.25, 0.13700000000000001, 0.51676686400747895},
    {3, 0, 1, 0, 2, 0.20999999999999999, -0.29999999999999999, 0.25, 0.13700000000000001, -0.007369330601267759},
    {3, 0, 1, 0, 0, -0.28000000000000003, 0.23999999999999999, -0.33000000000000002, 0, 1.1971318824905379},
    {3, 0, 1, 0, 1, -0.28000000000000003, 0.23999999999999999, -0.33000000000000002, 0, -0.57913434499601157},
    {3, 0, 1, 0, 2, -0.28000000000000003, 0.23999999999999999, -0.33000000000000002, 0, 0.073980679589154733},
    {3, 0, 1, 0, 0, -0.28000000000000003, 0.23999999999999999, -0.33000000000000002, 0.13700000000000001, 1.276504145709547},
    {3, 0, 1, 0, 1, -0.28000000000000003, 0.23999999999999999, -0.33000000000000002, 0.13700000000000001, -0.61753212250283118},
    {3, 0, 1, 0, 2, -0.28000000000000003, 0.23999999999999999, -0.33000000000000002, 0.13700000000000001, 0.07888574816126169},
    {3, 0, 1, 1, 0, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, 0.76350259782717167},
    {3, 0, 1, 1, 1, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, -0.37798973845150796},
    {3, 0, 1, 1, 2, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, 0.17689074049886516},
    {3, 0, 1, 1, 0, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0.13700000000000001, 0.81412436310591452},
    {3, 0, 1, 1, 1, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0.13700000000000001, -0.40305122203010996},
    {3, 0, 1, 1, 2, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0.13700000000000001, 0.18861895409106505},
    {3, 0, 1, 1, 0, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 0.53399364695281371},
    {3, 0, 1, 1, 1, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 0.5569551880396304},
    {3, 0, 1, 1, 2, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 0.28377868500468351},
    {3, 0, 1, 1, 0, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0.13700000000000001, 0.56939850494977906},
    {3, 0, 1, 1, 1, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0.13700000000000001, 0.59388244261604817},
    {3, 0, 1, 1, 2, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0.13700000000000001, 0.30259378533872217},
    {3, 0, 1, 1, 0, -0.12, 0.64000000000000001, 0.70999999999999996, 0, 1.8841114985326231},
    {3, 0, 1, 1, 1, -0.12, 0.64000000000000001, 0.70999999999999996, 0, -0.55904210591106729},
    {3, 0, 1, 1, 2, -0.12, 0.64000000000000001, 0.70999999999999996, 0, -0.13725360203551967},
    {3, 0, 1, 1, 0, -0.12, 0.64000000000000001, 0.70999999999999996, 0.13700000000000001, 2.0090318986846718},
    {3, 0, 1, 1, 1, -0.12, 0.64000000000000001, 0.70999999999999996, 0.13700000000000001, -0.59610772736003326},
    {3, 0, 1, 1, 2, -0.12, 0.64000000000000001, 0.70999999999999996, 0.13700000000000001, -0.14635379324073261},
    {3, 0, 1, 2, 0, 0.62, 0.57999999999999996, -0.58999999999999997, 0, -0.29675265340522983},
    {3, 0, 1, 2, 1, 0.62, 0.57999999999999996, -0.58999999999999997, 0, -0.24907483335285865},
    {3, 0, 1, 2, 2, 0.62, 0.57999999999999996, -0.58999999999999997, 0, 0.091177900131546488},
    {3, 0, 1, 2, 0, 0.62, 0.57999999999999996, -0.58999999999999997, 0.13700000000000001, -0.31642795406468366},
    {3, 0, 1, 2, 1, 0.62, 0.57999999999999996, -0.58999999999999997, 0.13700000000000001, -0.26558899818571297},
    {3, 0, 1, 2, 2, 0.62, 0.57999999999999996, -0.58999999999999997, 0.13700000000000001, 0.097223179181287978},
    {3, 0, 1, 2, 0, -1.1000000000000001, 0.40000000000000002, 0.82999999999999996, 0, -1.1341594390759469},
    {3, 0, 1, 2, 1, -1.1000000000000001, 0.40000000000000002, 0.82999999999999996, 0, -0.42892132043578057},
    {3, 0, 1, 2, 2, -1.1000000000000001, 0.40000000000000002, 0.82999999999999996, 0, 0.081184031421307035},
    {3, 0, 1, 2, 0, -1.1000000000000001, 0.40000000000000002, 0.82999999999999996, 0.13700000000000001, -1.2093565020289263},
    {3, 0, 1, 2, 1, -1.1000000000000001, 0.40000000000000002, 0.82999999999999996, 0.13700000000000001, -0.45735967083293733},
    {3, 0, 1, 2, 2, -1.1000000000000001, 0.40000000000000002, 0.82999999999999996, 0.13700000000000001, 0.086566696777897953},
    {3, 0, 1, 2, 0, 1.47, -1.4099999999999999, 1.3799999999999999, 0, 0.26511330167994568},
    {3, 0, 1, 2, 1, 1.47, -1.4099999999999999, 1.3799999999999999, 0, -0.60644999902243846},
    {3, 0, 1, 2, 2, 1.47, -1.4099999999999999, 1.3799999999999999, 0, -0.33934484665162917},
    {3, 0, 1, 2, 0, 1.47, -1.4099999999999999, 1.3799999999999999, 0.13700000000000001, 0.28269084937671546},
    {3, 0, 1, 2, 1, 1.47, -1.4099999999999999, 1.3799999999999999, 0.13700000000000001, -0.64665885959629188},
    {3, 0, 1, 2, 2, 1.47, -1.4099999999999999, 1.3799999999999999, 0.13700000000000001, -0.36184409580236782},
    {3, 1, 1, 0, 0, 0, 0, 0, 0, 1.3495574287564276},
    {3, 1, 1, 0, 1, 0, 0, 0, 0, -0.18568082204692038},
    {3, 1, 1, 0, 2, 0, 0, 0, 0, 0.076426919130048482},
    {3, 1, 1, 0, 0, 0, 0, 0, 0.13700000000000001, 1.4384479421479928},
    {3, 1, 1, 0, 1, 0, 0, 0, 0.13700000000000001, -0.18139699297220202},
    {3, 1, 1, 0, 2, 0, 0, 0, 0.13700000000000001, 0.081461577952849404},
    {3, 1, 1, 0, 0, 0.070000000000000007, 0.13, -0.050000000000000003, 0, 1.3928264995911162},
    {3, 1, 1, 0, 1, 0.070000000000000007, 0.13, -0.050000000000000003, 0, -0.3750298429511672},
    {3, 1, 1, 0, 2, 0.070000000000000007, 0.13, -0.050000000000000003, 0, 0.092048971401609642},
    {3, 1, 1, 0, 0, 0.070000000000000007, 0.13, -0.050000000000000003, 0.13700000000000001, 1.4581467540684232},
    {3, 1, 1, 0, 1, 0.070000000000000007, 0.13, -0.050000000000000003, 0.13700000000000001, -0.38469318810279418},
    {3, 1, 1, 0, 2, 0.070000000000000007, 0.13, -0.050000000000000003, 0.13700000000000001, 0.096355531594338018},
    {3, 1, 1, 0, 0, 0.20999999999999999, -0.29999999999999999, 0.25, 0, 0.98849769166779466},
    {3, 1, 1, 0, 1, 0.20999999999999999, -0.29999999999999999, 0.25, 0, 0.5858296666426418},
    {3, 1, 1, 0, 2, 0.20999999999999999, -0.29999999999999999, 0.25, 0, -0.004214813384448019},
    {3, 1, 1, 0, 0, 0.20999999999999999, -0.29999999999999999, 0.25, 0.13700000000000001, 1.0096309055510675},
    {3, 1, 1, 0, 1, 0.20999999999999999, -0.29999999999999999, 0.25, 0.13700000000000001, 0.62303146988426772},
    {3, 1, 1, 0, 2, 0.20999999999999999, -0.29999999999999999, 0.25, 0.13700000000000001, 0.014790445866356982},
    {3, 1, 1, 0, 0, -0.28000000000000003, 0.23999999999999999, -0.33000000000000002, 0, 1.0349903971792205},
    {3, 1, 1, 0, 1, -0.28000000000000003, 0.23999999999999999, -0.33000000000000002, 0, -0.69604874607356804},
    {3, 1, 1, 0, 2, -0.28000000000000003, 0.23999999999999999, -0.33000000000000002, 0, 0.094640906118872978},
    {3, 1, 1, 0, 0, -0.28000000000000003, 0.23999999999999999, -0.33000000000000002, 0.13700000000000001, 1.1453306408862627},
    {3, 1, 1, 0, 1, -0.28000000000000003, 0.23999999999999999, -0.33000000000000002, 0.13700000000000001, -0.71661439241470792},
    {3, 1, 1, 0, 2, -0.28000000000000003, 0.23999999999999999, -0.33000000000000002, 0.13700000000000001, 0.074140180222710614},
    {3, 1, 1, 1, 0, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, 0.58467426110856935},
    {3, 1, 1, 1, 1, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, -0.5850422260644832},
    {3, 1, 1, 1, 2, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, 0.2317771310800319},
    {3, 1, 1, 1, 0, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0.13700000000000001, 0.60309935281839167},
    {3, 1, 1, 1, 1, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0.13700000000000001, -0.57739786009994121},
    {3, 1, 1, 1, 2, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0.13700000000000001, 0.21231485239347581},
    {3, 1, 1, 1, 0, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 0.69537719875515835},
    {3, 1, 1, 1, 1, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 0.67403106366909693},
    {3, 1, 1, 1, 2, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 0.35595479610720943},
    {3, 1, 1, 1, 0, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0.13700000000000001, 0.63619017990506388},
    {3, 1, 1, 1, 1, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0.13700000000000001, 0.76031388173452197},
    {3, 1, 1, 1, 2, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0.13700000000000001, 0.36651498909749497},
    {3, 1, 1, 1, 0, -0.12, 0.64000000000000001, 0.70999999999999996, 0, 1.6306694605687575},
    {3, 1, 1, 1, 1, -0.12, 0.64000000000000001, 0.70999999999999996, 0, -0.61046332435174644},
    {3, 1, 1, 1, 2, -0.12, 0.64000000000000001, 0.70999999999999996, 0, -0.065594455409323946},
    {3, 1, 1, 1, 0, -0.12, 0.64000000000000001, 0.70999999999999996, 0.13700000000000001, 1.7665051841529726},
    {3, 1, 1, 1, 1, -0.12, 0.64000000000000001, 0.70999999999999996, 0.13700000000000001, -0.65034958173595114},
    {3, 1, 1, 1, 2, -0.12, 0.64000000000000001, 0.70999999999999996, 0.13700000000000001, -0.065129576085041968},
    {3, 1, 1, 2, 0, 0.62, 0.57999999999999996, -0.58999999999999997, 0, -0.27215857172957603},
    {3, 1, 1, 2, 1, 0.62, 0.57999999999999996, -0.58999999999999997, 0, -0.40749595567690133},
    {3, 1, 1, 2, 2, 0.62, 0.57999999999999996, -0.58999999999999997, 0, 0.11190386684206198},
    {3, 1, 1, 2, 0, 0.62, 0.57999999999999996, -0.58999999999999997, 0.13700000000000001, -0.33910839618674254},
    {3, 1, 1, 2, 1, 0.62, 0.57999999999999996, -0.58999999999999997, 0.13700000000000001, -0.40478903907886704},
    {3, 1, 1, 2, 2, 0.62, 0.57999999999999996, -0.58999999999999997, 0.13700000000000001, 0.088290855006521299},
    {3, 1, 1, 2, 0, -1.1000000000000001, 0.40000000000000002, 0.82999999999999996, 0, -1.0048895091914534},
    {3, 1, 1, 2, 1, -1.1000000000000001, 0.40000000000000002, 0.82999999999999996, 0, -0.61405498506541112},
    {3, 1, 1, 2, 2, -1.1000000000000001, 0.40000000000000002, 0.82999999999999996, 0, 0.1015787085617427},
    {3, 1, 1, 2, 0, -1.1000000000000001, 0.40000000000000002, 0.82999999999999996, 0.13700000000000001, -0.92462734882450415},
    {3, 1, 1, 2, 1, -1.1000000000000001, 0.40000000000000002, 0.82999999999999996, 0.13700000000000001, -0.66515185493085593},
    {3, 1, 1, 2, 2, -1.1000000000000001, 0.40000000000000002, 0.82999999999999996, 0.13700000000000001, 0.067915262668835422},
    {3, 1, 1, 2, 0, 1.47, -1.4099999999999999, 1.3799999999999999, 0, 0.21947413896995249},
    {3, 1, 1, 2, 1, 1.47, -1.4099999999999999, 1.3799999999999999, 0, -0.74034015194130232},
    {3, 1, 1, 2, 2, 1.47, -1.4099999999999999, 1.3799999999999999, 0, -0.42895979553550556},
    {3, 1, 1, 2, 0, 1.47, -1.4099999999999999, 1.3799999999999999, 0.13700000000000001, -0.21200044955089617},
    {3, 1, 1, 2, 1, 1.47, -1.4099999999999999, 1.3799999999999999, 0.13700000000000001, -0.65746262604851891},
    {3, 1, 1, 2, 2, 1.47, -1.4099999999999999, 1.3799999999999999, 0.13700000000000001, -0.93160419842984865},
    {4, 0, 1, 1, 0, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, -0.71363580789936953},
    {4, 0, 1, 1, 1, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, -0.10498823702630272},
    {4, 0, 1, 1, 2, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, -0.23251168071969047},
    {4, 0, 1, 1, 0, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0.13700000000000001, -0.76095130422485213},
    {4, 0, 1, 1, 1, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0.13700000000000001, -0.11194916932293068},
    {4, 0, 1, 1, 2, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0.13700000000000001, -0.24792767505874638},
    {4, 0, 1, 1, 0, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 0.065951484696455245},
    {4, 0, 1, 1, 1, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, -0.036163694120820801},
    {4, 0, 1, 1, 2, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 0.051776125577351777},
    {4, 0, 1, 1, 0, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0.13700000000000001, 0.070324201420130739},
    {4, 0, 1, 1, 1, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0.13700000000000001, -0.038561420128048925},
    {4, 0, 1, 1, 2, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0.13700000000000001, 0.055208987342954714},
    {4, 0, 1, 1, 0, -0.12, 0.64000000000000001, 0.70999999999999996, 0, -1.7394466258744459},
    {4, 0, 1, 1, 1, -0.12, 0.64000000000000001, 0.70999999999999996, 0, -0.053672911550404746},
    {4, 0, 1, 1, 2, -0.12, 0.64000000000000001, 0.70999999999999996, 0, 0.33153670355802262},
    {4, 0, 1, 1, 0, -0.12, 0.64000000000000001, 0.70999999999999996, 0.13700000000000001, -1.854775452600782},
    {4, 0, 1, 1, 1, -0.12, 0.64000000000000001, 0.70999999999999996, 0.13700000000000001, -0.057231534059435617},
    {4, 0, 1, 1, 2, -0.12, 0.64000000000000001, 0.70999999999999996, 0.13700000000000001, 0.35351825704135653},
    {4, 1, 1, 1, 0, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, -0.53112395158162207},
    {4, 1, 1, 1, 1, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, -0.052617373293260702},
    {4, 1, 1, 1, 2, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, -0.38562174471485244},
    {4, 1, 1, 1, 0, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0.13700000000000001, -0.56633854297580788},
    {4, 1, 1, 1, 1, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0.13700000000000001, -0.056106011482594503},
    {4, 1, 1, 1, 2, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0.13700000000000001, -0.41118924573303933},
    {4, 1, 1, 1, 0, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, -0.23062773523385358},
    {4, 1, 1, 1, 1, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 0.028097347671348114},
    {4, 1, 1, 1, 2, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, -0.17093299837587389},
    {4, 1, 1, 1, 0, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0.13700000000000001, -0.24591882017973468},
    {4, 1, 1, 1, 1, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0.13700000000000001, 0.029960258606847222},
    {4, 1, 1, 1, 2, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0.13700000000000001, -0.18226620162469101},
    {4, 1, 1, 1, 0, -0.12, 0.64000000000000001, 0.70999999999999996, 0, -1.5741423369200636},
    {4, 1, 1, 1, 1, -0.12, 0.64000000000000001, 0.70999999999999996, 0, -0.072077917555205651},
    {4, 1, 1, 1, 2, -0.12, 0.64000000000000001, 0.70999999999999996, 0, 0.26084103114163143},
    {4, 1, 1, 1, 0, -0.12, 0.64000000000000001, 0.70999999999999996, 0.13700000000000001, -1.6785111552078789},
    {4, 1, 1, 1, 1, -0.12, 0.64000000000000001, 0.70999999999999996, 0.13700000000000001, -0.076856829158969603},
    {4, 1, 1, 1, 2, -0.12, 0.64000000000000001, 0.70999999999999996, 0.13700000000000001, 0.27813531866742958},
    {5, 0, 1, 1, 0, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, 0.75339334077047471},
    {5, 0, 1, 1, 0, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0.13700000000000001, 0.88525145001178174},
    {5, 0, 1, 1, 0, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 0.64326565772854927},
    {5, 0, 1, 1, 0, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0.13700000000000001, 0.75494373351741872},
    {5, 0, 1, 1, 0, -0.12, 0.64000000000000001, 0.70999999999999996, 0, 1.0244615513637556},
    {5, 0, 1, 1, 0, -0.12, 0.64000000000000001, 0.70999999999999996, 0.13700000000000001, 1.162880675629395},
    {5, 1, 1, 1, 0, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, 0.87421090460085693},
    {5, 1, 1, 1, 0, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0.13700000000000001, 1.0140794624970557},
    {5, 1, 1, 1, 0, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 0.64513668653167511},
    {5, 1, 1, 1, 0, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0.13700000000000001, 0.75693881531155094},
    {5, 1, 1, 1, 0, -0.12, 0.64000000000000001, 0.70999999999999996, 0, 0.96702557235050846},
    {5, 1, 1, 1, 0, -0.12, 0.64000000000000001, 0.70999999999999996, 0.13700000000000001, 1.1016365751291513},
    {6, 0, 1, 1, 0, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, 0.92906057005891196},
    {6, 0, 1, 1, 0, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0.13700000000000001, 0.81711571685040629},
    {6, 0, 1, 1, 0, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 1.012022173634135},
    {6, 0, 1, 1, 0, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0.13700000000000001, 0.89008106740029469},
    {6, 0, 1, 1, 0, -0.12, 0.64000000000000001, 0.70999999999999996, 0, 1.0014368378048242},
    {6, 0, 1, 1, 0, -0.12, 0.64000000000000001, 0.70999999999999996, 0.13700000000000001, 0.88077118540441879},
    {6, 1, 1, 1, 0, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, 0.92906057005891196},
    {6, 1, 1, 1, 0, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0.13700000000000001, 0.81711571685040629},
    {6, 1, 1, 1, 0, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 1.012022173634135},
    {6, 1, 1, 1, 0, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0.13700000000000001, 0.89008106740029469},
    {6, 1, 1, 1, 0, -0.12, 0.64000000000000001, 0.70999999999999996, 0, 1.0014368378048242},
    {6, 1, 1, 1, 0, -0.12, 0.64000000000000001, 0.70999999999999996, 0.13700000000000001, 0.88077118540441879},
    {7, 0, 1, 0, 0, 0, 0, 0, 0, 0.10847717934605706},
    {7, 0, 1, 0, 1, 0, 0, 0, 0, -0.026629437666952697},
    {7, 0, 1, 0, 2, 0, 0, 0, 0, -0.069961495493920681},
    {7, 0, 1, 0, 3, 0, 0, 0, 0, 0.023242664657922539},
    {7, 0, 1, 0, 4, 0, 0, 0, 0, 0.025684480776004458},
    {7, 0, 1, 0, 5, 0, 0, 0, 0, 0.067280155996020399},
    {7, 0, 1, 0, 0, 0.070000000000000007, 0.13, -0.050000000000000003, 0, 0.11523144244951063},
    {7, 0, 1, 0, 1, 0.070000000000000007, 0.13, -0.050000000000000003, 0, -0.023633645169464729},
    {7, 0, 1, 0, 2, 0.070000000000000007, 0.13, -0.050000000000000003, 0, -0.06935962409530716},
    {7, 0, 1, 0, 3, 0.070000000000000007, 0.13, -0.050000000000000003, 0, 0.020961766457106453},
    {7, 0, 1, 0, 4, 0.070000000000000007, 0.13, -0.050000000000000003, 0, 0.021909047597716929},
    {7, 0, 1, 0, 5, 0.070000000000000007, 0.13, -0.050000000000000003, 0, 0.062806791093382927},
    {7, 0, 1, 0, 0, 0.20999999999999999, -0.29999999999999999, 0.25, 0, 0.15163956984135579},
    {7, 0, 1, 0, 1, 0.20999999999999999, -0.29999999999999999, 0.25, 0, 0.0070430784163413446},
    {7, 0, 1, 0, 2, 0.20999999999999999, -0.29999999999999999, 0.25, 0, -0.049845400484316883},
    {7, 0, 1, 0, 3, 0.20999999999999999, -0.29999999999999999, 0.25, 0, 0.014550834674099244},
    {7, 0, 1, 0, 4, 0.20999999999999999, -0.29999999999999999, 0.25, 0, 0.002225613029571403},
    {7, 0, 1, 0, 5, 0.20999999999999999, -0.29999999999999999, 0.25, 0, 0.032809595484544964},
    {7, 0, 1, 0, 0, -0.28000000000000003, 0.23999999999999999, -0.33000000000000002, 0, 0.15665607371437085},
    {7, 0, 1, 0, 1, -0.28000000000000003, 0.23999999999999999, -0.33000000000000002, 0, 0.017228986251390712},
    {7, 0, 1, 0, 2, -0.28000000000000003, 0.23999999999999999, -0.33000000000000002, 0, -0.040483130781091219},
    {7, 0, 1, 0, 3, -0.28000000000000003, 0.23999999999999999, -0.33000000000000002, 0, 0.016080970646107182},
    {7, 0, 1, 0, 4, -0.28000000000000003, 0.23999999999999999, -0.33000000000000002, 0, -0.00021262922430756586},
    {7, 0, 1, 0, 5, -0.28000000000000003, 0.23999999999999999, -0.33000000000000002, 0, 0.026262955639521978},
    {7, 0, 1, 1, 0, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, 0.15738472957834365},
    {7, 0, 1, 1, 1, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, 0.019450757721302565},
    {7, 0, 1, 1, 2, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, -0.038292011629080386},
    {7, 0, 1, 1, 3, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, 0.016603840783215282},
    {7, 0, 1, 1, 4, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, -0.00053267868819488625},
    {7, 0, 1, 1, 5, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, 0.025011429638441073},
    {7, 0, 1, 1, 0, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 0.097201186639712764},
    {7, 0, 1, 1, 1, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 0.065832905238739109},
    {7, 0, 1, 1, 2, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 0.038181559989613846},
    {7, 0, 1, 1, 3, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 0.066523826271050968},
    {7, 0, 1, 1, 4, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 0.036466280788684788},
    {7, 0, 1, 1, 5, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 0.035274987089236293},
    {7, 0, 1, 1, 0, -0.12, 0.64000000000000001, 0.70999999999999996, 0, 0.02837550659096431},
    {7, 0, 1, 1, 1, -0.12, 0.64000000000000001, 0.70999999999999996, 0, 0.028658266129057097},
    {7, 0, 1, 1, 2, -0.12, 0.64000000000000001, 0.70999999999999996, 0, 0.024603952904978122},
    {7, 0, 1, 1, 3, -0.12, 0.64000000000000001, 0.70999999999999996, 0, 0.087073767876423389},
    {7, 0, 1, 1, 4, -0.12, 0.64000000000000001, 0.70999999999999996, 0, 0.074632313391501143},
    {7, 0, 1, 1, 5, -0.12, 0.64000000000000001, 0.70999999999999996, 0, 0.083550725532612308},
    {7, 0, 1, 2, 0, 0.62, 0.57999999999999996, -0.58999999999999997, 0, 0.018470491549359005},
    {7, 0, 1, 2, 1, 0.62, 0.57999999999999996, -0.58999999999999997, 0, -0.00071602532766771004},
    {7, 0, 1, 2, 2, 0.62, 0.57999999999999996, -0.58999999999999997, 0, -0.0042544679056134146},
    {7, 0, 1, 2, 3, 0.62, 0.57999999999999996, -0.58999999999999997, 0, 0.080301233509256567},
    {7, 0, 1, 2, 4, 0.62, 0.57999999999999996, -0.58999999999999997, 0, 0.07902094457240956},
    {7, 0, 1, 2, 5, 0.62, 0.57999999999999996, -0.58999999999999997, 0, 0.10022827494138443},
    {7, 0, 1, 2, 0, -1.1000000000000001, 0.40000000000000002, 0.82999999999999996, 0, 0.12335562030724605},
    {7, 0, 1, 2, 1, -1.1000000000000001, 0.40000000000000002, 0.82999999999999996, 0, 0.065181733726939195},
    {7, 0, 1, 2, 2, -1.1000000000000001, 0.40000000000000002, 0.82999999999999996, 0, 0.026791679859805254},
    {7, 0, 1, 2, 3, -1.1000000000000001, 0.40000000000000002, 0.82999999999999996, 0, 0.052729524222980169},
    {7, 0, 1, 2, 4, -1.1000000000000001, 0.40000000000000002, 0.82999999999999996, 0, 0.021283709379620375},
    {7, 0, 1, 2, 5, -1.1000000000000001, 0.40000000000000002, 0.82999999999999996, 0, 0.022914855469773791},
    {7, 0, 1, 2, 0, 1.47, -1.4099999999999999, 1.3799999999999999, 0, 0.12295262253232631},
    {7, 0, 1, 2, 1, 1.47, -1.4099999999999999, 1.3799999999999999, 0, -0.01951008914522142},
    {7, 0, 1, 2, 2, 1.47, -1.4099999999999999, 1.3799999999999999, 0, -0.067888909530602407},
    {7, 0, 1, 2, 3, 1.47, -1.4099999999999999, 1.3799999999999999, 0, 0.018637399191462612},
    {7, 0, 1, 2, 4, 1.47, -1.4099999999999999, 1.3799999999999999, 0, 0.017625253568512397},
    {7, 0, 1, 2, 5, 1.47, -1.4099999999999999, 1.3799999999999999, 0, 0.057409978276211096},
    {7, 1, 1, 0, 0, 0, 0, 0, 0, 0.10847717934605706},
    {7, 1, 1, 0, 1, 0, 0, 0, 0, -0.026629437666952697},
    {7, 1, 1, 0, 2, 0, 0, 0, 0, -0.069961495493920681},
    {7, 1, 1, 0, 3, 0, 0, 0, 0, 0.023242664657922539},
    {7, 1, 1, 0, 4, 0, 0, 0, 0, 0.025684480776004458},
    {7, 1, 1, 0, 5, 0, 0, 0, 0, 0.067280155996020399},
    {7, 1, 1, 0, 0, 0.070000000000000007, 0.13, -0.050000000000000003, 0, 0.11523144244951063},
    {7, 1, 1, 0, 1, 0.070000000000000007, 0.13, -0.050000000000000003, 0, -0.023633645169464729},
    {7, 1, 1, 0, 2, 0.070000000000000007, 0.13, -0.050000000000000003, 0, -0.06935962409530716},
    {7, 1, 1, 0, 3, 0.070000000000000007, 0.13, -0.050000000000000003, 0, 0.020961766457106453},
    {7, 1, 1, 0, 4, 0.070000000000000007, 0.13, -0.050000000000000003, 0, 0.021909047597716929},
    {7, 1, 1, 0, 5, 0.070000000000000007, 0.13, -0.050000000000000003, 0, 0.062806791093382927},
    {7, 1, 1, 0, 0, 0.20999999999999999, -0.29999999999999999, 0.25, 0, 0.15163956984135579},
    {7, 1, 1, 0, 1, 0.20999999999999999, -0.29999999999999999, 0.25, 0, 0.0070430784163413446},
    {7, 1, 1, 0, 2, 0.20999999999999999, -0.29999999999999999, 0.25, 0, -0.049845400484316883},
    {7, 1, 1, 0, 3, 0.20999999999999999, -0.29999999999999999, 0.25, 0, 0.014550834674099244},
    {7, 1, 1, 0, 4, 0.20999999999999999, -0.29999999999999999, 0.25, 0, 0.002225613029571403},
    {7, 1, 1, 0, 5, 0.20999999999999999, -0.29999999999999999, 0.25, 0, 0.032809595484544964},
    {7, 1, 1, 0, 0, -0.28000000000000003, 0.23999999999999999, -0.33000000000000002, 0, 0.15665607371437085},
    {7, 1, 1, 0, 1, -0.28000000000000003, 0.23999999999999999, -0.33000000000000002, 0, 0.017228986251390712},
    {7, 1, 1, 0, 2, -0.28000000000000003, 0.23999999999999999, -0.33000000000000002, 0, -0.040483130781091219},
    {7, 1, 1, 0, 3, -0.28000000000000003, 0.23999999999999999, -0.33000000000000002, 0, 0.016080970646107182},
    {7, 1, 1, 0, 4, -0.28000000000000003, 0.23999999999999999, -0.33000000000000002, 0, -0.00021262922430756586},
    {7, 1, 1, 0, 5, -0.28000000000000003, 0.23999999999999999, -0.33000000000000002, 0, 0.026262955639521978},
    {7, 1, 1, 1, 0, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, 0.15738472957834365},
    {7, 1, 1, 1, 1, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, 0.019450757721302565},
    {7, 1, 1, 1, 2, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, -0.038292011629080386},
    {7, 1, 1, 1, 3, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, 0.016603840783215282},
    {7, 1, 1, 1, 4, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, -0.00053267868819488625},
    {7, 1, 1, 1, 5, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, 0.025011429638441073},
    {7, 1, 1, 1, 0, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 0.097201186639712764},
    {7, 1, 1, 1, 1, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 0.065832905238739109},
    {7, 1, 1, 1, 2, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 0.038181559989613846},
    {7, 1, 1, 1, 3, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 0.066523826271050968},
    {7, 1, 1, 1, 4, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 0.036466280788684788},
    {7, 1, 1, 1, 5, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 0.035274987089236293},
    {7, 1, 1, 1, 0, -0.12, 0.64000000000000001, 0.70999999999999996, 0, 0.02837550659096431},
    {7, 1, 1, 1, 1, -0.12, 0.64000000000000001, 0.70999999999999996, 0, 0.028658266129057097},
    {7, 1, 1, 1, 2, -0.12, 0.64000000000000001, 0.70999999999999996, 0, 0.024603952904978122},
    {7, 1, 1, 1, 3, -0.12, 0.64000000000000001, 0.70999999999999996, 0, 0.087073767876423389},
    {7, 1, 1, 1, 4, -0.12, 0.64000000000000001, 0.70999999999999996, 0, 0.074632313391501143},
    {7, 1, 1, 1, 5, -0.12, 0.64000000000000001, 0.70999999999999996, 0, 0.083550725532612308},
    {7, 1, 1, 2, 0, 0.62, 0.57999999999999996, -0.58999999999999997, 0, 0.018470491549359005},
    {7, 1, 1, 2, 1, 0.62, 0.57999999999999996, -0.58999999999999997, 0, -0.00071602532766771004},
    {7, 1, 1, 2, 2, 0.62, 0.57999999999999996, -0.58999999999999997, 0, -0.0042544679056134146},
    {7, 1, 1, 2, 3, 0.62, 0.57999999999999996, -0.58999999999999997, 0, 0.080301233509256567},
    {7, 1, 1, 2, 4, 0.62, 0.57999999999999996, -0.58999999999999997, 0, 0.07902094457240956},
    {7, 1, 1, 2, 5, 0.62, 0.57999999999999996, -0.58999999999999997, 0, 0.10022827494138443},
    {7, 1, 1, 2, 0, -1.1000000000000001, 0.40000000000000002, 0.82999999999999996, 0, 0.12335562030724605},
    {7, 1, 1, 2, 1, -1.1000000000000001, 0.40000000000000002, 0.82999999999999996, 0, 0.065181733726939195},
    {7, 1, 1, 2, 2, -1.1000000000000001, 0.40000000000000002, 0.82999999999999996, 0, 0.026791679859805254},
    {7, 1, 1, 2, 3, -1.1000000000000001, 0.40000000000000002, 0.82999999999999996, 0, 0.052729524222980169},
    {7, 1, 1, 2, 4, -1.1000000000000001, 0.40000000000000002, 0.82999999999999996, 0, 0.021283709379620375},
    {7, 1, 1, 2, 5, -1.1000000000000001, 0.40000000000000002, 0.82999999999999996, 0, 0.022914855469773791},
    {7, 1, 1, 2, 0, 1.47, -1.4099999999999999, 1.3799999999999999, 0, 0.12295262253232631},
    {7, 1, 1, 2, 1, 1.47, -1.4099999999999999, 1.3799999999999999, 0, -0.01951008914522142},
    {7, 1, 1, 2, 2, 1.47, -1.4099999999999999, 1.3799999999999999, 0, -0.067888909530602407},
    {7, 1, 1, 2, 3, 1.47, -1.4099999999999999, 1.3799999999999999, 0, 0.018637399191462612},
    {7, 1, 1, 2, 4, 1.47, -1.4099999999999999, 1.3799999999999999, 0, 0.017625253568512397},
    {7, 1, 1, 2, 5, 1.47, -1.4099999999999999, 1.3799999999999999, 0, 0.057409978276211096},
    {8, 0, 1, 0, 0, 0, 0, 0, 0, 0.41558195887766869},
    {8, 0, 1, 0, 1, 0, 0, 0, 0, -0.086253459128365462},
    {8, 0, 1, 0, 2, 0, 0, 0, 0, -0.23337816089643781},
    {8, 0, 1, 0, 3, 0, 0, 0, 0, 0.12876382195137942},
    {8, 0, 1, 0, 4, 0, 0, 0, 0, 0.089766961280588975},
    {8, 0, 1, 0, 5, 0, 0, 0, 0, 0.2756542191709519},
    {8, 0, 1, 0, 0, 0.070000000000000007, 0.13, -0.050000000000000003, 0, 0.43794706849175336},
    {8, 0, 1, 0, 1, 0.070000000000000007, 0.13, -0.050000000000000003, 0, -0.076333616421451658},
    {8, 0, 1, 0, 2, 0.070000000000000007, 0.13, -0.050000000000000003, 0, -0.23138520924539963},
    {8, 0, 1, 0, 3, 0.070000000000000007, 0.13, -0.050000000000000003, 0, 0.12121117890231953},
    {8, 0, 1, 0, 4, 0.070000000000000007, 0.13, -0.050000000000000003, 0, 0.077265526915398494},
    {8, 0, 1, 0, 5, 0.070000000000000007, 0.13, -0.050000000000000003, 0, 0.26084175260592712},
    {8, 0, 1, 0, 0, 0.20999999999999999, -0.29999999999999999, 0.25, 0, 0.55850378170978365},
    {8, 0, 1, 0, 1, 0.20999999999999999, -0.29999999999999999, 0.25, 0, 0.025244938498435995},
    {8, 0, 1, 0, 2, 0.20999999999999999, -0.29999999999999999, 0.25, 0, -0.16676857477192189},
    {8, 0, 1, 0, 3, 0.20999999999999999, -0.29999999999999999, 0.25, 0, 0.099982927965209564},
    {8, 0, 1, 0, 4, 0.20999999999999999, -0.29999999999999999, 0.25, 0, 0.012088591259287485},
    {8, 0, 1, 0, 5, 0.20999999999999999, -0.29999999999999999, 0.25, 0, 0.16151329032500675},
    {8, 0, 1, 0, 0, -0.28000000000000003, 0.23999999999999999, -0.33000000000000002, 0, 0.57511472168665467},
    {8, 0, 1, 0, 1, -0.28000000000000003, 0.23999999999999999, -0.33000000000000002, 0, 0.058973110137672306},
    {8, 0, 1, 0, 2, -0.28000000000000003, 0.23999999999999999, -0.33000000000000002, 0, -0.13576768171488326},
    {8, 0, 1, 0, 3, -0.28000000000000003, 0.23999999999999999, -0.33000000000000002, 0, 0.10504960336920936},
    {8, 0, 1, 0, 4, -0.28000000000000003, 0.23999999999999999, -0.33000000000000002, 0, 0.0040149414120061325},
    {8, 0, 1, 0, 5, -0.28000000000000003, 0.23999999999999999, -0.33000000000000002, 0, 0.13983567494413593},
    {8, 0, 1, 1, 0, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, 0.57752748944815402},
    {8, 0, 1, 1, 1, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, 0.066329969309565862},
    {8, 0, 1, 1, 2, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, -0.12851232028438384},
    {8, 0, 1, 1, 3, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, 0.10678096143910373},
    {8, 0, 1, 1, 4, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, 0.002955174975293151},
    {8, 0, 1, 1, 5, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, 0.1356915491127422},
    {8, 0, 1, 1, 0, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 0.3782442346844756},
    {8, 0, 1, 1, 1, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 0.21991323923485243},
    {8, 0, 1, 1, 2, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 0.12471142679738513},
    {8, 0, 1, 1, 3, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 0.27207892663061262},
    {8, 0, 1, 1, 4, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 0.12546828582588809},
    {8, 0, 1, 1, 5, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 0.16967683868491179},
    {8, 0, 1, 1, 0, -0.12, 0.64000000000000001, 0.70999999999999996, 0, 0.15034463187405031},
    {8, 0, 1, 1, 1, -0.12, 0.64000000000000001, 0.70999999999999996, 0, 0.096818407745839136},
    {8, 0, 1, 1, 2, -0.12, 0.64000000000000001, 0.70999999999999996, 0, 0.079752462940975447},
    {8, 0, 1, 1, 3, -0.12, 0.64000000000000001, 0.70999999999999996, 0, 0.3401250908868127},
    {8, 0, 1, 1, 4, -0.12, 0.64000000000000001, 0.70999999999999996, 0, 0.25184587722594226},
    {8, 0, 1, 1, 5, -0.12, 0.64000000000000001, 0.70999999999999996, 0, 0.32953027723913697},
    {8, 0, 1, 2, 0, 0.62, 0.57999999999999996, -0.58999999999999997, 0, 0.11754656882237711},
    {8, 0, 1, 2, 1, 0.62, 0.57999999999999996, -0.58999999999999997, 0, -0.00044745800490524366},
    {8, 0, 1, 2, 2, 0.62, 0.57999999999999996, -0.58999999999999997, 0, -0.015805221862307785},
    {8, 0, 1, 2, 3, 0.62, 0.57999999999999996, -0.58999999999999997, 0, 0.31769948039950535},
    {8, 0, 1, 2, 4, 0.62, 0.57999999999999996, -0.58999999999999997, 0, 0.26637776855345352},
    {8, 0, 1, 2, 5, 0.62, 0.57999999999999996, -0.58999999999999997, 0, 0.38475395077811753},
    {8, 0, 1, 2, 0, -1.1000000000000001, 0.40000000000000002, 0.82999999999999996, 0, 0.46484831967630769},
    {8, 0, 1, 2, 1, -1.1000000000000001, 0.40000000000000002, 0.82999999999999996, 0, 0.2177570421759123},
    {8, 0, 1, 2, 2, -1.1000000000000001, 0.40000000000000002, 0.82999999999999996, 0, 0.086996591930469255},
    {8, 0, 1, 2, 3, -1.1000000000000001, 0.40000000000000002, 0.82999999999999996, 0, 0.22640242978269606},
    {8, 0, 1, 2, 4, -1.1000000000000001, 0.40000000000000002, 0.82999999999999996, 0, 0.0751948705640854},
    {8, 0, 1, 2, 5, -1.1000000000000001, 0.40000000000000002, 0.82999999999999996, 0, 0.12874925054099623},
    {8, 0, 1, 2, 0, 1.47, -1.4099999999999999, 1.3799999999999999, 0, 0.46351388995803039},
    {8, 0, 1, 2, 1, 1.47, -1.4099999999999999, 1.3799999999999999, 0, -0.062679457400778446},
    {8, 0, 1, 2, 2, 1.47, -1.4099999999999999, 1.3799999999999999, 0, -0.2265152934682316},
    {8, 0, 1, 2, 3, 1.47, -1.4099999999999999, 1.3799999999999999, 0, 0.11351459855250548},
    {8, 0, 1, 2, 4, 1.47, -1.4099999999999999, 1.3799999999999999, 0, 0.063080778474323876},
    {8, 0, 1, 2, 5, 1.47, -1.4099999999999999, 1.3799999999999999, 0, 0.24297151148946414},
    {8, 1, 1, 0, 0, 0, 0, 0, 0, 0.41558195887766869},
    {8, 1, 1, 0, 1, 0, 0, 0, 0, -0.086253459128365462},
    {8, 1, 1, 0, 2, 0, 0, 0, 0, -0.23337816089643781},
    {8, 1, 1, 0, 3, 0, 0, 0, 0, 0.12876382195137942},
    {8, 1, 1, 0, 4, 0, 0, 0, 0, 0.089766961280588975},
    {8, 1, 1, 0, 5, 0, 0, 0, 0, 0.2756542191709519},
    {8, 1, 1, 0, 0, 0.070000000000000007, 0.13, -0.050000000000000003, 0, 0.43794706849175336},
    {8, 1, 1, 0, 1, 0.070000000000000007, 0.13, -0.050000000000000003, 0, -0.076333616421451658},
    {8, 1, 1, 0, 2, 0.070000000000000007, 0.13, -0.050000000000000003, 0, -0.23138520924539963},
    {8, 1, 1, 0, 3, 0.070000000000000007, 0.13, -0.050000000000000003, 0, 0.12121117890231953},
    {8, 1, 1, 0, 4, 0.070000000000000007, 0.13, -0.050000000000000003, 0, 0.077265526915398494},
    {8, 1, 1, 0, 5, 0.070000000000000007, 0.13, -0.050000000000000003, 0, 0.26084175260592712},
    {8, 1, 1, 0, 0, 0.20999999999999999, -0.29999999999999999, 0.25, 0, 0.55850378170978365},
    {8, 1, 1, 0, 1, 0.20999999999999999, -0.29999999999999999, 0.25, 0, 0.025244938498435995},
    {8, 1, 1, 0, 2, 0.20999999999999999, -0.29999999999999999, 0.25, 0, -0.16676857477192189},
    {8, 1, 1, 0, 3, 0.20999999999999999, -0.29999999999999999, 0.25, 0, 0.099982927965209564},
    {8, 1, 1, 0, 4, 0.20999999999999999, -0.29999999999999999, 0.25, 0, 0.012088591259287485},
    {8, 1, 1, 0, 5, 0.20999999999999999, -0.29999999999999999, 0.25, 0, 0.16151329032500675},
    {8, 1, 1, 0, 0, -0.28000000000000003, 0.23999999999999999, -0.33000000000000002, 0, 0.57511472168665467},
    {8, 1, 1, 0, 1, -0.28000000000000003, 0.23999999999999999, -0.33000000000000002, 0, 0.058973110137672306},
    {8, 1, 1, 0, 2, -0.28000000000000003, 0.23999999999999999, -0.33000000000000002, 0, -0.13576768171488326},
    {8, 1, 1, 0, 3, -0.28000000000000003, 0.23999999999999999, -0.33000000000000002, 0, 0.10504960336920936},
    {8, 1, 1, 0, 4, -0.28000000000000003, 0.23999999999999999, -0.33000000000000002, 0, 0.0040149414120061325},
    {8, 1, 1, 0, 5, -0.28000000000000003, 0.23999999999999999, -0.33000000000000002, 0, 0.13983567494413593},
    {8, 1, 1, 1, 0, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, 0.57752748944815402},
    {8, 1, 1, 1, 1, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, 0.066329969309565862},
    {8, 1, 1, 1, 2, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, -0.12851232028438384},
    {8, 1, 1, 1, 3, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, 0.10678096143910373},
    {8, 1, 1, 1, 4, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, 0.002955174975293151},
    {8, 1, 1, 1, 5, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, 0.1356915491127422},
    {8, 1, 1, 1, 0, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 0.3782442346844756},
    {8, 1, 1, 1, 1, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 0.21991323923485243},
    {8, 1, 1, 1, 2, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 0.12471142679738513},
    {8, 1, 1, 1, 3, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 0.27207892663061262},
    {8, 1, 1, 1, 4, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 0.12546828582588809},
    {8, 1, 1, 1, 5, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 0.16967683868491179},
    {8, 1, 1, 1, 0, -0.12, 0.64000000000000001, 0.70999999999999996, 0, 0.15034463187405031},
    {8, 1, 1, 1, 1, -0.12, 0.64000000000000001, 0.70999999999999996, 0, 0.096818407745839136},
    {8, 1, 1, 1, 2, -0.12, 0.64000000000000001, 0.70999999999999996, 0, 0.079752462940975447},
    {8, 1, 1, 1, 3, -0.12, 0.64000000000000001, 0.70999999999999996, 0, 0.3401250908868127},
    {8, 1, 1, 1, 4, -0.12, 0.64000000000000001, 0.70999999999999996, 0, 0.25184587722594226},
    {8, 1, 1, 1, 5, -0.12, 0.64000000000000001, 0.70999999999999996, 0, 0.32953027723913697},
    {8, 1, 1, 2, 0, 0.62, 0.57999999999999996, -0.58999999999999997, 0, 0.11754656882237711},
    {8, 1, 1, 2, 1, 0.62, 0.57999999999999996, -0.58999999999999997, 0, -0.00044745800490524366},
    {8, 1, 1, 2, 2, 0.62, 0.57999999999999996, -0.58999999999999997, 0, -0.015805221862307785},
    {8, 1, 1, 2, 3, 0.62, 0.57999999999999996, -0.58999999999999997, 0, 0.31769948039950535},
    {8, 1, 1, 2, 4, 0.62, 0.57999999999999996, -0.58999999999999997, 0, 0.26637776855345352},
    {8, 1, 1, 2, 5, 0.62, 0.57999999999999996, -0.58999999999999997, 0, 0.38475395077811753},
    {8, 1, 1, 2, 0, -1.1000000000000001, 0.40000000000000002, 0.82999999999999996, 0, 0.46484831967630769},
    {8, 1, 1, 2, 1, -1.1000000000000001, 0.40000000000000002, 0.82999999999999996, 0, 0.2177570421759123},
    {8, 1, 1, 2, 2, -1.1000000000000001, 0.40000000000000002, 0.82999999999999996, 0, 0.086996591930469255},
    {8, 1, 1, 2, 3, -1.1000000000000001, 0.40000000000000002, 0.82999999999999996, 0, 0.22640242978269606},
    {8, 1, 1, 2, 4, -1.1000000000000001, 0.40000000000000002, 0.82999999999999996, 0, 0.0751948705640854},
    {8, 1, 1, 2, 5, -1.1000000000000001, 0.40000000000000002, 0.82999999999999996, 0, 0.12874925054099623},
    {8, 1, 1, 2, 0, 1.47, -1.4099999999999999, 1.3799999999999999, 0, 0.46351388995803039},
    {8, 1, 1, 2, 1, 1.47, -1.4099999999999999, 1.3799999999999999, 0, -0.062679457400778446},
    {8, 1, 1, 2, 2, 1.47, -1.4099999999999999, 1.3799999999999999, 0, -0.2265152934682316},
    {8, 1, 1, 2, 3, 1.47, -1.4099999999999999, 1.3799999999999999, 0, 0.11351459855250548},
    {8, 1, 1, 2, 4, 1.47, -1.4099999999999999, 1.3799999999999999, 0, 0.063080778474323876},
    {8, 1, 1, 2, 5, 1.47, -1.4099999999999999, 1.3799999999999999, 0, 0.24297151148946414},
    {9, 0, 1, 0, 0, 0, 0, 0, 0, 0.62384031969540732},
    {9, 0, 1, 0, 0, 0.070000000000000007, 0.13, -0.050000000000000003, 0, 0.64138618790564406},
    {9, 0, 1, 0, 0, 0.20999999999999999, -0.29999999999999999, 0.25, 0, 0.67233608581589133},
    {9, 0, 1, 0, 0, -0.28000000000000003, 0.23999999999999999, -0.33000000000000002, 0, 0.56868159301619181},
    {9, 0, 1, 2, 0, 0.62, 0.57999999999999996, -0.58999999999999997, 0, 0.24458237802427851},
    {9, 0, 1, 2, 0, -1.1000000000000001, 0.40000000000000002, 0.82999999999999996, 0, 0.25220799276500722},
    {9, 0, 1, 2, 0, 1.47, -1.4099999999999999, 1.3799999999999999, 0, 0.25742286979843099},
    {9, 1, 1, 0, 0, 0, 0, 0, 0, 0.62384031969540732},
    {9, 1, 1, 0, 0, 0.070000000000000007, 0.13, -0.050000000000000003, 0, 0.64138618790564406},
    {9, 1, 1, 0, 0, 0.20999999999999999, -0.29999999999999999, 0.25, 0, 0.67233608581589133},
    {9, 1, 1, 0, 0, -0.28000000000000003, 0.23999999999999999, -0.33000000000000002, 0, 0.56868159301619181},
    {9, 1, 1, 2, 0, 0.62, 0.57999999999999996, -0.58999999999999997, 0, 0.24458237802427851},
    {9, 1, 1, 2, 0, -1.1000000000000001, 0.40000000000000002, 0.82999999999999996, 0, 0.25220799276500722},
    {9, 1, 1, 2, 0, 1.47, -1.4099999999999999, 1.3799999999999999, 0, 0.25742286979843099},
    {10, 0, 1, 1, 0, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, -0.36682944027802533},
    {10, 0, 1, 1, 0, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0.13700000000000001, -0.39115097353270362},
    {10, 0, 1, 1, 0, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, -0.15680035073312398},
    {10, 0, 1, 1, 0, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0.13700000000000001, -0.16719653088107092},
    {10, 0, 1, 1, 0, -0.12, 0.64000000000000001, 0.70999999999999996, 0, 0.98040474129199973},
    {10, 0, 1, 1, 0, -0.12, 0.64000000000000001, 0.70999999999999996, 0.13700000000000001, 1.0454075570428436},
    {10, 1, 1, 1, 0, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, -0.23817981956298079},
    {10, 1, 1, 1, 0, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0.13700000000000001, -0.32902575867419959},
    {10, 1, 1, 1, 0, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, -0.14535140070803546},
    {10, 1, 1, 1, 0, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0.13700000000000001, -0.1193871049971969},
    {10, 1, 1, 1, 0, -0.12, 0.64000000000000001, 0.70999999999999996, 0, 1.2533778679077823},
    {10, 1, 1, 1, 0, -0.12, 0.64000000000000001, 0.70999999999999996, 0.13700000000000001, 1.2918376546728996},
    {11, 0, 1, 0, 0, 0, 0, 0, 0, -0.2907185550267159},
    {11, 0, 1, 0, 0, 0, 0, 0, 0.13700000000000001, -0.30999378276872941},
    {11, 0, 1, 0, 0, 0.070000000000000007, 0.13, -0.050000000000000003, 0, -0.79594228434634873},
    {11, 0, 1, 0, 0, 0.070000000000000007, 0.13, -0.050000000000000003, 0.13700000000000001, -0.84871486640209159},
    {11, 0, 1, 0, 0, 0.20999999999999999, -0.29999999999999999, 0.25, 0, -1.8741163078799035},
    {11, 0, 1, 0, 0, 0.20999999999999999, -0.29999999999999999, 0.25, 0.13700000000000001, -1.9983740066913431},
    {11, 0, 1, 0, 0, -0.28000000000000003, 0.23999999999999999, -0.33000000000000002, 0, 1.6677918766445881},
    {11, 0, 1, 0, 0, -0.28000000000000003, 0.23999999999999999, -0.33000000000000002, 0.13700000000000001, 1.7783698486823563},
    {11, 0, 1, 1, 0, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, -0.35943045032403798},
    {11, 0, 1, 1, 0, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0.13700000000000001, -0.38326141559136895},
    {11, 0, 1, 1, 0, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, -1.451892829091354},
    {11, 0, 1, 1, 0, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0.13700000000000001, -1.548156257943222},
    {11, 0, 1, 1, 0, -0.12, 0.64000000000000001, 0.70999999999999996, 0, -1.3471504257263971},
    {11, 0, 1, 1, 0, -0.12, 0.64000000000000001, 0.70999999999999996, 0.13700000000000001, -1.4364692215502157},
    {11, 0, 1, 2, 0, 0.62, 0.57999999999999996, -0.58999999999999997, 0, -0.99711047355113613},
    {11, 0, 1, 2, 0, 0.62, 0.57999999999999996, -0.58999999999999997, 0.13700000000000001, -1.0632209131131343},
    {11, 0, 1, 2, 0, -1.1000000000000001, 0.40000000000000002, 0.82999999999999996, 0, -0.29888761772040673},
    {11, 0, 1, 2, 0, -1.1000000000000001, 0.40000000000000002, 0.82999999999999996, 0.13700000000000001, -0.31870447082873099},
    {11, 0, 1, 2, 0, 1.47, -1.4099999999999999, 1.3799999999999999, 0, 0.61018548637174286},
    {11, 0, 1, 2, 0, 1.47, -1.4099999999999999, 1.3799999999999999, 0.13700000000000001, 0.65064201730629501},
    {11, 1, 1, 0, 0, 0, 0, 0, 0, -0.65584932635366522},
    {11, 1, 1, 0, 0, 0, 0, 0, 0.13700000000000001, -0.82848152830189059},
    {11, 1, 1, 0, 0, 0.070000000000000007, 0.13, -0.050000000000000003, 0, -1.0087208439057602},
    {11, 1, 1, 0, 0, 0.070000000000000007, 0.13, -0.050000000000000003, 0.13700000000000001, -1.3620017313000909},
    {11, 1, 1, 0, 0, 0.20999999999999999, -0.29999999999999999, 0.25, 0, -1.9887015354450022},
    {11, 1, 1, 0, 0, 0.20999999999999999, -0.29999999999999999, 0.25, 0.13700000000000001, -1.9953010827410307},
    {11, 1, 1, 0, 0, -0.28000000000000003, 0.23999999999999999, -0.33000000000000002, 0, 1.2746127219168226},
    {11, 1, 1, 0, 0, -0.28000000000000003, 0.23999999999999999, -0.33000000000000002, 0.13700000000000001, 1.5460826838698802},
    {11, 1, 1, 1, 0, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, 0.40034280493899133},
    {11, 1, 1, 1, 0, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0.13700000000000001, 0.082543873488736666},
    {11, 1, 1, 1, 0, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, -1.52264472589912},
    {11, 1, 1, 1, 0, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0.13700000000000001, -1.4530587939700503},
    {11, 1, 1, 1, 0, -0.12, 0.64000000000000001, 0.70999999999999996, 0, -0.4606709550215593},
    {11, 1, 1, 1, 0, -0.12, 0.64000000000000001, 0.70999999999999996, 0.13700000000000001, -0.70959342045945795},
    {11, 1, 1, 2, 0, 0.62, 0.57999999999999996, -0.58999999999999997, 0, -0.79235412131931393},
    {11, 1, 1, 2, 0, 0.62, 0.57999999999999996, -0.58999999999999997, 0.13700000000000001, -0.8111704012601173},
    {11, 1, 1, 2, 0, -1.1000000000000001, 0.40000000000000002, 0.82999999999999996, 0, -0.25906848605643545},
    {11, 1, 1, 2, 0, -1.1000000000000001, 0.40000000000000002, 0.82999999999999996, 0.13700000000000001, -0.43634063502175019},
    {11, 1, 1, 2, 0, 1.47, -1.4099999999999999, 1.3799999999999999, 0, 0.68240521806754273},
    {11, 1, 1, 2, 0, 1.47, -1.4099999999999999, 1.3799999999999999, 0.13700000000000001, 1.3502997819366973},
    {12, 0, 1, 1, 0, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, 0.39313905126111881},
    {12, 0, 1, 1, 1, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, 0.010803739285366052},
    {12, 0, 1, 1, 2, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, -0.075233937603922996},
    {12, 0, 1, 1, 0, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0.13700000000000001, 0.41920496489583958},
    {12, 0, 1, 1, 1, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0.13700000000000001, 0.011520049034400242},
    {12, 0, 1, 1, 2, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0.13700000000000001, -0.080222099715250184},
    {12, 0, 1, 1, 0, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 0.35985228253011037},
    {12, 0, 1, 1, 1, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 0.30457326753686137},
    {12, 0, 1, 1, 2, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 0.18462599249514822},
    {12, 0, 1, 1, 0, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0.13700000000000001, 0.38371121612522902},
    {12, 0, 1, 1, 1, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0.13700000000000001, 0.32476709071874504},
    {12, 0, 1, 1, 2, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0.13700000000000001, 0.19686706892768702},
    {12, 0, 1, 1, 0, -0.12, 0.64000000000000001, 0.70999999999999996, 0, 0.21819417028950902},
    {12, 0, 1, 1, 1, -0.12, 0.64000000000000001, 0.70999999999999996, 0, -0.042294325586057231},
    {12, 0, 1, 1, 2, -0.12, 0.64000000000000001, 0.70999999999999996, 0, -0.035759134636857093},
    {12, 0, 1, 1, 0, -0.12, 0.64000000000000001, 0.70999999999999996, 0.13700000000000001, 0.23266088475127936},
    {12, 0, 1, 1, 1, -0.12, 0.64000000000000001, 0.70999999999999996, 0.13700000000000001, -0.045098524849469246},
    {12, 0, 1, 1, 2, -0.12, 0.64000000000000001, 0.70999999999999996, 0.13700000000000001, -0.038130037532681557},
    {12, 1, 1, 1, 0, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, 0.26907340835937665},
    {12, 1, 1, 1, 1, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, -0.023005003608601537},
    {12, 1, 1, 1, 2, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, -0.045416550086418858},
    {12, 1, 1, 1, 0, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0.13700000000000001, 0.28272259846925207},
    {12, 1, 1, 1, 1, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0.13700000000000001, -0.021024109532613626},
    {12, 1, 1, 1, 2, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0.13700000000000001, -0.050402677678742985},
    {12, 1, 1, 1, 0, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 0.45564240141446594},
    {12, 1, 1, 1, 1, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 0.38097333868615774},
    {12, 1, 1, 1, 2, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 0.23168828928928994},
    {12, 1, 1, 1, 0, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0.13700000000000001, 0.45354696351240958},
    {12, 1, 1, 1, 1, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0.13700000000000001, 0.39275803548856614},
    {12, 1, 1, 1, 2, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0.13700000000000001, 0.23692456916200097},
    {12, 1, 1, 1, 0, -0.12, 0.64000000000000001, 0.70999999999999996, 0, 0.18082699333233263},
    {12, 1, 1, 1, 1, -0.12, 0.64000000000000001, 0.70999999999999996, 0, -0.066274766109889194},
    {12, 1, 1, 1, 2, -0.12, 0.64000000000000001, 0.70999999999999996, 0, -0.045308124789024976},
    {12, 1, 1, 1, 0, -0.12, 0.64000000000000001, 0.70999999999999996, 0.13700000000000001, 0.19742451653014909},
    {12, 1, 1, 1, 1, -0.12, 0.64000000000000001, 0.70999999999999996, 0.13700000000000001, -0.066572606614170998},
    {12, 1, 1, 1, 2, -0.12, 0.64000000000000001, 0.70999999999999996, 0.13700000000000001, -0.044366888945385666},
    {13, 0, 1, 1, 0, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, 1.1585424265442228},
    {13, 0, 1, 1, 0, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 0.97638240328397263},
    {13, 0, 1, 1, 0, -0.12, 0.64000000000000001, 0.70999999999999996, 0, 0.99713250605483494},
    {13, 1, 1, 1, 0, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, 1.1585424265442228},
    {13, 1, 1, 1, 0, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 0.97638240328397263},
    {13, 1, 1, 1, 0, -0.12, 0.64000000000000001, 0.70999999999999996, 0, 0.99713250605483494},
};

inline const label nSelfTestSamples = 592;

} // End namespace shellTransmissionMMS
} // End namespace Foam

#endif
