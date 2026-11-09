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
    lead-field conductor transmission manufactured solution on the spherical
    shell (bulk conductor Gi + Ge in the heart, source Vm sealed at the
    spheres), the forward problem of the shell lead-field ECG test
    (manufacturedShellLeadFieldECGVerifier).

    Regenerate from generator/:
        python shell_mms.py --write ../../../src/verificationModels/ecgVerification/manufacturedShellLeadFieldECGReference.H

    sympy 1.14.0
\*---------------------------------------------------------------------------*/

#ifndef manufacturedShellLeadFieldECGReference_H
#define manufacturedShellLeadFieldECGReference_H

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
    const scalar theta0 = p.theta0;
    const scalar theta1 = p.theta1;
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
    const scalar w5 = hA1*cos(w4);
    const scalar w6 = fA1*y;
    const scalar w7 = (1.0/2.0)*y + 1;
    const scalar w8 = k*x;
    const scalar w9 = fA2*sin(w8);
    const scalar w10 = w4 + 3.0/10.0;
    const scalar w11 = fA3*exp((2.0/5.0)*z);
    const scalar w12 = w11*cos(w10);
    const scalar w13 = r1 + r2;
    const scalar w14 = (1.0/(r2));
    const scalar w15 = cos(gamma2);
    const scalar w16 = sin(gamma2);
    const scalar w17 = (w16*w16);
    const scalar w18 = (r1*r1);
    const scalar w19 = -(r2*r2);
    const scalar w20 = w2 + (y*y) + (z*z);
    const scalar w21 = (-w18 + w20)/(-w18 - w19);
    const scalar w22 = theta0 + theta1*w21;
    const scalar w23 = sin(w22);
    const scalar w24 = cos(w22);
    const scalar w25 = gil*w17*(w23*w23) + gin*(w15*w15) + git*w17*(w24*w24);
    const scalar w26 = w14*z;
    const scalar w27 = cos(gamma1);
    const scalar w28 = gin*w15*w16;
    const scalar w29 = sin(gamma1);
    const scalar w30 = w24*w29;
    const scalar w31 = w23*w27;
    const scalar w32 = w15*w31 + w30;
    const scalar w33 = w16*w23;
    const scalar w34 = w23*w29;
    const scalar w35 = w15*w24*w27 - w34;
    const scalar w36 = w16*w24;
    const scalar w37 = gil*w32*w33 + git*w35*w36 - w27*w28;
    const scalar w38 = w14*y;
    const scalar w39 = -w15*w34 + w24*w27;
    const scalar w40 = gil*w39;
    const scalar w41 = -w15*w30 - w31;
    const scalar w42 = git*w41;
    const scalar w43 = w28*w29 + w33*w40 + w36*w42;
    const scalar w44 = w14*x;
    const scalar w45 = gin*w17;
    const scalar w46 = gil*(w32*w32) + git*(w35*w35) + (w27*w27)*w45;
    const scalar w47 = -w27*w29*w45 + w32*w40 + w35*w42;
    const scalar w48 = gil*(w39*w39) + git*(w41*w41) + (w29*w29)*w45;
    const scalar w49 = St*(hA0 + hA2*w8*cos(w0)) - St*((2.0/5.0)*w12 + w6);
    const scalar w50 = -St*(fA0 + fA2*k*w7*cos(w8)) + St*(hA2*w1 - 3.0/5.0*w5*x);
    const scalar w51 = -St*hA1*k*w3*sin(w4) - St*(fA1*z - k*w11*sin(w10) + (1.0/2.0)*w9);
    const scalar w52 = w25*w49 + w37*w51 + w43*w50;
    const scalar w53 = w37*w49 + w46*w51 + w47*w50;
    const scalar w54 = w43*w49 + w47*w51 + w48*w50;
    const scalar w55 = (1.0/2.0)*w13;
    const scalar w56 = (1.0/(r1));
    const scalar w57 = w56*z;
    const scalar w58 = w56*y;
    const scalar w59 = w56*x;
    return St*(hA0*z + hA2*w1*x + w3*w5) - St*(fA0*x + w12 + w6*z + w7*w9) + w21*(w19 + w20)*(w14*w21*w55*(-w26*w52 - w38*w53 - w44*w54)/(w26*(w25*w26 + w37*w38 + w43*w44) + w38*(w26*w37 + w38*w46 + w44*w47) + w44*(w26*w43 + w38*w47 + w44*w48)) + w55*w56*(1 - w21)*(w52*w57 + w53*w58 + w54*w59)/(-w57*(-w25*w57 - w37*w58 - w43*w59) - w58*(-w37*w57 - w46*w58 - w47*w59) - w59*(-w43*w57 - w47*w58 - w48*w59)))/w13;
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
        const scalar w6 = sigmaEps*((1.0/2.0)*sin(k*z) + 1)*sin(w0 + 1.0/5.0)*cos(w2) + 1;
        const scalar w7 = cos(gamma1);
        const scalar w8 = cos(gamma2);
        const scalar w9 = sin(gamma2);
        const scalar w10 = w8*w9;
        const scalar w11 = gen*w10*w7;
        const scalar w12 = gin*w10;
        const scalar w13 = w12*w7;
        const scalar w14 = sin(gamma1);
        const scalar w15 = cos(theta0);
        const scalar w16 = w14*w15;
        const scalar w17 = sin(theta0);
        const scalar w18 = w17*w7;
        const scalar w19 = w16 + w18*w8;
        const scalar w20 = w17*w9;
        const scalar w21 = w19*w20;
        const scalar w22 = gel*w21;
        const scalar w23 = gil*w21;
        const scalar w24 = w14*w17;
        const scalar w25 = w15*w7*w8 - w24;
        const scalar w26 = w15*w9;
        const scalar w27 = w25*w26;
        const scalar w28 = get*w27;
        const scalar w29 = git*w27;
        const scalar w30 = (w8*w8);
        const scalar w31 = (w9*w9);
        const scalar w32 = (w17*w17)*w31;
        const scalar w33 = (w15*w15)*w31;
        const scalar w34 = (1.0/5.0)*z;
        const scalar w35 = gen*w14;
        const scalar w36 = w15*w7 - w24*w8;
        const scalar w37 = w20*w36;
        const scalar w38 = -w16*w8 - w18;
        const scalar w39 = w26*w38;
        const scalar w40 = gel*w37 + get*w39 + gil*w37 + git*w39 + w10*w35 + w12*w14;
        const scalar w41 = (1.0/5.0)*k*w1;
        const scalar w42 = (1.0/10.0)*w31*(w7*w7);
        const scalar w43 = (1.0/10.0)*(w19*w19);
        const scalar w44 = (1.0/10.0)*(w25*w25);
        const scalar w45 = w31*w7;
        const scalar w46 = w35*w45;
        const scalar w47 = gin*w14*w45;
        const scalar w48 = w19*w36;
        const scalar w49 = gel*w48;
        const scalar w50 = gil*w48;
        const scalar w51 = w25*w38;
        const scalar w52 = get*w51;
        const scalar w53 = git*w51;
        const scalar w54 = (w14*w14)*w31;
        const scalar w55 = (w36*w36);
        const scalar w56 = (w38*w38);
        return St*(fA0*x + fA1*y*z + fA2*w1*((1.0/2.0)*y + 1) + fA3*exp((2.0/5.0)*z)*cos(w2 + 3.0/10.0)) + Tt*(w3 + (1.0/10.0)*w4 + (1.0/5.0)*cos(w0) + (1.0/2.0)*w5*((r1*r1) - w4 - (x*x) - (y*y))*(-sigma0In*w6*((1.0/5.0)*k*w1*w5*x - w3*w5 - 1.0/5.0*w4*w5) - w5*x*(w34*w40 - w41*(gel*w55 + gen*w54 + get*w56 + gil*w55 + gin*w54 + git*w56) - 1.0/10.0*w46 - 1.0/10.0*w47 + (1.0/10.0)*w49 + (1.0/10.0)*w50 + (1.0/10.0)*w52 + (1.0/10.0)*w53) - w5*y*(gel*w43 + gen*w42 + get*w44 + gil*w43 + gin*w42 + git*w44 + w34*(-w11 - w13 + w22 + w23 + w28 + w29) - w41*(-w46 - w47 + w49 + w50 + w52 + w53)) - w5*z*(-1.0/10.0*w11 - 1.0/10.0*w13 + (1.0/10.0)*w22 + (1.0/10.0)*w23 + (1.0/10.0)*w28 + (1.0/10.0)*w29 + w34*(gel*w32 + gen*w30 + get*w33 + gil*w32 + gin*w30 + git*w33) - w40*w41))/(sigma0In*w6));
    }
    if (zone == bath)
    {
        const scalar w0 = k*x;
        const scalar w1 = sin(w0);
        const scalar w2 = k*y;
        const scalar w3 = (1.0/10.0)*y;
        const scalar w4 = (z*z);
        const scalar w5 = (1.0/(r2));
        const scalar w6 = sigmaEps*((1.0/2.0)*sin(k*z) + 1)*sin(w0 + 1.0/5.0)*cos(w2) + 1;
        const scalar w7 = (1.0/5.0)*w5;
        const scalar w8 = cos(gamma1);
        const scalar w9 = cos(gamma2);
        const scalar w10 = sin(gamma2);
        const scalar w11 = w10*w9;
        const scalar w12 = gen*w11*w8;
        const scalar w13 = gin*w11;
        const scalar w14 = w13*w8;
        const scalar w15 = sin(gamma1);
        const scalar w16 = theta0 + theta1;
        const scalar w17 = cos(w16);
        const scalar w18 = w15*w17;
        const scalar w19 = sin(w16);
        const scalar w20 = w19*w8;
        const scalar w21 = w18 + w20*w9;
        const scalar w22 = w10*w19;
        const scalar w23 = w21*w22;
        const scalar w24 = gel*w23;
        const scalar w25 = gil*w23;
        const scalar w26 = w15*w19;
        const scalar w27 = w17*w8*w9 - w26;
        const scalar w28 = w10*w17;
        const scalar w29 = w27*w28;
        const scalar w30 = get*w29;
        const scalar w31 = git*w29;
        const scalar w32 = (w9*w9);
        const scalar w33 = (w10*w10);
        const scalar w34 = (w19*w19)*w33;
        const scalar w35 = (w17*w17)*w33;
        const scalar w36 = (1.0/5.0)*z;
        const scalar w37 = gen*w15;
        const scalar w38 = w17*w8 - w26*w9;
        const scalar w39 = w22*w38;
        const scalar w40 = -w18*w9 - w20;
        const scalar w41 = w28*w40;
        const scalar w42 = gel*w39 + get*w41 + gil*w39 + git*w41 + w11*w37 + w13*w15;
        const scalar w43 = (1.0/5.0)*k*w1;
        const scalar w44 = (1.0/10.0)*w33*(w8*w8);
        const scalar w45 = (1.0/10.0)*(w21*w21);
        const scalar w46 = (1.0/10.0)*(w27*w27);
        const scalar w47 = w33*w8;
        const scalar w48 = w37*w47;
        const scalar w49 = gin*w15*w47;
        const scalar w50 = w21*w38;
        const scalar w51 = gel*w50;
        const scalar w52 = gil*w50;
        const scalar w53 = w27*w40;
        const scalar w54 = get*w53;
        const scalar w55 = git*w53;
        const scalar w56 = (w15*w15)*w33;
        const scalar w57 = (w38*w38);
        const scalar w58 = (w40*w40);
        return St*(fA0*x + fA1*y*z + fA2*w1*((1.0/2.0)*y + 1) + fA3*exp((2.0/5.0)*z)*cos(w2 + 3.0/10.0)) + Tt*(w3 + (1.0/10.0)*w4 + (1.0/5.0)*cos(w0) + (1.0/2.0)*w5*(-(r2*r2) + w4 + (x*x) + (y*y))*(-sigma0Out*w6*(-w0*w1*w7 + w3*w5 + w4*w7) + w5*x*(w36*w42 - w43*(gel*w57 + gen*w56 + get*w58 + gil*w57 + gin*w56 + git*w58) - 1.0/10.0*w48 - 1.0/10.0*w49 + (1.0/10.0)*w51 + (1.0/10.0)*w52 + (1.0/10.0)*w54 + (1.0/10.0)*w55) + w5*y*(gel*w45 + gen*w44 + get*w46 + gil*w45 + gin*w44 + git*w46 + w36*(-w12 - w14 + w24 + w25 + w30 + w31) - w43*(-w48 - w49 + w51 + w52 + w54 + w55)) + w5*z*(-1.0/10.0*w12 - 1.0/10.0*w14 + (1.0/10.0)*w24 + (1.0/10.0)*w25 + (1.0/10.0)*w30 + (1.0/10.0)*w31 + w36*(gel*w34 + gen*w32 + get*w35 + gil*w34 + gin*w32 + git*w35) - w42*w43))/(sigma0Out*w6));
    }
    const scalar w0 = (z*z);
    const scalar w1 = k*x;
    const scalar w2 = cos(w1);
    const scalar w3 = fA1*y;
    const scalar w4 = (1.0/2.0)*y + 1;
    const scalar w5 = fA2*sin(w1);
    const scalar w6 = k*y;
    const scalar w7 = w6 + 3.0/10.0;
    const scalar w8 = fA3*exp((2.0/5.0)*z);
    const scalar w9 = w8*cos(w7);
    const scalar w10 = r1 + r2;
    const scalar w11 = (1.0/(r2));
    const scalar w12 = cos(gamma2);
    const scalar w13 = (w12*w12);
    const scalar w14 = sin(gamma2);
    const scalar w15 = (w14*w14);
    const scalar w16 = (r1*r1);
    const scalar w17 = -(r2*r2);
    const scalar w18 = w0 + (x*x) + (y*y);
    const scalar w19 = (-w16 + w18)/(-w16 - w17);
    const scalar w20 = theta0 + theta1*w19;
    const scalar w21 = sin(w20);
    const scalar w22 = w15*(w21*w21);
    const scalar w23 = cos(w20);
    const scalar w24 = w15*(w23*w23);
    const scalar w25 = gel*w22 + gen*w13 + get*w24 + gil*w22 + gin*w13 + git*w24;
    const scalar w26 = w11*z;
    const scalar w27 = cos(gamma1);
    const scalar w28 = w12*w14;
    const scalar w29 = w27*w28;
    const scalar w30 = sin(gamma1);
    const scalar w31 = w23*w30;
    const scalar w32 = w21*w27;
    const scalar w33 = w12*w32 + w31;
    const scalar w34 = w14*w21;
    const scalar w35 = w33*w34;
    const scalar w36 = w21*w30;
    const scalar w37 = w12*w23*w27 - w36;
    const scalar w38 = w14*w23;
    const scalar w39 = w37*w38;
    const scalar w40 = gel*w35 - gen*w29 + get*w39 + gil*w35 - gin*w29 + git*w39;
    const scalar w41 = w11*y;
    const scalar w42 = gen*w30;
    const scalar w43 = gin*w30;
    const scalar w44 = -w12*w36 + w23*w27;
    const scalar w45 = w34*w44;
    const scalar w46 = -w12*w31 - w32;
    const scalar w47 = w38*w46;
    const scalar w48 = gel*w45 + get*w47 + gil*w45 + git*w47 + w28*w42 + w28*w43;
    const scalar w49 = w11*x;
    const scalar w50 = w15*(w27*w27);
    const scalar w51 = (w33*w33);
    const scalar w52 = (w37*w37);
    const scalar w53 = gel*w51 + gen*w50 + get*w52 + gil*w51 + gin*w50 + git*w52;
    const scalar w54 = w15*w27;
    const scalar w55 = w33*w44;
    const scalar w56 = w37*w46;
    const scalar w57 = gel*w55 + get*w56 + gil*w55 + git*w56 - w42*w54 - w43*w54;
    const scalar w58 = w15*(w30*w30);
    const scalar w59 = (w44*w44);
    const scalar w60 = (w46*w46);
    const scalar w61 = gel*w59 + gen*w58 + get*w60 + gil*w59 + gin*w58 + git*w60;
    const scalar w62 = sigmaEps*((1.0/2.0)*sin(k*z) + 1)*sin(w1 + 1.0/5.0)*cos(w6) + 1;
    const scalar w63 = -sigma0Out*w62;
    const scalar w64 = St*(w3 + (2.0/5.0)*w9);
    const scalar w65 = -w48;
    const scalar w66 = St*(fA0 + fA2*k*w2*w4);
    const scalar w67 = -w40;
    const scalar w68 = St*(fA1*z - k*w8*sin(w7) + (1.0/2.0)*w5);
    const scalar w69 = w65*w66 + w67*w68;
    const scalar w70 = -w57;
    const scalar w71 = w64*w67 + w66*w70;
    const scalar w72 = w64*w65 + w68*w70;
    const scalar w73 = (1.0/2.0)*w10;
    const scalar w74 = (1.0/(r1));
    const scalar w75 = w74*z;
    const scalar w76 = w74*y;
    const scalar w77 = w74*x;
    const scalar w78 = -sigma0In*w62;
    return St*(fA0*x + w3*z + w4*w5 + w9) + Tt*((1.0/10.0)*w0 + (1.0/5.0)*w2 + (1.0/10.0)*y) + w19*(w17 + w18)*(w11*w19*w73*(w26*(w64*(-w25 - w63) + w69) + w41*(w68*(-w53 - w63) + w71) + w49*(w66*(-w61 - w63) + w72))/(w26*(w25*w26 + w40*w41 + w48*w49) + w41*(w26*w40 + w41*w53 + w49*w57) + w49*(w26*w48 + w41*w57 + w49*w61)) + w73*w74*(1 - w19)*(-w75*(w64*(-w25 - w78) + w69) - w76*(w68*(-w53 - w78) + w71) - w77*(w66*(-w61 - w78) + w72))/(-w75*(-w25*w75 - w40*w76 - w48*w77) - w76*(-w40*w75 - w53*w76 - w57*w77) - w77*(-w48*w75 - w57*w76 - w61*w77)))/w10;
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

    const scalar w0 = (z*z);
    const scalar w1 = k*x;
    const scalar w2 = cos(w1);
    const scalar w3 = k*z;
    const scalar w4 = sin(w3);
    const scalar w5 = (x*x);
    const scalar w6 = 1 - 3.0/10.0*w5;
    const scalar w7 = k*y;
    const scalar w8 = cos(w7);
    const scalar w9 = hA1*w8;
    const scalar w10 = cos(gamma2);
    const scalar w11 = (w10*w10);
    const scalar w12 = sin(gamma2);
    const scalar w13 = (w12*w12);
    const scalar w14 = (r1*r1);
    const scalar w15 = -(r2*r2);
    const scalar w16 = w0 + w5 + (y*y);
    const scalar w17 = (-w14 + w16)/(-w14 - w15);
    const scalar w18 = theta0 + theta1*w17;
    const scalar w19 = sin(w18);
    const scalar w20 = w13*(w19*w19);
    const scalar w21 = cos(w18);
    const scalar w22 = w13*(w21*w21);
    const scalar w23 = gil*w20 + gin*w11 + git*w22;
    const scalar w24 = (1.0/(r2));
    const scalar w25 = w24*z;
    const scalar w26 = cos(gamma1);
    const scalar w27 = w10*w12;
    const scalar w28 = gin*w27;
    const scalar w29 = sin(gamma1);
    const scalar w30 = w21*w29;
    const scalar w31 = w19*w26;
    const scalar w32 = w10*w31 + w30;
    const scalar w33 = w12*w19;
    const scalar w34 = gil*w33;
    const scalar w35 = w19*w29;
    const scalar w36 = w10*w21*w26 - w35;
    const scalar w37 = w12*w21;
    const scalar w38 = git*w37;
    const scalar w39 = -w26*w28 + w32*w34 + w36*w38;
    const scalar w40 = w24*y;
    const scalar w41 = -w10*w35 + w21*w26;
    const scalar w42 = -w10*w30 - w31;
    const scalar w43 = w28*w29 + w34*w41 + w38*w42;
    const scalar w44 = w24*x;
    const scalar w45 = (w26*w26);
    const scalar w46 = gin*w13;
    const scalar w47 = (w32*w32);
    const scalar w48 = (w36*w36);
    const scalar w49 = gil*w47 + git*w48 + w45*w46;
    const scalar w50 = w26*w29;
    const scalar w51 = w32*w41;
    const scalar w52 = w36*w42;
    const scalar w53 = gil*w51 + git*w52 - w46*w50;
    const scalar w54 = (w29*w29);
    const scalar w55 = (w41*w41);
    const scalar w56 = (w42*w42);
    const scalar w57 = gil*w55 + git*w56 + w46*w54;
    const scalar w58 = w7 + 3.0/10.0;
    const scalar w59 = fA3*exp((2.0/5.0)*z);
    const scalar w60 = St*(fA1*y + (2.0/5.0)*w59*cos(w58));
    const scalar w61 = St*(hA0 + hA2*w1*cos(w3)) - w60;
    const scalar w62 = St*(fA0 + fA2*k*w2*((1.0/2.0)*y + 1));
    const scalar w63 = St*(hA2*w4 - 3.0/5.0*w9*x) - w62;
    const scalar w64 = St*(fA1*z + (1.0/2.0)*fA2*sin(w1) - k*w59*sin(w58));
    const scalar w65 = -St*hA1*k*w6*sin(w7) - w64;
    const scalar w66 = w23*w61 + w39*w65 + w43*w63;
    const scalar w67 = w39*w61 + w49*w65 + w53*w63;
    const scalar w68 = w43*w61 + w53*w65 + w57*w63;
    const scalar w69 = r1 + r2;
    const scalar w70 = (1.0/2.0)*w69;
    const scalar w71 = w17*w24*w70;
    const scalar w72 = (1.0/(r1));
    const scalar w73 = w72*z;
    const scalar w74 = w72*y;
    const scalar w75 = w72*x;
    const scalar w76 = w70*w72*(1 - w17);
    const scalar w77 = w17*(w15 + w16)/w69;
    const scalar w78 = gel*w20 + gen*w11 + get*w22 + w23;
    const scalar w79 = gel*w33;
    const scalar w80 = get*w37;
    const scalar w81 = gen*w27;
    const scalar w82 = -w26*w81 + w32*w79 + w36*w80 + w39;
    const scalar w83 = w29*w81 + w41*w79 + w42*w80 + w43;
    const scalar w84 = gen*w13;
    const scalar w85 = gel*w47 + get*w48 + w45*w84 + w49;
    const scalar w86 = gel*w51 + get*w52 - w50*w84 + w53;
    const scalar w87 = gel*w55 + get*w56 + w54*w84 + w57;
    const scalar w88 = sigmaEps*w8*((1.0/2.0)*w4 + 1)*sin(w1 + 1.0/5.0) + 1;
    const scalar w89 = -sigma0Out*w88;
    const scalar w90 = -w83;
    const scalar w91 = -w82;
    const scalar w92 = w62*w90 + w64*w91;
    const scalar w93 = -w86;
    const scalar w94 = w60*w91 + w62*w93;
    const scalar w95 = w60*w90 + w64*w93;
    const scalar w96 = -sigma0In*w88;
    return St*(hA0*z + hA2*w4*x + w6*w9) + Tt*((1.0/10.0)*w0 + (1.0/5.0)*w2 + (1.0/10.0)*y) + w77*(w71*(-w25*w66 - w40*w67 - w44*w68)/(w25*(w23*w25 + w39*w40 + w43*w44) + w40*(w25*w39 + w40*w49 + w44*w53) + w44*(w25*w43 + w40*w53 + w44*w57)) + w76*(w66*w73 + w67*w74 + w68*w75)/(-w73*(-w23*w73 - w39*w74 - w43*w75) - w74*(-w39*w73 - w49*w74 - w53*w75) - w75*(-w43*w73 - w53*w74 - w57*w75))) + w77*(w71*(w25*(w60*(-w78 - w89) + w92) + w40*(w64*(-w85 - w89) + w94) + w44*(w62*(-w87 - w89) + w95))/(w25*(w25*w78 + w40*w82 + w44*w83) + w40*(w25*w82 + w40*w85 + w44*w86) + w44*(w25*w83 + w40*w86 + w44*w87)) + w76*(-w73*(w60*(-w78 - w96) + w92) - w74*(w64*(-w85 - w96) + w94) - w75*(w62*(-w87 - w96) + w95))/(-w73*(-w73*w78 - w74*w82 - w75*w83) - w74*(-w73*w82 - w74*w85 - w75*w86) - w75*(-w73*w83 - w74*w86 - w75*w87)));
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
        const scalar w5 = (1.0/(sigma0In));
        const scalar w6 = (z*z);
        const scalar w7 = (r1*r1) - w6 - (x*x) - (y*y);
        const scalar w8 = w0 + 1.0/5.0;
        const scalar w9 = sin(w8);
        const scalar w10 = k*z;
        const scalar w11 = (1.0/2.0)*sin(w10) + 1;
        const scalar w12 = k*y;
        const scalar w13 = cos(w12);
        const scalar w14 = sigmaEps*w11*w13;
        const scalar w15 = w14*w9 + 1;
        const scalar w16 = (1.0/(w15));
        const scalar w17 = (1.0/5.0)*w4;
        const scalar w18 = (k*k)*w1;
        const scalar w19 = w17*w18;
        const scalar w20 = w19*x;
        const scalar w21 = sigma0In*w15;
        const scalar w22 = w4*y;
        const scalar w23 = (1.0/5.0)*k*w2*w4*x - w17*w6 - 1.0/10.0*w22;
        const scalar w24 = k*sigma0In*w23;
        const scalar w25 = w14*cos(w8);
        const scalar w26 = sin(gamma1);
        const scalar w27 = sin(gamma2);
        const scalar w28 = (w27*w27);
        const scalar w29 = (w26*w26)*w28;
        const scalar w30 = cos(gamma1);
        const scalar w31 = cos(theta0);
        const scalar w32 = cos(gamma2);
        const scalar w33 = sin(theta0);
        const scalar w34 = w26*w33;
        const scalar w35 = w30*w31 - w32*w34;
        const scalar w36 = (w35*w35);
        const scalar w37 = w30*w33;
        const scalar w38 = w26*w31;
        const scalar w39 = -w32*w38 - w37;
        const scalar w40 = (w39*w39);
        const scalar w41 = gel*w36 + gen*w29 + get*w40 + gil*w36 + gin*w29 + git*w40;
        const scalar w42 = gen*w26;
        const scalar w43 = w27*w32;
        const scalar w44 = w42*w43;
        const scalar w45 = gin*w26;
        const scalar w46 = w43*w45;
        const scalar w47 = w27*w33;
        const scalar w48 = w35*w47;
        const scalar w49 = gel*w48;
        const scalar w50 = gil*w48;
        const scalar w51 = w27*w31;
        const scalar w52 = w39*w51;
        const scalar w53 = get*w52;
        const scalar w54 = git*w52;
        const scalar w55 = w44 + w46 + w49 + w50 + w53 + w54;
        const scalar w56 = (1.0/5.0)*z;
        const scalar w57 = w55*w56;
        const scalar w58 = w28*w30;
        const scalar w59 = w42*w58;
        const scalar w60 = w45*w58;
        const scalar w61 = w32*w37 + w38;
        const scalar w62 = w35*w61;
        const scalar w63 = gel*w62;
        const scalar w64 = gil*w62;
        const scalar w65 = w30*w31*w32 - w34;
        const scalar w66 = w39*w65;
        const scalar w67 = get*w66;
        const scalar w68 = git*w66;
        const scalar w69 = -w59 - w60 + w63 + w64 + w67 + w68;
        const scalar w70 = w4*(-w3*w41 + w57 - 1.0/10.0*w59 - 1.0/10.0*w60 + (1.0/10.0)*w63 + (1.0/10.0)*w64 + (1.0/10.0)*w67 + (1.0/10.0)*w68);
        const scalar w71 = w4*x;
        const scalar w72 = w30*w43;
        const scalar w73 = gen*w72;
        const scalar w74 = gin*w72;
        const scalar w75 = w47*w61;
        const scalar w76 = gel*w75;
        const scalar w77 = gil*w75;
        const scalar w78 = w51*w65;
        const scalar w79 = get*w78;
        const scalar w80 = git*w78;
        const scalar w81 = (w32*w32);
        const scalar w82 = gen*w81;
        const scalar w83 = gin*w81;
        const scalar w84 = w28*(w33*w33);
        const scalar w85 = gel*w84;
        const scalar w86 = w28*(w31*w31);
        const scalar w87 = get*w86;
        const scalar w88 = gil*w84;
        const scalar w89 = git*w86;
        const scalar w90 = w4*(-w3*w55 + w56*(w82 + w83 + w85 + w87 + w88 + w89) - 1.0/10.0*w73 - 1.0/10.0*w74 + (1.0/10.0)*w76 + (1.0/10.0)*w77 + (1.0/10.0)*w79 + (1.0/10.0)*w80);
        const scalar w91 = (1.0/10.0)*w28*(w30*w30);
        const scalar w92 = (1.0/10.0)*(w61*w61);
        const scalar w93 = (1.0/10.0)*(w65*w65);
        const scalar w94 = w4*(gel*w92 + gen*w91 + get*w93 + gil*w92 + gin*w91 + git*w93 - w3*w69 + w56*(-w73 - w74 + w76 + w77 + w79 + w80));
        const scalar w95 = -w21*w23 - w70*x - w90*z - w94*y;
        const scalar w96 = w16*w5;
        const scalar w97 = w95*w96;
        const scalar w98 = w4*w7;
        const scalar w99 = (1.0/2.0)*w98;
        const scalar w100 = k*w5*w95/(w15*w15);
        const scalar w101 = w100*w99;
        const scalar w102 = w12 + 3.0/10.0;
        const scalar w103 = fA3*exp((2.0/5.0)*z);
        const scalar w104 = sigmaEps*w9;
        const scalar w105 = w104*w24;
        const scalar w106 = w11*sin(w12);
        const scalar w107 = w13*cos(w10);
        const scalar w108 = w4*z;
        return vector(St*(fA0 + fA2*k*w1*((1.0/2.0)*y + 1)) + Tt*(-w101*w25 + (1.0/2.0)*w16*w4*w5*w7*(w18*w4*w57 + w19*w69*y + w20*w41 - w21*(w20 + w3*w4) - w24*w25 - w70) - w3 - w71*w97), St*(fA1*z + (1.0/2.0)*fA2*w2 - k*w103*sin(w102)) + Tt*(w101*w104*w106 - w22*w97 + w96*w99*(w105*w106 + (1.0/10.0)*w21*w4 - w94) + 1.0/10.0), St*(fA1*y + (2.0/5.0)*w103*cos(w102)) + Tt*(-1.0/4.0*w100*w104*w107*w98 - w108*w97 + (1.0/2.0)*w16*w4*w5*w7*((2.0/5.0)*sigma0In*w15*w4*z - 1.0/2.0*w105*w107 - w108*((1.0/5.0)*w82 + (1.0/5.0)*w83 + (1.0/5.0)*w85 + (1.0/5.0)*w87 + (1.0/5.0)*w88 + (1.0/5.0)*w89) - w22*(-1.0/5.0*w73 - 1.0/5.0*w74 + (1.0/5.0)*w76 + (1.0/5.0)*w77 + (1.0/5.0)*w79 + (1.0/5.0)*w80) - w71*((1.0/5.0)*w44 + (1.0/5.0)*w46 + (1.0/5.0)*w49 + (1.0/5.0)*w50 + (1.0/5.0)*w53 + (1.0/5.0)*w54) - w90) + (1.0/5.0)*z));
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
        const scalar w7 = (z*z);
        const scalar w8 = -(r2*r2) + w7 + (x*x) + (y*y);
        const scalar w9 = w0 + 1.0/5.0;
        const scalar w10 = sin(w9);
        const scalar w11 = k*z;
        const scalar w12 = (1.0/2.0)*sin(w11) + 1;
        const scalar w13 = k*y;
        const scalar w14 = cos(w13);
        const scalar w15 = sigmaEps*w12*w14;
        const scalar w16 = w10*w15 + 1;
        const scalar w17 = (1.0/(w16));
        const scalar w18 = (1.0/5.0)*w5;
        const scalar w19 = (k*k)*w1;
        const scalar w20 = w18*w19;
        const scalar w21 = w20*x;
        const scalar w22 = sigma0Out*w16;
        const scalar w23 = w5*y;
        const scalar w24 = -w0*w3*w5 + w18*w7 + (1.0/10.0)*w23;
        const scalar w25 = k*sigma0Out*w24;
        const scalar w26 = w15*cos(w9);
        const scalar w27 = sin(gamma1);
        const scalar w28 = sin(gamma2);
        const scalar w29 = (w28*w28);
        const scalar w30 = (w27*w27)*w29;
        const scalar w31 = cos(gamma1);
        const scalar w32 = theta0 + theta1;
        const scalar w33 = cos(w32);
        const scalar w34 = cos(gamma2);
        const scalar w35 = sin(w32);
        const scalar w36 = w27*w35;
        const scalar w37 = w31*w33 - w34*w36;
        const scalar w38 = (w37*w37);
        const scalar w39 = w31*w35;
        const scalar w40 = w27*w33;
        const scalar w41 = -w34*w40 - w39;
        const scalar w42 = (w41*w41);
        const scalar w43 = gel*w38 + gen*w30 + get*w42 + gil*w38 + gin*w30 + git*w42;
        const scalar w44 = w28*w34;
        const scalar w45 = gen*w27;
        const scalar w46 = w44*w45;
        const scalar w47 = gin*w44;
        const scalar w48 = w27*w47;
        const scalar w49 = w28*w35;
        const scalar w50 = w37*w49;
        const scalar w51 = gel*w50;
        const scalar w52 = gil*w50;
        const scalar w53 = w28*w33;
        const scalar w54 = w41*w53;
        const scalar w55 = get*w54;
        const scalar w56 = git*w54;
        const scalar w57 = w46 + w48 + w51 + w52 + w55 + w56;
        const scalar w58 = (1.0/5.0)*z;
        const scalar w59 = w57*w58;
        const scalar w60 = w29*w31;
        const scalar w61 = w45*w60;
        const scalar w62 = gin*w27*w60;
        const scalar w63 = w34*w39 + w40;
        const scalar w64 = w37*w63;
        const scalar w65 = gel*w64;
        const scalar w66 = gil*w64;
        const scalar w67 = w31*w33*w34 - w36;
        const scalar w68 = w41*w67;
        const scalar w69 = get*w68;
        const scalar w70 = git*w68;
        const scalar w71 = -w61 - w62 + w65 + w66 + w69 + w70;
        const scalar w72 = -w4*w43 + w59 - 1.0/10.0*w61 - 1.0/10.0*w62 + (1.0/10.0)*w65 + (1.0/10.0)*w66 + (1.0/10.0)*w69 + (1.0/10.0)*w70;
        const scalar w73 = gen*w31*w44;
        const scalar w74 = w31*w47;
        const scalar w75 = w49*w63;
        const scalar w76 = gel*w75;
        const scalar w77 = gil*w75;
        const scalar w78 = w53*w67;
        const scalar w79 = get*w78;
        const scalar w80 = git*w78;
        const scalar w81 = (w34*w34);
        const scalar w82 = gen*w81;
        const scalar w83 = gin*w81;
        const scalar w84 = w29*(w35*w35);
        const scalar w85 = gel*w84;
        const scalar w86 = w29*(w33*w33);
        const scalar w87 = get*w86;
        const scalar w88 = gil*w84;
        const scalar w89 = git*w86;
        const scalar w90 = w5*(-w4*w57 + w58*(w82 + w83 + w85 + w87 + w88 + w89) - 1.0/10.0*w73 - 1.0/10.0*w74 + (1.0/10.0)*w76 + (1.0/10.0)*w77 + (1.0/10.0)*w79 + (1.0/10.0)*w80);
        const scalar w91 = (1.0/10.0)*w29*(w31*w31);
        const scalar w92 = (1.0/10.0)*(w63*w63);
        const scalar w93 = (1.0/10.0)*(w67*w67);
        const scalar w94 = w5*(gel*w92 + gen*w91 + get*w93 + gil*w92 + gin*w91 + git*w93 - w4*w71 + w58*(-w73 - w74 + w76 + w77 + w79 + w80));
        const scalar w95 = -w22*w24 + w5*w72*x + w90*z + w94*y;
        const scalar w96 = w5*w8;
        const scalar w97 = (1.0/2.0)*w96;
        const scalar w98 = k*w6*w95/(w16*w16);
        const scalar w99 = w97*w98;
        const scalar w100 = w13 + 3.0/10.0;
        const scalar w101 = (2.0/5.0)*z;
        const scalar w102 = fA3*exp(w101);
        const scalar w103 = w22*w5;
        const scalar w104 = sigmaEps*w10;
        const scalar w105 = w104*w25;
        const scalar w106 = w12*sin(w13);
        const scalar w107 = w17*w6;
        const scalar w108 = w107*w97;
        const scalar w109 = w107*w95;
        const scalar w110 = w14*cos(w11);
        const scalar w111 = w5*z;
        return vector(St*(fA0 + fA2*k*w1*((1.0/2.0)*y + 1)) + Tt*((1.0/2.0)*w17*w5*w6*w8*(-w19*w5*w59 - w20*w71*y - w21*w43 - w22*(-w21 - w4*w5) - w25*w26 + w5*w72) + w17*w5*w6*w95*x - w26*w99 - w4), St*(fA1*z + (1.0/2.0)*fA2*w2 - k*w102*sin(w100)) + Tt*(w104*w106*w99 + w108*(-1.0/10.0*w103 + w105*w106 + w94) + w109*w23 + 1.0/10.0), St*(fA1*y + (2.0/5.0)*w102*cos(w100)) + Tt*(-1.0/4.0*w104*w110*w96*w98 + w108*(-w101*w103 - 1.0/2.0*w105*w110 + w111*((1.0/5.0)*w82 + (1.0/5.0)*w83 + (1.0/5.0)*w85 + (1.0/5.0)*w87 + (1.0/5.0)*w88 + (1.0/5.0)*w89) + w23*(-1.0/5.0*w73 - 1.0/5.0*w74 + (1.0/5.0)*w76 + (1.0/5.0)*w77 + (1.0/5.0)*w79 + (1.0/5.0)*w80) + w5*x*((1.0/5.0)*w46 + (1.0/5.0)*w48 + (1.0/5.0)*w51 + (1.0/5.0)*w52 + (1.0/5.0)*w55 + (1.0/5.0)*w56) + w90) + w109*w111 + w58));
    }
    const scalar w0 = k*x;
    const scalar w1 = sin(w0);
    const scalar w2 = (1.0/5.0)*Tt;
    const scalar w3 = k*cos(w0);
    const scalar w4 = fA2*((1.0/2.0)*y + 1);
    const scalar w5 = St*(fA0 + w3*w4);
    const scalar w6 = r1 + r2;
    const scalar w7 = (1.0/(w6));
    const scalar w8 = (r1*r1);
    const scalar w9 = -(r2*r2);
    const scalar w10 = (1.0/((-w8 - w9)));
    const scalar w11 = (x*x) + (y*y) + (z*z);
    const scalar w12 = w10*(w11 - w8);
    const scalar w13 = w12*w7;
    const scalar w14 = (1.0/(r2));
    const scalar w15 = w0 + 1.0/5.0;
    const scalar w16 = sin(w15);
    const scalar w17 = k*z;
    const scalar w18 = (1.0/2.0)*sin(w17) + 1;
    const scalar w19 = k*y;
    const scalar w20 = cos(w19);
    const scalar w21 = sigmaEps*w18*w20;
    const scalar w22 = w16*w21 + 1;
    const scalar w23 = -sigma0Out*w22;
    const scalar w24 = cos(gamma2);
    const scalar w25 = (w24*w24);
    const scalar w26 = sin(gamma2);
    const scalar w27 = (w26*w26);
    const scalar w28 = theta0 + theta1*w12;
    const scalar w29 = sin(w28);
    const scalar w30 = w27*(w29*w29);
    const scalar w31 = cos(w28);
    const scalar w32 = w27*(w31*w31);
    const scalar w33 = gel*w30 + gen*w25 + get*w32 + gil*w30 + gin*w25 + git*w32;
    const scalar w34 = -w23 - w33;
    const scalar w35 = w19 + 3.0/10.0;
    const scalar w36 = cos(w35);
    const scalar w37 = fA3*exp((2.0/5.0)*z);
    const scalar w38 = (2.0/5.0)*w37;
    const scalar w39 = St*(fA1*y + w36*w38);
    const scalar w40 = sin(gamma1);
    const scalar w41 = gen*w40;
    const scalar w42 = w24*w26;
    const scalar w43 = gin*w40;
    const scalar w44 = cos(gamma1);
    const scalar w45 = w29*w40;
    const scalar w46 = w24*w45;
    const scalar w47 = w31*w44 - w46;
    const scalar w48 = w26*w29;
    const scalar w49 = w47*w48;
    const scalar w50 = w26*w31;
    const scalar w51 = w29*w44;
    const scalar w52 = w31*w40;
    const scalar w53 = w24*w52;
    const scalar w54 = -w51 - w53;
    const scalar w55 = git*w50;
    const scalar w56 = gel*w49 + get*w50*w54 + gil*w49 + w41*w42 + w42*w43 + w54*w55;
    const scalar w57 = -w56;
    const scalar w58 = w42*w44;
    const scalar w59 = w24*w51;
    const scalar w60 = w52 + w59;
    const scalar w61 = w48*w60;
    const scalar w62 = w24*w31*w44 - w45;
    const scalar w63 = get*w50;
    const scalar w64 = gel*w61 - gen*w58 + gil*w61 - gin*w58 + w55*w62 + w62*w63;
    const scalar w65 = -w64;
    const scalar w66 = (1.0/2.0)*fA2;
    const scalar w67 = k*sin(w35);
    const scalar w68 = St*(fA1*z + w1*w66 - w37*w67);
    const scalar w69 = w5*w57 + w65*w68;
    const scalar w70 = w14*(w34*w39 + w69);
    const scalar w71 = w27*(w44*w44);
    const scalar w72 = (w60*w60);
    const scalar w73 = (w62*w62);
    const scalar w74 = gel*w72 + gen*w71 + get*w73 + gil*w72 + gin*w71 + git*w73;
    const scalar w75 = -w23 - w74;
    const scalar w76 = w27*w44;
    const scalar w77 = w47*w60;
    const scalar w78 = w54*w62;
    const scalar w79 = gel*w77 + get*w78 + gil*w77 + git*w78 - w41*w76 - w43*w76;
    const scalar w80 = -w79;
    const scalar w81 = w39*w65 + w5*w80;
    const scalar w82 = w14*(w68*w75 + w81);
    const scalar w83 = w27*(w40*w40);
    const scalar w84 = (w47*w47);
    const scalar w85 = (w54*w54);
    const scalar w86 = gel*w84 + gen*w83 + get*w85 + gil*w84 + gin*w83 + git*w85;
    const scalar w87 = -w23 - w86;
    const scalar w88 = w39*w57 + w68*w80;
    const scalar w89 = w14*(w5*w87 + w88);
    const scalar w90 = w70*z + w82*y + w89*x;
    const scalar w91 = w14*w33;
    const scalar w92 = w14*w64;
    const scalar w93 = w14*w56;
    const scalar w94 = w14*(w91*z + w92*y + w93*x);
    const scalar w95 = w14*w74;
    const scalar w96 = w14*w79;
    const scalar w97 = w14*(w92*z + w95*y + w96*x);
    const scalar w98 = w14*w86;
    const scalar w99 = w14*(w93*z + w96*y + w98*x);
    const scalar w100 = w94*z + w97*y + w99*x;
    const scalar w101 = w6/w100;
    const scalar w102 = w101*w90;
    const scalar w103 = (1.0/2.0)*w12*w14;
    const scalar w104 = (1.0/(r1));
    const scalar w105 = -sigma0In*w22;
    const scalar w106 = -w105 - w33;
    const scalar w107 = w104*(w106*w39 + w69);
    const scalar w108 = -w105 - w74;
    const scalar w109 = w104*(w108*w68 + w81);
    const scalar w110 = -w105 - w86;
    const scalar w111 = w104*(w110*w5 + w88);
    const scalar w112 = -w107*z - w109*y - w111*x;
    const scalar w113 = w104*w33;
    const scalar w114 = w104*w64;
    const scalar w115 = w104*w56;
    const scalar w116 = w104*(-w113*z - w114*y - w115*x);
    const scalar w117 = w104*w74;
    const scalar w118 = w104*w79;
    const scalar w119 = w104*(-w114*z - w117*y - w118*x);
    const scalar w120 = w104*w86;
    const scalar w121 = w104*(-w115*z - w118*y - w120*x);
    const scalar w122 = -w116*z - w119*y - w121*x;
    const scalar w123 = w6/w122;
    const scalar w124 = (1.0/2.0)*w104*(1 - w12);
    const scalar w125 = w123*w124;
    const scalar w126 = 2*w102*w103 + 2*w112*w125;
    const scalar w127 = w126*w13;
    const scalar w128 = w10*x;
    const scalar w129 = w11 + w9;
    const scalar w130 = w126*w129*w7;
    const scalar w131 = w14*x;
    const scalar w132 = w10*w102;
    const scalar w133 = w104*x;
    const scalar w134 = w10*w112*w123;
    const scalar w135 = cos(w15);
    const scalar w136 = w135*w21;
    const scalar w137 = k*sigma0Out;
    const scalar w138 = theta1*w128;
    const scalar w139 = 4*w138;
    const scalar w140 = w27*w29*w31;
    const scalar w141 = w139*w140;
    const scalar w142 = gel*w141;
    const scalar w143 = gil*w141;
    const scalar w144 = get*w141;
    const scalar w145 = git*w141;
    const scalar w146 = -w142 - w143 + w144 + w145;
    const scalar w147 = 2*w138;
    const scalar w148 = w147*w50;
    const scalar w149 = w148*w47;
    const scalar w150 = -w147*w51 - w147*w53;
    const scalar w151 = w150*w48;
    const scalar w152 = w31*w44;
    const scalar w153 = -w147*w152 + w147*w46;
    const scalar w154 = gel*w149 + gel*w151 - 2*get*theta1*w10*w26*w29*w54*x + gil*w149 + gil*w151 - 2*git*theta1*w10*w26*w29*w54*x + w153*w55 + w153*w63;
    const scalar w155 = -w154;
    const scalar w156 = w148*w60;
    const scalar w157 = 2*theta1*w10*w24*w31*w44*x - w147*w45;
    const scalar w158 = w157*w48;
    const scalar w159 = -w147*w52 - w147*w59;
    const scalar w160 = gel*w156 + gel*w158 - 2*get*theta1*w10*w26*w29*w62*x + gil*w156 + gil*w158 - 2*git*theta1*w10*w26*w29*w62*x + w159*w55 + w159*w63;
    const scalar w161 = -w160;
    const scalar w162 = w3*w66;
    const scalar w163 = St*w162;
    const scalar w164 = (k*k);
    const scalar w165 = w1*w164*w4;
    const scalar w166 = St*w165;
    const scalar w167 = w155*w5 + w161*w68 + w163*w65 - w166*w57;
    const scalar w168 = w14*z;
    const scalar w169 = -k*sigma0Out*sigmaEps*w135*w18*w20;
    const scalar w170 = w60*(4*theta1*w10*w24*w31*w44*x - w139*w45);
    const scalar w171 = w62*(-w139*w52 - w139*w59);
    const scalar w172 = gel*w170 + get*w171 + gil*w170 + git*w171;
    const scalar w173 = w150*w60;
    const scalar w174 = w157*w47;
    const scalar w175 = w153*w62;
    const scalar w176 = w159*w54;
    const scalar w177 = gel*w173 + gel*w174 + get*w175 + get*w176 + gil*w173 + gil*w174 + git*w175 + git*w176;
    const scalar w178 = -w177;
    const scalar w179 = St*w80;
    const scalar w180 = w161*w39 - w165*w179 + w178*w5;
    const scalar w181 = w14*y;
    const scalar w182 = w47*(-w139*w51 - w139*w53);
    const scalar w183 = w54*(-w139*w152 + w139*w46);
    const scalar w184 = gel*w182 + get*w183 + gil*w182 + git*w183;
    const scalar w185 = w162*w179;
    const scalar w186 = w155*w39 + w178*w68 + w185;
    const scalar w187 = w101*w103;
    const scalar w188 = k*sigma0In;
    const scalar w189 = w104*z;
    const scalar w190 = -k*sigma0In*sigmaEps*w135*w18*w20;
    const scalar w191 = w104*y;
    const scalar w192 = St*w110;
    const scalar w193 = w142 + w143 - w144 - w145;
    const scalar w194 = w103*w6*w90/(w100*w100);
    const scalar w195 = w112*w124*w6/(w122*w122);
    const scalar w196 = w129*w13;
    const scalar w197 = w10*y;
    const scalar w198 = St*(fA1 - w38*w67);
    const scalar w199 = w18*sin(w19);
    const scalar w200 = sigmaEps*w16;
    const scalar w201 = w137*w200;
    const scalar w202 = w199*w201;
    const scalar w203 = theta1*w197;
    const scalar w204 = 4*w203;
    const scalar w205 = w140*w204;
    const scalar w206 = gel*w205 - get*w205 + gil*w205 - git*w205;
    const scalar w207 = 2*w203;
    const scalar w208 = w207*w50;
    const scalar w209 = w208*w47;
    const scalar w210 = -w207*w51 - w207*w53;
    const scalar w211 = w210*w48;
    const scalar w212 = -w152*w207 + w207*w46;
    const scalar w213 = gel*w209 + gel*w211 - 2*get*theta1*w10*w26*w29*w54*y + gil*w209 + gil*w211 - 2*git*theta1*w10*w26*w29*w54*y + w212*w55 + w212*w63;
    const scalar w214 = -w213;
    const scalar w215 = w208*w60;
    const scalar w216 = 2*theta1*w10*w24*w31*w44*y - w207*w45;
    const scalar w217 = w216*w48;
    const scalar w218 = -w207*w52 - w207*w59;
    const scalar w219 = gel*w215 + gel*w217 - 2*get*theta1*w10*w26*w29*w62*y + gil*w215 + gil*w217 - 2*git*theta1*w10*w26*w29*w62*y + w218*w55 + w218*w63;
    const scalar w220 = -w219;
    const scalar w221 = w36*w37;
    const scalar w222 = w164*w221;
    const scalar w223 = St*w222;
    const scalar w224 = w163*w57 + w214*w5 + w220*w68 - w223*w65;
    const scalar w225 = w60*(4*theta1*w10*w24*w31*w44*y - w204*w45);
    const scalar w226 = w62*(-w204*w52 - w204*w59);
    const scalar w227 = gel*w225 + get*w226 + gil*w225 + git*w226;
    const scalar w228 = w198*w65;
    const scalar w229 = w210*w60;
    const scalar w230 = w216*w47;
    const scalar w231 = w212*w62;
    const scalar w232 = w218*w54;
    const scalar w233 = gel*w229 + gel*w230 + get*w231 + get*w232 + gil*w229 + gil*w230 + git*w231 + git*w232;
    const scalar w234 = -w233;
    const scalar w235 = w185 + w220*w39 + w228 + w234*w5;
    const scalar w236 = w47*(-w204*w51 - w204*w53);
    const scalar w237 = w54*(-w152*w204 + w204*w46);
    const scalar w238 = gel*w236 + get*w237 + gil*w236 + git*w237;
    const scalar w239 = -w179*w222 + w198*w57 + w214*w39 + w234*w68;
    const scalar w240 = w188*w200;
    const scalar w241 = w199*w240;
    const scalar w242 = w10*z;
    const scalar w243 = cos(w17);
    const scalar w244 = (1.0/2.0)*w20*w243;
    const scalar w245 = theta1*w242;
    const scalar w246 = 4*w245;
    const scalar w247 = w140*w246;
    const scalar w248 = gel*w247;
    const scalar w249 = gil*w247;
    const scalar w250 = get*w247;
    const scalar w251 = git*w247;
    const scalar w252 = -w248 - w249 + w250 + w251;
    const scalar w253 = (4.0/25.0)*St*w221;
    const scalar w254 = 2*w245;
    const scalar w255 = w254*w50;
    const scalar w256 = w255*w47;
    const scalar w257 = -w254*w51 - w254*w53;
    const scalar w258 = w257*w48;
    const scalar w259 = -w152*w254 + w254*w46;
    const scalar w260 = gel*w256 + gel*w258 - 2*get*theta1*w10*w26*w29*w54*z + gil*w256 + gil*w258 - 2*git*theta1*w10*w26*w29*w54*z + w259*w55 + w259*w63;
    const scalar w261 = -w260;
    const scalar w262 = w255*w60;
    const scalar w263 = 2*theta1*w10*w24*w31*w44*z - w254*w45;
    const scalar w264 = w263*w48;
    const scalar w265 = -w254*w52 - w254*w59;
    const scalar w266 = gel*w262 + gel*w264 - 2*get*theta1*w10*w26*w29*w62*z + gil*w262 + gil*w264 - 2*git*theta1*w10*w26*w29*w62*z + w265*w55 + w265*w63;
    const scalar w267 = -w266;
    const scalar w268 = w228 + w261*w5 + w267*w68;
    const scalar w269 = -1.0/2.0*k*sigma0Out*sigmaEps*w16*w20*w243;
    const scalar w270 = w60*(4*theta1*w10*w24*w31*w44*z - w246*w45);
    const scalar w271 = w62*(-w246*w52 - w246*w59);
    const scalar w272 = gel*w270 + get*w271 + gil*w270 + git*w271;
    const scalar w273 = w257*w60;
    const scalar w274 = w263*w47;
    const scalar w275 = w259*w62;
    const scalar w276 = w265*w54;
    const scalar w277 = gel*w273 + gel*w274 + get*w275 + get*w276 + gil*w273 + gil*w274 + git*w275 + git*w276;
    const scalar w278 = -w277;
    const scalar w279 = w253*w65 + w267*w39 + w278*w5;
    const scalar w280 = w47*(-w246*w51 - w246*w53);
    const scalar w281 = w54*(-w152*w246 + w246*w46);
    const scalar w282 = gel*w280 + get*w281 + gil*w280 + git*w281;
    const scalar w283 = w198*w80 + w253*w57 + w261*w39 + w278*w68;
    const scalar w284 = -1.0/2.0*k*sigma0In*sigmaEps*w16*w20*w243;
    const scalar w285 = w248 + w249 - w250 - w251;
    return vector(-k*w1*w2 + w127*x + w128*w130 + w196*(w125*(-w111 - w133*(-w165*w192 + w186 + w5*(-w184 - w190)) - w189*(w167 + w39*(w136*w188 + w146)) - w191*(w108*w163 + w180 + w68*(-w172 - w190))) + w131*w132 - w133*w134 + w187*(w131*(-w166*w87 + w186 + w5*(-w169 - w184)) + w168*(w167 + w39*(w136*w137 + w146)) + w181*(w163*w75 + w180 + w68*(-w169 - w172)) + w89) + w194*(-w131*(w131*w184 + w154*w168 + w177*w181 + w98) - w168*(w131*w154 + w160*w181 + w168*w193 + w93) - w181*(w131*w177 + w160*w168 + w172*w181 + w96) - w99) + w195*(w121 + w133*(-w120 - w133*w184 - w154*w189 - w177*w191) + w189*(-w115 - w133*w154 - w160*w191 - w189*w193) + w191*(-w118 - w133*w177 - w160*w189 - w172*w191))) + w5, (1.0/10.0)*Tt + w127*y + w130*w197 + w196*(w125*(-w109 - w133*(w162*w192 + w239 + w5*(-w238 - w241)) - w189*(w106*w198 + w224 + w39*(-w206 - w241)) - w191*(-w108*w223 + w235 + w68*(-w227 - w241))) + w132*w181 - w134*w191 + w187*(w131*(w163*w87 + w239 + w5*(-w202 - w238)) + w168*(w198*w34 + w224 + w39*(-w202 - w206)) + w181*(-w223*w75 + w235 + w68*(-w202 - w227)) + w82) + w194*(-w131*(w131*w238 + w168*w213 + w181*w233 + w96) - w168*(w131*w213 + w168*w206 + w181*w219 + w92) - w181*(w131*w233 + w168*w219 + w181*w227 + w95) - w97) + w195*(w119 + w133*(-w118 - w133*w238 - w189*w213 - w191*w233) + w189*(-w114 - w133*w213 - w189*w206 - w191*w219) + w191*(-w117 - w133*w233 - w189*w219 - w191*w227))) + w68, w127*z + w130*w242 + w196*(w125*(-w107 - w133*(w283 + w5*(-w282 - w284)) - w189*(w106*w253 + w268 + w39*(w240*w244 + w252)) - w191*(w108*w198 + w279 + w68*(-w272 - w284))) + w132*w168 - w134*w189 + w187*(w131*(w283 + w5*(-w269 - w282)) + w168*(w253*w34 + w268 + w39*(w201*w244 + w252)) + w181*(w198*w75 + w279 + w68*(-w269 - w272)) + w70) + w194*(-w131*(w131*w282 + w168*w260 + w181*w277 + w93) - w168*(w131*w260 + w168*w285 + w181*w266 + w91) - w181*(w131*w277 + w168*w266 + w181*w272 + w92) - w94) + w195*(w116 + w133*(-w115 - w133*w282 - w189*w260 - w191*w277) + w189*(-w113 - w133*w260 - w189*w285 - w191*w266) + w191*(-w114 - w133*w277 - w189*w266 - w191*w272))) + w2*z + w39);
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
    const scalar theta0 = p.theta0;
    const scalar theta1 = p.theta1;
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

    const scalar w0 = (1.0/2.0)*y + 1;
    const scalar w1 = k*x;
    const scalar w2 = fA2*k*cos(w1);
    const scalar w3 = St*(fA0 + w0*w2);
    const scalar w4 = k*z;
    const scalar w5 = sin(w4);
    const scalar w6 = k*y;
    const scalar w7 = hA1*cos(w6);
    const scalar w8 = (3.0/5.0)*w7;
    const scalar w9 = hA2*w5 - w8*x;
    const scalar w10 = r1 + r2;
    const scalar w11 = (1.0/(w10));
    const scalar w12 = (r1*r1);
    const scalar w13 = -(r2*r2);
    const scalar w14 = (1.0/((-w12 - w13)));
    const scalar w15 = (x*x);
    const scalar w16 = w15 + (y*y) + (z*z);
    const scalar w17 = w14*(-w12 + w16);
    const scalar w18 = w11*w17;
    const scalar w19 = (1.0/(r2));
    const scalar w20 = hA2*cos(w4);
    const scalar w21 = w6 + 3.0/10.0;
    const scalar w22 = fA3*exp((2.0/5.0)*z);
    const scalar w23 = w22*cos(w21);
    const scalar w24 = St*(hA0 + w1*w20) - St*(fA1*y + (2.0/5.0)*w23);
    const scalar w25 = cos(gamma2);
    const scalar w26 = theta0 + theta1*w17;
    const scalar w27 = sin(w26);
    const scalar w28 = sin(gamma2);
    const scalar w29 = (w28*w28);
    const scalar w30 = gil*w29;
    const scalar w31 = cos(w26);
    const scalar w32 = git*w29;
    const scalar w33 = gin*(w25*w25) + (w27*w27)*w30 + (w31*w31)*w32;
    const scalar w34 = St*w9 - w3;
    const scalar w35 = sin(gamma1);
    const scalar w36 = gin*w25*w28;
    const scalar w37 = cos(gamma1);
    const scalar w38 = w27*w35;
    const scalar w39 = w25*w38;
    const scalar w40 = w31*w37 - w39;
    const scalar w41 = gil*w40;
    const scalar w42 = w27*w28;
    const scalar w43 = w27*w37;
    const scalar w44 = w31*w35;
    const scalar w45 = w25*w44;
    const scalar w46 = -w43 - w45;
    const scalar w47 = git*w46;
    const scalar w48 = w28*w31;
    const scalar w49 = w35*w36 + w41*w42 + w47*w48;
    const scalar w50 = 1 - 3.0/10.0*w15;
    const scalar w51 = St*k;
    const scalar w52 = sin(w6);
    const scalar w53 = hA1*w50*w51*w52;
    const scalar w54 = fA2*sin(w1);
    const scalar w55 = k*w22*sin(w21);
    const scalar w56 = St*(fA1*z + (1.0/2.0)*w54 - w55);
    const scalar w57 = -w53 - w56;
    const scalar w58 = w25*w43;
    const scalar w59 = w44 + w58;
    const scalar w60 = gil*w59;
    const scalar w61 = w25*w31*w37 - w38;
    const scalar w62 = git*w61;
    const scalar w63 = -w36*w37 + w42*w60 + w48*w62;
    const scalar w64 = w24*w33 + w34*w49 + w57*w63;
    const scalar w65 = w19*w64;
    const scalar w66 = gin*w29;
    const scalar w67 = gil*(w59*w59) + git*(w61*w61) + (w37*w37)*w66;
    const scalar w68 = -w35*w37*w66 + w41*w59 + w47*w61;
    const scalar w69 = w24*w63 + w34*w68 + w57*w67;
    const scalar w70 = w19*w69;
    const scalar w71 = gil*(w40*w40) + git*(w46*w46) + (w35*w35)*w66;
    const scalar w72 = w24*w49 + w34*w71 + w57*w68;
    const scalar w73 = w19*w72;
    const scalar w74 = -w65*z - w70*y - w73*x;
    const scalar w75 = w19*w33;
    const scalar w76 = w19*w63;
    const scalar w77 = w19*w49;
    const scalar w78 = w19*(w75*z + w76*y + w77*x);
    const scalar w79 = w19*w67;
    const scalar w80 = w19*w68;
    const scalar w81 = w19*(w76*z + w79*y + w80*x);
    const scalar w82 = w19*w71;
    const scalar w83 = w19*(w77*z + w80*y + w82*x);
    const scalar w84 = w78*z + w81*y + w83*x;
    const scalar w85 = w10/w84;
    const scalar w86 = w74*w85;
    const scalar w87 = (1.0/2.0)*w17*w19;
    const scalar w88 = (1.0/(r1));
    const scalar w89 = w64*w88;
    const scalar w90 = w69*w88;
    const scalar w91 = w72*w88;
    const scalar w92 = w89*z + w90*y + w91*x;
    const scalar w93 = w33*w88;
    const scalar w94 = w63*w88;
    const scalar w95 = w49*w88;
    const scalar w96 = w88*(-w93*z - w94*y - w95*x);
    const scalar w97 = w67*w88;
    const scalar w98 = w68*w88;
    const scalar w99 = w88*(-w94*z - w97*y - w98*x);
    const scalar w100 = w71*w88;
    const scalar w101 = w88*(-w100*x - w95*z - w98*y);
    const scalar w102 = -w101*x - w96*z - w99*y;
    const scalar w103 = w10/w102;
    const scalar w104 = (1.0/2.0)*w88*(1 - w17);
    const scalar w105 = w103*w104;
    const scalar w106 = w105*w92 + w86*w87;
    const scalar w107 = 2*x;
    const scalar w108 = w106*w107;
    const scalar w109 = w13 + w16;
    const scalar w110 = w109*w11*w14;
    const scalar w111 = w19*x;
    const scalar w112 = w14*w86;
    const scalar w113 = w88*x;
    const scalar w114 = w103*w14*w92;
    const scalar w115 = w20*w51;
    const scalar w116 = theta1*w14;
    const scalar w117 = 4*x;
    const scalar w118 = w116*w117;
    const scalar w119 = w27*w31;
    const scalar w120 = w118*w119;
    const scalar w121 = w120*w30 - w120*w32;
    const scalar w122 = (3.0/5.0)*St*hA1*k*w52*x - 1.0/2.0*St*w2;
    const scalar w123 = St*(k*k);
    const scalar w124 = -St*w8 + w0*w123*w54;
    const scalar w125 = w107*w116;
    const scalar w126 = w125*w48;
    const scalar w127 = w125*w42;
    const scalar w128 = w116*w43;
    const scalar w129 = -w107*w128 - w125*w45;
    const scalar w130 = gil*w42;
    const scalar w131 = w116*w31*w37;
    const scalar w132 = -w107*w131 + w125*w39;
    const scalar w133 = git*w48;
    const scalar w134 = w126*w41 - w127*w47 + w129*w130 + w132*w133;
    const scalar w135 = 2*theta1*w14*w25*w31*w37*x - w125*w38;
    const scalar w136 = -w125*w44 - w125*w58;
    const scalar w137 = w126*w60 - w127*w62 + w130*w135 + w133*w136;
    const scalar w138 = w115*w33 + w121*w24 + w122*w63 + w124*w49 + w134*w34 + w137*w57;
    const scalar w139 = w19*z;
    const scalar w140 = w60*(4*theta1*w14*w25*w31*w37*x - w118*w38) + w62*(-w118*w44 - w118*w58);
    const scalar w141 = w129*w60 + w132*w62 + w135*w41 + w136*w47;
    const scalar w142 = w115*w63 + w122*w67 + w124*w68 + w137*w24 + w140*w57 + w141*w34;
    const scalar w143 = w19*y;
    const scalar w144 = w115*w49;
    const scalar w145 = w122*w68;
    const scalar w146 = w41*(-w117*w128 - w118*w45) + w47*(-w117*w131 + w118*w39);
    const scalar w147 = w124*w71 + w134*w24 + w141*w57 + w144 + w145 + w146*w34;
    const scalar w148 = w85*w87;
    const scalar w149 = w88*z;
    const scalar w150 = w88*y;
    const scalar w151 = w10*w74*w87/(w84*w84);
    const scalar w152 = w10*w104*w92/(w102*w102);
    const scalar w153 = w109*w18;
    const scalar w154 = 2*y;
    const scalar w155 = w106*w154;
    const scalar w156 = St*(fA1 - 2.0/5.0*w55);
    const scalar w157 = 4*y;
    const scalar w158 = w116*w157;
    const scalar w159 = w119*w158;
    const scalar w160 = w159*w30 - w159*w32;
    const scalar w161 = w123*w23 - w123*w50*w7;
    const scalar w162 = w116*w154;
    const scalar w163 = w162*w48;
    const scalar w164 = w162*w42;
    const scalar w165 = -w128*w154 - w162*w45;
    const scalar w166 = -w131*w154 + w162*w39;
    const scalar w167 = w130*w165 + w133*w166 + w163*w41 - w164*w47;
    const scalar w168 = 2*theta1*w14*w25*w31*w37*y - w162*w38;
    const scalar w169 = -w162*w44 - w162*w58;
    const scalar w170 = w130*w168 + w133*w169 + w163*w60 - w164*w62;
    const scalar w171 = w122*w49 - w156*w33 + w160*w24 + w161*w63 + w167*w34 + w170*w57;
    const scalar w172 = -w156*w63;
    const scalar w173 = w60*(4*theta1*w14*w25*w31*w37*y - w158*w38) + w62*(-w158*w44 - w158*w58);
    const scalar w174 = w165*w60 + w166*w62 + w168*w41 + w169*w47;
    const scalar w175 = w145 + w161*w67 + w170*w24 + w172 + w173*w57 + w174*w34;
    const scalar w176 = w41*(-w128*w157 - w158*w45) + w47*(-w131*w157 + w158*w39);
    const scalar w177 = w122*w71 - w156*w49 + w161*w68 + w167*w24 + w174*w57 + w176*w34;
    const scalar w178 = -4.0/25.0*St*w23 - hA2*w123*w5*x;
    const scalar w179 = 4*z;
    const scalar w180 = w116*w179;
    const scalar w181 = w119*w180;
    const scalar w182 = w181*w30 - w181*w32;
    const scalar w183 = 2*z;
    const scalar w184 = w116*w183;
    const scalar w185 = w184*w48;
    const scalar w186 = w184*w42;
    const scalar w187 = -w128*w183 - w184*w45;
    const scalar w188 = -w131*w183 + w184*w39;
    const scalar w189 = w130*w187 + w133*w188 + w185*w41 - w186*w47;
    const scalar w190 = 2*theta1*w14*w25*w31*w37*z - w184*w38;
    const scalar w191 = -w184*w44 - w184*w58;
    const scalar w192 = w130*w190 + w133*w191 + w185*w60 - w186*w62;
    const scalar w193 = w144 + w172 + w178*w33 + w182*w24 + w189*w34 + w192*w57;
    const scalar w194 = w60*(4*theta1*w14*w25*w31*w37*z - w180*w38) + w62*(-w180*w44 - w180*w58);
    const scalar w195 = w187*w60 + w188*w62 + w190*w41 + w191*w47;
    const scalar w196 = w115*w68 - w156*w67 + w178*w63 + w192*w24 + w194*w57 + w195*w34;
    const scalar w197 = w41*(-w128*w179 - w180*w45) + w47*(-w131*w179 + w180*w39);
    const scalar w198 = w115*w71 - w156*w68 + w178*w49 + w189*w24 + w195*w57 + w197*w34;
    const scalar w199 = w106*w183;
    return vector(St*w9 + w108*w110 + w108*w18 + w153*(w105*(w113*w147 + w138*w149 + w142*w150 + w91) + w111*w112 - w113*w114 + w148*(-w111*w147 - w138*w139 - w142*w143 - w73) + w151*(-w111*(w111*w146 + w134*w139 + w141*w143 + w82) - w139*(w111*w134 + w121*w139 + w137*w143 + w77) - w143*(w111*w141 + w137*w139 + w140*w143 + w80) - w83) + w152*(w101 + w113*(-w100 - w113*w146 - w134*w149 - w141*w150) + w149*(-w113*w134 - w121*w149 - w137*w150 - w95) + w150*(-w113*w141 - w137*w149 - w140*w150 - w98))) - w3, w110*w155 + w153*(w105*(w113*w177 + w149*w171 + w150*w175 + w90) + w112*w143 - w114*w150 + w148*(-w111*w177 - w139*w171 - w143*w175 - w70) + w151*(-w111*(w111*w176 + w139*w167 + w143*w174 + w80) - w139*(w111*w167 + w139*w160 + w143*w170 + w76) - w143*(w111*w174 + w139*w170 + w143*w173 + w79) - w81) + w152*(w113*(-w113*w176 - w149*w167 - w150*w174 - w98) + w149*(-w113*w167 - w149*w160 - w150*w170 - w94) + w150*(-w113*w174 - w149*w170 - w150*w173 - w97) + w99)) + w155*w18 - w53 - w56, w110*w199 + w153*(w105*(w113*w198 + w149*w193 + w150*w196 + w89) + w112*w139 - w114*w149 + w148*(-w111*w198 - w139*w193 - w143*w196 - w65) + w151*(-w111*(w111*w197 + w139*w189 + w143*w195 + w77) - w139*(w111*w189 + w139*w182 + w143*w192 + w75) - w143*(w111*w195 + w139*w192 + w143*w194 + w76) - w78) + w152*(w113*(-w113*w197 - w149*w189 - w150*w195 - w95) + w149*(-w113*w189 - w149*w182 - w150*w192 - w93) + w150*(-w113*w195 - w149*w192 - w150*w194 - w94) + w96)) + w18*w199 + w24);
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
    const scalar theta0 = p.theta0;
    const scalar theta1 = p.theta1;
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
    const scalar w3 = k*y;
    const scalar w4 = sin(w0);
    const scalar w5 = (x*x);
    const scalar w6 = 1 - 3.0/10.0*w5;
    const scalar w7 = hA1*cos(w3);
    const scalar w8 = fA1*y;
    const scalar w9 = (1.0/2.0)*y + 1;
    const scalar w10 = fA2*sin(w2);
    const scalar w11 = w3 + 3.0/10.0;
    const scalar w12 = fA3*exp((2.0/5.0)*z);
    const scalar w13 = w12*cos(w11);
    const scalar w14 = r1 + r2;
    const scalar w15 = (1.0/(r2));
    const scalar w16 = cos(gamma2);
    const scalar w17 = sin(gamma2);
    const scalar w18 = (w17*w17);
    const scalar w19 = (r1*r1);
    const scalar w20 = -(r2*r2);
    const scalar w21 = w5 + (y*y) + (z*z);
    const scalar w22 = (-w19 + w21)/(-w19 - w20);
    const scalar w23 = theta0 + theta1*w22;
    const scalar w24 = sin(w23);
    const scalar w25 = cos(w23);
    const scalar w26 = gil*w18*(w24*w24) + gin*(w16*w16) + git*w18*(w25*w25);
    const scalar w27 = w15*z;
    const scalar w28 = cos(gamma1);
    const scalar w29 = gin*w16*w17;
    const scalar w30 = sin(gamma1);
    const scalar w31 = w25*w30;
    const scalar w32 = w24*w28;
    const scalar w33 = w16*w32 + w31;
    const scalar w34 = w17*w24;
    const scalar w35 = w24*w30;
    const scalar w36 = w16*w25*w28 - w35;
    const scalar w37 = w17*w25;
    const scalar w38 = gil*w33*w34 + git*w36*w37 - w28*w29;
    const scalar w39 = w15*y;
    const scalar w40 = -w16*w35 + w25*w28;
    const scalar w41 = gil*w40;
    const scalar w42 = -w16*w31 - w32;
    const scalar w43 = git*w42;
    const scalar w44 = w29*w30 + w34*w41 + w37*w43;
    const scalar w45 = w15*x;
    const scalar w46 = gin*w18;
    const scalar w47 = gil*(w33*w33) + git*(w36*w36) + (w28*w28)*w46;
    const scalar w48 = -w28*w30*w46 + w33*w41 + w36*w43;
    const scalar w49 = gil*(w40*w40) + git*(w42*w42) + (w30*w30)*w46;
    const scalar w50 = St*(hA0 + hA2*w1*w2) - St*((2.0/5.0)*w13 + w8);
    const scalar w51 = -St*(fA0 + fA2*k*w9*cos(w2)) + St*(hA2*w4 - 3.0/5.0*w7*x);
    const scalar w52 = -St*hA1*k*w6*sin(w3) - St*(fA1*z - k*w12*sin(w11) + (1.0/2.0)*w10);
    const scalar w53 = w26*w50 + w38*w52 + w44*w51;
    const scalar w54 = w38*w50 + w47*w52 + w48*w51;
    const scalar w55 = w44*w50 + w48*w52 + w49*w51;
    const scalar w56 = (1.0/2.0)*w14;
    const scalar w57 = (1.0/(r1));
    const scalar w58 = w57*z;
    const scalar w59 = w57*y;
    const scalar w60 = w57*x;
    return St*(hA0*z + hA2*w4*x + w6*w7) - St*(fA0*x + w10*w9 + w13 + w8*z) + (t + 1)*((1.0/4.0)*w1*sin(w2 + 1.0/5.0)*cos(w3 - 1.0/3.0) + 1) + w22*(w20 + w21)*(w15*w22*w56*(-w27*w53 - w39*w54 - w45*w55)/(w27*(w26*w27 + w38*w39 + w44*w45) + w39*(w27*w38 + w39*w47 + w45*w48) + w45*(w27*w44 + w39*w48 + w45*w49)) + w56*w57*(1 - w22)*(w53*w58 + w54*w59 + w55*w60)/(-w58*(-w26*w58 - w38*w59 - w44*w60) - w59*(-w38*w58 - w47*w59 - w48*w60) - w60*(-w44*w58 - w48*w59 - w49*w60)))/w14;
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
    const scalar w19 = w18*z;
    const scalar w20 = git*w5;
    const scalar w21 = w16*w20;
    const scalar w22 = w21*z;
    const scalar w23 = w19 - w22;
    const scalar w24 = (1.0/5.0)*Tt;
    const scalar w25 = k*x;
    const scalar w26 = k*z;
    const scalar w27 = cos(w26);
    const scalar w28 = hA2*w27;
    const scalar w29 = St*(hA0 + w25*w28);
    const scalar w30 = r1 + r2;
    const scalar w31 = (1.0/(w30));
    const scalar w32 = w10*w31;
    const scalar w33 = (1.0/(r2));
    const scalar w34 = cos(gamma2);
    const scalar w35 = (w34*w34);
    const scalar w36 = (w13*w13)*w15;
    const scalar w37 = gil*w36;
    const scalar w38 = (w12*w12)*w15;
    const scalar w39 = git*w38;
    const scalar w40 = gin*w35 + w37 + w39;
    const scalar w41 = w33*w40;
    const scalar w42 = cos(gamma1);
    const scalar w43 = w14*w34;
    const scalar w44 = gin*w43;
    const scalar w45 = sin(gamma1);
    const scalar w46 = w12*w45;
    const scalar w47 = w13*w42;
    const scalar w48 = w34*w47;
    const scalar w49 = w46 + w48;
    const scalar w50 = gil*w49;
    const scalar w51 = w13*w14;
    const scalar w52 = w50*w51;
    const scalar w53 = w13*w45;
    const scalar w54 = w12*w34*w42 - w53;
    const scalar w55 = git*w54;
    const scalar w56 = w12*w14;
    const scalar w57 = w55*w56;
    const scalar w58 = -w42*w44 + w52 + w57;
    const scalar w59 = w33*w58;
    const scalar w60 = w34*w53;
    const scalar w61 = w12*w42 - w60;
    const scalar w62 = gil*w61;
    const scalar w63 = w51*w62;
    const scalar w64 = w34*w46;
    const scalar w65 = -w47 - w64;
    const scalar w66 = git*w65;
    const scalar w67 = w56*w66;
    const scalar w68 = w44*w45 + w63 + w67;
    const scalar w69 = w33*w68;
    const scalar w70 = w33*(w41*z + w59*y + w69*x);
    const scalar w71 = (w42*w42);
    const scalar w72 = gin*w15;
    const scalar w73 = (w49*w49);
    const scalar w74 = (w54*w54);
    const scalar w75 = gil*w73 + git*w74 + w71*w72;
    const scalar w76 = w33*w75;
    const scalar w77 = w42*w45;
    const scalar w78 = w49*w62 + w54*w66 - w72*w77;
    const scalar w79 = w33*w78;
    const scalar w80 = w33*(w59*z + w76*y + w79*x);
    const scalar w81 = (w45*w45);
    const scalar w82 = (w61*w61);
    const scalar w83 = (w65*w65);
    const scalar w84 = gil*w82 + git*w83 + w72*w81;
    const scalar w85 = w33*w84;
    const scalar w86 = w33*(w69*z + w79*y + w85*x);
    const scalar w87 = w70*z + w80*y + w86*x;
    const scalar w88 = (1.0/(w87));
    const scalar w89 = k*y;
    const scalar w90 = w89 + 3.0/10.0;
    const scalar w91 = fA3*exp((2.0/5.0)*z);
    const scalar w92 = w91*cos(w90);
    const scalar w93 = (2.0/5.0)*w92;
    const scalar w94 = St*(fA1*y + w93);
    const scalar w95 = w29 - w94;
    const scalar w96 = (1.0/2.0)*y + 1;
    const scalar w97 = cos(w25);
    const scalar w98 = fA2*k*w97;
    const scalar w99 = St*(fA0 + w96*w98);
    const scalar w100 = sin(w26);
    const scalar w101 = cos(w89);
    const scalar w102 = hA1*w101;
    const scalar w103 = (3.0/5.0)*w102;
    const scalar w104 = w103*x;
    const scalar w105 = hA2*w100 - w104;
    const scalar w106 = St*w105 - w99;
    const scalar w107 = 1 - 3.0/10.0*w6;
    const scalar w108 = St*k;
    const scalar w109 = sin(w89);
    const scalar w110 = hA1*w109;
    const scalar w111 = w108*w110;
    const scalar w112 = w107*w111;
    const scalar w113 = sin(w25);
    const scalar w114 = fA2*w113;
    const scalar w115 = (1.0/2.0)*w114;
    const scalar w116 = w91*sin(w90);
    const scalar w117 = k*w116;
    const scalar w118 = St*(fA1*z + w115 - w117);
    const scalar w119 = -w112 - w118;
    const scalar w120 = w106*w68 + w119*w58 + w40*w95;
    const scalar w121 = w120*w33;
    const scalar w122 = w106*w78 + w119*w75 + w58*w95;
    const scalar w123 = w122*w33;
    const scalar w124 = w106*w84 + w119*w78 + w68*w95;
    const scalar w125 = w124*w33;
    const scalar w126 = -w121*z - w123*y - w125*x;
    const scalar w127 = w126*w88;
    const scalar w128 = (1.0/2.0)*w30;
    const scalar w129 = w10*w33;
    const scalar w130 = w128*w129;
    const scalar w131 = (1.0/(r1));
    const scalar w132 = w131*w40;
    const scalar w133 = w131*w58;
    const scalar w134 = w131*w68;
    const scalar w135 = w131*(-w132*z - w133*y - w134*x);
    const scalar w136 = w131*w75;
    const scalar w137 = w131*w78;
    const scalar w138 = w131*(-w133*z - w136*y - w137*x);
    const scalar w139 = w131*w84;
    const scalar w140 = w131*(-w134*z - w137*y - w139*x);
    const scalar w141 = -w135*z - w138*y - w140*x;
    const scalar w142 = (1.0/(w141));
    const scalar w143 = w120*w131;
    const scalar w144 = w122*w131;
    const scalar w145 = w124*w131;
    const scalar w146 = w143*z + w144*y + w145*x;
    const scalar w147 = w142*w146;
    const scalar w148 = w131*(1 - w10);
    const scalar w149 = w128*w148;
    const scalar w150 = w127*w130 + w147*w149;
    const scalar w151 = 2*w150;
    const scalar w152 = w151*w32;
    const scalar w153 = w1 + w9;
    const scalar w154 = w3*w31;
    const scalar w155 = w153*w154;
    const scalar w156 = w151*w155;
    const scalar w157 = gel*w36;
    const scalar w158 = get*w38;
    const scalar w159 = gen*w35 + w157 + w158 + w40;
    const scalar w160 = w159*w33;
    const scalar w161 = gel*w49;
    const scalar w162 = w161*w51;
    const scalar w163 = get*w54;
    const scalar w164 = w163*w56;
    const scalar w165 = gen*w43;
    const scalar w166 = w162 + w164 - w165*w42 + w58;
    const scalar w167 = w166*w33;
    const scalar w168 = gel*w61;
    const scalar w169 = w168*w51;
    const scalar w170 = get*w65;
    const scalar w171 = w170*w56;
    const scalar w172 = w165*w45 + w169 + w171 + w68;
    const scalar w173 = w172*w33;
    const scalar w174 = w33*(w160*z + w167*y + w173*x);
    const scalar w175 = gen*w15;
    const scalar w176 = gel*w73 + get*w74 + w175*w71 + w75;
    const scalar w177 = w176*w33;
    const scalar w178 = w168*w49 + w170*w54 - w175*w77 + w78;
    const scalar w179 = w178*w33;
    const scalar w180 = w33*(w167*z + w177*y + w179*x);
    const scalar w181 = gel*w82 + get*w83 + w175*w81 + w84;
    const scalar w182 = w181*w33;
    const scalar w183 = w33*(w173*z + w179*y + w182*x);
    const scalar w184 = w174*z + w180*y + w183*x;
    const scalar w185 = (1.0/(w184));
    const scalar w186 = w25 + 1.0/5.0;
    const scalar w187 = sin(w186);
    const scalar w188 = (1.0/2.0)*w100;
    const scalar w189 = w188 + 1;
    const scalar w190 = sigmaEps*w189;
    const scalar w191 = w101*w190;
    const scalar w192 = w187*w191;
    const scalar w193 = w192 + 1;
    const scalar w194 = -sigma0Out*w193;
    const scalar w195 = -w159 - w194;
    const scalar w196 = -w172;
    const scalar w197 = -w166;
    const scalar w198 = w118*w197 + w196*w99;
    const scalar w199 = w33*(w195*w94 + w198);
    const scalar w200 = -w176 - w194;
    const scalar w201 = -w178;
    const scalar w202 = w197*w94 + w201*w99;
    const scalar w203 = w33*(w118*w200 + w202);
    const scalar w204 = -w181 - w194;
    const scalar w205 = w118*w201 + w196*w94;
    const scalar w206 = w33*(w204*w99 + w205);
    const scalar w207 = w199*z + w203*y + w206*x;
    const scalar w208 = w185*w207;
    const scalar w209 = w131*w159;
    const scalar w210 = w131*w166;
    const scalar w211 = w131*w172;
    const scalar w212 = w131*(-w209*z - w210*y - w211*x);
    const scalar w213 = w131*w176;
    const scalar w214 = w131*w178;
    const scalar w215 = w131*(-w210*z - w213*y - w214*x);
    const scalar w216 = w131*w181;
    const scalar w217 = w131*(-w211*z - w214*y - w216*x);
    const scalar w218 = -w212*z - w215*y - w217*x;
    const scalar w219 = (1.0/(w218));
    const scalar w220 = -sigma0In*w193;
    const scalar w221 = -w159 - w220;
    const scalar w222 = w131*(w198 + w221*w94);
    const scalar w223 = -w176 - w220;
    const scalar w224 = w131*(w118*w223 + w202);
    const scalar w225 = -w181 - w220;
    const scalar w226 = w131*(w205 + w225*w99);
    const scalar w227 = -w222*z - w224*y - w226*x;
    const scalar w228 = w219*w227;
    const scalar w229 = w130*w208 + w149*w228;
    const scalar w230 = 2*w229;
    const scalar w231 = w230*w32;
    const scalar w232 = w155*w230;
    const scalar w233 = w33*z;
    const scalar w234 = w3*w30;
    const scalar w235 = w233*w234;
    const scalar w236 = w131*z;
    const scalar w237 = w234*w236;
    const scalar w238 = (k*k);
    const scalar w239 = St*w238;
    const scalar w240 = hA2*w100*w239;
    const scalar w241 = w240*x;
    const scalar w242 = St*w92;
    const scalar w243 = (4.0/25.0)*w242;
    const scalar w244 = -w241 - w243;
    const scalar w245 = w108*w28;
    const scalar w246 = w245*w68;
    const scalar w247 = St*(fA1 - 2.0/5.0*w117);
    const scalar w248 = -w247*w58;
    const scalar w249 = 2*w4;
    const scalar w250 = w249*w56;
    const scalar w251 = w250*w62;
    const scalar w252 = w249*w51;
    const scalar w253 = w252*w66;
    const scalar w254 = w249*w47;
    const scalar w255 = w249*w64;
    const scalar w256 = -w254*z - w255*z;
    const scalar w257 = gil*w256;
    const scalar w258 = w12*w42;
    const scalar w259 = w249*w258;
    const scalar w260 = w249*w60;
    const scalar w261 = -w259*z + w260*z;
    const scalar w262 = git*w261;
    const scalar w263 = w251*z - w253*z + w257*w51 + w262*w56;
    const scalar w264 = w250*w50;
    const scalar w265 = w252*w55;
    const scalar w266 = w249*w53;
    const scalar w267 = 2*theta1*w12*w3*w34*w42*z - w266*z;
    const scalar w268 = gil*w267;
    const scalar w269 = w249*w46;
    const scalar w270 = w249*w48;
    const scalar w271 = -w269*z - w270*z;
    const scalar w272 = git*w271;
    const scalar w273 = w264*z - w265*z + w268*w51 + w272*w56;
    const scalar w274 = w106*w263 + w119*w273 + w23*w95 + w244*w40 + w246 + w248;
    const scalar w275 = w5*w53;
    const scalar w276 = 4*theta1*w12*w3*w34*w42*z - w275*z;
    const scalar w277 = w46*w5;
    const scalar w278 = w48*w5;
    const scalar w279 = -w277*z - w278*z;
    const scalar w280 = w276*w50 + w279*w55;
    const scalar w281 = w256*w50 + w261*w55 + w267*w62 + w271*w66;
    const scalar w282 = w106*w281 + w119*w280 + w244*w58 + w245*w78 - w247*w75 + w273*w95;
    const scalar w283 = w282*w33;
    const scalar w284 = w47*w5;
    const scalar w285 = w5*w64;
    const scalar w286 = -w284*z - w285*z;
    const scalar w287 = w258*w5;
    const scalar w288 = w5*w60;
    const scalar w289 = -w287*z + w288*z;
    const scalar w290 = w286*w62 + w289*w66;
    const scalar w291 = w106*w290 + w119*w281 + w244*w68 + w245*w84 - w247*w78 + w263*w95;
    const scalar w292 = w291*w33;
    const scalar w293 = -w121 - w233*w274 - w283*y - w292*x;
    const scalar w294 = w293*w88;
    const scalar w295 = w131*w282;
    const scalar w296 = w131*w291;
    const scalar w297 = w143 + w236*w274 + w295*y + w296*x;
    const scalar w298 = w142*w149;
    const scalar w299 = w273*w33;
    const scalar w300 = w263*w33;
    const scalar w301 = w23*w233 + w299*y + w300*x + w41;
    const scalar w302 = w233*w301;
    const scalar w303 = w280*w33;
    const scalar w304 = w281*w33;
    const scalar w305 = w33*(w299*z + w303*y + w304*x + w59);
    const scalar w306 = w305*y;
    const scalar w307 = w290*w33;
    const scalar w308 = w33*(w300*z + w304*y + w307*x + w69);
    const scalar w309 = w308*x;
    const scalar w310 = -w302 - w306 - w309 - w70;
    const scalar w311 = (1.0/(w87*w87));
    const scalar w312 = w126*w311;
    const scalar w313 = w310*w312;
    const scalar w314 = w131*w273;
    const scalar w315 = w131*w263;
    const scalar w316 = -w132 - w23*w236 - w314*y - w315*x;
    const scalar w317 = w236*w316;
    const scalar w318 = w131*w280;
    const scalar w319 = w131*w281;
    const scalar w320 = w131*(-w133 - w314*z - w318*y - w319*x);
    const scalar w321 = w320*y;
    const scalar w322 = w131*w290;
    const scalar w323 = w131*(-w134 - w315*z - w319*y - w322*x);
    const scalar w324 = w323*x;
    const scalar w325 = w135 + w317 + w321 + w324;
    const scalar w326 = (1.0/(w141*w141));
    const scalar w327 = w146*w149;
    const scalar w328 = w326*w327;
    const scalar w329 = w127*w235 + w130*w294 + w130*w313 - w147*w237 + w297*w298 + w325*w328;
    const scalar w330 = w153*w32;
    const scalar w331 = k*sigma0Out;
    const scalar w332 = w187*w331;
    const scalar w333 = sigmaEps*w101;
    const scalar w334 = (1.0/2.0)*w27;
    const scalar w335 = w333*w334;
    const scalar w336 = gel*w17;
    const scalar w337 = w336*z;
    const scalar w338 = get*w17;
    const scalar w339 = w338*z;
    const scalar w340 = -w19 + w22 - w337 + w339;
    const scalar w341 = w332*w335 + w340;
    const scalar w342 = gel*w256;
    const scalar w343 = get*w261;
    const scalar w344 = w168*w250;
    const scalar w345 = -2*get*theta1*w13*w14*w3*w65*z + w263 + w342*w51 + w343*w56 + w344*z;
    const scalar w346 = -w345;
    const scalar w347 = w197*w247;
    const scalar w348 = gel*w267;
    const scalar w349 = get*w271;
    const scalar w350 = w161*w250;
    const scalar w351 = -2*get*theta1*w13*w14*w3*w54*z + w273 + w348*w51 + w349*w56 + w350*z;
    const scalar w352 = -w351;
    const scalar w353 = w118*w352 + w346*w99 + w347;
    const scalar w354 = w195*w243 + w341*w94 + w353;
    const scalar w355 = -1.0/2.0*k*sigma0Out*sigmaEps*w101*w187*w27;
    const scalar w356 = w161*w276 + w163*w279 + w280;
    const scalar w357 = -w355 - w356;
    const scalar w358 = w161*w256 + w163*w261 + w168*w267 + w170*w271 + w281;
    const scalar w359 = -w358;
    const scalar w360 = w197*w243 + w352*w94 + w359*w99;
    const scalar w361 = w33*(w118*w357 + w200*w247 + w360);
    const scalar w362 = w168*w286 + w170*w289 + w290;
    const scalar w363 = -w355 - w362;
    const scalar w364 = w118*w359 + w196*w243 + w201*w247 + w346*w94;
    const scalar w365 = w33*(w363*w99 + w364);
    const scalar w366 = w199 + w233*w354 + w361*y + w365*x;
    const scalar w367 = w185*w366;
    const scalar w368 = k*sigma0In;
    const scalar w369 = w187*w368;
    const scalar w370 = w335*w369 + w340;
    const scalar w371 = w221*w243 + w353 + w370*w94;
    const scalar w372 = -1.0/2.0*k*sigma0In*sigmaEps*w101*w187*w27;
    const scalar w373 = -w356 - w372;
    const scalar w374 = w131*(w118*w373 + w223*w247 + w360);
    const scalar w375 = -w362 - w372;
    const scalar w376 = w131*(w364 + w375*w99);
    const scalar w377 = -w222 - w236*w371 - w374*y - w376*x;
    const scalar w378 = w149*w219;
    const scalar w379 = w23 + w337 - w339;
    const scalar w380 = w33*w351;
    const scalar w381 = w33*w345;
    const scalar w382 = w160 + w233*w379 + w380*y + w381*x;
    const scalar w383 = w233*w382;
    const scalar w384 = w33*w356;
    const scalar w385 = w33*w358;
    const scalar w386 = w33*(w167 + w380*z + w384*y + w385*x);
    const scalar w387 = w386*y;
    const scalar w388 = w33*w362;
    const scalar w389 = w33*(w173 + w381*z + w385*y + w388*x);
    const scalar w390 = w389*x;
    const scalar w391 = -w174 - w383 - w387 - w390;
    const scalar w392 = (1.0/(w184*w184));
    const scalar w393 = w207*w392;
    const scalar w394 = w391*w393;
    const scalar w395 = w131*w351;
    const scalar w396 = w131*w345;
    const scalar w397 = -w209 - w236*w379 - w395*y - w396*x;
    const scalar w398 = w236*w397;
    const scalar w399 = w131*w356;
    const scalar w400 = w131*w358;
    const scalar w401 = w131*(-w210 - w395*z - w399*y - w400*x);
    const scalar w402 = w401*y;
    const scalar w403 = w131*w362;
    const scalar w404 = w131*(-w211 - w396*z - w400*y - w403*x);
    const scalar w405 = w404*x;
    const scalar w406 = w212 + w398 + w402 + w405;
    const scalar w407 = (1.0/(w218*w218));
    const scalar w408 = w149*w407;
    const scalar w409 = w227*w408;
    const scalar w410 = w130*w367 + w130*w394 + w208*w235 - w228*w237 + w377*w378 + w406*w409;
    const scalar w411 = w152*z + w156*z + w231*z + w232*z + w24*z + w29 + w329*w330 + w330*w410;
    const scalar w412 = 2*theta1*w12*w3*w34*w42*y - w266*y;
    const scalar w413 = gil*w412;
    const scalar w414 = -w269*y - w270*y;
    const scalar w415 = git*w414;
    const scalar w416 = w264*y - w265*y + w413*w51 + w415*w56;
    const scalar w417 = -w254*x - w255*x;
    const scalar w418 = gil*w51;
    const scalar w419 = -w259*x + w260*x;
    const scalar w420 = git*w56;
    const scalar w421 = w251*x - w253*x + w417*w418 + w419*w420;
    const scalar w422 = -w284*x - w285*x;
    const scalar w423 = -w287*x + w288*x;
    const scalar w424 = w422*w62 + w423*w66;
    const scalar w425 = w33*x;
    const scalar w426 = w234*w425;
    const scalar w427 = w131*x;
    const scalar w428 = w234*w427;
    const scalar w429 = w18*x;
    const scalar w430 = w21*x;
    const scalar w431 = w429 - w430;
    const scalar w432 = St*w98;
    const scalar w433 = (1.0/2.0)*w432;
    const scalar w434 = (3.0/5.0)*St*hA1*k*w109*x - w433;
    const scalar w435 = -St*w103;
    const scalar w436 = w114*w239*w96;
    const scalar w437 = w435 + w436;
    const scalar w438 = 2*theta1*w12*w3*w34*w42*x - w266*x;
    const scalar w439 = -w269*x - w270*x;
    const scalar w440 = w264*x - w265*x + w418*w438 + w420*w439;
    const scalar w441 = w106*w421 + w119*w440 + w245*w40 + w431*w95 + w434*w58 + w437*w68;
    const scalar w442 = w33*w441;
    const scalar w443 = 4*theta1*w12*w3*w34*w42*x - w275*x;
    const scalar w444 = -w277*x - w278*x;
    const scalar w445 = w443*w50 + w444*w55;
    const scalar w446 = w417*w50 + w419*w55 + w438*w62 + w439*w66;
    const scalar w447 = w106*w446 + w119*w445 + w245*w58 + w434*w75 + w437*w78 + w440*w95;
    const scalar w448 = w33*w447;
    const scalar w449 = w434*w78;
    const scalar w450 = w106*w424 + w119*w446 + w246 + w421*w95 + w437*w84 + w449;
    const scalar w451 = -w125 - w425*w450 - w442*z - w448*y;
    const scalar w452 = w451*w88;
    const scalar w453 = w131*w441;
    const scalar w454 = w131*w447;
    const scalar w455 = w145 + w427*w450 + w453*z + w454*y;
    const scalar w456 = w142*w455;
    const scalar w457 = w33*w431;
    const scalar w458 = w33*w440;
    const scalar w459 = w33*w421;
    const scalar w460 = w33*(w457*z + w458*y + w459*x + w69);
    const scalar w461 = w460*z;
    const scalar w462 = w33*w446;
    const scalar w463 = w424*w425 + w459*z + w462*y + w85;
    const scalar w464 = w425*w463;
    const scalar w465 = w33*w445;
    const scalar w466 = w33*(w458*z + w462*x + w465*y + w79);
    const scalar w467 = w466*y;
    const scalar w468 = -w461 - w464 - w467 - w86;
    const scalar w469 = w311*w468;
    const scalar w470 = w126*w469;
    const scalar w471 = w131*w431;
    const scalar w472 = w131*w440;
    const scalar w473 = w131*w421;
    const scalar w474 = w131*(-w134 - w471*z - w472*y - w473*x);
    const scalar w475 = w474*z;
    const scalar w476 = w131*w446;
    const scalar w477 = -w139 - w424*w427 - w473*z - w476*y;
    const scalar w478 = w427*w477;
    const scalar w479 = w131*w445;
    const scalar w480 = w131*(-w137 - w472*z - w476*x - w479*y);
    const scalar w481 = w480*y;
    const scalar w482 = w140 + w475 + w478 + w481;
    const scalar w483 = w326*w482;
    const scalar w484 = w127*w426 + w130*w452 + w130*w470 - w147*w428 + w149*w456 + w327*w483;
    const scalar w485 = cos(w186);
    const scalar w486 = w191*w485;
    const scalar w487 = w336*x;
    const scalar w488 = w338*x;
    const scalar w489 = -w429 + w430 - w487 + w488;
    const scalar w490 = w331*w486 + w489;
    const scalar w491 = gel*w51;
    const scalar w492 = get*w56;
    const scalar w493 = -2*get*theta1*w13*w14*w3*w65*x + w344*x + w417*w491 + w419*w492 + w421;
    const scalar w494 = -w493;
    const scalar w495 = -2*get*theta1*w13*w14*w3*w54*x + w350*x + w438*w491 + w439*w492 + w440;
    const scalar w496 = -w495;
    const scalar w497 = w118*w496 - w196*w436 + w197*w433 + w494*w99;
    const scalar w498 = w33*(w490*w94 + w497);
    const scalar w499 = -k*sigma0Out*sigmaEps*w101*w189*w485;
    const scalar w500 = w161*w443 + w163*w444 + w445;
    const scalar w501 = -w499 - w500;
    const scalar w502 = w161*w417 + w163*w419 + w168*w438 + w170*w439 + w446;
    const scalar w503 = -w502;
    const scalar w504 = -w201*w436 + w496*w94 + w503*w99;
    const scalar w505 = w33*(w118*w501 + w200*w433 + w504);
    const scalar w506 = w168*w422 + w170*w423 + w424;
    const scalar w507 = -w499 - w506;
    const scalar w508 = w201*w433;
    const scalar w509 = w118*w503 + w494*w94 + w508;
    const scalar w510 = -w204*w436 + w507*w99 + w509;
    const scalar w511 = w206 + w425*w510 + w498*z + w505*y;
    const scalar w512 = w185*w511;
    const scalar w513 = w368*w486 + w489;
    const scalar w514 = w131*(w497 + w513*w94);
    const scalar w515 = -k*sigma0In*sigmaEps*w101*w189*w485;
    const scalar w516 = -w500 - w515;
    const scalar w517 = w131*(w118*w516 + w223*w433 + w504);
    const scalar w518 = -w506 - w515;
    const scalar w519 = -w225*w436 + w509 + w518*w99;
    const scalar w520 = -w226 - w427*w519 - w514*z - w517*y;
    const scalar w521 = w431 + w487 - w488;
    const scalar w522 = w33*w521;
    const scalar w523 = w33*w495;
    const scalar w524 = w33*w493;
    const scalar w525 = w33*(w173 + w522*z + w523*y + w524*x);
    const scalar w526 = w525*z;
    const scalar w527 = w33*w502;
    const scalar w528 = w182 + w425*w506 + w524*z + w527*y;
    const scalar w529 = w425*w528;
    const scalar w530 = w33*w500;
    const scalar w531 = w33*(w179 + w523*z + w527*x + w530*y);
    const scalar w532 = w531*y;
    const scalar w533 = -w183 - w526 - w529 - w532;
    const scalar w534 = w393*w533;
    const scalar w535 = w131*w521;
    const scalar w536 = w131*w495;
    const scalar w537 = w131*w493;
    const scalar w538 = w131*(-w211 - w535*z - w536*y - w537*x);
    const scalar w539 = w538*z;
    const scalar w540 = w131*w502;
    const scalar w541 = -w216 - w427*w506 - w537*z - w540*y;
    const scalar w542 = w427*w541;
    const scalar w543 = w131*w500;
    const scalar w544 = w131*(-w214 - w536*z - w540*x - w543*y);
    const scalar w545 = w544*y;
    const scalar w546 = w217 + w539 + w542 + w545;
    const scalar w547 = w130*w512 + w130*w534 + w208*w426 - w228*w428 + w378*w520 + w409*w546;
    const scalar w548 = St*w105 - k*w113*w24 + w152*x + w156*x + w231*x + w232*x + w330*w484 + w330*w547;
    const scalar w549 = -w254*y - w255*y;
    const scalar w550 = -w259*y + w260*y;
    const scalar w551 = w412*w62 + w414*w66 + w50*w549 + w55*w550;
    const scalar w552 = 4*theta1*w12*w3*w34*w42*y - w275*y;
    const scalar w553 = -w277*y - w278*y;
    const scalar w554 = w50*w552 + w55*w553;
    const scalar w555 = w33*y;
    const scalar w556 = w234*w555;
    const scalar w557 = w131*y;
    const scalar w558 = w234*w557;
    const scalar w559 = w18*y - w21*y;
    const scalar w560 = w239*w92;
    const scalar w561 = -w102*w107*w239;
    const scalar w562 = w560 + w561;
    const scalar w563 = gil*w549;
    const scalar w564 = git*w550;
    const scalar w565 = w251*y - w253*y + w51*w563 + w56*w564;
    const scalar w566 = w106*w565 + w119*w416 - w247*w40 + w434*w68 + w559*w95 + w562*w58;
    const scalar w567 = w33*w566;
    const scalar w568 = w106*w551 + w119*w554 + w248 + w416*w95 + w449 + w562*w75;
    const scalar w569 = -w284*y - w285*y;
    const scalar w570 = -w287*y + w288*y;
    const scalar w571 = w569*w62 + w570*w66;
    const scalar w572 = w106*w571 + w119*w551 - w247*w68 + w434*w84 + w562*w78 + w565*w95;
    const scalar w573 = w33*w572;
    const scalar w574 = -w123 - w555*w568 - w567*z - w573*x;
    const scalar w575 = w574*w88;
    const scalar w576 = w131*w566;
    const scalar w577 = w131*w572;
    const scalar w578 = w144 + w557*w568 + w576*z + w577*x;
    const scalar w579 = w33*w559;
    const scalar w580 = w33*w416;
    const scalar w581 = w33*w565;
    const scalar w582 = w33*(w579*z + w580*y + w581*x + w59);
    const scalar w583 = w582*z;
    const scalar w584 = w33*w551;
    const scalar w585 = w554*w555 + w580*z + w584*x + w76;
    const scalar w586 = w555*w585;
    const scalar w587 = w33*w571;
    const scalar w588 = w33*(w581*z + w584*y + w587*x + w79);
    const scalar w589 = w588*x;
    const scalar w590 = -w583 - w586 - w589 - w80;
    const scalar w591 = w312*w590;
    const scalar w592 = w131*w559;
    const scalar w593 = w131*w416;
    const scalar w594 = w131*w565;
    const scalar w595 = w131*(-w133 - w592*z - w593*y - w594*x);
    const scalar w596 = w595*z;
    const scalar w597 = w131*w551;
    const scalar w598 = -w136 - w554*w557 - w593*z - w597*x;
    const scalar w599 = w557*w598;
    const scalar w600 = w131*w571;
    const scalar w601 = w131*(-w137 - w594*z - w597*y - w600*x);
    const scalar w602 = w601*x;
    const scalar w603 = w138 + w596 + w599 + w602;
    const scalar w604 = w127*w556 + w130*w575 + w130*w591 - w147*w558 + w298*w578 + w328*w603;
    const scalar w605 = w109*w190;
    const scalar w606 = w332*w605;
    const scalar w607 = w336*y - w338*y + w559;
    const scalar w608 = -w606 - w607;
    const scalar w609 = gel*w549;
    const scalar w610 = get*w550;
    const scalar w611 = -2*get*theta1*w13*w14*w3*w65*y + w344*y + w51*w609 + w56*w610 + w565;
    const scalar w612 = -w611;
    const scalar w613 = gel*w412;
    const scalar w614 = get*w414;
    const scalar w615 = -2*get*theta1*w13*w14*w3*w54*y + w350*y + w416 + w51*w613 + w56*w614;
    const scalar w616 = -w615;
    const scalar w617 = w118*w616 + w196*w433 - w197*w560 + w612*w99;
    const scalar w618 = w33*(w195*w247 + w608*w94 + w617);
    const scalar w619 = w161*w552 + w163*w553 + w554;
    const scalar w620 = -w606 - w619;
    const scalar w621 = w161*w549 + w163*w550 + w168*w412 + w170*w414 + w551;
    const scalar w622 = -w621;
    const scalar w623 = w347 + w508 + w616*w94 + w622*w99;
    const scalar w624 = w118*w620 - w200*w560 + w623;
    const scalar w625 = w168*w569 + w170*w570 + w571;
    const scalar w626 = -w606 - w625;
    const scalar w627 = w118*w622 + w196*w247 - w201*w560 + w612*w94;
    const scalar w628 = w33*(w204*w433 + w626*w99 + w627);
    const scalar w629 = w203 + w555*w624 + w618*z + w628*x;
    const scalar w630 = w185*w629;
    const scalar w631 = w369*w605;
    const scalar w632 = -w607 - w631;
    const scalar w633 = w131*(w221*w247 + w617 + w632*w94);
    const scalar w634 = -w619 - w631;
    const scalar w635 = w118*w634 - w223*w560 + w623;
    const scalar w636 = -w625 - w631;
    const scalar w637 = w131*(w225*w433 + w627 + w636*w99);
    const scalar w638 = -w224 - w557*w635 - w633*z - w637*x;
    const scalar w639 = w33*w607;
    const scalar w640 = w33*w615;
    const scalar w641 = w33*w611;
    const scalar w642 = w33*(w167 + w639*z + w640*y + w641*x);
    const scalar w643 = w642*z;
    const scalar w644 = w33*w621;
    const scalar w645 = w177 + w555*w619 + w640*z + w644*x;
    const scalar w646 = w555*w645;
    const scalar w647 = w33*w625;
    const scalar w648 = w33*(w179 + w641*z + w644*y + w647*x);
    const scalar w649 = w648*x;
    const scalar w650 = -w180 - w643 - w646 - w649;
    const scalar w651 = w393*w650;
    const scalar w652 = w131*w607;
    const scalar w653 = w131*w615;
    const scalar w654 = w131*w611;
    const scalar w655 = w131*(-w210 - w652*z - w653*y - w654*x);
    const scalar w656 = w655*z;
    const scalar w657 = w131*w621;
    const scalar w658 = -w213 - w557*w619 - w653*z - w657*x;
    const scalar w659 = w557*w658;
    const scalar w660 = w131*w625;
    const scalar w661 = w131*(-w214 - w654*z - w657*y - w660*x);
    const scalar w662 = w661*x;
    const scalar w663 = w215 + w656 + w659 + w662;
    const scalar w664 = w130*w630 + w130*w651 + w208*w556 - w228*w558 + w378*w638 + w409*w663;
    const scalar w665 = (1.0/10.0)*Tt - w112 + w152*y + w156*y + w231*y + w232*y + w330*w604 + w330*w664;
    const scalar w666 = 8*w154;
    const scalar w667 = w666*w8;
    const scalar w668 = 2*w131;
    const scalar w669 = (k*k*k);
    const scalar w670 = St*w669;
    const scalar w671 = (8.0/125.0)*w242;
    const scalar w672 = -w28*w670*x - w671;
    const scalar w673 = 2*w244;
    const scalar w674 = -w240*w68;
    const scalar w675 = (4.0/25.0)*St*w117;
    const scalar w676 = w58*w675;
    const scalar w677 = (theta1*theta1);
    const scalar w678 = (1.0/(w2*w2));
    const scalar w679 = w677*w678;
    const scalar w680 = w679*w8;
    const scalar w681 = 8*w680;
    const scalar w682 = w38*w681;
    const scalar w683 = gil*w682;
    const scalar w684 = w36*w681;
    const scalar w685 = git*w684;
    const scalar w686 = w37*w681;
    const scalar w687 = w39*w681;
    const scalar w688 = -w18 + w21;
    const scalar w689 = -w683 - w685 + w686 + w687 + w688;
    const scalar w690 = -w689;
    const scalar w691 = w245*w263;
    const scalar w692 = w247*w273;
    const scalar w693 = 4*w680;
    const scalar w694 = w254 + w255;
    const scalar w695 = 4*w13*w34*w45*w677*w678*w8 - w258*w693 - w694;
    const scalar w696 = -w259 + w260;
    const scalar w697 = w47*w693 + w64*w693 + w696;
    const scalar w698 = w5*z;
    const scalar w699 = w56*w698;
    const scalar w700 = w51*w698;
    const scalar w701 = -w251 + w253;
    const scalar w702 = -w257*w699 + w262*w700 - w418*w695 - w420*w697 + w63*w693 + w67*w693 + w701;
    const scalar w703 = -w702;
    const scalar w704 = -2*theta1*w12*w3*w34*w42 + w266;
    const scalar w705 = -w46*w693 - w48*w693 - w704;
    const scalar w706 = w258*w34;
    const scalar w707 = w269 + w270;
    const scalar w708 = 4*w13*w45*w677*w678*w8 - w693*w706 - w707;
    const scalar w709 = -w264 + w265;
    const scalar w710 = -w268*w699 + w272*w700 - w418*w705 - w420*w708 + w52*w693 + w57*w693 + w709;
    const scalar w711 = -w710;
    const scalar w712 = w106*w703 + w119*w711 + w23*w673 + w40*w672 + w674 + w676 + w690*w95 + 2*w691 - 2*w692;
    const scalar w713 = 2*w247;
    const scalar w714 = 2*w245;
    const scalar w715 = w267*w276;
    const scalar w716 = w271*w279;
    const scalar w717 = -4*theta1*w12*w3*w34*w42 + w275;
    const scalar w718 = -w46*w681 - w48*w681 - w717;
    const scalar w719 = w277 + w278;
    const scalar w720 = 8*w13*w45*w677*w678*w8 - w681*w706 - w719;
    const scalar w721 = gil*w715 + git*w716 + w50*w718 + w55*w720;
    const scalar w722 = 2*w267;
    const scalar w723 = 2*w271;
    const scalar w724 = w257*w722 + w262*w723 + w50*w695 + w55*w697 + w62*w705 + w66*w708;
    const scalar w725 = w106*w724 + w119*w721 - w240*w78 + w273*w673 - w280*w713 + w281*w714 + w58*w672 + w675*w75 + w711*w95;
    const scalar w726 = w256*w286;
    const scalar w727 = w261*w289;
    const scalar w728 = w284 + w285;
    const scalar w729 = 8*w13*w34*w45*w677*w678*w8 - w258*w681 - w728;
    const scalar w730 = -w287 + w288;
    const scalar w731 = w47*w681 + w64*w681 + w730;
    const scalar w732 = gil*w726 + git*w727 + w62*w729 + w66*w731;
    const scalar w733 = w106*w732 + w119*w724 - w240*w84 + w263*w673 - w281*w713 + w290*w714 + w672*w68 + w675*w78 + w703*w95;
    const scalar w734 = w148*w30;
    const scalar w735 = w326*w734;
    const scalar w736 = 2*w33;
    const scalar w737 = w130*w88;
    const scalar w738 = w142*w237;
    const scalar w739 = 2*w235;
    const scalar w740 = w129*w30;
    const scalar w741 = w311*w740;
    const scalar w742 = 2*w135 + 2*w317 + 2*w321 + 2*w324;
    const scalar w743 = w327/(w141*w141*w141);
    const scalar w744 = w325*w743;
    const scalar w745 = w130*w312;
    const scalar w746 = w325*w326;
    const scalar w747 = w146*w237;
    const scalar w748 = -2*w302 - 2*w306 - 2*w309 - 2*w70;
    const scalar w749 = w126*w130/(w87*w87*w87);
    const scalar w750 = w310*w749;
    const scalar w751 = w234*w33;
    const scalar w752 = w131*w234;
    const scalar w753 = w127*w751 - w147*w752;
    const scalar w754 = sigma0In*w238;
    const scalar w755 = w187*w754;
    const scalar w756 = w188*w333;
    const scalar w757 = w755*w756;
    const scalar w758 = w157*w681;
    const scalar w759 = w158*w681;
    const scalar w760 = gel*w682;
    const scalar w761 = get*w684;
    const scalar w762 = w18 - w21 + w336 - w338;
    const scalar w763 = w683 + w685 - w686 - w687 - w758 - w759 + w760 + w761 + w762;
    const scalar w764 = (8.0/25.0)*w242;
    const scalar w765 = w170*w252 - w344;
    const scalar w766 = w169*w693 + w171*w693 - w342*w699 + w343*w700 - w491*w695 - w492*w697 + w702 + w765;
    const scalar w767 = w163*w252 - w350;
    const scalar w768 = w162*w693 + w164*w693 - w348*w699 + w349*w700 - w491*w705 - w492*w708 + w710 + w767;
    const scalar w769 = w247*w352;
    const scalar w770 = -w197*w675;
    const scalar w771 = w118*w768 + w766*w99 + 2*w769 + w770;
    const scalar w772 = gel*w715 + get*w716 + w161*w718 + w163*w720 + w721;
    const scalar w773 = w161*w695 + w163*w697 + w168*w705 + w170*w708 + w342*w722 + w343*w723 + w724;
    const scalar w774 = -w773;
    const scalar w775 = w197*w671 + w352*w764 + w768*w94 + w774*w99;
    const scalar w776 = gel*w726 + get*w727 + w168*w729 + w170*w731 + w732;
    const scalar w777 = w118*w774 + w196*w671 - w201*w675 + w346*w764 + w359*w713 + w766*w94;
    const scalar w778 = w407*w734;
    const scalar w779 = -w336 + w338;
    const scalar w780 = -w689 - w758 - w759 + w760 + w761 - w779;
    const scalar w781 = -w768;
    const scalar w782 = -w766;
    const scalar w783 = sigma0Out*w238;
    const scalar w784 = w187*w783;
    const scalar w785 = w756*w784;
    const scalar w786 = w130*w185;
    const scalar w787 = w219*w237;
    const scalar w788 = w392*w740;
    const scalar w789 = 2*w212 + 2*w398 + 2*w402 + 2*w405;
    const scalar w790 = w149*w227/(w218*w218*w218);
    const scalar w791 = w406*w790;
    const scalar w792 = w130*w393;
    const scalar w793 = w227*w407;
    const scalar w794 = w237*w793;
    const scalar w795 = -2*w174 - 2*w383 - 2*w387 - 2*w390;
    const scalar w796 = w130*w207/(w184*w184*w184);
    const scalar w797 = w391*w796;
    const scalar w798 = w208*w751 - w228*w752;
    const scalar w799 = 4*z;
    const scalar w800 = w32*w329;
    const scalar w801 = w410*w799;
    const scalar w802 = w155*w329;
    const scalar w803 = w152 + w156 + w231 + w232;
    const scalar w804 = w6*w666;
    const scalar w805 = fA2*w670*w96*w97;
    const scalar w806 = w115*w239;
    const scalar w807 = (3.0/5.0)*w111 + w806;
    const scalar w808 = w6*w679;
    const scalar w809 = 8*w808;
    const scalar w810 = w38*w809;
    const scalar w811 = gil*w810;
    const scalar w812 = w36*w809;
    const scalar w813 = git*w812;
    const scalar w814 = w37*w809;
    const scalar w815 = w39*w809;
    const scalar w816 = w688 - w811 - w813 + w814 + w815;
    const scalar w817 = -w816;
    const scalar w818 = 2*w434;
    const scalar w819 = 2*w437;
    const scalar w820 = 4*w808;
    const scalar w821 = 4*w13*w34*w45*w6*w677*w678 - w258*w820 - w694;
    const scalar w822 = w47*w820 + w64*w820 + w696;
    const scalar w823 = gil*w417;
    const scalar w824 = w5*x;
    const scalar w825 = w56*w824;
    const scalar w826 = git*w419;
    const scalar w827 = w51*w824;
    const scalar w828 = -w418*w821 - w420*w822 + w63*w820 + w67*w820 + w701 - w823*w825 + w826*w827;
    const scalar w829 = -w828;
    const scalar w830 = -w46*w820 - w48*w820 - w704;
    const scalar w831 = 4*w13*w45*w6*w677*w678 - w706*w820 - w707;
    const scalar w832 = w438*w825;
    const scalar w833 = -gil*w832 + w20*w439*w51*x - w418*w830 - w420*w831 + w52*w820 + w57*w820 + w709;
    const scalar w834 = -w833;
    const scalar w835 = w106*w829 + w119*w834 + w421*w819 + w431*w714 + w440*w818 + w58*w807 + w68*w805 + w817*w95;
    const scalar w836 = w438*w443;
    const scalar w837 = w439*w444;
    const scalar w838 = -w46*w809 - w48*w809 - w717;
    const scalar w839 = 8*w13*w45*w6*w677*w678 - w706*w809 - w719;
    const scalar w840 = gil*w836 + git*w837 + w50*w838 + w55*w839;
    const scalar w841 = 2*w438;
    const scalar w842 = 2*w439;
    const scalar w843 = w50*w821 + w55*w822 + w62*w830 + w66*w831 + w823*w841 + w826*w842;
    const scalar w844 = w106*w843 + w119*w840 + w440*w714 + w445*w818 + w446*w819 + w75*w807 + w78*w805 + w834*w95;
    const scalar w845 = w78*w807;
    const scalar w846 = w245*w421;
    const scalar w847 = w434*w446;
    const scalar w848 = w417*w422;
    const scalar w849 = w419*w423;
    const scalar w850 = 8*w13*w34*w45*w6*w677*w678 - w258*w809 - w728;
    const scalar w851 = w47*w809 + w64*w809 + w730;
    const scalar w852 = gil*w848 + git*w849 + w62*w850 + w66*w851;
    const scalar w853 = w106*w852 + w119*w843 + w424*w819 + w805*w84 + w829*w95 + w845 + 2*w846 + 2*w847;
    const scalar w854 = 2*w428;
    const scalar w855 = w426*w88;
    const scalar w856 = 2*w140 + 2*w475 + 2*w478 + 2*w481;
    const scalar w857 = w482*w743;
    const scalar w858 = w146*w483;
    const scalar w859 = -2*w461 - 2*w464 - 2*w467 - 2*w86;
    const scalar w860 = w468*w749;
    const scalar w861 = w192*w238;
    const scalar w862 = sigma0In*w861;
    const scalar w863 = w762 + w862;
    const scalar w864 = w157*w809;
    const scalar w865 = w158*w809;
    const scalar w866 = gel*w810;
    const scalar w867 = get*w812;
    const scalar w868 = w811 + w813 - w814 - w815 - w864 - w865 + w866 + w867;
    const scalar w869 = gel*w417;
    const scalar w870 = get*w419;
    const scalar w871 = w169*w820 + w171*w820 - w491*w821 - w492*w822 + w765 - w825*w869 + w827*w870 + w828;
    const scalar w872 = -gel*w832 + get*w439*w827 + w162*w820 + w164*w820 - w491*w830 - w492*w831 + w767 + w833;
    const scalar w873 = 2*w436;
    const scalar w874 = w118*w872 - w196*w805 - w197*w806 + w432*w496 - w494*w873 + w871*w99;
    const scalar w875 = gel*w836 + get*w837 + w161*w838 + w163*w839 + w840;
    const scalar w876 = w161*w821 + w163*w822 + w168*w830 + w170*w831 + w841*w869 + w842*w870 + w843;
    const scalar w877 = -w876;
    const scalar w878 = -w201*w805 - w503*w873 + w872*w94 + w877*w99;
    const scalar w879 = gel*w848 + get*w849 + w168*w850 + w170*w851 + w852;
    const scalar w880 = -w201*w806;
    const scalar w881 = w118*w877 + w432*w503 + w871*w94 + w880;
    const scalar w882 = -w779 - w816 - w864 - w865 + w866 + w867;
    const scalar w883 = -w872;
    const scalar w884 = -w871;
    const scalar w885 = sigma0Out*w861;
    const scalar w886 = w762 + w885;
    const scalar w887 = w219*w428;
    const scalar w888 = w185*w426;
    const scalar w889 = 2*w217 + 2*w539 + 2*w542 + 2*w545;
    const scalar w890 = w546*w790;
    const scalar w891 = w546*w793;
    const scalar w892 = w393*w426;
    const scalar w893 = -2*w183 - 2*w526 - 2*w529 - 2*w532;
    const scalar w894 = w533*w796;
    const scalar w895 = 4*x;
    const scalar w896 = w32*w484;
    const scalar w897 = w547*w895;
    const scalar w898 = w155*w484;
    const scalar w899 = w666*w7;
    const scalar w900 = w239*w93;
    const scalar w901 = w104*w239;
    const scalar w902 = w116*w670;
    const scalar w903 = St*hA1*w107*w109*w669 - w902;
    const scalar w904 = w679*w7;
    const scalar w905 = 8*w904;
    const scalar w906 = w38*w905;
    const scalar w907 = gil*w906;
    const scalar w908 = w36*w905;
    const scalar w909 = git*w908;
    const scalar w910 = w37*w905;
    const scalar w911 = w39*w905;
    const scalar w912 = w688 - w907 - w909 + w910 + w911;
    const scalar w913 = -w912;
    const scalar w914 = 2*w562;
    const scalar w915 = 4*w904;
    const scalar w916 = 4*w13*w34*w45*w677*w678*w7 - w258*w915 - w694;
    const scalar w917 = w47*w915 + w64*w915 + w696;
    const scalar w918 = w5*y;
    const scalar w919 = w56*w918;
    const scalar w920 = w51*w918;
    const scalar w921 = -w418*w916 - w420*w917 - w563*w919 + w564*w920 + w63*w915 + w67*w915 + w701;
    const scalar w922 = -w921;
    const scalar w923 = -w46*w915 - w48*w915 - w704;
    const scalar w924 = 4*w13*w45*w677*w678*w7 - w706*w915 - w707;
    const scalar w925 = -w413*w919 + w415*w920 - w418*w923 - w420*w924 + w52*w915 + w57*w915 + w709;
    const scalar w926 = -w925;
    const scalar w927 = w106*w922 + w119*w926 + w40*w900 + w416*w914 - w559*w713 + w565*w818 + w58*w903 + w68*w901 + w913*w95;
    const scalar w928 = w58*w900;
    const scalar w929 = w78*w901;
    const scalar w930 = w247*w416;
    const scalar w931 = w434*w551;
    const scalar w932 = w412*w552;
    const scalar w933 = w414*w553;
    const scalar w934 = -w46*w905 - w48*w905 - w717;
    const scalar w935 = 8*w13*w45*w677*w678*w7 - w706*w905 - w719;
    const scalar w936 = gil*w932 + git*w933 + w50*w934 + w55*w935;
    const scalar w937 = 2*w412;
    const scalar w938 = 2*w414;
    const scalar w939 = w50*w916 + w55*w917 + w563*w937 + w564*w938 + w62*w923 + w66*w924;
    const scalar w940 = w106*w939 + w119*w936 + w554*w914 + w75*w903 + w926*w95 + w928 + w929 - 2*w930 + 2*w931;
    const scalar w941 = w549*w569;
    const scalar w942 = w550*w570;
    const scalar w943 = 8*w13*w34*w45*w677*w678*w7 - w258*w905 - w728;
    const scalar w944 = w47*w905 + w64*w905 + w730;
    const scalar w945 = gil*w941 + git*w942 + w62*w943 + w66*w944;
    const scalar w946 = w106*w945 + w119*w939 + w551*w914 - w565*w713 + w571*w818 + w68*w900 + w78*w903 + w84*w901 + w922*w95;
    const scalar w947 = w142*w558;
    const scalar w948 = 2*w556;
    const scalar w949 = 2*w138 + 2*w596 + 2*w599 + 2*w602;
    const scalar w950 = w603*w743;
    const scalar w951 = w326*w603;
    const scalar w952 = w146*w558;
    const scalar w953 = -2*w583 - 2*w586 - 2*w589 - 2*w80;
    const scalar w954 = w590*w749;
    const scalar w955 = w157*w905;
    const scalar w956 = w158*w905;
    const scalar w957 = gel*w906;
    const scalar w958 = get*w908;
    const scalar w959 = w907 + w909 - w910 - w911 - w955 - w956 + w957 + w958;
    const scalar w960 = w169*w915 + w171*w915 - w491*w916 - w492*w917 - w609*w919 + w610*w920 + w765 + w921;
    const scalar w961 = w162*w915 + w164*w915 - w491*w923 - w492*w924 - w613*w919 + w614*w920 + w767 + w925;
    const scalar w962 = 2*w560;
    const scalar w963 = w118*w961 + w197*w902 + w432*w612 - w616*w962 + w960*w99;
    const scalar w964 = gel*w932 + get*w933 + w161*w934 + w163*w935 + w936;
    const scalar w965 = w161*w916 + w163*w917 + w168*w923 + w170*w924 + w609*w937 + w610*w938 + w939;
    const scalar w966 = -w965;
    const scalar w967 = w247*w616;
    const scalar w968 = -w197*w900;
    const scalar w969 = w432*w622 + w94*w961 + w966*w99 + 2*w967 + w968;
    const scalar w970 = gel*w941 + get*w942 + w168*w943 + w170*w944 + w945;
    const scalar w971 = w118*w966 - w196*w900 + w201*w902 + w612*w713 - w622*w962 + w94*w960;
    const scalar w972 = -w779 - w912 - w955 - w956 + w957 + w958;
    const scalar w973 = -w961;
    const scalar w974 = -w960;
    const scalar w975 = w219*w558;
    const scalar w976 = 2*w215 + 2*w656 + 2*w659 + 2*w662;
    const scalar w977 = w663*w790;
    const scalar w978 = w558*w793;
    const scalar w979 = -2*w180 - 2*w643 - 2*w646 - 2*w649;
    const scalar w980 = w650*w796;
    const scalar w981 = 4*y;
    const scalar w982 = w32*w604;
    const scalar w983 = w664*w981;
    const scalar w984 = w155*w604;
    const scalar w985 = 8*x;
    const scalar w986 = w679*z;
    const scalar w987 = w985*w986;
    const scalar w988 = -w258*w987 + w60*w987;
    const scalar w989 = w47*w987 + w64*w987;
    const scalar w990 = w62*w988 + w66*w989;
    const scalar w991 = w257*w422 + w262*w423 + w990;
    const scalar w992 = w895*w986;
    const scalar w993 = -w258*w992 + w60*w992;
    const scalar w994 = w47*w992 + w64*w992;
    const scalar w995 = w250*x;
    const scalar w996 = w250*z;
    const scalar w997 = w252*x;
    const scalar w998 = w252*z;
    const scalar w999 = -w257*w995 + w262*w997 - w418*w993 - w420*w994 + w63*w992 + w67*w992 - w823*w996 + w826*w998;
    const scalar w1000 = -w999;
    const scalar w1001 = -w46*w992 - w48*w992;
    const scalar w1002 = w53*w992 - w706*w992;
    const scalar w1003 = w1001*w62 + w1002*w66 + w257*w438 + w262*w439 + w268*w417 + w272*w419 + w50*w993 + w55*w994;
    const scalar w1004 = w281*w434;
    const scalar w1005 = w1000*w95 + w1003*w119 + w1004 + w244*w421 + w245*w424 - w247*w446 + w290*w437 + w674 + w691;
    const scalar w1006 = w1005 + w106*w991;
    const scalar w1007 = -w46*w987 - w48*w987;
    const scalar w1008 = w53*w987 - w706*w987;
    const scalar w1009 = w1007*w50 + w1008*w55;
    const scalar w1010 = w1009 + w268*w443 + w272*w444;
    const scalar w1011 = w56*w841;
    const scalar w1012 = w4*z;
    const scalar w1013 = w1011*w1012;
    const scalar w1014 = w51*w842;
    const scalar w1015 = w1012*w1014;
    const scalar w1016 = -gil*w1013 + git*w1015 - w1001*w418 - w1002*w420 - w268*w995 + w272*w997 + w52*w992 + w57*w992;
    const scalar w1017 = -w1016;
    const scalar w1018 = w1003*w106 + w1017*w95 - w240*w58 + w244*w440 + w245*w273 + w245*w446 - w247*w445 + w280*w434 + w281*w437;
    const scalar w1019 = w1010*w119 + w1018;
    const scalar w1020 = w38*w987;
    const scalar w1021 = w36*w987;
    const scalar w1022 = -gil*w1020 - git*w1021 + w37*w987 + w39*w987;
    const scalar w1023 = -w1022;
    const scalar w1024 = -w247*w440;
    const scalar w1025 = w1000*w106 + w1017*w119 + w1023*w95 + w1024 + w23*w245 - w240*w40 + w244*w431 + w263*w437 + w273*w434 + w846;
    const scalar w1026 = w1025*w236 + w296 + w453;
    const scalar w1027 = w1000*w236 + w1003*w557 + w322 + w473;
    const scalar w1028 = w1003*w427 + w1017*w236 + w319 + w472;
    const scalar w1029 = w236*(-w1000*w427 - w1017*w557 - w1023*w236 - w315 - w471) + w323 + w474;
    const scalar w1030 = w1025*w233 + w292 + w442;
    const scalar w1031 = w1000*w233 + w1003*w555 + w307 + w459;
    const scalar w1032 = w1003*w425 + w1017*w233 + w304 + w458;
    const scalar w1033 = w233*(w1000*w425 + w1017*w555 + w1023*w233 + w300 + w457) + w308 + w460;
    const scalar w1034 = w149*w483;
    const scalar w1035 = w149*w455;
    const scalar w1036 = w142*w428;
    const scalar w1037 = w312*w426;
    const scalar w1038 = w130*w469;
    const scalar w1039 = w130*w311;
    const scalar w1040 = w1039*w451;
    const scalar w1041 = w146*w428;
    const scalar w1042 = w1034*w297 + w1035*w746 - w1036*w297 + w1037*w310 + w1038*w293 + w1040*w310 - w1041*w746 + w235*w452 + w235*w470 - w237*w456 - w237*w858 + w293*w855;
    const scalar w1043 = -1.0/2.0*sigma0In*sigmaEps*w101*w238*w27*w485;
    const scalar w1044 = w168*w988 + w170*w989;
    const scalar w1045 = w1044 + w342*w422 + w343*w423 + w991;
    const scalar w1046 = w169*w992 + w171*w992 - w342*w995 + w343*w997 - w491*w993 - w492*w994 - w869*w996 + w870*w998 + w999;
    const scalar w1047 = w1001*w168 + w1002*w170 + w1003 + w161*w993 + w163*w994 + w342*w438 + w343*w439 + w348*w417 + w349*w419;
    const scalar w1048 = -w1047;
    const scalar w1049 = w359*w433;
    const scalar w1050 = w1046*w94 + w1048*w118 + w1049 + w243*w494 + w247*w503;
    const scalar w1051 = w1050 - w375*w436;
    const scalar w1052 = w1007*w161 + w1008*w163;
    const scalar w1053 = w1010 + w1052 + w348*w443 + w349*w444;
    const scalar w1054 = -gel*w1013 + get*w1015 - w1001*w491 - w1002*w492 + w1016 + w162*w992 + w164*w992 - w348*w995 + w349*w997;
    const scalar w1055 = w1048*w99 + w1054*w94 + w243*w496 - w359*w436;
    const scalar w1056 = w1055 + w247*w516 + w373*w433;
    const scalar w1057 = w335*w485;
    const scalar w1058 = -gel*w1020 - get*w1021 + w1022 + w157*w987 + w158*w987;
    const scalar w1059 = w247*w496;
    const scalar w1060 = w1046*w99 + w1054*w118 + w1059 - w346*w436 + w352*w433;
    const scalar w1061 = w236*(w1060 + w243*w513 + w94*(w1057*w754 + w1058)) + w376 + w514;
    const scalar w1062 = -w1046;
    const scalar w1063 = w1047*w557 + w1062*w236 + w403 + w537;
    const scalar w1064 = -w1054;
    const scalar w1065 = w1047*w427 + w1064*w236 + w400 + w536;
    const scalar w1066 = -w1058;
    const scalar w1067 = w236*(-w1062*w427 - w1064*w557 - w1066*w236 - w396 - w535) + w404 + w538;
    const scalar w1068 = -1.0/2.0*sigma0Out*sigmaEps*w101*w238*w27*w485;
    const scalar w1069 = w1050 - w363*w436;
    const scalar w1070 = w1055 + w247*w501 + w357*w433;
    const scalar w1071 = w233*(w1060 + w243*w490 + w94*(w1057*w783 + w1058)) + w365 + w498;
    const scalar w1072 = w1047*w555 + w1062*w233 + w388 + w524;
    const scalar w1073 = w1047*w425 + w1064*w233 + w385 + w523;
    const scalar w1074 = w233*(w1062*w425 + w1064*w555 + w1066*w233 + w381 + w522) + w389 + w525;
    const scalar w1075 = w408*w546;
    const scalar w1076 = w408*w520;
    const scalar w1077 = w219*w520;
    const scalar w1078 = w130*w392;
    const scalar w1079 = w1078*w533;
    const scalar w1080 = w1078*w511;
    const scalar w1081 = w428*w793;
    const scalar w1082 = w1075*w377 + w1076*w406 - w1077*w237 + w1079*w366 + w1080*w391 - w1081*w406 + w235*w512 + w235*w534 - w237*w891 + w366*w888 - w377*w887 + w391*w892;
    const scalar w1083 = 2*x;
    const scalar w1084 = w1083*w410;
    const scalar w1085 = 2*z;
    const scalar w1086 = w1085*w547;
    const scalar w1087 = w666*x;
    const scalar w1088 = w1087*z;
    const scalar w1089 = w1083*w800 + w1083*w802 + w1084*w155 + w1084*w32 + w1085*w896 + w1085*w898 + w1086*w155 + w1086*w32 + w1088*w150 + w1088*w229 + w245;
    const scalar w1090 = w286*w417;
    const scalar w1091 = w289*w419;
    const scalar w1092 = gil*w1090 + git*w1091 + w990;
    const scalar w1093 = w1005 + w106*w1092;
    const scalar w1094 = w276*w438;
    const scalar w1095 = w279*w439;
    const scalar w1096 = gil*w1094 + git*w1095 + w1009;
    const scalar w1097 = w1018 + w1096*w119;
    const scalar w1098 = gel*w1090 + get*w1091 + w1044 + w1092;
    const scalar w1099 = gel*w1094 + get*w1095 + w1052 + w1096;
    const scalar w1100 = y*z;
    const scalar w1101 = 8*w1100*w679;
    const scalar w1102 = -w1101*w258 + w1101*w60;
    const scalar w1103 = w1101*w47 + w1101*w64;
    const scalar w1104 = w1102*w62 + w1103*w66;
    const scalar w1105 = w1104 + w257*w569 + w262*w570;
    const scalar w1106 = w981*w986;
    const scalar w1107 = -w1106*w258 + w1106*w60;
    const scalar w1108 = w1106*w47 + w1106*w64;
    const scalar w1109 = w250*y;
    const scalar w1110 = w252*y;
    const scalar w1111 = w1106*w63 + w1106*w67 - w1107*w418 - w1108*w420 - w1109*w257 + w1110*w262 - w563*w996 + w564*w998;
    const scalar w1112 = -w1111;
    const scalar w1113 = -w1106*w46 - w1106*w48;
    const scalar w1114 = w1106*w53 - w1106*w706;
    const scalar w1115 = w1107*w50 + w1108*w55 + w1113*w62 + w1114*w66 + w256*w413 + w261*w415 + w267*w563 + w271*w564;
    const scalar w1116 = w1112*w95 + w1115*w119 + w244*w565 + w245*w571 - w247*w263 - w247*w551 + w281*w562 + w290*w434 + w675*w68 + w78*w900;
    const scalar w1117 = w106*w1105 + w1116;
    const scalar w1118 = -w1101*w46 - w1101*w48;
    const scalar w1119 = w1101*w53 - w1101*w706;
    const scalar w1120 = w1118*w50 + w1119*w55;
    const scalar w1121 = w1120 + w268*w552 + w272*w553;
    const scalar w1122 = w1106*w52 + w1106*w57 - w1109*w268 + w1110*w272 - w1113*w418 - w1114*w420 - w413*w996 + w415*w998;
    const scalar w1123 = -w1122;
    const scalar w1124 = w1004 + w106*w1115 + w1123*w95 + w244*w416 + w245*w551 - w247*w554 + w280*w562 + w676 - w692 + w75*w900;
    const scalar w1125 = w1121*w119 + w1124;
    const scalar w1126 = w1101*w37;
    const scalar w1127 = w1101*w38;
    const scalar w1128 = gil*w1127;
    const scalar w1129 = w1101*w36;
    const scalar w1130 = git*w1129;
    const scalar w1131 = w1101*w39;
    const scalar w1132 = w1126 - w1128 - w1130 + w1131;
    const scalar w1133 = -w1132;
    const scalar w1134 = w245*w565;
    const scalar w1135 = w106*w1112 + w1123*w119 + w1133*w95 + w1134 - w23*w247 + w244*w559 + w263*w434 + w273*w562 + w40*w675 + w928 - w930;
    const scalar w1136 = w1135*w236 + w295 + w576;
    const scalar w1137 = w1112*w236 + w1115*w557 + w319 + w594;
    const scalar w1138 = w1115*w427 + w1123*w236 + w318 + w593;
    const scalar w1139 = w236*(-w1112*w427 - w1123*w557 - w1133*w236 - w314 - w592) + w320 + w595;
    const scalar w1140 = w1135*w233 + w283 + w567;
    const scalar w1141 = w1112*w233 + w1115*w555 + w304 + w581;
    const scalar w1142 = w1115*w425 + w1123*w233 + w303 + w580;
    const scalar w1143 = w233*(w1112*w425 + w1123*w555 + w1133*w233 + w299 + w579) + w305 + w582;
    const scalar w1144 = w1039*w293*w590 + w1039*w310*w574 + w149*w297*w951 + w149*w578*w746 + w235*w575 + w235*w591 + w294*w556 - w297*w947 + w313*w556 - w578*w738 - w746*w952 - w747*w951;
    const scalar w1145 = sigmaEps*w109*w334;
    const scalar w1146 = w1145*w755;
    const scalar w1147 = w1102*w168 + w1103*w170;
    const scalar w1148 = w1105 + w1147 + w342*w569 + w343*w570;
    const scalar w1149 = w1106*w169 + w1106*w171 - w1107*w491 - w1108*w492 - w1109*w342 + w1110*w343 + w1111 - w609*w996 + w610*w998;
    const scalar w1150 = w1107*w161 + w1108*w163 + w1113*w168 + w1114*w170 + w1115 + w256*w613 + w261*w614 + w267*w609 + w271*w610;
    const scalar w1151 = -w1150;
    const scalar w1152 = w1149*w94 + w1151*w118 - w196*w675 - w201*w900 + w243*w612 + w247*w346 + w247*w622 - w359*w560;
    const scalar w1153 = w1152 + w375*w433;
    const scalar w1154 = w1118*w161 + w1119*w163;
    const scalar w1155 = w1121 + w1154 + w348*w552 + w349*w553;
    const scalar w1156 = w1106*w162 + w1106*w164 - w1109*w348 + w1110*w349 - w1113*w491 - w1114*w492 + w1122 - w613*w996 + w614*w998;
    const scalar w1157 = w1049 + w1151*w99 + w1156*w94 + w243*w616 + w769 + w770;
    const scalar w1158 = w1157 - w223*w900 + w247*w634 - w373*w560;
    const scalar w1159 = w1101*w157;
    const scalar w1160 = w1101*w158;
    const scalar w1161 = gel*w1127;
    const scalar w1162 = get*w1129;
    const scalar w1163 = -w1126 + w1128 + w1130 - w1131 - w1159 - w1160 + w1161 + w1162;
    const scalar w1164 = w1149*w99 + w1156*w118 + w346*w433 - w352*w560 + w967 + w968;
    const scalar w1165 = w236*(w1164 - w221*w675 + w243*w632 + w247*w370 + w94*(-w1146 - w1163)) + w374 + w633;
    const scalar w1166 = -w1149;
    const scalar w1167 = w1150*w557 + w1166*w236 + w400 + w654;
    const scalar w1168 = -w1156;
    const scalar w1169 = w1150*w427 + w1168*w236 + w399 + w653;
    const scalar w1170 = -w1132 - w1159 - w1160 + w1161 + w1162;
    const scalar w1171 = w236*(-w1166*w427 - w1168*w557 - w1170*w236 - w395 - w652) + w401 + w655;
    const scalar w1172 = w1145*w784;
    const scalar w1173 = w1152 + w363*w433;
    const scalar w1174 = w1157 - w200*w900 + w247*w620 - w357*w560;
    const scalar w1175 = w233*(w1164 - w195*w675 + w243*w608 + w247*w341 + w94*(-w1163 - w1172)) + w361 + w618;
    const scalar w1176 = w1150*w555 + w1166*w233 + w385 + w641;
    const scalar w1177 = w1150*w425 + w1168*w233 + w384 + w640;
    const scalar w1178 = w233*(w1166*w425 + w1168*w555 + w1170*w233 + w380 + w639) + w386 + w642;
    const scalar w1179 = w1078*w366*w650 + w1078*w391*w629 + w235*w630 + w235*w651 + w367*w556 + w377*w408*w663 - w377*w975 + w394*w556 + w406*w408*w638 - w406*w978 - w638*w787 - w663*w794;
    const scalar w1180 = 2*y;
    const scalar w1181 = w1180*w410;
    const scalar w1182 = w1085*w664;
    const scalar w1183 = w1100*w666;
    const scalar w1184 = w1085*w982 + w1085*w984 + w1180*w800 + w1180*w802 + w1181*w155 + w1181*w32 + w1182*w155 + w1182*w32 + w1183*w150 + w1183*w229;
    const scalar w1185 = w1104 + w286*w563 + w289*w564;
    const scalar w1186 = w106*w1185 + w1116;
    const scalar w1187 = w1120 + w276*w413 + w279*w415;
    const scalar w1188 = w1124 + w1187*w119;
    const scalar w1189 = w1147 + w1185 + w286*w609 + w289*w610;
    const scalar w1190 = w1154 + w1187 + w276*w613 + w279*w614;
    const scalar w1191 = w679*y;
    const scalar w1192 = w1191*w985;
    const scalar w1193 = -w1192*w258 + w1192*w60;
    const scalar w1194 = w1192*w47 + w1192*w64;
    const scalar w1195 = w1193*w62 + w1194*w66;
    const scalar w1196 = w1195 + w422*w563 + w423*w564;
    const scalar w1197 = w1191*w895;
    const scalar w1198 = -w1197*w258 + w1197*w60;
    const scalar w1199 = w1197*w47 + w1197*w64;
    const scalar w1200 = -w1109*w823 + w1110*w826 + w1197*w63 + w1197*w67 - w1198*w418 - w1199*w420 - w563*w995 + w564*w997;
    const scalar w1201 = -w1200;
    const scalar w1202 = -w1197*w46 - w1197*w48;
    const scalar w1203 = w1197*w53 - w1197*w706;
    const scalar w1204 = w1198*w50 + w1199*w55 + w1202*w62 + w1203*w66 + w413*w417 + w415*w419 + w438*w563 + w439*w564;
    const scalar w1205 = w1134 + w119*w1204 + w1201*w95 - w247*w421 + w424*w434 + w437*w571 + w446*w562 + w807*w84 + w929 + w931;
    const scalar w1206 = w106*w1196 + w1205;
    const scalar w1207 = -w1192*w46 - w1192*w48;
    const scalar w1208 = w1192*w53 - w1192*w706;
    const scalar w1209 = w1207*w50 + w1208*w55;
    const scalar w1210 = w1209 + w413*w443 + w415*w444;
    const scalar w1211 = w4*y;
    const scalar w1212 = w1011*w1211;
    const scalar w1213 = w1014*w1211;
    const scalar w1214 = -gil*w1212 + git*w1213 + w1197*w52 + w1197*w57 - w1202*w418 - w1203*w420 - w413*w995 + w415*w997;
    const scalar w1215 = -w1214;
    const scalar w1216 = w1024 + w106*w1204 + w1215*w95 + w245*w416 + w434*w554 + w437*w551 + w445*w562 + w75*w901 + w845 + w847;
    const scalar w1217 = w119*w1210 + w1216;
    const scalar w1218 = w1192*w37;
    const scalar w1219 = w1192*w38;
    const scalar w1220 = gil*w1219;
    const scalar w1221 = w1192*w36;
    const scalar w1222 = git*w1221;
    const scalar w1223 = w1192*w39;
    const scalar w1224 = w1218 - w1220 - w1222 + w1223;
    const scalar w1225 = -w1224;
    const scalar w1226 = w106*w1201 + w119*w1215 + w1225*w95 + w245*w559 - w247*w431 + w416*w434 + w421*w434 + w437*w565 + w440*w562 + w58*w901 + w68*w807;
    const scalar w1227 = w1226*w236 + w454 + w577;
    const scalar w1228 = w1201*w236 + w1204*w557 + w476 + w600;
    const scalar w1229 = w1204*w427 + w1215*w236 + w479 + w597;
    const scalar w1230 = w236*(-w1201*w427 - w1215*w557 - w1225*w236 - w472 - w594) + w480 + w601;
    const scalar w1231 = w1226*w233 + w448 + w573;
    const scalar w1232 = w1201*w233 + w1204*w555 + w462 + w587;
    const scalar w1233 = w1204*w425 + w1215*w233 + w465 + w584;
    const scalar w1234 = w233*(w1201*w425 + w1215*w555 + w1225*w233 + w458 + w581) + w466 + w588;
    const scalar w1235 = w1034*w578 + w1035*w951 - w1036*w578 + w1037*w590 + w1038*w574 + w1040*w590 - w1041*w951 + w452*w556 - w456*w558 + w470*w556 - w558*w858 + w574*w855;
    const scalar w1236 = w485*w605;
    const scalar w1237 = w1236*w754;
    const scalar w1238 = w1193*w168 + w1194*w170;
    const scalar w1239 = w1196 + w1238 + w422*w609 + w423*w610;
    const scalar w1240 = -w1109*w869 + w1110*w870 + w1197*w169 + w1197*w171 - w1198*w491 - w1199*w492 + w1200 - w609*w995 + w610*w997;
    const scalar w1241 = w1198*w161 + w1199*w163 + w1202*w168 + w1203*w170 + w1204 + w417*w613 + w419*w614 + w438*w609 + w439*w610;
    const scalar w1242 = -w1241;
    const scalar w1243 = w118*w1242 + w1240*w94 + w247*w494 + w433*w622 - w503*w560;
    const scalar w1244 = w1243 - w225*w806 + w433*w518 - w436*w636;
    const scalar w1245 = w1207*w161 + w1208*w163;
    const scalar w1246 = w1210 + w1245 + w443*w613 + w444*w614;
    const scalar w1247 = -gel*w1212 + get*w1213 + w1197*w162 + w1197*w164 - w1202*w491 - w1203*w492 + w1214 - w613*w995 + w614*w997;
    const scalar w1248 = w1059 + w1242*w99 + w1247*w94 + w433*w503 - w436*w622 + w880;
    const scalar w1249 = w1248 + w433*w634 - w516*w560;
    const scalar w1250 = w1192*w157;
    const scalar w1251 = w1192*w158;
    const scalar w1252 = gel*w1219;
    const scalar w1253 = get*w1221;
    const scalar w1254 = -w1218 + w1220 + w1222 - w1223 - w1250 - w1251 + w1252 + w1253;
    const scalar w1255 = w118*w1247 + w1240*w99 - w196*w806 + w433*w494 + w433*w616 - w436*w612 - w496*w560;
    const scalar w1256 = w236*(w1255 + w247*w513 + w94*(-w1237 - w1254)) + w517 + w637;
    const scalar w1257 = -w1240;
    const scalar w1258 = w1241*w557 + w1257*w236 + w540 + w660;
    const scalar w1259 = -w1247;
    const scalar w1260 = w1241*w427 + w1259*w236 + w543 + w657;
    const scalar w1261 = -w1224 - w1250 - w1251 + w1252 + w1253;
    const scalar w1262 = w236*(-w1257*w427 - w1259*w557 - w1261*w236 - w536 - w654) + w544 + w661;
    const scalar w1263 = w1236*w783;
    const scalar w1264 = w1243 - w204*w806 + w433*w507 - w436*w626;
    const scalar w1265 = w1248 + w433*w620 - w501*w560;
    const scalar w1266 = w233*(w1255 + w247*w490 + w94*(-w1254 - w1263)) + w505 + w628;
    const scalar w1267 = w1241*w555 + w1257*w233 + w527 + w647;
    const scalar w1268 = w1241*w425 + w1259*w233 + w530 + w644;
    const scalar w1269 = w233*(w1257*w425 + w1259*w555 + w1261*w233 + w523 + w641) + w531 + w648;
    const scalar w1270 = w1075*w638 + w1076*w663 - w1077*w558 + w1079*w629 + w1080*w650 - w1081*w663 + w512*w556 + w534*w556 - w558*w891 + w629*w888 - w638*w887 + w650*w892;
    const scalar w1271 = w1083*w664;
    const scalar w1272 = w1180*w547;
    const scalar w1273 = w1087*y;
    const scalar w1274 = (3.0/5.0)*St*w110*w25 + w1083*w982 + w1083*w984 + w1180*w896 + w1180*w898 + w1271*w155 + w1271*w32 + w1272*w155 + w1272*w32 + w1273*w150 + w1273*w229;
    const scalar w1275 = w417*w569;
    const scalar w1276 = w419*w570;
    const scalar w1277 = gil*w1275 + git*w1276 + w1195;
    const scalar w1278 = w106*w1277 + w1205;
    const scalar w1279 = w438*w552;
    const scalar w1280 = w439*w553;
    const scalar w1281 = gil*w1279 + git*w1280 + w1209;
    const scalar w1282 = w119*w1281 + w1216;
    const scalar w1283 = gel*w1275 + get*w1276 + w1238 + w1277;
    const scalar w1284 = gel*w1279 + get*w1280 + w1245 + w1281;
    return w23*w411 + w263*w548 + w273*w665 + w40*(w150*w667 + w155*w801 + w229*w667 + w24 - w241 + w32*w801 + w330*(w293*w310*w741 + w294*w739 + w297*w325*w735 - 2*w297*w738 + w298*(w236*w712 + w274*w668 + w427*w733 + w557*w725) + w313*w739 + w328*(w236*(-w23*w668 - w236*w690 - w427*w703 - w557*w711) + w316*w668 + w427*(-w236*w703 - 2*w315 - w427*w732 - w557*w724) + w557*(-w236*w711 - 2*w314 - w427*w724 - w557*w721)) + w737*(-w233*w712 - w274*w736 - w425*w733 - w555*w725) + w742*w744 + w745*(-w233*(w23*w736 + w233*w690 + w425*w703 + w555*w711) - w301*w736 - w425*(w233*w703 + 2*w300 + w425*w732 + w555*w724) - w555*(w233*w711 + 2*w299 + w425*w724 + w555*w721)) - 2*w746*w747 + w748*w750 + w753) + w330*(w366*w391*w788 + w367*w739 + w377*w406*w778 - 2*w377*w787 + w378*(-w236*(w221*w671 + w370*w764 + w771 + w94*(-w757 - w763)) - w371*w668 - w427*(w777 + w99*(-w757 - w776)) - w557*(w118*(-w757 - w772) - w223*w675 + w373*w713 + w775)) + w394*w739 - 2*w406*w794 + w409*(w236*(-w236*w780 - w379*w668 - w427*w782 - w557*w781) + w397*w668 + w427*(-w236*w782 - 2*w396 - w427*w776 - w557*w773) + w557*(-w236*w781 - 2*w395 - w427*w773 - w557*w772)) + w786*(w233*(w195*w671 + w341*w764 + w771 + w94*(-w763 - w785)) + w354*w736 + w425*(w777 + w99*(-w776 - w785)) + w555*(w118*(-w772 - w785) - w200*w675 + w357*w713 + w775)) + w789*w791 + w792*(-w233*(w233*w780 + w379*w736 + w425*w782 + w555*w781) - w382*w736 - w425*(w233*w782 + 2*w381 + w425*w776 + w555*w773) - w555*(w233*w781 + 2*w380 + w425*w773 + w555*w772)) + w795*w797 + w798) + w799*w800 + w799*w802 + w803) + w411*w416 + w411*w421 + w424*w548 + w446*w665 + w548*w551 + w554*w665 + w58*(w1184 + w330*(w1144 + w298*(w1136 + w1186*w427 + w1188*w557) + w328*(w1139 + w427*(-w1137 - w1185*w427) + w557*(-w1138 - w1187*w557)) + w737*(-w1140 - w1186*w425 - w1188*w555) + w744*w949 + w745*(-w1143 - w425*(w1141 + w1185*w425) - w555*(w1142 + w1187*w555)) + w750*w953) + w330*(w1179 + w378*(-w1165 - w427*(w1153 + w99*(-w1146 - w1189)) - w557*(w1158 + w118*(-w1146 - w1190))) + w409*(w1171 + w427*(-w1167 - w1189*w427) + w557*(-w1169 - w1190*w557)) + w786*(w1175 + w425*(w1173 + w99*(-w1172 - w1189)) + w555*(w1174 + w118*(-w1172 - w1190))) + w791*w976 + w792*(-w1178 - w425*(w1176 + w1189*w425) - w555*(w1177 + w1190*w555)) + w797*w979)) + w58*(w1184 + w330*(w1144 + w298*(w1117*w427 + w1125*w557 + w1136) + w328*(w1139 + w427*(-w1105*w427 - w1137) + w557*(-w1121*w557 - w1138)) + w737*(-w1117*w425 - w1125*w555 - w1140) + w742*w950 + w745*(-w1143 - w425*(w1105*w425 + w1141) - w555*(w1121*w555 + w1142)) + w748*w954) + w330*(w1179 + w378*(-w1165 - w427*(w1153 + w99*(-w1146 - w1148)) - w557*(w1158 + w118*(-w1146 - w1155))) + w409*(w1171 + w427*(-w1148*w427 - w1167) + w557*(-w1155*w557 - w1169)) + w786*(w1175 + w425*(w1173 + w99*(-w1148 - w1172)) + w555*(w1174 + w118*(-w1155 - w1172))) + w789*w977 + w792*(-w1178 - w425*(w1148*w425 + w1176) - w555*(w1155*w555 + w1177)) + w795*w980)) + w68*(w1089 + w330*(w1042 + w298*(w1026 + w1093*w427 + w1097*w557) + w328*(w1029 + w427*(-w1027 - w1092*w427) + w557*(-w1028 - w1096*w557)) + w737*(-w1030 - w1093*w425 - w1097*w555) + w744*w856 + w745*(-w1033 - w425*(w1031 + w1092*w425) - w555*(w1032 + w1096*w555)) + w750*w859) + w330*(w1082 + w378*(-w1061 - w427*(w1051 + w99*(-w1043 - w1098)) - w557*(w1056 + w118*(-w1043 - w1099))) + w409*(w1067 + w427*(-w1063 - w1098*w427) + w557*(-w1065 - w1099*w557)) + w786*(w1071 + w425*(w1069 + w99*(-w1068 - w1098)) + w555*(w1070 + w118*(-w1068 - w1099))) + w791*w889 + w792*(-w1074 - w425*(w1072 + w1098*w425) - w555*(w1073 + w1099*w555)) + w797*w893)) + w68*(w1089 + w330*(w1042 + w298*(w1006*w427 + w1019*w557 + w1026) + w328*(w1029 + w427*(-w1027 - w427*w991) + w557*(-w1010*w557 - w1028)) + w737*(-w1006*w425 - w1019*w555 - w1030) + w742*w857 + w745*(-w1033 - w425*(w1031 + w425*w991) - w555*(w1010*w555 + w1032)) + w748*w860) + w330*(w1082 + w378*(-w1061 - w427*(w1051 + w99*(-w1043 - w1045)) - w557*(w1056 + w118*(-w1043 - w1053))) + w409*(w1067 + w427*(-w1045*w427 - w1063) + w557*(-w1053*w557 - w1065)) + w786*(w1071 + w425*(w1069 + w99*(-w1045 - w1068)) + w555*(w1070 + w118*(-w1053 - w1068))) + w789*w890 + w792*(-w1074 - w425*(w1045*w425 + w1072) - w555*(w1053*w555 + w1073)) + w795*w894)) + w75*(w150*w899 + w155*w983 + w229*w899 + w32*w983 + w330*(w298*(w236*w927 + w427*w946 + w557*w940 + w568*w668) + w328*(w236*(-w236*w913 - w427*w922 - w557*w926 - 2*w593) + w427*(-w236*w922 - w427*w945 - w557*w939 - 2*w597) + w557*(-w236*w926 - w427*w939 - w554*w668 - w557*w936) + w598*w668) + w574*w590*w741 + w575*w948 + w578*w603*w735 - 2*w578*w947 + w591*w948 + w737*(-w233*w927 - w425*w946 - w555*w940 - w568*w736) + w745*(-w233*(w233*w913 + w425*w922 + w555*w926 + 2*w580) - w425*(w233*w922 + w425*w945 + w555*w939 + 2*w584) - w555*(w233*w926 + w425*w939 + w554*w736 + w555*w936) - w585*w736) + w753 + w949*w950 - 2*w951*w952 + w953*w954) + w330*(w378*(-w236*(-w221*w900 + w632*w713 + w94*(-w863 - w959) + w963) - w427*(w432*w636 + w971 + w99*(-w862 - w970)) - w557*(w118*(-w862 - w964) + w223*w902 - w634*w962 + w969) - w635*w668) + w409*(w236*(-w236*w972 - w427*w974 - w557*w973 - 2*w653) + w427*(-w236*w974 - w427*w970 - w557*w965 - 2*w657) + w557*(-w236*w973 - w427*w965 - w557*w964 - w619*w668) + w658*w668) + w629*w650*w788 + w630*w948 + w638*w663*w778 - 2*w638*w975 + w651*w948 - 2*w663*w978 + w786*(w233*(-w195*w900 + w608*w713 + w94*(-w886 - w959) + w963) + w425*(w432*w626 + w971 + w99*(-w885 - w970)) + w555*(w118*(-w885 - w964) + w200*w902 - w620*w962 + w969) + w624*w736) + w792*(-w233*(w233*w972 + w425*w974 + w555*w973 + 2*w640) - w425*(w233*w974 + w425*w970 + w555*w965 + 2*w644) - w555*(w233*w973 + w425*w965 + w555*w964 + w619*w736) - w645*w736) + w798 + w976*w977 + w979*w980) + w561 + w803 + w981*w982 + w981*w984) + w78*(w1274 + w330*(w1235 + w298*(w1227 + w1278*w427 + w1282*w557) + w328*(w1230 + w427*(-w1228 - w1277*w427) + w557*(-w1229 - w1281*w557)) + w737*(-w1231 - w1278*w425 - w1282*w555) + w745*(-w1234 - w425*(w1232 + w1277*w425) - w555*(w1233 + w1281*w555)) + w856*w950 + w859*w954) + w330*(w1270 + w378*(-w1256 - w427*(w1244 + w99*(-w1237 - w1283)) - w557*(w118*(-w1237 - w1284) + w1249)) + w409*(w1262 + w427*(-w1258 - w1283*w427) + w557*(-w1260 - w1284*w557)) + w786*(w1266 + w425*(w1264 + w99*(-w1263 - w1283)) + w555*(w118*(-w1263 - w1284) + w1265)) + w792*(-w1269 - w425*(w1267 + w1283*w425) - w555*(w1268 + w1284*w555)) + w889*w977 + w893*w980)) + w78*(w1274 + w330*(w1235 + w298*(w1206*w427 + w1217*w557 + w1227) + w328*(w1230 + w427*(-w1196*w427 - w1228) + w557*(-w1210*w557 - w1229)) + w737*(-w1206*w425 - w1217*w555 - w1231) + w745*(-w1234 - w425*(w1196*w425 + w1232) - w555*(w1210*w555 + w1233)) + w857*w949 + w860*w953) + w330*(w1270 + w378*(-w1256 - w427*(w1244 + w99*(-w1237 - w1239)) - w557*(w118*(-w1237 - w1246) + w1249)) + w409*(w1262 + w427*(-w1239*w427 - w1258) + w557*(-w1246*w557 - w1260)) + w786*(w1266 + w425*(w1264 + w99*(-w1239 - w1263)) + w555*(w118*(-w1246 - w1263) + w1265)) + w792*(-w1269 - w425*(w1239*w425 + w1267) - w555*(w1246*w555 + w1268)) + w890*w976 + w894*w979)) + w84*(w150*w804 + w155*w897 + w229*w804 - w238*w24*w97 + w32*w897 + w330*(w298*(w236*w835 + w427*w853 + w450*w668 + w557*w844) + w328*(w236*(-w236*w817 - w427*w829 - 2*w473 - w557*w834) + w427*(-w236*w829 - w424*w668 - w427*w852 - w557*w843) + w477*w668 + w557*(-w236*w834 - w427*w843 - 2*w476 - w557*w840)) + 2*w426*w470 + w451*w469*w740 + 2*w451*w855 + w455*w483*w734 - w456*w854 + w737*(-w233*w835 - w425*w853 - w450*w736 - w555*w844) + w745*(-w233*(w233*w817 + w425*w829 + 2*w459 + w555*w834) - w425*(w233*w829 + w424*w736 + w425*w852 + w555*w843) - w463*w736 - w555*(w233*w834 + w425*w843 + 2*w462 + w555*w840)) + w753 - w854*w858 + w856*w857 + w859*w860) + w330*(w378*(-w236*(w874 + w94*(-w863 - w868)) - w427*(-w225*w805 - w518*w873 + w881 + w99*(-w862 - w879)) - w519*w668 - w557*(w118*(-w862 - w875) - w223*w806 + w432*w516 + w878)) + w409*(w236*(-w236*w882 - w427*w884 - 2*w537 - w557*w883) + w427*(-w236*w884 - w427*w879 - w506*w668 - w557*w876) + w541*w668 + w557*(-w236*w883 - w427*w876 - 2*w540 - w557*w875)) + w511*w533*w788 + 2*w511*w888 + w520*w546*w778 - 2*w520*w887 + 2*w533*w892 + w786*(w233*(w874 + w94*(-w868 - w886)) + w425*(-w204*w805 - w507*w873 + w881 + w99*(-w879 - w885)) + w510*w736 + w555*(w118*(-w875 - w885) - w200*w806 + w432*w501 + w878)) + w792*(-w233*(w233*w882 + w425*w884 + 2*w524 + w555*w883) - w425*(w233*w884 + w425*w879 + w506*w736 + w555*w876) - w528*w736 - w555*(w233*w883 + w425*w876 + 2*w527 + w555*w875)) + w798 - w854*w891 + w889*w890 + w893*w894) + w435 + w803 + w895*w896 + w895*w898);
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
        const scalar w7 = cos(w3);
        const scalar w8 = k*z;
        const scalar w9 = sin(w8);
        const scalar w10 = (1.0/2.0)*w9;
        const scalar w11 = w10 + 1;
        const scalar w12 = w0 + 1.0/5.0;
        const scalar w13 = sin(w12);
        const scalar w14 = sigmaEps*w11*w13;
        const scalar w15 = w14*w7;
        const scalar w16 = w15 + 1;
        const scalar w17 = sigma0In*w16;
        const scalar w18 = w17*w6;
        const scalar w19 = w6*y;
        const scalar w20 = (z*z);
        const scalar w21 = (1.0/5.0)*w6;
        const scalar w22 = (1.0/5.0)*k*w1*w6*x - 1.0/10.0*w19 - w20*w21;
        const scalar w23 = sigma0In*w22;
        const scalar w24 = k*w23;
        const scalar w25 = sin(w3);
        const scalar w26 = w14*w25;
        const scalar w27 = sin(gamma2);
        const scalar w28 = (w27*w27);
        const scalar w29 = cos(gamma1);
        const scalar w30 = (1.0/10.0)*w28*(w29*w29);
        const scalar w31 = sin(gamma1);
        const scalar w32 = cos(theta0);
        const scalar w33 = w31*w32;
        const scalar w34 = cos(gamma2);
        const scalar w35 = sin(theta0);
        const scalar w36 = w29*w35;
        const scalar w37 = w33 + w34*w36;
        const scalar w38 = (1.0/10.0)*(w37*w37);
        const scalar w39 = w31*w35;
        const scalar w40 = w29*w32*w34 - w39;
        const scalar w41 = (1.0/10.0)*(w40*w40);
        const scalar w42 = w27*w34;
        const scalar w43 = gen*w29*w42;
        const scalar w44 = gin*w42;
        const scalar w45 = w29*w44;
        const scalar w46 = w27*w35;
        const scalar w47 = w37*w46;
        const scalar w48 = gel*w47;
        const scalar w49 = gil*w47;
        const scalar w50 = w27*w32;
        const scalar w51 = w40*w50;
        const scalar w52 = get*w51;
        const scalar w53 = git*w51;
        const scalar w54 = (1.0/5.0)*z;
        const scalar w55 = gen*w31;
        const scalar w56 = w28*w29;
        const scalar w57 = w55*w56;
        const scalar w58 = gin*w31*w56;
        const scalar w59 = w29*w32 - w34*w39;
        const scalar w60 = w37*w59;
        const scalar w61 = gel*w60;
        const scalar w62 = gil*w60;
        const scalar w63 = -w33*w34 - w36;
        const scalar w64 = w40*w63;
        const scalar w65 = get*w64;
        const scalar w66 = git*w64;
        const scalar w67 = -w57 - w58 + w61 + w62 + w65 + w66;
        const scalar w68 = k*w1;
        const scalar w69 = (1.0/5.0)*w68;
        const scalar w70 = w6*(gel*w38 + gen*w30 + get*w41 + gil*w38 + gin*w30 + git*w41 + w54*(-w43 - w45 + w48 + w49 + w52 + w53) - w67*w69);
        const scalar w71 = (1.0/10.0)*w18 + w24*w26 - w70;
        const scalar w72 = (1.0/(sigma0In));
        const scalar w73 = (1.0/(w16));
        const scalar w74 = w72*w73;
        const scalar w75 = w6*w74;
        const scalar w76 = (r1*r1) - w20 - (x*x) - (y*y);
        const scalar w77 = (1.0/2.0)*w76;
        const scalar w78 = w75*w77;
        const scalar w79 = (w34*w34);
        const scalar w80 = gen*w79;
        const scalar w81 = gin*w79;
        const scalar w82 = w28*(w35*w35);
        const scalar w83 = gel*w82;
        const scalar w84 = w28*(w32*w32);
        const scalar w85 = get*w84;
        const scalar w86 = gil*w82;
        const scalar w87 = git*w84;
        const scalar w88 = w42*w55;
        const scalar w89 = w31*w44;
        const scalar w90 = w46*w59;
        const scalar w91 = gel*w90;
        const scalar w92 = gil*w90;
        const scalar w93 = w50*w63;
        const scalar w94 = get*w93;
        const scalar w95 = git*w93;
        const scalar w96 = w88 + w89 + w91 + w92 + w94 + w95;
        const scalar w97 = w6*(-1.0/10.0*w43 - 1.0/10.0*w45 + (1.0/10.0)*w48 + (1.0/10.0)*w49 + (1.0/10.0)*w52 + (1.0/10.0)*w53 + w54*(w80 + w81 + w83 + w85 + w86 + w87) - w69*w96);
        const scalar w98 = w54*w96;
        const scalar w99 = w28*(w31*w31);
        const scalar w100 = (w59*w59);
        const scalar w101 = (w63*w63);
        const scalar w102 = gel*w100 + gen*w99 + get*w101 + gil*w100 + gin*w99 + git*w101;
        const scalar w103 = w6*(-w102*w69 - 1.0/10.0*w57 - 1.0/10.0*w58 + (1.0/10.0)*w61 + (1.0/10.0)*w62 + (1.0/10.0)*w65 + (1.0/10.0)*w66 + w98);
        const scalar w104 = -w103*x - w17*w22 - w70*y - w97*z;
        const scalar w105 = w104*w74;
        const scalar w106 = (1.0/(w16*w16));
        const scalar w107 = w106*w72;
        const scalar w108 = k*w107;
        const scalar w109 = w108*w77;
        const scalar w110 = w104*w109*w6;
        const scalar w111 = k*sigma0In;
        const scalar w112 = w111*w26;
        const scalar w113 = w5*cos(w4);
        const scalar w114 = sigmaEps*w7;
        const scalar w115 = w114*w24;
        const scalar w116 = cos(w8);
        const scalar w117 = (1.0/2.0)*w116;
        const scalar w118 = w6*((1.0/5.0)*w80 + (1.0/5.0)*w81 + (1.0/5.0)*w83 + (1.0/5.0)*w85 + (1.0/5.0)*w86 + (1.0/5.0)*w87);
        const scalar w119 = w6*x;
        const scalar w120 = (2.0/5.0)*sigma0In*w16*w6*z - w115*w117*w13 - w118*z - w119*((1.0/5.0)*w88 + (1.0/5.0)*w89 + (1.0/5.0)*w91 + (1.0/5.0)*w92 + (1.0/5.0)*w94 + (1.0/5.0)*w95) - w19*(-1.0/5.0*w43 - 1.0/5.0*w45 + (1.0/5.0)*w48 + (1.0/5.0)*w49 + (1.0/5.0)*w52 + (1.0/5.0)*w53) - w97;
        const scalar w121 = w105*w6;
        const scalar w122 = w114*w13;
        const scalar w123 = w116*w122*w6;
        const scalar w124 = cos(w0);
        const scalar w125 = (1.0/2.0)*y + 1;
        const scalar w126 = (k*k);
        const scalar w127 = w124*w126;
        const scalar w128 = w127*w21;
        const scalar w129 = w128*x;
        const scalar w130 = w129 + w21*w68;
        const scalar w131 = cos(w12);
        const scalar w132 = w11*w131;
        const scalar w133 = w67*y;
        const scalar w134 = w102*w129 - w103 - w115*w132 + w127*w6*w98 + w128*w133 - w130*w17;
        const scalar w135 = w114*w132;
        const scalar w136 = w111*w135;
        const scalar w137 = St*w113;
        const scalar w138 = w126*w15;
        const scalar w139 = 2*w74;
        const scalar w140 = w6*w76;
        const scalar w141 = w108*w140;
        const scalar w142 = -w121;
        const scalar w143 = w104*w107;
        const scalar w144 = w143*w6;
        const scalar w145 = (w13*w13);
        const scalar w146 = (sigmaEps*sigmaEps);
        const scalar w147 = (w11*w11);
        const scalar w148 = (1.0/(w16*w16*w16));
        const scalar w149 = w104*w140;
        const scalar w150 = w123*w8;
        const scalar w151 = w122*w126;
        const scalar w152 = (1.0/4.0)*w149;
        const scalar w153 = (w7*w7);
        const scalar w154 = (k*k*k);
        const scalar w155 = w1*w154*w21;
        const scalar w156 = w155*x;
        return w111*w117*w122*(St*(fA1*y + (2.0/5.0)*w113) + Tt*(-1.0/4.0*w104*w108*w123*w76 + (1.0/2.0)*w120*w6*w72*w73*w76 - w121*z + (1.0/5.0)*z)) - w112*(St*(fA1*z - k*w5*sin(w4) + (1.0/2.0)*w2) + Tt*(-w105*w19 + w110*w26 + w71*w78 + 1.0/10.0)) + w136*(St*(fA0 + fA2*k*w124*w125) + Tt*(-w105*w119 - w110*w135 + (1.0/2.0)*w134*w6*w72*w73*w76 - w69)) + w17*(Tt*(w126*w145*w146*w147*w148*w149*(w25*w25)*w72 + w138*w144*w77 - w139*w19*w71 + w141*w26*w71 + w142 - 2*w144*w26*w3 + w78*(-w112*w21 + w138*w23)) - w126*w137) + w17*(Tt*(w107*w151*w152*w9 - w109*w120*w123 + (w116*w116)*w126*w145*w146*w148*w152*w153*w72 - 2*w120*w75*z + w142 + w143*w150 + w78*((2.0/5.0)*sigma0In*w150 + w10*w151*w23 - 2*w118 + (2.0/5.0)*w18) + 1.0/5.0) + (4.0/25.0)*w137) + w17*(-St*w125*w126*w2 + Tt*(2*k*sigmaEps*w104*w106*w11*w131*w6*w7*w72*x + (1.0/2.0)*sigmaEps*w104*w106*w11*w126*w13*w6*w7*w72*w76 + w104*w126*(w131*w131)*w146*w147*w148*w153*w6*w72*w76 - w119*w134*w139 - w121 - 1.0/5.0*w127 - w134*w135*w141 + (1.0/2.0)*w6*w72*w73*w76*(sigma0In*sigmaEps*w11*w126*w13*w22*w7 - w1*w154*w6*w98 + (2.0/5.0)*w102*w124*w126*w6 - w102*w156 - 2*w130*w136 - w133*w155 - w17*((2.0/5.0)*w124*w126*w6 - w156))));
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
        const scalar w8 = cos(w3);
        const scalar w9 = k*z;
        const scalar w10 = sin(w9);
        const scalar w11 = (1.0/2.0)*w10;
        const scalar w12 = w11 + 1;
        const scalar w13 = w0 + 1.0/5.0;
        const scalar w14 = sin(w13);
        const scalar w15 = sigmaEps*w12*w14;
        const scalar w16 = w15*w8;
        const scalar w17 = w16 + 1;
        const scalar w18 = sigma0Out*w17;
        const scalar w19 = w18*w7;
        const scalar w20 = w7*y;
        const scalar w21 = (z*z);
        const scalar w22 = (1.0/5.0)*w7;
        const scalar w23 = w1*w22;
        const scalar w24 = -w0*w23 + (1.0/10.0)*w20 + w21*w22;
        const scalar w25 = sigma0Out*w24;
        const scalar w26 = k*w25;
        const scalar w27 = sin(w3);
        const scalar w28 = w15*w27;
        const scalar w29 = sin(gamma2);
        const scalar w30 = (w29*w29);
        const scalar w31 = cos(gamma1);
        const scalar w32 = (1.0/10.0)*w30*(w31*w31);
        const scalar w33 = sin(gamma1);
        const scalar w34 = theta0 + theta1;
        const scalar w35 = cos(w34);
        const scalar w36 = w33*w35;
        const scalar w37 = cos(gamma2);
        const scalar w38 = sin(w34);
        const scalar w39 = w31*w38;
        const scalar w40 = w36 + w37*w39;
        const scalar w41 = (1.0/10.0)*(w40*w40);
        const scalar w42 = w33*w38;
        const scalar w43 = w31*w35*w37 - w42;
        const scalar w44 = (1.0/10.0)*(w43*w43);
        const scalar w45 = w29*w37;
        const scalar w46 = gen*w31*w45;
        const scalar w47 = gin*w45;
        const scalar w48 = w31*w47;
        const scalar w49 = w29*w38;
        const scalar w50 = w40*w49;
        const scalar w51 = gel*w50;
        const scalar w52 = gil*w50;
        const scalar w53 = w29*w35;
        const scalar w54 = w43*w53;
        const scalar w55 = get*w54;
        const scalar w56 = git*w54;
        const scalar w57 = (1.0/5.0)*z;
        const scalar w58 = gen*w33;
        const scalar w59 = w30*w31;
        const scalar w60 = w58*w59;
        const scalar w61 = gin*w33*w59;
        const scalar w62 = w31*w35 - w37*w42;
        const scalar w63 = w40*w62;
        const scalar w64 = gel*w63;
        const scalar w65 = gil*w63;
        const scalar w66 = -w36*w37 - w39;
        const scalar w67 = w43*w66;
        const scalar w68 = get*w67;
        const scalar w69 = git*w67;
        const scalar w70 = -w60 - w61 + w64 + w65 + w68 + w69;
        const scalar w71 = k*w1;
        const scalar w72 = (1.0/5.0)*w71;
        const scalar w73 = w7*(gel*w41 + gen*w32 + get*w44 + gil*w41 + gin*w32 + git*w44 + w57*(-w46 - w48 + w51 + w52 + w55 + w56) - w70*w72);
        const scalar w74 = -1.0/10.0*w19 + w26*w28 + w73;
        const scalar w75 = (1.0/(sigma0Out));
        const scalar w76 = (1.0/(w17));
        const scalar w77 = w75*w76;
        const scalar w78 = w7*w77;
        const scalar w79 = -(r2*r2) + w21 + (x*x) + (y*y);
        const scalar w80 = (1.0/2.0)*w79;
        const scalar w81 = w78*w80;
        const scalar w82 = (w37*w37);
        const scalar w83 = gen*w82;
        const scalar w84 = gin*w82;
        const scalar w85 = w30*(w38*w38);
        const scalar w86 = gel*w85;
        const scalar w87 = w30*(w35*w35);
        const scalar w88 = get*w87;
        const scalar w89 = gil*w85;
        const scalar w90 = git*w87;
        const scalar w91 = w45*w58;
        const scalar w92 = w33*w47;
        const scalar w93 = w49*w62;
        const scalar w94 = gel*w93;
        const scalar w95 = gil*w93;
        const scalar w96 = w53*w66;
        const scalar w97 = get*w96;
        const scalar w98 = git*w96;
        const scalar w99 = w91 + w92 + w94 + w95 + w97 + w98;
        const scalar w100 = w7*(-1.0/10.0*w46 - 1.0/10.0*w48 + (1.0/10.0)*w51 + (1.0/10.0)*w52 + (1.0/10.0)*w55 + (1.0/10.0)*w56 + w57*(w83 + w84 + w86 + w88 + w89 + w90) - w72*w99);
        const scalar w101 = w30*(w33*w33);
        const scalar w102 = (w62*w62);
        const scalar w103 = (w66*w66);
        const scalar w104 = gel*w102 + gen*w101 + get*w103 + gil*w102 + gin*w101 + git*w103;
        const scalar w105 = w57*w99;
        const scalar w106 = -w104*w72 + w105 - 1.0/10.0*w60 - 1.0/10.0*w61 + (1.0/10.0)*w64 + (1.0/10.0)*w65 + (1.0/10.0)*w68 + (1.0/10.0)*w69;
        const scalar w107 = w100*z + w106*w7*x - w18*w24 + w73*y;
        const scalar w108 = w107*w77;
        const scalar w109 = w75/(w17*w17);
        const scalar w110 = k*w109;
        const scalar w111 = w110*w80;
        const scalar w112 = w107*w111*w7;
        const scalar w113 = k*sigma0Out;
        const scalar w114 = w113*w28;
        const scalar w115 = w6*cos(w4);
        const scalar w116 = sigmaEps*w8;
        const scalar w117 = w116*w26;
        const scalar w118 = cos(w9);
        const scalar w119 = (1.0/2.0)*w118;
        const scalar w120 = w7*((1.0/5.0)*w83 + (1.0/5.0)*w84 + (1.0/5.0)*w86 + (1.0/5.0)*w88 + (1.0/5.0)*w89 + (1.0/5.0)*w90);
        const scalar w121 = w7*x;
        const scalar w122 = w100 - w117*w119*w14 + w120*z + w121*((1.0/5.0)*w91 + (1.0/5.0)*w92 + (1.0/5.0)*w94 + (1.0/5.0)*w95 + (1.0/5.0)*w97 + (1.0/5.0)*w98) - w19*w5 + w20*(-1.0/5.0*w46 - 1.0/5.0*w48 + (1.0/5.0)*w51 + (1.0/5.0)*w52 + (1.0/5.0)*w55 + (1.0/5.0)*w56);
        const scalar w123 = w108*w7;
        const scalar w124 = w116*w14;
        const scalar w125 = w118*w124*w7;
        const scalar w126 = cos(w0);
        const scalar w127 = (1.0/2.0)*y + 1;
        const scalar w128 = (k*k);
        const scalar w129 = w126*w128;
        const scalar w130 = w129*w22;
        const scalar w131 = w130*x;
        const scalar w132 = -w131 - w22*w71;
        const scalar w133 = cos(w13);
        const scalar w134 = w12*w133;
        const scalar w135 = w129*w7;
        const scalar w136 = w70*y;
        const scalar w137 = -w104*w131 - w105*w135 + w106*w7 - w117*w134 - w130*w136 - w132*w18;
        const scalar w138 = w116*w134;
        const scalar w139 = w113*w138;
        const scalar w140 = St*w115;
        const scalar w141 = w128*w16;
        const scalar w142 = w141*w25;
        const scalar w143 = 2*w77;
        const scalar w144 = w7*w79;
        const scalar w145 = w110*w144;
        const scalar w146 = (w14*w14);
        const scalar w147 = (sigmaEps*sigmaEps);
        const scalar w148 = (1.0/(w17*w17*w17));
        const scalar w149 = w107*w144;
        const scalar w150 = (w12*w12)*w128*w147*w148*w149*w75;
        const scalar w151 = w107*w109;
        const scalar w152 = w151*w7;
        const scalar w153 = 2*w152;
        const scalar w154 = w123 + w141*w152*w80;
        const scalar w155 = w125*w9;
        const scalar w156 = w124*w128;
        const scalar w157 = (1.0/4.0)*w149;
        const scalar w158 = (w8*w8);
        const scalar w159 = (2.0/5.0)*w135;
        const scalar w160 = (k*k*k);
        const scalar w161 = w160*w23;
        const scalar w162 = w161*x;
        return w113*w119*w124*(St*(fA1*y + (2.0/5.0)*w115) + Tt*(-1.0/4.0*w107*w110*w125*w79 + w122*w81 + w123*z + w57)) - w114*(St*(fA1*z - k*w6*sin(w4) + (1.0/2.0)*w2) + Tt*(w108*w20 + w112*w28 + w74*w81 + 1.0/10.0)) + w139*(St*(fA0 + fA2*k*w126*w127) + Tt*(w107*w7*w75*w76*x - w112*w138 + (1.0/2.0)*w137*w7*w75*w76*w79 - w72)) + w18*(Tt*(w143*w20*w74 + w145*w28*w74 + w146*w150*(w27*w27) + w153*w28*w3 + w154 + w81*(w114*w22 + w142)) - w128*w140) + w18*(Tt*(w10*w109*w156*w157 - w111*w122*w125 + (w118*w118)*w128*w146*w147*w148*w157*w158*w75 + 2*w122*w78*z + w123 - w151*w155 + w81*(-2.0/5.0*sigma0Out*w155 + w11*w156*w25 + 2*w120 - 2.0/5.0*w19) + 1.0/5.0) + (4.0/25.0)*w140) + w18*(-St*w127*w128*w2 + Tt*(-w0*w138*w153 + w121*w137*w143 - 1.0/5.0*w129 + (w133*w133)*w150*w158 - w137*w138*w145 + w154 + w81*(w1*w105*w160*w7 - w104*w159 + w104*w162 - 2*w132*w139 + w136*w161 + w142 - w18*(-w159 + w162))));
    }
    const scalar w0 = (r1*r1);
    const scalar w1 = -(r2*r2);
    const scalar w2 = -w0 - w1;
    const scalar w3 = (1.0/(w2));
    const scalar w4 = theta1*w3;
    const scalar w5 = 4*w4;
    const scalar w6 = gel*w5;
    const scalar w7 = (x*x);
    const scalar w8 = (y*y);
    const scalar w9 = (z*z);
    const scalar w10 = w7 + w8 + w9;
    const scalar w11 = w3*(-w0 + w10);
    const scalar w12 = theta0 + theta1*w11;
    const scalar w13 = cos(w12);
    const scalar w14 = sin(w12);
    const scalar w15 = sin(gamma2);
    const scalar w16 = (w15*w15);
    const scalar w17 = w13*w14*w16;
    const scalar w18 = w17*w6;
    const scalar w19 = w18*z;
    const scalar w20 = w17*w5;
    const scalar w21 = get*w20;
    const scalar w22 = w21*z;
    const scalar w23 = w19 - w22;
    const scalar w24 = k*y;
    const scalar w25 = w24 + 3.0/10.0;
    const scalar w26 = cos(w25);
    const scalar w27 = fA3*exp((2.0/5.0)*z);
    const scalar w28 = (2.0/5.0)*w27;
    const scalar w29 = w26*w28;
    const scalar w30 = St*(fA1*y + w29);
    const scalar w31 = (1.0/5.0)*Tt;
    const scalar w32 = (1.0/(r2));
    const scalar w33 = w32*z;
    const scalar w34 = r1 + r2;
    const scalar w35 = w3*w34;
    const scalar w36 = cos(gamma2);
    const scalar w37 = (w36*w36);
    const scalar w38 = (w14*w14)*w16;
    const scalar w39 = gil*w38;
    const scalar w40 = (w13*w13)*w16;
    const scalar w41 = git*w40;
    const scalar w42 = gin*w37 + w39 + w41;
    const scalar w43 = gel*w38;
    const scalar w44 = get*w40;
    const scalar w45 = gen*w37 + w43 + w44;
    const scalar w46 = w42 + w45;
    const scalar w47 = w32*w46;
    const scalar w48 = cos(gamma1);
    const scalar w49 = w15*w36;
    const scalar w50 = w48*w49;
    const scalar w51 = sin(gamma1);
    const scalar w52 = w13*w51;
    const scalar w53 = w14*w48;
    const scalar w54 = w36*w53;
    const scalar w55 = w52 + w54;
    const scalar w56 = gil*w55;
    const scalar w57 = w14*w15;
    const scalar w58 = w56*w57;
    const scalar w59 = w14*w51;
    const scalar w60 = w13*w36*w48 - w59;
    const scalar w61 = git*w60;
    const scalar w62 = w13*w15;
    const scalar w63 = w61*w62;
    const scalar w64 = -gin*w50 + w58 + w63;
    const scalar w65 = gel*w55;
    const scalar w66 = w57*w65;
    const scalar w67 = get*w60;
    const scalar w68 = w62*w67;
    const scalar w69 = -gen*w50 + w66 + w68;
    const scalar w70 = w64 + w69;
    const scalar w71 = w32*w70;
    const scalar w72 = gin*w51;
    const scalar w73 = w36*w59;
    const scalar w74 = w13*w48 - w73;
    const scalar w75 = gil*w74;
    const scalar w76 = w57*w75;
    const scalar w77 = w36*w52;
    const scalar w78 = -w53 - w77;
    const scalar w79 = git*w78;
    const scalar w80 = w62*w79;
    const scalar w81 = w49*w72 + w76 + w80;
    const scalar w82 = gen*w51;
    const scalar w83 = gel*w74;
    const scalar w84 = w57*w83;
    const scalar w85 = get*w78;
    const scalar w86 = w62*w85;
    const scalar w87 = w49*w82 + w84 + w86;
    const scalar w88 = w81 + w87;
    const scalar w89 = w32*w88;
    const scalar w90 = w32*(w47*z + w71*y + w89*x);
    const scalar w91 = w16*(w48*w48);
    const scalar w92 = (w55*w55);
    const scalar w93 = (w60*w60);
    const scalar w94 = gil*w92 + gin*w91 + git*w93;
    const scalar w95 = gel*w92 + gen*w91 + get*w93;
    const scalar w96 = w94 + w95;
    const scalar w97 = w32*w96;
    const scalar w98 = w16*w48;
    const scalar w99 = w55*w75 + w60*w79 - w72*w98;
    const scalar w100 = w55*w83 + w60*w85 - w82*w98;
    const scalar w101 = w100 + w99;
    const scalar w102 = w101*w32;
    const scalar w103 = w32*(w102*x + w71*z + w97*y);
    const scalar w104 = w16*(w51*w51);
    const scalar w105 = (w74*w74);
    const scalar w106 = (w78*w78);
    const scalar w107 = gil*w105 + gin*w104 + git*w106;
    const scalar w108 = gel*w105 + gen*w104 + get*w106;
    const scalar w109 = w107 + w108;
    const scalar w110 = w109*w32;
    const scalar w111 = w32*(w102*y + w110*x + w89*z);
    const scalar w112 = w103*y + w111*x + w90*z;
    const scalar w113 = (1.0/(w112));
    const scalar w114 = k*x;
    const scalar w115 = w114 + 1.0/5.0;
    const scalar w116 = sin(w115);
    const scalar w117 = cos(w24);
    const scalar w118 = k*z;
    const scalar w119 = sin(w118);
    const scalar w120 = (1.0/2.0)*w119;
    const scalar w121 = w120 + 1;
    const scalar w122 = sigmaEps*w121;
    const scalar w123 = w117*w122;
    const scalar w124 = w116*w123;
    const scalar w125 = w124 + 1;
    const scalar w126 = -sigma0Out*w125;
    const scalar w127 = -w126 - w46;
    const scalar w128 = -w88;
    const scalar w129 = (1.0/2.0)*y + 1;
    const scalar w130 = cos(w114);
    const scalar w131 = fA2*k*w130;
    const scalar w132 = St*(fA0 + w129*w131);
    const scalar w133 = -w70;
    const scalar w134 = sin(w114);
    const scalar w135 = fA2*w134;
    const scalar w136 = (1.0/2.0)*w135;
    const scalar w137 = sin(w25);
    const scalar w138 = k*w137;
    const scalar w139 = w138*w27;
    const scalar w140 = St*(fA1*z + w136 - w139);
    const scalar w141 = w128*w132 + w133*w140;
    const scalar w142 = w32*(w127*w30 + w141);
    const scalar w143 = -w126 - w96;
    const scalar w144 = -w101;
    const scalar w145 = w132*w144 + w133*w30;
    const scalar w146 = w32*(w140*w143 + w145);
    const scalar w147 = -w109 - w126;
    const scalar w148 = w128*w30 + w140*w144;
    const scalar w149 = w32*(w132*w147 + w148);
    const scalar w150 = w142*z + w146*y + w149*x;
    const scalar w151 = w113*w150;
    const scalar w152 = w151*w35;
    const scalar w153 = (1.0/(r1));
    const scalar w154 = w153*z;
    const scalar w155 = w153*w46;
    const scalar w156 = w153*w70;
    const scalar w157 = w153*w88;
    const scalar w158 = w153*(-w155*z - w156*y - w157*x);
    const scalar w159 = w153*w96;
    const scalar w160 = w101*w153;
    const scalar w161 = w153*(-w156*z - w159*y - w160*x);
    const scalar w162 = w109*w153;
    const scalar w163 = w153*(-w157*z - w160*y - w162*x);
    const scalar w164 = -w158*z - w161*y - w163*x;
    const scalar w165 = (1.0/(w164));
    const scalar w166 = -sigma0In*w125;
    const scalar w167 = -w166 - w46;
    const scalar w168 = w153*(w141 + w167*w30);
    const scalar w169 = -w166 - w96;
    const scalar w170 = w153*(w140*w169 + w145);
    const scalar w171 = -w109 - w166;
    const scalar w172 = w153*(w132*w171 + w148);
    const scalar w173 = -w168*z - w170*y - w172*x;
    const scalar w174 = w165*w173;
    const scalar w175 = w174*w35;
    const scalar w176 = k*sigma0Out;
    const scalar w177 = w116*w176;
    const scalar w178 = sigmaEps*w117;
    const scalar w179 = cos(w118);
    const scalar w180 = (1.0/2.0)*w179;
    const scalar w181 = w178*w180;
    const scalar w182 = gil*w20;
    const scalar w183 = w182*z;
    const scalar w184 = git*w20;
    const scalar w185 = w184*z;
    const scalar w186 = -w183 + w185 - w19 + w22;
    const scalar w187 = w177*w181 + w186;
    const scalar w188 = w26*w27;
    const scalar w189 = St*w188;
    const scalar w190 = (4.0/25.0)*w189;
    const scalar w191 = 2*w4;
    const scalar w192 = w191*w62;
    const scalar w193 = w192*w75;
    const scalar w194 = w191*w57;
    const scalar w195 = w194*w79;
    const scalar w196 = w191*w53;
    const scalar w197 = w191*w77;
    const scalar w198 = -w196*z - w197*z;
    const scalar w199 = w198*w57;
    const scalar w200 = w13*w48;
    const scalar w201 = w200*w4;
    const scalar w202 = 2*w201;
    const scalar w203 = w191*w73;
    const scalar w204 = -w202*z + w203*z;
    const scalar w205 = w204*w62;
    const scalar w206 = gil*w199 + git*w205 + w193*z - w195*z;
    const scalar w207 = w192*w83;
    const scalar w208 = w194*w85;
    const scalar w209 = gel*w199 + get*w205 + w207*z - w208*z;
    const scalar w210 = w206 + w209;
    const scalar w211 = -w210;
    const scalar w212 = St*(fA1 - w138*w28);
    const scalar w213 = w133*w212;
    const scalar w214 = w192*w56;
    const scalar w215 = w194*w61;
    const scalar w216 = w191*w59;
    const scalar w217 = 2*theta1*w13*w3*w36*w48*z - w216*z;
    const scalar w218 = w217*w57;
    const scalar w219 = w191*w52;
    const scalar w220 = w191*w54;
    const scalar w221 = -w219*z - w220*z;
    const scalar w222 = w221*w62;
    const scalar w223 = gil*w218 + git*w222 + w214*z - w215*z;
    const scalar w224 = w192*w65;
    const scalar w225 = w194*w67;
    const scalar w226 = gel*w218 + get*w222 + w224*z - w225*z;
    const scalar w227 = w223 + w226;
    const scalar w228 = -w227;
    const scalar w229 = w132*w211 + w140*w228 + w213;
    const scalar w230 = w127*w190 + w187*w30 + w229;
    const scalar w231 = -1.0/2.0*k*sigma0Out*sigmaEps*w116*w117*w179;
    const scalar w232 = w5*w59;
    const scalar w233 = 4*theta1*w13*w3*w36*w48*z - w232*z;
    const scalar w234 = w5*w52;
    const scalar w235 = w5*w54;
    const scalar w236 = -w234*z - w235*z;
    const scalar w237 = w233*w56 + w236*w61;
    const scalar w238 = w233*w65 + w236*w67 + w237;
    const scalar w239 = -w231 - w238;
    const scalar w240 = w198*w56 + w204*w61 + w217*w75 + w221*w79;
    const scalar w241 = w198*w65 + w204*w67 + w217*w83 + w221*w85 + w240;
    const scalar w242 = -w241;
    const scalar w243 = w132*w242 + w133*w190 + w228*w30;
    const scalar w244 = w32*(w140*w239 + w143*w212 + w243);
    const scalar w245 = 4*w53;
    const scalar w246 = w245*w4;
    const scalar w247 = w5*w77;
    const scalar w248 = -w246*z - w247*z;
    const scalar w249 = 4*w201;
    const scalar w250 = w5*w73;
    const scalar w251 = -w249*z + w250*z;
    const scalar w252 = w248*w75 + w251*w79;
    const scalar w253 = w248*w83 + w251*w85 + w252;
    const scalar w254 = -w231 - w253;
    const scalar w255 = w128*w190 + w140*w242 + w144*w212 + w211*w30;
    const scalar w256 = w32*(w132*w254 + w255);
    const scalar w257 = w142 + w230*w33 + w244*y + w256*x;
    const scalar w258 = (1.0/2.0)*w34;
    const scalar w259 = w11*w32;
    const scalar w260 = w258*w259;
    const scalar w261 = w113*w260;
    const scalar w262 = k*sigma0In;
    const scalar w263 = w116*w262;
    const scalar w264 = w181*w263 + w186;
    const scalar w265 = w167*w190 + w229 + w264*w30;
    const scalar w266 = -1.0/2.0*k*sigma0In*sigmaEps*w116*w117*w179;
    const scalar w267 = -w238 - w266;
    const scalar w268 = w153*(w140*w267 + w169*w212 + w243);
    const scalar w269 = -w253 - w266;
    const scalar w270 = w153*(w132*w269 + w255);
    const scalar w271 = -w154*w265 - w168 - w268*y - w270*x;
    const scalar w272 = w153*(1 - w11);
    const scalar w273 = w258*w272;
    const scalar w274 = w165*w273;
    const scalar w275 = (1.0/(w112*w112));
    const scalar w276 = w183 - w185;
    const scalar w277 = w23 + w276;
    const scalar w278 = w227*w32;
    const scalar w279 = w210*w32;
    const scalar w280 = w277*w33 + w278*y + w279*x + w47;
    const scalar w281 = w280*w33;
    const scalar w282 = w238*w32;
    const scalar w283 = w241*w32;
    const scalar w284 = w32*(w278*z + w282*y + w283*x + w71);
    const scalar w285 = w284*y;
    const scalar w286 = w253*w32;
    const scalar w287 = w32*(w279*z + w283*y + w286*x + w89);
    const scalar w288 = w287*x;
    const scalar w289 = -w281 - w285 - w288 - w90;
    const scalar w290 = w275*w289;
    const scalar w291 = w150*w260;
    const scalar w292 = w153*w227;
    const scalar w293 = w153*w210;
    const scalar w294 = -w154*w277 - w155 - w292*y - w293*x;
    const scalar w295 = w154*w294;
    const scalar w296 = w153*w238;
    const scalar w297 = w153*w241;
    const scalar w298 = w153*(-w156 - w292*z - w296*y - w297*x);
    const scalar w299 = w298*y;
    const scalar w300 = w153*w253;
    const scalar w301 = w153*(-w157 - w293*z - w297*y - w300*x);
    const scalar w302 = w301*x;
    const scalar w303 = w158 + w295 + w299 + w302;
    const scalar w304 = (1.0/(w164*w164));
    const scalar w305 = w173*w273;
    const scalar w306 = w304*w305;
    const scalar w307 = w152*w33 - w154*w175 + w257*w261 + w271*w274 + w290*w291 + w303*w306;
    const scalar w308 = w1 + w10;
    const scalar w309 = (1.0/(w34));
    const scalar w310 = w11*w309;
    const scalar w311 = w308*w310;
    const scalar w312 = w151*w260 + w174*w273;
    const scalar w313 = 2*w312;
    const scalar w314 = w310*w313;
    const scalar w315 = w3*w309;
    const scalar w316 = w308*w315;
    const scalar w317 = w313*w316;
    const scalar w318 = w307*w311 + w31*z + w314*z + w317*z;
    const scalar w319 = w30 + w318;
    const scalar w320 = 2*theta1*w13*w3*w36*w48*y - w216*y;
    const scalar w321 = w320*w57;
    const scalar w322 = -w219*y - w220*y;
    const scalar w323 = w322*w62;
    const scalar w324 = gel*w321 + get*w323 + w224*y - w225*y;
    const scalar w325 = -w196*x - w197*x;
    const scalar w326 = w325*w57;
    const scalar w327 = -w202*x + w203*x;
    const scalar w328 = w327*w62;
    const scalar w329 = gel*w326 + get*w328 + w207*x - w208*x;
    const scalar w330 = -w246*x - w247*x;
    const scalar w331 = -w249*x + w250*x;
    const scalar w332 = w330*w83 + w331*w85;
    const scalar w333 = w32*x;
    const scalar w334 = w153*x;
    const scalar w335 = cos(w115);
    const scalar w336 = w123*w335;
    const scalar w337 = w18*x;
    const scalar w338 = w182*x;
    const scalar w339 = w21*x;
    const scalar w340 = w184*x;
    const scalar w341 = -w337 - w338 + w339 + w340;
    const scalar w342 = w176*w336 + w341;
    const scalar w343 = gil*w326 + git*w328 + w193*x - w195*x;
    const scalar w344 = w329 + w343;
    const scalar w345 = -w344;
    const scalar w346 = 2*theta1*w13*w3*w36*w48*x - w216*x;
    const scalar w347 = w346*w57;
    const scalar w348 = -w219*x - w220*x;
    const scalar w349 = w348*w62;
    const scalar w350 = gil*w347 + git*w349 + w214*x - w215*x;
    const scalar w351 = gel*w347 - 2*get*theta1*w14*w15*w3*w60*x + get*w349 + w224*x + w350;
    const scalar w352 = -w351;
    const scalar w353 = St*w131;
    const scalar w354 = (1.0/2.0)*w353;
    const scalar w355 = (k*k);
    const scalar w356 = St*w355;
    const scalar w357 = w129*w135*w356;
    const scalar w358 = -w128*w357 + w132*w345 + w133*w354 + w140*w352;
    const scalar w359 = w32*(w30*w342 + w358);
    const scalar w360 = -k*sigma0Out*sigmaEps*w117*w121*w335;
    const scalar w361 = 4*theta1*w13*w3*w36*w48*x - w232*x;
    const scalar w362 = w361*w55;
    const scalar w363 = -w234*x - w235*x;
    const scalar w364 = w363*w60;
    const scalar w365 = gil*w362 + git*w364;
    const scalar w366 = gel*w362 + get*w364 + w365;
    const scalar w367 = -w360 - w366;
    const scalar w368 = w325*w55;
    const scalar w369 = w327*w60;
    const scalar w370 = gil*w368 + git*w369 + w346*w75 + w348*w79;
    const scalar w371 = gel*w368 + get*w369 + w346*w83 + w348*w85;
    const scalar w372 = w370 + w371;
    const scalar w373 = -w372;
    const scalar w374 = w132*w373 - w144*w357 + w30*w352;
    const scalar w375 = w32*(w140*w367 + w143*w354 + w374);
    const scalar w376 = w330*w75 + w331*w79;
    const scalar w377 = w332 + w376;
    const scalar w378 = -w360 - w377;
    const scalar w379 = w144*w354;
    const scalar w380 = w140*w373 + w30*w345 + w379;
    const scalar w381 = w132*w378 - w147*w357 + w380;
    const scalar w382 = w149 + w333*w381 + w359*z + w375*y;
    const scalar w383 = w262*w336 + w341;
    const scalar w384 = w153*(w30*w383 + w358);
    const scalar w385 = -k*sigma0In*sigmaEps*w117*w121*w335;
    const scalar w386 = -w366 - w385;
    const scalar w387 = w153*(w140*w386 + w169*w354 + w374);
    const scalar w388 = -w377 - w385;
    const scalar w389 = w132*w388 - w171*w357 + w380;
    const scalar w390 = -w172 - w334*w389 - w384*z - w387*y;
    const scalar w391 = w165*w390;
    const scalar w392 = w338 - w340;
    const scalar w393 = w337 - w339 + w392;
    const scalar w394 = w32*w393;
    const scalar w395 = w32*w351;
    const scalar w396 = w32*w344;
    const scalar w397 = w32*(w394*z + w395*y + w396*x + w89);
    const scalar w398 = w397*z;
    const scalar w399 = w32*w372;
    const scalar w400 = w110 + w333*w377 + w396*z + w399*y;
    const scalar w401 = w333*w400;
    const scalar w402 = w32*w366;
    const scalar w403 = w32*(w102 + w395*z + w399*x + w402*y);
    const scalar w404 = w403*y;
    const scalar w405 = -w111 - w398 - w401 - w404;
    const scalar w406 = w275*w405;
    const scalar w407 = w153*w393;
    const scalar w408 = w153*w351;
    const scalar w409 = w153*w344;
    const scalar w410 = w153*(-w157 - w407*z - w408*y - w409*x);
    const scalar w411 = w410*z;
    const scalar w412 = w153*w372;
    const scalar w413 = -w162 - w334*w377 - w409*z - w412*y;
    const scalar w414 = w334*w413;
    const scalar w415 = w153*w366;
    const scalar w416 = w153*(-w160 - w408*z - w412*x - w415*y);
    const scalar w417 = w416*y;
    const scalar w418 = w163 + w411 + w414 + w417;
    const scalar w419 = w304*w418;
    const scalar w420 = w152*w333 - w175*w334 + w261*w382 + w273*w391 + w291*w406 + w305*w419;
    const scalar w421 = -k*w134*w31 + w311*w420 + w314*x + w317*x;
    const scalar w422 = w132 + w421;
    const scalar w423 = -w196*y - w197*y;
    const scalar w424 = -w202*y + w203*y;
    const scalar w425 = w320*w83 + w322*w85 + w423*w65 + w424*w67;
    const scalar w426 = 4*theta1*w13*w3*w36*w48*y - w232*y;
    const scalar w427 = -w234*y - w235*y;
    const scalar w428 = w426*w65 + w427*w67;
    const scalar w429 = w32*y;
    const scalar w430 = w153*y;
    const scalar w431 = sin(w24);
    const scalar w432 = w122*w431;
    const scalar w433 = w177*w432;
    const scalar w434 = w182*y - w184*y;
    const scalar w435 = w18*y - w21*y + w434;
    const scalar w436 = -w433 - w435;
    const scalar w437 = w423*w57;
    const scalar w438 = w424*w62;
    const scalar w439 = gil*w437 + git*w438 + w193*y - w195*y;
    const scalar w440 = gel*w437 - 2*get*theta1*w14*w15*w3*w78*y + get*w438 + w207*y + w439;
    const scalar w441 = -w440;
    const scalar w442 = gil*w321 + git*w323 + w214*y - w215*y;
    const scalar w443 = w324 + w442;
    const scalar w444 = -w443;
    const scalar w445 = w188*w356;
    const scalar w446 = w128*w354 + w132*w441 - w133*w445 + w140*w444;
    const scalar w447 = w32*(w127*w212 + w30*w436 + w446);
    const scalar w448 = w426*w56 + w427*w61;
    const scalar w449 = w428 + w448;
    const scalar w450 = -w433 - w449;
    const scalar w451 = w320*w75 + w322*w79 + w423*w56 + w424*w61;
    const scalar w452 = w425 + w451;
    const scalar w453 = -w452;
    const scalar w454 = w132*w453 + w213 + w30*w444 + w379;
    const scalar w455 = w140*w450 - w143*w445 + w454;
    const scalar w456 = -w246*y - w247*y;
    const scalar w457 = -w249*y + w250*y;
    const scalar w458 = w456*w75 + w457*w79;
    const scalar w459 = w456*w83 + w457*w85 + w458;
    const scalar w460 = -w433 - w459;
    const scalar w461 = w128*w212 + w140*w453 - w144*w445 + w30*w441;
    const scalar w462 = w32*(w132*w460 + w147*w354 + w461);
    const scalar w463 = w146 + w429*w455 + w447*z + w462*x;
    const scalar w464 = w263*w432;
    const scalar w465 = -w435 - w464;
    const scalar w466 = w153*(w167*w212 + w30*w465 + w446);
    const scalar w467 = -w449 - w464;
    const scalar w468 = w140*w467 - w169*w445 + w454;
    const scalar w469 = -w459 - w464;
    const scalar w470 = w153*(w132*w469 + w171*w354 + w461);
    const scalar w471 = -w170 - w430*w468 - w466*z - w470*x;
    const scalar w472 = w32*w435;
    const scalar w473 = w32*w443;
    const scalar w474 = w32*w440;
    const scalar w475 = w32*(w472*z + w473*y + w474*x + w71);
    const scalar w476 = w475*z;
    const scalar w477 = w32*w452;
    const scalar w478 = w429*w449 + w473*z + w477*x + w97;
    const scalar w479 = w429*w478;
    const scalar w480 = w32*w459;
    const scalar w481 = w32*(w102 + w474*z + w477*y + w480*x);
    const scalar w482 = w481*x;
    const scalar w483 = -w103 - w476 - w479 - w482;
    const scalar w484 = w275*w483;
    const scalar w485 = w153*w435;
    const scalar w486 = w153*w443;
    const scalar w487 = w153*w440;
    const scalar w488 = w153*(-w156 - w485*z - w486*y - w487*x);
    const scalar w489 = w488*z;
    const scalar w490 = w153*w452;
    const scalar w491 = -w159 - w430*w449 - w486*z - w490*x;
    const scalar w492 = w430*w491;
    const scalar w493 = w153*w459;
    const scalar w494 = w153*(-w160 - w487*z - w490*y - w493*x);
    const scalar w495 = w494*x;
    const scalar w496 = w161 + w489 + w492 + w495;
    const scalar w497 = w152*w429 - w175*w430 + w261*w463 + w274*w471 + w291*w484 + w306*w496;
    const scalar w498 = (1.0/10.0)*Tt + w311*w497 + w314*y + w317*y;
    const scalar w499 = w140 + w498;
    const scalar w500 = hA2*w179;
    const scalar w501 = St*(hA0 + w114*w500);
    const scalar w502 = w32*w42;
    const scalar w503 = w32*w64;
    const scalar w504 = w32*w81;
    const scalar w505 = w32*(w502*z + w503*y + w504*x);
    const scalar w506 = w32*w94;
    const scalar w507 = w32*w99;
    const scalar w508 = w32*(w503*z + w506*y + w507*x);
    const scalar w509 = w107*w32;
    const scalar w510 = w32*(w504*z + w507*y + w509*x);
    const scalar w511 = w505*z + w508*y + w510*x;
    const scalar w512 = (1.0/(w511));
    const scalar w513 = -w30 + w501;
    const scalar w514 = hA1*w117;
    const scalar w515 = (3.0/5.0)*w514;
    const scalar w516 = w515*x;
    const scalar w517 = hA2*w119 - w516;
    const scalar w518 = St*w517 - w132;
    const scalar w519 = 1 - 3.0/10.0*w7;
    const scalar w520 = St*k;
    const scalar w521 = hA1*w431;
    const scalar w522 = w520*w521;
    const scalar w523 = w519*w522;
    const scalar w524 = -w140 - w523;
    const scalar w525 = w42*w513 + w518*w81 + w524*w64;
    const scalar w526 = w32*w525;
    const scalar w527 = w513*w64 + w518*w99 + w524*w94;
    const scalar w528 = w32*w527;
    const scalar w529 = w107*w518 + w513*w81 + w524*w99;
    const scalar w530 = w32*w529;
    const scalar w531 = -w526*z - w528*y - w530*x;
    const scalar w532 = w512*w531;
    const scalar w533 = w35*w532;
    const scalar w534 = w153*w42;
    const scalar w535 = w153*w64;
    const scalar w536 = w153*w81;
    const scalar w537 = w153*(-w534*z - w535*y - w536*x);
    const scalar w538 = w153*w94;
    const scalar w539 = w153*w99;
    const scalar w540 = w153*(-w535*z - w538*y - w539*x);
    const scalar w541 = w107*w153;
    const scalar w542 = w153*(-w536*z - w539*y - w541*x);
    const scalar w543 = -w537*z - w540*y - w542*x;
    const scalar w544 = (1.0/(w543));
    const scalar w545 = w153*w525;
    const scalar w546 = w153*w527;
    const scalar w547 = w153*w529;
    const scalar w548 = w545*z + w546*y + w547*x;
    const scalar w549 = w544*w548;
    const scalar w550 = w35*w549;
    const scalar w551 = hA2*w119*w356;
    const scalar w552 = w551*x;
    const scalar w553 = -w190 - w552;
    const scalar w554 = w500*w520;
    const scalar w555 = w554*w81;
    const scalar w556 = -w212*w64;
    const scalar w557 = w206*w518 + w223*w524 + w276*w513 + w42*w553 + w555 + w556;
    const scalar w558 = -w212*w94 + w223*w513 + w237*w524 + w240*w518 + w553*w64 + w554*w99;
    const scalar w559 = w32*w558;
    const scalar w560 = w107*w554 + w206*w513 - w212*w99 + w240*w524 + w252*w518 + w553*w81;
    const scalar w561 = w32*w560;
    const scalar w562 = -w33*w557 - w526 - w559*y - w561*x;
    const scalar w563 = w260*w512;
    const scalar w564 = w153*w558;
    const scalar w565 = w153*w560;
    const scalar w566 = w154*w557 + w545 + w564*y + w565*x;
    const scalar w567 = w273*w544;
    const scalar w568 = w223*w32;
    const scalar w569 = w206*w32;
    const scalar w570 = w276*w33 + w502 + w568*y + w569*x;
    const scalar w571 = w33*w570;
    const scalar w572 = w237*w32;
    const scalar w573 = w240*w32;
    const scalar w574 = w32*(w503 + w568*z + w572*y + w573*x);
    const scalar w575 = w574*y;
    const scalar w576 = w252*w32;
    const scalar w577 = w32*(w504 + w569*z + w573*y + w576*x);
    const scalar w578 = w577*x;
    const scalar w579 = -w505 - w571 - w575 - w578;
    const scalar w580 = (1.0/(w511*w511));
    const scalar w581 = w531*w580;
    const scalar w582 = w260*w581;
    const scalar w583 = w153*w223;
    const scalar w584 = w153*w206;
    const scalar w585 = -w154*w276 - w534 - w583*y - w584*x;
    const scalar w586 = w154*w585;
    const scalar w587 = w153*w237;
    const scalar w588 = w153*w240;
    const scalar w589 = w153*(-w535 - w583*z - w587*y - w588*x);
    const scalar w590 = w589*y;
    const scalar w591 = w153*w252;
    const scalar w592 = w153*(-w536 - w584*z - w588*y - w591*x);
    const scalar w593 = w592*x;
    const scalar w594 = w537 + w586 + w590 + w593;
    const scalar w595 = (1.0/(w543*w543));
    const scalar w596 = w273*w595;
    const scalar w597 = w548*w596;
    const scalar w598 = -w154*w550 + w33*w533 + w562*w563 + w566*w567 + w579*w582 + w594*w597;
    const scalar w599 = w260*w532 + w273*w549;
    const scalar w600 = 2*w599;
    const scalar w601 = w310*w600;
    const scalar w602 = w316*w600;
    const scalar w603 = w311*w598 + w318 + w501 + w601*z + w602*z;
    const scalar w604 = (3.0/5.0)*St*hA1*k*w431*x - w354;
    const scalar w605 = -St*w515;
    const scalar w606 = w357 + w605;
    const scalar w607 = w343*w518 + w350*w524 + w392*w513 + w42*w554 + w604*w64 + w606*w81;
    const scalar w608 = w32*w607;
    const scalar w609 = w350*w513 + w365*w524 + w370*w518 + w554*w64 + w604*w94 + w606*w99;
    const scalar w610 = w32*w609;
    const scalar w611 = w604*w99;
    const scalar w612 = w107*w606 + w343*w513 + w370*w524 + w376*w518 + w555 + w611;
    const scalar w613 = -w333*w612 - w530 - w608*z - w610*y;
    const scalar w614 = w153*w607;
    const scalar w615 = w153*w609;
    const scalar w616 = w334*w612 + w547 + w614*z + w615*y;
    const scalar w617 = w32*w392;
    const scalar w618 = w32*w350;
    const scalar w619 = w32*w343;
    const scalar w620 = w32*(w504 + w617*z + w618*y + w619*x);
    const scalar w621 = w620*z;
    const scalar w622 = w32*w370;
    const scalar w623 = w333*w376 + w509 + w619*z + w622*y;
    const scalar w624 = w333*w623;
    const scalar w625 = w32*w365;
    const scalar w626 = w32*(w507 + w618*z + w622*x + w625*y);
    const scalar w627 = w626*y;
    const scalar w628 = -w510 - w621 - w624 - w627;
    const scalar w629 = w153*w392;
    const scalar w630 = w153*w350;
    const scalar w631 = w153*w343;
    const scalar w632 = w153*(-w536 - w629*z - w630*y - w631*x);
    const scalar w633 = w632*z;
    const scalar w634 = w153*w370;
    const scalar w635 = -w334*w376 - w541 - w631*z - w634*y;
    const scalar w636 = w334*w635;
    const scalar w637 = w153*w365;
    const scalar w638 = w153*(-w539 - w630*z - w634*x - w637*y);
    const scalar w639 = w638*y;
    const scalar w640 = w542 + w633 + w636 + w639;
    const scalar w641 = w333*w533 - w334*w550 + w563*w613 + w567*w616 + w582*w628 + w597*w640;
    const scalar w642 = St*w517 + w311*w641 + w421 + w601*x + w602*x;
    const scalar w643 = -w356*w514*w519;
    const scalar w644 = w445 + w643;
    const scalar w645 = -w212*w42 + w434*w513 + w439*w518 + w442*w524 + w604*w81 + w64*w644;
    const scalar w646 = w32*w645;
    const scalar w647 = w442*w513 + w448*w524 + w451*w518 + w556 + w611 + w644*w94;
    const scalar w648 = w107*w604 - w212*w81 + w439*w513 + w451*w524 + w458*w518 + w644*w99;
    const scalar w649 = w32*w648;
    const scalar w650 = -w429*w647 - w528 - w646*z - w649*x;
    const scalar w651 = w153*w645;
    const scalar w652 = w153*w648;
    const scalar w653 = w430*w647 + w546 + w651*z + w652*x;
    const scalar w654 = w32*w434;
    const scalar w655 = w32*w442;
    const scalar w656 = w32*w439;
    const scalar w657 = w32*(w503 + w654*z + w655*y + w656*x);
    const scalar w658 = w657*z;
    const scalar w659 = w32*w451;
    const scalar w660 = w429*w448 + w506 + w655*z + w659*x;
    const scalar w661 = w429*w660;
    const scalar w662 = w32*w458;
    const scalar w663 = w32*(w507 + w656*z + w659*y + w662*x);
    const scalar w664 = w663*x;
    const scalar w665 = -w508 - w658 - w661 - w664;
    const scalar w666 = w153*w434;
    const scalar w667 = w153*w442;
    const scalar w668 = w153*w439;
    const scalar w669 = w153*(-w535 - w666*z - w667*y - w668*x);
    const scalar w670 = w669*z;
    const scalar w671 = w153*w451;
    const scalar w672 = -w430*w448 - w538 - w667*z - w671*x;
    const scalar w673 = w430*w672;
    const scalar w674 = w153*w458;
    const scalar w675 = w153*(-w539 - w668*z - w671*y - w674*x);
    const scalar w676 = w675*x;
    const scalar w677 = w540 + w670 + w673 + w676;
    const scalar w678 = w429*w533 - w430*w550 + w563*w650 + w567*w653 + w582*w665 + w597*w677;
    const scalar w679 = w311*w678 + w498 - w523 + w601*y + w602*y;
    const scalar w680 = 8*w315;
    const scalar w681 = w312*w680;
    const scalar w682 = 2*w153;
    const scalar w683 = sigma0In*w355;
    const scalar w684 = w116*w683;
    const scalar w685 = w120*w178;
    const scalar w686 = w684*w685;
    const scalar w687 = (theta1*theta1);
    const scalar w688 = (1.0/(w2*w2));
    const scalar w689 = w687*w688;
    const scalar w690 = w689*w9;
    const scalar w691 = 8*w690;
    const scalar w692 = w43*w691;
    const scalar w693 = w44*w691;
    const scalar w694 = w39*w691;
    const scalar w695 = w41*w691;
    const scalar w696 = w40*w691;
    const scalar w697 = gel*w696;
    const scalar w698 = w38*w691;
    const scalar w699 = get*w698;
    const scalar w700 = gil*w696;
    const scalar w701 = git*w698;
    const scalar w702 = w18 + w182 - w184 - w21;
    const scalar w703 = -w692 - w693 - w694 - w695 + w697 + w699 + w700 + w701 + w702;
    const scalar w704 = (8.0/25.0)*w189;
    const scalar w705 = (8.0/125.0)*w189;
    const scalar w706 = 4*w690;
    const scalar w707 = w196 + w197;
    const scalar w708 = 4*w14*w36*w51*w687*w688*w9 - w200*w706 - w707;
    const scalar w709 = w57*w708;
    const scalar w710 = -w202 + w203;
    const scalar w711 = w245*w690 + w706*w77 + w710;
    const scalar w712 = w62*w711;
    const scalar w713 = w62*z;
    const scalar w714 = w198*w713;
    const scalar w715 = w5*w57;
    const scalar w716 = w715*z;
    const scalar w717 = w204*w716;
    const scalar w718 = gil*w5;
    const scalar w719 = -w193 + w195;
    const scalar w720 = -gil*w709 - git*w712 + git*w717 + w706*w76 + w706*w80 - w714*w718 + w719;
    const scalar w721 = -w207 + w208;
    const scalar w722 = -gel*w709 - get*w712 + get*w717 - w6*w714 + w706*w84 + w706*w86 + w720 + w721;
    const scalar w723 = -2*theta1*w13*w3*w36*w48 + w216;
    const scalar w724 = -w52*w706 - w54*w706 - w723;
    const scalar w725 = w57*w724;
    const scalar w726 = w200*w36;
    const scalar w727 = w219 + w220;
    const scalar w728 = 4*w14*w51*w687*w688*w9 - w706*w726 - w727;
    const scalar w729 = w62*w728;
    const scalar w730 = w217*w713;
    const scalar w731 = w221*w716;
    const scalar w732 = -w214 + w215;
    const scalar w733 = -gil*w725 - git*w729 + git*w731 + w58*w706 + w63*w706 - w718*w730 + w732;
    const scalar w734 = -w224 + w225;
    const scalar w735 = -gel*w725 - get*w729 + get*w731 - w6*w730 + w66*w706 + w68*w706 + w733 + w734;
    const scalar w736 = w212*w228;
    const scalar w737 = (4.0/25.0)*St*w139;
    const scalar w738 = -w133*w737;
    const scalar w739 = w132*w722 + w140*w735 + 2*w736 + w738;
    const scalar w740 = -4*theta1*w13*w3*w36*w48 + w232;
    const scalar w741 = -w52*w691 - w54*w691 - w740;
    const scalar w742 = w217*w233;
    const scalar w743 = w234 + w235;
    const scalar w744 = 8*w14*w51*w687*w688*w9 - w691*w726 - w743;
    const scalar w745 = w221*w236;
    const scalar w746 = gil*w742 + git*w745 + w56*w741 + w61*w744;
    const scalar w747 = gel*w742 + get*w745 + w65*w741 + w67*w744 + w746;
    const scalar w748 = 2*w212;
    const scalar w749 = 2*gel;
    const scalar w750 = w198*w217;
    const scalar w751 = 2*get;
    const scalar w752 = w204*w221;
    const scalar w753 = 2*gil;
    const scalar w754 = 2*git;
    const scalar w755 = w56*w708 + w61*w711 + w724*w75 + w728*w79 + w750*w753 + w752*w754;
    const scalar w756 = w65*w708 + w67*w711 + w724*w83 + w728*w85 + w749*w750 + w751*w752 + w755;
    const scalar w757 = -w756;
    const scalar w758 = w132*w757 + w133*w705 + w228*w704 + w30*w735;
    const scalar w759 = w246 + w247;
    const scalar w760 = 8*w14*w36*w51*w687*w688*w9 - w200*w691 - w759;
    const scalar w761 = w198*w248;
    const scalar w762 = -w249 + w250;
    const scalar w763 = w53*w691 + w691*w77 + w762;
    const scalar w764 = w204*w251;
    const scalar w765 = gil*w761 + git*w764 + w75*w760 + w763*w79;
    const scalar w766 = gel*w761 + get*w764 + w760*w83 + w763*w85 + w765;
    const scalar w767 = w128*w705 + w140*w757 - w144*w737 + w211*w704 + w242*w748 + w30*w722;
    const scalar w768 = w272*w34;
    const scalar w769 = w304*w768;
    const scalar w770 = -w182 + w184;
    const scalar w771 = w694 + w695 - w700 - w701 + w770;
    const scalar w772 = -w18 + w21;
    const scalar w773 = -w692 - w693 + w697 + w699 - w771 - w772;
    const scalar w774 = -w735;
    const scalar w775 = -w722;
    const scalar w776 = 2*w32;
    const scalar w777 = sigma0Out*w355;
    const scalar w778 = w116*w777;
    const scalar w779 = w685*w778;
    const scalar w780 = w154*w35;
    const scalar w781 = w165*w780;
    const scalar w782 = w113*w35;
    const scalar w783 = w33*w782;
    const scalar w784 = w259*w34;
    const scalar w785 = w275*w784;
    const scalar w786 = 2*w158 + 2*w295 + 2*w299 + 2*w302;
    const scalar w787 = w305/(w164*w164*w164);
    const scalar w788 = w303*w787;
    const scalar w789 = w275*w291;
    const scalar w790 = w303*w304;
    const scalar w791 = w173*w780;
    const scalar w792 = w150*w35;
    const scalar w793 = w33*w792;
    const scalar w794 = -2*w281 - 2*w285 - 2*w288 - 2*w90;
    const scalar w795 = w291/(w112*w112*w112);
    const scalar w796 = w289*w795;
    const scalar w797 = w32*w35;
    const scalar w798 = w153*w35;
    const scalar w799 = w151*w797 - w174*w798;
    const scalar w800 = 4*w310;
    const scalar w801 = w307*z;
    const scalar w802 = 4*w316;
    const scalar w803 = w314 + w317;
    const scalar w804 = w31 + w311*(w257*w289*w785 + 2*w257*w783 + w261*(w230*w776 + w33*(w127*w705 + w187*w704 + w30*(-w703 - w779) + w739) + w333*(w132*(-w766 - w779) + w767) + w429*(w140*(-w747 - w779) - w143*w737 + w239*w748 + w758)) + w271*w303*w769 - 2*w271*w781 + w274*(-w154*(w167*w705 + w264*w704 + w30*(-w686 - w703) + w739) - w265*w682 - w334*(w132*(-w686 - w766) + w767) - w430*(w140*(-w686 - w747) - w169*w737 + w267*w748 + w758)) + 2*w290*w793 + w306*(w154*(-w154*w773 - w277*w682 - w334*w775 - w430*w774) + w294*w682 + w334*(-w154*w775 - 2*w293 - w334*w766 - w430*w756) + w430*(-w154*w774 - 2*w292 - w334*w756 - w430*w747)) + w786*w788 + w789*(-w280*w776 - w33*(w277*w776 + w33*w773 + w333*w775 + w429*w774) - w333*(2*w279 + w33*w775 + w333*w766 + w429*w756) - w429*(2*w278 + w33*w774 + w333*w756 + w429*w747)) - 2*w790*w791 + w794*w796 + w799) + w681*w9 + w800*w801 + w801*w802 + w803;
    const scalar w805 = w124*w355;
    const scalar w806 = sigma0In*w805;
    const scalar w807 = w702 + w806;
    const scalar w808 = w689*w7;
    const scalar w809 = 8*w808;
    const scalar w810 = w43*w809;
    const scalar w811 = w44*w809;
    const scalar w812 = w39*w809;
    const scalar w813 = w41*w809;
    const scalar w814 = w40*w809;
    const scalar w815 = gel*w814;
    const scalar w816 = w38*w809;
    const scalar w817 = get*w816;
    const scalar w818 = gil*w814;
    const scalar w819 = git*w816;
    const scalar w820 = -w810 - w811 - w812 - w813 + w815 + w817 + w818 + w819;
    const scalar w821 = 4*w808;
    const scalar w822 = 4*w14*w36*w51*w687*w688*w7 - w200*w821 - w707;
    const scalar w823 = w57*w822;
    const scalar w824 = w245*w808 + w710 + w77*w821;
    const scalar w825 = w62*w824;
    const scalar w826 = w62*x;
    const scalar w827 = w325*w826;
    const scalar w828 = w715*x;
    const scalar w829 = w327*w828;
    const scalar w830 = -gil*w823 - git*w825 + git*w829 - w718*w827 + w719 + w76*w821 + w80*w821;
    const scalar w831 = -gel*w823 - get*w825 + get*w829 - w6*w827 + w721 + w821*w84 + w821*w86 + w830;
    const scalar w832 = -w52*w821 - w54*w821 - w723;
    const scalar w833 = w57*w832;
    const scalar w834 = 4*w14*w51*w687*w688*w7 - w726*w821 - w727;
    const scalar w835 = w62*w834;
    const scalar w836 = w346*w826;
    const scalar w837 = w348*w828;
    const scalar w838 = -gil*w833 - git*w835 + git*w837 + w58*w821 + w63*w821 - w718*w836 + w732;
    const scalar w839 = -gel*w833 - get*w835 + get*w837 - w6*w836 + w66*w821 + w68*w821 + w734 + w838;
    const scalar w840 = w136*w356;
    const scalar w841 = (k*k*k);
    const scalar w842 = St*w841;
    const scalar w843 = fA2*w129*w130*w842;
    const scalar w844 = 2*w357;
    const scalar w845 = -w128*w843 + w132*w831 - w133*w840 + w140*w839 - w345*w844 + w352*w353;
    const scalar w846 = -w52*w809 - w54*w809 - w740;
    const scalar w847 = w346*w361;
    const scalar w848 = 8*w14*w51*w687*w688*w7 - w726*w809 - w743;
    const scalar w849 = w348*w363;
    const scalar w850 = gil*w847 + git*w849 + w56*w846 + w61*w848;
    const scalar w851 = gel*w847 + get*w849 + w65*w846 + w67*w848 + w850;
    const scalar w852 = w325*w346;
    const scalar w853 = w327*w348;
    const scalar w854 = w56*w822 + w61*w824 + w75*w832 + w753*w852 + w754*w853 + w79*w834;
    const scalar w855 = w65*w822 + w67*w824 + w749*w852 + w751*w853 + w83*w832 + w834*w85 + w854;
    const scalar w856 = -w855;
    const scalar w857 = w132*w856 - w144*w843 + w30*w839 - w373*w844;
    const scalar w858 = 8*w14*w36*w51*w687*w688*w7 - w200*w809 - w759;
    const scalar w859 = w325*w330;
    const scalar w860 = w53*w809 + w762 + w77*w809;
    const scalar w861 = w327*w331;
    const scalar w862 = gil*w859 + git*w861 + w75*w858 + w79*w860;
    const scalar w863 = gel*w859 + get*w861 + w83*w858 + w85*w860 + w862;
    const scalar w864 = -w144*w840;
    const scalar w865 = w140*w856 + w30*w831 + w353*w373 + w864;
    const scalar w866 = w770 + w812 + w813 - w818 - w819;
    const scalar w867 = -w772 - w810 - w811 + w815 + w817 - w866;
    const scalar w868 = -w839;
    const scalar w869 = -w831;
    const scalar w870 = sigma0Out*w805;
    const scalar w871 = w702 + w870;
    const scalar w872 = w334*w35;
    const scalar w873 = 2*w872;
    const scalar w874 = w333*w782;
    const scalar w875 = 2*w163 + 2*w411 + 2*w414 + 2*w417;
    const scalar w876 = w418*w787;
    const scalar w877 = w173*w419;
    const scalar w878 = w333*w792;
    const scalar w879 = -2*w111 - 2*w398 - 2*w401 - 2*w404;
    const scalar w880 = w405*w795;
    const scalar w881 = w420*x;
    const scalar w882 = -w130*w31*w355 + w311*(w261*(w33*(w30*(-w820 - w871) + w845) + w333*(w132*(-w863 - w870) - w147*w843 - w378*w844 + w865) + w381*w776 + w429*(w140*(-w851 - w870) - w143*w840 + w353*w367 + w857)) + w274*(-w154*(w30*(-w807 - w820) + w845) - w334*(w132*(-w806 - w863) - w171*w843 - w388*w844 + w865) - w389*w682 - w430*(w140*(-w806 - w851) - w169*w840 + w353*w386 + w857)) + w306*(w154*(-w154*w867 - w334*w869 - 2*w409 - w430*w868) + w334*(-w154*w869 - w334*w863 - w377*w682 - w430*w855) + w413*w682 + w430*(-w154*w868 - w334*w855 - 2*w412 - w430*w851)) + w382*w406*w784 + 2*w382*w874 + w390*w419*w768 - w391*w873 + 2*w406*w878 + w789*(-w33*(w33*w867 + w333*w869 + 2*w396 + w429*w868) - w333*(w33*w869 + w333*w863 + w377*w776 + w429*w855) - w400*w776 - w429*(w33*w868 + w333*w855 + 2*w399 + w429*w851)) + w799 - w873*w877 + w875*w876 + w879*w880) + w681*w7 + w800*w881 + w802*w881 + w803;
    const scalar w883 = w689*w8;
    const scalar w884 = 8*w883;
    const scalar w885 = w43*w884;
    const scalar w886 = w44*w884;
    const scalar w887 = w39*w884;
    const scalar w888 = w41*w884;
    const scalar w889 = w40*w884;
    const scalar w890 = gel*w889;
    const scalar w891 = w38*w884;
    const scalar w892 = get*w891;
    const scalar w893 = gil*w889;
    const scalar w894 = git*w891;
    const scalar w895 = -w885 - w886 - w887 - w888 + w890 + w892 + w893 + w894;
    const scalar w896 = w29*w356;
    const scalar w897 = 4*w883;
    const scalar w898 = 4*w14*w36*w51*w687*w688*w8 - w200*w897 - w707;
    const scalar w899 = w57*w898;
    const scalar w900 = w245*w883 + w710 + w77*w897;
    const scalar w901 = w62*w900;
    const scalar w902 = w62*y;
    const scalar w903 = w423*w902;
    const scalar w904 = w715*y;
    const scalar w905 = w424*w904;
    const scalar w906 = -gil*w899 - git*w901 + git*w905 - w718*w903 + w719 + w76*w897 + w80*w897;
    const scalar w907 = -gel*w899 - get*w901 + get*w905 - w6*w903 + w721 + w84*w897 + w86*w897 + w906;
    const scalar w908 = -w52*w897 - w54*w897 - w723;
    const scalar w909 = w57*w908;
    const scalar w910 = 4*w14*w51*w687*w688*w8 - w726*w897 - w727;
    const scalar w911 = w62*w910;
    const scalar w912 = w320*w902;
    const scalar w913 = w322*w904;
    const scalar w914 = -gil*w909 - git*w911 + git*w913 + w58*w897 + w63*w897 - w718*w912 + w732;
    const scalar w915 = -gel*w909 - get*w911 + get*w913 - w6*w912 + w66*w897 + w68*w897 + w734 + w914;
    const scalar w916 = w137*w27*w842;
    const scalar w917 = 2*w445;
    const scalar w918 = w132*w907 + w133*w916 + w140*w915 + w353*w441 - w444*w917;
    const scalar w919 = -w52*w884 - w54*w884 - w740;
    const scalar w920 = w320*w426;
    const scalar w921 = 8*w14*w51*w687*w688*w8 - w726*w884 - w743;
    const scalar w922 = w322*w427;
    const scalar w923 = gil*w920 + git*w922 + w56*w919 + w61*w921;
    const scalar w924 = gel*w920 + get*w922 + w65*w919 + w67*w921 + w923;
    const scalar w925 = w320*w423;
    const scalar w926 = w322*w424;
    const scalar w927 = w56*w898 + w61*w900 + w75*w908 + w753*w925 + w754*w926 + w79*w910;
    const scalar w928 = w65*w898 + w67*w900 + w749*w925 + w751*w926 + w83*w908 + w85*w910 + w927;
    const scalar w929 = -w928;
    const scalar w930 = w212*w444;
    const scalar w931 = -w133*w896;
    const scalar w932 = w132*w929 + w30*w915 + w353*w453 + 2*w930 + w931;
    const scalar w933 = 8*w14*w36*w51*w687*w688*w8 - w200*w884 - w759;
    const scalar w934 = w423*w456;
    const scalar w935 = w53*w884 + w762 + w77*w884;
    const scalar w936 = w424*w457;
    const scalar w937 = gil*w934 + git*w936 + w75*w933 + w79*w935;
    const scalar w938 = gel*w934 + get*w936 + w83*w933 + w85*w935 + w937;
    const scalar w939 = -w128*w896 + w140*w929 + w144*w916 + w30*w907 + w441*w748 - w453*w917;
    const scalar w940 = w770 + w887 + w888 - w893 - w894;
    const scalar w941 = -w772 - w885 - w886 + w890 + w892 - w940;
    const scalar w942 = -w915;
    const scalar w943 = -w907;
    const scalar w944 = w35*w430;
    const scalar w945 = w165*w944;
    const scalar w946 = w429*w782;
    const scalar w947 = 2*w161 + 2*w489 + 2*w492 + 2*w495;
    const scalar w948 = w496*w787;
    const scalar w949 = w304*w496;
    const scalar w950 = w173*w944;
    const scalar w951 = w429*w792;
    const scalar w952 = -2*w103 - 2*w476 - 2*w479 - 2*w482;
    const scalar w953 = w483*w795;
    const scalar w954 = w497*y;
    const scalar w955 = w311*(w261*(w33*(-w127*w896 + w30*(-w871 - w895) + w436*w748 + w918) + w333*(w132*(-w870 - w938) + w353*w460 + w939) + w429*(w140*(-w870 - w924) + w143*w916 - w450*w917 + w932) + w455*w776) + w274*(-w154*(-w167*w896 + w30*(-w807 - w895) + w465*w748 + w918) - w334*(w132*(-w806 - w938) + w353*w469 + w939) - w430*(w140*(-w806 - w924) + w169*w916 - w467*w917 + w932) - w468*w682) + w306*(w154*(-w154*w941 - w334*w943 - w430*w942 - 2*w486) + w334*(-w154*w943 - w334*w938 - w430*w928 - 2*w490) + w430*(-w154*w942 - w334*w928 - w430*w924 - w449*w682) + w491*w682) + w463*w483*w785 + 2*w463*w946 + w471*w496*w769 - 2*w471*w945 + 2*w484*w951 + w789*(-w33*(w33*w941 + w333*w943 + w429*w942 + 2*w473) - w333*(w33*w943 + w333*w938 + w429*w928 + 2*w477) - w429*(w33*w942 + w333*w928 + w429*w924 + w449*w776) - w478*w776) + w799 + w947*w948 - 2*w949*w950 + w952*w953) + w681*w8 + w800*w954 + w802*w954 + w803;
    const scalar w956 = -1.0/2.0*sigma0In*sigmaEps*w117*w179*w335*w355;
    const scalar w957 = w198*w330;
    const scalar w958 = w204*w331;
    const scalar w959 = w689*z;
    const scalar w960 = w959*x;
    const scalar w961 = 8*w960;
    const scalar w962 = -w200*w961 + w73*w961;
    const scalar w963 = w53*w961 + w77*w961;
    const scalar w964 = w75*w962 + w79*w963;
    const scalar w965 = gil*w957 + git*w958 + w964;
    const scalar w966 = w83*w962 + w85*w963;
    const scalar w967 = gel*w957 + get*w958 + w965 + w966;
    const scalar w968 = 4*w960;
    const scalar w969 = -w200*w968 + w73*w968;
    const scalar w970 = w57*w969;
    const scalar w971 = w245*w960 + w77*w968;
    const scalar w972 = w62*w971;
    const scalar w973 = w4*w749;
    const scalar w974 = w198*w826;
    const scalar w975 = w325*w713;
    const scalar w976 = w4*w57;
    const scalar w977 = w751*w976;
    const scalar w978 = w204*x;
    const scalar w979 = w327*z;
    const scalar w980 = w4*w753;
    const scalar w981 = w754*w976;
    const scalar w982 = -gil*w970 - git*w972 + w76*w968 + w80*w968 - w974*w980 - w975*w980 + w978*w981 + w979*w981;
    const scalar w983 = -gel*w970 - get*w972 + w84*w968 + w86*w968 - w973*w974 - w973*w975 + w977*w978 + w977*w979 + w982;
    const scalar w984 = -w52*w968 - w54*w968;
    const scalar w985 = w217*w325;
    const scalar w986 = w198*w346;
    const scalar w987 = w59*w968 - w726*w968;
    const scalar w988 = w221*w327;
    const scalar w989 = w204*w348;
    const scalar w990 = gil*w985 + gil*w986 + git*w988 + git*w989 + w56*w969 + w61*w971 + w75*w984 + w79*w987;
    const scalar w991 = gel*w985 + gel*w986 + get*w988 + get*w989 + w65*w969 + w67*w971 + w83*w984 + w85*w987 + w990;
    const scalar w992 = -w991;
    const scalar w993 = w242*w354;
    const scalar w994 = w140*w992 + w190*w345 + w212*w373 + w30*w983 + w993;
    const scalar w995 = -w269*w357 + w994;
    const scalar w996 = w217*w361;
    const scalar w997 = w221*w363;
    const scalar w998 = -w52*w961 - w54*w961;
    const scalar w999 = w59*w961 - w726*w961;
    const scalar w1000 = w56*w998 + w61*w999;
    const scalar w1001 = gil*w996 + git*w997 + w1000;
    const scalar w1002 = w65*w998 + w67*w999;
    const scalar w1003 = gel*w996 + get*w997 + w1001 + w1002;
    const scalar w1004 = w57*w984;
    const scalar w1005 = w62*w987;
    const scalar w1006 = w217*w826;
    const scalar w1007 = w346*w713;
    const scalar w1008 = w221*x;
    const scalar w1009 = w348*z;
    const scalar w1010 = -gil*w1004 - git*w1005 - w1006*w980 - w1007*w980 + w1008*w981 + w1009*w981 + w58*w968 + w63*w968;
    const scalar w1011 = -gel*w1004 - get*w1005 - w1006*w973 - w1007*w973 + w1008*w977 + w1009*w977 + w1010 + w66*w968 + w68*w968;
    const scalar w1012 = w1011*w30 + w132*w992 + w190*w352 - w242*w357;
    const scalar w1013 = w1012 + w212*w386 + w267*w354;
    const scalar w1014 = w181*w335;
    const scalar w1015 = w40*w961;
    const scalar w1016 = w38*w961;
    const scalar w1017 = -gil*w1015 - git*w1016 + w39*w961 + w41*w961;
    const scalar w1018 = -gel*w1015 - get*w1016 + w1017 + w43*w961 + w44*w961;
    const scalar w1019 = w212*w352;
    const scalar w1020 = w1011*w140 + w1019 + w132*w983 - w211*w357 + w228*w354;
    const scalar w1021 = w154*(w1020 + w190*w383 + w30*(w1014*w683 + w1018)) + w270 + w384;
    const scalar w1022 = -w983;
    const scalar w1023 = w1022*w154 + w300 + w409 + w430*w991;
    const scalar w1024 = -w1011;
    const scalar w1025 = w1024*w154 + w297 + w334*w991 + w408;
    const scalar w1026 = -w1018;
    const scalar w1027 = w154*(-w1022*w334 - w1024*w430 - w1026*w154 - w293 - w407) + w301 + w410;
    const scalar w1028 = -1.0/2.0*sigma0Out*sigmaEps*w117*w179*w335*w355;
    const scalar w1029 = -w254*w357 + w994;
    const scalar w1030 = w1012 + w212*w367 + w239*w354;
    const scalar w1031 = w256 + w33*(w1020 + w190*w342 + w30*(w1014*w777 + w1018)) + w359;
    const scalar w1032 = w1022*w33 + w286 + w396 + w429*w991;
    const scalar w1033 = w1024*w33 + w283 + w333*w991 + w395;
    const scalar w1034 = w287 + w33*(w1022*w333 + w1024*w429 + w1026*w33 + w279 + w394) + w397;
    const scalar w1035 = w382*w782;
    const scalar w1036 = w273*w419;
    const scalar w1037 = w273*w390;
    const scalar w1038 = w165*w872;
    const scalar w1039 = w406*w792;
    const scalar w1040 = w260*w406;
    const scalar w1041 = w260*w382;
    const scalar w1042 = w173*w872;
    const scalar w1043 = w1035*w33 + w1036*w271 + w1037*w790 - w1038*w271 + w1039*w33 + w1040*w257 + w1041*w290 - w1042*w790 + w257*w874 + w290*w878 - w391*w780 - w780*w877;
    const scalar w1044 = 2*w310;
    const scalar w1045 = w307*x;
    const scalar w1046 = 2*w316;
    const scalar w1047 = w420*z;
    const scalar w1048 = w681*z;
    const scalar w1049 = w1044*w1045 + w1044*w1047 + w1045*w1046 + w1046*w1047 + w1048*x;
    const scalar w1050 = w1049 + w311*(w1043 + w261*(w1031 + w333*(w1029 + w132*(-w1028 - w967)) + w429*(w1030 + w140*(-w1003 - w1028))) + w274*(-w1021 - w334*(w132*(-w956 - w967) + w995) - w430*(w1013 + w140*(-w1003 - w956))) + w306*(w1027 + w334*(-w1023 - w334*w967) + w430*(-w1003*w430 - w1025)) + w786*w876 + w789*(-w1034 - w333*(w1032 + w333*w967) - w429*(w1003*w429 + w1033)) + w794*w880);
    const scalar w1051 = w248*w325;
    const scalar w1052 = w251*w327;
    const scalar w1053 = gil*w1051 + git*w1052 + w964;
    const scalar w1054 = gel*w1051 + get*w1052 + w1053 + w966;
    const scalar w1055 = w233*w346;
    const scalar w1056 = w236*w348;
    const scalar w1057 = gil*w1055 + git*w1056 + w1000;
    const scalar w1058 = gel*w1055 + get*w1056 + w1002 + w1057;
    const scalar w1059 = w1049 + w311*(w1043 + w261*(w1031 + w333*(w1029 + w132*(-w1028 - w1054)) + w429*(w1030 + w140*(-w1028 - w1058))) + w274*(-w1021 - w334*(w132*(-w1054 - w956) + w995) - w430*(w1013 + w140*(-w1058 - w956))) + w306*(w1027 + w334*(-w1023 - w1054*w334) + w430*(-w1025 - w1058*w430)) + w788*w875 + w789*(-w1034 - w333*(w1032 + w1054*w333) - w429*(w1033 + w1058*w429)) + w796*w879);
    const scalar w1060 = sigmaEps*w180*w431;
    const scalar w1061 = w1060*w684;
    const scalar w1062 = w198*w456;
    const scalar w1063 = w204*w457;
    const scalar w1064 = w959*y;
    const scalar w1065 = 8*w1064;
    const scalar w1066 = -w1065*w200 + w1065*w73;
    const scalar w1067 = w1065*w53 + w1065*w77;
    const scalar w1068 = w1066*w75 + w1067*w79;
    const scalar w1069 = gil*w1062 + git*w1063 + w1068;
    const scalar w1070 = w1066*w83 + w1067*w85;
    const scalar w1071 = gel*w1062 + get*w1063 + w1069 + w1070;
    const scalar w1072 = 4*w1064;
    const scalar w1073 = -w1072*w200 + w1072*w73;
    const scalar w1074 = w1073*w57;
    const scalar w1075 = w1064*w245 + w1072*w77;
    const scalar w1076 = w1075*w62;
    const scalar w1077 = w198*w902;
    const scalar w1078 = w423*w713;
    const scalar w1079 = w204*y;
    const scalar w1080 = w424*z;
    const scalar w1081 = -gil*w1074 - git*w1076 + w1072*w76 + w1072*w80 - w1077*w980 - w1078*w980 + w1079*w981 + w1080*w981;
    const scalar w1082 = -gel*w1074 - get*w1076 + w1072*w84 + w1072*w86 - w1077*w973 - w1078*w973 + w1079*w977 + w1080*w977 + w1081;
    const scalar w1083 = -w1072*w52 - w1072*w54;
    const scalar w1084 = w217*w423;
    const scalar w1085 = w198*w320;
    const scalar w1086 = w1072*w59 - w1072*w726;
    const scalar w1087 = w221*w424;
    const scalar w1088 = w204*w322;
    const scalar w1089 = gil*w1084 + gil*w1085 + git*w1087 + git*w1088 + w1073*w56 + w1075*w61 + w1083*w75 + w1086*w79;
    const scalar w1090 = gel*w1084 + gel*w1085 + get*w1087 + get*w1088 + w1073*w65 + w1075*w67 + w1083*w83 + w1086*w85 + w1089;
    const scalar w1091 = -w1090;
    const scalar w1092 = w1082*w30 + w1091*w140 - w128*w737 - w144*w896 + w190*w441 + w211*w212 + w212*w453 - w242*w445;
    const scalar w1093 = w1092 + w269*w354;
    const scalar w1094 = w217*w426;
    const scalar w1095 = w221*w427;
    const scalar w1096 = -w1065*w52 - w1065*w54;
    const scalar w1097 = w1065*w59 - w1065*w726;
    const scalar w1098 = w1096*w56 + w1097*w61;
    const scalar w1099 = gil*w1094 + git*w1095 + w1098;
    const scalar w1100 = w1096*w65 + w1097*w67;
    const scalar w1101 = gel*w1094 + get*w1095 + w1099 + w1100;
    const scalar w1102 = w1083*w57;
    const scalar w1103 = w1086*w62;
    const scalar w1104 = w217*w902;
    const scalar w1105 = w320*w713;
    const scalar w1106 = w221*y;
    const scalar w1107 = w322*z;
    const scalar w1108 = -gil*w1102 - git*w1103 + w1072*w58 + w1072*w63 - w1104*w980 - w1105*w980 + w1106*w981 + w1107*w981;
    const scalar w1109 = -gel*w1102 - get*w1103 + w1072*w66 + w1072*w68 - w1104*w973 - w1105*w973 + w1106*w977 + w1107*w977 + w1108;
    const scalar w1110 = w1091*w132 + w1109*w30 + w190*w444 + w736 + w738 + w993;
    const scalar w1111 = w1110 - w169*w896 + w212*w467 - w267*w445;
    const scalar w1112 = w1065*w43;
    const scalar w1113 = w1065*w44;
    const scalar w1114 = w1065*w39;
    const scalar w1115 = w1065*w41;
    const scalar w1116 = w1065*w40;
    const scalar w1117 = gel*w1116;
    const scalar w1118 = w1065*w38;
    const scalar w1119 = get*w1118;
    const scalar w1120 = gil*w1116;
    const scalar w1121 = git*w1118;
    const scalar w1122 = -w1112 - w1113 - w1114 - w1115 + w1117 + w1119 + w1120 + w1121;
    const scalar w1123 = w1082*w132 + w1109*w140 + w211*w354 - w228*w445 + w930 + w931;
    const scalar w1124 = w154*(w1123 - w167*w737 + w190*w465 + w212*w264 + w30*(-w1061 - w1122)) + w268 + w466;
    const scalar w1125 = -w1082;
    const scalar w1126 = w1090*w430 + w1125*w154 + w297 + w487;
    const scalar w1127 = -w1109;
    const scalar w1128 = w1090*w334 + w1127*w154 + w296 + w486;
    const scalar w1129 = w1114 + w1115 - w1120 - w1121;
    const scalar w1130 = -w1112 - w1113 + w1117 + w1119 - w1129;
    const scalar w1131 = w154*(-w1125*w334 - w1127*w430 - w1130*w154 - w292 - w485) + w298 + w488;
    const scalar w1132 = w1060*w778;
    const scalar w1133 = w1092 + w254*w354;
    const scalar w1134 = w1110 - w143*w896 + w212*w450 - w239*w445;
    const scalar w1135 = w244 + w33*(w1123 - w127*w737 + w187*w212 + w190*w436 + w30*(-w1122 - w1132)) + w447;
    const scalar w1136 = w1090*w429 + w1125*w33 + w283 + w474;
    const scalar w1137 = w1090*w333 + w1127*w33 + w282 + w473;
    const scalar w1138 = w284 + w33*(w1125*w333 + w1127*w429 + w1130*w33 + w278 + w472) + w475;
    const scalar w1139 = w257*w260*w484 + w257*w946 + w260*w290*w463 + w271*w273*w949 - w271*w945 + w273*w471*w790 + w290*w951 + w463*w783 - w471*w781 + w484*w793 - w790*w950 - w791*w949;
    const scalar w1140 = w311*(w1139 + w261*(w1135 + w333*(w1133 + w132*(-w1071 - w1132)) + w429*(w1134 + w140*(-w1101 - w1132))) + w274*(-w1124 - w334*(w1093 + w132*(-w1061 - w1071)) - w430*(w1111 + w140*(-w1061 - w1101))) + w306*(w1131 + w334*(-w1071*w334 - w1126) + w430*(-w1101*w430 - w1128)) + w786*w948 + w789*(-w1138 - w333*(w1071*w333 + w1136) - w429*(w1101*w429 + w1137)) + w794*w953);
    const scalar w1141 = w307*y;
    const scalar w1142 = w497*z;
    const scalar w1143 = w1044*w1141 + w1044*w1142 + w1046*w1141 + w1046*w1142 + w1048*y;
    const scalar w1144 = w1143 + w212;
    const scalar w1145 = w248*w423;
    const scalar w1146 = w251*w424;
    const scalar w1147 = gil*w1145 + git*w1146 + w1068;
    const scalar w1148 = gel*w1145 + get*w1146 + w1070 + w1147;
    const scalar w1149 = w233*w320;
    const scalar w1150 = w236*w322;
    const scalar w1151 = gil*w1149 + git*w1150 + w1098;
    const scalar w1152 = gel*w1149 + get*w1150 + w1100 + w1151;
    const scalar w1153 = w311*(w1139 + w261*(w1135 + w333*(w1133 + w132*(-w1132 - w1148)) + w429*(w1134 + w140*(-w1132 - w1152))) + w274*(-w1124 - w334*(w1093 + w132*(-w1061 - w1148)) - w430*(w1111 + w140*(-w1061 - w1152))) + w306*(w1131 + w334*(-w1126 - w1148*w334) + w430*(-w1128 - w1152*w430)) + w788*w947 + w789*(-w1138 - w333*(w1136 + w1148*w333) - w429*(w1137 + w1152*w429)) + w796*w952);
    const scalar w1154 = w599*w680;
    const scalar w1155 = -w500*w842*x - w705;
    const scalar w1156 = 2*w553;
    const scalar w1157 = -w551*w81;
    const scalar w1158 = w64*w737;
    const scalar w1159 = -w771;
    const scalar w1160 = w206*w554;
    const scalar w1161 = w212*w223;
    const scalar w1162 = -w720;
    const scalar w1163 = -w733;
    const scalar w1164 = w1155*w42 + w1156*w276 + w1157 + w1158 + w1159*w513 + 2*w1160 - 2*w1161 + w1162*w518 + w1163*w524;
    const scalar w1165 = 2*w554;
    const scalar w1166 = w1155*w64 + w1156*w223 + w1163*w513 + w1165*w240 - w237*w748 + w518*w755 + w524*w746 - w551*w99 + w737*w94;
    const scalar w1167 = -w107*w551 + w1155*w81 + w1156*w206 + w1162*w513 + w1165*w252 - w240*w748 + w518*w765 + w524*w755 + w737*w99;
    const scalar w1168 = w595*w768;
    const scalar w1169 = w544*w780;
    const scalar w1170 = w35*w512;
    const scalar w1171 = w1170*w33;
    const scalar w1172 = w580*w784;
    const scalar w1173 = 2*w537 + 2*w586 + 2*w590 + 2*w593;
    const scalar w1174 = w273*w548/(w543*w543*w543);
    const scalar w1175 = w1174*w594;
    const scalar w1176 = w548*w595;
    const scalar w1177 = w1176*w780;
    const scalar w1178 = w35*w581;
    const scalar w1179 = w1178*w33;
    const scalar w1180 = -2*w505 - 2*w571 - 2*w575 - 2*w578;
    const scalar w1181 = w260*w531/(w511*w511*w511);
    const scalar w1182 = w1181*w579;
    const scalar w1183 = w532*w797 - w549*w798;
    const scalar w1184 = w598*z;
    const scalar w1185 = w601 + w602;
    const scalar w1186 = w335*w432;
    const scalar w1187 = w1186*w683;
    const scalar w1188 = w330*w423;
    const scalar w1189 = w331*w424;
    const scalar w1190 = x*y;
    const scalar w1191 = w1190*w689;
    const scalar w1192 = 8*w1191;
    const scalar w1193 = -w1192*w200 + w1192*w73;
    const scalar w1194 = w1192*w53 + w1192*w77;
    const scalar w1195 = w1193*w75 + w1194*w79;
    const scalar w1196 = gil*w1188 + git*w1189 + w1195;
    const scalar w1197 = w1193*w83 + w1194*w85;
    const scalar w1198 = gel*w1188 + get*w1189 + w1196 + w1197;
    const scalar w1199 = 4*w1191;
    const scalar w1200 = -w1199*w200 + w1199*w73;
    const scalar w1201 = w1200*w57;
    const scalar w1202 = w1191*w245 + w1199*w77;
    const scalar w1203 = w1202*w62;
    const scalar w1204 = w423*w826;
    const scalar w1205 = w325*w902;
    const scalar w1206 = w424*x;
    const scalar w1207 = w327*y;
    const scalar w1208 = -gil*w1201 - git*w1203 + w1199*w76 + w1199*w80 - w1204*w980 - w1205*w980 + w1206*w981 + w1207*w981;
    const scalar w1209 = -gel*w1201 - get*w1203 + w1199*w84 + w1199*w86 - w1204*w973 - w1205*w973 + w1206*w977 + w1207*w977 + w1208;
    const scalar w1210 = -w1199*w52 - w1199*w54;
    const scalar w1211 = w320*w325;
    const scalar w1212 = w346*w423;
    const scalar w1213 = w1199*w59 - w1199*w726;
    const scalar w1214 = w322*w327;
    const scalar w1215 = w348*w424;
    const scalar w1216 = gil*w1211 + gil*w1212 + git*w1214 + git*w1215 + w1200*w56 + w1202*w61 + w1210*w75 + w1213*w79;
    const scalar w1217 = gel*w1211 + gel*w1212 + get*w1214 + get*w1215 + w1200*w65 + w1202*w67 + w1210*w83 + w1213*w85 + w1216;
    const scalar w1218 = -w1217;
    const scalar w1219 = w1209*w30 + w1218*w140 + w212*w345 + w354*w453 - w373*w445;
    const scalar w1220 = w1219 - w171*w840 + w354*w388 - w357*w469;
    const scalar w1221 = w320*w361;
    const scalar w1222 = w322*w363;
    const scalar w1223 = -w1192*w52 - w1192*w54;
    const scalar w1224 = w1192*w59 - w1192*w726;
    const scalar w1225 = w1223*w56 + w1224*w61;
    const scalar w1226 = gil*w1221 + git*w1222 + w1225;
    const scalar w1227 = w1223*w65 + w1224*w67;
    const scalar w1228 = gel*w1221 + get*w1222 + w1226 + w1227;
    const scalar w1229 = w1210*w57;
    const scalar w1230 = w1213*w62;
    const scalar w1231 = w320*w826;
    const scalar w1232 = w346*w902;
    const scalar w1233 = w322*x;
    const scalar w1234 = w348*y;
    const scalar w1235 = -gil*w1229 - git*w1230 + w1199*w58 + w1199*w63 - w1231*w980 - w1232*w980 + w1233*w981 + w1234*w981;
    const scalar w1236 = -gel*w1229 - get*w1230 + w1199*w66 + w1199*w68 - w1231*w973 - w1232*w973 + w1233*w977 + w1234*w977 + w1235;
    const scalar w1237 = w1019 + w1218*w132 + w1236*w30 + w354*w373 - w357*w453 + w864;
    const scalar w1238 = w1237 + w354*w467 - w386*w445;
    const scalar w1239 = w1192*w43;
    const scalar w1240 = w1192*w44;
    const scalar w1241 = w1192*w39;
    const scalar w1242 = w1192*w41;
    const scalar w1243 = w1192*w40;
    const scalar w1244 = gel*w1243;
    const scalar w1245 = w1192*w38;
    const scalar w1246 = get*w1245;
    const scalar w1247 = gil*w1243;
    const scalar w1248 = git*w1245;
    const scalar w1249 = -w1239 - w1240 - w1241 - w1242 + w1244 + w1246 + w1247 + w1248;
    const scalar w1250 = w1209*w132 + w1236*w140 - w128*w840 + w345*w354 - w352*w445 + w354*w444 - w357*w441;
    const scalar w1251 = w154*(w1250 + w212*w383 + w30*(-w1187 - w1249)) + w387 + w470;
    const scalar w1252 = -w1209;
    const scalar w1253 = w1217*w430 + w1252*w154 + w412 + w493;
    const scalar w1254 = -w1236;
    const scalar w1255 = w1217*w334 + w1254*w154 + w415 + w490;
    const scalar w1256 = w1241 + w1242 - w1247 - w1248;
    const scalar w1257 = -w1239 - w1240 + w1244 + w1246 - w1256;
    const scalar w1258 = w154*(-w1252*w334 - w1254*w430 - w1257*w154 - w408 - w487) + w416 + w494;
    const scalar w1259 = w1186*w777;
    const scalar w1260 = w1219 - w147*w840 + w354*w378 - w357*w460;
    const scalar w1261 = w1237 + w354*w450 - w367*w445;
    const scalar w1262 = w33*(w1250 + w212*w342 + w30*(-w1249 - w1259)) + w375 + w462;
    const scalar w1263 = w1217*w429 + w1252*w33 + w399 + w480;
    const scalar w1264 = w1217*w333 + w1254*w33 + w402 + w477;
    const scalar w1265 = w33*(w1252*w333 + w1254*w429 + w1257*w33 + w395 + w474) + w403 + w481;
    const scalar w1266 = w1035*w429 + w1036*w471 + w1037*w949 - w1038*w471 + w1039*w429 + w1040*w463 + w1041*w484 - w1042*w949 - w391*w944 + w463*w874 + w484*w878 - w877*w944;
    const scalar w1267 = w311*(w1266 + w261*(w1262 + w333*(w1260 + w132*(-w1198 - w1259)) + w429*(w1261 + w140*(-w1228 - w1259))) + w274*(-w1251 - w334*(w1220 + w132*(-w1187 - w1198)) - w430*(w1238 + w140*(-w1187 - w1228))) + w306*(w1258 + w334*(-w1198*w334 - w1253) + w430*(-w1228*w430 - w1255)) + w789*(-w1265 - w333*(w1198*w333 + w1263) - w429*(w1228*w429 + w1264)) + w876*w947 + w880*w952);
    const scalar w1268 = w497*x;
    const scalar w1269 = w420*y;
    const scalar w1270 = w1044*w1268 + w1044*w1269 + w1046*w1268 + w1046*w1269 + w1190*w681;
    const scalar w1271 = w1270 + w354;
    const scalar w1272 = w325*w456;
    const scalar w1273 = w327*w457;
    const scalar w1274 = gil*w1272 + git*w1273 + w1195;
    const scalar w1275 = gel*w1272 + get*w1273 + w1197 + w1274;
    const scalar w1276 = w346*w426;
    const scalar w1277 = w348*w427;
    const scalar w1278 = gil*w1276 + git*w1277 + w1225;
    const scalar w1279 = gel*w1276 + get*w1277 + w1227 + w1278;
    const scalar w1280 = w311*(w1266 + w261*(w1262 + w333*(w1260 + w132*(-w1259 - w1275)) + w429*(w1261 + w140*(-w1259 - w1279))) + w274*(-w1251 - w334*(w1220 + w132*(-w1187 - w1275)) - w430*(w1238 + w140*(-w1187 - w1279))) + w306*(w1258 + w334*(-w1253 - w1275*w334) + w430*(-w1255 - w1279*w430)) + w789*(-w1265 - w333*(w1263 + w1275*w333) - w429*(w1264 + w1279*w429)) + w875*w948 + w879*w953);
    const scalar w1281 = (3.0/5.0)*w522 + w840;
    const scalar w1282 = -w866;
    const scalar w1283 = 2*w604;
    const scalar w1284 = 2*w606;
    const scalar w1285 = -w830;
    const scalar w1286 = -w838;
    const scalar w1287 = w1165*w392 + w1281*w64 + w1282*w513 + w1283*w350 + w1284*w343 + w1285*w518 + w1286*w524 + w81*w843;
    const scalar w1288 = w1165*w350 + w1281*w94 + w1283*w365 + w1284*w370 + w1286*w513 + w518*w854 + w524*w850 + w843*w99;
    const scalar w1289 = w1281*w99;
    const scalar w1290 = w343*w554;
    const scalar w1291 = w370*w604;
    const scalar w1292 = w107*w843 + w1284*w376 + w1285*w513 + w1289 + 2*w1290 + 2*w1291 + w518*w862 + w524*w854;
    const scalar w1293 = w544*w872;
    const scalar w1294 = w1170*w333;
    const scalar w1295 = 2*w542 + 2*w633 + 2*w636 + 2*w639;
    const scalar w1296 = w1174*w640;
    const scalar w1297 = w1176*w640;
    const scalar w1298 = w1178*w333;
    const scalar w1299 = -2*w510 - 2*w621 - 2*w624 - 2*w627;
    const scalar w1300 = w1181*w628;
    const scalar w1301 = w641*x;
    const scalar w1302 = w356*w516;
    const scalar w1303 = St*hA1*w431*w519*w841 - w916;
    const scalar w1304 = -w940;
    const scalar w1305 = 2*w644;
    const scalar w1306 = -w906;
    const scalar w1307 = -w914;
    const scalar w1308 = w1283*w439 + w1302*w81 + w1303*w64 + w1304*w513 + w1305*w442 + w1306*w518 + w1307*w524 + w42*w896 - w434*w748;
    const scalar w1309 = w64*w896;
    const scalar w1310 = w1302*w99;
    const scalar w1311 = w212*w442;
    const scalar w1312 = w451*w604;
    const scalar w1313 = w1303*w94 + w1305*w448 + w1307*w513 + w1309 + w1310 - 2*w1311 + 2*w1312 + w518*w927 + w524*w923;
    const scalar w1314 = w107*w1302 + w1283*w458 + w1303*w99 + w1305*w451 + w1306*w513 - w439*w748 + w518*w937 + w524*w927 + w81*w896;
    const scalar w1315 = w544*w944;
    const scalar w1316 = w1170*w429;
    const scalar w1317 = 2*w540 + 2*w670 + 2*w673 + 2*w676;
    const scalar w1318 = w1174*w677;
    const scalar w1319 = w1176*w944;
    const scalar w1320 = w1178*w429;
    const scalar w1321 = -2*w508 - 2*w658 - 2*w661 - 2*w664;
    const scalar w1322 = w1181*w665;
    const scalar w1323 = w678*y;
    const scalar w1324 = -w982;
    const scalar w1325 = w240*w604;
    const scalar w1326 = w1157 + w1160 + w1324*w513 + w1325 - w212*w370 + w252*w606 + w343*w553 + w376*w554 + w524*w990;
    const scalar w1327 = w1326 + w518*w965;
    const scalar w1328 = -w1010;
    const scalar w1329 = w1328*w513 - w212*w365 + w223*w554 + w237*w604 + w240*w606 + w350*w553 + w370*w554 + w518*w990 - w551*w64;
    const scalar w1330 = w1001*w524 + w1329;
    const scalar w1331 = -w1017;
    const scalar w1332 = -w212*w350;
    const scalar w1333 = w1290 + w1324*w518 + w1328*w524 + w1331*w513 + w1332 + w206*w606 + w223*w604 + w276*w554 + w392*w553 - w42*w551;
    const scalar w1334 = w1333*w154 + w565 + w614;
    const scalar w1335 = w1324*w154 + w430*w990 + w591 + w631;
    const scalar w1336 = w1328*w154 + w334*w990 + w588 + w630;
    const scalar w1337 = w154*(-w1324*w334 - w1328*w430 - w1331*w154 - w584 - w629) + w592 + w632;
    const scalar w1338 = w1333*w33 + w561 + w608;
    const scalar w1339 = w1324*w33 + w429*w990 + w576 + w619;
    const scalar w1340 = w1328*w33 + w333*w990 + w573 + w618;
    const scalar w1341 = w33*(w1324*w333 + w1328*w429 + w1331*w33 + w569 + w617) + w577 + w620;
    const scalar w1342 = w1170*w613;
    const scalar w1343 = w596*w640;
    const scalar w1344 = w596*w616;
    const scalar w1345 = w544*w616;
    const scalar w1346 = w1178*w628;
    const scalar w1347 = w260*w580;
    const scalar w1348 = w1347*w628;
    const scalar w1349 = w1347*w613;
    const scalar w1350 = w1176*w872;
    const scalar w1351 = -w1293*w566 + w1294*w562 - w1297*w780 + w1298*w579 + w1342*w33 + w1343*w566 + w1344*w594 - w1345*w780 + w1346*w33 + w1348*w562 + w1349*w579 - w1350*w594;
    const scalar w1352 = w598*x;
    const scalar w1353 = w641*z;
    const scalar w1354 = w1154*z;
    const scalar w1355 = w1044*w1352 + w1044*w1353 + w1046*w1352 + w1046*w1353 + w1354*x + w554;
    const scalar w1356 = w1053*w518 + w1326;
    const scalar w1357 = w1057*w524 + w1329;
    const scalar w1358 = -w1081;
    const scalar w1359 = w1089*w524 + w1358*w513 - w206*w212 - w212*w451 + w240*w644 + w252*w604 + w439*w553 + w458*w554 + w737*w81 + w896*w99;
    const scalar w1360 = w1069*w518 + w1359;
    const scalar w1361 = -w1108;
    const scalar w1362 = w1089*w518 + w1158 - w1161 + w1325 + w1361*w513 - w212*w448 + w237*w644 + w442*w553 + w451*w554 + w896*w94;
    const scalar w1363 = w1099*w524 + w1362;
    const scalar w1364 = -w1129;
    const scalar w1365 = w439*w554;
    const scalar w1366 = w1309 - w1311 + w1358*w518 + w1361*w524 + w1364*w513 + w1365 + w206*w604 - w212*w276 + w223*w644 + w42*w737 + w434*w553;
    const scalar w1367 = w1366*w154 + w564 + w651;
    const scalar w1368 = w1089*w430 + w1358*w154 + w588 + w668;
    const scalar w1369 = w1089*w334 + w1361*w154 + w587 + w667;
    const scalar w1370 = w154*(-w1358*w334 - w1361*w430 - w1364*w154 - w583 - w666) + w589 + w669;
    const scalar w1371 = w1366*w33 + w559 + w646;
    const scalar w1372 = w1089*w429 + w1358*w33 + w573 + w656;
    const scalar w1373 = w1089*w333 + w1361*w33 + w572 + w655;
    const scalar w1374 = w33*(w1358*w333 + w1361*w429 + w1364*w33 + w568 + w654) + w574 + w657;
    const scalar w1375 = -w1169*w653 + w1171*w650 - w1177*w677 + w1179*w665 - w1315*w566 + w1316*w562 - w1319*w594 + w1320*w579 + w1347*w562*w665 + w1347*w579*w650 + w566*w596*w677 + w594*w596*w653;
    const scalar w1376 = w598*y;
    const scalar w1377 = w678*z;
    const scalar w1378 = w1044*w1376 + w1044*w1377 + w1046*w1376 + w1046*w1377 + w1143 + w1354*y;
    const scalar w1379 = w1147*w518 + w1359;
    const scalar w1380 = w1151*w524 + w1362;
    const scalar w1381 = -w1208;
    const scalar w1382 = w107*w1281 + w1216*w524 + w1310 + w1312 + w1365 + w1381*w513 - w212*w343 + w370*w644 + w376*w604 + w458*w606;
    const scalar w1383 = w1196*w518 + w1382;
    const scalar w1384 = -w1235;
    const scalar w1385 = w1216*w518 + w1289 + w1291 + w1302*w94 + w1332 + w1384*w513 + w365*w644 + w442*w554 + w448*w604 + w451*w606;
    const scalar w1386 = w1226*w524 + w1385;
    const scalar w1387 = -w1256;
    const scalar w1388 = w1281*w81 + w1302*w64 + w1381*w518 + w1384*w524 + w1387*w513 - w212*w392 + w343*w604 + w350*w644 + w434*w554 + w439*w606 + w442*w604;
    const scalar w1389 = w1388*w154 + w615 + w652;
    const scalar w1390 = w1216*w430 + w1381*w154 + w634 + w674;
    const scalar w1391 = w1216*w334 + w1384*w154 + w637 + w671;
    const scalar w1392 = w154*(-w1381*w334 - w1384*w430 - w1387*w154 - w630 - w668) + w638 + w675;
    const scalar w1393 = w1388*w33 + w610 + w649;
    const scalar w1394 = w1216*w429 + w1381*w33 + w622 + w662;
    const scalar w1395 = w1216*w333 + w1384*w33 + w625 + w659;
    const scalar w1396 = w33*(w1381*w333 + w1384*w429 + w1387*w33 + w618 + w656) + w626 + w663;
    const scalar w1397 = -w1293*w653 + w1294*w650 - w1297*w944 + w1298*w665 + w1342*w429 + w1343*w653 + w1344*w677 - w1345*w944 + w1346*w429 + w1348*w650 + w1349*w665 - w1350*w677;
    const scalar w1398 = w678*x;
    const scalar w1399 = w641*y;
    const scalar w1400 = (3.0/5.0)*St*w114*w521 + w1044*w1398 + w1044*w1399 + w1046*w1398 + w1046*w1399 + w1154*w1190 + w1270;
    const scalar w1401 = w1274*w518 + w1382;
    const scalar w1402 = w1278*w524 + w1385;
    return w100*(w1267 + w1271) + w100*(w1271 + w1280) + w1050*w87 + w1059*w87 + w107*(w1154*w7 + w1185 + w1301*w800 + w1301*w802 + w311*(w1168*w616*w640 + w1172*w613*w628 + w1183 - 2*w1293*w616 + 2*w1294*w613 + w1295*w1296 - w1297*w873 + 2*w1298*w628 + w1299*w1300 + w563*(-w1287*w33 - w1288*w429 - w1292*w333 - w612*w776) + w567*(w1287*w154 + w1288*w430 + w1292*w334 + w612*w682) + w582*(-w33*(w1282*w33 + w1285*w333 + w1286*w429 + 2*w619) - w333*(w1285*w33 + w333*w862 + w376*w776 + w429*w854) - w429*(w1286*w33 + w333*w854 + w429*w850 + 2*w622) - w623*w776) + w597*(w154*(-w1282*w154 - w1285*w334 - w1286*w430 - 2*w631) + w334*(-w1285*w154 - w334*w862 - w376*w682 - w430*w854) + w430*(-w1286*w154 - w334*w854 - w430*w850 - 2*w634) + w635*w682)) + w605 + w882) + w108*(-w357 + w882) + w206*w642 + w209*w422 + w223*w679 + w226*w499 + w23*w319 + w276*w603 + w319*w324 + w319*w329 + w332*w422 + w343*w603 + w370*w679 + w371*w499 + w376*w642 + w42*(w1154*w9 + w1184*w800 + w1184*w802 + w1185 + w311*(w1168*w566*w594 - 2*w1169*w566 + 2*w1171*w562 + w1172*w562*w579 + w1173*w1175 - 2*w1177*w594 + 2*w1179*w579 + w1180*w1182 + w1183 + w563*(-w1164*w33 - w1166*w429 - w1167*w333 - w557*w776) + w567*(w1164*w154 + w1166*w430 + w1167*w334 + w557*w682) + w582*(-w33*(w1159*w33 + w1162*w333 + w1163*w429 + w276*w776) - w333*(w1162*w33 + w333*w765 + w429*w755 + 2*w569) - w429*(w1163*w33 + w333*w755 + w429*w746 + 2*w568) - w570*w776) + w597*(w154*(-w1159*w154 - w1162*w334 - w1163*w430 - w276*w682) + w334*(-w1162*w154 - w334*w765 - w430*w755 - 2*w584) + w430*(-w1163*w154 - w334*w755 - w430*w746 - 2*w583) + w585*w682)) - w552 + w804) + w422*w425 + w428*w499 + w442*w603 + w448*w679 + w45*(w190 + w804) + w451*w642 + w64*(w1140 + w1378 + w311*(w1173*w1318 + w1180*w1322 + w1375 + w563*(-w1360*w333 - w1363*w429 - w1371) + w567*(w1360*w334 + w1363*w430 + w1367) + w582*(-w1374 - w333*(w1069*w333 + w1372) - w429*(w1099*w429 + w1373)) + w597*(w1370 + w334*(-w1069*w334 - w1368) + w430*(-w1099*w430 - w1369)))) + w64*(w1153 + w1378 + w311*(w1175*w1317 + w1182*w1321 + w1375 + w563*(-w1371 - w1379*w333 - w1380*w429) + w567*(w1367 + w1379*w334 + w1380*w430) + w582*(-w1374 - w333*(w1147*w333 + w1372) - w429*(w1151*w429 + w1373)) + w597*(w1370 + w334*(-w1147*w334 - w1368) + w430*(-w1151*w430 - w1369)))) + w69*(w1140 + w1144) + w69*(w1144 + w1153) + w81*(w1050 + w1355 + w311*(w1173*w1296 + w1180*w1300 + w1351 + w563*(-w1327*w333 - w1330*w429 - w1338) + w567*(w1327*w334 + w1330*w430 + w1334) + w582*(-w1341 - w333*(w1339 + w333*w965) - w429*(w1001*w429 + w1340)) + w597*(w1337 + w334*(-w1335 - w334*w965) + w430*(-w1001*w430 - w1336)))) + w81*(w1059 + w1355 + w311*(w1175*w1295 + w1182*w1299 + w1351 + w563*(-w1338 - w1356*w333 - w1357*w429) + w567*(w1334 + w1356*w334 + w1357*w430) + w582*(-w1341 - w333*(w1053*w333 + w1339) - w429*(w1057*w429 + w1340)) + w597*(w1337 + w334*(-w1053*w334 - w1335) + w430*(-w1057*w430 - w1336)))) + w94*(w1154*w8 + w1185 + w1323*w800 + w1323*w802 + w311*(w1168*w653*w677 + w1172*w650*w665 + w1183 - 2*w1315*w653 + 2*w1316*w650 + w1317*w1318 - 2*w1319*w677 + 2*w1320*w665 + w1321*w1322 + w563*(-w1308*w33 - w1313*w429 - w1314*w333 - w647*w776) + w567*(w1308*w154 + w1313*w430 + w1314*w334 + w647*w682) + w582*(-w33*(w1304*w33 + w1306*w333 + w1307*w429 + 2*w655) - w333*(w1306*w33 + w333*w937 + w429*w927 + 2*w659) - w429*(w1307*w33 + w333*w927 + w429*w923 + w448*w776) - w660*w776) + w597*(w154*(-w1304*w154 - w1306*w334 - w1307*w430 - 2*w667) + w334*(-w1306*w154 - w334*w937 - w430*w927 - 2*w671) + w430*(-w1307*w154 - w334*w927 - w430*w923 - w448*w682) + w672*w682)) + w643 + w955) + w95*(-w445 + w955) + w99*(w1267 + w1400 + w311*(w1296*w1317 + w1300*w1321 + w1397 + w563*(-w1383*w333 - w1386*w429 - w1393) + w567*(w1383*w334 + w1386*w430 + w1389) + w582*(-w1396 - w333*(w1196*w333 + w1394) - w429*(w1226*w429 + w1395)) + w597*(w1392 + w334*(-w1196*w334 - w1390) + w430*(-w1226*w430 - w1391)))) + w99*(w1280 + w1400 + w311*(w1295*w1318 + w1299*w1322 + w1397 + w563*(-w1393 - w1401*w333 - w1402*w429) + w567*(w1389 + w1401*w334 + w1402*w430) + w582*(-w1396 - w333*(w1274*w333 + w1394) - w429*(w1278*w429 + w1395)) + w597*(w1392 + w334*(-w1274*w334 - w1390) + w430*(-w1278*w430 - w1391))));
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

    const scalar w0 = sin(gamma1);
    const scalar w1 = gen*w0;
    const scalar w2 = cos(gamma2);
    const scalar w3 = sin(gamma2);
    const scalar w4 = w2*w3;
    const scalar w5 = cos(gamma1);
    const scalar w6 = (r1*r1);
    const scalar w7 = -(r2*r2);
    const scalar w8 = (1.0/((-w6 - w7)));
    const scalar w9 = (x*x) + (y*y) + (z*z);
    const scalar w10 = w8*(-w6 + w9);
    const scalar w11 = theta0 + theta1*w10;
    const scalar w12 = cos(w11);
    const scalar w13 = sin(w11);
    const scalar w14 = w0*w13;
    const scalar w15 = w14*w2;
    const scalar w16 = w12*w5 - w15;
    const scalar w17 = w13*w3;
    const scalar w18 = w16*w17;
    const scalar w19 = w12*w3;
    const scalar w20 = w13*w5;
    const scalar w21 = w0*w12;
    const scalar w22 = w2*w21;
    const scalar w23 = -w20 - w22;
    const scalar w24 = gel*w18 + get*w19*w23 + w1*w4;
    const scalar w25 = (1.0/5.0)*Tt;
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
    const scalar w36 = k*x;
    const scalar w37 = w36 + 1.0/5.0;
    const scalar w38 = sin(w37);
    const scalar w39 = k*z;
    const scalar w40 = (1.0/2.0)*sin(w39) + 1;
    const scalar w41 = cos(w26);
    const scalar w42 = sigmaEps*w40*w41;
    const scalar w43 = w38*w42 + 1;
    const scalar w44 = -sigma0Out*w43;
    const scalar w45 = (w2*w2);
    const scalar w46 = (w3*w3);
    const scalar w47 = (w13*w13)*w46;
    const scalar w48 = (w12*w12)*w46;
    const scalar w49 = gel*w47 + gen*w45 + get*w48;
    const scalar w50 = gil*w47 + gin*w45 + git*w48 + w49;
    const scalar w51 = -w44 - w50;
    const scalar w52 = gin*w0;
    const scalar w53 = git*w19;
    const scalar w54 = gil*w18 + w23*w53 + w24 + w4*w52;
    const scalar w55 = -w54;
    const scalar w56 = k*cos(w36);
    const scalar w57 = fA2*((1.0/2.0)*y + 1);
    const scalar w58 = St*(fA0 + w56*w57);
    const scalar w59 = w2*w20;
    const scalar w60 = w21 + w59;
    const scalar w61 = w17*w60;
    const scalar w62 = w12*w2*w5 - w14;
    const scalar w63 = w4*w5;
    const scalar w64 = get*w19;
    const scalar w65 = gel*w61 - gen*w63 + w62*w64;
    const scalar w66 = gil*w61 - gin*w63 + w53*w62 + w65;
    const scalar w67 = -w66;
    const scalar w68 = sin(w36);
    const scalar w69 = (1.0/2.0)*fA2;
    const scalar w70 = k*sin(w27);
    const scalar w71 = St*(fA1*z - w29*w70 + w68*w69);
    const scalar w72 = w55*w58 + w67*w71;
    const scalar w73 = w35*(w31*w51 + w72);
    const scalar w74 = (w60*w60);
    const scalar w75 = (w62*w62);
    const scalar w76 = w46*(w5*w5);
    const scalar w77 = gel*w74 + gen*w76 + get*w75;
    const scalar w78 = gil*w74 + gin*w76 + git*w75 + w77;
    const scalar w79 = -w44 - w78;
    const scalar w80 = w16*w60;
    const scalar w81 = w23*w62;
    const scalar w82 = w46*w5;
    const scalar w83 = gel*w80 + get*w81 - w1*w82;
    const scalar w84 = gil*w80 + git*w81 - w52*w82 + w83;
    const scalar w85 = -w84;
    const scalar w86 = w31*w67 + w58*w85;
    const scalar w87 = w35*(w71*w79 + w86);
    const scalar w88 = (w16*w16);
    const scalar w89 = (w23*w23);
    const scalar w90 = (w0*w0)*w46;
    const scalar w91 = gel*w88 + gen*w90 + get*w89;
    const scalar w92 = gil*w88 + gin*w90 + git*w89 + w91;
    const scalar w93 = -w44 - w92;
    const scalar w94 = w31*w55 + w71*w85;
    const scalar w95 = w35*(w58*w93 + w94);
    const scalar w96 = w73*z + w87*y + w95*x;
    const scalar w97 = w35*w50;
    const scalar w98 = w35*w66;
    const scalar w99 = w35*w54;
    const scalar w100 = w35*(w97*z + w98*y + w99*x);
    const scalar w101 = w35*w78;
    const scalar w102 = w35*w84;
    const scalar w103 = w35*(w101*y + w102*x + w98*z);
    const scalar w104 = w35*w92;
    const scalar w105 = w35*(w102*y + w104*x + w99*z);
    const scalar w106 = w100*z + w103*y + w105*x;
    const scalar w107 = w32/w106;
    const scalar w108 = w107*w96;
    const scalar w109 = (1.0/2.0)*w10*w35;
    const scalar w110 = (1.0/(r1));
    const scalar w111 = -sigma0In*w43;
    const scalar w112 = -w111 - w50;
    const scalar w113 = w110*(w112*w31 + w72);
    const scalar w114 = -w111 - w78;
    const scalar w115 = w110*(w114*w71 + w86);
    const scalar w116 = -w111 - w92;
    const scalar w117 = w110*(w116*w58 + w94);
    const scalar w118 = -w113*z - w115*y - w117*x;
    const scalar w119 = w110*w50;
    const scalar w120 = w110*w66;
    const scalar w121 = w110*w54;
    const scalar w122 = w110*(-w119*z - w120*y - w121*x);
    const scalar w123 = w110*w78;
    const scalar w124 = w110*w84;
    const scalar w125 = w110*(-w120*z - w123*y - w124*x);
    const scalar w126 = w110*w92;
    const scalar w127 = w110*(-w121*z - w124*y - w126*x);
    const scalar w128 = -w122*z - w125*y - w127*x;
    const scalar w129 = w32/w128;
    const scalar w130 = (1.0/2.0)*w110*(1 - w10);
    const scalar w131 = w129*w130;
    const scalar w132 = 2*w108*w109 + 2*w118*w131;
    const scalar w133 = w132*w34;
    const scalar w134 = w7 + w9;
    const scalar w135 = w132*w134*w33;
    const scalar w136 = w8*z;
    const scalar w137 = w35*z;
    const scalar w138 = w108*w8;
    const scalar w139 = w110*z;
    const scalar w140 = w118*w129*w8;
    const scalar w141 = k*sigma0Out;
    const scalar w142 = sigmaEps*w38;
    const scalar w143 = w141*w142;
    const scalar w144 = cos(w39);
    const scalar w145 = (1.0/2.0)*w144*w41;
    const scalar w146 = theta1*w136;
    const scalar w147 = 4*w146;
    const scalar w148 = w12*w13*w46;
    const scalar w149 = w147*w148;
    const scalar w150 = gel*w149;
    const scalar w151 = gil*w149;
    const scalar w152 = get*w149;
    const scalar w153 = git*w149;
    const scalar w154 = -w150 - w151 + w152 + w153;
    const scalar w155 = w28*w29;
    const scalar w156 = (4.0/25.0)*St*w155;
    const scalar w157 = 2*w146;
    const scalar w158 = w157*w19;
    const scalar w159 = w158*w16;
    const scalar w160 = -w157*w20 - w157*w22;
    const scalar w161 = w160*w17;
    const scalar w162 = w12*w5;
    const scalar w163 = w15*w157 - w157*w162;
    const scalar w164 = gel*w159 + gel*w161 - 2*get*theta1*w13*w23*w3*w8*z + gil*w159 + gil*w161 - 2*git*theta1*w13*w23*w3*w8*z + w163*w53 + w163*w64;
    const scalar w165 = -w164;
    const scalar w166 = St*(fA1 - w30*w70);
    const scalar w167 = w166*w67;
    const scalar w168 = w158*w60;
    const scalar w169 = 2*theta1*w12*w2*w5*w8*z - w14*w157;
    const scalar w170 = w169*w17;
    const scalar w171 = -w157*w21 - w157*w59;
    const scalar w172 = gel*w168 + gel*w170 - 2*get*theta1*w13*w3*w62*w8*z + gil*w168 + gil*w170 - 2*git*theta1*w13*w3*w62*w8*z + w171*w53 + w171*w64;
    const scalar w173 = -w172;
    const scalar w174 = w165*w58 + w167 + w173*w71;
    const scalar w175 = -1.0/2.0*k*sigma0Out*sigmaEps*w144*w38*w41;
    const scalar w176 = w60*(4*theta1*w12*w2*w5*w8*z - w14*w147);
    const scalar w177 = w62*(-w147*w21 - w147*w59);
    const scalar w178 = gel*w176 + get*w177 + gil*w176 + git*w177;
    const scalar w179 = w160*w60;
    const scalar w180 = w16*w169;
    const scalar w181 = w163*w62;
    const scalar w182 = w171*w23;
    const scalar w183 = gel*w179 + gel*w180 + get*w181 + get*w182 + gil*w179 + gil*w180 + git*w181 + git*w182;
    const scalar w184 = -w183;
    const scalar w185 = w156*w67 + w173*w31 + w184*w58;
    const scalar w186 = w35*y;
    const scalar w187 = w16*(-w147*w20 - w147*w22);
    const scalar w188 = w23*(w147*w15 - w147*w162);
    const scalar w189 = gel*w187 + get*w188 + gil*w187 + git*w188;
    const scalar w190 = w156*w55 + w165*w31 + w166*w85 + w184*w71;
    const scalar w191 = w35*x;
    const scalar w192 = w107*w109;
    const scalar w193 = k*sigma0In;
    const scalar w194 = w142*w193;
    const scalar w195 = -1.0/2.0*k*sigma0In*sigmaEps*w144*w38*w41;
    const scalar w196 = w110*y;
    const scalar w197 = w110*x;
    const scalar w198 = w150 + w151 - w152 - w153;
    const scalar w199 = w109*w32*w96/(w106*w106);
    const scalar w200 = w118*w130*w32/(w128*w128);
    const scalar w201 = w134*w34;
    const scalar w202 = w133*z + w135*w136 + w201*(w131*(-w113 - w139*(w112*w156 + w174 + w31*(w145*w194 + w154)) - w196*(w114*w166 + w185 + w71*(-w178 - w195)) - w197*(w190 + w58*(-w189 - w195))) + w137*w138 - w139*w140 + w192*(w137*(w156*w51 + w174 + w31*(w143*w145 + w154)) + w186*(w166*w79 + w185 + w71*(-w175 - w178)) + w191*(w190 + w58*(-w175 - w189)) + w73) + w199*(-w100 - w137*(w137*w198 + w164*w191 + w172*w186 + w97) - w186*(w137*w172 + w178*w186 + w183*w191 + w98) - w191*(w137*w164 + w183*w186 + w189*w191 + w99)) + w200*(w122 + w139*(-w119 - w139*w198 - w164*w197 - w172*w196) + w196*(-w120 - w139*w172 - w178*w196 - w183*w197) + w197*(-w121 - w139*w164 - w183*w196 - w189*w197))) + w25*z + w31;
    const scalar w203 = w8*x;
    const scalar w204 = cos(w37);
    const scalar w205 = w204*w42;
    const scalar w206 = theta1*w203;
    const scalar w207 = 4*w206;
    const scalar w208 = w148*w207;
    const scalar w209 = gel*w208;
    const scalar w210 = gil*w208;
    const scalar w211 = get*w208;
    const scalar w212 = git*w208;
    const scalar w213 = -w209 - w210 + w211 + w212;
    const scalar w214 = 2*w206;
    const scalar w215 = w19*w214;
    const scalar w216 = w16*w215;
    const scalar w217 = -w20*w214 - w214*w22;
    const scalar w218 = w17*w217;
    const scalar w219 = w15*w214 - w162*w214;
    const scalar w220 = gel*w216 + gel*w218 - 2*get*theta1*w13*w23*w3*w8*x + gil*w216 + gil*w218 - 2*git*theta1*w13*w23*w3*w8*x + w219*w53 + w219*w64;
    const scalar w221 = -w220;
    const scalar w222 = w215*w60;
    const scalar w223 = 2*theta1*w12*w2*w5*w8*x - w14*w214;
    const scalar w224 = w17*w223;
    const scalar w225 = -w21*w214 - w214*w59;
    const scalar w226 = gel*w222 + gel*w224 - 2*get*theta1*w13*w3*w62*w8*x + gil*w222 + gil*w224 - 2*git*theta1*w13*w3*w62*w8*x + w225*w53 + w225*w64;
    const scalar w227 = -w226;
    const scalar w228 = w56*w69;
    const scalar w229 = St*w228;
    const scalar w230 = (k*k);
    const scalar w231 = w230*w57*w68;
    const scalar w232 = St*w231;
    const scalar w233 = w221*w58 + w227*w71 + w229*w67 - w232*w55;
    const scalar w234 = -k*sigma0Out*sigmaEps*w204*w40*w41;
    const scalar w235 = w60*(4*theta1*w12*w2*w5*w8*x - w14*w207);
    const scalar w236 = w62*(-w207*w21 - w207*w59);
    const scalar w237 = gel*w235 + get*w236 + gil*w235 + git*w236;
    const scalar w238 = w217*w60;
    const scalar w239 = w16*w223;
    const scalar w240 = w219*w62;
    const scalar w241 = w225*w23;
    const scalar w242 = gel*w238 + gel*w239 + get*w240 + get*w241 + gil*w238 + gil*w239 + git*w240 + git*w241;
    const scalar w243 = -w242;
    const scalar w244 = St*w85;
    const scalar w245 = w227*w31 - w231*w244 + w243*w58;
    const scalar w246 = w16*(-w20*w207 - w207*w22);
    const scalar w247 = w23*(w15*w207 - w162*w207);
    const scalar w248 = gel*w246 + get*w247 + gil*w246 + git*w247;
    const scalar w249 = w228*w244;
    const scalar w250 = w221*w31 + w243*w71 + w249;
    const scalar w251 = -k*sigma0In*sigmaEps*w204*w40*w41;
    const scalar w252 = St*w116;
    const scalar w253 = w209 + w210 - w211 - w212;
    const scalar w254 = -k*w25*w68 + w133*x + w135*w203 + w201*(w131*(-w117 - w139*(w233 + w31*(w193*w205 + w213)) - w196*(w114*w229 + w245 + w71*(-w237 - w251)) - w197*(-w231*w252 + w250 + w58*(-w248 - w251))) + w138*w191 - w140*w197 + w192*(w137*(w233 + w31*(w141*w205 + w213)) + w186*(w229*w79 + w245 + w71*(-w234 - w237)) + w191*(-w232*w93 + w250 + w58*(-w234 - w248)) + w95) + w199*(-w105 - w137*(w137*w253 + w186*w226 + w191*w220 + w99) - w186*(w102 + w137*w226 + w186*w237 + w191*w242) - w191*(w104 + w137*w220 + w186*w242 + w191*w248)) + w200*(w127 + w139*(-w121 - w139*w253 - w196*w226 - w197*w220) + w196*(-w124 - w139*w226 - w196*w237 - w197*w242) + w197*(-w126 - w139*w220 - w196*w242 - w197*w248))) + w58;
    const scalar w255 = w8*y;
    const scalar w256 = w40*sin(w26);
    const scalar w257 = w143*w256;
    const scalar w258 = theta1*w255;
    const scalar w259 = 4*w258;
    const scalar w260 = w148*w259;
    const scalar w261 = gel*w260 - get*w260 + gil*w260 - git*w260;
    const scalar w262 = 2*w258;
    const scalar w263 = w19*w262;
    const scalar w264 = w16*w263;
    const scalar w265 = -w20*w262 - w22*w262;
    const scalar w266 = w17*w265;
    const scalar w267 = w15*w262 - w162*w262;
    const scalar w268 = gel*w264 + gel*w266 - 2*get*theta1*w13*w23*w3*w8*y + gil*w264 + gil*w266 - 2*git*theta1*w13*w23*w3*w8*y + w267*w53 + w267*w64;
    const scalar w269 = -w268;
    const scalar w270 = w263*w60;
    const scalar w271 = 2*theta1*w12*w2*w5*w8*y - w14*w262;
    const scalar w272 = w17*w271;
    const scalar w273 = -w21*w262 - w262*w59;
    const scalar w274 = gel*w270 + gel*w272 - 2*get*theta1*w13*w3*w62*w8*y + gil*w270 + gil*w272 - 2*git*theta1*w13*w3*w62*w8*y + w273*w53 + w273*w64;
    const scalar w275 = -w274;
    const scalar w276 = w155*w230;
    const scalar w277 = St*w276;
    const scalar w278 = w229*w55 + w269*w58 + w275*w71 - w277*w67;
    const scalar w279 = w60*(4*theta1*w12*w2*w5*w8*y - w14*w259);
    const scalar w280 = w62*(-w21*w259 - w259*w59);
    const scalar w281 = gel*w279 + get*w280 + gil*w279 + git*w280;
    const scalar w282 = w265*w60;
    const scalar w283 = w16*w271;
    const scalar w284 = w267*w62;
    const scalar w285 = w23*w273;
    const scalar w286 = gel*w282 + gel*w283 + get*w284 + get*w285 + gil*w282 + gil*w283 + git*w284 + git*w285;
    const scalar w287 = -w286;
    const scalar w288 = w167 + w249 + w275*w31 + w287*w58;
    const scalar w289 = w16*(-w20*w259 - w22*w259);
    const scalar w290 = w23*(w15*w259 - w162*w259);
    const scalar w291 = gel*w289 + get*w290 + gil*w289 + git*w290;
    const scalar w292 = w166*w55 - w244*w276 + w269*w31 + w287*w71;
    const scalar w293 = w194*w256;
    const scalar w294 = (1.0/10.0)*Tt + w133*y + w135*w255 + w201*(w131*(-w115 - w139*(w112*w166 + w278 + w31*(-w261 - w293)) - w196*(-w114*w277 + w288 + w71*(-w281 - w293)) - w197*(w228*w252 + w292 + w58*(-w291 - w293))) + w138*w186 - w140*w196 + w192*(w137*(w166*w51 + w278 + w31*(-w257 - w261)) + w186*(-w277*w79 + w288 + w71*(-w257 - w281)) + w191*(w229*w93 + w292 + w58*(-w257 - w291)) + w87) + w199*(-w103 - w137*(w137*w261 + w186*w274 + w191*w268 + w98) - w186*(w101 + w137*w274 + w186*w281 + w191*w286) - w191*(w102 + w137*w268 + w186*w286 + w191*w291)) + w200*(w125 + w139*(-w120 - w139*w261 - w196*w274 - w197*w268) + w196*(-w123 - w139*w274 - w196*w281 - w197*w286) + w197*(-w124 - w139*w268 - w196*w286 - w197*w291))) + w71;
    return n.x()*(w202*w24 + w254*w91 + w294*w83) + n.y()*(w202*w65 + w254*w83 + w294*w77) + n.z()*(w202*w49 + w24*w254 + w294*w65);
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
    {0, 0, 1, 1, 0, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, -0.40157830615360224},
    {0, 0, 1, 1, 0, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0.13700000000000001, -0.42820375944347355},
    {0, 0, 1, 1, 0, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, -0.38328676616230672},
    {0, 0, 1, 1, 0, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0.13700000000000001, -0.40869945338345581},
    {0, 0, 1, 1, 0, -0.12, 0.64000000000000001, 0.70999999999999996, 0, 0.030214210695758515},
    {0, 0, 1, 1, 0, -0.12, 0.64000000000000001, 0.70999999999999996, 0.13700000000000001, 0.032217473927967435},
    {0, 1, 1, 1, 0, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, -0.28207538960781398},
    {0, 1, 1, 1, 0, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0.13700000000000001, -0.30077755801467115},
    {0, 1, 1, 1, 0, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, -0.39695250744998178},
    {0, 1, 1, 1, 0, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0.13700000000000001, -0.42327126093703937},
    {0, 1, 1, 1, 0, -0.12, 0.64000000000000001, 0.70999999999999996, 0, -0.026615297663650318},
    {0, 1, 1, 1, 0, -0.12, 0.64000000000000001, 0.70999999999999996, 0.13700000000000001, -0.028379945688408133},
    {1, 0, 1, 0, 0, 0, 0, 0, 0, 0.14330047336884089},
    {1, 0, 1, 0, 0, 0, 0, 0, 0.13700000000000001, 0.15280158436421207},
    {1, 0, 1, 0, 0, 0.070000000000000007, 0.13, -0.050000000000000003, 0, 0.22428334221294519},
    {1, 0, 1, 0, 0, 0.070000000000000007, 0.13, -0.050000000000000003, 0.13700000000000001, 0.23915378107948812},
    {1, 0, 1, 0, 0, 0.20999999999999999, -0.29999999999999999, 0.25, 0, 0.38911062314617195},
    {1, 0, 1, 0, 0, 0.20999999999999999, -0.29999999999999999, 0.25, 0.13700000000000001, 0.41490944385539713},
    {1, 0, 1, 0, 0, -0.28000000000000003, 0.23999999999999999, -0.33000000000000002, 0, -0.3800782337906668},
    {1, 0, 1, 0, 0, -0.28000000000000003, 0.23999999999999999, -0.33000000000000002, 0.13700000000000001, -0.40527818883111516},
    {1, 0, 1, 1, 0, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, 0.56053607746677603},
    {1, 0, 1, 1, 0, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0.13700000000000001, 0.59770075224921004},
    {1, 0, 1, 1, 0, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 0.46521585054762427},
    {1, 0, 1, 1, 0, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0.13700000000000001, 0.4960606016426341},
    {1, 0, 1, 1, 0, -0.12, 0.64000000000000001, 0.70999999999999996, 0, -0.2665929663506722},
    {1, 0, 1, 1, 0, -0.12, 0.64000000000000001, 0.70999999999999996, 0.13700000000000001, -0.28426861880552134},
    {1, 0, 1, 2, 0, 0.62, 0.57999999999999996, -0.58999999999999997, 0, 0.53526941704720399},
    {1, 0, 1, 2, 0, 0.62, 0.57999999999999996, -0.58999999999999997, 0.13700000000000001, 0.5707588611797656},
    {1, 0, 1, 2, 0, -1.1000000000000001, 0.40000000000000002, 0.82999999999999996, 0, -0.11231176414675761},
    {1, 0, 1, 2, 0, -1.1000000000000001, 0.40000000000000002, 0.82999999999999996, 0.13700000000000001, -0.11975826109235892},
    {1, 0, 1, 2, 0, 1.47, -1.4099999999999999, 1.3799999999999999, 0, -0.2089984284522364},
    {1, 0, 1, 2, 0, 1.47, -1.4099999999999999, 1.3799999999999999, 0.13700000000000001, -0.22285544664555251},
    {1, 1, 1, 0, 0, 0, 0, 0, 0, 0.19106729782512122},
    {1, 1, 1, 0, 0, 0, 0, 0, 0.13700000000000001, 0.22718853255285659},
    {1, 1, 1, 0, 0, 0.070000000000000007, 0.13, -0.050000000000000003, 0, 0.24606032173458534},
    {1, 1, 1, 0, 0, 0.070000000000000007, 0.13, -0.050000000000000003, 0.13700000000000001, 0.28762722880235386},
    {1, 1, 1, 0, 0, 0.20999999999999999, -0.29999999999999999, 0.25, 0, 0.39305281500064043},
    {1, 1, 1, 0, 0, 0.20999999999999999, -0.29999999999999999, 0.25, 0.13700000000000001, 0.43428199429034525},
    {1, 1, 1, 0, 0, -0.28000000000000003, 0.23999999999999999, -0.33000000000000002, 0, -0.30523892587092694},
    {1, 1, 1, 0, 0, -0.28000000000000003, 0.23999999999999999, -0.33000000000000002, 0.13700000000000001, -0.30641357511560441},
    {1, 1, 1, 1, 0, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, 0.50664392848004169},
    {1, 1, 1, 1, 0, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0.13700000000000001, 0.55400778437610265},
    {1, 1, 1, 1, 0, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 0.40125756172280685},
    {1, 1, 1, 1, 0, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0.13700000000000001, 0.42020116256981838},
    {1, 1, 1, 1, 0, -0.12, 0.64000000000000001, 0.70999999999999996, 0, -0.2613037586932872},
    {1, 1, 1, 1, 0, -0.12, 0.64000000000000001, 0.70999999999999996, 0.13700000000000001, -0.24340625915811884},
    {1, 1, 1, 2, 0, 0.62, 0.57999999999999996, -0.58999999999999997, 0, 0.40650375048370047},
    {1, 1, 1, 2, 0, 0.62, 0.57999999999999996, -0.58999999999999997, 0.13700000000000001, 0.43592103521066033},
    {1, 1, 1, 2, 0, -1.1000000000000001, 0.40000000000000002, 0.82999999999999996, 0, -0.058266090958165311},
    {1, 1, 1, 2, 0, -1.1000000000000001, 0.40000000000000002, 0.82999999999999996, 0.13700000000000001, -0.062801338521836461},
    {1, 1, 1, 2, 0, 1.47, -1.4099999999999999, 1.3799999999999999, 0, -0.4128912616840319},
    {1, 1, 1, 2, 0, 1.47, -1.4099999999999999, 1.3799999999999999, 0.13700000000000001, -0.57944112373787438},
    {2, 0, 1, 1, 0, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, 0.15895777131317385},
    {2, 0, 1, 1, 0, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0.13700000000000001, 0.16949699280573652},
    {2, 0, 1, 1, 0, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 0.081929084385317527},
    {2, 0, 1, 1, 0, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0.13700000000000001, 0.087361148259178331},
    {2, 0, 1, 1, 0, -0.12, 0.64000000000000001, 0.70999999999999996, 0, -0.2363787556549137},
    {2, 0, 1, 1, 0, -0.12, 0.64000000000000001, 0.70999999999999996, 0.13700000000000001, -0.25205114487755387},
    {2, 1, 1, 1, 0, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, 0.22456853887222766},
    {2, 1, 1, 1, 0, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0.13700000000000001, 0.2532302263614315},
    {2, 1, 1, 1, 0, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 0.0043050542728250455},
    {2, 1, 1, 1, 0, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0.13700000000000001, -0.0030700983672210191},
    {2, 1, 1, 1, 0, -0.12, 0.64000000000000001, 0.70999999999999996, 0, -0.28791905635693754},
    {2, 1, 1, 1, 0, -0.12, 0.64000000000000001, 0.70999999999999996, 0.13700000000000001, -0.27178620484652699},
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
    {3, 0, 1, 1, 0, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, 0.66960737097822465},
    {3, 0, 1, 1, 1, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, -0.45031346931242466},
    {3, 0, 1, 1, 2, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, 0.23867710579655918},
    {3, 0, 1, 1, 0, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0.13700000000000001, 0.71400369295413069},
    {3, 0, 1, 1, 1, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0.13700000000000001, -0.48017016241374927},
    {3, 0, 1, 1, 2, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0.13700000000000001, 0.25450188027856829},
    {3, 0, 1, 1, 0, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 0.5090940765514842},
    {3, 0, 1, 1, 1, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 0.56695031145887076},
    {3, 0, 1, 1, 2, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 0.25532034322493591},
    {3, 0, 1, 1, 0, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0.13700000000000001, 0.54284804270867759},
    {3, 0, 1, 1, 1, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0.13700000000000001, 0.60454026291818153},
    {3, 0, 1, 1, 2, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0.13700000000000001, 0.27224859798451756},
    {3, 0, 1, 1, 0, -0.12, 0.64000000000000001, 0.70999999999999996, 0, 1.8885056224295915},
    {3, 0, 1, 1, 1, -0.12, 0.64000000000000001, 0.70999999999999996, 0, -0.55568533482873494},
    {3, 0, 1, 1, 2, -0.12, 0.64000000000000001, 0.70999999999999996, 0, -0.13471499266133755},
    {3, 0, 1, 1, 0, -0.12, 0.64000000000000001, 0.70999999999999996, 0.13700000000000001, 2.0137173618765565},
    {3, 0, 1, 1, 1, -0.12, 0.64000000000000001, 0.70999999999999996, 0.13700000000000001, -0.59252839557089021},
    {3, 0, 1, 1, 2, -0.12, 0.64000000000000001, 0.70999999999999996, 0.13700000000000001, -0.14364686893449918},
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
    {3, 1, 1, 0, 0, 0, 0, 0, 0.13700000000000001, 1.438275864549414},
    {3, 1, 1, 0, 1, 0, 0, 0, 0.13700000000000001, -0.18155044163697748},
    {3, 1, 1, 0, 2, 0, 0, 0, 0.13700000000000001, 0.081465127080082134},
    {3, 1, 1, 0, 0, 0.070000000000000007, 0.13, -0.050000000000000003, 0, 1.3928264995911162},
    {3, 1, 1, 0, 1, 0.070000000000000007, 0.13, -0.050000000000000003, 0, -0.3750298429511672},
    {3, 1, 1, 0, 2, 0.070000000000000007, 0.13, -0.050000000000000003, 0, 0.092048971401609642},
    {3, 1, 1, 0, 0, 0.070000000000000007, 0.13, -0.050000000000000003, 0.13700000000000001, 1.4726150080717464},
    {3, 1, 1, 0, 1, 0.070000000000000007, 0.13, -0.050000000000000003, 0.13700000000000001, -0.38384720970944552},
    {3, 1, 1, 0, 2, 0.070000000000000007, 0.13, -0.050000000000000003, 0.13700000000000001, 0.094708130161579321},
    {3, 1, 1, 0, 0, 0.20999999999999999, -0.29999999999999999, 0.25, 0, 0.98849769166779466},
    {3, 1, 1, 0, 1, 0.20999999999999999, -0.29999999999999999, 0.25, 0, 0.5858296666426418},
    {3, 1, 1, 0, 2, 0.20999999999999999, -0.29999999999999999, 0.25, 0, -0.004214813384448019},
    {3, 1, 1, 0, 0, 0.20999999999999999, -0.29999999999999999, 0.25, 0.13700000000000001, 1.0120976271250841},
    {3, 1, 1, 0, 1, 0.20999999999999999, -0.29999999999999999, 0.25, 0.13700000000000001, 0.63120619063466932},
    {3, 1, 1, 0, 2, 0.20999999999999999, -0.29999999999999999, 0.25, 0.13700000000000001, 0.0065672306670670512},
    {3, 1, 1, 0, 0, -0.28000000000000003, 0.23999999999999999, -0.33000000000000002, 0, 1.0349903971792205},
    {3, 1, 1, 0, 1, -0.28000000000000003, 0.23999999999999999, -0.33000000000000002, 0, -0.69604874607356804},
    {3, 1, 1, 0, 2, -0.28000000000000003, 0.23999999999999999, -0.33000000000000002, 0, 0.094640906118872978},
    {3, 1, 1, 0, 0, -0.28000000000000003, 0.23999999999999999, -0.33000000000000002, 0.13700000000000001, 1.1626816756856369},
    {3, 1, 1, 0, 1, -0.28000000000000003, 0.23999999999999999, -0.33000000000000002, 0.13700000000000001, -0.73258493373469902},
    {3, 1, 1, 0, 2, -0.28000000000000003, 0.23999999999999999, -0.33000000000000002, 0.13700000000000001, 0.096120733571538122},
    {3, 1, 1, 1, 0, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, 0.52583512845376235},
    {3, 1, 1, 1, 1, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, -0.6298006093052172},
    {3, 1, 1, 1, 2, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, 0.27067015331365835},
    {3, 1, 1, 1, 0, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0.13700000000000001, 0.49403136225445721},
    {3, 1, 1, 1, 1, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0.13700000000000001, -0.65983111916557269},
    {3, 1, 1, 1, 2, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0.13700000000000001, 0.28322192155537124},
    {3, 1, 1, 1, 0, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 0.64383836401715366},
    {3, 1, 1, 1, 1, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 0.69225546401845417},
    {3, 1, 1, 1, 2, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 0.31558871419774759},
    {3, 1, 1, 1, 0, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0.13700000000000001, 0.61375322738381188},
    {3, 1, 1, 1, 1, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0.13700000000000001, 0.74987994370179711},
    {3, 1, 1, 1, 2, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0.13700000000000001, 0.34495599497973067},
    {3, 1, 1, 1, 0, -0.12, 0.64000000000000001, 0.70999999999999996, 0, 1.6342951052726546},
    {3, 1, 1, 1, 1, -0.12, 0.64000000000000001, 0.70999999999999996, 0, -0.60468187969956211},
    {3, 1, 1, 1, 2, -0.12, 0.64000000000000001, 0.70999999999999996, 0, -0.060907178348402778},
    {3, 1, 1, 1, 0, -0.12, 0.64000000000000001, 0.70999999999999996, 0.13700000000000001, 1.7697756072740281},
    {3, 1, 1, 1, 1, -0.12, 0.64000000000000001, 0.70999999999999996, 0.13700000000000001, -0.63304696702197683},
    {3, 1, 1, 1, 2, -0.12, 0.64000000000000001, 0.70999999999999996, 0.13700000000000001, -0.048293755785558438},
    {3, 1, 1, 2, 0, 0.62, 0.57999999999999996, -0.58999999999999997, 0, -0.27215857172957603},
    {3, 1, 1, 2, 1, 0.62, 0.57999999999999996, -0.58999999999999997, 0, -0.40749595567690133},
    {3, 1, 1, 2, 2, 0.62, 0.57999999999999996, -0.58999999999999997, 0, 0.11190386684206198},
    {3, 1, 1, 2, 0, 0.62, 0.57999999999999996, -0.58999999999999997, 0.13700000000000001, -0.3544668175061732},
    {3, 1, 1, 2, 1, 0.62, 0.57999999999999996, -0.58999999999999997, 0.13700000000000001, -0.42032193880788393},
    {3, 1, 1, 2, 2, 0.62, 0.57999999999999996, -0.58999999999999997, 0.13700000000000001, 0.10186342000378648},
    {3, 1, 1, 2, 0, -1.1000000000000001, 0.40000000000000002, 0.82999999999999996, 0, -1.0048895091914534},
    {3, 1, 1, 2, 1, -1.1000000000000001, 0.40000000000000002, 0.82999999999999996, 0, -0.61405498506541112},
    {3, 1, 1, 2, 2, -1.1000000000000001, 0.40000000000000002, 0.82999999999999996, 0, 0.1015787085617427},
    {3, 1, 1, 2, 0, -1.1000000000000001, 0.40000000000000002, 0.82999999999999996, 0.13700000000000001, -1.0127007316485963},
    {3, 1, 1, 2, 1, -1.1000000000000001, 0.40000000000000002, 0.82999999999999996, 0.13700000000000001, -0.61829044070026051},
    {3, 1, 1, 2, 2, -1.1000000000000001, 0.40000000000000002, 0.82999999999999996, 0.13700000000000001, 0.16922725170046721},
    {3, 1, 1, 2, 0, 1.47, -1.4099999999999999, 1.3799999999999999, 0, 0.21947413896995249},
    {3, 1, 1, 2, 1, 1.47, -1.4099999999999999, 1.3799999999999999, 0, -0.74034015194130232},
    {3, 1, 1, 2, 2, 1.47, -1.4099999999999999, 1.3799999999999999, 0, -0.42895979553550556},
    {3, 1, 1, 2, 0, 1.47, -1.4099999999999999, 1.3799999999999999, 0.13700000000000001, 0.12581340778232719},
    {3, 1, 1, 2, 1, 1.47, -1.4099999999999999, 1.3799999999999999, 0.13700000000000001, -0.58946194902399751},
    {3, 1, 1, 2, 2, 1.47, -1.4099999999999999, 1.3799999999999999, 0.13700000000000001, -0.39427529420984697},
    {4, 0, 1, 1, 0, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, -0.26976825153107703},
    {4, 0, 1, 1, 1, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, 0.23760497373473496},
    {4, 0, 1, 1, 2, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, -0.52499783849759163},
    {4, 0, 1, 1, 0, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0.13700000000000001, -0.2876544318106552},
    {4, 0, 1, 1, 1, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0.13700000000000001, 0.25335866369426058},
    {4, 0, 1, 1, 2, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0.13700000000000001, -0.55980625621339919},
    {4, 0, 1, 1, 0, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 0.21135458452409836},
    {4, 0, 1, 1, 1, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, -0.0097433039785737672},
    {4, 0, 1, 1, 2, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 0.21261016051042034},
    {4, 0, 1, 1, 0, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0.13700000000000001, 0.22536782062678284},
    {4, 0, 1, 1, 1, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0.13700000000000001, -0.010389304723622247},
    {4, 0, 1, 1, 2, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0.13700000000000001, 0.22670664383852371},
    {4, 0, 1, 1, 0, -0.12, 0.64000000000000001, 0.70999999999999996, 0, -1.7601714567981577},
    {4, 0, 1, 1, 1, -0.12, 0.64000000000000001, 0.70999999999999996, 0, -0.069931016343391966},
    {4, 0, 1, 1, 2, -0.12, 0.64000000000000001, 0.70999999999999996, 0, 0.31911572658136172},
    {4, 0, 1, 1, 0, -0.12, 0.64000000000000001, 0.70999999999999996, 0.13700000000000001, -1.8768743816997291},
    {4, 0, 1, 1, 1, -0.12, 0.64000000000000001, 0.70999999999999996, 0.13700000000000001, -0.074567584057913944},
    {4, 0, 1, 1, 2, -0.12, 0.64000000000000001, 0.70999999999999996, 0.13700000000000001, 0.34027374418828266},
    {4, 1, 1, 1, 0, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, -0.25297839695816865},
    {4, 1, 1, 1, 1, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, 0.15940453483905678},
    {4, 1, 1, 1, 2, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, -0.5697348901917012},
    {4, 1, 1, 1, 0, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0.13700000000000001, -0.2697513759471778},
    {4, 1, 1, 1, 1, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0.13700000000000001, 0.16997337765629697},
    {4, 1, 1, 1, 2, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0.13700000000000001, -0.6075094648486471},
    {4, 1, 1, 1, 0, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 0.076278639322644212},
    {4, 1, 1, 1, 1, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 0.031108048067524763},
    {4, 1, 1, 1, 2, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 0.058751561196862423},
    {4, 1, 1, 1, 0, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0.13700000000000001, 0.081336067269270326},
    {4, 1, 1, 1, 1, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0.13700000000000001, 0.033170574523931858},
    {4, 1, 1, 1, 2, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0.13700000000000001, 0.062646908441431323},
    {4, 1, 1, 1, 0, -0.12, 0.64000000000000001, 0.70999999999999996, 0, -1.5912761303373866},
    {4, 1, 1, 1, 1, -0.12, 0.64000000000000001, 0.70999999999999996, 0, -0.099827480129052279},
    {4, 1, 1, 1, 2, -0.12, 0.64000000000000001, 0.70999999999999996, 0, 0.23822375020724054},
    {4, 1, 1, 1, 0, -0.12, 0.64000000000000001, 0.70999999999999996, 0.13700000000000001, -1.6967809537562579},
    {4, 1, 1, 1, 1, -0.12, 0.64000000000000001, 0.70999999999999996, 0.13700000000000001, -0.10644624381347548},
    {4, 1, 1, 1, 2, -0.12, 0.64000000000000001, 0.70999999999999996, 0.13700000000000001, 0.25401846629744379},
    {5, 0, 1, 1, 0, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, 0.75696412039062055},
    {5, 0, 1, 1, 0, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0.13700000000000001, 0.88905897953730773},
    {5, 0, 1, 1, 0, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 0.59309563712166591},
    {5, 0, 1, 1, 0, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0.13700000000000001, 0.70144733915042112},
    {5, 0, 1, 1, 0, -0.12, 0.64000000000000001, 0.70999999999999996, 0, 1.0273467167505934},
    {5, 0, 1, 1, 0, -0.12, 0.64000000000000001, 0.70999999999999996, 0.13700000000000001, 1.1659571333123149},
    {5, 1, 1, 1, 0, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, 0.87646703693640882},
    {5, 1, 1, 1, 0, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0.13700000000000001, 1.0164851809661102},
    {5, 1, 1, 1, 0, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 0.57942989583399085},
    {5, 1, 1, 1, 0, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0.13700000000000001, 0.68687553159683756},
    {5, 1, 1, 1, 0, -0.12, 0.64000000000000001, 0.70999999999999996, 0, 0.9705172083911846},
    {5, 1, 1, 1, 0, -0.12, 0.64000000000000001, 0.70999999999999996, 0.13700000000000001, 1.1053597136959392},
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
    {10, 0, 1, 1, 0, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, -1.6699373995644406},
    {10, 0, 1, 1, 0, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0.13700000000000001, -1.7806576241079082},
    {10, 0, 1, 1, 0, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, -0.046218604809468281},
    {10, 0, 1, 1, 0, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0.13700000000000001, -0.049282991716381577},
    {10, 0, 1, 1, 0, -0.12, 0.64000000000000001, 0.70999999999999996, 0, 0.63866877882174622},
    {10, 0, 1, 1, 0, -0.12, 0.64000000000000001, 0.70999999999999996, 0.13700000000000001, 0.68101380961061875},
    {10, 1, 1, 1, 0, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, -1.1938167970429197},
    {10, 1, 1, 1, 0, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0.13700000000000001, -1.2908257280497368},
    {10, 1, 1, 1, 0, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, -0.13845766264163695},
    {10, 1, 1, 1, 0, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0.13700000000000001, -0.12665128067246037},
    {10, 1, 1, 1, 0, -0.12, 0.64000000000000001, 0.70999999999999996, 0, 0.89012444722757955},
    {10, 1, 1, 1, 0, -0.12, 0.64000000000000001, 0.70999999999999996, 0.13700000000000001, 0.94023192112699727},
    {11, 0, 1, 0, 0, 0, 0, 0, 0, -0.2907185550267159},
    {11, 0, 1, 0, 0, 0, 0, 0, 0.13700000000000001, -0.30999378276872941},
    {11, 0, 1, 0, 0, 0.070000000000000007, 0.13, -0.050000000000000003, 0, -0.79594228434634873},
    {11, 0, 1, 0, 0, 0.070000000000000007, 0.13, -0.050000000000000003, 0.13700000000000001, -0.84871486640209159},
    {11, 0, 1, 0, 0, 0.20999999999999999, -0.29999999999999999, 0.25, 0, -1.8741163078799035},
    {11, 0, 1, 0, 0, 0.20999999999999999, -0.29999999999999999, 0.25, 0.13700000000000001, -1.9983740066913431},
    {11, 0, 1, 0, 0, -0.28000000000000003, 0.23999999999999999, -0.33000000000000002, 0, 1.6677918766445881},
    {11, 0, 1, 0, 0, -0.28000000000000003, 0.23999999999999999, -0.33000000000000002, 0.13700000000000001, 1.7783698486823563},
    {11, 0, 1, 1, 0, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, -0.4130876421242608},
    {11, 0, 1, 1, 0, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0.13700000000000001, -0.44047618764941582},
    {11, 0, 1, 1, 0, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, -1.5471574581976846},
    {11, 0, 1, 1, 0, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0.13700000000000001, -1.6497371244896233},
    {11, 0, 1, 1, 0, -0.12, 0.64000000000000001, 0.70999999999999996, 0, -1.3400878126323703},
    {11, 0, 1, 1, 0, -0.12, 0.64000000000000001, 0.70999999999999996, 0.13700000000000001, -1.4289383429344762},
    {11, 0, 1, 2, 0, 0.62, 0.57999999999999996, -0.58999999999999997, 0, -0.99711047355113613},
    {11, 0, 1, 2, 0, 0.62, 0.57999999999999996, -0.58999999999999997, 0.13700000000000001, -1.0632209131131343},
    {11, 0, 1, 2, 0, -1.1000000000000001, 0.40000000000000002, 0.82999999999999996, 0, -0.29888761772040673},
    {11, 0, 1, 2, 0, -1.1000000000000001, 0.40000000000000002, 0.82999999999999996, 0.13700000000000001, -0.31870447082873099},
    {11, 0, 1, 2, 0, 1.47, -1.4099999999999999, 1.3799999999999999, 0, 0.61018548637174286},
    {11, 0, 1, 2, 0, 1.47, -1.4099999999999999, 1.3799999999999999, 0.13700000000000001, 0.65064201730629501},
    {11, 1, 1, 0, 0, 0, 0, 0, 0, -0.65584932635366522},
    {11, 1, 1, 0, 0, 0, 0, 0, 0.13700000000000001, -0.7926722344106214},
    {11, 1, 1, 0, 0, 0.070000000000000007, 0.13, -0.050000000000000003, 0, -1.0087208439057602},
    {11, 1, 1, 0, 0, 0.070000000000000007, 0.13, -0.050000000000000003, 0.13700000000000001, -1.2195098231242334},
    {11, 1, 1, 0, 0, 0.20999999999999999, -0.29999999999999999, 0.25, 0, -1.9887015354450022},
    {11, 1, 1, 0, 0, 0.20999999999999999, -0.29999999999999999, 0.25, 0.13700000000000001, -2.1476603351268682},
    {11, 1, 1, 0, 0, -0.28000000000000003, 0.23999999999999999, -0.33000000000000002, 0, 1.2746127219168226},
    {11, 1, 1, 0, 0, -0.28000000000000003, 0.23999999999999999, -0.33000000000000002, 0.13700000000000001, 1.2719088702434285},
    {11, 1, 1, 1, 0, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, 0.36964823017861548},
    {11, 1, 1, 1, 0, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0.13700000000000001, 0.31285564689219136},
    {11, 1, 1, 1, 0, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, -1.6892792253522444},
    {11, 1, 1, 1, 0, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0.13700000000000001, -1.7075222687877951},
    {11, 1, 1, 1, 0, -0.12, 0.64000000000000001, 0.70999999999999996, 0, -0.45254946290969061},
    {11, 1, 1, 1, 0, -0.12, 0.64000000000000001, 0.70999999999999996, 0.13700000000000001, -0.53186114694772513},
    {11, 1, 1, 2, 0, 0.62, 0.57999999999999996, -0.58999999999999997, 0, -0.79235412131931393},
    {11, 1, 1, 2, 0, 0.62, 0.57999999999999996, -0.58999999999999997, 0.13700000000000001, -0.82395744011753902},
    {11, 1, 1, 2, 0, -1.1000000000000001, 0.40000000000000002, 0.82999999999999996, 0, -0.25906848605643545},
    {11, 1, 1, 2, 0, -1.1000000000000001, 0.40000000000000002, 0.82999999999999996, 0.13700000000000001, -0.28337405276005823},
    {11, 1, 1, 2, 0, 1.47, -1.4099999999999999, 1.3799999999999999, 0, 0.68240521806754273},
    {11, 1, 1, 2, 0, 1.47, -1.4099999999999999, 1.3799999999999999, 0.13700000000000001, 1.0373137159745633},
    {12, 0, 1, 1, 0, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, 0.32617443661317841},
    {12, 0, 1, 1, 1, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, -0.0029645362254889234},
    {12, 0, 1, 1, 2, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, -0.05499708579620638},
    {12, 0, 1, 1, 0, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0.13700000000000001, 0.34780046096089917},
    {12, 0, 1, 1, 1, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0.13700000000000001, -0.0031610909685822804},
    {12, 0, 1, 1, 2, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0.13700000000000001, -0.058643503733897068},
    {12, 0, 1, 1, 0, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 0.34908314313972455},
    {12, 0, 1, 1, 1, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 0.29824636544524397},
    {12, 0, 1, 1, 2, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 0.17794608107838444},
    {12, 0, 1, 1, 0, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0.13700000000000001, 0.37222806102877282},
    {12, 0, 1, 1, 1, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0.13700000000000001, 0.31802070223175072},
    {12, 0, 1, 1, 2, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0.13700000000000001, 0.18974426588385537},
    {12, 0, 1, 1, 0, -0.12, 0.64000000000000001, 0.70999999999999996, 0, 0.21938226081060244},
    {12, 0, 1, 1, 1, -0.12, 0.64000000000000001, 0.70999999999999996, 0, -0.040087833132675409},
    {12, 0, 1, 1, 2, -0.12, 0.64000000000000001, 0.70999999999999996, 0, -0.033726754824853301},
    {12, 0, 1, 1, 0, -0.12, 0.64000000000000001, 0.70999999999999996, 0.13700000000000001, 0.23392774807505853},
    {12, 0, 1, 1, 1, -0.12, 0.64000000000000001, 0.70999999999999996, 0.13700000000000001, -0.042745737487095252},
    {12, 0, 1, 1, 2, -0.12, 0.64000000000000001, 0.70999999999999996, 0.13700000000000001, -0.035962906831691531},
    {12, 1, 1, 1, 0, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, 0.22712513707911991},
    {12, 1, 1, 1, 1, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, -0.031572208980683469},
    {12, 1, 1, 1, 2, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0, -0.032709811042959461},
    {12, 1, 1, 1, 0, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0.13700000000000001, 0.20515260817323677},
    {12, 1, 1, 1, 1, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0.13700000000000001, -0.036851345860603638},
    {12, 1, 1, 1, 2, 0.35999999999999999, 0.27000000000000002, -0.23000000000000001, 0.13700000000000001, -0.027008211789317157},
    {12, 1, 1, 1, 0, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 0.43512180955725049},
    {12, 1, 1, 1, 1, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 0.36953307877538549},
    {12, 1, 1, 1, 2, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0, 0.22069820277718361},
    {12, 1, 1, 1, 0, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0.13700000000000001, 0.44007710154986168},
    {12, 1, 1, 1, 1, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0.13700000000000001, 0.38228002788449866},
    {12, 1, 1, 1, 2, 0.55000000000000004, -0.46999999999999997, 0.35999999999999999, 0.13700000000000001, 0.22915923451358208},
    {12, 1, 1, 1, 0, -0.12, 0.64000000000000001, 0.70999999999999996, 0, 0.1823056617064365},
    {12, 1, 1, 1, 1, -0.12, 0.64000000000000001, 0.70999999999999996, 0, -0.06277685117161573},
    {12, 1, 1, 1, 2, -0.12, 0.64000000000000001, 0.70999999999999996, 0, -0.042018337984695682},
    {12, 1, 1, 1, 0, -0.12, 0.64000000000000001, 0.70999999999999996, 0.13700000000000001, 0.20093411683131879},
    {12, 1, 1, 1, 1, -0.12, 0.64000000000000001, 0.70999999999999996, 0.13700000000000001, -0.056130884120623434},
    {12, 1, 1, 1, 2, -0.12, 0.64000000000000001, 0.70999999999999996, 0.13700000000000001, -0.034200559934841897},
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
