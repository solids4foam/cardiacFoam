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
    bath-bidomain transmission manufactured solution on the oblique slab.

    Regenerate from generator/:
        python transmission_mms.py --write ../manufacturedTransmissionBathBidomainReference.H

    sympy 1.14.0
\*---------------------------------------------------------------------------*/

#ifndef manufacturedTransmissionBathBidomainReference_H
#define manufacturedTransmissionBathBidomainReference_H

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
namespace transmissionMMS
{

struct Params
{
    scalar alpha;
    scalar s1, s2, period;
    scalar gil, git, gin, gel, get, gen, theta0, theta1;
    scalar sigma0Low, sigma0High, sigmaEps;
    scalar fA[4], hA[3];
    scalar cmAmplitude, cmOmega;
    scalar P0, P1, C;
    label solution;
};

enum Zone : label { bathLow = 0, myocardium = 1, bathHigh = 2 };

inline Params defaultParams()
{
    Params p;
    p.alpha = degToRad(25.0);
    p.s1 = 0.5;
    p.s2 = 1.5;
    p.period = 1.0;
    p.gil = 0.17;
    p.git = 0.019;
    p.gin = 0.019;
    p.gel = 0.62;
    p.get = 0.12;
    p.gen = 0.12;
    p.theta0 = degToRad(-20.0);
    p.theta1 = degToRad(90.0);
    p.sigma0Low = 0.3;
    p.sigma0High = 0.45;
    p.sigmaEps = 0.2;
    p.cmAmplitude = 0.0;
    p.cmOmega = 0.0;
    p.P0 = 0.1;
    p.P1 = 0.4;
    p.C = 0.3;
    p.fA[0] = 0.3;
    p.fA[1] = 0.2;
    p.fA[2] = 0.4;
    p.fA[3] = 0.15;
    p.hA[0] = -0.2;
    p.hA[1] = 0.3;
    p.hA[2] = 0.1;
    p.solution = 1;
    return p;
}

inline Params alternateParams()
{
    Params p;
    p.alpha = degToRad(25.0);
    p.s1 = 0.5;
    p.s2 = 1.5;
    p.period = 1.0;
    p.gil = 0.17;
    p.git = 0.019;
    p.gin = 0.019;
    p.gel = 0.62;
    p.get = 0.12;
    p.gen = 0.12;
    p.theta0 = degToRad(-20.0);
    p.theta1 = degToRad(90.0);
    p.sigma0Low = 0.3;
    p.sigma0High = 0.45;
    p.sigmaEps = 0.2;
    p.cmAmplitude = 0.3;
    p.cmOmega = 20.0;
    p.P0 = 0.1;
    p.P1 = 0.4;
    p.C = 0.3;
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
    const scalar s = cos(p.alpha)*pt.x() + sin(p.alpha)*pt.y();
    if (s < p.s1)
    {
        return bathLow;
    }
    if (s > p.s2)
    {
        return bathHigh;
    }
    return myocardium;
}

inline void checkZone(const label zone, const bool heartAllowed)
{
    if (zone < bathLow || zone > bathHigh || (!heartAllowed && zone == myocardium))
    {
        FatalErrorInFunction
            << "Zone " << zone << " is not valid here" << exit(FatalError);
    }
}

inline scalar Vm(const point& pt, const scalar t, const Params& p)
{
    const scalar alpha = p.alpha;
    const scalar s1 = p.s1;
    const scalar s2 = p.s2;
    const scalar period = p.period;
    const scalar gil = p.gil;
    const scalar git = p.git;
    const scalar gel = p.gel;
    const scalar get = p.get;
    const scalar theta0 = p.theta0;
    const scalar theta1 = p.theta1;
    const scalar sigma0Low = p.sigma0Low;
    const scalar sigma0High = p.sigma0High;
    const scalar sigmaEps = p.sigmaEps;
    const scalar fA0 = p.fA[0];
    const scalar fA1 = p.fA[1];
    const scalar fA2 = p.fA[2];
    const scalar fA3 = p.fA[3];
    const scalar hA0 = p.hA[0];
    const scalar hA1 = p.hA[1];
    const scalar hA2 = p.hA[2];
    const scalar P0 = p.P0;
    const scalar P1 = p.P1;
    const scalar C = p.C;
    const scalar s = cos(alpha)*pt.x() + sin(alpha)*pt.y();
    const scalar tau = -sin(alpha)*pt.x() + cos(alpha)*pt.y();
    const scalar St = sqrt(1 + t);

    if (p.solution == 0)
    {
        return C*St - St*(P0 + P1*(s - s1));
    }

    const scalar w0 = (1.0/(period));
    const scalar w1 = constant::mathematical::pi*w0;
    const scalar w2 = tau*w1;
    const scalar w3 = 4*w2;
    const scalar w4 = hA2*sin(w3);
    const scalar w5 = (s*s);
    const scalar w6 = 1 - 3.0/10.0*w5;
    const scalar w7 = 2*w2;
    const scalar w8 = cos(w7);
    const scalar w9 = hA1*w8;
    const scalar w10 = (1.0/2.0)*s + 1;
    const scalar w11 = sin(w7);
    const scalar w12 = fA2*w11;
    const scalar w13 = w3 + 3.0/10.0;
    const scalar w14 = fA3*exp((2.0/5.0)*s);
    const scalar w15 = w14*cos(w13);
    const scalar w16 = -s2;
    const scalar w17 = (s - s1)/(-s1 - w16);
    const scalar w18 = -alpha + theta0 + theta1*w17;
    const scalar w19 = cos(w18);
    const scalar w20 = (w19*w19);
    const scalar w21 = sin(w18);
    const scalar w22 = (w21*w21);
    const scalar w23 = gil*w20 + git*w22;
    const scalar w24 = (1.0/(w23));
    const scalar w25 = 2*w1;
    const scalar w26 = St*w19*w21;
    const scalar w27 = St*w23*(hA0 - 3.0/5.0*s*w9 + w4) + w26*(gil - git)*(-hA1*w11*w25*w6 + 4*constant::mathematical::pi*hA2*s*w0*cos(w3));
    const scalar w28 = 1 - w17;
    const scalar w29 = w17*(s + w16);
    const scalar w30 = gel*w20 + get*w22;
    const scalar w31 = (1.0/(w30));
    const scalar w32 = w26*(gel - get)*(fA2*w10*w25*w8 - 4*w1*w14*sin(w13));
    const scalar w33 = 2*s;
    const scalar w34 = fA0 + fA1*w33 + (1.0/2.0)*w12 + (2.0/5.0)*w15;
    const scalar w35 = sigmaEps*w11*((1.0/2.0)*sin(constant::mathematical::pi*w33) + 1) + 1;
    return St*(hA0*s + s*w4 + w6*w9) - St*(fA0*s + fA1*w5 + w10*w12 + w15) + w29*(-w17*w24*w27 + w24*w27*w28) - w29*(w17*w31*(St*w34*(sigma0High*w35 - w30) - w32) + w28*w31*(-St*w34*(sigma0Low*w35 - w30) + w32));
}

inline scalar phiE(const point& pt, const scalar t, const label zone, const Params& p)
{
    const scalar alpha = p.alpha;
    const scalar s1 = p.s1;
    const scalar s2 = p.s2;
    const scalar period = p.period;
    const scalar gil = p.gil;
    const scalar git = p.git;
    const scalar gel = p.gel;
    const scalar get = p.get;
    const scalar theta0 = p.theta0;
    const scalar theta1 = p.theta1;
    const scalar sigma0Low = p.sigma0Low;
    const scalar sigma0High = p.sigma0High;
    const scalar sigmaEps = p.sigmaEps;
    const scalar fA0 = p.fA[0];
    const scalar fA1 = p.fA[1];
    const scalar fA2 = p.fA[2];
    const scalar fA3 = p.fA[3];
    const scalar cmAmplitude = p.cmAmplitude;
    const scalar cmOmega = p.cmOmega;
    const scalar P0 = p.P0;
    const scalar P1 = p.P1;
    const scalar s = cos(alpha)*pt.x() + sin(alpha)*pt.y();
    const scalar tau = -sin(alpha)*pt.x() + cos(alpha)*pt.y();
    const scalar St = sqrt(1 + t);
    const scalar Tt = cmAmplitude*sin(cmOmega*t);
    checkZone(zone, true);

    if (p.solution == 0)
    {
        if (zone == bathLow)
        {
            const scalar w0 = alpha - theta0;
            return St*(P0 + P1*(s - s1)*(gel*((cos(w0))*(cos(w0))) + get*((sin(w0))*(sin(w0))))/sigma0Low);
        }
        if (zone == bathHigh)
        {
            const scalar w0 = -s2;
            const scalar w1 = alpha - theta0;
            return St*(P0 + P1*(-s1 - w0) + P1*(s + w0)*(gel*((cos(w1))*(cos(w1))) + get*((sin(w1))*(sin(w1))))/sigma0High);
        }
        return St*(P0 + P1*(s - s1));
    }

    if (zone == bathLow)
    {
        const scalar w0 = 2*constant::mathematical::pi;
        const scalar w1 = (1.0/(period));
        const scalar w2 = tau*w1;
        const scalar w3 = w0*w2;
        const scalar w4 = sin(w3);
        const scalar w5 = alpha - theta0;
        const scalar w6 = cos(w5);
        const scalar w7 = (w6*w6);
        const scalar w8 = gel*w7;
        const scalar w9 = sin(w5);
        const scalar w10 = (w9*w9);
        const scalar w11 = get*w10;
        const scalar w12 = (2.0/5.0)*constant::mathematical::pi*w1*w4*w6*w9;
        const scalar w13 = gil*w7;
        const scalar w14 = git*w10;
        return St*(fA0*s + fA1*(s*s) + fA2*w4*((1.0/2.0)*s + 1) + fA3*exp((2.0/5.0)*s)*cos(4*constant::mathematical::pi*w2 + 3.0/10.0)) + Tt*((1.0/10.0)*s + (1.0/10.0 + (-1.0/10.0*w11 - w12*(gel - get) - 1.0/10.0*w8 - (w11 + w8)*(-w12*(gil - git) - 1.0/10.0*w13 - 1.0/10.0*w14)/(w13 + w14))/(sigma0Low*(sigmaEps*w4*((1.0/2.0)*sin(s*w0) + 1) + 1)))*(-s + s1) + (1.0/5.0)*cos(w3));
    }
    if (zone == bathHigh)
    {
        const scalar w0 = 2*constant::mathematical::pi;
        const scalar w1 = (1.0/(period));
        const scalar w2 = tau*w1;
        const scalar w3 = w0*w2;
        const scalar w4 = sin(w3);
        const scalar w5 = -alpha + theta0 + theta1;
        const scalar w6 = cos(w5);
        const scalar w7 = (w6*w6);
        const scalar w8 = gel*w7;
        const scalar w9 = sin(w5);
        const scalar w10 = (w9*w9);
        const scalar w11 = get*w10;
        const scalar w12 = (2.0/5.0)*constant::mathematical::pi*w1*w4*w6*w9;
        const scalar w13 = gil*w7;
        const scalar w14 = git*w10;
        return St*(fA0*s + fA1*(s*s) + fA2*w4*((1.0/2.0)*s + 1) + fA3*exp((2.0/5.0)*s)*cos(4*constant::mathematical::pi*w2 + 3.0/10.0)) + Tt*((1.0/10.0)*s + (-1.0/10.0 + ((1.0/10.0)*w11 - w12*(gel - get) + (1.0/10.0)*w8 - (w11 + w8)*(-w12*(gil - git) + (1.0/10.0)*w13 + (1.0/10.0)*w14)/(w13 + w14))/(sigma0High*(sigmaEps*w4*((1.0/2.0)*sin(s*w0) + 1) + 1)))*(s - s2) + (1.0/5.0)*cos(w3));
    }
    const scalar w0 = (1.0/2.0)*s + 1;
    const scalar w1 = constant::mathematical::pi/period;
    const scalar w2 = tau*w1;
    const scalar w3 = 2*w2;
    const scalar w4 = sin(w3);
    const scalar w5 = fA2*w4;
    const scalar w6 = 4*w2 + 3.0/10.0;
    const scalar w7 = fA3*exp((2.0/5.0)*s);
    const scalar w8 = w7*cos(w6);
    const scalar w9 = cos(w3);
    const scalar w10 = -s2;
    const scalar w11 = (s - s1)/(-s1 - w10);
    const scalar w12 = -alpha + theta0 + theta1*w11;
    const scalar w13 = cos(w12);
    const scalar w14 = (w13*w13);
    const scalar w15 = gil*w14;
    const scalar w16 = sin(w12);
    const scalar w17 = (w16*w16);
    const scalar w18 = git*w17;
    const scalar w19 = (1.0/((w15 + w18)));
    const scalar w20 = w13*w16;
    const scalar w21 = -2.0/5.0*w1*w20*w4*(gil - git) + (1.0/10.0)*w15 + (1.0/10.0)*w18;
    const scalar w22 = 1 - w11;
    const scalar w23 = w11*(s + w10);
    const scalar w24 = gel*w14 + get*w17;
    const scalar w25 = (1.0/(w24));
    const scalar w26 = St*w20*(gel - get)*(2*fA2*w0*w1*w9 - 4*w1*w7*sin(w6));
    const scalar w27 = 2*s;
    const scalar w28 = fA0 + fA1*w27 + (1.0/2.0)*w5 + (2.0/5.0)*w8;
    const scalar w29 = sigmaEps*w4*((1.0/2.0)*sin(constant::mathematical::pi*w27) + 1) + 1;
    return St*(fA0*s + fA1*(s*s) + w0*w5 + w8) + Tt*((1.0/10.0)*s + w23*(-w11*w19*w21 + w19*w21*w22) + (1.0/5.0)*w9) + w23*(w11*w25*(St*w28*(sigma0High*w29 - w24) - w26) + w22*w25*(-St*w28*(sigma0Low*w29 - w24) + w26));
}

inline scalar phiI(const point& pt, const scalar t, const Params& p)
{
    const scalar alpha = p.alpha;
    const scalar s1 = p.s1;
    const scalar s2 = p.s2;
    const scalar period = p.period;
    const scalar gil = p.gil;
    const scalar git = p.git;
    const scalar theta0 = p.theta0;
    const scalar theta1 = p.theta1;
    const scalar hA0 = p.hA[0];
    const scalar hA1 = p.hA[1];
    const scalar hA2 = p.hA[2];
    const scalar cmAmplitude = p.cmAmplitude;
    const scalar cmOmega = p.cmOmega;
    const scalar C = p.C;
    const scalar s = cos(alpha)*pt.x() + sin(alpha)*pt.y();
    const scalar tau = -sin(alpha)*pt.x() + cos(alpha)*pt.y();
    const scalar St = sqrt(1 + t);
    const scalar Tt = cmAmplitude*sin(cmOmega*t);

    if (p.solution == 0)
    {
        return C*St;
    }

    const scalar w0 = (1.0/(period));
    const scalar w1 = constant::mathematical::pi*w0;
    const scalar w2 = tau*w1;
    const scalar w3 = 4*w2;
    const scalar w4 = hA2*sin(w3);
    const scalar w5 = 1 - 3.0/10.0*(s*s);
    const scalar w6 = 2*w2;
    const scalar w7 = cos(w6);
    const scalar w8 = hA1*w7;
    const scalar w9 = -s2;
    const scalar w10 = (s - s1)/(-s1 - w9);
    const scalar w11 = -alpha + theta0 + theta1*w10;
    const scalar w12 = cos(w11);
    const scalar w13 = gil*(w12*w12);
    const scalar w14 = sin(w11);
    const scalar w15 = git*(w14*w14);
    const scalar w16 = w12*w14*(gil - git);
    const scalar w17 = w1*sin(w6);
    const scalar w18 = (1.0/10.0)*w13 + (1.0/10.0)*w15 - 2.0/5.0*w16*w17;
    const scalar w19 = w13 + w15;
    const scalar w20 = (1.0/(w19));
    const scalar w21 = w10*w20;
    const scalar w22 = w20*(1 - w10);
    const scalar w23 = w10*(s + w9);
    const scalar w24 = St*w16*(-2*hA1*w17*w5 + 4*constant::mathematical::pi*hA2*s*w0*cos(w3)) + St*w19*(hA0 - 3.0/5.0*s*w8 + w4);
    return St*(hA0*s + s*w4 + w5*w8) + Tt*((1.0/10.0)*s + w23*(-w18*w21 + w18*w22) + (1.0/5.0)*w7) + w23*(-w21*w24 + w22*w24);
}

inline vector gradPhiE(const point& pt, const scalar t, const label zone, const Params& p)
{
    const scalar alpha = p.alpha;
    const scalar s1 = p.s1;
    const scalar s2 = p.s2;
    const scalar period = p.period;
    const scalar gil = p.gil;
    const scalar git = p.git;
    const scalar gel = p.gel;
    const scalar get = p.get;
    const scalar theta0 = p.theta0;
    const scalar theta1 = p.theta1;
    const scalar sigma0Low = p.sigma0Low;
    const scalar sigma0High = p.sigma0High;
    const scalar sigmaEps = p.sigmaEps;
    const scalar fA0 = p.fA[0];
    const scalar fA1 = p.fA[1];
    const scalar fA2 = p.fA[2];
    const scalar fA3 = p.fA[3];
    const scalar cmAmplitude = p.cmAmplitude;
    const scalar cmOmega = p.cmOmega;
    const scalar P1 = p.P1;
    const scalar s = cos(alpha)*pt.x() + sin(alpha)*pt.y();
    const scalar tau = -sin(alpha)*pt.x() + cos(alpha)*pt.y();
    const scalar St = sqrt(1 + t);
    const scalar Tt = cmAmplitude*sin(cmOmega*t);
    checkZone(zone, true);

    if (p.solution == 0)
    {
        if (zone == bathLow)
        {
            const scalar w0 = alpha - theta0;
            const scalar w1 = P1*St*(gel*((cos(w0))*(cos(w0))) + get*((sin(w0))*(sin(w0))))/sigma0Low;
            return vector(w1*cos(alpha), w1*sin(alpha), 0);
        }
        if (zone == bathHigh)
        {
            const scalar w0 = alpha - theta0;
            const scalar w1 = P1*St*(gel*((cos(w0))*(cos(w0))) + get*((sin(w0))*(sin(w0))))/sigma0High;
            return vector(w1*cos(alpha), w1*sin(alpha), 0);
        }
        const scalar w0 = P1*St;
        return vector(w0*cos(alpha), w0*sin(alpha), 0);
    }

    if (zone == bathLow)
    {
        const scalar w0 = sin(alpha);
        const scalar w1 = constant::mathematical::pi/period;
        const scalar w2 = tau*w1;
        const scalar w3 = 2*w2;
        const scalar w4 = cos(w3);
        const scalar w5 = 2*w1*w4;
        const scalar w6 = 4*w2 + 3.0/10.0;
        const scalar w7 = fA3*exp((2.0/5.0)*s);
        const scalar w8 = sin(w3);
        const scalar w9 = (2.0/5.0)*w1*w8;
        const scalar w10 = -s + s1;
        const scalar w11 = gel - get;
        const scalar w12 = ((constant::mathematical::pi)*(constant::mathematical::pi));
        const scalar w13 = (1.0/(period*period));
        const scalar w14 = alpha - theta0;
        const scalar w15 = cos(w14);
        const scalar w16 = sin(w14);
        const scalar w17 = w15*w16;
        const scalar w18 = gil - git;
        const scalar w19 = (w15*w15);
        const scalar w20 = gel*w19;
        const scalar w21 = (w16*w16);
        const scalar w22 = get*w21;
        const scalar w23 = w20 + w22;
        const scalar w24 = gil*w19;
        const scalar w25 = git*w21;
        const scalar w26 = (1.0/((w24 + w25)));
        const scalar w27 = (1.0/(sigma0Low));
        const scalar w28 = 2*s;
        const scalar w29 = constant::mathematical::pi*w28;
        const scalar w30 = (1.0/2.0)*sin(w29) + 1;
        const scalar w31 = sigmaEps*w8;
        const scalar w32 = w30*w31 + 1;
        const scalar w33 = w27/w32;
        const scalar w34 = w17*w9;
        const scalar w35 = -w11*w34 - 1.0/10.0*w20 - 1.0/10.0*w22 - w23*w26*(-w18*w34 - 1.0/10.0*w24 - 1.0/10.0*w25);
        const scalar w36 = w27*w35/(w32*w32);
        const scalar w37 = St*(fA2*w5*((1.0/2.0)*s + 1) - 4*w1*w7*sin(w6)) + Tt*(w10*(-sigmaEps*w30*w36*w5 + w33*(-4.0/5.0*w11*w12*w13*w17*w4 + (4.0/5.0)*w12*w13*w15*w16*w18*w23*w26*w4)) - w9);
        const scalar w38 = cos(alpha);
        const scalar w39 = St*(fA0 + fA1*w28 + (1.0/2.0)*fA2*w8 + (2.0/5.0)*w7*cos(w6)) + Tt*(-constant::mathematical::pi*w10*w31*w36*cos(w29) - w33*w35);
        return vector(-w0*w37 + w38*w39, w0*w39 + w37*w38, 0);
    }
    if (zone == bathHigh)
    {
        const scalar w0 = sin(alpha);
        const scalar w1 = constant::mathematical::pi/period;
        const scalar w2 = tau*w1;
        const scalar w3 = 2*w2;
        const scalar w4 = cos(w3);
        const scalar w5 = 2*w1*w4;
        const scalar w6 = 4*w2 + 3.0/10.0;
        const scalar w7 = fA3*exp((2.0/5.0)*s);
        const scalar w8 = sin(w3);
        const scalar w9 = (2.0/5.0)*w1*w8;
        const scalar w10 = s - s2;
        const scalar w11 = gel - get;
        const scalar w12 = ((constant::mathematical::pi)*(constant::mathematical::pi));
        const scalar w13 = (1.0/(period*period));
        const scalar w14 = -alpha + theta0 + theta1;
        const scalar w15 = cos(w14);
        const scalar w16 = sin(w14);
        const scalar w17 = w15*w16;
        const scalar w18 = gil - git;
        const scalar w19 = (w15*w15);
        const scalar w20 = gel*w19;
        const scalar w21 = (w16*w16);
        const scalar w22 = get*w21;
        const scalar w23 = w20 + w22;
        const scalar w24 = gil*w19;
        const scalar w25 = git*w21;
        const scalar w26 = (1.0/((w24 + w25)));
        const scalar w27 = (1.0/(sigma0High));
        const scalar w28 = 2*s;
        const scalar w29 = constant::mathematical::pi*w28;
        const scalar w30 = (1.0/2.0)*sin(w29) + 1;
        const scalar w31 = sigmaEps*w8;
        const scalar w32 = w30*w31 + 1;
        const scalar w33 = (1.0/(w32));
        const scalar w34 = w17*w9;
        const scalar w35 = -w11*w34 + (1.0/10.0)*w20 + (1.0/10.0)*w22 - w23*w26*(-w18*w34 + (1.0/10.0)*w24 + (1.0/10.0)*w25);
        const scalar w36 = w27*w35/(w32*w32);
        const scalar w37 = St*(fA2*w5*((1.0/2.0)*s + 1) - 4*w1*w7*sin(w6)) + Tt*(w10*(-sigmaEps*w30*w36*w5 + w27*w33*(-4.0/5.0*w11*w12*w13*w17*w4 + (4.0/5.0)*w12*w13*w15*w16*w18*w23*w26*w4)) - w9);
        const scalar w38 = cos(alpha);
        const scalar w39 = St*(fA0 + fA1*w28 + (1.0/2.0)*fA2*w8 + (2.0/5.0)*w7*cos(w6)) + Tt*(-constant::mathematical::pi*w10*w31*w36*cos(w29) + w27*w33*w35);
        return vector(-w0*w37 + w38*w39, w0*w39 + w37*w38, 0);
    }
    const scalar w0 = sin(alpha);
    const scalar w1 = (1.0/2.0)*s + 1;
    const scalar w2 = 2*constant::mathematical::pi;
    const scalar w3 = (1.0/(period));
    const scalar w4 = tau*w3;
    const scalar w5 = w2*w4;
    const scalar w6 = cos(w5);
    const scalar w7 = fA2*w6;
    const scalar w8 = w2*w3;
    const scalar w9 = 4*constant::mathematical::pi;
    const scalar w10 = w4*w9 + 3.0/10.0;
    const scalar w11 = fA3*exp((2.0/5.0)*s);
    const scalar w12 = w11*sin(w10);
    const scalar w13 = St*(w1*w7*w8 - w12*w3*w9);
    const scalar w14 = sin(w5);
    const scalar w15 = constant::mathematical::pi*w3;
    const scalar w16 = (2.0/5.0)*w14*w15;
    const scalar w17 = -s2;
    const scalar w18 = (1.0/((-s1 - w17)));
    const scalar w19 = w18*(s - s1);
    const scalar w20 = -alpha + theta0 + theta1*w19;
    const scalar w21 = cos(w20);
    const scalar w22 = (w21*w21);
    const scalar w23 = gil*w22;
    const scalar w24 = sin(w20);
    const scalar w25 = (w24*w24);
    const scalar w26 = git*w25;
    const scalar w27 = w23 + w26;
    const scalar w28 = (1.0/(w27));
    const scalar w29 = w19*w28;
    const scalar w30 = gil - git;
    const scalar w31 = w21*w24;
    const scalar w32 = ((constant::mathematical::pi)*(constant::mathematical::pi))/(period*period);
    const scalar w33 = (4.0/5.0)*w30*w31*w32*w6;
    const scalar w34 = 1 - w19;
    const scalar w35 = w28*w34;
    const scalar w36 = s + w17;
    const scalar w37 = w19*w36;
    const scalar w38 = s*w2;
    const scalar w39 = (1.0/2.0)*sin(w38) + 1;
    const scalar w40 = 2*fA1;
    const scalar w41 = fA2*w14;
    const scalar w42 = w11*cos(w10);
    const scalar w43 = fA0 + s*w40 + (1.0/2.0)*w41 + (2.0/5.0)*w42;
    const scalar w44 = St*w43;
    const scalar w45 = sigmaEps*w39*w44*w6*w8;
    const scalar w46 = gel - get;
    const scalar w47 = w31*w46;
    const scalar w48 = -St*w47*(-4*w1*w32*w41 - 16*w32*w42);
    const scalar w49 = sigmaEps*w14;
    const scalar w50 = w39*w49 + 1;
    const scalar w51 = gel*w22 + get*w25;
    const scalar w52 = sigma0High*w50 - w51;
    const scalar w53 = St*(-8.0/5.0*w12*w15 + w15*w7);
    const scalar w54 = (1.0/(w51));
    const scalar w55 = w19*w54;
    const scalar w56 = sigma0Low*w50 - w51;
    const scalar w57 = w34*w54;
    const scalar w58 = Tt*(-w16 + w37*(w29*w33 - w33*w35)) + w13 + w37*(w55*(sigma0High*w45 + w48 + w52*w53) + w57*(-sigma0Low*w45 - w48 - w53*w56));
    const scalar w59 = cos(alpha);
    const scalar w60 = w13*w46;
    const scalar w61 = w31*w60;
    const scalar w62 = St*w43*w52 - w61;
    const scalar w63 = -w44*w56 + w61;
    const scalar w64 = w55*w62 + w57*w63;
    const scalar w65 = w18*w36;
    const scalar w66 = w16*w30;
    const scalar w67 = (1.0/10.0)*w23 + (1.0/10.0)*w26 - w31*w66;
    const scalar w68 = -w67;
    const scalar w69 = w29*w68 + w35*w67;
    const scalar w70 = w18*w28;
    const scalar w71 = theta1*w18;
    const scalar w72 = w31*w71;
    const scalar w73 = (1.0/5.0)*gil*w72 - 1.0/5.0*git*theta1*w18*w21*w24 - 2.0/5.0*constant::mathematical::pi*theta1*w14*w18*w25*w3*w30 + w22*w66*w71;
    const scalar w74 = 2*w72;
    const scalar w75 = (gil*w74 - git*w74)/(w27*w27);
    const scalar w76 = w18*w54;
    const scalar w77 = gel*w74 - get*w74;
    const scalar w78 = w77/(w51*w51);
    const scalar w79 = St*(w40 + (4.0/25.0)*w42);
    const scalar w80 = constant::mathematical::pi*w49*cos(w38);
    const scalar w81 = w60*w71;
    const scalar w82 = -w22*w81 + w25*w81 - w47*w53;
    const scalar w83 = Tt*(w19*w69 + w37*(w19*w68*w75 + w29*w73 + w34*w67*w75 - w35*w73 - w67*w70 + w68*w70) + w65*w69 + 1.0/10.0) + w19*w64 + w37*(w19*w62*w78 + w34*w63*w78 + w55*(w44*(sigma0High*w80 + w77) + w52*w79 + w82) + w57*(-w44*(sigma0Low*w80 + w77) - w56*w79 - w82) + w62*w76 - w63*w76) + w44 + w64*w65;
    return vector(-w0*w58 + w59*w83, w0*w83 + w58*w59, 0);
}

inline vector gradVm(const point& pt, const scalar t, const Params& p)
{
    const scalar alpha = p.alpha;
    const scalar s1 = p.s1;
    const scalar s2 = p.s2;
    const scalar period = p.period;
    const scalar gil = p.gil;
    const scalar git = p.git;
    const scalar gel = p.gel;
    const scalar get = p.get;
    const scalar theta0 = p.theta0;
    const scalar theta1 = p.theta1;
    const scalar sigma0Low = p.sigma0Low;
    const scalar sigma0High = p.sigma0High;
    const scalar sigmaEps = p.sigmaEps;
    const scalar fA0 = p.fA[0];
    const scalar fA1 = p.fA[1];
    const scalar fA2 = p.fA[2];
    const scalar fA3 = p.fA[3];
    const scalar hA0 = p.hA[0];
    const scalar hA1 = p.hA[1];
    const scalar hA2 = p.hA[2];
    const scalar P1 = p.P1;
    const scalar s = cos(alpha)*pt.x() + sin(alpha)*pt.y();
    const scalar tau = -sin(alpha)*pt.x() + cos(alpha)*pt.y();
    const scalar St = sqrt(1 + t);

    if (p.solution == 0)
    {
        const scalar w0 = P1*St;
        return vector(-w0*cos(alpha), -w0*sin(alpha), 0);
    }

    const scalar w0 = sin(alpha);
    const scalar w1 = (1.0/(period));
    const scalar w2 = tau*w1;
    const scalar w3 = 4*constant::mathematical::pi;
    const scalar w4 = w2*w3;
    const scalar w5 = cos(w4);
    const scalar w6 = 2*constant::mathematical::pi;
    const scalar w7 = w2*w6;
    const scalar w8 = sin(w7);
    const scalar w9 = w1*w6;
    const scalar w10 = hA1*(1 - 3.0/10.0*(s*s));
    const scalar w11 = 4*constant::mathematical::pi*hA2*s*w1*w5 - w10*w8*w9;
    const scalar w12 = (1.0/2.0)*s + 1;
    const scalar w13 = cos(w7);
    const scalar w14 = fA2*w13;
    const scalar w15 = w1*w3;
    const scalar w16 = w4 + 3.0/10.0;
    const scalar w17 = fA3*exp((2.0/5.0)*s);
    const scalar w18 = w17*sin(w16);
    const scalar w19 = St*(w12*w14*w9 - w15*w18);
    const scalar w20 = s - s1;
    const scalar w21 = -s2;
    const scalar w22 = s + w21;
    const scalar w23 = (1.0/((-s1 - w21)));
    const scalar w24 = hA1*s;
    const scalar w25 = constant::mathematical::pi*w1;
    const scalar w26 = hA2*w15*w5 + (6.0/5.0)*w24*w25*w8;
    const scalar w27 = w20*w23;
    const scalar w28 = -alpha + theta0 + theta1*w27;
    const scalar w29 = cos(w28);
    const scalar w30 = (w29*w29);
    const scalar w31 = sin(w28);
    const scalar w32 = (w31*w31);
    const scalar w33 = gil*w30 + git*w32;
    const scalar w34 = St*w33;
    const scalar w35 = hA2*sin(w4);
    const scalar w36 = ((constant::mathematical::pi)*(constant::mathematical::pi))/(period*period);
    const scalar w37 = 16*w36;
    const scalar w38 = 4*w36;
    const scalar w39 = gil - git;
    const scalar w40 = w29*w31;
    const scalar w41 = St*w39*w40;
    const scalar w42 = w26*w34 + w41*(-s*w35*w37 - w10*w13*w38);
    const scalar w43 = (1.0/(w33));
    const scalar w44 = w27*w43;
    const scalar w45 = 1 - w27;
    const scalar w46 = w43*w45;
    const scalar w47 = s*w6;
    const scalar w48 = (1.0/2.0)*sin(w47) + 1;
    const scalar w49 = 2*fA1;
    const scalar w50 = fA2*w8;
    const scalar w51 = w17*cos(w16);
    const scalar w52 = fA0 + s*w49 + (1.0/2.0)*w50 + (2.0/5.0)*w51;
    const scalar w53 = St*w52;
    const scalar w54 = sigmaEps*w13*w48*w53*w9;
    const scalar w55 = gel - get;
    const scalar w56 = w40*w55;
    const scalar w57 = -St*w56*(-w12*w38*w50 - w37*w51);
    const scalar w58 = sigmaEps*w8;
    const scalar w59 = w48*w58 + 1;
    const scalar w60 = gel*w30 + get*w32;
    const scalar w61 = sigma0High*w59 - w60;
    const scalar w62 = St*(w14*w25 - 8.0/5.0*w18*w25);
    const scalar w63 = (1.0/(w60));
    const scalar w64 = w27*w63;
    const scalar w65 = sigma0Low*w59 - w60;
    const scalar w66 = w45*w63;
    const scalar w67 = w22*w27;
    const scalar w68 = St*w11 - w19 + w20*w22*w23*(-w42*w44 + w42*w46) - w67*(w64*(sigma0High*w54 + w57 + w61*w62) + w66*(-sigma0Low*w54 - w57 - w62*w65));
    const scalar w69 = cos(alpha);
    const scalar w70 = (3.0/5.0)*w13;
    const scalar w71 = St*(hA0 - w24*w70 + w35);
    const scalar w72 = St*w11*w39;
    const scalar w73 = w33*w71 + w40*w72;
    const scalar w74 = -w73;
    const scalar w75 = w44*w74 + w46*w73;
    const scalar w76 = w22*w23;
    const scalar w77 = w19*w56;
    const scalar w78 = St*w52*w61 - w77;
    const scalar w79 = -w53*w65 + w77;
    const scalar w80 = w64*w78 + w66*w79;
    const scalar w81 = w23*w43;
    const scalar w82 = theta1*w23;
    const scalar w83 = 2*w40*w82;
    const scalar w84 = gil*w83 - 2*git*theta1*w23*w29*w31;
    const scalar w85 = w84/(w33*w33);
    const scalar w86 = w72*w82;
    const scalar w87 = -hA1*w34*w70 + w26*w41 + w30*w86 - w32*w86 - w71*w84;
    const scalar w88 = w23*w63;
    const scalar w89 = gel*w83 - get*w83;
    const scalar w90 = w89/(w60*w60);
    const scalar w91 = St*(w49 + (4.0/25.0)*w51);
    const scalar w92 = constant::mathematical::pi*w58*cos(w47);
    const scalar w93 = w19*w55*w82;
    const scalar w94 = -w30*w93 + w32*w93 - w56*w62;
    const scalar w95 = w27*w75 - w27*w80 - w53 + w67*(w27*w74*w85 - w44*w87 + w45*w73*w85 + w46*w87 - w73*w81 + w74*w81) - w67*(w27*w78*w90 + w45*w79*w90 + w64*(w53*(sigma0High*w92 + w89) + w61*w91 + w94) + w66*(-w53*(sigma0Low*w92 + w89) - w65*w91 - w94) + w78*w88 - w79*w88) + w71 + w75*w76 - w76*w80;
    return vector(-w0*w68 + w69*w95, w0*w95 + w68*w69, 0);
}

inline scalar u1(const point& pt, const scalar t, const Params& p)
{
    const scalar alpha = p.alpha;
    const scalar s1 = p.s1;
    const scalar s2 = p.s2;
    const scalar period = p.period;
    const scalar gil = p.gil;
    const scalar git = p.git;
    const scalar gel = p.gel;
    const scalar get = p.get;
    const scalar theta0 = p.theta0;
    const scalar theta1 = p.theta1;
    const scalar sigma0Low = p.sigma0Low;
    const scalar sigma0High = p.sigma0High;
    const scalar sigmaEps = p.sigmaEps;
    const scalar fA0 = p.fA[0];
    const scalar fA1 = p.fA[1];
    const scalar fA2 = p.fA[2];
    const scalar fA3 = p.fA[3];
    const scalar hA0 = p.hA[0];
    const scalar hA1 = p.hA[1];
    const scalar hA2 = p.hA[2];
    const scalar P0 = p.P0;
    const scalar P1 = p.P1;
    const scalar C = p.C;
    const scalar s = cos(alpha)*pt.x() + sin(alpha)*pt.y();
    const scalar tau = -sin(alpha)*pt.x() + cos(alpha)*pt.y();
    const scalar St = sqrt(1 + t);

    if (p.solution == 0)
    {
        return C*St - St*(P0 + P1*(s - s1)) + (t + 1)*((1.0/4.0)*sin(2*constant::mathematical::pi*tau/period)*cos(constant::mathematical::pi*s) + 1);
    }

    const scalar w0 = constant::mathematical::pi*s;
    const scalar w1 = (1.0/(period));
    const scalar w2 = constant::mathematical::pi*w1;
    const scalar w3 = tau*w2;
    const scalar w4 = 2*w3;
    const scalar w5 = sin(w4);
    const scalar w6 = 4*w3;
    const scalar w7 = hA2*sin(w6);
    const scalar w8 = (s*s);
    const scalar w9 = 1 - 3.0/10.0*w8;
    const scalar w10 = cos(w4);
    const scalar w11 = hA1*w10;
    const scalar w12 = (1.0/2.0)*s + 1;
    const scalar w13 = fA2*w5;
    const scalar w14 = w6 + 3.0/10.0;
    const scalar w15 = fA3*exp((2.0/5.0)*s);
    const scalar w16 = w15*cos(w14);
    const scalar w17 = -s2;
    const scalar w18 = (s - s1)/(-s1 - w17);
    const scalar w19 = -alpha + theta0 + theta1*w18;
    const scalar w20 = cos(w19);
    const scalar w21 = (w20*w20);
    const scalar w22 = sin(w19);
    const scalar w23 = (w22*w22);
    const scalar w24 = gil*w21 + git*w23;
    const scalar w25 = (1.0/(w24));
    const scalar w26 = 2*w2;
    const scalar w27 = St*w20*w22;
    const scalar w28 = St*w24*(hA0 - 3.0/5.0*s*w11 + w7) + w27*(gil - git)*(-hA1*w26*w5*w9 + 4*constant::mathematical::pi*hA2*s*w1*cos(w6));
    const scalar w29 = 1 - w18;
    const scalar w30 = w18*(s + w17);
    const scalar w31 = gel*w21 + get*w23;
    const scalar w32 = (1.0/(w31));
    const scalar w33 = w27*(gel - get)*(fA2*w10*w12*w26 - 4*w15*w2*sin(w14));
    const scalar w34 = fA0 + 2*fA1*s + (1.0/2.0)*w13 + (2.0/5.0)*w16;
    const scalar w35 = sigmaEps*w5*((1.0/2.0)*sin(2*w0) + 1) + 1;
    return St*(hA0*s + s*w7 + w11*w9) - St*(fA0*s + fA1*w8 + w12*w13 + w16) + w30*(-w18*w25*w28 + w25*w28*w29) - w30*(w18*w32*(St*w34*(sigma0High*w35 - w31) - w33) + w29*w32*(-St*w34*(sigma0Low*w35 - w31) + w33)) + (t + 1)*((1.0/4.0)*w5*cos(w0) + 1);
}

inline scalar u2(const point& pt, const scalar t, const Params& p)
{
    const scalar alpha = p.alpha;
    const scalar period = p.period;
    const scalar s = cos(alpha)*pt.x() + sin(alpha)*pt.y();
    const scalar tau = -sin(alpha)*pt.x() + cos(alpha)*pt.y();

    if (p.solution == 0)
    {
        return 1/((t + 1)*sqrt((1.0/4.0)*sin(2*constant::mathematical::pi*tau/period)*cos(constant::mathematical::pi*s) + 1));
    }

    return 1/((t + 1)*sqrt((1.0/4.0)*sin(2*constant::mathematical::pi*tau/period)*cos(constant::mathematical::pi*s) + 1));
}

inline symmTensor Gi(const point& pt, const Params& p)
{
    const scalar alpha = p.alpha;
    const scalar s1 = p.s1;
    const scalar s2 = p.s2;
    const scalar gil = p.gil;
    const scalar git = p.git;
    const scalar gin = p.gin;
    const scalar theta0 = p.theta0;
    const scalar theta1 = p.theta1;
    const scalar s = cos(alpha)*pt.x() + sin(alpha)*pt.y();
    const scalar thC = theta0 + theta1*min(max((s - s1)/(s2 - s1), scalar(0)), scalar(1));

    if (p.solution == 0)
    {
        const scalar w0 = cos(theta0);
        const scalar w1 = (w0*w0);
        const scalar w2 = sin(theta0);
        const scalar w3 = (w2*w2);
        return symmTensor(gil*w1 + git*w3, w0*w2*(gil - git), 0, gil*w3 + git*w1, 0, gin);
    }

    const scalar w0 = cos(thC);
    const scalar w1 = (w0*w0);
    const scalar w2 = sin(thC);
    const scalar w3 = (w2*w2);
    return symmTensor(gil*w1 + git*w3, w0*w2*(gil - git), 0, gil*w3 + git*w1, 0, gin);
}

inline symmTensor Ge(const point& pt, const Params& p)
{
    const scalar alpha = p.alpha;
    const scalar s1 = p.s1;
    const scalar s2 = p.s2;
    const scalar gel = p.gel;
    const scalar get = p.get;
    const scalar gen = p.gen;
    const scalar theta0 = p.theta0;
    const scalar theta1 = p.theta1;
    const scalar s = cos(alpha)*pt.x() + sin(alpha)*pt.y();
    const scalar thC = theta0 + theta1*min(max((s - s1)/(s2 - s1), scalar(0)), scalar(1));

    if (p.solution == 0)
    {
        const scalar w0 = cos(theta0);
        const scalar w1 = (w0*w0);
        const scalar w2 = sin(theta0);
        const scalar w3 = (w2*w2);
        return symmTensor(gel*w1 + get*w3, w0*w2*(gel - get), 0, gel*w3 + get*w1, 0, gen);
    }

    const scalar w0 = cos(thC);
    const scalar w1 = (w0*w0);
    const scalar w2 = sin(thC);
    const scalar w3 = (w2*w2);
    return symmTensor(gel*w1 + get*w3, w0*w2*(gel - get), 0, gel*w3 + get*w1, 0, gen);
}

inline scalar sigmaB(const point& pt, const label zone, const Params& p)
{
    const scalar alpha = p.alpha;
    const scalar period = p.period;
    const scalar sigma0Low = p.sigma0Low;
    const scalar sigma0High = p.sigma0High;
    const scalar sigmaEps = p.sigmaEps;
    const scalar s = cos(alpha)*pt.x() + sin(alpha)*pt.y();
    const scalar tau = -sin(alpha)*pt.x() + cos(alpha)*pt.y();
    checkZone(zone, false);

    if (p.solution == 0)
    {
        if (zone == bathLow)
        {
            return sigma0Low;
        }
        if (zone == bathHigh)
        {
            return sigma0High;
        }
        return 0;
    }

    if (zone == bathLow)
    {
        const scalar w0 = 2*constant::mathematical::pi;
        return sigma0Low*(sigmaEps*((1.0/2.0)*sin(s*w0) + 1)*sin(tau*w0/period) + 1);
    }
    if (zone == bathHigh)
    {
        const scalar w0 = 2*constant::mathematical::pi;
        return sigma0High*(sigmaEps*((1.0/2.0)*sin(s*w0) + 1)*sin(tau*w0/period) + 1);
    }
    return 0;
}

inline scalar divGiGradPhiI(const point& pt, const scalar t, const Params& p)
{
    const scalar alpha = p.alpha;
    const scalar s1 = p.s1;
    const scalar s2 = p.s2;
    const scalar period = p.period;
    const scalar gil = p.gil;
    const scalar git = p.git;
    const scalar theta0 = p.theta0;
    const scalar theta1 = p.theta1;
    const scalar hA0 = p.hA[0];
    const scalar hA1 = p.hA[1];
    const scalar hA2 = p.hA[2];
    const scalar cmAmplitude = p.cmAmplitude;
    const scalar cmOmega = p.cmOmega;
    const scalar s = cos(alpha)*pt.x() + sin(alpha)*pt.y();
    const scalar tau = -sin(alpha)*pt.x() + cos(alpha)*pt.y();
    const scalar St = sqrt(1 + t);
    const scalar Tt = cmAmplitude*sin(cmOmega*t);

    if (p.solution == 0)
    {
        return 0;
    }

    const scalar w0 = s - s1;
    const scalar w1 = -s2;
    const scalar w2 = -s1 - w1;
    const scalar w3 = (1.0/(w2));
    const scalar w4 = w0*w3;
    const scalar w5 = -alpha + theta0 + theta1*w4;
    const scalar w6 = cos(w5);
    const scalar w7 = (w6*w6);
    const scalar w8 = (1.0/(period));
    const scalar w9 = constant::mathematical::pi*w8;
    const scalar w10 = tau*w9;
    const scalar w11 = 4*w10;
    const scalar w12 = cos(w11);
    const scalar w13 = 1 - 3.0/10.0*(s*s);
    const scalar w14 = 2*w10;
    const scalar w15 = sin(w14);
    const scalar w16 = w15*w9;
    const scalar w17 = hA1*w16;
    const scalar w18 = St*(4*constant::mathematical::pi*hA2*s*w12*w8 - 2*w13*w17);
    const scalar w19 = (2.0/5.0)*w16;
    const scalar w20 = s + w1;
    const scalar w21 = gil*w7;
    const scalar w22 = sin(w5);
    const scalar w23 = (w22*w22);
    const scalar w24 = git*w23;
    const scalar w25 = w21 + w24;
    const scalar w26 = (1.0/(w25));
    const scalar w27 = w26*w4;
    const scalar w28 = cos(w14);
    const scalar w29 = ((constant::mathematical::pi)*(constant::mathematical::pi))/(period*period);
    const scalar w30 = w28*w29;
    const scalar w31 = (4.0/5.0)*w30;
    const scalar w32 = gil - git;
    const scalar w33 = w22*w6;
    const scalar w34 = w32*w33;
    const scalar w35 = w31*w34;
    const scalar w36 = 1 - w4;
    const scalar w37 = w26*w36;
    const scalar w38 = w27*w35 - w35*w37;
    const scalar w39 = w38*w4;
    const scalar w40 = hA2*w12;
    const scalar w41 = (6.0/5.0)*w17;
    const scalar w42 = s*w41 + 4*w40*w9;
    const scalar w43 = St*w42;
    const scalar w44 = hA2*sin(w11);
    const scalar w45 = 16*w29*w44;
    const scalar w46 = hA1*w28;
    const scalar w47 = w29*w46;
    const scalar w48 = St*(-s*w45 - 4*w13*w47);
    const scalar w49 = w25*w43 + w34*w48;
    const scalar w50 = -w49;
    const scalar w51 = w27*w50 + w37*w49;
    const scalar w52 = w4*w51;
    const scalar w53 = theta1*w3;
    const scalar w54 = w32*w53;
    const scalar w55 = w54*(Tt*(-w19 + w20*w39) + w18 + w20*w52);
    const scalar w56 = gil*w23;
    const scalar w57 = git*w7;
    const scalar w58 = ((constant::mathematical::pi)*(constant::mathematical::pi)*(constant::mathematical::pi));
    const scalar w59 = (1.0/(period*period*period));
    const scalar w60 = w58*w59;
    const scalar w61 = w15*w60;
    const scalar w62 = (8.0/5.0)*w34;
    const scalar w63 = w20*w4;
    const scalar w64 = (12.0/5.0)*s*w47 - w45;
    const scalar w65 = St*w25;
    const scalar w66 = St*w34;
    const scalar w67 = w64*w65 + w66*(8*hA1*w13*w61 - 64*s*w40*w60);
    const scalar w68 = w31*w54;
    const scalar w69 = w23*w68;
    const scalar w70 = w68*w7;
    const scalar w71 = w69 - w70;
    const scalar w72 = w26*w3;
    const scalar w73 = gil*w33;
    const scalar w74 = 2*w3;
    const scalar w75 = theta1*w74;
    const scalar w76 = git*w33;
    const scalar w77 = w73*w75 - w75*w76;
    const scalar w78 = (1.0/(w25*w25));
    const scalar w79 = w36*w78;
    const scalar w80 = w77*w79;
    const scalar w81 = w77*w78;
    const scalar w82 = w4*w81;
    const scalar w83 = w30*w62*w72 - w35*w80 + w35*w82;
    const scalar w84 = w20*w3;
    const scalar w85 = w38*w84 + w39;
    const scalar w86 = w48*w54;
    const scalar w87 = -w77;
    const scalar w88 = -w23*w86 + w41*w65 + w43*w87 + w64*w66 + w7*w86;
    const scalar w89 = w43 + w51*w84 + w52 + w63*(-w27*w88 + w37*w88 - w49*w72 + w49*w80 + w50*w72 + w50*w82);
    const scalar w90 = w23*w32;
    const scalar w91 = (1.0/(w2*w2));
    const scalar w92 = theta1*w0*w26*w31*w91;
    const scalar w93 = (3.0/5.0)*w46;
    const scalar w94 = hA0 - s*w93 + w44;
    const scalar w95 = St*w94;
    const scalar w96 = w18*w34;
    const scalar w97 = w25*w95 + w96;
    const scalar w98 = -w97;
    const scalar w99 = w27*w98 + w37*w97;
    const scalar w100 = w4*w98;
    const scalar w101 = w79*w97;
    const scalar w102 = w18*w54;
    const scalar w103 = -w102*w23 + w102*w7 + w34*w43 - w65*w93 + w87*w95;
    const scalar w104 = -w103;
    const scalar w105 = w100*w81 + w101*w77 + w103*w37 + w104*w27 - w72*w97 + w72*w98;
    const scalar w106 = -w19*w34 + (1.0/10.0)*w21 + (1.0/10.0)*w24;
    const scalar w107 = -w106;
    const scalar w108 = w106*w37 + w107*w27;
    const scalar w109 = -1.0/5.0*git*theta1*w22*w3*w6 - 2.0/5.0*constant::mathematical::pi*theta1*w15*w23*w3*w32*w8 + w19*w54*w7 + (1.0/5.0)*w53*w73;
    const scalar w110 = -w109;
    const scalar w111 = w107*w81;
    const scalar w112 = w106*w79;
    const scalar w113 = -w106*w72 + w107*w72 + w109*w27 + w110*w37 + w111*w4 + w112*w77;
    const scalar w114 = 2*w4;
    const scalar w115 = w20*w74;
    const scalar w116 = w74*w81;
    const scalar w117 = (theta1*theta1);
    const scalar w118 = w117*w91;
    const scalar w119 = 2*w118;
    const scalar w120 = -2*gil*w117*w7*w91 - 2*git*w117*w23*w91 + w119*w56 + w119*w57;
    const scalar w121 = -w120;
    const scalar w122 = w121*w78;
    const scalar w123 = 4*w53;
    const scalar w124 = w77*(w123*w73 - w123*w76)/(w25*w25*w25);
    const scalar w125 = w124*w36;
    const scalar w126 = w26*w74;
    const scalar w127 = w114*w81;
    const scalar w128 = 2*w80;
    const scalar w129 = -6.0/5.0*constant::mathematical::pi*St*hA1*w15*w22*w32*w6*w8 - 2*St*theta1*w3*w32*w42*w7 - St*w120*w94 + (6.0/5.0)*St*w46*w87 + 4*w118*w96 + w43*w75*w90;
    const scalar w130 = (1.0/5.0)*w118;
    const scalar w131 = w118*w16*w62 - w130*w21 - w130*w24 + w130*w56 + w130*w57;
    const scalar w132 = w107*w4;
    return -w23*w55 + w25*(-St*w93 + Tt*(w108*w74 + w113*w114 + w113*w115 + w63*(-w106*w116 + w106*w125 + w109*w126 + w109*w127 - w110*w126 + w110*w128 + w111*w74 + w112*w121 + w122*w132 + w124*w132 - w131*w27 + w131*w37)) + w105*w114 + w105*w115 + w63*(w100*w122 + w100*w124 + w101*w121 - w103*w126 + w103*w128 + w104*w126 + w104*w127 - w116*w97 + w116*w98 + w125*w97 + w129*w27 - w129*w37) + w74*w99) + w34*(Tt*(w63*(-w27*w71 + w37*w71 + w83) + w85) + w89) + w34*(Tt*(w63*(w32*w7*w92 + w37*w69 - w37*w70 + w83 - w90*w92) + w85) + w89) + w55*w7 + w87*(Tt*(w108*w4 + w108*w84 + w113*w63 + 1.0/10.0) + w105*w63 + w4*w99 + w84*w99 + w95) + (w56 + w57)*(Tt*(-w31 + w63*((8.0/5.0)*w15*w22*w26*w32*w36*w58*w59*w6 - w27*w61*w62)) + w48 + w63*(-w27*w67 + w37*w67));
}

inline scalar phiESource(const point& pt, const scalar t, const label zone, const Params& p)
{
    const scalar alpha = p.alpha;
    const scalar s1 = p.s1;
    const scalar s2 = p.s2;
    const scalar period = p.period;
    const scalar gil = p.gil;
    const scalar git = p.git;
    const scalar gel = p.gel;
    const scalar get = p.get;
    const scalar theta0 = p.theta0;
    const scalar theta1 = p.theta1;
    const scalar sigma0Low = p.sigma0Low;
    const scalar sigma0High = p.sigma0High;
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
    const scalar s = cos(alpha)*pt.x() + sin(alpha)*pt.y();
    const scalar tau = -sin(alpha)*pt.x() + cos(alpha)*pt.y();
    const scalar St = sqrt(1 + t);
    const scalar Tt = cmAmplitude*sin(cmOmega*t);
    checkZone(zone, true);

    if (p.solution == 0)
    {
        if (zone == bathLow)
        {
            return 0;
        }
        if (zone == bathHigh)
        {
            return 0;
        }
        return 0;
    }

    if (zone == bathLow)
    {
        const scalar w0 = (1.0/2.0)*s + 1;
        const scalar w1 = (1.0/(period));
        const scalar w2 = 2*constant::mathematical::pi;
        const scalar w3 = tau*w1;
        const scalar w4 = w2*w3;
        const scalar w5 = cos(w4);
        const scalar w6 = w1*w2*w5;
        const scalar w7 = 4*constant::mathematical::pi;
        const scalar w8 = w3*w7 + 3.0/10.0;
        const scalar w9 = fA3*exp((2.0/5.0)*s);
        const scalar w10 = w1*w7;
        const scalar w11 = sin(w4);
        const scalar w12 = (2.0/5.0)*constant::mathematical::pi*w1*w11;
        const scalar w13 = -s + s1;
        const scalar w14 = ((constant::mathematical::pi)*(constant::mathematical::pi));
        const scalar w15 = (1.0/(period*period));
        const scalar w16 = w14*w15;
        const scalar w17 = (4.0/5.0)*w16*w5;
        const scalar w18 = gel - get;
        const scalar w19 = alpha - theta0;
        const scalar w20 = cos(w19);
        const scalar w21 = sin(w19);
        const scalar w22 = w20*w21;
        const scalar w23 = w18*w22;
        const scalar w24 = gil - git;
        const scalar w25 = (w20*w20);
        const scalar w26 = gel*w25;
        const scalar w27 = (w21*w21);
        const scalar w28 = get*w27;
        const scalar w29 = w26 + w28;
        const scalar w30 = gil*w25;
        const scalar w31 = git*w27;
        const scalar w32 = (1.0/((w30 + w31)));
        const scalar w33 = (4.0/5.0)*w14*w15*w20*w21*w24*w29*w32*w5 - w17*w23;
        const scalar w34 = (1.0/(sigma0Low));
        const scalar w35 = s*w2;
        const scalar w36 = sin(w35);
        const scalar w37 = (1.0/2.0)*w36 + 1;
        const scalar w38 = sigmaEps*w11;
        const scalar w39 = w37*w38;
        const scalar w40 = w39 + 1;
        const scalar w41 = w34/w40;
        const scalar w42 = (1.0/(w40*w40));
        const scalar w43 = w12*w22;
        const scalar w44 = w29*w32;
        const scalar w45 = -w18*w43 - 1.0/10.0*w26 - 1.0/10.0*w28 - w44*(-w24*w43 - 1.0/10.0*w30 - 1.0/10.0*w31);
        const scalar w46 = w34*w45;
        const scalar w47 = w42*w46;
        const scalar w48 = sigmaEps*w37;
        const scalar w49 = w48*w6;
        const scalar w50 = 2*fA1;
        const scalar w51 = fA2*w11;
        const scalar w52 = w9*cos(w8);
        const scalar w53 = cos(w35);
        const scalar w54 = w38*w53;
        const scalar w55 = constant::mathematical::pi*w54;
        const scalar w56 = w13*w47;
        const scalar w57 = 2*w14;
        const scalar w58 = (sigmaEps*sigmaEps)*w46/(w40*w40*w40);
        const scalar w59 = sigma0Low*w40;
        const scalar w60 = 4*w16;
        const scalar w61 = (8.0/5.0)*((constant::mathematical::pi)*(constant::mathematical::pi)*(constant::mathematical::pi))*w11/(period*period*period);
        return sigma0Low*w49*(St*(fA2*w0*w6 - w10*w9*sin(w8)) + Tt*(-w12 + w13*(w33*w41 - w47*w49))) + sigma0Low*w55*(St*(fA0 + s*w50 + (1.0/2.0)*w51 + (2.0/5.0)*w52) + Tt*(-w41*w45 - w55*w56)) + w59*(St*(w50 + (4.0/25.0)*w52) + Tt*((w11*w11)*w13*(w53*w53)*w57*w58 + w2*w47*w54 + w36*w38*w56*w57)) + w59*(St*(-w0*w51*w60 - 16*w16*w52) + Tt*(w13*(-w10*w33*w34*w42*w48*w5 + 8*w16*(w37*w37)*(w5*w5)*w58 + w39*w47*w60 + w41*(-w22*w24*w44*w61 + w23*w61)) - w17));
    }
    if (zone == bathHigh)
    {
        const scalar w0 = (1.0/2.0)*s + 1;
        const scalar w1 = (1.0/(period));
        const scalar w2 = 2*constant::mathematical::pi;
        const scalar w3 = tau*w1;
        const scalar w4 = w2*w3;
        const scalar w5 = cos(w4);
        const scalar w6 = w1*w2*w5;
        const scalar w7 = 4*constant::mathematical::pi;
        const scalar w8 = w3*w7 + 3.0/10.0;
        const scalar w9 = fA3*exp((2.0/5.0)*s);
        const scalar w10 = w1*w7;
        const scalar w11 = sin(w4);
        const scalar w12 = (2.0/5.0)*constant::mathematical::pi*w1*w11;
        const scalar w13 = s - s2;
        const scalar w14 = ((constant::mathematical::pi)*(constant::mathematical::pi));
        const scalar w15 = (1.0/(period*period));
        const scalar w16 = w14*w15;
        const scalar w17 = (4.0/5.0)*w16*w5;
        const scalar w18 = gel - get;
        const scalar w19 = -alpha + theta0 + theta1;
        const scalar w20 = cos(w19);
        const scalar w21 = sin(w19);
        const scalar w22 = w20*w21;
        const scalar w23 = w18*w22;
        const scalar w24 = gil - git;
        const scalar w25 = (w20*w20);
        const scalar w26 = gel*w25;
        const scalar w27 = (w21*w21);
        const scalar w28 = get*w27;
        const scalar w29 = w26 + w28;
        const scalar w30 = gil*w25;
        const scalar w31 = git*w27;
        const scalar w32 = (1.0/((w30 + w31)));
        const scalar w33 = (4.0/5.0)*w14*w15*w20*w21*w24*w29*w32*w5 - w17*w23;
        const scalar w34 = (1.0/(sigma0High));
        const scalar w35 = s*w2;
        const scalar w36 = sin(w35);
        const scalar w37 = (1.0/2.0)*w36 + 1;
        const scalar w38 = sigmaEps*w11;
        const scalar w39 = w37*w38;
        const scalar w40 = w39 + 1;
        const scalar w41 = (1.0/(w40));
        const scalar w42 = w34*w41;
        const scalar w43 = (1.0/(w40*w40));
        const scalar w44 = w12*w22;
        const scalar w45 = w29*w32;
        const scalar w46 = -w18*w44 + (1.0/10.0)*w26 + (1.0/10.0)*w28 - w45*(-w24*w44 + (1.0/10.0)*w30 + (1.0/10.0)*w31);
        const scalar w47 = w34*w46;
        const scalar w48 = w43*w47;
        const scalar w49 = sigmaEps*w37;
        const scalar w50 = w49*w6;
        const scalar w51 = 2*fA1;
        const scalar w52 = fA2*w11;
        const scalar w53 = w9*cos(w8);
        const scalar w54 = cos(w35);
        const scalar w55 = w38*w54;
        const scalar w56 = constant::mathematical::pi*w55;
        const scalar w57 = w13*w48;
        const scalar w58 = 2*w14;
        const scalar w59 = (sigmaEps*sigmaEps)*w47/(w40*w40*w40);
        const scalar w60 = sigma0High*w40;
        const scalar w61 = 4*w16;
        const scalar w62 = (8.0/5.0)*((constant::mathematical::pi)*(constant::mathematical::pi)*(constant::mathematical::pi))*w11/(period*period*period);
        return sigma0High*w50*(St*(fA2*w0*w6 - w10*w9*sin(w8)) + Tt*(-w12 + w13*(w33*w42 - w48*w50))) + sigma0High*w56*(St*(fA0 + s*w51 + (1.0/2.0)*w52 + (2.0/5.0)*w53) + Tt*(w34*w41*w46 - w56*w57)) + w60*(St*(w51 + (4.0/25.0)*w53) + Tt*((w11*w11)*w13*(w54*w54)*w58*w59 - w2*w48*w55 + w36*w38*w57*w58)) + w60*(St*(-w0*w52*w61 - 16*w16*w53) + Tt*(w13*(-w10*w33*w34*w43*w49*w5 + 8*w16*(w37*w37)*(w5*w5)*w59 + w39*w48*w61 + w42*(-w22*w24*w45*w62 + w23*w62)) - w17));
    }
    const scalar w0 = s - s1;
    const scalar w1 = -s2;
    const scalar w2 = -s1 - w1;
    const scalar w3 = (1.0/(w2));
    const scalar w4 = w0*w3;
    const scalar w5 = -alpha + theta0 + theta1*w4;
    const scalar w6 = cos(w5);
    const scalar w7 = (w6*w6);
    const scalar w8 = theta1*w3;
    const scalar w9 = w7*w8;
    const scalar w10 = gil - git;
    const scalar w11 = (1.0/(period));
    const scalar w12 = tau*w11;
    const scalar w13 = 4*constant::mathematical::pi;
    const scalar w14 = w12*w13;
    const scalar w15 = cos(w14);
    const scalar w16 = 2*constant::mathematical::pi;
    const scalar w17 = w11*w16;
    const scalar w18 = w12*w16;
    const scalar w19 = sin(w18);
    const scalar w20 = hA1*(1 - 3.0/10.0*(s*s));
    const scalar w21 = w19*w20;
    const scalar w22 = St*(4*constant::mathematical::pi*hA2*s*w11*w15 - w17*w21);
    const scalar w23 = constant::mathematical::pi*w11;
    const scalar w24 = w19*w23;
    const scalar w25 = (2.0/5.0)*w24;
    const scalar w26 = s + w1;
    const scalar w27 = gil*w7;
    const scalar w28 = sin(w5);
    const scalar w29 = (w28*w28);
    const scalar w30 = git*w29;
    const scalar w31 = w27 + w30;
    const scalar w32 = (1.0/(w31));
    const scalar w33 = w32*w4;
    const scalar w34 = cos(w18);
    const scalar w35 = ((constant::mathematical::pi)*(constant::mathematical::pi));
    const scalar w36 = w35/(period*period);
    const scalar w37 = w34*w36;
    const scalar w38 = (4.0/5.0)*w37;
    const scalar w39 = w28*w6;
    const scalar w40 = w10*w39;
    const scalar w41 = w38*w40;
    const scalar w42 = 1 - w4;
    const scalar w43 = w32*w42;
    const scalar w44 = w33*w41 - w41*w43;
    const scalar w45 = w4*w44;
    const scalar w46 = Tt*(-w25 + w26*w45);
    const scalar w47 = w11*w13;
    const scalar w48 = hA2*w15;
    const scalar w49 = hA1*s;
    const scalar w50 = (6.0/5.0)*w24;
    const scalar w51 = w47*w48 + w49*w50;
    const scalar w52 = St*w51;
    const scalar w53 = hA2*sin(w14);
    const scalar w54 = 16*w36;
    const scalar w55 = w53*w54;
    const scalar w56 = St*(-s*w55 - 4*w20*w37);
    const scalar w57 = w31*w52 + w40*w56;
    const scalar w58 = -w57;
    const scalar w59 = w33*w58 + w43*w57;
    const scalar w60 = w4*w59;
    const scalar w61 = w10*(w22 + w26*w60 + w46);
    const scalar w62 = w29*w8;
    const scalar w63 = gil*w29;
    const scalar w64 = git*w7;
    const scalar w65 = ((constant::mathematical::pi)*(constant::mathematical::pi)*(constant::mathematical::pi));
    const scalar w66 = (1.0/(period*period*period));
    const scalar w67 = w65*w66;
    const scalar w68 = (8.0/5.0)*w40;
    const scalar w69 = w26*w4;
    const scalar w70 = Tt*(-w38 + w69*((8.0/5.0)*w10*w19*w28*w32*w42*w6*w65*w66 - w19*w33*w67*w68));
    const scalar w71 = (12.0/5.0)*w37*w49 - w55;
    const scalar w72 = St*w31;
    const scalar w73 = 8*w67;
    const scalar w74 = St*w40;
    const scalar w75 = w71*w72 + w74*(-64*s*w48*w67 + w21*w73);
    const scalar w76 = gel - get;
    const scalar w77 = (1.0/2.0)*s + 1;
    const scalar w78 = fA2*w34;
    const scalar w79 = w77*w78;
    const scalar w80 = w14 + 3.0/10.0;
    const scalar w81 = sin(w80);
    const scalar w82 = exp((2.0/5.0)*s);
    const scalar w83 = fA3*w82;
    const scalar w84 = w81*w83;
    const scalar w85 = St*(w17*w79 - w47*w84);
    const scalar w86 = s*w16;
    const scalar w87 = sin(w86);
    const scalar w88 = (1.0/2.0)*w87 + 1;
    const scalar w89 = sigmaEps*w34;
    const scalar w90 = w88*w89;
    const scalar w91 = 2*fA1;
    const scalar w92 = fA2*w19;
    const scalar w93 = w83*cos(w80);
    const scalar w94 = fA0 + s*w91 + (1.0/2.0)*w92 + (2.0/5.0)*w93;
    const scalar w95 = St*w94;
    const scalar w96 = sigma0High*w95;
    const scalar w97 = w36*w92;
    const scalar w98 = St*(-w54*w93 - 4*w77*w97);
    const scalar w99 = w39*w76;
    const scalar w100 = -w98*w99;
    const scalar w101 = sigmaEps*w19;
    const scalar w102 = w101*w88;
    const scalar w103 = w102 + 1;
    const scalar w104 = gel*w7 + get*w29;
    const scalar w105 = sigma0High*w103 - w104;
    const scalar w106 = w23*w84;
    const scalar w107 = -8.0/5.0*w106 + w23*w78;
    const scalar w108 = St*w107;
    const scalar w109 = w100 + w105*w108 + w17*w90*w96;
    const scalar w110 = (1.0/(w104));
    const scalar w111 = w110*w4;
    const scalar w112 = sigma0Low*w90;
    const scalar w113 = w112*w17;
    const scalar w114 = sigma0Low*w103 - w104;
    const scalar w115 = -w100 - w108*w114 - w113*w95;
    const scalar w116 = w110*w42;
    const scalar w117 = w109*w111 + w115*w116;
    const scalar w118 = w117*w4;
    const scalar w119 = w76*(w118*w26 + w46 + w85);
    const scalar w120 = gel*w29;
    const scalar w121 = get*w7;
    const scalar w122 = 4*w102*w36;
    const scalar w123 = St*w99;
    const scalar w124 = w123*(64*fA3*w65*w66*w81*w82 - w73*w79);
    const scalar w125 = -32.0/5.0*w36*w93 - 2*w97;
    const scalar w126 = sigma0Low*w95;
    const scalar w127 = St*w125;
    const scalar w128 = w10*w38;
    const scalar w129 = w128*w62;
    const scalar w130 = w128*w9;
    const scalar w131 = w129 - w130;
    const scalar w132 = w3*w32;
    const scalar w133 = gil*w39;
    const scalar w134 = 2*w3;
    const scalar w135 = theta1*w134;
    const scalar w136 = git*w39;
    const scalar w137 = w133*w135 - w135*w136;
    const scalar w138 = (1.0/(w31*w31));
    const scalar w139 = w138*w42;
    const scalar w140 = w137*w139;
    const scalar w141 = w137*w138;
    const scalar w142 = w141*w4;
    const scalar w143 = w132*w37*w68 - w140*w41 + w142*w41;
    const scalar w144 = w26*w3;
    const scalar w145 = w144*w44 + w45;
    const scalar w146 = Tt*(w145 + w69*(-w131*w33 + w131*w43 + w143));
    const scalar w147 = w10*w56;
    const scalar w148 = -w137;
    const scalar w149 = hA1*w50*w72 - w147*w62 + w147*w9 + w148*w52 + w71*w74;
    const scalar w150 = w144*w59 + w52 + w60 + w69*(-w132*w57 + w132*w58 + w140*w57 + w142*w58 - w149*w33 + w149*w43);
    const scalar w151 = (1.0/(w2*w2));
    const scalar w152 = theta1*w0*w128*w151*w32;
    const scalar w153 = Tt*(w145 + w69*(w129*w43 - w130*w43 + w143 - w152*w29 + w152*w7));
    const scalar w154 = hA1*w34;
    const scalar w155 = (3.0/5.0)*w154;
    const scalar w156 = hA0 - s*w155 + w53;
    const scalar w157 = St*w156;
    const scalar w158 = w22*w40;
    const scalar w159 = w157*w31 + w158;
    const scalar w160 = -w159;
    const scalar w161 = w159*w43 + w160*w33;
    const scalar w162 = w160*w4;
    const scalar w163 = w139*w159;
    const scalar w164 = w10*w22;
    const scalar w165 = w148*w157 - w155*w72 - w164*w62 + w164*w9 + w40*w52;
    const scalar w166 = -w165;
    const scalar w167 = -w132*w159 + w132*w160 + w137*w163 + w141*w162 + w165*w43 + w166*w33;
    const scalar w168 = -w25*w40 + (1.0/10.0)*w27 + (1.0/10.0)*w30;
    const scalar w169 = -w168;
    const scalar w170 = w168*w43 + w169*w33;
    const scalar w171 = -1.0/5.0*git*theta1*w28*w3*w6 - 2.0/5.0*constant::mathematical::pi*theta1*w10*w11*w19*w29*w3 + w10*w25*w9 + (1.0/5.0)*w133*w8;
    const scalar w172 = -w171;
    const scalar w173 = w141*w169;
    const scalar w174 = w139*w168;
    const scalar w175 = -w132*w168 + w132*w169 + w137*w174 + w171*w33 + w172*w43 + w173*w4;
    const scalar w176 = Tt*(w144*w170 + w170*w4 + w175*w69 + 1.0/10.0);
    const scalar w177 = w110*w3;
    const scalar w178 = (1.0/(w104*w104));
    const scalar w179 = w135*w39;
    const scalar w180 = gel*w179 - get*w179;
    const scalar w181 = w178*w180;
    const scalar w182 = w181*w42;
    const scalar w183 = cos(w86);
    const scalar w184 = constant::mathematical::pi*w183;
    const scalar w185 = sigma0High*w101;
    const scalar w186 = w180 + w184*w185;
    const scalar w187 = St*w105;
    const scalar w188 = (16.0/25.0)*w106;
    const scalar w189 = St*(w91 + (4.0/25.0)*w93);
    const scalar w190 = 2*w35;
    const scalar w191 = w11*w183*w190*w89;
    const scalar w192 = w76*w98;
    const scalar w193 = -w127*w99 + w192*w62 - w192*w9;
    const scalar w194 = sigma0Low*w101;
    const scalar w195 = w180 + w184*w194;
    const scalar w196 = w108 + w117*w144 + w118 + w69*(w109*w177 + w109*w181*w4 + w111*(sigma0High*w17*w189*w90 + w108*w186 - w187*w188 + w191*w96 + w193) - w115*w177 + w115*w182 + w116*((16.0/25.0)*constant::mathematical::pi*St*fA3*w11*w114*w81*w82 - w108*w195 - w113*w189 - w126*w191 - w193));
    const scalar w197 = w85*w99;
    const scalar w198 = St*w105*w94 - w197;
    const scalar w199 = -w114*w95 + w197;
    const scalar w200 = w111*w198 + w116*w199;
    const scalar w201 = w181*w198;
    const scalar w202 = w199*w42;
    const scalar w203 = w76*w85;
    const scalar w204 = -w108*w99 + w203*w62 - w203*w9;
    const scalar w205 = w105*w189 + w186*w95 + w204;
    const scalar w206 = -w114*w189 - w195*w95 - w204;
    const scalar w207 = w111*w205 + w116*w206 + w177*w198 - w177*w199 + w181*w202 + w201*w4;
    const scalar w208 = 2*w4;
    const scalar w209 = w134*w26;
    const scalar w210 = w134*w141;
    const scalar w211 = (theta1*theta1);
    const scalar w212 = w151*w211;
    const scalar w213 = 2*w212;
    const scalar w214 = -2*gil*w151*w211*w7 - 2*git*w151*w211*w29 + w213*w63 + w213*w64;
    const scalar w215 = -w214;
    const scalar w216 = w138*w215;
    const scalar w217 = 4*w8;
    const scalar w218 = w137*(w133*w217 - w136*w217)/(w31*w31*w31);
    const scalar w219 = w218*w42;
    const scalar w220 = w134*w32;
    const scalar w221 = w141*w208;
    const scalar w222 = 2*w140;
    const scalar w223 = w135*w29;
    const scalar w224 = 4*w212;
    const scalar w225 = -6.0/5.0*constant::mathematical::pi*St*hA1*w10*w11*w19*w28*w6 - 2*St*theta1*w10*w3*w51*w7 + (6.0/5.0)*St*w148*w154 - St*w156*w214 + w10*w223*w52 + w158*w224;
    const scalar w226 = (1.0/5.0)*w212;
    const scalar w227 = w212*w24*w68 - w226*w27 - w226*w30 + w226*w63 + w226*w64;
    const scalar w228 = w169*w4;
    const scalar w229 = Tt*(w134*w170 + w175*w208 + w175*w209 + w69*(w134*w173 - w168*w210 + w168*w219 + w171*w220 + w171*w221 - w172*w220 + w172*w222 + w174*w215 + w216*w228 + w218*w228 - w227*w33 + w227*w43));
    const scalar w230 = -2*gel*w151*w211*w7 - 2*get*w151*w211*w29 + w120*w213 + w121*w213;
    const scalar w231 = -w178*w230;
    const scalar w232 = w198*w4;
    const scalar w233 = w217*w39;
    const scalar w234 = w180*(gel*w233 - get*w233)/(w104*w104*w104);
    const scalar w235 = w110*w134;
    const scalar w236 = w190*w87;
    const scalar w237 = 2*w189;
    const scalar w238 = (8.0/125.0)*w93;
    const scalar w239 = w108*w76;
    const scalar w240 = w123*w188 - w135*w239*w7 + w197*w224 + w223*w239;
    return w104*(w134*w200 + w189 + w207*w208 + w207*w209 + w229 + w69*(w111*(w186*w237 + w187*w238 + w240 + w95*(-w185*w236 - w230)) + w116*(-St*w114*w238 - w195*w237 - w240 - w95*(-w194*w236 - w230)) - w134*w181*w199 + w134*w201 + w181*w205*w208 + 2*w182*w206 + w202*w231 + w202*w234 + w205*w235 - w206*w235 + w231*w232 + w232*w234)) - w119*w62 + w119*w9 + w148*(w144*w161 + w157 + w161*w4 + w167*w69 + w176) - w180*(w144*w200 + w176 + w200*w4 + w207*w69 + w95) + w31*(-St*w155 + w134*w161 + w167*w208 + w167*w209 + w229 + w69*(-w159*w210 + w159*w219 + w160*w210 + w162*w216 + w162*w218 + w163*w215 - w165*w220 + w165*w222 + w166*w220 + w166*w221 + w225*w33 - w225*w43)) + w40*(w146 + w150) + w40*(w150 + w153) - w61*w62 + w61*w9 + w99*(w146 + w196) + w99*(w153 + w196) + (w120 + w121)*(w69*(w111*(4*constant::mathematical::pi*St*sigma0High*sigmaEps*w107*w11*w34*w88 + St*w105*w125 - w122*w96 - w124) + w116*(-w108*w112*w47 - w114*w127 + w122*w126 + w124)) + w70 + w98) + (w63 + w64)*(w56 + w69*(-w33*w75 + w43*w75) + w70);
}

inline scalar interfaceFlux(const point& pt, const vector& n, const scalar t, const Params& p)
{
    const scalar alpha = p.alpha;
    const scalar s1 = p.s1;
    const scalar s2 = p.s2;
    const scalar period = p.period;
    const scalar gil = p.gil;
    const scalar git = p.git;
    const scalar gel = p.gel;
    const scalar get = p.get;
    const scalar theta0 = p.theta0;
    const scalar theta1 = p.theta1;
    const scalar sigma0Low = p.sigma0Low;
    const scalar sigma0High = p.sigma0High;
    const scalar sigmaEps = p.sigmaEps;
    const scalar fA0 = p.fA[0];
    const scalar fA1 = p.fA[1];
    const scalar fA2 = p.fA[2];
    const scalar fA3 = p.fA[3];
    const scalar cmAmplitude = p.cmAmplitude;
    const scalar cmOmega = p.cmOmega;
    const scalar P1 = p.P1;
    const scalar s = cos(alpha)*pt.x() + sin(alpha)*pt.y();
    const scalar tau = -sin(alpha)*pt.x() + cos(alpha)*pt.y();
    const scalar St = sqrt(1 + t);
    const scalar Tt = cmAmplitude*sin(cmOmega*t);

    if (p.solution == 0)
    {
        const scalar w0 = sin(alpha);
        const scalar w1 = alpha - theta0;
        const scalar w2 = cos(w1);
        const scalar w3 = sin(w1);
        const scalar w4 = P1*St;
        const scalar w5 = w2*w3*w4*(gel - get);
        const scalar w6 = cos(alpha);
        const scalar w7 = gel*(w2*w2) + get*(w3*w3);
        return n.x()*(w0*w5 + w4*w6*w7) + n.y()*(P1*St*w0*w7 - w5*w6) + n.z()*(0);
    }

    const scalar w0 = cos(alpha);
    const scalar w1 = (1.0/2.0)*s + 1;
    const scalar w2 = 2*constant::mathematical::pi;
    const scalar w3 = (1.0/(period));
    const scalar w4 = tau*w3;
    const scalar w5 = w2*w4;
    const scalar w6 = cos(w5);
    const scalar w7 = fA2*w6;
    const scalar w8 = w2*w3;
    const scalar w9 = 4*constant::mathematical::pi;
    const scalar w10 = w4*w9 + 3.0/10.0;
    const scalar w11 = fA3*exp((2.0/5.0)*s);
    const scalar w12 = w11*sin(w10);
    const scalar w13 = St*(w1*w7*w8 - w12*w3*w9);
    const scalar w14 = sin(w5);
    const scalar w15 = constant::mathematical::pi*w3;
    const scalar w16 = (2.0/5.0)*w14*w15;
    const scalar w17 = -s2;
    const scalar w18 = (1.0/((-s1 - w17)));
    const scalar w19 = w18*(s - s1);
    const scalar w20 = -alpha + theta0 + theta1*w19;
    const scalar w21 = cos(w20);
    const scalar w22 = (w21*w21);
    const scalar w23 = gil*w22;
    const scalar w24 = sin(w20);
    const scalar w25 = (w24*w24);
    const scalar w26 = git*w25;
    const scalar w27 = w23 + w26;
    const scalar w28 = (1.0/(w27));
    const scalar w29 = w19*w28;
    const scalar w30 = gil - git;
    const scalar w31 = w21*w24;
    const scalar w32 = ((constant::mathematical::pi)*(constant::mathematical::pi))/(period*period);
    const scalar w33 = (4.0/5.0)*w30*w31*w32*w6;
    const scalar w34 = 1 - w19;
    const scalar w35 = w28*w34;
    const scalar w36 = s + w17;
    const scalar w37 = w19*w36;
    const scalar w38 = s*w2;
    const scalar w39 = (1.0/2.0)*sin(w38) + 1;
    const scalar w40 = 2*fA1;
    const scalar w41 = fA2*w14;
    const scalar w42 = w11*cos(w10);
    const scalar w43 = fA0 + s*w40 + (1.0/2.0)*w41 + (2.0/5.0)*w42;
    const scalar w44 = St*w43;
    const scalar w45 = sigmaEps*w39*w44*w6*w8;
    const scalar w46 = gel - get;
    const scalar w47 = w31*w46;
    const scalar w48 = -St*w47*(-4*w1*w32*w41 - 16*w32*w42);
    const scalar w49 = sigmaEps*w14;
    const scalar w50 = w39*w49 + 1;
    const scalar w51 = gel*w22 + get*w25;
    const scalar w52 = sigma0High*w50 - w51;
    const scalar w53 = St*(-8.0/5.0*w12*w15 + w15*w7);
    const scalar w54 = (1.0/(w51));
    const scalar w55 = w19*w54;
    const scalar w56 = sigma0Low*w50 - w51;
    const scalar w57 = w34*w54;
    const scalar w58 = Tt*(-w16 + w37*(w29*w33 - w33*w35)) + w13 + w37*(w55*(sigma0High*w45 + w48 + w52*w53) + w57*(-sigma0Low*w45 - w48 - w53*w56));
    const scalar w59 = w13*w47;
    const scalar w60 = St*w43*w52 - w59;
    const scalar w61 = -w44*w56 + w59;
    const scalar w62 = w55*w60 + w57*w61;
    const scalar w63 = w18*w36;
    const scalar w64 = w16*w30;
    const scalar w65 = (1.0/10.0)*w23 + (1.0/10.0)*w26 - w31*w64;
    const scalar w66 = -w65;
    const scalar w67 = w29*w66 + w35*w65;
    const scalar w68 = w18*w28;
    const scalar w69 = theta1*w18;
    const scalar w70 = w22*w69;
    const scalar w71 = w31*w69;
    const scalar w72 = (1.0/5.0)*gil*w71 - 1.0/5.0*git*theta1*w18*w21*w24 - 2.0/5.0*constant::mathematical::pi*theta1*w14*w18*w25*w3*w30 + w64*w70;
    const scalar w73 = 2*w71;
    const scalar w74 = (gil*w73 - git*w73)/(w27*w27);
    const scalar w75 = w18*w54;
    const scalar w76 = gel*w73 - get*w73;
    const scalar w77 = w76/(w51*w51);
    const scalar w78 = St*(w40 + (4.0/25.0)*w42);
    const scalar w79 = constant::mathematical::pi*w49*cos(w38);
    const scalar w80 = w13*w46;
    const scalar w81 = w25*w69*w80 - w47*w53 - w70*w80;
    const scalar w82 = Tt*(w19*w67 + w37*(w19*w66*w74 + w29*w72 + w34*w65*w74 - w35*w72 - w65*w68 + w66*w68) + w63*w67 + 1.0/10.0) + w19*w62 + w37*(w19*w60*w77 + w34*w61*w77 + w55*(w44*(sigma0High*w79 + w76) + w52*w78 + w81) + w57*(-w44*(sigma0Low*w79 + w76) - w56*w78 - w81) + w60*w75 - w61*w75) + w44 + w62*w63;
    const scalar w83 = w47*w58 + w51*w82;
    const scalar w84 = sin(alpha);
    const scalar w85 = w47*w82 + w58*(gel*w25 + get*w22);
    return n.x()*(w0*w83 - w84*w85) + n.y()*(w0*w85 + w83*w84) + n.z()*(0);
}

inline scalar G(const point& pt, const Params& p)
{
    const scalar alpha = p.alpha;
    const scalar period = p.period;
    const scalar s = cos(alpha)*pt.x() + sin(alpha)*pt.y();
    const scalar tau = -sin(alpha)*pt.x() + cos(alpha)*pt.y();

    if (p.solution == 0)
    {
        return (1.0/4.0)*sin(2*constant::mathematical::pi*tau/period)*cos(constant::mathematical::pi*s) + 1;
    }

    return (1.0/4.0)*sin(2*constant::mathematical::pi*tau/period)*cos(constant::mathematical::pi*s) + 1;
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
    {0, 0, 1, 1, 0, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0, -0.69648824658912223},
    {0, 0, 1, 1, 0, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0.13700000000000001, -0.74266682494442293},
    {0, 0, 1, 1, 0, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0, 0.044923311931916378},
    {0, 0, 1, 1, 0, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0.13700000000000001, 0.047901818303254205},
    {0, 0, 1, 1, 0, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0, -0.54851681181692935},
    {0, 0, 1, 1, 0, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0.13700000000000001, -0.58488458499577889},
    {0, 1, 1, 1, 0, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0, -0.55921477879752535},
    {0, 1, 1, 1, 0, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0.13700000000000001, -0.59629184880783648},
    {0, 1, 1, 1, 0, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0, 0.043739180689507516},
    {0, 1, 1, 1, 0, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0.13700000000000001, 0.046639176766338115},
    {0, 1, 1, 1, 0, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0, -0.60949293945769867},
    {0, 1, 1, 1, 0, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0.13700000000000001, -0.64990355313220871},
    {0, 0, 0, 1, 0, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0, 0.19875999999999999},
    {0, 0, 0, 1, 0, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0.13700000000000001, 0.2119381896950146},
    {0, 0, 0, 1, 0, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0, 0.023999999999999976},
    {0, 0, 0, 1, 0, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0.13700000000000001, 0.025591248504127317},
    {0, 0, 0, 1, 0, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0, -0.19848000000000002},
    {0, 0, 0, 1, 0, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0.13700000000000001, -0.21163962512913315},
    {1, 0, 1, 0, 0, 0.0085011710662745762, 0.14740329063661345, 0.050000000000000003, 0, 0.26902279800550444},
    {1, 0, 1, 0, 0, 0.0085011710662745762, 0.14740329063661345, 0.050000000000000003, 0.13700000000000001, 0.28685955320977158},
    {1, 0, 1, 0, 0, 0.018932091702127847, 0.69292248910233978, 0.050000000000000003, 0, -0.24393600047238506},
    {1, 0, 1, 0, 0, 0.018932091702127847, 0.69292248910233978, 0.050000000000000003, 0.13700000000000001, -0.26010945029965549},
    {1, 0, 1, 0, 0, 0.62398034980701278, -0.16141813112137662, 0.050000000000000003, 0, -0.043466519229646457},
    {1, 0, 1, 0, 0, 0.62398034980701278, -0.16141813112137662, 0.050000000000000003, 0.13700000000000001, -0.046348437300638007},
    {1, 0, 1, 1, 0, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0, 0.54356931418490084},
    {1, 0, 1, 1, 0, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0.13700000000000001, 0.57960905827182785},
    {1, 0, 1, 1, 0, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0, -0.20526991488608773},
    {1, 0, 1, 1, 0, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0.13700000000000001, -0.2188797250946225},
    {1, 0, 1, 1, 0, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0, 0.17257307297304736},
    {1, 0, 1, 1, 0, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0.13700000000000001, 0.18401501648225654},
    {1, 0, 1, 2, 0, 1.1759566764140661, 1.033843888213875, 0.050000000000000003, 0, 1.4062157242894286},
    {1, 0, 1, 2, 0, 1.1759566764140661, 1.033843888213875, 0.050000000000000003, 0.13700000000000001, 1.4994506687792579},
    {1, 0, 1, 2, 0, 1.5793695941545203, 0.65923568248403053, 0.050000000000000003, 0, 1.0306033082696604},
    {1, 0, 1, 2, 0, 1.5793695941545203, 0.65923568248403053, 0.050000000000000003, 0.13700000000000001, 1.0989343904626931},
    {1, 0, 1, 2, 0, 1.4491150631448377, 1.4591317739555714, 0.050000000000000003, 0, 0.24470310988252914},
    {1, 0, 1, 2, 0, 1.4491150631448377, 1.4591317739555714, 0.050000000000000003, 0.13700000000000001, 0.2609274206140243},
    {1, 1, 1, 0, 0, 0.0085011710662745762, 0.14740329063661345, 0.050000000000000003, 0, 0.20979360028745969},
    {1, 1, 1, 0, 0, 0.0085011710662745762, 0.14740329063661345, 0.050000000000000003, 0.13700000000000001, 0.25158850520669418},
    {1, 1, 1, 0, 0, 0.018932091702127847, 0.69292248910233978, 0.050000000000000003, 0, -0.22839681761611358},
    {1, 1, 1, 0, 0, 0.018932091702127847, 0.69292248910233978, 0.050000000000000003, 0.13700000000000001, -0.25841807945686363},
    {1, 1, 1, 0, 0, 0.62398034980701278, -0.16141813112137662, 0.050000000000000003, 0, -0.014011123755265758},
    {1, 1, 1, 0, 0, 0.62398034980701278, -0.16141813112137662, 0.050000000000000003, 0.13700000000000001, -0.028915252178676677},
    {1, 1, 1, 1, 0, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0, 0.41783639283734264},
    {1, 1, 1, 1, 0, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0.13700000000000001, 0.4481080032827196},
    {1, 1, 1, 1, 0, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0, -0.17642720843364326},
    {1, 1, 1, 1, 0, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0.13700000000000001, -0.16847500273375579},
    {1, 1, 1, 1, 0, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0, 0.25617611640852739},
    {1, 1, 1, 1, 0, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0.13700000000000001, 0.30113264889778885},
    {1, 1, 1, 2, 0, 1.1759566764140661, 1.033843888213875, 0.050000000000000003, 0, 1.4937164394061146},
    {1, 1, 1, 2, 0, 1.1759566764140661, 1.033843888213875, 0.050000000000000003, 0.13700000000000001, 1.5885503793864972},
    {1, 1, 1, 2, 0, 1.5793695941545203, 0.65923568248403053, 0.050000000000000003, 0, 1.2137031850993021},
    {1, 1, 1, 2, 0, 1.5793695941545203, 0.65923568248403053, 0.050000000000000003, 0.13700000000000001, 1.3315882059655848},
    {1, 1, 1, 2, 0, 1.4491150631448377, 1.4591317739555714, 0.050000000000000003, 0, 0.32362123030662826},
    {1, 1, 1, 2, 0, 1.4491150631448377, 1.4591317739555714, 0.050000000000000003, 0.13700000000000001, 0.34949574309606146},
    {1, 0, 0, 0, 0, 0.0085011710662745762, 0.14740329063661345, 0.050000000000000003, 0, -0.11213333333333335},
    {1, 0, 0, 0, 0, 0.0085011710662745762, 0.14740329063661345, 0.050000000000000003, 0.13700000000000001, -0.11956799995539499},
    {1, 0, 0, 0, 0, 0.018932091702127847, 0.69292248910233978, 0.050000000000000003, 0, 0.0062666666666666565},
    {1, 0, 0, 0, 0, 0.018932091702127847, 0.69292248910233978, 0.050000000000000003, 0.13700000000000001, 0.00668215933163324},
    {1, 0, 0, 0, 0, 0.62398034980701278, -0.16141813112137662, 0.050000000000000003, 0, 0.098668000000000006},
    {1, 0, 0, 0, 0, 0.62398034980701278, -0.16141813112137662, 0.050000000000000003, 0.13700000000000001, 0.10520988780855153},
    {1, 0, 0, 1, 0, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0, 0.10124},
    {1, 0, 0, 1, 0, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0.13700000000000001, 0.10795241660657717},
    {1, 0, 0, 1, 0, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0, 0.27600000000000002},
    {1, 0, 0, 1, 0, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0.13700000000000001, 0.29429935779746447},
    {1, 0, 0, 1, 0, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0, 0.49848000000000003},
    {1, 0, 0, 1, 0, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0.13700000000000001, 0.53153023143072498},
    {1, 0, 0, 2, 0, 1.1759566764140661, 1.033843888213875, 0.050000000000000003, 0, 0.500888},
    {1, 0, 0, 2, 0, 1.1759566764140661, 1.033843888213875, 0.050000000000000003, 0.13700000000000001, 0.53409788669730573},
    {1, 0, 0, 2, 0, 1.5793695941545203, 0.65923568248403053, 0.050000000000000003, 0, 0.56906666666666661},
    {1, 0, 0, 2, 0, 1.5793695941545203, 0.65923568248403053, 0.050000000000000003, 0.13700000000000001, 0.6067969367534195},
    {1, 0, 0, 2, 0, 1.4491150631448377, 1.4591317739555714, 0.050000000000000003, 0, 0.64142222222222223},
    {1, 0, 0, 2, 0, 1.4491150631448377, 1.4591317739555714, 0.050000000000000003, 0.13700000000000001, 0.68394981187327009},
    {2, 0, 1, 1, 0, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0, -0.15291893240422141},
    {2, 0, 1, 1, 0, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0.13700000000000001, -0.16305776667259508},
    {2, 0, 1, 1, 0, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0, -0.16034660295417136},
    {2, 0, 1, 1, 0, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0.13700000000000001, -0.17097790679136829},
    {2, 0, 1, 1, 0, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0, -0.37594373884388205},
    {2, 0, 1, 1, 0, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0.13700000000000001, -0.40086956851352235},
    {2, 1, 1, 1, 0, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0, -0.14137838596018273},
    {2, 1, 1, 1, 0, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0.13700000000000001, -0.14818384552511685},
    {2, 1, 1, 1, 0, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0, -0.13268802774413577},
    {2, 1, 1, 1, 0, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0.13700000000000001, -0.12183582596741768},
    {2, 1, 1, 1, 0, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0, -0.35331682304917134},
    {2, 1, 1, 1, 0, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0.13700000000000001, -0.34877090423441992},
    {2, 0, 0, 1, 0, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0, 0.29999999999999999},
    {2, 0, 0, 1, 0, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0.13700000000000001, 0.31989060630159177},
    {2, 0, 0, 1, 0, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0, 0.29999999999999999},
    {2, 0, 0, 1, 0, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0.13700000000000001, 0.31989060630159177},
    {2, 0, 0, 1, 0, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0, 0.29999999999999999},
    {2, 0, 0, 1, 0, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0.13700000000000001, 0.31989060630159177},
    {3, 0, 1, 0, 0, 0.0085011710662745762, 0.14740329063661345, 0.050000000000000003, 0, 0.42291197422763804},
    {3, 0, 1, 0, 1, 0.0085011710662745762, 0.14740329063661345, 0.050000000000000003, 0, 0.1623337908527015},
    {3, 0, 1, 0, 2, 0.0085011710662745762, 0.14740329063661345, 0.050000000000000003, 0, 0},
    {3, 0, 1, 0, 0, 0.0085011710662745762, 0.14740329063661345, 0.050000000000000003, 0.13700000000000001, 0.45095189282627429},
    {3, 0, 1, 0, 1, 0.0085011710662745762, 0.14740329063661345, 0.050000000000000003, 0.13700000000000001, 0.17309684926368826},
    {3, 0, 1, 0, 2, 0.0085011710662745762, 0.14740329063661345, 0.050000000000000003, 0.13700000000000001, 0},
    {3, 0, 1, 0, 0, 0.018932091702127847, 0.69292248910233978, 0.050000000000000003, 0, 2.0165599397598282},
    {3, 0, 1, 0, 1, 0.018932091702127847, 0.69292248910233978, 0.050000000000000003, 0, -3.6829726340234608},
    {3, 0, 1, 0, 2, 0.018932091702127847, 0.69292248910233978, 0.050000000000000003, 0, 0},
    {3, 0, 1, 0, 0, 0.018932091702127847, 0.69292248910233978, 0.050000000000000003, 0.13700000000000001, 2.1502619392442428},
    {3, 0, 1, 0, 1, 0.018932091702127847, 0.69292248910233978, 0.050000000000000003, 0.13700000000000001, -3.9271611629664513},
    {3, 0, 1, 0, 2, 0.018932091702127847, 0.69292248910233978, 0.050000000000000003, 0.13700000000000001, 0},
    {3, 0, 1, 0, 0, 0.62398034980701278, -0.16141813112137662, 0.050000000000000003, 0, 2.4465496673353853},
    {3, 0, 1, 0, 1, 0.62398034980701278, -0.16141813112137662, 0.050000000000000003, 0, -4.2955312363463198},
    {3, 0, 1, 0, 2, 0.62398034980701278, -0.16141813112137662, 0.050000000000000003, 0, 0},
    {3, 0, 1, 0, 0, 0.62398034980701278, -0.16141813112137662, 0.050000000000000003, 0.13700000000000001, 2.6087608547695802},
    {3, 0, 1, 0, 1, 0.62398034980701278, -0.16141813112137662, 0.050000000000000003, 0.13700000000000001, -4.5803336386075015},
    {3, 0, 1, 0, 2, 0.62398034980701278, -0.16141813112137662, 0.050000000000000003, 0.13700000000000001, 0},
    {3, 0, 1, 1, 0, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0, 0.68561833680408646},
    {3, 0, 1, 1, 1, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0, 1.2282384544343237},
    {3, 0, 1, 1, 2, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0, 0},
    {3, 0, 1, 1, 0, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0.13700000000000001, 0.73107621817249402},
    {3, 0, 1, 1, 1, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0.13700000000000001, 1.3096731462397526},
    {3, 0, 1, 1, 2, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0.13700000000000001, 0},
    {3, 0, 1, 1, 0, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0, -1.3605071253159293},
    {3, 0, 1, 1, 1, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0, 3.5312603105828151},
    {3, 0, 1, 1, 2, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0, 0},
    {3, 0, 1, 1, 0, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0.13700000000000001, -1.4507114973164945},
    {3, 0, 1, 1, 1, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0.13700000000000001, 3.7653900058702803},
    {3, 0, 1, 1, 2, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0.13700000000000001, 0},
    {3, 0, 1, 1, 0, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0, -4.538465107529488},
    {3, 0, 1, 1, 1, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0, 3.5591328722911721},
    {3, 0, 1, 1, 2, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0, 0},
    {3, 0, 1, 1, 0, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0.13700000000000001, -4.8393745164207562},
    {3, 0, 1, 1, 1, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0.13700000000000001, 3.7951105747504963},
    {3, 0, 1, 1, 2, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0.13700000000000001, 0},
    {3, 0, 1, 2, 0, 1.1759566764140661, 1.033843888213875, 0.050000000000000003, 0, 2.0648490301699352},
    {3, 0, 1, 2, 1, 1.1759566764140661, 1.033843888213875, 0.050000000000000003, 0, -1.8889967613143217},
    {3, 0, 1, 2, 2, 1.1759566764140661, 1.033843888213875, 0.050000000000000003, 0, 0},
    {3, 0, 1, 2, 0, 1.1759566764140661, 1.033843888213875, 0.050000000000000003, 0.13700000000000001, 2.2017526939410481},
    {3, 0, 1, 2, 1, 1.1759566764140661, 1.033843888213875, 0.050000000000000003, 0.13700000000000001, -2.0142410642619391},
    {3, 0, 1, 2, 2, 1.1759566764140661, 1.033843888213875, 0.050000000000000003, 0.13700000000000001, 0},
    {3, 0, 1, 2, 0, 1.5793695941545203, 0.65923568248403053, 0.050000000000000003, 0, -1.7426984102242788},
    {3, 0, 1, 2, 1, 1.5793695941545203, 0.65923568248403053, 0.050000000000000003, 0, 6.0994758905037925},
    {3, 0, 1, 2, 2, 1.5793695941545203, 0.65923568248403053, 0.050000000000000003, 0, 0},
    {3, 0, 1, 2, 0, 1.5793695941545203, 0.65923568248403053, 0.050000000000000003, 0.13700000000000001, -1.8582428368248824},
    {3, 0, 1, 2, 1, 1.5793695941545203, 0.65923568248403053, 0.050000000000000003, 0.13700000000000001, 6.5038834691173317},
    {3, 0, 1, 2, 2, 1.5793695941545203, 0.65923568248403053, 0.050000000000000003, 0.13700000000000001, 0},
    {3, 0, 1, 2, 0, 1.4491150631448377, 1.4591317739555714, 0.050000000000000003, 0, 1.5467532962627251},
    {3, 0, 1, 2, 1, 1.4491150631448377, 1.4591317739555714, 0.050000000000000003, 0, -1.5397805971779162},
    {3, 0, 1, 2, 2, 1.4491150631448377, 1.4591317739555714, 0.050000000000000003, 0, 0},
    {3, 0, 1, 2, 0, 1.4491150631448377, 1.4591317739555714, 0.050000000000000003, 0.13700000000000001, 1.6493061658015624},
    {3, 0, 1, 2, 1, 1.4491150631448377, 1.4591317739555714, 0.050000000000000003, 0.13700000000000001, -1.6418711626755689},
    {3, 0, 1, 2, 2, 1.4491150631448377, 1.4591317739555714, 0.050000000000000003, 0.13700000000000001, 0},
    {3, 1, 1, 0, 0, 0.0085011710662745762, 0.14740329063661345, 0.050000000000000003, 0, 0.71017626320462468},
    {3, 1, 1, 0, 1, 0.0085011710662745762, 0.14740329063661345, 0.050000000000000003, 0, -0.61584863994888694},
    {3, 1, 1, 0, 2, 0.0085011710662745762, 0.14740329063661345, 0.050000000000000003, 0, 0},
    {3, 1, 1, 0, 0, 0.0085011710662745762, 0.14740329063661345, 0.050000000000000003, 0.13700000000000001, 0.77553525241044385},
    {3, 1, 1, 0, 1, 0.0085011710662745762, 0.14740329063661345, 0.050000000000000003, 0.13700000000000001, -0.73367526061145172},
    {3, 1, 1, 0, 2, 0.0085011710662745762, 0.14740329063661345, 0.050000000000000003, 0.13700000000000001, 0},
    {3, 1, 1, 0, 0, 0.018932091702127847, 0.69292248910233978, 0.050000000000000003, 0, 2.1904227369936722},
    {3, 1, 1, 0, 1, 0.018932091702127847, 0.69292248910233978, 0.050000000000000003, 0, -4.0728727147405959},
    {3, 1, 1, 0, 2, 0.018932091702127847, 0.69292248910233978, 0.050000000000000003, 0, 0},
    {3, 1, 1, 0, 0, 0.018932091702127847, 0.69292248910233978, 0.050000000000000003, 0.13700000000000001, 2.3239526693859802},
    {3, 1, 1, 0, 1, 0.018932091702127847, 0.69292248910233978, 0.050000000000000003, 0.13700000000000001, -4.2707233919734371},
    {3, 1, 1, 0, 2, 0.018932091702127847, 0.69292248910233978, 0.050000000000000003, 0.13700000000000001, 0},
    {3, 1, 1, 0, 0, 0.62398034980701278, -0.16141813112137662, 0.050000000000000003, 0, 2.6423709137167837},
    {3, 1, 1, 0, 1, 0.62398034980701278, -0.16141813112137662, 0.050000000000000003, 0, -4.6763662719326851},
    {3, 1, 1, 0, 2, 0.62398034980701278, -0.16141813112137662, 0.050000000000000003, 0, 0},
    {3, 1, 1, 0, 0, 0.62398034980701278, -0.16141813112137662, 0.050000000000000003, 0.13700000000000001, 2.7965762584233116},
    {3, 1, 1, 0, 1, 0.62398034980701278, -0.16141813112137662, 0.050000000000000003, 0.13700000000000001, -4.9095291409269546},
    {3, 1, 1, 0, 2, 0.62398034980701278, -0.16141813112137662, 0.050000000000000003, 0.13700000000000001, 0},
    {3, 1, 1, 1, 0, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0, 0.71446892878371004},
    {3, 1, 1, 1, 1, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0, 1.7644469608308277},
    {3, 1, 1, 1, 2, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0, 0},
    {3, 1, 1, 1, 0, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0.13700000000000001, 0.72142635762857343},
    {3, 1, 1, 1, 1, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0.13700000000000001, 1.7015882973974656},
    {3, 1, 1, 1, 2, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0.13700000000000001, 0},
    {3, 1, 1, 1, 0, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0, -1.609268616653271},
    {3, 1, 1, 1, 1, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0, 4.1012471321833672},
    {3, 1, 1, 1, 2, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0, 0},
    {3, 1, 1, 1, 0, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0.13700000000000001, -1.7679365936544926},
    {3, 1, 1, 1, 1, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0.13700000000000001, 4.4991772248549271},
    {3, 1, 1, 1, 2, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0.13700000000000001, 0},
    {3, 1, 1, 1, 0, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0, -5.2763074799309351},
    {3, 1, 1, 1, 1, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0, 4.1293766494178499},
    {3, 1, 1, 1, 2, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0, 0},
    {3, 1, 1, 1, 0, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0.13700000000000001, -5.774299500427448},
    {3, 1, 1, 1, 1, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0.13700000000000001, 4.4798996985125559},
    {3, 1, 1, 1, 2, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0.13700000000000001, 0},
    {3, 1, 1, 2, 0, 1.1759566764140661, 1.033843888213875, 0.050000000000000003, 0, 1.7484834986578239},
    {3, 1, 1, 2, 1, 1.1759566764140661, 1.033843888213875, 0.050000000000000003, 0, -0.9174867555593601},
    {3, 1, 1, 2, 2, 1.1759566764140661, 1.033843888213875, 0.050000000000000003, 0, 0},
    {3, 1, 1, 2, 0, 1.1759566764140661, 1.033843888213875, 0.050000000000000003, 0.13700000000000001, 1.8920702871362556},
    {3, 1, 1, 2, 1, 1.1759566764140661, 1.033843888213875, 0.050000000000000003, 0.13700000000000001, -1.0255010339454904},
    {3, 1, 1, 2, 2, 1.1759566764140661, 1.033843888213875, 0.050000000000000003, 0.13700000000000001, 0},
    {3, 1, 1, 2, 0, 1.5793695941545203, 0.65923568248403053, 0.050000000000000003, 0, -1.6587233795125786},
    {3, 1, 1, 2, 1, 1.5793695941545203, 0.65923568248403053, 0.050000000000000003, 0, 6.3093551245149468},
    {3, 1, 1, 2, 2, 1.5793695941545203, 0.65923568248403053, 0.050000000000000003, 0, 0},
    {3, 1, 1, 2, 0, 1.5793695941545203, 0.65923568248403053, 0.050000000000000003, 0.13700000000000001, -1.8094002493792447},
    {3, 1, 1, 2, 1, 1.5793695941545203, 0.65923568248403053, 0.050000000000000003, 0.13700000000000001, 6.799451753905192},
    {3, 1, 1, 2, 2, 1.5793695941545203, 0.65923568248403053, 0.050000000000000003, 0.13700000000000001, 0},
    {3, 1, 1, 2, 0, 1.4491150631448377, 1.4591317739555714, 0.050000000000000003, 0, 1.710656218697111},
    {3, 1, 1, 2, 1, 1.4491150631448377, 1.4591317739555714, 0.050000000000000003, 0, -1.5959257951092976},
    {3, 1, 1, 2, 2, 1.4491150631448377, 1.4591317739555714, 0.050000000000000003, 0, 0},
    {3, 1, 1, 2, 0, 1.4491150631448377, 1.4591317739555714, 0.050000000000000003, 0.13700000000000001, 1.7498686201150786},
    {3, 1, 1, 2, 1, 1.4491150631448377, 1.4591317739555714, 0.050000000000000003, 0.13700000000000001, -1.594271643885532},
    {3, 1, 1, 2, 2, 1.4491150631448377, 1.4591317739555714, 0.050000000000000003, 0.13700000000000001, 0},
    {3, 0, 0, 0, 0, 0.0085011710662745762, 0.14740329063661345, 0.050000000000000003, 0, 0.44711184160474737},
    {3, 0, 0, 0, 1, 0.0085011710662745762, 0.14740329063661345, 0.050000000000000003, 0, 0.20849167579207842},
    {3, 0, 0, 0, 2, 0.0085011710662745762, 0.14740329063661345, 0.050000000000000003, 0, 0},
    {3, 0, 0, 0, 0, 0.0085011710662745762, 0.14740329063661345, 0.050000000000000003, 0.13700000000000001, 0.4767562603185464},
    {3, 0, 0, 0, 1, 0.0085011710662745762, 0.14740329063661345, 0.050000000000000003, 0.13700000000000001, 0.22231509525987625},
    {3, 0, 0, 0, 2, 0.0085011710662745762, 0.14740329063661345, 0.050000000000000003, 0.13700000000000001, 0},
    {3, 0, 0, 0, 0, 0.018932091702127847, 0.69292248910233978, 0.050000000000000003, 0, 0.44711184160474737},
    {3, 0, 0, 0, 1, 0.018932091702127847, 0.69292248910233978, 0.050000000000000003, 0, 0.20849167579207842},
    {3, 0, 0, 0, 2, 0.018932091702127847, 0.69292248910233978, 0.050000000000000003, 0, 0},
    {3, 0, 0, 0, 0, 0.018932091702127847, 0.69292248910233978, 0.050000000000000003, 0.13700000000000001, 0.4767562603185464},
    {3, 0, 0, 0, 1, 0.018932091702127847, 0.69292248910233978, 0.050000000000000003, 0.13700000000000001, 0.22231509525987625},
    {3, 0, 0, 0, 2, 0.018932091702127847, 0.69292248910233978, 0.050000000000000003, 0.13700000000000001, 0},
    {3, 0, 0, 0, 0, 0.62398034980701278, -0.16141813112137662, 0.050000000000000003, 0, 0.44711184160474737},
    {3, 0, 0, 0, 1, 0.62398034980701278, -0.16141813112137662, 0.050000000000000003, 0, 0.20849167579207842},
    {3, 0, 0, 0, 2, 0.62398034980701278, -0.16141813112137662, 0.050000000000000003, 0, 0},
    {3, 0, 0, 0, 0, 0.62398034980701278, -0.16141813112137662, 0.050000000000000003, 0.13700000000000001, 0.4767562603185464},
    {3, 0, 0, 0, 1, 0.62398034980701278, -0.16141813112137662, 0.050000000000000003, 0.13700000000000001, 0.22231509525987625},
    {3, 0, 0, 0, 2, 0.62398034980701278, -0.16141813112137662, 0.050000000000000003, 0.13700000000000001, 0},
    {3, 0, 0, 1, 0, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0, 0.36252311481465999},
    {3, 0, 0, 1, 1, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0, 0.16904730469627979},
    {3, 0, 0, 1, 2, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0, 0},
    {3, 0, 0, 1, 0, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0.13700000000000001, 0.38655912998801056},
    {3, 0, 0, 1, 1, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0.13700000000000001, 0.18025548264314289},
    {3, 0, 0, 1, 2, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0.13700000000000001, 0},
    {3, 0, 0, 1, 0, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0, 0.36252311481465999},
    {3, 0, 0, 1, 1, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0, 0.16904730469627979},
    {3, 0, 0, 1, 2, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0, 0},
    {3, 0, 0, 1, 0, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0.13700000000000001, 0.38655912998801056},
    {3, 0, 0, 1, 1, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0.13700000000000001, 0.18025548264314289},
    {3, 0, 0, 1, 2, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0.13700000000000001, 0},
    {3, 0, 0, 1, 0, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0, 0.36252311481465999},
    {3, 0, 0, 1, 1, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0, 0.16904730469627979},
    {3, 0, 0, 1, 2, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0, 0},
    {3, 0, 0, 1, 0, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0.13700000000000001, 0.38655912998801056},
    {3, 0, 0, 1, 1, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0.13700000000000001, 0.18025548264314289},
    {3, 0, 0, 1, 2, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0.13700000000000001, 0},
    {3, 0, 0, 2, 0, 1.1759566764140661, 1.033843888213875, 0.050000000000000003, 0, 0.29807456106983155},
    {3, 0, 0, 2, 1, 1.1759566764140661, 1.033843888213875, 0.050000000000000003, 0, 0.13899445052805226},
    {3, 0, 0, 2, 2, 1.1759566764140661, 1.033843888213875, 0.050000000000000003, 0, 0},
    {3, 0, 0, 2, 0, 1.1759566764140661, 1.033843888213875, 0.050000000000000003, 0.13700000000000001, 0.3178375068790309},
    {3, 0, 0, 2, 1, 1.1759566764140661, 1.033843888213875, 0.050000000000000003, 0.13700000000000001, 0.14821006350658414},
    {3, 0, 0, 2, 2, 1.1759566764140661, 1.033843888213875, 0.050000000000000003, 0.13700000000000001, 0},
    {3, 0, 0, 2, 0, 1.5793695941545203, 0.65923568248403053, 0.050000000000000003, 0, 0.29807456106983155},
    {3, 0, 0, 2, 1, 1.5793695941545203, 0.65923568248403053, 0.050000000000000003, 0, 0.13899445052805226},
    {3, 0, 0, 2, 2, 1.5793695941545203, 0.65923568248403053, 0.050000000000000003, 0, 0},
    {3, 0, 0, 2, 0, 1.5793695941545203, 0.65923568248403053, 0.050000000000000003, 0.13700000000000001, 0.3178375068790309},
    {3, 0, 0, 2, 1, 1.5793695941545203, 0.65923568248403053, 0.050000000000000003, 0.13700000000000001, 0.14821006350658414},
    {3, 0, 0, 2, 2, 1.5793695941545203, 0.65923568248403053, 0.050000000000000003, 0.13700000000000001, 0},
    {3, 0, 0, 2, 0, 1.4491150631448377, 1.4591317739555714, 0.050000000000000003, 0, 0.29807456106983155},
    {3, 0, 0, 2, 1, 1.4491150631448377, 1.4591317739555714, 0.050000000000000003, 0, 0.13899445052805226},
    {3, 0, 0, 2, 2, 1.4491150631448377, 1.4591317739555714, 0.050000000000000003, 0, 0},
    {3, 0, 0, 2, 0, 1.4491150631448377, 1.4591317739555714, 0.050000000000000003, 0.13700000000000001, 0.3178375068790309},
    {3, 0, 0, 2, 1, 1.4491150631448377, 1.4591317739555714, 0.050000000000000003, 0.13700000000000001, 0.14821006350658414},
    {3, 0, 0, 2, 2, 1.4491150631448377, 1.4591317739555714, 0.050000000000000003, 0.13700000000000001, 0},
    {4, 0, 1, 1, 0, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0, -1.343208854517717},
    {4, 0, 1, 1, 1, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0, -4.1060665746851539},
    {4, 0, 1, 1, 2, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0, 0},
    {4, 0, 1, 1, 0, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0.13700000000000001, -1.4322663162044635},
    {4, 0, 1, 1, 1, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0.13700000000000001, -4.3783070869691132},
    {4, 0, 1, 1, 2, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0.13700000000000001, 0},
    {4, 0, 1, 1, 0, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0, 0.67174565924998675},
    {4, 0, 1, 1, 1, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0, -3.3002024474809875},
    {4, 0, 1, 1, 2, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0, 0},
    {4, 0, 1, 1, 0, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0.13700000000000001, 0.71628375405980249},
    {4, 0, 1, 1, 1, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0.13700000000000001, -3.5190125394756344},
    {4, 0, 1, 1, 2, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0.13700000000000001, 0},
    {4, 0, 1, 1, 0, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0, 5.2196172616052685},
    {4, 0, 1, 1, 1, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0, -3.9012196941208392},
    {4, 0, 1, 1, 2, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0, 0},
    {4, 0, 1, 1, 0, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0.13700000000000001, 5.5656884349238789},
    {4, 0, 1, 1, 1, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0.13700000000000001, -4.1598784442267522},
    {4, 0, 1, 1, 2, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0.13700000000000001, 0},
    {4, 1, 1, 1, 0, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0, -1.539028490119938},
    {4, 1, 1, 1, 1, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0, -5.3707399301151506},
    {4, 1, 1, 1, 2, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0, 0},
    {4, 1, 1, 1, 0, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0.13700000000000001, -1.6410691893996345},
    {4, 1, 1, 1, 1, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0.13700000000000001, -5.7268308417756808},
    {4, 1, 1, 1, 2, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0.13700000000000001, 0},
    {4, 1, 1, 1, 0, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0, 1.021282189600389},
    {4, 1, 1, 1, 1, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0, -4.028628565613591},
    {4, 1, 1, 1, 2, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0, 0},
    {4, 1, 1, 1, 0, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0.13700000000000001, 1.0889952627876187},
    {4, 1, 1, 1, 1, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0.13700000000000001, -4.2957347813934792},
    {4, 1, 1, 1, 2, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0.13700000000000001, 0},
    {4, 1, 1, 1, 0, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0, 6.5103729677110858},
    {4, 1, 1, 1, 1, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0, -4.7455785343229406},
    {4, 1, 1, 1, 2, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0, 0},
    {4, 1, 1, 1, 0, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0.13700000000000001, 6.9420238529686431},
    {4, 1, 1, 1, 1, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0.13700000000000001, -5.0602199819879496},
    {4, 1, 1, 1, 2, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0.13700000000000001, 0},
    {4, 0, 0, 1, 0, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0, -0.36252311481465999},
    {4, 0, 0, 1, 1, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0, -0.16904730469627979},
    {4, 0, 0, 1, 2, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0, 0},
    {4, 0, 0, 1, 0, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0.13700000000000001, -0.38655912998801056},
    {4, 0, 0, 1, 1, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0.13700000000000001, -0.18025548264314289},
    {4, 0, 0, 1, 2, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0.13700000000000001, 0},
    {4, 0, 0, 1, 0, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0, -0.36252311481465999},
    {4, 0, 0, 1, 1, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0, -0.16904730469627979},
    {4, 0, 0, 1, 2, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0, 0},
    {4, 0, 0, 1, 0, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0.13700000000000001, -0.38655912998801056},
    {4, 0, 0, 1, 1, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0.13700000000000001, -0.18025548264314289},
    {4, 0, 0, 1, 2, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0.13700000000000001, 0},
    {4, 0, 0, 1, 0, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0, -0.36252311481465999},
    {4, 0, 0, 1, 1, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0, -0.16904730469627979},
    {4, 0, 0, 1, 2, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0, 0},
    {4, 0, 0, 1, 0, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0.13700000000000001, -0.38655912998801056},
    {4, 0, 0, 1, 1, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0.13700000000000001, -0.18025548264314289},
    {4, 0, 0, 1, 2, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0.13700000000000001, 0},
    {5, 0, 1, 1, 0, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0, 0.3010962558957363},
    {5, 0, 1, 1, 0, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0.13700000000000001, 0.39158675438086121},
    {5, 0, 1, 1, 0, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0, 1.2732502089607542},
    {5, 0, 1, 1, 0, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0.13700000000000001, 1.444509500225043},
    {5, 0, 1, 1, 0, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0, 0.4541835921594759},
    {5, 0, 1, 1, 0, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0.13700000000000001, 0.55518577432539384},
    {5, 1, 1, 1, 0, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0, 0.43836972368733312},
    {5, 1, 1, 1, 0, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0.13700000000000001, 0.53796173051744767},
    {5, 1, 1, 1, 0, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0, 1.2720660777183455},
    {5, 1, 1, 1, 0, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0.13700000000000001, 1.4432468586881269},
    {5, 1, 1, 1, 0, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0, 0.39320746451870653},
    {5, 1, 1, 1, 0, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0.13700000000000001, 0.49016680618896402},
    {5, 0, 0, 1, 0, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0, 1.1963445024848585},
    {5, 0, 0, 1, 0, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0.13700000000000001, 1.3461917690202987},
    {5, 0, 0, 1, 0, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0, 1.2523268970288379},
    {5, 0, 0, 1, 0, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0.13700000000000001, 1.4221989304259162},
    {5, 0, 0, 1, 0, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0, 0.80422040397640526},
    {5, 0, 0, 1, 0, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0.13700000000000001, 0.9284307341920397},
    {6, 0, 1, 1, 0, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0, 1.0012099411567199},
    {6, 0, 1, 1, 0, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0.13700000000000001, 0.88057162810617406},
    {6, 0, 1, 1, 0, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0, 0.90228350740237806},
    {6, 0, 1, 1, 0, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0.13700000000000001, 0.79356509006365705},
    {6, 0, 1, 1, 0, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0, 0.998652526440717},
    {6, 0, 1, 1, 0, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0.13700000000000001, 0.87832236274469389},
    {6, 1, 1, 1, 0, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0, 1.0012099411567199},
    {6, 1, 1, 1, 0, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0.13700000000000001, 0.88057162810617406},
    {6, 1, 1, 1, 0, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0, 0.90228350740237806},
    {6, 1, 1, 1, 0, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0.13700000000000001, 0.79356509006365705},
    {6, 1, 1, 1, 0, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0, 0.998652526440717},
    {6, 1, 1, 1, 0, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0.13700000000000001, 0.87832236274469389},
    {6, 0, 0, 1, 0, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0, 1.0012099411567199},
    {6, 0, 0, 1, 0, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0.13700000000000001, 0.88057162810617406},
    {6, 0, 0, 1, 0, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0, 0.90228350740237806},
    {6, 0, 0, 1, 0, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0.13700000000000001, 0.79356509006365705},
    {6, 0, 0, 1, 0, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0, 0.998652526440717},
    {6, 0, 0, 1, 0, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0.13700000000000001, 0.87832236274469389},
    {7, 0, 1, 0, 0, 0.0085011710662745762, 0.14740329063661345, 0.050000000000000003, 0, 0.15233635545548285},
    {7, 0, 1, 0, 1, 0.0085011710662745762, 0.14740329063661345, 0.050000000000000003, 0, -0.048530464531333722},
    {7, 0, 1, 0, 2, 0.0085011710662745762, 0.14740329063661345, 0.050000000000000003, 0, 0},
    {7, 0, 1, 0, 3, 0.0085011710662745762, 0.14740329063661345, 0.050000000000000003, 0, 0.036663644544517156},
    {7, 0, 1, 0, 4, 0.0085011710662745762, 0.14740329063661345, 0.050000000000000003, 0, 0},
    {7, 0, 1, 0, 5, 0.0085011710662745762, 0.14740329063661345, 0.050000000000000003, 0, 0.019},
    {7, 0, 1, 0, 0, 0.018932091702127847, 0.69292248910233978, 0.050000000000000003, 0, 0.15233635545548285},
    {7, 0, 1, 0, 1, 0.018932091702127847, 0.69292248910233978, 0.050000000000000003, 0, -0.048530464531333722},
    {7, 0, 1, 0, 2, 0.018932091702127847, 0.69292248910233978, 0.050000000000000003, 0, 0},
    {7, 0, 1, 0, 3, 0.018932091702127847, 0.69292248910233978, 0.050000000000000003, 0, 0.036663644544517156},
    {7, 0, 1, 0, 4, 0.018932091702127847, 0.69292248910233978, 0.050000000000000003, 0, 0},
    {7, 0, 1, 0, 5, 0.018932091702127847, 0.69292248910233978, 0.050000000000000003, 0, 0.019},
    {7, 0, 1, 0, 0, 0.62398034980701278, -0.16141813112137662, 0.050000000000000003, 0, 0.15233635545548285},
    {7, 0, 1, 0, 1, 0.62398034980701278, -0.16141813112137662, 0.050000000000000003, 0, -0.048530464531333722},
    {7, 0, 1, 0, 2, 0.62398034980701278, -0.16141813112137662, 0.050000000000000003, 0, 0},
    {7, 0, 1, 0, 3, 0.62398034980701278, -0.16141813112137662, 0.050000000000000003, 0, 0.036663644544517156},
    {7, 0, 1, 0, 4, 0.62398034980701278, -0.16141813112137662, 0.050000000000000003, 0, 0},
    {7, 0, 1, 0, 5, 0.62398034980701278, -0.16141813112137662, 0.050000000000000003, 0, 0.019},
    {7, 0, 1, 1, 0, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0, 0.15280624035403112},
    {7, 0, 1, 1, 1, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0, -0.047964907336280282},
    {7, 0, 1, 1, 2, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0, 0},
    {7, 0, 1, 1, 3, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0, 0.036193759645968888},
    {7, 0, 1, 1, 4, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0, 0},
    {7, 0, 1, 1, 5, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0, 0.019},
    {7, 0, 1, 1, 0, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0, 0.15300830889716557},
    {7, 0, 1, 1, 1, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0, 0.047718212351196242},
    {7, 0, 1, 1, 2, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0, 0},
    {7, 0, 1, 1, 3, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0, 0.035991691102834442},
    {7, 0, 1, 1, 4, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0, 0},
    {7, 0, 1, 1, 5, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0, 0.019},
    {7, 0, 1, 1, 0, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0, 0.037247111281922887},
    {7, 0, 1, 1, 1, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0, 0.049217443385810675},
    {7, 0, 1, 1, 2, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0, 0},
    {7, 0, 1, 1, 3, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0, 0.15175288871807713},
    {7, 0, 1, 1, 4, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0, 0},
    {7, 0, 1, 1, 5, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0, 0.019},
    {7, 0, 1, 2, 0, 1.1759566764140661, 1.033843888213875, 0.050000000000000003, 0, 0.036663644544517163},
    {7, 0, 1, 2, 1, 1.1759566764140661, 1.033843888213875, 0.050000000000000003, 0, 0.048530464531333729},
    {7, 0, 1, 2, 2, 1.1759566764140661, 1.033843888213875, 0.050000000000000003, 0, 0},
    {7, 0, 1, 2, 3, 1.1759566764140661, 1.033843888213875, 0.050000000000000003, 0, 0.15233635545548285},
    {7, 0, 1, 2, 4, 1.1759566764140661, 1.033843888213875, 0.050000000000000003, 0, 0},
    {7, 0, 1, 2, 5, 1.1759566764140661, 1.033843888213875, 0.050000000000000003, 0, 0.019},
    {7, 0, 1, 2, 0, 1.5793695941545203, 0.65923568248403053, 0.050000000000000003, 0, 0.036663644544517163},
    {7, 0, 1, 2, 1, 1.5793695941545203, 0.65923568248403053, 0.050000000000000003, 0, 0.048530464531333729},
    {7, 0, 1, 2, 2, 1.5793695941545203, 0.65923568248403053, 0.050000000000000003, 0, 0},
    {7, 0, 1, 2, 3, 1.5793695941545203, 0.65923568248403053, 0.050000000000000003, 0, 0.15233635545548285},
    {7, 0, 1, 2, 4, 1.5793695941545203, 0.65923568248403053, 0.050000000000000003, 0, 0},
    {7, 0, 1, 2, 5, 1.5793695941545203, 0.65923568248403053, 0.050000000000000003, 0, 0.019},
    {7, 0, 1, 2, 0, 1.4491150631448377, 1.4591317739555714, 0.050000000000000003, 0, 0.036663644544517163},
    {7, 0, 1, 2, 1, 1.4491150631448377, 1.4591317739555714, 0.050000000000000003, 0, 0.048530464531333729},
    {7, 0, 1, 2, 2, 1.4491150631448377, 1.4591317739555714, 0.050000000000000003, 0, 0},
    {7, 0, 1, 2, 3, 1.4491150631448377, 1.4591317739555714, 0.050000000000000003, 0, 0.15233635545548285},
    {7, 0, 1, 2, 4, 1.4491150631448377, 1.4591317739555714, 0.050000000000000003, 0, 0},
    {7, 0, 1, 2, 5, 1.4491150631448377, 1.4591317739555714, 0.050000000000000003, 0, 0.019},
    {7, 1, 1, 0, 0, 0.0085011710662745762, 0.14740329063661345, 0.050000000000000003, 0, 0.15233635545548285},
    {7, 1, 1, 0, 1, 0.0085011710662745762, 0.14740329063661345, 0.050000000000000003, 0, -0.048530464531333722},
    {7, 1, 1, 0, 2, 0.0085011710662745762, 0.14740329063661345, 0.050000000000000003, 0, 0},
    {7, 1, 1, 0, 3, 0.0085011710662745762, 0.14740329063661345, 0.050000000000000003, 0, 0.036663644544517156},
    {7, 1, 1, 0, 4, 0.0085011710662745762, 0.14740329063661345, 0.050000000000000003, 0, 0},
    {7, 1, 1, 0, 5, 0.0085011710662745762, 0.14740329063661345, 0.050000000000000003, 0, 0.019},
    {7, 1, 1, 0, 0, 0.018932091702127847, 0.69292248910233978, 0.050000000000000003, 0, 0.15233635545548285},
    {7, 1, 1, 0, 1, 0.018932091702127847, 0.69292248910233978, 0.050000000000000003, 0, -0.048530464531333722},
    {7, 1, 1, 0, 2, 0.018932091702127847, 0.69292248910233978, 0.050000000000000003, 0, 0},
    {7, 1, 1, 0, 3, 0.018932091702127847, 0.69292248910233978, 0.050000000000000003, 0, 0.036663644544517156},
    {7, 1, 1, 0, 4, 0.018932091702127847, 0.69292248910233978, 0.050000000000000003, 0, 0},
    {7, 1, 1, 0, 5, 0.018932091702127847, 0.69292248910233978, 0.050000000000000003, 0, 0.019},
    {7, 1, 1, 0, 0, 0.62398034980701278, -0.16141813112137662, 0.050000000000000003, 0, 0.15233635545548285},
    {7, 1, 1, 0, 1, 0.62398034980701278, -0.16141813112137662, 0.050000000000000003, 0, -0.048530464531333722},
    {7, 1, 1, 0, 2, 0.62398034980701278, -0.16141813112137662, 0.050000000000000003, 0, 0},
    {7, 1, 1, 0, 3, 0.62398034980701278, -0.16141813112137662, 0.050000000000000003, 0, 0.036663644544517156},
    {7, 1, 1, 0, 4, 0.62398034980701278, -0.16141813112137662, 0.050000000000000003, 0, 0},
    {7, 1, 1, 0, 5, 0.62398034980701278, -0.16141813112137662, 0.050000000000000003, 0, 0.019},
    {7, 1, 1, 1, 0, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0, 0.15280624035403112},
    {7, 1, 1, 1, 1, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0, -0.047964907336280282},
    {7, 1, 1, 1, 2, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0, 0},
    {7, 1, 1, 1, 3, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0, 0.036193759645968888},
    {7, 1, 1, 1, 4, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0, 0},
    {7, 1, 1, 1, 5, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0, 0.019},
    {7, 1, 1, 1, 0, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0, 0.15300830889716557},
    {7, 1, 1, 1, 1, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0, 0.047718212351196242},
    {7, 1, 1, 1, 2, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0, 0},
    {7, 1, 1, 1, 3, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0, 0.035991691102834442},
    {7, 1, 1, 1, 4, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0, 0},
    {7, 1, 1, 1, 5, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0, 0.019},
    {7, 1, 1, 1, 0, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0, 0.037247111281922887},
    {7, 1, 1, 1, 1, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0, 0.049217443385810675},
    {7, 1, 1, 1, 2, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0, 0},
    {7, 1, 1, 1, 3, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0, 0.15175288871807713},
    {7, 1, 1, 1, 4, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0, 0},
    {7, 1, 1, 1, 5, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0, 0.019},
    {7, 1, 1, 2, 0, 1.1759566764140661, 1.033843888213875, 0.050000000000000003, 0, 0.036663644544517163},
    {7, 1, 1, 2, 1, 1.1759566764140661, 1.033843888213875, 0.050000000000000003, 0, 0.048530464531333729},
    {7, 1, 1, 2, 2, 1.1759566764140661, 1.033843888213875, 0.050000000000000003, 0, 0},
    {7, 1, 1, 2, 3, 1.1759566764140661, 1.033843888213875, 0.050000000000000003, 0, 0.15233635545548285},
    {7, 1, 1, 2, 4, 1.1759566764140661, 1.033843888213875, 0.050000000000000003, 0, 0},
    {7, 1, 1, 2, 5, 1.1759566764140661, 1.033843888213875, 0.050000000000000003, 0, 0.019},
    {7, 1, 1, 2, 0, 1.5793695941545203, 0.65923568248403053, 0.050000000000000003, 0, 0.036663644544517163},
    {7, 1, 1, 2, 1, 1.5793695941545203, 0.65923568248403053, 0.050000000000000003, 0, 0.048530464531333729},
    {7, 1, 1, 2, 2, 1.5793695941545203, 0.65923568248403053, 0.050000000000000003, 0, 0},
    {7, 1, 1, 2, 3, 1.5793695941545203, 0.65923568248403053, 0.050000000000000003, 0, 0.15233635545548285},
    {7, 1, 1, 2, 4, 1.5793695941545203, 0.65923568248403053, 0.050000000000000003, 0, 0},
    {7, 1, 1, 2, 5, 1.5793695941545203, 0.65923568248403053, 0.050000000000000003, 0, 0.019},
    {7, 1, 1, 2, 0, 1.4491150631448377, 1.4591317739555714, 0.050000000000000003, 0, 0.036663644544517163},
    {7, 1, 1, 2, 1, 1.4491150631448377, 1.4591317739555714, 0.050000000000000003, 0, 0.048530464531333729},
    {7, 1, 1, 2, 2, 1.4491150631448377, 1.4591317739555714, 0.050000000000000003, 0, 0},
    {7, 1, 1, 2, 3, 1.4491150631448377, 1.4591317739555714, 0.050000000000000003, 0, 0.15233635545548285},
    {7, 1, 1, 2, 4, 1.4491150631448377, 1.4591317739555714, 0.050000000000000003, 0, 0},
    {7, 1, 1, 2, 5, 1.4491150631448377, 1.4591317739555714, 0.050000000000000003, 0, 0.019},
    {7, 0, 0, 0, 0, 0.0085011710662745762, 0.14740329063661345, 0.050000000000000003, 0, 0.15233635545548285},
    {7, 0, 0, 0, 1, 0.0085011710662745762, 0.14740329063661345, 0.050000000000000003, 0, -0.048530464531333722},
    {7, 0, 0, 0, 2, 0.0085011710662745762, 0.14740329063661345, 0.050000000000000003, 0, 0},
    {7, 0, 0, 0, 3, 0.0085011710662745762, 0.14740329063661345, 0.050000000000000003, 0, 0.036663644544517156},
    {7, 0, 0, 0, 4, 0.0085011710662745762, 0.14740329063661345, 0.050000000000000003, 0, 0},
    {7, 0, 0, 0, 5, 0.0085011710662745762, 0.14740329063661345, 0.050000000000000003, 0, 0.019},
    {7, 0, 0, 0, 0, 0.018932091702127847, 0.69292248910233978, 0.050000000000000003, 0, 0.15233635545548285},
    {7, 0, 0, 0, 1, 0.018932091702127847, 0.69292248910233978, 0.050000000000000003, 0, -0.048530464531333722},
    {7, 0, 0, 0, 2, 0.018932091702127847, 0.69292248910233978, 0.050000000000000003, 0, 0},
    {7, 0, 0, 0, 3, 0.018932091702127847, 0.69292248910233978, 0.050000000000000003, 0, 0.036663644544517156},
    {7, 0, 0, 0, 4, 0.018932091702127847, 0.69292248910233978, 0.050000000000000003, 0, 0},
    {7, 0, 0, 0, 5, 0.018932091702127847, 0.69292248910233978, 0.050000000000000003, 0, 0.019},
    {7, 0, 0, 0, 0, 0.62398034980701278, -0.16141813112137662, 0.050000000000000003, 0, 0.15233635545548285},
    {7, 0, 0, 0, 1, 0.62398034980701278, -0.16141813112137662, 0.050000000000000003, 0, -0.048530464531333722},
    {7, 0, 0, 0, 2, 0.62398034980701278, -0.16141813112137662, 0.050000000000000003, 0, 0},
    {7, 0, 0, 0, 3, 0.62398034980701278, -0.16141813112137662, 0.050000000000000003, 0, 0.036663644544517156},
    {7, 0, 0, 0, 4, 0.62398034980701278, -0.16141813112137662, 0.050000000000000003, 0, 0},
    {7, 0, 0, 0, 5, 0.62398034980701278, -0.16141813112137662, 0.050000000000000003, 0, 0.019},
    {7, 0, 0, 1, 0, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0, 0.15233635545548285},
    {7, 0, 0, 1, 1, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0, -0.048530464531333722},
    {7, 0, 0, 1, 2, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0, 0},
    {7, 0, 0, 1, 3, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0, 0.036663644544517156},
    {7, 0, 0, 1, 4, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0, 0},
    {7, 0, 0, 1, 5, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0, 0.019},
    {7, 0, 0, 1, 0, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0, 0.15233635545548285},
    {7, 0, 0, 1, 1, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0, -0.048530464531333722},
    {7, 0, 0, 1, 2, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0, 0},
    {7, 0, 0, 1, 3, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0, 0.036663644544517156},
    {7, 0, 0, 1, 4, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0, 0},
    {7, 0, 0, 1, 5, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0, 0.019},
    {7, 0, 0, 1, 0, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0, 0.15233635545548285},
    {7, 0, 0, 1, 1, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0, -0.048530464531333722},
    {7, 0, 0, 1, 2, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0, 0},
    {7, 0, 0, 1, 3, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0, 0.036663644544517156},
    {7, 0, 0, 1, 4, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0, 0},
    {7, 0, 0, 1, 5, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0, 0.019},
    {7, 0, 0, 2, 0, 1.1759566764140661, 1.033843888213875, 0.050000000000000003, 0, 0.15233635545548285},
    {7, 0, 0, 2, 1, 1.1759566764140661, 1.033843888213875, 0.050000000000000003, 0, -0.048530464531333722},
    {7, 0, 0, 2, 2, 1.1759566764140661, 1.033843888213875, 0.050000000000000003, 0, 0},
    {7, 0, 0, 2, 3, 1.1759566764140661, 1.033843888213875, 0.050000000000000003, 0, 0.036663644544517156},
    {7, 0, 0, 2, 4, 1.1759566764140661, 1.033843888213875, 0.050000000000000003, 0, 0},
    {7, 0, 0, 2, 5, 1.1759566764140661, 1.033843888213875, 0.050000000000000003, 0, 0.019},
    {7, 0, 0, 2, 0, 1.5793695941545203, 0.65923568248403053, 0.050000000000000003, 0, 0.15233635545548285},
    {7, 0, 0, 2, 1, 1.5793695941545203, 0.65923568248403053, 0.050000000000000003, 0, -0.048530464531333722},
    {7, 0, 0, 2, 2, 1.5793695941545203, 0.65923568248403053, 0.050000000000000003, 0, 0},
    {7, 0, 0, 2, 3, 1.5793695941545203, 0.65923568248403053, 0.050000000000000003, 0, 0.036663644544517156},
    {7, 0, 0, 2, 4, 1.5793695941545203, 0.65923568248403053, 0.050000000000000003, 0, 0},
    {7, 0, 0, 2, 5, 1.5793695941545203, 0.65923568248403053, 0.050000000000000003, 0, 0.019},
    {7, 0, 0, 2, 0, 1.4491150631448377, 1.4591317739555714, 0.050000000000000003, 0, 0.15233635545548285},
    {7, 0, 0, 2, 1, 1.4491150631448377, 1.4591317739555714, 0.050000000000000003, 0, -0.048530464531333722},
    {7, 0, 0, 2, 2, 1.4491150631448377, 1.4591317739555714, 0.050000000000000003, 0, 0},
    {7, 0, 0, 2, 3, 1.4491150631448377, 1.4591317739555714, 0.050000000000000003, 0, 0.036663644544517156},
    {7, 0, 0, 2, 4, 1.4491150631448377, 1.4591317739555714, 0.050000000000000003, 0, 0},
    {7, 0, 0, 2, 5, 1.4491150631448377, 1.4591317739555714, 0.050000000000000003, 0, 0.019},
    {8, 0, 1, 0, 0, 0.0085011710662745762, 0.14740329063661345, 0.050000000000000003, 0, 0.56151111077974447},
    {8, 0, 1, 0, 1, 0.0085011710662745762, 0.14740329063661345, 0.050000000000000003, 0, -0.16069690242163481},
    {8, 0, 1, 0, 2, 0.0085011710662745762, 0.14740329063661345, 0.050000000000000003, 0, 0},
    {8, 0, 1, 0, 3, 0.0085011710662745762, 0.14740329063661345, 0.050000000000000003, 0, 0.17848888922025549},
    {8, 0, 1, 0, 4, 0.0085011710662745762, 0.14740329063661345, 0.050000000000000003, 0, 0},
    {8, 0, 1, 0, 5, 0.0085011710662745762, 0.14740329063661345, 0.050000000000000003, 0, 0.12},
    {8, 0, 1, 0, 0, 0.018932091702127847, 0.69292248910233978, 0.050000000000000003, 0, 0.56151111077974447},
    {8, 0, 1, 0, 1, 0.018932091702127847, 0.69292248910233978, 0.050000000000000003, 0, -0.16069690242163481},
    {8, 0, 1, 0, 2, 0.018932091702127847, 0.69292248910233978, 0.050000000000000003, 0, 0},
    {8, 0, 1, 0, 3, 0.018932091702127847, 0.69292248910233978, 0.050000000000000003, 0, 0.17848888922025549},
    {8, 0, 1, 0, 4, 0.018932091702127847, 0.69292248910233978, 0.050000000000000003, 0, 0},
    {8, 0, 1, 0, 5, 0.018932091702127847, 0.69292248910233978, 0.050000000000000003, 0, 0.12},
    {8, 0, 1, 0, 0, 0.62398034980701278, -0.16141813112137662, 0.050000000000000003, 0, 0.56151111077974447},
    {8, 0, 1, 0, 1, 0.62398034980701278, -0.16141813112137662, 0.050000000000000003, 0, -0.16069690242163481},
    {8, 0, 1, 0, 2, 0.62398034980701278, -0.16141813112137662, 0.050000000000000003, 0, 0},
    {8, 0, 1, 0, 3, 0.62398034980701278, -0.16141813112137662, 0.050000000000000003, 0, 0.17848888922025549},
    {8, 0, 1, 0, 4, 0.62398034980701278, -0.16141813112137662, 0.050000000000000003, 0, 0},
    {8, 0, 1, 0, 5, 0.62398034980701278, -0.16141813112137662, 0.050000000000000003, 0, 0.12},
    {8, 0, 1, 1, 0, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0, 0.56306702103983808},
    {8, 0, 1, 1, 1, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0, -0.15882419647774926},
    {8, 0, 1, 1, 2, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0, 0},
    {8, 0, 1, 1, 3, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0, 0.17693297896016186},
    {8, 0, 1, 1, 4, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0, 0},
    {8, 0, 1, 1, 5, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0, 0.12},
    {8, 0, 1, 1, 0, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0, 0.56373612217604485},
    {8, 0, 1, 1, 1, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0, 0.1580073256662127},
    {8, 0, 1, 1, 2, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0, 0},
    {8, 0, 1, 1, 3, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0, 0.17626387782395508},
    {8, 0, 1, 1, 4, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0, 0},
    {8, 0, 1, 1, 5, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0, 0.12},
    {8, 0, 1, 1, 0, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0, 0.18042089828451288},
    {8, 0, 1, 1, 1, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0, 0.16297166684043268},
    {8, 0, 1, 1, 2, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0, 0},
    {8, 0, 1, 1, 3, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0, 0.55957910171548708},
    {8, 0, 1, 1, 4, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0, 0},
    {8, 0, 1, 1, 5, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0, 0.12},
    {8, 0, 1, 2, 0, 1.1759566764140661, 1.033843888213875, 0.050000000000000003, 0, 0.17848888922025549},
    {8, 0, 1, 2, 1, 1.1759566764140661, 1.033843888213875, 0.050000000000000003, 0, 0.16069690242163484},
    {8, 0, 1, 2, 2, 1.1759566764140661, 1.033843888213875, 0.050000000000000003, 0, 0},
    {8, 0, 1, 2, 3, 1.1759566764140661, 1.033843888213875, 0.050000000000000003, 0, 0.56151111077974447},
    {8, 0, 1, 2, 4, 1.1759566764140661, 1.033843888213875, 0.050000000000000003, 0, 0},
    {8, 0, 1, 2, 5, 1.1759566764140661, 1.033843888213875, 0.050000000000000003, 0, 0.12},
    {8, 0, 1, 2, 0, 1.5793695941545203, 0.65923568248403053, 0.050000000000000003, 0, 0.17848888922025549},
    {8, 0, 1, 2, 1, 1.5793695941545203, 0.65923568248403053, 0.050000000000000003, 0, 0.16069690242163484},
    {8, 0, 1, 2, 2, 1.5793695941545203, 0.65923568248403053, 0.050000000000000003, 0, 0},
    {8, 0, 1, 2, 3, 1.5793695941545203, 0.65923568248403053, 0.050000000000000003, 0, 0.56151111077974447},
    {8, 0, 1, 2, 4, 1.5793695941545203, 0.65923568248403053, 0.050000000000000003, 0, 0},
    {8, 0, 1, 2, 5, 1.5793695941545203, 0.65923568248403053, 0.050000000000000003, 0, 0.12},
    {8, 0, 1, 2, 0, 1.4491150631448377, 1.4591317739555714, 0.050000000000000003, 0, 0.17848888922025549},
    {8, 0, 1, 2, 1, 1.4491150631448377, 1.4591317739555714, 0.050000000000000003, 0, 0.16069690242163484},
    {8, 0, 1, 2, 2, 1.4491150631448377, 1.4591317739555714, 0.050000000000000003, 0, 0},
    {8, 0, 1, 2, 3, 1.4491150631448377, 1.4591317739555714, 0.050000000000000003, 0, 0.56151111077974447},
    {8, 0, 1, 2, 4, 1.4491150631448377, 1.4591317739555714, 0.050000000000000003, 0, 0},
    {8, 0, 1, 2, 5, 1.4491150631448377, 1.4591317739555714, 0.050000000000000003, 0, 0.12},
    {8, 1, 1, 0, 0, 0.0085011710662745762, 0.14740329063661345, 0.050000000000000003, 0, 0.56151111077974447},
    {8, 1, 1, 0, 1, 0.0085011710662745762, 0.14740329063661345, 0.050000000000000003, 0, -0.16069690242163481},
    {8, 1, 1, 0, 2, 0.0085011710662745762, 0.14740329063661345, 0.050000000000000003, 0, 0},
    {8, 1, 1, 0, 3, 0.0085011710662745762, 0.14740329063661345, 0.050000000000000003, 0, 0.17848888922025549},
    {8, 1, 1, 0, 4, 0.0085011710662745762, 0.14740329063661345, 0.050000000000000003, 0, 0},
    {8, 1, 1, 0, 5, 0.0085011710662745762, 0.14740329063661345, 0.050000000000000003, 0, 0.12},
    {8, 1, 1, 0, 0, 0.018932091702127847, 0.69292248910233978, 0.050000000000000003, 0, 0.56151111077974447},
    {8, 1, 1, 0, 1, 0.018932091702127847, 0.69292248910233978, 0.050000000000000003, 0, -0.16069690242163481},
    {8, 1, 1, 0, 2, 0.018932091702127847, 0.69292248910233978, 0.050000000000000003, 0, 0},
    {8, 1, 1, 0, 3, 0.018932091702127847, 0.69292248910233978, 0.050000000000000003, 0, 0.17848888922025549},
    {8, 1, 1, 0, 4, 0.018932091702127847, 0.69292248910233978, 0.050000000000000003, 0, 0},
    {8, 1, 1, 0, 5, 0.018932091702127847, 0.69292248910233978, 0.050000000000000003, 0, 0.12},
    {8, 1, 1, 0, 0, 0.62398034980701278, -0.16141813112137662, 0.050000000000000003, 0, 0.56151111077974447},
    {8, 1, 1, 0, 1, 0.62398034980701278, -0.16141813112137662, 0.050000000000000003, 0, -0.16069690242163481},
    {8, 1, 1, 0, 2, 0.62398034980701278, -0.16141813112137662, 0.050000000000000003, 0, 0},
    {8, 1, 1, 0, 3, 0.62398034980701278, -0.16141813112137662, 0.050000000000000003, 0, 0.17848888922025549},
    {8, 1, 1, 0, 4, 0.62398034980701278, -0.16141813112137662, 0.050000000000000003, 0, 0},
    {8, 1, 1, 0, 5, 0.62398034980701278, -0.16141813112137662, 0.050000000000000003, 0, 0.12},
    {8, 1, 1, 1, 0, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0, 0.56306702103983808},
    {8, 1, 1, 1, 1, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0, -0.15882419647774926},
    {8, 1, 1, 1, 2, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0, 0},
    {8, 1, 1, 1, 3, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0, 0.17693297896016186},
    {8, 1, 1, 1, 4, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0, 0},
    {8, 1, 1, 1, 5, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0, 0.12},
    {8, 1, 1, 1, 0, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0, 0.56373612217604485},
    {8, 1, 1, 1, 1, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0, 0.1580073256662127},
    {8, 1, 1, 1, 2, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0, 0},
    {8, 1, 1, 1, 3, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0, 0.17626387782395508},
    {8, 1, 1, 1, 4, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0, 0},
    {8, 1, 1, 1, 5, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0, 0.12},
    {8, 1, 1, 1, 0, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0, 0.18042089828451288},
    {8, 1, 1, 1, 1, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0, 0.16297166684043268},
    {8, 1, 1, 1, 2, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0, 0},
    {8, 1, 1, 1, 3, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0, 0.55957910171548708},
    {8, 1, 1, 1, 4, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0, 0},
    {8, 1, 1, 1, 5, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0, 0.12},
    {8, 1, 1, 2, 0, 1.1759566764140661, 1.033843888213875, 0.050000000000000003, 0, 0.17848888922025549},
    {8, 1, 1, 2, 1, 1.1759566764140661, 1.033843888213875, 0.050000000000000003, 0, 0.16069690242163484},
    {8, 1, 1, 2, 2, 1.1759566764140661, 1.033843888213875, 0.050000000000000003, 0, 0},
    {8, 1, 1, 2, 3, 1.1759566764140661, 1.033843888213875, 0.050000000000000003, 0, 0.56151111077974447},
    {8, 1, 1, 2, 4, 1.1759566764140661, 1.033843888213875, 0.050000000000000003, 0, 0},
    {8, 1, 1, 2, 5, 1.1759566764140661, 1.033843888213875, 0.050000000000000003, 0, 0.12},
    {8, 1, 1, 2, 0, 1.5793695941545203, 0.65923568248403053, 0.050000000000000003, 0, 0.17848888922025549},
    {8, 1, 1, 2, 1, 1.5793695941545203, 0.65923568248403053, 0.050000000000000003, 0, 0.16069690242163484},
    {8, 1, 1, 2, 2, 1.5793695941545203, 0.65923568248403053, 0.050000000000000003, 0, 0},
    {8, 1, 1, 2, 3, 1.5793695941545203, 0.65923568248403053, 0.050000000000000003, 0, 0.56151111077974447},
    {8, 1, 1, 2, 4, 1.5793695941545203, 0.65923568248403053, 0.050000000000000003, 0, 0},
    {8, 1, 1, 2, 5, 1.5793695941545203, 0.65923568248403053, 0.050000000000000003, 0, 0.12},
    {8, 1, 1, 2, 0, 1.4491150631448377, 1.4591317739555714, 0.050000000000000003, 0, 0.17848888922025549},
    {8, 1, 1, 2, 1, 1.4491150631448377, 1.4591317739555714, 0.050000000000000003, 0, 0.16069690242163484},
    {8, 1, 1, 2, 2, 1.4491150631448377, 1.4591317739555714, 0.050000000000000003, 0, 0},
    {8, 1, 1, 2, 3, 1.4491150631448377, 1.4591317739555714, 0.050000000000000003, 0, 0.56151111077974447},
    {8, 1, 1, 2, 4, 1.4491150631448377, 1.4591317739555714, 0.050000000000000003, 0, 0},
    {8, 1, 1, 2, 5, 1.4491150631448377, 1.4591317739555714, 0.050000000000000003, 0, 0.12},
    {8, 0, 0, 0, 0, 0.0085011710662745762, 0.14740329063661345, 0.050000000000000003, 0, 0.56151111077974447},
    {8, 0, 0, 0, 1, 0.0085011710662745762, 0.14740329063661345, 0.050000000000000003, 0, -0.16069690242163481},
    {8, 0, 0, 0, 2, 0.0085011710662745762, 0.14740329063661345, 0.050000000000000003, 0, 0},
    {8, 0, 0, 0, 3, 0.0085011710662745762, 0.14740329063661345, 0.050000000000000003, 0, 0.17848888922025549},
    {8, 0, 0, 0, 4, 0.0085011710662745762, 0.14740329063661345, 0.050000000000000003, 0, 0},
    {8, 0, 0, 0, 5, 0.0085011710662745762, 0.14740329063661345, 0.050000000000000003, 0, 0.12},
    {8, 0, 0, 0, 0, 0.018932091702127847, 0.69292248910233978, 0.050000000000000003, 0, 0.56151111077974447},
    {8, 0, 0, 0, 1, 0.018932091702127847, 0.69292248910233978, 0.050000000000000003, 0, -0.16069690242163481},
    {8, 0, 0, 0, 2, 0.018932091702127847, 0.69292248910233978, 0.050000000000000003, 0, 0},
    {8, 0, 0, 0, 3, 0.018932091702127847, 0.69292248910233978, 0.050000000000000003, 0, 0.17848888922025549},
    {8, 0, 0, 0, 4, 0.018932091702127847, 0.69292248910233978, 0.050000000000000003, 0, 0},
    {8, 0, 0, 0, 5, 0.018932091702127847, 0.69292248910233978, 0.050000000000000003, 0, 0.12},
    {8, 0, 0, 0, 0, 0.62398034980701278, -0.16141813112137662, 0.050000000000000003, 0, 0.56151111077974447},
    {8, 0, 0, 0, 1, 0.62398034980701278, -0.16141813112137662, 0.050000000000000003, 0, -0.16069690242163481},
    {8, 0, 0, 0, 2, 0.62398034980701278, -0.16141813112137662, 0.050000000000000003, 0, 0},
    {8, 0, 0, 0, 3, 0.62398034980701278, -0.16141813112137662, 0.050000000000000003, 0, 0.17848888922025549},
    {8, 0, 0, 0, 4, 0.62398034980701278, -0.16141813112137662, 0.050000000000000003, 0, 0},
    {8, 0, 0, 0, 5, 0.62398034980701278, -0.16141813112137662, 0.050000000000000003, 0, 0.12},
    {8, 0, 0, 1, 0, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0, 0.56151111077974447},
    {8, 0, 0, 1, 1, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0, -0.16069690242163481},
    {8, 0, 0, 1, 2, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0, 0},
    {8, 0, 0, 1, 3, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0, 0.17848888922025549},
    {8, 0, 0, 1, 4, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0, 0},
    {8, 0, 0, 1, 5, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0, 0.12},
    {8, 0, 0, 1, 0, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0, 0.56151111077974447},
    {8, 0, 0, 1, 1, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0, -0.16069690242163481},
    {8, 0, 0, 1, 2, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0, 0},
    {8, 0, 0, 1, 3, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0, 0.17848888922025549},
    {8, 0, 0, 1, 4, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0, 0},
    {8, 0, 0, 1, 5, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0, 0.12},
    {8, 0, 0, 1, 0, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0, 0.56151111077974447},
    {8, 0, 0, 1, 1, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0, -0.16069690242163481},
    {8, 0, 0, 1, 2, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0, 0},
    {8, 0, 0, 1, 3, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0, 0.17848888922025549},
    {8, 0, 0, 1, 4, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0, 0},
    {8, 0, 0, 1, 5, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0, 0.12},
    {8, 0, 0, 2, 0, 1.1759566764140661, 1.033843888213875, 0.050000000000000003, 0, 0.56151111077974447},
    {8, 0, 0, 2, 1, 1.1759566764140661, 1.033843888213875, 0.050000000000000003, 0, -0.16069690242163481},
    {8, 0, 0, 2, 2, 1.1759566764140661, 1.033843888213875, 0.050000000000000003, 0, 0},
    {8, 0, 0, 2, 3, 1.1759566764140661, 1.033843888213875, 0.050000000000000003, 0, 0.17848888922025549},
    {8, 0, 0, 2, 4, 1.1759566764140661, 1.033843888213875, 0.050000000000000003, 0, 0},
    {8, 0, 0, 2, 5, 1.1759566764140661, 1.033843888213875, 0.050000000000000003, 0, 0.12},
    {8, 0, 0, 2, 0, 1.5793695941545203, 0.65923568248403053, 0.050000000000000003, 0, 0.56151111077974447},
    {8, 0, 0, 2, 1, 1.5793695941545203, 0.65923568248403053, 0.050000000000000003, 0, -0.16069690242163481},
    {8, 0, 0, 2, 2, 1.5793695941545203, 0.65923568248403053, 0.050000000000000003, 0, 0},
    {8, 0, 0, 2, 3, 1.5793695941545203, 0.65923568248403053, 0.050000000000000003, 0, 0.17848888922025549},
    {8, 0, 0, 2, 4, 1.5793695941545203, 0.65923568248403053, 0.050000000000000003, 0, 0},
    {8, 0, 0, 2, 5, 1.5793695941545203, 0.65923568248403053, 0.050000000000000003, 0, 0.12},
    {8, 0, 0, 2, 0, 1.4491150631448377, 1.4591317739555714, 0.050000000000000003, 0, 0.56151111077974447},
    {8, 0, 0, 2, 1, 1.4491150631448377, 1.4591317739555714, 0.050000000000000003, 0, -0.16069690242163481},
    {8, 0, 0, 2, 2, 1.4491150631448377, 1.4591317739555714, 0.050000000000000003, 0, 0},
    {8, 0, 0, 2, 3, 1.4491150631448377, 1.4591317739555714, 0.050000000000000003, 0, 0.17848888922025549},
    {8, 0, 0, 2, 4, 1.4491150631448377, 1.4591317739555714, 0.050000000000000003, 0, 0},
    {8, 0, 0, 2, 5, 1.4491150631448377, 1.4591317739555714, 0.050000000000000003, 0, 0.12},
    {9, 0, 1, 0, 0, 0.0085011710662745762, 0.14740329063661345, 0.050000000000000003, 0, 0.35304951001798424},
    {9, 0, 1, 0, 0, 0.018932091702127847, 0.69292248910233978, 0.050000000000000003, 0, 0.23983289956701823},
    {9, 0, 1, 0, 0, 0.62398034980701278, -0.16141813112137662, 0.050000000000000003, 0, 0.26757770275892379},
    {9, 0, 1, 2, 0, 1.1759566764140661, 1.033843888213875, 0.050000000000000003, 0, 0.48285019435548782},
    {9, 0, 1, 2, 0, 1.5793695941545203, 0.65923568248403053, 0.050000000000000003, 0, 0.43023798310666783},
    {9, 0, 1, 2, 0, 1.4491150631448377, 1.4591317739555714, 0.050000000000000003, 0, 0.38138563484594756},
    {9, 1, 1, 0, 0, 0.0085011710662745762, 0.14740329063661345, 0.050000000000000003, 0, 0.35304951001798424},
    {9, 1, 1, 0, 0, 0.018932091702127847, 0.69292248910233978, 0.050000000000000003, 0, 0.23983289956701823},
    {9, 1, 1, 0, 0, 0.62398034980701278, -0.16141813112137662, 0.050000000000000003, 0, 0.26757770275892379},
    {9, 1, 1, 2, 0, 1.1759566764140661, 1.033843888213875, 0.050000000000000003, 0, 0.48285019435548782},
    {9, 1, 1, 2, 0, 1.5793695941545203, 0.65923568248403053, 0.050000000000000003, 0, 0.43023798310666783},
    {9, 1, 1, 2, 0, 1.4491150631448377, 1.4591317739555714, 0.050000000000000003, 0, 0.38138563484594756},
    {9, 0, 0, 0, 0, 0.0085011710662745762, 0.14740329063661345, 0.050000000000000003, 0, 0.29999999999999999},
    {9, 0, 0, 0, 0, 0.018932091702127847, 0.69292248910233978, 0.050000000000000003, 0, 0.29999999999999999},
    {9, 0, 0, 0, 0, 0.62398034980701278, -0.16141813112137662, 0.050000000000000003, 0, 0.29999999999999999},
    {9, 0, 0, 2, 0, 1.1759566764140661, 1.033843888213875, 0.050000000000000003, 0, 0.45000000000000001},
    {9, 0, 0, 2, 0, 1.5793695941545203, 0.65923568248403053, 0.050000000000000003, 0, 0.45000000000000001},
    {9, 0, 0, 2, 0, 1.4491150631448377, 1.4591317739555714, 0.050000000000000003, 0, 0.45000000000000001},
    {10, 0, 1, 1, 0, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0, 1.233717755146869},
    {10, 0, 1, 1, 0, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0.13700000000000001, 1.3155157356632357},
    {10, 0, 1, 1, 0, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0, 0.40783428881189554},
    {10, 0, 1, 1, 0, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0.13700000000000001, 0.43487452639538587},
    {10, 0, 1, 1, 0, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0, 0.66691934289889576},
    {10, 0, 1, 1, 0, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0.13700000000000001, 0.71113744318062322},
    {10, 1, 1, 1, 0, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0, 1.5452079360158351},
    {10, 1, 1, 1, 0, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0.13700000000000001, 1.7462799123731461},
    {10, 1, 1, 1, 0, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0, 0.46059997025854699},
    {10, 1, 1, 1, 0, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0.13700000000000001, 0.54639443814442901},
    {10, 1, 1, 1, 0, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0, 0.92475062207851244},
    {10, 1, 1, 1, 0, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0.13700000000000001, 0.90167187140555738},
    {10, 0, 0, 1, 0, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0, 0},
    {10, 0, 0, 1, 0, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0.13700000000000001, 0},
    {10, 0, 0, 1, 0, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0, 0},
    {10, 0, 0, 1, 0, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0.13700000000000001, 0},
    {10, 0, 0, 1, 0, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0, 0},
    {10, 0, 0, 1, 0, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0.13700000000000001, 0},
    {11, 0, 1, 0, 0, 0.0085011710662745762, 0.14740329063661345, 0.050000000000000003, 0, -0.9695254596600601},
    {11, 0, 1, 0, 0, 0.0085011710662745762, 0.14740329063661345, 0.050000000000000003, 0.13700000000000001, -1.0338069570516204},
    {11, 0, 1, 0, 0, 0.018932091702127847, 0.69292248910233978, 0.050000000000000003, 0, 6.2994940505920987},
    {11, 0, 1, 0, 0, 0.018932091702127847, 0.69292248910233978, 0.050000000000000003, 0.13700000000000001, 6.7171632374572559},
    {11, 0, 1, 0, 0, 0.62398034980701278, -0.16141813112137662, 0.050000000000000003, 0, 3.4796317137789265},
    {11, 0, 1, 0, 0, 0.62398034980701278, -0.16141813112137662, 0.050000000000000003, 0.13700000000000001, 3.7103383287566256},
    {11, 0, 1, 1, 0, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0, -1.0408829139677929},
    {11, 0, 1, 1, 0, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0.13700000000000001, -1.1098955547937495},
    {11, 0, 1, 1, 0, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0, 11.25304135132841},
    {11, 0, 1, 1, 0, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0.13700000000000001, 11.999140735377761},
    {11, 0, 1, 1, 0, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0, -7.2683948207926585},
    {11, 0, 1, 1, 0, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0.13700000000000001, -7.7503040868757109},
    {11, 0, 1, 2, 0, 1.1759566764140661, 1.033843888213875, 0.050000000000000003, 0, -22.21738856609667},
    {11, 0, 1, 2, 0, 1.1759566764140661, 1.033843888213875, 0.050000000000000003, 0.13700000000000001, -23.690446329489056},
    {11, 0, 1, 2, 0, 1.5793695941545203, 0.65923568248403053, 0.050000000000000003, 0, -9.6593727521556421},
    {11, 0, 1, 2, 0, 1.5793695941545203, 0.65923568248403053, 0.050000000000000003, 0.13700000000000001, -10.299808687267147},
    {11, 0, 1, 2, 0, 1.4491150631448377, 1.4591317739555714, 0.050000000000000003, 0, 30.78656757827892},
    {11, 0, 1, 2, 0, 1.4491150631448377, 1.4591317739555714, 0.050000000000000003, 0.13700000000000001, 32.827779228535242},
    {11, 1, 1, 0, 0, 0.0085011710662745762, 0.14740329063661345, 0.050000000000000003, 0, 0.34073129482267428},
    {11, 1, 1, 0, 0, 0.0085011710662745762, 0.14740329063661345, 0.050000000000000003, 0.13700000000000001, 0.030074037928662679},
    {11, 1, 1, 0, 0, 0.018932091702127847, 0.69292248910233978, 0.050000000000000003, 0, 6.6238077093447743},
    {11, 1, 1, 0, 0, 0.018932091702127847, 0.69292248910233978, 0.050000000000000003, 0.13700000000000001, 7.2169328955611283},
    {11, 1, 1, 0, 0, 0.62398034980701278, -0.16141813112137662, 0.050000000000000003, 0, 2.9331553063058324},
    {11, 1, 1, 0, 0, 0.62398034980701278, -0.16141813112137662, 0.050000000000000003, 0.13700000000000001, 3.3104042034963848},
    {11, 1, 1, 1, 0, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0, -1.1720455086079518},
    {11, 1, 1, 1, 0, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0.13700000000000001, -0.74101852313581984},
    {11, 1, 1, 1, 0, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0, 12.651987868234501},
    {11, 1, 1, 1, 0, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0.13700000000000001, 13.719415445699161},
    {11, 1, 1, 1, 0, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0, -8.5434878351678858},
    {11, 1, 1, 1, 0, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0.13700000000000001, -9.5666691618098145},
    {11, 1, 1, 2, 0, 1.1759566764140661, 1.033843888213875, 0.050000000000000003, 0, -28.339695092249023},
    {11, 1, 1, 2, 0, 1.1759566764140661, 1.033843888213875, 0.050000000000000003, 0.13700000000000001, -29.774541448562172},
    {11, 1, 1, 2, 0, 1.5793695941545203, 0.65923568248403053, 0.050000000000000003, 0, -15.86816533466104},
    {11, 1, 1, 2, 0, 1.5793695941545203, 0.65923568248403053, 0.050000000000000003, 0.13700000000000001, -17.245302241668394},
    {11, 1, 1, 2, 0, 1.4491150631448377, 1.4591317739555714, 0.050000000000000003, 0, 35.751147443319219},
    {11, 1, 1, 2, 0, 1.4491150631448377, 1.4591317739555714, 0.050000000000000003, 0.13700000000000001, 38.311592804073875},
    {11, 0, 0, 0, 0, 0.0085011710662745762, 0.14740329063661345, 0.050000000000000003, 0, 0},
    {11, 0, 0, 0, 0, 0.0085011710662745762, 0.14740329063661345, 0.050000000000000003, 0.13700000000000001, 0},
    {11, 0, 0, 0, 0, 0.018932091702127847, 0.69292248910233978, 0.050000000000000003, 0, 0},
    {11, 0, 0, 0, 0, 0.018932091702127847, 0.69292248910233978, 0.050000000000000003, 0.13700000000000001, 0},
    {11, 0, 0, 0, 0, 0.62398034980701278, -0.16141813112137662, 0.050000000000000003, 0, 0},
    {11, 0, 0, 0, 0, 0.62398034980701278, -0.16141813112137662, 0.050000000000000003, 0.13700000000000001, 0},
    {11, 0, 0, 1, 0, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0, 0},
    {11, 0, 0, 1, 0, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0.13700000000000001, 0},
    {11, 0, 0, 1, 0, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0, 0},
    {11, 0, 0, 1, 0, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0.13700000000000001, 0},
    {11, 0, 0, 1, 0, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0, 0},
    {11, 0, 0, 1, 0, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0.13700000000000001, 0},
    {11, 0, 0, 2, 0, 1.1759566764140661, 1.033843888213875, 0.050000000000000003, 0, 0},
    {11, 0, 0, 2, 0, 1.1759566764140661, 1.033843888213875, 0.050000000000000003, 0.13700000000000001, 0},
    {11, 0, 0, 2, 0, 1.5793695941545203, 0.65923568248403053, 0.050000000000000003, 0, 0},
    {11, 0, 0, 2, 0, 1.5793695941545203, 0.65923568248403053, 0.050000000000000003, 0.13700000000000001, 0},
    {11, 0, 0, 2, 0, 1.4491150631448377, 1.4591317739555714, 0.050000000000000003, 0, 0},
    {11, 0, 0, 2, 0, 1.4491150631448377, 1.4591317739555714, 0.050000000000000003, 0.13700000000000001, 0},
    {12, 0, 1, 1, 0, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0, 0.19097508886596129},
    {12, 0, 1, 1, 1, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0, 0.10842310718317001},
    {12, 0, 1, 1, 2, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0, 0},
    {12, 0, 1, 1, 0, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0.13700000000000001, 0.20363712321944241},
    {12, 0, 1, 1, 1, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0.13700000000000001, 0.11561177831308909},
    {12, 0, 1, 1, 2, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0.13700000000000001, 0},
    {12, 0, 1, 1, 0, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0, -0.20900201311205002},
    {12, 0, 1, 1, 1, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0, 0.40746354352815417},
    {12, 0, 1, 1, 2, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0, 0},
    {12, 0, 1, 1, 0, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0.13700000000000001, -0.22285926897555641},
    {12, 0, 1, 1, 1, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0.13700000000000001, 0.43447919995005424},
    {12, 0, 1, 1, 2, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0.13700000000000001, 0},
    {12, 0, 1, 1, 0, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0, -0.23879613482951934},
    {12, 0, 1, 1, 1, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0, 1.2519751520915314},
    {12, 0, 1, 1, 2, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0, 0},
    {12, 0, 1, 1, 0, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0.13700000000000001, -0.25462880117697201},
    {12, 0, 1, 1, 1, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0.13700000000000001, 1.3349836349236253},
    {12, 0, 1, 1, 2, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0.13700000000000001, 0},
    {12, 1, 1, 1, 0, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0, 0.12205702057420495},
    {12, 1, 1, 1, 1, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0, 0.19871390347461138},
    {12, 1, 1, 1, 2, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0, 0},
    {12, 1, 1, 1, 0, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0.13700000000000001, 0.13595799601944777},
    {12, 1, 1, 1, 1, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0.13700000000000001, 0.18648712485405591},
    {12, 1, 1, 1, 2, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0.13700000000000001, 0},
    {12, 1, 1, 1, 0, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0, -0.25917575823920486},
    {12, 1, 1, 1, 1, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0, 0.4686254930370663},
    {12, 1, 1, 1, 2, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0, 0},
    {12, 1, 1, 1, 0, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0.13700000000000001, -0.28574675856225001},
    {12, 1, 1, 1, 1, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0.13700000000000001, 0.51369549155937011},
    {12, 1, 1, 1, 2, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0.13700000000000001, 0},
    {12, 1, 1, 1, 0, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0, -0.27898473958684572},
    {12, 1, 1, 1, 1, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0, 1.4508242513591609},
    {12, 1, 1, 1, 2, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0, 0},
    {12, 1, 1, 1, 0, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0.13700000000000001, -0.31170758168639101},
    {12, 1, 1, 1, 1, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0.13700000000000001, 1.5658110346485985},
    {12, 1, 1, 1, 2, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0.13700000000000001, 0},
    {12, 0, 0, 1, 0, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0, 0.17639537865549415},
    {12, 0, 0, 1, 1, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0, -0.028083275966041481},
    {12, 0, 0, 1, 2, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0, 0},
    {12, 0, 0, 1, 0, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0.13700000000000001, 0.18809074875634962},
    {12, 0, 0, 1, 1, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0.13700000000000001, -0.029945253919039769},
    {12, 0, 0, 1, 2, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0.13700000000000001, 0},
    {12, 0, 0, 1, 0, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0, 0.17639537865549415},
    {12, 0, 0, 1, 1, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0, -0.028083275966041481},
    {12, 0, 0, 1, 2, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0, 0},
    {12, 0, 0, 1, 0, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0.13700000000000001, 0.18809074875634962},
    {12, 0, 0, 1, 1, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0.13700000000000001, -0.029945253919039769},
    {12, 0, 0, 1, 2, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0.13700000000000001, 0},
    {12, 0, 0, 1, 0, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0, 0.17639537865549415},
    {12, 0, 0, 1, 1, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0, -0.028083275966041481},
    {12, 0, 0, 1, 2, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0, 0},
    {12, 0, 0, 1, 0, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0.13700000000000001, 0.18809074875634962},
    {12, 0, 0, 1, 1, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0.13700000000000001, -0.029945253919039769},
    {12, 0, 0, 1, 2, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0.13700000000000001, 0},
    {13, 0, 1, 1, 0, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0, 0.99758450248485853},
    {13, 0, 1, 1, 0, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0, 1.2283268970288379},
    {13, 0, 1, 1, 0, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0, 1.0027004039764054},
    {13, 1, 1, 1, 0, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0, 0.99758450248485853},
    {13, 1, 1, 1, 0, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0, 1.2283268970288379},
    {13, 1, 1, 1, 0, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0, 1.0027004039764054},
    {13, 0, 0, 1, 0, 0.34185651698814973, 0.4573223499816414, 0.050000000000000003, 0, 0.99758450248485853},
    {13, 0, 0, 1, 0, 0.5096085278044844, 1.131370473535944, 0.050000000000000003, 0, 1.2283268970288379},
    {13, 0, 0, 1, 0, 1.4320889980775615, 0.46918604154983751, 0.050000000000000003, 0, 1.0027004039764054},
};

inline const label nSelfTestSamples = 819;

} // End namespace transmissionMMS
} // End namespace Foam

#endif
