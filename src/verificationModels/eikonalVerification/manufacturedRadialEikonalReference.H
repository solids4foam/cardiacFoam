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
    Exact source-free eikonal solution around a held ball in a sealed sphere,
    in s = r/L with L = Ro - R, a = R/L, b = a + 1 and Pe = c0 L/sqrt(m):
    (w' + 2 w/s)/Pe = w - 1, w(b) = 0, where w = V dpsi/dr and V = c0 sqrt(m).
    The generator tutorials/manufacturedSolutions/generators/radial_eikonal_mms.py
    derives it and tests these samples.

Author
    Simao Nieto de Castro, UCD.
\*---------------------------------------------------------------------------*/

#ifndef manufacturedRadialEikonalReference_H
#define manufacturedRadialEikonalReference_H

#include "scalar.H"

#include <cmath>

namespace Foam
{
namespace radialEikonalMMS
{

struct SelfTestSample
{
    scalar Pe, a, s, slope, scaledTime;
};

constexpr int nSelfTestSamples = 7;

constexpr SelfTestSample selfTestSamples[nSelfTestSamples] =
{
    {2.0, 0.05, 0.07, 41.076130307536346, 1.1254141427520912},
    {5.0, 0.02, 0.03, 91.193019769818591, 1.3301088642241421},
    {21.0, 0.25, 0.6, 1.1713222637922323, 0.44395966643462248},
    {46.0, 0.25, 0.125, 1.4083175803402646, -0.15891755227387286},
    {46.0, 0.25, 1.0, 1.0444070520591153, 0.81310883375863241},
    {94.0, 0.1, 1.05, 1.01029166451527, 1.0019666300640726},
    {94.0, 0.5, 1.45, 1.0049087769924197, 0.9728434107394266},
};

inline scalar q(const scalar s, const scalar Pe)
{
    return s*s + 2*s/Pe + 2/(Pe*Pe);
}

//- w = V dpsi/dr, 1 for a plane front
inline scalar slope(const scalar s, const scalar Pe, const scalar a)
{
    const scalar b = a + 1;
    return (q(s, Pe) - q(b, Pe)*std::exp(-Pe*(b - s)))/(s*s);
}

//- Integral from s to b of exp(-Pe (b - t))/t^2 dt: Gauss-Legendre on panels
//  that follow the exponential layer at the wall and the growth of 1/t^2 at small t
inline scalar wallIntegral(const scalar s, const scalar Pe, const scalar a)
{
    static constexpr int nGauss = 8;
    static constexpr scalar x[nGauss] =
    {
        -0.96028985649753623, -0.79666647741362674,
        -0.52553240991632899, -0.18343464249564980,
         0.18343464249564980,  0.52553240991632899,
         0.79666647741362674,  0.96028985649753623
    };
    static constexpr scalar g[nGauss] =
    {
        0.10122853629037626, 0.22238103445337447,
        0.31370664587788729, 0.36268378337836198,
        0.36268378337836198, 0.31370664587788729,
        0.22238103445337447, 0.10122853629037626
    };

    const scalar b = a + 1;
    const scalar length = b - s;
    scalar tauLow = 0;
    scalar sum = 0;
    for (scalar scale = 2/Pe; tauLow < length; scale *= 2)
    {
        // Each panel is limited by the layer at the wall and by 1/t^2 at the other end
        const scalar step = scale < 0.5*(b - tauLow) ? scale : 0.5*(b - tauLow);
        const scalar tauHigh = tauLow + step < length ? tauLow + step : length;
        for (int i = 0; i < nGauss; ++i)
        {
            const scalar tau = 0.5*(tauHigh + tauLow + (tauHigh - tauLow)*x[i]);
            const scalar t = b - tau;
            sum += 0.5*(tauHigh - tauLow)*g[i]*std::exp(-Pe*tau)/(t*t);
        }
        tauLow = tauHigh;
    }
    return sum;
}

//- V (psi - psi(a))/L
inline scalar scaledTime(const scalar s, const scalar Pe, const scalar a)
{
    const scalar b = a + 1;
    const auto regular = [Pe](const scalar t)
    {
        return t + 2*std::log(t)/Pe - 2/(Pe*Pe*t);
    };
    return
        regular(s) - regular(a)
      - q(b, Pe)*(wallIntegral(a, Pe, a) - wallIntegral(s, Pe, a));
}

} // End namespace radialEikonalMMS
} // End namespace Foam

#endif

// ************************************************************************* //
