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

Class
    conormalManufacturedProfile

Description
    Spatial profile F = Re[v(x) exp(i(omegaY y + omegaZ z))] with
    v(x) = v0 + v1 x + v2 x^2, periodic in y and z.  conormalWallProfile
    returns the v with v0 = 1 that satisfies n.G.grad(F) = 0 on x = 0 and
    x = 1 for a constant tensor G.

Author
    Simao Nieto de Castro, UCD.
\*---------------------------------------------------------------------------*/

#ifndef conormalManufacturedProfile_H
#define conormalManufacturedProfile_H

#include "complex.H"
#include "tensor.H"

namespace Foam
{

struct conormalProfile
{
    complex v0;
    complex v1;
    complex v2;
    scalar omegaY;
    scalar omegaZ;
};


inline complex conormalWallFlux(const tensor& G, const conormalProfile& p)
{
    return complex(0, G.xy()*p.omegaY + G.xz()*p.omegaZ);
}


inline conormalProfile conormalWallProfile
(
    const tensor& G,
    const scalar omegaY,
    const scalar omegaZ
)
{
    conormalProfile p{complex(1, 0), complex(0, 0), complex(0, 0), omegaY, omegaZ};
    const complex ig = conormalWallFlux(G, p);
    p.v1 = -ig/G.xx();
    p.v2 = -(G.xx()*p.v1 + ig*(p.v0 + p.v1))/(complex(2*G.xx(), 0) + ig);
    return p;
}


inline complex conormalPhase(const conormalProfile& p, const vector& c)
{
    const scalar theta = p.omegaY*c.y() + p.omegaZ*c.z();
    return complex(Foam::cos(theta), Foam::sin(theta));
}


inline complex conormalV(const conormalProfile& p, const scalar x)
{
    return p.v0 + p.v1*x + p.v2*(x*x);
}


inline complex conormalDV(const conormalProfile& p, const scalar x)
{
    return p.v1 + p.v2*(2*x);
}


inline scalar conormalF(const conormalProfile& p, const vector& c)
{
    return (conormalV(p, c.x())*conormalPhase(p, c)).real();
}


inline vector conormalGradF(const conormalProfile& p, const vector& c)
{
    const complex e = conormalPhase(p, c);
    const complex v = conormalV(p, c.x());
    return vector
    (
        (conormalDV(p, c.x())*e).real(),
        (complex(0, p.omegaY)*v*e).real(),
        (complex(0, p.omegaZ)*v*e).real()
    );
}


inline scalar conormalDivGGradF
(
    const conormalProfile& p,
    const tensor& G,
    const vector& c
)
{
    const scalar oy = p.omegaY;
    const scalar oz = p.omegaZ;
    const complex mixed
    (
        0,
        (G.xy() + G.yx())*oy + (G.xz() + G.zx())*oz
    );
    const scalar tangential =
        G.yy()*oy*oy + (G.yz() + G.zy())*oy*oz + G.zz()*oz*oz;
    return
    (
        (
            p.v2*(2*G.xx())
          + mixed*conormalDV(p, c.x())
          - conormalV(p, c.x())*tangential
        )
       *conormalPhase(p, c)
    ).real();
}

}

#endif

// ************************************************************************* //
