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
    manufacturedUnequalConormalBidomainReference

Description
    Source-free conormal bidomain reference solution for unequal
    intracellular and extracellular conductivity tensors:
    phiE = s Re(E) + c(t), phiI = s Re(I) + c(t), Vm = s Re(I - E), with
    E = exp(iq(y - re x)), I = exp(iq(y - ri x)) W(x), s = sqrt(1 + t),
    div(Gi grad I) = -div(Ge grad E) and n.Ge.grad E = n.Gi.grad I = 0 on
    x = 0 and x = 1.

Author
    Simao Nieto de Castro, UCD.
\*---------------------------------------------------------------------------*/

#ifndef manufacturedUnequalConormalBidomainReference_H
#define manufacturedUnequalConormalBidomainReference_H

#include "complex.H"
#include "tensor.H"
#include "scalarField.H"
#include "vectorField.H"

namespace Foam
{

inline void computeUnequalConormalG
(
    scalarField& G,
    const vectorField& centres
)
{
    G.setSize(centres.size());

    forAll(centres, cellI)
    {
        const vector& c = centres[cellI];
        G[cellI] = 1.0 + c.x()*Foam::sqr(c.y())*Foam::pow(c.z(), 3);
    }
}


inline void computeUnequalConormalU
(
    scalarField& u1,
    scalarField& u2,
    scalarField& u3,
    const scalarField& G,
    const scalarField& V,
    const scalar t
)
{
    u1.setSize(G.size());
    u2.setSize(G.size());
    u3.setSize(G.size());

    const scalar inv1t = 1.0/(1.0 + t);

    forAll(G, cellI)
    {
        u1[cellI] = (1.0 + t)*G[cellI] + V[cellI];
        u2[cellI] = inv1t*Foam::pow(G[cellI], -0.5);
        u3[cellI] = 0.0;
    }
}


inline void computeUnequalConormalPhiI
(
    scalarField& phiI,
    const scalarField& V,
    const scalarField& phiE
)
{
    phiI.setSize(V.size());

    forAll(V, cellI)
    {
        phiI[cellI] = V[cellI] + phiE[cellI];
    }
}


struct unequalConormalProfile
{
    scalar q;
    scalar ri;
    scalar re;
    scalar delta;
    scalar alphaE;
    scalar lambda;
    complex A;
    complex B;
    complex C;
};


inline complex unequalConormalPhase(const scalar theta)
{
    return complex(Foam::cos(theta), Foam::sin(theta));
}


inline unequalConormalProfile unequalConormalWallProfile
(
    const tensor& Gi,
    const tensor& Ge,
    const scalar q
)
{
    unequalConormalProfile p;
    p.q = q;
    p.ri = Gi.xy()/Gi.xx();
    p.re = Ge.xy()/Ge.xx();
    p.delta = p.re - p.ri;
    const scalar alphaI = Gi.yy() - Gi.xy()*p.ri;
    p.alphaE = Ge.yy() - Ge.xy()*p.re;
    p.lambda = q*Foam::sqrt(alphaI/Gi.xx());
    p.C = complex(-p.alphaE/(Gi.xx()*sqr(p.delta) + alphaI), 0);
    const complex iqdC = complex(0, q*p.delta)*p.C;
    p.B = iqdC/p.lambda;
    p.A =
        iqdC*(unequalConormalPhase(-q*p.delta) - Foam::cosh(p.lambda))
       /(p.lambda*Foam::sinh(p.lambda));
    return p;
}


inline complex unequalConormalE(const unequalConormalProfile& p, const vector& c)
{
    return unequalConormalPhase(p.q*(c.y() - p.re*c.x()));
}


inline complex unequalConormalI(const unequalConormalProfile& p, const vector& c)
{
    const scalar x = c.x();
    const complex W =
        p.A*Foam::cosh(p.lambda*x)
      + p.B*Foam::sinh(p.lambda*x)
      + p.C*unequalConormalPhase(-p.q*p.delta*x);
    return unequalConormalPhase(p.q*(c.y() - p.ri*c.x()))*W;
}


inline void computeUnequalConormalV
(
    scalarField& V,
    const unequalConormalProfile& p,
    const vectorField& centres,
    const scalar t
)
{
    V.setSize(centres.size());
    const scalar s = Foam::sqrt(1.0 + t);

    forAll(centres, cellI)
    {
        V[cellI] =
            s*(unequalConormalI(p, centres[cellI])
          - unequalConormalE(p, centres[cellI])).real();
    }
}


inline scalar manufacturedUnequalConormalBidomainShift
(
    const unequalConormalProfile& p,
    const scalar t,
    const vector& referencePoint,
    const scalar phiEReferenceValue
)
{
    return
        phiEReferenceValue
      - Foam::sqrt(1.0 + t)*unequalConormalE(p, referencePoint).real();
}


inline void computeUnequalConormalPhiE
(
    scalarField& phiE,
    const unequalConormalProfile& p,
    const vectorField& centres,
    const scalar t,
    const scalar shift
)
{
    phiE.setSize(centres.size());
    const scalar s = Foam::sqrt(1.0 + t);

    forAll(centres, cellI)
    {
        phiE[cellI] = s*unequalConormalE(p, centres[cellI]).real() + shift;
    }
}


inline void computeUnequalConormalPdeSource
(
    scalarField& source,
    const unequalConormalProfile& p,
    const vectorField& centres,
    const scalar t,
    const scalar beta
)
{
    source.setSize(centres.size());
    const scalar s = Foam::sqrt(1.0 + t);
    const scalar q2AlphaE = sqr(p.q)*p.alphaE;

    forAll(centres, cellI)
    {
        const scalar E = unequalConormalE(p, centres[cellI]).real();
        const scalar I = unequalConormalI(p, centres[cellI]).real();
        source[cellI] = s*(beta*(I - E) - q2AlphaE*E);
    }
}

} // End namespace Foam

#endif

// ************************************************************************* //
