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
    manufacturedConormalMonodomainReference

Description
    Boundary-compatible manufactured reference solution for the monodomain
    equation with a constant, fully populated anisotropic conductivity
    tensor, periodic in y and z.

    The spatial factor F = conormalF(p, c) satisfies n.G.grad(F) = 0 on
    x = 0 and x = 1 for the case conductivity tensor G, while n.grad(F) and
    the tangential gradient of F are both non-zero there. The continuous PDE
    source is evaluated from the analytic divergence conormalDivGGradF; no
    discrete finite-volume operator is used to manufacture the source.

Author
    Simao Nieto de Castro, UCD.
\*---------------------------------------------------------------------------*/

#ifndef manufacturedConormalMonodomainReference_H
#define manufacturedConormalMonodomainReference_H

#include "conormalManufacturedProfile.H"
#include "scalarField.H"
#include "tensorField.H"
#include "vectorField.H"

namespace Foam
{

inline void computeConormalManufacturedF
(
    scalarField& F,
    const vectorField& centres,
    const conormalProfile& p
)
{
    F.setSize(centres.size());

    forAll(centres, cellI)
    {
        F[cellI] = conormalF(p, centres[cellI]);
    }
}


inline void computeConormalManufacturedV
(
    scalarField& V,
    const vectorField& centres,
    const scalar t,
    const conormalProfile& p
)
{
    scalarField F;
    computeConormalManufacturedF(F, centres, p);

    V.setSize(centres.size());
    const scalar sqrt1t = Foam::sqrt(1.0 + t);

    forAll(V, cellI)
    {
        V[cellI] = sqrt1t*F[cellI];
    }
}


inline void computeConormalManufacturedU
(
    scalarField& u1,
    scalarField& u2,
    scalarField& u3,
    const vectorField& centres,
    const scalar t,
    const conormalProfile& p
)
{
    scalarField F;
    computeConormalManufacturedF(F, centres, p);

    u1.setSize(centres.size());
    u2.setSize(centres.size());
    u3.setSize(centres.size());

    const scalar sqrt1t = Foam::sqrt(1.0 + t);
    const scalar inv1t = 1.0/(1.0 + t);

    forAll(centres, cellI)
    {
        const vector& c = centres[cellI];
        const scalar G =
            1.0 + c.x()*Foam::sqr(c.y())*Foam::pow(c.z(), 3);

        u1[cellI] = (1.0 + t)*G + sqrt1t*F[cellI];
        u2[cellI] = inv1t*Foam::pow(G, -0.5);
        u3[cellI] = 0.0;
    }
}


inline void computeConormalManufacturedPdeSource
(
    scalarField& source,
    const vectorField& centres,
    const tensorField& conductivity,
    const scalar t,
    const scalar beta,
    const conormalProfile& p
)
{
    if (conductivity.size() != centres.size())
    {
        FatalErrorInFunction
            << "Cannot evaluate the conormal manufactured source: "
            << "conductivity size " << conductivity.size()
            << " differs from the number of cell centres " << centres.size()
            << exit(FatalError);
    }

    source.setSize(centres.size());

    const scalar sqrt1t = Foam::sqrt(1.0 + t);

    forAll(centres, cellI)
    {
        const vector& c = centres[cellI];
        const tensor& K = conductivity[cellI];

        const scalar F = conormalF(p, c);
        const scalar divKGradF = conormalDivGGradF(p, K, c);

        // The solver uses +source on the PDE right-hand side.  The FDA
        // ionic model contributes -Cm*dV/dt + beta*V/chi, hence this source
        // must be beta*V - div(K grad(V)).
        source[cellI] = sqrt1t*(beta*F - divKGradF);
    }
}

} // End namespace Foam

#endif

// ************************************************************************* //
