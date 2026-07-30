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
    manufacturedAnisotropicMonodomainReference

Description
    Boundary-compatible manufactured reference solution for the monodomain
    equation with a constant, fully populated anisotropic conductivity tensor.

    The spatial factor

        F = sin^2(pi*x) sin^2(pi*y) sin^2(pi*z)

    has zero gradient on every face of the unit cube.  It therefore satisfies
    the homogeneous-flux boundary condition for any constant tensor.  The
    continuous PDE source is evaluated from the analytic Hessian; no discrete
    finite-volume operator is used to manufacture the source.

Author
    Simao Nieto de Castro. All rights reserved.
\*---------------------------------------------------------------------------*/

#ifndef manufacturedAnisotropicMonodomainReference_H
#define manufacturedAnisotropicMonodomainReference_H

#include "mathematicalConstants.H"
#include "scalarField.H"
#include "tensorField.H"
#include "vectorField.H"

namespace Foam
{

inline void computeAnisotropicManufacturedF
(
    scalarField& F,
    const vectorField& centres
)
{
    F.setSize(centres.size());
    const scalar pi = constant::mathematical::pi;

    forAll(centres, cellI)
    {
        const vector& c = centres[cellI];
        F[cellI] =
            Foam::sqr(Foam::sin(pi*c.x()))
           *Foam::sqr(Foam::sin(pi*c.y()))
           *Foam::sqr(Foam::sin(pi*c.z()));
    }
}


inline void computeAnisotropicManufacturedV
(
    scalarField& V,
    const vectorField& centres,
    const scalar t
)
{
    scalarField F;
    computeAnisotropicManufacturedF(F, centres);

    V.setSize(centres.size());
    const scalar sqrt1t = Foam::sqrt(1.0 + t);

    forAll(V, cellI)
    {
        V[cellI] = sqrt1t*F[cellI];
    }
}


inline void computeAnisotropicManufacturedU
(
    scalarField& u1,
    scalarField& u2,
    scalarField& u3,
    const vectorField& centres,
    const scalar t
)
{
    scalarField F;
    computeAnisotropicManufacturedF(F, centres);

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


inline void computeAnisotropicManufacturedPdeSource
(
    scalarField& source,
    const vectorField& centres,
    const tensorField& conductivity,
    const scalar t,
    const scalar beta
)
{
    if (conductivity.size() != centres.size())
    {
        FatalErrorInFunction
            << "Cannot evaluate the anisotropic manufactured source: "
            << "conductivity size " << conductivity.size()
            << " differs from the number of cell centres " << centres.size()
            << exit(FatalError);
    }

    source.setSize(centres.size());

    const scalar pi = constant::mathematical::pi;
    const scalar piSqr = Foam::sqr(pi);
    const scalar sqrt1t = Foam::sqrt(1.0 + t);

    forAll(centres, cellI)
    {
        const vector& c = centres[cellI];
        const tensor& K = conductivity[cellI];

        const scalar sx2 = Foam::sqr(Foam::sin(pi*c.x()));
        const scalar sy2 = Foam::sqr(Foam::sin(pi*c.y()));
        const scalar sz2 = Foam::sqr(Foam::sin(pi*c.z()));

        const scalar dx = pi*Foam::sin(2.0*pi*c.x());
        const scalar dy = pi*Foam::sin(2.0*pi*c.y());
        const scalar dz = pi*Foam::sin(2.0*pi*c.z());

        const scalar dxx = 2.0*piSqr*Foam::cos(2.0*pi*c.x());
        const scalar dyy = 2.0*piSqr*Foam::cos(2.0*pi*c.y());
        const scalar dzz = 2.0*piSqr*Foam::cos(2.0*pi*c.z());

        const scalar F = sx2*sy2*sz2;
        const scalar divKGradF =
            K.xx()*dxx*sy2*sz2
          + K.yy()*sx2*dyy*sz2
          + K.zz()*sx2*sy2*dzz
          + (K.xy() + K.yx())*dx*dy*sz2
          + (K.xz() + K.zx())*dx*sy2*dz
          + (K.yz() + K.zy())*sx2*dy*dz;

        // The solver uses +source on the PDE right-hand side.  The FDA
        // ionic model contributes -Cm*dV/dt + beta*V/chi, hence this source
        // must be beta*V - div(K grad(V)).
        source[cellI] = sqrt1t*(beta*F - divKGradF);
    }
}

} // End namespace Foam

#endif

// ************************************************************************* //
