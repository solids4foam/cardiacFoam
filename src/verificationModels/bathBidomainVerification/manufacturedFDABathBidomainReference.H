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
    manufacturedFDABathBidomainReference

Description
    Manufactured FDA bath bidomain reference solution.

Author
    Simao Nieto de Castro. All rights reserved.
\*---------------------------------------------------------------------------*/

#ifndef manufacturedFDABathBidomainReference_H
#define manufacturedFDABathBidomainReference_H

#include "mathematicalConstants.H"
#include "monodomainVerification/manufacturedFDAReference.H"
#include "scalarField.H"
#include "tensorField.H"
#include "volFields.H"

namespace Foam
{

inline scalar manufacturedFDABathSe(const volTensorField& Ge)
{
    scalar se = 0.0;

    if (!Ge.primitiveField().empty())
    {
        se = Ge.primitiveField()[0].xx();
    }

    reduce(se, maxOp<scalar>());

    if (se <= SMALL)
    {
        FatalErrorInFunction
            << "FDA bath-bidomain manufactured solution requires a positive "
            << "extracellular xx conductivity component s_e. Found " << se
            << exit(FatalError);
    }

    return se;
}


inline scalar manufacturedFDABathSigmaB(const scalar se)
{
    return 0.5*se;
}


inline scalar manufacturedFDABathC
(
    const scalar t,
    const scalar k,
    const scalar alpha,
    const scalar se
)
{
    return k*Foam::sqrt(1.0 + t) + alpha/manufacturedFDABathSigmaB(se);
}


inline void computeManufacturedFDABathF
(
    scalarField& F,
    const scalarField& x
)
{
    F.setSize(x.size());

    forAll(x, i)
    {
        F[i] = Foam::cos(constant::mathematical::pi*x[i]);
    }
}


inline void computeManufacturedFDABathV
(
    scalarField& V,
    const scalarField& x,
    const scalarField& y,
    const scalarField& z,
    const scalar t,
    const label dimension,
    const scalar alpha,
    const scalar se
)
{
    (void)y;
    (void)z;
    (void)dimension;

    scalarField F;
    computeManufacturedFDABathF(F, x);

    V.setSize(x.size());

    const scalar coeff = Foam::sqrt(1.0 + t);

    forAll(x, i)
    {
        V[i] = coeff*F[i] - alpha*x[i]/se;
    }
}


inline void computeManufacturedFDABathU
(
    scalarField& u1,
    scalarField& u2,
    scalarField& u3,
    const scalarField& x,
    const scalarField& y,
    const scalarField& z,
    const scalar t,
    const label dimension,
    const scalar alpha,
    const scalar se
)
{
    scalarField F, G;
    computeManufacturedFDABathF(F, x);
    computeG(G, x, y, z, dimension);

    u1.setSize(x.size());
    u2.setSize(x.size());
    u3.setSize(x.size());

    const scalar sqrt1t = Foam::sqrt(1.0 + t);
    const scalar inv1t = 1.0/(1.0 + t);

    forAll(x, i)
    {
        u1[i] = (1.0 + t)*G[i] + sqrt1t*F[i];
        u2[i] = inv1t*Foam::pow(G[i], -0.5);
        u3[i] = -alpha*x[i]/se;
    }
}


inline scalar manufacturedFDABathPhiEAt
(
    const scalar x,
    const scalar t,
    const scalar k,
    const scalar alpha,
    const scalar se
)
{
    const scalar sqrt1t = Foam::sqrt(1.0 + t);
    const scalar sigmaB = manufacturedFDABathSigmaB(se);
    const scalar C = manufacturedFDABathC(t, k, alpha, se);

    if (x <= 0.0)
    {
        return -k*sqrt1t + alpha*x/sigmaB + C;
    }
    else if (x <= 1.0)
    {
        return
            -k*sqrt1t*Foam::cos(constant::mathematical::pi*x)
          + alpha*x/se
          + C;
    }

    return
        -k*sqrt1t*Foam::cos(constant::mathematical::pi)
      + alpha/se
      + alpha*(x - 1.0)/sigmaB
      + C;
}


inline void computeManufacturedFDABathPhiE
(
    scalarField& phiE,
    const scalarField& x,
    const scalar t,
    const scalar k,
    const scalar alpha,
    const scalar se
)
{
    phiE.setSize(x.size());

    forAll(x, i)
    {
        phiE[i] = manufacturedFDABathPhiEAt(x[i], t, k, alpha, se);
    }
}


inline void computeManufacturedFDABathPhiI
(
    scalarField& phiI,
    const scalarField& x,
    const scalarField& y,
    const scalarField& z,
    const scalar t,
    const label dimension,
    const scalar k,
    const scalar alpha,
    const scalar se
)
{
    scalarField V;
    scalarField phiE;

    computeManufacturedFDABathV(V, x, y, z, t, dimension, alpha, se);
    computeManufacturedFDABathPhiE(phiE, x, t, k, alpha, se);

    phiI.setSize(x.size());

    forAll(x, i)
    {
        phiI[i] = V[i] + phiE[i];
    }
}

} // End namespace Foam

#endif

// ************************************************************************* //
