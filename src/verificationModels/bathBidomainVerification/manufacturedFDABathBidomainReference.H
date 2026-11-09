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
    Manufactured FDA bath bidomain reference solution: a heart slab
    0 <= x <= 1 between bath slabs, grounded on x = -1, with the current
    alpha into x = 2.

    The formulas between "// BEGIN GENERATED <name>" and
    "// END GENERATED <name>" are written from the sympy definitions in
    tutorials/manufacturedSolutions/generators/fda_bath_bidomain_mms.py.
    Regenerate them with --write and compare with --check; the rest of the
    file is hand-written.

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
    // BEGIN GENERATED manufacturedFDABathSigmaB
    return (1.0/2.0)*se;
    // END GENERATED manufacturedFDABathSigmaB
}


inline scalar manufacturedFDABathC
(
    const scalar t,
    const scalar k,
    const scalar alpha,
    const scalar se
)
{
    const scalar sigmaB = manufacturedFDABathSigmaB(se);

    // BEGIN GENERATED manufacturedFDABathC
    return alpha/sigmaB + k*sqrt(t + 1);
    // END GENERATED manufacturedFDABathC
}


inline void computeManufacturedFDABathF
(
    scalarField& F,
    const scalarField& x
)
{
    const scalar pi = constant::mathematical::pi;

    F.setSize(x.size());

    forAll(x, i)
    {
        // BEGIN GENERATED computeManufacturedFDABathF
        F[i] = cos(pi*x[i]);
        // END GENERATED computeManufacturedFDABathF
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

    forAll(x, i)
    {
        // BEGIN GENERATED computeManufacturedFDABathV
        V[i] = F[i]*sqrt(t + 1) - alpha*x[i]/se;
        // END GENERATED computeManufacturedFDABathV
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

    forAll(x, i)
    {
        // BEGIN GENERATED computeManufacturedFDABathU
        u1[i] = F[i]*sqrt(t + 1) + G[i]*(t + 1);
        u2[i] = 1/(sqrt(G[i])*(t + 1));
        u3[i] = -alpha*x[i]/se;
        // END GENERATED computeManufacturedFDABathU
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
    const scalar pi = constant::mathematical::pi;
    const scalar sigmaB = manufacturedFDABathSigmaB(se);
    const scalar C = manufacturedFDABathC(t, k, alpha, se);

    if (x <= 0)
    {
        // BEGIN GENERATED bathPhiE_leftBath
        return C + alpha*x/sigmaB - k*sqrt(t + 1);
        // END GENERATED bathPhiE_leftBath
    }
    if (x <= 1)
    {
        // BEGIN GENERATED bathPhiE_heart
        return C + alpha*x/se - k*sqrt(t + 1)*cos(pi*x);
        // END GENERATED bathPhiE_heart
    }

    // BEGIN GENERATED bathPhiE_rightBath
    return C + alpha*(x - 1)/sigmaB + alpha/se + k*sqrt(t + 1);
    // END GENERATED bathPhiE_rightBath
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
        // BEGIN GENERATED computeManufacturedFDABathPhiI
        phiI[i] = V[i] + phiE[i];
        // END GENERATED computeManufacturedFDABathPhiI
    }
}


//- V_m at point p
inline scalar manufacturedFDABathVmAt
(
    const point& p,
    const scalar t,
    const scalar alpha,
    const scalar se
)
{
    const scalar pi = constant::mathematical::pi;
    const scalar x = p.x();

    // BEGIN GENERATED manufacturedFDABathVmAt
    return -alpha*x/se + sqrt(t + 1)*cos(pi*x);
    // END GENERATED manufacturedFDABathVmAt
}


//- phi_i = V_m + phi_e at point p
inline scalar manufacturedFDABathPhiIAt
(
    const point& p,
    const scalar t,
    const scalar k,
    const scalar alpha,
    const scalar se
)
{
    const scalar V = manufacturedFDABathVmAt(p, t, alpha, se);
    const scalar phiE = manufacturedFDABathPhiEAt(p.x(), t, k, alpha, se);

    // BEGIN GENERATED manufacturedFDABathPhiIAt
    return V + phiE;
    // END GENERATED manufacturedFDABathPhiIAt
}


//- grad V_m at point p
inline vector manufacturedFDABathGradVm
(
    const point& p,
    const scalar t,
    const scalar alpha,
    const scalar se
)
{
    const scalar pi = constant::mathematical::pi;
    const scalar x = p.x();

    // BEGIN GENERATED manufacturedFDABathGradVm
    return vector(-alpha/se - pi*sqrt(t + 1)*sin(pi*x), 0, 0);
    // END GENERATED manufacturedFDABathGradVm
}


//- grad phi_e of the heart formula at point p
inline vector manufacturedFDABathGradPhiEHeart
(
    const point& p,
    const scalar t,
    const scalar k,
    const scalar alpha,
    const scalar se
)
{
    const scalar pi = constant::mathematical::pi;
    const scalar x = p.x();

    // BEGIN GENERATED manufacturedFDABathGradPhiEHeart
    return vector(alpha/se + pi*k*sqrt(t + 1)*sin(pi*x), 0, 0);
    // END GENERATED manufacturedFDABathGradPhiEHeart
}


//- grad phi_i of the heart formulas at point p
inline vector manufacturedFDABathGradPhiIHeart
(
    const point& p,
    const scalar t,
    const scalar k,
    const scalar alpha,
    const scalar se
)
{
    const scalar pi = constant::mathematical::pi;
    const scalar x = p.x();

    // BEGIN GENERATED manufacturedFDABathGradPhiIHeart
    return vector(pi*k*sqrt(t + 1)*sin(pi*x) - pi*sqrt(t + 1)*sin(pi*x), 0, 0);
    // END GENERATED manufacturedFDABathGradPhiIHeart
}


//- grad phi_e of the bath formula (uniform)
inline vector manufacturedFDABathGradPhiEBath
(
    const scalar alpha,
    const scalar se
)
{
    const scalar sigmaB = manufacturedFDABathSigmaB(se);

    // BEGIN GENERATED manufacturedFDABathGradPhiEBath
    return vector(alpha/sigmaB, 0, 0);
    // END GENERATED manufacturedFDABathGradPhiEBath
}

} // End namespace Foam

#endif

// ************************************************************************* //
