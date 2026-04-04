/*---------------------------------------------------------------------------*\
License
    This file is part of cardiacFoam.

    cardiacFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    cardiacFoam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
    Public License for more details.

    You should have received a copy of the GNU General Public License
    along with cardiacFoam.  If not, see <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#ifndef manufacturedFDAReference_H
#define manufacturedFDAReference_H

#include "List.H"
#include "Switch.H"
#include "mathematicalConstants.H"
#include "scalarField.H"
#include "tensor.H"
#include "vector.H"

namespace Foam
{

inline void computeF
(
    scalarField& F,
    const scalarField& x,
    const scalarField& y,
    const scalarField& z,
    const int dimension
)
{
    F.setSize(x.size());

    forAll(x, i)
    {
        if (dimension == 1)
        {
            F[i] = Foam::cos(constant::mathematical::pi*x[i]);
        }
        else if (dimension == 2)
        {
            F[i] = Foam::cos(constant::mathematical::pi*x[i])
                 * Foam::cos(2.0*constant::mathematical::pi*y[i]);
        }
        else
        {
            F[i] = Foam::cos(constant::mathematical::pi*x[i])
                 * Foam::cos(2.0*constant::mathematical::pi*y[i])
                 * Foam::cos(3.0*constant::mathematical::pi*z[i]);
        }
    }
}

inline void computeG
(
    scalarField& G,
    const scalarField& x,
    const scalarField& y,
    const scalarField& z,
    const int dimension
)
{
    G.setSize(x.size());

    forAll(x, i)
    {
        if (dimension == 1)
        {
            G[i] = 1.0 + x[i];
        }
        else if (dimension == 2)
        {
            G[i] = 1.0 + x[i]*Foam::pow(y[i], 2);
        }
        else
        {
            G[i] = 1.0 + x[i]*Foam::pow(y[i], 2)*Foam::pow(z[i], 3);
        }
    }
}

inline void computeManufacturedV
(
    scalarField& V,
    const scalarField& x,
    const scalarField& y,
    const scalarField& z,
    const scalar t,
    const int dimension
)
{
    scalarField Fx;
    computeF(Fx, x, y, z, dimension);

    V.setSize(x.size());

    const scalar coeff = Foam::sqrt(1.0 + t);

    forAll(x, i)
    {
        V[i] = coeff*Fx[i];
    }
}

inline void computeManufacturedU
(
    scalarField& u1,
    scalarField& u2,
    scalarField& u3,
    const scalarField& x,
    const scalarField& y,
    const scalarField& z,
    const scalar t,
    const int dimension
)
{
    scalarField Fx, Gx;
    computeF(Fx, x, y, z, dimension);
    computeG(Gx, x, y, z, dimension);

    u1.setSize(x.size());
    u2.setSize(x.size());
    u3.setSize(x.size());

    const scalar sqrt1t = Foam::sqrt(1.0 + t);
    const scalar inv1t = 1.0/(1.0 + t);

    forAll(x, i)
    {
        u1[i] = (1.0 + t)*Gx[i] + sqrt1t*Fx[i];
        u2[i] = inv1t*Foam::pow(Gx[i], -0.5);
        u3[i] = 0.0;
    }
}

inline void computeIion
(
    scalarField& Iion,
    const scalarField& u1,
    const scalarField& u2,
    const scalarField& u3,
    const scalarField& V,
    const scalar Cm,
    const scalar beta,
    const scalar chi
)
{
    Iion.setSize(V.size());

    forAll(V, i)
    {
        const scalar term1 =
            -(Cm/2.0)*(u1[i] + u3[i] - V[i])*Foam::pow(u2[i], 2)
           *(V[i] - u3[i]);
        const scalar term2 = (beta/chi)*(V[i] - u3[i]);
        Iion[i] = term1 + term2;
    }
}

inline bool pseudoECGManufacturedSupportsQuadratureOrder
(
    const label quadratureOrder
)
{
    return quadratureOrder > 0;
}

inline void pseudoECGManufacturedLegendrePolynomial
(
    const label order,
    const scalar x,
    scalar& value,
    scalar& derivative
)
{
    scalar pNm2 = 1.0;
    scalar pNm1 = x;

    if (order == 0)
    {
        value = pNm2;
        derivative = 0.0;
        return;
    }

    if (order == 1)
    {
        value = pNm1;
        derivative = 1.0;
        return;
    }

    scalar pN = pNm1;

    for (label n = 2; n <= order; ++n)
    {
        pN =
        (
            ((2.0*scalar(n) - 1.0)*x*pNm1 - (scalar(n) - 1.0)*pNm2)
          / scalar(n)
        );
        pNm2 = pNm1;
        pNm1 = pN;
    }

    value = pN;
    derivative = scalar(order)*(x*pN - pNm2)/(x*x - 1.0);
}

inline void pseudoECGManufacturedQuadratureRule
(
    const label quadratureOrder,
    List<scalar>& nodes,
    List<scalar>& weights
)
{
    if (!pseudoECGManufacturedSupportsQuadratureOrder(quadratureOrder))
    {
        FatalErrorInFunction
            << "Manufactured pseudo-ECG quadrature order must be positive. "
            << "Requested order " << quadratureOrder << "."
            << exit(FatalError);
    }

    nodes.setSize(quadratureOrder);
    weights.setSize(quadratureOrder);

    const label nRoots = (quadratureOrder + 1)/2;
    const scalar tolerance = 1e-14;
    const label maxIterations = 64;

    for (label i = 0; i < nRoots; ++i)
    {
        scalar root = Foam::cos
        (
            constant::mathematical::pi
          * (scalar(i) + 0.75)
          / (scalar(quadratureOrder) + 0.5)
        );
        scalar value = 0.0;
        scalar derivative = 0.0;

        label iteration = 0;
        for (; iteration < maxIterations; ++iteration)
        {
            pseudoECGManufacturedLegendrePolynomial
            (
                quadratureOrder,
                root,
                value,
                derivative
            );

            const scalar delta = value/derivative;
            root -= delta;

            if (Foam::mag(delta) <= tolerance)
            {
                break;
            }
        }

        if (iteration == maxIterations)
        {
            FatalErrorInFunction
                << "Failed to converge Gauss-Legendre quadrature node for "
                << "order " << quadratureOrder << "."
                << exit(FatalError);
        }

        pseudoECGManufacturedLegendrePolynomial
        (
            quadratureOrder,
            root,
            value,
            derivative
        );

        const scalar weight =
            2.0/((1.0 - root*root)*derivative*derivative);
        const label mirrorI = quadratureOrder - 1 - i;

        nodes[i] = -root;
        weights[i] = weight;
        nodes[mirrorI] = root;
        weights[mirrorI] = weight;
    }
}

inline scalar pseudoECGManufacturedTimeFactor(const scalar t)
{
    return Foam::sqrt(1.0 + t);
}

inline vector pseudoECGManufacturedGradF
(
    const scalar x,
    const scalar y,
    const scalar z,
    const label dimension
)
{
    const scalar pi = constant::mathematical::pi;

    if (dimension == 1)
    {
        return vector(-pi*Foam::sin(pi*x), 0.0, 0.0);
    }

    if (dimension == 2)
    {
        return vector
        (
           -pi*Foam::sin(pi*x)*Foam::cos(2.0*pi*y),
           -2.0*pi*Foam::cos(pi*x)*Foam::sin(2.0*pi*y),
            0.0
        );
    }

    return vector
    (
       -pi*Foam::sin(pi*x)*Foam::cos(2.0*pi*y)*Foam::cos(3.0*pi*z),
       -2.0*pi*Foam::cos(pi*x)*Foam::sin(2.0*pi*y)*Foam::cos(3.0*pi*z),
       -3.0*pi*Foam::cos(pi*x)*Foam::cos(2.0*pi*y)*Foam::sin(3.0*pi*z)
    );
}

inline scalar pseudoECGManufacturedIntegrand
(
    const scalar x,
    const scalar y,
    const scalar z,
    const vector& electrodePosition,
    const tensor& conductivity,
    const label dimension
)
{
    const vector sigmaGradF =
        conductivity & pseudoECGManufacturedGradF(x, y, z, dimension);

    if (dimension == 1)
    {
        const scalar dx = x - electrodePosition.x();
        const scalar absDx = Foam::mag(dx);

        if (absDx <= VSMALL)
        {
            return 0.0;
        }

        return sigmaGradF.x()*dx/Foam::pow(absDx, 3);
    }

    if (dimension == 2)
    {
        const scalar dx = x - electrodePosition.x();
        const scalar dy = y - electrodePosition.y();
        const scalar rSq = dx*dx + dy*dy;

        if (rSq <= VSMALL)
        {
            return 0.0;
        }

        return (sigmaGradF.x()*dx + sigmaGradF.y()*dy)/Foam::pow(rSq, 1.5);
    }

    const scalar dx = x - electrodePosition.x();
    const scalar dy = y - electrodePosition.y();
    const scalar dz = z - electrodePosition.z();
    const scalar rSq = dx*dx + dy*dy + dz*dz;

    if (rSq <= VSMALL)
    {
        return 0.0;
    }

    return
        (sigmaGradF.x()*dx + sigmaGradF.y()*dy + sigmaGradF.z()*dz)
       /Foam::pow(rSq, 1.5);
}

inline scalar computeManufacturedPseudoECGSpatialReference
(
    const vector& electrodePosition,
    const tensor& conductivity,
    const label dimension,
    const UList<scalar>& nodes,
    const UList<scalar>& weights
)
{
    const scalar jacobian = 0.5;
    scalar value = 0.0;

    if (dimension == 1)
    {
        forAll(nodes, i)
        {
            const scalar x = jacobian*(nodes[i] + 1.0);
            value += weights[i]*pseudoECGManufacturedIntegrand
            (
                x,
                0.0,
                0.0,
                electrodePosition,
                conductivity,
                dimension
            );
        }

        return jacobian*value;
    }

    if (dimension == 2)
    {
        forAll(nodes, i)
        {
            const scalar x = jacobian*(nodes[i] + 1.0);
            forAll(nodes, j)
            {
                const scalar y = jacobian*(nodes[j] + 1.0);
                value += weights[i]*weights[j]
                    * pseudoECGManufacturedIntegrand
                    (
                        x,
                        y,
                        0.0,
                        electrodePosition,
                        conductivity,
                        dimension
                    );
            }
        }

        return jacobian*jacobian*value;
    }

    forAll(nodes, i)
    {
        const scalar x = jacobian*(nodes[i] + 1.0);
        forAll(nodes, j)
        {
            const scalar y = jacobian*(nodes[j] + 1.0);
            forAll(nodes, k)
            {
                const scalar z = jacobian*(nodes[k] + 1.0);
                value += weights[i]*weights[j]*weights[k]
                    * pseudoECGManufacturedIntegrand
                    (
                        x,
                        y,
                        z,
                        electrodePosition,
                        conductivity,
                        dimension
                    );
            }
        }
    }

    return jacobian*jacobian*jacobian*value;
}

inline scalar computeManufacturedPseudoECGSpatialReference
(
    const vector& electrodePosition,
    const tensor& conductivity,
    const label dimension,
    const label quadratureOrder
)
{
    List<scalar> nodes, weights;
    pseudoECGManufacturedQuadratureRule(quadratureOrder, nodes, weights);

    return computeManufacturedPseudoECGSpatialReference
    (
        electrodePosition,
        conductivity,
        dimension,
        nodes,
        weights
    );
}

inline scalar computeManufacturedPseudoECGReference
(
    const scalar t,
    const vector& electrodePosition,
    const tensor& conductivity,
    const label dimension,
    const label quadratureOrder
)
{
    return
        pseudoECGManufacturedTimeFactor(t)
       * computeManufacturedPseudoECGSpatialReference
         (
             electrodePosition,
             conductivity,
             dimension,
             quadratureOrder
         );
}

} // End namespace Foam

#endif
