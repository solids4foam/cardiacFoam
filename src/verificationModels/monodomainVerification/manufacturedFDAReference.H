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
    manufacturedFDAReference

Description
    Manufactured FDA reference solution.

    The formulas between "// BEGIN GENERATED <name>" and
    "// END GENERATED <name>" are written from the sympy definitions in
    tutorials/manufacturedSolutions/generators: fda_monodomain_mms.py (F, G,
    V, u1-u3, Iion) and pseudo_ecg_mms.py (time factor, grad F, pseudo-ECG
    kernel).  Regenerate them with --write and compare with --check; the
    rest of the file is hand-written.

Author
    Simao Nieto de Castro. All rights reserved.
\*---------------------------------------------------------------------------*/

#ifndef manufacturedFDAReference_H
#define manufacturedFDAReference_H

#include "List.H"
#include "Switch.H"
#include "mathematicalConstants.H"
#include "scalarField.H"
#include "tensor.H"
#include "point.H"
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
    const scalar pi = constant::mathematical::pi;

    forAll(x, i)
    {
        if (dimension == 1)
        {
            // BEGIN GENERATED computeF_d1
            F[i] = cos(pi*x[i]);
            // END GENERATED computeF_d1
        }
        else if (dimension == 2)
        {
            // BEGIN GENERATED computeF_d2
            F[i] = cos(pi*x[i])*cos(2*pi*y[i]);
            // END GENERATED computeF_d2
        }
        else
        {
            // BEGIN GENERATED computeF_d3
            F[i] = cos(pi*x[i])*cos(2*pi*y[i])*cos(3*pi*z[i]);
            // END GENERATED computeF_d3
        }
    }
}

//- Spatial factor F at one point.
inline scalar manufacturedFDAFAt(const point& p, const label dimension)
{
    const scalar pi = constant::mathematical::pi;

    if (dimension <= 1)
    {
        // BEGIN GENERATED manufacturedFDAFAt_d1
        return cos(pi*p.x());
        // END GENERATED manufacturedFDAFAt_d1
    }

    if (dimension == 2)
    {
        // BEGIN GENERATED manufacturedFDAFAt_d2
        return cos(pi*p.x())*cos(2*pi*p.y());
        // END GENERATED manufacturedFDAFAt_d2
    }

    // BEGIN GENERATED manufacturedFDAFAt_d3
    return cos(pi*p.x())*cos(2*pi*p.y())*cos(3*pi*p.z());
    // END GENERATED manufacturedFDAFAt_d3
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
            // BEGIN GENERATED computeG_d1
            G[i] = x[i] + 1;
            // END GENERATED computeG_d1
        }
        else if (dimension == 2)
        {
            // BEGIN GENERATED computeG_d2
            G[i] = x[i]*(y[i]*y[i]) + 1;
            // END GENERATED computeG_d2
        }
        else
        {
            // BEGIN GENERATED computeG_d3
            G[i] = x[i]*(y[i]*y[i])*(z[i]*z[i]*z[i]) + 1;
            // END GENERATED computeG_d3
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

    forAll(x, i)
    {
        // BEGIN GENERATED computeManufacturedV
        V[i] = Fx[i]*sqrt(t + 1);
        // END GENERATED computeManufacturedV
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

    forAll(x, i)
    {
        // BEGIN GENERATED computeManufacturedU
        u1[i] = Fx[i]*sqrt(t + 1) + Gx[i]*(t + 1);
        u2[i] = 1/(sqrt(Gx[i])*(t + 1));
        u3[i] = 0;
        // END GENERATED computeManufacturedU
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
        // BEGIN GENERATED computeIion
        Iion[i] =
            -1.0/2.0*Cm*(u2[i]*u2[i])*(V[i] - u3[i])*(-V[i] + u1[i] + u3[i])
          + beta*(V[i] - u3[i])/chi;
        // END GENERATED computeIion
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
    // BEGIN GENERATED pseudoECGManufacturedTimeFactor
    return sqrt(t + 1);
    // END GENERATED pseudoECGManufacturedTimeFactor
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
        // BEGIN GENERATED pseudoECGManufacturedGradF_d1
        return vector(-pi*sin(pi*x), 0, 0);
        // END GENERATED pseudoECGManufacturedGradF_d1
    }

    if (dimension == 2)
    {
        // BEGIN GENERATED pseudoECGManufacturedGradF_d2
        return vector
        (
            -pi*sin(pi*x)*cos(2*pi*y),
            -2*pi*sin(2*pi*y)*cos(pi*x),
            0
        );
        // END GENERATED pseudoECGManufacturedGradF_d2
    }

    // BEGIN GENERATED pseudoECGManufacturedGradF_d3
    return vector
    (
        -pi*sin(pi*x)*cos(2*pi*y)*cos(3*pi*z),
        -2*pi*sin(2*pi*y)*cos(pi*x)*cos(3*pi*z),
        -3*pi*sin(3*pi*z)*cos(pi*x)*cos(2*pi*y)
    );
    // END GENERATED pseudoECGManufacturedGradF_d3
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

        // BEGIN GENERATED pseudoECGManufacturedIntegrand_d1
        return sigmaGradF.x()/(absDx*dx);
        // END GENERATED pseudoECGManufacturedIntegrand_d1
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

        // BEGIN GENERATED pseudoECGManufacturedIntegrand_d2
        return (dx*sigmaGradF.x() + dy*sigmaGradF.y())/pow(rSq, 1.5);
        // END GENERATED pseudoECGManufacturedIntegrand_d2
    }

    const scalar dx = x - electrodePosition.x();
    const scalar dy = y - electrodePosition.y();
    const scalar dz = z - electrodePosition.z();
    const scalar rSq = dx*dx + dy*dy + dz*dz;

    if (rSq <= VSMALL)
    {
        return 0.0;
    }

    // BEGIN GENERATED pseudoECGManufacturedIntegrand_d3
    return
        (dx*sigmaGradF.x() + dy*sigmaGradF.y() + dz*sigmaGradF.z())
       /pow(rSq, 1.5);
    // END GENERATED pseudoECGManufacturedIntegrand_d3
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
