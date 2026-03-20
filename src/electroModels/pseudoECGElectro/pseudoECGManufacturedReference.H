#ifndef pseudoECGManufacturedReference_H
#define pseudoECGManufacturedReference_H

#include "List.H"
#include "Switch.H"
#include "mathematicalConstants.H"
#include "tensor.H"
#include "vector.H"

namespace Foam
{

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
           *(scalar(i) + 0.75)
           /(scalar(quadratureOrder) + 0.5)
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
