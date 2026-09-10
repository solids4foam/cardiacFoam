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
    manufacturedEikonalReference

Description
    Fixed manufactured references for eikonal activation and eikonal ECG.

Author
    Simao Nieto de Castro. All rights reserved.
\*---------------------------------------------------------------------------*/

#ifndef manufacturedEikonalReference_H
#define manufacturedEikonalReference_H

#include "PstreamReduceOps.H"
#include "fvMesh.H"
#include "mathematicalConstants.H"
#include "monodomainVerification/manufacturedFDAReference.H"
#include "tensorField.H"
#include "volFields.H"

namespace Foam
{

inline vector manufacturedEikonalK(const label dimension)
{
    if (dimension == 1) return vector(1.0, 0.0, 0.0);
    if (dimension == 2) return vector(1.0, 1.0, 0.0);
    return vector(1.0, 1.0, 1.0);
}


inline scalar manufacturedEikonalTemplateValue(const scalar s)
{
    return Foam::sin(2.0*constant::mathematical::pi*s);
}


inline scalar manufacturedEikonalTemplateDerivative(const scalar s)
{
    return
        2.0*constant::mathematical::pi
       *Foam::cos(2.0*constant::mathematical::pi*s);
}


inline tensor manufacturedEikonalM
(
    const tensor& conductivity,
    const scalar chi,
    const scalar Cm
)
{
    return conductivity/(chi*Cm);
}


inline scalar manufacturedEikonalTau
(
    const vector& point,
    const vector& k
)
{
    return Foam::exp(k & point);
}


inline scalar manufacturedEikonalTau
(
    const scalar x,
    const scalar y,
    const scalar z,
    const vector& k
)
{
    return manufacturedEikonalTau(vector(x, y, z), k);
}


inline vector manufacturedEikonalGradTau
(
    const vector& point,
    const vector& k
)
{
    return k * manufacturedEikonalTau(point, k);
}


inline scalar manufacturedEikonalSourceTerm
(
    const vector& point,
    const tensor& conductivity,
    const scalar chi,
    const scalar Cm,
    const scalar c0,
    const vector& k
)
{
    const tensor M = manufacturedEikonalM(conductivity, chi, Cm);
    const scalar KM = k & (M & k);

    if (KM <= VSMALL)
    {
        FatalErrorInFunction
            << "Manufactured eikonal verification requires positive "
            << "conductivity metric."
            << exit(FatalError);
    }

    const scalar tau = manufacturedEikonalTau(point, k);
    return tau * (c0 * Foam::sqrt(KM) - KM) - 1.0;
}


inline void validateManufacturedEikonalUnitDomain
(
    const fvMesh& mesh,
    const label dimension,
    const scalar tolerance = 1e-8
)
{
    const pointField& points = mesh.points();

    scalar minX = GREAT;
    scalar minY = GREAT;
    scalar minZ = GREAT;
    scalar maxX = -GREAT;
    scalar maxY = -GREAT;
    scalar maxZ = -GREAT;

    forAll(points, pointI)
    {
        minX = min(minX, points[pointI].x());
        minY = min(minY, points[pointI].y());
        minZ = min(minZ, points[pointI].z());
        maxX = max(maxX, points[pointI].x());
        maxY = max(maxY, points[pointI].y());
        maxZ = max(maxZ, points[pointI].z());
    }

    reduce(minX, minOp<scalar>());
    reduce(minY, minOp<scalar>());
    reduce(minZ, minOp<scalar>());
    reduce(maxX, maxOp<scalar>());
    reduce(maxY, maxOp<scalar>());
    reduce(maxZ, maxOp<scalar>());

    const bool validX = mag(minX) <= tolerance && mag(maxX - 1.0) <= tolerance;
    const bool validY =
        dimension < 2
     || (mag(minY) <= tolerance && mag(maxY - 1.0) <= tolerance);
    const bool validZ =
        dimension < 3
     || (mag(minZ) <= tolerance && mag(maxZ - 1.0) <= tolerance);

    if (!validX || !validY || !validZ)
    {
        FatalErrorInFunction
            << "Manufactured eikonal verification requires the unit domain "
            << "[0,1]^d. Mesh bounds are min=(" << minX << " " << minY
            << " " << minZ << "), max=(" << maxX << " " << maxY
            << " " << maxZ << "), dimension=" << dimension << "."
            << exit(FatalError);
    }
}


inline tensor manufacturedEikonalConstantConductivity
(
    const volTensorField& conductivity,
    const scalar tolerance = 1e-10
)
{
    const tensorField& conductivityField = conductivity.primitiveField();

    if (conductivityField.empty())
    {
        FatalErrorInFunction
            << "Manufactured eikonal verification requires a non-empty "
            << "conductivity field."
            << exit(FatalError);
    }

    const tensor reference = conductivityField[0];
    scalar maxDiff = 0.0;

    forAll(conductivityField, cellI)
    {
        maxDiff =
            max(maxDiff, Foam::magSqr(conductivityField[cellI] - reference));
    }

    reduce(maxDiff, maxOp<scalar>());

    if (maxDiff > tolerance)
    {
        FatalErrorInFunction
            << "Manufactured eikonal verification requires constant "
            << "conductivity. max(|sigma - sigma0|^2)=" << maxDiff
            << " exceeds tolerance " << tolerance << "."
            << exit(FatalError);
    }

    return reference;
}


inline scalar manufacturedEikonalECGIntegrand
(
    const scalar x,
    const scalar y,
    const scalar z,
    const scalar time,
    const vector& electrodePosition,
    const tensor& conductivity,
    const vector& k,
    const label dimension
)
{
    const vector p(x, y, z);
    const scalar tau = manufacturedEikonalTau(p, k);
    const vector gradTau = manufacturedEikonalGradTau(p, k);

    const scalar dU =
        manufacturedEikonalTemplateDerivative(time - tau);
    const vector sigmaGradVm = conductivity & (-dU*gradTau);

    const scalar dx = x - electrodePosition.x();
    const scalar dy = y - electrodePosition.y();
    const scalar dz = z - electrodePosition.z();

    if (dimension == 1)
    {
        const scalar absDx = Foam::mag(dx);

        if (absDx <= VSMALL)
        {
            return 0.0;
        }

        return sigmaGradVm.x()*dx/Foam::pow(absDx, 3);
    }

    if (dimension == 2)
    {
        const scalar rSq = dx*dx + dy*dy;

        if (rSq <= VSMALL)
        {
            return 0.0;
        }

        return
            (sigmaGradVm.x()*dx + sigmaGradVm.y()*dy)
           /Foam::pow(rSq, 1.5);
    }

    const scalar rSq = dx*dx + dy*dy + dz*dz;

    if (rSq <= VSMALL)
    {
        return 0.0;
    }

    return
        (
            sigmaGradVm.x()*dx
          + sigmaGradVm.y()*dy
          + sigmaGradVm.z()*dz
        )/Foam::pow(rSq, 1.5);
}


inline scalar computeManufacturedEikonalECGReference
(
    const scalar time,
    const vector& electrodePosition,
    const tensor& conductivity,
    const vector& k,
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
            value += weights[i]*manufacturedEikonalECGIntegrand
            (
                x, 0.0, 0.0, time, electrodePosition, conductivity,
                k, dimension
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
                   *manufacturedEikonalECGIntegrand
                    (
                        x, y, 0.0, time, electrodePosition, conductivity,
                        k, dimension
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
            forAll(nodes, kk)
            {
                const scalar z = jacobian*(nodes[kk] + 1.0);
                value += weights[i]*weights[j]*weights[kk]
                   *manufacturedEikonalECGIntegrand
                    (
                        x, y, z, time, electrodePosition, conductivity,
                        k, dimension
                    );
            }
        }
    }

    return jacobian*jacobian*jacobian*value;
}

} // End namespace Foam

#endif

// ************************************************************************* //
