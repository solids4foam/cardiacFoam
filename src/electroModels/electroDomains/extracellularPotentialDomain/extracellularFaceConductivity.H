/*---------------------------------------------------------------------------*\
License
    This file is part of cardiacFoam.

    cardiacFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

Description
    Testable face-conductivity interpolation helpers for the global
    extracellular-potential equation.

\*---------------------------------------------------------------------------*/

#ifndef extracellularFaceConductivity_H
#define extracellularFaceConductivity_H

#include "tensor.H"

namespace Foam
{
namespace extracellularFaceConductivity
{

inline scalar unweightedHarmonic(const scalar a, const scalar b)
{
    const scalar denominator = a + b;

    return
        (mag(denominator) > SMALL && a*b > 0.0)
      ? 2.0*a*b/denominator
      : 0.5*denominator;
}


inline scalar distanceWeightedHarmonic
(
    const scalar a,
    const scalar b,
    const scalar distanceA,
    const scalar distanceB
)
{
    const scalar distance = distanceA + distanceB;

    if
    (
        distance > SMALL
     && distanceA >= 0.0
     && distanceB >= 0.0
     && a*b > 0.0
    )
    {
        const scalar resistance = distanceA/a + distanceB/b;

        if (mag(resistance) > SMALL)
        {
            return distance/resistance;
        }
    }

    return
        distance > SMALL
      ? (distanceB*a + distanceA*b)/distance
      : 0.5*(a + b);
}


inline tensor unweightedHarmonic(const tensor& a, const tensor& b)
{
    tensor result(tensor::zero);

    for (direction d = 0; d < tensor::nComponents; ++d)
    {
        result[d] = unweightedHarmonic(a[d], b[d]);
    }

    return result;
}


inline tensor distanceWeightedHarmonic
(
    const tensor& a,
    const tensor& b,
    const scalar distanceA,
    const scalar distanceB
)
{
    tensor result(tensor::zero);

    for (direction d = 0; d < tensor::nComponents; ++d)
    {
        result[d] = distanceWeightedHarmonic
        (
            a[d],
            b[d],
            distanceA,
            distanceB
        );
    }

    return result;
}


inline tensor linear
(
    const tensor& a,
    const tensor& b,
    const scalar weightA
)
{
    return weightA*a + (1.0 - weightA)*b;
}


inline tensor unweightedExtracellularFaceTensor
(
    const tensor& sigmaTotalA,
    const tensor& sigmaTotalB,
    const tensor& sigmaIntracellularA,
    const tensor& sigmaIntracellularB
)
{
    const bool aIsHeart = magSqr(sigmaIntracellularA) > SMALL;
    const bool bIsHeart = magSqr(sigmaIntracellularB) > SMALL;

    if (aIsHeart != bIsHeart)
    {
        const tensor heartSigmaE =
            aIsHeart
          ? sigmaTotalA - sigmaIntracellularA
          : sigmaTotalB - sigmaIntracellularB;
        const tensor bathSigma = aIsHeart ? sigmaTotalB : sigmaTotalA;

        return unweightedHarmonic(heartSigmaE, bathSigma);
    }

    return unweightedHarmonic(sigmaTotalA, sigmaTotalB);
}


inline tensor distanceWeightedExtracellularFaceTensor
(
    const tensor& sigmaTotalA,
    const tensor& sigmaTotalB,
    const tensor& sigmaIntracellularA,
    const tensor& sigmaIntracellularB,
    const scalar weightA
)
{
    const bool aIsHeart = magSqr(sigmaIntracellularA) > SMALL;
    const bool bIsHeart = magSqr(sigmaIntracellularB) > SMALL;

    if (aIsHeart != bIsHeart)
    {
        const tensor coefficientA =
            aIsHeart
          ? sigmaTotalA - sigmaIntracellularA
          : sigmaTotalA;
        const tensor coefficientB =
            bIsHeart
          ? sigmaTotalB - sigmaIntracellularB
          : sigmaTotalB;

        return distanceWeightedHarmonic
        (
            coefficientA,
            coefficientB,
            1.0 - weightA,
            weightA
        );
    }

    return linear(sigmaTotalA, sigmaTotalB, weightA);
}

} // End namespace extracellularFaceConductivity
} // End namespace Foam

#endif

// ************************************************************************* //
