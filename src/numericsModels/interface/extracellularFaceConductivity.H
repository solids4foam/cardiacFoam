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

Description
    Testable face-conductivity interpolation helpers for the global
    extracellular-potential equation.

Author
    Simao Nieto de Castro, UCD.
\*---------------------------------------------------------------------------*/

#ifndef extracellularFaceConductivity_H
#define extracellularFaceConductivity_H

#include "tensor.H"

namespace Foam
{
namespace extracellularFaceConductivity
{

//- Distance-weighted harmonic mean of two positive scalar conductivities at
//  distances distanceA and distanceB from the face; linear otherwise
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
        return distance/(distanceA/a + distanceB/b);
    }

    return
        distance > SMALL
      ? (distanceB*a + distanceA*b)/distance
      : 0.5*(a + b);
}


//- Face conductivity between two uniform materials a and b with unit face
//  normal n, from continuity of the potential and of the normal flux:
//  normal part the distance-weighted harmonic mean of n.a.n and n.b.n,
//  conormal part t = (kB ta + kA tb)/(kA + kB) with t = K.n - (n.K.n) n and
//  k = n.K.n/d, tangential block the mean of the two (it carries no flux)
inline tensor conormalHarmonic
(
    const tensor& a,
    const tensor& b,
    const vector& n,
    const scalar distanceA,
    const scalar distanceB
)
{
    const scalar aNN = n & a & n;
    const scalar bNN = n & b & n;
    const scalar kA = aNN/max(distanceA, SMALL);
    const scalar kB = bNN/max(distanceB, SMALL);
    const scalar kSum = kA + kB;
    const tensor P(tensor::I - sqr(n));
    const tensor tangential(P & (0.5*(a + b)) & P);

    if (kSum <= SMALL)
    {
        return tangential;
    }

    const vector t
    (
        (kB*((a & n) - aNN*n) + kA*((b & n) - bNN*n))/kSum
    );

    return
        kA*kB/kSum*(distanceA + distanceB)*sqr(n)
      + t*n
      + n*t
      + tangential;
}


//- Linear interpolation with owner weight weightA
inline tensor linear
(
    const tensor& a,
    const tensor& b,
    const scalar weightA
)
{
    return weightA*a + (1.0 - weightA)*b;
}

} // End namespace extracellularFaceConductivity
} // End namespace Foam

#endif

// ************************************************************************* //
