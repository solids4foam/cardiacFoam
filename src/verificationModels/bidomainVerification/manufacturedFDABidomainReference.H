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

#ifndef manufacturedFDABidomainReference_H
#define manufacturedFDABidomainReference_H

#include "monodomainVerification/manufacturedFDAReference.H"

namespace Foam
{

inline scalar manufacturedFDABidomainShift
(
    const scalar t,
    const label dimension,
    const scalar k,
    const vector& referencePoint,
    const scalar phiEReferenceValue
)
{
    scalarField x(1, referencePoint.x());
    scalarField y(1, referencePoint.y());
    scalarField z(1, referencePoint.z());
    scalarField F;

    computeF(F, x, y, z, dimension);

    return phiEReferenceValue + k*Foam::sqrt(1.0 + t)*F[0];
}

inline void computeManufacturedBidomainPhiE
(
    scalarField& phiE,
    const scalarField& x,
    const scalarField& y,
    const scalarField& z,
    const scalar t,
    const label dimension,
    const scalar k,
    const scalar shift
)
{
    scalarField F;
    computeF(F, x, y, z, dimension);

    phiE.setSize(x.size());

    const scalar coeff = -k*Foam::sqrt(1.0 + t);

    forAll(x, i)
    {
        phiE[i] = coeff*F[i] + shift;
    }
}

inline void computeManufacturedBidomainPhiI
(
    scalarField& phiI,
    const scalarField& x,
    const scalarField& y,
    const scalarField& z,
    const scalar t,
    const label dimension,
    const scalar k,
    const scalar shift
)
{
    scalarField V;
    scalarField phiE;

    computeManufacturedV(V, x, y, z, t, dimension);
    computeManufacturedBidomainPhiE(phiE, x, y, z, t, dimension, k, shift);

    phiI.setSize(x.size());

    forAll(x, i)
    {
        phiI[i] = V[i] + phiE[i];
    }
}

} // End namespace Foam

#endif

// ************************************************************************* //
