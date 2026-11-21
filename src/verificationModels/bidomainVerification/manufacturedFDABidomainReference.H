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
    manufacturedFDABidomainReference

Description
    Manufactured FDA bidomain reference solution: phi_e = -k V_m + c and
    phi_i = V_m + phi_e on the FDA monodomain field V_m = sqrt(1 + t) F_d
    (monodomainVerification/manufacturedFDAReference.H), with
    c = phi_e,ref + k V_m(x_ref, t).

    The formulas between "// BEGIN GENERATED <name>" and
    "// END GENERATED <name>" are written from the sympy definitions in
    tutorials/manufacturedSolutions/generators/fda_bidomain_mms.py.
    Regenerate them with --write and compare with --check; the rest of the
    file is hand-written.

Author
    Simao Nieto de Castro. All rights reserved.
\*---------------------------------------------------------------------------*/

#ifndef manufacturedFDABidomainReference_H
#define manufacturedFDABidomainReference_H

#include "monodomainVerification/manufacturedFDAReference.H"

namespace Foam
{

//- Constant c of phi_e = -k V_m + c that gives phi_e(x_ref, t) = phi_e,ref
inline scalar manufacturedFDABidomainShift
(
    const scalar t,
    const label dimension,
    const scalar k,
    const vector& referencePoint,
    const scalar phiEReferenceValue
)
{
    const scalar Fref = manufacturedFDAFAt(referencePoint, dimension);

    // BEGIN GENERATED manufacturedFDABidomainShift
    return Fref*k*sqrt(t + 1) + phiEReferenceValue;
    // END GENERATED manufacturedFDABidomainShift
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

    forAll(x, i)
    {
        // BEGIN GENERATED computeManufacturedBidomainPhiE
        phiE[i] = -F[i]*k*sqrt(t + 1) + shift;
        // END GENERATED computeManufacturedBidomainPhiE
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
        // BEGIN GENERATED computeManufacturedBidomainPhiI
        phiI[i] = V[i] + phiE[i];
        // END GENERATED computeManufacturedBidomainPhiI
    }
}


//- V_m at point p
inline scalar manufacturedFDABidomainVmAt
(
    const point& p,
    const scalar t,
    const label dimension
)
{
    const scalar F = manufacturedFDAFAt(p, dimension);

    // BEGIN GENERATED manufacturedFDABidomainVmAt
    return F*sqrt(t + 1);
    // END GENERATED manufacturedFDABidomainVmAt
}


//- phi_e at point p
inline scalar manufacturedFDABidomainPhiEAt
(
    const point& p,
    const scalar t,
    const label dimension,
    const scalar k,
    const scalar shift
)
{
    const scalar F = manufacturedFDAFAt(p, dimension);

    // BEGIN GENERATED manufacturedFDABidomainPhiEAt
    return -F*k*sqrt(t + 1) + shift;
    // END GENERATED manufacturedFDABidomainPhiEAt
}


//- phi_i at point p
inline scalar manufacturedFDABidomainPhiIAt
(
    const point& p,
    const scalar t,
    const label dimension,
    const scalar k,
    const scalar shift
)
{
    const scalar F = manufacturedFDAFAt(p, dimension);

    // BEGIN GENERATED manufacturedFDABidomainPhiIAt
    return -F*k*sqrt(t + 1) + F*sqrt(t + 1) + shift;
    // END GENERATED manufacturedFDABidomainPhiIAt
}


//- grad V_m at point p
inline vector manufacturedFDABidomainGradVm
(
    const point& p,
    const scalar t,
    const label dimension
)
{
    const vector gradF =
        pseudoECGManufacturedGradF(p.x(), p.y(), p.z(), dimension);

    // BEGIN GENERATED manufacturedFDABidomainGradVm
    return gradF*sqrt(t + 1);
    // END GENERATED manufacturedFDABidomainGradVm
}


//- grad phi_e at point p
inline vector manufacturedFDABidomainGradPhiE
(
    const point& p,
    const scalar t,
    const label dimension,
    const scalar k
)
{
    const vector gradF =
        pseudoECGManufacturedGradF(p.x(), p.y(), p.z(), dimension);

    // BEGIN GENERATED manufacturedFDABidomainGradPhiE
    return -gradF*k*sqrt(t + 1);
    // END GENERATED manufacturedFDABidomainGradPhiE
}


//- grad phi_i at point p
inline vector manufacturedFDABidomainGradPhiI
(
    const point& p,
    const scalar t,
    const label dimension,
    const scalar k
)
{
    const vector gradF =
        pseudoECGManufacturedGradF(p.x(), p.y(), p.z(), dimension);

    // BEGIN GENERATED manufacturedFDABidomainGradPhiI
    return gradF*(-k*sqrt(t + 1) + sqrt(t + 1));
    // END GENERATED manufacturedFDABidomainGradPhiI
}

} // End namespace Foam

#endif

// ************************************************************************* //
