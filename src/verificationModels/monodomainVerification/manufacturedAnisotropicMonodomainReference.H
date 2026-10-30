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
    manufacturedAnisotropicMonodomainReference

Description
    Boundary-compatible manufactured reference solution for the monodomain
    equation with a constant, fully populated anisotropic conductivity tensor.

    The spatial factor

        F = sin^2(pi*x) sin^2(pi*y) sin^2(pi*z)

    has zero gradient on every face of the unit cube.  It therefore satisfies
    the homogeneous-flux boundary condition for any constant tensor.  The
    continuous PDE source is evaluated from the analytic Hessian; no discrete
    finite-volume operator is used to manufacture the source.

    The formulas between "// BEGIN GENERATED <name>" and
    "// END GENERATED <name>" are written from the sympy definitions in
    tutorials/manufacturedSolutions/generators: rotated_monodomain_mms.py (F,
    V, u1-u3, source) and pseudo_ecg_mms.py (grad F, pseudo-ECG kernel).
    Regenerate them with --write and compare with --check; the rest of the
    file is hand-written.

Author
    Simao Nieto de Castro. All rights reserved.
\*---------------------------------------------------------------------------*/

#ifndef manufacturedAnisotropicMonodomainReference_H
#define manufacturedAnisotropicMonodomainReference_H

#include "mathematicalConstants.H"
#include "scalarField.H"
#include "tensorField.H"
#include "vectorField.H"

namespace Foam
{

inline void computeAnisotropicManufacturedF
(
    scalarField& F,
    const vectorField& centres
)
{
    F.setSize(centres.size());
    const scalar pi = constant::mathematical::pi;

    forAll(centres, cellI)
    {
        const vector& c = centres[cellI];
        // BEGIN GENERATED computeAnisotropicManufacturedF
        F[cellI] = sqr(sin(pi*c.x()))*sqr(sin(pi*c.y()))*sqr(sin(pi*c.z()));
        // END GENERATED computeAnisotropicManufacturedF
    }
}


inline void computeAnisotropicManufacturedV
(
    scalarField& V,
    const vectorField& centres,
    const scalar t
)
{
    scalarField F;
    computeAnisotropicManufacturedF(F, centres);

    V.setSize(centres.size());

    forAll(V, cellI)
    {
        // BEGIN GENERATED computeAnisotropicManufacturedV
        V[cellI] = F[cellI]*sqrt(t + 1);
        // END GENERATED computeAnisotropicManufacturedV
    }
}


inline void computeAnisotropicManufacturedU
(
    scalarField& u1,
    scalarField& u2,
    scalarField& u3,
    const vectorField& centres,
    const scalar t
)
{
    scalarField F;
    computeAnisotropicManufacturedF(F, centres);

    u1.setSize(centres.size());
    u2.setSize(centres.size());
    u3.setSize(centres.size());

    forAll(centres, cellI)
    {
        const vector& c = centres[cellI];
        // BEGIN GENERATED computeAnisotropicManufacturedU
        const scalar G = c.x()*(c.y()*c.y())*(c.z()*c.z()*c.z()) + 1;
        u1[cellI] = F[cellI]*sqrt(t + 1) + G*(t + 1);
        u2[cellI] = 1/(sqrt(G)*(t + 1));
        u3[cellI] = 0;
        // END GENERATED computeAnisotropicManufacturedU
    }
}


inline void computeAnisotropicManufacturedPdeSource
(
    scalarField& source,
    const vectorField& centres,
    const tensorField& conductivity,
    const scalar t,
    const scalar beta
)
{
    if (conductivity.size() != centres.size())
    {
        FatalErrorInFunction
            << "Cannot evaluate the anisotropic manufactured source: "
            << "conductivity size " << conductivity.size()
            << " differs from the number of cell centres " << centres.size()
            << exit(FatalError);
    }

    source.setSize(centres.size());

    const scalar pi = constant::mathematical::pi;

    forAll(centres, cellI)
    {
        const vector& c = centres[cellI];
        const tensor& K = conductivity[cellI];

        // The solver uses +source on the PDE right-hand side.  The FDA
        // ionic model contributes -Cm*dV/dt + beta*V/chi, hence this source
        // is beta*V - div(K grad(V)).
        // BEGIN GENERATED computeAnisotropicManufacturedPdeSource
        const scalar w0 = pi*c.x();
        const scalar w1 = sin(w0);
        const scalar w2 = (w1*w1);
        const scalar w3 = pi*c.y();
        const scalar w4 = sin(w3);
        const scalar w5 = (w4*w4);
        const scalar w6 = pi*c.z();
        const scalar w7 = sin(w6);
        const scalar w8 = (w7*w7);
        const scalar w9 = w5*w8;
        const scalar w10 = w2*w9;
        const scalar w11 = sqr(pi);
        const scalar w12 = 2*w11;
        const scalar w13 = w10*w12;
        const scalar w14 = cos(w0);
        const scalar w15 = cos(w3);
        const scalar w16 = w12*w2;
        const scalar w17 = cos(w6);
        const scalar w18 = 4*w1*w11*w14;
        const scalar w19 = w15*w4;
        const scalar w20 = w18*w19*w8;
        const scalar w21 = w17*w7;
        const scalar w22 = w18*w21*w5;
        const scalar w23 = 4*w11*w19*w2*w21;
        const scalar F = w10;
        const scalar divKGradF =
            K.xx()*w12*(w14*w14)*w9
          - K.xx()*w13
          + K.xy()*w20
          + K.xz()*w22
          + K.yx()*w20
          - K.yy()*w13
          + K.yy()*(w15*w15)*w16*w8
          + K.yz()*w23
          + K.zx()*w22
          + K.zy()*w23
          - K.zz()*w13
          + K.zz()*w16*(w17*w17)*w5;
        source[cellI] = sqrt(t + 1)*(F*beta - divKGradF);
        // END GENERATED computeAnisotropicManufacturedPdeSource
    }
}

inline vector pseudoECGManufacturedGradFAnisotropic
(
    const scalar x,
    const scalar y,
    const scalar z
)
{
    const scalar pi = constant::mathematical::pi;

    // BEGIN GENERATED pseudoECGManufacturedGradFAnisotropic
    const scalar w0 = pi*y;
    const scalar w1 = sin(w0);
    const scalar w2 = (w1*w1);
    const scalar w3 = pi*x;
    const scalar w4 = sin(w3);
    const scalar w5 = pi*z;
    const scalar w6 = sin(w5);
    const scalar w7 = 2*pi*(w6*w6);
    const scalar w8 = (w4*w4);
    return vector(w2*w4*w7*cos(w3), w1*w7*w8*cos(w0), 2*pi*w2*w6*w8*cos(w5));
    // END GENERATED pseudoECGManufacturedGradFAnisotropic
}


inline scalar pseudoECGManufacturedIntegrandAnisotropic
(
    const scalar x,
    const scalar y,
    const scalar z,
    const vector& electrodePosition,
    const tensor& conductivity
)
{
    const vector sigmaGradF =
        conductivity & pseudoECGManufacturedGradFAnisotropic(x, y, z);

    const scalar dx = x - electrodePosition.x();
    const scalar dy = y - electrodePosition.y();
    const scalar dz = z - electrodePosition.z();
    const scalar rSq = dx*dx + dy*dy + dz*dz;

    if (rSq <= VSMALL)
    {
        return 0.0;
    }

    // BEGIN GENERATED pseudoECGManufacturedIntegrandAnisotropic
    return
        (dx*sigmaGradF.x() + dy*sigmaGradF.y() + dz*sigmaGradF.z())
       /pow(rSq, 1.5);
    // END GENERATED pseudoECGManufacturedIntegrandAnisotropic
}


inline scalar computeManufacturedPseudoECGSpatialReferenceAnisotropic
(
    const vector& electrodePosition,
    const tensor& conductivity,
    const UList<scalar>& nodes,
    const UList<scalar>& weights
)
{
    const scalar jacobian = 0.5;
    scalar value = 0.0;

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
                    * pseudoECGManufacturedIntegrandAnisotropic
                    (
                        x,
                        y,
                        z,
                        electrodePosition,
                        conductivity
                    );
            }
        }
    }

    return jacobian*jacobian*jacobian*value;
}

} // End namespace Foam

#endif

// ************************************************************************* //
