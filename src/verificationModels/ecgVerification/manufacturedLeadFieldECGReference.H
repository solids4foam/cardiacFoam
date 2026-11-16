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
    Exact passive-conductor references of manufacturedLeadFieldECGVerifier.

    slab: heart 0 < x < 1 between torso slabs, insulated outer boundary,
    Gi and Ge uniform with no xy or xz components, and
        Vm  = sqrt(1+t) cos(pi x)                      (heart)
        phi = -k sqrt(1+t) cos(pi x) + C               (heart)
        phi = -k sqrt(1+t) + C   (x < 0),   k sqrt(1+t) + C   (x > 1)
    with k = Gi_xx/(Gi_xx + Ge_xx) and C = k sqrt(1+t); the verifier's gauge
    removes C.

    layered: the slab with each torso side split at x = xL < 0 and
    x = xR > 1 into two layers of isotropic conductivity, and the current
    density -alpha along x entering through the surface currents -alpha on
    the x = min face and alpha on the x = max face. The potential is the
    slab's plus the piecewise linear
        alpha R(x),   R(x) = integral from 0 to x of dx'/sigma(x')
    with sigma = (Gi + Ge)_xx in the heart and the layer conductivity in
    the torso: continuous, with current -alpha in every layer. Vm is the
    slab's.

    wallCurrent: the slab with Vm = sqrt(1+t) (cos(pi x) + x/2), so that
    n.Gi grad(Vm) is non-zero on both walls; phi = -k Vm + C in the heart
    and the wall values in the insulated torso slabs (zero total current).

    The formulas between "// BEGIN GENERATED <name>" and
    "// END GENERATED <name>" are written from the sympy definitions in
    tutorials/manufacturedSolutions/generators/slab_leadfield_mms.py.
    Regenerate them with --write and compare with --check; the rest of the
    file is hand-written.

Author
    Simao Nieto de Castro. All rights reserved.
\*---------------------------------------------------------------------------*/

#ifndef manufacturedLeadFieldECGReference_H
#define manufacturedLeadFieldECGReference_H

#include "mathematicalConstants.H"
#include "verificationUtils.H"
#include "volFields.H"

namespace Foam
{

//- Heart Vm of the slab at point p and time t
inline scalar manufacturedLeadFieldSlabVm(const point& p, const scalar t)
{
    const scalar pi = constant::mathematical::pi;
    const scalar x = p.x();

    // BEGIN GENERATED manufacturedLeadFieldSlabVm
    return sqrt(t + 1)*cos(pi*x);
    // END GENERATED manufacturedLeadFieldSlabVm
}


//- Potential of the slab at point p and time t, up to a constant
inline scalar manufacturedLeadFieldSlabPhi
(
    const point& p,
    const scalar t,
    const scalar k
)
{
    const scalar pi = constant::mathematical::pi;
    const scalar x = p.x();

    if (x <= 0)
    {
        // BEGIN GENERATED slabPhi_leftTorso
        return 0;
        // END GENERATED slabPhi_leftTorso
    }
    if (x <= 1)
    {
        // BEGIN GENERATED slabPhi_heart
        return -k*sqrt(t + 1)*cos(pi*x) + k*sqrt(t + 1);
        // END GENERATED slabPhi_heart
    }

    // BEGIN GENERATED slabPhi_rightTorso
    return 2*k*sqrt(t + 1);
    // END GENERATED slabPhi_rightTorso
}


//- Heart Vm of wallCurrent at point p and time t: the slab's Vm with a
//  wall slope, so that the wall layer n.Gi grad(Vm) is non-zero
inline scalar manufacturedLeadFieldWallVm(const point& p, const scalar t)
{
    const scalar pi = constant::mathematical::pi;
    const scalar x = p.x();

    // BEGIN GENERATED manufacturedLeadFieldWallVm
    return sqrt(t + 1)*((1.0/2.0)*x + cos(pi*x));
    // END GENERATED manufacturedLeadFieldWallVm
}


//- Potential of wallCurrent at point p and time t, up to a constant: -k Vm
//  in the heart, the wall values in the insulated torso slabs
inline scalar manufacturedLeadFieldWallPhi
(
    const point& p,
    const scalar t,
    const scalar k
)
{
    const scalar pi = constant::mathematical::pi;
    const scalar x = p.x();

    if (x <= 0)
    {
        // BEGIN GENERATED wallPhi_leftTorso
        return 0;
        // END GENERATED wallPhi_leftTorso
    }
    if (x <= 1)
    {
        // BEGIN GENERATED wallPhi_heart
        return -k*sqrt(t + 1)*((1.0/2.0)*x + cos(pi*x)) + k*sqrt(t + 1);
        // END GENERATED wallPhi_heart
    }

    // BEGIN GENERATED wallPhi_rightTorso
    return (3.0/2.0)*k*sqrt(t + 1);
    // END GENERATED wallPhi_rightTorso
}


//- k = Gi_xx/(Gi_xx + Ge_xx) of uniform Gi and Ge without xy or xz
//  components; the heart mesh must span 0 <= x <= 1
inline scalar manufacturedLeadFieldSlabK
(
    const volTensorField& Gi,
    const volTensorField& Ge
)
{
    const fvMesh& heart = Gi.mesh();
    const scalar xMin = gMin(heart.points().component(vector::X)());
    const scalar xMax = gMax(heart.points().component(vector::X)());
    const scalar geometryTol = 1e-12;

    if (!(mag(xMin) <= geometryTol && mag(xMax - 1.0) <= geometryTol))
    {
        FatalErrorInFunction
            << "The slab reference needs the heart on 0 <= x <= 1; the heart "
            << "mesh spans " << xMin << " <= x <= " << xMax << "."
            << exit(FatalError);
    }

    const tensor gi =
        verificationUtils::uniformTensor(Gi.primitiveField(), "slab reference");
    const tensor ge =
        verificationUtils::uniformTensor(Ge.primitiveField(), "slab reference");

    // x-coupling components only; yz and zy do not enter
    const scalar crossTol = 1e-12*(mag(gi) + mag(ge));
    const scalar cross =
        mag(gi.xy()) + mag(gi.xz()) + mag(gi.yx()) + mag(gi.zx())
      + mag(ge.xy()) + mag(ge.xz()) + mag(ge.yx()) + mag(ge.zx());

    if (!(cross <= crossTol) || !(gi.xx() > 0) || !(ge.xx() > 0))
    {
        FatalErrorInFunction
            << "The slab reference needs Gi and Ge with positive xx and no "
            << "xy or xz components; Gi " << gi << ", Ge " << ge << "."
            << exit(FatalError);
    }

    return gi.xx()/(gi.xx() + ge.xx());
}


//- Layer interfaces and the x conductivities of the layered slab
struct manufacturedLeadFieldLayers
{
    //- Interfaces x = xLeft < 0 and x = xRight > 1
    scalar xLeft;
    scalar xRight;

    //- Torso layers x < xLeft, xLeft < x < 0, 1 < x < xRight, x > xRight
    scalar sigmaOuterLeft;
    scalar sigmaInnerLeft;
    scalar sigmaInnerRight;
    scalar sigmaOuterRight;

    //- (Gi + Ge)_xx
    scalar sigmaHeart;
};


//- R(x) = integral from 0 to x of dx'/sigma(x') of the layered slab
inline scalar manufacturedLeadFieldLayeredResistance
(
    const scalar x,
    const manufacturedLeadFieldLayers& l
)
{
    if (x <= l.xLeft)
    {
        // BEGIN GENERATED resistance_outerLeft
        return (-l.xLeft + x)/l.sigmaOuterLeft + l.xLeft/l.sigmaInnerLeft;
        // END GENERATED resistance_outerLeft
    }
    if (x <= 0)
    {
        // BEGIN GENERATED resistance_innerLeft
        return x/l.sigmaInnerLeft;
        // END GENERATED resistance_innerLeft
    }
    if (x <= 1)
    {
        // BEGIN GENERATED resistance_heart
        return x/l.sigmaHeart;
        // END GENERATED resistance_heart
    }
    if (x <= l.xRight)
    {
        // BEGIN GENERATED resistance_innerRight
        return (x - 1)/l.sigmaInnerRight + (1.0/(l.sigmaHeart));
        // END GENERATED resistance_innerRight
    }

    // BEGIN GENERATED resistance_outerRight
    return
        (-l.xRight + x)/l.sigmaOuterRight
      + (l.xRight - 1)/l.sigmaInnerRight
      + (1.0/(l.sigmaHeart));
    // END GENERATED resistance_outerRight
}


//- Potential of the layered slab at point p and time t, up to a constant
inline scalar manufacturedLeadFieldLayeredPhi
(
    const point& p,
    const scalar t,
    const scalar k,
    const scalar alpha,
    const manufacturedLeadFieldLayers& l
)
{
    const scalar phiSlab = manufacturedLeadFieldSlabPhi(p, t, k);
    const scalar R = manufacturedLeadFieldLayeredResistance(p.x(), l);

    // BEGIN GENERATED manufacturedLeadFieldLayeredPhi
    return R*alpha + phiSlab;
    // END GENERATED manufacturedLeadFieldLayeredPhi
}

} // End namespace Foam

#endif

// ************************************************************************* //
