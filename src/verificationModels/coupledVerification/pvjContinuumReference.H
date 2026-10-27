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
    pvjContinuumReference

Description
    Continuum definition of the PVJ map, independent of pvjMapper. A
    terminal at t with support radius R and radial kernel w exchanges with
    the tissue through

        W      = int_{Omega cap B} w dx                (normalisation)
        <f>    = int_{Omega cap B} w f dx / W          (gathered value)
        S(x)   = I w(|x - t|) / W                      (deposited source)

    over the part of the ball B = {|x - t| <= R} inside the true tissue
    Omega. The kernels are written here from their definitions, not taken
    from pvjMapper, so a kernel defect there is not reproduced here.

    supportIntegral() integrates over Omega cap B in spherical coordinates
    about t with Gauss-Legendre in r and in mu = cos(theta) and the periodic
    trapezoidal rule in phi. Omega cap B is exact for three support kinds:
    interior (B inside Omega), plane (t on a flat wall, mu >= 0 about the
    inward normal) and cavitySphere (t on a sphere of radius rho with the
    tissue outside it, mu >= -r/(2 rho) about the outward radius).

    cellIntegrals() integrates w and w f over K cap B for each mesh cell K:
    K is split into tetrahedra (cell centre, face centre, face edge), tets
    cut by the sphere |x - t| = R are subdivided down to cutSize and then
    clipped by the linear interpolant of |x - t| - R, and each piece is
    integrated with a collapsed (Duffy) Gauss rule, 4x4x4 points on whole
    tets and 2x2x2 on clipped pieces. The cell average
    of the continuum source is then I*omega_K/(W*|K|).

SourceFiles
    pvjContinuumReference.C

\*---------------------------------------------------------------------------*/

#ifndef pvjContinuumReference_H
#define pvjContinuumReference_H

#include "fvMesh.H"
#include "scalarField.H"

#include <functional>

namespace Foam
{

class pvjContinuumReference
{
public:

    enum class supportKind
    {
        interior,
        plane,
        cavitySphere
    };

    //- Geometry of Omega cap B about one terminal
    struct support
    {
        point centre;
        //- Unit normal into the tissue (plane, cavitySphere)
        vector normal;
        //- Sphere radius (cavitySphere)
        scalar wallRadius;
        supportKind kind;
    };

    typedef std::function<scalar(const point&)> scalarFunction;

private:

    word kernel_;
    scalar R_;

    label nRadial_;
    label nPolar_;
    label nAzimuthal_;

    //- Largest tet edge for a tet cut by |x - t| = R before clipping
    scalar cutSize_;
    //- Largest tet edge for a tet inside the ball
    scalar insideSize_;
    //- Largest tet edge for a tet near t (the linear kernel's cone tip)
    scalar tipSize_;

    void integrateTet
    (
        const point& t,
        const point& a,
        const point& b,
        const point& c,
        const point& d,
        const scalarFunction& f,
        scalar& omega,
        scalar& omegaF
    ) const;

    //- Collapsed Gauss rule on one tet, kernel extended past R; the
    //  coarser rule on clipped pieces
    void quadratureTet
    (
        const point& t,
        const point& a,
        const point& b,
        const point& c,
        const point& d,
        const scalarFunction& f,
        scalar sign,
        bool clipped,
        scalar& omega,
        scalar& omegaF
    ) const;

    void clipTet
    (
        const point& t,
        const point& a,
        const point& b,
        const point& c,
        const point& d,
        const scalarFunction& f,
        scalar& omega,
        scalar& omegaF
    ) const;

public:

    pvjContinuumReference
    (
        const word& kernel,
        scalar R,
        const dictionary& dict
    );

    const word& kernel() const
    {
        return kernel_;
    }

    scalar radius() const
    {
        return R_;
    }

    //- Continuum kernel w(r), zero for r > R
    scalar weight(scalar r) const;

    //- The kernel formula continued past R (used inside clipped pieces)
    scalar extendedWeight(scalar r) const;

    //- int_{Omega cap B} w f dx; f = 1 when f is empty
    scalar supportIntegral
    (
        const support& s,
        const scalarFunction& f = scalarFunction()
    ) const;

    //- Local cells K meeting B, with int_{K cap B} w and int_{K cap B} w f
    void cellIntegrals
    (
        const fvMesh& mesh,
        const point& t,
        const scalarFunction& f,
        labelList& cells,
        scalarField& omega,
        scalarField& omegaF
    ) const;
};

} // End namespace Foam

#endif

// ************************************************************************* //
