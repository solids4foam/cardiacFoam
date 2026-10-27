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

\*---------------------------------------------------------------------------*/

#include "pvjContinuumReference.H"

#include "mathematicalConstants.H"
#include "monodomainVerification/manufacturedFDAReference.H"

namespace Foam
{

namespace
{

//- Gauss-Legendre nodes and weights on [0, 1]
void unitGaussRule(const label n, List<scalar>& x, List<scalar>& w)
{
    pseudoECGManufacturedQuadratureRule(n, x, w);
    forAll(x, i)
    {
        x[i] = 0.5*(x[i] + 1.0);
        w[i] *= 0.5;
    }
}


scalar tetDet(const point& a, const point& b, const point& c, const point& d)
{
    return mag((b - a) & ((c - a) ^ (d - a)));
}


scalar maxEdge(const point& a, const point& b, const point& c, const point& d)
{
    return max
    (
        max(max(mag(b - a), mag(c - a)), max(mag(d - a), mag(c - b))),
        max(mag(d - b), mag(d - c))
    );
}


//- Point on edge p-q where the linear interpolant of phi vanishes
point zeroOnEdge(const point& p, const point& q, scalar phiP, scalar phiQ)
{
    const scalar s = phiP/(phiP - phiQ);
    return p + s*(q - p);
}

} // End anonymous namespace


pvjContinuumReference::pvjContinuumReference
(
    const word& kernel,
    scalar R,
    const dictionary& dict
)
:
    kernel_(kernel),
    R_(R),
    nRadial_(dict.lookupOrDefault<label>("referenceRadialPoints", 48)),
    nPolar_(dict.lookupOrDefault<label>("referencePolarPoints", 48)),
    nAzimuthal_(dict.lookupOrDefault<label>("referenceAzimuthalPoints", 96)),
    // Kernels that jump at R (uniform, gaussian) need the finest cut; the
    // linear kernel vanishes at R and the cosine kernel's slope does too
    cutSize_
    (
        R*dict.lookupOrDefault<scalar>
        (
            "referenceCutSize",
            kernel == "cosine" ? 1.0/24.0
          : kernel == "linear" ? 1.0/48.0
          : 1.0/96.0
        )
    ),
    insideSize_(R*dict.lookupOrDefault<scalar>("referenceInsideSize", 0.25)),
    tipSize_(R*dict.lookupOrDefault<scalar>("referenceTipSize", 1.0/64.0))
{
    if
    (
        kernel_ != "uniform" && kernel_ != "linear"
     && kernel_ != "gaussian" && kernel_ != "cosine"
    )
    {
        FatalErrorInFunction
            << "No continuum definition of pvjKernel '" << kernel_
            << "'. Defined: uniform linear gaussian cosine."
            << exit(FatalError);
    }

    if (R_ <= 0)
    {
        FatalErrorInFunction
            << "Support radius must be positive; found " << R_ << "."
            << exit(FatalError);
    }
}


scalar pvjContinuumReference::extendedWeight(scalar r) const
{
    const scalar s = r/R_;

    if (kernel_ == "linear")
    {
        return 1.0 - s;
    }
    if (kernel_ == "gaussian")
    {
        return Foam::exp(-4.5*s*s);
    }
    if (kernel_ == "cosine")
    {
        return 0.5*(1.0 + Foam::cos(constant::mathematical::pi*s));
    }

    return 1.0;
}


scalar pvjContinuumReference::weight(scalar r) const
{
    return r <= R_ ? extendedWeight(r) : 0.0;
}


scalar pvjContinuumReference::supportIntegral
(
    const support& s,
    const scalarFunction& f
) const
{
    List<scalar> xr, wr, xm, wm;
    unitGaussRule(nRadial_, xr, wr);
    unitGaussRule(nPolar_, xm, wm);

    // Orthonormal frame with e3 along the normal into the tissue
    vector e3(0, 0, 1);
    if (s.kind != supportKind::interior)
    {
        e3 = s.normal/mag(s.normal);
    }
    vector trial = mag(e3.x()) < 0.9 ? vector(1, 0, 0) : vector(0, 1, 0);
    vector e1 = trial - (trial & e3)*e3;
    e1 /= mag(e1);
    const vector e2 = e3 ^ e1;

    const scalar dPhi = constant::mathematical::twoPi/scalar(nAzimuthal_);

    scalar sum = 0;

    forAll(xr, i)
    {
        const scalar r = R_*xr[i];
        const scalar wR = R_*wr[i]*r*r*weight(r);

        scalar muMin = -1.0;
        if (s.kind == supportKind::plane)
        {
            muMin = 0.0;
        }
        else if (s.kind == supportKind::cavitySphere)
        {
            muMin = -r/(2.0*s.wallRadius);
        }
        const scalar muSpan = 1.0 - muMin;

        scalar sumMu = 0;
        forAll(xm, j)
        {
            const scalar mu = muMin + muSpan*xm[j];
            const scalar sinTheta = Foam::sqrt(max(1.0 - mu*mu, 0.0));

            scalar sumPhi = 0;
            for (label k = 0; k < nAzimuthal_; ++k)
            {
                const scalar phi = (scalar(k) + 0.5)*dPhi;
                const point x =
                    s.centre
                  + r
                   *(
                        mu*e3
                      + sinTheta*(Foam::cos(phi)*e1 + Foam::sin(phi)*e2)
                    );
                sumPhi += f ? f(x) : 1.0;
            }

            sumMu += muSpan*wm[j]*dPhi*sumPhi;
        }

        sum += wR*sumMu;
    }

    return sum;
}


void pvjContinuumReference::quadratureTet
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
) const
{
    // 4 points per direction on whole tets, 2 on the small clipped pieces
    static List<scalar> x4, w4, x2, w2;
    if (x4.empty())
    {
        unitGaussRule(4, x4, w4);
        unitGaussRule(2, x2, w2);
    }
    const List<scalar>& x = clipped ? x2 : x4;
    const List<scalar>& w = clipped ? w2 : w4;

    const scalar det = tetDet(a, b, c, d);
    if (det <= VSMALL)
    {
        return;
    }

    const vector ab = b - a;
    const vector ac = c - a;
    const vector ad = d - a;

    scalar sumW = 0;
    scalar sumWF = 0;

    forAll(x, i)
    {
        const scalar u = x[i];
        forAll(x, j)
        {
            const scalar v = x[j];
            forAll(x, k)
            {
                const scalar l1 = u;
                const scalar l2 = (1.0 - u)*v;
                const scalar l3 = (1.0 - u)*(1.0 - v)*x[k];
                const scalar jac =
                    w[i]*w[j]*w[k]*(1.0 - u)*(1.0 - u)*(1.0 - v);

                const point p = a + l1*ab + l2*ac + l3*ad;
                const scalar wp = jac*extendedWeight(mag(p - t));

                sumW += wp;
                if (f)
                {
                    sumWF += wp*f(p);
                }
            }
        }
    }

    omega += sign*det*sumW;
    omegaF += sign*det*sumWF;
}


void pvjContinuumReference::clipTet
(
    const point& t,
    const point& a,
    const point& b,
    const point& c,
    const point& d,
    const scalarFunction& f,
    scalar& omega,
    scalar& omegaF
) const
{
    const point v[4] = {a, b, c, d};
    scalar phi[4];
    label nInside = 0;
    for (label i = 0; i < 4; ++i)
    {
        phi[i] = mag(v[i] - t) - R_;
        if (phi[i] < 0)
        {
            ++nInside;
        }
    }

    if (nInside == 0)
    {
        return;
    }

    if (nInside == 4)
    {
        quadratureTet(t, a, b, c, d, f, 1.0, true, omega, omegaF);
        return;
    }

    if (nInside == 1 || nInside == 3)
    {
        // The vertex alone on its side of the cut plane
        const bool loneInside = (nInside == 1);
        label lone = -1;
        for (label i = 0; i < 4; ++i)
        {
            if ((phi[i] < 0) == loneInside)
            {
                lone = i;
            }
        }

        point p[3];
        label n = 0;
        for (label i = 0; i < 4; ++i)
        {
            if (i != lone)
            {
                p[n++] = zeroOnEdge(v[lone], v[i], phi[lone], phi[i]);
            }
        }

        if (loneInside)
        {
            quadratureTet(t, v[lone], p[0], p[1], p[2], f, 1.0, true, omega, omegaF);
        }
        else
        {
            quadratureTet(t, a, b, c, d, f, 1.0, true, omega, omegaF);
            quadratureTet(t, v[lone], p[0], p[1], p[2], f, -1.0, true, omega, omegaF);
        }
        return;
    }

    // Two inside (i, j), two outside (k, l): a prism
    label in[2], out[2];
    label nIn = 0, nOut = 0;
    for (label m = 0; m < 4; ++m)
    {
        if (phi[m] < 0)
        {
            in[nIn++] = m;
        }
        else
        {
            out[nOut++] = m;
        }
    }

    const point& a0 = v[in[0]];
    const point a1 = zeroOnEdge(v[in[0]], v[out[0]], phi[in[0]], phi[out[0]]);
    const point a2 = zeroOnEdge(v[in[0]], v[out[1]], phi[in[0]], phi[out[1]]);
    const point& b0 = v[in[1]];
    const point b1 = zeroOnEdge(v[in[1]], v[out[0]], phi[in[1]], phi[out[0]]);
    const point b2 = zeroOnEdge(v[in[1]], v[out[1]], phi[in[1]], phi[out[1]]);

    quadratureTet(t, a0, a1, a2, b2, f, 1.0, true, omega, omegaF);
    quadratureTet(t, a0, a1, b1, b2, f, 1.0, true, omega, omegaF);
    quadratureTet(t, a0, b0, b1, b2, f, 1.0, true, omega, omegaF);
}


void pvjContinuumReference::integrateTet
(
    const point& t,
    const point& a,
    const point& b,
    const point& c,
    const point& d,
    const scalarFunction& f,
    scalar& omega,
    scalar& omegaF
) const
{
    const point g = 0.25*(a + b + c + d);
    const scalar rad = max
    (
        max(mag(a - g), mag(b - g)),
        max(mag(c - g), mag(d - g))
    );
    const scalar dg = mag(g - t);

    if (dg - rad >= R_)
    {
        return;
    }

    const scalar dMax = max
    (
        max(mag(a - t), mag(b - t)),
        max(mag(c - t), mag(d - t))
    );
    const scalar size = maxEdge(a, b, c, d);

    bool split = false;
    if (dMax <= R_)
    {
        split =
            size > insideSize_
         || (kernel_ == "linear" && dg < rad && size > tipSize_);
    }
    else
    {
        split = size > cutSize_;
    }

    if (!split)
    {
        if (dMax <= R_)
        {
            quadratureTet(t, a, b, c, d, f, 1.0, false, omega, omegaF);
        }
        else
        {
            clipTet(t, a, b, c, d, f, omega, omegaF);
        }
        return;
    }

    const point m01 = 0.5*(a + b);
    const point m02 = 0.5*(a + c);
    const point m03 = 0.5*(a + d);
    const point m12 = 0.5*(b + c);
    const point m13 = 0.5*(b + d);
    const point m23 = 0.5*(c + d);

    integrateTet(t, a, m01, m02, m03, f, omega, omegaF);
    integrateTet(t, b, m01, m12, m13, f, omega, omegaF);
    integrateTet(t, c, m02, m12, m23, f, omega, omegaF);
    integrateTet(t, d, m03, m13, m23, f, omega, omegaF);

    // Octahedron split about its shortest diagonal, which keeps the
    // children shape-regular so the longest edge halves every two levels
    const scalar d0 = mag(m01 - m23);
    const scalar d1 = mag(m02 - m13);
    const scalar d2 = mag(m03 - m12);

    // Diagonal (p, q) and the cycle of the other four midpoints about it
    const point* p = &m02;
    const point* q = &m13;
    const point* ring[4] = {&m01, &m03, &m23, &m12};
    if (d0 <= d1 && d0 <= d2)
    {
        p = &m01;
        q = &m23;
        ring[0] = &m02; ring[1] = &m03; ring[2] = &m13; ring[3] = &m12;
    }
    else if (d2 < d1)
    {
        p = &m03;
        q = &m12;
        ring[0] = &m01; ring[1] = &m02; ring[2] = &m23; ring[3] = &m13;
    }

    for (label k = 0; k < 4; ++k)
    {
        integrateTet
        (
            t, *p, *q, *ring[k], *ring[(k + 1) % 4], f, omega, omegaF
        );
    }
}


void pvjContinuumReference::cellIntegrals
(
    const fvMesh& mesh,
    const point& t,
    const scalarFunction& f,
    labelList& cells,
    scalarField& omega,
    scalarField& omegaF
) const
{
    const pointField& points = mesh.points();
    const faceList& faces = mesh.faces();
    const cellList& meshCells = mesh.cells();
    const vectorField& cellCentres = mesh.cellCentres();
    const vectorField& faceCentres = mesh.faceCentres();

    DynamicList<label> cellList;
    DynamicList<scalar> omegaList;
    DynamicList<scalar> omegaFList;

    forAll(meshCells, cellI)
    {
        const point& C = cellCentres[cellI];
        const labelList cellPoints = meshCells[cellI].labels(faces);

        scalar rad = 0;
        forAll(cellPoints, pI)
        {
            rad = max(rad, mag(points[cellPoints[pI]] - C));
        }

        if (mag(C - t) - rad >= R_)
        {
            continue;
        }

        scalar cellOmega = 0;
        scalar cellOmegaF = 0;

        const cell& c = meshCells[cellI];
        forAll(c, cfI)
        {
            const face& fc = faces[c[cfI]];
            const point& Cf = faceCentres[c[cfI]];

            forAll(fc, pI)
            {
                const point& p0 = points[fc[pI]];
                const point& p1 = points[fc.nextLabel(pI)];

                integrateTet(t, C, Cf, p0, p1, f, cellOmega, cellOmegaF);
            }
        }

        if (cellOmega != 0 || cellOmegaF != 0)
        {
            cellList.append(cellI);
            omegaList.append(cellOmega);
            omegaFList.append(cellOmegaF);
        }
    }

    cells.transfer(cellList);
    omega.transfer(omegaList);
    omegaF.transfer(omegaFList);
}

} // End namespace Foam

// ************************************************************************* //
