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
    manufacturedConormalBathBidomainReference

Description
    Manufactured conormal bath bidomain reference solution.

Author
    Simao Nieto de Castro. All rights reserved.
\*---------------------------------------------------------------------------*/

#ifndef manufacturedConormalBathBidomainReference_H
#define manufacturedConormalBathBidomainReference_H

#include "complex.H"
#include "conormalManufacturedProfile.H"
#include "mathematicalConstants.H"
#include "monodomainVerification/manufacturedFDAReference.H"
#include "scalarField.H"
#include "scalarMatrices.H"
#include "tensor.H"
#include "tensorField.H"
#include "vectorField.H"

namespace Foam
{

inline complex conormalBathComplexSqrt(const complex& z)
{
    const scalar r = Foam::sqrt(Foam::sqr(z.Re()) + Foam::sqr(z.Im()));
    const scalar re = Foam::sqrt(0.5*(r + z.Re()));
    const scalar im = Foam::sqrt(0.5*(r - z.Re()));
    return complex(re, z.Im() >= 0 ? im : -im);
}


inline complex conormalBathComplexExp(const complex& z)
{
    const scalar e = Foam::exp(z.Re());
    return complex(e*Foam::cos(z.Im()), e*Foam::sin(z.Im()));
}


struct conormalBathCoefficients
{
    complex lambda1;
    complex lambda2;
    complex v1;
    complex v2;
    complex h1;
    complex h2;
    complex bL;
    complex bR;
};


inline void conormalBathAddTerm
(
    scalarSquareMatrix& M,
    const label rowRe,
    const label rowIm,
    const label col,
    const complex& coeff
)
{
    M(rowRe, col)     += coeff.Re();
    M(rowRe, col + 1) += -coeff.Im();
    M(rowIm, col)     += coeff.Im();
    M(rowIm, col + 1) += coeff.Re();
}


inline void conormalBathAddRhs
(
    scalarField& rhs,
    const label rowRe,
    const label rowIm,
    const complex& value
)
{
    rhs[rowRe] += value.Re();
    rhs[rowIm] += value.Im();
}


inline void conormalBathBuildInterfaceEquations
(
    scalarSquareMatrix& M,
    scalarField& rhs,
    const label row0,
    const scalar X,
    const scalar xOut,
    const label bCol,
    const scalar G11,
    const scalar G12,
    const scalar kappa,
    const scalar k,
    const scalar sigmaB,
    const scalar omega,
    const complex& lambda1,
    const complex& lambda2
)
{
    const label colV1 = 0;
    const label colV2 = 2;
    const label colH1 = 4;
    const label colH2 = 6;

    const complex e1 = conormalBathComplexExp(lambda1*X);
    const complex e2 = conormalBathComplexExp(lambda2*X);

    const complex iOmegaG12(0, omega*G12);

    const complex CvCoeffV1(G11, omega*G12*X);
    const complex CvCoeffV2(2.0*G11*X, omega*G12*X*X);
    const complex CvConst = iOmegaG12;

    const complex ChCoeffH1 = e1*(complex(G11, 0)*lambda1 + iOmegaG12);
    const complex ChCoeffH2 = e2*(complex(G11, 0)*lambda2 + iOmegaG12);

    const scalar coshv = Foam::cosh(omega*(X - xOut));
    const scalar sinhv = Foam::sinh(omega*(X - xOut));

    // (i) sealed conormal condition: (1-k)*C(v) + C(h) = 0
    {
        const label rowRe = 2*(row0 + 0);
        const label rowIm = rowRe + 1;

        conormalBathAddTerm(M, rowRe, rowIm, colV1, (1.0 - k)*CvCoeffV1);
        conormalBathAddTerm(M, rowRe, rowIm, colV2, (1.0 - k)*CvCoeffV2);
        conormalBathAddTerm(M, rowRe, rowIm, colH1, ChCoeffH1);
        conormalBathAddTerm(M, rowRe, rowIm, colH2, ChCoeffH2);

        conormalBathAddRhs(rhs, rowRe, rowIm, -(1.0 - k)*CvConst);
    }

    // (ii) phiE continuity: -k*v(X) + h(X) - b*cosh(omega*(X-xOut)) = 0
    {
        const label rowRe = 2*(row0 + 1);
        const label rowIm = rowRe + 1;

        conormalBathAddTerm(M, rowRe, rowIm, colV1, complex(-k*X, 0));
        conormalBathAddTerm(M, rowRe, rowIm, colV2, complex(-k*X*X, 0));
        conormalBathAddTerm(M, rowRe, rowIm, colH1, e1);
        conormalBathAddTerm(M, rowRe, rowIm, colH2, e2);
        conormalBathAddTerm(M, rowRe, rowIm, bCol, complex(-coshv, 0));

        conormalBathAddRhs(rhs, rowRe, rowIm, complex(k, 0));
    }

    // (iii) flux matching:
    // kappa*[-k*C(v) + C(h)] - sigmaB*omega*b*sinh(omega*(X-xOut)) = 0
    {
        const label rowRe = 2*(row0 + 2);
        const label rowIm = rowRe + 1;

        conormalBathAddTerm(M, rowRe, rowIm, colV1, kappa*(-k)*CvCoeffV1);
        conormalBathAddTerm(M, rowRe, rowIm, colV2, kappa*(-k)*CvCoeffV2);
        conormalBathAddTerm(M, rowRe, rowIm, colH1, kappa*ChCoeffH1);
        conormalBathAddTerm(M, rowRe, rowIm, colH2, kappa*ChCoeffH2);
        conormalBathAddTerm
        (
            M, rowRe, rowIm, bCol, complex(-sigmaB*omega*sinhv, 0)
        );

        conormalBathAddRhs(rhs, rowRe, rowIm, kappa*k*CvConst);
    }
}


inline conormalBathCoefficients conormalBathSolveCoefficients
(
    const tensor& Gi,
    const scalar sigmaB
)
{
    const scalar omega = 2.0*constant::mathematical::pi;
    const scalar G11 = Gi.xx();
    const scalar G12 = Gi.xy();
    const scalar G22 = Gi.yy();
    const scalar kappa = Foam::sqrt(2.0) - 1.0;
    const scalar k = 1.0/(1.0 + kappa);

    const complex A(G11, 0);
    const complex B(0, 2.0*omega*G12);
    const complex C(-omega*omega*G22, 0);
    const complex disc = B*B - complex(4, 0)*A*C;
    const complex sq = conormalBathComplexSqrt(disc);

    const complex lambda1 = (-B + sq)/(2.0*A);
    const complex lambda2 = (-B - sq)/(2.0*A);

    scalarSquareMatrix M(12, 0.0);
    scalarField rhs(12, 0.0);

    conormalBathBuildInterfaceEquations
    (
        M, rhs, 0, 0.0, -1.0, 8, G11, G12, kappa, k, sigmaB, omega,
        lambda1, lambda2
    );
    conormalBathBuildInterfaceEquations
    (
        M, rhs, 3, 1.0, 2.0, 10, G11, G12, kappa, k, sigmaB, omega,
        lambda1, lambda2
    );

    Foam::solve(M, rhs);

    conormalBathCoefficients coeffs;
    coeffs.lambda1 = lambda1;
    coeffs.lambda2 = lambda2;
    coeffs.v1 = complex(rhs[0], rhs[1]);
    coeffs.v2 = complex(rhs[2], rhs[3]);
    coeffs.h1 = complex(rhs[4], rhs[5]);
    coeffs.h2 = complex(rhs[6], rhs[7]);
    coeffs.bL = complex(rhs[8], rhs[9]);
    coeffs.bR = complex(rhs[10], rhs[11]);

    return coeffs;
}


inline scalar conormalBathH
(
    const conormalBathCoefficients& coeffs,
    const vector& c
)
{
    const scalar theta = 2.0*constant::mathematical::pi*c.y();
    const complex phase(Foam::cos(theta), Foam::sin(theta));
    const complex h =
        coeffs.h1*conormalBathComplexExp(coeffs.lambda1*c.x())
      + coeffs.h2*conormalBathComplexExp(coeffs.lambda2*c.x());

    return (h*phase).Re();
}


inline scalar conormalBathBathPotential
(
    const complex& b,
    const scalar xOut,
    const vector& c
)
{
    const scalar omega = 2.0*constant::mathematical::pi;
    const scalar theta = omega*c.y();
    const complex phase(Foam::cos(theta), Foam::sin(theta));
    const scalar coshv = Foam::cosh(omega*(c.x() - xOut));

    return (b*coshv*phase).Re();
}


inline void computeConormalBathF
(
    scalarField& F,
    const vectorField& centres,
    const conormalProfile& p
)
{
    F.setSize(centres.size());

    forAll(centres, cellI)
    {
        F[cellI] = conormalF(p, centres[cellI]);
    }
}


inline void computeConormalBathV
(
    scalarField& V,
    const vectorField& centres,
    const scalar t,
    const conormalProfile& p
)
{
    scalarField F;
    computeConormalBathF(F, centres, p);

    V.setSize(centres.size());
    const scalar sqrt1t = Foam::sqrt(1.0 + t);

    forAll(centres, cellI)
    {
        V[cellI] = sqrt1t*F[cellI];
    }
}


inline void computeConormalBathU
(
    scalarField& u1,
    scalarField& u2,
    scalarField& u3,
    const vectorField& centres,
    const scalar t,
    const conormalProfile& p
)
{
    scalarField F;
    computeConormalBathF(F, centres, p);

    u1.setSize(centres.size());
    u2.setSize(centres.size());
    u3.setSize(centres.size());

    const scalar sqrt1t = Foam::sqrt(1.0 + t);
    const scalar inv1t = 1.0/(1.0 + t);

    forAll(centres, cellI)
    {
        const vector& c = centres[cellI];
        const scalar G = 1.0 + c.x()*Foam::sqr(c.y());

        u1[cellI] = (1.0 + t)*G + sqrt1t*F[cellI];
        u2[cellI] = inv1t*Foam::pow(G, -0.5);
        u3[cellI] = 0.0;
    }
}


inline void computeConormalBathPhiEHeart
(
    scalarField& phiE,
    const vectorField& centres,
    const scalar t,
    const scalar k,
    const conormalProfile& p,
    const conormalBathCoefficients& coeffs
)
{
    scalarField V;
    computeConormalBathV(V, centres, t, p);

    phiE.setSize(centres.size());
    const scalar sqrt1t = Foam::sqrt(1.0 + t);

    forAll(centres, cellI)
    {
        phiE[cellI] =
           -k*V[cellI] + sqrt1t*conormalBathH(coeffs, centres[cellI]);
    }
}


inline void computeConormalBathPhiIHeart
(
    scalarField& phiI,
    const vectorField& centres,
    const scalar t,
    const scalar k,
    const conormalProfile& p,
    const conormalBathCoefficients& coeffs
)
{
    scalarField V, phiE;
    computeConormalBathV(V, centres, t, p);
    computeConormalBathPhiEHeart(phiE, centres, t, k, p, coeffs);

    phiI.setSize(centres.size());

    forAll(centres, cellI)
    {
        phiI[cellI] = V[cellI] + phiE[cellI];
    }
}


inline void computeConormalBathPhiEWhole
(
    scalarField& phiE,
    const vectorField& centres,
    const scalar t,
    const scalar k,
    const conormalProfile& p,
    const conormalBathCoefficients& coeffs
)
{
    phiE.setSize(centres.size());
    const scalar sqrt1t = Foam::sqrt(1.0 + t);

    forAll(centres, cellI)
    {
        const vector& c = centres[cellI];

        if (c.x() <= 0.0)
        {
            phiE[cellI] =
                sqrt1t*conormalBathBathPotential(coeffs.bL, -1.0, c);
        }
        else if (c.x() >= 1.0)
        {
            phiE[cellI] =
                sqrt1t*conormalBathBathPotential(coeffs.bR, 2.0, c);
        }
        else
        {
            phiE[cellI] =
               -k*sqrt1t*conormalF(p, c)
              + sqrt1t*conormalBathH(coeffs, c);
        }
    }
}

} // End namespace Foam

#endif

// ************************************************************************* //
