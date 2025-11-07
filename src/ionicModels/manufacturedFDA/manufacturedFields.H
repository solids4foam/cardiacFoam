#ifndef MANUFACTUREDFIELDS_H
#define MANUFACTUREDFIELDS_H

#include "scalarField.H"
#include "mathematicalConstants.H"

namespace Foam
{

//- Compute F(x) = cos(pi * x)
inline void computeF(scalarField& F, const scalarField& x)
{
    F.setSize(x.size());
    forAll(x, i)
    {
        F[i] = Foam::cos(constant::mathematical::pi * x[i]);
    }
}

//- Compute G(x) = 1 + x
inline void computeG(scalarField& G, const scalarField& x)
{
    G.setSize(x.size());
    forAll(x, i)
    {
        G[i] = 1.0 + x[i];
    }
}

//- Compute manufactured voltage field: V = sqrt(1 + t) * F(x)
inline void computeManufacturedV
(
    scalarField& V,
    const scalarField& x,
    const scalar t
)
{
    scalarField Fx;
    computeF(Fx, x);

    V.setSize(x.size());
    const scalar coeff = Foam::sqrt(1.0 + t);
    forAll(x, i)
    {
        V[i] = coeff*Fx[i];
    }
}

//- Compute manufactured u-fields:
//  u1 = (1 + t)*G(x) + sqrt(1 + t)*F(x)
//  u2 = (1 + t)^(-1)*G(x)^(-1/2)
//  u3 = 0
inline void computeManufacturedU
(
    scalarField& u1,
    scalarField& u2,
    scalarField& u3,
    const scalarField& x,
    const scalar t
)
{
    scalarField Fx, Gx;
    computeF(Fx, x);
    computeG(Gx, x);

    u1.setSize(x.size());
    u2.setSize(x.size());
    u3.setSize(x.size());

    const scalar sqrt1t = Foam::sqrt(1.0 + t);
    const scalar inv1t  = 1.0/(1.0 + t);

    forAll(x, i)
    {
        u1[i] = (1.0 + t)*Gx[i] + sqrt1t*Fx[i];
        u2[i] = inv1t * Foam::pow(Gx[i], -0.5);
        u3[i] = 0.0;
    }
}

inline void computeIion
(
    scalarField& Iion,
    const scalarField& u1,
    const scalarField& u2,
    const scalarField& u3,
    const scalarField& V,
    const scalar Cm,
    const scalar beta,
    const scalar chi
)
{
    const label n = V.size();
    Iion.setSize(n);

    forAll(V, i)
    {
        const scalar term1 = -(Cm / 2.0) * (u1[i] + u3[i] - V[i]) * Foam::pow(u2[i], 2) * (V[i] - u3[i]);
        const scalar term2 = (beta / chi) * (V[i] - u3[i]);
        Iion[i] = term1 + term2;
    }
}

inline void computeIion_point
(
    scalar& Iion,
    const scalar u1,
    const scalar u2,
    const scalar u3,
    const scalar V,
    const scalar Cm,
    const scalar beta,
    const scalar chi
)
{
    const scalar term1 = -(Cm / 2.0) * (u1 + u3 - V) * Foam::pow(u2, 2) * (V - u3);
    const scalar term2 = (beta / chi) * (V - u3);
    Iion = term1 + term2;
}

inline void computeAndPrintErrors
(
    const scalarField& Vm,
    const scalarField& u1m,
    const scalarField& u2m,
    const scalarField& u3m,
    const scalarField& x,
    const scalar t
)
{
    // Manufactured (exact) fields
    scalarField Vex, u1ex, u2ex, u3ex;
    computeManufacturedV(Vex, x, t);
    computeManufacturedU(u1ex, u2ex, u3ex, x, t);

    auto computeNorms = [&](const scalarField& num, const scalarField& exact)
    {
        const label n = num.size();
        scalar sumAbs = 0.0, sumSq = 0.0, maxAbs = 0.0;
        forAll(num, i)
        {
            scalar diff = Foam::mag(num[i] - exact[i]);
            sumAbs += diff;
            sumSq  += diff*diff;
            if (diff > maxAbs) maxAbs = diff;
        }
        scalar L1 = sumAbs / n;
        scalar L2 = Foam::sqrt(sumSq / n);
        scalar Linf = maxAbs;
        return std::tuple<scalar, scalar, scalar>(L1, L2, Linf);
    };

    auto [L1_V, L2_V, Linf_V] = computeNorms(Vm, Vex);
    auto [L1_u1, L2_u1, Linf_u1] = computeNorms(u1m, u1ex);
    auto [L1_u2, L2_u2, Linf_u2] = computeNorms(u2m, u2ex);
    auto [L1_u3, L2_u3, Linf_u3] = computeNorms(u3m, u3ex);

    Info<< "\nManufactured-solution error summary (t = " << t << "):" << nl
        << "-------------------------------------------------" << nl
        << "Field     L1-error       L2-error       Linf-error" << nl
        << "Vm     " << L1_V   << "   " << L2_V   << "   " << Linf_V   << nl
        << "u1     " << L1_u1 << "   " << L2_u1 << "   " << Linf_u1 << nl
        << "u2     " << L1_u2 << "   " << L2_u2 << "   " << Linf_u2 << nl
        << "u3     " << L1_u3 << "   " << L2_u3 << "   " << Linf_u3 << nl
        << "-------------------------------------------------" << nl << endl;
}

//==============================================================================
// RHS of u-ODEs: f(u,V)
//==============================================================================

inline void compute_f_RHS
(
    scalarField& f1,
    scalarField& f2,
    scalarField& f3,
    const scalarField& u1,
    const scalarField& u2,
    const scalarField& u3,
    const scalarField& V
)
{
    const label n = V.size();
    f1.setSize(n);
    f2.setSize(n);
    f3.setSize(n);

    forAll(V, i)
    {
        const scalar d = (u1[i] + u3[i] - V[i]);     // (u1 + u3 - V)
        const scalar Vm_u3 = (V[i] - u3[i]);         // (V - u3)
        const scalar u2_2 = Foam::pow(u2[i], 2);
        const scalar u2_3 = u2_2 * u2[i];

        f1[i] = d*d*u2[i] + 0.5*d*u2_2*Vm_u3;
        f2[i] = -d*u2_3;
        f3[i] = 0.0;
    }
}

//- Compute f(u, V) for one integration point
inline void compute_f_RHS_point
(
    scalar& f1,
    scalar& f2,
    scalar& f3,
    const scalar u1,
    const scalar u2,
    const scalar u3,
    const scalar V
)
{
    const scalar d = (u1 + u3 - V);     // (u1 + u3 - V)
    const scalar Vm_u3 = (V - u3);      // (V - u3)
    const scalar u2_2 = Foam::pow(u2, 2);
    const scalar u2_3 = u2_2 * u2;

    f1 = d*d*u2 + 0.5*d*u2_2*Vm_u3;
    f2 = -d*u2_3;
    f3 = 0.0;
}


} // namespace Foam

#endif
