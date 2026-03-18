#ifndef TMANUFACTUREDEXACTFIELDS_H
#define TMANUFACTUREDEXACTFIELDS_H

#include "scalarField.H"
#include "mathematicalConstants.H"

namespace Foam
{

inline void computeF
(
    scalarField& F,
    const scalarField& x,
    const scalarField& y,
    const scalarField& z,
    const int dimension

)
{
    F.setSize(x.size());
    forAll(x, i)
    {
        if (dimension == 1)
        {
            F[i] = Foam::cos(constant::mathematical::pi * x[i]);
        }
        else if (dimension == 2)
        {
            F[i] = Foam::cos(constant::mathematical::pi * x[i])
                 * Foam::cos(2.0 * constant::mathematical::pi * y[i]);
        }
        else if (dimension == 3)
        {
            F[i] = Foam::cos(constant::mathematical::pi * x[i])
                 * Foam::cos(2.0 * constant::mathematical::pi * y[i])
                 * Foam::cos(3.0 * constant::mathematical::pi * z[i]);
        }
    }
}

inline void computeG
(
    scalarField& G,
    const scalarField& x,
    const scalarField& y,
    const scalarField& z,
    const int dimension
)
{
    G.setSize(x.size());
    forAll(x, i)
    {
        if (dimension == 1)
            G[i] = 1.0 + x[i];
        else if (dimension == 2)
            G[i] = 1.0 + x[i]*Foam::pow(y[i], 2);
        else
            G[i] = 1.0 + x[i]*Foam::pow(y[i], 2)*Foam::pow(z[i], 3);
    }
}
//- Compute manufactured voltage field: V = sqrt(1 + t) * F(x)
inline void computeManufacturedV
(
    scalarField& V,
    const scalarField& x,
    const scalarField& y,
    const scalarField& z,
    const scalar t,
    const int dimension
)
{
    scalarField Fx;
    computeF(Fx, x, y, z, dimension);

    V.setSize(x.size());
    const scalar coeff = Foam::sqrt(1.0 + t);

    forAll(x, i)
    {
        V[i] = coeff * Fx[i];
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
    const scalarField& y,
    const scalarField& z,
    const scalar t,
    const int dimension
)
{
    scalarField Fx, Gx;
    computeF(Fx, x, y, z, dimension);
    computeG(Gx, x, y, z, dimension);

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

} // namespace Foam

#endif
