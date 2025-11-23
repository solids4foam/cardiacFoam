#ifndef MANUFACTUREDFIELDS_H
#define MANUFACTUREDFIELDS_H

#include "scalarField.H"
#include "mathematicalConstants.H"
#include <sstream>

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



inline void computeAndPrintErrors
(
    const scalarField& Vm,
    const scalarField& u1m,
    const scalarField& u2m,
    const scalarField& x,
    const scalarField& y,
    const scalarField& z,
    const scalar t,
    const scalar dimension,
    const int N,
    const scalar dx,
    const scalar dt,
    const int nsteps,
    const bool solveExplicit //,
    //const std::string& filename = ""
)
{
    // Pout<< __FILE__<< __LINE__ << endl;
    // Manufactured (exact) fields
    scalarField Vex, u1ex, u2ex, u3ex;
    computeManufacturedV(Vex, x, y, z, t, dimension);
    computeManufacturedU(u1ex, u2ex, u3ex, x, y, z, t, dimension);

    auto computeNorms = [&](const scalarField& num, const scalarField& exact)
    {
        scalar sumAbs = 0.0, sumSq = 0.0, maxAbs = 0.0;
        forAll(num, i)
        {
            scalar diff = Foam::mag(num[i] - exact[i]);
            sumAbs += diff;
            sumSq  += diff*diff;
            if (diff > maxAbs) maxAbs = diff;
        }
        // Sync sumAbs
        reduce(maxAbs, maxOp<scalar>());
        reduce(sumAbs, sumOp<scalar>());
        reduce(sumSq, sumOp<scalar>());
        label n = num.size();
        reduce(n, sumOp<int>());

        scalar L1 = sumAbs / n;
        scalar L2 = Foam::sqrt(sumSq / n);
        scalar Linf = maxAbs;
        return std::tuple<scalar, scalar, scalar>(L1, L2, Linf);
    };

    auto [L1_V, L2_V, Linf_V] = computeNorms(Vm, Vex);
    auto [L1_u1, L2_u1, Linf_u1] = computeNorms(u1m, u1ex);
    auto [L1_u2, L2_u2, Linf_u2] = computeNorms(u2m, u2ex);

    // --- Print to console ---
    Info << "\nSimulation summary:\n";
    Info << "-------------------\n";
    Info << "Number of cells (N)   = " << N << nl;
    Info << "Solver type           = " << (solveExplicit ? "Explicit" : "Implicit") << nl;
    Info << "Grid spacing (dx)     = " << dx << nl;
    Info << "Time step (dt)        = " << dt << nl;
    Info << "Number of steps       = " << nsteps << nl;
    Info << "Final simulation time = " << t << nl;
    Info << "-------------------\n";
    Info << "\nManufactured-solution error summary (t = " << t << "):" << nl
         << "-------------------------------------------------\n"
         << "Field     L1-error       L2-error       Linf-error\n"
         << "Vm     " << L1_V   << "   " << L2_V   << "   " << Linf_V << nl
         << "u1     " << L1_u1 << "   " << L2_u1 << "   " << Linf_u1 << nl
         << "u2     " << L1_u2 << "   " << L2_u2 << "   " << Linf_u2 << nl
         << "-------------------------------------------------\n" << endl;


        

    // std::string outFilename = filename;
    // if (outFilename.empty())
    if (Pstream::master())
    {
        // Pstream:nProcs()
        // Pstream:myProcNum()
        // Pstream::parRun()

        // std::stringstream filenameSS;
        // filenameSS << N << "_cells_" << (solveExplicit ? "explicit" : "implicit") << ".dat";
        // outFilename = filenameSS.str();
        fileName outFilename
        (
            Foam::name(N) + "_cells_"
          + word(solveExplicit ? "explicit" : "implicit") + ".dat"
        );

        //        std::ofstream fout(outFilename, std::ios::out | std::ios::app);
        OFstream fout(outFilename);
        // if (fout.is_open())
        {
            fout << "Manufactured-solution error summary (t = " << t << "):\n";
            fout << "Field     L1-error       L2-error       Linf-error\n";
            fout << "Vm     " << L1_V   << "   " << L2_V   << "   " << Linf_V << "\n";
            fout << "u1     " << L1_u1 << "   " << L2_u1 << "   " << Linf_u1 << "\n";
            fout << "u2     " << L1_u2 << "   " << L2_u2 << "   " << Linf_u2 << "\n";
            fout << "-------------------------------------------------\n\n";

            fout << "\nSimulation summary:\n";
            fout << "-------------------\n";
            fout << "Number of cells (N)   = " << N << "\n";
            fout << "Solver type           = " << (solveExplicit ? "Explicit" : "Implicit") << "\n";
            fout << "Grid spacing (dx)     = " << dx << "\n";
            fout << "Time step (dt)        = " << dt << "\n";
            fout << "Number of steps       = " << nsteps << "\n";
            fout << "Final simulation time = " << t << "\n";
            fout << "-------------------\n\n";

            // fout.close();
        }
        // else
        // {
        //     Info << "Warning: Could not open file " << filename << " for writing.\n";
        // }
    }
}

} // namespace Foam

#endif

