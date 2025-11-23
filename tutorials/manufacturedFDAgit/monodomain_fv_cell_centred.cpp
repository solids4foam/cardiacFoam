#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <tuple>
#include <algorithm>
#include <iostream>
#include <iomanip>

using std::vector;
using std::size_t;

// -----------------------------
// Parameters
// -----------------------------
static constexpr double chi    = 3.0;
static constexpr double Cm     = 2.0;
static constexpr double beta   = -1.1;                  // so beta/chi = -1.1/3
static constexpr double sigma  = 1.1 / (M_PI*M_PI);     // 1.1 / pi^2
static constexpr double Tfinal = 0.2;

// -----------------------------
// Manufactured fields (Option 2):
// cell-centred FV with Neumann BCs
// -----------------------------
inline double F(double x) { return std::cos(M_PI * (x)); }
inline double G(double x)            { return 1.0 + x; }

inline double manufactured_V(double t, double x)
{
    return std::sqrt(1.0 + t) * F(x);
}

inline void manufactured_u(double t, double x,
                           double& u1, double& u2, double& u3)
{
    u1 = (1.0 + t)*G(x) + std::sqrt(1.0 + t)*F(x);
    u2 = std::pow(1.0 + t, -1.0) * std::pow(G(x), -0.5);
    u3 = 0.0;
}

// -----------------------------
// Ionic current and (CORRECT) RHS of u-ODEs
// (general form reduces to u3=0 here)
// -----------------------------
inline double I_ion(double u1, double u2, double u3, double V)
{
    // I_ion = -(Cm/2)*(u1 + u3 - V)*u2^2*(V - u3) + (beta/chi)*(V - u3)
    const double d      = (u1 + u3 - V);
    const double Vm_u3  = (V - u3);
    const double u2_2   = u2*u2;
    return -(Cm/2.0)*d*u2_2*Vm_u3 + (beta/chi)*Vm_u3;
}

inline void f_rhs(double u1, double u2, double u3, double V,
                  double& f1, double& f2, double& f3)
{
    // CORRECT manufactured ODEs:
    // f1 = (u1 - 0.5*V)*(u1 - V)*u2^2
    // f2 = -(u1 - V)*u2^3
    // f3 = 0
    const double d    = (u1 + u3 - V);  // (u1 - V) when u3=0
    const double u2_2 = u2*u2;
    const double u2_3 = u2_2*u2;

    f1 = (u1 - 0.5*V) * d * u2_2;
    f2 = - d * u2_3;
    f3 = 0.0;
}

// -----------------------------
// Utility: norms
// -----------------------------
struct Norms { double L1, L2, Linf; };

inline Norms error_norms(const vector<double>& a, const vector<double>& b)
{
    const size_t n = a.size();
    double s1 = 0.0, s2 = 0.0, sInf = 0.0;
    for (size_t i=0; i<n; ++i)
    {
        const double d = std::abs(a[i] - b[i]);
        s1  += d;
        s2  += d*d;
        sInf = std::max(sInf, d);
    }
    return { s1 / double(n), std::sqrt(s2 / double(n)), sInf };
}

// -----------------------------
// Laplacian with zero-flux BC (mirrored ghost)
// -----------------------------
inline void laplacian_neumann(const vector<double>& V, double dx, vector<double>& L)
{
    const size_t N = V.size();
    L.assign(N, 0.0);
    if (N <= 1) { if (N==1) L[0]=0.0; return; }

    for (size_t i=1; i+1<N; ++i)
        L[i] = (V[i+1] - 2.0*V[i] + V[i-1]) / (dx*dx);

    // boundaries (mirrored neighbours)
    // L[0]   = (V[1]   - 2.0*V[0]   + V[1])   / (dx*dx); // V[-1]=V[1]
    // L[N-1] = (V[N-2] - 2.0*V[N-1] + V[N-2]) / (dx*dx); // V[N]=V[N-2]
    L[0]   = (V[1]   - V[0])   / (dx*dx);
    L[N-1] = (V[N-2] - V[N-1]) / (dx*dx);
}

// -----------------------------
// One run at resolution N (cell-centred FV grid)
// -----------------------------
struct RunResult {
    int    N;
    double dx;
    double dt;
    int    steps;
    double errV_inf;
    double errU1_inf;
    double errU2_inf;
};

RunResult run_solver_ccfv(int N, double cfl = 0.1)
{
    // Cell-centred grid on [0,1]: centres x_i = (i+0.5)*dx
    const double dx = 1.0 / double(N);
    vector<double> x(N);
    for (int i=0; i<N; ++i) x[i] = (i + 0.5)*dx;

    // Explicit diffusion stability
    double dt = cfl * dx*dx / sigma;
    int nsteps = int(std::ceil(Tfinal / dt));
    dt = Tfinal / double(nsteps); // snap to hit Tfinal

    // Initial conditions
    vector<double> V(N), u1(N), u2(N), u3(N, 0.0);
    for (int i=0; i<N; ++i)
    {
        V[i] = manufactured_V(0.0, x[i]);
        double a,b,c;
        manufactured_u(0.0, x[i], a, b, c);
        u1[i] = a; u2[i] = b; u3[i] = c;
    }

    // Time loop: forward Euler
    vector<double> L(N);
    for (int step=0; step<nsteps; ++step)
    {
        // Laplacian(V)
        laplacian_neumann(V, dx, L);

        // dV/dt = (sigma*L - chi*Iion)/(chi*Cm)
        vector<double> dVdt(N), Iion(N);
        for (int i=0; i<N; ++i)
        {
            Iion[i] = I_ion(u1[i], u2[i], u3[i], V[i]);
            dVdt[i] = (sigma*L[i] - chi*Iion[i]) / (chi*Cm);
        }

        // Update V explicitly
        for (int i=0; i<N; ++i) V[i] += dt * dVdt[i];

        // Update u explicitly using the (new) V (matches the Python ordering)
        for (int i=0; i<N; ++i)
        {
            double rhs1, rhs2, rhs3;
            f_rhs(u1[i], u2[i], u3[i], V[i], rhs1, rhs2, rhs3);
            u1[i] += dt * rhs1;
            u2[i] += dt * rhs2;
            // u3 stays 0
        }
    }

    // Errors vs manufactured solution at Tfinal
    vector<double> Vex(N), u1e(N), u2e(N), u3e(N);
    for (int i=0; i<N; ++i)
    {
        Vex[i] = manufactured_V(Tfinal, x[i]);
        double a,b,c;
        manufactured_u(Tfinal, x[i], a, b, c);
        u1e[i]=a; u2e[i]=b; u3e[i]=c;
    }
    auto eV  = error_norms(V,   Vex);
    auto eU1 = error_norms(u1,  u1e);
    auto eU2 = error_norms(u2,  u2e);

    RunResult r;
    r.N = N; r.dx = dx; r.dt = dt; r.steps = nsteps;
    r.errV_inf  = eV.Linf;
    r.errU1_inf = eU1.Linf;
    r.errU2_inf = eU2.Linf;
    return r;
}

// -----------------------------
// Observed orders
// -----------------------------
inline double observed_order(double eCoarse, double eFine, double dxCoarse, double dxFine)
{
    return std::log(eCoarse/eFine) / std::log(dxCoarse/dxFine);
}

// -----------------------------
// main
// -----------------------------
int main()
{
    std::vector<int> Ns = {50, 100, 200, 400};

    std::vector<RunResult> results;
    results.reserve(Ns.size());
    for (int N : Ns) results.push_back(run_solver_ccfv(N));

    std::cout.setf(std::ios::fixed);
    std::cout << "   N      dx        dt       steps      errV_inf      pV"
              << "      errU1_inf     pU1      errU2_inf     pU2\n";

    for (size_t i=0; i<results.size(); ++i)
    {
        const auto& r = results[i];
        double pV  = NAN, pU1 = NAN, pU2 = NAN;
        if (i>0)
        {
            const auto& rc = results[i-1];
            pV  = observed_order(rc.errV_inf,  r.errV_inf,  rc.dx, r.dx);
            pU1 = observed_order(rc.errU1_inf, r.errU1_inf, rc.dx, r.dx);
            pU2 = observed_order(rc.errU2_inf, r.errU2_inf, rc.dx, r.dx);
        }

        std::cout << std::setw(4) << r.N
                  << std::setw(9) << std::setprecision(6) << r.dx
                  << std::setw(11) << std::setprecision(6) << r.dt
                  << std::setw(10) << r.steps
                  << std::setw(14) << std::setprecision(6) << r.errV_inf
                  << std::setw(9)  << std::setprecision(3) << (std::isnan(pV)?0.0:pV)
                  << std::setw(14) << std::setprecision(6) << r.errU1_inf
                  << std::setw(9)  << std::setprecision(3) << (std::isnan(pU1)?0.0:pU1)
                  << std::setw(14) << std::setprecision(6) << r.errU2_inf
                  << std::setw(9)  << std::setprecision(3) << (std::isnan(pU2)?0.0:pU2)
                  << "\n";
    }

    std::cout << "\nNotes:\n"
              << " - Cell-centred FV grid: x_i = (i+0.5)*dx on [0,1].\n"
              << " - Explicit Euler in time with dt ~ dx^2 for diffusion stability.\n";
    return 0;
}
