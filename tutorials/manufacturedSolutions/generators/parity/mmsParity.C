/*---------------------------------------------------------------------------*\
Application
    mmsParity

Description
    Evaluates the manufactured-solution reference functions of
    src/verificationModels at points read from stdin, so that generated NumPy
    references can be compared with the C++ the verifiers use.

Usage
    mmsParity <case> <d> [key=value ...]

    stdin   one point per line: x y z t
    stdout  one line per point, space-separated, 17 significant digits

    For pseudoECG, pseudoECGRotated and eikonalECG the point x y z is the
    electrode position and t the time.

    case              d        columns
    fda               1,2,3    V u1 u2 u3 Iion
    rotated           3        V u1 u2 u3 pdeSource
    pseudoECG         1,2,3    spatialReference reference(t)
    pseudoECGRotated  3        spatialReference reference(t)
    eikonal           1,2,3    tau gradTau_x gradTau_y gradTau_z sourceTerm
    eikonalECG        1,2,3    reference(t)
    graph             1,2,3    V u1 u2 u3
    coupled           3        V1D_0 Vm3D_0 I_0 V1D_1 Vm3D_1 I_1

    keys
        beta=<v>                    Beta of the FDA ionic closure
        order=<n>                   Gauss-Legendre order of the ECG cases
        conductivity=rotated|axis   tensor of eikonal and eikonalECG
        c0=<v>                      eikonal speed constant
        template=singleEvent|sine   eikonalECG template
        duration=<v> rise=<v>       singleEvent template parameters

    Case parameters are those of the tutorials under
    tutorials/manufacturedSolutions; each constant cites its source.
\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "List.H"
#include "tensor.H"
#include "vector.H"

#include "monodomainVerification/manufacturedFDAReference.H"
#include "monodomainVerification/manufacturedAnisotropicMonodomainReference.H"
#include "eikonalVerification/manufacturedEikonalReference.H"

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <map>
#include <sstream>
#include <string>

using namespace Foam;

namespace
{

// chi, Cm of every case: tutorials/manufacturedSolutions/monodomainPseudoECG/
// constant/electroProperties:25-26 (also eikonalECG and monodomain1D3D)
const scalar chi = 3.0;
const scalar Cm = 2.0;

// FDA diagonal conductivity:
// tutorials/manufacturedSolutions/monodomainPseudoECG/README.md:27 and
// monodomainPseudoECG/Allrun:178
tensor fdaConductivity()
{
    return tensor(0.111453302, 0, 0, 0, 0.121585420, 0, 0, 0, 0.030396355);
}

// Rotated conductivity, symmetric tensor (xx xy xz yy yz zz):
// tutorials/manufacturedSolutions/monodomainPseudoECG/constant/
// electroProperties:23 and eikonalECG/constant/electroProperties
tensor rotatedConductivity()
{
    const scalar xx = 0.106875, xy = -0.0084931, xz = -0.022561;
    const scalar yy = 0.116682, yz = -0.0130256, zz = 0.0398782;
    return tensor(xx, xy, xz, xy, yy, yz, xz, yz, zz);
}

// Axis conductivity of eikonalECG (-conductivity axis):
// tutorials/manufacturedSolutions/eikonalECG/README.md (Problem)
tensor axisConductivity()
{
    return tensor(0.111453302, 0, 0, 0, 0.121585420, 0, 0, 0, 0.030396355);
}

// Beta of monodomainFDAManufactured per dimension:
// src/ionicModels/verificationModels/monodomainFDAManufactured/
// monodomainFDAManufactured_2014.H:55-67
scalar defaultBeta(const label d)
{
    if (d == 1) return -1.1;
    if (d == 2) return -5.9;
    return -8.6;
}

// Terminals of the coupled case: monodomain1D3D/constant/purkinjeGraph.nodes003
// pvjLocations
const scalar terminalY = 0.166666666667;
const scalar terminalZ = 0.333333333333;
const scalar terminalX[2] = {0.0, 1.0};

// Junction resistance of the coupled case: monodomain1D3D/constant/
// electroProperties (rPvj 1.0)
const scalar rPvj = 1.0;


// Graph/tissue exact voltage at one point, the call the verifiers make:
// manufacturedGraphVerifier.C:preProcess, coupled1D3DMonodomainVerifier.C:
// computeExactTerminalVm (d = 1) and computeExactVmField (d = 3)
scalar exactV
(
    const scalar x,
    const scalar y,
    const scalar z,
    const scalar t,
    const label d
)
{
    scalarField X(1, x), Y(1, y), Z(1, z), V;
    computeManufacturedV(V, X, Y, Z, t, d);
    return V[0];
}


// Junction current at terminal k from the exact fields at the terminal point:
// (V1D - Vm3D)/R_pvj, coupled1D3DMonodomainVerifier.C:computeExactTerminalCurrent
// (the verifier replaces Vm3D by its volume-weighted mesh map M_k[Vm]; the
// point value is the continuum limit of that map)
scalar terminalCurrent(const scalar V1D, const scalar Vm3D)
{
    return (V1D - Vm3D)/rPvj;
}


struct Options
{
    std::map<std::string, std::string> kv;

    bool has(const std::string& key) const { return kv.count(key) > 0; }

    scalar number(const std::string& key, const scalar fallback) const
    {
        const auto it = kv.find(key);
        return it == kv.end() ? fallback : std::atof(it->second.c_str());
    }

    std::string word(const std::string& key, const std::string& fallback) const
    {
        const auto it = kv.find(key);
        return it == kv.end() ? fallback : it->second;
    }
};


void fail(const std::string& message)
{
    std::cerr << "mmsParity: " << message << std::endl;
    std::exit(1);
}


void emit(std::initializer_list<scalar> values)
{
    bool first = true;
    for (const scalar v : values)
    {
        std::printf(first ? "%.17g" : " %.17g", v);
        first = false;
    }
    std::printf("\n");
}


void fdaRow
(
    const scalar x, const scalar y, const scalar z, const scalar t,
    const label d, const Options& opt
)
{
    scalarField X(1, x), Y(1, y), Z(1, z), V, u1, u2, u3, Iion;
    computeManufacturedV(V, X, Y, Z, t, d);
    computeManufacturedU(u1, u2, u3, X, Y, Z, t, d);
    computeIion
    (
        Iion, u1, u2, u3, V, Cm, opt.number("beta", defaultBeta(d)), chi
    );
    emit({V[0], u1[0], u2[0], u3[0], Iion[0]});
}


void rotatedRow
(
    const scalar x, const scalar y, const scalar z, const scalar t,
    const Options& opt
)
{
    vectorField centres(1, vector(x, y, z));
    tensorField K(1, rotatedConductivity());
    scalarField V, u1, u2, u3, source;
    computeAnisotropicManufacturedV(V, centres, t);
    computeAnisotropicManufacturedU(u1, u2, u3, centres, t);
    computeAnisotropicManufacturedPdeSource
    (
        source, centres, K, t, opt.number("beta", defaultBeta(3))
    );
    emit({V[0], u1[0], u2[0], u3[0], source[0]});
}


void pseudoECGRow
(
    const vector& electrode, const scalar t, const label d,
    const bool rotated, const Options& opt
)
{
    const label order = label(opt.number("order", 96));
    List<scalar> nodes, weights;
    pseudoECGManufacturedQuadratureRule(order, nodes, weights);

    const scalar spatial =
        rotated
      ? computeManufacturedPseudoECGSpatialReferenceAnisotropic
        (
            electrode, rotatedConductivity(), nodes, weights
        )
      : computeManufacturedPseudoECGSpatialReference
        (
            electrode, fdaConductivity(), d, nodes, weights
        );

    emit({spatial, pseudoECGManufacturedTimeFactor(t)*spatial});
}


tensor eikonalConductivity(const Options& opt)
{
    const std::string name = opt.word("conductivity", "rotated");
    if (name == "rotated") return rotatedConductivity();
    if (name == "axis") return axisConductivity();
    fail("conductivity must be rotated or axis, got '" + name + "'");
    return tensor::zero;
}


void eikonalRow
(
    const scalar x, const scalar y, const scalar z, const label d,
    const Options& opt
)
{
    const vector p(x, y, z);
    const vector k = manufacturedEikonalK(d);
    const vector grad = manufacturedEikonalGradTau(p, k);
    emit
    ({
        manufacturedEikonalTau(p, k),
        grad.x(), grad.y(), grad.z(),
        manufacturedEikonalSourceTerm
        (
            p, eikonalConductivity(opt), chi, Cm, opt.number("c0", 1.0), k
        )
    });
}


void eikonalECGRow
(
    const vector& electrode, const scalar t, const label d, const Options& opt
)
{
    const std::string tmpl = opt.word("template", "singleEvent");
    if (tmpl != "singleEvent" && tmpl != "sine")
    {
        fail("template must be singleEvent or sine, got '" + tmpl + "'");
    }

    List<scalar> nodes, weights;
    pseudoECGManufacturedQuadratureRule
    (
        label(opt.number("order", 48)), nodes, weights
    );

    emit
    ({
        computeManufacturedEikonalECGReference
        (
            t, electrode, eikonalConductivity(opt), manufacturedEikonalK(d), d,
            tmpl == "sine", opt.number("duration", 64), opt.number("rise", 4),
            nodes, weights
        )
    });
}


void graphRow
(
    const scalar x, const scalar y, const scalar z, const scalar t,
    const label d
)
{
    scalarField X(1, x), Y(1, y), Z(1, z), V, u1, u2, u3;
    computeManufacturedV(V, X, Y, Z, t, d);
    computeManufacturedU(u1, u2, u3, X, Y, Z, t, d);
    emit({V[0], u1[0], u2[0], u3[0]});
}


void coupledRow(const scalar t)
{
    const scalar v1D0 = exactV(terminalX[0], 0, 0, t, 1);
    const scalar v3D0 = exactV(terminalX[0], terminalY, terminalZ, t, 3);
    const scalar v1D1 = exactV(terminalX[1], 0, 0, t, 1);
    const scalar v3D1 = exactV(terminalX[1], terminalY, terminalZ, t, 3);

    emit
    ({
        v1D0, v3D0, terminalCurrent(v1D0, v3D0),
        v1D1, v3D1, terminalCurrent(v1D1, v3D1)
    });
}

} // End anonymous namespace


int main(int argc, char* argv[])
{
    if (argc < 3)
    {
        fail("usage: mmsParity <case> <d> [key=value ...]");
    }

    const std::string name(argv[1]);
    const label d = std::atoi(argv[2]);

    if (d < 1 || d > 3)
    {
        fail("d must be 1, 2 or 3");
    }

    Options opt;
    for (int i = 3; i < argc; ++i)
    {
        const std::string arg(argv[i]);
        const auto eq = arg.find('=');
        if (eq == std::string::npos)
        {
            fail("expected key=value, got '" + arg + "'");
        }
        opt.kv[arg.substr(0, eq)] = arg.substr(eq + 1);
    }

    if ((name == "rotated" || name == "pseudoECGRotated" || name == "coupled")
     && d != 3)
    {
        fail(name + " exists for d = 3 only");
    }

    std::string line;
    while (std::getline(std::cin, line))
    {
        if (line.find_first_not_of(" \t\r") == std::string::npos)
        {
            continue;
        }

        std::istringstream in(line);
        scalar x, y, z, t;
        if (!(in >> x >> y >> z >> t))
        {
            fail("expected 'x y z t', got '" + line + "'");
        }

        if (name == "fda") fdaRow(x, y, z, t, d, opt);
        else if (name == "rotated") rotatedRow(x, y, z, t, opt);
        else if (name == "pseudoECG")
        {
            pseudoECGRow(vector(x, y, z), t, d, false, opt);
        }
        else if (name == "pseudoECGRotated")
        {
            pseudoECGRow(vector(x, y, z), t, d, true, opt);
        }
        else if (name == "eikonal") eikonalRow(x, y, z, d, opt);
        else if (name == "eikonalECG")
        {
            eikonalECGRow(vector(x, y, z), t, d, opt);
        }
        else if (name == "graph") graphRow(x, y, z, t, d);
        else if (name == "coupled") coupledRow(t);
        else fail("unknown case '" + name + "'");
    }

    return 0;
}


// ************************************************************************* //
