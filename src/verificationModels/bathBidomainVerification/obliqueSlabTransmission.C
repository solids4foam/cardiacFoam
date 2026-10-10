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

#include "obliqueSlabTransmission.H"

#include "FixedList.H"
#include "IStringStream.H"
#include "surfaceFields.H"
#include "syncTools.H"
#include "unitConversion.H"
#include "volFields.H"

#include <cmath>

namespace Foam
{

namespace
{

scalar sampleValue
(
    const transmissionMMS::SelfTestSample& s,
    const transmissionMMS::Params& p
)
{
    const point pt(s.x, s.y, s.z);
    const label c = s.component;

    switch (s.fn)
    {
        case 0: return transmissionMMS::Vm(pt, s.t, p);
        case 1: return transmissionMMS::phiE(pt, s.t, s.zone, p);
        case 2: return transmissionMMS::phiI(pt, s.t, p);
        case 3:
            return transmissionMMS::gradPhiE(pt, s.t, s.zone, p).component(c);
        case 4:
            return transmissionMMS::gradVm(pt, s.t, p).component(c);
        case 5: return transmissionMMS::u1(pt, s.t, p);
        case 6: return transmissionMMS::u2(pt, s.t, p);
        case 7: return transmissionMMS::Gi(pt, p).component(c);
        case 8: return transmissionMMS::Ge(pt, p).component(c);
        case 9: return transmissionMMS::sigmaB(pt, s.zone, p);
        case 10:
            return transmissionMMS::divGiGradPhiI(pt, s.t, p);
        case 11:
            return transmissionMMS::phiESource(pt, s.t, s.zone, p);
        case 12:
        {
            vector n(Zero);
            n.component(c) = 1;
            return transmissionMMS::interfaceFlux(pt, n, s.t, p);
        }
        case 13: return transmissionMMS::G(pt, p);
    }

    FatalErrorInFunction
        << "Self-test sample function index " << s.fn << " is not known."
        << exit(FatalError);

    return 0;
}

} // End anonymous namespace


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

vector obliqueSlabTransmission::normal() const
{
    return vector(Foam::cos(p_.alpha), Foam::sin(p_.alpha), 0);
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

transmissionMMS::Params
obliqueSlabTransmission::readParams
(
    const dictionary& dict
)
{
    transmissionMMS::Params p;

    const word solution(dict.get<word>("solution"));
    if (solution == "affine")
    {
        p.solution = 0;
    }
    else if (solution == "levelSet")
    {
        p.solution = 1;
    }
    else
    {
        FatalIOErrorInFunction(dict)
            << "Unknown solution '" << solution << "'. Valid solutions: "
            << "(affine levelSet)."
            << exit(FatalIOError);
    }

    const word geometry(dict.get<word>("geometry"));
    if (geometry != "obliqueSlab")
    {
        FatalIOErrorInFunction(dict)
            << "Unknown geometry '" << geometry << "'. Valid geometries: "
            << "(obliqueSlab)."
            << exit(FatalIOError);
    }

    const dictionary& geom = dict.subDict("geometryCoeffs");
    p.alpha = degToRad(geom.get<scalar>("alpha"));
    p.s1 = geom.get<scalar>("s1");
    p.s2 = geom.get<scalar>("s2");
    p.period = geom.get<scalar>("period");

    const dictionary& cond = dict.subDict("conductivityCoeffs");
    const vector gi(cond.get<vector>("gi"));
    const vector ge(cond.get<vector>("ge"));
    p.gil = gi.x();
    p.git = gi.y();
    p.gin = gi.z();
    p.gel = ge.x();
    p.get = ge.y();
    p.gen = ge.z();
    p.theta0 = degToRad(cond.get<scalar>("theta0"));
    p.theta1 = degToRad(cond.get<scalar>("theta1"));
    p.sigma0Low = cond.get<scalar>("sigma0Low");
    p.sigma0High = cond.get<scalar>("sigma0High");
    p.sigmaEps = cond.get<scalar>("sigmaEps");

    const dictionary& sol = dict.subDict("solutionCoeffs");
    const FixedList<scalar, 4> f(sol.get<FixedList<scalar, 4>>("f"));
    const FixedList<scalar, 3> h(sol.get<FixedList<scalar, 3>>("h"));
    forAll(f, i)
    {
        p.fA[i] = f[i];
    }
    forAll(h, i)
    {
        p.hA[i] = h[i];
    }

    p.cmAmplitude = 0;
    p.cmOmega = 0;
    if (const dictionary* cm = sol.findDict("commonMode"))
    {
        p.cmAmplitude = cm->getOrDefault<scalar>("amplitude", 0);
        p.cmOmega = cm->getOrDefault<scalar>("omega", 0);
    }

    const dictionary& affine = sol.subDict("affine");
    p.P0 = affine.get<scalar>("P0");
    p.P1 = affine.get<scalar>("P1");
    p.C = affine.get<scalar>("C");

    if (p.solution == 0 && p.cmAmplitude != 0)
    {
        FatalIOErrorInFunction(dict)
            << "The affine solution has no common mode, but "
            << "solutionCoeffs.commonMode.amplitude is " << p.cmAmplitude
            << ". Set it to 0 or use the levelSet solution."
            << exit(FatalIOError);
    }

    return p;
}


void obliqueSlabTransmission::selfTest()
{
    IStringStream defaultsStream
    (
        "solution levelSet;"
        "geometry obliqueSlab;"
        "geometryCoeffs { alpha 25; s1 0.5; s2 1.5; period 1; }"
        "conductivityCoeffs"
        "{"
        "    gi (0.17 0.019 0.019); ge (0.62 0.12 0.12);"
        "    theta0 -20; theta1 90;"
        "    sigma0Low 0.3; sigma0High 0.45; sigmaEps 0.2;"
        "}"
        "solutionCoeffs"
        "{"
        "    f (0.3 0.2 0.4 0.15); h (-0.2 0.3 0.1);"
        "    commonMode { amplitude 0; omega 0; }"
        "    affine { P0 0.1; P1 0.4; C 0.3; }"
        "}"
    );
    const dictionary defaultsDict(defaultsStream);

    const transmissionMMS::Params read = readParams(defaultsDict);
    const transmissionMMS::Params ref = transmissionMMS::defaultParams();

    wordList differing;
    auto compare = [&differing](const char* name, scalar a, scalar b)
    {
        if (a != b)
        {
            differing.append(word(name));
        }
    };

    compare("alpha", read.alpha, ref.alpha);
    compare("s1", read.s1, ref.s1);
    compare("s2", read.s2, ref.s2);
    compare("period", read.period, ref.period);
    compare("gil", read.gil, ref.gil);
    compare("git", read.git, ref.git);
    compare("gin", read.gin, ref.gin);
    compare("gel", read.gel, ref.gel);
    compare("get", read.get, ref.get);
    compare("gen", read.gen, ref.gen);
    compare("theta0", read.theta0, ref.theta0);
    compare("theta1", read.theta1, ref.theta1);
    compare("sigma0Low", read.sigma0Low, ref.sigma0Low);
    compare("sigma0High", read.sigma0High, ref.sigma0High);
    compare("sigmaEps", read.sigmaEps, ref.sigmaEps);
    for (label i = 0; i < 4; ++i)
    {
        compare("fA", read.fA[i], ref.fA[i]);
    }
    for (label i = 0; i < 3; ++i)
    {
        compare("hA", read.hA[i], ref.hA[i]);
    }
    compare("cmAmplitude", read.cmAmplitude, ref.cmAmplitude);
    compare("cmOmega", read.cmOmega, ref.cmOmega);
    compare("P0", read.P0, ref.P0);
    compare("P1", read.P1, ref.P1);
    compare("C", read.C, ref.C);
    compare("solution", read.solution, ref.solution);

    if (!differing.empty())
    {
        FatalErrorInFunction
            << "The verification-dictionary reader does not reproduce the "
            << "reference default parameters; differing fields: "
            << differing << "."
            << exit(FatalError);
    }

    const transmissionMMS::Params paramSets[2] =
    {
        transmissionMMS::defaultParams(),
        transmissionMMS::alternateParams()
    };

    for (label i = 0; i < transmissionMMS::nSelfTestSamples; ++i)
    {
        const transmissionMMS::SelfTestSample& s =
            transmissionMMS::selfTestSamples[i];

        if (s.param < 0 || s.param > 1)
        {
            FatalErrorInFunction
                << "Self-test sample " << i << " names parameter set "
                << s.param << "." << exit(FatalError);
        }

        transmissionMMS::Params pp = paramSets[s.param];
        pp.solution = s.solution;

        const scalar value = sampleValue(s, pp);
        const scalar err = mag(value - s.value);

        if (!(err <= 1e-11*max(scalar(1), mag(s.value))))
        {
            FatalErrorInFunction
                << "Transmission MMS self-test sample " << i
                << " (function " << s.fn << ", parameter set " << s.param
                << ", solution " << s.solution << ", zone " << s.zone
                << ", component " << s.component << ", point ("
                << s.x << ' ' << s.y << ' ' << s.z << "), t " << s.t
                << ") evaluates to " << value << " against the reference "
                << s.value << "."
                << exit(FatalError);
        }
    }

    Info<< "transmission MMS self-test: " << transmissionMMS::nSelfTestSamples
        << " samples OK" << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

obliqueSlabTransmission::obliqueSlabTransmission(const dictionary& dict)
:
    p_(readParams(dict))
{
    selfTest();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

word obliqueSlabTransmission::zoneName(const label zone) const
{
    static const word names[3] = {"bathLow", "myocardium", "bathHigh"};
    return names[zone];
}


label obliqueSlabTransmission::nMisplacedCells
(
    const fvMesh& mesh,
    const labelUList& zone
) const
{
    const vectorField& C = mesh.C().primitiveField();
    label nMismatched = 0;

    forAll(zone, cellI)
    {
        if (zone[cellI] >= 0 && transmissionMMS::zoneOf(C[cellI], p_) != zone[cellI])
        {
            ++nMismatched;
        }
    }

    return nMismatched;
}


void obliqueSlabTransmission::validateInterface
(
    const fvMesh& baseMesh,
    const labelUList& zone,
    const fvMesh& heartMesh
) const
{
    // Heart-bath faces lie on s = s1 and s = s2; heart cells fill s1 < s < s2 over one period

    const vector nHat(normal());
    const labelList nbrZone
    (
        syncTools::swapBoundaryCellList(baseMesh, zone)
    );

    scalar maxPlaneErr = 0;
    label nInterfaceFaces = 0;

    auto checkPlane = [&](const label bathZone, const point& Cf)
    {
        const scalar sExact =
            (bathZone == transmissionMMS::bathLow ? p_.s1 : p_.s2);
        maxPlaneErr = max(maxPlaneErr, mag((nHat & Cf) - sExact));
        ++nInterfaceFaces;
    };

    const labelUList& own = baseMesh.owner();
    const labelUList& nei = baseMesh.neighbour();

    forAll(nei, faceI)
    {
        const label zo = zone[own[faceI]];
        const label zn = zone[nei[faceI]];

        if (zo == zn)
        {
            continue;
        }
        if (zo == transmissionMMS::myocardium)
        {
            checkPlane(zn, baseMesh.Cf()[faceI]);
        }
        else if (zn == transmissionMMS::myocardium)
        {
            checkPlane(zo, baseMesh.Cf()[faceI]);
        }
    }

    forAll(baseMesh.boundary(), patchI)
    {
        const fvPatch& fp = baseMesh.boundary()[patchI];
        if (!fp.coupled())
        {
            continue;
        }

        const labelUList& faceCells = fp.faceCells();
        const label start = fp.patch().start() - baseMesh.nInternalFaces();

        forAll(faceCells, i)
        {
            const label zLocal = zone[faceCells[i]];
            const label zNbr = nbrZone[start + i];

            if (zLocal == zNbr)
            {
                continue;
            }
            if (zLocal == transmissionMMS::myocardium)
            {
                checkPlane(zNbr, baseMesh.Cf().boundaryField()[patchI][i]);
            }
            else if (zNbr == transmissionMMS::myocardium)
            {
                checkPlane(zLocal, baseMesh.Cf().boundaryField()[patchI][i]);
            }
        }
    }

    reduce(maxPlaneErr, maxOp<scalar>());
    reduce(nInterfaceFaces, sumOp<label>());

    const scalar planeTol = 1e-9*p_.period;

    if (nInterfaceFaces == 0 || maxPlaneErr > planeTol)
    {
        FatalErrorInFunction
            << "The heart-bath interface of mesh '" << baseMesh.name()
            << "' does not lie on the manufactured planes s = " << p_.s1
            << " and s = " << p_.s2 << " with normal (cos alpha, "
            << "sin alpha, 0): alpha " << radToDeg(p_.alpha)
            << " deg (geometryCoeffs) and period " << p_.period
            << " must match the mesh. Interface faces: " << nInterfaceFaces
            << "; maximum distance of a face centre from its plane: "
            << maxPlaneErr << " (tolerance " << planeTol << ")."
            << exit(FatalError);
    }

    const scalar thickness = baseMesh.bounds().span().z();
    const scalar heartArea = gSum(heartMesh.V())/thickness;
    const scalar expectedArea = (p_.s2 - p_.s1)*p_.period;

    if (mag(heartArea - expectedArea) > 1e-9*expectedArea)
    {
        FatalErrorInFunction
            << "The heart of mesh '" << heartMesh.name() << "' has area "
            << heartArea << " but (s2 - s1)*period = " << expectedArea
            << " (s1 " << p_.s1 << ", s2 " << p_.s2
            << ", period " << p_.period << ", alpha "
            << radToDeg(p_.alpha) << " deg in geometryCoeffs): the "
            << "geometry must match the mesh."
            << exit(FatalError);
    }
}


labelList obliqueSlabTransmission::decomposition
(
    const fvMesh& mesh,
    const labelUList& zone
) const
{
    const vectorField& C = mesh.C().primitiveField();

    // Manual decomposition into four domains so the interfaces carry processor and internal faces
    const vector tHat(-Foam::sin(p_.alpha), Foam::cos(p_.alpha), 0);
    const vector nHat(normal());
    const scalar sMid = 0.5*(p_.s1 + p_.s2);

    // Position along the interface in periods, shifted by 1e-6 off the cell centres
    auto tauMod = [&](const point& p)
    {
        const scalar tau = (tHat & p)/p_.period + 1e-6;
        return tau - std::floor(tau);
    };

    labelList domain(mesh.nCells(), -1);

    forAll(C, cellI)
    {
        if (zone[cellI] == transmissionMMS::myocardium)
        {
            domain[cellI] = tauMod(C[cellI]) < 0.5 ? 0 : 1;
        }
        else
        {
            domain[cellI] = (nHat & C[cellI]) < sMid ? 2 : 3;
        }
    }

    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();
    labelList movedDomain(domain);

    forAll(nei, faceI)
    {
        const label o = own[faceI];
        const label n = nei[faceI];
        const bool oHeart = zone[o] == transmissionMMS::myocardium;
        const bool nHeart = zone[n] == transmissionMMS::myocardium;

        if (oHeart == nHeart)
        {
            continue;
        }

        const label heartCell = oHeart ? o : n;
        const label bathCell = oHeart ? n : o;

        if (domain[heartCell] == 0 && tauMod(C[bathCell]) < 0.25)
        {
            movedDomain[bathCell] = 0;
        }
    }

    return movedDomain;
}

} // End namespace Foam

// ************************************************************************* //
