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

#include "sphericalShellTransmission.H"

#include "FixedList.H"
#include "IStringStream.H"
#include "syncTools.H"
#include "unitConversion.H"
#include "volFields.H"

namespace Foam
{

namespace
{

scalar sampleValue
(
    const shellTransmissionMMS::SelfTestSample& s,
    const shellTransmissionMMS::Params& p
)
{
    const point pt(s.x, s.y, s.z);
    const label c = s.component;

    switch (s.fn)
    {
        case 0: return shellTransmissionMMS::Vm(pt, s.t, p);
        case 1: return shellTransmissionMMS::phiE(pt, s.t, s.zone, p);
        case 2: return shellTransmissionMMS::phiI(pt, s.t, p);
        case 3:
            return shellTransmissionMMS::gradPhiE(pt, s.t, s.zone, p).component(c);
        case 4:
            return shellTransmissionMMS::gradVm(pt, s.t, p).component(c);
        case 5: return shellTransmissionMMS::u1(pt, s.t, p);
        case 6: return shellTransmissionMMS::u2(pt, s.t, p);
        case 7: return shellTransmissionMMS::Gi(pt, p).component(c);
        case 8: return shellTransmissionMMS::Ge(pt, p).component(c);
        case 9: return shellTransmissionMMS::sigmaB(pt, s.zone, p);
        case 10:
            return shellTransmissionMMS::divGiGradPhiI(pt, s.t, p);
        case 11:
            return shellTransmissionMMS::phiESource(pt, s.t, s.zone, p);
        case 12:
        {
            vector n(Zero);
            n.component(c) = 1;
            return shellTransmissionMMS::interfaceFlux(pt, n, s.t, p);
        }
        case 13: return shellTransmissionMMS::G(pt, p);
    }

    FatalErrorInFunction
        << "Self-test sample function index " << s.fn << " is not known."
        << exit(FatalError);

    return 0;
}

} // End anonymous namespace


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

scalar sphericalShellTransmission::rhoTolerance() const
{
    return 1e-8*sqr(p_.r2);
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

shellTransmissionMMS::Params sphericalShellTransmission::readParams
(
    const dictionary& dict
)
{
    shellTransmissionMMS::Params p;

    const word solution(dict.get<word>("solution"));
    if (solution != "levelSet")
    {
        FatalIOErrorInFunction(dict)
            << "Unknown solution '" << solution << "' for the sphericalShell "
            << "geometry. Valid solutions: (levelSet)."
            << exit(FatalIOError);
    }
    p.solution = 1;

    const dictionary& geom = dict.subDict("geometryCoeffs");
    p.r1 = geom.get<scalar>("r1");
    p.r2 = geom.get<scalar>("r2");
    if (!(p.r1 > 0 && p.r2 > p.r1))
    {
        FatalIOErrorInFunction(geom)
            << "The spheres need 0 < r1 < r2; found r1 " << p.r1 << ", r2 "
            << p.r2 << "." << exit(FatalIOError);
    }

    const dictionary& cond = dict.subDict("conductivityCoeffs");
    const vector gi(cond.get<vector>("gi"));
    const vector ge(cond.get<vector>("ge"));
    p.gil = gi.x();
    p.git = gi.y();
    p.gin = gi.z();
    p.gel = ge.x();
    p.get = ge.y();
    p.gen = ge.z();
    const FixedList<scalar, 2> frame(cond.get<FixedList<scalar, 2>>("frame"));
    p.gamma1 = degToRad(frame[0]);
    p.gamma2 = degToRad(frame[1]);
    p.theta0 = degToRad(cond.get<scalar>("theta0"));
    p.theta1 = degToRad(cond.get<scalar>("theta1"));
    p.sigma0In = cond.get<scalar>("sigma0In");
    p.sigma0Out = cond.get<scalar>("sigma0Out");
    p.sigmaEps = cond.get<scalar>("sigmaEps");

    const dictionary& sol = dict.subDict("solutionCoeffs");
    p.k = constant::mathematical::twoPi/sol.get<scalar>("wavelength");
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

    return p;
}


void sphericalShellTransmission::selfTest()
{
    IStringStream defaultsStream
    (
        "solution levelSet;"
        "geometry sphericalShell;"
        "geometryCoeffs { r1 0.5; r2 1; }"
        "conductivityCoeffs"
        "{"
        "    gi (0.17 0.019 0.01); ge (0.62 0.12 0.08);"
        "    frame (20 50); theta0 -20; theta1 90;"
        "    sigma0In 0.6; sigma0Out 0.25; sigmaEps 0.2;"
        "}"
        "solutionCoeffs"
        "{"
        "    wavelength 2;"
        "    f (0.3 0.2 0.4 0.15); h (-0.1 0.3 0.25);"
        "    commonMode { amplitude 0; omega 0; }"
        "}"
    );
    const dictionary defaultsDict(defaultsStream);

    const shellTransmissionMMS::Params read = readParams(defaultsDict);
    const shellTransmissionMMS::Params ref = shellTransmissionMMS::defaultParams();

    wordList differing;
    auto compare = [&differing](const char* name, scalar a, scalar b)
    {
        if (a != b)
        {
            differing.append(word(name));
        }
    };

    compare("r1", read.r1, ref.r1);
    compare("r2", read.r2, ref.r2);
    compare("gamma1", read.gamma1, ref.gamma1);
    compare("gamma2", read.gamma2, ref.gamma2);
    compare("k", read.k, ref.k);
    compare("gil", read.gil, ref.gil);
    compare("git", read.git, ref.git);
    compare("gin", read.gin, ref.gin);
    compare("gel", read.gel, ref.gel);
    compare("get", read.get, ref.get);
    compare("gen", read.gen, ref.gen);
    compare("theta0", read.theta0, ref.theta0);
    compare("theta1", read.theta1, ref.theta1);
    compare("sigma0In", read.sigma0In, ref.sigma0In);
    compare("sigma0Out", read.sigma0Out, ref.sigma0Out);
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
    compare("solution", read.solution, ref.solution);

    if (!differing.empty())
    {
        FatalErrorInFunction
            << "The verification-dictionary reader does not reproduce the "
            << "reference default parameters; differing fields: "
            << differing << "."
            << exit(FatalError);
    }

    const shellTransmissionMMS::Params paramSets[2] =
    {
        shellTransmissionMMS::defaultParams(),
        shellTransmissionMMS::alternateParams()
    };

    for (label i = 0; i < shellTransmissionMMS::nSelfTestSamples; ++i)
    {
        const shellTransmissionMMS::SelfTestSample& s =
            shellTransmissionMMS::selfTestSamples[i];

        if (s.param < 0 || s.param > 1)
        {
            FatalErrorInFunction
                << "Self-test sample " << i << " names parameter set "
                << s.param << "." << exit(FatalError);
        }

        const scalar value = sampleValue(s, paramSets[s.param]);
        const scalar err = mag(value - s.value);

        if (!(err <= 1e-11*max(scalar(1), mag(s.value))))
        {
            FatalErrorInFunction
                << "Shell transmission MMS self-test sample " << i
                << " (function " << s.fn << ", parameter set " << s.param
                << ", zone " << s.zone << ", component " << s.component
                << ", point (" << s.x << ' ' << s.y << ' ' << s.z << "), t "
                << s.t << ") evaluates to " << value << " against the "
                << "reference " << s.value << "."
                << exit(FatalError);
        }
    }

    Info<< "shell transmission MMS self-test: "
        << shellTransmissionMMS::nSelfTestSamples << " samples OK" << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

sphericalShellTransmission::sphericalShellTransmission(const dictionary& dict)
:
    p_(readParams(dict))
{
    selfTest();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

word sphericalShellTransmission::zoneName(const label zone) const
{
    static const word names[3] = {"bloodPool", "myocardium", "bath"};
    return names[zone];
}


vector sphericalShellTransmission::interfaceNormal(const point& pt) const
{
    return pt/max(mag(pt), VSMALL);
}


label sphericalShellTransmission::nMisplacedCells
(
    const fvMesh& mesh,
    const labelUList& zone
) const
{
    const scalar tol = rhoTolerance();
    const scalar lo[3] = {-tol, sqr(p_.r1) - tol, sqr(p_.r2) - tol};
    const scalar hi[3] = {sqr(p_.r1) + tol, sqr(p_.r2) + tol, GREAT};
    const pointField& points = mesh.points();
    label nMisplaced = 0;

    forAll(zone, cellI)
    {
        const label z = zone[cellI];
        if (z < 0)
        {
            continue;
        }

        for (const label pointI : mesh.cellPoints(cellI))
        {
            const scalar rho = magSqr(points[pointI]);
            if (rho < lo[z] || rho > hi[z])
            {
                ++nMisplaced;
                break;
            }
        }
    }

    return nMisplaced;
}


void sphericalShellTransmission::validateInterface
(
    const fvMesh& baseMesh,
    const labelUList& zone,
    const fvMesh&
) const
{
    // Heart-bath face vertices lie on the sphere of the bath zone; the two baths never touch

    const labelList nbrZone
    (
        syncTools::swapBoundaryCellList(baseMesh, zone)
    );
    const pointField& points = baseMesh.points();
    const faceList& faces = baseMesh.faces();

    scalar maxVertexErr = 0;
    scalar maxCentreDistance = 0;
    label nInterfaceFaces[2] = {0, 0};
    label nBathBathFaces = 0;

    // A processor face is counted on the rank of its heart (or first bath) cell
    auto check = [&](const label zA, const label zB, const label faceI)
    {
        if (zA == zB)
        {
            return;
        }
        if (zA != heartZone && zB != heartZone)
        {
            ++nBathBathFaces;
            return;
        }

        const label bathZone = (zA == heartZone ? zB : zA);
        const scalar r = (bathZone == shellTransmissionMMS::bloodPool ? p_.r1 : p_.r2);
        for (const label pointI : faces[faceI])
        {
            maxVertexErr = max(maxVertexErr, mag(magSqr(points[pointI]) - sqr(r)));
        }
        maxCentreDistance =
            max(maxCentreDistance, mag(mag(baseMesh.faceCentres()[faceI]) - r));
        ++nInterfaceFaces[bathZone == shellTransmissionMMS::bloodPool ? 0 : 1];
    };

    const labelUList& own = baseMesh.owner();
    const labelUList& nei = baseMesh.neighbour();

    forAll(nei, faceI)
    {
        check(zone[own[faceI]], zone[nei[faceI]], faceI);
    }

    forAll(baseMesh.boundary(), patchI)
    {
        const fvPatch& fp = baseMesh.boundary()[patchI];
        if (!fp.coupled())
        {
            continue;
        }

        const labelUList& faceCells = fp.faceCells();
        const label start = fp.patch().start();

        forAll(faceCells, i)
        {
            const label zLocal = zone[faceCells[i]];
            const label zNbr = nbrZone[start - baseMesh.nInternalFaces() + i];
            if (zLocal == heartZone || (zNbr != heartZone && zLocal < zNbr))
            {
                check(zLocal, zNbr, start + i);
            }
        }
    }

    reduce(maxVertexErr, maxOp<scalar>());
    reduce(maxCentreDistance, maxOp<scalar>());
    reduce(nInterfaceFaces[0], sumOp<label>());
    reduce(nInterfaceFaces[1], sumOp<label>());
    reduce(nBathBathFaces, sumOp<label>());

    const scalar tol = rhoTolerance();

    if
    (
        nInterfaceFaces[0] == 0
     || nInterfaceFaces[1] == 0
     || nBathBathFaces
     || maxVertexErr > tol
    )
    {
        FatalErrorInFunction
            << "The heart-bath interface of mesh '" << baseMesh.name()
            << "' does not lie on the manufactured spheres r1 = " << p_.r1
            << " and r2 = " << p_.r2 << " (geometryCoeffs). Interface faces "
            << "on r1: " << nInterfaceFaces[0] << ", on r2: "
            << nInterfaceFaces[1] << "; faces between the two bath zones: "
            << nBathBathFaces << "; maximum | |x|^2 - r^2 | of an interface "
            << "vertex: " << maxVertexErr << " (tolerance " << tol << ")."
            << exit(FatalError);
    }

    Info<< "sphericalShell: maximum distance of an interface face centre "
        << "from its sphere " << maxCentreDistance << endl;
}


labelList sphericalShellTransmission::decomposition
(
    const fvMesh& mesh,
    const labelUList& zone
) const
{
    const vectorField& C = mesh.C().primitiveField();

    // Heart halves on ranks 0 and 1, bloodPool on 2, bath on 3; part of the bath next to rank 0 joins rank 0
    const vector heartSplit(vector(1, 0.3, 0.2)/mag(vector(1, 0.3, 0.2)));
    const vector bathSplit(vector(0.2, 1, -0.3)/mag(vector(0.2, 1, -0.3)));

    labelList domain(mesh.nCells(), -1);

    forAll(C, cellI)
    {
        if (zone[cellI] == heartZone)
        {
            domain[cellI] = (heartSplit & C[cellI]) > 0 ? 0 : 1;
        }
        else
        {
            domain[cellI] = zone[cellI] == shellTransmissionMMS::bloodPool ? 2 : 3;
        }
    }

    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();
    labelList movedDomain(domain);

    forAll(nei, faceI)
    {
        const label o = own[faceI];
        const label n = nei[faceI];
        const bool oHeart = zone[o] == heartZone;
        const bool nHeart = zone[n] == heartZone;

        if (oHeart == nHeart)
        {
            continue;
        }

        const label heartCell = oHeart ? o : n;
        const label bathCell = oHeart ? n : o;

        if (domain[heartCell] == 0 && (bathSplit & C[bathCell]) > 0)
        {
            movedDomain[bathCell] = 0;
        }
    }

    return movedDomain;
}

} // End namespace Foam

// ************************************************************************* //
