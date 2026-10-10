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

#include "manufacturedTransmissionBathBidomainVerifier.H"

#include "FixedList.H"
#include "IStringStream.H"
#include "OFstream.H"
#include "OSspecific.H"
#include "dimVoltage.H"
#include "fixedValueFvPatchFields.H"
#include "fvcGrad.H"
#include "ionicModel.H"
#include "labelIOList.H"
#include "syncTools.H"
#include "verificationUtils.H"
#include "addToRunTimeSelectionTable.H"

#include <cmath>
#include <limits>

namespace Foam
{
using namespace verificationUtils;

defineTypeNameAndDebug(manufacturedTransmissionBathBidomainVerifier, 0);
addToRunTimeSelectionTable
(
    electroVerificationModel,
    manufacturedTransmissionBathBidomainVerifier,
    dictionary
);

namespace
{

//- Cell-zone names indexed by transmissionMMS::Zone.
const word zoneNames[3] = {"bathLow", "myocardium", "bathHigh"};

const word exposedPatchName("exposed");

const dimensionSet conductivityDim(-1, -3, 3, 0, 0, 2, 0);

struct errorNorms
{
    scalar L1;
    scalar L2;
    scalar Linf;
};


//- Weighted L1, L2 and maximum of the errors over the masked entries, reduced over ranks
errorNorms weightedNorms
(
    const scalarField& weights,
    const scalarField& e,
    const boolList& mask
)
{
    scalar sumAbs = 0;
    scalar sumSq = 0;
    scalar sumW = 0;
    scalar maxAbs = 0;

    forAll(e, i)
    {
        if (mask[i])
        {
            const scalar a = mag(e[i]);
            sumAbs += weights[i]*a;
            sumSq += weights[i]*a*a;
            sumW += weights[i];
            maxAbs = max(maxAbs, a);
        }
    }

    reduce(sumAbs, sumOp<scalar>());
    reduce(sumSq, sumOp<scalar>());
    reduce(sumW, sumOp<scalar>());
    reduce(maxAbs, maxOp<scalar>());

    if (sumW <= VSMALL)
    {
        const scalar nan = std::numeric_limits<scalar>::quiet_NaN();
        return {nan, nan, nan};
    }

    return {sumAbs/sumW, Foam::sqrt(sumSq/sumW), maxAbs};
}


errorNorms tupleNorms(const Tuple2<Tuple2<scalar, scalar>, scalar>& n)
{
    return {n.first().first(), n.first().second(), n.second()};
}


errorNorms singleValue(const scalar v)
{
    return {v, v, v};
}


scalar relativeError(const tensor& computed, const tensor& exact)
{
    return mag(computed - exact)/max(mag(exact), VSMALL);
}


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

const dictionary&
manufacturedTransmissionBathBidomainVerifier::verificationDict() const
{
    return dict().subDict("verificationModel");
}


transmissionMMS::Params
manufacturedTransmissionBathBidomainVerifier::readParams
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


void manufacturedTransmissionBathBidomainVerifier::selfTest()
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


labelList manufacturedTransmissionBathBidomainVerifier::buildZones
(
    const fvMesh& mesh,
    const transmissionMMS::Params& p
)
{
    labelList zone(mesh.nCells(), -1);

    label nMissingZones = 0;
    label nMultiplyZoned = 0;

    for (label z = 0; z < 3; ++z)
    {
        const label zoneI = mesh.cellZones().findZoneID(zoneNames[z]);
        if (zoneI < 0)
        {
            ++nMissingZones;
            continue;
        }

        for (const label cellI : mesh.cellZones()[zoneI])
        {
            if (zone[cellI] >= 0)
            {
                ++nMultiplyZoned;
            }
            zone[cellI] = z;
        }
    }

    label nUnzoned = 0;
    label nMismatched = 0;
    const vectorField& C = mesh.C().primitiveField();

    forAll(zone, cellI)
    {
        if (zone[cellI] < 0)
        {
            ++nUnzoned;
        }
        else if (transmissionMMS::zoneOf(C[cellI], p) != zone[cellI])
        {
            ++nMismatched;
        }
    }

    reduce(nMissingZones, maxOp<label>());
    reduce(nMultiplyZoned, sumOp<label>());
    reduce(nUnzoned, sumOp<label>());
    reduce(nMismatched, sumOp<label>());

    if (nMissingZones || nMultiplyZoned || nUnzoned || nMismatched)
    {
        FatalErrorInFunction
            << "Mesh '" << mesh.name() << "' does not match the "
            << "manufactured zones: the cellZones " << zoneNames[0] << ", "
            << zoneNames[1] << " and " << zoneNames[2]
            << " must exist and partition the cells. Missing zones on some "
            << "rank: " << nMissingZones << "; cells in several zones: "
            << nMultiplyZoned << "; cells in no zone: " << nUnzoned
            << "; cells whose centre lies in another zone: " << nMismatched
            << "."
            << exit(FatalError);
    }

    return zone;
}


scalar manufacturedTransmissionBathBidomainVerifier::ionicConstant
(
    const ionicModel& model,
    const word& name
)
{
    const wordList names = model.constantVariableNames();
    const scalarField values = model.constantVariableValues();

    if (names.size() != values.size())
    {
        FatalErrorInFunction
            << "Ionic model " << model.type()
            << " exposes inconsistent constant metadata: "
            << names.size() << " names but " << values.size() << " values."
            << exit(FatalError);
    }

    const label i = names.find(name);
    if (i < 0)
    {
        FatalErrorInFunction
            << typeName << " requires the manufactured ionic constant '"
            << name << "', but ionic model " << model.type()
            << " exposes constants " << names << "."
            << exit(FatalError);
    }

    return values[i];
}


void manufacturedTransmissionBathBidomainVerifier::validateConfiguration
(
    const ionicModel& model
) const
{
    const word family = model.verificationFamily();

    if
    (
        (
            family != "bidomainFDAManufactured"
         && family != "bathBidomainFDAManufactured"
        )
     || model.geometricDimension() != 2
    )
    {
        FatalErrorInFunction
            << typeName << " requires the 2D bidomainFDAManufactured or "
            << "bathBidomainFDAManufactured ionic model. Received ionic "
            << "model " << model.type() << " with verification family '"
            << family << "' and geometric dimension "
            << model.geometricDimension() << "."
            << exit(FatalError);
    }
}


void manufacturedTransmissionBathBidomainVerifier::initialiseIonicConstants
(
    const ionicModel& model,
    const label nLocalCells
)
{
    // Take the ionic constants from ranks with heart cells
    const bool hasCells = nLocalCells > 0;
    FixedList<scalar, 4> sums
    ({
        scalar(hasCells),
        hasCells ? ionicConstant(model, "Beta") : scalar(0),
        hasCells ? ionicConstant(model, "Chi") : scalar(0),
        hasCells ? ionicConstant(model, "Cm") : scalar(0)
    });
    for (scalar& v : sums)
    {
        reduce(v, sumOp<scalar>());
    }

    if (sums[0] < 0.5)
    {
        FatalErrorInFunction
            << typeName << " found no heart cells on any rank."
            << exit(FatalError);
    }

    const scalar Beta = sums[1]/sums[0];
    const scalar Chi = sums[2]/sums[0];
    const scalar Cm = sums[3]/sums[0];

    const dimensionedScalar chi("chi", dimArea/dimVolume, dict());
    const dimensionedScalar cm
    (
        "cm",
        dimCurrent*dimTime/(dimVoltage*dimArea),
        dict()
    );

    const scalar chiErr = mag(chi.value() - Chi)/max(mag(Chi), VSMALL);
    const scalar cmErr = mag(cm.value() - Cm)/max(mag(Cm), VSMALL);

    const scalar ratio = (chi.value()*cm.value())/(Chi*Cm);

    if (!betaInitialised_ && (chiErr > 1e-12 || cmErr > 1e-12))
    {
        WarningInFunction
            << "chi and cm differ from the ionic model constants Chi and Cm"
            << " (chi = " << chi.value() << ", Chi = " << Chi
            << ", cm = " << cm.value() << ", Cm = " << Cm << "). The"
            << " manufactured Vm source carries the ratio (chi*cm)/(Chi*Cm) = "
            << ratio << "." << endl;
    }

    beta_ = Beta;
    capacitanceRatio_ = ratio;
    betaInitialised_ = true;
}


void manufacturedTransmissionBathBidomainVerifier::validateBoundBath
(
    const volScalarField& Vm
) const
{
    const fvMesh& phiEMesh = phiEPtr_->mesh();
    const fvMesh& heartMesh = Vm.mesh();

    if (phiEMesh.nSolutionD() != 2)
    {
        FatalErrorInFunction
            << typeName << " requires a 2D extracellular-potential mesh; "
            << "mesh '" << phiEMesh.name() << "' solves in "
            << phiEMesh.nSolutionD() << " directions."
            << exit(FatalError);
    }

    const word exposedName
    (
        dict().getOrDefault<word>("exposedFacesPatch", word::null)
    );
    if (exposedName != exposedPatchName)
    {
        FatalIOErrorInFunction(dict())
            << typeName << " requires exposedFacesPatch " << exposedPatchName
            << "; found '" << exposedName << "'."
            << exit(FatalIOError);
    }

    const label exposedI =
        phiEMesh.boundaryMesh().findPatchID(exposedPatchName);
    const label missing = returnReduce(label(exposedI < 0), maxOp<label>());
    const label exposedSize = returnReduce
    (
        exposedI < 0 ? label(0) : phiEMesh.boundaryMesh()[exposedI].size(),
        sumOp<label>()
    );

    if (missing || exposedSize != 0)
    {
        FatalErrorInFunction
            << typeName << " requires an empty patch '" << exposedPatchName
            << "' on mesh '" << phiEMesh.name() << "'. Missing on some "
            << "rank: " << missing << "; faces: " << exposedSize << "."
            << exit(FatalError);
    }

    const volTensorField& GiField =
        heartMesh.lookupObject<volTensorField>("conductivityIntracellular");
    const volTensorField& GeField =
        heartMesh.lookupObject<volTensorField>("conductivityExtracellular");
    const volTensorField& sigmaTotal =
        phiEMesh.lookupObject<volTensorField>("sigmaTotal");

    const vectorField& heartC = heartMesh.C().primitiveField();
    scalar maxGiErr = 0;
    scalar maxGeErr = 0;

    forAll(heartC, cellI)
    {
        const point& c = heartC[cellI];
        maxGiErr = max
        (
            maxGiErr,
            relativeError
            (
                GiField.primitiveField()[cellI],
                tensor(transmissionMMS::Gi(c, params_))
            )
        );
        maxGeErr = max
        (
            maxGeErr,
            relativeError
            (
                GeField.primitiveField()[cellI],
                tensor(transmissionMMS::Ge(c, params_))
            )
        );
    }

    const vectorField& baseC = phiEMesh.C().primitiveField();
    scalar maxSigmaErr = 0;

    forAll(baseC, cellI)
    {
        const label z = zoneOfCell_[cellI];
        if (z != transmissionMMS::myocardium)
        {
            maxSigmaErr = max
            (
                maxSigmaErr,
                relativeError
                (
                    sigmaTotal.primitiveField()[cellI],
                    transmissionMMS::sigmaB(baseC[cellI], z, params_)
                   *tensor::I
                )
            );
        }
    }

    reduce(maxGiErr, maxOp<scalar>());
    reduce(maxGeErr, maxOp<scalar>());
    reduce(maxSigmaErr, maxOp<scalar>());

    if (maxGiErr > 1e-10 || maxGeErr > 1e-10 || maxSigmaErr > 1e-10)
    {
        FatalErrorInFunction
            << "The conductivities differ from the manufactured ones "
            << "(maximum relative error, tolerance 1e-10): "
            << "conductivityIntracellular " << maxGiErr
            << ", conductivityExtracellular " << maxGeErr
            << ", bath sigmaTotal " << maxSigmaErr
            << ". Write them with the manufactured-field utility."
            << exit(FatalError);
    }

    // Heart-bath faces lie on s = s1 and s = s2; heart cells fill s1 < s < s2 over one period

    const vector nHat(Foam::cos(params_.alpha), Foam::sin(params_.alpha), 0);
    const labelList nbrZone
    (
        syncTools::swapBoundaryCellList(phiEMesh, zoneOfCell_)
    );

    scalar maxPlaneErr = 0;
    label nInterfaceFaces = 0;

    auto checkPlane = [&](const label bathZone, const point& Cf)
    {
        const scalar sExact =
            (bathZone == transmissionMMS::bathLow ? params_.s1 : params_.s2);
        maxPlaneErr = max(maxPlaneErr, mag((nHat & Cf) - sExact));
        ++nInterfaceFaces;
    };

    const labelUList& own = phiEMesh.owner();
    const labelUList& nei = phiEMesh.neighbour();

    forAll(nei, faceI)
    {
        const label zo = zoneOfCell_[own[faceI]];
        const label zn = zoneOfCell_[nei[faceI]];

        if (zo == zn)
        {
            continue;
        }
        if (zo == transmissionMMS::myocardium)
        {
            checkPlane(zn, phiEMesh.Cf()[faceI]);
        }
        else if (zn == transmissionMMS::myocardium)
        {
            checkPlane(zo, phiEMesh.Cf()[faceI]);
        }
    }

    forAll(phiEMesh.boundary(), patchI)
    {
        const fvPatch& fp = phiEMesh.boundary()[patchI];
        if (!fp.coupled())
        {
            continue;
        }

        const labelUList& faceCells = fp.faceCells();
        const label start = fp.patch().start() - phiEMesh.nInternalFaces();

        forAll(faceCells, i)
        {
            const label zLocal = zoneOfCell_[faceCells[i]];
            const label zNbr = nbrZone[start + i];

            if (zLocal == zNbr)
            {
                continue;
            }
            if (zLocal == transmissionMMS::myocardium)
            {
                checkPlane(zNbr, phiEMesh.Cf().boundaryField()[patchI][i]);
            }
            else if (zNbr == transmissionMMS::myocardium)
            {
                checkPlane(zLocal, phiEMesh.Cf().boundaryField()[patchI][i]);
            }
        }
    }

    reduce(maxPlaneErr, maxOp<scalar>());
    reduce(nInterfaceFaces, sumOp<label>());

    const scalar planeTol = 1e-9*params_.period;

    if (nInterfaceFaces == 0 || maxPlaneErr > planeTol)
    {
        FatalErrorInFunction
            << "The heart-bath interface of mesh '" << phiEMesh.name()
            << "' does not lie on the manufactured planes s = " << params_.s1
            << " and s = " << params_.s2 << " with normal (cos alpha, "
            << "sin alpha, 0): alpha " << radToDeg(params_.alpha)
            << " deg (geometryCoeffs) and period " << params_.period
            << " must match the mesh. Interface faces: " << nInterfaceFaces
            << "; maximum distance of a face centre from its plane: "
            << maxPlaneErr << " (tolerance " << planeTol << ")."
            << exit(FatalError);
    }

    const scalar thickness = phiEMesh.bounds().span().z();
    const scalar heartArea = gSum(heartMesh.V())/thickness;
    const scalar expectedArea = (params_.s2 - params_.s1)*params_.period;

    if (mag(heartArea - expectedArea) > 1e-9*expectedArea)
    {
        FatalErrorInFunction
            << "The heart of mesh '" << heartMesh.name() << "' has area "
            << heartArea << " but (s2 - s1)*period = " << expectedArea
            << " (s1 " << params_.s1 << ", s2 " << params_.s2
            << ", period " << params_.period << ", alpha "
            << radToDeg(params_.alpha) << " deg in geometryCoeffs): the "
            << "geometry must match the mesh."
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

manufacturedTransmissionBathBidomainVerifier::
manufacturedTransmissionBathBidomainVerifier
(
    const dictionary& dict
)
:
    electroVerificationModel(dict),
    phiEPtr_(nullptr),
    phiEHeartCellMapPtr_(nullptr),
    errorsReported_(false),
    params_(),
    solutionName_(),
    nPerDirection_(0),
    dirichletPatches_(),
    beta_(0),
    capacitanceRatio_(1),
    betaInitialised_(false),
    zoneOfCell_(),
    zoneMeshPtr_(nullptr)
{
    selfTest();

    const dictionary& cfg = verificationDict();

    params_ = readParams(cfg);
    solutionName_ = cfg.get<word>("solution");

    nPerDirection_ = cfg.get<label>("N");
    if (nPerDirection_ <= 0)
    {
        FatalIOErrorInFunction(cfg)
            << "N must be positive; found " << nPerDirection_ << "."
            << exit(FatalIOError);
    }

    const dictionary& bc = cfg.subDict("boundaryData");
    dirichletPatches_ = bc.get<wordList>("dirichlet");
    if (dirichletPatches_.empty())
    {
        FatalIOErrorInFunction(bc)
            << "boundaryData.dirichlet names no patch."
            << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void manufacturedTransmissionBathBidomainVerifier::bindBidomainField
(
    volScalarField& phiE
)
{
    phiEPtr_ = &phiE;
    phiEHeartCellMapPtr_ = nullptr;
}


void manufacturedTransmissionBathBidomainVerifier::bindBidomainField
(
    volScalarField& phiE,
    const labelUList& heartCellMap
)
{
    phiEPtr_ = &phiE;
    phiEHeartCellMapPtr_ = &heartCellMap;
}


void manufacturedTransmissionBathBidomainVerifier::unbindBidomainField()
{
    phiEPtr_ = nullptr;
    phiEHeartCellMapPtr_ = nullptr;
}


wordList manufacturedTransmissionBathBidomainVerifier::preProcessFieldNames
(
    const ionicModel&
) const
{
    return wordList({"u1", "u2", "u3"});
}


wordList
manufacturedTransmissionBathBidomainVerifier::requiredPostProcessFieldNames
(
    const ionicModel&
) const
{
    return wordList({"u1", "u2", "u3"});
}


bool manufacturedTransmissionBathBidomainVerifier::shouldPostProcess
(
    const ionicModel&,
    const volScalarField& Vm
) const
{
    return !errorsReported_ && shouldReportManufacturedErrors(Vm);
}


void manufacturedTransmissionBathBidomainVerifier::preProcess
(
    ionicModel& model,
    volScalarField& Vm,
    PtrList<volScalarField>& fields
)
{
    const wordList names = preProcessFieldNames(model);

    if (fields.size() != names.size())
    {
        FatalErrorInFunction
            << type() << " preProcess expected " << names.size()
            << " fields " << names << " but received " << fields.size()
            << exit(FatalError);
    }

    validateConfiguration(model);
    initialiseIonicConstants(model, Vm.mesh().nCells());

    if (phiEHeartCellMapPtr_)
    {
        if (!phiEPtr_ || &phiEPtr_->mesh() == &Vm.mesh())
        {
            FatalErrorInFunction
                << typeName << " requires a bath: the extracellular "
                << "potential must be bound on the whole heart-and-bath mesh."
                << exit(FatalError);
        }

        zoneOfCell_ = buildZones(phiEPtr_->mesh(), params_);
        zoneMeshPtr_ = &phiEPtr_->mesh();

        validateBoundBath(Vm);
    }

    const label u1Field = requireFieldIndex(names, "u1", "preProcess");
    const label u2Field = requireFieldIndex(names, "u2", "preProcess");
    const label u3Field = requireFieldIndex(names, "u3", "preProcess");

    scalarField& VmI = Vm.primitiveFieldRef();
    scalarField& u1I = fields[u1Field].primitiveFieldRef();
    scalarField& u2I = fields[u2Field].primitiveFieldRef();
    scalarField& u3I = fields[u3Field].primitiveFieldRef();

    const vectorField& C = Vm.mesh().C().primitiveField();
    const scalar t = Vm.mesh().time().value();

    forAll(C, cellI)
    {
        VmI[cellI] = transmissionMMS::Vm(C[cellI], t, params_);
        u1I[cellI] = transmissionMMS::u1(C[cellI], t, params_);
        u2I[cellI] = transmissionMMS::u2(C[cellI], t, params_);
        u3I[cellI] = 0;
    }

    Vm.correctBoundaryConditions();
    fields[u1Field].correctBoundaryConditions();
    fields[u2Field].correctBoundaryConditions();
    fields[u3Field].correctBoundaryConditions();

    if (phiEPtr_)
    {
        const fvMesh& phiEMesh = phiEPtr_->mesh();
        const bool zoned = (zoneMeshPtr_ == &phiEMesh);

        if (!zoned && &phiEMesh != &Vm.mesh())
        {
            FatalErrorInFunction
                << "The bound extracellular potential lives on mesh '"
                << phiEMesh.name() << "', which is neither the myocardium "
                << "mesh nor a zoned heart-and-bath mesh."
                << exit(FatalError);
        }

        const vectorField& phiEC = phiEMesh.C().primitiveField();
        scalarField& phiEI = phiEPtr_->primitiveFieldRef();

        forAll(phiEC, cellI)
        {
            const label z =
                zoned ? zoneOfCell_[cellI] : label(transmissionMMS::myocardium);
            phiEI[cellI] = transmissionMMS::phiE(phiEC[cellI], t, z, params_);
        }

        phiEPtr_->correctBoundaryConditions();
    }

    model.importFields(Vm, names, fields);
}


void manufacturedTransmissionBathBidomainVerifier::addManufacturedPdeSource
(
    volScalarField& sourceField,
    const volTensorField& conductivity,
    const scalar evaluationTime
) const
{
    if (!betaInitialised_)
    {
        FatalErrorInFunction
            << "The manufactured source was requested before preProcess "
            << "read the ionic constant Beta."
            << exit(FatalError);
    }

    if (&sourceField.mesh() != &conductivity.mesh())
    {
        FatalErrorInFunction
            << "Source field and conductivity field belong to different "
            << "finite-volume meshes."
            << exit(FatalError);
    }

    const vectorField& C = sourceField.mesh().C().primitiveField();
    scalarField& source = sourceField.primitiveFieldRef();

    forAll(C, cellI)
    {
        source[cellI] +=
            capacitanceRatio_*beta_
           *transmissionMMS::Vm(C[cellI], evaluationTime, params_)
          - transmissionMMS::divGiGradPhiI(C[cellI], evaluationTime, params_);
    }

    sourceField.correctBoundaryConditions();
}


void manufacturedTransmissionBathBidomainVerifier::addManufacturedPhiESource
(
    volScalarField& rhs,
    const scalar evaluationTime
) const
{
    if (zoneMeshPtr_ != &rhs.mesh())
    {
        FatalErrorInFunction
            << "The extracellular-potential source was requested on mesh '"
            << rhs.mesh().name() << "' before preProcess built the "
            << "manufactured zones of that mesh."
            << exit(FatalError);
    }

    const vectorField& C = rhs.mesh().C().primitiveField();
    scalarField& r = rhs.primitiveFieldRef();

    forAll(C, cellI)
    {
        r[cellI] += transmissionMMS::phiESource
        (
            C[cellI], evaluationTime, zoneOfCell_[cellI], params_
        );
    }
}


void manufacturedTransmissionBathBidomainVerifier::applyPhiEBoundaryData
(
    volScalarField& phiE,
    const scalar evaluationTime
) const
{
    const fvMesh& mesh = phiE.mesh();

    for (const word& patchName : dirichletPatches_)
    {
        const label patchI = mesh.boundaryMesh().findPatchID(patchName);

        if (patchI < 0)
        {
            FatalErrorInFunction
                << "Dirichlet patch '" << patchName << "' is not a patch of "
                << "mesh '" << mesh.name() << "'."
                << exit(FatalError);
        }

        fvPatchScalarField& pf = phiE.boundaryFieldRef()[patchI];

        if (!isA<fixedValueFvPatchScalarField>(pf))
        {
            FatalErrorInFunction
                << "Dirichlet patch '" << patchName << "' of "
                << phiE.name() << " is " << pf.type()
                << ", not fixedValue."
                << exit(FatalError);
        }

        const vectorField& Cf = mesh.Cf().boundaryField()[patchI];
        scalarField values(Cf.size());

        forAll(Cf, faceI)
        {
            values[faceI] = transmissionMMS::phiE
            (
                Cf[faceI],
                evaluationTime,
                transmissionMMS::zoneOf(Cf[faceI], params_),
                params_
            );
        }

        pf == values;
    }
}


wordList manufacturedTransmissionBathBidomainVerifier::boundaryDataPatches()
const
{
    return dirichletPatches_;
}


bool manufacturedTransmissionBathBidomainVerifier::
requiresPhiECompatibilityProjection() const
{
    return false;
}


bool manufacturedTransmissionBathBidomainVerifier::
requiresInterfaceFluxDiagnostic() const
{
    return true;
}


void manufacturedTransmissionBathBidomainVerifier::writeManufacturedFields
(
    const fvMesh& mesh
) const
{
    if (Pstream::parRun())
    {
        FatalErrorInFunction
            << typeName << " writes the manufactured fields and the cell "
            << "decomposition of the undecomposed case; run it in serial."
            << exit(FatalError);
    }

    const labelList zone = buildZones(mesh, params_);
    const vectorField& C = mesh.C().primitiveField();
    const word startTime("0");

    const unsigned int oldPrecision = IOstream::defaultPrecision(17);

    volSymmTensorField GiField
    (
        IOobject
        (
            "ConductivityIntracellular",
            startTime,
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedSymmTensor(conductivityDim, Zero),
        "calculated"
    );
    volSymmTensorField GeField
    (
        IOobject
        (
            "ConductivityExtracellular",
            startTime,
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedSymmTensor(conductivityDim, Zero),
        "calculated"
    );
    volScalarField sigmaBField
    (
        IOobject
        (
            "bodyAndOrgansConductivity",
            startTime,
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(conductivityDim, Zero),
        "calculated"
    );

    symmTensorField& GiI = GiField.primitiveFieldRef();
    symmTensorField& GeI = GeField.primitiveFieldRef();
    scalarField& sigmaBI = sigmaBField.primitiveFieldRef();

    forAll(C, cellI)
    {
        GiI[cellI] = transmissionMMS::Gi(C[cellI], params_);
        GeI[cellI] = transmissionMMS::Ge(C[cellI], params_);
        sigmaBI[cellI] =
            zone[cellI] == transmissionMMS::myocardium
          ? 1.0
          : transmissionMMS::sigmaB(C[cellI], zone[cellI], params_);
    }

    forAll(mesh.boundary(), patchI)
    {
        if (polyPatch::constraintType(mesh.boundaryMesh()[patchI].type()))
        {
            continue;
        }

        GiField.boundaryFieldRef()[patchI] ==
            GiField.boundaryField()[patchI].patchInternalField();
        GeField.boundaryFieldRef()[patchI] ==
            GeField.boundaryField()[patchI].patchInternalField();
        sigmaBField.boundaryFieldRef()[patchI] ==
            sigmaBField.boundaryField()[patchI].patchInternalField();
    }

    GiField.write();
    GeField.write();
    sigmaBField.write();

    IOstream::defaultPrecision(oldPrecision);

    // Manual decomposition into four domains so the interfaces carry processor and internal faces
    const vector tHat(-Foam::sin(params_.alpha), Foam::cos(params_.alpha), 0);
    const vector nHat(Foam::cos(params_.alpha), Foam::sin(params_.alpha), 0);
    const scalar sMid = 0.5*(params_.s1 + params_.s2);

    // Position along the interface in periods, shifted by 1e-6 off the cell centres
    auto tauMod = [&](const point& p)
    {
        const scalar tau = (tHat & p)/params_.period + 1e-6;
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

    labelIOList decomposition
    (
        IOobject
        (
            "cellDecomposition",
            mesh.time().constant(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        movedDomain
    );
    decomposition.write();

    Info<< typeName << ": wrote ConductivityIntracellular, "
        << "ConductivityExtracellular and bodyAndOrgansConductivity in "
        << startTime << " and cellDecomposition in "
        << mesh.time().constant() << "." << endl;
}


void manufacturedTransmissionBathBidomainVerifier::postProcess
(
    const ionicModel& model,
    const volScalarField& Vm,
    const PtrList<volScalarField>& fields
)
{
    if (!shouldPostProcess(model, Vm))
    {
        return;
    }

    if (!phiEPtr_ || !phiEHeartCellMapPtr_ || zoneMeshPtr_ != &phiEPtr_->mesh())
    {
        FatalErrorInFunction
            << type() << " requires the extracellular potential bound on the "
            << "zoned heart-and-bath mesh before postProcess()."
            << exit(FatalError);
    }

    const wordList requiredNames = requiredPostProcessFieldNames(model);
    if (fields.size() != requiredNames.size())
    {
        FatalErrorInFunction
            << type() << " expected " << requiredNames.size()
            << " postProcess fields " << requiredNames
            << " but received " << fields.size()
            << exit(FatalError);
    }

    const label u1Field = requireFieldIndex(requiredNames, "u1", "postProcess");
    const label u2Field = requireFieldIndex(requiredNames, "u2", "postProcess");

    const fvMesh& heartMesh = Vm.mesh();
    const Time& time = heartMesh.time();
    const scalar t = time.value();

    const volScalarField& phiE = *phiEPtr_;
    const fvMesh& mesh = phiE.mesh();
    const labelUList& heartCellMap = *phiEHeartCellMapPtr_;
    const vector nHat(Foam::cos(params_.alpha), Foam::sin(params_.alpha), 0);
    const tensor tangential(tensor::I - sqr(nHat));

    // Heart mesh: Vm, phiI, u1, u2

    const vectorField& heartC = heartMesh.C().primitiveField();
    const scalarField& VmI = Vm.primitiveField();
    const scalarField& phiEI = phiE.primitiveField();

    if (heartCellMap.size() != VmI.size())
    {
        FatalErrorInFunction
            << "Heart cell map size " << heartCellMap.size()
            << " does not match the Vm cell count " << VmI.size() << "."
            << exit(FatalError);
    }

    scalarField VmExact(heartC.size());
    scalarField phiIExact(heartC.size());
    scalarField phiIValues(heartC.size());
    scalarField u1Exact(heartC.size());
    scalarField u2Exact(heartC.size());

    forAll(heartC, cellI)
    {
        VmExact[cellI] = transmissionMMS::Vm(heartC[cellI], t, params_);
        phiIExact[cellI] = transmissionMMS::phiI(heartC[cellI], t, params_);
        u1Exact[cellI] = transmissionMMS::u1(heartC[cellI], t, params_);
        u2Exact[cellI] = transmissionMMS::u2(heartC[cellI], t, params_);
        phiIValues[cellI] = VmI[cellI] + phiEI[heartCellMap[cellI]];
    }

    const errorNorms VmNorms =
        tupleNorms(computeNorms(heartMesh, VmI, VmExact));
    const errorNorms phiINorms =
        tupleNorms(computeNorms(heartMesh, phiIValues, phiIExact));
    const errorNorms u1Norms = tupleNorms
    (
        computeNorms(heartMesh, fields[u1Field].primitiveField(), u1Exact)
    );
    const errorNorms u2Norms = tupleNorms
    (
        computeNorms(heartMesh, fields[u2Field].primitiveField(), u2Exact)
    );

    // Heart trace gradient on the cells of the exposed (interface) patch

    boolList heartTraceCell(heartMesh.nCells(), false);
    const label heartExposedI =
        heartMesh.boundaryMesh().findPatchID(exposedPatchName);
    if (heartExposedI >= 0)
    {
        const polyPatch& exposed = heartMesh.boundaryMesh()[heartExposedI];
        for (const label cellI : exposed.faceCells())
        {
            heartTraceCell[cellI] = true;
        }
    }
    if (!returnReduceOr(heartExposedI >= 0))
    {
        FatalErrorInFunction
            << "The myocardium mesh has no patch '" << exposedPatchName
            << "' holding the heart-bath interface faces."
            << exit(FatalError);
    }

    const tmp<volVectorField> tgradVm(fvc::grad(Vm));
    const vectorField& gradVmI = tgradVm().primitiveField();
    scalarField heartTraceErr(heartC.size(), 0);
    forAll(heartC, cellI)
    {
        if (heartTraceCell[cellI])
        {
            heartTraceErr[cellI] = mag
            (
                gradVmI[cellI]
              - transmissionMMS::gradVm(heartC[cellI], t, params_)
            );
        }
    }
    const errorNorms heartTraceNorms =
        weightedNorms(heartMesh.V(), heartTraceErr, heartTraceCell);

    // Base mesh: interface faces and bands

    const vectorField& C = mesh.C().primitiveField();
    const scalarField& V = mesh.V();
    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();
    const surfaceScalarField& interfaceFlux =
        mesh.lookupObject<surfaceScalarField>("phiEInterfaceFlux");

    const labelList nbrZone
    (
        syncTools::swapBoundaryCellList(mesh, zoneOfCell_)
    );

    boolList interfaceCell(mesh.nCells(), false);
    boolList boundaryCell(mesh.nCells(), false);

    DynamicList<scalar> faceErr;
    DynamicList<scalar> faceArea;
    scalar currentComputed[2] = {0, 0};
    scalar currentExact[2] = {0, 0};

    auto addInterfaceFace = [&]
    (
        const label bathZone,
        const point& Cf,
        const vector& Sf,
        const scalar sign,
        const scalar flux
    )
    {
        const scalar magSf = mag(Sf);
        const vector nOut = sign*Sf/magSf;
        const scalar exact =
            transmissionMMS::interfaceFlux(Cf, nOut, t, params_)*magSf;

        faceErr.append((sign*flux - exact)/magSf);
        faceArea.append(magSf);

        const label j = (bathZone == transmissionMMS::bathLow ? 0 : 1);
        currentComputed[j] += sign*flux;
        currentExact[j] += exact;
    };

    forAll(nei, faceI)
    {
        const label o = own[faceI];
        const label n = nei[faceI];
        const label zo = zoneOfCell_[o];
        const label zn = zoneOfCell_[n];

        if
        (
            zo == zn
         || (
                zo != transmissionMMS::myocardium
             && zn != transmissionMMS::myocardium
            )
        )
        {
            continue;
        }

        interfaceCell[o] = true;
        interfaceCell[n] = true;

        const bool ownerHeart = (zo == transmissionMMS::myocardium);
        addInterfaceFace
        (
            ownerHeart ? zn : zo,
            mesh.Cf()[faceI],
            mesh.Sf()[faceI],
            ownerHeart ? 1.0 : -1.0,
            interfaceFlux[faceI]
        );
    }

    forAll(mesh.boundary(), patchI)
    {
        const fvPatch& fp = mesh.boundary()[patchI];
        if (!fp.coupled())
        {
            continue;
        }

        const labelUList& faceCells = fp.faceCells();
        const label start = fp.patch().start() - mesh.nInternalFaces();

        forAll(faceCells, i)
        {
            const label cellI = faceCells[i];
            const label zLocal = zoneOfCell_[cellI];
            const label zNbr = nbrZone[start + i];

            if
            (
                zLocal == zNbr
             || (
                    zLocal != transmissionMMS::myocardium
                 && zNbr != transmissionMMS::myocardium
                )
            )
            {
                continue;
            }

            interfaceCell[cellI] = true;

            if (zLocal == transmissionMMS::myocardium)
            {
                addInterfaceFace
                (
                    zNbr,
                    mesh.Cf().boundaryField()[patchI][i],
                    mesh.Sf().boundaryField()[patchI][i],
                    1.0,
                    interfaceFlux.boundaryField()[patchI][i]
                );
            }
        }
    }

    const errorNorms interfaceFluxNorms = weightedNorms
    (
        scalarField(faceArea),
        scalarField(faceErr),
        boolList(faceErr.size(), true)
    );

    reduce(currentComputed[0], sumOp<scalar>());
    reduce(currentComputed[1], sumOp<scalar>());
    reduce(currentExact[0], sumOp<scalar>());
    reduce(currentExact[1], sumOp<scalar>());

    const errorNorms current1 =
        singleValue(mag(currentComputed[0] - currentExact[0]));
    const errorNorms current2 =
        singleValue(mag(currentComputed[1] - currentExact[1]));

    // Outer Dirichlet patches: band cells and current

    const volTensorField& sigmaTotal =
        mesh.lookupObject<volTensorField>("sigmaTotal");
    scalar outerComputed = 0;
    scalar outerExact = 0;

    for (const word& patchName : dirichletPatches_)
    {
        const label patchI = mesh.boundaryMesh().findPatchID(patchName);
        if (patchI < 0)
        {
            FatalErrorInFunction
                << "Dirichlet patch '" << patchName << "' is not a patch of "
                << "mesh '" << mesh.name() << "'."
                << exit(FatalError);
        }

        const fvPatch& fp = mesh.boundary()[patchI];
        for (const label cellI : fp.faceCells())
        {
            boundaryCell[cellI] = true;
        }

        const vectorField nf(fp.nf());
        const scalarField& magSf = fp.magSf();
        const vectorField& Cf = fp.Cf();
        const scalarField snGrad(phiE.boundaryField()[patchI].snGrad());
        const tensorField sigmaOwner
        (
            sigmaTotal.boundaryField()[patchI].patchInternalField()
        );

        forAll(Cf, i)
        {
            outerComputed +=
                (nf[i] & sigmaOwner[i] & nf[i])*snGrad[i]*magSf[i];

            const label z = transmissionMMS::zoneOf(Cf[i], params_);
            outerExact +=
                transmissionMMS::sigmaB(Cf[i], z, params_)
               *(nf[i] & transmissionMMS::gradPhiE(Cf[i], t, z, params_))
               *magSf[i];
        }
    }

    reduce(outerComputed, sumOp<scalar>());
    reduce(outerExact, sumOp<scalar>());
    const errorNorms outerCurrent =
        singleValue(mag(outerComputed - outerExact));

    // Base-mesh phiE and tangential gradient

    const tmp<volVectorField> tgradPhiE(fvc::grad(phiE));
    const vectorField& gradPhiEI = tgradPhiE().primitiveField();

    scalarField phiEErr(C.size());
    scalarField gradTangentialErr(C.size(), 0);
    boolList heartCell(C.size(), false);
    boolList bathCell(C.size(), false);

    forAll(C, cellI)
    {
        const label z = zoneOfCell_[cellI];
        phiEErr[cellI] =
            phiEI[cellI] - transmissionMMS::phiE(C[cellI], t, z, params_);
        heartCell[cellI] = (z == transmissionMMS::myocardium);
        bathCell[cellI] = !heartCell[cellI];

        if (interfaceCell[cellI])
        {
            gradTangentialErr[cellI] = mag
            (
                tangential
              & (
                    gradPhiEI[cellI]
                  - transmissionMMS::gradPhiE(C[cellI], t, z, params_)
                )
            );
        }
    }

    const errorNorms phiEHeartNorms = weightedNorms(V, phiEErr, heartCell);
    const errorNorms phiEBathNorms = weightedNorms(V, phiEErr, bathCell);
    const errorNorms phiEInterfaceNorms =
        weightedNorms(V, phiEErr, interfaceCell);
    const errorNorms phiEBoundaryNorms =
        weightedNorms(V, phiEErr, boundaryCell);
    const errorNorms gradTangentialNorms =
        weightedNorms(V, gradTangentialErr, interfaceCell);

    // Header data

    const scalar thickness = mesh.bounds().span().z();
    const scalar heartArea = gSum(heartMesh.V())/thickness;
    const label nHeartCells = returnReduce(heartMesh.nCells(), sumOp<label>());
    const scalar h = Foam::sqrt(heartArea/max(nHeartCells, label(1)));
    const scalar dt = time.deltaTValue();
    const label nSteps = max(label(0), time.timeIndex());

    const fileName outputDir(time.globalPath()/"postProcessing");
    const fileName outputFile =
        outputDir
      / (dimensionName(2) + "_" + Foam::name(nPerDirection_) + "_cells.dat");

    const List<Tuple2<word, errorNorms>> rows
    ({
        {"Vm", VmNorms},
        {"phiI", phiINorms},
        {"phiE_heart", phiEHeartNorms},
        {"phiE_bath", phiEBathNorms},
        {"phiE_interfaceBand", phiEInterfaceNorms},
        {"phiE_boundaryBand", phiEBoundaryNorms},
        {"u1", u1Norms},
        {"u2", u2Norms},
        {"interfaceFlux", interfaceFluxNorms},
        {"interfaceCurrent_1", current1},
        {"interfaceCurrent_2", current2},
        {"gradTangentialError", gradTangentialNorms},
        {"gradHeartTraceError", heartTraceNorms},
        {"outerCurrent", outerCurrent}
    });

    if (Pstream::master())
    {
        mkDir(outputDir);

        Info<< nl
            << "Transmission bath-bidomain manufactured-solution error "
            << "summary (t = " << t << ", solution " << solutionName_
            << "):" << nl
            << "field L1 L2 Linf" << nl;
        for (const auto& row : rows)
        {
            Info<< row.first() << ' ' << row.second().L1 << ' '
                << row.second().L2 << ' ' << row.second().Linf << nl;
        }
        Info<< endl;

        OFstream os(outputFile);
        os.precision(12);
        os  << "# Bath-bidomain manufactured solution error summary\n"
            << "# dimension " << dimensionName(2) << "\n"
            << "# cellsPerDirection " << nPerDirection_ << "\n"
            << "# h " << h << "\n"
            << "# steps " << nSteps << "\n"
            << "# dt " << dt << "\n"
            << "# time " << t << "\n"
            << "# solution " << solutionName_ << "\n"
            << "# field L1 L2 Linf\n";
        for (const auto& row : rows)
        {
            os  << row.first() << ' ' << row.second().L1 << ' '
                << row.second().L2 << ' ' << row.second().Linf << "\n";
        }
    }

    errorsReported_ = true;
}

} // End namespace Foam

// ************************************************************************* //
