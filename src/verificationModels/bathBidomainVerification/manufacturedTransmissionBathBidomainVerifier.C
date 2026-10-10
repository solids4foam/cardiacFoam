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
#include "OFstream.H"
#include "OSspecific.H"
#include "dimVoltage.H"
#include "fixedGradientFvPatchFields.H"
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


} // End anonymous namespace


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const dictionary&
manufacturedTransmissionBathBidomainVerifier::verificationDict() const
{
    return dict().subDict("verificationModel");
}


labelList manufacturedTransmissionBathBidomainVerifier::buildZones
(
    const fvMesh& mesh
) const
{
    labelList zone(mesh.nCells(), -1);

    label nMissingZones = 0;
    label nMultiplyZoned = 0;

    for (label z = 0; z < 3; ++z)
    {
        const label zoneI = mesh.cellZones().findZoneID(solution_->zoneName(z));
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
    label nMismatched = solution_->nMisplacedCells(mesh, zone);

    forAll(zone, cellI)
    {
        if (zone[cellI] < 0)
        {
            ++nUnzoned;
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
            << "manufactured zones: the cellZones " << solution_->zoneName(0)
            << ", " << solution_->zoneName(1) << " and "
            << solution_->zoneName(2)
            << " must exist and partition the cells. Missing zones on some "
            << "rank: " << nMissingZones << "; cells in several zones: "
            << nMultiplyZoned << "; cells in no zone: " << nUnzoned
            << "; cells placed outside their zone (" << solution_->geometry()
            << "): " << nMismatched << "."
            << exit(FatalError);
    }

    return zone;
}


wordList manufacturedTransmissionBathBidomainVerifier::outerPatches() const
{
    wordList names(dirichletPatches_);
    names.append(neumannPatches_);
    return names;
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
     || model.geometricDimension() != solution_->nSolutionD()
    )
    {
        FatalErrorInFunction
            << typeName << " requires the bidomainFDAManufactured or "
            << "bathBidomainFDAManufactured ionic model with geometric "
            << "dimension " << solution_->nSolutionD() << " for the "
            << solution_->geometry() << " geometry. Received ionic "
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

    if (phiEMesh.nSolutionD() != solution_->nSolutionD())
    {
        FatalErrorInFunction
            << typeName << " requires a " << solution_->nSolutionD()
            << "D extracellular-potential mesh for the "
            << solution_->geometry() << " geometry; mesh '"
            << phiEMesh.name() << "' solves in " << phiEMesh.nSolutionD()
            << " directions."
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
                tensor(solution_->Gi(c))
            )
        );
        maxGeErr = max
        (
            maxGeErr,
            relativeError
            (
                GeField.primitiveField()[cellI],
                tensor(solution_->Ge(c))
            )
        );
    }

    const vectorField& baseC = phiEMesh.C().primitiveField();
    scalar maxSigmaErr = 0;

    forAll(baseC, cellI)
    {
        const label z = zoneOfCell_[cellI];
        if (z != transmissionManufacturedSolution::heartZone)
        {
            maxSigmaErr = max
            (
                maxSigmaErr,
                relativeError
                (
                    sigmaTotal.primitiveField()[cellI],
                    solution_->sigmaB(baseC[cellI], z)*tensor::I
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

    solution_->validateInterface(phiEMesh, zoneOfCell_, heartMesh);
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
    solution_(),
    solutionName_(),
    nPerDirection_(0),
    dirichletPatches_(),
    neumannPatches_(),
    beta_(0),
    capacitanceRatio_(1),
    betaInitialised_(false),
    zoneOfCell_(),
    zoneMeshPtr_(nullptr)
{
    const dictionary& cfg = verificationDict();
    solution_ = transmissionManufacturedSolution::New(cfg);
    solutionName_ = cfg.get<word>("solution");

    nPerDirection_ = cfg.get<label>("N");
    if (nPerDirection_ <= 0)
    {
        FatalIOErrorInFunction(cfg)
            << "N must be positive; found " << nPerDirection_ << "."
            << exit(FatalIOError);
    }

    const dictionary& bc = cfg.subDict("boundaryData");
    dirichletPatches_ = bc.getOrDefault<wordList>("dirichlet", wordList());
    neumannPatches_ = bc.getOrDefault<wordList>("neumann", wordList());
    if (dirichletPatches_.empty() && neumannPatches_.empty())
    {
        FatalIOErrorInFunction(bc)
            << "boundaryData names no dirichlet and no neumann patch."
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

        zoneOfCell_ = buildZones(phiEPtr_->mesh());
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
        VmI[cellI] = solution_->Vm(C[cellI], t);
        u1I[cellI] = solution_->u1(C[cellI], t);
        u2I[cellI] = solution_->u2(C[cellI], t);
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
                zoned
              ? zoneOfCell_[cellI]
              : label(transmissionManufacturedSolution::heartZone);
            phiEI[cellI] = solution_->phiE(phiEC[cellI], t, z);
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
           *solution_->Vm(C[cellI], evaluationTime)
          - solution_->divGiGradPhiI(C[cellI], evaluationTime);
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
        r[cellI] +=
            solution_->phiESource(C[cellI], evaluationTime, zoneOfCell_[cellI]);
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
            values[faceI] = solution_->phiE
            (
                Cf[faceI],
                evaluationTime,
                solution_->zoneOf(Cf[faceI])
            );
        }

        pf == values;
    }

    // Neumann data: the exact current sigmaB(Cf) dn(phiE) over the owner conductivity the solver applies
    const volTensorField& sigmaTotal =
        mesh.lookupObject<volTensorField>("sigmaTotal");

    for (const word& patchName : neumannPatches_)
    {
        const label patchI = mesh.boundaryMesh().findPatchID(patchName);

        if (patchI < 0)
        {
            FatalErrorInFunction
                << "Neumann patch '" << patchName << "' is not a patch of "
                << "mesh '" << mesh.name() << "'."
                << exit(FatalError);
        }

        fvPatchScalarField& pf = phiE.boundaryFieldRef()[patchI];

        if (!isA<fixedGradientFvPatchScalarField>(pf))
        {
            FatalErrorInFunction
                << "Neumann patch '" << patchName << "' of "
                << phiE.name() << " is " << pf.type()
                << ", not fixedGradient."
                << exit(FatalError);
        }

        const fvPatch& fp = mesh.boundary()[patchI];
        const vectorField& Cf = fp.Cf();
        const vectorField nf(fp.nf());
        const tensorField sigmaOwner
        (
            sigmaTotal.boundaryField()[patchI].patchInternalField()
        );
        scalarField& gradient =
            refCast<fixedGradientFvPatchScalarField>(pf).gradient();

        forAll(Cf, faceI)
        {
            const label z = solution_->zoneOf(Cf[faceI]);
            gradient[faceI] =
                solution_->sigmaB(Cf[faceI], z)
               *(nf[faceI] & solution_->gradPhiE(Cf[faceI], evaluationTime, z))
               /(nf[faceI] & sigmaOwner[faceI] & nf[faceI]);
        }
    }
}


wordList manufacturedTransmissionBathBidomainVerifier::boundaryDataPatches()
const
{
    return outerPatches();
}


bool manufacturedTransmissionBathBidomainVerifier::
requiresPhiECompatibilityProjection() const
{
    return dirichletPatches_.empty();
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

    const labelList zone = buildZones(mesh);
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
        GiI[cellI] = solution_->Gi(C[cellI]);
        GeI[cellI] = solution_->Ge(C[cellI]);
        sigmaBI[cellI] =
            zone[cellI] == transmissionManufacturedSolution::heartZone
          ? 1.0
          : solution_->sigmaB(C[cellI], zone[cellI]);
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
        solution_->decomposition(mesh, zone)
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
    const label heartZone = transmissionManufacturedSolution::heartZone;

    const vectorField& C = mesh.C().primitiveField();
    const scalarField& V = mesh.V();
    const scalarField& phiEI = phiE.primitiveField();

    // Without a Dirichlet patch phiE is fixed up to a constant: remove the volume-weighted mean error
    scalar phiEOffset = 0;
    if (dirichletPatches_.empty())
    {
        scalar sumErr = 0;
        forAll(C, cellI)
        {
            sumErr +=
                V[cellI]
               *(phiEI[cellI] - solution_->phiE(C[cellI], t, zoneOfCell_[cellI]));
        }
        reduce(sumErr, sumOp<scalar>());
        phiEOffset = sumErr/gSum(V);
    }

    // Heart mesh: Vm, phiI, u1, u2

    const vectorField& heartC = heartMesh.C().primitiveField();
    const scalarField& VmI = Vm.primitiveField();

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
        VmExact[cellI] = solution_->Vm(heartC[cellI], t);
        phiIExact[cellI] = solution_->phiI(heartC[cellI], t);
        u1Exact[cellI] = solution_->u1(heartC[cellI], t);
        u2Exact[cellI] = solution_->u2(heartC[cellI], t);
        phiIValues[cellI] =
            VmI[cellI] + phiEI[heartCellMap[cellI]] - phiEOffset;
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
                gradVmI[cellI] - solution_->gradVm(heartC[cellI], t)
            );
        }
    }
    const errorNorms heartTraceNorms =
        weightedNorms(heartMesh.V(), heartTraceErr, heartTraceCell);

    // Base mesh: interface faces and bands

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
    DynamicList<scalar> faceErrBath;
    DynamicList<scalar> faceJump;
    DynamicList<scalar> faceArea;
    scalar currentComputed[2] = {0, 0};
    scalar currentExact[2] = {0, 0};
    const label nParts = solution_->nInterfaceParts();
    scalarField partComputed(2*nParts, 0);
    scalarField partExact(2*nParts, 0);

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
        const scalar exact = solution_->interfaceFlux(Cf, nOut, t)*magSf;
        const scalar exactBath =
            solution_->sigmaB(Cf, bathZone)
           *(nOut & solution_->gradPhiE(Cf, t, bathZone))*magSf;

        faceErr.append((sign*flux - exact)/magSf);
        faceErrBath.append((sign*flux - exactBath)/magSf);
        faceJump.append((exact - exactBath)/magSf);
        faceArea.append(magSf);

        const label j = (bathZone == 0 ? 0 : 1);
        currentComputed[j] += sign*flux;
        currentExact[j] += exact;

        if (nParts)
        {
            const label part = j*nParts + solution_->interfacePart(Cf);
            partComputed[part] += sign*flux;
            partExact[part] += exact;
        }
    };

    forAll(nei, faceI)
    {
        const label o = own[faceI];
        const label n = nei[faceI];
        const label zo = zoneOfCell_[o];
        const label zn = zoneOfCell_[n];

        if (zo == zn || (zo != heartZone && zn != heartZone))
        {
            continue;
        }

        interfaceCell[o] = true;
        interfaceCell[n] = true;

        const bool ownerHeart = (zo == heartZone);
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

            if (zLocal == zNbr || (zLocal != heartZone && zNbr != heartZone))
            {
                continue;
            }

            interfaceCell[cellI] = true;

            if (zLocal == heartZone)
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

    const boolList allFaces(faceErr.size(), true);
    const errorNorms interfaceFluxNorms =
        weightedNorms(scalarField(faceArea), scalarField(faceErr), allFaces);
    const errorNorms interfaceFluxBathNorms =
        weightedNorms(scalarField(faceArea), scalarField(faceErrBath), allFaces);
    const errorNorms interfaceFluxJumpNorms =
        weightedNorms(scalarField(faceArea), scalarField(faceJump), allFaces);

    reduce(currentComputed[0], sumOp<scalar>());
    reduce(currentComputed[1], sumOp<scalar>());
    reduce(currentExact[0], sumOp<scalar>());
    reduce(currentExact[1], sumOp<scalar>());

    const errorNorms current1 =
        singleValue(mag(currentComputed[0] - currentExact[0]));
    const errorNorms current2 =
        singleValue(mag(currentComputed[1] - currentExact[1]));

    // Partial currents over the interface parts: L1 and L2 are means over the parts, Linf the largest
    Pstream::listCombineReduce(partComputed, plusEqOp<scalar>());
    Pstream::listCombineReduce(partExact, plusEqOp<scalar>());
    errorNorms partNorms[2];
    for (label j = 0; j < 2; ++j)
    {
        scalar sumAbs = 0;
        scalar sumSq = 0;
        scalar maxAbs = 0;
        for (label part = 0; part < nParts; ++part)
        {
            const scalar d =
                mag(partComputed[j*nParts + part] - partExact[j*nParts + part]);
            sumAbs += d;
            sumSq += d*d;
            maxAbs = max(maxAbs, d);
        }
        const label n = max(nParts, label(1));
        partNorms[j] = {sumAbs/n, Foam::sqrt(sumSq/n), maxAbs};
    }

    // Outer patches: band cells, and the current through the Dirichlet patches

    const volTensorField& sigmaTotal =
        mesh.lookupObject<volTensorField>("sigmaTotal");
    scalar outerComputed = 0;
    scalar outerExact = 0;

    for (const word& patchName : outerPatches())
    {
        const label patchI = mesh.boundaryMesh().findPatchID(patchName);
        if (patchI < 0)
        {
            FatalErrorInFunction
                << "Outer patch '" << patchName << "' is not a patch of "
                << "mesh '" << mesh.name() << "'."
                << exit(FatalError);
        }

        const fvPatch& fp = mesh.boundary()[patchI];
        for (const label cellI : fp.faceCells())
        {
            boundaryCell[cellI] = true;
        }

        if (!dirichletPatches_.found(patchName))
        {
            continue;
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

            const label z = solution_->zoneOf(Cf[i]);
            outerExact +=
                solution_->sigmaB(Cf[i], z)
               *(nf[i] & solution_->gradPhiE(Cf[i], t, z))
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
    const vectorField& gradJumpAwareI =
        mesh.lookupObject<volVectorField>("phiEInterfaceGradient")
       .primitiveField();

    scalarField phiEErr(C.size());
    scalarField gradTangentialErr(C.size(), 0);
    scalarField gradJumpAwareTangentialErr(C.size(), 0);
    boolList heartCell(C.size(), false);
    boolList bathCell(C.size(), false);

    forAll(C, cellI)
    {
        const label z = zoneOfCell_[cellI];
        phiEErr[cellI] =
            phiEI[cellI] - phiEOffset - solution_->phiE(C[cellI], t, z);
        heartCell[cellI] = (z == heartZone);
        bathCell[cellI] = !heartCell[cellI];

        if (interfaceCell[cellI])
        {
            const vector gradExact(solution_->gradPhiE(C[cellI], t, z));
            const tensor tangential
            (
                tensor::I - sqr(solution_->interfaceNormal(C[cellI]))
            );
            gradTangentialErr[cellI] =
                mag(tangential & (gradPhiEI[cellI] - gradExact));
            gradJumpAwareTangentialErr[cellI] =
                mag(tangential & (gradJumpAwareI[cellI] - gradExact));
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
    const errorNorms gradJumpAwareTangentialNorms =
        weightedNorms(V, gradJumpAwareTangentialErr, interfaceCell);

    // Header data

    const label nD = solution_->nSolutionD();
    const label nHeartCells = returnReduce(heartMesh.nCells(), sumOp<label>());
    const scalar heartMeasure =
        nD == 2
      ? gSum(heartMesh.V())/mesh.bounds().span().z()
      : gSum(heartMesh.V());
    const scalar cellMeasure = heartMeasure/max(nHeartCells, label(1));
    const scalar h =
        nD == 2 ? Foam::sqrt(cellMeasure) : Foam::cbrt(cellMeasure);
    const scalar dt = time.deltaTValue();
    const label nSteps = max(label(0), time.timeIndex());

    const fileName outputDir(time.globalPath()/"postProcessing");
    const fileName outputFile =
        outputDir
      / (dimensionName(nD) + "_" + Foam::name(nPerDirection_) + "_cells.dat");

    List<Tuple2<word, errorNorms>> rows
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
        {"interfaceFlux_bathSide", interfaceFluxBathNorms},
        {"interfaceFluxJumpScale", interfaceFluxJumpNorms},
        {"interfaceCurrent_1", current1},
        {"interfaceCurrent_2", current2},
        {"gradTangentialError", gradTangentialNorms},
        {"gradJumpAwareTangentialError", gradJumpAwareTangentialNorms},
        {"gradHeartTraceError", heartTraceNorms},
        {"outerCurrent", outerCurrent}
    });
    if (nParts)
    {
        rows.append({"interfacePartCurrent_1", partNorms[0]});
        rows.append({"interfacePartCurrent_2", partNorms[1]});
    }

    if (Pstream::master())
    {
        mkDir(outputDir);

        Info<< nl
            << "Transmission bath-bidomain manufactured-solution error "
            << "summary (t = " << t << ", " << solution_->geometry()
            << ", solution " << solutionName_ << "):" << nl
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
            << "# dimension " << dimensionName(nD) << "\n"
            << "# geometry " << solution_->geometry() << "\n"
            << "# cellsPerDirection " << nPerDirection_ << "\n"
            << "# h " << h << "\n"
            << "# steps " << nSteps << "\n"
            << "# dt " << dt << "\n"
            << "# time " << t << "\n"
            << "# solution " << solutionName_ << "\n"
            << "# phiEGaugeOffset " << phiEOffset << "\n"
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
