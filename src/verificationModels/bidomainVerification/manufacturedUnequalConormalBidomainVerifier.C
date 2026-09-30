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

#include "manufacturedUnequalConormalBidomainVerifier.H"

#include "OFstream.H"
#include "OSspecific.H"
#include "PstreamReduceOps.H"
#include "boundBox.H"
#include "dimensionedSymmTensor.H"
#include "ionicModel.H"
#include "mathematicalConstants.H"
#include "verificationUtils.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
using namespace verificationUtils;

defineTypeNameAndDebug(manufacturedUnequalConormalBidomainVerifier, 0);
addToRunTimeSelectionTable
(
    electroVerificationModel,
    manufacturedUnequalConormalBidomainVerifier,
    dictionary
);

manufacturedUnequalConormalBidomainVerifier::manufacturedUnequalConormalBidomainVerifier
(
    const dictionary& dict
)
:
    electroVerificationModel(dict),
    phiEPtr_(nullptr),
    errorsReported_(false),
    phiEReferenceValue_(0.0),
    phiEReferencePoint_(point::zero),
    Gi_(tensor::zero),
    Ge_(tensor::zero),
    profile_(),
    beta_(0.0),
    betaInitialised_(false),
    conductivityValidated_(false)
{
    const dictionary& coeffDict = this->dict();
    phiEReferenceValue_ = coeffDict.lookupOrDefault<scalar>("phiEReferenceValue", 0.0);
    phiEReferencePoint_ = coeffDict.get<point>("phiERefPoint");

    Gi_ = tensor(dimensionedSymmTensor("conductivityIntracellular", coeffDict).value());
    Ge_ = tensor(dimensionedSymmTensor("conductivityExtracellular", coeffDict).value());

    profile_ = unequalConormalWallProfile
    (
        Gi_,
        Ge_,
        2.0*constant::mathematical::pi
    );
}


void manufacturedUnequalConormalBidomainVerifier::bindBidomainField(volScalarField& phiE)
{
    phiEPtr_ = &phiE;
}


void manufacturedUnequalConormalBidomainVerifier::initialiseBeta(const ionicModel& model)
{
    const wordList constantNames = model.constantVariableNames();
    const scalarField constantValues = model.constantVariableValues();

    if (constantNames.size() != constantValues.size())
    {
        FatalErrorInFunction
            << "Ionic model " << model.type()
            << " exposes inconsistent constant metadata: "
            << constantNames.size() << " names but "
            << constantValues.size() << " values."
            << exit(FatalError);
    }

    label betaIndex = -1;
    forAll(constantNames, i)
    {
        if (constantNames[i] == "Beta")
        {
            betaIndex = i;
            break;
        }
    }

    if (betaIndex < 0)
    {
        FatalErrorInFunction
            << "manufacturedUnequalConormalBidomainVerifier requires the "
            << "manufactured ionic constant 'Beta', but ionic model "
            << model.type() << " exposes constants " << constantNames << "."
            << exit(FatalError);
    }

    beta_ = constantValues[betaIndex];
    betaInitialised_ = true;
}


void manufacturedUnequalConormalBidomainVerifier::validateConfiguration
(
    const ionicModel& model,
    const fvMesh& mesh
) const
{
    if
    (
        model.verificationFamily() != "bidomainFDAManufactured"
     || model.geometricDimension() != 3
    )
    {
        FatalErrorInFunction
            << "manufacturedUnequalConormalBidomainVerifier requires the 3D "
            << "bidomainFDAManufactured ionic model. Received ionic model "
            << model.type() << " with verification family '"
            << model.verificationFamily() << "' and geometric dimension "
            << model.geometricDimension() << "."
            << exit(FatalError);
    }

    const boundBox bounds(mesh.points(), true);
    const vector unitMin(vector::zero);
    const vector unitMax(vector::one);
    const scalar cubeTolerance = 1e-10;

    if
    (
        mag(bounds.min() - unitMin) > cubeTolerance
     || mag(bounds.max() - unitMax) > cubeTolerance
    )
    {
        FatalErrorInFunction
            << "manufacturedUnequalConormalBidomainVerifier is defined on the "
            << "unit cube [0,1]^3. Mesh bounds are " << bounds << "."
            << exit(FatalError);
    }
}


void manufacturedUnequalConormalBidomainVerifier::validateConductivity
(
    const fvMesh& mesh
) const
{
    const wordList names
    ({
        "conductivityIntracellular",
        "conductivityExtracellular"
    });
    const List<tensor> configured({Gi_, Ge_});

    forAll(names, i)
    {
        const tensor meanSigma = uniformTensor
        (
            mesh.lookupObject<volTensorField>(names[i]).primitiveField(),
            typeName
        );
        const scalar scale = max(scalar(1), mag(meanSigma));

        if (mag(meanSigma - configured[i]) > 1e-10*scale)
        {
            FatalErrorInFunction
                << typeName << " requires the solver's " << names[i]
                << " to match the configured tensor. Solver mean tensor is "
                << meanSigma << ", configuration tensor is " << configured[i]
                << "."
                << exit(FatalError);
        }
    }
}


wordList manufacturedUnequalConormalBidomainVerifier::preProcessFieldNames
(
    const ionicModel&
) const
{
    return wordList({"u1", "u2", "u3"});
}


wordList manufacturedUnequalConormalBidomainVerifier::requiredPostProcessFieldNames
(
    const ionicModel&
) const
{
    return wordList({"u1", "u2"});
}


bool manufacturedUnequalConormalBidomainVerifier::shouldPostProcess
(
    const ionicModel&,
    const volScalarField& Vm
) const
{
    return !errorsReported_ && shouldReportManufacturedErrors(Vm);
}


void manufacturedUnequalConormalBidomainVerifier::preProcess
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
            << "manufacturedUnequalConormalBidomainVerifier preProcess expected "
            << names.size() << " preProcess fields " << names
            << " but received " << fields.size()
            << exit(FatalError);
    }

    validateConfiguration(model, Vm.mesh());
    initialiseBeta(model);

    const label u1Field = requireFieldIndex(names, "u1", "preProcess");
    const label u2Field = requireFieldIndex(names, "u2", "preProcess");
    const label u3Field = requireFieldIndex(names, "u3", "preProcess");

    volScalarField& u1m = fields[u1Field];
    volScalarField& u2m = fields[u2Field];
    volScalarField& u3m = fields[u3Field];

    const vectorField& centres = Vm.mesh().C().primitiveField();
    const scalar t = Vm.mesh().time().value();

    scalarField G;
    computeUnequalConormalG(G, centres);

    scalarField& VmI = Vm.primitiveFieldRef();
    computeUnequalConormalV(VmI, profile_, centres, t);

    scalarField& u1I = u1m.primitiveFieldRef();
    scalarField& u2I = u2m.primitiveFieldRef();
    scalarField& u3I = u3m.primitiveFieldRef();
    computeUnequalConormalU(u1I, u2I, u3I, G, VmI, t);

    Vm.correctBoundaryConditions();
    u1m.correctBoundaryConditions();
    u2m.correctBoundaryConditions();
    u3m.correctBoundaryConditions();

    if (phiEPtr_)
    {
        const scalar shift = manufacturedUnequalConormalBidomainShift
        (
            profile_,
            t,
            phiEReferencePoint_,
            phiEReferenceValue_
        );

        scalarField& phiEI = phiEPtr_->primitiveFieldRef();
        computeUnequalConormalPhiE(phiEI, profile_, centres, t, shift);
        phiEPtr_->correctBoundaryConditions();
    }

    model.importFields(Vm, names, fields);
}


void manufacturedUnequalConormalBidomainVerifier::addManufacturedPdeSource
(
    volScalarField& sourceField,
    const volTensorField& conductivity,
    const scalar evaluationTime
) const
{
    if (!betaInitialised_)
    {
        FatalErrorInFunction
            << "The manufactured source was requested before ionic constant "
            << "Beta was initialised by preProcess."
            << exit(FatalError);
    }

    if (&sourceField.mesh() != &conductivity.mesh())
    {
        FatalErrorInFunction
            << "Source field and conductivity field belong to different "
            << "finite-volume meshes."
            << exit(FatalError);
    }

    if (!conductivityValidated_)
    {
        validateConductivity(sourceField.mesh());
        conductivityValidated_ = true;
    }

    scalarField source;
    computeUnequalConormalPdeSource
    (
        source,
        profile_,
        sourceField.mesh().C().primitiveField(),
        evaluationTime,
        beta_
    );

    sourceField.primitiveFieldRef() += source;
    sourceField.correctBoundaryConditions();
}


void manufacturedUnequalConormalBidomainVerifier::postProcess
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

    if (!phiEPtr_)
    {
        FatalErrorInFunction
            << "manufacturedUnequalConormalBidomainVerifier requires a bound phiE "
            << "field before postProcess()."
            << exit(FatalError);
    }

    const wordList requiredNames = requiredPostProcessFieldNames(model);
    if (fields.size() != requiredNames.size())
    {
        FatalErrorInFunction
            << "manufacturedUnequalConormalBidomainVerifier expected "
            << requiredNames.size() << " postProcess fields " << requiredNames
            << " but received " << fields.size()
            << exit(FatalError);
    }

    const label u1Field =
        requireFieldIndex(requiredNames, "u1", "postProcess");
    const label u2Field =
        requireFieldIndex(requiredNames, "u2", "postProcess");

    const fvMesh& mesh = Vm.mesh();
    const Time& time = mesh.time();
    const vectorField& centres = mesh.C().primitiveField();
    const scalar t = time.value();

    scalarField G;
    computeUnequalConormalG(G, centres);

    const scalar shift = manufacturedUnequalConormalBidomainShift
    (
        profile_,
        t,
        phiEReferencePoint_,
        phiEReferenceValue_
    );

    scalarField VmExact, phiEExact, phiIExact, u1Exact, u2Exact, u3Exact;
    computeUnequalConormalV(VmExact, profile_, centres, t);
    computeUnequalConormalPhiE(phiEExact, profile_, centres, t, shift);
    computeUnequalConormalPhiI(phiIExact, VmExact, phiEExact);
    computeUnequalConormalU(u1Exact, u2Exact, u3Exact, G, VmExact, t);

    const scalarField& VmValues = Vm.primitiveField();
    const scalarField& phiEValues = phiEPtr_->primitiveField();
    const scalarField& u1Values = fields[u1Field].primitiveField();
    const scalarField& u2Values = fields[u2Field].primitiveField();

    scalarField phiIValues(VmValues.size(), 0.0);
    forAll(phiIValues, i)
    {
        phiIValues[i] = VmValues[i] + phiEValues[i];
    }

    const scalar phiEGaugeOffset =
        computeVolumeMeanError(mesh, phiEValues, phiEExact);

    scalarField phiEExactGauge(phiEExact);
    scalarField phiIExactGauge(phiIExact);
    phiEExactGauge += phiEGaugeOffset;
    phiIExactGauge += phiEGaugeOffset;

    const auto VmNorms = computeNorms(mesh, VmValues, VmExact);
    const auto phiENorms = computeNorms(mesh, phiEValues, phiEExact);
    const auto phiEGaugeNorms =
        computeNorms(mesh, phiEValues, phiEExactGauge);
    const auto phiINorms = computeNorms(mesh, phiIValues, phiIExact);
    const auto phiIGaugeNorms =
        computeNorms(mesh, phiIValues, phiIExactGauge);
    const auto u1Norms = computeNorms(mesh, u1Values, u1Exact);
    const auto u2Norms = computeNorms(mesh, u2Values, u2Exact);

    const label totalCells = globalManufacturedCellCount(mesh);
    const label nPerDirection = structuredCellsPerDirection(totalCells, 3);
    const scalar dx = structuredManufacturedDx(nPerDirection);
    const scalar dt = time.deltaTValue();
    const label nSteps = max(label(0), time.timeIndex());

    const fileName outputDir(time.globalPath()/"postProcessing");
    mkDir(outputDir);
    const fileName outputFile =
        outputDir
      / (
            dimensionName(3) + "_" + Foam::name(nPerDirection) + "_cells.dat"
        );

    if (Pstream::master())
    {
        Info<< nl
            << "Bidomain manufactured-solution error summary (t = " << t
            << "):" << nl
            << "-------------------------------------------------" << nl
            << "Field     L1-error       L2-error       Linf-error" << nl
            << "Vm        " << VmNorms.first().first() << "   "
            << VmNorms.first().second() << "   " << VmNorms.second() << nl
            << "phiE      " << phiENorms.first().first() << "   "
            << phiENorms.first().second() << "   " << phiENorms.second() << nl
            << "phiE_gauge " << phiEGaugeNorms.first().first() << "   "
            << phiEGaugeNorms.first().second() << "   "
            << phiEGaugeNorms.second() << nl
            << "phiI      " << phiINorms.first().first() << "   "
            << phiINorms.first().second() << "   " << phiINorms.second() << nl
            << "phiI_gauge " << phiIGaugeNorms.first().first() << "   "
            << phiIGaugeNorms.first().second() << "   "
            << phiIGaugeNorms.second() << nl
            << "u1        " << u1Norms.first().first() << "   "
            << u1Norms.first().second() << "   " << u1Norms.second() << nl
            << "u2        " << u2Norms.first().first() << "   "
            << u2Norms.first().second() << "   " << u2Norms.second() << nl
            << "-------------------------------------------------" << endl;

        OFstream out(outputFile);
        out << "Bidomain manufactured-solution error summary (t = "
            << t << "):\n";
        out << "Field     L1-error       L2-error       Linf-error\n";
        out << "Vm        " << VmNorms.first().first() << "   "
            << VmNorms.first().second() << "   " << VmNorms.second() << "\n";
        out << "phiE      " << phiENorms.first().first() << "   "
            << phiENorms.first().second() << "   " << phiENorms.second() << "\n";
        out << "phiE_gauge " << phiEGaugeNorms.first().first() << "   "
            << phiEGaugeNorms.first().second() << "   "
            << phiEGaugeNorms.second() << "\n";
        out << "phiI      " << phiINorms.first().first() << "   "
            << phiINorms.first().second() << "   " << phiINorms.second() << "\n";
        out << "phiI_gauge " << phiIGaugeNorms.first().first() << "   "
            << phiIGaugeNorms.first().second() << "   "
            << phiIGaugeNorms.second() << "\n";
        out << "u1        " << u1Norms.first().first() << "   "
            << u1Norms.first().second() << "   " << u1Norms.second() << "\n";
        out << "u2        " << u2Norms.first().first() << "   "
            << u2Norms.first().second() << "   " << u2Norms.second() << "\n\n";

        out << "Number of cells (N)   = " << nPerDirection << "\n";

        out << "Grid spacing (dx)     = " << dx << "\n";
        out << "Time step (dt)        = " << dt << "\n";
        out << "Number of steps       = " << nSteps << "\n";
        out << "Final simulation time = " << t << "\n";
        out << "Beta                  = " << beta_ << "\n";
        out << "phiE reference point  = " << phiEReferencePoint_ << "\n";
        out << "phiE reference value  = " << phiEReferenceValue_ << "\n";
        out << "phiE gauge offset     = " << phiEGaugeOffset << "\n";
    }

    errorsReported_ = true;
}

} // End namespace Foam

// ************************************************************************* //
