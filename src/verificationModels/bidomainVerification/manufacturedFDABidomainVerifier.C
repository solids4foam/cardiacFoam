/*---------------------------------------------------------------------------*\
License
    This file is part of cardiacFoam.

    cardiacFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    cardiacFoam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
    or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License along
    with cardiacFoam.  If not, see <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "manufacturedFDABidomainVerifier.H"

#include "OFstream.H"
#include "OSspecific.H"
#include "PstreamReduceOps.H"
#include "bidomainVerification/manufacturedFDABidomainReference.H"
#include "ionicModel.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug(manufacturedFDABidomainVerifier, 0);
addToRunTimeSelectionTable
(
    electroVerificationModel,
    manufacturedFDABidomainVerifier,
    dictionary
);

namespace
{

label requireFieldIndex(const wordList& names, const word& name, const char* phase)
{
    forAll(names, i)
    {
        if (names[i] == name)
        {
            return i;
        }
    }

    FatalErrorInFunction
        << "Required " << phase << " field '" << name
        << "' is missing from configured hook fields " << names
        << exit(FatalError);

    return -1;
}


bool shouldReportManufacturedErrors(const volScalarField& Vm)
{
    const Time& time = Vm.mesh().time();
    const scalar t = time.value();
    const scalar dt = time.deltaTValue();
    const scalar endTime = time.endTime().value();

    return t + 0.5*dt >= endTime;
}


label globalManufacturedCellCount(const fvMesh& mesh)
{
    label totalCells = mesh.nCells();
    reduce(totalCells, sumOp<label>());
    return totalCells;
}


label structuredCellsPerDirection(const label totalCells, const label dimension)
{
    if (totalCells <= 0 || dimension <= 0)
    {
        return 0;
    }

    if (dimension == 1)
    {
        return totalCells;
    }

    return max
    (
        label(1),
        label(Foam::pow(scalar(totalCells), 1.0/scalar(dimension)) + 0.5)
    );
}


scalar structuredManufacturedDx(const label nPerDirection)
{
    if (nPerDirection <= 0)
    {
        return 0.0;
    }

    return 1.0/scalar(nPerDirection);
}


Tuple2<Tuple2<scalar, scalar>, scalar> computeNorms
(
    const scalarField& numeric,
    const scalarField& exact
)
{
    scalar sumAbs = 0.0;
    scalar sumSq = 0.0;
    scalar maxAbs = 0.0;

    forAll(numeric, i)
    {
        const scalar diff = Foam::mag(numeric[i] - exact[i]);
        sumAbs += diff;
        sumSq += diff*diff;
        maxAbs = max(maxAbs, diff);
    }

    reduce(sumAbs, sumOp<scalar>());
    reduce(sumSq, sumOp<scalar>());
    reduce(maxAbs, maxOp<scalar>());

    label n = numeric.size();
    reduce(n, sumOp<label>());

    return Tuple2<Tuple2<scalar, scalar>, scalar>
    (
        Tuple2<scalar, scalar>(sumAbs/scalar(n), Foam::sqrt(sumSq/scalar(n))),
        maxAbs
    );
}


word dimensionName(const label dimension)
{
    return
        dimension == 1 ? "1D"
      : dimension == 2 ? "2D"
      : dimension == 3 ? "3D"
                       : "unknown";
}

} // End anonymous namespace


manufacturedFDABidomainVerifier::manufacturedFDABidomainVerifier
(
    const dictionary& dict
)
:
    electroVerificationModel(dict),
    phiEPtr_(nullptr),
    enabled_(true),
    useExplicitAlgorithm_(false),
    errorsReported_(false),
    k_(0.5),
    phiEReferenceCell_(0),
    phiEReferenceValue_(0.0),
    outputFileName_()
{
    const dictionary& cfg = verificationDict();

    enabled_ = cfg.lookupOrDefault<Switch>("enabled", true);
    useExplicitAlgorithm_ =
        cfg.lookupOrDefault<word>("solutionAlgorithm", "implicit") == "explicit";
    k_ = cfg.lookupOrDefault<scalar>("k", 0.5);
    phiEReferenceCell_ = cfg.lookupOrDefault<label>("phiEReferenceCell", 0);
    phiEReferenceValue_ = cfg.lookupOrDefault<scalar>("phiEReferenceValue", 0.0);
    outputFileName_ = cfg.lookupOrDefault<fileName>("outputFile", fileName());
}


const dictionary& manufacturedFDABidomainVerifier::verificationDict() const
{
    if (dict().found("manufacturedFDABidomainVerifierCoeffs"))
    {
        return dict().subDict("manufacturedFDABidomainVerifierCoeffs");
    }

    return dict();
}


wordList manufacturedFDABidomainVerifier::readConfiguredNames
(
    const word& subDictName,
    const word& entryName,
    const wordList& defaults
) const
{
    const dictionary& cfg = verificationDict();

    if (!cfg.found(subDictName))
    {
        return defaults;
    }

    const dictionary& hookDict = cfg.subDict(subDictName);
    return hookDict.lookupOrDefault<wordList>(entryName, defaults);
}


void manufacturedFDABidomainVerifier::bindBidomainField(volScalarField& phiE)
{
    phiEPtr_ = &phiE;
}


wordList manufacturedFDABidomainVerifier::preProcessFieldNames
(
    const ionicModel&
) const
{
    return readConfiguredNames
    (
        "solverHookFields",
        "preProcess",
        wordList({"u1", "u2", "u3"})
    );
}


wordList manufacturedFDABidomainVerifier::requiredPostProcessFieldNames
(
    const ionicModel&
) const
{
    return readConfiguredNames
    (
        "solverHookFields",
        "postProcess",
        wordList({"u1", "u2"})
    );
}


bool manufacturedFDABidomainVerifier::shouldPostProcess
(
    const ionicModel&,
    const volScalarField& Vm
) const
{
    return enabled_ && !errorsReported_ && shouldReportManufacturedErrors(Vm);
}


void manufacturedFDABidomainVerifier::preProcess
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
            << "manufacturedFDABidomainVerifier preProcess expected "
            << names.size() << " preProcess fields " << names
            << " but received " << fields.size()
            << exit(FatalError);
    }

    const label u1Field = requireFieldIndex(names, "u1", "preProcess");
    const label u2Field = requireFieldIndex(names, "u2", "preProcess");
    const label u3Field = requireFieldIndex(names, "u3", "preProcess");

    volScalarField& u1m = fields[u1Field];
    volScalarField& u2m = fields[u2Field];
    volScalarField& u3m = fields[u3Field];

    const vectorField& centres = Vm.mesh().C().primitiveField();
    scalarField X(centres.component(vector::X));
    scalarField Y(centres.component(vector::Y));
    scalarField Z(centres.component(vector::Z));
    const scalar t = Vm.mesh().time().value();
    const label dimension = model.geometricDimension();

    scalarField& VmI = Vm.primitiveFieldRef();
    scalarField& u1I = u1m.primitiveFieldRef();
    scalarField& u2I = u2m.primitiveFieldRef();
    scalarField& u3I = u3m.primitiveFieldRef();

    computeManufacturedV(VmI, X, Y, Z, t, dimension);
    computeManufacturedU(u1I, u2I, u3I, X, Y, Z, t, dimension);

    Vm.correctBoundaryConditions();
    u1m.correctBoundaryConditions();
    u2m.correctBoundaryConditions();
    u3m.correctBoundaryConditions();

    if (phiEPtr_)
    {
        scalarField& phiEI = phiEPtr_->primitiveFieldRef();

        label referenceCell = phiEReferenceCell_;
        if (referenceCell < 0 || referenceCell >= centres.size())
        {
            referenceCell = 0;
        }

        const scalar shift = manufacturedFDABidomainShift
        (
            t,
            dimension,
            k_,
            centres[referenceCell],
            phiEReferenceValue_
        );

        computeManufacturedBidomainPhiE(phiEI, X, Y, Z, t, dimension, k_, shift);
        phiEPtr_->correctBoundaryConditions();
    }

    model.importFields(Vm, names, fields);
}


void manufacturedFDABidomainVerifier::postProcess
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
            << "manufacturedFDABidomainVerifier requires a bound phiE field "
            << "before postProcess()."
            << exit(FatalError);
    }

    const wordList requiredNames = requiredPostProcessFieldNames(model);
    if (fields.size() != requiredNames.size())
    {
        FatalErrorInFunction
            << "manufacturedFDABidomainVerifier expected "
            << requiredNames.size() << " postProcess fields " << requiredNames
            << " but received " << fields.size()
            << exit(FatalError);
    }

    const label dimension = model.geometricDimension();
    if (dimension < 1 || dimension > 3)
    {
        FatalErrorInFunction
            << "manufacturedFDABidomainVerifier requires a valid geometric "
            << "dimension in [1,3], but ionic model '" << model.type()
            << "' reported " << dimension << "."
            << exit(FatalError);
    }

    const label u1Field =
        requireFieldIndex(requiredNames, "u1", "postProcess");
    const label u2Field =
        requireFieldIndex(requiredNames, "u2", "postProcess");

    const fvMesh& mesh = Vm.mesh();
    const Time& time = mesh.time();
    const vectorField& centres = mesh.C().primitiveField();
    scalarField X(centres.component(vector::X));
    scalarField Y(centres.component(vector::Y));
    scalarField Z(centres.component(vector::Z));

    const scalar t = time.value();

    label referenceCell = phiEReferenceCell_;
    if (referenceCell < 0 || referenceCell >= centres.size())
    {
        referenceCell = 0;
    }

    const scalar shift = manufacturedFDABidomainShift
    (
        t,
        dimension,
        k_,
        centres[referenceCell],
        phiEReferenceValue_
    );

    scalarField VmExact;
    scalarField phiEExact;
    scalarField phiIExact;
    scalarField u1Exact;
    scalarField u2Exact;
    scalarField u3Exact;

    computeManufacturedV(VmExact, X, Y, Z, t, dimension);
    computeManufacturedBidomainPhiE(phiEExact, X, Y, Z, t, dimension, k_, shift);
    computeManufacturedBidomainPhiI(phiIExact, X, Y, Z, t, dimension, k_, shift);
    computeManufacturedU(u1Exact, u2Exact, u3Exact, X, Y, Z, t, dimension);

    const scalarField& VmValues = Vm.primitiveField();
    const scalarField& phiEValues = phiEPtr_->primitiveField();
    const scalarField& u1Values = fields[u1Field].primitiveField();
    const scalarField& u2Values = fields[u2Field].primitiveField();

    scalarField phiIValues(VmValues.size(), 0.0);
    forAll(phiIValues, i)
    {
        phiIValues[i] = VmValues[i] + phiEValues[i];
    }

    const auto VmNorms = computeNorms(VmValues, VmExact);
    const auto phiENorms = computeNorms(phiEValues, phiEExact);
    const auto phiINorms = computeNorms(phiIValues, phiIExact);
    const auto u1Norms = computeNorms(u1Values, u1Exact);
    const auto u2Norms = computeNorms(u2Values, u2Exact);

    const label totalCells = globalManufacturedCellCount(mesh);
    const label nPerDirection = structuredCellsPerDirection(totalCells, dimension);
    const scalar dx = structuredManufacturedDx(nPerDirection);
    const scalar dt = time.deltaTValue();
    const label nSteps = max(label(0), time.timeIndex());

    fileName outputFile = outputFileName_;
    if (outputFile.empty())
    {
        const fileName outputDir(time.path()/"postProcessing");
        mkDir(outputDir);
        outputFile =
            outputDir
          / (
                "bidomain_"
              + dimensionName(dimension)
              + "_"
              + Foam::name(nPerDirection)
              + "_cells_"
              + word(useExplicitAlgorithm_ ? "explicit" : "implicit")
              + ".dat"
            );
    }

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
            << "phiI      " << phiINorms.first().first() << "   "
            << phiINorms.first().second() << "   " << phiINorms.second() << nl
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
        out << "phiI      " << phiINorms.first().first() << "   "
            << phiINorms.first().second() << "   " << phiINorms.second() << "\n";
        out << "u1        " << u1Norms.first().first() << "   "
            << u1Norms.first().second() << "   " << u1Norms.second() << "\n";
        out << "u2        " << u2Norms.first().first() << "   "
            << u2Norms.first().second() << "   " << u2Norms.second() << "\n\n";
        out << "Number of cells (N)   = " << nPerDirection << "\n";
        out << "Solver type           = "
            << (useExplicitAlgorithm_ ? "Explicit" : "Implicit") << "\n";
        out << "Grid spacing (dx)     = " << dx << "\n";
        out << "Time step (dt)        = " << dt << "\n";
        out << "Number of steps       = " << nSteps << "\n";
        out << "Final simulation time = " << t << "\n";
        out << "k                     = " << k_ << "\n";
        out << "phiE reference cell   = " << phiEReferenceCell_ << "\n";
        out << "phiE reference value  = " << phiEReferenceValue_ << "\n";
    }

    errorsReported_ = true;
}

} // End namespace Foam

// ************************************************************************* //
