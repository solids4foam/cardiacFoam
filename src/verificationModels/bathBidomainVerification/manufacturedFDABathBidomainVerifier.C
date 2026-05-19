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

#include "manufacturedFDABathBidomainVerifier.H"

#include "OFstream.H"
#include "OSspecific.H"
#include "bathBidomainVerification/manufacturedFDABathBidomainReference.H"
#include "ionicModel.H"
#include "addToRunTimeSelectionTable.H"

#include <limits>

namespace Foam
{

defineTypeNameAndDebug(manufacturedFDABathBidomainVerifier, 0);
addToRunTimeSelectionTable
(
    electroVerificationModel,
    manufacturedFDABathBidomainVerifier,
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


Tuple2<Tuple2<scalar, scalar>, scalar> computeNorms
(
    const scalarField& numeric,
    const scalarField& exact
)
{
    if (numeric.size() != exact.size())
    {
        FatalErrorInFunction
            << "Cannot compute manufactured norms for fields with different "
            << "sizes: numeric=" << numeric.size()
            << ", exact=" << exact.size() << "."
            << exit(FatalError);
    }

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

    const scalar denom = max(scalar(1), scalar(n));

    return Tuple2<Tuple2<scalar, scalar>, scalar>
    (
        Tuple2<scalar, scalar>(sumAbs/denom, Foam::sqrt(sumSq/denom)),
        maxAbs
    );
}


Tuple2<Tuple2<scalar, scalar>, scalar> nanNorms()
{
    const scalar nan = std::numeric_limits<scalar>::quiet_NaN();
    return Tuple2<Tuple2<scalar, scalar>, scalar>
    (
        Tuple2<scalar, scalar>(nan, nan),
        nan
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


// Return the verificationModel sub-dictionary, which holds all
// bath-bidomain verifier parameters (k, alpha, groundElectrode,
// enabled, outputFile).
const dictionary&
manufacturedFDABathBidomainVerifier::verificationDict() const
{
    return dict().subDict("verificationModel");
}


manufacturedFDABathBidomainVerifier::manufacturedFDABathBidomainVerifier
(
    const dictionary& dict
)
:
    electroVerificationModel(dict),
    phiEPtr_(nullptr),
    phiEHeartCellMapPtr_(nullptr),
    enabled_(true),
    useExplicitAlgorithm_(false),
    errorsReported_(false),
    k_(1.0/Foam::sqrt(2.0)),
    alpha_(0.01),
    groundElectrode_(true),
    outputFileName_()
{
    const dictionary& cfg = verificationDict();

    enabled_ = cfg.lookupOrDefault<Switch>("enabled", true);
    // solutionAlgorithm is a solver-level key (lives in the parent
    // bidomainSolverCoeffs dict, not in verificationModel).
    useExplicitAlgorithm_ =
        dict.lookupOrDefault<word>("solutionAlgorithm", "implicit")
     == "explicit";
    k_ = cfg.lookupOrDefault<scalar>("k", 1.0/Foam::sqrt(2.0));
    alpha_ = cfg.lookupOrDefault<scalar>("alpha", 0.01);
    groundElectrode_ = cfg.lookupOrDefault<Switch>("groundElectrode", true);
    outputFileName_ = cfg.lookupOrDefault<fileName>("outputFile", fileName());

    if (!groundElectrode_)
    {
        FatalErrorInFunction
            << type()
            << " currently implements the FDA ground-electrode variant only. "
            << "Set verificationModel { groundElectrode yes; }."
            << exit(FatalError);
    }
}


void manufacturedFDABathBidomainVerifier::bindBidomainField(volScalarField& phiE)
{
    phiEPtr_ = &phiE;
    phiEHeartCellMapPtr_ = nullptr;
}


void manufacturedFDABathBidomainVerifier::bindBidomainField
(
    volScalarField& phiE,
    const labelUList& heartCellMap
)
{
    phiEPtr_ = &phiE;
    phiEHeartCellMapPtr_ = &heartCellMap;
}


void manufacturedFDABathBidomainVerifier::unbindBidomainField()
{
    phiEPtr_ = nullptr;
    phiEHeartCellMapPtr_ = nullptr;
}


wordList manufacturedFDABathBidomainVerifier::preProcessFieldNames
(
    const ionicModel&
) const
{
    return wordList({"u1", "u2", "u3"});
}


wordList manufacturedFDABathBidomainVerifier::requiredPostProcessFieldNames
(
    const ionicModel&
) const
{
    return wordList({"u1", "u2", "u3"});
}


bool manufacturedFDABathBidomainVerifier::shouldPostProcess
(
    const ionicModel&,
    const volScalarField& Vm
) const
{
    return enabled_ && !errorsReported_ && shouldReportManufacturedErrors(Vm);
}


void manufacturedFDABathBidomainVerifier::preProcess
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

    const label u1Field = requireFieldIndex(names, "u1", "preProcess");
    const label u2Field = requireFieldIndex(names, "u2", "preProcess");
    const label u3Field = requireFieldIndex(names, "u3", "preProcess");

    volScalarField& u1m = fields[u1Field];
    volScalarField& u2m = fields[u2Field];
    volScalarField& u3m = fields[u3Field];

    const volTensorField& Ge =
        Vm.mesh().lookupObject<volTensorField>("conductivityExtracellular");
    const scalar se = manufacturedFDABathSe(Ge);

    const vectorField& centres = Vm.mesh().C().primitiveField();
    scalarField X(centres.component(vector::X));
    scalarField Y(centres.component(vector::Y));
    scalarField Z(centres.component(vector::Z));
    const scalar t = Vm.mesh().time().value();
    const label dimension = model.geometricDimension();

    computeManufacturedFDABathV
    (
        Vm.primitiveFieldRef(), X, Y, Z, t, dimension, alpha_, se
    );
    computeManufacturedFDABathU
    (
        u1m.primitiveFieldRef(),
        u2m.primitiveFieldRef(),
        u3m.primitiveFieldRef(),
        X,
        Y,
        Z,
        t,
        dimension,
        alpha_,
        se
    );

    Vm.correctBoundaryConditions();
    u1m.correctBoundaryConditions();
    u2m.correctBoundaryConditions();
    u3m.correctBoundaryConditions();

    if (phiEPtr_)
    {
        const fvMesh& phiEMesh = phiEPtr_->mesh();
        const vectorField& phiECentres = phiEMesh.C().primitiveField();
        scalarField phiEX(phiECentres.component(vector::X));

        computeManufacturedFDABathPhiE
        (
            phiEPtr_->primitiveFieldRef(), phiEX, t, k_, alpha_, se
        );
        phiEPtr_->correctBoundaryConditions();
    }

    model.importFields(Vm, names, fields);
}


void manufacturedFDABathBidomainVerifier::postProcess
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
            << type() << " requires a bound phiE field before postProcess()."
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

    const label dimension = model.geometricDimension();
    if (dimension < 1 || dimension > 3)
    {
        FatalErrorInFunction
            << type() << " requires a valid geometric dimension in [1,3], "
            << "but ionic model '" << model.type() << "' reported "
            << dimension << "."
            << exit(FatalError);
    }

    const label u1Field =
        requireFieldIndex(requiredNames, "u1", "postProcess");
    const label u2Field =
        requireFieldIndex(requiredNames, "u2", "postProcess");
    const label u3Field =
        requireFieldIndex(requiredNames, "u3", "postProcess");

    const fvMesh& mesh = Vm.mesh();
    const Time& time = mesh.time();
    const volTensorField& Ge =
        mesh.lookupObject<volTensorField>("conductivityExtracellular");
    const scalar se = manufacturedFDABathSe(Ge);

    const vectorField& centres = mesh.C().primitiveField();
    scalarField X(centres.component(vector::X));
    scalarField Y(centres.component(vector::Y));
    scalarField Z(centres.component(vector::Z));
    const scalar t = time.value();

    scalarField VmExact;
    scalarField phiEHeartExact;
    scalarField phiIExact;
    scalarField u1Exact;
    scalarField u2Exact;
    scalarField u3Exact;

    computeManufacturedFDABathV
    (
        VmExact, X, Y, Z, t, dimension, alpha_, se
    );
    computeManufacturedFDABathPhiE
    (
        phiEHeartExact, X, t, k_, alpha_, se
    );
    computeManufacturedFDABathPhiI
    (
        phiIExact, X, Y, Z, t, dimension, k_, alpha_, se
    );
    computeManufacturedFDABathU
    (
        u1Exact, u2Exact, u3Exact, X, Y, Z, t, dimension, alpha_, se
    );

    const scalarField& VmValues = Vm.primitiveField();
    const scalarField& phiEValues = phiEPtr_->primitiveField();
    const scalarField& u1Values = fields[u1Field].primitiveField();
    const scalarField& u2Values = fields[u2Field].primitiveField();
    const scalarField& u3Values = fields[u3Field].primitiveField();

    const auto VmNorms = computeNorms(VmValues, VmExact);

    const fvMesh& phiEMesh = phiEPtr_->mesh();
    const vectorField& phiECentres = phiEMesh.C().primitiveField();
    scalarField phiEX(phiECentres.component(vector::X));
    scalarField phiEExact;
    computeManufacturedFDABathPhiE(phiEExact, phiEX, t, k_, alpha_, se);

    const auto phiENorms = computeNorms(phiEValues, phiEExact);

    auto phiINorms = nanNorms();
    if (phiEValues.size() == VmValues.size() && &phiEMesh == &mesh)
    {
        scalarField phiIValues(VmValues.size(), 0.0);
        forAll(phiIValues, i)
        {
            phiIValues[i] = VmValues[i] + phiEValues[i];
        }

        phiINorms = computeNorms(phiIValues, phiIExact);
    }
    else if (phiEHeartCellMapPtr_)
    {
        const labelUList& heartCellMap = *phiEHeartCellMapPtr_;

        if (heartCellMap.size() != VmValues.size())
        {
            FatalErrorInFunction
                << "Cannot compute manufactured phiI norms: mapped global "
                << "phiE cell map size " << heartCellMap.size()
                << " does not match Vm cell count " << VmValues.size() << "."
                << exit(FatalError);
        }

        scalarField phiIValues(VmValues.size(), 0.0);
        forAll(phiIValues, i)
        {
            const label baseCellI = heartCellMap[i];
            phiIValues[i] = VmValues[i] + phiEValues[baseCellI];
        }

        phiINorms = computeNorms(phiIValues, phiIExact);
    }
    const auto u1Norms = computeNorms(u1Values, u1Exact);
    const auto u2Norms = computeNorms(u2Values, u2Exact);
    const auto u3Norms = computeNorms(u3Values, u3Exact);

    const label totalCells = globalManufacturedCellCount(mesh);
    const label nPerDirection = structuredCellsPerDirection(totalCells, dimension);
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
                "bathBidomain_"
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
            << "Bath-bidomain manufactured-solution error summary (t = "
            << t << "):" << nl
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
            << "u3        " << u3Norms.first().first() << "   "
            << u3Norms.first().second() << "   " << u3Norms.second() << nl
            << "-------------------------------------------------" << endl;

        OFstream os(outputFile);
        os  << "# Bath-bidomain manufactured solution error summary\n"
            << "# dimension " << dimensionName(dimension) << "\n"
            << "# cellsPerDirection " << nPerDirection << "\n"
            << "# steps " << nSteps << "\n"
            << "# dt " << dt << "\n"
            << "# time " << t << "\n"
            << "# k " << k_ << "\n"
            << "# alpha " << alpha_ << "\n"
            << "# se " << se << "\n"
            << "# sigmaB " << manufacturedFDABathSigmaB(se) << "\n"
            << "# groundElectrode " << groundElectrode_ << "\n"
            << "# phiEMesh " << phiEMesh.name() << "\n"
            << "# field L1 L2 Linf\n"
            << "Vm " << VmNorms.first().first() << " "
            << VmNorms.first().second() << " " << VmNorms.second() << "\n"
            << "phiE " << phiENorms.first().first() << " "
            << phiENorms.first().second() << " " << phiENorms.second() << "\n"
            << "phiI " << phiINorms.first().first() << " "
            << phiINorms.first().second() << " " << phiINorms.second() << "\n"
            << "u1 " << u1Norms.first().first() << " "
            << u1Norms.first().second() << " " << u1Norms.second() << "\n"
            << "u2 " << u2Norms.first().first() << " "
            << u2Norms.first().second() << " " << u2Norms.second() << "\n"
            << "u3 " << u3Norms.first().first() << " "
            << u3Norms.first().second() << " " << u3Norms.second() << "\n";
    }

    errorsReported_ = true;
}

} // End namespace Foam

// ************************************************************************* //
