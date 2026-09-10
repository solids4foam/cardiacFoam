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

#include "manufacturedFDABidomainVerifier.H"

#include "OFstream.H"
#include "OSspecific.H"
#include "bidomainVerification/manufacturedFDABidomainReference.H"
#include "ionicModel.H"
#include "verificationUtils.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
using namespace verificationUtils;

defineTypeNameAndDebug(manufacturedFDABidomainVerifier, 0);
addToRunTimeSelectionTable
(
    electroVerificationModel,
    manufacturedFDABidomainVerifier,
    dictionary
);

manufacturedFDABidomainVerifier::manufacturedFDABidomainVerifier
(
    const dictionary& dict
)
:
    electroVerificationModel(dict),
    phiEPtr_(nullptr),
    errorsReported_(false),
    k_(1.0/Foam::sqrt(2.0)),
    phiEReferenceValue_(0.0),
    phiEReferencePoint_(point::zero)
{
    const dictionary& coeffDict = this->dict();
    const dictionary& cfg = verificationDict();

    k_ = cfg.lookupOrDefault<scalar>("k", 1.0/Foam::sqrt(2.0));
    // phiEReferenceValue and phiERefPoint are physics parameters that belong
    // in the bidomainSolverCoeffs dict (the parent dict), so we read them from coeffDict.
    phiEReferenceValue_ = coeffDict.lookupOrDefault<scalar>("phiEReferenceValue", 0.0);
    phiEReferencePoint_ = coeffDict.get<point>("phiERefPoint");
}


const dictionary& manufacturedFDABidomainVerifier::verificationDict() const
{
    return dict().subDict("verificationModel");
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
    return wordList({"u1", "u2", "u3"});
}


wordList manufacturedFDABidomainVerifier::requiredPostProcessFieldNames
(
    const ionicModel&
) const
{
    return wordList({"u1", "u2"});
}


bool manufacturedFDABidomainVerifier::shouldPostProcess
(
    const ionicModel&,
    const volScalarField& Vm
) const
{
    return !errorsReported_ && shouldReportManufacturedErrors(Vm);
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

        const scalar shift = manufacturedFDABidomainShift
        (
            t,
            dimension,
            k_,
            phiEReferencePoint_,
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

    const scalar shift = manufacturedFDABidomainShift
    (
        t,
        dimension,
        k_,
        phiEReferencePoint_,
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
    const label nPerDirection = structuredCellsPerDirection(totalCells, dimension);
    const scalar dx = structuredManufacturedDx(nPerDirection);
    const scalar dt = time.deltaTValue();
    const label nSteps = max(label(0), time.timeIndex());

    const fileName outputDir(time.globalPath()/"postProcessing");
    mkDir(outputDir);
    const fileName outputFile =
        outputDir
      / (
            dimensionName(dimension) + "_" + Foam::name(nPerDirection) + "_cells.dat"
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
        out << "k                     = " << k_ << "\n";
        out << "phiE reference point  = " << phiEReferencePoint_ << "\n";
        out << "phiE reference value  = " << phiEReferenceValue_ << "\n";
        out << "phiE gauge offset     = " << phiEGaugeOffset << "\n";
    }

    errorsReported_ = true;
}

} // End namespace Foam

// ************************************************************************* //
