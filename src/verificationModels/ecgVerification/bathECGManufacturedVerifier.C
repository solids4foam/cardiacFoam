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

#include "bathECGManufacturedVerifier.H"

#include "DynamicList.H"
#include "PstreamReduceOps.H"
#include "addToRunTimeSelectionTable.H"
#include "bathBidomainVerification/manufacturedFDABathBidomainReference.H"
#include "ecgModelIO.H"

namespace Foam
{

defineTypeNameAndDebug(bathECGManufacturedVerifier, 0);
addToRunTimeSelectionTable
(
    ecgVerificationModel,
    bathECGManufacturedVerifier,
    dictionary
);

namespace
{

Tuple2<Tuple2<scalar, scalar>, scalar> computeNorms
(
    const scalarField& numeric,
    const scalarField& exact
)
{
    if (numeric.size() != exact.size())
    {
        FatalErrorInFunction
            << "Cannot compute bath ECG manufactured norms for fields with "
            << "different sizes: numeric=" << numeric.size()
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


scalarField phiEOnVerifierMesh
(
    const volScalarField& phiE,
    const electroStateProvider& stateProvider,
    const fvMesh& verifierMesh
)
{
    if (phiE.mesh().nCells() == verifierMesh.nCells())
    {
        return phiE.primitiveField();
    }

    const labelUList* cellMapPtr = stateProvider.subsetCellMapPtr();

    if (!cellMapPtr)
    {
        FatalErrorInFunction
            << "Cannot map phiE from mesh '" << phiE.mesh().name()
            << "' to ECG verifier mesh '" << verifierMesh.name()
            << "' because no subset cell map is available."
            << exit(FatalError);
    }

    const labelUList& cellMap = *cellMapPtr;

    if (cellMap.size() != verifierMesh.nCells())
    {
        FatalErrorInFunction
            << "ECG verifier subset map size " << cellMap.size()
            << " does not match verifier mesh cell count "
            << verifierMesh.nCells() << "."
            << exit(FatalError);
    }

    scalarField mapped(verifierMesh.nCells(), 0.0);
    const scalarField& basePhiE = phiE.primitiveField();

    forAll(mapped, cellI)
    {
        mapped[cellI] = basePhiE[cellMap[cellI]];
    }

    return mapped;
}


void validateGroundedFDABathBoundarySetup
(
    const dictionary& dict,
    const scalar alpha
)
{
    if (!dict.found("groundPatches"))
    {
        FatalErrorInFunction
            << "bathECGManufacturedVerifier supports only the FDA grounded "
            << "bath variant. Configure groundPatches { xMin 0; }."
            << exit(FatalError);
    }

    const dictionary& groundDict = dict.subDict("groundPatches");
    const wordList groundNames(groundDict.toc());

    if (groundNames.size() != 1 || !groundDict.found("xMin"))
    {
        FatalErrorInFunction
            << "bathECGManufacturedVerifier supports only the FDA grounded "
            << "bath variant with exactly one ground patch: "
            << "groundPatches { xMin 0; }."
            << exit(FatalError);
    }

    const scalar groundValue = groundDict.get<scalar>("xMin");

    if (mag(groundValue) > SMALL)
    {
        FatalErrorInFunction
            << "bathECGManufacturedVerifier requires phiE = 0 on xMin. "
            << "Found groundPatches { xMin " << groundValue << "; }."
            << exit(FatalError);
    }

    if (!dict.found("surfaceCurrentPatches"))
    {
        FatalErrorInFunction
            << "bathECGManufacturedVerifier supports only the FDA grounded "
            << "bath variant. Configure surfaceCurrentPatches { xMax alpha; }."
            << exit(FatalError);
    }

    const dictionary& currentDict = dict.subDict("surfaceCurrentPatches");
    const wordList currentNames(currentDict.toc());

    if (currentNames.size() != 1 || !currentDict.found("xMax"))
    {
        FatalErrorInFunction
            << "bathECGManufacturedVerifier supports only the FDA grounded "
            << "bath variant with exactly one surface-current patch: "
            << "surfaceCurrentPatches { xMax alpha; }."
            << exit(FatalError);
    }

    const scalar currentValue = currentDict.get<scalar>("xMax");

    if (mag(currentValue - alpha) > SMALL)
    {
        FatalErrorInFunction
            << "bathECGManufacturedVerifier requires +alpha on xMax. "
            << "Expected " << alpha << " but found " << currentValue << "."
            << exit(FatalError);
    }
}

}


bathECGManufacturedVerifier::bathECGManufacturedVerifier
(
    const electroStateProvider& stateProvider,
    const dictionary& dict,
    const wordList& electrodeNames,
    const List<vector>& electrodePositions
)
:
    ecgVerificationModel(stateProvider, electrodeNames, electrodePositions),
    outputPtr_(),
    enabled_(true),
    k_(1.0/Foam::sqrt(2.0)),
    alpha_(0.01),
    groundElectrode_(true),
    sampleCount_(0),
    fieldErrorL1Sum_(0.0),
    fieldErrorL2Sum_(0.0),
    fieldErrorLinf_(0.0),
    electrodeErrorL1Sum_(electrodePositions_.size(), 0.0),
    electrodeErrorL2Sum_(electrodePositions_.size(), 0.0),
    electrodeErrorLinf_(electrodePositions_.size(), 0.0),
    summaryWritten_(false)
{
    read(dict);
}


ecgVerificationModel::Requirements
bathECGManufacturedVerifier::requirements() const
{
    Requirements needs;
    needs.needPhiE = true;
    return needs;
}


void bathECGManufacturedVerifier::resizeStatistics()
{
    electrodeErrorL1Sum_.setSize(electrodePositions_.size(), 0.0);
    electrodeErrorL2Sum_.setSize(electrodePositions_.size(), 0.0);
    electrodeErrorLinf_.setSize(electrodePositions_.size(), 0.0);
}


void bathECGManufacturedVerifier::initialiseOutput()
{
    const fileName outDir(mesh_.time().path() / "postProcessing");
    wordList columns;

    columns.append("field_L1");
    columns.append("field_L2");
    columns.append("field_Linf");

    forAll(electrodeNames_, electrodeI)
    {
        const word& name = electrodeNames_[electrodeI];
        columns.append("numeric_" + name);
        columns.append("exact_" + name);
        columns.append("error_" + name);
    }

    outputPtr_ =
        ecgModelIO::openTimeSeries(outDir, "manufacturedBathECG.dat", columns);
}


bool bathECGManufacturedVerifier::read(const dictionary& dict)
{
    const dictionary& manufactured = dict.subDict("manufacturedBidomain");

    enabled_ = manufactured.lookupOrDefault<Switch>("enabled", true);
    k_ = manufactured.lookupOrDefault<scalar>("k", 1.0/Foam::sqrt(2.0));
    alpha_ = manufactured.lookupOrDefault<scalar>("alpha", 0.01);
    groundElectrode_ =
        manufactured.lookupOrDefault<Switch>("groundElectrode", true);

    if (!groundElectrode_)
    {
        FatalErrorInFunction
            << type()
            << " currently implements the FDA ground-electrode variant only. "
            << "Set manufacturedBidomain { groundElectrode yes; }."
            << exit(FatalError);
    }

    if (enabled_)
    {
        validateGroundedFDABathBoundarySetup(dict, alpha_);
    }

    resizeStatistics();

    if (enabled_ && !outputPtr_.valid())
    {
        initialiseOutput();
    }

    Info<< (enabled_ ? "Enabled" : "Disabled")
        << " bath ECG manufactured verification: k=" << k_
        << ", alpha=" << alpha_
        << ", groundElectrode=" << groundElectrode_ << "." << endl;

    return true;
}


void bathECGManufacturedVerifier::record(const List<scalar>& numericValues)
{
    if (!enabled_)
    {
        return;
    }

    const volScalarField& phiE = requirePhiE();
    const volTensorField* GePtr = stateProvider_.extracellularConductivityPtr();

    if (!GePtr)
    {
        FatalErrorInFunction
            << type()
            << " requires an upstream extracellular conductivity field."
            << exit(FatalError);
    }

    const scalar se = manufacturedFDABathSe(*GePtr);
    const scalar t = mesh_.time().value();

    const vectorField& centres = mesh_.C().primitiveField();
    scalarField X(centres.component(vector::X));
    scalarField phiEExact;
    computeManufacturedFDABathPhiE(phiEExact, X, t, k_, alpha_, se);

    const scalarField phiEValues =
        phiEOnVerifierMesh(phiE, stateProvider_, mesh_);
    const auto fieldNorms = computeNorms(phiEValues, phiEExact);

    ++sampleCount_;
    fieldErrorL1Sum_ += fieldNorms.first().first();
    fieldErrorL2Sum_ += sqr(fieldNorms.first().second());
    fieldErrorLinf_ = max(fieldErrorLinf_, fieldNorms.second());

    DynamicList<scalar> row;
    row.reserve(3 + 3*electrodePositions_.size());
    row.append(fieldNorms.first().first());
    row.append(fieldNorms.first().second());
    row.append(fieldNorms.second());

    forAll(electrodePositions_, electrodeI)
    {
        const scalar exact = manufacturedFDABathPhiEAt
        (
            electrodePositions_[electrodeI].x(),
            t,
            k_,
            alpha_,
            se
        );
        const scalar error = Foam::mag(numericValues[electrodeI] - exact);

        electrodeErrorL1Sum_[electrodeI] += error;
        electrodeErrorL2Sum_[electrodeI] += error*error;
        electrodeErrorLinf_[electrodeI] =
            max(electrodeErrorLinf_[electrodeI], error);

        row.append(numericValues[electrodeI]);
        row.append(exact);
        row.append(error);
    }

    List<scalar> values(row.size(), scalar(0));
    forAll(values, valueI)
    {
        values[valueI] = row[valueI];
    }

    ecgModelIO::writeRow(outputPtr_.ref(), t, values);

    const Time& time = mesh_.time();
    if (time.value() + 0.5*time.deltaTValue() >= time.endTime().value())
    {
        writeSummary();
    }
}


void bathECGManufacturedVerifier::writeSummary()
{
    if (!enabled_ || summaryWritten_)
    {
        return;
    }

    summaryWritten_ = true;

    if (!Pstream::master())
    {
        return;
    }

    const scalar count = max(scalar(1), scalar(sampleCount_));
    const fileName outputFile
    (
        mesh_.time().path() / "postProcessing" / "manufacturedBathECGSummary.dat"
    );
    OFstream os(outputFile);

    os << "Manufactured bath ECG summary\n";
    os << "samples " << sampleCount_ << "\n";
    os << "k " << k_ << "\n";
    os << "alpha " << alpha_ << "\n";
    os << "groundElectrode " << groundElectrode_ << "\n";
    os << "field_L1 " << fieldErrorL1Sum_/count << "\n";
    os << "field_L2 " << Foam::sqrt(fieldErrorL2Sum_/count) << "\n";
    os << "field_Linf " << fieldErrorLinf_ << "\n";
    os << "Electrode L1 L2 Linf\n";

    forAll(electrodeNames_, electrodeI)
    {
        os << electrodeNames_[electrodeI] << " "
           << electrodeErrorL1Sum_[electrodeI]/count << " "
           << Foam::sqrt(electrodeErrorL2Sum_[electrodeI]/count) << " "
           << electrodeErrorLinf_[electrodeI] << "\n";
    }
}


void bathECGManufacturedVerifier::end()
{
    writeSummary();
}

} // End namespace Foam

// ************************************************************************* //
