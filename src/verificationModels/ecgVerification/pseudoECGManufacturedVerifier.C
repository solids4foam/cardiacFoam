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

#include "pseudoECGManufacturedVerifier.H"

#include "DynamicList.H"
#include "PstreamReduceOps.H"
#include "addToRunTimeSelectionTable.H"
#include "ecgModelIO.H"
#include "monodomainVerification/manufacturedFDAReference.H"

namespace Foam
{

defineTypeNameAndDebug(pseudoECGManufacturedVerifier, 0);
addToRunTimeSelectionTable
(
    ecgVerificationModel,
    pseudoECGManufacturedVerifier,
    dictionary
);


const dictionary& pseudoECGManufacturedVerifier::verificationDict
(
    const dictionary& dict
)
{
    if (dict.found("pseudoECGManufacturedVerifierCoeffs"))
    {
        return dict.subDict("pseudoECGManufacturedVerifierCoeffs");
    }

    if (dict.found("manufactured"))
    {
        return dict.subDict("manufactured");
    }

    return dict;
}


pseudoECGManufacturedVerifier::pseudoECGManufacturedVerifier
(
    const electroStateProvider& stateProvider,
    const dictionary& dict,
    const wordList& electrodeNames,
    const List<vector>& electrodePositions
)
:
    ecgVerificationModel(stateProvider, electrodeNames, electrodePositions),
    outputPtr_(),
    enabled_(false),
    dimension_(max(label(1), min(mesh_.nGeometricD(), label(3)))),
    referenceQuadratureOrder_(96),
    checkQuadratureOrders_(),
    referenceSpatialValues_(electrodePositions_.size(), scalar(0)),
    checkSpatialValues_(),
    cachedConductivity_(tensor::zero),
    spatialCacheValid_(false),
    sampleCount_(0),
    referenceErrorL1Sum_(electrodePositions_.size(), scalar(0)),
    referenceErrorL2Sum_(electrodePositions_.size(), scalar(0)),
    referenceErrorLinf_(electrodePositions_.size(), scalar(0)),
    checkErrorL1Sum_(),
    checkErrorL2Sum_(),
    checkErrorLinf_(),
    referenceDeltaL1Sum_(),
    referenceDeltaL2Sum_(),
    referenceDeltaLinf_(),
    summaryWritten_(false)
{
    read(dict);
}


ecgVerificationModel::Requirements
pseudoECGManufacturedVerifier::requirements() const
{
    Requirements needs;
    needs.needConductivity = true;
    return needs;
}


void pseudoECGManufacturedVerifier::resizeCheckStorage()
{
    const label nChecks = checkQuadratureOrders_.size();
    const label nElectrodes = electrodePositions_.size();

    checkSpatialValues_.setSize(nChecks);
    checkErrorL1Sum_.setSize(nChecks);
    checkErrorL2Sum_.setSize(nChecks);
    checkErrorLinf_.setSize(nChecks);
    referenceDeltaL1Sum_.setSize(nChecks);
    referenceDeltaL2Sum_.setSize(nChecks);
    referenceDeltaLinf_.setSize(nChecks);

    for (label checkI = 0; checkI < nChecks; ++checkI)
    {
        checkSpatialValues_[checkI].setSize(nElectrodes, scalar(0));
        checkErrorL1Sum_[checkI].setSize(nElectrodes, scalar(0));
        checkErrorL2Sum_[checkI].setSize(nElectrodes, scalar(0));
        checkErrorLinf_[checkI].setSize(nElectrodes, scalar(0));
        referenceDeltaL1Sum_[checkI].setSize(nElectrodes, scalar(0));
        referenceDeltaL2Sum_[checkI].setSize(nElectrodes, scalar(0));
        referenceDeltaLinf_[checkI].setSize(nElectrodes, scalar(0));
    }
}


void pseudoECGManufacturedVerifier::initialiseOutput()
{
    const fileName outDir(mesh_.time().path() / "postProcessing");
    wordList columns;

    forAll(electrodeNames_, electrodeI)
    {
        const word& name = electrodeNames_[electrodeI];
        columns.append("numeric_" + name);

        forAll(checkQuadratureOrders_, checkI)
        {
            columns.append
            (
                "refQ" + Foam::name(checkQuadratureOrders_[checkI]) + "_" + name
            );
        }

        columns.append
        (
            "refQ" + Foam::name(referenceQuadratureOrder_) + "_" + name
        );

        forAll(checkQuadratureOrders_, checkI)
        {
            columns.append
            (
                "errQ" + Foam::name(checkQuadratureOrders_[checkI]) + "_" + name
            );
        }

        columns.append
        (
            "errQ" + Foam::name(referenceQuadratureOrder_) + "_" + name
        );

        forAll(checkQuadratureOrders_, checkI)
        {
            columns.append
            (
                "deltaQuadratureQ"
              + Foam::name(checkQuadratureOrders_[checkI])
              + "_Q"
              + Foam::name(referenceQuadratureOrder_)
              + "_"
              + name
            );
        }
    }

    outputPtr_ =
        ecgModelIO::openTimeSeries(outDir, "manufacturedPseudoECG.dat", columns);
}


void pseudoECGManufacturedVerifier::invalidateReferenceCache()
{
    referenceSpatialValues_.clear();
    checkSpatialValues_.clear();
    spatialCacheValid_ = false;
    cachedConductivity_ = tensor::zero;
}


void pseudoECGManufacturedVerifier::rebuildReferenceCache
(
    const tensor& referenceConductivity
)
{
    List<scalar> referenceNodes, referenceWeights;
    pseudoECGManufacturedQuadratureRule
    (
        referenceQuadratureOrder_,
        referenceNodes,
        referenceWeights
    );

    List<List<scalar>> checkNodes(checkQuadratureOrders_.size());
    List<List<scalar>> checkWeights(checkQuadratureOrders_.size());

    forAll(checkQuadratureOrders_, checkI)
    {
        pseudoECGManufacturedQuadratureRule
        (
            checkQuadratureOrders_[checkI],
            checkNodes[checkI],
            checkWeights[checkI]
        );
    }

    referenceSpatialValues_.setSize(electrodePositions_.size());
    checkSpatialValues_.setSize(checkQuadratureOrders_.size());

    forAll(checkSpatialValues_, checkI)
    {
        checkSpatialValues_[checkI].setSize(electrodePositions_.size());
    }

    forAll(electrodePositions_, electrodeI)
    {
        forAll(checkQuadratureOrders_, checkI)
        {
            checkSpatialValues_[checkI][electrodeI] =
                computeManufacturedPseudoECGSpatialReference
                (
                    electrodePositions_[electrodeI],
                    referenceConductivity,
                    dimension_,
                    checkNodes[checkI],
                    checkWeights[checkI]
                );
        }

        referenceSpatialValues_[electrodeI] =
            computeManufacturedPseudoECGSpatialReference
            (
                electrodePositions_[electrodeI],
                referenceConductivity,
                dimension_,
                referenceNodes,
                referenceWeights
            );
    }

    cachedConductivity_ = referenceConductivity;
    spatialCacheValid_ = true;
}


void pseudoECGManufacturedVerifier::updateStatistics
(
    const List<scalar>& numericValues,
    const List<List<scalar>>& checkReferenceValues,
    const List<scalar>& referenceValues
)
{
    ++sampleCount_;

    forAll(numericValues, electrodeI)
    {
        const scalar referenceError =
            Foam::mag(numericValues[electrodeI] - referenceValues[electrodeI]);

        referenceErrorL1Sum_[electrodeI] += referenceError;
        referenceErrorL2Sum_[electrodeI] += referenceError*referenceError;
        referenceErrorLinf_[electrodeI] =
            max(referenceErrorLinf_[electrodeI], referenceError);

        forAll(checkQuadratureOrders_, checkI)
        {
            const scalar checkError = Foam::mag
            (
                numericValues[electrodeI] - checkReferenceValues[checkI][electrodeI]
            );
            const scalar referenceDelta = Foam::mag
            (
                referenceValues[electrodeI] - checkReferenceValues[checkI][electrodeI]
            );

            checkErrorL1Sum_[checkI][electrodeI] += checkError;
            checkErrorL2Sum_[checkI][electrodeI] += checkError*checkError;
            checkErrorLinf_[checkI][electrodeI] =
                max(checkErrorLinf_[checkI][electrodeI], checkError);

            referenceDeltaL1Sum_[checkI][electrodeI] += referenceDelta;
            referenceDeltaL2Sum_[checkI][electrodeI] +=
                referenceDelta*referenceDelta;
            referenceDeltaLinf_[checkI][electrodeI] =
                max(referenceDeltaLinf_[checkI][electrodeI], referenceDelta);
        }
    }
}


void pseudoECGManufacturedVerifier::writeSummary()
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

    const fileName outputFile
    (
        mesh_.time().path() / "postProcessing" / "manufacturedPseudoECGSummary.dat"
    );
    OFstream os(outputFile);

    os << "Manufactured pseudo-ECG summary\n";
    os << "samples " << sampleCount_ << "\n";
    os << "dimension " << dimension_ << "D\n";
    os << "qChecks";
    forAll(checkQuadratureOrders_, checkI)
    {
        os << " " << checkQuadratureOrders_[checkI];
    }
    os << "\n";
    os << "qReference " << referenceQuadratureOrder_ << "\n";
    os << "Electrode  L1_err_ref  L2_err_ref  Linf_err_ref";
    forAll(checkQuadratureOrders_, checkI)
    {
        const label qCheck = checkQuadratureOrders_[checkI];
        os << "  L1_err_q" << qCheck
           << "  L2_err_q" << qCheck
           << "  Linf_err_q" << qCheck
           << "  L1_delta_q" << qCheck << "_ref"
           << "  L2_delta_q" << qCheck << "_ref"
           << "  Linf_delta_q" << qCheck << "_ref";
    }
    os << "\n";

    const scalar count = max(scalar(1), scalar(sampleCount_));

    forAll(electrodeNames_, electrodeI)
    {
        os << electrodeNames_[electrodeI] << " "
           << referenceErrorL1Sum_[electrodeI]/count << " "
           << Foam::sqrt(referenceErrorL2Sum_[electrodeI]/count) << " "
           << referenceErrorLinf_[electrodeI];

        forAll(checkQuadratureOrders_, checkI)
        {
            os << " " << checkErrorL1Sum_[checkI][electrodeI]/count
               << " "
               << Foam::sqrt(checkErrorL2Sum_[checkI][electrodeI]/count)
               << " " << checkErrorLinf_[checkI][electrodeI]
               << " "
               << referenceDeltaL1Sum_[checkI][electrodeI]/count
               << " "
               << Foam::sqrt(referenceDeltaL2Sum_[checkI][electrodeI]/count)
               << " " << referenceDeltaLinf_[checkI][electrodeI];
        }

        os << "\n";
    }
}


bool pseudoECGManufacturedVerifier::read(const dictionary& dict)
{
    const Switch wasEnabled = enabled_;
    const label previousDimension = dimension_;
    const label previousReferenceQuadratureOrder = referenceQuadratureOrder_;
    const List<label> previousCheckQuadratureOrders = checkQuadratureOrders_;

    enabled_ = false;
    dimension_ = max(label(1), min(mesh_.nGeometricD(), label(3)));
    referenceQuadratureOrder_ = 96;
    checkQuadratureOrders_.setSize(4);
    checkQuadratureOrders_[0] = 6;
    checkQuadratureOrders_[1] = 12;
    checkQuadratureOrders_[2] = 24;
    checkQuadratureOrders_[3] = 48;

    const dictionary& cfg = verificationDict(dict);
    enabled_ = cfg.lookupOrDefault<Switch>("enabled", true);

    if (cfg.found("dimension"))
    {
        const word dimensionName(cfg.lookup("dimension"));

        if (dimensionName == "1D")
        {
            dimension_ = 1;
        }
        else if (dimensionName == "2D")
        {
            dimension_ = 2;
        }
        else if (dimensionName == "3D")
        {
            dimension_ = 3;
        }
        else
        {
            FatalErrorInFunction
                << "Unsupported manufactured pseudo-ECG dimension '"
                << dimensionName << "'. Expected one of 1D, 2D, or 3D."
                << exit(FatalError);
        }
    }

    referenceQuadratureOrder_ =
        cfg.lookupOrDefault<label>("referenceQuadratureOrder", 96);

    if (cfg.found("checkQuadratureOrders"))
    {
        cfg.lookup("checkQuadratureOrders") >> checkQuadratureOrders_;
    }
    else
    {
        checkQuadratureOrders_.setSize(1);
        checkQuadratureOrders_[0] =
            cfg.lookupOrDefault<label>("checkQuadratureOrder", 6);
    }

    if (checkQuadratureOrders_.empty())
    {
        FatalErrorInFunction
            << "Manufactured pseudo-ECG requires at least one entry in "
               "checkQuadratureOrders."
            << exit(FatalError);
    }

    for (label i = 0; i < checkQuadratureOrders_.size(); ++i)
    {
        for (label j = i + 1; j < checkQuadratureOrders_.size(); ++j)
        {
            if (checkQuadratureOrders_[j] < checkQuadratureOrders_[i])
            {
                Swap(checkQuadratureOrders_[i], checkQuadratureOrders_[j]);
            }
        }
    }

    label uniqueCount = 0;
    forAll(checkQuadratureOrders_, checkI)
    {
        const label current = checkQuadratureOrders_[checkI];
        if (checkI == 0 || current != checkQuadratureOrders_[uniqueCount - 1])
        {
            checkQuadratureOrders_[uniqueCount++] = current;
        }
    }
    checkQuadratureOrders_.setSize(uniqueCount);

    if (!pseudoECGManufacturedSupportsQuadratureOrder(referenceQuadratureOrder_))
    {
        FatalErrorInFunction
            << "Manufactured pseudo-ECG quadrature orders must be positive. "
            << "Requested qRef=" << referenceQuadratureOrder_ << "."
            << exit(FatalError);
    }

    forAll(checkQuadratureOrders_, checkI)
    {
        if (!pseudoECGManufacturedSupportsQuadratureOrder
            (
                checkQuadratureOrders_[checkI]
            ))
        {
            FatalErrorInFunction
                << "Manufactured pseudo-ECG quadrature orders must be positive. "
                << "Invalid qCheck=" << checkQuadratureOrders_[checkI] << "."
                << exit(FatalError);
        }
    }

    resizeCheckStorage();

    const bool configurationChanged =
        enabled_ != wasEnabled
     || dimension_ != previousDimension
     || referenceQuadratureOrder_ != previousReferenceQuadratureOrder
     || checkQuadratureOrders_ != previousCheckQuadratureOrders;

    if (outputPtr_.valid() && configurationChanged && enabled_)
    {
        FatalErrorInFunction
            << "Changing manufactured pseudo-ECG dimension or quadrature "
               "orders during read() is not supported because the output "
               "layout and semantics are fixed at startup."
            << exit(FatalError);
    }

    if (!enabled_)
    {
        invalidateReferenceCache();
    }
    else
    {
        invalidateReferenceCache();

        if (!outputPtr_.valid())
        {
            initialiseOutput();
        }
    }

    if (enabled_)
    {
        Info << "Enabled manufactured pseudo-ECG references: qRef="
             << referenceQuadratureOrder_ << ", qChecks=(";
        forAll(checkQuadratureOrders_, checkI)
        {
            if (checkI > 0)
            {
                Info << " ";
            }
            Info << checkQuadratureOrders_[checkI];
        }
        Info << "), dimension=" << dimension_ << "." << endl;
    }
    else
    {
        Info << "Manufactured pseudo-ECG verification disabled." << endl;
    }

    return true;
}


void pseudoECGManufacturedVerifier::record(const List<scalar>& numericValues)
{
    if (!enabled_)
    {
        return;
    }

    const tensorField& conductivityField =
        requireConductivity().primitiveField();

    if (conductivityField.empty())
    {
        return;
    }

    const scalar timeValue = mesh_.time().value();
    const tensor referenceConductivity = conductivityField[0];

    if
    (
        !spatialCacheValid_
     || Foam::magSqr(referenceConductivity - cachedConductivity_) > SMALL
    )
    {
        rebuildReferenceCache(referenceConductivity);
    }

    const scalar timeFactor = pseudoECGManufacturedTimeFactor(timeValue);
    List<List<scalar>> checkReferenceValues(checkQuadratureOrders_.size());

    forAll(checkReferenceValues, checkI)
    {
        checkReferenceValues[checkI].setSize
        (
            electrodePositions_.size(), scalar(0)
        );
    }

    List<scalar> referenceValues(electrodePositions_.size(), scalar(0));
    DynamicList<scalar> row;
    row.reserve
    (
        electrodePositions_.size()*(2 + 3*checkQuadratureOrders_.size())
    );

    forAll(electrodePositions_, electrodeI)
    {
        referenceValues[electrodeI] =
            timeFactor*referenceSpatialValues_[electrodeI];
        row.append(numericValues[electrodeI]);

        forAll(checkQuadratureOrders_, checkI)
        {
            checkReferenceValues[checkI][electrodeI] =
                timeFactor*checkSpatialValues_[checkI][electrodeI];
            row.append(checkReferenceValues[checkI][electrodeI]);
        }

        row.append(referenceValues[electrodeI]);

        forAll(checkQuadratureOrders_, checkI)
        {
            row.append
            (
                Foam::mag
                (
                    numericValues[electrodeI]
                  - checkReferenceValues[checkI][electrodeI]
                )
            );
        }

        row.append(Foam::mag(numericValues[electrodeI] - referenceValues[electrodeI]));

        forAll(checkQuadratureOrders_, checkI)
        {
            row.append
            (
                Foam::mag
                (
                    referenceValues[electrodeI]
                  - checkReferenceValues[checkI][electrodeI]
                )
            );
        }
    }

    updateStatistics(numericValues, checkReferenceValues, referenceValues);

    List<scalar> values(row.size(), scalar(0));
    forAll(values, valueI)
    {
        values[valueI] = row[valueI];
    }
    ecgModelIO::writeRow(outputPtr_.ref(), timeValue, values);

    const Time& time = mesh_.time();
    if (time.value() + 0.5*time.deltaTValue() >= time.endTime().value())
    {
        writeSummary();
    }
}


void pseudoECGManufacturedVerifier::end()
{
    writeSummary();
}

} // End namespace Foam

// ************************************************************************* //
