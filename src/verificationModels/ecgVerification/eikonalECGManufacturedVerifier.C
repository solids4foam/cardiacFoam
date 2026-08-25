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

#include "eikonalECGManufacturedVerifier.H"

#include "DynamicList.H"
#include "PstreamReduceOps.H"
#include "addToRunTimeSelectionTable.H"
#include "ecgModelIO.H"
#include "eikonalVerification/manufacturedEikonalReference.H"

namespace Foam
{

defineTypeNameAndDebug(eikonalECGManufacturedVerifier, 0);
addToRunTimeSelectionTable
(
    ecgVerificationModel,
    eikonalECGManufacturedVerifier,
    dictionary
);


eikonalECGManufacturedVerifier::eikonalECGManufacturedVerifier
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
    k_(vector::zero),
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
eikonalECGManufacturedVerifier::requirements() const
{
    Requirements needs;
    needs.needActivationTime = true;
    needs.needConductivity = true;
    needs.needChi = true;
    needs.needCm = true;
    needs.needC0 = true;
    return needs;
}


const dictionary& eikonalECGManufacturedVerifier::manufacturedDict
(
    const dictionary& dict
) const
{
    return dict.subDict("verificationModel");
}


void eikonalECGManufacturedVerifier::resizeStorage()
{
    const label nChecks = checkQuadratureOrders_.size();
    const label nElectrodes = electrodePositions_.size();

    referenceErrorL1Sum_.setSize(nElectrodes, scalar(0));
    referenceErrorL2Sum_.setSize(nElectrodes, scalar(0));
    referenceErrorLinf_.setSize(nElectrodes, scalar(0));
    checkErrorL1Sum_.setSize(nChecks);
    checkErrorL2Sum_.setSize(nChecks);
    checkErrorLinf_.setSize(nChecks);
    referenceDeltaL1Sum_.setSize(nChecks);
    referenceDeltaL2Sum_.setSize(nChecks);
    referenceDeltaLinf_.setSize(nChecks);

    for (label checkI = 0; checkI < nChecks; ++checkI)
    {
        checkErrorL1Sum_[checkI].setSize(nElectrodes, scalar(0));
        checkErrorL2Sum_[checkI].setSize(nElectrodes, scalar(0));
        checkErrorLinf_[checkI].setSize(nElectrodes, scalar(0));
        referenceDeltaL1Sum_[checkI].setSize(nElectrodes, scalar(0));
        referenceDeltaL2Sum_[checkI].setSize(nElectrodes, scalar(0));
        referenceDeltaLinf_[checkI].setSize(nElectrodes, scalar(0));
    }
}


void eikonalECGManufacturedVerifier::initialiseOutput()
{
    const fileName outDir(mesh_.time().globalPath() / "postProcessing");
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
        ecgModelIO::openTimeSeries(outDir, "manufacturedEikonalECG.dat", columns);
}


bool eikonalECGManufacturedVerifier::read(const dictionary& dict)
{
    const dictionary& cfg = manufacturedDict(dict);

    enabled_ = cfg.lookupOrDefault<Switch>("enabled", true);
    const label meshDimension = max(label(1), min(mesh_.nGeometricD(), label(3)));
    dimension_ = meshDimension;

    if (cfg.found("dimension"))
    {
        const word dimensionName(cfg.lookup("dimension"));
        label requestedDimension = 0;

        if (dimensionName == "1D")
        {
            requestedDimension = 1;
        }
        else if (dimensionName == "2D")
        {
            requestedDimension = 2;
        }
        else if (dimensionName == "3D")
        {
            requestedDimension = 3;
        }
        else
        {
            FatalErrorInFunction
                << "Unsupported manufactured eikonal ECG dimension '"
                << dimensionName << "'. Expected one of 1D, 2D, or 3D."
                << exit(FatalError);
        }

        if (requestedDimension != meshDimension)
        {
            FatalErrorInFunction
                << "Manufactured eikonal ECG dimension '" << dimensionName
                << "' does not match mesh geometric dimension "
                << meshDimension << "D. Use a matching dimension or omit the "
                << "dimension entry."
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
        checkQuadratureOrders_[0] = 6;
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
            << "Manufactured eikonal ECG quadrature orders must be positive. "
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
                << "Manufactured eikonal ECG quadrature orders must be "
                << "positive. Invalid qCheck="
                << checkQuadratureOrders_[checkI] << "."
                << exit(FatalError);
        }
    }

    if (enabled_)
    {
        validateManufacturedEikonalUnitDomain(mesh_, dimension_);
        (void)manufacturedEikonalConstantConductivity(requireConductivity());
        (void)requireActivationTime();
    }

    resizeStorage();

    if (enabled_ && !outputPtr_.valid())
    {
        initialiseOutput();
    }

    Info<< (enabled_ ? "Enabled" : "Disabled")
        << " manufactured eikonal ECG verification: qRef="
        << referenceQuadratureOrder_ << ", dimension=" << dimension_
        << "." << endl;

    return true;
}


void eikonalECGManufacturedVerifier::updateStatistics
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
                numericValues[electrodeI]
              - checkReferenceValues[checkI][electrodeI]
            );
            const scalar referenceDelta = Foam::mag
            (
                referenceValues[electrodeI]
              - checkReferenceValues[checkI][electrodeI]
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


void eikonalECGManufacturedVerifier::record
(
    const List<scalar>& numericValues
)
{
    record(mesh_.time().value(), numericValues);
}


void eikonalECGManufacturedVerifier::record
(
    scalar sampleTime,
    const List<scalar>& numericValues
)
{
    if (!enabled_)
    {
        return;
    }

    const tensor conductivity =
        manufacturedEikonalConstantConductivity(requireConductivity());

    k_ = manufacturedEikonalK(dimension_);

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

    List<List<scalar>> checkReferenceValues(checkQuadratureOrders_.size());
    forAll(checkReferenceValues, checkI)
    {
        checkReferenceValues[checkI].setSize
        (
            electrodePositions_.size(),
            scalar(0)
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
        forAll(checkQuadratureOrders_, checkI)
        {
            checkReferenceValues[checkI][electrodeI] =
                computeManufacturedEikonalECGReference
                (
                    sampleTime,
                    electrodePositions_[electrodeI],
                    conductivity,
                    k_,
                    dimension_,
                    checkNodes[checkI],
                    checkWeights[checkI]
                );
        }

        referenceValues[electrodeI] =
            computeManufacturedEikonalECGReference
            (
                sampleTime,
                electrodePositions_[electrodeI],
                conductivity,
                k_,
                dimension_,
                referenceNodes,
                referenceWeights
            );

        row.append(numericValues[electrodeI]);

        forAll(checkQuadratureOrders_, checkI)
        {
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

        row.append
        (
            Foam::mag(numericValues[electrodeI] - referenceValues[electrodeI])
        );

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
    ecgModelIO::writeRow(outputPtr_.ref(), sampleTime, values);
}


void eikonalECGManufacturedVerifier::writeSummary()
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
        mesh_.time().globalPath()
      / "postProcessing"
      / "manufacturedEikonalECGSummary.dat"
    );
    OFstream os(outputFile);

    os << "Manufactured eikonal ECG summary\n";
    os << "samples " << sampleCount_ << "\n";
    os << "dimension " << dimension_ << "D\n";
    os << "qChecks";
    forAll(checkQuadratureOrders_, checkI)
    {
        os << " " << checkQuadratureOrders_[checkI];
    }
    os << "\n";
    os << "qReference " << referenceQuadratureOrder_ << "\n";
    os << "k " << k_ << "\n";
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


void eikonalECGManufacturedVerifier::end()
{
    writeSummary();
}

} // End namespace Foam

// ************************************************************************* //
