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

#include "manufacturedEikonalECGVerifier.H"

#include "DynamicList.H"
#include "PstreamReduceOps.H"
#include "addToRunTimeSelectionTable.H"
#include "ecgModelIO.H"
#include "eikonalVerification/manufacturedEikonalReference.H"
#include "eikonalVerification/manufacturedShellEikonalReference.H"

namespace Foam
{

namespace
{

//- Nodes of the spherical rule of order q on the shell r1 < |x| < r2: q
//  Gauss-Legendre nodes in r (Jacobian r^2), q in mu = cos(theta) and 2q
//  trapezoid nodes in phi. At each node the exact psi and, per electrode,
//  the weight times -(G grad psi).(x - e)/|x - e|^3, so the reference at
//  time t is the sum over the nodes of kernel*U'(t - psi). The radial nodes
//  are dealt round the ranks; the caller sums over them.
void shellEikonalECGNodes
(
    const List<vector>& electrodePositions,
    const shellEikonalMMS::Params& p,
    const label quadratureOrder,
    scalarField& psi,
    List<scalarField>& kernel
)
{
    List<scalar> nodes, weights;
    pseudoECGManufacturedQuadratureRule(quadratureOrder, nodes, weights);

    const label nPhi = 2*quadratureOrder;
    const scalar dPhi = constant::mathematical::twoPi/scalar(nPhi);
    const scalar rHalf = 0.5*(p.r2 - p.r1);
    const scalar rMid = 0.5*(p.r2 + p.r1);

    label nLocal = 0;
    forAll(nodes, i)
    {
        if (i % Pstream::nProcs() == Pstream::myProcNo())
        {
            nLocal += quadratureOrder*nPhi;
        }
    }

    psi.setSize(nLocal);
    kernel.setSize(electrodePositions.size());
    forAll(kernel, electrodeI)
    {
        kernel[electrodeI].setSize(nLocal);
    }

    label nodeI = 0;

    forAll(nodes, i)
    {
        if (i % Pstream::nProcs() != Pstream::myProcNo())
        {
            continue;
        }

        const scalar r = rMid + rHalf*nodes[i];
        const scalar wr = rHalf*weights[i]*r*r;

        forAll(nodes, j)
        {
            const scalar mu = nodes[j];
            const scalar sinTheta = Foam::sqrt(max(scalar(0), 1 - mu*mu));
            const scalar w = wr*weights[j]*dPhi;

            for (label k = 0; k < nPhi; ++k)
            {
                const scalar phi = scalar(k)*dPhi;
                const point pt
                (
                    r*sinTheta*Foam::cos(phi),
                    r*sinTheta*Foam::sin(phi),
                    r*mu
                );

                // grad Vm = -U'(t - psi) grad psi
                const vector GGradPsi
                (
                    shellEikonalMMS::sigma(pt, p)
                  & shellEikonalMMS::gradPsi(pt, p)
                );

                psi[nodeI] = shellEikonalMMS::psi(pt, p);

                forAll(electrodePositions, electrodeI)
                {
                    const vector d(pt - electrodePositions[electrodeI]);
                    const scalar dMag = mag(d);
                    kernel[electrodeI][nodeI] =
                        -w*(GGradPsi & d)/(dMag*dMag*dMag);
                }

                ++nodeI;
            }
        }
    }
}

} // End anonymous namespace


defineTypeNameAndDebug(manufacturedEikonalECGVerifier, 0);
addToRunTimeSelectionTable
(
    ecgVerificationModel,
    manufacturedEikonalECGVerifier,
    dictionary
);


manufacturedEikonalECGVerifier::manufacturedEikonalECGVerifier
(
    const electroStateProvider& stateProvider,
    const word& domainName,
    const dictionary& dict,
    const wordList& electrodeNames,
    const List<vector>& electrodePositions
)
:
    ecgVerificationModel(stateProvider, domainName, electrodeNames, electrodePositions),
    outputPtr_(),
    enabled_(false),
    dimension_(max(label(1), min(mesh_.nGeometricD(), label(3)))),
    domain_("cube"),
    referenceQuadratureOrder_(96),
    checkQuadratureOrders_(),
    k_(vector::zero),
    sineTemplate_(false),
    templateDuration_(0),
    templateRise_(0),
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
    summaryWritten_(false),
    shellCacheValid_(false),
    shellPsi_(),
    shellKernel_()
{
    read(dict);
}


ecgVerificationModel::Requirements
manufacturedEikonalECGVerifier::requirements() const
{
    Requirements needs;
    needs.needActivationTime = true;
    needs.needConductivity = true;
    needs.needChi = true;
    needs.needCm = true;
    needs.needC0 = true;
    return needs;
}


const dictionary& manufacturedEikonalECGVerifier::manufacturedDict
(
    const dictionary& dict
) const
{
    return dict.subDict("verificationModel");
}


void manufacturedEikonalECGVerifier::resizeStorage()
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


void manufacturedEikonalECGVerifier::initialiseOutput()
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
        ecgModelIO::openTimeSeries
        (
            outDir, "manufacturedEikonalECG_" + domainName_ + ".dat", columns
        );
}


bool manufacturedEikonalECGVerifier::read(const dictionary& dict)
{
    const dictionary& cfg = manufacturedDict(dict);

    enabled_ = cfg.lookupOrDefault<Switch>("enabled", true);
    const label meshDimension = max(label(1), min(mesh_.nGeometricD(), label(3)));
    dimension_ = meshDimension;

    if (cfg.found("dimension"))
    {
        const word dimensionName(cfg.get<word>("dimension"));
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

    domain_ = cfg.lookupOrDefault<word>("domain", "cube");

    if (domain_ != "cube" && domain_ != "shell")
    {
        FatalErrorInFunction
            << "Unsupported manufactured eikonal ECG domain '" << domain_
            << "'. Expected cube or shell." << exit(FatalError);
    }

    if (domain_ == "shell" && dimension_ != 3)
    {
        FatalErrorInFunction
            << "Manufactured eikonal ECG 'domain shell' is the 3D shell "
               "field; the dimension resolved to " << dimension_ << "D."
            << exit(FatalError);
    }

    // Electrodes, orders or the template may have changed
    shellCacheValid_ = false;

    k_ = cfg.lookupOrDefault<vector>("k", manufacturedEikonalK(dimension_));
    eikonalECG_templates::readManufacturedTemplate
    (
        cfg, sineTemplate_, templateDuration_, templateRise_
    );
    referenceQuadratureOrder_ =
        cfg.lookupOrDefault<label>("referenceQuadratureOrder", 96);

    if (cfg.found("checkQuadratureOrders"))
    {
        checkQuadratureOrders_ = cfg.get<List<label>>("checkQuadratureOrders");
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

    if (enabled_ && domain_ == "shell")
    {
        // The conductivity, chi and Cm are checked at the first sample
        checkShellElectrodes();
        (void)requireActivationTime();
    }
    else if (enabled_)
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
        << referenceQuadratureOrder_ << ", dimension=" << dimension_;
    if (domain_ == "shell")
    {
        Info<< ", domain=shell";
    }
    Info<< "." << endl;

    return true;
}


void manufacturedEikonalECGVerifier::checkShellElectrodes() const
{
    const scalar r2 = shellEikonalMMS::defaultParams().r2;

    forAll(electrodePositions_, electrodeI)
    {
        if (mag(electrodePositions_[electrodeI]) <= r2*(1 + 1e-6))
        {
            FatalErrorInFunction
                << "Manufactured eikonal ECG 'domain shell' needs every "
                   "electrode outside the outer sphere |x| = " << r2
                << " so the kernel is smooth on the shell; electrode "
                << electrodeNames_[electrodeI] << " at "
                << electrodePositions_[electrodeI] << " has |x| = "
                << mag(electrodePositions_[electrodeI]) << "."
                << exit(FatalError);
        }
    }
}


void manufacturedEikonalECGVerifier::checkShellConfiguration() const
{
    const shellEikonalMMS::Params p(shellEikonalMMS::defaultParams());

    // The header fixes chi and Cm through the stimulus metric
    const scalar chi = requireChi().value();
    const scalar Cm = requireCm().value();

    if (mag(chi - p.chi) > 1e-12 || mag(Cm - p.cm) > 1e-12)
    {
        FatalErrorInFunction
            << "Manufactured eikonal ECG 'domain shell' requires chi = "
            << p.chi << " and Cm = " << p.cm << "; the solver has chi = "
            << chi << ", Cm = " << Cm << "." << exit(FatalError);
    }

    const tensorField& conductivity = requireConductivity().primitiveField();
    const vectorField& C = mesh_.C().primitiveField();
    scalar worst = 0;

    forAll(conductivity, cellI)
    {
        const tensor G(shellEikonalMMS::sigma(C[cellI], p) & tensor(I));
        worst = max(worst, mag(conductivity[cellI] - G)/mag(G));
    }
    reduce(worst, maxOp<scalar>());

    if (worst > 1e-9)
    {
        FatalErrorInFunction
            << "Manufactured eikonal ECG 'domain shell' integrates the shell "
               "conductivity G, but the ECG conductivity differs from G by a "
               "relative " << worst << "; write it with setManufacturedFields."
            << exit(FatalError);
    }
}


void manufacturedEikonalECGVerifier::rebuildShellCache()
{
    checkShellElectrodes();
    checkShellConfiguration();

    // psi scales as 1/c0; the header's family is exact for any c0
    shellEikonalMMS::Params p(shellEikonalMMS::defaultParams());
    p.c0 = requireC0().value();

    const label nChecks = checkQuadratureOrders_.size();
    shellPsi_.setSize(nChecks + 1);
    shellKernel_.setSize(nChecks + 1);

    label nNodes = 0;

    for (label orderI = 0; orderI <= nChecks; ++orderI)
    {
        const label q =
        (
            orderI < nChecks
          ? checkQuadratureOrders_[orderI]
          : referenceQuadratureOrder_
        );

        shellEikonalECGNodes
        (
            electrodePositions_, p, q, shellPsi_[orderI], shellKernel_[orderI]
        );

        nNodes += shellPsi_[orderI].size();
    }

    reduce(nNodes, sumOp<label>());

    Info<< "Manufactured shell eikonal ECG reference: " << nNodes
        << " quadrature nodes over qRef=" << referenceQuadratureOrder_
        << " and " << nChecks << " check orders, c0=" << p.c0 << "."
        << endl;

    shellCacheValid_ = true;
}


void manufacturedEikonalECGVerifier::shellReferences
(
    const scalar sampleTime,
    List<List<scalar>>& checkReferenceValues,
    List<scalar>& referenceValues
) const
{
    const label nChecks = checkQuadratureOrders_.size();
    const label nElectrodes = electrodePositions_.size();

    // This rank's part of every order, summed over the ranks in one reduce
    List<scalar> sums((nChecks + 1)*nElectrodes, scalar(0));

    for (label orderI = 0; orderI <= nChecks; ++orderI)
    {
        const scalarField& psi = shellPsi_[orderI];
        const List<scalarField>& kernel = shellKernel_[orderI];

        forAll(psi, nodeI)
        {
            const scalar dU =
                eikonalECG_templates::evaluateManufacturedTemplateDerivative
                (
                    sampleTime - psi[nodeI],
                    sineTemplate_,
                    templateDuration_,
                    templateRise_
                );

            forAll(kernel, electrodeI)
            {
                sums[orderI*nElectrodes + electrodeI] +=
                    kernel[electrodeI][nodeI]*dU;
            }
        }
    }

    Pstream::listCombineReduce(sums, plusEqOp<scalar>());

    forAll(checkReferenceValues, checkI)
    {
        forAll(checkReferenceValues[checkI], electrodeI)
        {
            checkReferenceValues[checkI][electrodeI] =
                sums[checkI*nElectrodes + electrodeI];
        }
    }

    forAll(referenceValues, electrodeI)
    {
        referenceValues[electrodeI] = sums[nChecks*nElectrodes + electrodeI];
    }
}


void manufacturedEikonalECGVerifier::updateStatistics
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


void manufacturedEikonalECGVerifier::record
(
    const List<scalar>& numericValues
)
{
    record(mesh_.time().value(), numericValues);
}


void manufacturedEikonalECGVerifier::record
(
    scalar sampleTime,
    const List<scalar>& numericValues
)
{
    if (!enabled_)
    {
        return;
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

    if (domain_ == "shell")
    {
        // G(x) and psi are fixed: build the nodes once
        if (!shellCacheValid_)
        {
            rebuildShellCache();
        }

        shellReferences(sampleTime, checkReferenceValues, referenceValues);
    }
    else
    {
        const tensor conductivity =
            manufacturedEikonalConstantConductivity(requireConductivity());

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
                        sineTemplate_,
                        templateDuration_,
                        templateRise_,
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
                    sineTemplate_,
                    templateDuration_,
                    templateRise_,
                    referenceNodes,
                    referenceWeights
                );
        }
    }

    DynamicList<scalar> row;
    row.reserve
    (
        electrodePositions_.size()*(2 + 3*checkQuadratureOrders_.size())
    );

    forAll(electrodePositions_, electrodeI)
    {
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


void manufacturedEikonalECGVerifier::writeSummary()
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
      / ("manufacturedEikonalECGSummary_" + domainName_ + ".dat")
    );
    OFstream os(outputFile);

    os << "Manufactured eikonal ECG summary\n";
    os << "samples " << sampleCount_ << "\n";
    os << "dimension " << dimension_ << "D\n";
    if (domain_ == "shell")
    {
        os << "domain shell\n";
    }
    os << "qChecks";
    forAll(checkQuadratureOrders_, checkI)
    {
        os << " " << checkQuadratureOrders_[checkI];
    }
    os << "\n";
    os << "qReference " << referenceQuadratureOrder_ << "\n";
    if (domain_ != "shell")
    {
        os << "k " << k_ << "\n";
    }
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


void manufacturedEikonalECGVerifier::end()
{
    writeSummary();
}

} // End namespace Foam

// ************************************************************************* //
