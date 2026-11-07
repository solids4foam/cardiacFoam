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

#include "leadFieldECGVerifierBase.H"

#include "ecgModelIO.H"
#include "ecgSolver.H"
#include "verificationUtils.H"
#include "fvMatrices.H"
#include "PstreamReduceOps.H"
#include <limits>

namespace Foam
{

const wordList leadFieldECGVerifierBase::quantityNames_
({
    "leadField-direct",
    "leadFieldStatic-directStatic",
    "leadFieldStatic-exactPoint",
    "leadFieldStatic-exactCentroid",
    "directStatic-exactPoint",
    "directStatic-exactCentroid",
    "localisation"
});


const wordList leadFieldECGVerifierBase::fieldRegionNames_
({
    "phiStatic",
    "phiStatic_heart",
    "phiStatic_torso",
    "phiStatic_interfaceBand",
    "phiStatic_surface"
});


leadFieldECGVerifierBase::leadFieldECGVerifierBase
(
    const electroStateProvider& stateProvider,
    const word& domainName,
    const wordList& electrodeNames,
    const List<vector>& electrodePositions
)
:
    ecgVerificationModel
    (
        stateProvider, domainName, electrodeNames, electrodePositions
    ),
    leadFieldPtr_(nullptr),
    electrodeCells_(),
    electrodeCentroids_(),
    referenceCentroid_(point::zero),
    directPhi_(),
    directStaticPhi_(),
    heartPhi_(),
    outputPtr_(),
    sampleCount_(0),
    exactScale_(0),
    electrodeNorms_
    (
        quantityNames_.size(),
        List<errorNorms>(electrodePositions.size())
    ),
    maxNorms_(quantityNames_.size()),
    relativeGapNorms_(),
    fieldErrors_(fieldRegionNames_.size(), FixedList<scalar, 3>(Zero)),
    summaryWritten_(false),
    heartCells_(),
    interfaceBand_(),
    surfaceCells_()
{}


ecgVerificationModel::Requirements
leadFieldECGVerifierBase::requirements() const
{
    Requirements needs;
    needs.needVm = true;
    return needs;
}


void leadFieldECGVerifierBase::initialise()
{
    // The ECG solver has sampled once before the first record, so its lead
    // fields are loaded
    leadFieldPtr_ = solver().leadField();
    if (!leadFieldPtr_ || !leadFieldPtr_->loaded())
    {
        FatalErrorInFunction
            << "ECG domain " << domainName_ << ": " << type()
            << " verifies the lead-field ECG and needs ecgSolver torsoECG "
            << "with torsoCoupling leadField, whose lead fields are loaded "
            << "at its first solve"
            << exit(FatalError);
    }

    const leadFieldConductor& lfc = leadFieldPtr_->conductor();
    const fvMesh& baseMesh = stateProvider_.baseMesh();

    // Region masks of the phiStatic rows and of the reference checks
    heartCells_ = boolList(baseMesh.nCells(), false);
    for (const label cellI : lfc.conductor().heartCellToBaseCell())
    {
        heartCells_[cellI] = true;
    }
    interfaceBand_ =
        verificationUtils::interfaceBandCells(baseMesh, heartCells_);
    surfaceCells_ = verificationUtils::boundaryAdjacentCells(baseMesh);

    initialiseReference(lfc);

    // Electrode cells, their centroids and the centroid of the pinned cell
    electrodeCells_ = lfc.electrodeCells();
    electrodeCentroids_.setSize(electrodeCells_.size(), point::zero);
    forAll(electrodeCells_, e)
    {
        if (electrodeCells_[e] >= 0)
        {
            electrodeCentroids_[e] = baseMesh.C()[electrodeCells_[e]];
        }
    }
    for (point& c : electrodeCentroids_)
    {
        reduce(c, sumOp<point>());
    }

    if (lfc.conductor().torso().pinsReference())
    {
        const label refCell = lfc.conductor().referenceCell();
        referenceCentroid_ =
            refCell >= 0 ? baseMesh.C()[refCell] : point::zero;
        reduce(referenceCentroid_, sumOp<point>());
    }

    directPhi_.reset(lfc.forwardPotential().ptr());
    directStaticPhi_.reset(lfc.forwardPotential().ptr());
    heartPhi_.reset(lfc.heartPotential().ptr());

    if (Pstream::master())
    {
        wordList columns;
        for (const word& name : electrodeNames_)
        {
            for
            (
                const word& column
              : {
                    "leadField", "direct", "exactPoint", "exactCentroid",
                    "leadFieldStatic", "directStatic"
                }
            )
            {
                columns.append(word(column) + "_" + name);
            }
        }
        outputPtr_ = ecgModelIO::openTimeSeries
        (
            mesh_.time().globalPath()/"postProcessing",
            "manufacturedLeadFieldECG_" + domainName_ + ".dat",
            columns
        );
        outputPtr_->precision(std::numeric_limits<scalar>::max_digits10);
    }

    Info<< type() << ": reference " << referenceName();
    reportParameters(Info);
    Info<< ", electrode cell centroids " << electrodeCentroids_ << endl;
}


scalar leadFieldECGVerifierBase::exactPhi
(
    const point& p,
    const scalar t,
    const label zone
) const
{
    const passiveConductor::torsoSettings& torso =
        leadFieldPtr_->conductor().conductor().torso();

    const scalar phi = referencePhi(p, t, zone);
    const scalar gauge = referencePhi(referenceCentroid_, t, -1);

    // A pinned conductor fixes the constant at the pinned cell; with ground
    // patches the reference is absolute
    return
        torso.pinsReference()
      ? phi - gauge + torso.phiEReferenceValue
      : phi;
}


void leadFieldECGVerifierBase::solveDirect
(
    const volScalarField& Vm,
    volScalarField& phi,
    scalarField& values,
    const volScalarField* extraSource
)
{
    const passiveConductor& conductor = leadFieldPtr_->conductor().conductor();
    const dictionary& solverDict =
        phi.mesh().solverDict(passiveConductor::potentialName);
    const leadFieldConductor::fixedPointControls fixedPoint
    (
        leadFieldConductor::fixedPointControls::read(phi.mesh())
    );
    const scalar directTolerance = fixedPoint.tolerance;
    const label directMaxPasses = fixedPoint.maxPasses;
    tmp<volScalarField> tb(conductor.heartSource(Vm));
    if (extraSource)
    {
        tb.ref().primitiveFieldRef() += extraSource->primitiveField();
    }
    volScalarField& heartPhi = heartPhi_.ref();

    const scalarField& phiI = phi.primitiveField();
    scalarField phiPrev(phiI.size());
    scalar change = GREAT;
    label pass = 0;

    while
    (
        pass < leadFieldConductor::fixedPointControls::minPasses
     || change >= directTolerance
    )
    {
        if (pass == directMaxPasses)
        {
            FatalErrorInFunction
                << "Direct forward solve not converged after " << pass
                << " passes; last relative change " << change
                << " (tolerance " << directTolerance << ")."
                << exit(FatalError);
        }

        phiPrev = phiI;

        // The heart companion takes the base values before each assembly
        conductor.fillHeart(phi, heartPhi);

        tmp<fvScalarMatrix> tEqn(conductor.assemble(phi, heartPhi, &tb()));
        conductor.applyReference
        (
            tEqn.ref(), conductor.torso().phiEReferenceValue
        );
        solve(tEqn.ref(), solverDict);

        ++pass;
        change = gMax(mag(phiI - phiPrev))/max(gMax(mag(phiI)), VSMALL);
    }

    values.setSize(electrodeCells_.size());
    forAll(electrodeCells_, e)
    {
        values[e] = electrodeCells_[e] >= 0 ? phiI[electrodeCells_[e]] : 0;
    }
    Pstream::listCombineReduce(values, plusEqOp<scalar>());
}


void leadFieldECGVerifierBase::recordFieldError(const scalar t)
{
    const fvMesh& baseMesh = stateProvider_.baseMesh();
    const scalarField& phi = directStaticPhi_().primitiveField();
    const vectorField& C = baseMesh.C().primitiveField();
    const scalarField& V = baseMesh.V().field();

    // Per region: sum |e| V, sum e^2 V, sum V, max |e|
    List<FixedList<scalar, 4>> sums(fieldRegionNames_.size(), Zero);
    scalar maxE = -1;
    point maxAt = point::zero;
    forAll(C, cellI)
    {
        const scalar e =
            mag(phi[cellI] - exactPhi(C[cellI], t, cellZone(cellI)));
        const bool inHeart = heartCells_[cellI];
        const bool member[5] =
        {
            true, inHeart, !inHeart, interfaceBand_[cellI], surfaceCells_[cellI]
        };
        for (label r = 0; r < 5; ++r)
        {
            if (member[r])
            {
                sums[r][0] += e*V[cellI];
                sums[r][1] += sqr(e)*V[cellI];
                sums[r][2] += V[cellI];
                sums[r][3] = max(sums[r][3], e);
            }
        }
        if (e > maxE)
        {
            maxE = e;
            maxAt = C[cellI];
        }
    }

    forAll(sums, r)
    {
        for (label k = 0; k < 3; ++k)
        {
            reduce(sums[r][k], sumOp<scalar>());
        }
        reduce(sums[r][3], maxOp<scalar>());
        const scalar vol = max(sums[r][2], VSMALL);
        fieldErrors_[r][0] = max(fieldErrors_[r][0], sums[r][0]/vol);
        fieldErrors_[r][1] =
            max(fieldErrors_[r][1], Foam::sqrt(sums[r][1]/vol));
        fieldErrors_[r][2] = max(fieldErrors_[r][2], sums[r][3]);
    }

    const scalar globalMax = returnReduce(maxE, maxOp<scalar>());
    if (maxE == globalMax)
    {
        Pout<< type() << ": largest phiStatic error " << globalMax
            << " at " << maxAt << endl;
    }
}


void leadFieldECGVerifierBase::record
(
    const List<scalar>& numericValues
)
{
    if (!initialised())
    {
        initialise();
    }

    const volScalarField& Vm = requireVm();
    const scalar t = mesh_.time().value();
    const label nElectrodes = electrodeNames_.size();

    // Direct forward on Vm_h
    scalarField direct;
    solveDirect(Vm, directPhi_.ref(), direct);

    // Lead field and direct forward on the exact Vm at the cell centres
    volScalarField VmExact
    (
        IOobject
        (
            Vm.name(), Vm.instance(), Vm.mesh(),
            IOobject::NO_READ, IOobject::NO_WRITE, IOobject::NO_REGISTER
        ),
        Vm
    );
    const vectorField& centres = Vm.mesh().C().primitiveField();
    forAll(centres, cellI)
    {
        VmExact.primitiveFieldRef()[cellI] = exactVm(centres[cellI], t);
    }
    VmExact.correctBoundaryConditions();

    scalarField leadFieldStatic;
    leadFieldPtr_->evaluate(VmExact, leadFieldStatic);
    const tmp<volScalarField> tS
    (
        staticSource(VmExact, t, leadFieldStatic)
    );
    scalarField directStatic;
    solveDirect
    (
        VmExact,
        directStaticPhi_.ref(),
        directStatic,
        tS.valid() ? &tS() : nullptr
    );

    scalarField exactPoint(nElectrodes);
    scalarField exactCentroid(nElectrodes);
    forAll(exactPoint, e)
    {
        exactPoint[e] = exactPhi(electrodePositions_[e], t);
        exactCentroid[e] = exactPhi(electrodeCentroids_[e], t);
    }

    // Errors per quantity, in quantityNames_ order
    List<scalarField> errors(quantityNames_.size());
    errors[0] = scalarField(numericValues) - direct;
    errors[1] = leadFieldStatic - directStatic;
    errors[2] = leadFieldStatic - exactPoint;
    errors[3] = leadFieldStatic - exactCentroid;
    errors[4] = directStatic - exactPoint;
    errors[5] = directStatic - exactCentroid;
    errors[6] = exactPoint - exactCentroid;

    ++sampleCount_;
    forAll(errors, q)
    {
        forAll(errors[q], e)
        {
            electrodeNorms_[q][e].add(errors[q][e]);
        }
        maxNorms_[q].add(max(mag(errors[q])));
    }
    recordFieldError(t);

    relativeGapNorms_.add
    (
        max(mag(errors[0]))/max(max(mag(direct)), VSMALL)
    );
    exactScale_ = max(exactScale_, max(mag(exactPoint)));

    if (Pstream::master())
    {
        List<scalar> row;
        forAll(electrodeNames_, e)
        {
            row.append(numericValues[e]);
            row.append(direct[e]);
            row.append(exactPoint[e]);
            row.append(exactCentroid[e]);
            row.append(leadFieldStatic[e]);
            row.append(directStatic[e]);
        }
        ecgModelIO::writeRow(outputPtr_.ref(), t, row);
    }

    if (verificationUtils::shouldReportManufacturedErrors(Vm))
    {
        writeSummary();
    }
}


void leadFieldECGVerifierBase::writeSummary()
{
    if (summaryWritten_ || sampleCount_ == 0)
    {
        return;
    }
    summaryWritten_ = true;

    // Reduced on every rank before the master writes
    const label nCells =
        verificationUtils::globalManufacturedCellCount
        (
            stateProvider_.baseMesh()
        );
    const label nHeartCells =
        verificationUtils::globalManufacturedCellCount(mesh_);

    if (!Pstream::master())
    {
        return;
    }

    const fileName outputFile
    (
        mesh_.time().globalPath()/"postProcessing"
      / ("manufacturedLeadFieldECGSummary_" + domainName_ + ".dat")
    );
    OFstream os(outputFile);
    os.precision(std::numeric_limits<scalar>::max_digits10);

    const scalar n = sampleCount_;
    const auto writeRow = [&os, n](const word& name, const errorNorms& norms)
    {
        os  << name << " " << norms.l1/n << " " << Foam::sqrt(norms.l2/n)
            << " " << norms.linf << "\n";
    };

    os  << "# Manufactured lead-field ECG summary\n"
        << "# reference " << referenceName() << "\n";
    writeParameters(os);
    os  << "# cells " << nCells << "\n"
        << "# heartCells " << nHeartCells << "\n"
        << "# samples " << sampleCount_ << "\n"
        << "# time " << mesh_.time().value() << "\n"
        << "# exactScale " << exactScale_ << "\n"
        << "# Errors over the samples: <quantity> is the maximum over the "
        << "electrodes, <quantity>_<electrode> one electrode;\n"
        << "# leadField-direct_relative is that maximum divided by the "
        << "largest |direct|\n"
        << "# row L1 L2 Linf\n";

    writeRow("leadField-direct_relative", relativeGapNorms_);
    forAll(fieldRegionNames_, r)
    {
        os  << fieldRegionNames_[r] << " " << fieldErrors_[r][0] << " "
            << fieldErrors_[r][1] << " " << fieldErrors_[r][2] << "\n";
    }
    forAll(quantityNames_, q)
    {
        writeRow(quantityNames_[q], maxNorms_[q]);
        forAll(electrodeNames_, e)
        {
            writeRow
            (
                quantityNames_[q] + "_" + electrodeNames_[e],
                electrodeNorms_[q][e]
            );
        }
    }
}


void leadFieldECGVerifierBase::end()
{
    writeSummary();
}

} // End namespace Foam

// ************************************************************************* //
