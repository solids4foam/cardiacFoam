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

#include "manufacturedShellBidomainVerifier.H"

#include "OFstream.H"
#include "OSspecific.H"
#include "PstreamReduceOps.H"
#include "ionicModel.H"
#include "verificationUtils.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

using namespace verificationUtils;

defineTypeNameAndDebug(manufacturedShellBidomainVerifier, 0);
addToRunTimeSelectionTable
(
    electroVerificationModel,
    manufacturedShellBidomainVerifier,
    dictionary
);

namespace
{

const dimensionSet conductivityDim(-1, -3, 3, 0, 0, 2, 0);

scalar sampleValue
(
    const shellBidomainMMS::SelfTestSample& s,
    const shellBidomainMMS::Params& p
)
{
    const point pt(s.x, s.y, s.z);
    switch (s.fn)
    {
        case 0: return shellBidomainMMS::Fi(pt, p);
        case 1: return shellBidomainMMS::gradFi(pt, p)[s.component];
        case 2: return shellBidomainMMS::Fe(pt, p);
        case 3: return shellBidomainMMS::gradFe(pt, p)[s.component];
        case 4: return shellBidomainMMS::divGiGradFi(pt, p);
        case 5: return shellBidomainMMS::divGeGradFe(pt, p);
        case 6: return shellBidomainMMS::Gi(pt, p)[s.component];
        case 7: return shellBidomainMMS::Ge(pt, p)[s.component];
        case 8: return shellBidomainMMS::g(pt, p);
    }
    FatalErrorInFunction
        << "Unknown self-test function " << s.fn << exit(FatalError);
    return 0;
}

} // End anonymous namespace


void manufacturedShellBidomainVerifier::selfTest()
{
    const shellBidomainMMS::Params paramSets[2] =
    {
        shellBidomainMMS::defaultParams(),
        shellBidomainMMS::alternateParams()
    };

    for (label i = 0; i < shellBidomainMMS::nSelfTestSamples; ++i)
    {
        const shellBidomainMMS::SelfTestSample& s =
            shellBidomainMMS::selfTestSamples[i];
        const scalar value = sampleValue(s, paramSets[s.param]);

        if (!(mag(value - s.value) <= 1e-11*max(scalar(1), mag(s.value))))
        {
            FatalErrorInFunction
                << "Shell bidomain MMS self-test sample " << i
                << " (function " << s.fn << ", parameter set " << s.param
                << ", component " << s.component << "): " << value
                << " instead of " << s.value << exit(FatalError);
        }
    }

    Info<< "shell bidomain MMS self-test: "
        << shellBidomainMMS::nSelfTestSamples << " samples OK" << endl;
}


manufacturedShellBidomainVerifier::manufacturedShellBidomainVerifier
(
    const dictionary& dict
)
:
    electroVerificationModel(dict),
    p_(shellBidomainMMS::defaultParams()),
    nPerDirection_(dict.subDict("verificationModel").get<label>("N")),
    phiEPtr_(nullptr),
    errorsReported_(false),
    phiEReferenceValue_(dict.lookupOrDefault<scalar>("phiEReferenceValue", 0.0)),
    phiEReferencePoint_(dict.get<point>("phiERefPoint")),
    beta_(0),
    betaInitialised_(false)
{
    if (nPerDirection_ <= 0)
    {
        FatalIOErrorInFunction(dict.subDict("verificationModel"))
            << "N must be positive; found " << nPerDirection_ << "."
            << exit(FatalIOError);
    }
    selfTest();
}


void manufacturedShellBidomainVerifier::bindBidomainField(volScalarField& phiE)
{
    phiEPtr_ = &phiE;
}


void manufacturedShellBidomainVerifier::initialiseBeta(const ionicModel& model)
{
    const wordList names = model.constantVariableNames();
    const scalarField values = model.constantVariableValues();
    const label i = names.find("Beta");

    if (i < 0 || names.size() != values.size())
    {
        FatalErrorInFunction
            << typeName << " requires the manufactured ionic constant 'Beta'; "
            << "ionic model " << model.type() << " exposes " << names << "."
            << exit(FatalError);
    }

    beta_ = values[i];
    betaInitialised_ = true;
}


void manufacturedShellBidomainVerifier::validateConfiguration
(
    const ionicModel& model,
    const fvMesh& mesh
) const
{
    if
    (
        model.verificationFamily() != "monodomainFDAManufactured"
     || model.geometricDimension() != 3
    )
    {
        FatalErrorInFunction
            << typeName << " requires the 3D monodomainFDAManufactured ionic "
            << "model. Received " << model.type() << " (family '"
            << model.verificationFamily() << "', dimension "
            << model.geometricDimension() << ")." << exit(FatalError);
    }

    if (!dict().get<Switch>("sealedHeartBoundary"))
    {
        FatalErrorInFunction
            << typeName << " is sealed on both spheres and needs "
            << "sealedHeartBoundary true." << exit(FatalError);
    }

    const scalar rRef = mag(phiEReferencePoint_);
    if (rRef <= p_.r1 || rRef >= p_.r2)
    {
        FatalErrorInFunction
            << typeName << ": phiERefPoint " << phiEReferencePoint_
            << " is not inside the shell " << p_.r1 << " < |x| < " << p_.r2
            << "." << exit(FatalError);
    }

    // Mesh points lie in the closed shell, so the walls are the two spheres
    const scalar tol = 1e-8;
    label nOutside = 0;
    for (const point& pt : mesh.points())
    {
        const scalar r = mag(pt);
        if (r < p_.r1 - tol || r > p_.r2 + tol)
        {
            ++nOutside;
        }
    }
    if (returnReduce(nOutside, sumOp<label>()))
    {
        FatalErrorInFunction
            << typeName << ": " << returnReduce(nOutside, sumOp<label>())
            << " mesh points lie outside the shell " << p_.r1 << " <= |x| <= "
            << p_.r2 << "." << exit(FatalError);
    }
}


void manufacturedShellBidomainVerifier::validateConductivity
(
    const fvMesh& mesh
) const
{
    const vectorField& C = mesh.C().primitiveField();

    const word names[2] = {"conductivityIntracellular", "conductivityExtracellular"};

    for (label i = 0; i < 2; ++i)
    {
        const word& name = names[i];
        const volTensorField& conductivity =
            mesh.lookupObject<volTensorField>(name);
        scalar worst = 0;
        forAll(C, cellI)
        {
            const tensor G
            (
                (
                    i == 0
                  ? shellBidomainMMS::Gi(C[cellI], p_)
                  : shellBidomainMMS::Ge(C[cellI], p_)
                ) & tensor(I)
            );
            worst = max(worst, mag(conductivity[cellI] - G)/mag(G));
        }
        reduce(worst, maxOp<scalar>());
        if (worst > 1e-9)
        {
            FatalErrorInFunction
                << typeName << " needs " << name << " written by "
                << "setManufacturedFields (conductivitySource field); the "
                << "largest relative difference from the manufactured tensor is "
                << worst << "." << exit(FatalError);
        }
    }
}


scalar manufacturedShellBidomainVerifier::phiEShift(const scalar t) const
{
    return
        phiEReferenceValue_
      - Foam::sqrt(1 + t)*shellBidomainMMS::Fe(phiEReferencePoint_, p_);
}


wordList manufacturedShellBidomainVerifier::preProcessFieldNames
(
    const ionicModel&
) const
{
    return wordList({"u1", "u2", "u3"});
}


wordList manufacturedShellBidomainVerifier::requiredPostProcessFieldNames
(
    const ionicModel&
) const
{
    return wordList({"u1", "u2"});
}


bool manufacturedShellBidomainVerifier::shouldPostProcess
(
    const ionicModel&,
    const volScalarField& Vm
) const
{
    return !errorsReported_ && shouldReportManufacturedErrors(Vm);
}


void manufacturedShellBidomainVerifier::preProcess
(
    ionicModel& model,
    volScalarField& Vm,
    PtrList<volScalarField>& fields
)
{
    const wordList names = preProcessFieldNames(model);
    validateConfiguration(model, Vm.mesh());
    initialiseBeta(model);
    validateConductivity(Vm.mesh());

    const label u1Field = requireFieldIndex(names, "u1", "preProcess");
    const label u2Field = requireFieldIndex(names, "u2", "preProcess");
    const label u3Field = requireFieldIndex(names, "u3", "preProcess");

    const vectorField& C = Vm.mesh().C().primitiveField();
    const scalar t = Vm.mesh().time().value();
    const scalar S = Foam::sqrt(1 + t);

    scalarField& VmI = Vm.primitiveFieldRef();
    scalarField& u1I = fields[u1Field].primitiveFieldRef();
    scalarField& u2I = fields[u2Field].primitiveFieldRef();
    scalarField& u3I = fields[u3Field].primitiveFieldRef();

    forAll(C, cellI)
    {
        const scalar g = shellBidomainMMS::g(C[cellI], p_);
        VmI[cellI] =
            S
           *(
                shellBidomainMMS::Fi(C[cellI], p_)
              - shellBidomainMMS::Fe(C[cellI], p_)
            );
        u1I[cellI] = (1 + t)*g + VmI[cellI];
        u2I[cellI] = 1/((1 + t)*Foam::sqrt(g));
        u3I[cellI] = 0;
    }

    Vm.correctBoundaryConditions();
    fields[u1Field].correctBoundaryConditions();
    fields[u2Field].correctBoundaryConditions();
    fields[u3Field].correctBoundaryConditions();

    if (phiEPtr_)
    {
        const scalar shift = phiEShift(t);
        scalarField& phiEI = phiEPtr_->primitiveFieldRef();
        forAll(C, cellI)
        {
            phiEI[cellI] = S*shellBidomainMMS::Fe(C[cellI], p_) + shift;
        }
        phiEPtr_->correctBoundaryConditions();
    }

    model.importFields(Vm, names, fields);
}


void manufacturedShellBidomainVerifier::addManufacturedPdeSource
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
            << "read the ionic constant Beta." << exit(FatalError);
    }

    const vectorField& C = sourceField.mesh().C().primitiveField();
    scalarField& source = sourceField.primitiveFieldRef();
    const scalar S = Foam::sqrt(1 + evaluationTime);

    // The FDA ionic model contributes -Cm dV/dt + beta V/chi, so the source is beta V - div(Gi grad phiI)
    forAll(C, cellI)
    {
        source[cellI] +=
            S*(beta_
              *(
                    shellBidomainMMS::Fi(C[cellI], p_)
                  - shellBidomainMMS::Fe(C[cellI], p_)
                )
          - shellBidomainMMS::divGiGradFi(C[cellI], p_));
    }

    sourceField.correctBoundaryConditions();
}


void manufacturedShellBidomainVerifier::addManufacturedPhiESource
(
    volScalarField& rhs,
    const scalar evaluationTime
) const
{
    const vectorField& C = rhs.mesh().C().primitiveField();
    scalarField& r = rhs.primitiveFieldRef();
    const scalar S = Foam::sqrt(1 + evaluationTime);

    forAll(C, cellI)
    {
        r[cellI] +=
            S
           *(
                shellBidomainMMS::divGiGradFi(C[cellI], p_)
              + shellBidomainMMS::divGeGradFe(C[cellI], p_)
            );
    }
}


void manufacturedShellBidomainVerifier::writeManufacturedFields
(
    const fvMesh& mesh
) const
{
    if (Pstream::parRun())
    {
        FatalErrorInFunction
            << typeName << " writes the manufactured conductivities of the "
            << "undecomposed case; run it in serial." << exit(FatalError);
    }

    const unsigned int oldPrecision = IOstream::defaultPrecision(17);
    const vectorField& C = mesh.C().primitiveField();

    for (label i = 0; i < 2; ++i)
    {
        volSymmTensorField G
        (
            IOobject
            (
                i == 0 ? "ConductivityIntracellular" : "ConductivityExtracellular",
                "0",
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedSymmTensor(conductivityDim, Zero),
            "calculated"
        );

        forAll(C, cellI)
        {
            G[cellI] =
                i == 0
              ? shellBidomainMMS::Gi(C[cellI], p_)
              : shellBidomainMMS::Ge(C[cellI], p_);
        }
        forAll(mesh.boundary(), patchI)
        {
            if (polyPatch::constraintType(mesh.boundaryMesh()[patchI].type()))
            {
                continue;
            }
            G.boundaryFieldRef()[patchI] ==
                G.boundaryField()[patchI].patchInternalField();
        }
        G.write();

        Info<< typeName << ": wrote " << G.name() << " in 0." << endl;
    }

    IOstream::defaultPrecision(oldPrecision);
}


void manufacturedShellBidomainVerifier::postProcess
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
            << typeName << " requires a bound phiE field before postProcess()."
            << exit(FatalError);
    }

    const wordList names = requiredPostProcessFieldNames(model);
    const label u1Idx = requireFieldIndex(names, "u1", "postProcess");
    const label u2Idx = requireFieldIndex(names, "u2", "postProcess");

    const fvMesh& mesh = Vm.mesh();
    const Time& time = mesh.time();
    const vectorField& C = mesh.C().primitiveField();
    const scalar t = time.value();
    const scalar S = Foam::sqrt(1 + t);
    const scalar shift = phiEShift(t);

    scalarField VmExact(C.size()), phiEExact(C.size());
    scalarField u1Exact(C.size()), u2Exact(C.size());
    vectorField gradVmExact(C.size()), gradPhiEExact(C.size());
    forAll(C, cellI)
    {
        const scalar g = shellBidomainMMS::g(C[cellI], p_);
        VmExact[cellI] =
            S
           *(
                shellBidomainMMS::Fi(C[cellI], p_)
              - shellBidomainMMS::Fe(C[cellI], p_)
            );
        phiEExact[cellI] = S*shellBidomainMMS::Fe(C[cellI], p_) + shift;
        u1Exact[cellI] = (1 + t)*g + VmExact[cellI];
        u2Exact[cellI] = 1/((1 + t)*Foam::sqrt(g));
        gradVmExact[cellI] =
            S
           *(
                shellBidomainMMS::gradFi(C[cellI], p_)
              - shellBidomainMMS::gradFe(C[cellI], p_)
            );
        gradPhiEExact[cellI] = S*shellBidomainMMS::gradFe(C[cellI], p_);
    }

    // phiE is fixed up to a constant: the exact field is shifted by the volume mean of the error
    const scalarField& phiEValues = phiEPtr_->primitiveField();
    const scalar phiEGaugeOffset =
        computeVolumeMeanError(mesh, phiEValues, phiEExact);
    scalarField phiEExactGauge(phiEExact);
    phiEExactGauge += phiEGaugeOffset;

    const auto VmNorms = computeNorms(mesh, Vm.primitiveField(), VmExact);
    const auto u1Norms = computeNorms(mesh, fields[u1Idx].primitiveField(), u1Exact);
    const auto u2Norms = computeNorms(mesh, fields[u2Idx].primitiveField(), u2Exact);
    const auto phiEGaugeNorms = computeNorms(mesh, phiEValues, phiEExactGauge);

    const tmp<volVectorField> tgradVm(fvc::grad(Vm));
    const tmp<volVectorField> tgradPhiE(fvc::grad(*phiEPtr_));
    const tmp<volVectorField> tgradReconVm
    (
        exactFieldGradient
        (
            mesh,
            Vm.name(),
            [&](const point& c)
            {
                return
                    S
                   *(
                        shellBidomainMMS::Fi(c, p_)
                      - shellBidomainMMS::Fe(c, p_)
                    );
            }
        )
    );
    const tmp<volVectorField> tgradReconPhiE
    (
        exactFieldGradient
        (
            mesh,
            phiEPtr_->name(),
            [&](const point& c)
            {
                return
                    S*shellBidomainMMS::Fe(c, p_) + shift + phiEGaugeOffset;
            }
        )
    );

    auto splitRows = bulkBoundaryRows
    (
        mesh, "Vm", Vm.primitiveField(), VmExact,
        tgradVm().primitiveField(), tgradReconVm().primitiveField(), gradVmExact
    );
    splitRows.append
    (
        bulkBoundaryRows
        (
            mesh, "phiE_gauge", phiEValues, phiEExactGauge,
            tgradPhiE().primitiveField(), tgradReconPhiE().primitiveField(),
            gradPhiEExact, "phiE"
        )
    );

    const label nCells = returnReduce(mesh.nCells(), sumOp<label>());
    const scalar h = Foam::cbrt(gSum(mesh.V())/max(nCells, label(1)));
    const word cellKind(meshCellKind(mesh));
    const fileName outputDir(time.globalPath()/"postProcessing");
    const fileName outputFile
    (
        outputDir/("3D_" + cellKind + "_" + Foam::name(nPerDirection_) + "_cells.dat")
    );

    if (Pstream::master())
    {
        mkDir(outputDir);
        OFstream os(outputFile);
        os.precision(12);
        os  << "# Bidomain sealed-shell manufactured solution error summary\n"
            << "# mesh " << cellKind << "\n"
            << "# cellsPerDirection " << nPerDirection_ << "\n"
            << "# h " << h << "\n"
            << "# steps " << max(label(0), time.timeIndex()) << "\n"
            << "# dt " << time.deltaTValue() << "\n"
            << "# time " << t << "\n"
            << "# Beta " << beta_ << "\n"
            << "# field L1 L2 Linf\n";
        writeNormRow(os, "Vm", VmNorms);
        writeNormRow(os, "u1", u1Norms);
        writeNormRow(os, "u2", u2Norms);
        writeNormRow(os, "phiE_gauge", phiEGaugeNorms);
        for (const auto& row : splitRows)
        {
            writeNormRow(os, row.first(), row.second());
        }
        Info<< nl << "Bidomain sealed-shell manufactured-solution error summary"
            << " (t = " << t << ") written to " << outputFile << endl;
    }

    errorsReported_ = true;
}

} // End namespace Foam

// ************************************************************************* //
