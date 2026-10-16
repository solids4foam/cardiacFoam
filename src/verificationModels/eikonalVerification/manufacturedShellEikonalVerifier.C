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

#include "manufacturedShellEikonalVerifier.H"

#include "OFstream.H"
#include "OSspecific.H"
#include "PstreamReduceOps.H"
#include "verificationUtils.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

using namespace verificationUtils;

defineTypeNameAndDebug(manufacturedShellEikonalVerifier, 0);
addToRunTimeSelectionTable
(
    eikonalVerificationModel,
    manufacturedShellEikonalVerifier,
    dictionary
);

namespace
{

const dimensionSet conductivityDim(-1, -3, 3, 0, 0, 2, 0);

point stimulusPoint(const shellEikonalMMS::Params& p)
{
    return point(p.xs[0], p.xs[1], p.xs[2]);
}

scalar sampleValue
(
    const shellEikonalMMS::SelfTestSample& s,
    const shellEikonalMMS::Params& p
)
{
    const point pt(s.x, s.y, s.z);
    switch (s.fn)
    {
        case 0: return shellEikonalMMS::psi(pt, p);
        case 1: return shellEikonalMMS::gradPsi(pt, p)[s.component];
        case 2: return shellEikonalMMS::source(pt, p);
        case 3: return shellEikonalMMS::sigma(pt, p)[s.component];
    }
    FatalErrorInFunction
        << "Unknown self-test function " << s.fn << exit(FatalError);
    return 0;
}

} // End anonymous namespace


void manufacturedShellEikonalVerifier::selfTest()
{
    const shellEikonalMMS::Params paramSets[2] =
    {
        shellEikonalMMS::defaultParams(),
        shellEikonalMMS::alternateParams()
    };

    for (label i = 0; i < shellEikonalMMS::nSelfTestSamples; ++i)
    {
        const shellEikonalMMS::SelfTestSample& s =
            shellEikonalMMS::selfTestSamples[i];
        const scalar value = sampleValue(s, paramSets[s.param]);

        if (!(mag(value - s.value) <= 1e-11*max(scalar(1), mag(s.value))))
        {
            FatalErrorInFunction
                << "Shell eikonal MMS self-test sample " << i
                << " (function " << s.fn << ", parameter set " << s.param
                << ", component " << s.component << "): " << value
                << " instead of " << s.value << exit(FatalError);
        }
    }

    Info<< "shell eikonal MMS self-test: "
        << shellEikonalMMS::nSelfTestSamples << " samples OK" << endl;
}


manufacturedShellEikonalVerifier::manufacturedShellEikonalVerifier
(
    const dictionary& electroProperties,
    const fvMesh& mesh,
    const volTensorField& conductivity,
    const dimensionedScalar& chi,
    const dimensionedScalar& Cm,
    const dimensionedScalar& c0,
    const Switch& eikonalAdvectionDiffusionApproach
)
:
    eikonalVerificationModel(electroProperties),
    mesh_(mesh),
    conductivity_(conductivity),
    p_(shellEikonalMMS::defaultParams()),
    nPerDirection_
    (
        electroProperties.subDict("verificationModel").get<label>("N")
    ),
    initialiseExact_
    (
        electroProperties.subDict("verificationModel")
       .getOrDefault<Switch>("initialiseExact", false)
    ),
    advectionDiffusion_(eikonalAdvectionDiffusionApproach),
    activationTimePtr_(nullptr),
    nPinned_(0),
    errorsReported_(false)
{
    if (nPerDirection_ <= 0)
    {
        FatalIOErrorInFunction(electroProperties.subDict("verificationModel"))
            << "N must be positive; found " << nPerDirection_ << "."
            << exit(FatalIOError);
    }

    p_.rPin =
        electroProperties.subDict("verificationModel")
       .getOrDefault<scalar>("pinnedRadius", p_.rPin);

    selfTest();

    // The family is exact for any front-speed constant
    p_.c0 = c0.value();
    validateConfiguration(chi, Cm);

    Info<< "Enabled manufactured shell eikonal verification: stimulus point "
        << stimulusPoint(p_) << ", pinned radius " << p_.rPin
        << (initialiseExact_ ? ", starting from the exact field" : "")
        << endl;
}


void manufacturedShellEikonalVerifier::validateConfiguration
(
    const dimensionedScalar& chi,
    const dimensionedScalar& Cm
) const
{
    // The reference header fixes chi and Cm through the stimulus metric
    if
    (
        mag(chi.value() - p_.chi) > 1e-12
     || mag(Cm.value() - p_.cm) > 1e-12
    )
    {
        FatalErrorInFunction
            << typeName << " requires chi = " << p_.chi << " and cm = "
            << p_.cm << "; the solver has chi = " << chi.value()
            << ", cm = " << Cm.value() << "." << exit(FatalError);
    }

    // Mesh points lie in the closed shell, so the walls are the two spheres
    const scalar tol = 1e-8;
    label nOutside = 0;
    for (const point& pt : mesh_.points())
    {
        const scalar r = mag(pt);
        if (r < p_.r1 - tol || r > p_.r2 + tol)
        {
            ++nOutside;
        }
    }
    reduce(nOutside, sumOp<label>());
    if (nOutside)
    {
        FatalErrorInFunction
            << typeName << ": " << nOutside
            << " mesh points lie outside the shell " << p_.r1 << " <= |x| <= "
            << p_.r2 << "." << exit(FatalError);
    }
}


void manufacturedShellEikonalVerifier::validateConductivity() const
{
    const vectorField& C = mesh_.C().primitiveField();
    scalar worst = 0;
    forAll(C, cellI)
    {
        const tensor G(shellEikonalMMS::sigma(C[cellI], p_) & tensor(I));
        worst = max(worst, mag(conductivity_[cellI] - G)/mag(G));
    }
    reduce(worst, maxOp<scalar>());
    if (worst > 1e-9)
    {
        FatalErrorInFunction
            << typeName << " needs the conductivity written by "
            << "setManufacturedFields (conductivitySource field); the largest "
            << "relative difference from G is " << worst << "." << exit(FatalError);
    }
}


void manufacturedShellEikonalVerifier::pinCell
(
    scalarField& values,
    const label cellI,
    const point& centre
) const
{
    const scalar psi = shellEikonalMMS::psi(centre, p_);

    // The solver treats a negative activation time as not constrained
    if (psi < 0)
    {
        FatalErrorInFunction
            << "The exact activation time " << psi << " at the pinned cell "
            << "centre " << centre << " is negative, so the solver would not "
            << "constrain the cell." << exit(FatalError);
    }
    values[cellI] = psi;
}


void manufacturedShellEikonalVerifier::applyConstraints
(
    volScalarField& activationTime
) const
{
    validateConductivity();

    const vectorField& C = mesh_.C().primitiveField();
    const point xs(stimulusPoint(p_));

    scalarField& values = activationTime.primitiveFieldRef();
    values = -1.0;

    nPinned_ = 0;
    scalar nearest = GREAT;
    label nearestCell = -1;
    forAll(C, cellI)
    {
        const scalar distance = mag(C[cellI] - xs);
        if (distance < nearest)
        {
            nearest = distance;
            nearestCell = cellI;
        }
        if (distance < p_.rPin)
        {
            pinCell(values, cellI, C[cellI]);
            ++nPinned_;
        }
    }

    // A mesh coarser than rPin pins the cell nearest the stimulus point (the gauge)
    label nPinned = nPinned_;
    reduce(nPinned, sumOp<label>());
    if (nPinned == 0)
    {
        scalar globalNearest = nearest;
        reduce(globalNearest, minOp<scalar>());
        if (nearestCell >= 0 && nearest == globalNearest)
        {
            pinCell(values, nearestCell, C[nearestCell]);
            ++nPinned_;
        }
        nPinned = nPinned_;
        reduce(nPinned, sumOp<label>());
        if (nPinned == 0)
        {
            FatalErrorInFunction
                << "No cell can be pinned near the stimulus point " << xs
                << exit(FatalError);
        }
    }
    nPinned_ = nPinned;

    activationTimePtr_ = &activationTime;
    activationTime.correctBoundaryConditions();
}


tmp<volScalarField> manufacturedShellEikonalVerifier::sourceTerm() const
{
    const vectorField& C = mesh_.C().primitiveField();

    tmp<volScalarField> tSmms
    (
        new volScalarField
        (
            IOobject
            (
                "Smms",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("Smms", dimless, 0.0)
        )
    );

    scalarField& Smms = tSmms.ref().primitiveFieldRef();
    forAll(Smms, cellI)
    {
        Smms[cellI] = shellEikonalMMS::source(C[cellI], p_);
    }

    // Called after the constraints and initial guess are set: the unconstrained cells start exact
    if (initialiseExact_ && activationTimePtr_)
    {
        const point xs(stimulusPoint(p_));
        scalarField& values = activationTimePtr_->primitiveFieldRef();
        forAll(C, cellI)
        {
            if (mag(C[cellI] - xs) >= p_.rPin)
            {
                values[cellI] = shellEikonalMMS::psi(C[cellI], p_);
            }
        }
        activationTimePtr_->correctBoundaryConditions();
    }

    return tSmms;
}


void manufacturedShellEikonalVerifier::postProcess
(
    const volScalarField& activationTime
)
{
    if (!shouldPostProcess())
    {
        return;
    }

    const vectorField& C = mesh_.C().primitiveField();
    const Time& time = mesh_.time();

    scalarField exact(C.size());
    vectorField gradExact(C.size());
    forAll(C, cellI)
    {
        exact[cellI] = shellEikonalMMS::psi(C[cellI], p_);
        gradExact[cellI] = shellEikonalMMS::gradPsi(C[cellI], p_);
    }

    const auto norms =
        computeNorms(mesh_, activationTime.primitiveField(), exact);

    const tmp<volVectorField> tgrad(fvc::grad(activationTime));
    const tmp<volVectorField> tgradRecon
    (
        exactFieldGradient
        (
            mesh_,
            activationTime.name(),
            [&](const point& c) { return shellEikonalMMS::psi(c, p_); }
        )
    );
    const auto splitRows = bulkBoundaryRows
    (
        mesh_, "activationTime", activationTime.primitiveField(), exact,
        tgrad().primitiveField(), tgradRecon().primitiveField(), gradExact
    );

    const label nCells = returnReduce(mesh_.nCells(), sumOp<label>());
    const scalar h = Foam::cbrt(gSum(mesh_.V())/max(nCells, label(1)));
    const word cellKind(meshCellKind(mesh_));
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
        os  << "# Eikonal sealed-shell manufactured solution error summary\n"
            << "# mesh " << cellKind << "\n"
            << "# cellsPerDirection " << nPerDirection_ << "\n"
            << "# h " << h << "\n"
            << "# time " << time.value() << "\n"
            << "# pinnedCells " << nPinned_ << "\n"
            << "# pinnedRadius " << p_.rPin << "\n"
            << "# chi " << p_.chi << "\n"
            << "# Cm " << p_.cm << "\n"
            << "# c0 " << p_.c0 << "\n"
            << "# eikonalAdvectionDiffusionApproach " << advectionDiffusion_ << "\n"
            << "# field L1 L2 Linf\n";
        writeNormRow(os, "activationTime", norms);
        for (const auto& row : splitRows)
        {
            writeNormRow(os, row.first(), row.second());
        }
        Info<< nl << "Eikonal sealed-shell manufactured-solution error summary"
            << " written to " << outputFile << nl;
        writeNormRow(Info, "activationTime", norms);
        for (const auto& row : splitRows)
        {
            writeNormRow(Info, row.first(), row.second());
        }
        Info<< "pinned cells " << nPinned_ << endl;
    }

    errorsReported_ = true;
}


void manufacturedShellEikonalVerifier::writeManufacturedFields
(
    const fvMesh& mesh
) const
{
    if (Pstream::parRun())
    {
        FatalErrorInFunction
            << typeName << " writes the manufactured conductivity of the "
            << "undecomposed case; run it in serial." << exit(FatalError);
    }

    const unsigned int oldPrecision = IOstream::defaultPrecision(17);

    volSymmTensorField G
    (
        IOobject
        (
            "Conductivity",
            "0",
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedSymmTensor(conductivityDim, Zero),
        "calculated"
    );

    const vectorField& C = mesh.C().primitiveField();
    forAll(C, cellI)
    {
        G[cellI] = shellEikonalMMS::sigma(C[cellI], p_);
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

    IOstream::defaultPrecision(oldPrecision);

    Info<< typeName << ": wrote Conductivity in 0." << endl;
}

} // End namespace Foam

// ************************************************************************* //
