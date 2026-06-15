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

#include "manufacturedEikonalVerifier.H"

#include "OFstream.H"
#include "OSspecific.H"
#include "PstreamReduceOps.H"
#include "eikonalVerification/manufacturedEikonalReference.H"
#include "addToRunTimeSelectionTable.H"

#include "verificationUtils.H"

namespace Foam
{
using namespace verificationUtils;

defineTypeNameAndDebug(manufacturedEikonalVerifier, 0);

addToRunTimeSelectionTable
(
    eikonalVerificationModel,
    manufacturedEikonalVerifier,
    dictionary
);

namespace
{


} // End anonymous namespace


manufacturedEikonalVerifier::manufacturedEikonalVerifier
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
    chi_(chi.value()),
    Cm_(Cm.value()),
    c0_(c0.value()),
    enabled_(true),
    dimension_(max(label(1), min(mesh.nGeometricD(), label(3)))),
    errorsReported_(false)
{
    const dictionary& cfg = electroProperties.subDict("verificationModel");
    enabled_ = cfg.lookupOrDefault<Switch>("enabled", true);

    if (!enabled_)
    {
        return;
    }

    validateManufacturedEikonalUnitDomain(mesh_, dimension_);
    (void)manufacturedEikonalConstantConductivity(conductivity_);

    Info<< "Enabled manufactured eikonal activation-time verification."
        << endl;
}



tmp<volScalarField> manufacturedEikonalVerifier::sourceTerm() const
{
    const tensor conductivity =
        manufacturedEikonalConstantConductivity(conductivity_);
    const vector k = manufacturedEikonalK(dimension_);

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

    volScalarField& Smms = tSmms.ref();
    const vectorField& centres = mesh_.C().primitiveField();
    scalarField& SmmsValues = Smms.primitiveFieldRef();

    forAll(SmmsValues, cellI)
    {
        SmmsValues[cellI] = manufacturedEikonalSourceTerm
        (
            centres[cellI],
            conductivity,
            chi_,
            Cm_,
            c0_,
            k
        );
    }

    return tSmms;
}


void manufacturedEikonalVerifier::applyConstraints
(
    volScalarField& activationTime
) const
{
    if (!enabled_)
    {
        return;
    }

    scalar minX = GREAT, minY = GREAT, minZ = GREAT;
    const vectorField& centres = mesh_.C().primitiveField();

    forAll(centres, cellI)
    {
        minX = min(minX, centres[cellI].x());
        minY = min(minY, centres[cellI].y());
        minZ = min(minZ, centres[cellI].z());
    }

    reduce(minX, minOp<scalar>());
    reduce(minY, minOp<scalar>());
    reduce(minZ, minOp<scalar>());

    const vector k = manufacturedEikonalK(dimension_);

    scalarField& activationValues = activationTime.primitiveFieldRef();
    activationValues = -1.0;

    const scalar tolerance = 1e-8;
    forAll(activationValues, cellI)
    {
        const vector& c = centres[cellI];
        bool constrain = false;

        if (dimension_ >= 1 && Foam::mag(c.x() - minX) <= tolerance) constrain = true;
        if (dimension_ >= 2 && Foam::mag(c.y() - minY) <= tolerance) constrain = true;
        if (dimension_ >= 3 && Foam::mag(c.z() - minZ) <= tolerance) constrain = true;

        if (constrain)
        {
            activationValues[cellI] = manufacturedEikonalTau(c, k);
        }
    }

    volScalarField::Boundary& boundary = activationTime.boundaryFieldRef();
    forAll(boundary, patchI)
    {
        fvPatchScalarField& patchField = boundary[patchI];

        if (patchField.empty() || patchField.type() == "empty")
        {
            continue;
        }

        const vectorField& faceCentres = patchField.patch().Cf();

        forAll(patchField, faceI)
        {
            patchField[faceI] = manufacturedEikonalTau(faceCentres[faceI], k);
        }
    }

    activationTime.correctBoundaryConditions();
}


void manufacturedEikonalVerifier::postProcess
(
    const volScalarField& activationTime
)
{
    using namespace verificationUtils;

    if (!shouldPostProcess())
    {
        return;
    }

    const vector k = manufacturedEikonalK(dimension_);

    const vectorField& centres = mesh_.C().primitiveField();
    scalarField exact(centres.size(), 0.0);

    forAll(exact, cellI)
    {
        exact[cellI] = manufacturedEikonalTau(centres[cellI], k);
    }

    const auto norms = computeNorms(activationTime.primitiveField(), exact);

    const fileName outputDir(mesh_.time().globalPath()/"postProcessing");
    mkDir(outputDir);
    const fileName outputFile(outputDir/"manufacturedEikonalActivationTime.dat");

    label globalCells = exact.size();
    reduce(globalCells, sumOp<label>());

    if (Pstream::master())
    {
        Info<< nl
            << "Eikonal manufactured activation-time error summary:" << nl
            << "-------------------------------------------------" << nl
            << "Field           L1-error       L2-error       Linf-error" << nl
            << "activationTime  " << norms.first().first() << "   "
            << norms.first().second() << "   " << norms.second() << nl
            << "-------------------------------------------------" << endl;

        OFstream os(outputFile);
        os  << "# Eikonal manufactured activation-time summary\n"
            << "Number of cells = " << globalCells << "\n"
            << "# dimension " << dimensionName(dimension_) << "\n"
            << "# chi " << chi_ << "\n"
            << "# Cm " << Cm_ << "\n"
            << "# c0 " << c0_ << "\n"
            << "# k " << k << "\n"
            << "# field L1 L2 Linf\n"
            << "activationTime " << norms.first().first() << " "
            << norms.first().second() << " " << norms.second() << "\n";
    }

    errorsReported_ = true;
}

} // End namespace Foam

// ************************************************************************* //
