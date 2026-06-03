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

namespace Foam
{

namespace
{

Tuple2<Tuple2<scalar, scalar>, scalar> computeNorms
(
    const scalarField& numeric,
    const scalarField& exact
)
{
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


word dimensionName(const label dimension)
{
    return
        dimension == 1 ? "1D"
      : dimension == 2 ? "2D"
      : dimension == 3 ? "3D"
                       : "unknown";
}

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


autoPtr<manufacturedEikonalVerifier> manufacturedEikonalVerifier::New
(
    const dictionary& electroProperties,
    const fvMesh& mesh,
    const volTensorField& conductivity,
    const dimensionedScalar& chi,
    const dimensionedScalar& Cm,
    const dimensionedScalar& c0,
    const Switch& eikonalAdvectionDiffusionApproach
)
{
    const dictionary* verificationDictPtr =
        electroProperties.findDict("verificationModel");

    if (!verificationDictPtr)
    {
        return autoPtr<manufacturedEikonalVerifier>(nullptr);
    }

    const word modelType =
        verificationDictPtr->lookupOrDefault<word>("type", word::null);

    if (modelType.empty() || modelType == "none")
    {
        return autoPtr<manufacturedEikonalVerifier>(nullptr);
    }

    if (modelType != "manufacturedEikonalVerifier")
    {
        FatalErrorInFunction
            << "Unknown eikonal verificationModel type " << modelType
            << ". Valid type is manufacturedEikonalVerifier."
            << exit(FatalError);
    }

    return autoPtr<manufacturedEikonalVerifier>
    (
        new manufacturedEikonalVerifier
        (
            electroProperties,
            mesh,
            conductivity,
            chi,
            Cm,
            c0,
            eikonalAdvectionDiffusionApproach
        )
    );
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

    scalar minX =  GREAT;
    scalar maxX = -GREAT;
    const vectorField& centres = mesh_.C().primitiveField();

    forAll(centres, cellI)
    {
        minX = min(minX, centres[cellI].x());
        maxX = max(maxX, centres[cellI].x());
    }

    reduce(minX, minOp<scalar>());
    reduce(maxX, maxOp<scalar>());

    const tensor conductivity =
        manufacturedEikonalConstantConductivity(conductivity_);
    const vector gradTau =
        manufacturedEikonalGradTau(conductivity, chi_, Cm_, c0_);

    scalarField& activationValues = activationTime.primitiveFieldRef();
    activationValues = -1.0;

    const scalar tolerance = 1e-8;
    forAll(activationValues, cellI)
    {
        const scalar x = centres[cellI].x();
        if (   Foam::mag(x - minX) <= tolerance
            || Foam::mag(x - maxX) <= tolerance)
        {
            activationValues[cellI] =
                manufacturedEikonalTau(centres[cellI], gradTau);
        }
    }

    volScalarField::Boundary& boundary = activationTime.boundaryFieldRef();
    forAll(boundary, patchI)
    {
        fvPatchScalarField& patchField = boundary[patchI];

        if (patchField.empty())
        {
            continue;
        }

        const vectorField& faceCentres = patchField.patch().Cf();

        forAll(patchField, faceI)
        {
            patchField[faceI] =
                manufacturedEikonalTau(faceCentres[faceI], gradTau);
        }
    }

    activationTime.correctBoundaryConditions();
}


void manufacturedEikonalVerifier::postProcess
(
    const volScalarField& activationTime
)
{
    if (!shouldPostProcess())
    {
        return;
    }

    const tensor conductivity =
        manufacturedEikonalConstantConductivity(conductivity_);
    const vector gradTau =
        manufacturedEikonalGradTau(conductivity, chi_, Cm_, c0_);

    const vectorField& centres = mesh_.C().primitiveField();
    scalarField exact(centres.size(), 0.0);

    forAll(exact, cellI)
    {
        exact[cellI] = manufacturedEikonalTau(centres[cellI], gradTau);
    }

    const auto norms = computeNorms(activationTime.primitiveField(), exact);

    const fileName outputDir(mesh_.time().globalPath()/"postProcessing");
    mkDir(outputDir);
    const fileName outputFile(outputDir/"manufacturedEikonalActivationTime.dat");

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
            << "# dimension " << dimensionName(dimension_) << "\n"
            << "# chi " << chi_ << "\n"
            << "# Cm " << Cm_ << "\n"
            << "# c0 " << c0_ << "\n"
            << "# gradTau " << gradTau << "\n"
            << "# field L1 L2 Linf\n"
            << "activationTime " << norms.first().first() << " "
            << norms.first().second() << " " << norms.second() << "\n";
    }

    errorsReported_ = true;
}

} // End namespace Foam

// ************************************************************************* //
