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

#include "manufacturedConormalBathBidomainVerifier.H"

#include "OFstream.H"
#include "OSspecific.H"
#include "boundBox.H"
#include "ionicModel.H"
#include "verificationUtils.H"
#include "addToRunTimeSelectionTable.H"

#include <limits>

namespace Foam
{
using namespace verificationUtils;

defineTypeNameAndDebug(manufacturedConormalBathBidomainVerifier, 0);
addToRunTimeSelectionTable
(
    electroVerificationModel,
    manufacturedConormalBathBidomainVerifier,
    dictionary
);

namespace
{

Tuple2<Tuple2<scalar, scalar>, scalar> nanNorms()
{
    const scalar nan = std::numeric_limits<scalar>::quiet_NaN();
    return Tuple2<Tuple2<scalar, scalar>, scalar>
    (
        Tuple2<scalar, scalar>(nan, nan),
        nan
    );
}

} // End anonymous namespace


Foam::scalar Foam::manufacturedConormalBathBidomainVerifier::volumeMean
(
    const fvMesh& m,
    const scalarField& f
)
{
    const scalarField& V = m.V();
    const scalar sumFV = gSum(f*V);
    const scalar sumV = gSum(V);
    return sumV > SMALL ? sumFV/sumV : 0.0;
}


const dictionary&
manufacturedConormalBathBidomainVerifier::verificationDict() const
{
    return dict().subDict("verificationModel");
}


manufacturedConormalBathBidomainVerifier::
manufacturedConormalBathBidomainVerifier
(
    const dictionary& dict
)
:
    electroVerificationModel(dict),
    phiEPtr_(nullptr),
    phiEHeartCellMapPtr_(nullptr),
    errorsReported_(false),
    k_(1.0/Foam::sqrt(2.0)),
    beta_(0.0),
    betaInitialised_(false),
    conductivityValidated_(false),
    coeffsInitialised_(false),
    profile_{complex(1, 0), complex(0, 0), complex(0, 0), 0.0, 0.0},
    coeffs_()
{
    const dictionary& cfg = verificationDict();
    k_ = cfg.lookupOrDefault<scalar>("k", 1.0/Foam::sqrt(2.0));
}


void manufacturedConormalBathBidomainVerifier::initialiseBeta
(
    const ionicModel& model
)
{
    const wordList constantNames = model.constantVariableNames();
    const scalarField constantValues = model.constantVariableValues();

    if (constantNames.size() != constantValues.size())
    {
        FatalErrorInFunction
            << "Ionic model " << model.type()
            << " exposes inconsistent constant metadata: "
            << constantNames.size() << " names but "
            << constantValues.size() << " values."
            << exit(FatalError);
    }

    label betaIndex = -1;
    forAll(constantNames, i)
    {
        if (constantNames[i] == "Beta")
        {
            betaIndex = i;
            break;
        }
    }

    if (betaIndex < 0)
    {
        FatalErrorInFunction
            << "manufacturedConormalBathBidomainVerifier requires the "
            << "manufactured ionic constant 'Beta', but ionic model "
            << model.type() << " exposes constants " << constantNames << "."
            << exit(FatalError);
    }

    beta_ = constantValues[betaIndex];
    betaInitialised_ = true;
}


void manufacturedConormalBathBidomainVerifier::initialiseCoefficients
(
    const tensor& Gi,
    const scalar sigmaB
)
{
    coeffs_ = conormalBathSolveCoefficients(Gi, sigmaB);
    profile_ = conormalProfile
    {
        complex(1, 0),
        coeffs_.v1,
        coeffs_.v2,
        2.0*constant::mathematical::pi,
        0.0
    };
    coeffsInitialised_ = true;
}


void manufacturedConormalBathBidomainVerifier::validateConfiguration
(
    const ionicModel& model
) const
{
    if
    (
        model.verificationFamily() != "bidomainFDAManufactured"
     || model.geometricDimension() != 2
    )
    {
        FatalErrorInFunction
            << "manufacturedConormalBathBidomainVerifier requires the 2D "
            << "bidomainFDAManufactured ionic model. Received ionic model "
            << model.type() << " with verification family '"
            << model.verificationFamily() << "' and geometric dimension "
            << model.geometricDimension() << "."
            << exit(FatalError);
    }
}


void manufacturedConormalBathBidomainVerifier::validateConductivity
(
    const volTensorField& Gi,
    const volTensorField& Ge
) const
{
    const tensor meanI = uniformTensor(Gi.primitiveField(), typeName);
    const tensor meanE = uniformTensor(Ge.primitiveField(), typeName);

    const scalar kappa = Foam::sqrt(2.0) - 1.0;
    const scalar mismatch = mag(meanE - kappa*meanI);

    if (mismatch > 1e-12)
    {
        FatalErrorInFunction
            << "manufacturedConormalBathBidomainVerifier requires "
            << "conductivityExtracellular = kappa*conductivityIntracellular "
            << "with kappa = sqrt(2)-1. ||Ge - kappa*Gi|| = " << mismatch
            << exit(FatalError);
    }
}


void manufacturedConormalBathBidomainVerifier::bindBidomainField
(
    volScalarField& phiE
)
{
    phiEPtr_ = &phiE;
    phiEHeartCellMapPtr_ = nullptr;
}


void manufacturedConormalBathBidomainVerifier::bindBidomainField
(
    volScalarField& phiE,
    const labelUList& heartCellMap
)
{
    phiEPtr_ = &phiE;
    phiEHeartCellMapPtr_ = &heartCellMap;
}


void manufacturedConormalBathBidomainVerifier::unbindBidomainField()
{
    phiEPtr_ = nullptr;
    phiEHeartCellMapPtr_ = nullptr;
}


wordList manufacturedConormalBathBidomainVerifier::preProcessFieldNames
(
    const ionicModel&
) const
{
    return wordList({"u1", "u2", "u3"});
}


wordList
manufacturedConormalBathBidomainVerifier::requiredPostProcessFieldNames
(
    const ionicModel&
) const
{
    return wordList({"u1", "u2", "u3"});
}


bool manufacturedConormalBathBidomainVerifier::shouldPostProcess
(
    const ionicModel&,
    const volScalarField& Vm
) const
{
    return !errorsReported_ && shouldReportManufacturedErrors(Vm);
}


void manufacturedConormalBathBidomainVerifier::preProcess
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
            << type() << " preProcess expected " << names.size()
            << " fields " << names << " but received " << fields.size()
            << exit(FatalError);
    }

    validateConfiguration(model);
    initialiseBeta(model);

    const label u1Field = requireFieldIndex(names, "u1", "preProcess");
    const label u2Field = requireFieldIndex(names, "u2", "preProcess");
    const label u3Field = requireFieldIndex(names, "u3", "preProcess");

    volScalarField& u1m = fields[u1Field];
    volScalarField& u2m = fields[u2Field];
    volScalarField& u3m = fields[u3Field];

    const volTensorField& Gi =
        Vm.mesh().lookupObject<volTensorField>("conductivityIntracellular");
    const volTensorField& Ge =
        Vm.mesh().lookupObject<volTensorField>("conductivityExtracellular");

    if (!conductivityValidated_)
    {
        validateConductivity(Gi, Ge);
        conductivityValidated_ = true;
    }

    if (!coeffsInitialised_)
    {
        const dictionary& cfg = verificationDict();
        const scalar sigmaB =
            cfg.lookupOrDefault<scalar>("sigmaB", 0.0230827345);

        scalar Gxx = 0, Gxy = 0, Gyy = 0, Gzz = 0;
        if (!Gi.primitiveField().empty())
        {
            const tensor& t = Gi.primitiveField()[0];
            Gxx = t.xx();
            Gxy = t.xy();
            Gyy = t.yy();
            Gzz = t.zz();
        }
        reduce(Gxx, maxOp<scalar>());
        reduce(Gxy, maxOp<scalar>());
        reduce(Gyy, maxOp<scalar>());
        reduce(Gzz, maxOp<scalar>());

        initialiseCoefficients(tensor(Gxx, Gxy, 0, Gxy, Gyy, 0, 0, 0, Gzz), sigmaB);
    }

    const vectorField& centres = Vm.mesh().C().primitiveField();
    const scalar t = Vm.mesh().time().value();

    computeConormalBathV(Vm.primitiveFieldRef(), centres, t, profile_);
    computeConormalBathU
    (
        u1m.primitiveFieldRef(),
        u2m.primitiveFieldRef(),
        u3m.primitiveFieldRef(),
        centres,
        t,
        profile_
    );

    Vm.correctBoundaryConditions();
    u1m.correctBoundaryConditions();
    u2m.correctBoundaryConditions();
    u3m.correctBoundaryConditions();

    if (phiEPtr_)
    {
        const vectorField& phiECentres = phiEPtr_->mesh().C().primitiveField();
        computeConormalBathPhiEWhole
        (
            phiEPtr_->primitiveFieldRef(),
            phiECentres,
            t,
            k_,
            profile_,
            coeffs_
        );
        phiEPtr_->correctBoundaryConditions();
    }

    model.importFields(Vm, names, fields);
}


void manufacturedConormalBathBidomainVerifier::addManufacturedPdeSource
(
    volScalarField& sourceField,
    const volTensorField& conductivity,
    const scalar evaluationTime
) const
{
    if (!betaInitialised_ || !coeffsInitialised_)
    {
        FatalErrorInFunction
            << "The manufactured source was requested before preProcess "
            << "initialised the manufactured coefficients."
            << exit(FatalError);
    }

    if (&sourceField.mesh() != &conductivity.mesh())
    {
        FatalErrorInFunction
            << "Source field and conductivity field belong to different "
            << "finite-volume meshes."
            << exit(FatalError);
    }

    const scalar kappa = Foam::sqrt(2.0) - 1.0;
    const scalar k = 1.0/(1.0 + kappa);
    const scalar sqrt1t = Foam::sqrt(1.0 + evaluationTime);

    const vectorField& centres = sourceField.mesh().C().primitiveField();
    const tensorField& Gi = conductivity.primitiveField();

    scalarField source(centres.size(), 0.0);

    forAll(centres, cellI)
    {
        const scalar F = conormalF(profile_, centres[cellI]);
        const scalar divGGradF =
            conormalDivGGradF(profile_, Gi[cellI], centres[cellI]);

        source[cellI] =
            sqrt1t*(beta_*F - (1.0 - k)*divGGradF);
    }

    sourceField.primitiveFieldRef() += source;
    sourceField.correctBoundaryConditions();
}


void manufacturedConormalBathBidomainVerifier::postProcess
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
            << type() << " requires a bound phiE field before postProcess()."
            << exit(FatalError);
    }

    const wordList requiredNames = requiredPostProcessFieldNames(model);
    if (fields.size() != requiredNames.size())
    {
        FatalErrorInFunction
            << type() << " expected " << requiredNames.size()
            << " postProcess fields " << requiredNames
            << " but received " << fields.size()
            << exit(FatalError);
    }

    const label u1Field = requireFieldIndex(requiredNames, "u1", "postProcess");
    const label u2Field = requireFieldIndex(requiredNames, "u2", "postProcess");
    const label u3Field = requireFieldIndex(requiredNames, "u3", "postProcess");

    const fvMesh& mesh = Vm.mesh();
    const Time& time = mesh.time();
    const scalar t = time.value();

    const vectorField& centres = mesh.C().primitiveField();

    scalarField VmExact;
    scalarField phiEHeartExact;
    scalarField phiIExact;
    scalarField u1Exact;
    scalarField u2Exact;
    scalarField u3Exact;

    computeConormalBathV(VmExact, centres, t, profile_);
    computeConormalBathPhiEHeart
    (
        phiEHeartExact, centres, t, k_, profile_, coeffs_
    );
    computeConormalBathPhiIHeart
    (
        phiIExact, centres, t, k_, profile_, coeffs_
    );
    computeConormalBathU
    (
        u1Exact, u2Exact, u3Exact, centres, t, profile_
    );

    const scalarField& VmValues = Vm.primitiveField();
    const scalarField& phiEValues = phiEPtr_->primitiveField();
    const scalarField& u1Values = fields[u1Field].primitiveField();
    const scalarField& u2Values = fields[u2Field].primitiveField();
    const scalarField& u3Values = fields[u3Field].primitiveField();

    const auto VmNorms = computeNorms(mesh, VmValues, VmExact);

    const fvMesh& phiEMesh = phiEPtr_->mesh();
    const vectorField& phiECentres = phiEMesh.C().primitiveField();
    scalarField phiEExact;
    computeConormalBathPhiEWhole
    (
        phiEExact, phiECentres, t, k_, profile_, coeffs_
    );

    scalarField phiEComputed(phiEValues);
    phiEComputed -= volumeMean(phiEMesh, phiEComputed);
    phiEExact -= volumeMean(phiEMesh, phiEExact);

    const auto phiENorms = computeNorms(phiEMesh, phiEComputed, phiEExact);

    auto phiINorms = nanNorms();
    if (phiEValues.size() == VmValues.size() && &phiEMesh == &mesh)
    {
        scalarField phiIValues(VmValues.size(), 0.0);
        forAll(phiIValues, i)
        {
            phiIValues[i] = VmValues[i] + phiEValues[i];
        }

        phiIValues -= volumeMean(mesh, phiIValues);
        phiIExact -= volumeMean(mesh, phiIExact);

        phiINorms = computeNorms(mesh, phiIValues, phiIExact);
    }
    else if (phiEHeartCellMapPtr_)
    {
        const labelUList& heartCellMap = *phiEHeartCellMapPtr_;

        if (heartCellMap.size() != VmValues.size())
        {
            FatalErrorInFunction
                << "Cannot compute manufactured phiI norms: mapped global "
                << "phiE cell map size " << heartCellMap.size()
                << " does not match Vm cell count " << VmValues.size() << "."
                << exit(FatalError);
        }

        scalarField phiIValues(VmValues.size(), 0.0);
        forAll(phiIValues, i)
        {
            const label baseCellI = heartCellMap[i];
            if (baseCellI < 0 || baseCellI >= phiEValues.size())
            {
                FatalErrorInFunction
                    << "Cannot compute manufactured phiI norms: mapped phiE "
                    << "cell index " << baseCellI << " is outside local phiE "
                    << "field size " << phiEValues.size()
                    << " for myocardium cell " << i << "."
                    << exit(FatalError);
            }

            phiIValues[i] = VmValues[i] + phiEValues[baseCellI];
        }

        phiIValues -= volumeMean(mesh, phiIValues);
        phiIExact -= volumeMean(mesh, phiIExact);

        phiINorms = computeNorms(mesh, phiIValues, phiIExact);
    }

    const auto u1Norms = computeNorms(mesh, u1Values, u1Exact);
    const auto u2Norms = computeNorms(mesh, u2Values, u2Exact);
    const auto u3Norms = computeNorms(mesh, u3Values, u3Exact);

    const label totalCells = globalManufacturedCellCount(mesh);
    const label nPerDirection = structuredCellsPerDirection(totalCells, 2);
    const scalar dt = time.deltaTValue();
    const label nSteps = max(label(0), time.timeIndex());

    const fileName outputDir(time.globalPath()/"postProcessing");
    mkDir(outputDir);
    const fileName outputFile =
        outputDir
      / (
            dimensionName(2) + "_" + Foam::name(nPerDirection) + "_cells.dat"
        );

    if (Pstream::master())
    {
        Info<< nl
            << "Conormal bath-bidomain manufactured-solution error summary "
            << "(t = " << t << "):" << nl
            << "-------------------------------------------------" << nl
            << "Field     L1-error       L2-error       Linf-error" << nl
            << "Vm        " << VmNorms.first().first() << "   "
            << VmNorms.first().second() << "   " << VmNorms.second() << nl
            << "phiE      " << phiENorms.first().first() << "   "
            << phiENorms.first().second() << "   " << phiENorms.second() << nl
            << "phiI      " << phiINorms.first().first() << "   "
            << phiINorms.first().second() << "   " << phiINorms.second() << nl
            << "u1        " << u1Norms.first().first() << "   "
            << u1Norms.first().second() << "   " << u1Norms.second() << nl
            << "u2        " << u2Norms.first().first() << "   "
            << u2Norms.first().second() << "   " << u2Norms.second() << nl
            << "u3        " << u3Norms.first().first() << "   "
            << u3Norms.first().second() << "   " << u3Norms.second() << nl
            << "-------------------------------------------------" << endl;

        OFstream os(outputFile);
        os  << "# Bath-bidomain manufactured solution error summary\n"
            << "# dimension " << dimensionName(2) << "\n"
            << "# cellsPerDirection " << nPerDirection << "\n"
            << "# steps " << nSteps << "\n"
            << "# dt " << dt << "\n"
            << "# time " << t << "\n"
            << "# k " << k_ << "\n"
            << "# phiEMesh " << phiEMesh.name() << "\n"
            << "# field L1 L2 Linf\n"
            << "Vm " << VmNorms.first().first() << " "
            << VmNorms.first().second() << " " << VmNorms.second() << "\n"
            << "phiE " << phiENorms.first().first() << " "
            << phiENorms.first().second() << " " << phiENorms.second() << "\n"
            << "phiI " << phiINorms.first().first() << " "
            << phiINorms.first().second() << " " << phiINorms.second() << "\n"
            << "u1 " << u1Norms.first().first() << " "
            << u1Norms.first().second() << " " << u1Norms.second() << "\n"
            << "u2 " << u2Norms.first().first() << " "
            << u2Norms.first().second() << " " << u2Norms.second() << "\n"
            << "u3 " << u3Norms.first().first() << " "
            << u3Norms.first().second() << " " << u3Norms.second() << "\n";
    }

    errorsReported_ = true;
}

} // End namespace Foam

// ************************************************************************* //
