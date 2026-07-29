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

#include "manufacturedElectromechanicsVerifier.H"

#include "OFstream.H"
#include "OSspecific.H"
#include "PstreamReduceOps.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcGrad.H"
#include "manufacturedElectromechanicsReference.H"
#include "fvc.H"
#include "fvm.H"
#include "verificationUtils.H"

#include <tuple>

namespace Foam
{

using namespace verificationUtils;

defineTypeNameAndDebug(manufacturedElectromechanicsVerifier, 0);
addToRunTimeSelectionTable
(
    electromechanicalVerificationModel,
    manufacturedElectromechanicsVerifier,
    dictionary
);

namespace
{

void computeNumericalLambda
(
    scalarField& lambda,
    const volVectorField& D
)
{
    const fvMesh& mesh = D.mesh();

    if (!mesh.foundObject<volVectorField>("f0"))
    {
        FatalErrorInFunction
            << "manufacturedElectromechanicsVerifier requires 'f0' on the "
            << mesh.name() << " mesh to evaluate lambda."
            << exit(FatalError);
    }

    const volVectorField& f0 = mesh.lookupObject<volVectorField>("f0");
    const volTensorField gradD(fvc::grad(D));

    lambda.setSize(D.size());

    forAll(lambda, cellI)
    {
        const tensor F(I + gradD[cellI].T());
        lambda[cellI] = mag(F & f0[cellI]);
    }
}

} // End anonymous namespace


manufacturedElectromechanicsVerifier::manufacturedElectromechanicsVerifier
(
    const dictionary& dict
)
:
    electromechanicalVerificationModel(dict),
    amplitude_(dict.lookupOrDefault<vector>("amplitude", vector(0.02, 0.02, 0.02))),
    Tmax_(1.0),
    V0_(1.0),
    gamma_(1.0),
    TaScale_(1.0),
    initializeFields_(dict.lookupOrDefault<Switch>("initializeFields", true)),
    errorsReported_(false)
{
    const dictionary& parentDict = dict.parent();
    if (parentDict.found("constants"))
    {
        const dictionary& constDict = parentDict.subDict("constants");
        Tmax_ = constDict.lookupOrDefault<scalar>("Tmax", dict.lookupOrDefault<scalar>("Tmax", 1.0));
        V0_ = constDict.lookupOrDefault<scalar>("V0", dict.lookupOrDefault<scalar>("V0", 1.0));
        gamma_ = constDict.lookupOrDefault<scalar>("gamma", dict.lookupOrDefault<scalar>("gamma", 1.0));
    }
    else
    {
        Tmax_ = dict.lookupOrDefault<scalar>("Tmax", 1.0);
        V0_ = dict.lookupOrDefault<scalar>("V0", 1.0);
        gamma_ = dict.lookupOrDefault<scalar>("gamma", 1.0);
    }
    TaScale_ = parentDict.lookupOrDefault<scalar>("TaScale", dict.lookupOrDefault<scalar>("TaScale", 1e3));

    const word dimStr = parentDict.lookupOrDefault<word>("dimension", dict.lookupOrDefault<word>("dimension", "3D"));
    dimension_ = (dimStr == "1D") ? 1 : ((dimStr == "2D") ? 2 : 3);


    if (V0_ <= SMALL)
    {
        FatalErrorInFunction
            << "manufacturedElectromechanicsVerifier requires V0 > 0. "
            << "Current value: " << V0_
            << exit(FatalError);
    }
}


void manufacturedElectromechanicsVerifier::setExactFields
(
    volScalarField& Vm,
    volVectorField& D
) const
{
    const scalar t = Vm.mesh().time().value();

    scalarField VmExact;
    vectorField DExact;

    computeManufacturedElectromechanicsVm
    (
        VmExact,
        Vm.mesh().C().primitiveField(),
        t,
        dimension_
    );

    computeManufacturedElectromechanicsD
    (
        DExact,
        D.mesh().C().primitiveField(),
        t,
        amplitude_
    );

    Vm.primitiveFieldRef() = VmExact;
    D.primitiveFieldRef() = DExact;

    Vm.correctBoundaryConditions();
    D.correctBoundaryConditions();
}


void manufacturedElectromechanicsVerifier::initialize
(
    volScalarField& Vm,
    volVectorField& D
)
{
    if (!initializeFields_)
    {
        return;
    }

    Info<< "Initializing manufactured electromechanics Vm and D fields."
        << nl << endl;

    setExactFields(Vm, D);

    // Initialise D.oldTime() to ensure the solids4foam d2dt2 solver computes
    // the correct initial velocity and acceleration, avoiding a massive initial shock
    if (D.nOldTimes() == 0)
    {
        D.storeOldTime();
    }

    const scalar t = Vm.mesh().time().value();
    const scalar dt = Vm.mesh().time().deltaTValue();
    vectorField DOldExact;
    computeManufacturedElectromechanicsD
    (
        DOldExact,
        D.mesh().C().primitiveField(),
        t - dt,
        amplitude_
    );

    D.oldTime().primitiveFieldRef() = DOldExact;
    D.oldTime().correctBoundaryConditions();
}


void manufacturedElectromechanicsVerifier::preSolve
(
    volScalarField&,
    volVectorField&
)
{
    // Deliberately empty. This hook runs once per solve, so writing the
    // manufactured fields here would overwrite the numerical solution with the
    // exact one on every step and drive the reported errors to round-off. The
    // manufactured fields are imposed once, as an initial condition, by
    // initializeFields; anything beyond that would invalidate the measurement
    // this class exists to make.
}


bool manufacturedElectromechanicsVerifier::shouldPostProcess
(
    const volScalarField& Vm,
    const volVectorField&
) const
{
    return !errorsReported_ && shouldReportManufacturedErrors(Vm);
}


void manufacturedElectromechanicsVerifier::postProcess
(
    const volScalarField& Vm,
    const volVectorField& D,
    const volScalarField& Ta
)
{
    if (!shouldPostProcess(Vm, D))
    {
        return;
    }

    const scalar t = Vm.mesh().time().value();

    scalarField VmExact;
    vectorField DExact;
    scalarField lambdaExact;
    scalarField lambdaNum;
    scalarField TaExact;

    computeManufacturedElectromechanicsVm
    (
        VmExact,
        Vm.mesh().C().primitiveField(),
        t,
        dimension_
    );

    computeManufacturedElectromechanicsD
    (
        DExact,
        D.mesh().C().primitiveField(),
        t,
        amplitude_
    );

    computeManufacturedElectromechanicsTa
    (
        TaExact,
        Ta.mesh().C().primitiveField(),
        t,
        amplitude_,
        Tmax_*TaScale_,
        V0_,
        gamma_,
        dimension_
    );

    computeManufacturedElectromechanicsLambda
    (
        lambdaExact,
        D.mesh().C().primitiveField(),
        t,
        amplitude_
    );

    computeNumericalLambda(lambdaNum, D);

    const auto VmNorms = computeNorms(Vm.primitiveField(), VmExact);
    const auto DNorms = computeNorms(D.primitiveField(), DExact);
    const auto lambdaNorms = computeNorms(lambdaNum, lambdaExact);
    const auto TaNorms = computeNorms(Ta.primitiveField(), TaExact);

    if (Pstream::master())
    {
        const fileName outputDir(Vm.mesh().time().globalPath()/"postProcessing");
        mkDir(outputDir);

        OFstream out(outputDir/"manufacturedElectromechanicsSummary.dat");

        Info<< nl
            << "Manufactured electromechanics error summary (t = "
            << t << "):" << nl
            << "Field     L1-error       L2-error       Linf-error" << nl
            << "Vm     " << VmNorms.first().first() << "   " << VmNorms.first().second() << "   " << VmNorms.second() << nl
            << "D      " << DNorms.first().first() << "   " << DNorms.first().second() << "   " << DNorms.second() << nl
            << "lambda " << lambdaNorms.first().first() << "   " << lambdaNorms.first().second() << "   " << lambdaNorms.second() << nl
            << "Ta     " << TaNorms.first().first() << "   " << TaNorms.first().second() << "   " << TaNorms.second() << nl
            << endl;

        out << "Manufactured electromechanics error summary (t = "
            << t << "):\n";
        out << "Field     L1-error       L2-error       Linf-error\n";
        out << "Vm     " << VmNorms.first().first() << "   " << VmNorms.first().second() << "   " << VmNorms.second() << "\n";
        out << "D      " << DNorms.first().first() << "   " << DNorms.first().second() << "   " << DNorms.second() << "\n";
        out << "lambda " << lambdaNorms.first().first() << "   " << lambdaNorms.first().second() << "   " << lambdaNorms.second() << "\n";
        out << "Ta     " << TaNorms.first().first() << "   " << TaNorms.first().second() << "   " << TaNorms.second() << "\n";
    }

    errorsReported_ = true;
}

} // End namespace Foam

// ************************************************************************* //
