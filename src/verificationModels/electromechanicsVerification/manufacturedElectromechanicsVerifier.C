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
#include "manufacturedElectromechanicsReference.H"

#include <tuple>

namespace Foam
{

defineTypeNameAndDebug(manufacturedElectromechanicsVerifier, 0);
addToRunTimeSelectionTable
(
    electromechanicalVerificationModel,
    manufacturedElectromechanicsVerifier,
    dictionary
);

namespace
{

template<class FieldType1, class FieldType2>
std::tuple<scalar, scalar, scalar> errorNorms
(
    const FieldType1& num,
    const FieldType2& exact
)
{
    scalar sumAbs = 0.0;
    scalar sumSq = 0.0;
    scalar maxAbs = 0.0;

    forAll(num, i)
    {
        const scalar diff = Foam::mag(num[i] - exact[i]);
        sumAbs += diff;
        sumSq += diff*diff;
        maxAbs = max(maxAbs, diff);
    }

    reduce(maxAbs, maxOp<scalar>());
    reduce(sumAbs, sumOp<scalar>());
    reduce(sumSq, sumOp<scalar>());

    label n = num.size();
    reduce(n, sumOp<label>());

    return std::tuple<scalar, scalar, scalar>
    (
        sumAbs/scalar(n),
        Foam::sqrt(sumSq/scalar(n)),
        maxAbs
    );
}


bool finalTimeReached(const Time& time)
{
    const scalar t = time.value();
    const scalar dt = time.deltaTValue();
    const scalar endTime = time.endTime().value();

    return t + 0.5*dt >= endTime;
}

} // End anonymous namespace


manufacturedElectromechanicsVerifier::manufacturedElectromechanicsVerifier
(
    const dictionary& dict
)
:
    electromechanicalVerificationModel(dict),
    uAmplitude_(dict.lookupOrDefault<scalar>("uAmplitude", 0.02)),
    vAmplitude_(dict.lookupOrDefault<scalar>("vAmplitude", 0.02)),
    wAmplitude_(dict.lookupOrDefault<scalar>("wAmplitude", 0.02)),
    Tmax_(dict.lookupOrDefault<scalar>("Tmax", 1.0)),
    V0_(dict.lookupOrDefault<scalar>("V0", 1.0)),
    gamma_(dict.lookupOrDefault<scalar>("gamma", 1.0)),
    TaScale_(dict.lookupOrDefault<scalar>("TaScale", 1e3)),
    initializeFields_(dict.lookupOrDefault<Switch>("initializeFields", true)),
    enforceExactFields_(dict.lookupOrDefault<Switch>("enforceExactFields", false)),
    errorsReported_(false)
{
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
        t
    );

    computeManufacturedElectromechanicsD
    (
        DExact,
        D.mesh().C().primitiveField(),
        t,
        uAmplitude_,
        vAmplitude_,
        wAmplitude_
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
}


void manufacturedElectromechanicsVerifier::preSolve
(
    volScalarField& Vm,
    volVectorField& D
)
{
    if (enforceExactFields_)
    {
        setExactFields(Vm, D);
    }
}


bool manufacturedElectromechanicsVerifier::shouldPostProcess
(
    const volScalarField& Vm,
    const volVectorField&
) const
{
    return !errorsReported_ && finalTimeReached(Vm.mesh().time());
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
    scalarField TaExact;

    computeManufacturedElectromechanicsVm
    (
        VmExact,
        Vm.mesh().C().primitiveField(),
        t
    );

    computeManufacturedElectromechanicsD
    (
        DExact,
        D.mesh().C().primitiveField(),
        t,
        uAmplitude_,
        vAmplitude_,
        wAmplitude_
    );

    computeManufacturedElectromechanicsTa
    (
        TaExact,
        Ta.mesh().C().primitiveField(),
        t,
        uAmplitude_,
        wAmplitude_,
        Tmax_*TaScale_,
        V0_,
        gamma_
    );

    const auto [L1Vm, L2Vm, LinfVm] =
        errorNorms(Vm.primitiveField(), VmExact);
    const auto [L1D, L2D, LinfD] =
        errorNorms(D.primitiveField(), DExact);
    const auto [L1Ta, L2Ta, LinfTa] =
        errorNorms(Ta.primitiveField(), TaExact);

    if (Pstream::master())
    {
        const fileName outputDir(Vm.mesh().time().path()/"postProcessing");
        mkDir(outputDir);

        OFstream out(outputDir/"manufacturedElectromechanicsSummary.dat");

        Info<< nl
            << "Manufactured electromechanics error summary (t = "
            << t << "):" << nl
            << "Field     L1-error       L2-error       Linf-error" << nl
            << "Vm     " << L1Vm << "   " << L2Vm << "   " << LinfVm << nl
            << "D      " << L1D << "   " << L2D << "   " << LinfD << nl
            << "Ta     " << L1Ta << "   " << L2Ta << "   " << LinfTa << nl
            << endl;

        out << "Manufactured electromechanics error summary (t = "
            << t << "):\n";
        out << "Field     L1-error       L2-error       Linf-error\n";
        out << "Vm     " << L1Vm << "   " << L2Vm << "   " << LinfVm << "\n";
        out << "D      " << L1D << "   " << L2D << "   " << LinfD << "\n";
        out << "Ta     " << L1Ta << "   " << L2Ta << "   " << LinfTa << "\n";
    }

    errorsReported_ = true;
}

} // End namespace Foam

// ************************************************************************* //
