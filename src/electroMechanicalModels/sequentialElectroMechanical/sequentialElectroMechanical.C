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

#include "sequentialElectroMechanical.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace electroMechanicalModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(sequentialElectroMechanical, 0);
addToRunTimeSelectionTable
(
    electroMechanicalModel, sequentialElectroMechanical, dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

sequentialElectroMechanical::sequentialElectroMechanical
(
    Time& runTime,
    const word& region
)
:
    electroMechanicalModel(typeName, runTime, region),
    Ta_
    (
        IOobject
        (
            "Ta",
            runTime.timeName(),
            solid().mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        solid().mesh(),
        dimensionedScalar("zero", dimPressure, 0.0),
        "zeroGradient"
    ),
    TaScale_
    (
        electroMechanicalProperties().lookupOrDefault<scalar>("TaScale", 1e3)
    ),
    lambdaField_(electro().mesh().nCells(), 1.0),
    activeTensionModel_
    (
        activeTensionModel::New
        (
            electroMechanicalProperties(),
            electro().mesh().nCells()
        )
    ),
    verificationModelPtr_(),
    activeTensionRequirements_(activeTensionModel_->requirements())
{
    const ElectromechanicalSignalProvider* prov = electro().provider();

    if (prov)
    {
        activeTensionModel_->setElectromechanicalSignalProvider(*prov);
    }

    activeTensionModel_->validateProvider();

    if (solid().mesh().nCells() != electro().mesh().nCells())
    {
        FatalErrorInFunction
            << "sequentialElectroMechanical requires conforming meshes "
            << "(same cell count). Solid has "
            << solid().mesh().nCells() << " cells, electro has "
            << electro().mesh().nCells() << " cells."
            << abort(FatalError);
    }

    if (!solid().mesh().foundObject<volVectorField>("f0"))
    {
        new volVectorField
        (
            IOobject
            (
                "f0",
                runTime.timeName(),
                solid().mesh(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            solid().mesh()
        );

        Info<< "    Registered f0 in solid objectRegistry." << nl << endl;
    }

    if (activeTensionRequirements_.needsLambda)
    {
        if (!solid().mesh().foundObject<volVectorField>("D"))
        {
            FatalErrorInFunction
                << "Active tension model '" << activeTensionModel_->type()
                << "' requires fibre stretch (lambda) but field D "
                << "is not in the solid objectRegistry."
                << abort(FatalError);
        }
        if (!solid().mesh().foundObject<volVectorField>("f0"))
        {
            FatalErrorInFunction
                << "Active tension model '" << activeTensionModel_->type()
                << "' requires fibre stretch (lambda) but field f0 "
                << "is not in the solid objectRegistry."
                << abort(FatalError);
        }
    }

    const bool activeTensionRestarted =
        activeTensionModel_->readRestartState(solid().mesh());

    if (activeTensionRestarted)
    {
        activeTensionModel_->refreshRestartState(solid().mesh());
        scalarField restartTa(Ta_.primitiveField());
        if (activeTensionModel_->restartTension(restartTa))
        {
            if (TaScale_ != 1.0)
            {
                restartTa *= TaScale_;
            }
            Ta_.primitiveFieldRef() = restartTa;
            Ta_.correctBoundaryConditions();
        }
    }

    if
    (
        prov
     && activeTensionRequirements_.needCai
     && !activeTensionRestarted
    )
    {
        const scalar restingCai = prov->signal(0, CouplingSignal::CAI);
        activeTensionModel_->preconditionToRestingState(restingCai);
    }

    if
    (
        electromechanicalVerificationModel::configured
        (
            electroMechanicalProperties()
        )
    )
    {
        verificationModelPtr_ =
            electromechanicalVerificationModel::New
            (
                electroMechanicalProperties()
            );
    }

    if (verificationModelPtr_.valid())
    {
        verificationModelPtr_->initialize
        (
            const_cast<volScalarField&>(electro().Vm()),
            solid().D()
        );
    }

    Info<< "    Active tension model: "
        << activeTensionModel_->type() << nl
        << "    TaScale (model units -> Pa): " << TaScale_ << nl
        << "    Integration points: " << electro().mesh().nCells() << nl
        << endl;
}


void sequentialElectroMechanical::writeFields(const Time& runTime)
{
    electroMechanicalModel::writeFields(runTime);
    electro().writeRestartState();
    activeTensionModel_->writeRestartState(solid().mesh());
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void sequentialElectroMechanical::updateLambda()
{
    const fvMesh& solidMesh = solid().mesh();

    const bool hasD  = solidMesh.foundObject<volVectorField>("D");
    const bool hasF0 = solidMesh.foundObject<volVectorField>("f0");

    if (!hasD || !hasF0)
    {
        if (activeTensionRequirements_.needsLambda)
        {
            FatalErrorInFunction
                << "Active tension model '" << activeTensionModel_->type()
                << "' requires fibre stretch (lambda) but field "
                << (!hasD ? "D" : "f0")
                << " disappeared from the solid objectRegistry at t="
                << runTime().value() << "."
                << abort(FatalError);
        }
        lambdaField_ = 1.0;
        return;
    }

    const volVectorField& D  = solidMesh.lookupObject<volVectorField>("D");
    const volVectorField& f0 = solidMesh.lookupObject<volVectorField>("f0");

    // Deformation gradient F = I + grad(D)^T (total Lagrangian convention,
    // matching solids4foam mechanicalLaw). The fibre stretch follows from
    // lambda^2 = f0 & C & f0 = (F & f0) & (F & f0), i.e. lambda = mag(F & f0).
    const volTensorField gradD(fvc::grad(D));

    forAll(lambdaField_, cellI)
    {
        const tensor F(I + gradD[cellI].T());
        lambdaField_[cellI] = mag(F & f0[cellI]);
    }
}


bool sequentialElectroMechanical::evolve()
{
    Info<< "Evolving " << type() << endl;

    if (verificationModelPtr_.valid())
    {
        verificationModelPtr_->preSolve
        (
            const_cast<volScalarField&>(electro().Vm()),
            solid().D()
        );
    }

    electro().evolve();

    // Update the fibre stretch from the (lagged) solid deformation before
    // evaluating the active tension.
    updateLambda();

    const scalar t  = runTime().value();
    const scalar dt = runTime().deltaT().value();

    scalarField& TaI = Ta_.primitiveFieldRef();

    activeTensionModel_->calculateTension(t, dt, lambdaField_, TaI);

    if (TaScale_ != 1.0)  // skip no-op multiply; 1.0 is exactly representable
    {
        TaI *= TaScale_;
    }

    Ta_.correctBoundaryConditions();

    solid().evolve();
    solid().updateTotalFields();

    if
    (
        verificationModelPtr_.valid()
     && verificationModelPtr_->shouldPostProcess(electro().Vm(), solid().D())
    )
    {
        verificationModelPtr_->postProcess(electro().Vm(), solid().D(), Ta_);
    }

    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace electroMechanicalModels

} // End namespace Foam

// ************************************************************************* //
