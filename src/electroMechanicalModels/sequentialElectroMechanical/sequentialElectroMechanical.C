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
    firstTimeStep_(true)
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

    Info<< "    Active tension model: "
        << activeTensionModel_->type() << nl
        << "    TaScale (model units -> Pa): " << TaScale_ << nl
        << "    Integration points: " << electro().mesh().nCells() << nl
        << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool sequentialElectroMechanical::evolve()
{
    Info<< "Evolving " << type() << endl;

    // Phase 2 readiness probe: on the first time step report whether the
    // solid objectRegistry holds the fields needed to compute lambda.
    if (firstTimeStep_)
    {
        const bool hasD  =
            solid().mesh().foundObject<volVectorField>("D");
        const bool hasF0 =
            solid().mesh().foundObject<volVectorField>("f0");

        Info<< nl
            << "  [Phase2 probe] solid objectRegistry fields for lambda:" << nl
            << "    D   (displacement)      : "
            << (hasD  ? "FOUND"  : "NOT FOUND") << nl
            << "    f0  (fibre direction)   : "
            << (hasF0 ? "FOUND"  : "NOT FOUND") << nl;

        if (hasD && hasF0)
        {
            Info<< "    -> Both present: Phase 2 lambda feedback is ready." << nl;
        }
        else
        {
            Info<< "    -> Missing fields: lambda will stay at 1.0 "
                << "until Phase 2 is wired." << nl;
        }
        Info<< endl;

        firstTimeStep_ = false;
    }

    electro().evolve();

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

    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace electroMechanicalModels

} // End namespace Foam

// ************************************************************************* //
