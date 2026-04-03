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
    kTa_
    (
        "kTa",
        dimPressure,
        electroMechanicalProperties()
    ),
    CaiThreshold_
    (
        "CaiThreshold",
        dimless,
        electroMechanicalProperties()
    )
{
    Info<< "    Active tension coupling parameters:" << nl
        << "        kTa = " << kTa_.value() << " Pa/(Cai unit)" << nl
        << "        CaiThreshold = " << CaiThreshold_.value() << nl
        << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool sequentialElectroMechanical::evolve()
{
    Info<< "Evolving " << type() << endl;

    // Evolve the electro model
    electro().evolve();

    // Extract intracellular calcium from the electro model
    const tmp<volScalarField> tCai = electro().couplingField("Cai");
    const scalarField& Cai = tCai().primitiveField();

    // Compute active tension using a simple linear model:
    //   Ta = kTa * max(Cai - CaiThreshold, 0)
    // This is a placeholder that will be replaced by a dedicated
    // runtime-selectable active tension model in the future.
    scalarField& TaI = Ta_.primitiveFieldRef();
    forAll(TaI, cellI)
    {
        TaI[cellI] =
            kTa_.value()
           *max(Cai[cellI] - CaiThreshold_.value(), scalar(0));
    }
    Ta_.correctBoundaryConditions();

    // Evolve the solid model
    solid().evolve();

    // Update total fields at the end of the time-step
    solid().updateTotalFields();

    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace electroMechanicalModels

} // End namespace Foam

// ************************************************************************* //
