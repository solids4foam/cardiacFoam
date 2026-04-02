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
    electroMechanicalModel(typeName, runTime, region)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool sequentialElectroMechanical::evolve()
{
    Info<< "Evolving " << type() << endl;

    // Evolve the electro model
    electro().evolve();

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
