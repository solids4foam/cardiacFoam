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

#include "electrophysiologyModel.H"
#include "addToRunTimeSelectionTable.H"
#include "electrophysicsSystemBuilder.H"
#include "electroVolumeFieldDomain.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameWithName(electrophysiologyModel, "electrophysiologyModel");
defineDebugSwitch(electrophysiologyModel, 0);

//  Register under the dict keys that users write in electroProperties.
//  Both "monodomainSolver" and "bidomainSolver" map to this single class;
//  electroModel::New passes the selected name, which the builder uses to
//  build the myocardium domain and its kernel.
addNamedToRunTimeSelectionTable
(
    electroModel,
    electrophysiologyModel,
    dictionary,
    monodomainSolver
);

addNamedToRunTimeSelectionTable
(
    electroModel,
    electrophysiologyModel,
    dictionary,
    bidomainSolver
);

addNamedToRunTimeSelectionTable
(
    electroModel,
    electrophysiologyModel,
    dictionary,
    eikonalSolver
);

} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::electrophysiologyModel::electrophysiologyModel
(
    const word& solverType,
    Time& runTime,
    const word& region
)
:
    electroModel(solverType, runTime, region)
{
    // The torso coupling, read once: with a bath the bath domain owns phiE
    const torsoCoupling::type coupling = torsoCoupling::read(*this);

    // Assemble MyocardiumDomain (FVM kernel + ionic model + optional
    // verification) or EikonalMyocardiumDomain for the solver type
    electrophysicsSystemBuilder::configureMyocardiumDomain
    (
        domainSystem_,
        mesh(),
        myocardiumSolverType(),
        electroProperties(),
        runTime.deltaTValue(),
        coupling
    );

    readRestartState();

    electrophysicsSystemBuilder::configureAdvanceScheme
    (
        domainSystem_,
        electroProperties()
    );

    electrophysicsSystemBuilder::configureBathPotentialDomain
    (
        domainSystem_,
        mesh(),
        *this,
        coupling
    );

    electrophysicsSystemBuilder::configureConductionDomains
    (
        domainSystem_,
        mesh(),
        electroProperties(),
        runTime.deltaTValue()
    );

    configureECGDomains();

#   include "printElectrophysiologySummary.H"
}


bool Foam::electrophysiologyModel::readRestartState()
{
    return
        domainSystem_.hasMyocardium()
     && domainSystem_.myocardium().readRestartState();
}


void Foam::electrophysiologyModel::writeRestartState() const
{
    if (domainSystem_.hasMyocardium())
    {
        domainSystem_.myocardium().writeRestartState();
    }
}


void Foam::electrophysiologyModel::writeFields(const Time& runTime)
{
    electroModel::writeFields(runTime);
    writeRestartState();
}


void Foam::electrophysiologyModel::setDeltaT(Time& runTime)
{
    if
    (
        domainSystem_.hasMyocardium()
     && domainSystem_.myocardium().applyModelTimeControls(runTime)
    )
    {
        return;
    }

    electroModel::setDeltaT(runTime);
}


// ************************************************************************* //
