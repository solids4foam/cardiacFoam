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

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameWithName(electrophysiologyModel, "electrophysiologyModel");
defineDebugSwitch(electrophysiologyModel, 0);

//  Register under the dict keys that users write in electroProperties.
//  Both "monodomainSolver" and "bidomainSolver" map to this single class;
//  the actual domain/kernel is selected by the builder from the Coeffs subdict.
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

Foam::word Foam::electrophysiologyModel::readSolverType
(
    Time& runTime,
    const word& region
)
{
    // Pre-read electroProperties without registering so the base class
    // constructor can subsequently register and own it properly.
    IOdictionary props
    (
        IOobject
        (
            "electroProperties",
            bool(region == dynamicFvMesh::defaultRegion)
          ? fileName(runTime.caseConstant())
          : fileName(runTime.caseConstant()/region),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false  // do not register
        )
    );

    return word(props.lookup("myocardiumSolver"));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::electrophysiologyModel::electrophysiologyModel
(
    Time& runTime,
    const word& region
)
:
    // Pass the user-selected solver type (e.g. "monodomainSolver") to the
    // base class so it reads the correct <type>Coeffs sub-dictionary.
    electroModel(readSolverType(runTime, region), runTime, region),
    ionicModelPtr_(),
    verificationModelPtr_(),
    outFields_(),
    postProcessFields_()
{
    // Assemble MyocardiumDomain (FVM kernel + ionic model + optional
    // verification) or EikonalMyocardiumDomain. The builder reads the
    // solver type from the active Coeffs sub-dict name.
    electrophysicsSystemBuilder::configureMyocardiumDomain
    (
        domainSystem_,
        mesh(),
        electroProperties(),
        outFields_,
        wordList(),         // postProcessFieldNames (filled by verification)
        postProcessFields_,
        ionicModelPtr_,
        verificationModelPtr_,
        runTime.deltaTValue()
    );

    electrophysicsSystemBuilder::configureAdvanceScheme
    (
        domainSystem_,
        electroProperties()
    );

    electrophysicsSystemBuilder::configureConductionDomains
    (
        domainSystem_,
        mesh(),
        electroProperties(),
        runTime.deltaTValue()
    );

    configureECGDomains();
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
