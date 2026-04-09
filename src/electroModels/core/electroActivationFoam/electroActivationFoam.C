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

#include "electroActivationFoam.H"
#include "addToRunTimeSelectionTable.H"
#include "electrophysicsSystemBuilder.H"
#include "ionicModel.H"
#include "electroVerificationModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameWithName(electroActivationFoam, "electroActivationFoam");
defineDebugSwitch(electroActivationFoam, 0);

//  Register under the dict keys that users write in electroProperties.
//  Both "monodomainSolver" and "bidomainSolver" map to this single class;
//  the actual FVM kernel is selected by the builder from the Coeffs subdict.
addNamedToRunTimeSelectionTable
(
    electroModel,
    electroActivationFoam,
    dictionary,
    monodomainSolver
);

addNamedToRunTimeSelectionTable
(
    electroModel,
    electroActivationFoam,
    dictionary,
    bidomainSolver
);

} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::word Foam::electroActivationFoam::readSolverType
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

Foam::electroActivationFoam::electroActivationFoam
(
    Time& runTime,
    const word& region
)
:
    // Pass the user-selected solver type (e.g. "monodomainSolver") to the
    // base class so it reads the correct <type>Coeffs sub-dictionary.
    electroModel(readSolverType(runTime, region), runTime, region),
    ionicModelPtr_
    (
        ionicModel::New
        (
            electroProperties(),
            MyocardiumDomain::configuredCellCount(mesh(), electroProperties()),
            runTime.deltaTValue()
        )
    ),
    verificationModelPtr_
    (
        electroVerificationModel::New(electroProperties())
    ),
    outFields_(),
    postProcessFields_()
{
    // Assemble MyocardiumDomain (FVM kernel + ionic model + optional
    // verification).  The builder reads the solver type from the Coeffs
    // sub-dict name, so MonodomainSolver/BidomainSolver are selected here.
    electrophysicsSystemBuilder::configureMyocardiumDomain
    (
        domainSystem_,
        mesh(),
        electroProperties(),
        outFields_,
        wordList(),         // postProcessFieldNames (filled by verification)
        postProcessFields_,
        ionicModelPtr_(),
        verificationModelPtr_.get()
    );

    electrophysicsSystemBuilder::configureAdvanceScheme
    (
        domainSystem_,
        electroProperties()
    );

    electrophysicsSystemBuilder::configureConductionSystemDomain
    (
        domainSystem_,
        mesh(),
        electroProperties(),
        runTime.deltaTValue()
    );

    configureECGSystem();
}


// ************************************************************************* //
