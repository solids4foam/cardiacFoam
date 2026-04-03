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

#include "biDomainSolver.H"
#include "myocardiumDomainConfig.H"
#include "electrophysicsSystemBuilder.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
namespace electroModels
{

defineTypeNameAndDebug(BiDomainSolver, 0);
addToRunTimeSelectionTable(electroModel, BiDomainSolver, dictionary);

BiDomainSolver::BiDomainSolver(Time& runTime, const word& region)
:
    BiDomainSolver(typeName, runTime, region)
{}

BiDomainSolver::BiDomainSolver
(
    const word& type,
    Time& runTime,
    const word& region
)
:
    electroModel(type, runTime, region),
    phiE_
    (
        IOobject
        (
            "phiE",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedScalar("phiE", dimVoltage, 0.0),
        "zeroGradient"
    ),
    ionicModelPtr_
    (
        ionicModel::New
        (
            ionicProperties(),
            mesh().nCells(),
            runTime.deltaTValue()
        )
    ),
    verificationModelPtr_
    (
        electroVerificationModel::New
        (
            electroProperties()
        )
    ),
    outFields_(),
    postProcessFieldNames_
    (
        verificationModelPtr_.valid()
      ? verificationModelPtr_->requiredPostProcessFieldNames(*ionicModelPtr_)
      : wordList()
    ),
    postProcessFields_()
{
    configureSystem();
}

void BiDomainSolver::configureSystem()
{
    myocardiumDomainConfig config
    {
        outFields_,
        postProcessFieldNames_,
        postProcessFields_,
        *ionicModelPtr_,
        verificationModelPtr_.valid() ? &(*verificationModelPtr_) : nullptr
    };

    electrophysicsSystemBuilder::configureBidomainSystem
    (
        domainSystem_,
        mesh(),
        electroProperties(),
        phiE_,
        config,
        runTime().deltaTValue()
    );

    configureECGSystem();
}

} // End namespace electroModels
} // End namespace Foam

// ************************************************************************* //
