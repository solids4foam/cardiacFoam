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

#include "monoDomainSolver.H"
#include "myocardiumDomainConfig.H"
#include "electrophysicsSystemBuilder.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
namespace electroModels
{

defineTypeNameAndDebug(MonoDomainSolver, 0);
addToRunTimeSelectionTable(electroModel, MonoDomainSolver, dictionary);

MonoDomainSolver::MonoDomainSolver(Time& runTime, const word& region)
:
    MonoDomainSolver(typeName, runTime, region)
{}


MonoDomainSolver::MonoDomainSolver
(
    const word& type,
    Time& runTime,
    const word& region
)
:
    electroModel(type, runTime, region),
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

void MonoDomainSolver::configureSystem()
{
    myocardiumDomainConfig config
    {
        outFields_,
        postProcessFieldNames_,
        postProcessFields_,
        *ionicModelPtr_,
        verificationModelPtr_.valid() ? &(*verificationModelPtr_) : nullptr
    };

    electrophysicsSystemBuilder::configureMonodomainSystem
    (
        domainSystem_,
        mesh(),
        electroProperties(),
        config,
        runTime().deltaTValue()
    );

    configureECGSystem();
}

} // End namespace electroModels
} // End namespace Foam

// ************************************************************************* //
