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

#include "ecgDomain.H"
#include "Switch.H"
#include "ecgModelIO.H"
#include "ecgVerificationModel.H"
#include "ecgElectrodes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(ecgDomain, 0);


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void ecgDomain::readElectrodes
(
    const dictionary& dict,
    const dictionary* inheritedElectrodePositionsPtr
)
{
    const ecgElectrodes electrodes
    (
        readEcgElectrodes(dict, inheritedElectrodePositionsPtr)
    );
    electrodeNames_ = electrodes.names;
    electrodePositions_ = electrodes.positions;

    Info<< "ECG electrodes (" << electrodeNames_.size() << "):" << nl;
    forAll(electrodeNames_, i)
    {
        Info<< "  " << electrodeNames_[i]
            << "  @  " << electrodePositions_[i] << nl;
    }
    Info<< endl;
}


const volScalarField& ecgDomain::Vm() const
{
    const volScalarField* VmPtr = stateProvider_.VmPtr();

    if (!VmPtr)
    {
        FatalErrorInFunction
            << "ECG model requires a transmembrane voltage field, "
            << "but the selected electroStateProvider does not expose Vm."
            << exit(FatalError);
    }

    return *VmPtr;
}


const volScalarField& ecgDomain::activationTime() const
{
    const volScalarField* activationTimePtr = stateProvider_.activationTimePtr();

    if (!activationTimePtr)
    {
        FatalErrorInFunction
            << "ECG model requires an activation-time field, "
            << "but the selected electroStateProvider does not expose one."
            << exit(FatalError);
    }

    return *activationTimePtr;
}


const volTensorField& ecgDomain::conductivity() const
{
    const volTensorField* conductivityPtr = stateProvider_.conductivityPtr();

    if (!conductivityPtr)
    {
        FatalErrorInFunction
            << "ECG model requires a conductivity tensor field, "
            << "but the selected electroStateProvider does not expose one."
            << exit(FatalError);
    }

    return *conductivityPtr;
}


const volScalarField& ecgDomain::phiE() const
{
    const volScalarField* phiEPtr = stateProvider_.phiEPtr();

    if (!phiEPtr)
    {
        FatalErrorInFunction
            << "ECG model requires an extracellular potential field, "
            << "but the selected electroStateProvider does not expose phiE."
            << exit(FatalError);
    }

    return *phiEPtr;
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

ecgDomain::ecgDomain
(
    const electroStateProvider& stateProvider,
    autoPtr<ecgSolver>&& solver,
    const dictionary& dict,
    const word& domainName,
    const dictionary* inheritedElectrodePositionsPtr
)
:
    stateProvider_(stateProvider),
    mesh_(stateProvider.mesh()),
    domainName_(domainName),
    solverType_(solver->type()),
    outputPtr_(),
    solverPtr_(std::move(solver)),
    verificationModelPtr_(),
    numericValues_(),
    writeEveryTimeStep_(dict.lookupOrDefault<Switch>("writeEveryTimeStep", false)),
    electrodeNames_(),
    electrodePositions_()
{
    Info<< domainName << " ECG solver: " << solverType_
        << " on mesh '" << mesh_.name() << "'" << nl
        << "  state source: electrophysiology provider"
        << nl << endl;

    solverPtr_->checkState(*this);

    readElectrodes(dict, inheritedElectrodePositionsPtr);

    const fileName outDir(mesh_.time().globalPath() / "postProcessing");
    if (!solverPtr_->writesOwnTimeSeries())
    {
        outputPtr_ =
            ecgModelIO::openTimeSeries
            (
                outDir,
                solverType_ + ".dat",
                electrodeNames_
            );
    }

    verificationModelPtr_ = ecgVerificationModel::New
    (
        *this,
        domainName_,
        dict,
        electrodeNames_,
        electrodePositions_
    );
    if (verificationModelPtr_.valid())
    {
        verificationModelPtr_->bindSolver(*solverPtr_);
        verificationModelPtr_->validateProvider();
    }

    numericValues_.setSize(electrodeNames_.size(), 0.0);
}


ecgDomain::~ecgDomain() = default;


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void ecgDomain::evolve
(
    scalar t0,
    scalar dt
)
{
    solverPtr_->solve(*this, t0, dt, numericValues_);

    if (!solverPtr_->writesOwnTimeSeries() && verificationModelPtr_.valid())
    {
        verificationModelPtr_->record(numericValues_);
    }

    if
    (
        !solverPtr_->writesOwnTimeSeries()
     && (writeEveryTimeStep_ || mesh_.time().outputTime())
    )
    {
        ecgModelIO::writeRow
        (
            outputPtr_.ref(), mesh_.time().value(), numericValues_
        );
    }
}


void ecgDomain::recordVerification
(
    scalar sampleTime,
    const scalarField& values
)
{
    if (verificationModelPtr_.valid())
    {
        verificationModelPtr_->record(sampleTime, values);
    }
}


void ecgDomain::end()
{
    if (verificationModelPtr_.valid())
    {
        verificationModelPtr_->end();
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
