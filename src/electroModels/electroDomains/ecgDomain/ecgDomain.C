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
#include "ecgModelIO.H"
#include "ecgVerificationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(ecgDomain, 0);


// * * * * * * * * * * * * Private Helpers * * * * * * * * * * * * * * * * * //

namespace
{

word selectedECGSolverType(const dictionary& dict)
{
    return dict.lookupOrDefault<word>("ecgSolver", "pseudoECG");
}


word outputFileName(const word& solverType)
{
    return
        solverType == "torsoECG"
      ? word("torsoECG.dat")
      : word("pseudoECG.dat");
}


dictionary withInheritedManufacturedBidomain
(
    const dictionary& dict,
    const dictionary* inheritedManufacturedBidomainPtr,
    const dictionary* inheritedPotentialDomainPtr
)
{
    dictionary merged(dict);

    if
    (
        inheritedManufacturedBidomainPtr
     && !merged.found("manufacturedBidomain")
    )
    {
        merged.add("manufacturedBidomain", *inheritedManufacturedBidomainPtr);
    }

    if (inheritedPotentialDomainPtr)
    {
        if
        (
            !merged.found("groundPatches")
         && inheritedPotentialDomainPtr->found("groundPatches")
        )
        {
            merged.add
            (
                "groundPatches",
                inheritedPotentialDomainPtr->subDict("groundPatches")
            );
        }

        if
        (
            !merged.found("surfaceCurrentPatches")
         && inheritedPotentialDomainPtr->found("surfaceCurrentPatches")
        )
        {
            merged.add
            (
                "surfaceCurrentPatches",
                inheritedPotentialDomainPtr->subDict("surfaceCurrentPatches")
            );
        }
    }

    return merged;
}


void finalizeVerificationModel(autoPtr<ecgVerificationModel>& verifierPtr)
{
    if (verifierPtr.valid())
    {
        verifierPtr->end();
        verifierPtr.clear();
    }
}

}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void ecgDomain::readElectrodes(const dictionary& dict)
{
    electrodeNames_.clear();
    electrodePositions_.clear();

    const dictionary* eDictPtr = dict.findDict("electrodePositions");

    if (!eDictPtr)
    {
        eDictPtr = inheritedElectrodePositionsPtr_;
    }

    if (!eDictPtr)
    {
        return;
    }

    const dictionary& eDict = *eDictPtr;
    const wordList names(eDict.toc());

    electrodeNames_.setSize(names.size());
    electrodePositions_.setSize(names.size());

    forAll(names, i)
    {
        electrodeNames_[i]     = names[i];
        electrodePositions_[i] = eDict.get<vector>(names[i]);
    }

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
    const dictionary& dict,
    const word& domainName,
    const dictionary* inheritedElectrodePositionsPtr,
    const dictionary* inheritedManufacturedBidomainPtr,
    const dictionary* inheritedPotentialDomainPtr
)
:
    stateProvider_(stateProvider),
    mesh_(stateProvider.mesh()),
    solverType_(selectedECGSolverType(dict)),
    outputPtr_(),
    solverPtr_(ecgSolver::New(dict)),
    verificationModelPtr_(),
    numericValues_(),
    inheritedElectrodePositionsPtr_(inheritedElectrodePositionsPtr),
    inheritedManufacturedBidomainPtr_(inheritedManufacturedBidomainPtr),
    inheritedPotentialDomainPtr_(inheritedPotentialDomainPtr),
    electrodeNames_(),
    electrodePositions_()
{
    Info<< domainName << " ECG solver: " << solverType_
        << " on mesh '" << mesh_.name() << "'" << nl
        << "  state source: electrophysiology provider"
        << nl << endl;

    if (solverType_ == "pseudoECG")
    {
        (void)Vm();
        (void)conductivity();
    }
    else if (solverType_ == "torsoECG")
    {
        (void)phiE();
    }

    readElectrodes(dict);

    const fileName outDir(mesh_.time().globalPath() / "postProcessing");
    outputPtr_ =
        ecgModelIO::openTimeSeries
        (
            outDir,
            outputFileName(solverType_),
            electrodeNames_
        );

    const dictionary verificationDict =
        withInheritedManufacturedBidomain
        (
            dict,
            inheritedManufacturedBidomainPtr_,
            inheritedPotentialDomainPtr_
        );

    verificationModelPtr_ = ecgVerificationModel::New
    (
        *this,
        verificationDict,
        electrodeNames_,
        electrodePositions_
    );

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

    if (verificationModelPtr_.valid())
    {
        verificationModelPtr_->record(numericValues_);
    }

    if (mesh_.time().outputTime())
    {
        ecgModelIO::writeRow
        (
            outputPtr_.ref(), mesh_.time().value(), numericValues_
        );
    }
}


bool ecgDomain::read(const dictionary& dict)
{
    const wordList previousElectrodeNames(electrodeNames_);

    readElectrodes(dict);

    numericValues_.setSize(electrodeNames_.size(), 0.0);

    if (electrodeNames_.size() != previousElectrodeNames.size())
    {
        FatalErrorInFunction
            << "Changing the number of ECG electrodes during read() is not "
               "supported because the output column layout is fixed at startup."
            << exit(FatalError);
    }

    forAll(previousElectrodeNames, electrodeI)
    {
        if (electrodeNames_[electrodeI] != previousElectrodeNames[electrodeI])
        {
            FatalErrorInFunction
                << "Changing ECG electrode names during read() is not "
                   "supported because the output column layout is fixed at "
                   "startup."
                << exit(FatalError);
        }
    }

    const dictionary verificationDict =
        withInheritedManufacturedBidomain
        (
            dict,
            inheritedManufacturedBidomainPtr_,
            inheritedPotentialDomainPtr_
        );

    const word requestedType
    (
        ecgVerificationModel::selectedType(verificationDict)
    );

    if (requestedType.empty())
    {
        finalizeVerificationModel(verificationModelPtr_);
        return true;
    }

    if
    (
        verificationModelPtr_.valid()
     && verificationModelPtr_->type() == requestedType
    )
    {
        verificationModelPtr_->updateElectrodes
        (
            electrodeNames_,
            electrodePositions_
        );

        return verificationModelPtr_->read(verificationDict);
    }

    finalizeVerificationModel(verificationModelPtr_);
    verificationModelPtr_ = ecgVerificationModel::New
    (
        *this,
        verificationDict,
        electrodeNames_,
        electrodePositions_
    );

    return true;
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
