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
#include "Pstream.H"
#include <cmath>

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


dictionary withInheritedPotentialDomain
(
    const dictionary& dict,
    const dictionary* inheritedPotentialDomainPtr
)
{
    dictionary merged(dict);

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
    const dictionary& dict,
    const word& domainName,
    const dictionary* inheritedElectrodePositionsPtr,
    const dictionary* inheritedPotentialDomainPtr
)
:
    stateProvider_(stateProvider),
    mesh_(stateProvider.mesh()),
    domainName_(domainName),
    solverType_(selectedECGSolverType(dict)),
    channelNames_(),
    channelPositions_(),
    sampledLeadField_(false),
    sampleStart_(0), sampleEnd_(0), sampleDelta_(0), nextSample_(0),
    outputPtr_(),
    solverPtr_(ecgSolver::New(dict)),
    verificationModelPtr_(),
    numericValues_(),
    inheritedElectrodePositionsPtr_(inheritedElectrodePositionsPtr),
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
    else if (solverType_ == "eikonalECG")
    {
        (void)activationTime();
        (void)conductivity();
    }
    else if (solverType_ == "torsoECG")
    {
        (void)phiE();
    }

    readElectrodes(dict);
    solverPtr_->initialise(*this);
    channelNames_ = solverPtr_->channelNames(*this);
    channelPositions_ = solverPtr_->channelPositions(*this);

    if (solverType_ == "leadFieldECG" && dict.found("sampling"))
    {
        const dictionary& sampling = dict.subDict("sampling");
        sampledLeadField_ = true;
        sampleStart_ = sampling.get<scalar>("start");
        sampleEnd_ = sampling.get<scalar>("end");
        sampleDelta_ = sampling.get<scalar>("deltaT");
        if (!std::isfinite(sampleStart_) || !std::isfinite(sampleEnd_)
         || !std::isfinite(sampleDelta_) || sampleDelta_ <= 0 || sampleEnd_ < sampleStart_)
        {
            FatalIOErrorInFunction(sampling) << "Invalid lead-field sampling interval."
                << exit(FatalIOError);
        }
        nextSample_ = sampleStart_ + max(scalar(0), std::ceil
            ((mesh_.time().value()-sampleStart_)/sampleDelta_))*sampleDelta_;
    }
    const fileName outDir(mesh_.time().globalPath() / "postProcessing");
    if (!solverPtr_->writesOwnTimeSeries() && Pstream::master())
    {
        outputPtr_ = solverType_ == "leadFieldECG"
          ? ecgModelIO::openRestartTimeSeries(outDir, word(solverType_+".dat"),
                channelNames_, mesh_.time().value())
          : ecgModelIO::openTimeSeries(outDir, word(solverType_+".dat"), channelNames_);
    }

    const dictionary verificationDict =
        withInheritedPotentialDomain
        (
            dict,
            inheritedPotentialDomainPtr_
        );

    verificationModelPtr_ = ecgVerificationModel::New
    (
        *this,
        domainName_,
        verificationDict,
        channelNames_,
        channelPositions_
    );

    numericValues_.setSize(channelNames_.size(), 0.0);
}



autoPtr<ecgDomain> ecgDomain::New
(
    const electroStateProvider& myocardium,
    const electroStateProvider* potential,
    const word& myocardiumType,
    const dictionary& dict,
    const word& name,
    const dictionary* inheritedElectrodes,
    const dictionary* inheritedPotential
)
{
    const word type(selectedECGSolverType(dict));
    const electroStateProvider* provider = &myocardium;
    if (type == "leadFieldECG")
    {
        if (myocardiumType != "monodomainSolver" && myocardiumType != "bidomainSolver")
        {
            FatalIOErrorInFunction(dict)
                << "leadFieldECG currently requires monodomainSolver or bidomainSolver."
                << exit(FatalIOError);
        }
    }
    else if (type == "eikonalECG")
    {
        if (myocardiumType != "eikonalSolver")
        {
            FatalIOErrorInFunction(dict) << "eikonalECG requires eikonalSolver."
                << exit(FatalIOError);
        }
    }
    else if (type == "torsoECG")
    {
        if (myocardiumType != "bidomainSolver" || !potential)
        {
            FatalIOErrorInFunction(dict)
                << "torsoECG requires bidomainSolver and bathPotentialDomain."
                << exit(FatalIOError);
        }
        provider = potential;
    }
    else if (type != "pseudoECG")
    {
        FatalIOErrorInFunction(dict) << "Unknown ecgSolver '" << type << "'."
            << exit(FatalIOError);
    }
    return autoPtr<ecgDomain>
    (
        new ecgDomain(*provider, dict, name, inheritedElectrodes, inheritedPotential)
    );
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

    bool emit = mesh_.time().outputTime();
    if (sampledLeadField_)
    {
        if (sampleDelta_+SMALL < dt)
        {
            FatalErrorInFunction << "leadFieldECG sampling.deltaT must be at least the solver timestep."
                << exit(FatalError);
        }
        const scalar time=mesh_.time().value();
        emit = time+SMALL >= sampleStart_ && time <= sampleEnd_+SMALL
            && time+SMALL >= nextSample_;
        if (emit) while (nextSample_ <= time+SMALL) nextSample_ += sampleDelta_;
    }
    if (!solverPtr_->writesOwnTimeSeries() && emit && Pstream::master())
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

    const List<vector> previousPositions(electrodePositions_);
    readElectrodes(dict);
    if (solverType_ == "leadFieldECG" && previousPositions != electrodePositions_)
    {
        FatalErrorInFunction << "Static leadFieldECG electrodes cannot change "
            << "during read(). Regenerate the fields and restart."
            << exit(FatalError);
    }
    channelNames_ = solverPtr_->channelNames(*this);
    channelPositions_ = solverPtr_->channelPositions(*this);

    numericValues_.setSize(channelNames_.size(), 0.0);

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
        withInheritedPotentialDomain
        (
            dict,
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
            channelNames_,
            channelPositions_
        );

        return verificationModelPtr_->read(verificationDict);
    }

    finalizeVerificationModel(verificationModelPtr_);
    verificationModelPtr_ = ecgVerificationModel::New
    (
        *this,
        domainName_,
        verificationDict,
        channelNames_,
        channelPositions_
    );

    return true;
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
