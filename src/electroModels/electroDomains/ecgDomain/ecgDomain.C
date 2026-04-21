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

defineTypeNameAndDebug(ECGDomain, 0);


// * * * * * * * * * * * * Private Helpers * * * * * * * * * * * * * * * * * //

namespace
{

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

void ECGDomain::readElectrodes(const dictionary& dict)
{
    electrodeNames_.clear();
    electrodePositions_.clear();

    if (!dict.found("electrodePositions"))
    {
        return;
    }

    const dictionary& eDict = dict.subDict("electrodePositions");
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


const volScalarField& ECGDomain::Vm() const
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


const volTensorField& ECGDomain::conductivity() const
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


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

ECGDomain::ECGDomain
(
    const electroStateProvider& stateProvider,
    const dictionary& dict,
    const word& domainName
)
:
    outputPtr_(),
    solverPtr_(ECGSolver::New(dict)),
    verificationModelPtr_(),
    numericValues_(),
    stateProvider_(stateProvider),
    mesh_(stateProvider.mesh()),
    electrodeNames_(),
    electrodePositions_()
{
    Info<< domainName << " state source: electrophysiology provider"
        << nl << endl;

    // Validate state access
    (void)Vm();
    (void)conductivity();

    readElectrodes(dict);

    // Open output file
    const fileName outDir(mesh_.time().path() / "postProcessing");
    outputPtr_ =
        ecgModelIO::openTimeSeries(outDir, "pseudoECG.dat", electrodeNames_);

    // Optional verification model
    verificationModelPtr_ = ecgVerificationModel::New
    (
        stateProvider,
        dict,
        electrodeNames_,
        electrodePositions_
    );

    numericValues_.setSize(electrodeNames_.size(), 0.0);
}


ECGDomain::~ECGDomain() = default;


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void ECGDomain::evolve
(
    scalar t0,
    scalar dt
)
{
    (void)t0;
    (void)dt;

    solverPtr_->compute(*this, numericValues_);

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


bool ECGDomain::read(const dictionary& dict)
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

    const word requestedType(ecgVerificationModel::selectedType(dict));

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

        return verificationModelPtr_->read(dict);
    }

    finalizeVerificationModel(verificationModelPtr_);
    verificationModelPtr_ = ecgVerificationModel::New
    (
        stateProvider_,
        dict,
        electrodeNames_,
        electrodePositions_
    );

    return true;
}


void ECGDomain::end()
{
    if (verificationModelPtr_.valid())
    {
        verificationModelPtr_->end();
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
