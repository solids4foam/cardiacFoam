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
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(ECGDomain, 0);
defineRunTimeSelectionTable(ECGDomain, dictionary);


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void ECGDomain::readElectrodes(const dictionary& dict)
{
    electrodeNames_.clear();
    electrodePositions_.clear();

    if (!dict.found("electrodes"))
    {
        return;
    }

    const dictionary& eDict = dict.subDict("electrodes");
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
    stateProvider_(stateProvider),
    mesh_(stateProvider.mesh()),
    electrodeNames_(),
    electrodePositions_()
{
    Info<< domainName << " state source: electrophysiology provider"
        << nl << endl;

    // Validate that the current diagnostic provider exposes the fields used by
    // the existing ECG implementations. Future torso-domain ECG models may
    // derive from a different base and use different state channels.
    (void)Vm();
    (void)conductivity();

    readElectrodes(dict);
}


// * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * * //

autoPtr<ECGDomain> ECGDomain::New
(
    const electroStateProvider& stateProvider,
    const dictionary& dict,
    const word& domainName
)
{
    const word modelType
    (
        dict.lookupOrDefault<word>("ECGDomain", word::null)
    );

    if (modelType.empty())
    {
        FatalErrorInFunction
            << "No ECGDomain specified in dictionary '"
            << dict.dictName() << "'." << nl
            << "Expected key ECGDomain."
            << exit(FatalError);
    }

    Info<< "Selecting ECGDomain " << modelType << nl;

    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalErrorInFunction
            << "Unknown ECGDomain type " << modelType << nl
            << "Valid types:" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<ECGDomain>
    (
        ctorPtr
        (
            stateProvider,
            dict.subDict(modelType + "Coeffs"),
            domainName
        )
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool ECGDomain::read(const dictionary& dict)
{
    readElectrodes(dict);

    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
