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

#include "ecgModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(ecgModel, 0);
defineRunTimeSelectionTable(ecgModel, dictionary);


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void ecgModel::readElectrodes(const dictionary& dict)
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
        electrodeNames_[i] = names[i];
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


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

ecgModel::ecgModel
(
    const volScalarField& Vm,
    const dictionary& dict,
    const volTensorField& conductivity
)
:
    Vm_(Vm),
    mesh_(Vm.mesh()),
    conductivity_(conductivity),
    electrodeNames_(),
    electrodePositions_()
{
    Info<< "ECG conductivity source: monodomain conductivity field"
        << nl << endl;

    readElectrodes(dict);
}


// * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * * //

autoPtr<ecgModel> ecgModel::New
(
    const volScalarField& Vm,
    const dictionary& dict,
    const volTensorField& conductivity
)
{
    const word modelType(dict.lookup("ecgModel"));

    Info<< "Selecting ecgModel " << modelType << nl;

    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalErrorInFunction
            << "Unknown ecgModel type " << modelType << nl
            << "Valid types:" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<ecgModel>
    (
        ctorPtr(Vm, dict.subDict(modelType + "Coeffs"), conductivity)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool ecgModel::read(const dictionary& dict)
{
    readElectrodes(dict);

    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
