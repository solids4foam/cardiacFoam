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


// * * * * * * * * * * * * * Private Static Members * * * * * * * * * * * * //

tmp<volTensorField> ecgModel::initialiseGi
(
    const dictionary& dict,
    const fvMesh& mesh
)
{
    tmp<volTensorField> tresult
    (
        new volTensorField
        (
            IOobject
            (
                "Gi",
                mesh.time().timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedTensor
            (
                "zero",
                pow3(dimTime)*sqr(dimCurrent)/(dimMass*dimVolume),
                tensor::zero
            )
        )
    );
    volTensorField& result = tresult.ref();

    if (!result.headerOk())
    {
        Info<< "\nGi not found on disk, using Gi from dict" << nl << endl;

        result =
            dimensionedTensor
            (
                dimensionedSymmTensor
                (
                    "Gi",
                    pow3(dimTime)*sqr(dimCurrent)/(dimMass*dimVolume),
                    dict
                ) & tensor(I)
            );

        if (result.size() > 0)
        {
            Info<< "Gi tensor (cell 0): " << result[0] << nl;
        }
    }
    else
    {
        Info<< "Gi field read from " << mesh.time().timeName() << nl << endl;
    }

    return tresult;
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void ecgModel::readElectrodes(const dictionary& dict)
{
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


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

ecgModel::ecgModel
(
    const volScalarField& Vm,
    const dictionary& dict,
    const volTensorField* conductivityPtr
)
:
    Vm_(Vm),
    mesh_(Vm.mesh()),
    GiOwnedPtr_(),
    GiPtr_(nullptr),
    sigmaT_
    (
        "sigmaT",
        pow3(dimTime)*sqr(dimCurrent)/(dimMass*dimVolume),
        dict
    ),
    electrodeNames_(),
    electrodePositions_()
{
    if (conductivityPtr)
    {
        GiPtr_ = conductivityPtr;
        Info<< "ECG conductivity source: monodomain conductivity field"
            << nl << endl;
    }
    else
    {
        GiOwnedPtr_.reset(new volTensorField(initialiseGi(dict, Vm.mesh())));
        GiPtr_ = GiOwnedPtr_.operator->();
    }

    readElectrodes(dict);
}


// * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * * //

autoPtr<ecgModel> ecgModel::New
(
    const volScalarField& Vm,
    const dictionary& dict,
    const volTensorField* conductivityPtr
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
        ctorPtr(Vm, dict.subDict(modelType + "Coeffs"), conductivityPtr)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool ecgModel::read(const dictionary& dict)
{
    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
