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

#include "electroMechanicalModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(electroMechanicalModel, 0);
    defineRunTimeSelectionTable(electroMechanicalModel, dictionary);
    addToRunTimeSelectionTable
    (
        physicsModel, electroMechanicalModel, physicsModel
    );

    namespace
    {
        word resolvedRegion
        (
            const Time& runTime,
            const dictionary& electroMechanicalProperties,
            const word& modelName,
            const word& defaultRegion
        )
        {
            const word controlDictRegion
            (
                runTime.controlDict().subOrEmptyDict(modelName)
                    .lookupOrDefault<word>("region", defaultRegion)
            );

            return electroMechanicalProperties.lookupOrDefault<word>
            (
                modelName + "Region",
                controlDictRegion
            );
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::electroMechanicalModel::electroMechanicalModel
(
    const word& type,
    Time& runTime,
    const word& region
)
:
    physicsModel(type, runTime, region),
    IOdictionary
    (
        IOobject
        (
            "electroMechanicalProperties",
            runTime.constant(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    electroMechanicalProperties_(subDict(type + "Coeffs")),
    electro_(),
    solid_()
{
    const word electroRegion
    (
        resolvedRegion(runTime, electroMechanicalProperties_, "electro", "electro")
    );

    const word solidRegion
    (
        resolvedRegion(runTime, electroMechanicalProperties_, "solid", "solid")
    );

    Info<< "Electro-mechanical model: electro region = " << electroRegion
        << ", solid region = " << solidRegion << endl;

    electro_ = electroModel::New(runTime, electroRegion);
    solid_ = solidModel::New(runTime, solidRegion);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::electroMechanicalModel::~electroMechanicalModel()
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::electroMechanicalModel>
Foam::electroMechanicalModel::New
(
    Time& runTime,
    const word& region
)
{
    // NB: dictionary must be unregistered to avoid adding to the database

    IOdictionary props
    (
        IOobject
        (
            "electroMechanicalProperties",
            runTime.constant(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false  // Do not register
        )
    );

    const word modelType(props.lookup("electroMechanicalModel"));

    Info<< "Selecting electroMechanicalModel " << modelType << endl;

#if (OPENFOAM >= 2112)
    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            props,
            "electroMechanicalModel",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

#else
    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "electroMechanicalModel::New(Time&, const word&)"
        )   << "Unknown electroMechanicalModel type " << modelType
            << endl << endl
            << "Valid electroMechanicalModel types are :" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    auto* ctorPtr = cstrIter();
#endif

    return autoPtr<electroMechanicalModel>(ctorPtr(runTime, region));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::electroMechanicalModel::setDeltaT(Time& runTime)
{
    // The electro model typically requires the smallest time-step
    electro().setDeltaT(runTime);
}


void Foam::electroMechanicalModel::writeFields(const Time& runTime)
{
    // The solid calls runTime.write() which writes both region fields
    solid().writeFields(runTime);
}


void Foam::electroMechanicalModel::end()
{
    this->IOobject::rename(this->IOobject::name() + ".withDefaultValues");
    this->regIOobject::write();
    electro().end();
    solid().end();
}


// ************************************************************************* //
