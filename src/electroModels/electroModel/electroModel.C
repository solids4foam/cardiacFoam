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

#include "electroModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

    defineTypeNameAndDebug(electroModel, 0);
    defineRunTimeSelectionTable(electroModel, dictionary);
    addToRunTimeSelectionTable(physicsModel, electroModel, physicsModel);

    const Enum<electroModel::solutionAlgorithm>
    electroModel::solutionAlgorithmNames_
    ({
        {
            electroModel::solutionAlgorithm::IMPLICIT,
            "implicit"
        },
        {
            electroModel::solutionAlgorithm::EXPLICIT,
            "explicit"
        },
    });
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::electroModel::makePimpleControl() const
{
    if (pimplePtr_)
    {
        FatalErrorIn("void Foam::electroModel::makePimpleControl() const")
            << "pointer already set" << abort(FatalError);
    }

    pimplePtr_.reset
    (
        new pimpleControl
        (
            const_cast<fvMesh&>
            (
                refCast<const fvMesh>(mesh())
            )
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::electroModel::electroModel
(
    const word& type,
    Time& runTime,
    const word& region
)
:
    physicsModel(type, runTime),
    IOdictionary
    (
        bool(region == dynamicFvMesh::defaultRegion)
      ? IOobject
        (
            "electroProperties",
            runTime.caseConstant(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
      : IOobject
        (
            "electroProperties",
            runTime.caseConstant(),
            region,
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    meshPtr_
    (
        dynamicFvMesh::New
        (
            IOobject
            (
                region,
                runTime.timeName(),
                runTime,
                IOobject::MUST_READ
            )
        )
    ),
    electroProperties_(subDict(type + "Coeffs")),
    pimplePtr_(),
    solutionAlgorithm_
    (
        solutionAlgorithmNames_.get("solutionAlgorithm", electroProperties_)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::electroModel::~electroModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::electroModel> Foam::electroModel::New
(
    Time& runTime,
    const word& region
)
{
    IOdictionary props
    (
        IOobject
        (
            "electroProperties",
            bool(region == dynamicFvMesh::defaultRegion)
          ? fileName(runTime.caseConstant())
          : fileName(runTime.caseConstant()/region),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    const word modelType(props.lookup("electroModel"));

    Info<< nl << "Selecting electroModel " << modelType << endl;

#if (OPENFOAM >= 2112)
    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            props,
            "electroModel",
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
            "electroModel::New(Time&, const word&)"
        )   << "Unknown electroModel type " << modelType
            << endl << endl
            << "Valid electroModel types are :" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    auto* ctorPtr = cstrIter();
#endif

    return autoPtr<electroModel>(ctorPtr(runTime, region));
}


Foam::pimpleControl& Foam::electroModel::pimple()
{
    if (!pimplePtr_)
    {
        makePimpleControl();
    }

    return pimplePtr_();
}


void Foam::electroModel::writeFields(const Time& runTime)
{
    physicsModel::writeFields(runTime);
}


const Foam::dictionary& Foam::electroModel::ionicProperties() const
{
    return electroProperties_;
}


bool Foam::electroModel::read()
{
    if (regIOobject::read())
    {
        electroProperties_ = subDict(type() + "Coeffs");
        return true;
    }

    return false;
}


void Foam::electroModel::end()
{
    this->IOobject::rename(this->IOobject::name() + ".withDefaultValues");
    this->regIOobject::write();
}

// ************************************************************************* //
