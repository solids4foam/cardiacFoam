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
#include "electrophysicsSystemBuilder.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    const dimensionSet dimVoltage(dimMass*dimArea/(pow3(dimTime)*dimCurrent));

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
        // If region == "region0" then read from the main case
        // Otherwise, read from the region/sub-mesh directory e.g.
        // constant/electro
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
            region, // using 'local' property of IOobject
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
    ),
    advanceTimings_(),
    setDeltaT_(true),
    domainSystem_()
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
    // NB: dictionary must be unregistered to avoid adding to the database

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
            false  // Do not register
        )
    );

    const word modelType(props.lookup("myocardiumSolver"));

    Info<< nl << "Selecting myocardiumSolver " << modelType << endl;

    // electroModel is now a single concrete entry point (not subclassed
    // per-solver). The myocardiumSolver selection happens inside
    // configureMyocardiumDomain() via the reactionDiffusionSolver factory.
    // singleCellSolver is the only remaining electroModel subclass.

#if (OPENFOAM >= 2112)
    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (ctorPtr)
    {
        // Legacy path: singleCellSolver or any surviving subclass
        return autoPtr<electroModel>(ctorPtr(runTime, region));
    }
#else
    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter != dictionaryConstructorTablePtr_->end())
    {
        return autoPtr<electroModel>(cstrIter()(runTime, region));
    }
#endif

    // Default path: plain electroModel constructed in-place, will call
    // configureMyocardiumDomain() via its own constructor.
    return autoPtr<electroModel>(new electroModel(typeName, runTime, region));
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


void Foam::electroModel::configureECGSystem()
{
    electrophysicsSystemBuilder::configureECGDomain
    (
        domainSystem_,
        *this,
        electroProperties_
    );
}


void Foam::electroModel::setDeltaT(Time& runTime)
{
    if (solutionAlg() == solutionAlgorithm::EXPLICIT && setDeltaT_)
    {
        setDeltaT_ = false;

        const scalar maxCo =
            runTime.controlDict().lookupOrDefault<scalar>("maxCo", 0.1);

        const scalar newDeltaT =
            domainSystem_.myocardium().suggestExplicitDeltaT(maxCo);

        Info << "Setting deltaT = " << newDeltaT << ", maxCo = " << maxCo
             << endl;

        runTime.setDeltaT(newDeltaT);
    }
    
    physicsModel::setDeltaT(runTime);
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

        electrophysicsSystemBuilder::configureAdvanceScheme
        (
            domainSystem_,
            electroProperties_
        );

        electrophysicsSystemBuilder::configureConductionSystemDomain
        (
            domainSystem_,
            mesh(),
            electroProperties_,
            runTime().deltaTValue()
        );

        configureECGSystem();

        return true;
    }
    else
    {
        return false;
    }
}


bool Foam::electroModel::evolve()
{
    const scalar dt = runTime().deltaTValue();
    const scalar t0 = runTime().value() - dt;

    pimpleControl* pimplePtr = nullptr;

    if (solutionAlg() == solutionAlgorithm::IMPLICIT)
    {
        pimplePtr = &pimple();
    }
    else if (solutionAlg() != solutionAlgorithm::EXPLICIT)
    {
        FatalErrorInFunction
            << "Unrecognised solution algorithm. Available options are "
            << solutionAlgorithmNames_[solutionAlgorithm::IMPLICIT]
            << ", "
            << solutionAlgorithmNames_[solutionAlgorithm::EXPLICIT]
            << endl;
    }

    const bool converged =
        domainSystem_.advance(t0, dt, pimplePtr, advanceTimings_);

    if (domainSystem_.hasMyocardium())
    {
        if (domainSystem_.myocardium().shouldPostProcess())
        {
            domainSystem_.myocardium().postProcess();
        }
        else if (runTime().outputTime())
        {
            domainSystem_.myocardium().exportStates();
        }
    }

    domainSystem_.writeUpstreamCouplingModel();
    domainSystem_.writeUpstreamDomain();
    domainSystem_.writeDownstreamDomain();

    return converged;
}


void Foam::electroModel::end()
{
    domainSystem_.endDownstreamCouplingModel();
    domainSystem_.endDownstreamDomain();
    domainSystem_.endUpstreamCouplingModel();
    domainSystem_.endUpstreamDomain();

    this->IOobject::rename(this->IOobject::name()+".withDefaultValues");
    this->regIOobject::write();
}
// ************************************************************************* //
