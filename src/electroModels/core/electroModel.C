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
#include "dimVoltage.H"
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
    // Pre-read electroProperties without registering to the database.
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
            false  // do not register
        )
    );

    // Single entry point: 'myocardiumSolver' is the canonical selector key.
    // Direct electroModel subclasses (e.g. singleCellSolver) register directly
    // in this table.  Myocardium spatial solvers (monodomainSolver,
    // bidomainSolver, eikonalSolver) are all handled by electrophysiologyModel,
    // which registers itself under those three names.
    const word modelType(props.lookup("myocardiumSolver"));

    Info<< nl << "Selecting electroModel entry for myocardiumSolver "
        << modelType << endl;

    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (ctorPtr)
    {
        // Direct electroModel subclass or orchestration wrapper
        // (e.g. singleCellSolver, electrophysiologyModel).
        return autoPtr<electroModel>(ctorPtr(runTime, region));
    }

    // Type not found in either table -- give a clear error listing all valid
    // type names so the user knows exactly what to fix in their electroProperties.
    if (!dictionaryConstructorTablePtr_)
    {
        FatalErrorInFunction
            << "No electroModel implementations are registered in the "
            << "runtime table while selecting myocardiumSolver '"
            << modelType << "'."
            << exit(FatalError);
    }

    FatalErrorInLookup
    (
        "myocardiumSolver",
        modelType,
        *dictionaryConstructorTablePtr_
    ) << exit(FatalError);

    return autoPtr<electroModel>();
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

    if (domainSystem_.hasMyocardium())
    {
        domainSystem_.myocardium().write();
    }
}


void Foam::electroModel::configureECGDomains()
{
    electrophysicsSystemBuilder::configureECGDomains
    (
        domainSystem_,
        domainSystem_.myocardium(),
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

        electrophysicsSystemBuilder::configureConductionDomains
        (
            domainSystem_,
            mesh(),
            electroProperties_,
            runTime().deltaTValue()
        );

        configureECGDomains();

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

    domainSystem_.writeConductionCouplings();
    domainSystem_.writeConductionDomains();
    domainSystem_.writeECGDomains();

    return converged;
}


Foam::tmp<Foam::volScalarField>
Foam::electroModel::couplingField(const word& fieldName) const
{
    // Map field name to CouplingSignal enum
    CouplingSignal sig;
    if (fieldName == "Vm")
    {
        sig = CouplingSignal::VM;
    }
    else if (fieldName == "Cai")
    {
        sig = CouplingSignal::CAI;
    }
    else
    {
        FatalErrorInFunction
            << "Unknown coupling field \"" << fieldName << nl
            << "Valid options are: Vm, Cai."
            << abort(FatalError);
        // suppress compiler warning — unreachable
        sig = CouplingSignal::VM;
    }

    const ElectromechanicalSignalProvider* p = provider();

    if (!p)
    {
        FatalErrorInFunction
            << "No ElectromechanicalSignalProvider available for electroModel "
            << type() << nl
            << "The electro solver does not expose a coupling signal."
            << abort(FatalError);
    }

    if (!p->hasSignal(sig))
    {
        FatalErrorInFunction
            << "Coupling signal \"" << fieldName
            << "\" is not available from the ionic model." << nl
            << "Check that the ionic model has the corresponding state variable."
            << abort(FatalError);
    }

    const fvMesh& msh = refCast<const fvMesh>(mesh());

    tmp<volScalarField> tResult
    (
        new volScalarField
        (
            IOobject
            (
                fieldName,
                msh.time().timeName(),
                msh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            msh,
            dimensionedScalar(fieldName, dimless, scalar(0))
        )
    );

    scalarField& resultI = tResult.ref().primitiveFieldRef();
    forAll(resultI, cellI)
    {
        resultI[cellI] = p->signal(cellI, sig);
    }

    tResult.ref().correctBoundaryConditions();

    return tResult;
}


void Foam::electroModel::end()
{
    domainSystem_.endECGCouplings();
    domainSystem_.endECGDomains();
    domainSystem_.endConductionCouplings();
    domainSystem_.endConductionDomains();

    this->IOobject::rename(this->IOobject::name()+".withDefaultValues");
    this->regIOobject::write();
}
// ************************************************************************* //
