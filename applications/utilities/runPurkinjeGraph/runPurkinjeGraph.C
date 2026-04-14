/*---------------------------------------------------------------------------*\
License
    This file is part of cardiacFoam.

Application
    runPurkinjeGraph

Description
    Advance a graph-backed Purkinje conductionSystemDomain without myocardium
    coupling. This is a diagnostic utility for checking 1D graph evolution.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "IOdictionary.H"
#include "conductionSystemDomain.H"

using namespace Foam;

namespace
{

const dictionary& selectElectroModelCoeffs
(
    const IOdictionary& electroProperties
)
{
    const word solverType
    (
        electroProperties.lookupOrDefault<word>
        (
            "myocardiumSolver",
            "monodomainSolver"
        )
    );

    const word coeffsName(solverType + "Coeffs");

    if (electroProperties.found(coeffsName))
    {
        return electroProperties.subDict(coeffsName);
    }

    return electroProperties;
}


const dictionary& selectConductionDomainDict
(
    const dictionary& electroCoeffs,
    const word& requestedName,
    word& selectedName
)
{
    const dictionary& domains = electroCoeffs.subDict("conductionNetworkDomains");

    if (requestedName.size())
    {
        selectedName = requestedName;
        return domains.subDict(selectedName);
    }

    forAllConstIter(dictionary, domains, iter)
    {
        const entry& e = iter();

        if (e.isDict())
        {
            selectedName = e.keyword();
            return e.dict();
        }
    }

    FatalErrorInFunction
        << "No dictionary entries were found in conductionNetworkDomains."
        << exit(FatalError);

    return domains;
}

} // End anonymous namespace


int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Advance a graph-backed Purkinje conductionSystemDomain without "
        "constructing myocardium coupling."
    );

    argList::noParallel();
    argList::addOption
    (
        "conductionDomain",
        "word",
        "Name of the conductionNetworkDomains entry to run"
    );
    argList::addOption
    (
        "nSteps",
        "label",
        "Number of graph-only time steps to run"
    );
    argList::addOption
    (
        "deltaT",
        "scalar",
        "Override time-step size for this graph-only diagnostic"
    );

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    IOdictionary electroProperties
    (
        IOobject
        (
            "electroProperties",
            runTime.caseConstant(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    const dictionary& electroCoeffs =
        selectElectroModelCoeffs(electroProperties);

    word selectedDomainName;
    const dictionary& domainDict =
        selectConductionDomainDict
        (
            electroCoeffs,
            args.getOrDefault<word>("conductionDomain", word::null),
            selectedDomainName
        );

    const scalar deltaT =
        args.getOrDefault<scalar>("deltaT", runTime.deltaTValue());

    const label nSteps = args.getOrDefault<label>("nSteps", 10);

    if (nSteps < 0)
    {
        FatalErrorInFunction
            << "nSteps must be non-negative. Received " << nSteps << "."
            << exit(FatalError);
    }

    runTime.setDeltaT(deltaT);

    Info<< nl
        << "Running Purkinje graph domain '" << selectedDomainName << "'"
        << " without myocardium coupling" << nl
        << "  nSteps: " << nSteps << nl
        << "  deltaT: " << runTime.deltaTValue() << nl
        << endl;

    autoPtr<ConductionSystemDomain> conductionDomain
    (
        ConductionSystemDomain::New
        (
            mesh,
            domainDict,
            runTime.deltaTValue()
        )
    );

    conductionDomain->write();

    for (label stepI = 0; stepI < nSteps; ++stepI)
    {
        const scalar t0 = runTime.value();
        const scalar dt = runTime.deltaTValue();

        conductionDomain->advance(t0, dt);

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        if (runTime.outputTime())
        {
            conductionDomain->write();
        }

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    conductionDomain->end();

    Info<< nl << "End" << nl << endl;

    return 0;
}


// ************************************************************************* //
