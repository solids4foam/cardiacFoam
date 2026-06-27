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

Application
    runPurkinjeGraph

Description
    Advance a graph-backed Purkinje conductionSystemDomain without myocardium
    coupling. This is a diagnostic utility for checking 1D graph evolution.

Author
    Simao Nieto de Castro. All rights reserved.
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

    Info<< nl
        << "Running Purkinje graph domain '" << selectedDomainName << "'"
        << " without myocardium coupling" << nl
        << "  deltaT:  " << runTime.deltaTValue() << nl
        << "  endTime: " << runTime.endTime().value() << nl
        << endl;

    autoPtr<conductionSystemDomain> conductionDomain
    (
        conductionSystemDomain::New
        (
            mesh,
            selectedDomainName,
            domainDict,
            runTime.deltaTValue()
        )
    );

    conductionDomain->write();

    while (runTime.run())
    {
        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        conductionDomain->advance
        (
            runTime.value() - runTime.deltaTValue(),
            runTime.deltaTValue()
        );

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
