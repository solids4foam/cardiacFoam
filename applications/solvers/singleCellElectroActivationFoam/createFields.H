Info<< "Reading ionic model and tissue properties\n" << endl;
IOdictionary cardiacProperties
(
    IOobject
    (
        "cardiacProperties",
        runTime.constant(),
        runTime,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);
Info<< "Reading numerical solver properties\n" << endl;
IOdictionary numericalsolverProperties
(
    IOobject
    (
        "numericalsolverProperties",
        runTime.constant(),
        runTime,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

Info<< "Reading stimulus protocol\n" << endl;
IOdictionary stimulusProtocol
(
    IOobject
    (
        "stimulusProtocol",
        runTime.constant(),
        runTime,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

//Merge the dictionaries and create a memory object to read all properties.
dictionary solutionVariablesMemory(cardiacProperties);
solutionVariablesMemory.merge(numericalsolverProperties);
solutionVariablesMemory.merge(stimulusProtocol);






