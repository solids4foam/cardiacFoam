// Define voltage dimensions for convenience
const dimensionSet dimVoltage(dimMass*dimArea/(pow3(dimTime)*dimCurrent));

Foam::solverStartupLog::readingDictionary("stimulus protocol");
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

const scalar stimStart = readScalar(stimulusProtocol.lookup("stim_start"));
const scalar stimDuration = readScalar(stimulusProtocol.lookup("stim_duration"));
const scalar stimAmplitude = readScalar(stimulusProtocol.lookup("stim_amplitude"));
const scalar stimPeriodS1 = readScalar(stimulusProtocol.lookup("stim_period_S1"));
const label nStim1 = readLabel(stimulusProtocol.lookup("nstim1"));
const scalar stimPeriodS2 = readScalar(stimulusProtocol.lookup("stim_period_S2"));
const label nStim2 = readLabel(stimulusProtocol.lookup("nstim2"));

if (stimDuration <= 0)
{
    FatalErrorInFunction
        << "stim_duration must be > 0. Current value: " << stimDuration
        << abort(FatalError);
}

if (nStim1 < 0 || nStim2 < 0)
{
    FatalErrorInFunction
        << "nstim1 and nstim2 must be >= 0. Current values: "
        << "nstim1=" << nStim1 << ", nstim2=" << nStim2
        << abort(FatalError);
}

Info<< "Stimulus protocol:" << nl
    << "  stim_start      = " << stimStart << "[ms]"<< nl
    << "  stim_duration   = " << stimDuration << "[ms]"<< nl
    << "  stim_amplitude  = " << stimAmplitude << "[ms]"<< nl
    << "  stim_period_S1  = " << stimPeriodS1 << "[ms]"<< nl
    << "  nstim1          = " << nStim1  <<  nl
    << "  stim_period_S2  = " << stimPeriodS2 << "[ms]"<< nl
    << "  nstim2          = " << nStim2 << nl
    << "  protocol suffix = " << stimulusIO::protocolSuffix(stimulusProtocol)
    << nl << endl;

Foam::solverStartupLog::sectionDivider();

Foam::solverStartupLog::readingDictionary("numerical solver properties");
IOdictionary timeIntegrationProperties
(
    IOobject
    (
        "timeIntegrationProperties",
        runTime.constant(),
        runTime,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

const Switch solveExplicit(timeIntegrationProperties.lookup("solveExplicit"));
word odeSolverName;
timeIntegrationProperties.lookup("solver") >> odeSolverName;
const scalar initialODEStep =
    readScalar(timeIntegrationProperties.lookup("initialODEStep"));
const label maxSteps = readLabel(timeIntegrationProperties.lookup("maxSteps"));
const scalar absTol =
    timeIntegrationProperties.lookupOrDefault<scalar>("absTol", 1e-15);
const scalar relTol =
    timeIntegrationProperties.lookupOrDefault<scalar>("relTol", 1e-4);

Info<< "Numerical settings:" << nl
    << "  solveExplicit   = " << solveExplicit << nl
    << "  solver          = " << odeSolverName << nl
    << "  initialODEStep  = " << initialODEStep << nl
    << "  maxSteps        = " << maxSteps << nl
    << "  absTol          = " << absTol << nl
    << "  relTol          = " << relTol << nl
    << endl;

Foam::solverStartupLog::sectionDivider();

Foam::solverStartupLog::readingDictionary("tissue properties");
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

word ionicModelName;
cardiacProperties.lookup("ionicModel") >> ionicModelName;
word tissueName;
cardiacProperties.lookup("tissue") >> tissueName;

Info<< "Cell/tissue settings:" << nl
    << "  ionicModel      = " << ionicModelName << nl
    << "  tissue          = " << tissueName << nl;

if (cardiacProperties.found("chi"))
{
    const dimensionedScalar chi("chi", dimArea/dimVolume, cardiacProperties);
    Info<< "  chi             = " << chi.value() << " " << chi.dimensions()
        << nl;
}

if (cardiacProperties.found("Cm"))
{
    const dimensionedScalar Cm
    (
        "Cm",
        dimCurrent*dimTime/(dimVoltage*dimArea),
        cardiacProperties
    );
    Info<< "  Cm              = " << Cm.value() << " " << Cm.dimensions()
        << nl;
}

Info<< endl;
Foam::solverStartupLog::sectionDivider();

// Merge the dictionaries and create a memory object to read all properties.
dictionary solutionVariablesMemory(cardiacProperties);
solutionVariablesMemory.merge(timeIntegrationProperties);
solutionVariablesMemory.merge(stimulusProtocol);
