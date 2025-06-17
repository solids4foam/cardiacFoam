// Transmembrane potential
volScalarField Vm
(
    IOobject
    (
        "Vm",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("Vm", dimMass*dimArea/(pow3(dimTime)*dimCurrent), -85.23),
    "zeroGradient"
);

// Create commonly used dimensionSets for convenience
// Read electoActivationProperties
Info<< "Reading electroActivationProperties\n" << endl;
IOdictionary electroActivationProperties
(
    IOobject
    (
        "electroActivationProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

Info<< "Reading conductivity\n" << endl;
const dimensionedSymmTensor conductivity
(
    "conductivity",
    pow3(dimTime)*sqr(dimCurrent)/(dimMass*dimVolume),
    electroActivationProperties
);

Info<< "Reading cellular surface-to-volume ratio beta\n" << endl;
const dimensionedScalar beta
(
    "beta", dimArea/dimVolume, electroActivationProperties
);

Info<< "Reading membrane capacitance Cm\n" << endl;
const dimensionSet dimVoltage(dimMass*dimArea/(pow3(dimTime)*dimCurrent));
const dimensionedScalar Cm
(
    "Cm",
    dimCurrent*dimTime/(dimVoltage*dimArea),
    electroActivationProperties
);

// Info<< "Reading diffusion coefficient D\n" << endl;
// const dimensionedScalar D
// (
//     "D", dimArea/dimTime, electroActivationProperties
// );

// Find cells in the stimulus volume
Info<< "Reading stimulus location \n" << endl;
boundBox bb
(
    point(electroActivationProperties.lookup("stimulusLocationMin")),
    point(electroActivationProperties.lookup("stimulusLocationMax"))
);

labelHashSet stimCellSet;
forAll(mesh.C(), cellI)
{
    if (bb.contains(mesh.C()[cellI]))
    {
        stimCellSet.insert(cellI);
    }
}
const labelList stimulusCellIDs = stimCellSet.toc();
Info<< "Number of cells in the stimulus volume: "
    << stimulusCellIDs.size() << endl;

Info<< "Reading stimulus duration\n" << endl;
const dimensionedScalar stimulusDuration
(
    "stimulusDuration", dimTime, electroActivationProperties
);

Info<< "Reading stimulus intensity\n" << endl;
// Per unit area as we multiply by the cellular surface-to-volume ratio
const dimensionedScalar stimulusIntensity
(
    "stimulusIntensity",
    // dimCurrent/dimArea,
    dimCurrent/dimVolume,
    electroActivationProperties
);

// Ionic source term
volScalarField ionicCurrent
(
    IOobject
    (
        "ionicCurrent",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless/dimTime, 0.0),
    "zeroGradient"
);

// // Scale factor for Bueno-Orovio
// const dimensionedScalar BuenoOrovioScaleFactor
// (
//     "BuenoOrovioScaleFactor", dimVoltage, 84.0
// );

// External stimulus current
volScalarField externalStimulusCurrent
(
    IOobject
    (
        "externalStimulusCurrent",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    // dimensionedScalar("sourceTerm", dimCurrent/dimArea, 0.0),
    dimensionedScalar("sourceTerm", dimCurrent/dimVolume, 0.0),
    "zeroGradient"
);

// Read ionic model dict
const dictionary& ionicModelCoeffs =
    electroActivationProperties.subDict("ionicModelCoeffs");

// Create ionicModelCellML object
ionicModelCellML ionicModel(ionicModelCoeffs, mesh.nCells());

// Field to record the activation time, which is the time at which each cell
// reaches or exceeds a transmembrane voltage of 0.0
volScalarField activationTime
(
    IOobject
    (
        "activationTime",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimTime, 0.0),
    "zeroGradient"
);
boolList calculateActivationTime(mesh.nCells(), true);

// Activation velocity
volVectorField activationVelocity
(
    IOobject
    (
        "activationVelocity",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    fvc::grad(1.0/(activationTime + dimensionedScalar("SMALL", dimTime, SMALL)))
);

// Scale factor for Bueno-Orovio
const dimensionedScalar BuenoOrovioScaleFactor
(
    "BuenoOrovioScaleFactor", dimVoltage, 1 //84.0
);
