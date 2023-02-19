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


// Normalised transmembrane potential
// Info<< "Reading field u\n" << endl;
// volScalarField u
// (
//     IOobject
//     (
//         "u",
//         runTime.timeName(),
//         mesh,
//         IOobject::MUST_READ,
//         IOobject::AUTO_WRITE
//     ),
//     mesh
// );

// Create commonly used dimensionSets for convenience
const dimensionSet dimVoltage(dimMass*dimArea/(pow3(dimTime)*dimCurrent));

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
const dimensionedScalar Cm
(
    "Cm",
    // pow4(dimTime)*sqr(dimCurrent)/(dimMass*sqr(dimArea)),
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

// Find the closest cell to the stimulus location
// label stimulusLocationCellID = -1;
// {
//     scalar minDist = GREAT;
//     const vectorField& CI = mesh.C();
//     forAll(CI, cellI)
//     {
//         const scalar dist = mag(CI[cellI] - stimulusLocation);

//         if (dist < minDist)
//         {
//             minDist = dist;
//             stimulusLocationCellID = cellI;
//         }
//     }

//     // Ensure only one processor has the closest point
//     const scalar globalMinDist = returnReduce(minDist, minOp<scalar>());
//     if (mag(globalMinDist - minDist) > SMALL)
//     {
//         // This proc does not have the global minDist
//         stimulusLocationCellID = -1;
//     }

//     // Check that at least one processor found a cell
//     if (returnReduce(stimulusLocationCellID, maxOp<label>()) == -1)
//     {
//         FatalError
//             << "Problem finding the stimulus location cellID"
//             << abort(FatalError);
//     }
// }

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
    // dimensionedScalar("zero", dimCurrent/dimVolume, 0.0),
    // dimensionedScalar("zero", dimCurrent/dimArea, 0.0),
    // dimensionedScalar
    // (
    //     "zero", dimMass*dimArea/(sqr(sqr(dimTime))*dimCurrent), 0.0
    // ),
    dimensionedScalar("zero", dimless/dimTime, 0.0),
    "zeroGradient"
);

// Scale factor for Bueno-Orovio
const dimensionedScalar BuenoOrovioScaleFactor
(
    "BuenoOrovioScaleFactor", dimVoltage, 84.0
);

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

// Create ionic model ODE system
BuenoOrovioIonicModel ionicModel(ionicModelCoeffs, mesh.nCells());

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

// Active tension
volScalarField activeTension
(
    IOobject
    (
        "activeTension",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimPressure, 0.0),
    "zeroGradient"
);

// Read active tension model dict
const dictionary& activeTensionModelCoeffs =
    electroActivationProperties.subDict("activeTensionModelCoeffs");

// Create active tension model ODE system
NashPanfilovActiveTensionModel activeTensionModel
(
    activeTensionModelCoeffs, mesh.nCells()
);
