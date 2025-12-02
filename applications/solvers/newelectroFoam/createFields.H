const dimensionSet dimVoltage(dimMass*dimArea/(pow3(dimTime)*dimCurrent));

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
    dimensionedScalar("Vm", dimVoltage, -0.084),
    "zeroGradient"
);

 volScalarField u1
 (
     IOobject
     (
         "u1",
         runTime.timeName(),
         mesh,
         IOobject::READ_IF_PRESENT,
         IOobject::AUTO_WRITE
     ),
     mesh,
     dimensionedScalar("u1", dimless, 0.0),
       "zeroGradient"
 );

 volScalarField u2
 (
     IOobject
     (
         "u2",
         runTime.timeName(),
         mesh,
         IOobject::READ_IF_PRESENT,
         IOobject::AUTO_WRITE
     ),
     mesh,
     dimensionedScalar("u2", dimless, 0.0),
       "zeroGradient"
 );

 volScalarField u3
 (
     IOobject
     (
         "u3",
         runTime.timeName(),
         mesh,
         IOobject::READ_IF_PRESENT,
         IOobject::AUTO_WRITE
     ),
     mesh,
     dimensionedScalar("u3", dimless, 0.0),
       "zeroGradient"
 );
 

// ************************************************************************* //

// Create dictionaries for cardiac properties, stimulus conditions and solver parameters


Info<< "Reading stimulus protocol\n" << endl;
IOdictionary stimulusProtocol
(
IOobject
(
    "stimulusProtocol",
    runTime.constant(),
    mesh,
    IOobject::MUST_READ_IF_MODIFIED,
    IOobject::NO_WRITE
)
);

//read external stimulus from a square bounding box
// Find cells in the stimulus volume
boundBox bb
(
point(stimulusProtocol.lookup("stimulusLocationMin")),
point(stimulusProtocol.lookup("stimulusLocationMax"))
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

const dimensionedScalar stimulusDuration
(
"stimulusDuration", dimTime, stimulusProtocol
);

const dimensionedScalar stimulusIntensity
(
"stimulusIntensity", dimCurrent/dimVolume, stimulusProtocol
);

// Only print stimulus info if this is a 3D mesh
if (mesh.nGeometricD() == 3)
{
    Info<< "---------------------------------------------" << nl;
    Info<< "3D mesh detected — stimulus parameters:"       << nl;
    Info<< "Stimulus bounding box: " << bb                  << nl;
    Info<< "Number of cells in stimulus region: "
        << stimulusCellIDs.size()                           << nl;

    Info<< "Stimulus duration: "
        << stimulusDuration.value() << " "
        << stimulusDuration.dimensions()                    << nl;

    Info<< "Stimulus intensity: "
        << stimulusIntensity.value() << " "
        << stimulusIntensity.dimensions()                   << nl;
    Info<< "---------------------------------------------" << endl;
}


Info<< "Reading numerical solver properties\n" << endl;
IOdictionary numericalsolverProperties
(
IOobject
(
    "numericalsolverProperties",
    runTime.constant(),
    mesh,
    IOobject::MUST_READ_IF_MODIFIED,
    IOobject::NO_WRITE
)
);
Info<< "---------------------------------------------" << nl;

Info<< "Reading tissue properties\n" << endl;
IOdictionary cardiacProperties
(
IOobject
(
    "cardiacProperties",
    runTime.constant(),
    mesh,
    IOobject::MUST_READ_IF_MODIFIED,
    IOobject::NO_WRITE
)
);

const dimensionedScalar chi
(
"chi", dimArea/dimVolume, cardiacProperties
);
Info<< "Reading surface-to-volume ratio chi = " << chi.value() <<"  "<< chi.dimensions() << endl;

const dimensionedScalar Cm
(
"Cm",
dimCurrent*dimTime/(dimVoltage*dimArea),
cardiacProperties
);
Info<< "Reading membrane capacitance Cm = " << Cm.value()  <<"  " << Cm.dimensions()  << endl;


Info<< "Reading conductivity D:" << endl;

volTensorField conductivity
(
IOobject
(
    "conductivity",
    runTime.timeName(),
    mesh,
    IOobject::READ_IF_PRESENT,
    IOobject::NO_WRITE
),
mesh,
dimensionedTensor
(
    "zero",
    pow3(dimTime)*sqr(dimCurrent)/(dimMass*dimVolume),
    tensor::zero
)
);

// Check if we actually read it, otherwise assign constant value
if (!conductivity.headerOk())
{
Info<< "\n conductivity not found on disk, using conductivity from "
    << cardiacProperties.name() << nl << endl;

conductivity =
    dimensionedTensor
    (
        dimensionedSymmTensor
        (
            "conductivity",
            pow3(dimTime)*sqr(dimCurrent)/(dimMass*dimVolume),
            cardiacProperties
        ) & tensor(I)
    );

    Info<< "Conductivity tensor (cell 0): " << conductivity[0] << nl;
}
else
{
Info<< "conductivity field read from " << runTime.timeName() << nl << endl;
}
Info<< "---------------------------------------------" << nl;

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


//Merge the dictionaries and create a memory object to read all properties.
dictionary solutionVariablesMemory(cardiacProperties);
solutionVariablesMemory.merge(numericalsolverProperties);
solutionVariablesMemory.merge(stimulusProtocol);

// Create ionicModelCellML object
autoPtr<ionicModel> ionicModel =
ionicModel::New(solutionVariablesMemory, mesh.nCells(), runTime.deltaTValue());



//Read objects for the monodomain equation: physical variables and material properties.

volScalarField Iion
(
IOobject
(
"ionicCurrent",
runTime.timeName(),
mesh,
IOobject::NO_READ,
IOobject::AUTO_WRITE
),
mesh,
dimensionedScalar("zero", dimVoltage/ dimTime, 0.0),
"zeroGradient"
);


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

// ************************************************************************* //
// Allocate volScalarFields for export
List<volScalarField*> outFields;
wordList names = ionicModel->exportedFieldNames();

outFields.setSize(names.size());

forAll(names, i)
{
    outFields[i] = new volScalarField
    (
        IOobject
        (
            names[i],
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionSet(0,0,0,0,0,0,0),
        zeroGradientFvPatchScalarField::typeName
    );
}






















