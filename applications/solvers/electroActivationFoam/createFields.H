// Define voltage dimesions for convenience
const dimensionSet dimVoltage(dimMass*dimArea/(pow3(dimTime)*dimCurrent));

// Trans-membrane voltage
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

// Ionic current
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
    dimensionedScalar("zero", dimCurrent/dimVolume, 0.0),
    "zeroGradient"
);

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

// Read external stimulus from one or more bounding boxes
List<boundBox> stimulusBoxes;
List<scalar> stimulusStartTimes;

if
(
    stimulusProtocol.found("stimulusLocationMinList")
 && stimulusProtocol.found("stimulusLocationMaxList")
)
{
    List<point> mins;
    List<point> maxs;
    stimulusProtocol.lookup("stimulusLocationMinList") >> mins;
    stimulusProtocol.lookup("stimulusLocationMaxList") >> maxs;

    if (mins.size() != maxs.size())
    {
        FatalErrorInFunction
            << "stimulusLocationMinList and stimulusLocationMaxList must "
            << "have the same size."
            << abort(FatalError);
    }

    stimulusBoxes.setSize(mins.size());
    forAll(mins, i)
    {
        stimulusBoxes[i] = boundBox(mins[i], maxs[i]);
    }
}
else
{
    stimulusBoxes.setSize(1);
    stimulusBoxes[0] = boundBox
    (
        point(stimulusProtocol.lookup("stimulusLocationMin")),
        point(stimulusProtocol.lookup("stimulusLocationMax"))
    );
}

if (stimulusProtocol.found("stimulusStartTimeList"))
{
    stimulusProtocol.lookup("stimulusStartTimeList") >> stimulusStartTimes;
    if (stimulusStartTimes.size() != stimulusBoxes.size())
    {
        FatalErrorInFunction
            << "stimulusStartTimeList must have the same size as "
            << "the stimulus box list."
            << abort(FatalError);
    }
}
else
{
    stimulusStartTimes.setSize(stimulusBoxes.size());
    const scalar startTime =
        stimulusProtocol.lookupOrDefault<scalar>("stimulusStartTime", 0.0);
    forAll(stimulusStartTimes, i)
    {
        stimulusStartTimes[i] = startTime;
    }
}

List<labelHashSet> stimCellSets(stimulusBoxes.size());
labelHashSet stimCellSet;

forAll(mesh.C(), cellI)
{
    const point& c = mesh.C()[cellI];
    forAll(stimulusBoxes, bI)
    {
        if (stimulusBoxes[bI].contains(c))
        {
            stimCellSets[bI].insert(cellI);
            stimCellSet.insert(cellI);
            break;
        }
    }
}
List<labelList> stimulusCellIDsList(stimulusBoxes.size());
forAll(stimulusBoxes, bI)
{
    stimulusCellIDsList[bI] = stimCellSets[bI].toc();
}

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
    Info<< "---------------------------------------------" << nl
        << "3D mesh detected — stimulus parameters:"       << nl
        << "Stimulus boxes: " << stimulusBoxes.size()      << nl
        << "Number of cells in stimulus region: "
        << stimCellSet.size()                              << nl
        << "Stimulus duration: "
        << stimulusDuration.value() << " "
        << stimulusDuration.dimensions()                   << nl
        << "Stimulus start time(s): " << stimulusStartTimes << nl
        << "Stimulus intensity: "
        << stimulusIntensity.value() << " "
        << stimulusIntensity.dimensions()                  << nl
        << "---------------------------------------------" << endl;
}

Info<< "Reading numerical solver properties\n" << endl;
IOdictionary timeIntegrationProperties
(
    IOobject
    (
        "timeIntegrationProperties",
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

Info<< "Reading surface-to-volume ratio chi = " << chi.value()
    <<"  "<< chi.dimensions() << endl;
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
    Info<< "\nconductivity not found on disk, using conductivity from "
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


// Merge the dictionaries and create a memory object to read all properties
dictionary solutionVariablesMemory(cardiacProperties);
solutionVariablesMemory.merge(timeIntegrationProperties);
solutionVariablesMemory.merge(stimulusProtocol);

// Create ionicModelCellML object
autoPtr<ionicModel> ionicModel =
    ionicModel::New
    (
        solutionVariablesMemory, mesh.nCells(), runTime.deltaTValue()
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

// Allocate volScalarFields for export
PtrList<volScalarField> outFields;
const wordList names = ionicModel->exportedFieldNames();
outFields.setSize(names.size());

forAll(names, i)
{
    outFields.set
    (
        i,
        new volScalarField
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
        )
    );
}
