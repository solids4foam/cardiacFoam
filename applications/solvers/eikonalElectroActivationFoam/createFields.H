// Define voltage dimesions for convenience
const dimensionSet dimVoltage(dimMass*dimArea/(pow3(dimTime)*dimCurrent));

// Activation time
volScalarField psi
(
    IOobject
    (
        "psi",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("psi", dimTime, 0.0),
    "zeroGradient"
);

// Gradient of activation time
volVectorField gradPsi
(
    IOobject
    (
        "grad(psi)",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::grad(psi)
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

// Read external stimulus from a bounding box
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

// Calculate the M tensor
const volTensorField M("M", conductivity/(chi*Cm));
Info<< "M tensor field:" << nl
    << "    max(M) = " << max(M).value() << nl
    << "    min(M) = " << min(M).value() << nl
    << "    average(M) = " << average(M).value() << endl;

// M & gradPsi term
volVectorField w("w", M & gradPsi);

// Define a small constant to add to the nonlinear term G
const dimensionedScalar smallG("smallG", dimTime, SMALL);

// Calculate the nonlinear term
volScalarField G("G", sqrt((gradPsi & w) + smallG));

// Read reference velocity
const dimensionedScalar c0("c0", cardiacProperties);

// Create unit dimensioned scalar
const dimensionedScalar one("one", dimless, 1.0);

// Term related to the linearisation of c0*G
volVectorField a("a", w/G);

// Term related to the linearisation of c0*G
volVectorField u("u", c0*a);

// Update the phiU term
surfaceScalarField phiU("phiU", (fvc::interpolate(u) & mesh.Sf()));

Info<< "Reading timeIntegrationProperties\n" << endl;
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

// Lookup solution algorithm
const Switch eikonalAdvectionDiffusionApproach
(
    timeIntegrationProperties.lookup("eikonalAdvectionDiffusionApproach")
);
