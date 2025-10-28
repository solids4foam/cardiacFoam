const dimensionSet dimVoltage(dimMass*dimArea/(pow3(dimTime)*dimCurrent));

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
    dimensionedScalar("Vm", dimVoltage, -0.084),
    "zeroGradient"
);

Info<< "Setting initial conditions for V" << endl;
computeManufacturedV
(
    Vm.primitiveFieldRef(),
    mesh.C().component(vector::X),
    runTime.value()
);
Vm.correctBoundaryConditions();

// --- Create fields
// PC: why do we need u1, u2 and u3 at the solver level?
Info<< "Creating fields for ionic model gating variables\n" << endl;

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

Info<< "Setting initial conditions for u1, u2 and u3" << endl;
computeManufacturedU
(
    u1.primitiveFieldRef(),
    u2.primitiveFieldRef(),
    u3.primitiveFieldRef(),
    mesh.C().component(vector::X),
    runTime.value()
);
u1.correctBoundaryConditions();
u2.correctBoundaryConditions();
u3.correctBoundaryConditions();


// ************************************************************************* //

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

Info<< "Reading electric Conductivity Tensor" << endl;
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
    Info<< "conductivity not found on disk, using conductivity from "
        << electroActivationProperties.name() << nl << endl;

    conductivity =
        dimensionedTensor
        (
            dimensionedSymmTensor
            (
                "conductivity",
                pow3(dimTime)*sqr(dimCurrent)/(dimMass*dimVolume),
                electroActivationProperties
            ) & tensor(I)
        );
}
else
{
    Info<< "conductivity field read from " << runTime.timeName() << nl << endl;
}

Info<< "Reading cellular surface-to-volume ratio \n" << endl;
const dimensionedScalar chi
(
    "Chi", dimArea/dimVolume, electroActivationProperties
);

Info<< "Reading membrane capacitance Cm\n" << endl;
const dimensionedScalar Cm
(
    "Cm",
    dimCurrent*dimTime/(dimVoltage*dimArea),
    electroActivationProperties
);


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
    dimensionedScalar("zero",dimCurrent/dimArea, 0.0),
    "zeroGradient"
);


// Read ionic model dict
const dictionary& ionicModelCoeffs =
    electroActivationProperties.subDict("ionicModelCoeffs");

// Create ionicModelCellML object
autoPtr<ionicModelFDA> ionicModelFDA =
    ionicModelFDA::New(ionicModelCoeffs, mesh.nCells(), runTime.deltaTValue());





