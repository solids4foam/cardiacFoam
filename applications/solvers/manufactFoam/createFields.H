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
        dimensionedScalar("Vm", dimVoltage, 0),
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

// Create commonly used dimensionSets for convenience
// Read electoActivationProperties
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
        dimensionedScalar("zero", dimCurrent/dimArea, 0.0),
        "zeroGradient"
    );


Info<< "Reading cellular surface-to-volume ratio \n" << endl;
const dimensionedScalar chi
(
    "chi", dimArea/dimVolume, electroActivationProperties
);

Info<< "chi = " << chi.value() << chi.dimensions() << endl;

Info<< "Reading membrane capacitance Cm\n" << endl;
const dimensionedScalar Cm
(
    "Cm",
    dimCurrent*dimTime/(dimVoltage*dimArea),
    electroActivationProperties
);
Info<< "Cm = " << Cm.value()  << Cm.dimensions()  << endl;

// Read ionic model dict
const dictionary& ionicModelCoeffs =
electroActivationProperties.subDict("ionicModelCoeffs");

// Create ionicModelCellML object
autoPtr<ionicModelFDA> ionicModelFDA =
ionicModelFDA::New(ionicModelCoeffs, mesh.nCells(), runTime.deltaTValue());





