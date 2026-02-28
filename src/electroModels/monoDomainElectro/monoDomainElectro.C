/*---------------------------------------------------------------------------*\
License
    This file is part of solids4foam.

    solids4foam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    solids4foam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with solids4foam.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "monoDomainElectro.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"
#include "fvm.H"
#include "IOmanip.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace electroModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(monoDomainElectro, 0);
addToRunTimeSelectionTable(electroModel, monoDomainElectro, dictionary);

// * * * * * * * * * * * * * * * Private Members * * * * * * * * * * * * * * //

tmp<volTensorField> monoDomainElectro::initialiseConductivity() const
{
    // Prepare the result
    tmp<volTensorField> tresult
    (
        new volTensorField
        (
            IOobject
            (
                "conductivity",
                runTime().timeName(),
                mesh(),
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedTensor
            (
                "zero",
                pow3(dimTime)*sqr(dimCurrent)/(dimMass*dimVolume),
                tensor::zero
            )
        )
    );
    volTensorField& result = tresult.ref();

    // Check if we actually read it, otherwise assign constant value
    if (!result.headerOk())
    {
        Info<< "\nconductivity not found on disk, using conductivity from "
            << cardiacProperties_.name() << nl << endl;

        result =
            dimensionedTensor
            (
                dimensionedSymmTensor
                (
                    "conductivity",
                    pow3(dimTime)*sqr(dimCurrent)/(dimMass*dimVolume),
                    cardiacProperties_
                ) & tensor(I)
            );

        if (result.size() > 0)
        {
            Info<< "Conductivity tensor (cell 0): " << result[0] << nl;
        }
    }
    else
    {
        Info<< "conductivity field read from " << runTime().timeName() << nl
            << endl;
    }

    return tresult;
}


void monoDomainElectro::updateExternalStimulusCurrent
(
    volScalarField& externalStimulusCurrent,
    const List<stimulus>& externalStimulus,
    const scalar t0
) const
{
    scalarField& externalStimulusCurrentI = externalStimulusCurrent;
    externalStimulusCurrentI = 0.0;

    forAll(externalStimulus, stimI)
    {
        const stimulus& curStim = externalStimulus[stimI];

        const scalar tStart = curStim.startTime_;
        if (t0 < tStart || t0 > (tStart + curStim.duration_))
        {
            continue;
        }

        const labelList& stimulusCellIDs = curStim.cellIDs_;
        forAll(stimulusCellIDs, cI)
        {
            const label id = stimulusCellIDs[cI];
            externalStimulusCurrentI[id] = curStim.intensity_;
        }
    }

    externalStimulusCurrent.correctBoundaryConditions();
}


bool monoDomainElectro::evolveExplicit()
{
    if (time().timeIndex() == 1)
    {
        Info<< "Solving the mono-domain equation for Vm using an explicit "
            << "approach" << nl
            << setw(20) << "Simulation Time"
            << setw(20) << "Clock Time"
            << setw(20) << "Min Vm"
            << setw(20) << "Max Vm" << endl;
    }

    physicsModel::printInfo() = bool
    (
        time().timeIndex() % infoFrequency_ == 0
     || mag(time().value() - time().endTime().value()) < SMALL
    );

    if (physicsModel::printInfo())
    {
        Info<< setw(20) <<time().value()
            << setw(20) << time().elapsedClockTime()
            << setw(20) << min(Vm_).value()
            << setw(20) << max(Vm_).value() << endl;

        physicsModel::printInfo() = false;
    }

    // Disable OpenFOAM linear solver output
    const label debugOrg = SolverPerformance<scalar>::debug;
    SolverPerformance<scalar>::debug = 0;

    // Current time step
    const scalar dt = runTime().deltaTValue();

    // Old time
    const scalar t0 = runTime().value() - dt;

    // 1) Set the external stimulus
    updateExternalStimulusCurrent
    (
        externalStimulusCurrent_, externalStimulus_, t0
    );

    // 2) Update the ionic current using OLD Vm
    ionicModelPtr_->calculateCurrent
    (
        t0,
        dt,
        Vm_.oldTime(),
        Iion_,
        states_
    );
    Iion_.correctBoundaryConditions();

    // 3) Solve the Vm equation
    solve
    (
        chi_*Cm_*fvm::ddt(Vm_)
     == fvc::laplacian(conductivity_, Vm_)
      - chi_*Cm_*Iion_
      + externalStimulusCurrent_
    );

    // 4) Advance ionic model in time (ODE solve) with NEW Vm, once per timestep
    // PHILIP -> we can merge 2 and 4, i.e. update current once per time step
    ionicModelPtr_->solveODE
    (
        t0,
        dt,
        Vm_,    // Vm_new
        Iion_,
        states_
    );

    // 5) Update post-processing fields
    if (runTime().outputTime())
    {
        // Extract states for visualisation
        ionicModelPtr_->exportStates(states_, outFields_);

        // Update activationTime field
        const scalarField& VmI = Vm_;
        const scalarField& VmOldI = Vm_.oldTime();
        boolList& calculateActivationTimeI = calculateActivationTime_;
        scalarField& activationTimeI = activationTime_.primitiveFieldRef();
        const scalar oldTime = runTime().value() - runTime().deltaTValue();
        const scalar deltaT = runTime().deltaTValue();
        forAll(activationTimeI, cellI)
        {
            if (calculateActivationTimeI[cellI])
            {
                if (VmI[cellI] > SMALL)
                {
                    calculateActivationTimeI[cellI] = false;

                    // Linearly interpolate for more accuracy
                    const scalar w =
                        (0.0 - VmOldI[cellI])/(VmI[cellI] - VmOldI[cellI]);

                    activationTimeI[cellI] = oldTime + w*deltaT;
                }
            }
        }

        activationTime_.correctBoundaryConditions();
    }

    // Re-enable OpenFOAM linear solver output
    SolverPerformance<scalar>::debug = debugOrg;

    return true;
}


bool monoDomainElectro::evolveImplicit()
{
    if (time().timeIndex() == 1)
    {
        Info<< "Solving the mono-domain equation for Vm using an implicit "
            << "approach" << endl;
    }

    // Current time step
    const scalar dt = runTime().deltaTValue();

    // Old time
    const scalar t0 = runTime().value() - dt;

    // 1) Set the external stimulus
    updateExternalStimulusCurrent
    (
        externalStimulusCurrent_, externalStimulus_, t0
    );

    // 2) Update the ionic current using OLD Vm
    ionicModelPtr_->calculateCurrent
    (
        t0,
        dt,
        Vm_.oldTime(),
        Iion_,
        states_
    );
    Iion_.correctBoundaryConditions();

    // 3) Solve the Vm equation implicitly - no ODE update inside
    while (pimple().loop())
    {
        solve
        (
            chi_*Cm_*fvm::ddt(Vm_)
         == fvm::laplacian(conductivity_, Vm_)
          - chi_*Cm_*Iion_
          + externalStimulusCurrent_
        );
    }

    // 4) Advance ionic model in time (ODE solve) with NEW Vm, once per timestep
    // PHILIP -> we can merge 2 and 4, i.e. update current once per time step
    ionicModelPtr_->solveODE
    (
        t0,
        dt,
        Vm_,    // Vm_new
        Iion_,
        states_
    );

    // 5) Update post-processing fields
    if (runTime().outputTime())
    {
        // Extract states for visualisation
        ionicModelPtr_->exportStates(states_, outFields_);

        // Update activationTime field
        const scalarField& VmI = Vm_;
        const scalarField& VmOldI = Vm_.oldTime();
        boolList& calculateActivationTimeI = calculateActivationTime_;
        scalarField& activationTimeI = activationTime_.primitiveFieldRef();
        const scalar oldTime = runTime().value() - runTime().deltaTValue();
        const scalar deltaT = runTime().deltaTValue();
        forAll(activationTimeI, cellI)
        {
            if (calculateActivationTimeI[cellI])
            {
                if (VmI[cellI] > SMALL)
                {
                    calculateActivationTimeI[cellI] = false;

                    // Linearly interpolate for more accuracy
                    const scalar w =
                        (0.0 - VmOldI[cellI])/(VmI[cellI] - VmOldI[cellI]);

                    activationTimeI[cellI] = oldTime + w*deltaT;
                }
            }
        }

        activationTime_.correctBoundaryConditions();
    }

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

monoDomainElectro::monoDomainElectro
(
    Time& runTime,
    const word& region
)
:
    electroModel(typeName, runTime, region),
    Vm_
    (
        IOobject
        (
            "Vm",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedScalar("Vm", dimVoltage, -0.084),
        "zeroGradient"
    ),
    Iion_
    (
        IOobject
        (
            "ionicCurrent",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedScalar("zero", dimVoltage/ dimTime, 0.0),
        "zeroGradient"
    ),
    externalStimulusCurrent_
    (
        IOobject
        (
            "externalStimulusCurrent",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedScalar("zero", dimCurrent/dimVolume, 0.0),
        "zeroGradient"
    ),
    externalStimulus_(),
    // cardiacProperties_(electroProperties().subDict("cardiacProperties")),
    cardiacProperties_(electroProperties()),
    chi_("chi", dimArea/dimVolume, cardiacProperties_),
    Cm_("Cm",dimCurrent*dimTime/(dimVoltage*dimArea), cardiacProperties_),
    conductivity_(initialiseConductivity()),
    activationTime_
    (
        IOobject
        (
            "activationTime",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedScalar("zero", dimTime, 0.0),
        "zeroGradient"
    ),
    calculateActivationTime_(mesh().nCells(), true),
    activationVelocity_
    (
        IOobject
        (
            "activationVelocity",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvc::grad
        (
            1.0/(activationTime_ + dimensionedScalar("SMALL", dimTime, SMALL))
        )
    ),
    outFields_(),
    ionicModelPtr_
    (
        ionicModel::New
        (
            electroProperties(),
            mesh().nCells(),
            runTime.deltaTValue()
        )
    ),
    states_
    (
        mesh().nCells(), scalarField(ionicModelPtr_->nEqns(), 0.0)
    ),
    setDeltaT_(true),
    infoFrequency_
    (
        electroProperties().lookupOrDefault<label>("infoFrequency", 1)
    )
{
    // Initialise the external stimulii
    const dictionary& stimulusProtocol =
        // electroProperties().subDict("stimulusProtocol");
        electroProperties();

    // Read external stimulus from one or more bounding boxes
    List<boundBox> stimulusBoxes;
    List<scalar> stimulusStartTimes;

    if
    (
        stimulusProtocol.found("stimulusLocationMinList")
     && stimulusProtocol.found("stimulusLocationMaxList")
    )
    {
        const List<point> mins
        (
            stimulusProtocol.lookup("stimulusLocationMinList")
        );
        const List<point> maxs
        (
            stimulusProtocol.lookup("stimulusLocationMaxList")
        );

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
        stimulusProtocol.lookup("stimulusStartTimeList")
            >> stimulusStartTimes;

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
    const vectorField& CI = mesh().C();
    forAll(CI, cellI)
    {
        const point& c = CI[cellI];
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

    const dimensionedScalar stimulusDuration
    (
        "stimulusDuration", dimTime, stimulusProtocol
    );

    const dimensionedScalar stimulusIntensity
    (
        "stimulusIntensity", dimCurrent/dimVolume, stimulusProtocol
    );

    // Set the number of stimulii
    externalStimulus_.setSize(stimulusBoxes.size());

    // Initialise each stimulii
    forAll(stimulusBoxes, bI)
    {
        externalStimulus_[bI].cellIDs_ = stimCellSets[bI].toc();
        externalStimulus_[bI].startTime_ = stimulusStartTimes[bI];
        externalStimulus_[bI].duration_ = stimulusDuration.value();
        externalStimulus_[bI].intensity_ = stimulusIntensity.value();
    }

    // Only print stimulus info if this is a 3D mesh
    if (mesh().nGeometricD() == 3)
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

    // Write parameters
    Info<< "Surface-to-volume ratio chi = " << chi_ << nl
        << "Membrane capacitance Cm = " << Cm_ << nl
        << "M tensor field:" << endl;

    // Set output fields
    const wordList names = ionicModelPtr_->exportedFieldNames();
    outFields_.setSize(names.size());
    forAll(names, i)
    {
        outFields_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    names[i],
                    runTime.timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh(),
                dimless,
                "zeroGradient"
            )
        );
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void monoDomainElectro::setDeltaT(Time& runTime)
{
    // Update the time step
    if (solutionAlg() == solutionAlgorithm::EXPLICIT && setDeltaT_)
    {
        setDeltaT_ = false;

        // Unit face normals
        surfaceVectorField n("n", mesh().Sf());
        n /= mesh().magSf();

        // Interpolate the conductivity to the faces and take the normal
        // component and divide by chi and Cm
        const scalarField Df
        (
            (n & (n & fvc::interpolate(conductivity_)))/(chi_*Cm_)
        );

        // Lookup the desired Courant number
        const scalar maxCo =
            runTime.controlDict().lookupOrDefault<scalar>("maxCo", 0.1);

        // Use the face delta coeffs as measures of the local cell size
        const scalarField dx(1.0/mesh().deltaCoeffs());

        // Compute stable dt for every face and take the minium
        scalar newDeltaT = maxCo*gMin(sqr(dx)/Df);

        Info<< "Setting deltaT = " << newDeltaT
            << ", maxCo = " << maxCo << endl;

        runTime.setDeltaT(newDeltaT);
    }
}


bool monoDomainElectro::evolve()
{
    if (solutionAlg() == solutionAlgorithm::IMPLICIT)
    {
        return evolveImplicit();
    }
    else if (solutionAlg() == solutionAlgorithm::EXPLICIT)
    {
        return evolveExplicit();
    }
    else
    {
        FatalErrorInFunction
            << "Unrecognised solution algorithm. Available options are "
            << electroModel::solutionAlgorithmNames_
               [
                   electroModel::solutionAlgorithm::IMPLICIT
               ]
            << ", "
            << electroModel::solutionAlgorithmNames_
               [
                   electroModel::solutionAlgorithm::EXPLICIT
               ]
            << endl;
    }

    // Keep compiler happy
    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace electroModels
} // End namespace Foam

// ************************************************************************* //
