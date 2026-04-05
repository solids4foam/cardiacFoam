/*---------------------------------------------------------------------------*\
License
    This file is part of cardiacFoam.

    cardiacFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the Free
    Software Foundation, either version 3 of the License, or (at your option)
    any later version.

    cardiacFoam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
    Public License for more details.

    You should have received a copy of the GNU General Public License
    along with cardiacFoam.  If not, see <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "myocardiumDomain.H"
#include "fvc.H"

namespace Foam
{

namespace
{



} // End anonymous namespace


MyocardiumDomain::MyocardiumDomain
(
    const fvMesh& mesh,
    const dictionary& electroProperties,
    PtrList<volScalarField>& outFields,
    const wordList& postProcessFieldNames,
    PtrList<volScalarField>& postProcessFields,
    ionicModel& ionicModel,
    electroVerificationModel* verificationModelPtr,
    autoPtr<myocardiumSolver> diffusionSolverPtr
)
:
    electroDomainInterface(),
    diffusionSolverPtr_(diffusionSolverPtr),
    Vm_
    (
        IOobject
        (
            "Vm",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("Vm", dimVoltage, -0.084),
        "zeroGradient"
    ),
    sourceField_
    (
        IOobject
        (
            "externalStimulusCurrent",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimCurrent/dimVolume, 0.0),
        "zeroGradient"
    ),
    Iion_
    (
        IOobject
        (
            "ionicCurrent",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimVoltage/dimTime, 0.0),
        "zeroGradient"
    ),
    activationTime_
    (
        IOobject
        (
            "activationTime",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimTime, 0.0),
        "zeroGradient"
    ),
    calculateActivationTime_(mesh.nCells(), true),
    outFields_(outFields),
    preProcessFieldNames_(),
    preProcessFields_(),
    postProcessFieldNames_(postProcessFieldNames),
    postProcessFields_(postProcessFields),
    ionicModel_(ionicModel),
    verificationModelPtr_(verificationModelPtr),
    electroProperties_(electroProperties),
    chi_("chi", dimArea/dimVolume, electroProperties_),
    Cm_("Cm", dimCurrent*dimTime/(dimVoltage*dimArea), electroProperties_),
    monodomainStimulus_(stimulusIO::loadMonodomainStimulusProtocol
    (
        electroProperties_
    )),
    useExplicitAlgorithm_
    (
        electroProperties_.lookupOrDefault<word>
        (
            "solutionAlgorithm", "implicit"
        ) == "explicit"
    )
{
    if (diffusionSolverPtr_->phiEPtr())
    {
        bindBidomainField(*const_cast<volScalarField*>(diffusionSolverPtr_->phiEPtr()));
    }
    else
    {
        validateNoIonicStimulusInMonodomain();
    }

    initialiseProcessing();
}


void MyocardiumDomain::updateExternalStimulusCurrent
(
    volScalarField& externalStimulusCurrent,
    const MonodomainStimulusProtocol& externalStimulus,
    scalar t0
) const
{
    scalarField& externalStimulusCurrentI =
        externalStimulusCurrent.primitiveFieldRef();
    externalStimulusCurrentI = 0.0;

    const vectorField& centres = mesh().C().primitiveField();

    forAll(externalStimulus.boxes, stimI)
    {
        const boundBox& stimBox = externalStimulus.boxes[stimI];
        const scalar tStart = externalStimulus.startTimes[stimI];
        const scalar tEnd = tStart + externalStimulus.durations[stimI];

        if (t0 < tStart || t0 > tEnd)
        {
            continue;
        }

        forAll(centres, cellI)
        {
            if (stimBox.contains(centres[cellI]))
            {
                externalStimulusCurrentI[cellI] +=
                    externalStimulus.intensities[stimI];
            }
        }
    }

    externalStimulusCurrent.correctBoundaryConditions();
}


void MyocardiumDomain::updateActivationTime
(
    volScalarField& activationTime,
    boolList& calculateActivationTime,
    const volScalarField& Vm
) const
{
    const scalarField& VmI = Vm.primitiveField();
    const scalarField& VmOldI = Vm.oldTime().primitiveField();
    scalarField& activationTimeI = activationTime.primitiveFieldRef();

    const scalar oldTime =
        mesh().time().value() - mesh().time().deltaTValue();
    const scalar deltaT = mesh().time().deltaTValue();

    forAll(activationTimeI, cellI)
    {
        if (calculateActivationTime[cellI] && VmI[cellI] > SMALL)
        {
            calculateActivationTime[cellI] = false;

            const scalar w =
                (0.0 - VmOldI[cellI])/(VmI[cellI] - VmOldI[cellI]);

            activationTimeI[cellI] = oldTime + w*deltaT;
        }
    }

    activationTime.correctBoundaryConditions();
}


void MyocardiumDomain::validateNoIonicStimulusInMonodomain() const
{
    const Switch allowIonicStimulusInMonodomain =
        electroProperties_.lookupOrDefault<Switch>
        (
            "allowIonicStimulusInMonodomain",
            false
        );

    const StimulusProtocol& ionicStim = ionicModel_.stimulusProtocol();
    if
    (
        stimulusIO::hasActiveStimulus(ionicStim)
     && !allowIonicStimulusInMonodomain
    )
    {
        FatalErrorInFunction
            << "Detected active ionic-model stimulus protocol while using "
            << "a reaction-diffusion myocardium domain." << nl
            << "This can unintentionally combine ionic and PDE external "
            << "stimulation." << nl
            << "Either disable ionic protocol keys (stim_*, nstim*) or set "
            << "allowIonicStimulusInMonodomain true to override."
            << exit(FatalError);
    }
}


void MyocardiumDomain::bindBidomainField(volScalarField& phiE)
{
    if (verificationModelPtr_)
    {
        verificationModelPtr_->bindBidomainField(phiE);
    }
}


void MyocardiumDomain::initialiseProcessing()
{
    if (verificationModelPtr_)
    {
        preProcessFieldNames_ =
            verificationModelPtr_->preProcessFieldNames(ionicModel_);

        if (postProcessFieldNames_.empty())
        {
            postProcessFieldNames_ =
                verificationModelPtr_
                    ->requiredPostProcessFieldNames(ionicModel_);
        }
    }

    const wordList exportedNames = ionicModel_.exportedFieldNames();

    outFields_.setSize(exportedNames.size());
    forAll(exportedNames, i)
    {
        outFields_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    exportedNames[i],
                    mesh().time().timeName(),
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

    if (verificationModelPtr_)
    {
        electroVerificationModel::allocateFields
        (
            preProcessFieldNames_,
            mesh(),
            "preProcess_",
            preProcessFields_
        );

        electroVerificationModel::allocateFields
        (
            postProcessFieldNames_,
            mesh(),
            "postProcess_",
            postProcessFields_
        );

        Info << "Using verification model "
             << verificationModelPtr_->type() << "." << endl;

        if (!preProcessFieldNames_.empty())
        {
            Info << "Running preProcess with fields " << preProcessFieldNames_
                 << "." << endl;
        }

        verificationModelPtr_->preProcess(ionicModel_, Vm_, preProcessFields_);
    }
}


void MyocardiumDomain::advance(scalar t0, scalar dt)
{
    advance(t0, dt, nullptr);
}


void MyocardiumDomain::advance
(
    scalar t0,
    scalar dt,
    pimpleControl* pimplePtr
)
{
    updateExternalStimulusCurrent(sourceField_, monodomainStimulus_, t0);

    ionicModel_.solveODE(t0, dt, Vm_, Iion_);
    Iion_.correctBoundaryConditions();

    if (pimplePtr)
    {
        diffusionSolverPtr_->solveDiffusionImplicit(*this, dt, *pimplePtr);
    }
    else if (useExplicitAlgorithm_)
    {
        diffusionSolverPtr_->solveDiffusionExplicit(*this, dt);
    }
    else
    {
        FatalErrorInFunction
            << "advance() requires a pimpleControl instance "
               "for implicit algorithms."
            << exit(FatalError);
    }

    updateActivationTime(activationTime_, calculateActivationTime_, Vm_);
}


scalar MyocardiumDomain::suggestExplicitDeltaT(scalar maxCo) const
{
    surfaceVectorField n("n", mesh().Sf());
    n /= mesh().magSf();

    const scalarField Df
    (
        (
            n
          & (n & fvc::interpolate
                (diffusionSolverPtr_->explicitConductivityField()))
        )
      /(chi_*Cm_)
    );

    const scalarField dx(1.0/mesh().deltaCoeffs());

    return maxCo*gMin(sqr(dx)/Df);
}


bool MyocardiumDomain::shouldPostProcess() const
{
    if (verificationModelPtr_)
    {
        return verificationModelPtr_->shouldPostProcess(ionicModel_, Vm_);
    }
    return false;
}


void MyocardiumDomain::exportStates()
{
    if (!outFields_.empty())
    {
        ionicModel_.exportStates(outFields_);
    }
}


void MyocardiumDomain::exportPostProcessFields()
{
    if (!postProcessFields_.empty())
    {
        ionicModel_.exportFields(postProcessFieldNames_, postProcessFields_);
    }
}


void MyocardiumDomain::postProcess()
{
    if (!shouldPostProcess())
    {
        return;
    }

    exportStates();
    exportPostProcessFields();
    if (verificationModelPtr_)
    {
        verificationModelPtr_->postProcess
        (
            ionicModel_, Vm_, postProcessFields_
        );
    }
}

} // End namespace Foam

// ************************************************************************* //
