/*---------------------------------------------------------------------------*\
License
    This file is part of cardiacFoam.

    cardiacFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    cardiacFoam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with cardiacFoam.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "myocardiumDomain.H"
#include "fvc.H"
#include "nonOrthogonalCorrectorLoop.H"

namespace Foam
{

namespace
{

autoPtr<fvMeshSubset> createMyocardiumMeshSubset
(
    const fvMesh& supportMesh,
    const dictionary& electroProperties
)
{
    if (!electroProperties.found("cellZone"))
    {
        return autoPtr<fvMeshSubset>(nullptr);
    }

    const word cellZoneName(electroProperties.lookup("cellZone"));
    const label zoneId = supportMesh.cellZones().findZoneID(cellZoneName);

    if (zoneId < 0)
    {
        FatalErrorInFunction
            << "Cannot find myocardium cellZone '" << cellZoneName
            << "' on mesh '" << supportMesh.name() << "'."
            << exit(FatalError);
    }

    autoPtr<fvMeshSubset> subsetPtr(new fvMeshSubset(supportMesh));
    subsetPtr->setCellSubset(supportMesh.cellZones()[zoneId]);
    return subsetPtr;
}


const fvMesh& resolveMyocardiumMesh
(
    const fvMesh& supportMesh,
    const autoPtr<fvMeshSubset>& subsetPtr
)
{
    return
    (
        subsetPtr.valid() && subsetPtr->hasSubMesh()
      ? subsetPtr->subMesh()
      : supportMesh
    );
}


void writeMappedCellField
(
    const volScalarField& subField,
    const fvMesh& supportMesh,
    const labelUList& cellMap
)
{
    volScalarField fullField
    (
        IOobject
        (
            subField.name(),
            supportMesh.time().timeName(),
            supportMesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        supportMesh,
        dimensionedScalar("zero", subField.dimensions(), 0.0),
        "zeroGradient"
    );

    scalarField& fullValues = fullField.primitiveFieldRef();
    const scalarField& subValues = subField.primitiveField();

    forAll(cellMap, subCelli)
    {
        fullValues[cellMap[subCelli]] = subValues[subCelli];
    }

    fullField.write();
}


} // End anonymous namespace


autoPtr<myocardiumDomain> myocardiumDomain::New
(
    const fvMesh& supportMesh,
    const dictionary& electroProperties,
    PtrList<volScalarField>& outFields,
    const wordList& postProcessFieldNames,
    PtrList<volScalarField>& postProcessFields,
    ionicModel& ionicModel,
    electroVerificationModel* verificationModelPtr
)
{
    const word& cn = electroProperties.dictName();
    const word solverType =
        cn.endsWith("Coeffs") ? word(cn.substr(0, cn.size() - 6)) : cn;

    autoPtr<fvMeshSubset> meshSubsetPtr
    (
        createMyocardiumMeshSubset(supportMesh, electroProperties)
    );
    const fvMesh& myocardiumMesh =
        resolveMyocardiumMesh(supportMesh, meshSubsetPtr);

    return autoPtr<myocardiumDomain>
    (
        new myocardiumDomain
        (
            supportMesh,
            electroProperties,
            outFields,
            postProcessFieldNames,
            postProcessFields,
            ionicModel,
            verificationModelPtr,
            myocardiumSolver::New
            (
                myocardiumMesh,
                supportMesh,
                meshSubsetPtr.valid() ? &meshSubsetPtr() : nullptr,
                solverType,
                electroProperties
            ),
            meshSubsetPtr
        )
    );
}


label myocardiumDomain::configuredCellCount
(
    const fvMesh& mesh,
    const dictionary& electroProperties
)
{
    if (!electroProperties.found("cellZone"))
    {
        return mesh.nCells();
    }

    const word cellZoneName(electroProperties.lookup("cellZone"));
    const label zoneId = mesh.cellZones().findZoneID(cellZoneName);

    if (zoneId < 0)
    {
        FatalErrorInFunction
            << "Cannot find myocardium cellZone '" << cellZoneName
            << "' on mesh '" << mesh.name() << "'."
            << exit(FatalError);
    }

    return mesh.cellZones()[zoneId].size();
}


myocardiumDomain::myocardiumDomain
(
    const fvMesh& supportMesh,
    const dictionary& electroProperties,
    PtrList<volScalarField>& outFields,
    const wordList& postProcessFieldNames,
    PtrList<volScalarField>& postProcessFields,
    ionicModel& ionicModel,
    electroVerificationModel* verificationModelPtr,
    autoPtr<myocardiumSolver> diffusionSolverPtr,
    autoPtr<fvMeshSubset> meshSubsetPtr
)
:
    meshSubsetPtr_(meshSubsetPtr),
    supportMesh_(supportMesh),
    diffusionSolverPtr_(diffusionSolverPtr),
    Vm_
    (
        IOobject
        (
            "Vm",
            resolveMyocardiumMesh(supportMesh_, meshSubsetPtr_).time().timeName(),
            resolveMyocardiumMesh(supportMesh_, meshSubsetPtr_),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        resolveMyocardiumMesh(supportMesh_, meshSubsetPtr_),
        dimensionedScalar("Vm", dimVoltage, -0.084),
        "zeroGradient"
    ),
    gradVm_
    (
        IOobject
        (
            "grad(" + Vm_.name() + ")",
            resolveMyocardiumMesh(supportMesh_, meshSubsetPtr_).time().timeName(),
            resolveMyocardiumMesh(supportMesh_, meshSubsetPtr_),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        resolveMyocardiumMesh(supportMesh_, meshSubsetPtr_),
        dimensionedVector("0", dimVoltage/dimLength, vector::zero)
    ),
    sourceField_
    (
        IOobject
        (
            "externalStimulusCurrent",
            resolveMyocardiumMesh(supportMesh_, meshSubsetPtr_).time().timeName(),
            resolveMyocardiumMesh(supportMesh_, meshSubsetPtr_),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        resolveMyocardiumMesh(supportMesh_, meshSubsetPtr_),
        dimensionedScalar("zero", dimCurrent/dimVolume, 0.0),
        "zeroGradient"
    ),
    implicitSourceCoeff_
    (
        IOobject
        (
            "implicitSourceCoeff",
            resolveMyocardiumMesh(supportMesh_, meshSubsetPtr_).time().timeName(),
            resolveMyocardiumMesh(supportMesh_, meshSubsetPtr_),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        resolveMyocardiumMesh(supportMesh_, meshSubsetPtr_),
        dimensionedScalar
        (
            "zero",
            dimCurrent/(dimVolume*dimVoltage),
            0.0
        ),
        "zeroGradient"
    ),
    Iion_
    (
        IOobject
        (
            "ionicCurrent",
            resolveMyocardiumMesh(supportMesh_, meshSubsetPtr_).time().timeName(),
            resolveMyocardiumMesh(supportMesh_, meshSubsetPtr_),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        resolveMyocardiumMesh(supportMesh_, meshSubsetPtr_),
        dimensionedScalar("zero", dimVoltage/dimTime, 0.0),
        "zeroGradient"
    ),
    IionOld_
    (
        IOobject
        (
            "ionicCurrentOld",
            resolveMyocardiumMesh(supportMesh_, meshSubsetPtr_).time().timeName(),
            resolveMyocardiumMesh(supportMesh_, meshSubsetPtr_),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        resolveMyocardiumMesh(supportMesh_, meshSubsetPtr_),
        dimensionedScalar("zero", dimVoltage/dimTime, 0.0),
        "zeroGradient"
    ),
    IionOldOld_
    (
        IOobject
        (
            "ionicCurrentOldOld",
            resolveMyocardiumMesh(supportMesh_, meshSubsetPtr_).time().timeName(),
            resolveMyocardiumMesh(supportMesh_, meshSubsetPtr_),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        resolveMyocardiumMesh(supportMesh_, meshSubsetPtr_),
        dimensionedScalar("zero", dimVoltage/dimTime, 0.0),
        "zeroGradient"
    ),
    VmPrev_
    (
        IOobject
        (
            "VmPrevious",
            resolveMyocardiumMesh(supportMesh_, meshSubsetPtr_).time().timeName(),
            resolveMyocardiumMesh(supportMesh_, meshSubsetPtr_),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Vm_
    ),
    VmRate_
    (
        IOobject
        (
            "VmRate",
            resolveMyocardiumMesh(supportMesh_, meshSubsetPtr_).time().timeName(),
            resolveMyocardiumMesh(supportMesh_, meshSubsetPtr_),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        resolveMyocardiumMesh(supportMesh_, meshSubsetPtr_),
        dimensionedScalar("zero", dimVoltage/dimTime, 0.0),
        "zeroGradient"
    ),
    activationTime_
    (
        IOobject
        (
            "activationTime",
            resolveMyocardiumMesh(supportMesh_, meshSubsetPtr_).time().timeName(),
            resolveMyocardiumMesh(supportMesh_, meshSubsetPtr_),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        resolveMyocardiumMesh(supportMesh_, meshSubsetPtr_),
        dimensionedScalar("unactivated", dimTime, -1.0),
        "zeroGradient"
    ),
    outFields_(outFields),
    preProcessFieldNames_(),
    preProcessFields_(),
    postProcessFieldNames_(postProcessFieldNames),
    postProcessFields_(postProcessFields),
    ionicModel_(ionicModel),
    verificationModelPtr_(verificationModelPtr),
    electroProperties_(electroProperties),
    chi_("chi", dimArea/dimVolume, electroProperties_),
    Cm_("cm", dimCurrent*dimTime/(dimVoltage*dimArea), electroProperties_),
    externalStimulus_(stimulusIO::loadExternalStimulusProtocol
    (
        electroProperties_
    )),
    useExplicitAlgorithm_
    (
        electroProperties_.lookupOrDefault<word>
        (
            "solutionAlgorithm", "implicit"
        ) == "explicit"
    ),
    timeCouplingScheme_
    (
        electroProperties_.lookupOrDefault<word>("timeCouplingScheme", "godunov")
    ),
    activationThreshold_
    (
        electroProperties_.lookupOrDefault<scalar>("activationThreshold", 0.0)
    ),
    setDeltaT_(true)
{
    if (timeCouplingScheme_ != "godunov" && timeCouplingScheme_ != "sbdf2")
    {
        FatalErrorInFunction
            << "timeCouplingScheme must be 'godunov' or 'sbdf2'; got "
            << timeCouplingScheme_
            << exit(FatalError);
    }

    if
    (
        usesSbdf2Scheme()
     && (
            !ionicModel_.supportsVmRateCoupling()
         || !ionicModel_.supportsIonicCurrentEvaluation()
        )
    )
    {
        FatalErrorInFunction
            << "timeCouplingScheme 'sbdf2' requires an ionic model that "
            << "supports both second-order couplings, but "
            << ionicModel_.type()
            << " reports supportsVmRateCoupling="
            << Switch(ionicModel_.supportsVmRateCoupling())
            << " supportsIonicCurrentEvaluation="
            << Switch(ionicModel_.supportsIonicCurrentEvaluation())
            << ". Use timeCouplingScheme 'godunov' with this model."
            << exit(FatalError);
    }

    if (diffusionSolverPtr_->phiEPtr())
    {
        bindBidomainField(*const_cast<volScalarField*>(diffusionSolverPtr_->phiEPtr()));
    }
    else
    {
        validateNoIonicStimulusInMonodomain();
    }

    if (verificationModelPtr_)
    {
        verificationModelPtr_->bindSourceField(sourceField_);
    }

    initialiseProcessing();
}


void myocardiumDomain::updateExternalStimulusCurrent
(
    volScalarField& externalStimulusCurrent,
    const ExternalStimulusProtocol& externalStimulus,
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


void myocardiumDomain::updateGradVm()
{
    gradVm_ = fvc::grad(Vm_);
}


void myocardiumDomain::updateActivationTime()
{
    const scalarField& VmI = Vm_.primitiveField();
    const scalarField& VmOldI = Vm_.oldTime().primitiveField();
    scalarField& activationTimeI = activationTime_.primitiveFieldRef();

    const scalar oldTime =
        mesh().time().value() - mesh().time().deltaTValue();
    const scalar deltaT = mesh().time().deltaTValue();

    forAll(activationTimeI, cellI)
    {
        if
        (
            VmOldI[cellI] <= activationThreshold_
         && VmI[cellI] > activationThreshold_
        )
        {
            const scalar w =
                (activationThreshold_ - VmOldI[cellI])
               /(VmI[cellI] - VmOldI[cellI]);

            activationTimeI[cellI] = oldTime + w*deltaT;
        }
    }

    activationTime_.correctBoundaryConditions();
}


void myocardiumDomain::validateNoIonicStimulusInMonodomain() const
{
    const StimulusProtocol& ionicStim = ionicModel_.stimulusProtocol();
    if
    (
        stimulusIO::hasActiveStimulus(ionicStim)
    )
    {
        FatalErrorInFunction
            << "Detected active ionic-model stimulus protocol while using "
            << "a reaction-diffusion myocardium domain." << nl
            << "This can unintentionally combine ionic and PDE external "
            << "stimulation." << nl
            << "Disable the ionic protocol keys (stim_*, nstim*) and use "
            << "externalStimulus for PDE-scale stimulation instead."
            << exit(FatalError);
    }
}


void myocardiumDomain::bindBidomainField(volScalarField& phiE)
{
    if (verificationModelPtr_)
    {
        verificationModelPtr_->bindBidomainField(phiE);
    }
}


void myocardiumDomain::initialiseProcessing()
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

        verificationModelPtr_->preProcess(ionicModel_, Vm_, preProcessFields_);
    }
}


void myocardiumDomain::advance(scalar t0, scalar dt)
{
    advance(t0, dt, nullptr);
}


void myocardiumDomain::prepareTimeStep(scalar t0, scalar dt)
{
    // Reset the source field and apply the 3D external box stimulus.
    // Called by the advance scheme BEFORE any domain coupling deposits
    // current into sourceField_.  Keeping this here ensures that coupling
    // current added by preparePrimaryCoupling survives into the FVM solve.
    updateExternalStimulusCurrent(sourceField_, externalStimulus_, t0);

    IionOldOld_ = IionOld_;
    IionOld_ = Iion_;

    if (usesSbdf2Scheme())
    {
        const scalar deltaT0 = mesh().time().deltaT0Value();

        if (usesSbdf2VmRate() && deltaT0 > VSMALL)
        {
            VmRate_ =
                (Vm_ - VmPrev_)
              / dimensionedScalar("deltaT0", dimTime, deltaT0);
        }
        else
        {
            VmRate_ =
                dimensionedScalar("zero", VmRate_.dimensions(), 0.0);
        }

        VmPrev_ = Vm_;
    }

    if (verificationModelPtr_)
    {
        const volTensorField* conductivityPtr =
            diffusionSolverPtr_->conductivityPtr();

        if (conductivityPtr)
        {
            verificationModelPtr_->addManufacturedPdeSource
            (
                sourceField_,
                *conductivityPtr,
                t0 + dt
            );
        }
    }

    implicitSourceCoeff_ = dimensionedScalar
    (
        "zero",
        implicitSourceCoeff_.dimensions(),
        0.0
    );
    implicitSourceCoeff_.correctBoundaryConditions();
}


void myocardiumDomain::solveIonicCurrent(scalar t0, scalar dt)
{
    ionicModel_.solveODE(t0, dt, Vm_, Iion_);

    Iion_.correctBoundaryConditions();
}


void myocardiumDomain::refreshIonicCurrent(scalar t)
{
    ionicModel_.evaluateIonicCurrent(t, Vm_, Iion_);
    Iion_.correctBoundaryConditions();
}


void myocardiumDomain::advance
(
    scalar t0,
    scalar dt,
    pimpleControl* pimplePtr
)
{
    // sourceField_ was already set by prepareTimeStep (external stimulus)
    // and then augmented by the PVJ coupler (preparePrimaryCoupling).
    // Do NOT reset it here.

    if (usesSbdf2VmRate() && ionicModel_.supportsVmRateCoupling())
    {
        ionicModel_.setVmRate(VmRate_);
    }

    solveIonicCurrent(t0, dt);
    ionicModel_.clearVmRate();

    updateGradVm();

    if (useExplicitAlgorithm_)
    {
        diffusionSolverPtr_->solveDiffusionExplicit(*this, dt);
    }
    else if (pimplePtr)
    {
        diffusionSolverPtr_->solveDiffusionImplicit(*this, dt, *pimplePtr);
    }
    else
    {
        FatalErrorInFunction
            << "advance() requires a pimpleControl instance "
               "for implicit algorithms."
            << exit(FatalError);
    }

    // Move Iion onto its own time level now that Vm(t0+dt) is known
    if (usesSbdf2Scheme())
    {
        refreshIonicCurrent(t0 + dt);
    }

    updateActivationTime();
}


void myocardiumDomain::solveReactionStep(scalar t0, scalar dt)
{
    solveIonicCurrent(t0, dt);
}


void myocardiumDomain::solveDiffusionStep
(
    scalar t0,
    scalar dt,
    pimpleControl* pimplePtr
)
{
    (void)t0;

    updateGradVm();

    if (useExplicitAlgorithm_)
    {
        diffusionSolverPtr_->solveDiffusionExplicit(*this, dt);
    }
    else if (pimplePtr)
    {
        diffusionSolverPtr_->solveDiffusionImplicit(*this, dt, *pimplePtr);
    }
    else
    {
        FatalErrorInFunction
            << "solveDiffusionStep() requires a pimpleControl instance "
               "for implicit algorithms."
            << exit(FatalError);
    }
}


// Called twice per step by the bath-PDE predictor/corrector (the outer
// Vm(phiE^n) -> phiE^n+1 -> Vm(phiE^n+1) sweep owned by the advance scheme).
// Each call here is only the inner, single-stage non-orthogonal reassembly
// for that one Vm solve; it does not itself predict or correct anything.
void myocardiumDomain::solveDiffusionStepOnce
(
    scalar t0,
    scalar dt,
    pimpleControl* pimplePtr
)
{
    (void)t0;

    if (useExplicitAlgorithm_)
    {
        FatalErrorInFunction
            << "The bath-PDE predictor/corrector requires "
            << "solutionAlgorithm implicit."
            << exit(FatalError);
    }
    else if (pimplePtr)
    {
        correctNonOrthogonalLoop
        (
            *pimplePtr,
            [&]()
            {
                updateGradVm();
                diffusionSolverPtr_->solveDiffusionImplicit(*this, dt);
            }
        );
    }
    else
    {
        FatalErrorInFunction
            << "solveDiffusionStepOnce() requires a pimpleControl instance "
               "for implicit algorithms."
            << exit(FatalError);
    }
}


void myocardiumDomain::finalizeDiffusionStep()
{
    updateActivationTime();
}


void myocardiumDomain::bindExternalPhiE
(
    const volScalarField& phiE,
    const labelUList& heartCellMap
)
{
    diffusionSolverPtr_->bindExternalPhiE(phiE, heartCellMap);

    if (verificationModelPtr_)
    {
        verificationModelPtr_->bindBidomainField
        (
            const_cast<volScalarField&>(phiE),
            heartCellMap
        );

        verificationModelPtr_->preProcess(ionicModel_, Vm_, preProcessFields_);
    }

}


void myocardiumDomain::unbindExternalPhiE()
{
    diffusionSolverPtr_->unbindExternalPhiE();

    if (verificationModelPtr_)
    {
        if (const volScalarField* phiEPtr = diffusionSolverPtr_->phiEPtr())
        {
            verificationModelPtr_->bindBidomainField
            (
                *const_cast<volScalarField*>(phiEPtr)
            );
        }
        else
        {
            verificationModelPtr_->unbindBidomainField();
        }
    }

}


scalar myocardiumDomain::suggestExplicitDeltaT(scalar maxCo) const
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


bool myocardiumDomain::applyModelTimeControls(Time& runTime)
{
    if (!useExplicitAlgorithm_ || !setDeltaT_)
        return false;

    setDeltaT_ = false;

    const scalar maxCo =
        runTime.controlDict().lookupOrDefault<scalar>("maxCo", 0.1);

    const scalar stableMaxDeltaT = suggestExplicitDeltaT(maxCo);

    if (runTime.deltaTValue() > stableMaxDeltaT)
    {
        Info << "Capping deltaT from " << runTime.deltaTValue()
             << " to " << stableMaxDeltaT
             << " (explicit stability limit, maxCo = " << maxCo << ")"
             << endl;
        runTime.setDeltaT(stableMaxDeltaT);
        return true;
    }

    return false;
}


bool myocardiumDomain::shouldPostProcess() const
{
    if (verificationModelPtr_)
    {
        return verificationModelPtr_->shouldPostProcess(ionicModel_, Vm_);
    }
    return false;
}


void myocardiumDomain::exportStates()
{
    if (!outFields_.empty())
    {
        ionicModel_.exportStates(outFields_);
    }
}


void myocardiumDomain::exportPostProcessFields()
{
    if (!postProcessFields_.empty())
    {
        ionicModel_.exportFields(postProcessFieldNames_, postProcessFields_);
    }
}


void myocardiumDomain::write()
{
    if (meshSubsetPtr_.valid() && meshSubsetPtr_->hasSubMesh())
    {
        const labelUList& cellMap = meshSubsetPtr_->cellMap();

        writeMappedCellField(Vm_, supportMesh_, cellMap);
        writeMappedCellField(sourceField_, supportMesh_, cellMap);
        writeMappedCellField(implicitSourceCoeff_, supportMesh_, cellMap);
        writeMappedCellField(Iion_, supportMesh_, cellMap);
        writeMappedCellField(activationTime_, supportMesh_, cellMap);

        if (const volScalarField* phiEPtr = diffusionSolverPtr_->phiEPtr())
        {
            writeMappedCellField(*phiEPtr, supportMesh_, cellMap);
        }

        if (const volScalarField* phiIPtr = diffusionSolverPtr_->phiIPtr())
        {
            writeMappedCellField(*phiIPtr, supportMesh_, cellMap);
        }

        forAll(outFields_, i)
        {
            writeMappedCellField(outFields_[i], supportMesh_, cellMap);
        }

        forAll(preProcessFields_, i)
        {
            writeMappedCellField(preProcessFields_[i], supportMesh_, cellMap);
        }

        forAll(postProcessFields_, i)
        {
            writeMappedCellField(postProcessFields_[i], supportMesh_, cellMap);
        }

        return;
    }

    Vm_.write();
    sourceField_.write();
    Iion_.write();
    activationTime_.write();

    if (const volScalarField* phiEPtr = diffusionSolverPtr_->phiEPtr())
    {
        phiEPtr->write();
    }

    if (const volScalarField* phiIPtr = diffusionSolverPtr_->phiIPtr())
    {
        phiIPtr->write();
    }

    forAll(outFields_, i)
    {
        outFields_[i].write();
    }

    forAll(preProcessFields_, i)
    {
        preProcessFields_[i].write();
    }

    forAll(postProcessFields_, i)
    {
        postProcessFields_[i].write();
    }
}


void myocardiumDomain::postProcess()
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
