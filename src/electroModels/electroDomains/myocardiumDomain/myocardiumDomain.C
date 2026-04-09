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


autoPtr<MyocardiumDomain> MyocardiumDomain::New
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

    return autoPtr<MyocardiumDomain>
    (
        new MyocardiumDomain
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
                solverType,
                electroProperties
            ),
            meshSubsetPtr
        )
    );
}


label MyocardiumDomain::configuredCellCount
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


MyocardiumDomain::MyocardiumDomain
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
    electroDomainInterface(),
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
        dimensionedScalar("zero", dimTime, 0.0),
        "zeroGradient"
    ),
    calculateActivationTime_
    (
        resolveMyocardiumMesh(supportMesh_, meshSubsetPtr_).nCells(),
        true
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
    )
{
    if (meshSubsetPtr_.valid() && meshSubsetPtr_->hasSubMesh())
    {
        Info<< "Constructed MyocardiumDomain on submesh '"
            << mesh().name() << "' from cellZone '"
            << electroProperties_.lookupOrDefault<word>("cellZone", word::null)
            << "'." << nl << endl;
    }

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
    updateExternalStimulusCurrent(sourceField_, externalStimulus_, t0);

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


void MyocardiumDomain::write()
{
    if (meshSubsetPtr_.valid() && meshSubsetPtr_->hasSubMesh())
    {
        const labelUList& cellMap = meshSubsetPtr_->cellMap();

        writeMappedCellField(Vm_, supportMesh_, cellMap);
        writeMappedCellField(sourceField_, supportMesh_, cellMap);
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
