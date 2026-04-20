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

#ifndef MyocardiumDomain_H
#define MyocardiumDomain_H

#include "myocardiumDomainInterface.H"
#include "electromechanicalSignalProvider.H"
#include "ionicModel.H"
#include "electroVerificationModel.H"
#include "dimVoltage.H"
#include "electroVolumeFieldDomain.H"
#include "fvMeshSubset.H"
#include "myocardiumSolver.H"
#include "stimulusIO.H"
#include "Switch.H"

namespace Foam
{

class MyocardiumDomain
:
    public myocardiumDomainInterface,
    public electroVolumeFieldDomain
{
protected:

    autoPtr<fvMeshSubset> meshSubsetPtr_;
    const fvMesh& supportMesh_;
    autoPtr<myocardiumSolver> diffusionSolverPtr_;
    volScalarField Vm_;
    volScalarField sourceField_;
    volScalarField Iion_;
    volScalarField activationTime_;
    boolList calculateActivationTime_;
    PtrList<volScalarField>& outFields_;
    wordList preProcessFieldNames_;
    PtrList<volScalarField> preProcessFields_;
    wordList postProcessFieldNames_;
    PtrList<volScalarField>& postProcessFields_;
    ionicModel& ionicModel_;
    electroVerificationModel* verificationModelPtr_;
    const dictionary& electroProperties_;
    const dimensionedScalar chi_;
    const dimensionedScalar Cm_;
    ExternalStimulusProtocol externalStimulus_;
    Switch useExplicitAlgorithm_;
    Switch reportSetup_;

    void updateExternalStimulusCurrent
    (
        volScalarField& externalStimulusCurrent,
        const ExternalStimulusProtocol& externalStimulus,
        scalar t0
    ) const;

    void updateActivationTime
    (
        volScalarField& activationTime,
        boolList& calculateActivationTime,
        const volScalarField& Vm
    ) const;

    void validateNoIonicStimulusInMonodomain() const;

    void bindBidomainField(volScalarField& phiE);

    void initialiseProcessing();

public:

    using myocardiumDomainInterface::advance;

    //- Factory: derives the solver type from the Coeffs subdict name and
    //  delegates to the injecting constructor below.
    //  ("monodomainSolverCoeffs" → selects "monodomainSolver", etc.)
    static autoPtr<MyocardiumDomain> New
    (
        const fvMesh& mesh,
        const dictionary& electroProperties,
        PtrList<volScalarField>& outFields,
        const wordList& postProcessFieldNames,
        PtrList<volScalarField>& postProcessFields,
        ionicModel& ionicModel,
        electroVerificationModel* verificationModelPtr
    );

    static label configuredCellCount
    (
        const fvMesh& mesh,
        const dictionary& electroProperties
    );

    MyocardiumDomain
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
    );

    virtual ~MyocardiumDomain() = default;

    virtual const fvMesh& mesh() const override
    {
        return Vm_.mesh();
    }

    virtual const fvMesh& baseMesh() const override
    {
        return supportMesh_;
    }

    virtual const labelUList* subsetFaceMapPtr() const override
    {
        return
        (
            meshSubsetPtr_.valid() && meshSubsetPtr_->hasSubMesh()
          ? &meshSubsetPtr_->faceMap()
          : nullptr
        );
    }

    virtual const labelUList* subsetCellMapPtr() const override
    {
        return
        (
            meshSubsetPtr_.valid() && meshSubsetPtr_->hasSubMesh()
          ? &meshSubsetPtr_->cellMap()
          : nullptr
        );
    }

    virtual const Time& time() const override
    {
        return Vm_.mesh().time();
    }

    virtual const volScalarField& Vm() const override
    {
        return Vm_;
    }

    virtual const volScalarField* VmPtr() const override
    {
        return &Vm_;
    }

    virtual const volScalarField* phiEPtr() const override
    {
        return diffusionSolverPtr_->phiEPtr();
    }

    virtual const volTensorField* conductivityPtr() const override
    {
        return diffusionSolverPtr_->conductivityPtr();
    }

    virtual const volTensorField* intracellularConductivityPtr() const override
    {
        return diffusionSolverPtr_->intracellularConductivityPtr();
    }

    virtual const volTensorField* extracellularConductivityPtr() const override
    {
        return diffusionSolverPtr_->extracellularConductivityPtr();
    }

    virtual volScalarField& sourceField() override
    {
        return sourceField_;
    }

    virtual const volScalarField& activationTime() const override
    {
        return activationTime_;
    }

    virtual volScalarField& activationTimeField() override
    {
        return activationTime_;
    }

    virtual const ElectromechanicalSignalProvider* provider() const override
    {
        return &ionicModel_;
    }

    virtual volScalarField& VmRef() override
    {
        return Vm_;
    }

    virtual const volScalarField& Iion() const override
    {
        return Iion_;
    }

    virtual const dimensionedScalar& chi() const override
    {
        return chi_;
    }

    virtual const dimensionedScalar& Cm() const override
    {
        return Cm_;
    }

    virtual void prepareTimeStep(scalar t0, scalar dt) override;

    virtual void advance(scalar t0, scalar dt) override;

    virtual void advance
    (
        scalar t0,
        scalar dt,
        pimpleControl* pimplePtr
    ) override;

    virtual void write() override;

    virtual scalar suggestExplicitDeltaT(scalar maxCo) const override;

    virtual bool shouldPostProcess() const override;

    virtual void exportStates() override;

    virtual void exportPostProcessFields();

    virtual void postProcess() override;
};

} // End namespace Foam

#endif

// ************************************************************************* //
