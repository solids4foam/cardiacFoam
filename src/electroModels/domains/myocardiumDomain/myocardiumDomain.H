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

#include "electroModel/myocardialDomainModel.H"
#include "electroModel/myocardiumDomainConfig.H"
#include "ionicModel.H"
#include "electroVerificationModel.H"
#include "dimVoltage.H"
#include "stimulusIO.H"
#include "reactionDiffusionSolver.H"

namespace Foam
{

class MyocardiumDomain
:
    public myocardialDomainModel
{
protected:

    autoPtr<reactionDiffusionSolver> diffusionSolverPtr_;
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
    MonodomainStimulusProtocol monodomainStimulus_;
    Switch useExplicitAlgorithm_;

    void updateExternalStimulusCurrent
    (
        volScalarField& externalStimulusCurrent,
        const MonodomainStimulusProtocol& externalStimulus,
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

    MyocardiumDomain
    (
        const fvMesh& mesh,
        const dictionary& electroProperties,
        myocardiumDomainConfig& config,
        autoPtr<reactionDiffusionSolver> diffusionSolverPtr,
        volScalarField* phiEPtr = nullptr
    );

    virtual ~MyocardiumDomain() = default;

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

    virtual volScalarField& sourceField() override
    {
        return sourceField_;
    }

    virtual const ElectromechanicalSignalProvider* provider() const override
    {
        return &ionicModel_;
    }

    volScalarField& VmRef()
    {
        return Vm_;
    }

    const volScalarField& Iion() const
    {
        return Iion_;
    }

    const dimensionedScalar& chi() const
    {
        return chi_;
    }

    const dimensionedScalar& Cm() const
    {
        return Cm_;
    }

    virtual void prepareTimeStep(scalar t0, scalar dt) override {}

    virtual void advance(scalar t0, scalar dt) override;

    virtual void advance
    (
        scalar t0,
        scalar dt,
        pimpleControl* pimplePtr
    ) override;

    virtual scalar suggestExplicitDeltaT(scalar maxCo) const override;

    virtual bool shouldPostProcess() const override;

    virtual void exportStates() override;

    virtual void exportPostProcessFields() override;

    virtual void postProcess() override;
};

} // End namespace Foam

#endif

// ************************************************************************* //
