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

#ifndef eikonalMyocardiumDomain_H
#define eikonalMyocardiumDomain_H

#include "myocardiumDomainInterface.H"
#include "dimVoltage.H"
#include "fvMeshSubset.H"
#include "fvc.H"
#include "fvm.H"

namespace Foam
{

class EikonalMyocardiumDomain
:
    public myocardiumDomainInterface
{
    autoPtr<fvMeshSubset> meshSubsetPtr_;
    const fvMesh& supportMesh_;
    volScalarField activationTime_;
    volScalarField Vm_;
    volScalarField sourceField_;
    volVectorField gradActivationTime_;
    labelList stimulusCellIDs_;
    const dictionary& electroProperties_;
    const dimensionedScalar chi_;
    const dimensionedScalar Cm_;
    volTensorField conductivity_;
    const volTensorField M_;
    volVectorField w_;
    volScalarField G_;
    const dimensionedScalar c0_;
    volVectorField a_;
    volVectorField u_;
    surfaceScalarField phiU_;
    volScalarField divPhiU_;
    const Switch eikonalAdvectionDiffusionApproach_;

    tmp<volTensorField> initialiseConductivity() const;

    EikonalMyocardiumDomain(const EikonalMyocardiumDomain&) = delete;
    void operator=(const EikonalMyocardiumDomain&) = delete;

public:

    using myocardiumDomainInterface::advance;

    EikonalMyocardiumDomain
    (
        const fvMesh& supportMesh,
        const dictionary& electroProperties
    );

    virtual ~EikonalMyocardiumDomain() = default;

    virtual const fvMesh& mesh() const override
    {
        return activationTime_.mesh();
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
        return mesh().time();
    }

    virtual const volScalarField& Vm() const override
    {
        return Vm_;
    }

    virtual const volScalarField& activationTime() const override
    {
        return activationTime_;
    }

    virtual volScalarField& activationTimeField() override
    {
        return activationTime_;
    }

    virtual volScalarField& sourceField() override
    {
        return sourceField_;
    }

    virtual const volScalarField* VmPtr() const override
    {
        return &Vm_;
    }

    virtual const volTensorField* conductivityPtr() const override
    {
        return &conductivity_;
    }

    virtual const volTensorField* intracellularConductivityPtr() const override
    {
        return &conductivity_;
    }

    virtual const volTensorField* extracellularConductivityPtr() const override
    {
        return &conductivity_;
    }

    virtual void advance(scalar t0, scalar dt) override;

    virtual void advance
    (
        scalar t0,
        scalar dt,
        pimpleControl* pimplePtr
    ) override;

    virtual void write() override;

    virtual scalar suggestExplicitDeltaT(scalar maxCo) const override;

    virtual bool shouldPostProcess() const override
    {
        return false;
    }

    virtual void exportStates() override
    {}

    virtual void postProcess() override
    {}

    virtual bool applyModelTimeControls(Time& runTime) const override;
};

} // End namespace Foam

#endif

// ************************************************************************* //
