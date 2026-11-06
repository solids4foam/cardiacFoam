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

Class
    extracellularPotentialDomain

Description
    Owns the unified extracellular potential field phiE on the base
    mesh covering both heart and bath cellZones.

Author
    Simao Nieto de Castro. All rights reserved.
\*---------------------------------------------------------------------------*/

#ifndef extracellularPotentialDomain_H
#define extracellularPotentialDomain_H

#include "volFields.H"
#include "surfaceFields.H"
#include "fvMatricesFwd.H"
#include "fvMesh.H"
#include "fvMeshSubset.H"
#include "dictionary.H"
#include "autoPtr.H"
#include "Switch.H"
#include "electroStateDomain.H"
#include "myocardiumDomainInterface.H"
#include "passiveConductor.H"
#include "overrideTypeName.H"

namespace Foam
{

class extracellularPotentialDomain
:
    public electroStateDomain
{
    const fvMesh& baseMesh_;

    // Heart-side data sources (non-owning).
    myocardiumDomainInterface& heartDomain_;

    // Heart Vm, read from the myocardium
    const volScalarField& heartVm_;

    //- Passive heart+bath operator built from the top-level torso block,
    //  with its conductivities, boundary data, reference and heart -> base
    //  map
    autoPtr<passiveConductor> conductorPtr_;

    //- The heart phiE (heart mesh, named phiE): the restriction of the
    //  global phiE to the heart cells, with the exposed-face trace when
    //  bathHeartPhiETrace is global. Owned and written here only; the
    //  myocardium solver reads it in its Vm equation, and the operator uses
    //  it as its heart companion. Its patch types and wall states are set
    //  once, at construction.
    autoPtr<volScalarField> heartPhiEPtr_;

    // Global potential and the V_m field projected to the base mesh.
    autoPtr<volScalarField> phiEPtr_;
    autoPtr<volScalarField> VmGlobalPtr_;

    //- Shared number of additional corrected-Laplacian assemblies, read from
    //  PIMPLE/nNonOrthogonalCorrectors (required) and used for phiE and Vm.
    //  Zero means one equation solve.
    label nNonOrthogonalCorrectors_;

    //- Time at which manufactured data are evaluated
    scalar manufacturedEvaluationTime_;

    //- Face flux of the phiE equation on heart-bath faces (A, along Sf, extracellular part)
    autoPtr<surfaceScalarField> interfaceFluxPtr_;

    //- Cell gradient of phiE used by the face corrections of the phiE operator (jump-aware under conormalHarmonic)
    autoPtr<volVectorField> interfaceGradientPtr_;

    //- Write the heart phiE from the global phiE: heart cells, and the
    //  exposed-face trace under bathHeartPhiETrace global
    void updateHeartPhiE();

    //- Assemble and solve the global phiE equation once.
    void solvePhiEOnce();

    //- Verification model bound to the heart domain, or nullptr.
    const electroVerificationModel* verifier() const;

    //- Copy the heart-bath face values of F into phiEInterfaceFlux.
    void storeInterfaceFlux(const surfaceScalarField& F);

public:

    OverrideTypeName("extracellularPotentialDomain");

    //- Construct on the base mesh from the bidomain heart domain and the
    //  top-level electroProperties: its torso block and the
    //  sealedHeartBoundary of the myocardium Coeffs
    extracellularPotentialDomain
    (
        const fvMesh& baseMesh,
        myocardiumDomainInterface& heartDomain,
        const dictionary& electroProperties
    );

    virtual ~extracellularPotentialDomain() = default;

    // electroDomainInterface
    virtual const Time& time() const override { return baseMesh_.time(); }
    virtual void advance(scalar t0, scalar dt) override;
    virtual void prepareTimeStep(scalar t0, scalar dt) override;
    virtual void write() override;

    // electroStateProvider — exposes the global phiE.
    virtual const fvMesh& mesh() const override { return baseMesh_; }
    virtual const fvMesh& baseMesh() const override { return baseMesh_; }
    virtual const volScalarField* phiEPtr() const override
    {
        return phiEPtr_.valid() ? &phiEPtr_() : nullptr;
    }

    // Direct accessors for downstream consumers (heart solver and ECG probes).
    const volScalarField& phiE() const { return phiEPtr_(); }
    const volTensorField& sigmaTotal() const
    {
        return conductorPtr_->sigmaTotal();
    }

    virtual const volTensorField* conductivityPtr() const override
    {
        return conductorPtr_.valid() ? &conductorPtr_->sigmaTotal() : nullptr;
    }

    virtual const volTensorField* intracellularConductivityPtr()
        const override
    {
        return heartDomain_.intracellularConductivityPtr();
    }

    virtual const volTensorField* extracellularConductivityPtr()
        const override
    {
        return heartDomain_.extracellularConductivityPtr();
    }

    // Scatter heart V_m -> VmGlobalPtr_ (zero in bath cells).
    void scatterHeartVm();

    // Solve one timestep of the global phiE equation.
    void solvePhiE();
};

} // End namespace Foam

#endif

// ************************************************************************* //
