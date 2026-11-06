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
    bidomainSolver

Description
    Bidomain FVM reaction-diffusion kernel.

Author
    Simao Nieto de Castro. All rights reserved.
\*---------------------------------------------------------------------------*/

#ifndef bidomainSolver_H
#define bidomainSolver_H

#include "myocardiumSolver.H"
#include "conductivityFieldIO.H"
#include "overrideTypeName.H"
#include "Enum.H"
#include "Tuple2.H"
#include "fvc.H"
#include "fvm.H"

namespace Foam
{

class bidomainSolver
:
    public myocardiumSolver
{
public:

    //- Heart phiE at the heart-torso wall with a bath (bathHeartPhiETrace)
    enum class heartPhiETrace
    {
        //- zeroGradient heart phiE companion
        zeroGradient,

        //- The global phiE on the exposed heart faces (sealed heart only)
        global
    };

    static const Enum<heartPhiETrace> heartPhiETraceNames;

    //- bathHeartPhiETrace of the bidomainSolver Coeffs dict (required):
    //  the single reader of the key. global requires sealedHeartBoundary
    //  true.
    static heartPhiETrace readBathHeartPhiETrace(const dictionary& dict);


private:

    //- The heart phiE this solver owns and solves without a bath; null with
    //  a bath (torsoCoupling bath), whose domain owns the heart phiE
    autoPtr<volScalarField> ownPhiEPtr_;

    //- The heart phiE the Vm equation reads: the owned one, or the bath's
    //  heart phiE once bound by bindExternalPhiE. This solver never writes
    //  a bound phiE.
    const volScalarField* phiEPtr_;

    volScalarField phiI_;
    volTensorField Gi_;
    volTensorField Ge_;
    volTensorField GiPlusGe_;
    Switch sealedHeartBoundary_;
    //- phiERefPoint and the value pinned at its cell; null when the bath
    //  owns phiE (torsoCoupling bath)
    autoPtr<Tuple2<point, scalar>> phiEReferencePtr_;
    const electroVerificationModel* verifierPtr_;

    tmp<volTensorField> initialiseConductivityTensor
    (
        const fvMesh& mesh,
        const fvMesh& supportMesh,
        const fvMeshSubset* meshSubsetPtr,
        const conductivityFieldSpec& spec,
        const dictionary& dict
    ) const;

    //- The heart phiE trace at the heart-torso wall: readBathHeartPhiETrace
    //  with a torso (torsoCoupling bath or leadField). Without a torso the
    //  key is rejected and the heart wall is the outer, zeroGradient
    //  boundary.
    static heartPhiETrace readHeartPhiETrace(const dictionary& dict);

    label referenceCell() const;

    //- The heart phiE the Vm equation reads; fatal with a bath before the
    //  bath has bound its heart phiE
    const volScalarField& phiE() const;

    //- Solve the owned heart phiE; nothing with a bath, whose domain
    //  solves the global phiE and updates the bound heart phiE
    void solvePhiEImplicitOnce(electroVolumeFieldDomain& domain);

    void solveVmImplicitOnce(electroVolumeFieldDomain& domain);

public:

    OverrideTypeName("bidomainSolver");

    bidomainSolver
    (
        const fvMesh& mesh,
        const fvMesh& supportMesh,
        const fvMeshSubset* meshSubsetPtr,
        const dictionary& electroProperties
    );

    //- The heart phiE: owned, or bound from the bath; null with a bath
    //  before binding
    virtual const volScalarField* phiEPtr() const override
    {
        return phiEPtr_;
    }

    virtual bool hasExtracellularPotential() const override
    {
        return true;
    }

    //- The owned phiE (no bath), or null with a bath
    virtual volScalarField* ownedPhiEPtr() override
    {
        return ownPhiEPtr_.get();
    }

    virtual const volScalarField* phiIPtr() const override
    {
        return &phiI_;
    }

    virtual bool sealedHeartBoundary() const override
    {
        return sealedHeartBoundary_;
    }

    virtual const volTensorField* conductivityPtr() const override
    {
        return &Gi_;
    }

    virtual const volTensorField* intracellularConductivityPtr() const override
    {
        return &Gi_;
    }

    virtual const volTensorField* extracellularConductivityPtr() const override
    {
        return &Ge_;
    }

    virtual const volTensorField& explicitConductivityField() const override
    {
        return Gi_;
    }

    //- Read the bath's heart phiE (torsoCoupling bath) in the Vm
    //  equation; fatal for a solver that owns its phiE
    virtual void bindExternalPhiE(const volScalarField& heartPhiE) override;

    virtual void unbindExternalPhiE() override;

    virtual void bindVerifier(const electroVerificationModel* verifier) override
    {
        verifierPtr_ = verifier;
    }

    virtual void solveDiffusionExplicit
    (
        electroVolumeFieldDomain& domain,
        scalar dt
    ) override;

    virtual void solveDiffusionImplicit
    (
        electroVolumeFieldDomain& domain,
        scalar dt
    ) override;
};

} // End namespace Foam

#endif

// ************************************************************************* //
