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

/*---------------------------------------------------------------------------*\
  NOT_COMPILED
  This file is not included in Make/files and is not assembled by any active
  orchestration path. It is a stub for a future bath-coupling workflow.
  See electroCouplers/heartBathInterfaceCoupler for the matching coupler.
\*---------------------------------------------------------------------------*/

#ifndef BathDomain_H
#define BathDomain_H

#include "autoPtr.H"
#include "fvMesh.H"
#include "electroDomainCouplingEndpoints.H"
#include "electroDomainInterface.H"
#include "electroStateProvider.H"
#include "bathECGSolver.H"
#include "fvMeshSubset.H"
#include "volFields.H"
#include "dimVoltage.H"

namespace Foam
{

class BathDomain
:
    public electroDomainInterface,
    public bathCouplingEndpoint,
    public electroStateProvider
{
    autoPtr<BathECGSolver> solverPtr_;
    const electroStateProvider& stateProvider_;
    autoPtr<fvMesh> ownedMeshPtr_;
    autoPtr<fvMeshSubset> meshSubsetPtr_;
    const fvMesh& supportMesh_;
    const fvMesh& mesh_;
    word interfacePatchName_;
    label interfacePatchIndex_;
    volScalarField phiE_;
    volScalarField interfaceSource_;
    volScalarField sigmaBath_;

    BathDomain(const BathDomain&) = delete;
    void operator=(const BathDomain&) = delete;

public:

    TypeName("bathDomain");

    BathDomain
    (
        const electroStateProvider& stateProvider,
        const fvMesh& mesh,
        const dictionary& dict,
        const word& domainName = "bath"
    );

    virtual ~BathDomain() = default;

    virtual const fvMesh& mesh() const override
    {
        return mesh_;
    }

    virtual const fvMesh& baseMesh() const override
    {
        return supportMesh_;
    }

    virtual const Time& time() const override
    {
        return mesh_.time();
    }

    virtual const volScalarField* phiEPtr() const override
    {
        return &phiE_;
    }

    virtual volScalarField& interfaceSourceField() override
    {
        return interfaceSource_;
    }

    virtual volScalarField& bathPotentialField() override
    {
        return phiE_;
    }

    virtual label interfacePatchIndex() const override
    {
        return interfacePatchIndex_;
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

    volScalarField& phiERef()
    {
        return phiE_;
    }

    const volScalarField& phiE() const
    {
        return phiE_;
    }

    volScalarField& sigmaBathRef()
    {
        return sigmaBath_;
    }

    const volScalarField& sigmaBath() const
    {
        return sigmaBath_;
    }

    const electroStateProvider& stateProvider() const
    {
        return stateProvider_;
    }

    bool usesSubsetMesh() const
    {
        return meshSubsetPtr_.valid() && meshSubsetPtr_->hasSubMesh();
    }

    virtual void advance
    (
        scalar t0,
        scalar dt
    ) override;
};

} // End namespace Foam

#endif

// ************************************************************************* //
