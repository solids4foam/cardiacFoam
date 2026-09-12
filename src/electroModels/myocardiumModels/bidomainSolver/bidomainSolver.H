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
#include "fvc.H"
#include "fvm.H"

namespace Foam
{

class bidomainSolver
:
    public myocardiumSolver
{
    volScalarField phiE_;
    // Registered grad(phiE), looked up by name by "corrected" boundary conditions.
    volVectorField gradPhiE_;
    volScalarField phiI_;
    volTensorField Gi_;
    volTensorField Ge_;
    volTensorField GiPlusGe_;
    scalar phiEReferenceValue_;
    point phiEReferencePoint_;
    bool hasPhiEReferencePoint_;
    const volScalarField* externalPhiEBasePtr_;
    const labelUList* externalPhiECellMapPtr_;

    tmp<volTensorField> initialiseConductivityTensor
    (
        const fvMesh& mesh,
        const fvMesh& supportMesh,
        const fvMeshSubset* meshSubsetPtr,
        const conductivityFieldSpec& spec,
        const dictionary& dict
    ) const;

    label referenceCell() const;

    bool externalPhiEBound() const;

    void restrictExternalPhiE();

    //- Refresh gradPhiE_ before a phiE diffusion solve.
    void updateGradPhiE();

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

    virtual const volScalarField* phiEPtr() const override
    {
        return &phiE_;
    }

    virtual const volScalarField* phiIPtr() const override
    {
        return &phiI_;
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

    const volTensorField& intracellularConductivity() const
    {
        return Gi_;
    }

    const volTensorField& extracellularConductivity() const
    {
        return Ge_;
    }

    const volTensorField& totalConductivity() const
    {
        return GiPlusGe_;
    }

    virtual void bindExternalPhiE
    (
        const volScalarField& phiE,
        const labelUList& heartCellMap
    ) override;

    virtual void unbindExternalPhiE() override;

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
