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
#include "fvMesh.H"
#include "fvMeshSubset.H"
#include "dictionary.H"
#include "autoPtr.H"
#include "Switch.H"
#include "electroStateDomain.H"
#include "myocardiumDomainInterface.H"
#include "overrideTypeName.H"

namespace Foam
{

class bidomainSolver;

class extracellularPotentialDomain
:
    public electroStateDomain
{
    const fvMesh& baseMesh_;

    // Heart-side data sources (non-owning).
    myocardiumDomainInterface& heartDomain_;

    // Conductivity coefficients on the BASE mesh.
    autoPtr<volTensorField> sigmaTotalPtr_;     // sigma_i+sigma_e in heart, sigma_b in bath
    autoPtr<volTensorField> sigmaIglobalPtr_;   // sigma_i in heart, 0 in bath

    // Precomputed face coefficient for the global phiE operator.
    autoPtr<surfaceTensorField> sigmaTotalfPtr_;
    autoPtr<surfaceTensorField> sigmaExtracellularfPtr_;

    // Face interpolation selected by bathPotentialDomain. The legacy
    // unweighted harmonic behavior remains the default.
    word interfaceConductivityInterpolation_;
    word intracellularAssembly_;

    // Global potential and the V_m field projected to the base mesh.
    autoPtr<volScalarField> phiEPtr_;
    autoPtr<volScalarField> VmGlobalPtr_;

    // Heart -> base scatter map (built once).
    labelList heartCellToBaseCell_;

    // Reference pinning.
    point phiEReferencePoint_;
    scalar phiEReferenceValue_;
    bool hasPhiEReferencePoint_;

    // CellZone names.

    wordList bathCellZoneNames_;

    word bathConductivityFieldName_;
    wordList surfaceCurrentPatchNames_;
    scalarList surfaceCurrentPatchValues_;

    // True if any phiE patch is Dirichlet (fixedValue); when true, the
    // linear system is already non-singular and we MUST NOT inject an
    // additional setReference pin (it would over-determine the problem
    // with a value that need not match the analytical phiE).
    bool hasDirichletPatch_;

    Switch reportSetup_;

    // Build the heart-cell -> base-cell map from the heart subset.
    void buildHeartScatterMap();

    // Assemble sigmaTotal_ and sigmaIglobal_ from heart conductivities and
    // the bath cellZone scalar conductivities.
    void assembleConductivities();

    // Build sigmaTotalf on internal and coupled faces.
    void buildSigmaTotalSurface();

    // Locate the reference cell on the base mesh (parallel-safe).
    label referenceCell() const;

public:

    OverrideTypeName("extracellularPotentialDomain");

    extracellularPotentialDomain
    (
        const fvMesh& baseMesh,
        myocardiumDomainInterface& heartDomain,
        const dictionary& dict
    );

    virtual ~extracellularPotentialDomain() = default;

    // electroDomainInterface
    virtual const Time& time() const override { return baseMesh_.time(); }
    virtual void prepareTimeStep(scalar t0, scalar dt) override;
    virtual void advance(scalar t0, scalar dt) override;
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
    const volTensorField& sigmaTotal() const { return sigmaTotalPtr_(); }

    virtual const volTensorField* conductivityPtr() const override
    {
        return sigmaTotalPtr_.valid() ? &sigmaTotalPtr_() : nullptr;
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
