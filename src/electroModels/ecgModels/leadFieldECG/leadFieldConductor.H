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
    leadFieldConductor

Description
    The heart+torso passiveConductor of a torsoECG domain with torsoCoupling
    leadField, built from electroProperties, with its electrodes and
    signature.

    - Coupling: torsoCoupling must be leadField.
    - Heart: the myocardium cellZone sub-mesh (createMyocardiumMeshSubset).
    - Myocardium: bidomainSolver or monodomainSolver; eikonalSolver is not
      supported.
    - Gi, Ge: bidomainSolver: those of the running myocardium, or read from
      its Coeffs when built from the case alone; monodomainSolver reads them
      from the ECG domain entry with the same keys. For monodomainSolver
      the Coeffs conductivity must equal Gi (Gi+Ge)^-1 Ge in every cell.
    - Torso: the top-level torso block of electroProperties
      (passiveConductor::torsoSettings).
    - Sealing: sealedHeartBoundary of the Coeffs, with a heart field Vm
      that is zeroGradient on every non-coupled patch.
    - Electrodes: those of the ECG domain (readEcgElectrodes); each lies in
      a torso cell of the mesh (ecgElectrodeCells); an electrode outside the
      mesh or in a myocardium cell is fatal.
    - bidomainSolver: bathHeartPhiETrace must be zeroGradient, the heart
      phiE companion the conductor assembles.

    Built from the case alone, the heart sub-mesh and the sealing Vm are
    its own. Built from a running myocardium, they are the myocardium's
    sub-mesh and Vm.

    The conductor's base-mesh fields and the Gi and Ge it reads are not
    registered.

SourceFiles
    leadFieldConductor.C

Author
    Simao Nieto de Castro. All rights reserved.
\*---------------------------------------------------------------------------*/

#ifndef leadFieldConductor_H
#define leadFieldConductor_H

#include "passiveConductor.H"
#include "bidomainSolver.H"
#include "ecgElectrodes.H"
#include "electroStateProvider.H"
#include "fvMeshSubset.H"
#include "dictionary.H"

namespace Foam
{

class leadFieldConductor
{
    const fvMesh& baseMesh_;

    // Myocardium solver type (Coeffs name without "Coeffs").
    const word solverType_;

    // The myocardium is a bidomainSolver (else a monodomainSolver).
    const bool bidomain_;

    // bidomainSolver heart phiE trace at the heart-torso wall.
    bidomainSolver::heartPhiETrace bathHeartPhiETrace_;

    Switch sealedHeartBoundary_;

    // Myocardium cellZone of the base mesh.
    word myocardiumZone_;

    // Electrodes of the ECG domain.
    const ecgElectrodes electrodes_;

    // Heart sub-mesh: owned, or the running myocardium's.
    autoPtr<fvMeshSubset> ownSubsetPtr_;
    const fvMeshSubset* subsetPtr_;

    // Heart conductivities on the heart sub-mesh: read and owned, or the
    // running bidomain myocardium's.
    autoPtr<volTensorField> ownGiPtr_;
    autoPtr<volTensorField> ownGePtr_;
    const volTensorField* GiPtr_;
    const volTensorField* GePtr_;

    // Heart field whose patch types decide sealing: owned, or the
    // running myocardium's Vm.
    autoPtr<volScalarField> ownSealingFieldPtr_;
    const volScalarField* sealingFieldPtr_;

    // Base cell of each electrode on its owning processor, -1 elsewhere.
    labelList electrodeCells_;

    autoPtr<passiveConductor> conductorPtr_;

    // Members from electroProperties, the ECG domain entry name and the
    // electrodes; checks the coupling, the myocardium solver type and its
    // cellZone. The heart, conductivities and conductor are set by the
    // public constructors.
    leadFieldConductor
    (
        const fvMesh& baseMesh,
        const dictionary& electroProperties,
        const word& ecgDomainName,
        const ecgElectrodes& electrodes
    );

    // Check the coupling, the myocardium solver type and its cellZone.
    void checkMyocardium
    (
        const dictionary& electroProperties,
        const word& ecgDomainName
    ) const;

    // Set Gi and Ge (from heart, a running bidomain myocardium, when given)
    // and build the conductor and the electrode cells.
    void build
    (
        const dictionary& electroProperties,
        const word& ecgDomainName,
        const electroStateProvider* heart
    );

    // Set Gi and Ge; check the monodomain conductivity against them.
    void setHeartConductivities
    (
        const dictionary& electroCoeffs,
        const dictionary& ecgDict,
        const electroStateProvider* heart
    );

    // Locate the electrode cells; fatal for an electrode in a myocardium
    // cell.
    void locateElectrodes(const dictionary& ecgDict);

    // Zero potential named phiE on mesh (unregistered) with the patch
    // types, whose wall-gradient states (seal) are published unless the
    // forward's potential of the same name on mesh has published them
    tmp<volScalarField> newPotential
    (
        const fvMesh& mesh,
        const dimensionSet& dims,
        const wordList& patchTypes,
        const bool seal
    ) const;


public:

    //- The fixed point of the lead fields (buildLeadFields) and of the
    //  verifier's direct forward solve: the fvSolution dictionary
    //  leadFieldFixedPoint
    struct fixedPointControls
    {
        //- Relative change of the unknown between passes below which the
        //  fixed point stops (tolerance, required)
        scalar tolerance;

        //- Passes after which the fixed point is fatal (maxPasses,
        //  required)
        label maxPasses;

        //- Passes before the change can stop the fixed point: the change
        //  then compares two solves of the fixed point, not a solve with
        //  its start value
        static constexpr label minPasses = 2;

        //- Read leadFieldFixedPoint from the fvSolution of mesh
        static fixedPointControls read(const fvMesh& mesh);
    };


    //- Construct from the base mesh, the top-level electroProperties and
    //  the name of the ECG domain entry in its ecgDomains, whose electrodes
    //  are read (readEcgElectrodes)
    leadFieldConductor
    (
        const fvMesh& baseMesh,
        const dictionary& electroProperties,
        const word& ecgDomainName
    );

    //- Construct on a running myocardium (its sub-mesh, its Vm, which
    //  decides sealing, and the Gi and Ge of a bidomain myocardium) from
    //  the top-level electroProperties, the ECG domain entry name and the
    //  electrodes of the ECG domain
    leadFieldConductor
    (
        const electroStateProvider& heart,
        const dictionary& electroProperties,
        const word& ecgDomainName,
        const wordList& electrodeNames,
        const List<vector>& electrodePositions
    );

    leadFieldConductor(const leadFieldConductor&) = delete;
    void operator=(const leadFieldConductor&) = delete;

    ~leadFieldConductor() = default;


    const passiveConductor& conductor() const
    {
        return conductorPtr_();
    }

    const fvMesh& baseMesh() const
    {
        return baseMesh_;
    }

    const fvMesh& heartMesh() const
    {
        return subsetPtr_->subMesh();
    }

    const volTensorField& Gi() const
    {
        return *GiPtr_;
    }

    const volTensorField& Ge() const
    {
        return *GePtr_;
    }

    const wordList& electrodeNames() const
    {
        return electrodes_.names;
    }

    //- Base (torso) cell of each electrode on its owning processor, -1
    //  elsewhere (ecgElectrodeCells)
    const labelList& electrodeCells() const
    {
        return electrodeCells_;
    }

    //- Decomposition-independent signature: the operator (adjoint), global
    //  cell counts, the global sum of base-mesh cell volumes and of volume
    //  times cell centre, the global cell count of the myocardium and of
    //  each torso cellZone, sums of Gi, Ge and the torso conductivity,
    //  myocardium solver, interface option (and nInterfaceGradientCorrectors
    //  with conormalHarmonic), seal flag, sealing patch types,
    //  bidomain heart trace, geometry scheme, the base and heart fvSchemes
    //  entries that can resolve the schemes of phiE (default, patterns and
    //  keys containing phiE), electrodes, ground and surface-current
    //  patches and the reference
    dictionary signature() const;

    //- True if built, a signature written by buildLeadFields, equals
    //  signature(): words and labels exactly, other numbers to a relative
    //  1e-10 of the largest number of the entry. The differing entries are
    //  written to os.
    bool matches(const dictionary& built, Ostream& os) const;

    //- Zero potential phiE on the base mesh (unregistered) with the
    //  forward's patch types, ground values and surface-current gradients
    tmp<volScalarField> forwardPotential() const;

    //- Zero lead-field potential phiE on the base mesh (unregistered):
    //  fixedValue 0 on ground patches, zeroGradient on the other
    //  non-coupled patches
    tmp<volScalarField> leadFieldPotential() const;

    //- Zero heart companion phiE on the heart mesh (unregistered,
    //  zeroGradient), sealed as the forward's
    tmp<volScalarField> heartPotential() const;

    //- Constant right-hand side of the forward on the base mesh: the
    //  operator assembled on a zero potential that carries the forward's
    //  ground values, surface currents and reference value
    tmp<scalarField> constantSource() const;
};

} // End namespace Foam

#endif

// ************************************************************************* //
