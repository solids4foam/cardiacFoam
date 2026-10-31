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
    The heart+torso passiveConductor of a lead-field ECG domain, built from
    the case dictionaries, with its electrodes, reference and signature.

    - Heart: the myocardium cellZone sub-mesh (createMyocardiumMeshSubset).
    - Myocardium: bidomainSolver or monodomainSolver; eikonalSolver is not
      supported.
    - Gi, Ge: bidomainSolver reads them from its Coeffs; monodomainSolver
      reads them from the ECG domain entry with the same keys.
      For monodomainSolver the Coeffs conductivity must equal
      Gi (Gi+Ge)^-1 Ge in every cell.
    - Torso: the bathPotentialDomain block of the Coeffs, or the torso
      sub-dictionary of the ECG domain entry; exactly one of them.
    - Sealing: sealedHeartBoundary of the Coeffs, with a heart field Vm
      that is zeroGradient on every non-coupled patch.
    - Electrodes: each lies in a torso cell; an electrode whose cell is a
      myocardium cell is fatal.
    - bidomainSolver: bathHeartPhiETrace must be zeroGradient, the heart
      phiE companion the conductor assembles.

    Built from the case alone, the heart sub-mesh and the sealing Vm are
    its own. Built from a running myocardium, they are the myocardium's
    sub-mesh and Vm.

    The conductor's base-mesh fields and its Gi and Ge are not registered.

SourceFiles
    leadFieldConductor.C

Author
    Simao Nieto de Castro. All rights reserved.
\*---------------------------------------------------------------------------*/

#ifndef leadFieldConductor_H
#define leadFieldConductor_H

#include "passiveConductor.H"
#include "fvMeshSubset.H"
#include "dictionary.H"
#include "scalarList.H"

namespace Foam
{

class leadFieldConductor
{
    const fvMesh& baseMesh_;

    // Myocardium solver type (Coeffs name without "Coeffs").
    const word solverType_;

    // bidomainSolver heart phiE trace at the heart-torso wall.
    word bathHeartPhiETrace_;

    // Myocardium cellZone of the base mesh.
    word myocardiumZone_;

    // Heart sub-mesh: owned, or the running myocardium's.
    autoPtr<fvMeshSubset> ownSubsetPtr_;
    const fvMeshSubset* subsetPtr_;

    // Heart conductivities on the heart sub-mesh.
    autoPtr<volTensorField> GiPtr_;
    autoPtr<volTensorField> GePtr_;

    // Torso conductivity on the base mesh and its cellZones.
    autoPtr<volScalarField> torsoConductivityPtr_;
    wordList torsoCellZones_;

    word interfaceConductivityInterpolation_;
    Switch sealedHeartBoundary_;

    // Heart field whose patch types decide sealing: owned, or the
    // running myocardium's Vm.
    autoPtr<volScalarField> ownSealingFieldPtr_;
    const volScalarField* sealingFieldPtr_;

    wordList electrodeNames_;
    List<vector> electrodePositions_;

    // Base cell of each electrode on its owning processor, -1 elsewhere.
    labelList electrodeCells_;

    wordList groundPatches_;
    scalarList groundValues_;
    wordList surfaceCurrentPatches_;
    scalarList surfaceCurrentValues_;

    // Reference pin, used when there is no ground patch.
    bool hasReferencePoint_;
    point referencePoint_;
    scalar referenceValue_;
    label referenceCell_;

    autoPtr<passiveConductor> conductorPtr_;

    // Check the myocardium solver type and its cellZone.
    void checkMyocardium
    (
        const dictionary& electroCoeffs,
        const word& ecgDomainName
    ) const;

    // Read conductivities, torso and electrodes; build the conductor.
    void build(const dictionary& electroCoeffs, const word& ecgDomainName);

    // Read Gi and Ge; check the monodomain conductivity against them.
    void readHeartConductivities
    (
        const dictionary& electroCoeffs,
        const dictionary& ecgDict
    );

    // Read the torso, interface and reference keys.
    void readTorso(const dictionary& torsoDict, const bool bathBlock);

    void readElectrodes
    (
        const dictionary& electroCoeffs,
        const dictionary& ecgDict
    );

    // Locate the electrode cells; fatal for an electrode in a myocardium
    // cell.
    void locateElectrodes(const dictionary& ecgDict);


public:

    //- Construct from the base mesh, the myocardium solver Coeffs and the
    //  name of the ECG domain entry in its ecgDomains
    leadFieldConductor
    (
        const fvMesh& baseMesh,
        const dictionary& electroCoeffs,
        const word& ecgDomainName
    );

    //- Construct on the sub-mesh of a running myocardium, whose Vm decides
    //  sealing, from its solver Coeffs and the ECG domain entry name
    leadFieldConductor
    (
        const fvMeshSubset& heartSubset,
        const volScalarField& Vm,
        const dictionary& electroCoeffs,
        const word& ecgDomainName
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
        return GiPtr_();
    }

    const volTensorField& Ge() const
    {
        return GePtr_();
    }

    //- Heart field Vm whose patch types decide sealing
    const volScalarField& sealingField() const
    {
        return *sealingFieldPtr_;
    }

    const List<vector>& electrodePositions() const
    {
        return electrodePositions_;
    }

    const wordList& electrodeNames() const
    {
        return electrodeNames_;
    }

    //- Base (torso) cell of each electrode on its owning processor, -1
    //  elsewhere (ecgElectrodeCells)
    const labelList& electrodeCells() const
    {
        return electrodeCells_;
    }

    //- Pinned base cell on its owning processor; -1 elsewhere and when
    //  ground patches fix the constant
    label referenceCell() const
    {
        return referenceCell_;
    }

    scalar referenceValue() const
    {
        return referenceValue_;
    }

    const wordList& groundPatches() const
    {
        return groundPatches_;
    }

    //- Value of each ground patch, in groundPatches() order
    const scalarList& groundValues() const
    {
        return groundValues_;
    }

    const wordList& surfaceCurrentPatches() const
    {
        return surfaceCurrentPatches_;
    }

    //- Current of each surface-current patch, in surfaceCurrentPatches()
    //  order
    const scalarList& surfaceCurrentValues() const
    {
        return surfaceCurrentValues_;
    }

    //- Decomposition-independent signature: the operator (adjoint), global
    //  cell counts, the global sum of base-mesh cell volumes and of volume
    //  times cell centre, the global cell count of the myocardium and of
    //  each torso cellZone, sums of Gi, Ge and the torso conductivity,
    //  myocardium solver, interface option, seal flag, sealing patch types,
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
