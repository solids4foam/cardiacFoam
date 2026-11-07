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
    passiveConductor

Description
    Passive heart+torso conductor on the base mesh, built from the torso
    block of electroProperties (torsoSettings): the cell and face
    conductivities, the operator for a base-mesh potential and its heart
    companion, the boundary data of the potential (patch types, ground
    values, surface currents), the reference pin and the heart intracellular
    source.

    Its two users are the bath (extracellularPotentialDomain, which solves
    phiE with a bidomain heart every step) and the lead fields of torsoECG
    (leadFieldConductor, buildLeadFields and leadFieldECG).

    The extracellular/torso part is assembled on the base mesh. The
    intracellular part is assembled on the heart sub-mesh and injected
    through the heart-to-base cell and face maps.

    Heart-torso faces carry Ge on the heart side when the intracellular
    potential is insulated at the heart wall (bidomain), and Gi + Ge when the
    source Vm is sealed there (a monodomain source): the current crossing the
    wall is Ge grad(phi).n in the first case and (Gi + Ge) grad(phi).n in the
    second.

    Torso faces between cells of different conductivity take the
    distance-weighted harmonic mean of the two conductivities, which carries
    the exact flux of a potential that is linear on each side; torso faces
    between equal conductivities take the linear one.

    Ground patches fix the constant of the potential; without them the cell
    containing phiERefPoint is pinned to phiEReferenceValue.

    The potential is named phiE on the base and heart meshes. With
    interfaceConductivityInterpolation conormalHarmonic, the solver entry of
    phiE in fvSolution sets the number of jump-aware gradient passes
    (interfaceGradient) of the explicit corrections:

        solvers
        {
            phiE
            {
                ...
                nInterfaceGradientCorrectors 2;
            }
        }

    The conductor publishes the gradient fields whose schemes enter its
    operator and the stencil of its explicit part (leadFieldAdjoint).

Author
    Simao Nieto de Castro. All rights reserved.
\*---------------------------------------------------------------------------*/

#ifndef passiveConductor_H
#define passiveConductor_H

#include "volFields.H"
#include "surfaceFields.H"
#include "fvMatricesFwd.H"
#include "fvMesh.H"
#include "autoPtr.H"
#include "Switch.H"
#include "Enum.H"
#include "dictionary.H"
#include "point.H"
#include "wordList.H"
#include "scalarList.H"
#include "probedOperator.H"

namespace Foam
{

class passiveConductor
{
public:

    //- Face conductivity of the heart-torso faces
    enum class interfaceInterpolation
    {
        //- Distance-weighted harmonic mean; standard cell gradient in the
        //  cross term
        distanceWeightedHarmonic,

        //- Conormal harmonic mean; jump-aware cell gradient in the cross
        //  term (interfaceGradient)
        conormalHarmonic
    };

    static const Enum<interfaceInterpolation> interfaceInterpolationNames;


    //- Explicit face terms of the assembled operator, for flux diagnostics
    struct faceCorrection
    {
        //- Cell gradient of psi used by the face corrections
        //  (jump-aware under conormalHarmonic)
        tmp<volVectorField> gradient;

        //- Explicit face flux added to the operator (conormalHarmonic only)
        tmp<surfaceScalarField> crossFlux;
    };


    //- The keys of the top-level torso block of electroProperties: the
    //  single reader of these keys
    struct torsoSettings
    {
        //- Torso cellZones (cellZones, required)
        wordList cellZones;

        //- Torso scalar conductivity field name (conductivityField,
        //  required)
        word conductivityField;

        //- Face conductivity of the heart-torso faces
        //  (interfaceConductivityInterpolation, required)
        interfaceInterpolation interface;

        //- Dirichlet patches and their values
        wordList groundPatches;
        scalarList groundValues;

        //- Neumann (surface-current) patches and their currents
        wordList surfaceCurrentPatches;
        scalarList surfaceCurrentValues;

        //- phiERefPoint and the value pinned at its cell, both required
        //  without ground patches; not read (zero) with ground patches
        point phiERefPoint;
        scalar phiEReferenceValue;

        //- Read from dict. Every ground and surface-current patch must be
        //  a patch of mesh and no patch may be in both lists; without
        //  ground patches phiERefPoint and phiEReferenceValue are required.
        static torsoSettings read
        (
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- True when no ground patch fixes the constant of the potential,
        //  which is then pinned at phiERefPoint
        bool pinsReference() const
        {
            return groundPatches.empty();
        }

        //- Patch types of the potential on mesh: fixedValue on the ground
        //  patches, fixedGradient on the surface-current patches when
        //  surfaceCurrents (zeroGradient otherwise), zeroGradient on the
        //  other patches
        wordList patchTypes
        (
            const fvMesh& mesh,
            const bool surfaceCurrents
        ) const;

        //- Set the ground values on the ground patches of psi
        void setGroundValues(volScalarField& psi) const;
    };


private:

    const fvMesh& baseMesh_;
    const fvMesh& heartMesh_;

    // Heart fields (non-owning).
    const volTensorField& Gi_;
    const volScalarField& sealingField_;

    // Torso block: zones, conductivity, interface, boundary data and
    // reference.
    const torsoSettings torso_;

    const Switch sealedHeartBoundary_;

    // Jump-aware gradient passes of the corrections (conormalHarmonic),
    // 0 otherwise
    label nInterfaceGradientCorrectors_;

    // Heart-side wall conductivity is Gi + Ge (sealed source), else Ge.
    const bool totalWallConductivity_;

    // Register the base-mesh fields with the base mesh.
    const bool registerFields_;

    // Heart -> base cell and face maps.
    labelList heartCellToBaseCell_;
    labelList heartFaceToBaseFace_;

    // Conductivity coefficients on the BASE mesh.
    autoPtr<volTensorField> sigmaTotalPtr_;     // sigma_i+sigma_e in heart, sigma_b in torso
    autoPtr<volTensorField> sigmaIglobalPtr_;   // sigma_i in heart, 0 in torso

    // Extracellular/torso face coefficient of the base-mesh operator.
    autoPtr<surfaceTensorField> sigmaExtracellularfPtr_;

    // Pinned base cell on its owning processor; -1 elsewhere and when
    // ground patches fix the constant.
    label referenceCell_;

    // Build and check the heart -> base cell and face maps.
    void buildHeartScatterMap
    (
        const labelUList& heartCellMap,
        const labelUList& heartFaceMap
    );

    // Assemble sigmaTotal and sigmaIglobal from heart conductivities and
    // the torso cellZone scalar conductivities.
    void assembleConductivities
    (
        const volTensorField& Ge,
        const volScalarField& torsoConductivity
    );

    // Cell gradient of psi of the conormalHarmonic corrections: the
    // standard gradient refined by the jump-aware passes
    tmp<volVectorField> jumpAwareGradient(const volScalarField& psi) const;

    // Torso operator with conormalHarmonic heart-torso faces: normal
    // Laplacian, stabilisation and cross term, the source, then the
    // non-orthogonal correction of the jump-aware gradient
    tmp<fvScalarMatrix> conormalOperator
    (
        volScalarField& psi,
        const volScalarField* source,
        faceCorrection* correction
    ) const;

    // Torso operator with distanceWeightedHarmonic heart-torso faces:
    // tensor Laplacian and stabilisation, then the source
    tmp<fvScalarMatrix> harmonicOperator
    (
        volScalarField& psi,
        const volScalarField* source,
        faceCorrection* correction
    ) const;

    // Add the heart intracellular block, assembled on heartPsi, to the
    // base-mesh matrix through the heart -> base maps
    void addHeartBlock(fvScalarMatrix& eqn, volScalarField& heartPsi) const;


public:

    //- Name of the potential on the base and heart meshes
    static const word potentialName;

    //- The gradient fields whose schemes enter the operator, on the base
    //  and heart meshes: grad(phiE) (snGrad corrections and the cell
    //  gradient of the corrections) and grad(phiE.component(0)) (the
    //  corrections of the tensor-conductivity Laplacians)
    static wordList gradientNames();

    //- Extracellular/torso face conductivity on all faces of the mesh of
    //  sigmaTotal, from sigmaTotal (Gi + Ge in the heart cells, the torso
    //  conductivity elsewhere), sigmaIntracellular (Gi in the heart cells,
    //  0 elsewhere) and the heart cells: the linear mean of Ge on heart
    //  faces, the torso face value on torso faces, the conormal harmonic
    //  value of Ge (Gi + Ge when totalWallConductivity) and the torso
    //  conductivity on heart-torso faces, and sigmaTotal - sigmaIntracellular
    //  of the cell on physical boundary faces
    static tmp<surfaceTensorField> faceConductivity
    (
        const volTensorField& sigmaTotal,
        const volTensorField& sigmaIntracellular,
        const labelUList& heartCells,
        const bool totalWallConductivity,
        const bool registerObject
    );


    //- Construct from the base mesh, the heart sub-mesh with its cell and
    //  face maps to the base mesh, Gi and Ge on the heart mesh, the torso
    //  block (its conductivity field is read from the base mesh), the seal
    //  flag and the heart field whose patch types decide sealing. With
    //  conormalHarmonic, nInterfaceGradientCorrectors (required, at least
    //  1) is read from the phiE entry of the solvers of fvSolution of the
    //  base mesh.
    //  totalWallConductivity selects Gi + Ge (true) or Ge (false) on the
    //  heart side of heart-torso faces. registerFields false keeps the
    //  base-mesh fields out of the base-mesh registry (their names stay the
    //  same).
    passiveConductor
    (
        const fvMesh& baseMesh,
        const fvMesh& heartMesh,
        const labelUList& heartCellMap,
        const labelUList& heartFaceMap,
        const volTensorField& Gi,
        const volTensorField& Ge,
        const torsoSettings& torso,
        const Switch sealedHeartBoundary,
        const volScalarField& sealingField,
        const bool totalWallConductivity,
        const bool registerFields = true
    );

    passiveConductor(const passiveConductor&) = delete;
    void operator=(const passiveConductor&) = delete;

    ~passiveConductor() = default;


    //- Assemble A psi for the base unknown psi and its heart companion
    //  heartPsi. A given source enters as == source after the base-mesh
    //  Laplacian (and the conormalHarmonic cross term), before the
    //  conormalHarmonic non-orthogonal term and the heart block. A given
    //  correction receives the explicit face terms. psi and heartPsi are
    //  named potentialName.
    tmp<fvScalarMatrix> assemble
    (
        volScalarField& psi,
        volScalarField& heartPsi,
        const volScalarField* source = nullptr,
        faceCorrection* correction = nullptr
    ) const;

    //- Set the heart companion heartPsi to the base values of psi in the
    //  heart cells and correct its boundary conditions
    void fillHeart(const volScalarField& psi, volScalarField& heartPsi) const;

    //- Set the gradient I/(n.sigmaTotal.n) of the surface current I on the
    //  surface-current patches of psi (fixedGradient) and evaluate them;
    //  fatal for a non-positive normal conductivity
    void setSurfaceCurrents(volScalarField& psi) const;

    //- Pin the reference cell to value (no-op when ground patches fix the
    //  constant)
    void applyReference(fvScalarMatrix& eqn, const scalar value) const;

    //- Stencil of the explicit part of the operator as probedOperator
    //  hops, outermost first: the divergence and face interpolation (face),
    //  one face hop per jump-aware gradient pass, and the gradient stencil
    //  (face when every gradientNames() scheme on the base and heart meshes
    //  is leastSquares or Gauss linear, point otherwise)
    List<probedOperator::hop> explicitStencil() const;

    //- Heart source b = -laplacian(insulated Gi, Vm), scattered to the
    //  base mesh (zero in torso cells)
    tmp<volScalarField> heartSource(const volScalarField& Vm) const;

    //- The torso block the conductor is built from
    const torsoSettings& torso() const
    {
        return torso_;
    }

    //- Pinned base cell on its owning processor; -1 elsewhere and when
    //  ground patches fix the constant
    label referenceCell() const
    {
        return referenceCell_;
    }

    const labelList& heartCellToBaseCell() const
    {
        return heartCellToBaseCell_;
    }

    //- Gi+Ge in heart cells, sigma_b in torso cells
    const volTensorField& sigmaTotal() const
    {
        return sigmaTotalPtr_();
    }

    //- Gi in heart cells, 0 in torso cells
    const volTensorField& sigmaIntracellular() const
    {
        return sigmaIglobalPtr_();
    }

    //- Jump-aware gradient passes of the corrections (conormalHarmonic),
    //  0 otherwise
    label nInterfaceGradientCorrectors() const
    {
        return nInterfaceGradientCorrectors_;
    }

    //- True when heart-torso faces carry Gi + Ge on the heart side (a sealed
    //  source), false when they carry Ge (bidomain)
    bool totalWallConductivity() const
    {
        return totalWallConductivity_;
    }
};

} // End namespace Foam

#endif

// ************************************************************************* //
