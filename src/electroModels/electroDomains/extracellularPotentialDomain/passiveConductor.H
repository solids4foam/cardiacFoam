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
    Passive heart+torso conductor on the base mesh: the cell and face
    conductivities, the operator for a base-mesh potential and its heart
    companion, the reference pin and the heart intracellular source.

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
#include "dictionary.H"
#include "point.H"
#include "wordList.H"
#include "scalarList.H"

namespace Foam
{

class passiveConductor
{
public:

    //- Explicit face terms of the assembled operator, for flux diagnostics
    struct faceCorrection
    {
        //- Cell gradient of psi used by the face corrections
        //  (jump-aware under conormalHarmonic)
        tmp<volVectorField> gradient;

        //- Explicit face flux added to the operator (conormalHarmonic only)
        tmp<surfaceScalarField> crossFlux;
    };


    //- The torso/bath keys of a bathPotentialDomain block (or of an ECG
    //  domain torso sub-dictionary): the single reader of these keys.
    //  The interface option is checked by the passiveConductor constructor.
    struct torsoSettings
    {
        //- Torso cellZones (required)
        wordList bathCellZones;

        //- Torso scalar conductivity field name
        //  (default bodyAndOrgansConductivity)
        word bathConductivityField;

        //- Heart-torso faces carry the flux-continuous face conductivity;
        //  distanceWeightedHarmonic (default) uses the standard cell
        //  gradient for its cross term, conormalHarmonic the jump-aware one
        word interfaceConductivityInterpolation;

        //- Dirichlet patches and their values
        wordList groundPatches;
        scalarList groundValues;

        //- Neumann (surface-current) patches and their currents
        wordList surfaceCurrentPatches;
        scalarList surfaceCurrentValues;

        //- phiERefPoint is present; the point and the pinned value
        //  (default 0) are read whether or not ground patches exist
        bool hasPhiERefPoint;
        point phiERefPoint;
        scalar phiEReferenceValue;

        //- Read from dict. Every ground and surface-current patch must be
        //  a patch of mesh and no patch may be in both lists.
        static torsoSettings read
        (
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- FATAL unless phiERefPoint is present; owner names the caller
        void requirePhiERefPoint(const word& owner) const;
    };

    //- Read the torso scalar conductivity field fieldName on mesh
    //  (MUST_READ at its latest instance)
    static autoPtr<volScalarField> readTorsoConductivity
    (
        const fvMesh& mesh,
        const word& fieldName,
        const bool registerObject
    );


private:

    const fvMesh& baseMesh_;
    const fvMesh& heartMesh_;

    // Heart fields (non-owning).
    const volTensorField& Gi_;
    const volScalarField& sealingField_;

    // Face interpolation on heart-torso faces.
    const word interfaceConductivityInterpolation_;
    const Switch sealedHeartBoundary_;

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
        const volScalarField& torsoConductivity,
        const wordList& torsoCellZones
    );

    // Build the extracellular/torso face conductivity on all faces.
    void buildSigmaExtracellularSurface();


public:

    //- Construct from the base mesh, the heart sub-mesh with its cell and
    //  face maps to the base mesh, Gi and Ge on the heart mesh, the torso
    //  conductivity on the base mesh and its cellZones, the heart-torso face
    //  interpolation, the seal flag and the heart field whose patch types
    //  decide sealing. totalWallConductivity selects Gi + Ge (true) or Ge
    //  (false) on the heart side of heart-torso faces. registerFields false
    //  keeps the base-mesh fields out of the base-mesh registry (their names
    //  stay the same).
    passiveConductor
    (
        const fvMesh& baseMesh,
        const fvMesh& heartMesh,
        const labelUList& heartCellMap,
        const labelUList& heartFaceMap,
        const volTensorField& Gi,
        const volTensorField& Ge,
        const volScalarField& torsoConductivity,
        const wordList& torsoCellZones,
        const word& interfaceConductivityInterpolation,
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
    //  correction receives the explicit face terms.
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

    //- Pin refCell to value when psi has no fixedValue patch
    void applyReference
    (
        fvScalarMatrix& eqn,
        const volScalarField& psi,
        const label refCell,
        const scalar value
    ) const;

    //- Base cell containing p, owned by exactly one partition (-1 elsewhere)
    label referenceCell(const point& p) const;

    //- Heart source b = -laplacian(insulated Gi, Vm), scattered to the
    //  base mesh (zero in torso cells)
    tmp<volScalarField> heartSource(const volScalarField& Vm) const;

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

    //- Extracellular/torso face conductivity of the base-mesh operator
    const surfaceTensorField& sigmaExtracellularFace() const
    {
        return sigmaExtracellularfPtr_();
    }

    const word& interfaceConductivityInterpolation() const
    {
        return interfaceConductivityInterpolation_;
    }

    Switch sealedHeartBoundary() const
    {
        return sealedHeartBoundary_;
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
