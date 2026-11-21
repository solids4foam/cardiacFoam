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
    manufacturedShellBathBidomainVerifier

Description
    Manufactured bath bidomain verifier on the bath shell: the heart shell
    r1 < |x| < r2 (cellZone myocardium) coupled every step (torsoCoupling
    bath) to the bath ball bloodPool and the outer bath r2 < |x| < r3
    (cellZones bloodPool and bath). The intracellular space is sealed on both
    spheres and the extracellular current is continuous through Ge on the
    heart side; Gi and Ge are full, unequal and rotating through the wall,
    the bath conductivity varies in each bath zone. The outer sphere is
    insulated or grounded at c0. Fields from
    manufacturedShellBathBidomainReference.H
    (bath_shell_bidomain_mms.py).

    Keys of verificationModel:
        N               24;          // cells per unit length, for the file name
        outerClosure    insulated;   // insulated or grounded
        parameterSet    default;     // default, or alternate (common mode on)
        electrodePositions           // required; surface* are the body surface
        {
            surface1    (1.2228 0.8015 -0.1527);
            pool        (0.1131 -0.0733 0.2117);
        }

    The Vm equation takes the applied current beta Vm - div(Gi grad phiI) and
    the phiE equation the manufactured source of each zone at the cell
    centres; with the outer sphere insulated their volume-weighted mean over
    heart and bath is removed, so the pure-Neumann system is compatible on
    the mesh, and phiE is compared after removing the volume mean of its
    error. With the outer sphere grounded the source is used as it is and
    phiE is compared as it is. At each electrode the solved phiE of the cell
    that holds it is compared with the exact phiE at the cell centre (rows
    phiE*_electrode_<name>, and phiE*_bodySurface over the surface*
    electrodes); the exact value at the point is written beside it.

    The verifier needs sealedHeartBoundary true, the 3D
    monodomainFDAManufactured ionic model, the conductivities written by
    setManufacturedFields, and the outer patch insulated (no fixedValue
    patch of phiE) or grounded at c0 (every non-coupled patch of phiE
    fixedValue c0).

Author
    Simao Nieto de Castro, UCD.

SourceFiles
    manufacturedShellBathBidomainVerifier.C
\*---------------------------------------------------------------------------*/

#ifndef manufacturedShellBathBidomainVerifier_H
#define manufacturedShellBathBidomainVerifier_H

#include "electroVerificationModel.H"
#include "overrideTypeName.H"
#include "Enum.H"
#include "bathBidomainVerification/manufacturedShellBathBidomainReference.H"

namespace Foam
{

class manufacturedShellBathBidomainVerifier
:
    public electroVerificationModel
{
public:

    //- Data of the outer sphere r3
    enum class outerClosure
    {
        insulated,
        grounded
    };

    static const Enum<outerClosure> outerClosureNames;

private:

    shellBathBidomainMMS::Params p_;
    outerClosure closure_;
    word parameterSet_;
    label nPerDirection_;
    volScalarField* phiEPtr_;
    const labelUList* heartCellMapPtr_;
    mutable bool errorsReported_;
    scalar beta_;
    bool betaInitialised_;

    //- Electrodes: names and positions; those named surface* are the body
    //  surface
    wordList electrodeNames_;
    List<point> electrodePositions_;

    //- Zone of every base cell: the heart cells from the heart cell map, the
    //  bath cells from the cellZones bloodPool and bath
    mutable labelList zones_;

    //- Check the header functions against their self-test samples
    static void selfTest();

    void initialiseBeta(const ionicModel& model);

    //- Zones of the base mesh, checked against the radii of the spheres
    const labelList& zones(const fvMesh& baseMesh) const;

    void validateConfiguration(const ionicModel& model) const;

    //- Mesh points in the radial range of their zone, non-coupled boundary
    //  faces on r3, and the phiE patches against the outer closure
    void validateGeometry(const fvMesh& baseMesh) const;

    //- Gi and Ge on the heart mesh and sigma_b on the bath cells must be the
    //  manufactured ones
    void validateConductivity
    (
        const fvMesh& heartMesh,
        const fvMesh& baseMesh
    ) const;

    //- Exact phiE of the base cells at time t
    tmp<scalarField> exactPhiE(const fvMesh& baseMesh, const scalar t) const;

public:

    OverrideTypeName("manufacturedShellBathBidomainVerifier");

    explicit manufacturedShellBathBidomainVerifier(const dictionary& dict);

    virtual void bindBidomainField(volScalarField& phiE) override;

    virtual void bindBidomainField
    (
        volScalarField& phiE,
        const labelUList& heartCellMap
    ) override;

    virtual void unbindBidomainField() override;

    virtual wordList preProcessFieldNames(const ionicModel& model) const override;

    virtual wordList requiredPostProcessFieldNames
    (const ionicModel& model) const override;

    virtual bool shouldPostProcess
    (const ionicModel& model, const volScalarField& Vm) const override;

    virtual void preProcess
    (
        ionicModel& model,
        volScalarField& Vm,
        PtrList<volScalarField>& fields
    ) override;

    virtual void addManufacturedPdeSource
    (
        volScalarField& sourceField,
        const volTensorField& conductivity,
        const scalar evaluationTime
    ) const override;

    virtual void addManufacturedPhiESource
    (
        volScalarField& rhs,
        const scalar evaluationTime
    ) const override;

    //- The phiE rows read the gradient the bath operator uses
    virtual bool requiresInterfaceGradient() const override
    {
        return true;
    }

    virtual void writeManufacturedFields(const fvMesh& mesh) const override;

    virtual void postProcess
    (
        const ionicModel& model,
        const volScalarField& Vm,
        const PtrList<volScalarField>& fields
    ) override;
};

} // End namespace Foam

#endif

// ************************************************************************* //
