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
    manufacturedTransmissionBathBidomainVerifier

Description
    Manufactured transmission verifier for the heart-in-bath bidomain model
    on the oblique slab: a myocardium slab between two bath slabs, with
    planar heart-bath interfaces oblique to the mesh axes, rotated
    anisotropic heart conductivities and spatially varying bath
    conductivities.

    The exact fields, conductivities and sources come from
    manufacturedTransmissionBathBidomainReference.H. The verifier drives the
    myocardium source, the global extracellular-potential source and the
    Dirichlet data of the outer bath patches, checks the mesh zones, the
    exposed-face patch and the cell conductivities against the manufactured
    solution, and reports field, interface-flux, interface-current,
    gradient and outer-current errors.

    The face conductivity the solver assembles on the heart-bath interface
    is internal to the solver, so its positivity is not checked here.

Author
    Simao Nieto de Castro. All rights reserved.
\*---------------------------------------------------------------------------*/

#ifndef manufacturedTransmissionBathBidomainVerifier_H
#define manufacturedTransmissionBathBidomainVerifier_H

#include "electroVerificationModel.H"
#include "manufacturedTransmissionBathBidomainReference.H"
#include "overrideTypeName.H"

namespace Foam
{

class manufacturedTransmissionBathBidomainVerifier
:
    public electroVerificationModel
{
    volScalarField* phiEPtr_;
    const labelUList* phiEHeartCellMapPtr_;
    mutable bool errorsReported_;
    transmissionMMS::Params params_;
    word solutionName_;
    label nPerDirection_;
    wordList dirichletPatches_;
    scalar beta_;
    scalar capacitanceRatio_;
    bool betaInitialised_;
    labelList zoneOfCell_;
    const fvMesh* zoneMeshPtr_;

    const dictionary& verificationDict() const;

    //- Exact-solution parameters from a verification dictionary.
    static transmissionMMS::Params readParams(const dictionary& dict);

    //- Check the reader against the reference defaults and samples
    static void selfTest();

    //- Zone index of every cell from the bathLow, myocardium and bathHigh cellZones
    static labelList buildZones
    (
        const fvMesh& mesh,
        const transmissionMMS::Params& p
    );

    static scalar ionicConstant(const ionicModel& model, const word& name);

    void validateConfiguration(const ionicModel& model) const;

    //- Beta, Chi and Cm of the ionic model, taken from ranks holding heart cells
    void initialiseIonicConstants
    (
        const ionicModel& model,
        const label nLocalCells
    );

    void validateBoundBath(const volScalarField& Vm) const;

public:

    OverrideTypeName("manufacturedTransmissionBathBidomainVerifier");

    explicit manufacturedTransmissionBathBidomainVerifier
    (
        const dictionary& dict
    );

    virtual void bindBidomainField(volScalarField& phiE) override;

    virtual void bindBidomainField
    (
        volScalarField& phiE,
        const labelUList& heartCellMap
    ) override;

    virtual void unbindBidomainField() override;

    virtual wordList preProcessFieldNames
    (const ionicModel& model) const override;

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

    virtual void applyPhiEBoundaryData
    (
        volScalarField& phiE,
        const scalar evaluationTime
    ) const override;

    virtual wordList boundaryDataPatches() const override;

    virtual bool requiresPhiECompatibilityProjection() const override;

    virtual bool requiresInterfaceFluxDiagnostic() const override;

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
