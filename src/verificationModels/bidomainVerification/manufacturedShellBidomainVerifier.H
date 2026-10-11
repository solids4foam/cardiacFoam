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
    manufacturedShellBidomainVerifier

Description
    Manufactured heart-only bidomain verifier on the spherical shell
    r1 < |x| < r2, sealed on both spheres: phiI = sqrt(1 + t) Fi and
    phiE = sqrt(1 + t) Fe with n.Gi grad(Fi) = 0 and n.Ge grad(Fe) = 0 on each
    sphere, Gi and Ge full, unequal and rotating through the wall (written
    cell-wise as 0/ConductivityIntracellular and 0/ConductivityExtracellular
    by setManufacturedFields), Vm = phiI - phiE and the FDA ionic closure.
    The phiE equation takes the manufactured source and the compatibility
    projection of its pure-Neumann system. Fields from
    manufacturedShellBidomainReference.H.

Author
    Simao Nieto de Castro, UCD.

SourceFiles
    manufacturedShellBidomainVerifier.C
\*---------------------------------------------------------------------------*/

#ifndef manufacturedShellBidomainVerifier_H
#define manufacturedShellBidomainVerifier_H

#include "electroVerificationModel.H"
#include "overrideTypeName.H"
#include "bidomainVerification/manufacturedShellBidomainReference.H"

namespace Foam
{

class manufacturedShellBidomainVerifier
:
    public electroVerificationModel
{
    shellBidomainMMS::Params p_;
    label nPerDirection_;
    volScalarField* phiEPtr_;
    mutable bool errorsReported_;
    scalar phiEReferenceValue_;
    point phiEReferencePoint_;
    scalar beta_;
    bool betaInitialised_;

    //- Check the header functions against their self-test samples
    static void selfTest();

    void initialiseBeta(const ionicModel& model);

    void validateConfiguration(const ionicModel& model, const fvMesh& mesh) const;

    //- The conductivities must be the manufactured Gi and Ge at every cell centre
    void validateConductivity(const fvMesh& mesh) const;

    //- Constant that gives the exact phiE the value phiEReferenceValue at phiERefPoint
    scalar phiEShift(const scalar t) const;

public:

    OverrideTypeName("manufacturedShellBidomainVerifier");

    explicit manufacturedShellBidomainVerifier(const dictionary& dict);

    virtual void bindBidomainField(volScalarField& phiE) override;

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

    virtual bool requiresPhiECompatibilityProjection() const override
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
