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
    manufacturedShellMonodomainVerifier

Description
    Manufactured monodomain verifier on the spherical shell r1 < |x| < r2,
    sealed on both spheres: V = sqrt(1 + t) F with n.G grad(F) = 0 on each
    sphere, G full and rotating through the wall (written cell-wise as
    0/Conductivity by setManufacturedFields), and the FDA ionic closure.
    Fields from manufacturedShellMonodomainReference.H.

Author
    Simao Nieto de Castro, UCD.

SourceFiles
    manufacturedShellMonodomainVerifier.C
\*---------------------------------------------------------------------------*/

#ifndef manufacturedShellMonodomainVerifier_H
#define manufacturedShellMonodomainVerifier_H

#include "electroVerificationModel.H"
#include "overrideTypeName.H"
#include "monodomainVerification/manufacturedShellMonodomainReference.H"

namespace Foam
{

class manufacturedShellMonodomainVerifier
:
    public electroVerificationModel
{
    shellMonodomainMMS::Params p_;
    label nPerDirection_;
    mutable bool errorsReported_;
    scalar beta_;
    bool betaInitialised_;

    //- Check the header functions against their self-test samples
    static void selfTest();

    void initialiseBeta(const ionicModel& model);

    void validateConfiguration(const ionicModel& model, const fvMesh& mesh) const;

    //- The conductivity must be the manufactured G at every cell centre
    void validateConductivity(const volTensorField& conductivity) const;

public:

    OverrideTypeName("manufacturedShellMonodomainVerifier");

    explicit manufacturedShellMonodomainVerifier(const dictionary& dict);

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
