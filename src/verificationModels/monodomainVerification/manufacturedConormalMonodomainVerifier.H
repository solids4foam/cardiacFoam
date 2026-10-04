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
    manufacturedConormalMonodomainVerifier

Description
    Manufactured monodomain verifier for a constant, fully populated
    anisotropic conductivity tensor on the unit cube, periodic in y and z,
    with a wall profile satisfying the true conormal flux condition rather
    than a zero-gradient one.

Author
    Simao Nieto de Castro, UCD.
\*---------------------------------------------------------------------------*/

#ifndef manufacturedConormalMonodomainVerifier_H
#define manufacturedConormalMonodomainVerifier_H

#include "conormalManufacturedProfile.H"
#include "electroVerificationModel.H"
#include "overrideTypeName.H"

namespace Foam
{

class manufacturedConormalMonodomainVerifier
:
    public electroVerificationModel
{
    mutable bool errorsReported_;
    scalar beta_;
    bool betaInitialised_;
    mutable bool profileInitialised_;
    mutable conormalProfile profile_;

    void initialiseBeta(const ionicModel& model);

    void validateConfiguration
    (
        const ionicModel& model,
        const fvMesh& mesh
    ) const;

    void initialiseProfile(const volTensorField& conductivity) const;

public:

    OverrideTypeName("manufacturedConormalMonodomainVerifier");

    explicit manufacturedConormalMonodomainVerifier(const dictionary& dict);

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
