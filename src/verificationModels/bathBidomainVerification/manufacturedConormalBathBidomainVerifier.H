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
    manufacturedConormalBathBidomainVerifier

Description
    Manufactured conormal bath bidomain verifier.

Author
    Simao Nieto de Castro. All rights reserved.
\*---------------------------------------------------------------------------*/

#ifndef manufacturedConormalBathBidomainVerifier_H
#define manufacturedConormalBathBidomainVerifier_H

#include "conormalManufacturedProfile.H"
#include "electroVerificationModel.H"
#include "manufacturedConormalBathBidomainReference.H"
#include "overrideTypeName.H"

namespace Foam
{

class manufacturedConormalBathBidomainVerifier
:
    public electroVerificationModel
{
    volScalarField* phiEPtr_;
    const labelUList* phiEHeartCellMapPtr_;
    mutable bool errorsReported_;
    scalar k_;
    scalar beta_;
    bool betaInitialised_;
    mutable bool conductivityValidated_;
    bool coeffsInitialised_;
    conormalProfile profile_;
    conormalBathCoefficients coeffs_;

    static scalar volumeMean(const fvMesh& m, const scalarField& f);

    const dictionary& verificationDict() const;

    void initialiseBeta(const ionicModel& model);

    void initialiseCoefficients(const tensor& Gi, const scalar sigmaB);

    void validateConfiguration(const ionicModel& model) const;

    void validateConductivity
    (
        const volTensorField& conductivityIntracellular,
        const volTensorField& conductivityExtracellular
    ) const;

public:

    OverrideTypeName("manufacturedConormalBathBidomainVerifier");

    explicit manufacturedConormalBathBidomainVerifier(const dictionary& dict);

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
