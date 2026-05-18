/*---------------------------------------------------------------------------*\
License
    This file is part of cardiacFoam.

    cardiacFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    cardiacFoam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
    or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License along
    with cardiacFoam.  If not, see <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#ifndef manufacturedFDABathBidomainVerifier_H
#define manufacturedFDABathBidomainVerifier_H

#include "electroVerificationModel.H"
#include "overrideTypeName.H"

namespace Foam
{

class manufacturedFDABathBidomainVerifier
:
    public electroVerificationModel
{
    volScalarField* phiEPtr_;
    const labelUList* phiEHeartCellMapPtr_;
    Switch enabled_;
    bool useExplicitAlgorithm_;
    mutable bool errorsReported_;
    scalar k_;
    scalar alpha_;
    Switch groundElectrode_;
    fileName outputFileName_;
    word phiEErrorCellZoneName_;

    const dictionary& verificationDict() const;

public:

    OverrideTypeName("manufacturedFDABathBidomainVerifier");

    explicit manufacturedFDABathBidomainVerifier(const dictionary& dict);

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
