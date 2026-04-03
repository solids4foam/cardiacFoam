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

#ifndef manufacturedFDAMonodomainVerifier_H
#define manufacturedFDAMonodomainVerifier_H

#include "electroVerificationModel.H"
#include "overrideTypeName.H"

namespace Foam
{

class manufacturedFDAMonodomainVerifier : public electroVerificationModel
{
    bool useExplicitAlgorithm_;
    mutable bool errorsReported_;

    wordList readConfiguredNames
    (
        const word& subDictName,
        const word& entryName,
        const wordList& defaults
    ) const;

public:

    OverrideTypeName("manufacturedFDAMonodomainVerifier");

    explicit manufacturedFDAMonodomainVerifier(const dictionary& dict);

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
