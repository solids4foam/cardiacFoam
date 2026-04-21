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

\*---------------------------------------------------------------------------*/

#ifndef electroVerificationModel_H
#define electroVerificationModel_H

#include "runTimeSelectionTables.H"
#include "dictionary.H"
#include "volFields.H"

namespace Foam
{

class ionicModel;

class electroVerificationModel
{
    dictionary dict_;

public:

    TypeName("electroVerificationModel");

    declareRunTimeSelectionTable
    (
        autoPtr,
        electroVerificationModel,
        dictionary,
        (const dictionary& dict),
        (dict)
    );

    static autoPtr<electroVerificationModel> New(const dictionary& dict);

    static word selectedType(const dictionary& dict);

    explicit electroVerificationModel(const dictionary& dict);

    static void allocateFields
    (
        const wordList& fieldNames,
        const fvMesh& mesh,
        const word& prefix,
        PtrList<volScalarField>& fields
    );

    virtual ~electroVerificationModel() = default;

    virtual void preProcess
    (
        ionicModel& model,
        volScalarField& Vm,
        PtrList<volScalarField>& fields
    ) = 0;

    virtual wordList preProcessFieldNames(const ionicModel& model) const = 0;

    virtual wordList requiredPostProcessFieldNames
    (const ionicModel& model) const = 0;

    virtual bool shouldPostProcess
    (const ionicModel& model, const volScalarField& Vm) const = 0;

    virtual void bindBidomainField(volScalarField& phiE) {}

    virtual void postProcess
    (
        const ionicModel& model,
        const volScalarField& Vm,
        const PtrList<volScalarField>& fields
    ) = 0;

protected:

    const dictionary& dict() const { return dict_; }
};

} // End namespace Foam

#endif

// ************************************************************************* //
