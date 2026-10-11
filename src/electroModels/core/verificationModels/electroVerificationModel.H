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
    electroVerificationModel

Description
    Base class for electrophysiology verification models.

Author
    Simao Nieto de Castro. All rights reserved.
\*---------------------------------------------------------------------------*/

#ifndef electroVerificationModel_H
#define electroVerificationModel_H

#include "runTimeSelectionTables.H"
#include "dictionary.H"
#include "volFields.H"

namespace Foam
{

class ionicModel;
class conductionSystemDomain;

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

    virtual void bindSourceField(const volScalarField& sourceField)
    {
        (void)sourceField;
    }

    //- Add a continuous manufactured PDE source after the ordinary external
    //  stimulus has been reset and before domain coupling deposits current.
    virtual void addManufacturedPdeSource
    (
        volScalarField& sourceField,
        const volTensorField& conductivity,
        const scalar evaluationTime
    ) const
    {
        (void)sourceField;
        (void)conductivity;
        (void)evaluationTime;
    }

    //- Add the manufactured source of the phiE equation, L(phiE) == rhs, on the mesh phiE is solved on
    virtual void addManufacturedPhiESource
    (
        volScalarField& rhs,
        const scalar evaluationTime
    ) const
    {
        (void)rhs;
        (void)evaluationTime;
    }

    //- Set fixedValue values and fixedGradient gradients on the patches of boundaryDataPatches(); a gradient g gives n . sigmaTotal(owner cell) . n * g = prescribed flux
    virtual void applyPhiEBoundaryData
    (
        volScalarField& phiE,
        const scalar evaluationTime
    ) const
    {
        (void)phiE;
        (void)evaluationTime;
    }

    //- Extracellular-potential patches set by applyPhiEBoundaryData.
    virtual wordList boundaryDataPatches() const
    {
        return wordList();
    }

    //- Remove the discrete imbalance of a pure-Neumann phiE system before the reference cell is pinned
    virtual bool requiresPhiECompatibilityProjection() const
    {
        return false;
    }

    //- Store the face flux balanced by the phiE equation on heart-bath faces in phiEInterfaceFlux (A, owner->neighbour, extracellular part) and the gradient of its face corrections in phiEInterfaceGradient
    virtual bool requiresInterfaceFluxDiagnostic() const
    {
        return false;
    }

    //- Write the fields the manufactured solution prescribes.
    virtual void writeManufacturedFields(const fvMesh& mesh) const
    {
        (void)mesh;

        FatalErrorInFunction
            << "The selected verification model does not write "
            << "manufactured fields."
            << exit(FatalError);
    }

    virtual void bindBidomainField
    (
        volScalarField& phiE,
        const labelUList& heartCellMap
    )
    {
        (void)heartCellMap;
        bindBidomainField(phiE);
    }

    virtual void unbindBidomainField() {}

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
