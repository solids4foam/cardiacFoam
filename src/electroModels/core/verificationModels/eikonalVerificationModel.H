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
    eikonalVerificationModel

Description
    Base class for eikonal activation-time verification models.

Author
    Simao Nieto de Castro. All rights reserved.
\*---------------------------------------------------------------------------*/

#ifndef eikonalVerificationModel_H
#define eikonalVerificationModel_H

#include "runTimeSelectionTables.H"
#include "dictionary.H"
#include "volFields.H"
#include "Switch.H"
#include "dimensionedScalar.H"

namespace Foam
{

class eikonalVerificationModel
{
    dictionary dict_;

public:

    TypeName("eikonalVerificationModel");

    declareRunTimeSelectionTable
    (
        autoPtr,
        eikonalVerificationModel,
        dictionary,
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const volTensorField& conductivity,
            const dimensionedScalar& chi,
            const dimensionedScalar& Cm,
            const dimensionedScalar& c0,
            const Switch& eikonalAdvectionDiffusionApproach
        ),
        (dict, mesh, conductivity, chi, Cm, c0, eikonalAdvectionDiffusionApproach)
    );

    static autoPtr<eikonalVerificationModel> New
    (
        const dictionary& dict,
        const fvMesh& mesh,
        const volTensorField& conductivity,
        const dimensionedScalar& chi,
        const dimensionedScalar& Cm,
        const dimensionedScalar& c0,
        const Switch& eikonalAdvectionDiffusionApproach
    );

    static word selectedType(const dictionary& dict);

    explicit eikonalVerificationModel(const dictionary& dict);

    virtual ~eikonalVerificationModel() = default;

    virtual bool enabled() const = 0;

    virtual void applyConstraints(volScalarField& activationTime) const = 0;

    virtual tmp<volScalarField> sourceTerm() const = 0;

    virtual bool shouldPostProcess() const = 0;

    virtual void postProcess(const volScalarField& activationTime) = 0;

protected:

    const dictionary& dict() const { return dict_; }
};

} // End namespace Foam

#endif

// ************************************************************************* //
