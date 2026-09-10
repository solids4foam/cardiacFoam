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
    electromechanicalVerificationModel

Description
    Base class for coupled electromechanics verification hooks.

SourceFiles
    electromechanicalVerificationModel.C

Author
    Simao Nieto de Castro, UCD.
\*---------------------------------------------------------------------------*/

#ifndef electromechanicalVerificationModel_H
#define electromechanicalVerificationModel_H

#include "autoPtr.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"

namespace Foam
{

class electromechanicalVerificationModel
{
    dictionary dict_;

public:

    TypeName("electromechanicalVerificationModel");

    declareRunTimeSelectionTable
    (
        autoPtr,
        electromechanicalVerificationModel,
        dictionary,
        (const dictionary& dict),
        (dict)
    );

    explicit electromechanicalVerificationModel(const dictionary& dict);

    virtual ~electromechanicalVerificationModel() = default;

    static bool configured(const dictionary& dict);

    static autoPtr<electromechanicalVerificationModel> New
    (
        const dictionary& dict
    );

    virtual void initialize
    (
        volScalarField& Vm,
        volVectorField& D
    ) = 0;

    virtual void preSolve
    (
        volScalarField& Vm,
        volVectorField& D
    ) = 0;

    virtual bool shouldPostProcess
    (
        const volScalarField& Vm,
        const volVectorField& D
    ) const = 0;

    virtual void postProcess
    (
        const volScalarField& Vm,
        const volVectorField& D,
        const volScalarField& Ta
    ) = 0;

protected:

    const dictionary& dict() const
    {
        return dict_;
    }
};

} // End namespace Foam

#endif

// ************************************************************************* //
