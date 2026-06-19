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
    graphVerificationModel

Description
    Base class for 1D graph conduction system verification models.
    Provides preProcess and postProcess hooks to compute analytical errors
    on the Purkinje graph nodes.

Author
    Simao Nieto de Castro. All rights reserved.
\*---------------------------------------------------------------------------*/

#ifndef graphVerificationModel_H
#define graphVerificationModel_H

#include "runTimeSelectionTables.H"
#include "dictionary.H"
#include "scalarField.H"
#include "pointField.H"
#include "Time.H"

namespace Foam
{

class ionicModel;

class graphVerificationModel
{
    dictionary dict_;

public:

    TypeName("graphVerificationModel");

    declareRunTimeSelectionTable
    (
        autoPtr,
        graphVerificationModel,
        dictionary,
        (
            const dictionary& dict
        ),
        (dict)
    );

    static autoPtr<graphVerificationModel> New(const dictionary& dict);
    static word selectedType(const dictionary& dict);

    explicit graphVerificationModel(const dictionary& dict);

    virtual ~graphVerificationModel() = default;

    virtual void preProcess
    (
        const Time& runTime,
        ionicModel* model,
        scalarField& Vm1D,
        const pointField& nodeLocations,
        label localStartNode
    ) = 0;

    virtual void postProcess
    (
        const Time& runTime,
        const ionicModel* model,
        const scalarField& Vm1D,
        const pointField& nodeLocations,
        label localStartNode
    ) = 0;

protected:

    const dictionary& dict() const { return dict_; }
};

} // End namespace Foam

#endif

// ************************************************************************* //
