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
    couplingVerificationModel

Description
    Base class for electrophysiology coupled domain verification models.

Author
    Simao Nieto de Castro. All rights reserved.
\*---------------------------------------------------------------------------*/

#ifndef couplingVerificationModel_H
#define couplingVerificationModel_H

#include "runTimeSelectionTables.H"
#include "dictionary.H"

namespace Foam
{

class tissueCouplingEndpoint;
class electroDomainInterface;

class couplingVerificationModel
{
    dictionary dict_;

public:

    TypeName("couplingVerificationModel");

    declareRunTimeSelectionTable
    (
        autoPtr,
        couplingVerificationModel,
        dictionary,
        (const dictionary& dict),
        (dict)
    );

    static autoPtr<couplingVerificationModel> New(const dictionary& dict);

    static word selectedType(const dictionary& dict);

    explicit couplingVerificationModel(const dictionary& dict);

    virtual ~couplingVerificationModel() = default;

    virtual void preProcess
    (
        tissueCouplingEndpoint& primaryDomain,
        electroDomainInterface& secondaryDomain
    ) = 0;

    virtual void updateManufacturedSource
    (
        tissueCouplingEndpoint& primaryDomain,
        electroDomainInterface& secondaryDomain,
        scalar primaryTime,
        scalar secondaryTime,
        bool implicitCoupling,
        bool bidirectionalCoupling,
        const word& phaseName
    )
    {
        (void)primaryDomain;
        (void)secondaryDomain;
        (void)primaryTime;
        (void)secondaryTime;
        (void)implicitCoupling;
        (void)bidirectionalCoupling;
        (void)phaseName;
    }

    virtual void postProcess
    (
        tissueCouplingEndpoint& primaryDomain,
        electroDomainInterface& secondaryDomain
    ) = 0;

protected:

    const dictionary& dict() const { return dict_; }
};

} // End namespace Foam

#endif

// ************************************************************************* //
