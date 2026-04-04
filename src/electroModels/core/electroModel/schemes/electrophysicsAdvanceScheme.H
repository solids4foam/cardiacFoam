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

#ifndef electrophysicsAdvanceScheme_H
#define electrophysicsAdvanceScheme_H

#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "pimpleControl.H"

namespace Foam
{

class electrophysicsSystem;

struct electrophysicsAdvanceTimings
{
    scalar couplingTime;
    //- Time spent advancing upstream domains (e.g. Purkinje)
    scalar upstreamDomainTime;
    //- Time spent advancing the primary myocardium domain
    scalar primaryDomainTime;
    //- Time spent preparing downstream couplings (e.g. ECG)
    scalar downstreamCouplingTime;
    //- Time spent advancing downstream domains (e.g. ECG)
    scalar downstreamDomainTime;

    electrophysicsAdvanceTimings()
    :
        couplingTime(0.0),
        upstreamDomainTime(0.0),
        primaryDomainTime(0.0),
        downstreamCouplingTime(0.0),
        downstreamDomainTime(0.0)
    {}
};


class electrophysicsAdvanceScheme
{
public:

    TypeName("electrophysicsAdvanceScheme");

    declareRunTimeSelectionTable
    (
        autoPtr,
        electrophysicsAdvanceScheme,
        dictionary,
        (
            const dictionary& dict
        ),
        (dict)
    );

    electrophysicsAdvanceScheme() = default;

    static autoPtr<electrophysicsAdvanceScheme> New(const dictionary& dict);

    virtual ~electrophysicsAdvanceScheme() = default;

    virtual bool advance
    (
        scalar t0,
        scalar dt,
        electrophysicsSystem& system,
        pimpleControl* pimplePtr,
        electrophysicsAdvanceTimings& timings
    ) = 0;
};

} // End namespace Foam

#endif

// ************************************************************************* //
