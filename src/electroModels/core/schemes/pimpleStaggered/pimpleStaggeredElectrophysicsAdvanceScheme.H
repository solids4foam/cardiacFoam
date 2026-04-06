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

    You should have received a copy of the GNU General Public License
    along with cardiacFoam.  If not, see <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#ifndef pimpleStaggeredElectrophysicsAdvanceScheme_H
#define pimpleStaggeredElectrophysicsAdvanceScheme_H

#include "../electrophysicsAdvanceScheme.H"
#include "overrideTypeName.H"

namespace Foam
{

class pimpleStaggeredElectrophysicsAdvanceScheme
:
    public electrophysicsAdvanceScheme
{
public:

    OverrideTypeName("pimpleStaggeredElectrophysicsAdvanceScheme");

    explicit pimpleStaggeredElectrophysicsAdvanceScheme(const dictionary&)
    {}

    virtual bool advance
    (
        scalar t0,
        scalar dt,
        electrophysicsSystem& system,
        pimpleControl* pimplePtr,
        electrophysicsAdvanceTimings& timings
    ) override;
};

} // End namespace Foam

#endif

// ************************************************************************* //
