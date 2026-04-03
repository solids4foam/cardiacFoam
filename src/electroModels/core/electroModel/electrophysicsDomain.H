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

#ifndef electrophysicsDomain_H
#define electrophysicsDomain_H

#include "Time.H"

namespace Foam
{

class electrophysicsDomain
{
public:

    virtual ~electrophysicsDomain() = default;

    //- Time database driving this domain.
    //  Not every electro domain is naturally an fvMesh-backed region
    //  (for example, a graph-based Purkinje network), so the base contract
    //  only requires access to the shared timestep state.
    virtual const Time& time() const = 0;

    //- Optional setup hook before a system-level timestep advance.
    virtual void prepareTimeStep
    (
        scalar t0,
        scalar dt
    )
    {}

    //- Advance this electrophysiology domain by one timestep.
    virtual void advance
    (
        scalar t0,
        scalar dt
    ) = 0;

    //- Optional output and teardown hooks used by the owning electroModel.
    virtual void write() {}

    virtual void end() {}
};

} // End namespace Foam

#endif

// ************************************************************************* //
