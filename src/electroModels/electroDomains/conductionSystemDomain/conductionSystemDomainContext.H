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

#ifndef ConductionSystemDomainContext_H
#define ConductionSystemDomainContext_H

#include "Time.H"
#include "fvMesh.H"
#include "error.H"

namespace Foam
{

//- Execution context for conduction-system domains.
//  The abstract conduction-system layer only requires time integration.
//  Legacy implementations may additionally consume a support mesh for PVJ
//  lookup or output placement during the transition to mesh-independent
//  graph-backed domains.
class ConductionSystemDomainContext
{
    const Time& time_;
    const fvMesh* supportMeshPtr_;

public:

    explicit ConductionSystemDomainContext
    (
        const Time& time,
        const fvMesh* supportMeshPtr = nullptr
    )
    :
        time_(time),
        supportMeshPtr_(supportMeshPtr)
    {}

    const Time& time() const
    {
        return time_;
    }

    bool hasSupportMesh() const
    {
        return supportMeshPtr_ != nullptr;
    }

    const fvMesh& supportMesh() const
    {
        if (!supportMeshPtr_)
        {
            FatalErrorInFunction
                << "This conduction-system domain context does not provide a "
                << "legacy support mesh."
                << exit(FatalError);
        }

        return *supportMeshPtr_;
    }
};

} // End namespace Foam

#endif

// ************************************************************************* //
