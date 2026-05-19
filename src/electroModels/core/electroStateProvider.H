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
    electroStateProvider

Description
    Abstract read-only state interface that enables one-way coupling between
    electrophysiology domains while maintaining clean separation of concerns.

Author
    Simao Nieto de Castro. All rights reserved.
\*---------------------------------------------------------------------------*/

#ifndef electroStateProvider_H
#define electroStateProvider_H

#include "volFields.H"

namespace Foam
{

class electroStateProvider
{
public:

    virtual ~electroStateProvider() = default;

    //- The mesh on which this state provider is defined.
    virtual const fvMesh& mesh() const = 0;

    //- Base mesh from which this provider's mesh may have been derived.
    virtual const fvMesh& baseMesh() const
    {
        return mesh();
    }

    //- Optional mapping from a subset mesh back to the base mesh.
    virtual const labelUList* subsetFaceMapPtr() const
    {
        return nullptr;
    }

    virtual const labelUList* subsetCellMapPtr() const
    {
        return nullptr;
    }

    //- Optional state channels exposed by the provider.
    virtual const volScalarField* VmPtr() const
    {
        return nullptr;
    }

    virtual const volScalarField* phiEPtr() const
    {
        return nullptr;
    }

    virtual const volTensorField* conductivityPtr() const
    {
        return nullptr;
    }

    //- Optional split conductivity channels for models that distinguish
    //  between intracellular and extracellular media.
    //  By default these fall back to the legacy conductivityPtr() channel.
    virtual const volTensorField* intracellularConductivityPtr() const
    {
        return conductivityPtr();
    }

    virtual const volTensorField* extracellularConductivityPtr() const
    {
        return conductivityPtr();
    }
};

} // End namespace Foam

#endif

// ************************************************************************* //
