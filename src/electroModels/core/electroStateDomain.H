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

#ifndef electroStateDomain_H
#define electroStateDomain_H

#include "autoPtr.H"
#include "dictionary.H"
#include "electroDomainInterface.H"
#include "electroStateProvider.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

class fvMesh;
class myocardiumDomainInterface;

//- Electro domain that both advances in time and exposes read-only state.
class electroStateDomain
:
    public electroDomainInterface,
    public electroStateProvider
{
public:

    TypeName("electroStateDomain");

    declareRunTimeSelectionTable
    (
        autoPtr,
        electroStateDomain,
        dictionary,
        (
            const fvMesh& baseMesh,
            myocardiumDomainInterface& myocardium,
            const dictionary& dict
        ),
        (baseMesh, myocardium, dict)
    );

    static autoPtr<electroStateDomain> New
    (
        const fvMesh& baseMesh,
        myocardiumDomainInterface& myocardium,
        const dictionary& dict
    );

    virtual ~electroStateDomain() = default;
};

} // End namespace Foam

#endif

// ************************************************************************* //
