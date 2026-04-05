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

#ifndef myocardiumSolver_H
#define myocardiumSolver_H

#include "pimpleControl.H"
#include "volFields.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "dictionary.H"
#include "fvMesh.H"

namespace Foam
{

class MyocardiumDomain;

class myocardiumSolver
{
public:

    TypeName("myocardiumSolver");

    declareRunTimeSelectionTable
    (
        autoPtr,
        myocardiumSolver,
        dictionary,
        (
            const fvMesh& mesh,
            const dictionary& electroProperties
        ),
        (mesh, electroProperties)
    );

    //- Select from the "myocardiumSolver" key in electroProperties.
    //  For bidomain, call the overload that accepts a phiE reference.
    static autoPtr<myocardiumSolver> New
    (
        const fvMesh& mesh,
        const dictionary& electroProperties
    );

    //- Bidomain overload: accepts the pre-built phiE field.
    static autoPtr<myocardiumSolver> NewBidomain
    (
        const fvMesh& mesh,
        const dictionary& electroProperties,
        volScalarField& phiE
    );

    virtual ~myocardiumSolver() = default;

    virtual const volScalarField* phiEPtr() const
    {
        return nullptr;
    }

    virtual const volTensorField* conductivityPtr() const = 0;

    virtual const volTensorField& explicitConductivityField() const = 0;

    virtual void solveDiffusionExplicit
    (
        MyocardiumDomain& domain,
        scalar dt
    ) = 0;

    virtual void solveDiffusionImplicit
    (
        MyocardiumDomain& domain,
        scalar dt,
        pimpleControl& pimple
    ) = 0;
};

} // End namespace Foam

#endif

// ************************************************************************* //
