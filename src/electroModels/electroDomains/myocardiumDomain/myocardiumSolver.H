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
#include "electroVolumeFieldDomain.H"

namespace Foam
{

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

    //- Construct from explicit solver type name and its Coeffs subdictionary.
    //  The solverType is resolved by the caller (electroModel::New / builder);
    //  no dict key lookup is performed here.
    static autoPtr<myocardiumSolver> New
    (
        const fvMesh& mesh,
        const word& solverType,
        const dictionary& coeffs
    );

    virtual ~myocardiumSolver() = default;

    virtual const volScalarField* phiEPtr() const
    {
        return nullptr;
    }

    virtual const volScalarField* phiIPtr() const
    {
        return nullptr;
    }

    virtual const volTensorField* conductivityPtr() const = 0;

    virtual const volTensorField* intracellularConductivityPtr() const
    {
        return conductivityPtr();
    }

    virtual const volTensorField* extracellularConductivityPtr() const
    {
        return conductivityPtr();
    }

    virtual const volTensorField& explicitConductivityField() const = 0;

    virtual void solveDiffusionExplicit
    (
        electroVolumeFieldDomain& domain,
        scalar dt
    ) = 0;

    virtual void solveDiffusionImplicit
    (
        electroVolumeFieldDomain& domain,
        scalar dt
    ) = 0;

    virtual void solveDiffusionImplicit
    (
        electroVolumeFieldDomain& domain,
        scalar dt,
        pimpleControl& pimple
    ) = 0;
};

} // End namespace Foam

#endif

// ************************************************************************* //
