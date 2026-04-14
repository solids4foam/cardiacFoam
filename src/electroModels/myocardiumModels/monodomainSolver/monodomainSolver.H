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

#ifndef MonodomainSolver_H
#define MonodomainSolver_H

#include "myocardiumSolver.H"
#include "overrideTypeName.H"
#include "fvc.H"
#include "fvm.H"

namespace Foam
{

//- Monodomain FVM reaction-diffusion kernel.
//  Registered in the myocardiumSolver runtime table under
//  the key "monodomainSolver".
class MonodomainSolver
:
    public myocardiumSolver
{
    volTensorField conductivity_;

    tmp<volTensorField> initialiseConductivity
    (
        const fvMesh& mesh,
        const dictionary& electroProperties
    ) const;

public:

    OverrideTypeName("monodomainSolver");

    MonodomainSolver
    (
        const fvMesh& mesh,
        const dictionary& electroProperties
    );

    virtual const volTensorField* conductivityPtr() const override
    {
        return &conductivity_;
    }

    virtual const volTensorField& explicitConductivityField() const override
    {
        return conductivity_;
    }

    virtual void solveDiffusionExplicit
    (
        electroVolumeFieldDomain& domain,
        scalar dt
    ) override;

    virtual void solveDiffusionImplicit
    (
        electroVolumeFieldDomain& domain,
        scalar dt
    ) override;

    virtual void solveDiffusionImplicit
    (
        electroVolumeFieldDomain& domain,
        scalar dt,
        pimpleControl& pimple
    ) override;
};

} // End namespace Foam

#endif

// ************************************************************************* //
