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

#ifndef electroVolumeFieldDomain_H
#define electroVolumeFieldDomain_H

#include "dimensionedScalar.H"
#include "fvMesh.H"
#include "volFields.H"

namespace Foam
{

//- Minimal field contract for volume electrophysiology diffusion solvers.
class electroVolumeFieldDomain
{
public:

    virtual ~electroVolumeFieldDomain() = default;

    virtual const fvMesh& mesh() const = 0;

    virtual volScalarField& VmRef() = 0;

    virtual const volScalarField& Vm() const = 0;

    virtual const volScalarField& Iion() const = 0;

    virtual volScalarField& sourceField() = 0;

    virtual const dimensionedScalar& chi() const = 0;

    virtual const dimensionedScalar& Cm() const = 0;
};

} // End namespace Foam

#endif

// ************************************************************************* //
