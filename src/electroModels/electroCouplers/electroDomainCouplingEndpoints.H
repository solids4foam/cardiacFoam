/*---------------------------------------------------------------------------*\
License
    This file is part of cardiacFoam.

    cardiacFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    cardiacFoam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
    Public License for more details.

    You should have received a copy of the GNU General Public License
    along with cardiacFoam.  If not, see <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#ifndef electroDomainCouplingEndpoints_H
#define electroDomainCouplingEndpoints_H

#include "fvMesh.H"
#include "labelList.H"
#include "pointField.H"
#include "volFields.H"

namespace Foam
{

//- Primary tissue domain interface used by coupling models.
class tissueCouplingEndpoint
{
public:

    virtual ~tissueCouplingEndpoint() = default;

    //- Mesh associated with the domain.
    virtual const fvMesh& mesh() const = 0;

    //- Source term field exposed to coupling models.
    virtual volScalarField& sourceField() = 0;

    //- Transmembrane voltage at the domain.
    virtual const volScalarField& Vm() const = 0;
};


//- Terminal network interface used by conduction-network domains.
class networkCouplingEndpoint
{
public:

    virtual ~networkCouplingEndpoint() = default;

    virtual const labelList& terminalNodes() const = 0;

    virtual const pointField& terminalLocations() const = 0;

    virtual void terminalVm(scalarField& values) const = 0;

    //- Receive terminal exchange terms prepared by the coupling model.
    //  Positive terminalCurrent removes charge from the network and injects
    //  it into the coupled myocardium. terminalSource carries the matching
    //  myocardium-side volumetric source representation for diagnostics.
    virtual void setTerminalCoupling
    (
        const scalarField& terminalCurrent,
        const scalarField& terminalSource
    ) = 0;
};

} // End namespace Foam

#endif

// ************************************************************************* //
