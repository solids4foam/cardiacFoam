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
    eikonalMonodomainPvjCoupler

Description
    Eikonal-to-Monodomain PVJ coupler. Uses voltage template for anterograde
    coupling and gathers activation times for retrograde coupling.

Author
    Simao Nieto de Castro. All rights reserved.
\*---------------------------------------------------------------------------*/

#ifndef eikonalMonodomainPvjCoupler_H
#define eikonalMonodomainPvjCoupler_H

#include "pvjCoupler.H"
#include "overrideTypeName.H"

namespace Foam
{

class eikonalMonodomainPvjCoupler
:
    public pvjCoupler
{
    mutable scalarField terminalActivationBuffer_;
    scalarField R_pvj_;

public:

    OverrideTypeName("eikonalMonodomainPvjCoupler");

    eikonalMonodomainPvjCoupler
    (
        tissueCouplingEndpoint& primaryDomain,
        electroDomainInterface& secondaryDomain,
        const dictionary& dict
    );

    virtual void prepareSecondaryCoupling(scalar t0, scalar dt) override;

    virtual void preparePrimaryCoupling(scalar t0, scalar dt) override;
};

} // End namespace Foam

#endif

// ************************************************************************* //
