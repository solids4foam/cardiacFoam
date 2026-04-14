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

/*---------------------------------------------------------------------------*\
  NOT_COMPILED
  This file is not included in Make/files. The heart-bath coupler is not
  currently wired into the active electrophysicsSystemBuilder path.
  Requires BathDomain to be activated first.
\*---------------------------------------------------------------------------*/

#ifndef heartBathInterfaceCoupler_H
#define heartBathInterfaceCoupler_H

#include "Switch.H"
#include "electroDomainCoupler.H"
#include "electroStateProvider.H"
#include "overrideTypeName.H"

namespace Foam
{

class heartBathInterfaceCoupler
:
    public ElectroDomainCoupler
{
    bathCouplingEndpoint& bathDomain_;
    const electroStateProvider& heartStateProvider_;
    Switch copyHeartPhiEToBath_;
    Switch reportCouplingExtrema_;

    static bathCouplingEndpoint& requireBathDomain
    (
        electroDomainInterface& secondaryDomain
    );

    static const electroStateProvider& requireHeartStateProvider
    (
        tissueCouplingEndpoint& primaryDomain
    );

public:

    OverrideTypeName("heartBathInterfaceCoupler");

    heartBathInterfaceCoupler
    (
        tissueCouplingEndpoint& primaryDomain,
        electroDomainInterface& secondaryDomain,
        const dictionary& dict
    );

    virtual void preparePostPrimaryCoupling(scalar t0, scalar dt) override;
};

} // End namespace Foam

#endif

// ************************************************************************* //
