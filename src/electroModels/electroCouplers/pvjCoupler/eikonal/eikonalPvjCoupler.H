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
    eikonalPvjCoupler

Description
    Eikonal PVJ coupler. seedProfile sets the seed time of the junction cells:
    flat (the terminal time), front (plus the metric distance over c0) or
    curvedFront (plus the radial source-free solution from the metric radius
    2/c0).

Author
    Simao Nieto de Castro. All rights reserved.
\*---------------------------------------------------------------------------*/

#ifndef eikonalPvjCoupler_H
#define eikonalPvjCoupler_H

#include "pvjCoupler.H"
#include "overrideTypeName.H"

namespace Foam
{

class eikonalPvjCoupler
:
    public pvjCoupler
{
    mutable scalarField terminalActivationBuffer_;
    const word seedProfile_;
    bool seedOffsetsSet_;

public:

    OverrideTypeName("eikonalPvjCoupler");

    //- The seedProfile of a coupler dictionary, checked against the valid
    //  profiles; front when absent
    static word seedProfile(const dictionary& dict);

    eikonalPvjCoupler
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
