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

#ifndef pvjCoupler_H
#define pvjCoupler_H

#include "electroDomainCoupler.H"
#include "pvjMapper.H"
#include "Switch.H"

namespace Foam
{

class PVJCoupler
:
    public ElectroDomainCoupler
{
public:

    enum CouplingMode
    {
        unidirectional,
        bidirectional
    };

protected:

    networkCouplingEndpoint& networkTerminalDomain_;
    PVJMapper mapper_;
    scalar pvjRadius_;
    CouplingMode couplingMode_;
    Switch reportSetup_;
    mutable scalarField terminalCurrentBuffer_;
    mutable scalarField terminalSourceBuffer_;

    static networkCouplingEndpoint& requireTerminalNetworkDomain
    (
        electroDomainInterface& secondaryDomain
    );

    static CouplingMode parseCouplingMode(const word& modeName);

    static word couplingModeName(CouplingMode mode);

    PVJCoupler
    (
        tissueCouplingEndpoint& primaryDomain,
        electroDomainInterface& secondaryDomain,
        const dictionary& dict
    );

    void clearTerminalCouplingBuffers() const;
};

} // End namespace Foam

#endif

// ************************************************************************* //
