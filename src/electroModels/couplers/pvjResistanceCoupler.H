/*---------------------------------------------------------------------------*\
License
    This file is part of cardiacFoam.

    cardiacFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    cardiacFoam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
    or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License along
    with cardiacFoam.  If not, see <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#ifndef PVJResistanceCoupler_H
#define PVJResistanceCoupler_H

#include "electroDomainCoupler.H"
#include "pvjMapper.H"
#include "overrideTypeName.H"

namespace Foam
{

class PVJResistanceCoupler
:
    public ElectroDomainCoupler
{
public:

    enum CouplingMode
    {
        networkToMyocardium,
        bidirectional
    };

private:

    networkCouplingEndpoint& networkTerminalDomain_;
    PVJMapper mapper_;
    scalar R_pvj_;
    scalar pvjRadius_;
    CouplingMode couplingMode_;
    mutable scalarField tissueVmBuffer_;
    mutable scalarField networkVmBuffer_;
    mutable scalarField terminalCurrentBuffer_;
    mutable scalarField terminalSourceBuffer_;

    static networkCouplingEndpoint& requireTerminalNetworkDomain
    (
        electrophysicsDomain& secondaryDomain
    );

    static CouplingMode parseCouplingMode(const word& modeName);

    static word couplingModeName(CouplingMode mode);

    void couplingCurrentAtPvjs
    (
        const scalarField& networkVm,
        const scalarField& tissueVm,
        scalarField& current
    ) const;

    void evaluateCoupling() const;

public:

    OverrideTypeName("PVJResistanceCoupler");

    PVJResistanceCoupler
    (
        tissueCouplingEndpoint& primaryDomain,
        electrophysicsDomain& secondaryDomain,
        const dictionary& dict
    );

    virtual void prepareSecondaryCoupling(scalar t0, scalar dt) override;

    virtual void preparePrimaryCoupling(scalar t0, scalar dt) override;
};

} // End namespace Foam

#endif

// ************************************************************************* //
