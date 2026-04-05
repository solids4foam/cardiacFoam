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

#include "pvjResistanceCoupler.H"

#include "addToRunTimeSelectionTable.H"
#include "electroDomainInterface.H"

namespace Foam
{

namespace
{

networkCouplingEndpoint& requireTerminalNetworkDomain
(
    electroDomainInterface& secondaryDomain
)
{
    networkCouplingEndpoint* terminalDomain =
        dynamic_cast<networkCouplingEndpoint*>(&secondaryDomain);

    if (!terminalDomain)
    {
        FatalErrorInFunction
            << "Configured network domain does not implement "
            << "networkCouplingEndpoint."
            << exit(FatalError);
    }

    return *terminalDomain;
}

} // End anonymous namespace


defineTypeNameAndDebug(PVJResistanceCoupler, 0);
addToRunTimeSelectionTable
(
    ElectroDomainCoupler,
    PVJResistanceCoupler,
    dictionary
);


networkCouplingEndpoint& PVJResistanceCoupler::requireTerminalNetworkDomain
(
    electroDomainInterface& secondaryDomain
)
{
    return Foam::requireTerminalNetworkDomain(secondaryDomain);
}


PVJResistanceCoupler::CouplingMode
PVJResistanceCoupler::parseCouplingMode(const word& modeName)
{
    if
    (
        modeName.empty()
     || modeName == "networkToMyocardium"
     || modeName == "oneWay"
     || modeName == "oneWayToMyocardium"
    )
    {
        return networkToMyocardium;
    }

    if (modeName == "bidirectional")
    {
        return bidirectional;
    }

    FatalErrorInFunction
        << "Unknown couplingMode '" << modeName << "'."
        << exit(FatalError);

    return networkToMyocardium;
}


word PVJResistanceCoupler::couplingModeName(CouplingMode mode)
{
    return mode == bidirectional ? "bidirectional" : "networkToMyocardium";
}


void PVJResistanceCoupler::couplingCurrentAtPvjs
(
    const scalarField& networkVm,
    const scalarField& tissueVm,
    scalarField& current
) const
{
    if (networkVm.size() != tissueVm.size())
    {
        FatalErrorInFunction
            << "PVJ coupling expected matching field sizes but received "
            << networkVm.size() << " and " << tissueVm.size()
            << exit(FatalError);
    }

    current.setSize(networkVm.size());
    forAll(current, i)
    {
        current[i] = (networkVm[i] - tissueVm[i])/R_pvj_;
    }
}


void PVJResistanceCoupler::evaluateCoupling
() const
{
    mapper_.gatherVm3DPvjs(primaryDomain_.Vm(), tissueVmBuffer_);
    networkTerminalDomain_.terminalVm(networkVmBuffer_);
    couplingCurrentAtPvjs
    (
        networkVmBuffer_,
        tissueVmBuffer_,
        terminalCurrentBuffer_
    );
    mapper_.volumetricSource(terminalCurrentBuffer_, terminalSourceBuffer_);
}


PVJResistanceCoupler::PVJResistanceCoupler
(
    tissueCouplingEndpoint& primaryDomain,
    electroDomainInterface& secondaryDomain,
    const dictionary& dict
)
:
    ElectroDomainCoupler(primaryDomain, secondaryDomain),
    networkTerminalDomain_(requireTerminalNetworkDomain(secondaryDomain)),
    mapper_
    (
        primaryDomain.mesh(),
        networkTerminalDomain_.terminalLocations(),
        dict.lookupOrDefault<scalar>("pvjRadius", 0.5e-3)
    ),
    R_pvj_(dict.get<scalar>("R_pvj")),
    pvjRadius_(dict.lookupOrDefault<scalar>("pvjRadius", 0.5e-3)),
    couplingMode_
    (
        parseCouplingMode
        (
            dict.lookupOrDefault<word>("couplingMode", word::null)
        )
    ),
    tissueVmBuffer_(),
    networkVmBuffer_(),
    terminalCurrentBuffer_(),
    terminalSourceBuffer_()
{
    Info << "PVJ coupling model: R_pvj=" << R_pvj_
         << ", pvjRadius=" << pvjRadius_
         << ", couplingMode=" << couplingModeName(couplingMode_)
         << endl;
}


void PVJResistanceCoupler::prepareSecondaryCoupling(scalar t0, scalar dt)
{
    (void)t0;
    (void)dt;

    evaluateCoupling();

    if (couplingMode_ == networkToMyocardium)
    {
        terminalCurrentBuffer_ = 0.0;
        terminalSourceBuffer_ = 0.0;
    }

    networkTerminalDomain_.setTerminalCoupling
    (
        terminalCurrentBuffer_,
        terminalSourceBuffer_
    );
}


void PVJResistanceCoupler::preparePrimaryCoupling(scalar t0, scalar dt)
{
    (void)t0;
    (void)dt;

    evaluateCoupling();

    mapper_.depositCoupling
    (
        terminalCurrentBuffer_,
        primaryDomain_.sourceField()
    );

    // Retain the latest physical exchange for diagnostics/output even when the
    // network-side advance is configured one-way. The next pre-network hook
    // will zero these terms again before evolve() when needed.
    networkTerminalDomain_.setTerminalCoupling
    (
        terminalCurrentBuffer_,
        terminalSourceBuffer_
    );
}

} // End namespace Foam

// ************************************************************************* //
