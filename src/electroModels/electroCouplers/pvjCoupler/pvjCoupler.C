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

#include "pvjCoupler.H"
#include "electroDomainInterface.H"

namespace Foam
{

networkCouplingEndpoint& pvjCoupler::requireTerminalNetworkDomain
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


pvjCoupler::CouplingMode pvjCoupler::parseCouplingMode(const word& modeName)
{
    if (modeName == "unidirectional")
    {
        return unidirectional;
    }

    if (modeName == "bidirectional")
    {
        return bidirectional;
    }

    FatalErrorInFunction
        << "Unknown couplingMode '" << modeName << "'. "
        << "Valid options are 'unidirectional' and 'bidirectional'."
        << exit(FatalError);

    return unidirectional;
}


word pvjCoupler::couplingModeName(CouplingMode mode)
{
    return mode == bidirectional ? "bidirectional" : "unidirectional";
}


pvjCoupler::pvjCoupler
(
    tissueCouplingEndpoint& primaryDomain,
    electroDomainInterface& secondaryDomain,
    const dictionary& dict
)
:
    electroDomainCoupler(primaryDomain, secondaryDomain),
    networkTerminalDomain_(requireTerminalNetworkDomain(secondaryDomain)),
    mapper_
    (
        primaryDomain.mesh(),
        networkTerminalDomain_.terminalLocations(),
        dict.lookupOrDefault<scalar>("pvjRadius", 0.5e-3),
        dict.lookupOrDefault<word>("pvjKernel", "uniform")
    ),
    pvjRadius_(dict.lookupOrDefault<scalar>("pvjRadius", 0.5e-3)),
    couplingMode_(parseCouplingMode(dict.get<word>("couplingMode"))),
    terminalCurrentBuffer_
    (
        networkTerminalDomain_.terminalNodes().size(),
        0.0
    ),
    terminalSourceBuffer_
    (
        networkTerminalDomain_.terminalNodes().size(),
        0.0
    )
{}


void pvjCoupler::clearTerminalCouplingBuffers() const
{
    terminalCurrentBuffer_ = 0.0;
    terminalSourceBuffer_ = 0.0;
}

} // End namespace Foam

// ************************************************************************* //
