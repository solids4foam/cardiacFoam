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
#include "electroVolumeFieldDomain.H"

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


word pvjCoupler::readCouplingScheme(const dictionary& dict)
{
    const word scheme
    (
        dict.lookupOrDefault<word>("pvjCouplingScheme", "explicit")
    );

    if (scheme != "explicit" && scheme != "implicit")
    {
        FatalIOErrorInFunction(dict)
            << "Unknown pvjCouplingScheme '" << scheme
            << "'. Valid options are 'explicit' and 'implicit'."
            << exit(FatalIOError);
    }

    return scheme;
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
    ),
    lastObservedTissueActivation_
    (
        networkTerminalDomain_.terminalNodes().size(),
        -1.0
    ),
    stabilityCheckedDeltaT_(-1.0)
{}


void pvjCoupler::clearTerminalCouplingBuffers() const
{
    terminalCurrentBuffer_ = 0.0;
    terminalSourceBuffer_ = 0.0;
}


void pvjCoupler::observeTerminalActivations()
{
    scalarField observedTissueTimes;
    scalarField latestTissueTimes;
    mapper_.gatherActivationTimes
    (
        primaryDomain_.activationTime(),
        lastObservedTissueActivation_,
        observedTissueTimes,
        latestTissueTimes
    );

    lastObservedTissueActivation_ = latestTissueTimes;
    networkTerminalDomain_.setTerminalActivationObservations(observedTissueTimes);
}



void pvjCoupler::checkExplicitCouplingStability
(
    const scalarField& resistance,
    const scalar dt
) const
{
    if (dt == stabilityCheckedDeltaT_)
    {
        return;
    }
    stabilityCheckedDeltaT_ = dt;

    const electroVolumeFieldDomain* tissue =
        dynamic_cast<const electroVolumeFieldDomain*>(&primaryDomain_);

    if (!tissue)
    {
        FatalErrorInFunction
            << "pvjCouplingScheme explicit needs a tissue domain with chi "
            << "and cm to bound its time step."
            << exit(FatalError);
    }

    // The explicit term relaxes the sphere average towards the network
    // voltage at rate a/dt, a = dt*sum(w^2 V)/(R chi Cm V_s^2). With the
    // term on the old level, Euler is stable for a < 2 and BDF2 (backward)
    // for a < 4; above 2, BDF2 rings with a decaying step-to-step sign flip.
    const word VmName(tissue->Vm().name());
    ITstream& ddtIs = mesh_.ddtScheme("ddt(" + VmName + ")");
    const word ddtSchemeName(ddtIs);
    ddtIs.rewind();

    const bool bdf2 = (ddtSchemeName == "backward");
    const scalar aMax = bdf2 ? 4.0 : 2.0;
    const scalar chiCm = (tissue->chi()*tissue->Cm()).value();
    const scalarField rates(mapper_.sphereAverageDecayRates());

    label nRinging = 0;
    scalar aRinging = 0.0;

    forAll(rates, i)
    {
        const scalar a = dt*rates[i]/(resistance[i]*chiCm);

        if (a >= aMax)
        {
            FatalErrorInFunction
                << "pvjCouplingScheme explicit is unstable at PVJ " << i
                << ": dt*sum(w^2 V)/(R chi Cm V_s^2) = " << a
                << ", bound " << aMax << " for ddt scheme '"
                << ddtSchemeName << "'." << nl
                << "  R      = " << resistance[i] << " Ohm (stable above "
                << dt*rates[i]/(aMax*chiCm) << " Ohm)" << nl
                << "  deltaT = " << dt << " s (stable below "
                << aMax*resistance[i]*chiCm/rates[i] << " s)" << nl
                << "Use pvjCouplingScheme implicit, a larger resistance or "
                << "a smaller deltaT."
                << exit(FatalError);
        }

        if (bdf2 && a > 2.0)
        {
            ++nRinging;
            aRinging = max(aRinging, a);
        }
    }

    if (nRinging)
    {
        WarningInFunction
            << "pvjCouplingScheme explicit: " << nRinging << " PVJ(s) have "
            << "dt*sum(w^2 V)/(R chi Cm V_s^2) in (2, 4) (largest "
            << aRinging << "), where backward is stable but the junction "
            << "voltage oscillates from step to step. pvjCouplingScheme "
            << "implicit does not." << endl;
    }
}

} // End namespace Foam

// ************************************************************************* //
