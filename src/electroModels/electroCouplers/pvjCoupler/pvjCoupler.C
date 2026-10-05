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
#include "purkinjeModelIO.H"

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


scalarField pvjCoupler::readResistances(const dictionary& dict) const
{
    if (const scalarField* graphResistances =
            networkTerminalDomain_.terminalResistances())
    {
        return *graphResistances;
    }

    const scalar resistance(dict.get<scalar>("rPvj"));

    if (resistance <= 0)
    {
        FatalIOErrorInFunction(dict)
            << "rPvj must be positive [Ohm]; got " << resistance << "."
            << exit(FatalIOError);
    }

    return scalarField
    (
        networkTerminalDomain_.terminalNodes().size(),
        resistance
    );
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
        IOobject
        (
            IOobject::groupName
            (
                "lastObservedTissueActivation",
                dict.dictName()
            ),
            primaryDomain.mesh().time().timeName(),
            primaryDomain.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        0
    ),
    stabilityCheckedDeltaT_(-1.0)
{
    const label nTerminals = networkTerminalDomain_.terminalNodes().size();

    if (lastObservedTissueActivation_.filePath().empty())
    {
        lastObservedTissueActivation_.setSize(nTerminals, -1.0);
    }
    else if (lastObservedTissueActivation_.size() != nTerminals)
    {
        FatalErrorInFunction
            << lastObservedTissueActivation_.name() << " has "
            << lastObservedTissueActivation_.size() << " entries for "
            << nTerminals << " PVJs."
            << exit(FatalError);
    }
}


void pvjCoupler::write()
{
    if (couplingMode_ == bidirectional && mesh_.time().outputTime())
    {
        purkinjeModelIO::writeGlobalField(lastObservedTissueActivation_);
    }
}


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

    // The explicit term relaxes the sphere averages towards the network
    // voltages at rates a/dt, a = dt*lambda/(chi Cm), lambda the largest
    // eigenvalue of the junction operator: one junction's sum(w^2 V)/(R
    // V_s^2) when no spheres overlap, more where they share cells. With the
    // term on the old level, Euler is stable for a < 2 and BDF2 (backward)
    // for a < 4; above 2, BDF2 rings with a decaying step-to-step sign flip.
    const word VmName(tissue->Vm().name());
    ITstream& ddtIs = mesh_.ddtScheme("ddt(" + VmName + ")");
    const word ddtSchemeName(ddtIs);
    ddtIs.rewind();

    const bool bdf2 = (ddtSchemeName == "backward");
    const scalar aMax = bdf2 ? 4.0 : 2.0;
    const scalar chiCm = (tissue->chi()*tissue->Cm()).value();

    label peak = -1;
    const scalar rate =
        mapper_.sphereAverageDecayRate(1.0/resistance, peak)/chiCm;
    const scalar a = dt*rate;

    if (a >= aMax)
    {
        FatalErrorInFunction
            << "pvjCouplingScheme explicit is unstable: dt*lambda/(chi Cm) = "
            << a << ", bound " << aMax << " for ddt scheme '"
            << ddtSchemeName << "', with lambda the largest eigenvalue of "
            << "the junction term (overlapping spheres included); the "
            << "unstable mode is centred on PVJ " << peak << "." << nl
            << "  stable with every resistance scaled by more than "
            << a/aMax << " (R at PVJ " << peak << " = " << resistance[peak]
            << " Ohm)" << nl
            << "  or deltaT below " << aMax/rate << " s (now " << dt << ")"
            << nl
            << "Use pvjCouplingScheme implicit, larger resistances or "
            << "a smaller deltaT."
            << exit(FatalError);
    }

    if (bdf2 && a > 2.0)
    {
        WarningInFunction
            << "pvjCouplingScheme explicit: dt*lambda/(chi Cm) = " << a
            << " lies in (2, 4) (mode centred on PVJ " << peak << "), "
            << "where backward is stable but junction voltages oscillate "
            << "from step to step. pvjCouplingScheme implicit does not."
            << endl;
    }
}

} // End namespace Foam

// ************************************************************************* //
