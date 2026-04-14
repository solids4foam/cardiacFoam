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

#include "eikonalPvjCoupler.H"

#include "addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug(EikonalPvjCoupler, 0);
addToRunTimeSelectionTable
(
    ElectroDomainCoupler,
    EikonalPvjCoupler,
    dictionary
);


void EikonalPvjCoupler::ensureSupportedMode() const
{
    // Mode checks removed to allow bidirectional eikonal coupling
}


EikonalPvjCoupler::EikonalPvjCoupler
(
    tissueCouplingEndpoint& primaryDomain,
    electroDomainInterface& secondaryDomain,
    const dictionary& dict
)
:
    PVJCoupler(primaryDomain, secondaryDomain, dict),
    terminalActivationBuffer_()
{
    Info << "Eikonal PVJ coupling model: pvjRadius=" << pvjRadius_
         << ", couplingMode=" << couplingModeName(couplingMode_)
         << endl;
}


void EikonalPvjCoupler::prepareSecondaryCoupling(scalar t0, scalar dt)
{
    (void)t0;
    (void)dt;

    if (couplingMode_ == bidirectional)
    {
        scalarField observedTissueTimes;
        mapper_.gatherActivationTimes
        (
            primaryDomain_.activationTime(),
            observedTissueTimes
        );

        networkTerminalDomain_.setTerminalActivationTime(observedTissueTimes);
    }

    clearTerminalCouplingBuffers();

    networkTerminalDomain_.setTerminalCoupling
    (
        terminalCurrentBuffer_,
        terminalSourceBuffer_
    );
}


void EikonalPvjCoupler::preparePrimaryCoupling(scalar t0, scalar dt)
{
    (void)t0;
    (void)dt;

    ensureSupportedMode();

    networkTerminalDomain_.terminalActivationTime(terminalActivationBuffer_);
    mapper_.depositActivationTimes
    (
        terminalActivationBuffer_,
        primaryDomain_.activationTimeField()
    );

    networkTerminalDomain_.setTerminalCoupling
    (
        terminalCurrentBuffer_,
        terminalSourceBuffer_
    );
}

} // End namespace Foam

// ************************************************************************* //
