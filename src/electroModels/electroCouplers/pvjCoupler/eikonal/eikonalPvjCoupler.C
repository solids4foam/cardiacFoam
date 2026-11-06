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

#include "eikonalPvjCoupler.H"
#include "dictionaryKeyCheck.H"

#include "addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug(eikonalPvjCoupler, 0);
addToRunTimeSelectionTable
(
    electroDomainCoupler,
    eikonalPvjCoupler,
    dictionary
);


word eikonalPvjCoupler::seedProfile(const dictionary& dict)
{
    return dictionaryKeyCheck::getWord
    (
        dict,
        "seedProfile",
        "front",
        wordList({"flat", "front", "curvedFront"})
    );
}


eikonalPvjCoupler::eikonalPvjCoupler
(
    tissueCouplingEndpoint& primaryDomain,
    electroDomainInterface& secondaryDomain,
    const dictionary& dict
)
:
    pvjCoupler(primaryDomain, secondaryDomain, dict),
    terminalActivationBuffer_(),
    seedProfile_(seedProfile(dict)),
    seedOffsetsSet_(false)
{
    if
    (
        seedProfile_ != "flat"
     && (
            !primaryDomain_.activationMetricPtr()
         || !primaryDomain_.frontCoefficientPtr()
        )
    )
    {
        FatalIOErrorInFunction(dict)
            << "seedProfile " << seedProfile_ << " needs an eikonal myocardium. "
            << "Set seedProfile flat for another tissue model."
            << exit(FatalIOError);
    }
}


void eikonalPvjCoupler::prepareSecondaryCoupling(scalar t0, scalar dt)
{
    (void)t0;
    (void)dt;

    if (couplingMode_ == bidirectional)
    {
        observeTerminalActivations();
    }

    clearTerminalCouplingBuffers();

    networkTerminalDomain_.setTerminalCoupling
    (
        terminalCurrentBuffer_,
        terminalSourceBuffer_
    );
}


void eikonalPvjCoupler::preparePrimaryCoupling(scalar t0, scalar dt)
{
    (void)t0;
    (void)dt;

    networkTerminalDomain_.terminalActivationTime(terminalActivationBuffer_);

    // Junction delay
    for (scalar& t : terminalActivationBuffer_)
    {
        if (t >= 0)
        {
            t += pvjDelay_;
        }
    }

    if (seedProfile_ != "flat" && !seedOffsetsSet_)
    {
        mapper_.setSeedOffsets
        (
            *primaryDomain_.activationMetricPtr(),
            primaryDomain_.frontCoefficientPtr()->value(),
            seedProfile_ == "curvedFront"
        );
        seedOffsetsSet_ = true;
    }

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
