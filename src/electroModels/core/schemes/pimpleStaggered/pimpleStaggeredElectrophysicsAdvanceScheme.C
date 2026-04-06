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

    You should have received a copy of the GNU General Public License
    along with cardiacFoam.  If not, see <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "pimpleStaggeredElectrophysicsAdvanceScheme.H"
#include "electrophysicsSystem.H"
#include "addToRunTimeSelectionTable.H"
#include "clockTime.H"

namespace Foam
{

defineTypeNameAndDebug(pimpleStaggeredElectrophysicsAdvanceScheme, 0);
addToRunTimeSelectionTable
(
    electrophysicsAdvanceScheme,
    pimpleStaggeredElectrophysicsAdvanceScheme,
    dictionary
);

bool pimpleStaggeredElectrophysicsAdvanceScheme::advance
(
    scalar t0,
    scalar dt,
    electrophysicsSystem& system,
    pimpleControl* pimplePtr,
    electrophysicsAdvanceTimings& timings
)
{
    if (!pimplePtr)
    {
        FatalErrorInFunction
            << "pimpleStaggeredElectrophysicsAdvanceScheme requires "
            << "solutionAlgorithm = implicit and a valid pimpleControl. "
            << "Set solutionAlgorithm implicit; in the electroModel Coeffs dict."
            << exit(FatalError);
    }

    timings = electrophysicsAdvanceTimings();

    MyocardiumDomain& myocardium = system.myocardium();

    clockTime timer;

    myocardium.prepareTimeStep(t0, dt);

    while (pimplePtr->loop())
    {
        system.prepareUpstreamCouplings(t0, dt);
        timings.couplingTime += timer.timeIncrement();

        system.advanceUpstreamDomains(t0, dt);
        timings.upstreamDomainTime += timer.timeIncrement();

        system.preparePrimaryCouplings(t0, dt);
        timings.couplingTime += timer.timeIncrement();

        myocardium.advance(t0, dt, pimplePtr);
        timings.primaryDomainTime += timer.timeIncrement();
    }

    system.prepareDownstreamCouplings(t0, dt);
    timings.downstreamCouplingTime = timer.timeIncrement();

    system.advanceDownstreamDomains(t0, dt);
    timings.downstreamDomainTime = timer.timeIncrement();

    return true;
}

} // End namespace Foam

// ************************************************************************* //
