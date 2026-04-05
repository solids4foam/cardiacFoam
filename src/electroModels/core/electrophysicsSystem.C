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

#include "electrophysicsSystem.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::electrophysicsSystem::prepareUpstreamCouplings(scalar t0, scalar dt)
{
    forAll(upstreamCouplingModels_, i)
    {
        upstreamCouplingModels_[i].prepareSecondaryCoupling(t0, dt);
    }
}

void Foam::electrophysicsSystem::advanceUpstreamDomains(scalar t0, scalar dt)
{
    forAll(upstreamDomains_, i)
    {
        upstreamDomains_[i].advance(t0, dt);
    }
}

void Foam::electrophysicsSystem::preparePrimaryCouplings(scalar t0, scalar dt)
{
    forAll(upstreamCouplingModels_, i)
    {
        upstreamCouplingModels_[i].preparePrimaryCoupling(t0, dt);
    }
}

void Foam::electrophysicsSystem::prepareDownstreamCouplings
(
    scalar t0,
    scalar dt
)
{
    forAll(downstreamCouplingModels_, i)
    {
        downstreamCouplingModels_[i].preparePostPrimaryCoupling(t0, dt);
    }
}

void Foam::electrophysicsSystem::advanceDownstreamDomains(scalar t0, scalar dt)
{
    forAll(downstreamDomains_, i)
    {
        downstreamDomains_[i].advance(t0, dt);
    }
}

void Foam::electrophysicsSystem::writeUpstreamDomain()
{
    forAll(upstreamDomains_, i)
    {
        upstreamDomains_[i].write();
    }
}


void Foam::electrophysicsSystem::writeUpstreamCouplingModel()
{
    forAll(upstreamCouplingModels_, i)
    {
        upstreamCouplingModels_[i].write();
    }
}


void Foam::electrophysicsSystem::endUpstreamDomain()
{
    forAll(upstreamDomains_, i)
    {
        upstreamDomains_[i].end();
    }
}


void Foam::electrophysicsSystem::writeDownstreamDomain()
{
    forAll(downstreamDomains_, i)
    {
        downstreamDomains_[i].write();
    }
}


void Foam::electrophysicsSystem::endUpstreamCouplingModel()
{
    forAll(upstreamCouplingModels_, i)
    {
        upstreamCouplingModels_[i].end();
    }
}


void Foam::electrophysicsSystem::endDownstreamDomain()
{
    forAll(downstreamDomains_, i)
    {
        downstreamDomains_[i].end();
    }
}


void Foam::electrophysicsSystem::endDownstreamCouplingModel()
{
    forAll(downstreamCouplingModels_, i)
    {
        downstreamCouplingModels_[i].end();
    }
}


// ************************************************************************* //
