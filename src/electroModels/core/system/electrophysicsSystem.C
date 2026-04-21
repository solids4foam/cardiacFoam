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

void Foam::electrophysicsSystem::prepareConductionCouplings(scalar t0, scalar dt)
{
    forAll(conductionCouplingModels_, i)
    {
        conductionCouplingModels_[i].prepareSecondaryCoupling(t0, dt);
    }
}


void Foam::electrophysicsSystem::advanceConductionDomains(scalar t0, scalar dt)
{
    forAll(conductionDomains_, i)
    {
        conductionDomains_[i].advance(t0, dt);
    }
}


void Foam::electrophysicsSystem::prepareMyocardiumCouplings(scalar t0, scalar dt)
{
    forAll(conductionCouplingModels_, i)
    {
        conductionCouplingModels_[i].preparePrimaryCoupling(t0, dt);
    }
}


void Foam::electrophysicsSystem::prepareECGCouplings
(
    scalar t0,
    scalar dt
)
{
    forAll(ecgCouplingModels_, i)
    {
        ecgCouplingModels_[i].preparePostPrimaryCoupling(t0, dt);
    }
}


void Foam::electrophysicsSystem::advanceECGDomains
(
    scalar t0,
    scalar dt
)
{
    forAll(ecgDomains_, i)
    {
        ecgDomains_[i].advance(t0, dt);
    }
}


void Foam::electrophysicsSystem::writeConductionDomains()
{
    forAll(conductionDomains_, i)
    {
        conductionDomains_[i].write();
    }
}


void Foam::electrophysicsSystem::writeConductionCouplings()
{
    forAll(conductionCouplingModels_, i)
    {
        conductionCouplingModels_[i].write();
    }
}


void Foam::electrophysicsSystem::endConductionDomains()
{
    forAll(conductionDomains_, i)
    {
        conductionDomains_[i].end();
    }
}


void Foam::electrophysicsSystem::writeECGDomains()
{
    forAll(ecgDomains_, i)
    {
        ecgDomains_[i].write();
    }
}


void Foam::electrophysicsSystem::endConductionCouplings()
{
    forAll(conductionCouplingModels_, i)
    {
        conductionCouplingModels_[i].end();
    }
}


void Foam::electrophysicsSystem::endECGDomains()
{
    forAll(ecgDomains_, i)
    {
        ecgDomains_[i].end();
    }
}


void Foam::electrophysicsSystem::endECGCouplings()
{
    forAll(ecgCouplingModels_, i)
    {
        ecgCouplingModels_[i].end();
    }
}


// ************************************************************************* //
