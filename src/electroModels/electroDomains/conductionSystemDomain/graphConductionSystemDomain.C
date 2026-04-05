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

#include "graphConductionSystemDomain.H"

namespace Foam
{

GraphConductionSystemDomain::GraphConductionSystemDomain
(
    const ConductionSystemDomainContext& context,
    const dictionary& dict
)
:
    ConductionSystemDomain(context, dict),
    solverPtr_(GraphConductionSystemSolver::New(dict))
{}


void GraphConductionSystemDomain::reportAdvanceDiagnostics
(
    scalar t0,
    scalar dt
) const
{
    (void)t0;
    (void)dt;
}


void GraphConductionSystemDomain::advance
(
    scalar t0,
    scalar dt
)
{
    solverPtr_->advance(*this, t0, dt);
}

} // End namespace Foam

// ************************************************************************* //
