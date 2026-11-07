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

Class
    electrophysicsSystemBuilder

Description
    Builder functions for assembling the electrophysics system components.

Author
    Simao Nieto de Castro. All rights reserved.
\*---------------------------------------------------------------------------*/

#ifndef electrophysicsSystemBuilder_H
#define electrophysicsSystemBuilder_H

#include "electrophysicsSystem.H"
#include "torsoCoupling.H"
#include "electroStateProvider.H"
#include "volFields.H"

namespace Foam
{

class ionicModel;
class electroVerificationModel;

namespace electrophysicsSystemBuilder
{

void configureMyocardiumDomain
(
    electrophysicsSystem&    system,
    const fvMesh&            mesh,
    const word&              solverType,
    const dictionary&        electroProperties,
    scalar                   initialDeltaT,
    const torsoCoupling::type coupling
);

void configureAdvanceScheme
(
    electrophysicsSystem& system,
    const dictionary&     electroProperties
);

void configureConductionDomains
(
    electrophysicsSystem& system,
    const fvMesh&         mesh,
    const dictionary&     electroProperties,
    scalar                initialDeltaT
);

//- The extracellularPotentialDomain of torsoCoupling bath, from the torso
//  block of the top-level electroProperties; none for any other coupling
void configureBathPotentialDomain
(
    electrophysicsSystem& system,
    const fvMesh&         mesh,
    const dictionary&     electroProperties,
    const torsoCoupling::type coupling
);

//- The ECG domains of the ecgDomains block of the top-level
//  electroProperties
void configureECGDomains
(
    electrophysicsSystem&       system,
    const myocardiumDomainInterface& myocardium,
    const electroStateProvider* potentialStateProviderPtr,
    const dictionary&           electroProperties
);

} // End namespace electrophysicsSystemBuilder
} // End namespace Foam

#endif

// ************************************************************************* //
