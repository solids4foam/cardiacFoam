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

#ifndef electrophysicsSystemBuilder_H
#define electrophysicsSystemBuilder_H

#include "electrophysicsSystem.H"
#include "electroStateProvider.H"
#include "volFields.H"

namespace Foam
{

class ionicModel;
class electroVerificationModel;

namespace electrophysicsSystemBuilder
{

// Terminology used by the orchestration layer:
//   primary domain        = myocardium
//   conduction domains    = pre-myocardium graph / Purkinje domains
//   ECG domains           = post-myocardium ECG domains

void configureMyocardiumDomain
(
    electrophysicsSystem&    system,
    const fvMesh&            mesh,
    const dictionary&        electroProperties,
    PtrList<volScalarField>& outFields,
    const wordList&          postProcessFieldNames,
    PtrList<volScalarField>& postProcessFields,
    autoPtr<ionicModel>&     ionicModelPtr,
    autoPtr<electroVerificationModel>& verificationModelPtr,
    scalar                   initialDeltaT
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

void configureECGDomains
(
    electrophysicsSystem&       system,
    const electroStateProvider& stateProvider,
    const dictionary&           electroProperties
);

} // End namespace electrophysicsSystemBuilder
} // End namespace Foam

#endif

// ************************************************************************* //
