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
    Foam::myocardiumDomainConfig

Description
    Bundles I/O and processing configuration into a single struct, addressing
    constructor parameter explosion in myocardium domain initialization.

Motivation
    The MyocardiumDomain constructor previously required
    12 separate parameters, making the interface difficult to understand and
    maintain. This struct groups 5 tightly-related I/O parameters:

    1. outFields           — Output field references (e.g., Vm, Iion)
    2. postProcessFieldNames — Names of fields to post-process
    3. postProcessFields   — Post-processed field values
    4. ionicModel          — The cellular model (AlievPanfilov, TNNP, etc.)
    5. prePostProcessor    — I/O handler for serialization

Design Pattern: Parameter Object
    Reduces method signature complexity by grouping related parameters.
    Makes the interface self-documenting (parameters are named, grouped
    logically). Simplifies future extensions (add new config fields without
    changing signature).

Scope
    This config object lives in electrophysicsSystemBuilder, which assembles
    and passes it to the myocardium domain during system construction.
    It is NOT exposed to ECG or other post-primary domains (one-way coupling).

Related
    - MyocardiumDomain uses this struct in its constructor
    - electrophysicsSystemBuilder creates and populates this config
    - ARCHITECTURE.md explains the overall domain separation design

\*---------------------------------------------------------------------------*/

#ifndef myocardiumDomainConfig_H
#define myocardiumDomainConfig_H

#include "PtrList.H"
#include "volFields.H"
#include "wordList.H"

namespace Foam
{

class ionicModel;
class electroVerificationModel;

struct myocardiumDomainConfig
{
    PtrList<volScalarField>& outFields;
    const wordList& postProcessFieldNames;
    PtrList<volScalarField>& postProcessFields;
    ionicModel& ionicModel;
    electroVerificationModel* verificationModelPtr;
};

} // End namespace Foam

#endif

// ************************************************************************* //
