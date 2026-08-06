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

Description
    Spatial ionic heterogeneity configuration functions for ionic models.
    Dispatches transmuralBands, namedRegions, cellZoneRegions, and
    apex-to-base setup.

Author
    Simao Nieto de Castro, UCD.
\*---------------------------------------------------------------------------*/

#ifndef ionicHeterogeneityOrchestrator_H
#define ionicHeterogeneityOrchestrator_H

#include "ionicHeterogeneity.H"
#include "PtrList.H"
#include "scalarField.H"
#include "dictionary.H"
#include "word.H"

namespace Foam
{

class ionicModel;

namespace ionicHeterogeneityOrchestrator
{

//- Configure region-based heterogeneity from the selected mode.
void configureTransmuralBandHeterogeneity
(
    const ionicModel& model,
    const scalarField& transmuralDistance,
    const dictionary& heterogeneityDict,
    PtrList<scalarField>& heterogeneousConstants,
    PtrList<scalarField>* heterogeneousInitialStates = nullptr
);

//- Apply apex-to-base exponential scaling to selected constants.
void configureApexBaseBandsHeterogeneityImpl
(
    const ionicModel& model,
    const scalarField& apexDist,
    const dictionary& dict,
    PtrList<scalarField>& heterogeneousConstants
);

//- Configure heterogeneity from named field regions.
void configureNamedRegionHeterogeneity
(
    const ionicModel& model,
    const scalarField& fieldValues,
    const dictionary& heterogeneityDict,
    PtrList<scalarField>& heterogeneousConstants,
    PtrList<scalarField>* heterogeneousInitialStates
);

//- Configure heterogeneity from mesh cell zones.
void configureCellZoneRegionHeterogeneity
(
    const ionicModel& model,
    const scalarField& regionIndices,
    const dictionary& heterogeneityDict,
    PtrList<scalarField>& heterogeneousConstants,
    PtrList<scalarField>* heterogeneousInitialStates
);

//- Blend region constants and initial states across the supplied regions.
void blendNamedRegions
(
    const ionicModel& model,
    const scalarField& fieldValues,
    const List<ionicHeterogeneity::NamedFieldRegion>& regions,
    scalar transitionWidth,
    const word& smoothing,
    const word& transitionMode,
    PtrList<scalarField>& heterogeneousConstants,
    PtrList<scalarField>* heterogeneousInitialStates
);

//- Return constants for a region and its baseline tissue.
scalarField constantsForRegion
(
    const ionicModel& model,
    const word& regionName,
    const word& baseline,
    const wordList& knownRegionNames
);

//- Return initial states for a region baseline.
scalarField initialStatesForRegion
(
    const ionicModel& model,
    const word& regionName,
    const word& baseline,
    const wordList& knownRegionNames
);

} // End namespace ionicHeterogeneityOrchestrator
} // End namespace Foam

#endif

// ************************************************************************* //
