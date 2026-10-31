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

Namespace
    Foam::ionicTemplates

Description
    Single-cell templates of Vm and Cai, one per ionic-heterogeneity region,
    made by pacing single cells of the case's ionicModel. A cell activated at
    time tau follows its region templates at t - tau.

SourceFiles
    ionicTemplates.C

Author
    Simao Nieto de Castro, UCD.
\*---------------------------------------------------------------------------*/

#ifndef ionicTemplates_H
#define ionicTemplates_H

#include "scalarField.H"
#include "dictionary.H"
#include "List.H"
#include "fvMesh.H"

namespace Foam
{
namespace ionicTemplates
{

//- Time zero is the steepest upstroke of the final S1 beat.
struct DynamicTemplate
{
    scalarField times;      // s; 0 = steepest upstroke
    scalarField valuesMv;   // mV, direct ionicModel::signal(VM) output
    scalarField valuesCai;  // mM, direct ionicModel::signal(CAI) output
};

//- nBeats overrides ionicModelConfig.singleCellStimulus.nstim1; neither
//  dict is mutated. Returns one template per anchor: one entry per region
//  for mode namedRegions (sorted by rangeMin, matching
//  ionicHeterogeneity::parseNamedFieldRegions()'s own order) or mode
//  cellZoneRegions (in the order 'regions' declares them). Fatals on any
//  other heterogeneityDict.mode.
//- One template per region; regionStates seeds each region's cell from its single-cell state file
List<DynamicTemplate> generatePersonalizedTemplates
(
    const dictionary& ionicModelConfig,
    const dictionary& heterogeneityDict,
    label nBeats,
    scalar captureDuration,
    scalar dt,
    const dictionary* regionStates = nullptr,
    const fileName& stateRoot = fileName::null
);

//- Keys and controls of a personalizedTemplates dictionary
void checkControls(const dictionary& templatesDict);

//- nBeats, 1 by default when regionStates is given
label templateBeats(const dictionary& templatesDict);

//- Templates from a personalizedTemplates dictionary
List<DynamicTemplate> generate
(
    const dictionary& templatesDict,
    const dictionary& heterogeneityDict,
    const Time& runTime
);

//- ionicHeterogeneity of electroProperties or of its myocardium solver Coeffs
const dictionary* findHeterogeneityDict(const dictionary& electroProperties);

//- Per-cell weight of each region, index-aligned with the templates
List<scalarField> regionWeights
(
    const fvMesh& mesh,
    const dictionary& heterogeneityDict
);

//- Linear interpolation, clamped to the end values
inline scalar evaluateTemplate
(
    scalar localTime,
    const scalar* times,
    const scalar* values,
    int numSamples
)
{
    if (numSamples < 2) return 0.0;
    if (localTime <= times[0]) return values[0];
    if (localTime >= times[numSamples - 1]) return values[numSamples - 1];

    // O(1) Direct Index Lookup
    // Assumes times array has a perfectly uniform deltaT
    const scalar dt = times[1] - times[0];
    label lo = static_cast<label>((localTime - times[0]) / dt);

    // Safeguard against floating point inaccuracies
    if (lo < 0) lo = 0;
    if (lo >= numSamples - 1) lo = numSamples - 2;

    const scalar t0 = times[lo];
    const scalar t1 = times[lo + 1];

    // Fallback if the array isn't perfectly uniform at this specific index
    if (localTime < t0 || localTime > t1)
    {
        // Rare edge case: standard binary search fallback
        label hi = numSamples - 1;
        lo = 0;
        while (hi - lo > 1)
        {
            const label mid = (lo + hi)/2;
            if (times[mid] <= localTime) lo = mid;
            else hi = mid;
        }
    }

    const scalar exactT0 = times[lo];
    const scalar exactT1 = times[lo + 1];
    const scalar alpha = (localTime - exactT0)/(exactT1 - exactT0);

    return (1.0 - alpha)*values[lo] + alpha*values[lo + 1];
}

} // End namespace ionicTemplates
} // End namespace Foam

#endif

// ************************************************************************* //
