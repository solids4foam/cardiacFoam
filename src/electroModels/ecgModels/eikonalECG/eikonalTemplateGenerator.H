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
    Foam::eikonalECG_templates

Description
    Standalone generator for personalized eikonalECG action-potential
    templates. eikonalECG reconstructs a pseudo-ECG as
    Vm(x,t) = U(t - psi(x)), where U is one of several fixed, compiled-in
    templates (endo/mid/epi) stored in tissueTemplates.H, or one of N
    dynamically-generated ones. This component generates templates at
    runtime by pacing single cells of a real ionicModel with a case's own
    ionic parameterisation (heterogeneity constants and S1 stimulus
    protocol): 3 anchors (endo/mid/epi) for ionicHeterogeneity.mode
    transmuralBands, or N anchors (one per named region) for mode
    namedRegions.

    This file does not wire the generated templates into eikonalECG; it is
    an independent, separately testable component.

SourceFiles
    eikonalTemplateGenerator.C

Author
    Simao Nieto de Castro, UCD.
\*---------------------------------------------------------------------------*/

#ifndef eikonalTemplateGenerator_H
#define eikonalTemplateGenerator_H

#include "scalarField.H"
#include "dictionary.H"
#include "List.H"

namespace Foam
{
namespace eikonalECG_templates
{

//- Local time zero is the final S1 stimulus onset.
struct DynamicTemplate
{
    scalarField times;      // s, starts exactly at final S1 onset
    scalarField valuesMv;   // mV, direct ionicModel::signal(VM) output
};

//- nBeats overrides ionicModelConfig.singleCellStimulus.nstim1; neither
//  dict is mutated. Returns one template per anchor: 3 entries
//  ([endo, mid, epi], in that order) for mode transmuralBands (the
//  default), or one entry per region (sorted by rangeMin, matching
//  ionicHeterogeneity::parseNamedFieldRegions()'s own order) for mode
//  namedRegions. Fatals on any other heterogeneityDict.mode.
List<DynamicTemplate> generatePersonalizedTemplates
(
    const dictionary& ionicModelConfig,
    const dictionary& heterogeneityDict,
    label nBeats,
    scalar captureDuration,
    scalar dt
);

} // End namespace eikonalECG_templates
} // End namespace Foam

#endif

// ************************************************************************* //
