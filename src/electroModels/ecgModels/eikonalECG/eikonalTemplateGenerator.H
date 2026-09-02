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
    Vm(x,t) = U(t - psi(x)), where U is one of three fixed, compiled-in
    templates (endo/mid/epi) stored in tissueTemplates.H. This component
    instead generates those three templates at runtime by pacing three
    single cells of a real ionicModel with a case's own ionic
    parameterisation (heterogeneity constants and S1 stimulus protocol).

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

namespace Foam
{
namespace eikonalECG_templates
{

//- A single generated action-potential trace, sampled at a fixed dt
//  starting at the final S1 stimulus onset (local time zero).
struct DynamicTemplate
{
    scalarField times;      // s, starts exactly at final S1 onset
    scalarField valuesMv;   // mV, direct ionicModel::signal(VM) output
};

//- The three transmural anchor templates eikonalECG selects between.
struct TemplateTriplet
{
    DynamicTemplate endo;
    DynamicTemplate mid;
    DynamicTemplate epi;
};

//- Generate endo/mid/epi action-potential templates by pacing three
//  single cells (transmural coordinates 0.0, 0.5*(endoMInterface +
//  mEpiInterface), 1.0) of the ionic model named in ionicModelConfig.
//
//  ionicModelConfig is the active solver's own ionicModel coefficients
//  dictionary (must contain "ionicModel" and "singleCellStimulus").
//  heterogeneityDict is the active solver's transmuralBands
//  ionicHeterogeneity dictionary (endoMInterface, mEpiInterface,
//  transitionWidth, smoothing, transitionMode). Neither is mutated.
//
//  nBeats S1 pulses are paced (singleCellStimulus's own nstim1 is
//  overridden locally); captureDuration/dt of Vm are then recorded
//  starting at the final S1 onset. Every returned trace is validated for
//  finiteness, strictly increasing times, minimum amplitude, and
//  end-of-capture recovery to rest before this function returns.
TemplateTriplet generatePersonalizedTemplates
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
