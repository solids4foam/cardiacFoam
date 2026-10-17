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

#ifndef eikonalTemplateSource_H
#define eikonalTemplateSource_H

#include "Switch.H"
#include "dictionary.H"
#include "volFields.H"
#include "eikonalTemplateGenerator.H"

namespace Foam
{
class ecgDomain;

// Shared template source for passive eikonal ECG observers. The legacy regional
// path differentiates activation time while holding local template weights fixed.
class eikonalTemplateSource
{
    Switch useManufacturedTemplate_;
    word templateType_;
    bool weightsCalculated_;
    List<scalarField> regionWeights_;
    bool personalizedTemplatesEnabled_;
    bool personalizedTemplatesGenerated_;
    dictionary personalizedTemplatesDict_;
    List<eikonalECG_templates::DynamicTemplate> personalizedTemplates_;

    scalar manufacturedTemplateDerivative(scalar localTime) const;
    const dictionary* findHeterogeneityDict(const dictionary&) const;
    void calculateTransmuralWeights(const ecgDomain&);
    void generatePersonalizedTemplates(const ecgDomain&);

public:
    explicit eikonalTemplateSource(const dictionary&);
    bool homogeneous() const { return templateType_ != "regional"; }
    const word& templateType() const { return templateType_; }
    void initialise(const ecgDomain&);
    // Chain rule for a uniform template; local-weight approximation for legacy
    // regional templates. No claim about a universal discrete convergence order.
    void reconstructGradVm
    (
        scalar sampleTime,
        const volScalarField& activationTime,
        const volVectorField& gradActivationTime,
        volVectorField& gradVm
    ) const;
};
}
#endif
