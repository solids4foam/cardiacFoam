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
    Manufactured voltage template U(s), s = t - psi, shared by eikonalECG and
    manufacturedEikonalECGVerifier: one action-potential-like event (default)
    or sin(2 pi s), selected by the verificationModel entry "template".

Author
    Simao Nieto de Castro, UCD.
\*---------------------------------------------------------------------------*/

#ifndef manufacturedTemplate_H
#define manufacturedTemplate_H

#include "dictionary.H"
#include "Enum.H"
#include "mathematicalConstants.H"

namespace Foam
{
namespace eikonalECG_templates
{

inline void readManufacturedTemplate
(
    const dictionary& verificationDict,
    bool& sine,
    scalar& duration,
    scalar& rise
)
{
    enum class templateKind
    {
        singleEvent,
        sine
    };

    // Unlike a plain word comparison, a misspelt value fails and lists the
    // valid names instead of silently selecting singleEvent.
    static const Enum<templateKind> templateKinds
    ({
        { templateKind::singleEvent, "singleEvent" },
        { templateKind::sine,        "sine"        }
    });

    sine =
        templateKinds.getOrDefault
        (
            "template",
            verificationDict,
            templateKind::singleEvent
        )
     == templateKind::sine;
    duration = verificationDict.lookupOrDefault<scalar>("templateDuration", 64);
    rise = verificationDict.lookupOrDefault<scalar>("templateRise", 4);
}


inline scalar evaluateManufacturedTemplateDerivative
(
    const scalar s,
    const bool sine,
    const scalar duration,
    const scalar rise
)
{
    if (sine)
    {
        return
            2.0*constant::mathematical::pi
           *Foam::cos(2.0*constant::mathematical::pi*s);
    }

    const scalar w = rise/(2.0*Foam::atanh(0.8));
    const scalar up = Foam::tanh(s/w);
    const scalar down = Foam::tanh((s - duration)/w);
    return
        0.25/w
       *((1.0 - sqr(up))*(1.0 - down) - (1.0 + up)*(1.0 - sqr(down)));
}

} // End namespace eikonalECG_templates
} // End namespace Foam

#endif
