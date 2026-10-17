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
    Manufactured voltage template U(s), s = t - psi, of eikonalECG and
    manufacturedEikonalECGVerifier, selected by the verificationModel entry
    "template":

        singleEvent (default):
            U = Su(s) (1 - Su(s - templateDuration)), Su = (1 + tanh(s/w))/2,
            w such that Su rises from 10 % to 90 % within templateRise.
        sine:
            U = sin(2 pi s).

    The default duration exceeds the range of psi = exp(x + y + z), so
    U'(t - psi) keeps one sign over the domain; the sine's does not, and the
    ECG integral cancels.

Author
    Simao Nieto de Castro, UCD.
\*---------------------------------------------------------------------------*/

#ifndef manufacturedTemplate_H
#define manufacturedTemplate_H

#include "dictionary.H"
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
    sine =
        verificationDict.getOrDefault<word>("template", "singleEvent")
     == "sine";
    duration = verificationDict.getOrDefault<scalar>("templateDuration", 64);
    rise = verificationDict.getOrDefault<scalar>("templateRise", 4);
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
