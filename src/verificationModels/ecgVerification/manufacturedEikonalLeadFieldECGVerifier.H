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

#ifndef manufacturedEikonalLeadFieldECGVerifier_H
#define manufacturedEikonalLeadFieldECGVerifier_H

#include "ecgVerificationModel.H"
#include "overrideTypeName.H"
#include "OFstream.H"

namespace Foam
{
// Independent boundary-integral reference on the 3D unit cube for
// tau=exp(x+y+z), U(s)=sin(2*pi*s), L_i=scale_i*(1,2,3) m^-2.
class manufacturedEikonalLeadFieldECGVerifier : public ecgVerificationModel
{
    scalarField scales_;
    List<scalar> nodes_, weights_, checkNodes_, checkWeights_;
    scalar maxError_, maxDelta_;
    label samples_;
    autoPtr<OFstream> output_;
    scalar reference(scalar,const UList<scalar>&,const UList<scalar>&) const;
public:
    OverrideTypeName("manufacturedEikonalLeadFieldECGVerifier");
    manufacturedEikonalLeadFieldECGVerifier(const electroStateProvider&,const word&,
        const dictionary&,const wordList&,const List<vector>&);
    virtual Requirements requirements() const override;
    virtual bool read(const dictionary&) override;
    virtual void record(const List<scalar>&) override;
    virtual void record(scalar,const List<scalar>&) override;
    virtual void end() override;
};
}
#endif
