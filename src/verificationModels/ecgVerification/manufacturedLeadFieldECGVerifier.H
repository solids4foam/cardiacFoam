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

#ifndef manufacturedLeadFieldECGVerifier_H
#define manufacturedLeadFieldECGVerifier_H
#include "ecgVerificationModel.H"
#include "overrideTypeName.H"
#include "OFstream.H"
namespace Foam
{
class manufacturedLeadFieldECGVerifier : public ecgVerificationModel
{
    scalarField scales_;
    scalar maxError_;
    label samples_;
    autoPtr<OFstream> output_;
public:
    OverrideTypeName("manufacturedLeadFieldECGVerifier");
    manufacturedLeadFieldECGVerifier(const electroStateProvider&, const word&,
        const dictionary&, const wordList&, const List<vector>&);
    virtual Requirements requirements() const override;
    virtual bool read(const dictionary&) override;
    virtual void record(const List<scalar>&) override;
    virtual void record(scalar, const List<scalar>&) override;
    virtual void end() override;
};
}
#endif
