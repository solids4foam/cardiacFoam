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

#ifndef eikonalLeadFieldECG_H
#define eikonalLeadFieldECG_H

#include "ecgSolver.H"
#include "eikonalTemplateSource.H"
#include "leadFieldData.H"
#include "OFstream.H"

namespace Foam
{
// Passive reciprocal recording from activation time and one spatially uniform
// voltage template. Owns the sampling loop, independently of solver time steps.
class eikonalLeadFieldECG : public ecgSolver
{
    eikonalTemplateSource templates_;
    scalar start_, end_, delta_;
    bool heartOnly_, sampled_;
    autoPtr<leadFieldData> data_;
    autoPtr<OFstream> output_;
    scalarField lastValues_;
public:
    TypeName("eikonalLeadFieldECG");
    explicit eikonalLeadFieldECG(const dictionary&);
    virtual void initialise(const ecgDomain&) override;
    virtual wordList channelNames(const ecgDomain&) const override;
    virtual List<vector> channelPositions(const ecgDomain&) const override;
    virtual bool writesOwnTimeSeries() const override { return true; }
    virtual void solve(ecgDomain&, scalar, scalar, scalarField&) override;
};
}
#endif
