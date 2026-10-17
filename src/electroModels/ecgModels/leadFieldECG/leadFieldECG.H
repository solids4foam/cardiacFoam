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

#ifndef leadFieldECG_H
#define leadFieldECG_H
#include "ecgSolver.H"
#include "leadFieldData.H"
#include "overrideTypeName.H"

namespace Foam
{
class leadFieldECG : public ecgSolver
{
    autoPtr<leadFieldData> data_;
    bool heartOnly_;
public:
    OverrideTypeName("leadFieldECG");
    explicit leadFieldECG(const dictionary&);
    virtual void initialise(const ecgDomain&) override;
    virtual wordList channelNames(const ecgDomain&) const override;
    virtual List<vector> channelPositions(const ecgDomain&) const override;
    virtual void solve(ecgDomain&, scalar, scalar, scalarField&) override;
};
}
#endif
