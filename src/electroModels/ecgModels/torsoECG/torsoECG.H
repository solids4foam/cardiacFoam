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

Class
    torsoECG

Description
    Torso ECG.

Author
    Simao Nieto de Castro. All rights reserved.
\*---------------------------------------------------------------------------*/

#ifndef torsoECG_H
#define torsoECG_H

#include "ecgSolver.H"
#include "labelList.H"
#include "overrideTypeName.H"
#include "Switch.H"

namespace Foam
{

class torsoECG
:
    public ecgSolver
{
    labelList electrodeCells_;
    bool cellsBuilt_;

    void buildElectrodeCells(const ecgDomain& domain);

public:

    OverrideTypeName("torsoECG");

    explicit torsoECG(const dictionary& dict);

    virtual void solve
    (
        ecgDomain& domain,
        scalar t0,
        scalar dt,
        scalarField& values
    ) override;
};

} // End namespace Foam

#endif

// ************************************************************************* //
