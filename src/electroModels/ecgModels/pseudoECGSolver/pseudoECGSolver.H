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
    pseudoECGSolver

Description
    Pseudo-ECG solver.

Author
    Simao Nieto de Castro. All rights reserved.
\*---------------------------------------------------------------------------*/

#ifndef pseudoECGSolver_H
#define pseudoECGSolver_H

#include "ecgSolver.H"
#include "overrideTypeName.H"

namespace Foam
{

class ecgDomain;

class pseudoECGSolver
:
    public ecgSolver
{
    scalar sigmaE_;
    bool reportedConductivitySource_;
    bool leadVectorsCalculated_;
    List<List<vector>> leadVectors_;

    void calculateLeadVectors(const ecgDomain& domain);

public:

    OverrideTypeName("pseudoECG");

    explicit pseudoECGSolver(const dictionary& dict)
    :
        sigmaE_(dict.lookupOrDefault<scalar>("sigmaExtracellular", 0.0)),
        reportedConductivitySource_(false),
        leadVectorsCalculated_(false)
    {}

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
