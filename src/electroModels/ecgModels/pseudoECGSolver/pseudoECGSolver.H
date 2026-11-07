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
#include "infiniteMediumLeadField.H"
#include "OFstream.H"
#include "autoPtr.H"
#include "overrideTypeName.H"

namespace Foam
{

class ecgDomain;

class pseudoECGSolver
:
    public ecgSolver
{
    infiniteMediumLeadField leadField_;
    bool reportedConductivitySource_;

    //- True when a 'sampling' sub-dict is present: the solver then owns its
    //  own time-series file, written at deltaT_ independently of writeInterval.
    bool hasOwnSampling_;
    scalar startTime_;
    scalar endTime_;
    scalar deltaT_;
    scalar nextSampleTime_;
    autoPtr<OFstream> outputPtr_;

public:

    OverrideTypeName("pseudoECG");

    explicit pseudoECGSolver(const dictionary& dict);

    virtual bool writesOwnTimeSeries() const override
    {
        return hasOwnSampling_;
    }

    //- Reads Vm and the conductivity
    virtual void checkState(const ecgDomain& domain) const override;

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
